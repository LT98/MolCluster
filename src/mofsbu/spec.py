"""A build is DATA, not a function call.

Ground rule 8's precondition: you cannot queue keyword arguments.  A run is described by
a serialisable spec, the spec is stored, and a worker somewhere executes it — in this
process on a laptop, across cores on the workstation, and later under MPI or on a GPU
without any of the callers changing.  Every interface (notebook, CLI, a future form) is
an editor of this object and nothing more.

The spec is also the reproducibility record: hand someone the JSON and they get your run.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

SPEC_VERSION = 8

# What to do with each structure once it is constructed.  Whether a mode can RUN is a
# property of the machine, not of the spec: `energy.relax.mode_status()` asks the
# backends, and the planner refuses with that reason.  `dft_go` is still a settled name
# with no body (ground rule 7) — there is no external code wired up.
RUN_MODES = ("construct", "ml_go", "xtb_go", "dft_go")

# `ml_go` names a RUNG, not a theory.  Two MACE foundation models serve it and they are
# not interchangeable: MACE-MP-0 (Materials Project) cannot see formal charge or spin,
# MACE-OMOL-0 (OMol25, wB97M-V/def2-TZVPD) takes both as inputs.  A spec may name one;
# `None` means "whatever this machine declares in MOFSBU_ML_MODEL", and either way the
# resolved model is written into the `methods` row of every number, so a stored energy
# never loses the name of the model that produced it.
ML_MODELS = ("mace-mp-0", "mace-omol-0")

# How many co-ligands (water or any named solvent/co-ligand) a product carries (D26).
# `fill` puts one on every vertex the ligands leave free; `range` also builds the counts up
# to `co_ligand_window` below that, the rest of the polyhedron left empty, and with
# `pathways` makes gaining one co-ligand a ladder step.
CO_LIGAND_COUNTS = ("fill", "range")


#: Separators a range may be written with.  `~` is the one the builder page offers;
#: `-` and `..` are accepted because people type them and both are unambiguous here —
#: every quantity this parser reads (a coordination number, a count of ligand copies) is
#: a positive integer, so a leading `-` can only ever be a range separator.
_RANGE_SEPARATORS = ("~", "..", "-")

#: How many values one range may expand to.  A guard rail, not a policy: `1~3` is the
#: point of the notation and `1~500` is a typo that would otherwise fill a queue before
#: anyone noticed.  It refuses rather than truncating, so nothing is silently dropped.
MAX_RANGE_SPAN = 64


def int_series(value: Any, *, what: str = "value") -> tuple[int, ...]:
    """Expand a list of positive integers written as numbers, ranges, or both.

    `[4, 6]`, `"4,6"`, `"1~3"` and `"1~3, 6"` are all read; the last two are why this
    exists — a run that sweeps one to three ligand copies is one experiment, and writing
    it as a range is how it gets asked for rather than typed out.

    The EXPANDED tuple is what a spec stores.  A spec is the reproducibility record, so it
    carries the integers a run actually enumerated; keeping `"1~3"` in the file would make
    the digest depend on how the request was phrased and would leave every reader of the
    JSON re-implementing this parser.

    Duplicates collapse (first occurrence wins) and order is otherwise preserved: the
    order decides which combination is queued first, and re-sorting a hand-written spec
    would quietly reorder its run.
    """
    if isinstance(value, (str, int)) and not isinstance(value, bool):
        items: list[Any] = [value]
    else:
        try:
            items = list(value)
        except TypeError:
            raise ValueError(f"{what}: expected numbers or a range, got {value!r}") from None

    out: list[int] = []
    for item in items:
        for n in _expand_token(item, what):
            if n not in out:
                out.append(n)
    return tuple(out)


def _expand_token(item: Any, what: str) -> list[int]:
    if isinstance(item, bool):                       # bool is an int; it is not a count
        raise ValueError(f"{what}: {item!r} is not a number")
    if isinstance(item, int):
        return [_positive(item, what)]
    if not isinstance(item, str):
        raise ValueError(f"{what}: expected a number or a range, got {item!r}")
    out: list[int] = []
    for token in item.replace(",", " ").split():
        lo_hi = _split_range(token, what)
        if lo_hi is None:
            out.append(_positive(_int(token, what), what))
            continue
        lo, hi = lo_hi
        if hi < lo:
            raise ValueError(
                f"{what}: range {token!r} counts down ({lo} to {hi}); write it low to high")
        if hi - lo + 1 > MAX_RANGE_SPAN:
            raise ValueError(
                f"{what}: range {token!r} is {hi - lo + 1} values, past the {MAX_RANGE_SPAN} "
                f"a single range may expand to. List the ones you mean.")
        out.extend(range(lo, hi + 1))
    return out


def _split_range(token: str, what: str) -> tuple[int, int] | None:
    for sep in _RANGE_SEPARATORS:
        head, found, tail = token.partition(sep)
        if found and head and tail:
            return _positive(_int(head, what), what), _positive(_int(tail, what), what)
    return None


def _int(text: str, what: str) -> int:
    try:
        return int(text.strip())
    except ValueError:
        raise ValueError(f"{what}: {text.strip()!r} is not a whole number") from None


def _positive(n: int, what: str) -> int:
    if n < 1:
        raise ValueError(f"{what}: {n} is not a count; these start at 1")
    return n


@dataclass(frozen=True)
class MoleculeSpec:
    name: str
    smiles: str
    multiplicity: int = 1
    max_deprotonations: int | None = None


@dataclass(frozen=True)
class MetalSpec:
    """A coordination centre.  No `multiplicity` field on purpose.

    There used to be one, defaulting to 1 and unconnected to `oxidation_state` or
    `spin_class` — nothing kept it in sync, so a Cu(II) `MetalSpec` silently carried a
    singlet (d9 has one unpaired electron; no singlet is reachable) until MACE-OMOL-0
    started checking.  The complex's multiplicity is now DERIVED from `spin_class` at
    build time (`runner.execute`, via `energy.backends.spin_class_multiplicity`), the
    same way `high_spin_multiplicity` already was for a bare ion — stated once, as a
    spin class, not copied around as a number that can go stale.
    """

    symbol: str
    oxidation_state: int = 2
    spin_class: str = "ls"


@dataclass(frozen=True)
class PocketPredicate:
    """Structural selection of chelate pockets.  Never a chemical name."""

    ring_size: int | None = None
    n_anionic: int | None = None
    min_anionic: int | None = None
    rigid: bool | None = None

    def as_kwargs(self) -> dict[str, Any]:
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass(frozen=True)
class BuildSpec:
    """Everything a run needs.  Serialise it, store it, execute it elsewhere."""

    molecules: tuple[MoleculeSpec, ...]
    metals: tuple[MetalSpec, ...] = ()
    degree: int = 1
    coordination: tuple[int, ...] = (4, 6)
    geometries: tuple[str, ...] | None = None
    ligands_per_metal: tuple[int, ...] = (1,)
    co_ligand: str | None = "O"
    binding: tuple[str, ...] = ("chelate", "mono")
    pocket: PocketPredicate = field(default_factory=PocketPredicate)
    seed: int = 0
    # With no co-ligand, a bidentate ligand cannot reach CN 4 on its own.  Setting this
    # builds the coordinatively unsaturated product instead of skipping the combination —
    # off by default, because quietly building something smaller than you asked for is
    # exactly the failure this project keeps running into.
    allow_unsaturated: bool = False
    # How many DIFFERENT ligand kinds may share one coordination sphere.
    #
    # 1 is homoleptic — n copies of one thing, plus the co-ligand — and it is the default
    # because it is what every spec written before this field meant.  Raising it is how
    # you ask for [Mg(dtBK)(Cl)]: a "kind" is one (molecule, protomer, donor set, binding
    # mode), so two protomers of one molecule are two kinds and a partially-deprotonated
    # set on one centre is expressible.
    #
    # It is a cap and not a target: a spec asking for 2 still builds the homoleptic
    # compositions too.  The number matters because the enumeration is over MULTISETS of
    # kinds — a ligand with several protomers and pockets contributes a lot of kinds — so
    # this is the knob that decides how large a run is, and it is set deliberately rather
    # than discovered when the queue has 40,000 tasks in it.
    max_distinct_ligands: int = 1
    # Record how the products of this run reach each other, not just that they exist.
    #
    # A run that sweeps a RANGE of ligand copies builds a ladder — M(L), M(L)2, M(L)3 —
    # whose rungs differ by one addition, and that relationship is chemistry the registry
    # can hold (a `reactions` edge, D8) rather than something a reader has to infer from
    # two formulas.  Setting this plans the intermediate each rung is reached FROM and
    # builds the step with `assembly.join`, so the edge is a construction that was
    # actually performed rather than an assertion about two rows.
    #
    # Off by default: it adds the coordinatively unsaturated intermediates to the run, and
    # a spec written before this field existed did not ask for them.
    pathways: bool = False
    # `fill` or `range` (CO_LIGAND_COUNTS, D26).  `range` is the default for a new spec;
    # every spec older than v8 migrates to `fill`, which is what it planned.  Under
    # `range` the lower-hydration products keep the requested CN and geometry with the
    # uncovered vertices empty — they are asked for by this field, so `allow_unsaturated`
    # (which governs the case with NO co-ligand) does not gate them.
    co_ligand_counts: str = "range"
    # Under `range`, how many vertices a rung may leave empty: a product carries between
    # `full - window` and `full` co-ligands, and a pathway rung outside that is not
    # planned, so a ladder's root is the lowest co-ligand state inside it (D26).  Inert
    # under `fill`.
    co_ligand_window: int = 2
    run_mode: str = "construct"          # construct | ml_go | xtb_go | dft_go (no body)
    # Which ML potential `ml_go` means.  None = the machine's declared default.
    ml_model: str | None = None          # mace-mp-0 | mace-omol-0 | None
    note: str = ""

    def __post_init__(self) -> None:
        if self.run_mode not in RUN_MODES:
            raise ValueError(f"run_mode must be one of {RUN_MODES}, got {self.run_mode!r}")
        if not self.molecules:
            raise ValueError("a spec needs at least one molecule")
        if self.co_ligand_counts not in CO_LIGAND_COUNTS:
            raise ValueError(
                f"co_ligand_counts must be one of {CO_LIGAND_COUNTS}, got "
                f"{self.co_ligand_counts!r}: 'fill' puts a co-ligand on every free vertex, "
                f"'range' also builds each lower count")
        if (isinstance(self.co_ligand_window, bool)
                or not isinstance(self.co_ligand_window, int) or self.co_ligand_window < 0):
            raise ValueError(
                f"co_ligand_window must be a whole number >= 0 (how many vertices a rung "
                f"may leave empty), got {self.co_ligand_window!r}")
        if self.max_distinct_ligands < 1:
            raise ValueError(
                f"max_distinct_ligands must be at least 1, got {self.max_distinct_ligands}; "
                f"1 means homoleptic (one ligand kind per centre)")
        # Ranges are expanded at CONSTRUCTION, so every way of making a spec — this
        # constructor, `from_dict`, `dataclasses.replace` — reads `"1~3"` identically and
        # the object always holds the integers a run enumerated.  Normalising in one of
        # those paths only is how a spec built in a notebook would mean something
        # different from the same spec typed into the page.
        for name in ("coordination", "ligands_per_metal"):
            object.__setattr__(self, name, int_series(getattr(self, name), what=name))
        if not self.coordination:
            raise ValueError("coordination is empty; a metal run needs at least one CN")
        if not self.ligands_per_metal:
            raise ValueError("ligands_per_metal is empty; say how many copies to place")
        if self.ml_model is not None:
            from mofsbu.config import resolve_ml_backend

            resolve_ml_backend(self.ml_model)   # raises on a typo rather than defaulting
        # metals are DELIBERATELY optional: a purely molecular construction (a COF, an
        # organic cage) is a first-class case, not a degenerate one.

    # -- serialisation --------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {"spec_version": SPEC_VERSION, **asdict(self)}

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    @property
    def digest(self) -> str:
        """Content hash of the spec — two identical runs share it."""
        canonical = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode()).hexdigest()

    @classmethod
    def looks_like_spec(cls, d: dict[str, Any]) -> bool:
        """Is this JSON a build spec at all?  Used to keep other files out of listings."""
        return isinstance(d, dict) and "spec_version" in d and "molecules" in d

    @staticmethod
    def _migrate(d: dict[str, Any], version: int) -> dict[str, Any]:
        """Bring an older spec forward.

        Specs are saved files, so changing the dataclass without a migration silently
        breaks every spec already on disk — which is what happened between v1 and v2.
        An unknown FUTURE version still refuses: guessing at a newer schema is worse than
        saying no.
        """
        if version > SPEC_VERSION:
            raise ValueError(
                f"spec version {version} is newer than this build understands "
                f"({SPEC_VERSION}); update mofsbu rather than editing the file")
        if version <= 1:
            d.pop("relax_to", None)             # v1 field, replaced by run_mode
            d.setdefault("run_mode", "construct")
            d.setdefault("allow_unsaturated", False)
        if version <= 3:
            # v3 -> v4 adds `ml_model`.  It defaults to None, which resolves to the
            # machine's declared model — for a v3 spec that is MACE-MP-0 unless the
            # environment says otherwise, which is exactly what a v3 run did.  The
            # field is added rather than back-filled with "mace-mp-0" because a v3 spec
            # never expressed a choice and writing one in would invent provenance.
            d.setdefault("ml_model", None)
        if version <= 7:
            # v7 -> v8 adds `co_ligand_counts` (D26).  `fill` is what every earlier spec
            # planned, so an old spec re-run queues the same tasks; the new default,
            # `range`, applies only to a spec that is written at v8.
            d.setdefault("co_ligand_counts", "fill")
            d.setdefault("co_ligand_window", 2)             # inert under fill
        if version <= 6:
            # v6 -> v7 adds `pathways`.  False is what every earlier spec did — the rungs
            # of a ligand-count sweep were built independently and nothing recorded that
            # one is the other plus a ligand — so an old spec re-run plans the same tasks.
            d.setdefault("pathways", False)
        if version <= 5:
            # v5 -> v6 adds `max_distinct_ligands`.  It defaults to 1, which is exactly
            # what every earlier spec did — one molecule per coordination sphere — so an
            # old spec re-run produces the same structures it produced before.  Defaulting
            # it to anything else would silently multiply the size of every stored run.
            d.setdefault("max_distinct_ligands", 1)
        if version <= 4:
            # v4 -> v5 drops `MetalSpec.multiplicity`.  It was never derived from
            # `oxidation_state`/`spin_class`, so an old spec's stored value cannot be
            # trusted to be the reachable one for a spin-aware backend — dropping it
            # forces every metal centre back through `spin_class_multiplicity` instead
            # of resurrecting a number that may be exactly the bug being fixed.
            for m in d.get("metals", ()):
                m.pop("multiplicity", None)
        # v2 -> v3 adds no field: `xtb_go` joins RUN_MODES, so every v2 spec is already a
        # valid v3 one and there is nothing to rewrite.  The version still moves, because
        # a v3 spec saying `run_mode: xtb_go` must be REFUSED by a v2 build with "this is
        # newer than I understand" rather than with a bare "not a valid run_mode".
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> BuildSpec:
        d = dict(d)
        version = int(d.pop("spec_version", SPEC_VERSION))
        d = cls._migrate(d, version)
        unknown = set(d) - {f for f in cls.__dataclass_fields__}
        if unknown:
            raise ValueError(f"unknown spec field(s): {sorted(unknown)}")
        return cls(
            molecules=tuple(MoleculeSpec(**m) for m in d.pop("molecules", ())),
            metals=tuple(MetalSpec(**m) for m in d.pop("metals", ())),
            pocket=PocketPredicate(**(d.pop("pocket", None) or {})),
            # Passed through as written — a list, a range, or both.  `__post_init__`
            # expands it, so a hand-written `"1~3"` and the page's list become the same
            # spec through one parser rather than one per entry point.
            coordination=d.pop("coordination", (4, 6)),
            geometries=tuple(g) if (g := d.pop("geometries", None)) else None,
            ligands_per_metal=d.pop("ligands_per_metal", (1,)),
            binding=tuple(d.pop("binding", ("chelate", "mono"))),
            **d,
        )

    @classmethod
    def from_json(cls, text: str) -> BuildSpec:
        return cls.from_dict(json.loads(text))

    @classmethod
    def load(cls, path: str | Path) -> BuildSpec:
        return cls.from_json(Path(path).read_text(encoding="utf-8"))

    def save(self, path: str | Path) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(self.to_json() + "\n", encoding="utf-8")
        return p
