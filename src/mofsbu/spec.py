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

SPEC_VERSION = 6

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
    run_mode: str = "construct"          # construct | ml_go | xtb_go | dft_go (no body)
    # Which ML potential `ml_go` means.  None = the machine's declared default.
    ml_model: str | None = None          # mace-mp-0 | mace-omol-0 | None
    note: str = ""

    def __post_init__(self) -> None:
        if self.run_mode not in RUN_MODES:
            raise ValueError(f"run_mode must be one of {RUN_MODES}, got {self.run_mode!r}")
        if not self.molecules:
            raise ValueError("a spec needs at least one molecule")
        if self.max_distinct_ligands < 1:
            raise ValueError(
                f"max_distinct_ligands must be at least 1, got {self.max_distinct_ligands}; "
                f"1 means homoleptic (one ligand kind per centre)")
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
            coordination=tuple(d.pop("coordination", (4, 6))),
            geometries=tuple(g) if (g := d.pop("geometries", None)) else None,
            ligands_per_metal=tuple(d.pop("ligands_per_metal", (1,))),
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
