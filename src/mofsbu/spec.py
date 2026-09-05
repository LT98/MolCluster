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

SPEC_VERSION = 3

# What to do with each structure once it is constructed.  Whether a mode can RUN is a
# property of the machine, not of the spec: `energy.relax.mode_status()` asks the
# backends, and the planner refuses with that reason.  `dft_go` is still a settled name
# with no body (ground rule 7) — there is no external code wired up.
RUN_MODES = ("construct", "ml_go", "xtb_go", "dft_go")


@dataclass(frozen=True)
class MoleculeSpec:
    name: str
    smiles: str
    multiplicity: int = 1
    max_deprotonations: int | None = None


@dataclass(frozen=True)
class MetalSpec:
    symbol: str
    oxidation_state: int = 2
    spin_class: str = "ls"
    multiplicity: int = 1


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
    run_mode: str = "construct"          # construct | ml_go | xtb_go | dft_go (no body)
    note: str = ""

    def __post_init__(self) -> None:
        if self.run_mode not in RUN_MODES:
            raise ValueError(f"run_mode must be one of {RUN_MODES}, got {self.run_mode!r}")
        if not self.molecules:
            raise ValueError("a spec needs at least one molecule")
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
