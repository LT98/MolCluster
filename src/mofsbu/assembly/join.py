"""Joining building blocks — the recursive assembly operation (D1, D13).  NOT IMPLEMENTED.

Ground rule 7: the signatures are the design decision and they are settled here, so
callers can be written against them now.  The bodies are scheduled work (M5).  Every one
of these raises; none returns a plausible-looking value, because a stub that quietly
answers is indistinguishable from a working function until the answers are wrong.

What has to exist before these can be written:

* `sites.frames` — DONE.  A join needs both partners' frames and a torsion well, and
  those are stored rather than re-derived.
* `compatible()` — a frame-alignment feasibility test between two open sites, the same
  question `sites.model.chelate_pockets` answers within one molecule, generalised to two.
* `join()` — align, bond, recompute the product's open sites by inheriting through the
  atom map, and emit the choice-vector that regenerates it.
* `geometry.placer.place_multicentre` — inter-centre constraints (M-M distance, bridge
  bite angle).  The plan's declared headline cost (M6).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from mofsbu._types import AtomMap, MofsbuError
from mofsbu.graph._types import TypedGraph
from mofsbu.sites.model import Site


class NotBuiltYet(NotImplementedError, MofsbuError):
    """Raised by a settled interface whose body is still scheduled work."""


@dataclass(frozen=True)
class BuildingBlock:
    """A molecule, a metal, or an assembled fragment — one interface for all three (D1)."""

    graph: TypedGraph
    sites: tuple[Site, ...]
    geometry: Any | None = None
    structure_id: int | None = None

    def open_sites(self) -> tuple[Site, ...]:
        raise NotBuiltYet(
            "BuildingBlock.open_sites needs site_state (open vs occupied) — M4, second half"
        )


@dataclass(frozen=True)
class Compatibility:
    feasible: bool
    strain: float
    reason: str


@dataclass(frozen=True)
class JoinResult:
    block: BuildingBlock
    atom_map: AtomMap
    choice_vector: dict[str, Any]
    strain: float


def compatible(a: Site, b: Site, *, partner: Any | None = None) -> Compatibility:
    """Can these two open sites be joined at acceptable strain?

    Frame-alignment feasibility, not merely "both are open" — the same convergence
    question `chelate_pockets` answers inside one molecule, across two blocks.
    """
    raise NotBuiltYet("assembly.compatible — M5")


def join(a: BuildingBlock, b: BuildingBlock, site_a: Site, site_b: Site, *,
         mode: str = "mono", torsion_well: int = 0, seed: int = 0) -> JoinResult:
    """Join two blocks at a pair of sites and return the product.

    Must: align by frame, add the typed bond, carry an atom map, inherit the parents'
    sites minus those consumed, and emit the choice-vector that regenerates the result.
    """
    raise NotBuiltYet("assembly.join — M5 (the recursive assembly operation)")


def grow(seed_block: BuildingBlock, partners: tuple[BuildingBlock, ...], *,
         degree: int, max_products: int = 1000, seed: int = 0) -> list[JoinResult]:
    """Repeatedly join partners onto a block, `degree` additions deep.

    This is what "iterate cluster formation to a certain degree" means: degree 1 is one
    metal's coordination sphere (which `geometry.placer.place_mononuclear` already does),
    degree 2+ is polynuclear growth, which needs `join` and the multi-centre placer.
    """
    raise NotBuiltYet(
        f"assembly.grow(degree={degree}) — M5/M6.  Degree 1 is served today by "
        "geometry.placer.place_mononuclear; polynuclear growth needs assembly.join and "
        "geometry.placer.place_multicentre, neither of which is built."
    )
