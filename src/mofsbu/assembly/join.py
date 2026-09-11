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
    #: Per-geometry state from `sites.state.refresh_state`, keyed by `(atom_idx, slot)`.
    #: The slot is in the key because a metal carries several vacancies on ONE atom, so
    #: an atom-keyed dict silently keeps whichever of them was inserted last.  Absent
    #: means "not computed for this block", which is not the same as "nothing is open".
    state: dict[tuple[int, int], Any] | None = None

    @staticmethod
    def state_key(site: Any) -> tuple[int, int]:
        """The key a site's state is stored under.  One definition, used by both sides."""
        return (site.atom_idx, getattr(site, "slot", 0))

    def open_sites(self) -> tuple[Site, ...]:
        """The sites a join may actually use — donors AND vacant metal vertices.

        Open means: not already dative-bonded to a metal, and not sterically walled off
        (`sites.state.SiteStatus`).  Both halves matter — "both ends are unoccupied" is
        the test §6.3 explicitly says is not sufficient, and a site the metal cannot
        reach is not a site a join can use however free its valence looks.

        The two kinds come back in one tuple on purpose.  A join needs a donor on one
        block and somewhere on the other block to put it, and on a metal that somewhere is
        a vacancy — so `compatible(a, b)` can take two `Site`s and ask one frame-alignment
        question, rather than needing a separate accessor and a separate predicate for the
        metal's side of every bond.

        A block with no state raises rather than returning every site.  Treating unknown
        as open is how an assembly step would confidently join onto a buried donor.
        """
        if self.state is None:
            raise NotBuiltYet(
                "this BuildingBlock carries no site state, so 'open' is unknown. Run "
                "sites.state.refresh_state on its geometry (or load it with "
                "registry.get_site_state) — returning every site would silently treat "
                "occupied and buried donors as available.")
        from mofsbu.sites.state import SiteStatus

        return tuple(s for s in self.sites
                     if getattr(self.state.get(self.state_key(s)), "status", None)
                     is SiteStatus.OPEN)

    def open_vacancies(self) -> tuple[Site, ...]:
        """Just the metal's usable empty vertices — where an incoming ligand can go."""
        return tuple(s for s in self.open_sites() if s.is_vacancy)

    def open_donors(self) -> tuple[Site, ...]:
        """Just the usable donor atoms — what this block can offer another one."""
        return tuple(s for s in self.open_sites() if not s.is_vacancy)


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
