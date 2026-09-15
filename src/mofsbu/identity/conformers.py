"""L3: which conformers are the same one.  Provenance-primary, geometry-verifier (D11).

The label on a conformer is **the choice vector that built it**, not a hash of its
coordinates.  That is the whole of D11 and it is what makes near-degenerate cis and trans
survive a clustering pass: they carry different choice vectors, so no threshold ever gets
the chance to merge them.

Geometry has exactly two jobs, and neither of them is naming anything:

* **collapse Kind-A duplicates** — one choice vector sampled twice, differing only by
  stochastic noise.  These are the same conformer and only the geometry can say so.
* **reconcile divergence and convergence** — one choice leading to two minima, or several
  choices falling into one.  Both really happen and neither is visible in the provenance.

So the order is fixed and it matters: **choice-vector dedup runs first**, and geometric
clustering only ever looks at what survives it.  Running the clustering first would let a
threshold decide identity, which is the mistake D11 was written to prevent, and it is also
the declared escape hatch for the conformer-explosion risk.

**On theta_geom.**  The threshold is a named constant here and its value is NOT yet
calibrated — see `docs/WORKPLAN_M5.md` for the measurement and why.  The short version:
over the M5 fixture set at RAW fidelity the Kind-A population has *exactly zero* spread,
because a join is a deterministic function of its choice vector and absorbs the ligand's
embedding noise entirely.  There is no stochastic-duplicate peak to put a threshold above,
so a number chosen now would be calibrated against a distribution with no noise in it.
`DEFAULT_THETA_GEOM` is therefore a placeholder wearing its provenance: it is used only
where a caller asks for clustering and does not supply one, and `cluster` takes the value
as an argument precisely so the calibration can arrive without this module changing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

import numpy as np

from mofsbu.assembly.choice import digest_of
from mofsbu.graph._types import EdgeType, TypedGraph
from mofsbu.graph.canon import canonical_index_map
from mofsbu.geometry._linalg import kabsch

#: Placeholder, not a calibration.  Two geometries closer than this over the rigid core are
#: the same conformer.  See the module docstring: the M5 fixture set cannot set this number
#: because its Kind-A spread is identically zero, and the fixture set that can is M7's
#: relaxed one.  Passed explicitly wherever it matters so the eventual value lands in one
#: place and nothing has to be rewritten to accept it.
DEFAULT_THETA_GEOM = 0.25        # angstrom, rigid core + coordination sphere

#: Exact-match tolerance.  Distinct from theta_geom and not a weaker version of it: this
#: one asks "did these two come out of the same arithmetic", which is a question about
#: floating point, and the answer does not depend on any chemistry.
IDENTICAL_A = 1e-9


def rotatable(g: TypedGraph, u: int, v: int) -> bool:
    """Can this bond twist?  Single, acyclic, non-aromatic — `sites.model`'s definition.

    Deliberately the same rule the pocket search uses, so "rigid" means one thing in the
    codebase.  It is order-blind (D15 keeps bond order out of identity), which means a
    delocalised carboxylate's C-O reads as rotatable when chemically it is not — see
    `docs/ISSUES.md`. The consequence here is a core slightly smaller than it should be,
    which makes this test more permissive, never less.
    """
    import networkx as nx

    if g.edge_type(u, v) is not EdgeType.COVALENT:
        return False
    nxg = g.to_networkx()
    return not any(u in cycle and v in cycle for cycle in nx.cycle_basis(nxg))


def core_atoms(g: TypedGraph) -> list[int]:
    """The rigid core plus the coordination sphere — what an L3 comparison looks at.

    Section 4.2 is explicit that whole-molecule RMSD is a bad sole trigger: a distant
    floppy tail gives a huge number and means nothing, while a small twist that parks a
    site open or blocked barely moves it.  So the comparison is restricted to the part
    that cannot twist away — the metals, everything dative-bonded to them, and each
    ligand's rigid body — reached by walking out from the metals and refusing to cross a
    rotatable bond.
    """
    keep = set(g.metals())
    frontier = list(keep)
    while frontier:
        atom = frontier.pop()
        for nb in g.neighbors(atom):
            if nb in keep:
                continue
            if g.edge_type(atom, nb) is not EdgeType.COVALENT or not rotatable(g, atom, nb):
                keep.add(nb)
                frontier.append(nb)
    return sorted(keep)


def core_rmsd(ga: TypedGraph, ca: Any, gb: TypedGraph, cb: Any) -> float | None:
    """Best-fit RMSD over the two structures' rigid cores, or None if they differ.

    Correspondence is by CANONICAL index, so the comparison does not assume the two were
    built in the same atom order.  One caveat worth knowing: canonical labelling resolves
    symmetry-equivalent atoms arbitrarily but deterministically, so a highly symmetric
    ligand can report a larger RMSD than an optimal atom matching would.  It is stable and
    it never merges two things that differ — it only ever fails to merge two that do not.
    """
    ca, cb = np.asarray(ca, dtype=float), np.asarray(cb, dtype=float)
    ia, ib = canonical_index_map(ga), canonical_index_map(gb)
    core_a = {ia[i]: i for i in core_atoms(ga)}
    core_b = {ib[i]: i for i in core_atoms(gb)}
    if set(core_a) != set(core_b) or len(core_a) < 3:
        return None
    shared = sorted(core_a)
    pa = np.array([ca[core_a[c]] for c in shared])
    pb = np.array([cb[core_b[c]] for c in shared])
    rot, trans = kabsch(pa, pb)
    moved = (rot @ pa.T).T + trans
    return float(np.sqrt(np.mean(np.sum((moved - pb) ** 2, axis=1))))


@dataclass(frozen=True)
class Conformer:
    """One candidate: where it came from, and optionally what it looks like."""

    choice_vector: dict[str, Any]
    graph: TypedGraph | None = None
    geometry: Any | None = None
    energy: float | None = None

    @property
    def label(self) -> str:
        """The provenance label — the answer, before geometry is consulted at all."""
        return digest_of(self.choice_vector)


@dataclass(frozen=True)
class Cluster:
    """A set of candidates judged to be one conformer, and why."""

    representative: Conformer
    members: tuple[Conformer, ...] = ()
    merged_by: str = "choice-vector"        # choice-vector | geometry
    reasons: tuple[str, ...] = ()
    #: Set only when one choice vector diverged into several wells, so the shared
    #: provenance label can still address them separately.
    suffix: str = ""

    @property
    def label(self) -> str:
        return self.representative.label + self.suffix

    @property
    def size(self) -> int:
        return len(self.members)


def l3_conformer_id(choice_vector: Any | None = None, geom: Any | None = None) -> str:
    """The conformer's identity: its choice vector's digest.  Geometry is not consulted.

    Provenance-primary means exactly this — the label is a function of how the thing was
    built, and it does not move when the coordinates do.  A conformer with no choice vector
    has no L3 identity and gets `""`, which is correct rather than degraded: something that
    arrived without provenance cannot claim a provenance-based key.  `geom` stays in the
    signature because `cluster` is the geometry half and callers reach for this name first.
    """
    if not choice_vector:
        return ""
    if isinstance(choice_vector, str):
        return choice_vector
    return digest_of(choice_vector)


def cluster(candidates: Sequence[Conformer], *, theta_geom: float = DEFAULT_THETA_GEOM,
            energy_window: float | None = None) -> list[Cluster]:
    """Group candidates into conformers: choice vectors first, geometry second.

    The two passes do different jobs and the order is the design:

    1. **By choice vector.**  Same provenance, same conformer — no geometry consulted, so
       no threshold can split a branch that a builder deliberately made.
    2. **By geometry, across the survivors.**  This is the only pass that can merge two
       DIFFERENT choice vectors, and it merges them only when their cores coincide within
       `theta_geom`: the convergence case D11 names, several choices falling into one well.

    `energy_window` (kcal/mol) gates the second pass: a geometry far off the
    representative's energy is not a duplicate of it, it is a worse structure that happens
    to look similar, and merging the two would hide a bad geometry inside a good one's
    identity.  Candidates without energies are never gated out — an absent number is not
    evidence of anything (D18's rule, applied here).
    """
    by_label: dict[str, list[Conformer]] = {}
    for candidate in candidates:
        by_label.setdefault(candidate.label, []).append(candidate)

    clusters: list[Cluster] = []
    for label, group in by_label.items():
        # Divergence: one choice vector can land in more than one well, and when it does
        # those really are different conformers that happen to share a provenance label.
        # The suffix is an ordinal within the group in input order — a disambiguator, not
        # a claim about which well is which.
        wells = _split_by_geometry(group, theta_geom)
        for k, well in enumerate(wells):
            suffix = "" if len(wells) == 1 else f"#{k}"
            reasons = []
            if len(well) > 1:
                reasons.append(f"{len(well)} sample(s) of one choice vector")
            if len(wells) > 1:
                reasons.append(
                    f"one choice vector diverged into {len(wells)} wells more than "
                    f"{theta_geom} A apart; this is {k}")
            clusters.append(Cluster(well[0], tuple(well), "choice-vector",
                                    tuple(reasons), suffix))

    merged: list[Cluster] = []
    for current in clusters:
        target = _converges_onto(current, merged, theta_geom, energy_window)
        if target is None:
            merged.append(current)
            continue
        i = merged.index(target)
        rmsd = core_rmsd(target.representative.graph, target.representative.geometry,
                         current.representative.graph, current.representative.geometry)
        merged[i] = Cluster(
            target.representative, target.members + current.members, "geometry",
            target.reasons + (
                f"{current.label[:16]} converged onto this one at "
                f"{rmsd:.4f} A over the rigid core (theta_geom {theta_geom})",),
            target.suffix)
    return merged


def _split_by_geometry(group: Sequence[Conformer], theta_geom: float,
                       ) -> list[list[Conformer]]:
    """Split one choice vector's samples into the wells they actually fell into.

    Members with no geometry cannot be split and stay with the first well — an absent
    measurement is not evidence that something diverged.
    """
    wells: list[list[Conformer]] = []
    for candidate in group:
        if candidate.graph is None or candidate.geometry is None:
            if not wells:
                wells.append([])
            wells[0].append(candidate)
            continue
        for well in wells:
            head = well[0]
            if head.graph is None or head.geometry is None:
                continue
            rmsd = core_rmsd(head.graph, head.geometry, candidate.graph, candidate.geometry)
            if rmsd is not None and rmsd <= theta_geom:
                well.append(candidate)
                break
        else:
            wells.append([candidate])
    return wells or [list(group)]


def _converges_onto(current: Cluster, existing: Iterable[Cluster], theta_geom: float,
                    energy_window: float | None) -> Cluster | None:
    """The already-kept cluster this one falls into, if any."""
    a = current.representative
    if a.graph is None or a.geometry is None:
        return None                     # no geometry: provenance is the only answer there is
    for other in existing:
        b = other.representative
        if b.graph is None or b.geometry is None:
            continue
        if energy_window is not None and a.energy is not None and b.energy is not None:
            if abs(a.energy - b.energy) > energy_window:
                continue                # similar shape, different energy: not a duplicate
        rmsd = core_rmsd(b.graph, b.geometry, a.graph, a.geometry)
        if rmsd is not None and rmsd <= theta_geom:
            return other
    return None
