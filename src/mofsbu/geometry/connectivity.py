"""Does a geometry still have the bonds its structure's graph says it has?  (gate E1)

A relaxation moves atoms under a potential that knows nothing about the typed graph the
geometry is filed under: a water drifts off its metal, a carboxylate O swings on, a proton
crosses an H-bond.  The result is then a realisation of a DIFFERENT species, and an energy
read from it under the old identity is wrong in a way nothing downstream can see.

This module only measures.  Bonding is judged from distances against the same ideal
lengths construction uses (`geometry.distances`): a graph edge is broken past `BROKEN`
times its ideal length, a non-edge is formed inside `FORMED` times its.  The thresholds are
the ones the E-MH study measured with (`docs/WORKPLAN_energy.md` §3a, metric (a)).
`registry.api.put_geometry` records the verdict; energy selection reads it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mofsbu._types import MofsbuError
from mofsbu.geometry.distances import (COVALENT_RADII, DEFAULT_COVALENT,
                                       metal_donor_distance, metal_metal_distance)
from mofsbu.graph._types import EdgeType, TypedGraph

#: A graph edge longer than this multiple of its ideal length is broken.
BROKEN = 1.30
#: A non-edge shorter than this multiple of the ideal length it would have is formed.
FORMED = 1.10


@dataclass(frozen=True)
class BondChange:
    """One pair whose bonding disagrees with the graph: measured against what it was judged by."""

    i: int
    j: int
    sym_i: str
    sym_j: str
    etype: EdgeType
    distance: float
    ideal: float
    limit: float

    def to_dict(self) -> dict:
        return {"i": self.i, "j": self.j, "sym_i": self.sym_i, "sym_j": self.sym_j,
                "etype": self.etype.value, "distance": round(self.distance, 3),
                "ideal": round(self.ideal, 3), "limit": round(self.limit, 3)}


@dataclass(frozen=True)
class ConnectivityReport:
    """`broken` graph edges and `formed` non-edges; `changed` if either is non-empty."""

    broken: tuple[BondChange, ...]
    formed: tuple[BondChange, ...]

    @property
    def changed(self) -> bool:
        return bool(self.broken or self.formed)

    def describe(self) -> str:
        if not self.changed:
            return "connectivity intact"
        parts = [f"{c.sym_i}{c.i}-{c.sym_j}{c.j} broken at {c.distance:.2f} A "
                 f"(limit {c.limit:.2f})" for c in self.broken]
        parts += [f"{c.sym_i}{c.i}-{c.sym_j}{c.j} formed at {c.distance:.2f} A "
                  f"(limit {c.limit:.2f})" for c in self.formed]
        return "connectivity changed: " + "; ".join(parts)

    def to_dict(self) -> dict:
        return {"changed": self.changed, "broken": [c.to_dict() for c in self.broken],
                "formed": [c.to_dict() for c in self.formed],
                "thresholds": {"broken": BROKEN, "formed": FORMED}}


def _radius(sym: str) -> float:
    return COVALENT_RADII.get(sym, DEFAULT_COVALENT)


def _ideal(sym_i: str, sym_j: str, etype: EdgeType, i_is_metal: bool) -> float:
    if etype is EdgeType.DATIVE:
        metal, donor = (sym_i, sym_j) if i_is_metal else (sym_j, sym_i)
        return metal_donor_distance(metal, donor).value
    if etype is EdgeType.METAL_METAL:
        return metal_metal_distance(sym_i, sym_j).value
    return _radius(sym_i) + _radius(sym_j)


def connectivity_check(symbols, coords, graph: TypedGraph) -> ConnectivityReport:
    """Compare the bonding implied by `coords` with `graph`.

    Atom `k` of `coords` must be atom `k` of the graph; a different element at any index is
    refused, because then the comparison is between two orderings, not two bondings.
    A METAL-METAL edge is declared, never inferred (D25): a declared one can break, but two
    metals coming close forms nothing.  A metal-H contact forms nothing either (agostic
    contacts are not dative bonds in this graph).
    """
    symbols = [str(s) for s in symbols]
    xyz = np.asarray(coords, dtype=float).reshape(-1, 3)
    n = len(graph)
    if len(symbols) != n or len(xyz) != n:
        raise MofsbuError(
            f"geometry has {len(xyz)} atoms and {len(symbols)} symbols; the graph has {n}")
    wrong = [k for k in range(n) if graph.label(k).element != symbols[k]]
    if wrong:
        k = wrong[0]
        raise MofsbuError(
            f"atom {k} is {symbols[k]} in the geometry and {graph.label(k).element} in the "
            f"graph ({len(wrong)} mismatches); a connectivity check needs one atom order")

    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    metal = [graph.is_metal(k) for k in range(n)]
    edges = {(u, v): et for u, v, et in graph.edges()}

    broken = []
    for (u, v), et in edges.items():
        ideal = _ideal(symbols[u], symbols[v], et, metal[u])
        if d[u, v] > BROKEN * ideal:
            broken.append(BondChange(u, v, symbols[u], symbols[v], et, float(d[u, v]),
                                     ideal, BROKEN * ideal))
    formed = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in edges or (metal[i] and metal[j]):
                continue
            if metal[i] != metal[j]:
                if "H" in (symbols[i], symbols[j]):
                    continue
                et = EdgeType.DATIVE
            else:
                et = EdgeType.COVALENT
            ideal = _ideal(symbols[i], symbols[j], et, metal[i])
            if d[i, j] < FORMED * ideal:
                formed.append(BondChange(i, j, symbols[i], symbols[j], et, float(d[i, j]),
                                         ideal, FORMED * ideal))
    return ConnectivityReport(tuple(broken), tuple(formed))
