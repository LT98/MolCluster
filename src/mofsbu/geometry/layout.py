"""Placeholder 3D coordinates for structures that have no real geometry yet.

This is NOT chemistry.  It is a spring embedding of the typed graph, refined so every
bond sits at its ideal length and every non-bonded pair well clear of one, so that a
structure can be looked at before the multi-centre placer (M6) exists.  Its guarantees are
the ones that matter for correctness: the atoms, their count and their order are exactly
the graph's, and its bonding is the graph's by `geometry.connectivity` (E1), so a layout
can never describe a molecule the graph does not.

Anything produced here is recorded at `Fidelity.RAW` with method `spring-embed-3d`, so
it is impossible to mistake for a relaxed structure in the registry or the viewer.
"""
from __future__ import annotations

import networkx as nx
import numpy as np

from mofsbu.geometry.connectivity import FORMED, _ideal
from mofsbu.graph import EdgeType, TypedGraph

TARGET_BOND = 1.5      # angstrom, roughly a C-C single bond
#: A non-bonded pair is pushed out to this multiple of the length it would have if bonded,
#: comfortably past `connectivity.FORMED`.
CLEAR = 1.30
#: A metal and an H, or two metals, have no bonded length to scale; kept this far apart.
CLEAR_FLOOR = 2.0


def embed(g: TypedGraph, *, seed: int = 0) -> np.ndarray:
    """(n_atoms, 3) coordinates in the graph's own atom order."""
    nxg = g.to_networkx()
    if len(g) == 1:
        return np.zeros((1, 3))
    pos = nx.spring_layout(nxg, dim=3, seed=seed, iterations=200)
    coords = np.array([pos[i] for i in g.nodes()], dtype=float)

    bonded = [(u, v) for u, v, _et in g.edges()]
    if bonded:
        idx = {atom: k for k, atom in enumerate(g.nodes())}
        lengths = [np.linalg.norm(coords[idx[u]] - coords[idx[v]]) for u, v in bonded]
        median = float(np.median(lengths)) or 1.0
        coords *= TARGET_BOND / median
        coords = _refine(g, coords)
    return coords - coords.mean(axis=0)


def _refine(g: TypedGraph, coords: np.ndarray, *, iterations: int = 4000,
            step: float = 0.1) -> np.ndarray:
    """Gradient descent on bond lengths (exact) and non-bonded clearances (lower bounds)."""
    n = len(g)
    sym = [g.label(k).element for k in g.nodes()]
    metal = [g.is_metal(k) for k in g.nodes()]
    edges = {(u, v): et for u, v, et in g.edges()}
    ii, jj = np.triu_indices(n, k=1)
    target = np.empty(len(ii))
    exact = np.zeros(len(ii), dtype=bool)
    for k, (i, j) in enumerate(zip(ii, jj)):
        et = edges.get((int(i), int(j)))
        if et is not None:
            target[k], exact[k] = _ideal(sym[i], sym[j], et, metal[i]), True
        elif metal[i] and metal[j] or (metal[i] != metal[j] and "H" in (sym[i], sym[j])):
            target[k] = CLEAR_FLOOR
        else:
            et = EdgeType.DATIVE if metal[i] != metal[j] else EdgeType.COVALENT
            target[k] = max(CLEAR, FORMED) * _ideal(sym[i], sym[j], et, metal[i])
    x = coords.copy()
    for _ in range(iterations):
        diff = x[ii] - x[jj]
        d = np.linalg.norm(diff, axis=1) + 1e-9
        err = d - target
        active = exact | (err < 0)
        if not np.any(np.abs(err[active]) > 1e-3):
            break
        f = np.where(active, err / d, 0.0)[:, None] * diff
        grad = np.zeros_like(x)
        np.add.at(grad, ii, f)
        np.subtract.at(grad, jj, f)
        x -= step * grad
    return x


def to_xyz(g: TypedGraph, *, seed: int = 0, comment: str | None = None) -> str:
    coords = embed(g, seed=seed)
    lines = [str(len(g)), comment or f"{g.name or 'structure'} (spring-embed-3d placeholder)"]
    for k, atom in enumerate(g.nodes()):
        x, y, z = coords[k]
        lines.append(f"{g.label(atom).element:<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"
