"""Placeholder 3D coordinates for structures that have no real geometry yet.

This is NOT chemistry.  It is a spring embedding of the typed graph, scaled to
plausible bond lengths, so that a structure can be looked at before the multi-centre
placer (M6) exists.  Its one guarantee is the one that matters for correctness: the
atoms, their count and their order are exactly the graph's, so a layout can never
describe a molecule the graph does not.

Anything produced here is recorded at `Fidelity.RAW` with method `spring-embed-3d`, so
it is impossible to mistake for a relaxed structure in the registry or the viewer.
"""
from __future__ import annotations

import networkx as nx
import numpy as np

from mofsbu.graph import TypedGraph

TARGET_BOND = 1.5      # angstrom, roughly a C-C single bond


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
    return coords - coords.mean(axis=0)


def to_xyz(g: TypedGraph, *, seed: int = 0, comment: str | None = None) -> str:
    coords = embed(g, seed=seed)
    lines = [str(len(g)), comment or f"{g.name or 'structure'} (spring-embed-3d placeholder)"]
    for k, atom in enumerate(g.nodes()):
        x, y, z = coords[k]
        lines.append(f"{g.label(atom).element:<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"
