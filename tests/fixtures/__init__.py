"""Test fixtures — the example structures now live in the package.

Kept as a shim so tests keep importing `fixtures`, while scripts and the demo
registry draw the same graphs from `mofsbu.examples`.
"""
from __future__ import annotations

from mofsbu.examples import *  # noqa: F401,F403
from mofsbu.examples import ALL  # noqa: F401


def realised_xyz(g, shift: float = 0.0, comment: str = "test") -> str:
    """XYZ text whose bonding is `g`'s (the refined layout), translated by `shift` along x.

    Anything stored above RAW is checked against its graph (E1), so a test geometry has to
    realise the graph it is filed under; `shift` gives distinct coordinates (and blob) with
    the same bonding.
    """
    from mofsbu.geometry.layout import embed

    coords = embed(g)
    lines = [str(len(g)), comment]
    for k, atom in enumerate(g.nodes()):
        x, y, z = coords[k]
        lines.append(f"{g.label(atom).element} {x + shift:.5f} {y:.5f} {z:.5f}")
    return "\n".join(lines) + "\n"
