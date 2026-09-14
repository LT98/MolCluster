"""Pure rotation and alignment math.  No chemistry, no state, no dependencies but numpy.

The salvage ledger (`PLAN_implementation.md` §5) named this file as the destination for
`ebu_core._kabsch_align`, `_rot`, `_rotate_around_axis` and `_arbitrary_perpendicular`
— "port as-is (pure math)".  It was not created at the time, so each consumer ported the
piece it needed: `geometry.placer` grew six of these, `sites.frames` grew three written
differently, and `assembly.join` grew two more and imported two of the placer's PRIVATE
names across a module boundary.  Three copies of `unit` is not a maintenance problem
worth a milestone, but a fourth consumer reaching into a second module's underscore
namespace is where it starts to cost something.

Names here are public even though the module is private: `_linalg` is internal to the
package, and underscoring every function inside an already-underscored module only makes
call sites noisier.  Consumers alias them back to their existing local spellings, so the
call sites they appear in are unchanged.

**The two rotation forms are both here, and neither is written in terms of the other.**
That is deliberate and it is the one thing to know before editing this file.
`axis_rotation` builds the matrix; `rotate_vector` applies Rodrigues' formula directly.
They are the same rotation in exact arithmetic and they are NOT the same in floating
point — over 20000 random (vector, axis, angle) triples, 94% of results differ, by up to
1.8e-15.  That is far below anything chemically meaningful and far below the 1e-6 the
frame-reproducibility gate works at, but `sites.frames` uses the vector form to build
every stored frame, and rewriting it as `axis_rotation(axis, theta) @ v` would move the
last bits of every frame in the corpus for no gain.  A recipe version exists to be bumped
when the recipe changes; it should not have to absorb a refactor that changed nothing
anyone asked about.  So both expressions live here, once each, and the duplication that
was worth removing — six near-identical `unit` and `perpendicular` bodies — is gone.
"""
from __future__ import annotations

import numpy as np

#: Below this norm a vector has no direction to report, and the fallback is returned.
EPS = 1e-9

#: A vector this close to the x axis gets its perpendicular seeded from y instead, so the
#: cross product never collapses.  0.9 rather than 1.0 keeps it away from the degenerate
#: case with margin.
_AXIS_SEED_LIMIT = 0.9


def unit(v: np.ndarray) -> np.ndarray:
    """`v` normalised.  A zero vector has no direction, so z is returned rather than NaN."""
    n = float(np.linalg.norm(v))
    return v / n if n > EPS else np.array([0.0, 0.0, 1.0])


def perpendicular(v: np.ndarray) -> np.ndarray:
    """Some unit vector perpendicular to `v`.  Arbitrary, but deterministic.

    Deterministic matters more than the particular choice: it is what a frame's torsion
    zero is measured from, and a reference direction that moved between runs would make
    the same construct emit different conformers.
    """
    trial = (np.array([1.0, 0.0, 0.0]) if abs(v[0]) < _AXIS_SEED_LIMIT
             else np.array([0.0, 1.0, 0.0]))
    return unit(np.cross(v, trial))


def axis_rotation(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rotation MATRIX of `theta` radians about `axis` (right-handed)."""
    axis = unit(axis)
    x, y, z = axis
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ])


def rotate_vector(v: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    """`v` rotated by `theta` about `axis`, via Rodrigues directly.

    Equivalent to `axis_rotation(axis, theta) @ v` in exact arithmetic and not bit-for-bit
    in floating point — see the module docstring.  This is the form every stored site frame
    was built with, so it stays the form they are rebuilt with.
    """
    axis = unit(axis)
    return (v * np.cos(theta)
            + np.cross(axis, v) * np.sin(theta)
            + axis * np.dot(axis, v) * (1 - np.cos(theta)))


def rotation_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector `a` onto unit vector `b`.

    The antiparallel case is the one worth having written down: the cross product vanishes
    there, so the axis is chosen from an arbitrary perpendicular and the rotation is a
    half turn about it.  Without that branch the formula divides by zero exactly when a
    ligand has to be flipped end for end, which is a common case and not an edge one.
    """
    a, b = unit(a), unit(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if float(np.linalg.norm(v)) < EPS:
        if c > 0:
            return np.eye(3)
        trial = (np.array([1.0, 0.0, 0.0]) if abs(a[0]) < _AXIS_SEED_LIMIT
                 else np.array([0.0, 1.0, 0.0]))
        return axis_rotation(unit(np.cross(a, trial)), np.pi)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1 / (1 + c))


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation + translation least-squares fitting point set `p` onto `q`.

    The `det` correction is what keeps it a ROTATION: without it the SVD is free to return
    a reflection, which fits the points just as well and turns a molecule into its mirror
    image.
    """
    cp, cq = p.mean(axis=0), q.mean(axis=0)
    h = (p - cp).T @ (q - cq)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, cq - rot @ cp


def angle_between(a: np.ndarray, b: np.ndarray) -> float:
    """Angle between two vectors, in DEGREES.  Clamped, because `arccos` is unforgiving."""
    cos = float(np.dot(unit(a), unit(b)))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))
