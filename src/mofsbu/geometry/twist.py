"""Twist a chelating ligand so its two donors can reach two vertices of one metal (#15).

The placer moves each ligand as a rigid body, in the shape it was embedded in. A free
ethylenediamine embeds anti (N...N 3.77 A) and a free 2,2'-bipyridine trans (3.65 A): a
metal cannot span either, and docking both donors pushes the backbone through the metal.
The fix is the ligand's own flexibility. This module turns the single bonds on the path
between the two donors until the donor-donor distance fits the two target vertices.

Only the torsions move. Bond lengths and angles are the embed's own, so any bending a
ring really needs is left to relaxation, as with the rest of a rigid-ligand model.

The search keeps the ligand's character:
* a bond between two sp3 atoms turns in 120 deg steps, so a staggered chain stays
  staggered (en: anti to gauche);
* a bond touching an sp2 atom turns in 180 deg steps, so a conjugated link stays planar
  (bipy: trans to cis, acac: extended to U).
A 30 deg grid on every bond is the fallback when neither finds a clash-free fit.
"""
from __future__ import annotations

import itertools
from collections.abc import Sequence

import numpy as np
from rdkit import Chem

__all__ = ["twist_bonds", "twist_to_bite"]

#: Fallback grid, in steps per full turn.
FINE_STEPS = 12

#: Closest two atoms 1-5 or further apart may come after a twist, in A: heavy-heavy,
#: heavy-H, H-H.  Coarse on purpose: this rejects a ligand folded into itself, not a
#: contact a relaxation would tidy up.
_MIN_CONTACT = {2: 2.6, 1: 2.0, 0: 1.6}


def twist_bonds(mol: Chem.Mol, donors: Sequence[int]) -> list[tuple[int, int]]:
    """The single, acyclic bonds on the shortest path between two donors.

    A bond touching a donor is left out: turning it moves that donor's own substituents,
    not the distance to the other donor.
    """
    path = Chem.GetShortestPath(mol, int(donors[0]), int(donors[1]))
    out = []
    for a, b in zip(path[1:-2], path[2:-1], strict=True):
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing():
            out.append((a, b))
    return out


def _side(mol: Chem.Mol, i: int, j: int) -> list[int]:
    """The atoms on j's side of the bond i-j."""
    seen, todo = {j}, [j]
    while todo:
        a = todo.pop()
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            k = n.GetIdx()
            if k not in seen and not (a == j and k == i):
                seen.add(k)
                todo.append(k)
    return sorted(seen)


def _rotate(xyz: np.ndarray, i: int, j: int, side: list[int], theta: float) -> np.ndarray:
    axis = xyz[j] - xyz[i]
    axis = axis / np.linalg.norm(axis)
    k = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]],
                  [-axis[1], axis[0], 0.0]])
    rot = np.eye(3) + np.sin(theta) * k + (1.0 - np.cos(theta)) * k @ k
    out = xyz.copy()
    out[side] = (rot @ (xyz[side] - xyz[j]).T).T + xyz[j]
    return out


def _folded(xyz: np.ndarray, n_heavy: np.ndarray, far: np.ndarray) -> bool:
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    limit = np.vectorize(_MIN_CONTACT.get)(n_heavy[:, None] + n_heavy[None])
    return bool((far & (d < limit)).any())


def _steps(mol: Chem.Mol, bond: tuple[int, int]) -> int:
    sp3 = Chem.HybridizationType.SP3
    both_sp3 = all(mol.GetAtomWithIdx(a).GetHybridization() == sp3 for a in bond)
    return 3 if both_sp3 else 2


def twist_to_bite(mol: Chem.Mol, xyz: np.ndarray, donors: Sequence[int], target: float,
                  *, window: tuple[float, float]) -> tuple[np.ndarray, list[dict]]:
    """Coordinates with the donor-donor distance moved toward `target`, and the twists.

    Returns the input unchanged and `[]` when there is nothing to turn, or when no twist
    both lands inside `window` (donor-donor distances, A) and leaves the ligand unfolded.
    Deterministic: the same ligand and target always give the same twists.
    """
    bonds = twist_bonds(mol, donors)
    if not bonds:
        return xyz, []
    sides = [_side(mol, a, b) for a, b in bonds]
    n_heavy = np.array([int(a.GetAtomicNum() > 1) for a in mol.GetAtoms()])
    far = Chem.GetDistanceMatrix(mol) >= 4
    d0, d1 = int(donors[0]), int(donors[1])
    for grids in ([_steps(mol, b) for b in bonds], [FINE_STEPS] * len(bonds)):
        best = None
        for ks in itertools.product(*(range(n) for n in grids)):
            x = xyz
            for (a, b), side, k, n in zip(bonds, sides, ks, grids, strict=True):
                if k:
                    x = _rotate(x, a, b, side, 2.0 * np.pi * k / n)
            dist = float(np.linalg.norm(x[d0] - x[d1]))
            miss = abs(dist - target)
            if not window[0] <= dist <= window[1] or (best and miss >= best[0]):
                continue
            if _folded(x, n_heavy, far):
                continue
            best = (miss, x, ks)
        if best is not None:
            twists = [{"bond": [int(a), int(b)], "deg": round(360.0 * k / n, 1)}
                      for (a, b), k, n in zip(bonds, best[2], grids, strict=True) if k]
            return (best[1], twists) if twists else (xyz, [])
    return xyz, []
