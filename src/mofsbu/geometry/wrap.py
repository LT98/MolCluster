"""A conformer of a ligand wrapped around one metal by three or more donors (#15).

A rigid ligand cannot wrap a metal: dien, tren and EDTA embed as open chains, and turning
their bonds alone cannot close a hexadentate wrap, which also needs bond angles to give.
This is distance geometry with the metal in the molecule: a stand-in metal atom bonded to
every donor, bounds that put each donor at its own metal distance and at the separation its
vertices have from each other, and the vertices other ligands will take as fixed points
nothing may come near.  The embed chooses the torsions and the bends together; a short MMFF
minimisation with the donors pinned on their vertices then tidies the rest.

Distances cannot tell a conformer from its mirror image, so an embed that only fits the
targets as a reflection is discarded and the next seed tried: the arrangement asked for is
one hand (Delta or Lambda), and reflecting the ligand would also invert its own stereocentres.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from rdkit import Chem, DistanceGeometry
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

__all__ = ["wrap_conformer"]

#: Bound widths, A: metal-donor, donor-donor and donor-vertex.
TOL_MD, TOL_DD = 0.04, 0.08
#: Closest a non-donor atom may come to the metal or to another ligand's vertex, A.
CLEAR_HEAVY, CLEAR_H = 2.6, 2.0
#: Seeds tried, in order; the first fit inside `GOOD_FIT` stops the search.
SEEDS = tuple(range(1, 41))
GOOD_FIT = 0.1
#: Largest out-of-plane RMS, A, an aromatic ring may have after the wrap.  A planar
#: ligand forced into an arrangement it cannot reach (fac terpy) puckers by ~0.2.
MAX_PUCKER = 0.05
#: How far, A, the final minimisation may move an atom from where the embed put it.
PIN_SLACK = 0.25
#: Stand-in metal for the embed's atom typing only; nothing of it is kept.
_STAND_IN = 26


def wrap_conformer(mol: Chem.Mol, donors: Sequence[int], targets: np.ndarray, *,
                   vertices: np.ndarray | None = None) -> np.ndarray | None:
    """Ligand coordinates with each donor on its target, metal at the origin, or None.

    `targets` are the donors' positions in donor order; `vertices` the positions other
    ligands will take.  None if no seed embeds in the right hand, or if the only fit
    bends an aromatic ring.  Deterministic.
    """
    donors = [int(d) for d in donors]
    vertices = np.zeros((0, 3)) if vertices is None else vertices
    n = mol.GetNumAtoms()
    rw = Chem.RWMol(Chem.Mol(mol))
    rw.RemoveAllConformers()
    metal = rw.AddAtom(Chem.Atom(_STAND_IN))
    for d in donors:
        rw.AddBond(metal, d, Chem.BondType.ZERO)
    marks = [rw.AddAtom(Chem.Atom(2)) for _ in vertices]       # He: a fixed point
    for m in marks:
        rw.AddBond(metal, m, Chem.BondType.ZERO)   # one fragment, or the bounds are ignored
    big = rw.GetMol()
    big.UpdatePropertyCache(strict=False)
    bm = rdDistGeom.GetMoleculeBoundsMatrix(big)

    def bound(i: int, j: int, lo: float, hi: float) -> None:
        a, b = min(i, j), max(i, j)
        bm[a][b], bm[b][a] = hi, lo

    fixed = [(metal, np.zeros(3))] + list(zip(donors, targets, strict=True)) \
        + list(zip(marks, vertices, strict=True))
    for x, (i, p) in enumerate(fixed):
        for j, q in fixed[x + 1:]:
            r = float(np.linalg.norm(p - q))
            tol = TOL_MD if metal in (i, j) else TOL_DD
            bound(i, j, max(r - tol, 0.1), r + tol)
    for a in range(n):
        if a in donors:
            continue
        lo = CLEAR_HEAVY if mol.GetAtomWithIdx(a).GetAtomicNum() > 1 else CLEAR_H
        for point in [metal, *marks]:
            i, j = min(a, point), max(a, point)
            bm[j][i] = max(bm[j][i], lo)
            bm[i][j] = max(bm[i][j], lo + 0.1)
    if not DistanceGeometry.DoTriangleSmoothing(bm):
        return None

    best = None
    for seed in SEEDS:
        params = rdDistGeom.EmbedParameters()
        params.randomSeed = seed
        params.useRandomCoords = True
        params.useExpTorsionAnglePrefs = False
        params.useBasicKnowledge = True
        params.SetBoundsMat(bm)
        if rdDistGeom.EmbedMolecule(big, params) < 0:
            continue
        x = np.array(big.GetConformer().GetPositions())
        x = x - x[metal]
        u, _s, vt = np.linalg.svd(x[donors].T @ targets)
        if np.linalg.det(vt.T @ u.T) < 0:
            continue                                   # the other hand
        y = x @ (vt.T @ u.T).T
        miss = float(np.linalg.norm(y[donors] - targets, axis=1).max())
        if best is None or miss < best[0]:
            best = (miss, y[:n])
        if miss < GOOD_FIT:
            break
    if best is None:
        return None
    xyz = _pin(mol, best[1], donors, targets)
    return None if _puckered(mol, xyz) else xyz


def _puckered(mol: Chem.Mol, xyz: np.ndarray) -> bool:
    """Has the wrap bent an aromatic ring out of its plane?  Then it is not a conformer."""
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            p = xyz[list(ring)] - xyz[list(ring)].mean(axis=0)
            if np.linalg.svd(p)[1][-1] / np.sqrt(len(ring)) > MAX_PUCKER:
                return True
    return False


def _pin(mol: Chem.Mol, xyz: np.ndarray, donors: list[int],
         targets: np.ndarray) -> np.ndarray:
    """Put the donors exactly on their targets and let MMFF settle the rest around them."""
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    xyz = xyz.copy()
    xyz[donors] = targets
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, p.tolist())
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(out)
    if props is None:
        return xyz
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(out, props)
    for d in donors:
        ff.AddFixedPoint(d)
    # MMFF cannot see the metal, so every other atom is held near where the embed put it,
    # clear of the metal and the other vertices; it may only tidy its bonds.
    for i in range(out.GetNumAtoms()):
        if i not in donors:
            ff.MMFFAddPositionConstraint(i, PIN_SLACK, 100.0)
    ff.Minimize(maxIts=200)
    return np.array(out.GetConformer().GetPositions())
