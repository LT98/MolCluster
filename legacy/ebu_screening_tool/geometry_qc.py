"""
geometry_qc.py -- fast post-placement geometry validation.

A cheap safety net for the EBU pipeline: not every metal/ligand/CN/geometry
combination that enumerate_sbus() is asked to build is actually physically
achievable (e.g. two bulky ligands forced onto adjacent sites of a tight
5-coordinate center). Rather than silently exporting a distorted or clashing
structure, run three fast checks after placement and report a verdict with
reasons attached:

  1. Non-bonded clash  -- any pair of atoms (excluding 1-2 and 1-3 bonded
     neighbors) closer than `clash_scale` times the sum of their vdW radii.
  2. M-donor bond sanity -- every intended metal-donor contact should land
     within a reasonable window of the requested d_ml; if placement (DG
     wrap, rigid fallback, or the UFF relax step) left one far off target,
     that EBU did not actually assemble as requested.
  3. Rigid-ring planarity -- aromatic/other rigid rings should stay flat.
     A ring that was planar in the free ligand conformer but is measurably
     puckered after placement indicates the placement step bent something
     that should not bend (the "extreme torsion" failure mode).

All three are O(n^2) at most on EBU-sized systems (tens to ~100 atoms), so
this is cheap enough to run on every candidate in a combinatorial sweep.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

PT = Chem.GetPeriodicTable()

# Fallback vdW radii (Angstrom) for elements PT might not have a good value
# for in some RDKit builds; covers everything likely to appear in an EBU.
_VDW_FALLBACK = {
    'H': 1.10, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
    'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'B': 1.92,
}


def _vdw_radius(symbol: str) -> float:
    try:
        r = PT.GetRvdw(PT.GetAtomicNumber(symbol))
        if r and r > 0:
            return r
    except Exception:
        pass
    return _VDW_FALLBACK.get(symbol, 1.7)


@dataclass
class QCResult:
    is_valid: bool
    issues: List[str] = field(default_factory=list)
    min_clash_ratio: Optional[float] = None   # smallest (distance / vdW_sum) seen; <1 means clash

    def __str__(self):
        if self.is_valid:
            return 'VALID'
        return 'INVALID: ' + '; '.join(self.issues)


def _bonded_graph_distance(mol: Chem.Mol, max_depth: int = 2):
    """Return dict atom_idx -> set of atom_idx within `max_depth` bonds (BFS)."""
    n = mol.GetNumAtoms()
    close = [set() for _ in range(n)]
    adj = [[nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors()] for i in range(n)]
    for start in range(n):
        frontier = {start}
        seen = {start}
        for _ in range(max_depth):
            nxt = set()
            for a in frontier:
                for nb in adj[a]:
                    if nb not in seen:
                        nxt.add(nb)
                        seen.add(nb)
            close[start] |= nxt
            frontier = nxt
    return close


def check_clashes(mol3d: Chem.Mol, clash_scale: float = 0.65, clash_scale_h: float = 0.50) -> QCResult:
    """
    Flag any non-bonded (>=1-4 relationship) atom pair closer than
    clash_scale * (vdW_i + vdW_j). clash_scale ~0.6-0.7 is a common
    permissive threshold (allows normal steric contact, catches real overlap).
    """
    conf = mol3d.GetConformer()
    n = mol3d.GetNumAtoms()
    close = _bonded_graph_distance(mol3d, max_depth=2)
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(n)])
    syms = [mol3d.GetAtomWithIdx(i).GetSymbol() for i in range(n)]
    radii = np.array([_vdw_radius(s) for s in syms])

    worst_ratio = np.inf
    issues = []
    for i in range(n):
        for j in range(i + 1, n):
            if j in close[i]:
                continue
            d = np.linalg.norm(pos[i] - pos[j])
            vdw_sum = radii[i] + radii[j]
            ratio = d / vdw_sum
            if ratio < worst_ratio:
                worst_ratio = ratio
            # H atoms have soft, small vdW shells and the least reliable
            # positions in an unrelaxed placement, so a near-threshold M...H
            # or H...H graze is not a real clash (it relaxes out). Use a more
            # permissive cutoff for any H-involving pair; heavy-heavy overlaps
            # (the genuine steric failures) keep the stricter clash_scale.
            eff_scale = clash_scale_h if (syms[i] == 'H' or syms[j] == 'H') else clash_scale
            if ratio < eff_scale:
                issues.append(f'clash {syms[i]}{i}-{syms[j]}{j}: {d:.2f} Ang '
                              f'({ratio:.2f}x vdW sum {vdw_sum:.2f})')
    issues.sort(key=lambda s: float(s.split('(')[1].split('x')[0]))
    return QCResult(is_valid=(len(issues) == 0), issues=issues[:8],
                     min_clash_ratio=float(worst_ratio) if np.isfinite(worst_ratio) else None)


def check_ml_bonds(mol3d: Chem.Mol, metal_idx: int, donor_idxs: List[int],
                    d_ml: float, tol: float = 0.6) -> QCResult:
    """Every intended metal-donor contact should be within `tol` Ang of d_ml."""
    conf = mol3d.GetConformer()
    m_pos = np.array(conf.GetAtomPosition(metal_idx))
    issues = []
    for d_i in donor_idxs:
        dist = np.linalg.norm(np.array(conf.GetAtomPosition(d_i)) - m_pos)
        if abs(dist - d_ml) > tol:
            sym = mol3d.GetAtomWithIdx(d_i).GetSymbol()
            issues.append(f'M-donor {sym}{d_i}: {dist:.2f} Ang (target {d_ml:.2f} +/- {tol})')
    return QCResult(is_valid=(len(issues) == 0), issues=issues)


def check_ring_planarity(free_mol: Chem.Mol, placed_mol: Chem.Mol,
                          atom_map: Optional[List[int]] = None,
                          max_deviation: float = 0.25) -> QCResult:
    """
    Compare each ring's out-of-plane RMSD between the free ligand conformer
    and its placed conformer. A ring that was flat and is no longer flat
    signals the placement step bent something rigid (the ATF-style
    "extreme torsion" failure).

    atom_map: optional index mapping free_mol atom idx -> placed_mol atom idx
              (identity mapping assumed if None, i.e. same atom order).
    """
    ri = free_mol.GetRingInfo()
    if ri.NumRings() == 0:
        return QCResult(is_valid=True)
    free_conf = free_mol.GetConformer()
    placed_conf = placed_mol.GetConformer()
    issues = []
    for ring in ri.AtomRings():
        if len(ring) < 5:
            continue
        def plane_rmsd(conf, idxs):
            pts = np.array([list(conf.GetAtomPosition(i)) for i in idxs])
            centroid = pts.mean(axis=0)
            _, _, vt = np.linalg.svd(pts - centroid)
            normal = vt[2]
            dev = (pts - centroid) @ normal
            return float(np.sqrt((dev ** 2).mean()))

        free_dev = plane_rmsd(free_conf, ring)
        mapped = ring if atom_map is None else [atom_map[i] for i in ring]
        placed_dev = plane_rmsd(placed_conf, mapped)
        if placed_dev - free_dev > max_deviation:
            issues.append(f'ring {list(ring)}: planarity RMSD {free_dev:.2f} -> {placed_dev:.2f} Ang '
                          f'(delta {placed_dev - free_dev:.2f} > {max_deviation})')
    return QCResult(is_valid=(len(issues) == 0), issues=issues)


def qc_ebu(mol3d: Chem.Mol, metal_symbol: str, ligands, d_ml: float,
           clash_scale: float = 0.65, ring_max_deviation: float = 0.15,
           clash_scale_h: float = 0.50) -> QCResult:
    """
    Run the full checklist on a placed EBU. `ligands` is the list of
    LigandBuilder instances used to build it (same order/offsets as
    EBUBuilder.build), used to (a) find the metal/donor indices and
    (b) compare each ligand's rings against its own free conformer.
    """
    an = PT.GetAtomicNumber(metal_symbol)
    metal_idx = next(a.GetIdx() for a in mol3d.GetAtoms() if a.GetAtomicNum() == an)

    donor_idxs = []
    offset = metal_idx + 1
    ring_issues = []
    for lig in ligands:
        for d in lig.donor_indices:
            donor_idxs.append(offset + d)
        ring_result = check_ring_planarity(lig.mol, mol3d,
                                            atom_map=list(range(offset, offset + lig.mol.GetNumAtoms())),
                                            max_deviation=ring_max_deviation)
        if not ring_result.is_valid:
            ring_issues += [f'{lig.name}: {msg}' for msg in ring_result.issues]
        offset += lig.mol.GetNumAtoms()

    clash = check_clashes(mol3d, clash_scale=clash_scale, clash_scale_h=clash_scale_h)
    mlbonds = check_ml_bonds(mol3d, metal_idx, donor_idxs, d_ml)

    all_issues = mlbonds.issues + ring_issues + clash.issues
    return QCResult(is_valid=(len(all_issues) == 0), issues=all_issues,
                     min_clash_ratio=clash.min_clash_ratio)
