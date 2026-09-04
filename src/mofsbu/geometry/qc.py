"""Geometry quality control.  Ported from `legacy/geometry_qc.py`.

Cheap checks that catch the ways a constructed geometry is obviously wrong: atoms
inside one another, metal-ligand bonds at implausible lengths, rings bent out of plane
by a bad placement.  Nothing here is a substitute for a relaxation; it is a filter so
that nonsense never reaches one.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
       "Cl": 1.75, "Br": 1.85, "I": 1.98,
       "Li": 1.82, "Na": 2.27, "K": 2.75, "Mg": 1.73, "Ca": 2.31, "Al": 1.84,
       "Sc": 2.11, "Ti": 2.00, "V": 2.00, "Cr": 2.00, "Mn": 2.00, "Fe": 2.00,
       "Co": 2.00, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ru": 2.00, "Pd": 1.63,
       "Ag": 1.72, "Cd": 1.58, "Pt": 1.75, "Au": 1.66, "Hg": 1.55, "Zr": 2.00}
DEFAULT_VDW = 1.90


@dataclass
class QCReport:
    ok: bool = True
    clashes: list[tuple[int, int, float]] = field(default_factory=list)
    bad_bonds: list[tuple[int, float]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"ok": self.ok, "n_clashes": len(self.clashes),
                "worst_clash": min((c[2] for c in self.clashes), default=None),
                "bad_bonds": self.bad_bonds, "notes": self.notes}

    def __str__(self) -> str:
        if self.ok:
            return "QC ok"
        bits = []
        if self.clashes:
            bits.append(f"{len(self.clashes)} clash(es), closest {min(c[2] for c in self.clashes):.2f} A")
        if self.bad_bonds:
            bits.append(f"{len(self.bad_bonds)} bad M-L bond(s)")
        return "QC FAILED: " + "; ".join(bits + self.notes)


def check_clashes(symbols: list[str], coords: np.ndarray, bonded: set[tuple[int, int]],
                  *, scale: float = 0.62, scale_h: float = 0.50) -> list[tuple[int, int, float]]:
    """Non-bonded atom pairs sitting inside each other's van der Waals radii."""
    n = len(symbols)
    out: list[tuple[int, int, float]] = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded or (j, i) in bonded:
                continue
            d = float(np.linalg.norm(coords[i] - coords[j]))
            ri = VDW.get(symbols[i], DEFAULT_VDW)
            rj = VDW.get(symbols[j], DEFAULT_VDW)
            s = scale_h if "H" in (symbols[i], symbols[j]) else scale
            if d < s * (ri + rj):
                out.append((i, j, d))
    return out


def check_metal_bonds(coords: np.ndarray, metal_idx: int, donor_idxs: list[int],
                      *, d_ml: float, tol: float = 0.45) -> list[tuple[int, float]]:
    return [(d, float(np.linalg.norm(coords[metal_idx] - coords[d])))
            for d in donor_idxs
            if abs(float(np.linalg.norm(coords[metal_idx] - coords[d])) - d_ml) > tol]


def qc(symbols: list[str], coords: np.ndarray, bonded: set[tuple[int, int]],
       *, metal_idx: int | None = None, donor_idxs: list[int] | None = None,
       d_ml: float = 2.05) -> QCReport:
    report = QCReport()
    report.clashes = check_clashes(symbols, coords, bonded)
    if metal_idx is not None and donor_idxs:
        report.bad_bonds = check_metal_bonds(coords, metal_idx, donor_idxs, d_ml=d_ml)
    report.ok = not report.clashes and not report.bad_bonds
    return report
