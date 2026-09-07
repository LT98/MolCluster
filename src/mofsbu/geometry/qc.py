"""Geometry quality control.  Ported from `legacy/geometry_qc.py`.

Cheap checks that catch the ways a constructed geometry is obviously wrong: atoms
inside one another, metal-ligand bonds at implausible lengths, rings bent out of plane
by a bad placement.  Nothing here is a substitute for a relaxation; it is a filter so
that nonsense never reaches one.

**A refusal has to say what it refused.**  This module used to report
"QC FAILED: 2 clash(es), closest 1.40 A" and nothing else — no elements, no atom
indices, no threshold, no indication of which ligand each atom belonged to.  That single
string is what a whole run's worth of rejected candidates collapsed into, and it made
two completely different situations look identical: a ligand genuinely too bulky for the
site, and an iodide placed at oxygen's bond length by a table that did not know about
iodine.  Every finding here therefore carries the atoms, the measured value, the limit
it was measured against, and the ligand each atom came from.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np

VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
       "Cl": 1.75, "Br": 1.85, "I": 1.98, "Se": 1.90, "As": 1.85, "Te": 2.06, "B": 1.92,
       "Si": 2.10,
       "Li": 1.82, "Na": 2.27, "K": 2.75, "Mg": 1.73, "Ca": 2.31, "Al": 1.84,
       "Sc": 2.11, "Ti": 2.00, "V": 2.00, "Cr": 2.00, "Mn": 2.00, "Fe": 2.00,
       "Co": 2.00, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Ru": 2.00, "Pd": 1.63,
       "Ag": 1.72, "Cd": 1.58, "Pt": 1.75, "Au": 1.66, "Hg": 1.55, "Zr": 2.00}
DEFAULT_VDW = 1.90


class Clash(NamedTuple):
    """Two non-bonded atoms inside each other's van der Waals radii.

    The first three fields are `(i, j, distance)` in that order deliberately: that was
    the whole of the old tuple, so callers that index `c[2]` for the distance keep
    working.  The rest is what was missing.
    """

    i: int
    j: int
    distance: float
    limit: float                 # what it was measured against
    sym_i: str = ""
    sym_j: str = ""
    owner_i: str = ""            # which ligand atom i belongs to ("" = unknown)
    owner_j: str = ""

    @property
    def overlap(self) -> float:
        """How far inside the limit, angstrom.  The severity, not just the fact."""
        return self.limit - self.distance

    def describe(self) -> str:
        a = f"{self.sym_i}{self.i}" + (f"[{self.owner_i}]" if self.owner_i else "")
        b = f"{self.sym_j}{self.j}" + (f"[{self.owner_j}]" if self.owner_j else "")
        return (f"{a}...{b} {self.distance:.2f} A "
                f"(limit {self.limit:.2f}, over by {self.overlap:.2f})")

    def to_dict(self) -> dict:
        return {"i": self.i, "j": self.j, "sym_i": self.sym_i, "sym_j": self.sym_j,
                "owner_i": self.owner_i, "owner_j": self.owner_j,
                "distance": round(self.distance, 3), "limit": round(self.limit, 3),
                "overlap": round(self.overlap, 3)}


class BadBond(NamedTuple):
    """A metal-donor bond that is not the length it was placed at.

    `(atom, distance)` first, for the same compatibility reason as `Clash`.
    """

    atom: int
    distance: float
    target: float = 0.0
    tol: float = 0.0
    sym: str = ""
    owner: str = ""
    source: str = ""             # where `target` came from — see geometry.distances

    def describe(self) -> str:
        who = f"{self.sym}{self.atom}" + (f"[{self.owner}]" if self.owner else "")
        src = f", {self.source}" if self.source else ""
        return (f"M-{who} {self.distance:.2f} A vs target {self.target:.2f} "
                f"+/- {self.tol:.2f}{src}")

    def to_dict(self) -> dict:
        return {"atom": self.atom, "sym": self.sym, "owner": self.owner,
                "distance": round(self.distance, 3), "target": round(self.target, 3),
                "tol": self.tol, "source": self.source}


@dataclass
class QCReport:
    ok: bool = True
    clashes: list[Clash] = field(default_factory=list)
    bad_bonds: list[BadBond] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def code(self) -> str:
        """A stable machine-readable reason, so a run's rejections can be GROUPED.

        Free text cannot be grouped: "2 clash(es), closest 1.40 A" and
        "2 clash(es), closest 1.41 A" are the same finding and different strings.
        """
        if self.ok:
            return "ok"
        if self.clashes and self.bad_bonds:
            return "qc_clash_and_bond"
        if self.clashes:
            return "qc_clash"
        if self.bad_bonds:
            return "qc_bond_length"
        return "qc_other"

    @property
    def worst_clash(self) -> Clash | None:
        return max(self.clashes, key=lambda c: c.overlap, default=None)

    def elements_involved(self) -> list[str]:
        """Which elements the failure is about.  This is the field that would have said
        `I` the first time a halide co-ligand was refused."""
        out: set[str] = set()
        for c in self.clashes:
            out.update(s for s in (c.sym_i, c.sym_j) if s)
        out.update(b.sym for b in self.bad_bonds if b.sym)
        return sorted(out)

    def to_dict(self) -> dict:
        worst = self.worst_clash
        return {"ok": self.ok, "code": self.code,
                "n_clashes": len(self.clashes),
                "worst_clash": worst.distance if worst else None,
                "worst_overlap": round(worst.overlap, 3) if worst else None,
                "clashes": [c.to_dict() for c in self.clashes[:20]],
                "bad_bonds": [b.to_dict() for b in self.bad_bonds[:20]],
                "elements": self.elements_involved(),
                "notes": self.notes}

    def __str__(self) -> str:
        if self.ok:
            return "QC ok"
        bits = []
        if self.clashes:
            worst = self.worst_clash
            bits.append(f"{len(self.clashes)} clash(es); worst {worst.describe()}")
        if self.bad_bonds:
            bits.append(f"{len(self.bad_bonds)} bad M-L bond(s): "
                        + "; ".join(b.describe() for b in self.bad_bonds[:3]))
        return "QC FAILED: " + "; ".join(bits + self.notes)


def check_clashes(symbols: list[str], coords: np.ndarray, bonded: set[tuple[int, int]],
                  *, scale: float = 0.62, scale_h: float = 0.50,
                  owners: list[str] | None = None) -> list[Clash]:
    """Non-bonded atom pairs sitting inside each other's van der Waals radii."""
    n = len(symbols)
    out: list[Clash] = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded or (j, i) in bonded:
                continue
            d = float(np.linalg.norm(coords[i] - coords[j]))
            ri = VDW.get(symbols[i], DEFAULT_VDW)
            rj = VDW.get(symbols[j], DEFAULT_VDW)
            s = scale_h if "H" in (symbols[i], symbols[j]) else scale
            limit = s * (ri + rj)
            if d < limit:
                out.append(Clash(i, j, d, limit, symbols[i], symbols[j],
                                 owners[i] if owners else "",
                                 owners[j] if owners else ""))
    return out


def check_metal_bonds(coords: np.ndarray, metal_idx: int, donor_idxs: list[int],
                      *, d_ml: float | list[float], tol: float = 0.45,
                      symbols: list[str] | None = None,
                      owners: list[str] | None = None,
                      sources: list[str] | None = None) -> list[BadBond]:
    """Metal-donor bonds that are not the length they were placed at.

    `d_ml` may be one number for every donor, or one per donor.  The per-donor form is
    the point: a centre carrying a water and an iodide has two different target lengths,
    and checking both against a single number is how the iodide came to look like the
    problem.
    """
    targets = ([float(d_ml)] * len(donor_idxs)
               if isinstance(d_ml, (int, float))
               else [float(x) for x in d_ml])
    if len(targets) != len(donor_idxs):
        raise ValueError(f"{len(targets)} target distances for {len(donor_idxs)} donors")
    out: list[BadBond] = []
    for k, (atom, target) in enumerate(zip(donor_idxs, targets)):
        dist = float(np.linalg.norm(coords[metal_idx] - coords[atom]))
        if abs(dist - target) > tol:
            out.append(BadBond(atom, dist, target, tol,
                               symbols[atom] if symbols else "",
                               owners[atom] if owners else "",
                               sources[k] if sources else ""))
    return out


def qc(symbols: list[str], coords: np.ndarray, bonded: set[tuple[int, int]],
       *, metal_idx: int | None = None, donor_idxs: list[int] | None = None,
       d_ml: float | list[float] = 2.05,
       owners: list[str] | None = None,
       sources: list[str] | None = None) -> QCReport:
    report = QCReport()
    report.clashes = check_clashes(symbols, coords, bonded, owners=owners)
    if metal_idx is not None and donor_idxs:
        report.bad_bonds = check_metal_bonds(coords, metal_idx, donor_idxs, d_ml=d_ml,
                                             symbols=symbols, owners=owners,
                                             sources=sources)
    report.ok = not report.clashes and not report.bad_bonds
    return report
