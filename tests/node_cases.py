"""The polynuclear-node targets M6 is measured against, read from the curated TSV.

A helper module rather than a test module, on the same footing as `build_routes.py`:
the rows are ground truth about the same nodes at two layers — the graph they must hash
to, and the geometry they must come out at — and a second copy of this reader would let
those two drift apart.

`data/reference/node_cases.tsv` is the interface for adding a case.  Nothing here needs
touching to add one.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

CASES_TSV = (Path(__file__).resolve().parents[1]
             / "data" / "reference" / "node_cases.tsv")

#: The bridging carboxylate, as one geometry for every row.  It is a property of the
#: GROUP and not of the node, so it is here rather than repeated down a column where
#: seven copies could disagree.
CARBOXYLATE_CO_A = 1.25
CARBOXYLATE_OCO_DEG = 125.0


def carboxylate_o_o() -> float:
    """O...O across a carboxylate — how far apart the two bridging donors are.

    This is the span the bridge has to reach with, and it does not change with the metal:
    a paddlewheel at 2.09 A and a mu4-oxo node at 3.17 A are both bridged by this same
    2.22 A bite, tilted differently against the M-M axis.
    """
    return 2.0 * CARBOXYLATE_CO_A * np.sin(np.radians(CARBOXYLATE_OCO_DEG / 2.0))


def _num(value: str) -> float | None:
    """A column that may legitimately be absent.  `-` is absent, and absent is not zero."""
    return None if value.strip() in ("", "-") else float(value)


@dataclass(frozen=True)
class NodeCase:
    """One polynuclear node: what it is, and the numbers a built one is measured against."""

    name: str
    fixture: str | None          # the `mofsbu.examples` graph it must hash to, or None
    metal: str
    n_metals: int
    bridge: str
    cn: int
    local_geometry: str
    mm_bond: bool
    n_bridges: int
    n_vacancies: int
    d_mm: float
    d_mm_lo: float
    d_mm_hi: float
    d_m_o: float
    d_m_mu: float | None
    mu_angle: float | None
    source: str
    note: str

    @property
    def has_central_oxo(self) -> bool:
        return self.d_m_mu is not None

    def mm_from_oxo(self) -> float:
        """M...M implied by the central oxo alone: 2 d sin(theta/2).

        For a symmetric mu-oxo node this is not an approximation, it is the definition —
        which is what makes it a check on the TABLE rather than on the placer.
        """
        if not self.has_central_oxo:
            raise ValueError(f"{self.name} has no central oxo")
        return 2.0 * self.d_m_mu * np.sin(np.radians(self.mu_angle / 2.0))

    def bridge_tilt_deg(self) -> float:
        """The M-M-O angle a syn-syn carboxylate needs to span this M...M distance.

        Planar construction: both donors lie in a plane containing the M-M axis, each at
        `d_m_o` from its own metal, and their separation along the axis is the
        carboxylate's own O...O bite.  Solving for the angle off the axis gives

            cos(tilt) = (d_mm - d_oo) / (2 d_m_o)

        which is undefined exactly when the bridge cannot reach at all.  That is the
        point of computing it: an unreachable row is a typo in the table, and it would
        otherwise be discovered as a pile of clashes out of the placer.
        """
        cos = (self.d_mm - carboxylate_o_o()) / (2.0 * self.d_m_o)
        if not -1.0 <= cos <= 1.0:
            raise ValueError(
                f"{self.name}: a syn-syn carboxylate cannot span M...M {self.d_mm} A "
                f"at M-O {self.d_m_o} A — the bite is {carboxylate_o_o():.2f} A")
        return float(np.degrees(np.arccos(cos)))


def load_node_cases(path: Path = CASES_TSV) -> list[NodeCase]:
    rows = [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")]
    out: list[NodeCase] = []
    for r in csv.DictReader(rows, delimiter="\t"):
        out.append(NodeCase(
            name=r["name"],
            fixture=None if r["fixture"].strip() == "-" else r["fixture"].strip(),
            metal=r["metal"], n_metals=int(r["n_metals"]), bridge=r["bridge"],
            cn=int(r["cn"]), local_geometry=r["local_geometry"],
            mm_bond=r["mm_bond"].strip() == "yes",
            n_bridges=int(r["n_bridges"]), n_vacancies=int(r["n_vacancies"]),
            d_mm=float(r["d_mm"]), d_mm_lo=float(r["d_mm_lo"]),
            d_mm_hi=float(r["d_mm_hi"]), d_m_o=float(r["d_m_o"]),
            d_m_mu=_num(r["d_m_mu"]), mu_angle=_num(r["mu_angle"]),
            source=r["source"], note=r["note"]))
    return out


def case(name: str) -> NodeCase:
    for c in load_node_cases():
        if c.name == name:
            return c
    raise KeyError(f"no node case {name!r} in {CASES_TSV}")
