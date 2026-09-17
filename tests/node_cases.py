"""The M6 battery, read from the curated TSV.

A helper module rather than a test module, on the same footing as `build_routes.py`: the
rows are ground truth about the same nodes at two layers — the graph they must hash to,
and the geometry they must come out at — and a second copy of this reader would let those
two drift apart.

`data/reference/node_cases.tsv` is the interface for adding a case.  Nothing here needs
touching to add one.  Read `WORKPLAN_M6.md` §6 for what the battery is for.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

CASES_TSV = (Path(__file__).resolve().parents[1]
             / "data" / "reference" / "node_cases.tsv")

#: The bridging carboxylate's own geometry, as one group rather than seven columns that
#: could disagree.  Used only by `bridge_reach`, which is a REACHABILITY check on the
#: table — the numbers the model actually produces are measured and live in `d_mm_bridge`.
CARBOXYLATE_CO_A = 1.25
CARBOXYLATE_OCO_DEG = 125.0


def carboxylate_o_o() -> float:
    """O...O across a carboxylate — how far apart the two bridging donors sit."""
    return 2.0 * CARBOXYLATE_CO_A * np.sin(np.radians(CARBOXYLATE_OCO_DEG / 2.0))


def _num(value: str) -> float | None:
    """A column that may legitimately be absent.  `-` is absent, and absent is not zero."""
    return None if value.strip() in ("", "-") else float(value)


def _txt(value: str) -> str | None:
    return None if value.strip() in ("", "-") else value.strip()


@dataclass(frozen=True)
class NodeCase:
    """One row of the M6 battery: what the node is, and what a built one is measured against."""

    name: str
    fixture: str | None          # the `mofsbu.examples` graph it must hash to, or None
    mechanism: str               # "A", "B" or "A+B" — WORKPLAN_M6 §1
    metal: str
    n_metals: int
    bridge: str
    cn: int
    local_geometry: str
    mm_bond: bool
    n_bridges: int
    n_vacancies: int
    d_mm: float | None
    d_mm_lo: float | None
    d_mm_hi: float | None
    d_m_o: float | None
    mu_geometry: str | None
    d_m_mu: float | None
    mu_angle: float | None
    d_mm_bridge: float | None
    blocked_on: str | None
    source: str
    note: str

    @property
    def has_bridging_centre(self) -> bool:
        """Is there a single atom bridging the metals — mechanism B (WORKPLAN_M6 §4)?"""
        return "B" in self.mechanism

    @property
    def buildable_today(self) -> bool:
        return self.blocked_on is None

    @property
    def determined(self) -> bool:
        """Does the bridging centre alone fix M...M?  Then the skeleton is built, not searched."""
        return self.d_m_mu is not None and self.mu_angle is not None

    def mm_from_bridging_centre(self) -> float:
        """M...M implied by the bridging atom's own local geometry: 2 d sin(theta/2).

        For a symmetric bridge this is not an approximation, it is a definition — which is
        what makes it a check on the TABLE rather than on the placer.  It is also exactly
        the computation WORKPLAN_M6 §4 reports as landing the µ3 and µ4 skeletons on the
        literature values, so a row that fails it disagrees with a measurement.
        """
        if not self.determined:
            raise ValueError(f"{self.name}: no determined bridging centre")
        return 2.0 * self.d_m_mu * np.sin(np.radians(self.mu_angle / 2.0))

    def collision(self) -> float | None:
        """How far apart the two determinants are, where the row has both (§5).

        The bridging centre fixes M...M and so does the multi-atom bridge; in an
        oxo-centred cluster they disagree, and that disagreement — not a constraint solve —
        is M6's real work.  `None` where only one determinant exists, which is the
        paddlewheel case and the reason it never needs `place_multicentre`.
        """
        if self.d_mm_bridge is None or not self.determined:
            return None
        return abs(self.mm_from_bridging_centre() - self.d_mm_bridge)

    def bridge_reach_deg(self) -> float:
        """The M-M-O angle a syn-syn carboxylate needs to span this row's M...M.

        Planar construction: both donors lie in a plane containing the M-M axis, each at
        `d_m_o` from its own metal, and their separation along the axis is the
        carboxylate's own O...O bite.  Solving for the angle off the axis gives

            cos(reach) = (d_mm - d_oo) / (2 d_m_o)

        which is undefined exactly when the bridge cannot reach at all.  That is the point
        of computing it: an unreachable row is a typo in the table, and it would otherwise
        be discovered as a pile of clashes out of a placer that is working correctly.
        """
        if self.d_mm is None or self.d_m_o is None:
            raise ValueError(f"{self.name}: no M...M or M-O to measure reach against")
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
            name=r["name"], fixture=_txt(r["fixture"]), mechanism=r["mechanism"].strip(),
            metal=r["metal"], n_metals=int(r["n_metals"]), bridge=r["bridge"],
            cn=int(r["cn"]), local_geometry=r["local_geometry"],
            mm_bond=r["mm_bond"].strip() == "yes",
            n_bridges=int(r["n_bridges"]), n_vacancies=int(r["n_vacancies"]),
            d_mm=_num(r["d_mm"]), d_mm_lo=_num(r["d_mm_lo"]), d_mm_hi=_num(r["d_mm_hi"]),
            d_m_o=_num(r["d_m_o"]), mu_geometry=_txt(r["mu_geometry"]),
            d_m_mu=_num(r["d_m_mu"]), mu_angle=_num(r["mu_angle"]),
            d_mm_bridge=_num(r["d_mm_bridge"]), blocked_on=_txt(r["blocked_on"]),
            source=r["source"], note=r["note"]))
    return out


def case(name: str) -> NodeCase:
    for c in load_node_cases():
        if c.name == name:
            return c
    raise KeyError(f"no node case {name!r} in {CASES_TSV}")
