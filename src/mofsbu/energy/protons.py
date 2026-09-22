"""Which structures are the same molecule with a proton removed, and what that costs.

The registry stores every protomer of a ligand as a first-class row and records nothing
connecting them: `runner`'s `ligand` branch calls `put_structure` with tags and no
provenance, so tHQ and tHQ(-1H) are two unrelated structures that happen to share a name.
The relationship is recoverable — one is the other minus an H and minus a charge — and
once it is written down as an edge the existing reference scheme prices it.

Why a couple rather than a bare proton: `AH -> A- + H+` cannot be weighed. A proton has
no electrons, so no gas-phase method has an energy for it and `check_balance` reads charge
off a graph with no atoms. `AH + H2O -> A- + H3O+` balances, and because no metal-donor
bond changes across the arrow it is **isodesmic** — which every assembly edge in this
project fails to be. These are the first numbers here that `strict=True` accepts.
"""
from __future__ import annotations

import re
from typing import Any

from mofsbu._types import Fidelity, MofsbuError
from mofsbu.energy.reference import put_balanced_reaction

__all__ = ["deprotonation_pairs", "link_protomers", "PROTON_DONOR", "PROTON_ACCEPTOR"]

#: The couple, by the `TypedGraph.name` `examples` gives them.
PROTON_DONOR = "water"          # takes the proton away: H2O -> H3O+
PROTON_ACCEPTOR = "hydronium"

_FORMULA = re.compile(r"([A-Z][a-z]?)(\d*)")


def _counts(formula: str) -> dict[str, int]:
    """Element counts from a Hill string, for NARROWING ONLY — the graph decides."""
    out: dict[str, int] = {}
    for element, n in _FORMULA.findall(formula or ""):
        if element:
            out[element] = out.get(element, 0) + (int(n) if n else 1)
    return out


def deprotonation_pairs(reg: Any) -> list[tuple[int, int, int]]:
    """`(protonated_id, deprotonated_id, n)` for every pair differing by n protons.

    Candidates only: `put_balanced_reaction` weighs each one on the typed graphs and
    refuses what does not balance, so nothing here is trusted beyond picking rows to try.
    """
    rows = [r for r in reg.conn.execute(
        "SELECT id, formula, net_charge FROM structures WHERE hidden = 0")]
    by_key: dict[tuple[str, int], list[int]] = {}
    for r in rows:
        counts = _counts(r["formula"])
        # Removing a proton takes one H and one unit of charge together, so
        # `charge - H` is what two protomers of the same molecule share.
        key = ("".join(f"{e}{counts[e]}" for e in sorted(counts) if e != "H"),
               int(r["net_charge"] or 0) - counts.get("H", 0))
        by_key.setdefault(key, []).append(int(r["id"]))

    out: list[tuple[int, int, int]] = []
    counts_by_id = {int(r["id"]): _counts(r["formula"]) for r in rows}
    charge_by_id = {int(r["id"]): int(r["net_charge"] or 0) for r in rows}
    for group in by_key.values():
        # Within a group the heavy-atom skeleton and the proton-corrected charge match, so
        # any two members differ only in how many H they carry.
        for a in group:
            for b in group:
                n = counts_by_id[a].get("H", 0) - counts_by_id[b].get("H", 0)
                if n > 0 and charge_by_id[a] - charge_by_id[b] == n:
                    out.append((a, b, n))
    return out


def link_protomers(reg: Any, *, water_id: int, hydronium_id: int,
                   limit: int | None = None) -> list[dict[str, Any]]:
    """Write one balanced deprotonation edge per pair, skipping what already exists.

    Returns a row per pair with the reaction id, or the reason it was refused. A pair
    whose equation does not balance on the graphs is reported, not raised: the point is
    to link what can be linked, and to say what could not.
    """
    existing = _already_linked(reg)
    # The couple is itself a protomer pair, and linking it through itself states
    # `H3O+ + H2O -> H2O + H3O+` — balanced, isodesmic, and exactly zero. True and empty.
    couple = {int(water_id), int(hydronium_id)}
    written: list[dict[str, Any]] = []
    for protonated, deprotonated, n in deprotonation_pairs(reg):
        if (protonated, deprotonated) in existing or {protonated, deprotonated} & couple:
            continue
        if limit is not None and len(written) >= limit:
            break
        entry: dict[str, Any] = {"protonated": protonated,
                                 "deprotonated": deprotonated, "n": n,
                                 "reaction_id": None, "why_not": None}
        try:
            entry["reaction_id"] = put_balanced_reaction(
                reg, deprotonated,
                reagents=[(protonated, 1), (water_id, n)],
                leaving=[(hydronium_id, n)],
                kind="deprotonation",
                note=f"deprotonation x{n} (H2O/H3O+ couple)")
            existing.add((protonated, deprotonated))
        except MofsbuError as exc:
            entry["why_not"] = str(exc).splitlines()[0]
        written.append(entry)
    return written


def _already_linked(reg: Any) -> set[tuple[int, int]]:
    """Pairs a deprotonation edge already connects, so re-running writes nothing new."""
    out: set[tuple[int, int]] = set()
    for row in reg.conn.execute(
            "SELECT id, product_structure_id FROM reactions WHERE kind = 'deprotonation'"):
        for r in reg.conn.execute(
                "SELECT structure_id FROM reaction_reagents "
                "WHERE reaction_id = ? AND role = 'reagent'", (row["id"],)):
            out.add((int(r["structure_id"]), int(row["product_structure_id"])))
    return out
