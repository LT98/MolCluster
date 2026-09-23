"""What one provenance edge costs, and why it sometimes cannot be said.

`energy.reference` answers "what is dE for this equation, or which rule stops me" by
raising. A reader — the viewer's provenance panel — needs the same answer as *data*, one
entry per edge, because a panel listing twelve routes cannot be a single exception.

Two things this module is careful about:

* **A refusal is a result.** `check_balance` and `check_reference_quality` return rather
  than raise, so the diagnosis is assembled before any number is attempted and survives
  the attempt failing. An edge that cannot be priced still says what is wrong with it.
* **The rung is checked after the fact, not asked for.** `reference._energy_row` matches
  a requested fidelity EXACTLY, so asking for `Fidelity.ML` would refuse an xTB or DFT
  energy — the opposite of a floor. The equation is resolved at whatever rung it has and
  `min_fidelity` is applied to the result, which is `min()` over the terms.

Shaped as a list of steps from the start: a route is one step today, and the caller's
contract does not change when `pathways` learns to compose several (workplan M7 S3/S4).
"""
from __future__ import annotations

import re
from typing import Any

from mofsbu._types import Fidelity, MofsbuError, ReferenceSchemeError
from mofsbu.energy.reference import (
    Term, balance_for_terms, check_balance, check_reference_quality, energy_for_terms,
    quality_for_terms, reaction_balanced_energy, reaction_terms,
)
from mofsbu.identity import hill_formula

__all__ = ["price_reaction", "price_incoming_routes", "decompositions"]


#: A refusal this module reports rather than raises.  `OSError` is in here because a
#: graph is read from the blob store, and a row whose blob is gone is a condition the
#: reader already knows how to show (`ui/app.py:_blob_ok`) — not a reason to fail a page.
REPORTABLE = (MofsbuError, OSError)


def _terms(reg: Any, reaction_id: int) -> list[dict[str, Any]]:
    """Every species in the edge. `display_label` is `None` when the row carries none —
    a reader shows the id rather than inventing a name."""
    try:
        terms = reaction_terms(reg, reaction_id)
        labels = _display_labels(reg, [t.structure_id for t in terms])
        return [{"structure_id": t.structure_id, "stoich": t.stoich,
                 "role": t.role, "side": t.side, "label": t.label,
                 "display_label": labels.get(t.structure_id)}
                for t in terms]
    except REPORTABLE:
        return []


def _display_labels(reg: Any, structure_ids: list[int]) -> dict[int, str | None]:
    """One query for the whole edge — a panel prices dozens of edges per page load."""
    if not structure_ids:
        return {}
    marks = ",".join("?" * len(structure_ids))
    rows = reg.conn.execute(
        f"SELECT id, display_label FROM structures WHERE id IN ({marks})",  # noqa: S608
        structure_ids).fetchall()
    return {int(r["id"]): r["display_label"] for r in rows}


def _step(reaction_id: int | None) -> dict[str, Any]:
    """An unpriced step. Every key the priced shape carries, so the caller never branches."""
    return {"reaction_id": reaction_id, "dE": None, "fidelity": None,
            "fidelity_name": None, "method": None, "isodesmic": None, "caveats": [],
            "all_converged": None, "terms": [], "can_price": False, "why_not": None}


def price_reaction(reg: Any, reaction_id: int, *,
                   min_fidelity: Fidelity = Fidelity.ML,
                   solvent: str | None = None) -> dict[str, Any]:
    """dE for one edge, or the stated reason there is none.

    `reg` needs only `.conn` and `.store` — `registry.ReadOnlyRegistry` is enough, and is
    what lets a read-only viewer call this without opening a writable connection.
    """
    row = reg.conn.execute(
        "SELECT id, kind, depth, note FROM reactions WHERE id = ?", (reaction_id,)).fetchone()
    if row is None:
        raise ReferenceSchemeError(f"no reaction {reaction_id}")

    step = _step(reaction_id)
    step["terms"] = _terms(reg, reaction_id)

    # The diagnosis first, so it is attached whether or not a number follows.
    try:
        quality = check_reference_quality(reg, reaction_id, medium=solvent)
        step["isodesmic"] = quality.isodesmic
        step["caveats"] = [i.code for i in quality.issues]
    except REPORTABLE:
        pass

    try:
        balance = check_balance(reg, reaction_id)
        if not balance.balanced:
            step["why_not"] = balance.describe()
        else:
            energy = reaction_balanced_energy(
                reg, reaction_id, fidelity=None, solvent=solvent, strict=False)
            if energy.fidelity < min_fidelity:
                step["why_not"] = (
                    f"best available is {energy.fidelity.name}; needs "
                    f"{min_fidelity.name} or better")
            else:
                step.update(
                    dE=energy.dE, fidelity=int(energy.fidelity),
                    fidelity_name=energy.fidelity.name,
                    method=energy.method.describe(),
                    isodesmic=energy.quality.isodesmic,
                    caveats=[i.code for i in energy.quality.issues],
                    all_converged=energy.all_converged, can_price=True)
    except REPORTABLE as exc:
        step["why_not"] = str(exc)

    return {"reaction_id": reaction_id, "kind": row["kind"], "depth": row["depth"],
            "note": row["note"] or "", "steps": [step],
            "can_price": step["can_price"],
            "why_not": step["why_not"],
            "total_dE": step["dE"] if step["can_price"] else None}


def price_incoming_routes(reg: Any, structure_id: int, *,
                          min_fidelity: Fidelity = Fidelity.ML,
                          solvent: str | None = None) -> list[dict[str, Any]]:
    """Every edge arriving at this structure, priced or refused, in `reactions` order."""
    ids = [r["id"] for r in reg.conn.execute(
        "SELECT id FROM reactions WHERE product_structure_id = ? ORDER BY id",
        (structure_id,))]
    return [price_reaction(reg, rid, min_fidelity=min_fidelity, solvent=solvent)
            for rid in ids]


# ── what this structure could have been made from ────────────────────────────
#
# A `place` edge records a construction, not a reaction: the runner built the whole
# coordination sphere at once, so the edge has no reagents and its balance report is the
# product's entire composition.  That leaves a real question unanswered — what pieces
# ALREADY IN THE REGISTRY add up to this thing — and the answer is arithmetic, not a
# record.  It is offered as `inferred` and never mixed with what was recorded.

_FORMULA = re.compile(r"([A-Z][a-z]?)(\d*)")


def _counts(formula: str) -> dict[str, int]:
    """Element counts from a Hill string, for NARROWING ONLY.

    `check_balance` is explicit that a formula parsed back into counts is a second
    encoding waiting to disagree with the graph, so nothing here is trusted: this picks
    candidates out of an indexed column, and `balance_for_terms` then weighs every
    candidate on the typed graphs before it is reported.
    """
    out: dict[str, int] = {}
    for element, n in _FORMULA.findall(formula or ""):
        if element:
            out[element] = out.get(element, 0) + (int(n) if n else 1)
    return out


def decompositions(reg: Any, structure_id: int, *,
                   min_fidelity: Fidelity = Fidelity.ML,
                   solvent: str | None = None,
                   limit: int = 25) -> list[dict[str, Any]]:
    """Pairs of registry structures whose atoms and charge add up to this one.

    The complement of `price_incoming_routes`: that says how this structure WAS reached,
    this says how it COULD be, from species that already exist.  Two-part splits only —
    the ladder that would compose more is M8's.
    """
    target = reg.conn.execute(
        "SELECT id, formula, net_charge FROM structures WHERE id = ?",
        (structure_id,)).fetchone()
    if target is None:
        raise ReferenceSchemeError(f"no structure {structure_id}")
    want = _counts(target["formula"])
    want_q = int(target["net_charge"] or 0)

    rows = [r for r in reg.conn.execute(
        "SELECT id, formula, net_charge, display_label FROM structures "
        "WHERE hidden = 0 AND id != ?", (structure_id,))]
    index: dict[tuple[str, int], list[Any]] = {}
    for r in rows:
        index.setdefault((r["formula"] or "", int(r["net_charge"] or 0)), []).append(r)

    recorded = _recorded_pairs(reg, structure_id)
    seen: set[tuple[int, int]] = set()
    out: list[dict[str, Any]] = []
    for a in rows:
        rest = dict(want)
        for element, n in _counts(a["formula"]).items():
            rest[element] = rest.get(element, 0) - n
        if any(n < 0 for n in rest.values()):
            continue
        rest = {e: n for e, n in rest.items() if n}
        if not rest:                      # the whole product: not a split
            continue
        for b in index.get((hill_formula(rest), want_q - int(a["net_charge"] or 0)), ()):
            key = (min(a["id"], b["id"]), max(a["id"], b["id"]))
            if key in seen:
                continue
            seen.add(key)
            out.append(_derived(reg, target, a, b, recorded.get(key),
                                min_fidelity=min_fidelity, solvent=solvent))
            if len(out) >= limit:
                return out
    return out


def _recorded_pairs(reg: Any, structure_id: int) -> dict[tuple[int, int], int]:
    """Which two-reagent splits are already written down, so a derived one can say so."""
    out: dict[tuple[int, int], int] = {}
    for row in reg.conn.execute(
            "SELECT id FROM reactions WHERE product_structure_id = ?", (structure_id,)):
        ids = [r["structure_id"] for r in reg.conn.execute(
            "SELECT structure_id FROM reaction_reagents WHERE reaction_id = ?",
            (row["id"],))]
        if len(ids) == 2:
            out.setdefault((min(ids), max(ids)), int(row["id"]))
    return out


def _derived(reg: Any, target: Any, a: Any, b: Any, recorded: int | None, *,
             min_fidelity: Fidelity, solvent: str | None) -> dict[str, Any]:
    terms = (Term(structure_id=int(target["id"]), stoich=1, role="product",
                  side="product"),
             Term(structure_id=int(a["id"]), stoich=1, role="reagent", side="reagent"),
             Term(structure_id=int(b["id"]), stoich=1, role="reagent", side="reagent"))
    step = _step(None)                    # derived: there is no row to point at
    step["terms"] = [{"structure_id": t.structure_id, "stoich": t.stoich,
                      "role": t.role, "side": t.side, "label": t.label} for t in terms]
    subject = f"{a['display_label']} + {b['display_label']}"
    try:
        quality = quality_for_terms(reg, terms, medium=solvent)
        step["isodesmic"] = quality.isodesmic
        step["caveats"] = [i.code for i in quality.issues]
    except REPORTABLE:
        pass
    try:
        balance = balance_for_terms(reg, terms)
        if not balance.balanced:
            # The formula index proposed it; the graphs are what decide.
            step["why_not"] = balance.describe()
        else:
            energy = energy_for_terms(reg, terms, fidelity=None, solvent=solvent,
                                      strict=False, subject=subject)
            if energy.fidelity < min_fidelity:
                step["why_not"] = (f"best available is {energy.fidelity.name}; needs "
                                   f"{min_fidelity.name} or better")
            else:
                step.update(dE=energy.dE, fidelity=int(energy.fidelity),
                            fidelity_name=energy.fidelity.name,
                            method=energy.method.describe(),
                            isodesmic=energy.quality.isodesmic,
                            caveats=[i.code for i in energy.quality.issues],
                            all_converged=energy.all_converged, can_price=True)
    except REPORTABLE as exc:
        step["why_not"] = str(exc)

    return {"origin": "inferred", "recorded_as": recorded,
            "parts": [{"structure_id": int(a["id"]), "stoich": 1,
                       "label": a["display_label"]},
                      {"structure_id": int(b["id"]), "stoich": 1,
                       "label": b["display_label"]}],
            "steps": [step], "can_price": step["can_price"],
            "why_not": step["why_not"],
            "total_dE": step["dE"] if step["can_price"] else None}
