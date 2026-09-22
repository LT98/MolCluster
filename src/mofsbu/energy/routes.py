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

from typing import Any

from mofsbu._types import Fidelity, MofsbuError, ReferenceSchemeError
from mofsbu.energy.reference import (
    check_balance, check_reference_quality, reaction_balanced_energy, reaction_terms,
)

__all__ = ["price_reaction", "price_incoming_routes"]


#: A refusal this module reports rather than raises.  `OSError` is in here because a
#: graph is read from the blob store, and a row whose blob is gone is a condition the
#: reader already knows how to show (`ui/app.py:_blob_ok`) — not a reason to fail a page.
REPORTABLE = (MofsbuError, OSError)


def _terms(reg: Any, reaction_id: int) -> list[dict[str, Any]]:
    try:
        return [{"structure_id": t.structure_id, "stoich": t.stoich,
                 "role": t.role, "side": t.side, "label": t.label}
                for t in reaction_terms(reg, reaction_id)]
    except REPORTABLE:
        return []


def _step(reaction_id: int) -> dict[str, Any]:
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
        quality = check_reference_quality(reg, reaction_id)
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
