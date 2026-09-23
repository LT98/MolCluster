"""What a whole route costs, with the target at zero.

`energy.routes` prices one edge. A reader comparing routes needs them composed: a
number at every node, on an axis both routes share. That is arithmetic over stored
energies and it is done here rather than in the browser, because nothing in
`scripts/check.sh` can execute a line of JavaScript and a cumulative sum no test can
reach is exactly the kind of number this project refuses to show.

**The anchor is the target.** `y(target) = 0` and `y` decreases leftwards by each step's
dE, so a route's precursors sit above a target they fall to. Every route ends at the one
node they all share, which is what pins them to a common point; their left ends are
wherever each one bottoms out, and those are not required to agree.

**What a y value is an energy OF.** Each step balances, so at every node

    comp(node) + comp(pieces not yet consumed) == comp(target) + comp(what has left)

— the node plus its spectators has exactly the target's atoms. `y` is that system's
energy relative to the target's, which is why it can be read down a column.

**When two routes may be compared.** Only where the accounting matches: same target, and
the same things shed on the way to it. The identity above then forces the unconsumed
pieces to match too, so the shed multiset alone decides it, and it is reported per node
as `basis`. An assembly step sheds nothing and a deprotonation sheds `n H3O+`, so what
this actually catches is a protonated route set beside a deprotonated one — where the
vertical gap between the curves is not a comparison of anything.

A step is `recorded` (a `reactions` row) or `inferred` (a split `decompositions` derived
from composition). The two are never merged: the distinction is what the whole reference
scheme rests on.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from mofsbu._types import Fidelity, ReferenceSchemeError
from mofsbu.energy.routes import decompositions, price_reaction

__all__ = ["price_path", "DERIVED"]

#: A leg walked along a split nobody recorded: `"d"` followed by the id of the OTHER
#: part.  Never a reaction id, and never stored.  Both parts are named because one is
#: not enough to identify the split: the complement is fixed by composition only up to
#: its (formula, charge), and the registry holds several protomers that share those —
#: measured, on structure 41, where one part appears in four different splits.
DERIVED = "d"

#: How far down the derived splits to look for the one a leg names.  `decompositions`
#: stops early by default because a panel lists them; resolving a named leg must not.
_DERIVED_SEARCH = 200


def price_path(reg: Any, nodes: Sequence[int], via: Sequence[int | str], *,
               min_fidelity: Fidelity = Fidelity.ML,
               solvent: str | None = None) -> dict[str, Any]:
    """Price a chain of steps read backwards from `nodes[0]`, the target.

    `nodes` is the walk: the target, then what each step was reached from. `via[i]` is the
    edge joining `nodes[i+1]` to `nodes[i]` — a reaction id, or `DERIVED` for a split that
    was derived rather than recorded.

    Raises for a malformed walk (lengths that disagree, a node that is not there, an edge
    that does not join the two nodes it is given for). Everything a reference-scheme rule
    refuses comes back as data, per step and on the path, because a route drawn on a chart
    cannot be a single exception any more than a panel of them could.
    """
    nodes = [int(n) for n in nodes]
    if not nodes:
        raise ReferenceSchemeError("a path needs at least the structure it ends at")
    if len(via) != len(nodes) - 1:
        raise ReferenceSchemeError(
            f"{len(nodes)} nodes need {len(nodes) - 1} edges between them, got {len(via)}")
    labels = _labels(reg, nodes)
    missing = [n for n in nodes if n not in labels]
    if missing:
        raise ReferenceSchemeError(f"no structure {missing[0]}")

    steps = [_price_leg(reg, nodes[i], nodes[i + 1], via[i],
                        min_fidelity=min_fidelity, solvent=solvent)
             for i in range(len(via))]

    # Walk out from the target accumulating both sides of the identity above.  `y` stops
    # at the first step that has no number and stays stopped: a gap is not a zero (D18).
    out_nodes = [{"structure_id": nodes[0], "display_label": labels[nodes[0]],
                  "y": 0.0, "shed": [], "free_pieces": [], "basis": _basis([])}]
    shed: dict[int, int] = {}
    free: dict[int, int] = {}
    y: float | None = 0.0
    for i, step in enumerate(steps):
        for t in step["shed"]:
            shed[t["structure_id"]] = shed.get(t["structure_id"], 0) + t["stoich"]
        for t in step["added"] + step["solvent"]:
            free[t["structure_id"]] = free.get(t["structure_id"], 0) + t["stoich"]
        y = None if (y is None or not step["can_price"]) else y - step["dE"]
        out_nodes.append({
            "structure_id": nodes[i + 1], "display_label": labels[nodes[i + 1]],
            "y": y,
            "shed": _multiset(reg, shed), "free_pieces": _multiset(reg, free),
            "basis": _basis(sorted(shed.items()))})

    priced = [s for s in steps if s["can_price"]]
    unpriced = [s for s in steps if not s["can_price"]]
    caveats: list[str] = []
    for s in steps:
        for c in s["caveats"]:
            if c not in caveats:
                caveats.append(c)
    # A route is only as good as its worst number, so the rung it reports is the lowest
    # one on it — not the best, which would describe a step rather than the path.
    rungs = [s["fidelity"] for s in priced if s["fidelity"] is not None]
    worst = max(priced, key=lambda s: s["dE"], default=None)
    return {
        "nodes": out_nodes, "steps": steps,
        "can_price": bool(steps) and not unpriced,
        "why_not": unpriced[0]["why_not"] if unpriced else None,
        "n_unpriced": len(unpriced),
        "total_dE": sum(s["dE"] for s in priced) if steps and not unpriced else None,
        "basis": out_nodes[-1]["basis"],
        "caveats": caveats,
        "isodesmic": all(s["isodesmic"] is True for s in steps) if steps else None,
        "fidelity": min(rungs) if rungs else None,
        "fidelity_name": Fidelity(min(rungs)).name if rungs else None,
        "worst_step": None if worst is None else
                      {"index": steps.index(worst), "dE": worst["dE"]},
    }


def _price_leg(reg: Any, product_id: int, source_id: int, via: int | str, *,
               min_fidelity: Fidelity, solvent: str | None) -> dict[str, Any]:
    """One step, flattened: the number, the diagnosis, and what moved either way."""
    if str(via).startswith(DERIVED):
        other = str(via)[len(DERIVED):]
        priced = _derived_leg(reg, product_id, source_id,
                              int(other) if other else None,
                              min_fidelity=min_fidelity, solvent=solvent)
        origin, reaction_id = "inferred", None
        recorded_as = priced.get("recorded_as")
    else:
        priced = price_reaction(reg, int(via), min_fidelity=min_fidelity, solvent=solvent)
        origin, reaction_id, recorded_as = "recorded", int(via), None

    step = priced["steps"][0]
    terms = step["terms"]
    # The edge has to be the one that joins these two nodes, or the sum is of a different
    # route than the one being drawn.  Checked rather than assumed: a link can be stale.
    if origin == "recorded":
        _check_joins(reg, int(via), product_id, source_id, terms)

    # `added` and `solvent` are both consumed, so both count in the arithmetic — but they
    # are not the same claim.  Water carrying a proton away is not a piece being built in,
    # and a step that listed it as "added" would read as chemistry it is not doing.
    added = [_term(t) for t in terms
             if t["side"] == "reagent" and (t["role"] or "reagent") == "reagent"
             and int(t["structure_id"]) != source_id]
    solvent = [_term(t) for t in terms
               if t["side"] == "reagent" and t["role"] == "solvent"]
    return {
        "via": reaction_id if reaction_id is not None else DERIVED,
        "origin": origin, "reaction_id": reaction_id, "recorded_as": recorded_as,
        "product": product_id, "source": source_id,
        "dE": priced["total_dE"], "can_price": priced["can_price"],
        "why_not": priced["why_not"],
        "fidelity": step["fidelity"], "fidelity_name": step["fidelity_name"],
        "isodesmic": step["isodesmic"], "caveats": list(step["caveats"]),
        "all_converged": step["all_converged"],
        "added": added, "solvent": solvent,
        "shed": [_term(t) for t in terms if t["role"] == "leaving"],
        "terms": terms,
    }


def _check_joins(reg: Any, reaction_id: int, product_id: int, source_id: int,
                 terms: list[dict[str, Any]]) -> None:
    row = reg.conn.execute("SELECT product_structure_id FROM reactions WHERE id = ?",
                           (reaction_id,)).fetchone()
    if row is None or int(row["product_structure_id"]) != product_id:
        raise ReferenceSchemeError(
            f"reaction {reaction_id} does not produce structure {product_id}")
    if not any(t["side"] == "reagent" and int(t["structure_id"]) == source_id
               for t in terms):
        raise ReferenceSchemeError(
            f"reaction {reaction_id} does not consume structure {source_id}")


def _derived_leg(reg: Any, product_id: int, source_id: int, other_id: int | None, *,
                 min_fidelity: Fidelity, solvent: str | None) -> dict[str, Any]:
    """The split of `product_id` into `source_id` and `other_id`, re-derived.

    Nothing derived is ever written down, so the split is recovered from the two parts
    the leg names rather than looked up.  `other_id` may be omitted, and then this
    refuses an ambiguous leg instead of choosing: several splits can share a part, and
    picking one silently would put a number on the chart that the link did not ask for.
    """
    found = [d for d in decompositions(reg, product_id, min_fidelity=min_fidelity,
                                       solvent=solvent, limit=_DERIVED_SEARCH)
             if any(int(p["structure_id"]) == source_id for p in d["parts"])
             and (other_id is None
                  or any(int(p["structure_id"]) == other_id for p in d["parts"]))]
    if not found:
        raise ReferenceSchemeError(
            f"structure {product_id} does not split into {source_id}"
            + (f" and {other_id}" if other_id is not None else ""))
    if len(found) > 1:
        others = sorted({int(p["structure_id"]) for d in found for p in d["parts"]
                         if int(p["structure_id"]) != source_id})
        raise ReferenceSchemeError(
            f"structure {product_id} splits into {source_id} in {len(found)} ways "
            f"(with {', '.join(map(str, others))}); say which")
    return found[0]


def _term(t: dict[str, Any]) -> dict[str, Any]:
    return {"structure_id": int(t["structure_id"]), "stoich": int(t["stoich"] or 1),
            "display_label": t.get("display_label"), "role": t["role"]}


def _multiset(reg: Any, counts: dict[int, int]) -> list[dict[str, Any]]:
    if not counts:
        return []
    labels = _labels(reg, list(counts))
    return [{"structure_id": sid, "stoich": n, "display_label": labels.get(sid)}
            for sid, n in sorted(counts.items())]


def _basis(items: Sequence[tuple[int, int]]) -> str:
    """What the y values at this node are measured against, as one comparable token.

    Two nodes with the same basis carry energies of systems with the same atoms, so the
    gap between them is a number about chemistry.  With different bases it is not, and
    the page says so rather than letting it be read.
    """
    return ",".join(f"{sid}x{n}" for sid, n in items if n) or "—"


def _labels(reg: Any, ids: Sequence[int]) -> dict[int, str | None]:
    ids = list(dict.fromkeys(int(i) for i in ids))
    if not ids:
        return {}
    marks = ",".join("?" * len(ids))
    return {int(r["id"]): r["display_label"] for r in reg.conn.execute(
        f"SELECT id, display_label FROM structures WHERE id IN ({marks})", ids)}  # noqa: S608
