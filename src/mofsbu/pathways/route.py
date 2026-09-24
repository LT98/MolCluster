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

**Two directions.** A leg is walked one of two ways. *made from* follows an edge against
its arrow, from its product to a reagent. *consumed by* follows it with the arrow, from a
reagent to its product — so in the route's own direction (towards the target) the edge
runs in reverse: it contributes `-dE`, what the edge added is released, and what it shed
is taken up. Walking back to a shared parent and forward again composes a ligand
exchange out of recorded assemblies.

The spectators are one signed tally, so a piece shed on one leg and taken up on another
cancels. A route's reference quality is judged on its **net equation** — every leg's
terms summed in the route's direction, with whatever appears on both sides cancelled — and
not on its steps, whose own diagnostics stay on them: an exchange can be isodesmic as a
whole when every leg of it changes the coordination count.

A **pivot** is a node where the walk changes direction. It is bookkeeping for the
composition, not an intermediate the route claims to pass through, and its `y` is not a
barrier (docs/gui/05_graph.md).

A step is `recorded` (a `reactions` row) or `inferred` (a split `decompositions` derived
from composition). The two are never merged: the distinction is what the whole reference
scheme rests on. A ROUTE is `recorded` only when its net equation is one row as written,
`inferred` when any leg is, and otherwise `composed`, with the rows behind it as `witness`.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from mofsbu._types import Fidelity, ReferenceSchemeError
from mofsbu.energy.reference import (
    Term, balance_for_terms, energy_for_terms, quality_for_terms,
)
from mofsbu.energy.routes import REPORTABLE, decompositions, price_reaction

__all__ = ["price_path", "proton_sinks", "proton_couple", "DERIVED", "CONSUMED", "MADE_FROM",
           "CONSUMED_BY"]

#: A leg walked along a split nobody recorded: `"d"` followed by the id of the OTHER
#: part.  Never a reaction id, and never stored.  Both parts are named because one is
#: not enough to identify the split: the complement is fixed by composition only up to
#: its (formula, charge), and the registry holds several protomers that share those —
#: measured, on structure 41, where one part appears in four different splits.
DERIVED = "d"

#: A leg walked WITH a recorded edge's arrow: `"c"` followed by the reaction id.  The node
#: the walk stands on is a `reagent` of that reaction and the next node is its product.
#: A bare reaction id keeps meaning "made from", so every existing link reads as before.
CONSUMED = "c"

MADE_FROM = "made_from"
CONSUMED_BY = "consumed_by"

#: How far down the derived splits to look for the one a leg names.  `decompositions`
#: stops early by default because a panel lists them; resolving a named leg must not.
_DERIVED_SEARCH = 200

#: The legs summed and the net equation priced whole must agree to this, or they drew
#: on different energies for one species.
_AGREE_EV = 1e-6


def price_path(reg: Any, nodes: Sequence[int], via: Sequence[int | str], *,
               min_fidelity: Fidelity = Fidelity.ML,
               solvent: str | None = None,
               proton_sink: int | None = None) -> dict[str, Any]:
    """Price a chain of steps read backwards from `nodes[0]`, the target.

    `nodes` is the walk: the target, then each node the walk stepped to. `via[i]` is the
    edge joining `nodes[i]` and `nodes[i+1]`: a reaction id (made from — the edge produces
    `nodes[i]` from `nodes[i+1]`), `c<reaction id>` (consumed by — the edge produces
    `nodes[i+1]` from `nodes[i]`), or `d<other part>` for a split that was derived.

    Raises for a malformed walk (lengths that disagree, a node that is not there, an edge
    that does not join the two nodes it is given for, a leg that walks straight back along
    the one before it). Everything a reference-scheme rule refuses comes back as data, per
    step and on the path, because a route drawn on a chart cannot be a single exception
    any more than a panel of them could.

    `proton_sink` names a free base A⁻ (one of `proton_sinks(reg)`) that takes each proton
    the route releases instead of water: every step that sheds n H₃O⁺ also runs
    n × (H₃O⁺ + A⁻ → H₂O + HA), priced in the same medium.  Exact by Hess's law, and the
    step, the spectators and the net equation all carry it, so nothing is hidden.
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
    for i in range(1, len(steps)):
        a, b = steps[i - 1], steps[i]
        if (a["reaction_id"] is not None and a["reaction_id"] == b["reaction_id"]
                and nodes[i + 1] == nodes[i - 1]):
            raise ReferenceSchemeError(
                f"leg {i + 1} walks reaction {b['reaction_id']} straight back to structure "
                f"{nodes[i + 1]}; that is stepping back along the trail, not a hop")
    sink = None
    if proton_sink is not None:
        sink = _apply_proton_sink(reg, steps, int(proton_sink), solvent=solvent)

    # Walk out from the target with one signed tally of spectators: positive is a piece
    # not yet consumed, negative is something that has left.  `y` stops at the first
    # step that has no number and stays stopped: a gap is not a zero (D18).
    out_nodes = [_node(reg, nodes[0], labels, 0.0, {})]
    spectators: dict[int, int] = {}
    net: dict[int, int] = {}
    y: float | None = 0.0
    for i, step in enumerate(steps):
        for sid, n in _route_counts(step).items():
            net[sid] = net.get(sid, 0) + n
            spectators[sid] = spectators.get(sid, 0) - n
        # The two chain nodes are the walk, not spectators.
        spectators[nodes[i + 1]] = spectators.get(nodes[i + 1], 0) - 1
        spectators[nodes[i]] = spectators.get(nodes[i], 0) + 1
        y = None if (y is None or not step["can_price"]) else y - step["dE"]
        out_nodes.append(_node(reg, nodes[i + 1], labels, y, spectators))

    pivots = []
    for i in range(1, len(steps)):
        if steps[i - 1]["direction"] != steps[i]["direction"]:
            kind = "parent" if steps[i - 1]["direction"] == MADE_FROM else "product"
            out_nodes[i].update(pivot=True, pivot_kind=kind)
            pivots.append(i)

    priced = [s for s in steps if s["can_price"]]
    unpriced = [s for s in steps if not s["can_price"]]
    total = sum(s["dE"] for s in priced) if steps and not unpriced else None
    equation = (_net_equation(reg, nodes[0], net, steps, total,
                              min_fidelity=min_fidelity, solvent=solvent)
                if steps else None)
    # A route is only as good as its worst number, so the rung it reports is the lowest
    # one on it — not the best, which would describe a step rather than the path.
    rungs = [s["fidelity"] for s in priced if s["fidelity"] is not None]
    worst = max(priced, key=lambda s: s["dE"], default=None)
    worst_why = None
    if pivots:
        worst, worst_why = None, (
            "the route changes direction at a pivot, so the steps into and out of it are "
            "halves of one exchange and neither alone is a step to single out")
    if not steps:
        origin = None
    elif any(s["origin"] == "inferred" for s in steps):
        origin = "inferred"
    elif len(steps) == 1 and steps[0]["direction"] == MADE_FROM:
        origin = "recorded"
    else:
        origin = "composed"
    return {
        "nodes": out_nodes, "steps": steps,
        "can_price": bool(steps) and not unpriced,
        "why_not": unpriced[0]["why_not"] if unpriced else None,
        "n_unpriced": len(unpriced),
        "total_dE": total,
        "basis": out_nodes[-1]["basis"],
        "net_equation": equation,
        "proton_sink": sink,
        "medium": solvent,
        "caveats": list(equation["caveats"]) if equation else [],
        "isodesmic": equation["isodesmic"] if equation else None,
        "origin": origin,
        "witness": [{"reaction_id": s["reaction_id"], "direction": s["direction"]}
                    for s in steps if s["reaction_id"] is not None],
        "pivots": pivots,
        "fidelity": min(rungs) if rungs else None,
        "fidelity_name": Fidelity(min(rungs)).name if rungs else None,
        "worst_step": None if worst is None else
                      {"index": steps.index(worst), "dE": worst["dE"]},
        "worst_step_why": worst_why,
    }


def _node(reg: Any, sid: int, labels: dict[int, str | None], y: float | None,
          spectators: dict[int, int]) -> dict[str, Any]:
    shed = {k: -n for k, n in spectators.items() if n < 0}
    free = {k: n for k, n in spectators.items() if n > 0}
    return {"structure_id": sid, "display_label": labels[sid], "y": y,
            "shed": _multiset(reg, shed), "free_pieces": _multiset(reg, free),
            "basis": _basis(sorted(shed.items())),
            "pivot": False, "pivot_kind": None}


def proton_couple(reg: Any) -> tuple[int, int] | None:
    """`(water_id, hydronium_id)` as the recorded deprotonation edges use them, or None."""
    row = reg.conn.execute(
        "SELECT MAX(CASE WHEN rr.role = 'solvent' THEN rr.structure_id END) AS water, "
        "       MAX(CASE WHEN rr.role = 'leaving' THEN rr.structure_id END) AS hydronium "
        "FROM reactions r JOIN reaction_reagents rr ON rr.reaction_id = r.id "
        "WHERE r.kind = 'deprotonation'").fetchone()
    if row is None or row["water"] is None or row["hydronium"] is None:
        return None
    return int(row["water"]), int(row["hydronium"])


def proton_sinks(reg: Any) -> list[dict[str, Any]]:
    """Free bases a released proton can be handed to instead of water.

    Every metal-free deprotonation edge whose acid is metal-free too — so only a base the
    registry actually holds, with its conjugate acid, at an energy the run computed.  For
    a chloride run that is the ligand's own anion and nothing else: HCl is not a species.
    """
    couple = proton_couple(reg)
    if couple is None:
        return []
    rows = reg.conn.execute(
        "SELECT DISTINCT r.product_structure_id AS base, rr.structure_id AS acid "
        "FROM reactions r JOIN reaction_reagents rr ON rr.reaction_id = r.id "
        "WHERE r.kind = 'deprotonation' AND rr.role = 'reagent' "
        "AND NOT EXISTS (SELECT 1 FROM structure_metals sm "
        "                WHERE sm.structure_id IN (r.product_structure_id, rr.structure_id)) "
        "ORDER BY base").fetchall()
    labels = _labels(reg, [int(r[k]) for r in rows for k in ("base", "acid")])
    return [{"base": int(r["base"]), "acid": int(r["acid"]),
             "base_label": labels.get(int(r["base"])), "acid_label": labels.get(int(r["acid"]))}
            for r in rows if int(r["base"]) not in couple]


def _apply_proton_sink(reg: Any, steps: list[dict[str, Any]], base: int, *,
                       solvent: str | None) -> dict[str, Any]:
    """Re-point every proton a route releases from water to `base`, in place.

    A step that releases n H₃O⁺ (route direction) also runs n × (H₃O⁺ + A⁻ → H₂O + HA):
    its terms gain that equation, so the spectator tally, the basis and the net equation
    follow, and its dE gains n × the conversion's dE.  A conversion that cannot be priced
    leaves every step that needed it unpriced, with the reason — never a silent water.
    """
    sinks = {s["base"]: s for s in proton_sinks(reg)}
    if base not in sinks:
        raise ReferenceSchemeError(
            f"structure {base} is not a proton sink in this registry; the free bases with a "
            f"conjugate acid are {sorted(sinks) or 'none'}")
    water, hydronium = proton_couple(reg)
    acid = sinks[base]["acid"]
    conversion = (Term(acid, 1, "product", "product"), Term(water, 1, "leaving", "product"),
                  Term(hydronium, 1, "reagent", "reagent"), Term(base, 1, "reagent", "reagent"))
    out: dict[str, Any] = {**sinks[base], "water": water, "hydronium": hydronium,
                           "dE_per_proton": None, "why_not": None, "protons": 0}
    try:
        out["dE_per_proton"] = energy_for_terms(
            reg, conversion, fidelity=None, solvent=solvent, strict=False,
            subject=f"H3O+ + {sinks[base]['base_label']} -> H2O + {sinks[base]['acid_label']}"
        ).dE
    except REPORTABLE as exc:
        out["why_not"] = str(exc).splitlines()[0]
    for step in steps:
        n = _route_counts(step).get(hydronium, 0)       # released (+) or taken up (-)
        if not n:
            continue
        out["protons"] += n
        wanted = {hydronium: -n, base: -n, water: n, acid: n}
        for sid, c in wanted.items():
            side = "product" if c * step["sign"] > 0 else "reagent"
            step["terms"].append({"structure_id": sid, "stoich": abs(c), "side": side,
                                  "role": "leaving" if side == "product" else "solvent",
                                  "display_label": None})
        step["proton_sink"] = {"base": base, "acid": acid, "protons": n}
        if out["dE_per_proton"] is None:
            step.update(can_price=False, dE=None,
                        why_not=f"the proton sink's conversion cannot be priced: "
                                f"{out['why_not']}")
        elif step["can_price"]:
            step["dE"] += n * out["dE_per_proton"]
    return out


def _route_counts(step: dict[str, Any]) -> dict[int, int]:
    """The leg's equation in the ROUTE's direction, as signed counts (products positive)."""
    out: dict[int, int] = {}
    for t in step["terms"]:
        n = step["sign"] * (1 if t["side"] == "product" else -1) * int(t["stoich"] or 1)
        out[int(t["structure_id"])] = out.get(int(t["structure_id"]), 0) + n
    return out


def _net_equation(reg: Any, target: int, net: dict[int, int],
                  steps: list[dict[str, Any]], total: float | None, *,
                  min_fidelity: Fidelity, solvent: str | None) -> dict[str, Any]:
    """The route as one equation, judged by the rules a single edge is judged by."""
    out: dict[str, Any] = {"terms": [], "balanced": None, "balance": None,
                           "isodesmic": None, "caveats": [], "dative_delta": None,
                           "dE": None, "fidelity": None, "fidelity_name": None,
                           "agrees_with_steps": None, "why_not": None}
    bare = [i for i, s in enumerate(steps) if not s["terms"]]
    if bare:
        out["why_not"] = (f"step {bare[0] + 1} has no species recorded, so the route has "
                          f"no net equation")
        return out
    kept = {sid: n for sid, n in net.items() if n}
    if not kept:
        out["why_not"] = "every species cancels: the walk returns to where it started"
        return out
    # The target first: `energy_for_terms` pins the theory to its first term, as a
    # stored reaction's product pins it.
    order = sorted(kept, key=lambda sid: (sid != target, kept[sid] < 0, sid))
    terms = tuple(
        Term(structure_id=sid, stoich=abs(kept[sid]),
             role=("product" if sid == target else "leaving") if kept[sid] > 0
             else "reagent",
             side="product" if kept[sid] > 0 else "reagent")
        for sid in order)
    labels = _labels(reg, order)
    out["terms"] = [{"structure_id": t.structure_id, "stoich": t.stoich, "role": t.role,
                     "side": t.side, "display_label": labels.get(t.structure_id)}
                    for t in terms]
    try:
        quality = quality_for_terms(reg, terms)
        out.update(isodesmic=quality.isodesmic, dative_delta=quality.dative_delta,
                   caveats=[i.code for i in quality.issues])
        balance = balance_for_terms(reg, terms)
        out.update(balanced=balance.balanced, balance=balance.describe())
        if not balance.balanced:
            out["why_not"] = balance.describe()
            return out
        energy = energy_for_terms(reg, terms, fidelity=None, solvent=solvent,
                                  strict=False, subject="the route's net equation")
        if energy.fidelity < min_fidelity:
            out["why_not"] = (f"best available is {energy.fidelity.name}; needs "
                              f"{min_fidelity.name} or better")
            return out
        out.update(dE=energy.dE, fidelity=int(energy.fidelity),
                   fidelity_name=energy.fidelity.name)
        if total is not None:
            out["agrees_with_steps"] = abs(energy.dE - total) <= _AGREE_EV
    except REPORTABLE as exc:
        out["why_not"] = str(exc)
    return out


def _price_leg(reg: Any, here_id: int, next_id: int, via: int | str, *,
               min_fidelity: Fidelity, solvent: str | None) -> dict[str, Any]:
    """One step, flattened: the number, the diagnosis, and what moved either way.

    `dE`, `added`, `solvent` and `shed` are in the route's direction (towards the
    target); `reaction_dE` is the edge's own, as recorded.
    """
    token = str(via).strip()
    direction = MADE_FROM
    if token.startswith(DERIVED):
        other = token[len(DERIVED):]
        priced = _derived_leg(reg, here_id, next_id,
                              int(other) if other else None,
                              min_fidelity=min_fidelity, solvent=solvent)
        origin, reaction_id = "inferred", None
        recorded_as = priced.get("recorded_as")
    else:
        if token.startswith(CONSUMED):
            direction, token = CONSUMED_BY, token[len(CONSUMED):]
        try:
            reaction_id = int(token)
        except ValueError:
            raise ReferenceSchemeError(
                f"{via!r} is not an edge: give a reaction id, c<reaction id> or "
                f"d<other part>") from None
        priced = price_reaction(reg, reaction_id, min_fidelity=min_fidelity,
                                solvent=solvent)
        origin, recorded_as = "recorded", None

    step = priced["steps"][0]
    terms = step["terms"]
    # The edge has to be the one that joins these two nodes, or the sum is of a different
    # route than the one being drawn.  Checked rather than assumed: a link can be stale.
    if origin == "recorded" and direction == MADE_FROM:
        _check_joins(reg, reaction_id, here_id, next_id, terms)
    elif origin == "recorded":
        _check_joins(reg, reaction_id, next_id, here_id, terms, role="reagent")

    # `added` and `solvent` are both consumed, so both count in the arithmetic — but they
    # are not the same claim.  Water carrying a proton away is not a piece being built in,
    # and a step that listed it as "added" would read as chemistry it is not doing.
    source = next_id if direction == MADE_FROM else here_id
    put_in = [_term(t) for t in terms
              if t["side"] == "reagent" and (t["role"] or "reagent") == "reagent"
              and int(t["structure_id"]) != source]
    carrier = [_term(t) for t in terms if t["side"] == "reagent" and t["role"] == "solvent"]
    let_go = [_term(t) for t in terms if t["role"] == "leaving"]
    sign = 1 if direction == MADE_FROM else -1
    dE = priced["total_dE"]
    if reaction_id is None:
        via_out: int | str = DERIVED
    else:
        via_out = reaction_id if sign > 0 else f"{CONSUMED}{reaction_id}"
    return {
        "via": via_out, "direction": direction, "sign": sign,
        "origin": origin, "reaction_id": reaction_id, "recorded_as": recorded_as,
        "product": here_id if sign > 0 else next_id, "source": source,
        "dE": None if dE is None else sign * dE, "reaction_dE": dE,
        "can_price": priced["can_price"],
        "why_not": priced["why_not"],
        "fidelity": step["fidelity"], "fidelity_name": step["fidelity_name"],
        "isodesmic": step["isodesmic"], "caveats": list(step["caveats"]),
        "all_converged": step["all_converged"],
        # Run against its arrow, an edge takes up what it shed and releases the rest.
        "added": put_in if sign > 0 else let_go,
        "solvent": carrier if sign > 0 else [],
        "shed": let_go if sign > 0 else put_in + carrier,
        "terms": terms,
    }


def _check_joins(reg: Any, reaction_id: int, product_id: int, source_id: int,
                 terms: list[dict[str, Any]], *, role: str | None = None) -> None:
    row = reg.conn.execute("SELECT product_structure_id FROM reactions WHERE id = ?",
                           (reaction_id,)).fetchone()
    if row is None or int(row["product_structure_id"]) != product_id:
        raise ReferenceSchemeError(
            f"reaction {reaction_id} does not produce structure {product_id}")
    if not any(t["side"] == "reagent" and int(t["structure_id"]) == source_id
               for t in terms):
        raise ReferenceSchemeError(
            f"reaction {reaction_id} does not consume structure {source_id}")
    # A consumed-by leg follows a `reagent` term only (registry.api.CONSUMING_ROLE).
    if role is not None and not any(
            t["side"] == "reagent" and int(t["structure_id"]) == source_id
            and (t["role"] or "reagent") == role for t in terms):
        raise ReferenceSchemeError(
            f"reaction {reaction_id} takes structure {source_id} only as a proton "
            f"carrier, not as a {role}; a consumed-by hop follows reagents only")


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
