"""A relaxation that changed connectivity is filed under the species it became (D31).

E1 flags a relaxed geometry whose bonding no longer matches its structure's graph and keeps
it out of energy selection.  When the change leaves one piece — a donor closing onto the
metal, a proton crossing an H-bond — the result is filed as its own structure (a different
graph, so a new L1) with the relaxed coordinates and energy as built, and a `relaxation`
edge from the species it was built as.  The flagged geometry stays where it is, as history.

What `derive_species` derives, and what it refuses with its reason:

* a DATIVE bond formed or broken between a donor and a metal — applied to the graph;
* a proton transfer — an X-H broken and a Y-H formed on the same H, X and Y non-metal,
  non-carbon atoms.  A localised atom's formal charge moves by one; a delocalised one keeps
  label 0 because its charge lives on the graph (D15), the same rule `from_rdkit` uses;
* anything else (a covalent skeleton broken, a C-H, a metal-metal bond) is refused — it
  reads as a failed relaxation, not a species.

When the result falls apart, every metal must stay in one piece, and every other piece must
be a neutral closed-shell metal-free leaver (no formal charge, nothing that could carry
delocalised charge, an even electron count: water, ammonia, an alcohol).  Anything else
would need a per-piece charge or spin nobody can state, so it is refused.

`file_derived` files no separation: whether a piece is held in the second shell or keeps
leaving is pending the second-shell probe (D31), so a separation comes back as a refusal.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

import networkx as nx

from mofsbu._types import Fidelity, MofsbuError
from mofsbu.graph import EdgeType, TypedGraph


class DeriveRefused(MofsbuError):
    """The connectivity change is not one this module files, and why."""


@dataclass(frozen=True)
class Piece:
    graph: TypedGraph
    rows: tuple[int, ...]           # piece atom k is row rows[k] of the relaxed geometry
    leaving: bool


@dataclass(frozen=True)
class Derived:
    pieces: tuple[Piece, ...]       # pieces[0] holds every metal (or is the whole graph)
    summary: str


def _orders(d: dict) -> dict[tuple[int, int], float]:
    return {(min(b["i"], b["j"]), max(b["i"], b["j"])): float(b.get("order", 1.0))
            for b in d["bonds"] if b["t"] == EdgeType.COVALENT.value}


def _delocalised(g: TypedGraph, orders: dict[tuple[int, int], float], a: int) -> bool:
    """`from_mol._is_delocalised`, read off the graph's covalent bond orders."""
    def order(u: int, v: int) -> float:
        return orders.get((min(u, v), max(u, v)), 0.0)

    for b in g.neighbors(a, EdgeType.COVALENT):
        if order(a, b) > 1.0:
            return True
        if any(order(b, c) > 1.0 for c in g.neighbors(b, EdgeType.COVALENT)):
            return True
    return False


def _hetero(g: TypedGraph, a: int) -> bool:
    return not g.is_metal(a) and g.label(a).element not in ("C", "H")


def derive_species(g: TypedGraph, verdict: dict[str, Any]) -> Derived:
    """The species a relaxation became: `g` numbered like the geometry's rows, E1's verdict."""
    if not verdict.get("changed"):
        raise DeriveRefused("the verdict records no connectivity change")
    remove: list[tuple[int, int]] = []
    add: list[tuple[int, int, EdgeType]] = []
    lost: dict[int, int] = {}                 # H -> heavy atom it left
    gained: dict[int, int] = {}               # H -> heavy atom it joined
    for side, store in (("broken", lost), ("formed", gained)):
        for c in verdict.get(side, []):
            i, j, et = int(c["i"]), int(c["j"]), c["etype"]
            if et == EdgeType.DATIVE.value:
                (remove if side == "broken" else add).append(
                    (i, j) if side == "broken" else (i, j, EdgeType.DATIVE))
                continue
            if et == EdgeType.COVALENT.value and (c["sym_i"] == "H") != (c["sym_j"] == "H"):
                h, x = (i, j) if c["sym_i"] == "H" else (j, i)
                if _hetero(g, x):
                    store[h] = x
                    continue
            raise DeriveRefused(
                f"{c['sym_i']}{i}-{c['sym_j']}{j} ({et}) {side}: the skeleton changed, which "
                f"reads as a failed relaxation rather than a species")
    if set(lost) != set(gained):
        raise DeriveRefused(f"hydrogens {sorted(set(lost) ^ set(gained))} changed partner on "
                            f"one side only; that is not a proton transfer")

    d = g.to_dict()
    orders = _orders(d)
    charge_shift = {a: 0 for a in g.nodes()}
    for h, x in lost.items():
        y = gained[h]
        remove.append((h, x))
        add.append((h, y, EdgeType.COVALENT))
        if not _delocalised(g, orders, x):
            charge_shift[x] -= 1
        if not _delocalised(g, orders, y):
            charge_shift[y] += 1
    gone = {(min(u, v), max(u, v)) for u, v in remove}
    d["bonds"] = [b for b in d["bonds"]
                  if (min(b["i"], b["j"]), max(b["i"], b["j"])) not in gone]
    d["bonds"] += [{"i": u, "j": v, "t": et.value, "order": 1.0} for u, v, et in add]
    for atom in d["atoms"]:
        atom["q"] = int(atom.get("q", 0)) + charge_shift[atom["i"]]
    new = TypedGraph.from_dict(d)

    parts = [sorted(c) for c in nx.connected_components(new.to_networkx())]
    summary = (f"{len([r for r in remove if g.edge_type(*r) is EdgeType.DATIVE])} dative "
               f"broken, {len([a for a in add if a[2] is EdgeType.DATIVE])} dative formed, "
               f"{len(lost)} proton(s) moved")
    if len(parts) == 1:
        return Derived((Piece(new, tuple(new.nodes()), False),), summary)

    holding = [p for p in parts if any(new.is_metal(a) for a in p)]
    if len(holding) > 1:
        raise DeriveRefused("the metals ended up in different pieces; a per-metal charge and "
                            "spin cannot be derived")
    if not holding:
        raise DeriveRefused("a metal-free species fell apart; its pieces' charges cannot be "
                            "derived")
    main = holding[0]
    new_orders = _orders(new.to_dict())
    pieces = [Piece(new.subgraph(main).with_charge(new.net_charge())
                    .with_multiplicity(new.multiplicity), tuple(main), False)]
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    for part in parts:
        if part is main:
            continue
        labels = [new.label(a) for a in part]
        if any(lab.formal_charge for lab in labels):
            raise DeriveRefused(f"leaving piece {_formula(labels)} carries a formal charge; "
                                f"only a neutral closed-shell leaver is filed")
        if any(_delocalised(new, new_orders, a) for a in part):
            raise DeriveRefused(f"leaving piece {_formula(labels)} could carry delocalised "
                                f"charge, which the graph does not apportion")
        if sum(table.GetAtomicNumber(lab.element) for lab in labels) % 2:
            raise DeriveRefused(f"leaving piece {_formula(labels)} has an odd electron "
                                f"count; its spin is not derivable")
        pieces.append(Piece(new.subgraph(part).with_charge(0).with_multiplicity(1),
                            tuple(part), True))
    for p in pieces:
        p.graph.name = g.name if not p.leaving else _formula([p.graph.label(a)
                                                              for a in p.graph.nodes()])
    return Derived(tuple(pieces), summary + f", split into {len(pieces)} pieces")


def _formula(labels) -> str:
    counts: dict[str, int] = {}
    for lab in labels:
        counts[lab.element] = counts.get(lab.element, 0) + 1
    return "".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(counts.items()))


def file_derived(reg: Any, geometry_id: int) -> dict[str, Any]:
    """File the species a connectivity-changed geometry became.  Writes; does not commit.

    Idempotent: a geometry already filed (another geometry `derived_from` it) is skipped.
    A refusal, a separation included, is returned with its reason, never raised.
    """
    from mofsbu.energy.reference import _method_of, check_balance
    from mofsbu.registry.api import (
        Provenance,
        RegistryError,
        geometry_graph,
        put_geometry,
        put_reaction,
        put_structure,
    )

    row = reg.conn.execute(
        "SELECT structure_id, fidelity, method_id, energy, converged, coords_hash, "
        "connectivity_changed, connectivity_json FROM geometries WHERE id=?",
        (geometry_id,)).fetchone()
    if row is None:
        raise RegistryError(f"no geometry {geometry_id}")
    if row["connectivity_changed"] != 1:
        return {"skipped": "connectivity not changed"}
    done = reg.conn.execute("SELECT id FROM geometries WHERE derived_from=?",
                            (geometry_id,)).fetchall()
    if done:
        return {"skipped": "already filed", "geometries": [int(r["id"]) for r in done]}
    try:
        derived = derive_species(geometry_graph(reg, geometry_id),
                                 json.loads(row["connectivity_json"]))
    except DeriveRefused as exc:
        return {"refused": str(exc)}

    if len(derived.pieces) > 1:
        return {"refused": "separation — second shell or leaving, pending",
                "summary": derived.summary}

    # All three rows or none: a failure after the first write must not leave a structure
    # that the "already filed" check would then accept.
    began = not reg.conn.in_transaction
    if began:
        reg.conn.execute("BEGIN")
    reg.conn.execute("SAVEPOINT file_derived")
    try:
        put = put_structure(reg, derived.pieces[0].graph, l2="", tags=["derived"])
        geom = put_geometry(
            reg, put.id, reg.store.get(row["coords_hash"]).decode(),
            fidelity=Fidelity(int(row["fidelity"])),
            method=_method_of(reg, row["method_id"]), energy=row["energy"],
            converged=None if row["converged"] is None else bool(row["converged"]),
            canonical_order=put.canonical_order, derived_from=geometry_id)
        rid = put_reaction(reg, put.id, Provenance(
            kind="relaxation", reagent_ids=(int(row["structure_id"]),),
            note=f"geometry {geometry_id}: {derived.summary}"))
        balance = check_balance(reg, rid)
        if not balance.balanced:
            raise MofsbuError(f"relaxation edge {rid} from geometry {geometry_id} is "
                              f"{balance.describe()}; a derived species must conserve atoms")
    except BaseException:
        reg.conn.execute("ROLLBACK TO SAVEPOINT file_derived")
        reg.conn.execute("RELEASE SAVEPOINT file_derived")
        if began:
            reg.conn.execute("ROLLBACK")
        raise
    reg.conn.execute("RELEASE SAVEPOINT file_derived")
    return {"filed": [put.id], "geometries": [geom.id], "created": [put.created],
            "edge": rid, "summary": derived.summary}
