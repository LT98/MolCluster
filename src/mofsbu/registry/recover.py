"""B28 recovery: link geometries stored before D30 to their structure's graph.

A geometry written before `geometries.canonical_order_json` existed may be in any build's
numbering.  Its coordinates are left exactly as they are; what is recovered is the link
(entry c = the row holding the atom at canonical position c):

* a geometry with no `relaxed_from` (a construct, an embed) is matched to its graph — the
  stored numbering if it realises the graph that way, else the numbering under which every
  graph bond is a contact (`BROKEN` times its ideal length, as E1 judges it).  More contacts
  than bonds are allowed: a construct may carry a marginal clash;
* a start pulled apart on purpose (choice vector with `perturb`) has bonds broken by design,
  so it is linked only if it is in the stored numbering with its covalent skeleton intact and
  no metal bond longer than the recorded pull-back explains;
* a relaxed geometry takes its parent's link: the optimiser keeps the rows it was handed.

Site states of a geometry that turns out not to be in the stored numbering were keyed
through the wrong map (B28) and are deleted — absent beats stale (ground rule 9); they are
written again the next time that geometry's state is.  A geometry that cannot be matched
keeps its NULL link, and readers refuse it.

Not called by `migrate`: it is the database owner's decision
(`scripts/recover_atom_links.py`, dry run by default).
"""
from __future__ import annotations

import json
from typing import Any

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism import GraphMatcher

from mofsbu.geometry.connectivity import BROKEN, _ideal, connectivity_check
from mofsbu.graph import EdgeType, TypedGraph
from mofsbu.registry.api import _refresh_open_sites, _stored_order, get_graph
from mofsbu.registry.db import Registry


def _parse(xyz_text: str) -> tuple[list[str], np.ndarray]:
    lines = xyz_text.strip().splitlines()
    rows = [ln.split() for ln in lines[2:2 + int(lines[0])]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def _contacts(symbols: list[str], xyz: np.ndarray, graph: TypedGraph) -> nx.Graph:
    """Element-labelled pairs close enough to be any bond the graph could hold."""
    metal = {graph.label(k).element: graph.is_metal(k) for k in graph.nodes()}
    out = nx.Graph()
    out.add_nodes_from((k, {"el": s}) for k, s in enumerate(symbols))
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            mi, mj = metal.get(symbols[i], False), metal.get(symbols[j], False)
            if mi != mj and "H" in (symbols[i], symbols[j]):
                continue
            et = (EdgeType.METAL_METAL if mi and mj else
                  EdgeType.DATIVE if mi != mj else EdgeType.COVALENT)
            if d[i, j] < BROKEN * _ideal(symbols[i], symbols[j], et, mi):
                out.add_edge(i, j)
    return out


def match_to_graph(symbols: list[str], xyz: np.ndarray, graph: TypedGraph,
                   ) -> list[int] | None:
    """Row k -> graph atom `[k]` under which every graph bond is a contact, or None."""
    n = len(graph)
    if len(symbols) != n:
        return None
    if symbols == [graph.label(k).element for k in graph.nodes()]:
        if not connectivity_check(symbols, xyz, graph).broken:
            return list(range(n))
    typed = nx.Graph()
    typed.add_nodes_from((k, {"el": graph.label(k).element}) for k in graph.nodes())
    typed.add_edges_from((u, v) for u, v, _et in graph.edges())
    gm = GraphMatcher(_contacts(symbols, xyz, graph), typed,
                      node_match=lambda a, b: a["el"] == b["el"])
    if not gm.subgraph_is_monomorphic():
        return None
    return [int(gm.mapping[k]) for k in range(n)]


def _pulled_back_in_stored_order(symbols: list[str], xyz: np.ndarray, graph: TypedGraph,
                                 pulled_by: float) -> bool:
    """A `perturb` start: covalent skeleton intact, metal bonds within the pull-back."""
    if symbols != [graph.label(k).element for k in graph.nodes()]:
        return False
    for u, v, et in graph.edges():
        d = float(np.linalg.norm(xyz[u] - xyz[v]))
        ideal = _ideal(symbols[u], symbols[v], et, graph.is_metal(u))
        slack = 0.0 if et is EdgeType.COVALENT else abs(pulled_by)
        if d > BROKEN * ideal + slack:
            return False
    return True


def recover_atom_links(reg: Registry, *, dry_run: bool = True) -> dict[str, Any]:
    """Recover (unless `dry_run`) the link of every geometry that has none.  Counts."""
    rows = reg.conn.execute(
        "SELECT id, structure_id, coords_hash, relaxed_from, choice_vector_json "
        "FROM geometries WHERE canonical_order_json IS NULL ORDER BY id").fetchall()
    linked = {int(r["id"]): [int(x) for x in json.loads(r["canonical_order_json"])]
              for r in reg.conn.execute(
                  "SELECT id, canonical_order_json FROM geometries "
                  "WHERE canonical_order_json IS NOT NULL")}
    graphs: dict[int, TypedGraph] = {}
    orders: dict[int, list[int]] = {}
    links: dict[int, list[int] | None] = {}
    out: dict[str, Any] = {"geometries": len(rows), "stored_numbering": 0,
                           "other_numbering": 0, "unmatched": [], "parent_unmatched": 0,
                           "states_dropped": 0, "structures": 0}
    touched: set[int] = set()
    for r in rows:
        gid, sid = int(r["id"]), int(r["structure_id"])
        graph = graphs.setdefault(sid, get_graph(reg, sid))
        stored = orders.setdefault(sid, _stored_order(reg, sid))
        parent = None if r["relaxed_from"] is None else int(r["relaxed_from"])
        if parent is not None and (parent in links or parent in linked):
            link = links[parent] if parent in links else linked[parent]
            if link is None:
                out["parent_unmatched"] += 1
                links[gid] = None
                continue
        else:
            symbols, xyz = _parse(reg.store.get(r["coords_hash"]).decode())
            pert = json.loads(r["choice_vector_json"] or "{}").get("perturb")
            if isinstance(pert, dict) and "by" in pert:
                row_to_atom = (list(range(len(graph)))
                               if _pulled_back_in_stored_order(symbols, xyz, graph,
                                                               float(pert["by"]))
                               else None)
            else:
                row_to_atom = match_to_graph(symbols, xyz, graph)
            if row_to_atom is None:
                out["unmatched"].append(gid)
                links[gid] = None
                continue
            row_of_atom = {atom: row for row, atom in enumerate(row_to_atom)}
            link = [row_of_atom[a] for a in stored]
        links[gid] = link
        in_stored = link == stored
        out["stored_numbering" if in_stored else "other_numbering"] += 1
        if dry_run:
            continue
        reg.conn.execute("UPDATE geometries SET canonical_order_json=? WHERE id=?",
                         (json.dumps(link), gid))
        if not in_stored:
            cur = reg.conn.execute("DELETE FROM site_state WHERE geometry_id=?", (gid,))
            out["states_dropped"] += cur.rowcount
            touched.add(sid)
    if not dry_run:
        for sid in touched:
            _refresh_open_sites(reg, sid)
        reg.conn.commit()
    out["structures"] = len(touched)
    return out
