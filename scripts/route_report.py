"""Report what relaxation did to every species a spec's routes stored (`spec.routes`).

    python scripts/route_report.py --db data/test/m6_s5_ground_up/registry.db
    python scripts/route_report.py --db … --store … --json report.json

Reads the registry read-only.  For every `route` task of every run it finds the species
the route stored (built products and `perturb` starts), every relaxation of each one, and
reports, per relaxed geometry:

  * the `watch` distances the route declared, at the start and after relaxation;
  * which built species of the same composition the relaxed connectivity matches (graph
    isomorphism of element-labelled contact graphs — Cu–O inside 2.5 Å, covalent pairs
    inside 1.15 × the radii sum), so "route B reached route A's product" is a query;
  * ΔE against each relaxed built species of the same composition under the same model;
  * the coordination number of every labelled metal, start → final;
  * a water audit: every co-ligand atom set (`<metal>.cap` from a `cap` step, `<load>.H2O`
    from a `load`) is checked for a co-ligand that lost or gained an H, left or changed its
    metal, or gained an H-bond (H···O < 2.1 Å) to an oxygen outside the set that it did not
    have at the start.  Empty is benign.

A failed relax task is reported with its error, not dropped.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.geometry.distances import COVALENT_RADII, DEFAULT_COVALENT   # noqa: E402
from mofsbu.registry import BlobStore, geometry_xyz                       # noqa: E402
from mofsbu.registry.db import ReadOnlyRegistry                           # noqa: E402
from mofsbu.runner import read_xyz                                        # noqa: E402

#: Metal–O counted as coordinated inside this (Å); H···O an H-bond inside HBOND.
M_O, OH_COV, HBOND = 2.5, 1.25, 2.1
METALS = {"Cu", "Ni", "Zn", "Fe", "Co", "Mn", "Cr", "Mo", "Rh", "Zr"}


def contact_graph(sym, xyz) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from((i, {"el": s}) for i, s in enumerate(sym))
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            mi, mj = sym[i] in METALS, sym[j] in METALS
            if mi or mj:
                ok = (mi != mj) and "O" in (sym[i], sym[j]) and d[i, j] < M_O
            else:
                ok = d[i, j] < 1.15 * (COVALENT_RADII.get(sym[i], DEFAULT_COVALENT)
                                       + COVALENT_RADII.get(sym[j], DEFAULT_COVALENT))
            if ok:
                g.add_edge(i, j)
    return g


def same_node(a: nx.Graph, b: nx.Graph) -> bool:
    return nx.is_isomorphic(a, b, node_match=lambda x, y: x["el"] == y["el"])


def _cn(sym, xyz, m) -> int:
    return int(sum(1 for o, s in enumerate(sym) if s == "O" and np.linalg.norm(xyz[m] - xyz[o]) < M_O))


def water_audit(sym, x0, x1, rows) -> list[str]:
    cap = sorted({r for lab, rr in rows.items()
                  if lab.endswith(".cap") or lab.endswith(".H2O") for r in rr})
    cap_o = [i for i in cap if sym[i] == "O"]
    other_o = [i for i, s in enumerate(sym) if s == "O" and i not in cap]
    metals = [i for i, s in enumerate(sym) if s in METALS]
    hs = [i for i, s in enumerate(sym) if s == "H"]

    def state(x):
        d = np.linalg.norm(x[:, None] - x[None], axis=-1)
        out = {}
        for o in cap_o:
            mine = [h for h in hs if d[o, h] < OH_COV]
            out[o] = ([m for m in metals if d[o, m] < M_O], len(mine),
                      {f for h in mine for f in other_o if d[h, f] < HBOND})
        return out

    a, b = state(x0), state(x1)
    issues = []
    for o in cap_o:
        (m0, h0, hb0), (m1, h1, hb1) = a[o], b[o]
        if h1 != h0:
            issues.append(f"cap O{o}: {h0}→{h1} H")
        if m1 != m0:
            issues.append(f"cap O{o}: metal {m0}→{m1}")
        if hb1 - hb0:
            issues.append(f"cap O{o}: new H-bond to O{sorted(hb1 - hb0)}")
    return issues


def _geom(reg, gid):
    sym, xyz = read_xyz(geometry_xyz(reg, gid))
    return sym, np.asarray(xyz, dtype=float)


def _model(con, method_id) -> str:
    r = con.execute("SELECT code, method, extras_json FROM methods WHERE id=?",
                    (method_id,)).fetchone()
    return f"{r['code']}:{r['method']}" if r else "?"


def report(con, reg) -> list[dict]:
    out = []
    for run in con.execute("SELECT id, spec_json FROM runs ORDER BY id").fetchall():
        spec = json.loads(run["spec_json"])
        for t in con.execute("SELECT id, status, error, detail_json, "
                             "json_extract(payload_json, '$.name') AS route FROM tasks WHERE "
                             "run_id=? AND kind='route' ORDER BY id", (run["id"],)):
            if t["status"] != "done":
                out.append({"run": run["id"], "ml_model": spec.get("ml_model"),
                            "species": f"route {t['route']}", "route_task": t["id"],
                            "error": _first_line(t["error"])})
                continue
            species = json.loads(t["detail_json"])["species"]
            built = {sp["name"]: sp for sp in species if sp["kind"] == "built"}
            out += _species_rows(con, reg, run["id"], spec.get("ml_model"), species, built)
    # Cross-route: compare each relaxed geometry with every built species of its formula.
    built_all = {}
    for row in out:
        if row.get("kind") == "built":
            built_all.setdefault((row["run"], row["formula"]), {})[row["species"]] = row
    for row in out:
        if "relaxed" not in row:
            continue
        peers = built_all.get((row["run"], row["formula"]), {})
        g = contact_graph(row["_sym"], row["_x1"])
        row["matches_built"] = sorted(n for n, p in peers.items()
                                      if same_node(g, contact_graph(p["_sym"], p["_x0"])))
        row["dE_vs_built"] = {n: round(row["energy"] - p["energy"], 3) for n, p in peers.items()
                              if p.get("energy") is not None and row.get("energy") is not None}
    for row in out:
        for k in ("_sym", "_x0", "_x1"):
            row.pop(k, None)
    return out


def _first_line(text: str | None) -> str:
    lines = (text or "").strip().splitlines()
    return lines[0][:300] if lines else ""


def built_report(con, reg) -> list[dict]:
    """Construct-only: per route variant, what it built — QC, bridges, and which node.

    Every stored species that declares `watch` distances is a route's product; each is
    compared, by contact-graph isomorphism, with every other such product of the same
    formula in its run, so "these variants reach one node" is read off the table.  A route
    that was refused is listed with its refusal.
    """
    from mofsbu.registry import get_graph

    rows, by_formula = [], {}
    for t in con.execute("SELECT run_id, status, error, detail_json, "
                         "json_extract(payload_json, '$.name') AS route FROM tasks "
                         "WHERE kind='route' ORDER BY id"):
        if t["status"] != "done":
            rows.append({"route": t["route"], "refused": _first_line(t["error"])})
            continue
        for sp in json.loads(t["detail_json"])["species"]:
            if sp["kind"] != "built" or not sp.get("watch"):
                continue
            sym, x = _geom(reg, sp["geometry_id"])
            r = sp["rows"]
            g = get_graph(reg, sp["structure_id"])
            qc = json.loads(con.execute("SELECT qc_json FROM geometries WHERE id=?",
                                        (sp["geometry_id"],)).fetchone()[0] or "{}")
            worst = (qc.get("worst_overlap") or {}).get("overlap") \
                if isinstance(qc.get("worst_overlap"), dict) else qc.get("worst_overlap")
            row = {"route": t["route"], "species": sp["name"],
                   "structure_id": sp["structure_id"],
                   "qc": "missing" if not qc else "ok" if qc.get("ok") else
                         "marginal" if qc.get("marginal") else "clash",
                   "worst_overlap": worst,
                   "mu2": sum(g.fragment_bridge_class(f).value == "mu2"
                              for f in g.ligand_fragments()),
                   "watch": {k: round(float(np.linalg.norm(x[r[a]].mean(axis=0)
                                                           - x[r[b]].mean(axis=0))), 3)
                             for k, (a, b) in sp["watch"].items()},
                   "_graph": contact_graph(sym, x)}
            formula = "".join(f"{e}{sym.count(e)}" for e in sorted(set(sym)))
            by_formula.setdefault((t["run_id"], formula), []).append(row)
            rows.append(row)
    for peers in by_formula.values():
        for row in peers:
            row["same_node_as"] = sorted(f"{p['route']}/{p['species']}" for p in peers
                                         if p is not row and same_node(row["_graph"], p["_graph"]))
    for row in rows:
        row.pop("_graph", None)
    return rows


def _species_rows(con, reg, run_id, ml_model, species, built):
    rows = []
    for sp in species:
        sym, x0 = _geom(reg, sp["geometry_id"])
        formula = "".join(f"{e}{sym.count(e)}" for e in sorted(set(sym)))
        base = {"run": run_id, "ml_model": ml_model, "species": sp["name"], "kind": sp["kind"],
                "structure_id": sp["structure_id"], "start_geometry": sp["geometry_id"],
                "formula": formula, "_sym": sym, "_x0": x0}
        relax = con.execute(
            "SELECT t.status, t.error, t.geometry_id, t.detail_json FROM tasks t WHERE "
            "t.run_id=? AND t.kind='relax' AND json_extract(t.payload_json,'$.geometry_id')=?"
            " ORDER BY t.id DESC LIMIT 1", (run_id, sp["geometry_id"])).fetchone()
        gid = None
        if relax is not None and relax["status"] == "done":
            gid = relax["geometry_id"]
        elif relax is None:
            reused = con.execute(
                "SELECT g.id FROM geometries g JOIN methods m ON m.id=g.method_id WHERE "
                "g.relaxed_from=? ORDER BY g.id", (sp["geometry_id"],)).fetchall()
            gid = reused[-1]["id"] if reused else None
        if gid is None:
            if relax is None:
                why = "no relax task"
            else:
                first = (relax["error"] or "").strip().splitlines()
                why = f"relax {relax['status']}" + (f": {first[0][:200]}" if first else "")
            rows.append(base | {"error": why})
            continue
        g = con.execute("SELECT energy, converged, method_id FROM geometries WHERE id=?",
                        (gid,)).fetchone()
        _, x1 = _geom(reg, gid)
        r = sp["rows"]
        dist = lambda x, a, b: round(float(np.linalg.norm(  # noqa: E731
            x[r[a]].mean(axis=0) - x[r[b]].mean(axis=0))), 3)
        metals = {lab: rr[0] for lab, rr in r.items()
                  if "." not in lab and len(rr) == 1 and sym[rr[0]] in METALS}
        rows.append(base | {
            "relaxed": gid, "method": _model(con, g["method_id"]), "energy": g["energy"],
            "converged": bool(g["converged"]), "_x1": x1,
            "watch": {k: [dist(x0, a, b), dist(x1, a, b)] for k, (a, b) in sp["watch"].items()},
            "cn": {lab: [_cn(sym, x0, m), _cn(sym, x1, m)] for lab, m in metals.items()},
            "water_issues": water_audit(sym, x0, x1, r),
            "perturbation": sp.get("perturbation")})
    return rows


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, required=True)
    p.add_argument("--store", type=Path, default=None, help="default: <db dir>/store")
    p.add_argument("--json", type=Path, default=None, help="also write the rows as JSON")
    p.add_argument("--built", action="store_true",
                   help="construct-only: QC, bridges and node identity of every route product")
    a = p.parse_args(argv)

    con = sqlite3.connect(f"file:{a.db.resolve()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    reg = ReadOnlyRegistry(conn=con, store=BlobStore(a.store or a.db.parent / "store"))
    if a.built:
        rows = built_report(con, reg)
        for row in rows:
            if "refused" in row:
                print(f"{row['route']:48s} REFUSED: {row['refused']}")
            else:
                print(f"{row['route']:48s} {row['species']:12s} qc {row['qc']:8s} "
                      f"mu2 {row['mu2']} watch {row['watch']} "
                      f"same node as {len(row['same_node_as'])} other product(s)")
        if a.json:
            a.json.write_text(json.dumps(rows, indent=1))
        return 0
    rows = report(con, reg)
    for row in rows:
        head = f"{(row.get('ml_model') or '?'):15s} {row['species']:22s}"
        if "error" in row:
            print(f"{head} FAILED: {row['error']}")
            continue
        print(f"{head} conv {row['converged']!s:5s} matches {row['matches_built']} "
              f"watch {row['watch']} cn {row['cn']} dE {row['dE_vs_built']} "
              f"water {row['water_issues'] or 'benign'}")
    if a.json:
        a.json.write_text(json.dumps(rows, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
