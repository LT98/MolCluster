"""Relax one sample with two ML backends from the same start and compare the minima (E-MH).

    python scripts/compare_relaxers.py --db data/salt_study_k1.db --out mh_vs_omol.json
    python scripts/compare_relaxers.py --db … --per-class 15 --conformers 3 --fmax 0.20

The simple form of `docs/WORKPLAN_energy.md` §3a study E-MH: reads the registry read-only,
writes nothing to it, and reports per relaxation —

  (a) bonds of the stored graph broken / non-bonds formed, judged from distances;
  (b) heavy-atom RMSD and the largest metal–donor distance change between the two minima;
  (c) the reference model's single point on the candidate minimum minus its own minimum;
  (d) within-identity conformer ranking agreement (top-1);
  (f) convergence, steps and wall time per relaxation.

A second pass starts the candidate from the reference minimum, which separates "different
basin" from "different start".  Metric (e), re-pricing the route equations, is not done here.
"""
from __future__ import annotations

import argparse
import json
import random
import sqlite3
import statistics
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.energy.backends import get_backend                           # noqa: E402
from mofsbu.geometry._linalg import kabsch                               # noqa: E402
from mofsbu.geometry.distances import COVALENT_RADII, DEFAULT_COVALENT, metal_donor_distance  # noqa: E402
from mofsbu.graph._types import EdgeType                                 # noqa: E402
from mofsbu.registry import BlobStore, get_graph, geometry_xyz           # noqa: E402
from mofsbu.registry.db import ReadOnlyRegistry                          # noqa: E402
from mofsbu.runner import read_xyz                                       # noqa: E402

#: Demo-route nodes (WORKPLAN_energy §3a step 2): always sampled when present.
ROUTE_NODES = (132, 94, 28, 68, 328, 301)

#: A bond is broken past this multiple of its ideal length, formed inside the second.
BROKEN, FORMED = 1.30, 1.10


def _rc(sym: str) -> float:
    return COVALENT_RADII.get(sym, DEFAULT_COVALENT)


def bond_changes(graph, symbols, xyz) -> dict:
    """Graph edges whose length left the bonded range, and non-edges that entered it."""
    n = len(symbols)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    edges = {(u, v): et for u, v, et in graph.edges()}
    metals = set(graph.metals())

    def ideal(i, j, dative):
        if dative:
            m, x = (i, j) if i in metals else (j, i)
            return metal_donor_distance(symbols[m], symbols[x]).value
        return _rc(symbols[i]) + _rc(symbols[j])

    broken, formed = [], []
    for (u, v), et in edges.items():
        if d[u, v] > BROKEN * ideal(u, v, et is EdgeType.DATIVE):
            broken.append([u, v, round(float(d[u, v]), 3)])
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in edges or (i in metals and j in metals):
                continue
            dative = (i in metals) != (j in metals)
            if dative and "H" in (symbols[i], symbols[j]):
                continue
            if d[i, j] < FORMED * ideal(i, j, dative):
                formed.append([i, j, round(float(d[i, j]), 3)])
    return {"broken": broken, "formed": formed}


def heavy_rmsd(symbols, a, b) -> float:
    idx = [i for i, s in enumerate(symbols) if s != "H"]
    if len(idx) < 3:
        idx = list(range(len(symbols)))
    p, q = a[idx], b[idx]
    if len(idx) < 3:
        return float(np.sqrt((((p - p.mean(0)) - (q - q.mean(0))) ** 2).sum(1).mean()))
    rot, t = kabsch(p, q)
    return float(np.sqrt((((p @ rot.T) + t - q) ** 2).sum(1).mean()))


def metal_donor_shift(graph, a, b) -> float | None:
    shifts = [abs(np.linalg.norm(a[u] - a[v]) - np.linalg.norm(b[u] - b[v]))
              for u, v, et in graph.edges() if et is EdgeType.DATIVE]
    return float(max(shifts)) if shifts else None


def sample(con: sqlite3.Connection, reference_method: str, per_class: int, conformers: int,
           seed: int) -> list[dict]:
    """Start geometries: each is the `relaxed_from` of a stored reference-model relaxation."""
    rows = con.execute(
        """SELECT s.id sid, s.n_metals, s.net_charge, s.n_atoms, s.formula,
                  r.id start_gid
           FROM geometries g JOIN methods m ON m.id = g.method_id
           JOIN geometries r ON r.id = g.relaxed_from
           JOIN structures s ON s.id = g.structure_id
           WHERE m.method = ? AND s.hidden = 0 AND r.fidelity <= 1
           ORDER BY s.id, r.id""", (reference_method,)).fetchall()
    by_sid: dict[int, dict] = {}
    for r in rows:
        e = by_sid.setdefault(r["sid"], {"sid": r["sid"], "n_metals": r["n_metals"],
                                         "charge": r["net_charge"], "n_atoms": r["n_atoms"],
                                         "formula": r["formula"], "starts": []})
        if r["start_gid"] not in e["starts"]:
            e["starts"].append(r["start_gid"])
    rng = random.Random(seed)
    classes: dict[str, list[dict]] = {"free": [], "metal_neutral": [], "metal_charged": []}
    for e in by_sid.values():
        cls = "free" if e["n_metals"] == 0 else (
            "metal_neutral" if e["charge"] == 0 else "metal_charged")
        e["class"] = cls
        classes[cls].append(e)
    picked: list[dict] = []
    for cls, pool in classes.items():
        forced = [e for e in pool if e["sid"] in ROUTE_NODES]
        rest = [e for e in pool if e["sid"] not in ROUTE_NODES]
        rng.shuffle(rest)
        take = pool if cls == "free" else (forced + rest)[:max(per_class, len(forced))]
        picked.extend(take)
    for e in picked:
        e["starts"] = e["starts"][:conformers]
    return picked


def relax(backend, symbols, xyz, graph, fmax, steps) -> dict:
    t0 = time.perf_counter()
    try:
        r = backend.relax(symbols, xyz, charge=int(graph.charge),
                          multiplicity=int(graph.multiplicity), fmax=fmax, steps=steps)
    except Exception as exc:                                             # noqa: BLE001
        return {"error": f"{type(exc).__name__}: {exc}"}
    return {"energy": r.energy, "converged": r.converged, "steps": r.n_steps,
            "fmax": r.fmax, "seconds": time.perf_counter() - t0,
            "xyz": np.asarray(r.positions, dtype=float)}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, required=True)
    p.add_argument("--store", type=Path, default=None,
                   help="blob store (default: <db dir>/store)")
    p.add_argument("--out", type=Path, required=True,
                   help="with several candidates, one file each: <out stem>.<key>.json")
    p.add_argument("--reference", default="mace_omol", help="backend key of today's relaxer")
    p.add_argument("--candidate", nargs="+", default=["mace_mh"],
                   help="one or more backend keys; the reference relaxes each start once")
    p.add_argument("--head", default=None, help="MACE-MH head (default: the backend's)")
    p.add_argument("--per-class", type=int, default=15)
    p.add_argument("--conformers", type=int, default=3)
    p.add_argument("--fmax", type=float, default=0.20, help="production default: 0.20 eV/Å")
    p.add_argument("--steps", type=int, default=500)
    p.add_argument("--device", default=None)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--resume", action="store_true",
                   help="keep the records already in the output files and skip their starts")
    a = p.parse_args(argv)

    con = sqlite3.connect(f"file:{a.db.resolve()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    reg = ReadOnlyRegistry(conn=con, store=BlobStore(a.store or a.db.parent / "store"))

    kw = {} if a.device is None else {"device": a.device}
    ref = get_backend(a.reference, **kw)
    cands = {key: get_backend(key, **kw, **({"head": a.head} if a.head else {}))
             for key in a.candidate}
    ref_method = ref.method
    picked = sample(con, ref_method, a.per_class, a.conformers, a.seed)
    print(f"{len(picked)} structures, {sum(len(e['starts']) for e in picked)} starts; "
          f"reference {ref.method}, candidates "
          f"{[c.method + ' ' + getattr(c, 'head', '') for c in cands.values()]}", flush=True)

    records: dict[str, list[dict]] = {key: [] for key in cands}
    done: set[tuple[int, int]] = set()
    if a.resume:
        for key in cands:
            out = _out(a, key, len(cands))
            if out.is_file():
                records[key] = json.loads(out.read_text())["records"]
        # a start counts as done only when every candidate has it
        done = set.intersection(*({(r["sid"], r["start_geometry"]) for r in recs}
                                  for recs in records.values()))
        for key in cands:
            records[key] = [r for r in records[key]
                            if (r["sid"], r["start_geometry"]) in done]
        print(f"resuming: {len(done)} starts already done", flush=True)
    for e in picked:
        graph = get_graph(reg, e["sid"])
        for gid in e["starts"]:
            if (e["sid"], gid) in done:
                continue
            symbols, coords = read_xyz(geometry_xyz(reg, gid))
            # A structure reached by two routes keeps one graph (D22); a geometry from the
            # other route may order its atoms differently, so per-atom metrics are skipped.
            ordered = [graph.label(i).element for i in graph.nodes()] == list(symbols)
            x0 = np.asarray(coords, dtype=float)
            rr = relax(ref, symbols, x0, graph, a.fmax, a.steps)
            for key, cand in cands.items():
                rec = compare_one(e, gid, graph, symbols, x0, rr, ref, cand, ordered, a)
                records[key].append(rec)
                print(json.dumps({"candidate": key} | {k: rec.get(k) for k in (
                    "sid", "class", "start_geometry", "heavy_rmsd", "cross_rmsd",
                    "dE_ref_at_candidate_min")}), flush=True)
            write(a, ref, cands, records)       # after every start: a killed run keeps its data

    for key, summary in write(a, ref, cands, records).items():
        print(key, json.dumps(summary, indent=1))
    return 0


def _out(a, key: str, n: int) -> Path:
    return a.out if n == 1 else a.out.with_name(f"{a.out.stem}.{key}.json")


def write(a, ref, cands, records) -> dict[str, dict]:
    summaries = {}
    for key, cand in cands.items():
        out = _out(a, key, len(cands))
        summaries[key] = summarise(records[key])
        out.write_text(json.dumps({
            "reference": ref.method_spec(charge=0, multiplicity=1).describe(),
            "candidate": cand.method_spec(charge=0, multiplicity=1).describe(),
            "candidate_head": getattr(cand, "head", None),
            "reference_key": a.reference, "candidate_key": key,
            "fmax": a.fmax, "steps": a.steps, "seed": a.seed,
            "summary": summaries[key], "records": records[key]}, indent=1))
    return summaries


def compare_one(e, gid, graph, symbols, x0, rr, ref, cand, ordered, a) -> dict:
    """One start: the candidate's relaxation against the reference's, metrics (a)–(c), (f)."""
    rc = relax(cand, symbols, x0, graph, a.fmax, a.steps)
    rec = {k: e[k] for k in ("sid", "class", "charge", "n_atoms", "formula")}
    rec["start_geometry"] = gid
    rec["reference"] = {k: v for k, v in rr.items() if k != "xyz"}
    rec["candidate"] = {k: v for k, v in rc.items() if k != "xyz"}
    if "xyz" in rr and "xyz" in rc:
        rec["heavy_rmsd"] = heavy_rmsd(symbols, rc["xyz"], rr["xyz"])
        if ordered:
            rec["reference"]["bonds"] = bond_changes(graph, symbols, rr["xyz"])
            rec["candidate"]["bonds"] = bond_changes(graph, symbols, rc["xyz"])
            rec["metal_donor_shift"] = metal_donor_shift(graph, rc["xyz"], rr["xyz"])
        sp = ref.single_point(symbols, rc["xyz"], charge=int(graph.charge),
                              multiplicity=int(graph.multiplicity))
        rec["ref_sp_on_candidate"] = sp.energy
        rec["dE_ref_at_candidate_min"] = sp.energy - rr["energy"]
        cross = relax(cand, symbols, rr["xyz"], graph, a.fmax, a.steps)
        if "xyz" in cross:
            rec["cross_rmsd"] = heavy_rmsd(symbols, cross["xyz"], rr["xyz"])
    return rec


def _med(xs):
    xs = [x for x in xs if x is not None]
    return round(statistics.median(xs), 4) if xs else None


def summarise(records: list[dict]) -> dict:
    out: dict = {}
    for cls in ("all", "free", "metal_neutral", "metal_charged"):
        rs = [r for r in records if cls == "all" or r["class"] == cls]
        ok = [r for r in rs if "heavy_rmsd" in r]
        if not rs:
            continue

        def changed(r, side):
            b = r[side].get("bonds")
            return bool(b and (b["broken"] or b["formed"]))

        def n_changed(side):
            return sum(1 for r in ok if changed(r, side))

        out[cls] = {
            "n": len(rs), "failed_ref": sum("error" in r["reference"] for r in rs),
            "failed_cand": sum("error" in r["candidate"] for r in rs),
            "unconverged_ref": sum(not r["reference"].get("converged", True) for r in rs),
            "unconverged_cand": sum(not r["candidate"].get("converged", True) for r in rs),
            "bond_change_ref": n_changed("reference"),
            "bond_change_cand": n_changed("candidate"),
            "bond_change_cand_only": sum(
                1 for r in ok if changed(r, "candidate") and not changed(r, "reference")),
            "bonds_unchecked": sum(1 for r in ok if "bonds" not in r["reference"]),
            "median_heavy_rmsd": _med([r["heavy_rmsd"] for r in ok]),
            "max_heavy_rmsd": round(max((r["heavy_rmsd"] for r in ok), default=0), 4),
            "median_cross_rmsd": _med([r.get("cross_rmsd") for r in ok]),
            "median_metal_donor_shift": _med([r.get("metal_donor_shift") for r in ok]),
            "n_dE_over_0.05": sum(1 for r in ok if r["dE_ref_at_candidate_min"] > 0.05),
            "median_dE_ref_at_cand_min": _med([r["dE_ref_at_candidate_min"] for r in ok]),
            "max_dE_ref_at_cand_min": round(max((r["dE_ref_at_candidate_min"] for r in ok),
                                                default=0), 4),
            "median_seconds_ref": _med([r["reference"].get("seconds") for r in rs]),
            "median_seconds_cand": _med([r["candidate"].get("seconds") for r in rs]),
            "median_steps_ref": _med([r["reference"].get("steps") for r in rs]),
            "median_steps_cand": _med([r["candidate"].get("steps") for r in rs]),
        }
    # (d) top-1 agreement within identity, for structures with ≥2 successful starts
    by_sid: dict[int, list[dict]] = {}
    for r in records:
        if "heavy_rmsd" in r:
            by_sid.setdefault(r["sid"], []).append(r)
    multi = [rs for rs in by_sid.values() if len(rs) >= 2]

    def top(rs, key):
        return min(rs, key=key)["start_geometry"]

    out["ranking"] = {
        "structures": len(multi),
        "top1_cand_energy": sum(top(rs, lambda r: r["reference"]["energy"])
                                == top(rs, lambda r: r["candidate"]["energy"]) for rs in multi),
        "top1_ref_sp_at_cand_min": sum(top(rs, lambda r: r["reference"]["energy"])
                                       == top(rs, lambda r: r["ref_sp_on_candidate"])
                                       for rs in multi),
    }
    return out


if __name__ == "__main__":
    raise SystemExit(main())
