"""Relax every conformer with a small model, refine only its top-k with the large one: how much
is lost against relaxing every conformer with the large model?

    python scripts/funnel_eval.py --db data/salt_study_k1.db --out funnel.json \\
        --small mace_polar_s mace_polar_m --large mace_polar_l --structures 12 --starts 10

Per structure (an identity with many raw starts): every start relaxed directly with `--large`
(the baseline), and with each `--small`; the small model's k lowest minima are then relaxed
with `--large` from where the small model left them.  Reported for k = 1..--k-max:

  regret   min E_large over the funnel's refined minima − min E_large over the baseline
           (> 0: the funnel missed a lower minimum; < 0: it found one the baseline did not)
  seconds  all small relaxations + k large refinements, against all large relaxations
  stored   geometries kept: one per start + k, against one per start

All energies compared are the large model's, so regret is in one theory.  Reads the registry
read-only and writes nothing to it.

`--trajectories DIR` (funnel F0, `docs/WORKPLAN_energy.md` §3a) also writes one `.npz` per
relaxation: energy and fmax at every optimiser point, coordinates every `--every` points plus
the last, and the run's identity (structure, start, stage, model).  F1 replays these offline.
`--sample-from raw` takes starts from raw constructs, so a registry with no relaxations yet works.
`--shard I/N` splits the structures over N processes, one card each, and `--merge` joins them
(`scripts/_shard.py`).
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
sys.path.insert(0, str(ROOT / "scripts"))

import _shard                                                            # noqa: E402
from compare_relaxers import relax                                       # noqa: E402

from mofsbu.energy.backends import get_backend                           # noqa: E402
from mofsbu.registry import BlobStore, geometry_xyz, get_graph           # noqa: E402
from mofsbu.registry.db import ReadOnlyRegistry                          # noqa: E402
from mofsbu.runner import read_xyz                                       # noqa: E402


def pick(con, method: str, n: int, starts: int, seed: int,
         classes: tuple[str, ...] = ("neutral", "charged")) -> list[dict]:
    """Ni structures with at least `starts` raw starts, `n` split evenly over `classes`.

    A start is a raw construct that `method` relaxed; with `method == "raw"`, any raw construct.
    """
    if method == "raw":
        rows = con.execute(
            """SELECT s.id sid, s.net_charge, g.id gid FROM geometries g
               JOIN structures s ON s.id = g.structure_id
               WHERE g.fidelity <= 1 AND s.hidden = 0 AND s.n_metals > 0
               ORDER BY s.id, g.id""").fetchall()
    else:
        rows = con.execute(
            """SELECT s.id sid, s.net_charge, r.id gid FROM geometries g
               JOIN methods m ON m.id = g.method_id JOIN geometries r ON r.id = g.relaxed_from
               JOIN structures s ON s.id = g.structure_id
               WHERE m.method = ? AND r.fidelity <= 1 AND s.hidden = 0 AND s.n_metals > 0
               ORDER BY s.id, r.id""", (method,)).fetchall()
    by: dict[int, dict] = {}
    for r in rows:
        e = by.setdefault(r["sid"], {"sid": r["sid"], "charge": r["net_charge"], "starts": []})
        if r["gid"] not in e["starts"]:
            e["starts"].append(r["gid"])
    rng = random.Random(seed)
    out = []
    for charged in (False, True):
        pool = [e for e in by.values() if len(e["starts"]) >= starts
                and (e["charge"] != 0) == charged]
        rng.shuffle(pool)
        if ("charged" if charged else "neutral") in classes:
            out.extend(pool[:n // len(classes)])
    for e in out:
        e["starts"] = rng.sample(e["starts"], starts)
    return out


class Trajectories:
    """Writes each relaxation's path to `<dir>/<sid>_<start>_<stage>_<model>.npz` (F0)."""

    def __init__(self, out: Path | None, every: int) -> None:
        self.out, self.every = out, every
        if out is not None:
            out.mkdir(parents=True, exist_ok=True)

    def relax(self, backend, key: str, stage: str, sid: int, gid: int, symbols, xyz, graph,
              fmax: float, steps: int) -> dict:
        if self.out is None:
            return relax(backend, symbols, xyz, graph, fmax, steps)
        energy: list[float] = []
        force: list[float] = []
        frames: list[tuple[int, np.ndarray]] = []

        def on_step(e: float, f: float, x) -> None:
            if len(energy) % self.every == 0:
                frames.append((len(energy), np.asarray(x, dtype=float)))
            energy.append(e)
            force.append(f)
            last[0] = np.asarray(x, dtype=float)

        last: list = [None]
        r = relax(backend, symbols, xyz, graph, fmax, steps, on_step=on_step)
        if energy and frames[-1][0] != len(energy) - 1:
            frames.append((len(energy) - 1, last[0]))
        meta = {"sid": sid, "start_geometry": gid, "stage": stage, "model": key,
                "method": backend.method, "charge": int(graph.charge),
                "multiplicity": int(graph.multiplicity), "fmax_target": fmax,
                "steps_max": steps, "every": self.every, "symbols": list(symbols),
                "converged": r.get("converged"), "error": r.get("error")}
        np.savez_compressed(
            self.out / f"{sid}_{gid}_{stage}_{key}.npz",
            energy=np.array(energy), fmax=np.array(force),
            frame_index=np.array([i for i, _ in frames], dtype=int),
            coords=np.array([x for _, x in frames]).reshape(len(frames), len(symbols), 3),
            start=np.asarray(xyz, dtype=float), meta=json.dumps(meta))
        return r


def summarise(results: list[dict], smalls: list[str], k_max: int, starts: int) -> dict:
    """Regret, time and storage per small model and k, over every structure in `results`."""
    summary: dict = {}
    base_s = sum(r["baseline"]["seconds"] for r in results)
    for key in smalls:
        summary[key] = {}
        for k in range(1, k_max + 1):
            regrets = [r[key]["per_k"][k]["regret"] for r in results
                       if r[key]["per_k"][k]["regret"] is not None]
            summary[key][f"k={k}"] = {
                "median_regret": round(statistics.median(regrets), 4) if regrets else None,
                "worst_regret": round(max(regrets), 4) if regrets else None,
                "missed_over_0.05eV": sum(x > 0.05 for x in regrets),
                "time_vs_all_large": round(sum(r[key]["per_k"][k]["seconds"] for r in results)
                                           / base_s, 3) if base_s else None,
                "stored_vs_all_large": round((starts + k) / starts, 2)}
        summary[key]["rank_of_large_best"] = [r[key]["rank_of_large_best"] for r in results]
    return summary


def merge(a: argparse.Namespace) -> int:
    """Join `--shard` outputs into `--out` (`_shard`), with the summary over all of them."""
    head, results = _shard.load(a.merge, "results",
                                ("large", "small", "fmax", "starts", "k_max", "seed"))
    for r in results:                          # JSON turned the k keys into strings
        for key in head["small"]:
            r[key]["per_k"] = {int(k): v for k, v in r[key]["per_k"].items()}
    head = {k: v for k, v in head.items() if k not in ("summary", "results")}
    summary = summarise(results, head["small"], head["k_max"], head["starts"])
    a.out.write_text(json.dumps(head | {"summary": summary, "results": results}, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, default=None, help="required unless --merge")
    p.add_argument("--store", type=Path, default=None)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--small", nargs="+", default=["mace_polar_s", "mace_polar_m"])
    p.add_argument("--large", default="mace_polar_l")
    p.add_argument("--structures", type=int, default=12)
    p.add_argument("--starts", type=int, default=10)
    p.add_argument("--k-max", type=int, default=3)
    p.add_argument("--fmax", type=float, default=0.05)
    p.add_argument("--steps", type=int, default=1000)
    p.add_argument("--sample-from", default="MACE-OMOL-0",
                   help="method whose stored relaxations define the starts; 'raw' for any "
                        "raw construct")
    p.add_argument("--trajectories", type=Path, default=None,
                   help="directory for one .npz per relaxation (F0 trajectory capture)")
    p.add_argument("--every", type=int, default=10,
                   help="with --trajectories: keep coordinates every N optimiser points")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--classes", nargs="+", default=["neutral", "charged"],
                   choices=["neutral", "charged"])
    p.add_argument("--skip", type=int, nargs="*", default=[],
                   help="structure ids already measured by an earlier run")
    _shard.add_args(p)
    a = p.parse_args(argv)
    if a.merge:
        return merge(a)
    if a.db is None:
        p.error("--db is required")
    shard = _shard.start(a)

    con = sqlite3.connect(f"file:{a.db.resolve()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    reg = ReadOnlyRegistry(conn=con, store=BlobStore(a.store or a.db.parent / "store"))
    large = get_backend(a.large)
    smalls = {k: get_backend(k) for k in a.small}
    picked = _shard.take(
        pick(con, a.sample_from, a.structures, a.starts, a.seed, tuple(a.classes)), shard)
    header = {"large": large.method, "small": list(smalls), "fmax": a.fmax,
              "starts": a.starts, "k_max": a.k_max, "seed": a.seed}
    traj = Trajectories(a.trajectories, a.every)
    print(f"{len(picked)} structures x {a.starts} starts; large {large.method}, small "
          f"{[b.method for b in smalls.values()]}", flush=True)

    results = []
    for e in picked:
        if e["sid"] in a.skip:
            continue
        graph = get_graph(reg, e["sid"])
        starts = []
        for gid in e["starts"]:
            symbols, coords = read_xyz(geometry_xyz(reg, gid))
            starts.append((gid, symbols, np.asarray(coords, dtype=float)))
        base = [traj.relax(large, a.large, "base", e["sid"], g, s, x, graph, a.fmax, a.steps)
                for g, s, x in starts]
        ok = [r for r in base if "energy" in r]
        if not ok:
            continue
        best = min(r["energy"] for r in ok)
        row = {"sid": e["sid"], "charge": e["charge"], "n_atoms": len(starts[0][1]),
               "baseline": {"seconds": sum(r.get("seconds", 0) for r in base),
                            "stored": len(base), "failed": len(base) - len(ok)}}
        for key, small in smalls.items():
            mins = [traj.relax(small, key, "small", e["sid"], g, s, x, graph, a.fmax, a.steps)
                    for g, s, x in starts]
            t_small = sum(r.get("seconds", 0) for r in mins)
            order = sorted((i for i, r in enumerate(mins) if "energy" in r),
                           key=lambda i: mins[i]["energy"])
            refined: list[float] = []
            t_refine = 0.0
            per_k = {}
            for k in range(1, a.k_max + 1):
                if k <= len(order):
                    i = order[k - 1]
                    r = traj.relax(large, a.large, f"refine-{key}", e["sid"], starts[i][0],
                                   starts[i][1], mins[i]["xyz"], graph, a.fmax, a.steps)
                    t_refine += r.get("seconds", 0)
                    if "energy" in r:
                        refined.append(r["energy"])
                per_k[k] = {"regret": (min(refined) - best) if refined else None,
                            "seconds": t_small + t_refine, "stored": len(starts) + k}
            # where the large model's own best lands in the small model's ranking
            best_i = min((i for i, r in enumerate(base) if "energy" in r),
                         key=lambda i: base[i]["energy"])
            row[key] = {"per_k": per_k, "seconds_small": t_small,
                        "rank_of_large_best": order.index(best_i) + 1 if best_i in order else None}
        results.append(row)
        print(json.dumps({"sid": row["sid"], "charge": row["charge"],
                          "baseline_s": round(row["baseline"]["seconds"], 1)}
                         | {k: {kk: {"regret": None if v["regret"] is None else round(v["regret"], 3),
                                     "s": round(v["seconds"], 1)}
                                for kk, v in row[k]["per_k"].items()} for k in smalls}), flush=True)
        a.out.write_text(json.dumps(header | {"results": results}, indent=1))

    summary = summarise(results, list(smalls), a.k_max, a.starts)
    a.out.write_text(json.dumps(header | {"summary": summary, "results": results}, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
