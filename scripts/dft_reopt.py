"""Stage 5a: how far is each ML relaxer's minimum from the DFT minimum it relaxes into?

    MOFSBU_DEVICE=cuda python scripts/dft_reopt.py --db data/salt_study_k1.db \\
        --out reopt.json --sids 132 68 569 28 94 718 678 804 --extra 1 --resume

`docs/WORKPLAN_energy.md` stage 5a.  For each structure, its first stored start geometry
is relaxed by each ML backend (same start, same settings), and each ML minimum is then
re-optimised by DFT (geomeTRIC, default ωB97M-V/def2-TZVPD).  Reported per relaxer:

  * ΔE = E_DFT(ML minimum) − E_DFT(DFT minimum reached from it): the accelerator claim
    holds if its mean is under 0.05 eV (proposal);
  * heavy-atom RMSD between the ML minimum and that DFT minimum;
  * ⟨S²⟩, metal spin population, SCF cycles and the largest DFT force at every point;
  * E1's connectivity verdict on the ML and DFT minima;
  * the DFT optimiser's energy trail, and peak GPU memory per structure.

Reads the registry read-only and writes nothing to it.  The output is rewritten after every
structure; `--resume` keeps the structures already in it.  Needs the GPU to itself.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from compare_relaxers import bond_changes, heavy_rmsd, relax, sample  # noqa: E402
from dft_arbitrate import dft_point  # noqa: E402
from mofsbu.energy.backends import get_backend  # noqa: E402
from mofsbu.energy.dft import DFTBackend, contraction_engine  # noqa: E402
from mofsbu.registry import BlobStore, RegistryError, geometry_graph, geometry_xyz, get_graph  # noqa: E402
from mofsbu.registry.db import ReadOnlyRegistry  # noqa: E402
from mofsbu.runner import read_xyz  # noqa: E402

#: 5a acceptance threshold (proposal): mean ΔE(ML minimum − DFT minimum) below this, eV.
ACCEPT_EV = 0.05


def pick(pool: list[dict], sids: list[int], extra: int, near_atoms: int) -> list[dict]:
    """The named structures, plus per metal class the `extra` structures nearest
    `near_atoms` atoms (the memory figure §3a lacks is for a ~40-atom Ni complex)."""
    by_sid = {e["sid"]: e for e in pool}
    chosen = [by_sid[s] for s in sids if s in by_sid]
    missing = [s for s in sids if s not in by_sid]
    if missing:
        print(f"not in the registry's reference relaxations, skipped: {missing}", flush=True)
    for cls in ("metal_neutral", "metal_charged"):
        rest = [e for e in pool if e["class"] == cls and e["sid"] not in sids]
        chosen += sorted(rest, key=lambda e: (abs(e["n_atoms"] - near_atoms), e["sid"]))[:extra]
    return chosen


class GpuPeak:
    """Peak device memory held by CuPy (GPU4PySCF) and torch (MACE) within one structure."""

    def __init__(self) -> None:
        self.cupy = 0

    def reset(self) -> None:
        import cupy
        import torch

        cupy.get_default_memory_pool().free_all_blocks()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()
        self.cupy = 0

    def sample(self) -> None:
        import cupy

        # the pool keeps what it has allocated, so its size after a stage is that stage's peak
        self.cupy = max(self.cupy, cupy.get_default_memory_pool().total_bytes())

    def report(self) -> dict:
        import torch

        gib = 1024 ** 3
        return {"cupy_gib": round(self.cupy / gib, 2),
                "torch_gib": round(torch.cuda.max_memory_allocated() / gib, 2)
                if torch.cuda.is_available() else None}


def point(dft: DFTBackend, symbols, xyz, q: int, m: int, metal: int | None) -> dict:
    out = dft_point(dft, symbols, xyz, q, m)
    pops = out.pop("spin_populations")
    out["metal_spin_population"] = None if pops is None or metal is None else pops[metal]
    return out


def one_relaxer(dft, backend, symbols, x0, graph, ordered, metal, a, peak) -> dict:
    q, m = int(graph.charge), int(graph.multiplicity)
    ml = relax(backend, symbols, x0, graph, a.fmax, a.steps)
    if "error" in ml:
        return {"error": ml["error"]}
    rec = {"ml": {k: ml[k] for k in ("energy", "converged", "steps", "fmax", "seconds")}}
    if ordered:
        rec["ml"]["connectivity"] = bond_changes(graph, symbols, ml["xyz"])
    rec["ml_xyz"] = np.asarray(ml["xyz"]).round(6).tolist()
    release_torch_cache()                 # MACE's cached blocks are DFT's working memory
    try:
        rec["dft_at_ml_min"] = point(dft, symbols, ml["xyz"], q, m, metal)
    except Exception as exc:                                             # noqa: BLE001
        rec["dft_at_ml_min"] = {"error": f"{type(exc).__name__}: {exc}"}
        return rec
    peak.sample()

    trail: list[list[float]] = []
    t0 = time.perf_counter()

    def on_step(e: float, f: float) -> None:
        trail.append([round(e, 6), round(f, 4)])
        print(f"    step {len(trail):3d}  E {e:.4f} eV  fmax {f:.3f} eV/Å", flush=True)

    try:
        r = dft.relax(symbols, ml["xyz"], charge=q, multiplicity=m, fmax=a.dft_fmax,
                      steps=a.dft_steps, on_step=on_step)
    except Exception as exc:                                             # noqa: BLE001
        rec["dft_relax"] = {"error": f"{type(exc).__name__}: {exc}", "trail": trail}
        return rec
    peak.sample()
    xd = np.asarray(r.positions, dtype=float)
    rec["dft_relax"] = {"converged": r.converged, "steps": r.n_steps,
                        "seconds": round(time.perf_counter() - t0, 1), "trail": trail}
    rec["dft_xyz"] = xd.round(6).tolist()
    try:
        rec["dft_min"] = point(dft, symbols, xd, q, m, metal)
    except Exception as exc:                                             # noqa: BLE001
        rec["dft_min"] = {"error": f"{type(exc).__name__}: {exc}"}
        return rec
    peak.sample()
    if ordered:
        rec["dft_min"]["connectivity"] = bond_changes(graph, symbols, xd)
    rec["heavy_rmsd_ml_dft"] = heavy_rmsd(symbols, np.asarray(ml["xyz"]), xd)
    rec["dE_ml_min_minus_dft_min"] = rec["dft_at_ml_min"]["energy"] - rec["dft_min"]["energy"]
    return rec


def release_torch_cache() -> None:
    import torch

    if torch.cuda.is_available():
        torch.cuda.empty_cache()


#: Two DFT points whose ⟨S²⟩ differ by more than this are on different SCF solutions, and
#: their energy difference is not a geometry error.
S2_TOL = 0.01


def _bonding(conn: dict | None) -> tuple | None:
    if conn is None:
        return None
    return (sorted((i, j) for i, j, _ in conn["broken"]),
            sorted((i, j) for i, j, _ in conn["formed"]))


def excluded(rec: dict) -> str | None:
    """Why a record's ΔE is not a geometry error at all, or None.

    Only cases where the two DFT energies are not comparable: no DFT minimum, an
    unconverged SCF or optimisation, or two different SCF solutions.  A DFT minimum with
    other bonding than the ML minimum is NOT excluded — that is the ML geometry being in
    the wrong basin, which is what 5a measures (`changed_bonding` marks it).
    """
    if "dE_ml_min_minus_dft_min" not in rec:
        return "no_dft_minimum"
    if not (rec["dft_at_ml_min"]["converged"] and rec["dft_min"]["converged"]):
        return "scf_not_converged"
    if not rec["dft_relax"]["converged"]:
        return "dft_relax_not_converged"
    s2a, s2b = rec["dft_at_ml_min"].get("s2"), rec["dft_min"].get("s2")
    if s2a is not None and s2b is not None and abs(s2a - s2b) > S2_TOL:
        return "scf_state_differs"
    return None


def changed_bonding(rec: dict) -> bool:
    """DFT re-optimisation broke or formed a bond relative to the ML minimum."""
    return (_bonding(rec["ml"].get("connectivity"))
            != _bonding(rec.get("dft_min", {}).get("connectivity")))


def clean(rec: dict) -> bool:
    return excluded(rec) is None


def summarise(rows: list[dict], keys: list[str]) -> dict:
    """The 5a acceptance figure is the mean ΔE over every comparable record (proposal: < ACCEPT_EV).
    The same-bonding mean is reported beside it, never instead of it."""
    out = {}
    for key in keys:
        recs = [r["relaxers"][key] for r in rows if key in r["relaxers"]]
        ok = [x for x in recs if clean(x)]
        de = [x["dE_ml_min_minus_dft_min"] for x in ok]
        same = [x["dE_ml_min_minus_dft_min"] for x in ok if not changed_bonding(x)]
        out[key] = {
            "n": len(recs), "n_comparable": len(ok),
            "excluded": dict(Counter(why for x in recs if (why := excluded(x)))),
            "n_changed_bonding": len(de) - len(same),
            "mean_dE": statistics.fmean(de) if de else None,
            "mean_dE_same_bonding": statistics.fmean(same) if same else None,
            "median_dE": statistics.median(de) if de else None,
            "max_dE": max(de) if de else None,
            "median_rmsd": statistics.median(x["heavy_rmsd_ml_dft"] for x in ok) if ok else None,
            "within_acceptance": (statistics.fmean(de) < ACCEPT_EV) if de else None,
        }
    return out


def write(a, dft, backends, rows) -> dict:
    keys = list(backends)
    summary = summarise(rows, keys)
    a.out.write_text(json.dumps({
        "dft": dft.method, "contraction_engine": contraction_engine(),
        "relaxers": {k: b.method for k, b in backends.items()},
        "fmax": a.fmax, "steps": a.steps, "dft_fmax": a.dft_fmax, "dft_steps": a.dft_steps,
        "acceptance_eV": ACCEPT_EV, "summary": summary, "records": rows}, indent=1))
    return summary


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, required=True)
    p.add_argument("--store", type=Path, default=None,
                   help="blob store (default: <db dir>/store)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--relaxers", nargs="+", default=["mace_omol", "mace_polar_l"],
                   help="ML backend keys; the first one's stored relaxations supply the starts")
    p.add_argument("--sids", type=int, nargs="*", default=[132, 68, 569, 28, 94, 718, 678, 804])
    p.add_argument("--extra", type=int, default=1, help="added structures per metal class")
    p.add_argument("--near-atoms", type=int, default=40)
    p.add_argument("--fmax", type=float, default=0.05, help="ML relaxation, eV/Å")
    p.add_argument("--steps", type=int, default=500)
    p.add_argument("--dft-fmax", type=float, default=0.05, help="DFT re-optimisation, eV/Å")
    p.add_argument("--dft-steps", type=int, default=150)
    p.add_argument("--xc", default="wb97m-v")
    p.add_argument("--basis", default="def2-tzvpd")
    p.add_argument("--resume", action="store_true",
                   help="keep the structures already in --out and skip them")
    a = p.parse_args(argv)

    con = sqlite3.connect(f"file:{a.db.resolve()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    reg = ReadOnlyRegistry(conn=con, store=BlobStore(a.store or a.db.parent / "store"))

    backends = {key: get_backend(key) for key in a.relaxers}
    dft = DFTBackend(xc=a.xc, basis=a.basis)
    if not dft.available():
        raise SystemExit(f"DFT backend unavailable: {dft.install_hint()}")
    ref_method = next(iter(backends.values())).method
    pool = sample(con, ref_method, per_class=10 ** 6, conformers=1, seed=0)
    chosen = pick(pool, a.sids, a.extra, a.near_atoms)

    rows: list[dict] = []
    if a.resume and a.out.is_file():
        rows = json.loads(a.out.read_text())["records"]
        # a structure whose DFT raised (out of memory, say) is run again, not kept
        rows = [r for r in rows if set(r["relaxers"]) >= set(backends)
                and not any(excluded(x) == "no_dft_minimum" for x in r["relaxers"].values())]
    done = {r["sid"] for r in rows}
    print(f"{len(chosen)} structures ({len(done)} done); {dft.method}; relaxers "
          f"{[b.method for b in backends.values()]}", flush=True)

    peak = GpuPeak()
    for e in chosen:
        if e["sid"] in done:
            continue
        gid = e["starts"][0]
        symbols, coords = read_xyz(geometry_xyz(reg, gid))
        # connectivity is judged only when the graph is numbered like the geometry (D30)
        try:
            graph, ordered = geometry_graph(reg, gid), True
        except (RegistryError, sqlite3.OperationalError):
            graph, ordered = get_graph(reg, e["sid"]), False
        metals = graph.metals()
        metal = metals[0] if metals else None
        x0 = np.asarray(coords, dtype=float)
        print(f"sid {e['sid']} ({e['class']}, {e['formula']}, {len(symbols)} atoms, "
              f"q={graph.charge} m={graph.multiplicity})", flush=True)
        peak.reset()
        row = {"sid": e["sid"], "class": e["class"], "formula": e["formula"],
               "n_atoms": len(symbols), "start_geometry": gid, "atom_order_checked": ordered,
               "charge": int(graph.charge), "multiplicity": int(graph.multiplicity),
               "symbols": list(symbols), "relaxers": {}}
        for key, backend in backends.items():
            print(f"  {backend.method}", flush=True)
            rec = one_relaxer(dft, backend, symbols, x0, graph, ordered, metal, a, peak)
            row["relaxers"][key] = rec
            print("  " + json.dumps({k: rec.get(k) for k in (
                "heavy_rmsd_ml_dft", "dE_ml_min_minus_dft_min", "error")}), flush=True)
        dmins = [r["dft_min"]["energy"] for r in row["relaxers"].values() if "dft_min" in r]
        if len(dmins) == len(backends) == 2:
            row["dE_dft_min_second_minus_first"] = dmins[1] - dmins[0]
        row["gpu_peak"] = peak.report()
        rows.append(row)
        write(a, dft, backends, rows)       # after every structure: a killed run keeps its data

    print(json.dumps(write(a, dft, backends, rows), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
