"""Which relaxer's minimum does DFT prefer?  Stage 5a in miniature (E-MH).

    python scripts/dft_arbitrate.py --db data/salt_study_k1.db --study mh_vs_omol.json \\
        --out arbitration.json --sids 132 68 28 94 --extra 2 --basis def2-svp

For each chosen structure, the start geometry of the study's first record is relaxed again
with both ML backends (same settings as the study), and both minima get a DFT single point
and gradient.  Reported per pair: ΔE_DFT = E(candidate minimum) − E(reference minimum), and
each minimum's largest DFT force.  A minimum that is lower in DFT and nearer a DFT stationary
point is the better geometry; neither number is a reaction energy.  Reads the registry
read-only and writes nothing to it.  Needs a declared CUDA device and the GPU to itself.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import statistics
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from compare_relaxers import heavy_rmsd, relax                           # noqa: E402

from mofsbu.energy.backends import get_backend                           # noqa: E402
from mofsbu.energy.dft import BOHR_ANGSTROM, HARTREE_EV, DFTBackend, _analysis  # noqa: E402
from mofsbu.registry import BlobStore, geometry_xyz, get_graph           # noqa: E402
from mofsbu.registry.db import ReadOnlyRegistry                          # noqa: E402
from mofsbu.runner import read_xyz                                       # noqa: E402


def dft_point(b: DFTBackend, symbols, xyz, charge: int, multiplicity: int) -> dict:
    """Energy (eV) and largest force (eV/Å) at one geometry."""
    t0 = time.perf_counter()
    mol, mf = b._prepare(symbols, xyz, charge, multiplicity, None)
    e = float(mf.kernel(dm0=b._guess(mol, multiplicity))) * HARTREE_EV
    g = np.asarray(mf.nuc_grad_method().kernel())
    g = g.get() if hasattr(g, "get") else g
    fmax = float(np.sqrt((g ** 2).sum(axis=1)).max()) * HARTREE_EV / BOHR_ANGSTROM
    # the electronic state, so two SCF solutions at one geometry are told apart
    state = _analysis(mf, multiplicity)
    return {"energy": e, "fmax": fmax, "converged": bool(mf.converged),
            "s2": state.get("s2"), "spin_populations": state.get("spin_populations"),
            "scf_cycles": state.get("scf_cycles"),
            "seconds": round(time.perf_counter() - t0, 1)}


def pick(records: list[dict], sids: list[int], extra: int) -> list[dict]:
    """The named structures' first record, plus per metal class the `extra` records whose
    reference-scored disagreement is nearest that class's median."""
    ok = [r for r in records if "dE_ref_at_candidate_min" in r]
    chosen = {}
    for sid in sids:
        hit = next((r for r in ok if r["sid"] == sid), None)
        if hit is not None:
            chosen[sid] = hit
    for cls in ("metal_neutral", "metal_charged"):
        pool = [r for r in ok if r["class"] == cls and r["sid"] not in chosen]
        if not pool:
            continue
        med = statistics.median(r["dE_ref_at_candidate_min"] for r in pool)
        for r in sorted(pool, key=lambda r: abs(r["dE_ref_at_candidate_min"] - med))[:extra]:
            if r["sid"] not in chosen:
                chosen[r["sid"]] = r
    return list(chosen.values())


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, required=True)
    p.add_argument("--store", type=Path, default=None)
    p.add_argument("--study", type=Path, required=True, help="compare_relaxers.py output")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--sids", type=int, nargs="*", default=[132, 68, 28, 94])
    p.add_argument("--extra", type=int, default=1, help="median-disagreement picks per class")
    p.add_argument("--xc", default="wb97m-v")
    p.add_argument("--basis", default="def2-svp")
    a = p.parse_args(argv)

    study = json.loads(a.study.read_text())
    con = sqlite3.connect(f"file:{a.db.resolve()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    reg = ReadOnlyRegistry(conn=con, store=BlobStore(a.store or a.db.parent / "store"))

    ref = get_backend(study.get("reference_key", "mace_omol"))
    cand = get_backend(study.get("candidate_key", "mace_mh"),
                       **({"head": study["candidate_head"]} if study.get("candidate_head") else {}))
    dft = DFTBackend(xc=a.xc, basis=a.basis)
    if not dft.available():
        raise SystemExit(f"DFT backend unavailable: {dft.install_hint()}")
    chosen = pick(study["records"], a.sids, a.extra)
    print(f"{len(chosen)} structures; {dft.method}; relaxers {ref.method} vs {cand.method}",
          flush=True)

    rows = []
    for r in chosen:
        graph = get_graph(reg, r["sid"])
        q, m = int(graph.charge), int(graph.multiplicity)
        symbols, coords = read_xyz(geometry_xyz(reg, r["start_geometry"]))
        x0 = np.asarray(coords, dtype=float)
        rr = relax(ref, symbols, x0, graph, study["fmax"], study["steps"])
        rc = relax(cand, symbols, x0, graph, study["fmax"], study["steps"])
        row = {k: r[k] for k in ("sid", "class", "formula", "start_geometry")}
        row["charge"], row["multiplicity"] = q, m
        row["heavy_rmsd"] = heavy_rmsd(symbols, rc["xyz"], rr["xyz"])
        row["ref_min"] = dft_point(dft, symbols, rr["xyz"], q, m)
        row["cand_min"] = dft_point(dft, symbols, rc["xyz"], q, m)
        metal = graph.metals()[0] if graph.metals() else None
        for side in ("ref_min", "cand_min"):
            pops = row[side].pop("spin_populations")
            row[side]["metal_spin_population"] = None if pops is None or metal is None \
                else pops[metal]
        row["ref_xyz"] = np.asarray(rr["xyz"]).round(6).tolist()
        row["cand_xyz"] = np.asarray(rc["xyz"]).round(6).tolist()
        row["symbols"] = list(symbols)
        row["dE_dft_cand_minus_ref"] = row["cand_min"]["energy"] - row["ref_min"]["energy"]
        rows.append(row)
        print(json.dumps({k: row[k] for k in ("sid", "class", "heavy_rmsd",
                                              "dE_dft_cand_minus_ref")}
                         | {"fmax_ref": row["ref_min"]["fmax"],
                            "fmax_cand": row["cand_min"]["fmax"]}), flush=True)

    summary = {
        "n": len(rows),
        "dft_prefers_candidate": sum(r["dE_dft_cand_minus_ref"] < 0 for r in rows),
        "median_dE_dft": statistics.median(r["dE_dft_cand_minus_ref"] for r in rows)
        if rows else None,
        "median_fmax_ref": statistics.median(r["ref_min"]["fmax"] for r in rows)
        if rows else None,
        "median_fmax_cand": statistics.median(r["cand_min"]["fmax"] for r in rows)
        if rows else None,
    }
    a.out.write_text(json.dumps({"dft": dft.method, "reference": ref.method,
                                 "candidate": cand.method, "fmax": study["fmax"],
                                 "summary": summary, "records": rows}, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
