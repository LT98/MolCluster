"""Add a continuum correction to every converged ML energy in a registry that lacks one.

    python scripts/solvate.py --db data/salt_study_k1.db
    python scripts/solvate.py --db data/salt_study_k1.db --medium alpb:water --limit 20
    python scripts/solvate.py --db data/campaigns/x/x.db --medium alpb:dmf --workers 32

xTB (tblite) single points with and without ALPB on each stored geometry; the difference is
stored as that geometry's dG_solv (WORKPLAN_solvation C17).  Resumable — a geometry already
corrected in the medium is skipped.  A `pathways` run does this itself when it finishes, for
each of its spec's `solvation_media` (`runner.finalise_run`); this is for a registry built
before that, or for a medium the run did not ask for.

`--workers N` (a declaration, ground rule 9) splits the geometries over N processes, each
with its share of the cores (`runner.finalise_env`), committing after every geometry so no
worker holds the write lock across an xTB call.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.config import data_root, store_root                         # noqa: E402
from mofsbu.energy.solvation import MEDIA, correct_registry              # noqa: E402
from mofsbu.registry import BlobStore, Registry                         # noqa: E402


def _shard(db: str, store: str, medium: str, mine: list[int],
           env: dict[str, str]) -> dict:
    """One worker: the geometries it was handed, each committed as soon as it is done.

    The parent takes the list once and slices it: a worker that queried for itself would
    see a list already shortened by the others' commits, and skip what fell between (#98).
    """
    import os

    os.environ.update(env)                    # before tblite starts its OpenMP threads
    from mofsbu._types import MofsbuError
    from mofsbu.energy.solvation import correct

    written, refused = 0, {}
    with Registry(Path(db), BlobStore(Path(store))) as reg:
        for gid in mine:
            try:
                if correct(reg, gid, medium).created:
                    written += 1
            except (MofsbuError, ValueError, RuntimeError) as exc:
                why = f"{type(exc).__name__}: {str(exc).splitlines()[0][:120]}"
                refused[why] = refused.get(why, 0) + 1
            reg.conn.commit()
    return {"candidates": len(mine), "written": written, "refused": refused}


def correct_in_parallel(db: Path, store: Path, medium: str, workers: int) -> dict:
    import multiprocessing as mp

    from mofsbu.energy.solvation import uncorrected
    from mofsbu.runner import finalise_env

    with Registry(db, BlobStore(store)) as reg:
        todo = uncorrected(reg, medium)
    env = finalise_env(workers)
    with mp.get_context("spawn").Pool(workers) as pool:
        parts = pool.starmap(_shard, [(str(db), str(store), medium, todo[i::workers], env)
                                      for i in range(workers)])
    refused: dict[str, int] = {}
    for part in parts:
        for why, k in part["refused"].items():
            refused[why] = refused.get(why, 0) + k
    return {"medium": medium, "candidates": sum(p["candidates"] for p in parts),
            "written": sum(p["written"] for p in parts), "refused": refused}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, default=data_root() / "registry.db")
    p.add_argument("--store", type=Path, default=None)
    p.add_argument("--medium", default="alpb:water", choices=sorted(MEDIA))
    p.add_argument("--limit", type=int, default=None, help="serial only: stop after N")
    p.add_argument("--workers", type=int, default=1,
                   help="processes to split the geometries over (default 1, in-process)")
    a = p.parse_args(argv)
    if a.workers > 1 and a.limit is not None:
        p.error("--limit is for a serial trial; drop it or use --workers 1")

    store = a.store or store_root()
    started = time.perf_counter()
    with Registry(a.db, BlobStore(store)) as reg:
        reg.migrate()
        if a.workers <= 1:
            out = correct_registry(reg, a.medium, limit=a.limit)
    if a.workers > 1:
        out = correct_in_parallel(a.db, store, a.medium, a.workers)
    took = time.perf_counter() - started
    print(f"{a.db}: {out['written']} of {out['candidates']} geometries corrected in "
          f"{out['medium']} ({took:.1f} s)")
    for reason, n in sorted(out["refused"].items(), key=lambda kv: -kv[1]):
        print(f"  [{n}] {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
