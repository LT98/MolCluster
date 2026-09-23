"""Add a continuum correction to every converged ML energy in a registry that lacks one.

    python scripts/solvate.py --db data/salt_study_k1.db
    python scripts/solvate.py --db data/salt_study_k1.db --medium alpb:water --limit 20

xTB (tblite) single points with and without ALPB on each stored geometry; the difference is
stored as that geometry's dG_solv (WORKPLAN_solvation C17).  Resumable — a geometry already
corrected in the medium is skipped.  A `pathways` run does this itself when it finishes
(`runner.finalise_run`); this is for a registry built before that.
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


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, default=data_root() / "registry.db")
    p.add_argument("--store", type=Path, default=None)
    p.add_argument("--medium", default="alpb:water", choices=sorted(MEDIA))
    p.add_argument("--limit", type=int, default=None)
    a = p.parse_args(argv)

    started = time.perf_counter()
    with Registry(a.db, BlobStore(a.store or store_root())) as reg:
        reg.migrate()
        out = correct_registry(reg, a.medium, limit=a.limit)
    took = time.perf_counter() - started
    print(f"{a.db}: {out['written']} of {out['candidates']} geometries corrected in "
          f"{out['medium']} ({took:.1f} s)")
    for reason, n in sorted(out["refused"].items(), key=lambda kv: -kv[1]):
        print(f"  [{n}] {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
