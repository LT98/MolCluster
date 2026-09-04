"""Run a build spec.  This is the input interface: edit a spec file, run it.

    python scripts/run_spec.py data/reference/spec_zn_thq.json
    python scripts/run_spec.py my_spec.json --db data/my_run.db --status

The spec is data, so the same file can be run here, queued on the workstation, or
re-run later to reproduce a result.  Nothing about the interface changes when the
execution does — see `mofsbu/runner.py`.

Parallelism is opt-in (ground rule 8).  This machine runs ONE in-process worker unless
you declare it a workstation:

    set MOFSBU_PROFILE=workstation        (windows)
    export MOFSBU_PROFILE=workstation     (linux)
    export MOFSBU_WORKERS=8               (or say exactly how many)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.config import data_root, machine_profile, max_workers      # noqa: E402
from mofsbu.registry import BlobStore, Registry, find                  # noqa: E402
from mofsbu.registry.jobs import iter_tasks                            # noqa: E402
from mofsbu.runner import run                                          # noqa: E402
from mofsbu.spec import BuildSpec                                      # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("spec", type=Path, help="the build spec (JSON)")
    p.add_argument("--db", type=Path, default=data_root() / "runs.db")
    p.add_argument("--store", type=Path, default=data_root() / "store")
    p.add_argument("--workers", type=int, default=None, help="override the worker count")
    p.add_argument("--status", action="store_true", help="list every task and its outcome")
    a = p.parse_args(argv)

    spec = BuildSpec.load(a.spec)
    print(f"spec {a.spec}  digest {spec.digest[:12]}  degree {spec.degree}")
    print(f"  machine profile: {machine_profile()}  workers: "
          f"{a.workers if a.workers is not None else max_workers()}")

    with Registry(a.db, BlobStore(a.store)) as reg:
        reg.migrate("run_spec")
        summary = run(reg, spec, workers=a.workers)
        print(f"\n  run {summary['run_id']}: {summary['tasks']} tasks -> {summary['counts']}")
        if a.status:
            for task in iter_tasks(reg, summary["run_id"]):
                mark = {"done": "ok  ", "rejected": "no  ", "failed": "FAIL"}.get(
                    task["status"], task["status"])
                first = ((task["error"] or "").strip().splitlines() or [""])[0]
                print(f"    {mark} {task['kind']:8s} #{task['id']:<4d} {first[:90]}")
        print()
        for row in find(reg, limit=200):
            print("   ", row["display_label"])
    return 0 if summary["status"] == "done" else 1


if __name__ == "__main__":
    raise SystemExit(main())
