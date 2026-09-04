"""Check a registry against itself: `python scripts/verify_registry.py [--db ...]`.

Exits non-zero if any row contradicts its own graph, so it can gate a commit or a
bulk import.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.config import data_root, registry_path      # noqa: E402
from mofsbu.registry import Registry                    # noqa: E402
from mofsbu.registry.verify import summarise, verify    # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--db", type=Path, default=None)
    p.add_argument("--verbose", "-v", action="store_true")
    a = p.parse_args(argv)

    db = a.db or (registry_path() if registry_path().exists()
                  else data_root() / "demo_registry.db")
    if not db.exists():
        p.error(f"no registry at {db}")

    with Registry(db) as reg:
        problems = verify(reg)
        print(f"{db}  ({reg.count('structures')} structures, {reg.count('geometries')} geometries)")
        print(summarise(problems))
        if a.verbose:
            for problem in problems:
                print("   ", problem)
    return 1 if any(p.severity == "error" for p in problems) else 0


if __name__ == "__main__":
    raise SystemExit(main())
