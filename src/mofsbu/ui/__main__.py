"""``python -m mofsbu.ui`` — serve the read-only registry viewer.

    python -m mofsbu.ui                      # finds the registry on its own
    python -m mofsbu.ui --db some/other.db --port 8001

Paths are resolved against the data root and the repo, not just the working directory,
so the command behaves the same wherever you run it from.  Requires the package to be
importable — `pip install -e .` once per machine, or use `python scripts/viewer.py`,
which needs no install.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from mofsbu.config import REPO_ROOT, data_root, registry_path, store_root


def _candidates(given: Path | None, *defaults: Path) -> list[Path]:
    """Where to look, in order.  A relative path is tried against cwd, then the data
    root, then the repo root — a relative path typed from `src/` should still find
    `data/` at the top of the checkout rather than failing with a bare 'not found'."""
    if given is None:
        return list(defaults)
    if given.is_absolute():
        return [given]
    return [given, data_root() / given.name, REPO_ROOT / given]


def _first_existing(paths: list[Path]) -> Path | None:
    return next((p for p in paths if p.exists()), None)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="python -m mofsbu.ui",
                                description=__doc__.splitlines()[0])
    p.add_argument("--db", type=Path, default=None,
                   help="registry SQLite file, opened read-only "
                        "(default: the real registry if present, else the demo one)")
    p.add_argument("--store", type=Path, default=None,
                   help="blob store root holding the .xyz payloads")
    p.add_argument("--host", default="127.0.0.1", help="bind address (default: localhost only)")
    p.add_argument("--port", type=int, default=8000)
    p.add_argument("--reload", action="store_true", help="uvicorn autoreload (development)")
    a = p.parse_args(argv)

    db_candidates = _candidates(a.db, registry_path(), data_root() / "demo_registry.db")
    db = _first_existing(db_candidates)
    if db is None:
        looked = "\n  ".join(str(c) for c in db_candidates)
        p.error(
            "no registry database found. Looked in:\n  " + looked +
            "\n\nSeed the demo one with:  python scripts/seed_demo_registry.py"
        )

    store = _first_existing(_candidates(a.store, store_root(), db.parent / "store")) \
        or store_root()

    import uvicorn

    from mofsbu.ui.app import create_app

    print(f"mofsbu viewer (read-only)  db={db}  store={store}")
    print(f"  → http://{a.host}:{a.port}/")
    uvicorn.run(create_app(db, store), host=a.host, port=a.port, reload=a.reload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
