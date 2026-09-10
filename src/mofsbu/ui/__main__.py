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
import os
import sys
from pathlib import Path

from mofsbu.config import (
    REPO_ROOT, compute_device, data_root, machine_profile, max_workers, registry_path,
    store_root,
)


def discover_databases(root: Path | None = None) -> list[Path]:
    """Every SQLite file under the data root, newest first."""

    base = Path(root or data_root())
    if not base.exists():
        return []
    found = [p for p in base.glob("*.db") if p.is_file()]
    return sorted(found, key=lambda p: p.stat().st_mtime, reverse=True)


def _candidates(given: Path | None, *defaults: Path) -> list[Path]:
    """Where to look, in order.

    A relative path is tried against cwd, then the data root, then the repo root — a
    relative path typed from `src/` should still find `data/` at the top of the checkout
    rather than failing with a bare 'not found'.
    """
    if given is None:
        return list(defaults)
    if given.is_absolute():
        return [given]
    return [given, data_root() / given.name, REPO_ROOT / given]


def resolve_database(given: Path | None) -> tuple[Path, list[Path], str]:
    """Pick a database.  Returns (chosen, others_available, how_it_was_chosen)."""
    if given is not None:
        for candidate in _candidates(given):
            if candidate.exists():
                return candidate, [], "given"
        # An explicitly named database that does not exist is still honoured: the user
        # said which file they meant, and it will be created on first write.
        return _candidates(given)[0], discover_databases(), "given (will be created)"

    discovered = discover_databases()
    preferred = registry_path()
    if preferred.exists():
        return preferred, [p for p in discovered if p != preferred], "the default registry"
    if discovered:
        return discovered[0], discovered[1:], "newest in the data folder"
    return preferred, [], "nothing found yet — will be created on first write"


def confirm_compute_settings() -> None:
    """Ask, once at startup, how this server should use CUDA and in-process workers.

    Ground rule 9 still holds — the device and worker count are DECLARED
    (`MOFSBU_DEVICE`, `MOFSBU_WORKERS`/`MOFSBU_PROFILE`), never auto-detected — this
    just puts that declaration in front of a human at the moment it matters, because
    `viewer.py` is the process that actually runs any `ml_go`/`xtb_go` task submitted
    through `/builder`, using whatever was already (maybe silently) set.  Skipped
    outright without a TTY, so a systemd/docker/nohup launch never blocks on stdin.
    """
    if not sys.stdin.isatty():
        return

    device = compute_device()
    hint = ""
    try:
        import torch

        if torch.cuda.is_available():
            hint = f"  (torch sees a GPU: {torch.cuda.get_device_name(0)})"
    except Exception:                                                    # noqa: BLE001
        pass
    print(f"\ncompute device: currently {device!r}{hint}")
    answer = input("  use CUDA for this server's ml_go/xtb_go runs? [y/N]: ").strip().lower()
    if answer in ("y", "yes"):
        chosen = input("  device [cuda]: ").strip() or "cuda"
        os.environ["MOFSBU_DEVICE"] = chosen
    print(f"worker profile: currently {machine_profile()!r}, max_workers={max_workers()}")
    answer = input("  allow more than one in-process worker (workstation profile)? "
                   "[y/N]: ").strip().lower()
    if answer in ("y", "yes"):
        os.environ["MOFSBU_PROFILE"] = "workstation"
        count = input(f"  how many workers? [{max_workers()}]: ").strip()
        if count:
            os.environ["MOFSBU_WORKERS"] = count
    print(f"→ device={compute_device()}  profile={machine_profile()}  "
          f"max_workers={max_workers()}\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="python -m mofsbu.ui",
                                description=__doc__.splitlines()[0])
    p.add_argument("--db", type=Path, default=None,
                   help="registry SQLite file, opened read-only "
                        "(default: the real registry if present, else the demo one)")
    p.add_argument("--store", type=Path, default=None,
                   help="blob store root holding the .xyz payloads")
    p.add_argument("--list-db", action="store_true",
                   help="list the databases that can be seen, then exit")
    p.add_argument("--host", default="127.0.0.1", help="bind address (default: localhost only)")
    p.add_argument("--port", type=int, default=8000)
    p.add_argument("--reload", action="store_true", help="uvicorn autoreload (development)")
    p.add_argument("--no-prompt", action="store_true",
                   help="skip the CUDA/worker confirmation; use MOFSBU_DEVICE and "
                        "MOFSBU_WORKERS/MOFSBU_PROFILE exactly as already declared")
    a = p.parse_args(argv)

    if a.list_db:
        found = discover_databases()
        print(f"data folder: {data_root()}")
        if not found:
            print("  no databases yet — one is created when you submit a run at /builder")
        for path in found:
            size = path.stat().st_size / 1e6
            print(f"  {path.name:28s} {size:8.2f} MB   {path}")
        return 0

    # `--reload` re-execs this process on every code change; asking again each time
    # would be a prompt loop, not a confirmation.
    if not a.no_prompt and not a.reload:
        confirm_compute_settings()

    db, others, why = resolve_database(a.db)
    store = a.store or store_root()

    print(f"mofsbu  db={db}  ({why})")
    if others:
        print("  also available: " + ", ".join(o.name for o in others)
              + "   — pass --db to choose")
    if not db.exists():
        print("  this file does not exist yet; open /builder and submit a run to create it")
    print(f"  store={store}")
    print(f"  → http://{a.host}:{a.port}/         registry viewer (read-only)")
    print(f"  → http://{a.host}:{a.port}/builder  spec builder")

    import uvicorn

    from mofsbu.ui.app import create_app

    print(f"mofsbu viewer (read-only)  db={db}  store={store}")
    print(f"  → http://{a.host}:{a.port}/")
    uvicorn.run(create_app(db, store), host=a.host, port=a.port, reload=a.reload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
