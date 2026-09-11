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
    """Ask, at startup, how this server should use CUDA and in-process workers.

    No longer the default, and that is the point.  The prompt made the device VISIBLE,
    which was the original problem, but it did not make it REACHABLE: answering it
    requires a terminal, which most of the people who use the page do not have open.
    `/builder` now carries the same control (`POST /api/compute`), so this is the
    second-best way to answer the same question and is opt-in behind `--prompt-device`.

    Ground rule 9 is unchanged either way: the device and worker count are DECLARED
    (`MOFSBU_DEVICE`, `MOFSBU_WORKERS`/`MOFSBU_PROFILE`), never auto-detected.  Skipped
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
    p.add_argument("--prompt-device", action="store_true",
                   help="ask about CUDA and workers at startup. Off by default: the "
                        "same control is on /builder, which does not need a terminal")
    p.add_argument("--no-prompt", action="store_true",
                   help=argparse.SUPPRESS)     # kept so old command lines keep working
    p.add_argument("--device", default=None,
                   help="declare the compute device for this server (cpu, cuda, cuda:1, "
                        "mps). Changeable afterwards on /builder")
    p.add_argument("--workers", type=int, default=None,
                   help="declare how many in-process workers a run may use (default 1)")
    p.add_argument("--log-level", default="warning",
                   choices=["critical", "error", "warning", "info", "debug"],
                   help="uvicorn log level. Defaults to 'warning', so an open run "
                        "inspector polling every few seconds does not fill the console "
                        "with access lines; pass 'info' to get them back")
    p.add_argument("--open", action="store_true",
                   help="open the viewer in a browser once the server is up")
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

    # A flag is a declaration too, and it is the one a launcher script can make on the
    # user's behalf.  Applied before anything reads the environment.
    if a.device or a.workers:
        from mofsbu.ui.active import declare_compute

        try:
            declare_compute(a.device, a.workers)
        except ValueError as exc:
            print(f"error: {exc}")
            return 2

    # `--reload` re-execs this process on every code change; asking again each time
    # would be a prompt loop, not a confirmation.
    if a.prompt_device and not a.reload:
        confirm_compute_settings()

    db, others, why = resolve_database(a.db)
    store = a.store or store_root()

    from mofsbu.ui.active import ActiveDatabase, compute_state

    active = ActiveDatabase(db)
    state = compute_state()

    print(f"mofsbu  db={db}  ({why})")
    if others:
        print("  also available: " + ", ".join(o.name for o in others)
              + "   — choose one on /builder, or pass --db")
    if not db.exists():
        print("  this file does not exist yet; open /builder and submit a run to create it")
    print(f"  store={store}")
    print(f"  device={state['device']}  workers={state['workers']}"
          f"  ({len(state['devices'])} device(s) visible — change on /builder)")
    print(f"  → http://{a.host}:{a.port}/         registry viewer (read-only)")
    print(f"  → http://{a.host}:{a.port}/builder  spec builder")
    print(f"  → http://{a.host}:{a.port}/runs     run inspector")

    import uvicorn

    from mofsbu.ui.app import create_app

    if a.open:
        _open_browser_when_up(a.host, a.port)
    uvicorn.run(create_app(db, store, active), host=a.host, port=a.port,
                reload=a.reload, log_level=a.log_level)
    return 0


def _open_browser_when_up(host: str, port: int, timeout: float = 30.0) -> None:
    """Open a browser once the port answers, on a daemon thread.

    Polling the port rather than sleeping a fixed second: importing torch or rdkit can
    make startup take a while on a cold cache, and a browser tab that opens onto a
    connection-refused page is worse than one that opens two seconds later.
    """
    import socket as _socket
    import threading
    import time
    import webbrowser

    target = "127.0.0.1" if host in ("0.0.0.0", "::") else host

    def wait() -> None:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            with _socket.socket() as s:
                s.settimeout(0.5)
                if s.connect_ex((target, port)) == 0:
                    webbrowser.open(f"http://{target}:{port}/")
                    return
            time.sleep(0.3)

    threading.Thread(target=wait, daemon=True, name="mofsbu-open-browser").start()


if __name__ == "__main__":
    raise SystemExit(main())
