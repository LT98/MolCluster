"""Which registry the running server is pointed at, and which device it will use.

Both used to be decided once, on the command line, and then baked into the process:
``create_app(db, store)`` closed over a path and ``confirm_compute_settings()`` asked
about CUDA at startup.  That put two ordinary choices behind a shell, which is the one
tool most of the people using this page do not have.

The holder below is the smallest thing that makes them changeable at runtime without
inventing a new process model:

* **One mutable path, read at request time.**  The viewer, the builder and the run
  inspector all ask this object where the registry is, so switching it moves all three
  together.  A page that showed runs from one database and structures from another
  would be worse than no switch at all.
* **Read-only stays read-only.**  This object holds a *path*, never a connection.  The
  viewer keeps opening its own ``mode=ro`` connections from it (`ui/app.py`) and the
  builder keeps opening its own writable ones (`ui/builder.py`).  Switching the path
  does not widen anybody's access.
* **A run keeps the database it was submitted against.**  `submit_run` snapshots
  ``path`` before starting its thread, so switching mid-run cannot redirect a build
  that is already executing into a different file.

The device follows ground rule 9 unchanged: **declared, never detected.**  Listing the
GPUs a machine has and letting a person pick one is still a declaration; defaulting to
CUDA because a card exists is not, and is exactly what the startup prompt was added to
prevent.
"""
from __future__ import annotations

import os
import threading
from pathlib import Path

from mofsbu.config import compute_device, machine_profile, max_workers


class ActiveDatabase:
    """The registry path the whole app currently points at."""

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)
        self._lock = threading.Lock()

    @property
    def path(self) -> Path:
        # Reading a reference is atomic under the GIL; the lock exists for the swap.
        return self._path

    def switch(self, path: Path | str) -> Path:
        with self._lock:
            self._path = Path(path)
        return self._path

    def __fspath__(self) -> str:                       # so Path(active) just works
        return str(self._path)

    def __str__(self) -> str:
        return str(self._path)


def _validate_name(name: str) -> str:
    """A database name is a path COMPONENT, not a path.

    Same rule `save_spec` applies to spec filenames, for the same reason: the value
    arrives from a browser and is about to be joined onto the data root.
    """
    name = (name or "").strip()
    if not name:
        raise ValueError("give the database a name")
    if "/" in name or "\\" in name or ".." in name or name.startswith("."):
        raise ValueError(f"{name!r} must be a plain file name, not a path")
    if not name.endswith(".db"):
        name += ".db"
    stem = name[:-3]
    if not stem or not all(c.isalnum() or c in "-_." for c in stem):
        raise ValueError(f"{name!r}: use letters, digits, '-', '_' and '.' only")
    return name


def available_devices() -> list[dict[str, object]]:
    """Devices this machine could be told to use — enumerated, never chosen.

    ``cpu`` is always present and always first.  Everything else is reported only if
    torch is installed AND can see it, with the card's own name attached, because
    "cuda:1" on a two-GPU workstation is not a self-explanatory string.  A machine with
    no torch reports cpu alone rather than offering a device that would fail at the
    first structure of an hour-long run.
    """
    devices: list[dict[str, object]] = [
        {"id": "cpu", "label": "CPU", "note": "always available"},
    ]
    try:
        import torch
    except Exception:                                                    # noqa: BLE001
        devices[0]["note"] = "torch is not installed, so the ML rung cannot run at all"
        return devices
    try:
        if torch.cuda.is_available():
            for i in range(torch.cuda.device_count()):
                try:
                    name = torch.cuda.get_device_name(i)
                except Exception:                                        # noqa: BLE001
                    name = "CUDA device"
                devices.append({"id": f"cuda:{i}", "label": f"{name} (cuda:{i})",
                                "note": "NVIDIA GPU"})
    except Exception:                                                    # noqa: BLE001
        pass
    try:
        if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
            devices.append({"id": "mps", "label": "Apple GPU (mps)", "note": "Metal"})
    except Exception:                                                    # noqa: BLE001
        pass
    return devices


def compute_state() -> dict[str, object]:
    """What is declared right now, and what could be declared instead."""
    try:
        declared = compute_device()
    except ValueError as exc:
        # A bad MOFSBU_DEVICE must not take the page down with it; say so instead.
        return {"devices": available_devices(), "device": None, "device_error": str(exc),
                "profile": machine_profile(), "workers": max_workers(),
                "cpu_count": os.cpu_count() or 1}
    return {
        "devices": available_devices(),
        "device": declared,
        "device_error": None,
        "profile": machine_profile(),
        "workers": max_workers(),
        "cpu_count": os.cpu_count() or 1,
    }


def declare_compute(device: str | None = None, workers: int | None = None) -> dict[str, object]:
    """Record a device / worker choice for this process.

    Writes the same environment variables the command line and the startup prompt
    write, because `compute_device()` and `max_workers()` read the environment at
    EXECUTION time and `submit_run` executes on a thread in this very process.  That is
    what makes a page control possible here without a new process model — and also why
    it is process-wide: two runs going at once share one device declaration.  The run
    row records the device each was submitted under, so a stored result never loses it.
    """
    if device is not None:
        device = str(device).strip().lower()
        previous = os.environ.get("MOFSBU_DEVICE")
        os.environ["MOFSBU_DEVICE"] = device
        try:
            compute_device()                     # shape check, same as the CLI path
        except ValueError:
            if previous is None:
                os.environ.pop("MOFSBU_DEVICE", None)
            else:
                os.environ["MOFSBU_DEVICE"] = previous
            raise
    if workers is not None:
        n = max(1, int(workers))
        os.environ["MOFSBU_WORKERS"] = str(n)
        # The profile is what `parallel_enabled()` and the CLI report, so keep the two
        # from disagreeing: asking for more than one worker IS declaring a workstation.
        os.environ["MOFSBU_PROFILE"] = "workstation" if n > 1 else "laptop"
    return compute_state()
