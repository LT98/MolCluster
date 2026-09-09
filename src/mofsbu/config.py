"""Where the regenerable data lives.

The registry and blob store default to `<repo>/data/`, which is git-ignored.  If the
repo sits inside a synced folder (OneDrive, Dropbox), point `MOFSBU_DATA` somewhere
local instead: a sync client copying a live SQLite file and its -wal/-shm sidecars is a
real corruption risk, and none of this data needs syncing — it is regenerable by
construction, which is why it is out of git in the first place.

    export MOFSBU_DATA=~/.local/share/mofsbu      # linux / wsl
    setx MOFSBU_DATA %LOCALAPPDATA%\\mofsbu        # windows
"""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def data_root() -> Path:
    root = Path(os.environ.get("MOFSBU_DATA") or REPO_ROOT / "data")
    root.mkdir(parents=True, exist_ok=True)
    return root


def registry_path() -> Path:
    return data_root() / "registry.db"


def store_root() -> Path:
    return data_root() / "store"


# ── machine profile: the laptop must never be swamped by a build ──────────────
# Ground rule 8.  Parallelism is opt-in.  An unconfigured machine is treated as the
# laptop and runs ONE in-process worker; the workstation is declared, not detected, so a
# build cannot decide on its own to take every core on the machine you are typing on.

LAPTOP = "laptop"
WORKSTATION = "workstation"


def machine_profile() -> str:
    """`MOFSBU_PROFILE` = laptop | workstation.  Defaults to laptop."""
    value = os.environ.get("MOFSBU_PROFILE", LAPTOP).strip().lower()
    return WORKSTATION if value == WORKSTATION else LAPTOP


def max_workers() -> int:
    """How many workers may run.  1 unless explicitly raised.

    `MOFSBU_WORKERS` overrides everything (including on the laptop, if you mean it).
    Otherwise a declared workstation uses cpu_count - 1, leaving a core for the machine
    to stay responsive; anything else uses 1.
    """
    explicit = os.environ.get("MOFSBU_WORKERS")
    if explicit:
        try:
            return max(1, int(explicit))
        except ValueError:
            pass
    if machine_profile() == WORKSTATION:
        return max(1, (os.cpu_count() or 2) - 1)
    return 1


def parallel_enabled() -> bool:
    """True only when more than one worker is permitted.

    When False the runner executes in-process: no subprocesses, no process pool, nothing
    that can leave orphans behind if a laptop is closed mid-build.
    """
    return max_workers() > 1


# ── compute device: declared, never detected ─────────────────────────────────
# Ground rule 9's stance applied to the accelerator.  `MACEBackend` used to default to
# `device="cpu"` with no way to say otherwise, so a workstation with a GPU ran an MLIP on
# its CPU and the only symptom was that the card never warmed up.  Detection was
# considered and rejected for the same reason parallelism is opt-in: a build should not
# decide on its own to seize hardware you are using for something else.

CPU = "cpu"


def compute_device() -> str:
    """`MOFSBU_DEVICE` = cpu | cuda | cuda:N | mps.  Defaults to cpu.

    Returned verbatim after a shape check so a typo fails at the backend with the string
    you typed, rather than silently becoming "cpu" and looking like slow hardware.
    """
    value = (os.environ.get("MOFSBU_DEVICE") or CPU).strip().lower()
    if value != CPU and not value.startswith(("cuda", "mps", "xpu")):
        raise ValueError(
            f"MOFSBU_DEVICE={value!r} is not a device torch would recognise; "
            f"expected cpu, cuda, cuda:<n>, or mps")
    return value


def device_note() -> str:
    """One line for a run's first line of output, so an idle GPU is visible immediately."""
    device = compute_device()
    if device == CPU:
        return "device=cpu  (set MOFSBU_DEVICE=cuda to use a GPU)"
    return f"device={device}"
