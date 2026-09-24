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


def name_process(name: str) -> None:
    """Name this process as a system monitor shows it (`top`, `ps -o comm`, GNOME).

    Otherwise every viewer and worker is just `python`, which is what an orphaned viewer
    looks like after an SSH session drops.  Linux only (`/proc/self/comm`, 15 characters);
    anywhere else, or if the kernel refuses, nothing happens — a name is a convenience.
    """
    try:
        with open("/proc/self/comm", "w") as f:
            f.write(name[:15])
    except OSError:
        pass


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


def relax_workers() -> int:
    """How many of the declared workers drain the RELAX queue, and only that queue.

    0 means no split: every worker claims any task, which is the right shape when
    relaxation is CPU work like construction and the two compete for the same cores.

    On a declared accelerator it is 1.  The card is one device — a second process feeding
    it divides its memory rather than multiplying its throughput — while construction is
    CPU work that the same pool was serialising behind it.  Splitting the queue is what
    lets the two run at once: the GPU relaxes what has been built while the remaining
    workers build the next batch.

    `MOFSBU_RELAX_WORKERS` overrides, including back to 0 (share one pool) or up to N for
    a multi-GPU box.  Declared, never detected — the device it keys off is itself declared.
    """
    explicit = os.environ.get("MOFSBU_RELAX_WORKERS")
    if explicit:
        try:
            return max(0, int(explicit))
        except ValueError:
            pass
    try:
        return 0 if compute_device() == CPU else 1
    except ValueError:
        # A malformed MOFSBU_DEVICE is reported where it is set, not here; until then
        # this run has no accelerator to feed.
        return 0


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

# ── which ML potential serves the ML rung: declared, never guessed ───────────
# There is more than one MACE foundation model and they are NOT interchangeable.
# MACE-MP-0 is trained on Materials Project relaxations and is blind to formal charge
# and spin; MACE-OMOL-0 is trained on OMol25 (wB97M-V/def2-TZVPD, eV) and takes total
# charge and spin multiplicity as INPUTS, which is what makes it usable for the charged
# species the reference scheme is built out of.  Their energies share a rung on the
# ladder and nothing else: they are different theories and must never be subtracted from
# one another.  Which one runs is therefore declared, like the device and the worker
# count, and the choice is recorded in the `methods` row of every number it produces.

MACE_MP = "mace"
MACE_OMOL = "mace_omol"
ML_BACKENDS = (MACE_MP, MACE_OMOL)

_ML_ALIASES = {
    "mace": MACE_MP, "mp": MACE_MP, "mace-mp": MACE_MP, "mace_mp": MACE_MP,
    "mace-mp-0": MACE_MP, "mace_mp_0": MACE_MP, "macemp0": MACE_MP,
    "mace_omol": MACE_OMOL, "omol": MACE_OMOL, "mace-omol": MACE_OMOL,
    "mace-omol-0": MACE_OMOL, "mace_omol_0": MACE_OMOL, "maceomol0": MACE_OMOL,
}


def resolve_ml_backend(name: str | None) -> str:
    """Map a user-facing model name onto a backend key, or refuse with the list.

    Accepts the spellings people actually type ("MACE-OMOL-0", "omol", "mace_omol").
    Refuses anything else rather than falling back to MP-0: a typo that silently
    selected the charge-blind model would produce numbers that look fine and mean
    something different.
    """
    if name is None:
        return ml_backend()
    key = _ML_ALIASES.get(str(name).strip().lower().replace(" ", ""))
    if key is None:
        raise ValueError(
            f"unknown ML model {name!r}; expected one of "
            f"{sorted({'mace-mp-0', 'mace-omol-0'})} "
            f"(backend keys {list(ML_BACKENDS)})")
    return key


def ml_backend() -> str:
    """`MOFSBU_ML_MODEL` = mace-mp-0 | mace-omol-0.  Defaults to MACE-MP-0.

    MP-0 is the default only because it is what every stored ML number in this project
    was produced with; it is not the better model for charged coordination complexes.
    """
    value = os.environ.get("MOFSBU_ML_MODEL")
    if not value:
        return MACE_MP
    return resolve_ml_backend(value)


def ml_note() -> str:
    """One line for a run's header, next to `device_note()`."""
    key = ml_backend()
    label = "MACE-MP-0 (charge- and spin-blind)" if key == MACE_MP else \
            "MACE-OMOL-0 (charge- and spin-aware)"
    return f"ml_model={label}  (set MOFSBU_ML_MODEL to change)"
