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
import re
import sys
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

    With a pool of cards (`MOFSBU_DEVICE=cuda:0,cuda:1`) it is one per card; the runner
    pins relax worker i to card i mod n.

    `MOFSBU_RELAX_WORKERS` overrides, including back to 0 (share one pool).  Declared, never
    detected — the devices it keys off are themselves declared.
    """
    explicit = os.environ.get("MOFSBU_RELAX_WORKERS")
    if explicit:
        try:
            return max(0, int(explicit))
        except ValueError:
            pass
    try:
        devices = compute_devices()
        return 0 if devices == [CPU] else len(devices)
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
_DEVICE = re.compile(r"cpu|mps|(?:cuda|xpu)(?::\d+)?")


_CARD = re.compile(r"cuda:\d+")


def compute_devices() -> list[str]:
    """`MOFSBU_DEVICE` = cpu | cuda | cuda:N | mps | xpu, or a pool of cards `cuda:0,cuda:1`.

    Defaults to cpu.  Returned verbatim after a shape check so a typo fails here with the
    string you typed, rather than silently becoming "cpu" and looking like slow hardware.
    A pool lists distinct indexed CUDA cards only: one process uses one card, and the
    runner gives each relax worker its own (`runner.worker_env`).
    """
    value = (os.environ.get("MOFSBU_DEVICE") or CPU).strip().lower()
    parts = [p.strip() for p in value.split(",")]
    if len(parts) == 1 and _DEVICE.fullmatch(parts[0]):
        return parts
    if len(parts) > 1 and all(_CARD.fullmatch(p) for p in parts) and len(set(parts)) == len(parts):
        return parts
    raise ValueError(
        f"MOFSBU_DEVICE={value!r} is not a device torch would recognise; expected cpu, "
        f"cuda, cuda:<n>, mps or xpu, or a list of distinct cards such as 'cuda:0,cuda:1'")


def compute_device() -> str:
    """The one device this process computes on: the first of a pool."""
    return compute_devices()[0]


def physical_device() -> str:
    """The device this process runs on, as the machine numbers it, for records.

    A pinned worker sees only its card, as `cuda:0`; `MOFSBU_DEVICE_PHYSICAL` says which.
    """
    return os.environ.get("MOFSBU_DEVICE_PHYSICAL") or compute_device()


def pin_env(device: str) -> dict[str, str]:
    """The environment that confines a process to one indexed card, or {} for anything else.

    The card is the only one the process sees, as `cuda:0`; `MOFSBU_DEVICE_PHYSICAL` keeps
    its real name.  GPU4PySCF works on every card it can see (#70), so this is how DFT is
    kept on one.
    """
    card = re.fullmatch(r"cuda:(\d+)", device)
    if card is None:
        return {}
    return {"CUDA_VISIBLE_DEVICES": card.group(1), "MOFSBU_DEVICE": "cuda:0",
            "MOFSBU_DEVICE_PHYSICAL": device}


def pin_process(index: int = 0) -> str:
    """Confine this process to card `index` mod n of the declared pool; return its name.

    Call before anything touches CUDA.  A process already confined (`CUDA_VISIBLE_DEVICES`
    set) is left as it is.  Refused once torch has initialised CUDA, because by then the
    other cards are already open.
    """
    if os.environ.get("CUDA_VISIBLE_DEVICES") is not None:
        return physical_device()
    devices = compute_devices()
    env = pin_env(devices[index % len(devices)])
    if not env:
        return physical_device()
    torch = sys.modules.get("torch")
    if torch is not None and torch.cuda.is_initialized():
        raise RuntimeError("pin_process() after CUDA was initialised: call it first")
    os.environ.update(env)
    return physical_device()


def visible_cards() -> int | None:
    """How many CUDA cards this process can see, without initialising CUDA; None if unknown."""
    value = os.environ.get("CUDA_VISIBLE_DEVICES")
    if value is not None:
        return len([v for v in value.split(",") if v.strip()])
    try:
        return len(os.listdir("/proc/driver/nvidia/gpus"))
    except OSError:
        return None


def device_note() -> str:
    """One line for a run's first line of output, so an idle GPU is visible immediately."""
    devices = compute_devices()
    if devices == [CPU]:
        return "device=cpu  (set MOFSBU_DEVICE=cuda to use a GPU)"
    if len(devices) > 1:
        return (f"devices={','.join(devices)}  (one card per relax worker; work in this "
                f"process uses {devices[0]})")
    return f"device={physical_device()}"


# ── headroom floors: declared, never detected ────────────────────────────────
# A worker that would start a task with less than this free leaves it pending instead
# (`headroom.shortfall`, asked before every claim).  Unset or 0: nothing is checked.


def _floor_gib(var: str) -> float:
    value = os.environ.get(var, "").strip()
    if not value:
        return 0.0
    try:
        floor = float(value)
    except ValueError:
        raise ValueError(f"{var}={value!r} is not a number of GiB") from None
    if not floor >= 0:
        raise ValueError(f"{var}={value!r} must be >= 0 GiB")
    return floor


def min_free_ram_gib() -> float:
    """`MOFSBU_MIN_FREE_RAM_GIB`: claim nothing while less system RAM than this is free."""
    return _floor_gib("MOFSBU_MIN_FREE_RAM_GIB")


def min_free_vram_gib() -> float:
    """`MOFSBU_MIN_FREE_VRAM_GIB`: a worker feeding `compute_device()` claims nothing while
    less device memory than this is free."""
    return _floor_gib("MOFSBU_MIN_FREE_VRAM_GIB")


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
MACE_MH = "mace_mh"
MACE_POLAR = {size: f"mace_polar_{size}" for size in "sml"}
ML_BACKENDS = (MACE_MP, MACE_OMOL, MACE_MH, *MACE_POLAR.values())

_ML_ALIASES = {
    "mace": MACE_MP, "mp": MACE_MP, "mace-mp": MACE_MP, "mace_mp": MACE_MP,
    "mace-mp-0": MACE_MP, "mace_mp_0": MACE_MP, "macemp0": MACE_MP,
    "mace_omol": MACE_OMOL, "omol": MACE_OMOL, "mace-omol": MACE_OMOL,
    "mace-omol-0": MACE_OMOL, "mace_omol_0": MACE_OMOL, "maceomol0": MACE_OMOL,
    "mace_mh": MACE_MH, "mh": MACE_MH, "mh-1": MACE_MH, "mace-mh": MACE_MH,
    "mace-mh-1": MACE_MH, "mace_mh_1": MACE_MH, "macemh1": MACE_MH,
    **{alias: key for size, key in MACE_POLAR.items()
       for alias in (key, f"mace-polar-1-{size}", f"mace_polar_1_{size}", f"polar-1-{size}")},
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
            f"{sorted({'mace-mp-0', 'mace-omol-0', 'mace-mh-1', 'mace-polar-1-s', 'mace-polar-1-m', 'mace-polar-1-l'})} "
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
    label = {MACE_MP: "MACE-MP-0 (charge- and spin-blind)",
             MACE_OMOL: "MACE-OMOL-0 (charge- and spin-aware)",
             MACE_MH: "MACE-MH-1 (charge- and spin-blind)",
             **{k: f"MACE-POLAR-1-{s.upper()} (charge- and spin-aware)"
                for s, k in MACE_POLAR.items()}}[key]
    return f"ml_model={label}  (set MOFSBU_ML_MODEL to change)"
