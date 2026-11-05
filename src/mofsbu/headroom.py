"""Is there room on this machine to start one more task?  Asked before every claim.

A task that is never claimed cannot fail for want of memory.  `runner.work` asks
`shortfall` before it claims, and while there is a reason it leaves the queue alone: the
task stays `pending`, for this worker later or for another worker now.  The floors are
declared (`config.min_free_ram_gib`, `config.min_free_vram_gib`); unset, nothing is read.

Both reads are cheap enough to make before every claim (~8 µs for /proc/meminfo, ~4 µs
for cudaMemGetInfo once the process has a CUDA context), so nothing is cached.  They are a
snapshot, not a reservation: two workers can both see room for one task.
"""
from __future__ import annotations

import sys

from mofsbu._types import MofsbuError
from mofsbu.config import compute_device, min_free_ram_gib, min_free_vram_gib

GIB = 1024 ** 3


def free_ram_gib() -> float | None:
    """`MemAvailable` from /proc/meminfo, or None where there is no such file."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024 / GIB
    except OSError:
        return None
    return None


def free_vram_gib(device: str) -> float | None:
    """Free memory on a CUDA `device`, plus what this process's own allocators hold and are
    not using (it is free to this worker's next task).  None for a non-CUDA device."""
    if not device.startswith("cuda"):
        return None
    import torch

    if not torch.cuda.is_available():
        return None
    index = torch.device(device).index
    index = torch.cuda.current_device() if index is None else index
    free, _total = torch.cuda.mem_get_info(index)
    free += torch.cuda.memory_reserved(index) - torch.cuda.memory_allocated(index)
    cupy = sys.modules.get("cupy")          # GPU4PySCF's pool, only if already in use
    if cupy is not None:
        free += cupy.get_default_memory_pool().free_bytes()
    return free / GIB


def shortfall(*, gpu: bool) -> str | None:
    """Why this worker should not claim now, or None.  `gpu`: it feeds `compute_device()`.

    A floor that is declared but cannot be read here is refused, not skipped: a check that
    silently stops checking is worse than none.
    """
    ram_floor = min_free_ram_gib()
    if ram_floor:
        ram = free_ram_gib()
        if ram is None:
            raise MofsbuError("MOFSBU_MIN_FREE_RAM_GIB is set but this machine has no "
                              "/proc/meminfo to read free memory from")
        if ram < ram_floor:
            return f"{ram:.1f} GiB RAM free < MOFSBU_MIN_FREE_RAM_GIB={ram_floor:g}"
    vram_floor = min_free_vram_gib() if gpu else 0.0
    if vram_floor:
        device = compute_device()
        vram = free_vram_gib(device)
        if vram is None:
            raise MofsbuError(f"MOFSBU_MIN_FREE_VRAM_GIB is set but MOFSBU_DEVICE={device!r} "
                              "is not a CUDA device this process can see")
        if vram < vram_floor:
            return (f"{vram:.1f} GiB free on {device} < "
                    f"MOFSBU_MIN_FREE_VRAM_GIB={vram_floor:g}")
    return None
