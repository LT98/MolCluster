"""Geometry optimisation — the caller-facing half of the energy layer (M7).

The signature was settled in M3 and is unchanged: the runner, the spec and the UI were
all written against it while the body raised.  What arrives now is the body, plus one
thing the stub could not have: `available_modes` reports what this MACHINE can run,
asked of the backends themselves, so the UI's disabled options track the environment
instead of a hard-coded dict that was true when it was typed.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from mofsbu._types import EnergyBackendUnavailable, Fidelity
from mofsbu.energy.backends import RelaxResult, backend_for

# Which fidelity each run mode is asking for.  One place, because the spec, the planner
# and the capability read-out all need the same answer.
MODE_FIDELITY: dict[str, Fidelity] = {
    "construct": Fidelity.RAW,
    "ml_go": Fidelity.ML,
    "xtb_go": Fidelity.XTB,
    "dft_go": Fidelity.DFT,
}


def relax_geometry(coords: Any, symbols: Sequence[str], *, charge: int, multiplicity: int,
                   target: Fidelity, solvent: str | None = None,
                   max_steps: int = 500, fmax: float = 0.20) -> RelaxResult:
    """Optimise a geometry to the requested fidelity, returning coords + energy.

    Raises `NotBuiltYet` if no backend serves that rung (DFT), and
    `EnergyBackendUnavailable` if one does but is not installed here — the laptop and
    the workstation give different answers to the second question and the same answer
    to the first, which is why they are different exceptions.
    """
    backend = backend_for(target)
    if not backend.available():
        raise EnergyBackendUnavailable(
            f"{target.name} relaxation needs the {backend.name} backend: "
            f"{backend.install_hint()}")
    return backend.relax(list(symbols), coords, charge=charge, multiplicity=multiplicity,
                         solvent=solvent, fmax=fmax, steps=max_steps)


def single_point(coords: Any, symbols: Sequence[str], *, charge: int, multiplicity: int,
                 target: Fidelity, solvent: str | None = None):
    """One energy, no relaxation.  The other half of what the reference scheme needs."""
    backend = backend_for(target)
    if not backend.available():
        raise EnergyBackendUnavailable(
            f"{target.name} single point needs the {backend.name} backend: "
            f"{backend.install_hint()}")
    return backend.single_point(list(symbols), coords, charge=charge,
                                multiplicity=multiplicity, solvent=solvent)


def available_modes() -> dict[str, bool]:
    """Which run modes can actually execute, asked of the backends, not asserted.

    `construct` is always true — it is the placer, not a backend.  Everything else is
    true only if the stack behind it imports on this machine, so the same build shows
    `xtb_go` enabled on the workstation and disabled on a laptop without xTB, with no
    code change and no stale string in the markup.
    """
    modes = {"construct": True}
    for mode, fidelity in MODE_FIDELITY.items():
        if mode == "construct":
            continue
        try:
            modes[mode] = backend_for(fidelity).available()
        except Exception:                       # no backend serves that rung at all
            modes[mode] = False
    return modes


def mode_status() -> dict[str, dict[str, Any]]:
    """`available_modes` plus the reason, for a UI that has to explain a disabled option."""
    out: dict[str, dict[str, Any]] = {
        "construct": {"available": True, "backend": "placer", "note": ""}}
    for mode, fidelity in MODE_FIDELITY.items():
        if mode == "construct":
            continue
        try:
            backend = backend_for(fidelity)
        except Exception as exc:                                        # noqa: BLE001
            out[mode] = {"available": False, "backend": None, "note": str(exc)}
            continue
        ok = backend.available()
        out[mode] = {"available": ok, "backend": backend.name,
                     "note": "" if ok else f"not installed here — {backend.install_hint()}"}
    return out
