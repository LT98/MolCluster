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
# Whether the RUNNER actually performs a relaxation.  It now does: a completed build task
# emits a `relax` follow-up task, and `runner._execute_relax` stores the result as a
# second geometry row chained by `relaxed_from`.
#
# The flag stays because the question it answers is a real one and was answered wrongly
# for a whole milestone: rev 16 shipped the planner's pre-flight WITHOUT the executor, so
# a 500-structure `ml_go` run passed every check, left the GPU idle, and wrote 281 raw
# constructs under a label saying they had been optimised.  `mode_status()` reports
# `backend_available` and `wired` separately for the same reason — "MACE imports" and
# "MACE will run" are different claims.
RELAXATION_IS_EXECUTED = True

MODE_FIDELITY: dict[str, Fidelity] = {
    "construct": Fidelity.RAW,
    "ml_go": Fidelity.ML,
    "xtb_go": Fidelity.XTB,
    "dft_go": Fidelity.DFT,
}


def relax_geometry(coords: Any, symbols: Sequence[str], *, charge: int, multiplicity: int,
                   target: Fidelity, solvent: str | None = None,
                   max_steps: int = 500, fmax: float = 0.20,
                   ml_model: str | None = None) -> RelaxResult:
    """Optimise a geometry to the requested fidelity, returning coords + energy.

    Raises `NotBuiltYet` if no backend serves that rung (DFT), and
    `EnergyBackendUnavailable` if one does but is not installed here — the laptop and
    the workstation give different answers to the second question and the same answer
    to the first, which is why they are different exceptions.
    """
    backend = backend_for(target, ml_model=ml_model)
    if not backend.available():
        raise EnergyBackendUnavailable(
            f"{target.name} relaxation needs the {backend.name} backend: "
            f"{backend.install_hint()}")
    return backend.relax(list(symbols), coords, charge=charge, multiplicity=multiplicity,
                         solvent=solvent, fmax=fmax, steps=max_steps)


def single_point(coords: Any, symbols: Sequence[str], *, charge: int, multiplicity: int,
                 target: Fidelity, solvent: str | None = None,
                 ml_model: str | None = None):
    """One energy, no relaxation.  The other half of what the reference scheme needs."""
    backend = backend_for(target, ml_model=ml_model)
    if not backend.available():
        raise EnergyBackendUnavailable(
            f"{target.name} single point needs the {backend.name} backend: "
            f"{backend.install_hint()}")
    return backend.single_point(list(symbols), coords, charge=charge,
                                multiplicity=multiplicity, solvent=solvent)


def available_modes(ml_model: str | None = None) -> dict[str, bool]:
    """Which run modes can actually execute, asked of the backends, not asserted.

    `construct` is always true — it is the placer, not a backend.  Everything else is
    true only if the stack behind it imports on this machine, so the same build shows
    `xtb_go` enabled on the workstation and disabled on a laptop without xTB, with no
    code change and no stale string in the markup.
    """
    return {mode: entry["available"] for mode, entry in mode_status(ml_model).items()}


def mode_status(ml_model: str | None = None) -> dict[str, dict[str, Any]]:
    """`available_modes` plus the reason, for a UI that has to explain a disabled option.

    `ml_model` is the model a SPEC asks for.  It matters here and not only at run time:
    a machine with mace-torch 0.3.6 can run `ml_go` with MACE-MP-0 and cannot run it
    with MACE-OMOL-0, and a planner that answered "ml_go is available" without being
    told which model would queue 500 tasks that all fail at the first import.
    """
    out: dict[str, dict[str, Any]] = {
        "construct": {"available": True, "backend": "placer", "backend_available": True,
                      "wired": True, "note": "", "method": "raw-construct",
                      "charge_aware": True}}
    for mode, fidelity in MODE_FIDELITY.items():
        if mode == "construct":
            continue
        try:
            backend = backend_for(fidelity, ml_model=ml_model)
        except Exception as exc:                                        # noqa: BLE001
            out[mode] = {"available": False, "backend": None, "backend_available": False,
                         "wired": RELAXATION_IS_EXECUTED, "note": str(exc),
                         "method": None, "charge_aware": None}
            continue
        # Two independent questions, kept apart on purpose: is the stack installed on
        # this machine, and does the pipeline call it?  Collapsing them is how "MACE is
        # importable" got reported to the user as "MACE will run".
        installed = backend.available()
        note = ""
        if not RELAXATION_IS_EXECUTED:
            note = ("the runner does not execute relaxations yet: `runner.execute` stores "
                    "the raw construct and never calls the backend, so this mode would "
                    "silently return unoptimised geometries")
        elif not installed:
            note = f"not installed here — {backend.install_hint()}"
        out[mode] = {"available": installed and RELAXATION_IS_EXECUTED,
                     "backend": backend.name, "backend_available": installed,
                     "wired": RELAXATION_IS_EXECUTED, "note": note,
                     # Which THEORY, not just which rung.  Two ML models sit on the same
                     # rung and produce numbers that must never be subtracted from one
                     # another, so the name travels with the capability read-out.
                     "method": getattr(backend, "method", None),
                     "charge_aware": bool(getattr(backend, "charge_aware", True))}
    return out


def ml_model_status() -> dict[str, dict[str, Any]]:
    """Both ML foundation models, side by side, for the UI's model picker.

    Not folded into `mode_status()`: that answers "can this run mode execute", and this
    answers "which theory would it be", which is the question the archived MP-0 numbers
    could not answer about themselves.
    """
    from mofsbu.config import ml_backend
    from mofsbu.energy.backends import ML_BACKENDS, get_backend

    declared = ml_backend()
    out: dict[str, dict[str, Any]] = {}
    for key in ML_BACKENDS:
        backend = get_backend(key)
        out[key] = {"method": backend.method, "available": backend.available(),
                    "declared": key == declared,
                    "charge_aware": backend.charge_aware,
                    "spin_aware": backend.spin_aware,
                    "training_set": getattr(backend, "training_set", ""),
                    "note": "" if backend.available() else backend.install_hint()}
    return out
