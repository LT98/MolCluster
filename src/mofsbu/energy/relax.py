"""Geometry optimisation.  NOT IMPLEMENTED (M7).

Ground rule 7: the signature is settled so the runner, the spec and the UI can all refer
to it today.  The body is scheduled work — it needs the energy backends (`tblite` for
xTB, MACE for the ML potential, an external code for DFT) behind one protocol, and the
reaction-balanced reference scheme, before any number it returns should be believed.
"""
from __future__ import annotations

from typing import Any

from mofsbu.assembly.join import NotBuiltYet
from mofsbu._types import Fidelity


def relax_geometry(coords: Any, symbols: list[str], *, charge: int, multiplicity: int,
                   target: Fidelity, solvent: str | None = None,
                   max_steps: int = 500) -> Any:
    """Optimise a geometry to the requested fidelity, returning coords + energy."""
    raise NotBuiltYet(
        f"energy.relax.relax_geometry(target={target!r}) - M7. "
        "Needs the backend protocol (xTB / MACE / DFT) and the reference scheme."
    )


def available_modes() -> dict[str, bool]:
    """Which run modes can actually execute.  The UI reads this rather than guessing."""
    return {"construct": True, "ml_go": False, "dft_go": False}
