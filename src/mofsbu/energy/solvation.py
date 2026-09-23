"""Continuum corrections on stored energies — the producer (WORKPLAN_solvation S2, C17).

MACE-OMOL-0 has no continuum, so an ML energy acquires a medium the way C17 describes: one
correction method, evaluated with and without the continuum on the SAME stored coordinates,
and the difference `dG_solv` stored beside the energy (`registry.put_solvation_correction`).
`energy.reference` then prices an equation in that medium from `energy + dG_solv` per term,
and refuses an equation whose terms do not all carry one.

Only ALPB through tblite is wired: it is what this environment runs (ddX is absent from the
installed build), and a medium names its model (`'alpb:water'`), never a solvent alone.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any

from mofsbu._types import Fidelity, MofsbuError, split_medium

__all__ = ["correct", "correct_registry", "MEDIA"]

#: Media this producer can compute: token -> the xTB backend's solvent key.
MEDIA = {"alpb:water": "water"}


def _backend_solvent(medium: str) -> str:
    model, solvent = split_medium(medium)
    if medium not in MEDIA:
        raise MofsbuError(
            f"medium {medium!r} ({model} in {solvent}) has no producer here; this build can "
            f"compute {sorted(MEDIA)}")
    return MEDIA[medium]


def correct(reg: Any, geometry_id: int, medium: str = "alpb:water") -> Any:
    """Compute and store the continuum correction for one geometry.  Returns the `Put`."""
    from mofsbu.energy.backends import get_backend
    from mofsbu.registry.api import put_solvation_correction

    solvent = _backend_solvent(medium)
    row = reg.conn.execute(
        "SELECT g.coords_hash, s.net_charge, s.multiplicity FROM geometries g "
        "JOIN structures s ON s.id = g.structure_id WHERE g.id=?", (geometry_id,)).fetchone()
    if row is None:
        raise MofsbuError(f"no geometry {geometry_id}")
    lines = [ln.split() for ln in reg.store.get(row["coords_hash"]).decode().splitlines()[2:]
             if ln.strip()]
    symbols = [ln[0] for ln in lines]
    positions = [[float(x) for x in ln[1:4]] for ln in lines]
    charge, multiplicity = int(row["net_charge"]), int(row["multiplicity"])

    xtb = get_backend("xtb")
    gas = xtb.single_point(symbols, positions, charge=charge, multiplicity=multiplicity)
    solv = xtb.single_point(symbols, positions, charge=charge, multiplicity=multiplicity,
                            solvent=solvent)
    return put_solvation_correction(reg, geometry_id, method=replace(gas.method, solvent=medium),
                                    e_gas=gas.energy, e_solv=solv.energy)


def correct_registry(reg: Any, medium: str = "alpb:water", *,
                     min_fidelity: Fidelity = Fidelity.ML, commit_every: int = 50,
                     limit: int | None = None) -> dict[str, Any]:
    """Correct every converged gas-phase energy at or above `min_fidelity` that lacks one.

    Resumable: a geometry that already has a correction in `medium` is skipped.  Every
    geometry rather than one per structure, because `energy.reference` takes the correction
    from the geometry whose energy it uses and never borrows one from a sibling.  A refusal
    is counted with its reason, not raised: one odd structure must not stop the rest.
    """
    _backend_solvent(medium)
    todo = [int(r["id"]) for r in reg.conn.execute(
        "SELECT g.id FROM geometries g JOIN methods m ON m.id = g.method_id "
        "WHERE g.energy IS NOT NULL AND g.converged = 1 AND m.solvent IS NULL "
        "AND g.fidelity >= ? AND NOT EXISTS ("
        "  SELECT 1 FROM solvation_corrections sc JOIN methods cm ON cm.id = sc.method_id "
        "  WHERE sc.geometry_id = g.id AND cm.solvent = ?) ORDER BY g.id",
        (int(min_fidelity), medium))]
    if limit is not None:
        todo = todo[:limit]
    written, refused = 0, {}
    for i, gid in enumerate(todo, 1):
        try:
            if correct(reg, gid, medium).created:
                written += 1
        except (MofsbuError, ValueError, RuntimeError) as exc:
            why = f"{type(exc).__name__}: {str(exc).splitlines()[0][:120]}"
            refused[why] = refused.get(why, 0) + 1
        if i % commit_every == 0:
            reg.conn.commit()
    reg.conn.commit()
    return {"medium": medium, "candidates": len(todo), "written": written, "refused": refused}
