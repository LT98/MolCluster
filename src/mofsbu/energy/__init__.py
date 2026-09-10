"""energy: backends behind one protocol, and a reference scheme that refuses bad subtractions.

The two halves are deliberately separate.  `backends` produces numbers and records what
produced them; `reference` decides which numbers may be subtracted from which.  The
project's first-listed known risk (`DESIGN §11`) is entirely about the second half — the
gas-phase energies were never the problem, the arithmetic done with them was.
"""
from __future__ import annotations

from mofsbu._types import EnergyBackendUnavailable, Fidelity, MethodSpec, ReferenceSchemeError
from mofsbu.energy.backends import (
    ML_BACKENDS, EnergyBackend, EnergyResult, MACEBackend, MACEOmolBackend, NullBackend,
    RelaxResult, XTBBackend, available_backends, backend_for, combined_multiplicity,
    d_electrons, get_backend, high_spin_multiplicity, minimal_multiplicity, ml_backend_key,
    check_spin, spin_class_multiplicity,
)
from mofsbu.energy.reference import (
    BalanceReport, QualityIssue, ReactionEnergy, ReferenceQuality, Term, check_balance,
    check_reference_quality, put_balanced_reaction, reaction_balanced_energy,
    reaction_terms, store_reaction_energy,
)
from mofsbu.energy.relax import (
    MODE_FIDELITY, available_modes, ml_model_status, mode_status, relax_geometry,
    single_point,
)

__all__ = [
    "BalanceReport", "EnergyBackend", "EnergyBackendUnavailable", "EnergyResult",
    "Fidelity", "MACEBackend", "MACEOmolBackend", "ML_BACKENDS", "MODE_FIDELITY",
    "MethodSpec", "NullBackend",
    "QualityIssue", "ReactionEnergy", "ReferenceQuality", "ReferenceSchemeError",
    "RelaxResult", "Term", "XTBBackend", "available_backends", "available_modes",
    "backend_for", "check_balance", "check_reference_quality", "check_spin",
    "combined_multiplicity", "d_electrons", "get_backend", "high_spin_multiplicity",
    "minimal_multiplicity", "ml_backend_key", "ml_model_status", "mode_status",
    "put_balanced_reaction", "reaction_balanced_energy", "reaction_terms",
    "relax_geometry", "single_point", "spin_class_multiplicity", "store_reaction_energy",
]
