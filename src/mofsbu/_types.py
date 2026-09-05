"""Shared type aliases, the fidelity ladder, and the exception hierarchy.

Named `_types` and not `types` on purpose: a module called `types.py` shadows the
standard library's `types` whenever its own directory lands first on `sys.path` —
which happens the moment anyone runs Python from inside this folder.  The failure
is spectacular and misleading (`functools` fails to import `GenericAlias` during
interpreter startup), so the name is simply avoided.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Any, TypeAlias

AtomMap: TypeAlias = dict[int, int]   # parent atom index -> child atom index


class Fidelity(IntEnum):
    """Ordered ladder.  A property of a GEOMETRY, never of a structure (D4)."""

    RAW = 0        # as constructed
    FF = 1         # force-field relaxed
    ML = 2         # MACE / other MLIP
    XTB = 3        # GFN2-xTB
    DFT = 4


@dataclass(frozen=True)
class MethodSpec:
    """What produced a number.  Ground rule 3: no bare floats anywhere.

    It lives here, next to `Fidelity`, because the two are the same rule seen from two
    sides — the ladder says how good a number is, the spec says what made it — and both
    the registry (which stores them) and `energy` (which produces them) need the type
    without depending on each other.
    """

    code: str                       # tblite | mace | heuristic | legacy | construct | null
    code_version: str
    method: str                     # GFN2-xTB | MACE-MP-0 | raw-construct
    solvent: str | None = None
    charge: int | None = None
    multiplicity: int | None = None
    extras: dict[str, Any] = field(default_factory=dict)

    def same_theory(self, other: MethodSpec) -> bool:
        """Are two numbers comparable — same code, version, method and medium?

        Charge and multiplicity are deliberately NOT compared: two species in one
        balanced reaction differ in both, and that is exactly the comparison the
        reference scheme has to make.  Everything else has to match or the difference
        of two energies is a difference of two theories.
        """
        return (self.code, self.code_version, self.method, self.solvent, self.extras) == (
            other.code, other.code_version, other.method, other.solvent, other.extras)

    def describe(self) -> str:
        medium = self.solvent or "gas"
        return f"{self.method}/{self.code}-{self.code_version} ({medium})"


class MofsbuError(Exception):
    """Base for every error this package raises deliberately."""


class AmbiguousSpecError(MofsbuError):
    """Something needed a value that was never specified.

    Ground rule 5: ambiguity branches explicitly. Never default silently.
    """


class GraphValidationError(MofsbuError):
    """A typed-graph invariant was violated."""


class CanonicalisationError(MofsbuError):
    """Canonical labelling failed or exceeded its search budget."""


class EnergyBackendUnavailable(MofsbuError):
    """A backend was asked for and its stack is not installed on this machine.

    Distinct from `NotBuiltYet`: the code exists, the environment does not.  The two
    are not interchangeable — one is scheduled work, the other is `conda install`.
    """


class ReferenceSchemeError(MofsbuError):
    """A reaction energy was requested for an equation that cannot carry one.

    Unbalanced atoms, unbalanced charge, species computed at different levels of
    theory, or a charge-blind backend asked about a charged reaction.  Ground rule 8:
    this raises rather than returning the number anyway, because the number would look
    exactly like a good one.
    """
