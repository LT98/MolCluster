"""Shared type aliases, the fidelity ladder, and the exception hierarchy.

Named `_types` and not `types` on purpose: a module called `types.py` shadows the
standard library's `types` whenever its own directory lands first on `sys.path` —
which happens the moment anyone runs Python from inside this folder.  The failure
is spectacular and misleading (`functools` fails to import `GenericAlias` during
interpreter startup), so the name is simply avoided.
"""
from __future__ import annotations

from enum import IntEnum
from typing import TypeAlias

AtomMap: TypeAlias = dict[int, int]   # parent atom index -> child atom index


class Fidelity(IntEnum):
    """Ordered ladder.  A property of a GEOMETRY, never of a structure (D4)."""

    RAW = 0        # as constructed
    FF = 1         # force-field relaxed
    ML = 2         # MACE / other MLIP
    XTB = 3        # GFN2-xTB
    DFT = 4


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
