"""Molar mass from a written formula, including hydrates: `Cu(C2H3O2)2·H2O`, `NiCl2·6H2O`."""
from __future__ import annotations

import re

from rdkit import Chem

_HYDRATE_SEP = re.compile(r"\s*[·•*]\s*")
_TOKEN = re.compile(r"([A-Z][a-z]?)|(\()|(\))|(\d+(?:\.\d+)?)")
_TABLE = Chem.GetPeriodicTable()
_WEIGHTS = {_TABLE.GetElementSymbol(z): _TABLE.GetAtomicWeight(z) for z in range(1, 119)}


def _part_mass(part: str) -> float:
    m = re.match(r"(\d+(?:\.\d+)?)?(.*)$", part)
    factor = float(m.group(1)) if m.group(1) else 1.0
    body = m.group(2)
    if not body:
        raise ValueError(f"no formula after the count in {part!r}")
    stack: list[float] = [0.0]
    last: float | None = None          # mass of the element or group a count applies to
    pos = 0
    while pos < len(body):
        tok = _TOKEN.match(body, pos)
        if tok is None:
            raise ValueError(f"cannot read {body[pos:]!r} in formula {part!r}")
        pos = tok.end()
        element, open_, close, count = tok.groups()
        if element:
            if element not in _WEIGHTS:
                raise ValueError(f"unknown element {element!r} in {part!r}")
            last = _WEIGHTS[element]
            stack[-1] += last
        elif open_:
            stack.append(0.0)
            last = None
        elif close:
            if len(stack) == 1:
                raise ValueError(f"unbalanced ')' in {part!r}")
            last = stack.pop()
            stack[-1] += last
        else:
            if last is None:
                raise ValueError(f"count with nothing to apply to in {part!r}")
            stack[-1] += last * (float(count) - 1)
            last = None
    if len(stack) != 1:
        raise ValueError(f"unbalanced '(' in {part!r}")
    return factor * stack[0]


def molar_mass(formula: str) -> float:
    """g/mol.  Hydrate parts are separated by `·`, `•` or `*`; a leading number multiplies."""
    parts = [p for p in _HYDRATE_SEP.split(formula.strip()) if p]
    if not parts:
        raise ValueError("empty formula")
    return sum(_part_mass(p) for p in parts)
