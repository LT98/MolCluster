"""Load the log folder, find what is still to record, and compute amounts.

Folder layout (under `data_root()/experiments` unless `--root` says otherwise)::

    reagents.yaml      code -> reagent
    phases.yaml        label -> crystalline phase
    entries/<id>.yaml  one experiment each
    raw/               instrument files the entries point at
    INDEX.md           written by `expt index`
    .schema/           written by every command; editor validation

Two kinds of finding: *pending* (blank, still to record) and *problems* (a value that is
wrong in kind, or points at something that does not exist).  Neither is ever filled in.
"""
from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .chem import molar_mass
from .spec import (
    OUTCOME_SECTIONS,
    PHASE,
    REAGENT,
    SENTINELS,
    Block,
    Items,
    Leaf,
    entry_spec,
)

RATIO_TOLERANCE = 0.05
PRODUCT_PREFIX = "product:"


@dataclass
class Findings:
    pending: list[str] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)


@dataclass
class Log:
    root: Path
    reagents: dict = field(default_factory=dict)
    phases: dict = field(default_factory=dict)
    entries: dict[str, dict] = field(default_factory=dict)
    files: dict[str, Path] = field(default_factory=dict)
    load_errors: dict[str, str] = field(default_factory=dict)

    def entry_path(self, entry_id: str) -> Path:
        return self.root / "entries" / f"{entry_id}.yaml"


def _read_yaml(path: Path, log: Log) -> object:
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        log.load_errors[path.name] = str(exc).replace("\n", " ")
        return None


def load_log(root: Path) -> Log:
    log = Log(root=root)
    for name in ("reagents", "phases"):
        path = root / f"{name}.yaml"
        if path.exists():
            data = _read_yaml(path, log) or {}
            if isinstance(data, dict):
                setattr(log, name, data)
            else:
                log.load_errors[path.name] = "expected a mapping of code -> fields"
    for path in sorted((root / "entries").glob("*.yaml")):
        data = _read_yaml(path, log)
        if data is None:
            log.load_errors.setdefault(path.name, "empty file")
            continue
        if not isinstance(data, dict):
            log.load_errors[path.name] = "expected a mapping of field -> value"
            continue
        log.entries[path.stem] = data
        log.files[path.stem] = path
    return log


# ---- values -----------------------------------------------------------------------------

def is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def recorded(value: object) -> bool:
    """True when anything under `value` was written down (a sentinel counts)."""
    if value is None:
        return False
    if isinstance(value, dict):
        return any(recorded(v) for v in value.values())
    if isinstance(value, list):
        return any(recorded(v) for v in value)
    return True


def _leaf_problem(leaf: Leaf, value: object, log: Log) -> str | None:
    kind = leaf.kind
    if kind == "num" and not is_number(value):
        return f"expected a number, got {value!r}"
    if kind == "int" and not (isinstance(value, int) and not isinstance(value, bool)):
        return f"expected a whole number, got {value!r}"
    if kind == "temp" and not (is_number(value) or value == "RT"):
        return f"expected °C or RT, got {value!r}"
    if kind in ("str", "path", "reagent", "phase", "entry") and (
            isinstance(value, bool) or not isinstance(value, (str, int, float))):
        return f"expected text, got {value!r} (put it in quotes)"
    if kind == "text" and not isinstance(value, str):
        return f"expected text, got {value!r}"
    if kind == "date" and not isinstance(value, dt.date):
        return f"expected a date YYYY-MM-DD, got {value!r}"
    if kind == "enum" and value not in leaf.choices:
        return f"{value!r} is not one of {', '.join(leaf.choices)}"
    if kind == "path" and not (log.root / str(value)).exists():
        return f"file {value} not found under {log.root.name}/"
    if kind == "reagent":
        text = str(value)
        if text.startswith(PRODUCT_PREFIX):
            if text[len(PRODUCT_PREFIX):] not in log.entries:
                return f"no entry {text[len(PRODUCT_PREFIX):]!r} to take a product from"
        elif text not in log.reagents:
            return f"reagent {text!r} is not in reagents.yaml"
    if kind == "phase" and str(value) not in log.phases:
        return f"phase {value!r} is not in phases.yaml"
    if kind == "entry" and str(value) not in log.entries:
        return f"no entry {value!r}"
    return None


def _item_name(node: Items, item: object, index: int) -> str:
    if isinstance(item, dict) and node.key and item.get(node.key) not in (None, ""):
        return str(item[node.key])
    return str(index + 1)


def walk(node, value: object, path: str, log: Log, out: Findings, skip: tuple = ()) -> None:
    if isinstance(value, str) and value in SENTINELS:
        return
    if not recorded(value):
        if node.required:
            out.pending.append(path if isinstance(node, Leaf) else f"{path} (nothing recorded)")
        return
    if isinstance(node, Leaf):
        problem = _leaf_problem(node, value, log)
        if problem:
            out.problems.append(f"{path}: {problem}")
        return
    if isinstance(node, Block):
        if not isinstance(value, dict):
            out.problems.append(f"{path}: expected fields, got {value!r}")
            return
        for key, child in node.fields.items():
            if any(child is s for s in skip):
                continue
            walk(child, value.get(key), f"{path}.{key}" if path else key, log, out, skip)
        for key in value:
            if key not in node.fields:
                out.problems.append(f"{path + '.' if path else ''}{key}: not a field the log knows")
        return
    if not isinstance(value, list):
        out.problems.append(f"{path}: expected a list, got {value!r}")
        return
    for i, item in enumerate(value):
        walk(node.item, item, f"{path}[{_item_name(node, item, i)}]", log, out, skip)


# ---- amounts ----------------------------------------------------------------------------

def reagent_mw(log: Log, code: object) -> float | None:
    info = log.reagents.get(code) if isinstance(code, str) else None
    if not isinstance(info, dict):
        return None
    if is_number(info.get("mw")):
        return float(info["mw"])
    formula = info.get("formula")
    if isinstance(formula, str) and formula not in SENTINELS:
        try:
            return molar_mass(formula)
        except ValueError:
            return None
    return None


def line_mmol(log: Log, line: dict) -> float | None:
    if is_number(line.get("mmol")):
        return float(line["mmol"])
    mw = reagent_mw(log, line.get("code"))
    if mw is None:
        return None
    if is_number(line.get("mass_mg")):
        return float(line["mass_mg"]) / mw
    info = log.reagents.get(line.get("code"), {})
    if is_number(line.get("volume_ml")) and is_number(info.get("density_g_ml")):
        return float(line["volume_ml"]) * float(info["density_g_ml"]) * 1000.0 / mw
    return None


def reagent_lines(entry: dict) -> list[dict]:
    lines = entry.get("reagents")
    return [r for r in lines if isinstance(r, dict)] if isinstance(lines, list) else []


def eta_ul_per_mg(entry: dict) -> float | None:
    """Liquid-to-solid ratio η for a milled entry, when every solid mass is known."""
    cond = entry.get("conditions")
    if entry.get("method") != "ball_mill" or not isinstance(cond, dict):
        return None
    solids = [r for r in reagent_lines(entry) if r.get("role") != "solvent"]
    if not is_number(cond.get("liquid_ul")) or not solids:
        return None
    if not all(is_number(r.get("mass_mg")) for r in solids):
        return None
    total = sum(float(r["mass_mg"]) for r in solids)
    return float(cond["liquid_ul"]) / total if total else None


def _check_amounts(log: Log, entry: dict, out: Findings) -> None:
    measured = []
    for i, line in enumerate(reagent_lines(entry)):
        name = f"reagents[{line.get('code') or i + 1}]"
        solvent = line.get("role") == "solvent"
        if solvent and f"{name}.ratio" in out.pending:
            out.pending.remove(f"{name}.ratio")
        if not any(recorded(line.get(k)) for k in ("mass_mg", "volume_ml", "mmol")):
            out.pending.append(f"{name}.amount (mass_mg, volume_ml or mmol)")
        mmol = line_mmol(log, line)
        if not solvent and mmol is not None and is_number(line.get("ratio")) and line["ratio"]:
            measured.append((name, mmol, float(line["ratio"])))
    if len(measured) < 2:
        return
    _, ref_mmol, ref_ratio = measured[0]
    for name, mmol, ratio in measured[1:]:
        stated = ratio / ref_ratio
        actual = mmol / ref_mmol
        if abs(actual / stated - 1) > RATIO_TOLERANCE:
            out.problems.append(
                f"{name}: ratio says {stated:.3g} eq of the first reagent, amounts give "
                f"{actual:.3g}")


# ---- checks -----------------------------------------------------------------------------

def check_entry(log: Log, entry_id: str) -> Findings:
    entry = log.entries[entry_id]
    out = Findings()
    spec = entry_spec(entry.get("method"))
    skip: tuple = ()
    if entry.get("status") in ("planned", "abandoned"):
        skip = tuple(spec.fields[k] for k in OUTCOME_SECTIONS)
    walk(spec, entry, "", log, out, skip)
    if entry.get("id") != entry_id:
        out.problems.append(f"id: {entry.get('id')!r} does not match the file name {entry_id}")
    _check_amounts(log, entry, out)
    return out


def check_reagents(log: Log) -> dict[str, Findings]:
    found = {}
    for code, info in log.reagents.items():
        out = Findings()
        walk(REAGENT, info, "", log, out)
        formula = info.get("formula") if isinstance(info, dict) else None
        if isinstance(formula, str) and formula not in SENTINELS:
            try:
                molar_mass(formula)
            except ValueError as exc:
                out.problems.append(f"formula: {exc}")
        found[str(code)] = out
    return found


def check_phases(log: Log) -> dict[str, Findings]:
    found = {}
    for label, info in log.phases.items():
        out = Findings()
        walk(PHASE, info, "", log, out)
        found[str(label)] = out
    return found
