"""Store the log as JSON, find what is still to record, and compute amounts.

Folder layout (`data_root()/experiments` unless the app is given another)::

    reagents.json      code -> reagent
    phases.json        label -> crystalline phase
    entries/<id>.json  one experiment each; entries/.deleted/ keeps removed ones
    raw/               instrument files the entries point at
    INDEX.md           rewritten on every save, for reading without the page

Two kinds of finding: *pending* (blank, still to record) and *problems* (a value that is
wrong in kind, or points at something that does not exist).  Neither is ever filled in.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

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
ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
DATE_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}$")


class StaleVersion(Exception):
    """The file changed since the caller read it."""


@dataclass
class Findings:
    pending: list[str] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)


# ---- files ------------------------------------------------------------------------------

def version_of(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16] if path.exists() else None


def read_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def write_json(path: Path, data: object, expect: str | None) -> str:
    """Write atomically, refusing if the file is not the version the caller read.

    `expect=None` means the caller believes there is no file yet."""
    if version_of(path) != expect:
        raise StaleVersion(f"{path.name} changed since it was read")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    return version_of(path)


@dataclass
class Log:
    root: Path
    reagents: dict = field(default_factory=dict)
    phases: dict = field(default_factory=dict)
    entries: dict[str, dict] = field(default_factory=dict)
    load_errors: dict[str, str] = field(default_factory=dict)

    def entry_path(self, entry_id: str) -> Path:
        if not ID_PATTERN.match(entry_id):
            raise ValueError(f"{entry_id!r} is not a valid id: letters, digits, . _ - only")
        return self.root / "entries" / f"{entry_id}.json"

    @property
    def reagents_path(self) -> Path:
        return self.root / "reagents.json"

    @property
    def phases_path(self) -> Path:
        return self.root / "phases.json"


def load_log(root: Path) -> Log:
    log = Log(root=root)
    for name in ("reagents", "phases"):
        path = root / f"{name}.json"
        try:
            data = read_json(path) or {}
        except json.JSONDecodeError as exc:
            log.load_errors[path.name] = str(exc)
            continue
        if isinstance(data, dict):
            setattr(log, name, data)
        else:
            log.load_errors[path.name] = "expected an object of code -> fields"
    for path in sorted((root / "entries").glob("*.json")):
        try:
            data = read_json(path)
        except json.JSONDecodeError as exc:
            log.load_errors[path.name] = str(exc)
            continue
        if isinstance(data, dict):
            log.entries[path.stem] = data
        else:
            log.load_errors[path.name] = "expected an object of field -> value"
    return log


def remove_entry(log: Log, entry_id: str) -> Path:
    """Move an entry to entries/.deleted/ — kept, out of the way, never destroyed."""
    src = log.entry_path(entry_id)
    if not src.exists():
        raise FileNotFoundError(entry_id)
    stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    dst = log.root / "entries" / ".deleted" / f"{entry_id}.{stamp}.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.replace(src, dst)
    return dst


# ---- values -----------------------------------------------------------------------------

def is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def recorded(value: object) -> bool:
    """True when anything under `value` was written down (a sentinel counts)."""
    if value is None or value == "":
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
    if kind in ("str", "text", "path", "reagent", "phase", "entry", "date") and not isinstance(
            value, str):
        return f"expected text, got {value!r}"
    if kind == "date" and not (DATE_PATTERN.match(value) and _is_date(value)):
        return f"expected a date YYYY-MM-DD, got {value!r}"
    if kind == "enum" and value not in leaf.choices:
        return f"{value!r} is not one of {', '.join(leaf.choices)}"
    if kind == "path" and not (log.root / value).exists():
        return f"file {value} not found in the log folder"
    if kind == "reagent":
        if value.startswith(PRODUCT_PREFIX):
            if value[len(PRODUCT_PREFIX):] not in log.entries:
                return f"no entry {value[len(PRODUCT_PREFIX):]!r} to take a product from"
        elif value not in log.reagents:
            return f"reagent {value!r} is not in the reagent table"
    if kind == "phase" and value not in log.phases:
        return f"phase {value!r} is not in the phase table"
    if kind == "entry" and value not in log.entries:
        return f"no entry {value!r}"
    return None


def _is_date(text: str) -> bool:
    try:
        dt.date.fromisoformat(text)
    except ValueError:
        return False
    return True


def item_name(node: Items, item: object, index: int) -> str:
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
        walk(node.item, item, f"{path}[{item_name(node, item, i)}]", log, out, skip)


# ---- amounts ----------------------------------------------------------------------------

def reagent_mw(log: Log, code: object) -> float | None:
    info = log.reagents.get(code) if isinstance(code, str) else None
    if not isinstance(info, dict):
        return None
    if is_number(info.get("mw")):
        return float(info["mw"])
    formula = info.get("formula")
    if isinstance(formula, str) and formula and formula not in SENTINELS:
        try:
            return molar_mass(formula)
        except ValueError:
            return None
    return None


def component_mmol(log: Log, comp: dict) -> tuple[float | None, str]:
    """mmol and where it came from: entered, mass, volume (density, assay) or molarity."""
    if is_number(comp.get("mmol")):
        return float(comp["mmol"]), "entered"
    info = log.reagents.get(comp.get("code"), {})
    info = info if isinstance(info, dict) else {}
    volume = comp.get("volume_ml")
    if is_number(volume) and is_number(info.get("molarity_m")):
        return float(volume) * float(info["molarity_m"]), "volume × molarity"
    mw = reagent_mw(log, comp.get("code"))
    if mw is None:
        return None, ""
    if is_number(comp.get("mass_mg")):
        mass, how = float(comp["mass_mg"]), "mass"
    elif is_number(volume) and is_number(info.get("density_g_ml")):
        mass, how = float(volume) * float(info["density_g_ml"]) * 1000.0, "volume × density"
    else:
        return None, ""
    if is_number(info.get("wt_pct")):
        mass, how = mass * float(info["wt_pct"]) / 100.0, how + " × assay"
    return mass / mw, how


def preparations(entry: dict) -> list[dict]:
    preps = entry.get("preparations")
    return [p for p in preps if isinstance(p, dict)] if isinstance(preps, list) else []


def components(prep: dict) -> list[dict]:
    comps = prep.get("components")
    return [c for c in comps if isinstance(c, dict)] if isinstance(comps, list) else []


def all_components(entry: dict) -> list[tuple[str, dict]]:
    return [(str(p.get("label") or i + 1), c)
            for i, p in enumerate(preparations(entry)) for c in components(p)]


def eta_ul_per_mg(entry: dict) -> float | None:
    """Liquid-to-solid ratio η for a milled entry, when every solid mass is known."""
    app = entry.get("apparatus")
    if entry.get("method") != "ball_mill" or not isinstance(app, dict):
        return None
    solids = [c for _, c in all_components(entry) if c.get("role") != "solvent"]
    if not is_number(app.get("liquid_ul")) or not solids:
        return None
    if not all(is_number(c.get("mass_mg")) for c in solids):
        return None
    total = sum(float(c["mass_mg"]) for c in solids)
    return float(app["liquid_ul"]) / total if total else None


def computed(log: Log, entry: dict) -> dict:
    """Everything derived from what was entered.  Nothing here is ever stored."""
    rows, preps = [], []
    totals: dict[str, float] = {}
    total_ml = 0.0
    for i, prep in enumerate(preparations(entry)):
        label = str(prep.get("label") or i + 1)
        solvent_ml = sum(float(c["volume_ml"]) for c in components(prep)
                         if c.get("role") == "solvent" and is_number(c.get("volume_ml")))
        total_ml += solvent_ml
        conc = {}
        for c in components(prep):
            mmol, how = component_mmol(log, c)
            rows.append({"preparation": label, "code": c.get("code"), "role": c.get("role"),
                         "mmol": mmol, "from": how, "ratio": c.get("ratio")})
            if c.get("role") != "solvent" and mmol is not None:
                totals[str(c.get("code"))] = totals.get(str(c.get("code")), 0.0) + mmol
                if solvent_ml:
                    conc[str(c.get("code"))] = mmol / solvent_ml          # mmol/mL = mol/L
        preps.append({"label": label, "solvent_ml": solvent_ml or None, "conc_m": conc})
    ref = next((r for r in rows if r["role"] != "solvent" and r["mmol"]), None)
    for r in rows:
        r["eq"] = (r["mmol"] / ref["mmol"] if ref and r["mmol"] is not None
                   and r["role"] != "solvent" else None)
    final = {code: {"mmol": m, "conc_m": (m / total_ml if total_ml else None)}
             for code, m in totals.items()}
    return {"components": rows, "preparations": preps,
            "mixture": {"solvent_ml": total_ml or None, "solutes": final},
            "eta_ul_per_mg": eta_ul_per_mg(entry)}


def _check_amounts(log: Log, entry: dict, out: Findings) -> None:
    measured = []
    for p_index, prep in enumerate(preparations(entry)):
        plabel = str(prep.get("label") or p_index + 1)
        for c_index, comp in enumerate(components(prep)):
            name = f"preparations[{plabel}].components[{comp.get('code') or c_index + 1}]"
            solvent = comp.get("role") == "solvent"
            if solvent and f"{name}.ratio" in out.pending:
                out.pending.remove(f"{name}.ratio")
            if not any(recorded(comp.get(k)) for k in ("mass_mg", "volume_ml", "mmol")):
                out.pending.append(f"{name}.amount (mass_mg, volume_ml or mmol)")
            mmol, _ = component_mmol(log, comp)
            if not solvent and mmol is not None and is_number(comp.get("ratio")) \
                    and comp["ratio"]:
                measured.append((name, mmol, float(comp["ratio"])))
    if len(measured) < 2:
        return
    _, ref_mmol, ref_ratio = measured[0]
    for name, mmol, ratio in measured[1:]:
        stated, actual = ratio / ref_ratio, mmol / ref_mmol
        if abs(actual / stated - 1) > RATIO_TOLERANCE:
            out.problems.append(f"{name}: ratio says {stated:.3g} eq of the first reagent, "
                                f"amounts give {actual:.3g}")


# ---- checks -----------------------------------------------------------------------------

def check_entry(log: Log, entry: dict, entry_id: str | None = None) -> Findings:
    out = Findings()
    spec = entry_spec(entry.get("method"))
    skip: tuple = ()
    if entry.get("status") in ("planned", "abandoned"):
        skip = tuple(spec.fields[k] for k in OUTCOME_SECTIONS)
    walk(spec, entry, "", log, out, skip)
    if entry_id is not None and entry.get("id") != entry_id:
        out.problems.append(f"id: {entry.get('id')!r} does not match the file name {entry_id}")
    _check_amounts(log, entry, out)
    return out


def check_table(log: Log, spec: Block, table: dict) -> dict[str, Findings]:
    found = {}
    for code, info in table.items():
        out = Findings()
        walk(spec, info, "", log, out)
        formula = info.get("formula") if isinstance(info, dict) else None
        if spec is REAGENT and isinstance(formula, str) and formula and formula not in SENTINELS:
            try:
                molar_mass(formula)
            except ValueError as exc:
                out.problems.append(f"formula: {exc}")
        found[str(code)] = out
    return found


def check_reagents(log: Log) -> dict[str, Findings]:
    return check_table(log, REAGENT, log.reagents)


def check_phases(log: Log) -> dict[str, Findings]:
    return check_table(log, PHASE, log.phases)


def copy_recipe(source: dict, new_id: str, source_id: str) -> dict:
    """A new planned entry with `source`'s recipe and none of its outcome."""
    entry = json.loads(json.dumps(source))
    for key in ("title", "date", "notebook_ref", "notes", *OUTCOME_SECTIONS):
        entry.pop(key, None)
    entry.update(id=new_id, derived_from=source_id, status="planned")
    return entry
