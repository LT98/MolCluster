"""Text views of the log: INDEX.md, one entry as markdown, and two entries side by side."""
from __future__ import annotations

import datetime as dt
from collections.abc import Callable

from .model import (
    Log,
    all_components,
    check_entry,
    check_phases,
    check_reagents,
    computed,
    is_number,
    item_name,
    recorded,
)
from .spec import Block, Items, Leaf, entry_spec

WARN = "⚠"


def fmt(value: object) -> str:
    if value is None:
        return ""
    if is_number(value):
        return f"{value:.4g}" if isinstance(value, float) else str(value)
    return str(value).replace("|", "\\|").replace("\n", " ")


def flatten(node, value: object, path: str,
            name: Callable[[Items, object, int, str], str] | None = None) -> list[tuple]:
    """(path, value, leaf) for every recorded leaf.  A sentinel block comes back as one row."""
    if not recorded(value):
        return []
    if isinstance(node, Leaf) or isinstance(value, str):
        return [(path, value, node if isinstance(node, Leaf) else None)]
    rows: list[tuple] = []
    if isinstance(node, Block) and isinstance(value, dict):
        for key, child in node.fields.items():
            rows += flatten(child, value.get(key), f"{path}.{key}" if path else key, name)
        rows += [(f"{path}.{k}" if path else k, v, None) for k, v in value.items()
                 if k not in node.fields and recorded(v)]
    elif isinstance(node, Items) and isinstance(value, list):
        for i, item in enumerate(value):
            label = name(node, item, i, path) if name else item_name(node, item, i)
            rows += flatten(node.item, item, f"{path}[{label}]", name)
    return rows


def _phase(entry: dict) -> str:
    ch = entry.get("characterization")
    pxrd = ch.get("pxrd") if isinstance(ch, dict) else None
    return fmt(pxrd.get("phase")) if isinstance(pxrd, dict) else fmt(pxrd)


def _roles(entry: dict, role: str) -> str:
    seen = []
    for _, c in all_components(entry):
        if c.get("role") == role and c.get("code"):
            text = f"{c['code']} ({fmt(c['ratio'])})" if is_number(c.get("ratio")) else c["code"]
            if text not in seen:
                seen.append(text)
    return ", ".join(seen)


def _cell_mean(entry: dict, key: str, extra: str = "") -> str:
    cells = entry.get("cells")
    if not isinstance(cells, list):
        return fmt(cells)
    vals, extras = [], set()
    for c in cells:
        r = c.get("results") if isinstance(c, dict) else None
        if isinstance(r, dict) and is_number(r.get(key)):
            vals.append(float(r[key]))
            if extra and r.get(extra) is not None:
                extras.add(fmt(r[extra]))
    if not vals:
        return ""
    text = f"{sum(vals) / len(vals):.4g}" + (f" (n={len(vals)})" if len(vals) > 1 else "")
    return text + (f" @ {'/'.join(sorted(extras))}" if extras else "")


def summary_row(log: Log, entry_id: str) -> dict:
    """What the page's list and INDEX.md show for one entry."""
    e = log.entries[entry_id]
    f = check_entry(log, e, entry_id)
    method = fmt(e.get("method"))
    start = _roles(e, "substrate")
    return {
        "id": entry_id, "title": e.get("title"), "status": e.get("status"),
        "date": e.get("date"), "method": method + (f" of {start}" if start else ""),
        "metal": _roles(e, "metal_source"), "ligand": _roles(e, "ligand"),
        "additive": _roles(e, "additive"), "phase": _phase(e),
        "capacity": _cell_mean(e, "capacity_mah_g", "capacity_rate_ma_g"),
        "retention": _cell_mean(e, "retention_pct", "retention_cycles"),
        "pending": len(f.pending), "problems": len(f.problems),
    }


def index(log: Log) -> str:
    out = ["# Experiment log", "",
           f"Rewritten on every save ({dt.date.today().isoformat()}). Edit through the "
           "`/experiments` page, not this file.", "",
           "| ID | Status | Date | Method | Metal | Ligand | Additive | Phase "
           "| Capacity mAh/g | Retention % @ cycles | To record | Problems |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i in log.entries:
        s = summary_row(log, i)
        flag = f"{WARN} {s['problems']}" if s["problems"] else ""
        out.append(
            f"| {i} | {fmt(s['status'])} | {fmt(s['date'])} | {s['method']} | {s['metal']} "
            f"| {s['ligand']} | {s['additive']} | {s['phase']} | {s['capacity']} "
            f"| {s['retention']} | {s['pending'] or ''} | {flag} |")
    if log.phases:
        out += ["", "## Phases", "", "| Phase | Description | Defining 2θ | Entries |",
                "|---|---|---|---|"]
        for label, info in log.phases.items():
            info = info if isinstance(info, dict) else {}
            members = ", ".join(i for i, e in log.entries.items() if _phase(e) == str(label))
            out.append(f"| {label} | {fmt(info.get('description'))} "
                       f"| {fmt(info.get('peaks_2theta'))} | {members} |")
    problems = [(i, p) for i, e in log.entries.items() for p in check_entry(log, e, i).problems]
    for name, table in (("reagent", check_reagents(log)), ("phase", check_phases(log))):
        problems += [(f"{name} {code}", p) for code, f in table.items() for p in f.problems]
    problems += list(log.load_errors.items())
    if problems:
        out += ["", f"## {WARN} Problems", ""] + [f"- **{i}** {p}" for i, p in problems]
    return "\n".join(out) + "\n"


def card(log: Log, entry_id: str) -> str:
    """One entry as markdown: recorded values, computed amounts, then what is missing."""
    e = log.entries[entry_id]
    found = check_entry(log, e, entry_id)
    calc = computed(log, e)
    out = [f"# {entry_id}" + (f" — {e['title']}" if e.get("title") else ""), "",
           " · ".join(f"{k}: {fmt(e[k])}" for k in ("status", "date", "method", "operator",
                                                     "notebook_ref", "derived_from")
                      if recorded(e.get(k))), "", "## Amounts (computed)", "",
           "| Preparation | Code | Role | mmol | from | eq | stated |",
           "|---|---|---|---|---|---|---|"]
    for r in calc["components"]:
        out.append(f"| {fmt(r['preparation'])} | {fmt(r['code'])} | {fmt(r['role'])} "
                   f"| {fmt(r['mmol'])} | {r['from']} | {fmt(r['eq'])} | {fmt(r['ratio'])} |")
    for p in calc["preparations"]:
        if p["conc_m"]:
            conc = ", ".join(f"{c} {m:.3g} M" for c, m in p["conc_m"].items())
            out.append(f"\n- **{p['label']}**: {fmt(p['solvent_ml'])} mL solvent; {conc}")
    mix = calc["mixture"]
    if mix["solvent_ml"]:
        conc = ", ".join(f"{c} {v['conc_m']:.3g} M" for c, v in mix["solutes"].items()
                         if v["conc_m"] is not None)
        out.append(f"- **reaction mixture**: {fmt(mix['solvent_ml'])} mL solvent; {conc}")
    if calc["eta_ul_per_mg"] is not None:
        out.append(f"- **η**: {calc['eta_ul_per_mg']:.3g} µL/mg")
    out += ["", "## Recorded", ""]
    out += [f"- **{p}**: {fmt(v)}{' ' + leaf.unit if leaf and leaf.unit and is_number(v) else ''}"
            for p, v, leaf in flatten(entry_spec(e.get("method")), e, "")
            if not p.startswith(("id", "notes"))]
    if e.get("notes"):
        out += ["", "## Notes", "", str(e["notes"]).rstrip()]
    if found.problems:
        out += ["", f"## Problems ({len(found.problems)})", ""] + [f"- {p}" for p in found.problems]
    if found.pending:
        out += ["", f"## Still to record ({len(found.pending)})", ""]
        out += [f"- {p}" for p in found.pending]
    return "\n".join(out) + "\n"


def diff(log: Log, a: str, b: str) -> list[dict]:
    """Fields that differ.  Preparations pair by position and components by role within
    them, so a swapped additive is one row rather than a removal and an addition."""
    rows: dict[str, list] = {}
    for side, entry_id in ((0, a), (1, b)):
        entry = log.entries[entry_id]
        seen: dict[str, int] = {}

        def name(node, item, i, path, _seen=seen):
            if path.endswith(".components") and isinstance(item, dict) and item.get("role"):
                key = f"{path}:{item['role']}"
                _seen[key] = _seen.get(key, 0) + 1
                return item["role"] + ("" if _seen[key] == 1 else f"#{_seen[key]}")
            if path in ("preparations", "steps", "cells", "workup.washes"):
                return str(i + 1)
            return item_name(node, item, i)

        for path, value, _ in flatten(entry_spec(entry.get("method")), entry, "", name):
            rows.setdefault(path, [None, None])[side] = value
    return [{"field": p, "a": va, "b": vb} for p, (va, vb) in rows.items()
            if va != vb and p not in ("id", "derived_from")]
