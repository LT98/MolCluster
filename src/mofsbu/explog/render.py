"""Human-readable views: one entry (`show`), the whole log (`index`), two entries (`diff`)."""
from __future__ import annotations

import datetime as dt
from collections.abc import Callable

from .model import (
    Findings,
    Log,
    check_entry,
    check_phases,
    check_reagents,
    eta_ul_per_mg,
    is_number,
    line_mmol,
    reagent_lines,
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


def _unit(leaf: Leaf) -> str:
    return f" {leaf.unit}" if leaf.unit and leaf.kind != "temp" else (
        " °C" if leaf.kind == "temp" else "")


def flatten(node, value: object, path: str,
            name: Callable[[Items, dict, int, str], str] | None = None) -> list[tuple]:
    """(path, value, leaf) for every recorded leaf.  Sentinel blocks come back as one row."""
    if value is None:
        return []
    if isinstance(node, Leaf) or isinstance(value, str) and not isinstance(node, Leaf):
        return [(path, value, node if isinstance(node, Leaf) else None)]
    rows: list[tuple] = []
    if isinstance(node, Block) and isinstance(value, dict):
        for key, child in node.fields.items():
            rows += flatten(child, value.get(key), f"{path}.{key}" if path else key, name)
        rows += [(f"{path}.{k}" if path else k, v, None) for k, v in value.items()
                 if k not in node.fields and v is not None]
    elif isinstance(node, Items) and isinstance(value, list):
        for i, item in enumerate(value):
            label = name(node, item, i, path) if name else (
                str(item.get(node.key)) if isinstance(item, dict) and item.get(node.key)
                not in (None, "") else str(i + 1))
            rows += flatten(node.item, item, f"{path}[{label}]", name)
    return rows


def _section(entry: dict, key: str) -> list[str]:
    spec = entry_spec(entry.get("method")).fields[key]
    rows = flatten(spec, entry.get(key), "")
    return [f"- **{p}**: {fmt(v)}{_unit(leaf) if leaf and is_number(v) else ''}"
            for p, v, leaf in rows]


def _cells_table(entry: dict) -> list[str]:
    cells = entry.get("cells")
    if not isinstance(cells, list):
        return [f"- {fmt(cells)}"] if cells else []
    cells = [c for c in cells if isinstance(c, dict) and recorded(c)]
    if not cells:
        return []
    out = ["| Cell | Basis | Capacity mAh/g @ mA/g | ICE % | Retention % @ cycles | CE % "
           "| Rate recovery % |", "|---|---|---|---|---|---|---|"]
    for c in cells:
        if not isinstance(c, dict):
            continue
        r = c.get("results") if isinstance(c.get("results"), dict) else {}
        out.append(
            f"| {fmt(c.get('id'))} | {fmt(c.get('capacity_basis'))} "
            f"| {fmt(r.get('capacity_mah_g'))} @ {fmt(r.get('capacity_rate_ma_g'))} "
            f"| {fmt(r.get('ice_pct'))} | {fmt(r.get('retention_pct'))} @ "
            f"{fmt(r.get('retention_cycles'))} | {fmt(r.get('ce_pct'))} "
            f"| {fmt(r.get('rate_recovery_pct'))} |")
        ladder = r.get("rate_ladder")
        if isinstance(ladder, list) and any(isinstance(s, dict) and s.get("rate_ma_g")
                                            for s in ladder):
            steps = ", ".join(f"{fmt(s.get('rate_ma_g'))}→{fmt(s.get('capacity_mah_g'))}"
                              for s in ladder if isinstance(s, dict))
            out.append(f"|  | rate ladder (mA/g→mAh/g): {steps} ||||||")
    return out


def _findings(found: Findings) -> list[str]:
    out = []
    if found.problems:
        out += ["", f"## Problems ({len(found.problems)})", ""]
        out += [f"- {p}" for p in found.problems]
    if found.pending:
        out += ["", f"## Still to record ({len(found.pending)})", ""]
        out += [f"- {p}" for p in found.pending]
    return out


def show(log: Log, entry_id: str) -> str:
    entry = log.entries[entry_id]
    found = check_entry(log, entry_id)
    head = f"# {entry_id}" + (f" — {entry['title']}" if entry.get("title") else "")
    meta = " · ".join(f"{k}: {fmt(entry[k])}" for k in
                      ("status", "date", "method", "operator", "notebook_ref", "derived_from")
                      if entry.get(k) not in (None, ""))
    out = [head, "", meta, "", "## Reagents", "",
           "| Code | Role | Mass mg | Volume mL | mmol | Stated eq | Eq from amounts | Note |",
           "|---|---|---|---|---|---|---|---|"]
    lines = reagent_lines(entry)
    mmols = [line_mmol(log, r) for r in lines]
    ref = next((m for r, m in zip(lines, mmols, strict=True)
                if m and r.get("role") != "solvent"), None)
    for r, m in zip(lines, mmols, strict=True):
        calc = "" if is_number(r.get("mmol")) or m is None else " (calc)"
        eq = f"{m / ref:.3g}" if m and ref and r.get("role") != "solvent" else ""
        out.append(f"| {fmt(r.get('code'))} | {fmt(r.get('role'))} | {fmt(r.get('mass_mg'))} "
                   f"| {fmt(r.get('volume_ml'))} | {fmt(m)}{calc} | {fmt(r.get('ratio'))} "
                   f"| {eq} | {fmt(r.get('note'))} |")
    eta = eta_ul_per_mg(entry)
    for key, title in (("conditions", "Conditions"), ("workup", "Workup"),
                       ("product", "Product"), ("characterization", "Characterization")):
        body = _section(entry, key)
        if key == "conditions" and eta is not None:
            body.append(f"- **η (computed)**: {eta:.3g} µL/mg")
        if body:
            out += ["", f"## {title}", ""] + body
    cells = _cells_table(entry)
    if cells:
        out += ["", "## Cells", ""] + cells
    if entry.get("notes"):
        out += ["", "## Notes", "", str(entry["notes"]).rstrip()]
    return "\n".join(out + _findings(found)) + "\n"


def _roles(entry: dict, role: str) -> str:
    return ", ".join(f"{r.get('code')} ({fmt(r.get('ratio'))})" if is_number(r.get("ratio"))
                     else str(r.get("code")) for r in reagent_lines(entry)
                     if r.get("role") == role and r.get("code"))


def _cell_summary(entry: dict, field: str, extra: str = "") -> str:
    cells = entry.get("cells")
    if not isinstance(cells, list):
        return fmt(cells)
    vals, extras = [], set()
    for c in cells:
        r = c.get("results") if isinstance(c, dict) else None
        if isinstance(r, dict) and is_number(r.get(field)):
            vals.append(float(r[field]))
            if extra and r.get(extra) is not None:
                extras.add(fmt(r[extra]))
    if not vals:
        return ""
    text = f"{sum(vals) / len(vals):.4g}" + (f" (n={len(vals)})" if len(vals) > 1 else "")
    return text + (f" @ {'/'.join(sorted(extras))}" if extras else "")


def _method(entry: dict) -> str:
    start = _roles(entry, "substrate")
    return fmt(entry.get("method")) + (f" of {start}" if start else "")


def _phase(entry: dict) -> str:
    ch = entry.get("characterization")
    pxrd = ch.get("pxrd") if isinstance(ch, dict) else None
    return fmt(pxrd.get("phase")) if isinstance(pxrd, dict) else fmt(pxrd)


def index(log: Log) -> str:
    out = ["# Experiment log", "",
           f"Written by `expt index` on {dt.date.today().isoformat()} — edit `entries/*.yaml`, "
           "not this file.", "",
           "| ID | Status | Date | Method | Metal | Ligand | Additive | Phase "
           "| Capacity mAh/g | Retention % @ cycles | To record | Problems |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    findings = {i: check_entry(log, i) for i in log.entries}
    for i, e in log.entries.items():
        f = findings[i]
        out.append(
            f"| {i} | {fmt(e.get('status'))} | {fmt(e.get('date'))} | {_method(e)} "
            f"| {_roles(e, 'metal_source')} | {_roles(e, 'ligand')} | {_roles(e, 'additive')} "
            f"| {_phase(e)} | {_cell_summary(e, 'capacity_mah_g', 'capacity_rate_ma_g')} "
            f"| {_cell_summary(e, 'retention_pct', 'retention_cycles')} "
            f"| {len(f.pending) or ''} "
            f"| {f'{WARN} {len(f.problems)}' if f.problems else ''} |")
    if log.phases:
        out += ["", "## Phases", "", "| Phase | Description | Defining 2θ | Entries |",
                "|---|---|---|---|"]
        for label, info in log.phases.items():
            info = info if isinstance(info, dict) else {}
            members = ", ".join(i for i, e in log.entries.items() if _phase(e) == str(label))
            out.append(f"| {label} | {fmt(info.get('description'))} "
                       f"| {fmt(info.get('peaks_2theta'))} | {members} |")
    problems = [(i, p) for i, f in findings.items() for p in f.problems]
    for name, table in (("reagents.yaml", check_reagents(log)), ("phases.yaml", check_phases(log))):
        problems += [(f"{name} {code}", p) for code, f in table.items() for p in f.problems]
    problems += [(name, err) for name, err in log.load_errors.items()]
    if problems:
        out += ["", f"## {WARN} Problems", ""] + [f"- **{i}** {p}" for i, p in problems]
    out += ["", "## Still to record", ""]
    for i, f in findings.items():
        if f.pending:
            out.append(f"<details><summary><b>{i}</b> — {len(f.pending)}</summary>\n")
            out += [f"- {p}" for p in f.pending] + ["", "</details>", ""]
    reagent_gaps = {c: f.pending for c, f in check_reagents(log).items() if f.pending}
    if reagent_gaps:
        out.append(f"<details><summary><b>reagents.yaml</b> — "
                   f"{sum(map(len, reagent_gaps.values()))}</summary>\n")
        out += [f"- {c}: {', '.join(g)}" for c, g in reagent_gaps.items()] + ["", "</details>"]
    return "\n".join(out) + "\n"


def _by_role(node: Items, item: dict, i: int, path: str) -> str:
    key = "role" if path == "reagents" else node.key
    if isinstance(item, dict) and item.get(key) not in (None, ""):
        return str(item[key])
    return str(i + 1)


def diff(log: Log, a: str, b: str) -> str:
    """Fields that differ.  Reagents are paired by role, so a swapped additive shows as one row."""
    rows = {}
    for side, entry_id in ((0, a), (1, b)):
        entry = log.entries[entry_id]
        seen: dict[str, int] = {}

        def name(node, item, i, path, _seen=seen):
            label = _by_role(node, item, i, path)
            if path == "reagents":
                _seen[label] = _seen.get(label, 0) + 1
                label = label if _seen[label] == 1 else f"{label}#{_seen[label]}"
            return label

        for path, value, _ in flatten(entry_spec(entry.get("method")), entry, "", name):
            rows.setdefault(path, [None, None])[side] = value
    changed = [(p, va, vb) for p, (va, vb) in rows.items()
               if va != vb and p not in ("id", "derived_from")]
    out = [f"# {a} → {b}", "", f"| Field | {a} | {b} |", "|---|---|---|"]
    out += [f"| {p} | {fmt(va)} | {fmt(vb)} |" for p, va, vb in changed]
    if not changed:
        out.append("| (no differences) | | |")
    return "\n".join(out) + "\n"
