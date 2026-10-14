"""Write entries, reagents and phases as commented YAML in spec order.

Blank blocks are written out field by field, so a file always shows what is still to record.
Keys the spec does not know are kept, marked, and reported by `check` — never dropped.
"""
from __future__ import annotations

import yaml

from .spec import (
    NA,
    PHASE,
    REAGENT,
    SENTINELS,
    UNKNOWN,
    Block,
    Leaf,
    entry_spec,
)

COMMENT_COL = 44

HEADER = f"""\
# Blank = not recorded yet (flagged by `expt check`).  {NA} = does not apply.
# {UNKNOWN} = cannot be recovered.  Amounts: enter what was weighed; mmol is computed.
"""


def _modeline(schema: str) -> str:
    return f"# yaml-language-server: $schema=../.schema/{schema}.schema.json\n"


def scalar(value: object) -> str:
    if value is None:
        return ""
    text = yaml.safe_dump(value, default_flow_style=True, allow_unicode=True, width=10_000)
    return text.removesuffix("...\n").strip()


def _comment(leaf: Leaf) -> str:
    parts = []
    if leaf.choices:
        parts.append(" | ".join(leaf.choices))
    if leaf.unit:
        parts.append(leaf.unit)
    if leaf.help:
        parts.append(leaf.help)
    if not leaf.required:
        parts.append("optional")
    return " — ".join(parts)


def _line(indent: int, text: str, comment: str = "") -> str:
    line = " " * indent + text
    if comment:
        line = line.ljust(COMMENT_COL - 1) + " # " + comment
    return line


def _emit(node, key: str, value: object, indent: int, out: list[str]) -> None:
    if isinstance(node, Leaf):
        if isinstance(value, str) and "\n" in value:
            out.append(_line(indent, f"{key}: |", _comment(node)))
            out.extend(" " * (indent + 2) + row for row in value.rstrip("\n").split("\n"))
        else:
            text = f"{key}: {scalar(value)}".rstrip()
            out.append(_line(indent, text, _comment(node)))
        return
    if isinstance(value, str) and value in SENTINELS:
        out.append(_line(indent, f"{key}: {value}", node.help))
        return
    if isinstance(node, Block):
        out.append(_line(indent, f"{key}:", node.help))
        _emit_fields(node, value if isinstance(value, dict) else {}, indent + 2, out)
        return
    items = value if isinstance(value, list) and value else [dict(s) for s in node.skeleton]
    out.append(_line(indent, f"{key}:", node.help))
    for item in items:
        body: list[str] = []
        _emit_fields(node.item, item if isinstance(item, dict) else {}, indent + 2, body)
        body[0] = " " * indent + "- " + body[0][indent + 2:]
        out.extend(body)


def _emit_fields(block: Block, values: dict, indent: int, out: list[str]) -> None:
    for key, node in block.fields.items():
        _emit(node, key, values.get(key), indent, out)
    for key in values:
        if key not in block.fields:
            out.append(_line(indent, f"{key}: {scalar(values[key])}".rstrip(), "not in the spec"))


def render_entry(values: dict) -> str:
    out: list[str] = []
    _emit_fields(entry_spec(values.get("method")), values, 0, out)
    return _modeline("entry") + HEADER + "\n" + "\n".join(out) + "\n"


def _render_keyed(spec: Block, table: dict, schema: str, intro: str) -> str:
    out: list[str] = []
    for code, values in table.items():
        out.append(f"{code}:")
        _emit_fields(spec, values if isinstance(values, dict) else {}, 2, out)
        out.append("")
    return _modeline(schema) + HEADER + intro + "\n" + "\n".join(out)


def render_reagents(table: dict) -> str:
    return _render_keyed(REAGENT, table, "reagents",
                         "# One block per reagent code; entries refer to these codes.\n")


def render_phases(table: dict) -> str:
    return _render_keyed(PHASE, table, "phases",
                         "# One block per crystalline phase; pxrd.phase in an entry names one.\n")
