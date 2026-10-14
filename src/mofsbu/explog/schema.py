"""JSON Schemas for editor support (VS Code + YAML extension): completion, dropdowns, red marks.

Generated from the spec and the current reagent and phase tables, so reagent codes and phase
labels appear as choices.  `expt check` remains the authority; the schema is a typing aid.
"""
from __future__ import annotations

import json
from pathlib import Path

from .model import Log
from .spec import CONDITIONS, METHODS, PHASE, REAGENT, SENTINELS, Block, Leaf, entry_spec

_BLANK = [{"type": "null"}, {"enum": list(SENTINELS)}]


def _describe(leaf: Leaf) -> str:
    return " — ".join(x for x in (leaf.unit, leaf.help) if x)


def _leaf(leaf: Leaf, log: Log) -> dict:
    kinds = {
        "num": {"type": "number"},
        "int": {"type": "integer"},
        "temp": {"anyOf": [{"type": "number"}, {"const": "RT"}]},
        "bool": {"type": "boolean"},
        "date": {"type": "string", "format": "date"},
        "enum": {"enum": list(leaf.choices)},
        "reagent": {"anyOf": [{"enum": sorted(map(str, log.reagents))},
                              {"type": "string", "pattern": "^product:"}]},
        "phase": {"enum": sorted(map(str, log.phases))},
        "entry": {"enum": sorted(log.entries)},
    }
    out = {"anyOf": [kinds.get(leaf.kind, {"type": ["string", "number"]})] + _BLANK}
    if _describe(leaf):
        out["description"] = _describe(leaf)
    return out


def _node(node, log: Log) -> dict:
    if isinstance(node, Leaf):
        return _leaf(node, log)
    if isinstance(node, Block):
        obj = {"type": "object", "additionalProperties": False,
               "properties": {k: _node(v, log) for k, v in node.fields.items()}}
        return {"anyOf": [obj] + _BLANK, **({"description": node.help} if node.help else {})}
    return {"anyOf": [{"type": "array", "items": _node(node.item, log)}] + _BLANK}


def entry_schema(log: Log) -> dict:
    base = _node(entry_spec("other"), log)["anyOf"][0]
    base["properties"]["conditions"] = {}
    base["allOf"] = [
        {"if": {"properties": {"method": {"const": m}}, "required": ["method"]},
         "then": {"properties": {"conditions": _node(CONDITIONS[m], log)}}}
        for m in METHODS]
    return {"$schema": "http://json-schema.org/draft-07/schema#", "title": "Experiment entry",
            **base}


def _keyed(spec: Block, log: Log, title: str) -> dict:
    return {"$schema": "http://json-schema.org/draft-07/schema#", "title": title,
            "type": "object", "additionalProperties": _node(spec, log)}


def write_schemas(log: Log) -> list[Path]:
    """Write `.schema/*.schema.json` under the log root; return the files that changed."""
    folder = log.root / ".schema"
    folder.mkdir(parents=True, exist_ok=True)
    changed = []
    for name, schema in (("entry", entry_schema(log)),
                         ("reagents", _keyed(REAGENT, log, "Reagents")),
                         ("phases", _keyed(PHASE, log, "Phases"))):
        path = folder / f"{name}.schema.json"
        text = json.dumps(schema, indent=1, ensure_ascii=False) + "\n"
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            path.write_text(text, encoding="utf-8")
            changed.append(path)
    return changed

