"""The experiment-log page and its API: wet-lab records, kept apart from the registry.

Writes only under its own folder (`<data root>/experiments`): one JSON file per entry,
the reagent and phase tables, and INDEX.md.  It never opens the registry.  Every write
names the version it was based on and is refused if the file changed since (two tabs, two
people), and a removed entry is moved aside rather than deleted.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from fastapi import APIRouter, Body, HTTPException
from fastapi.responses import HTMLResponse, PlainTextResponse

from mofsbu.explog.model import (
    ID_PATTERN,
    Log,
    StaleVersion,
    check_entry,
    check_phases,
    check_reagents,
    computed,
    copy_recipe,
    load_log,
    reagent_mw,
    remove_entry,
    version_of,
    write_json,
)
from mofsbu.explog.render import card, diff, index, summary_row
from mofsbu.explog.spec import METHODS, describe_all

STATIC = Path(__file__).parent / "static"
JSON_BODY = Body(...)


def prune(value: Any) -> Any:
    """Drop blanks so a file holds only what was recorded; blank and absent mean the same."""
    if isinstance(value, dict):
        kept = {k: prune(v) for k, v in value.items()}
        return {k: v for k, v in kept.items() if v not in (None, "", {}, [])}
    if isinstance(value, list):
        kept = [prune(v) for v in value]
        return [v for v in kept if v not in (None, "", {}, [])]
    if isinstance(value, str):
        return value.strip()
    return value


def _findings(found) -> dict[str, list[str]]:
    return {"pending": found.pending, "problems": found.problems}


def build_explog_router(root: Path) -> APIRouter:
    router = APIRouter()

    def log() -> Log:
        return load_log(root)

    def entry_path(lg: Log, entry_id: str) -> Path:
        try:
            return lg.entry_path(entry_id)
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc

    def refresh_index() -> None:
        lg = log()
        (root / "INDEX.md").write_text(index(lg), encoding="utf-8")

    def save(path: Path, data: Any, version: str | None) -> str:
        try:
            new = write_json(path, data, version)
        except StaleVersion as exc:
            raise HTTPException(409, f"{exc}; reload to see the newer version") from exc
        refresh_index()
        return new

    @router.get("/experiments", response_class=HTMLResponse)
    def page() -> str:
        return (STATIC / "experiments.html").read_text(encoding="utf-8")

    @router.get("/api/explog/spec")
    def spec() -> dict[str, Any]:
        return describe_all()

    @router.get("/api/explog/entries")
    def entries() -> dict[str, Any]:
        lg = log()
        return {"root": str(root), "entries": [summary_row(lg, i) for i in lg.entries],
                "load_errors": lg.load_errors}

    @router.get("/api/explog/entries/{entry_id}")
    def get_entry(entry_id: str) -> dict[str, Any]:
        lg = log()
        path = entry_path(lg, entry_id)
        if entry_id not in lg.entries:
            raise HTTPException(404, f"no entry {entry_id}")
        e = lg.entries[entry_id]
        return {"entry": e, "version": version_of(path),
                "findings": _findings(check_entry(lg, e, entry_id)), "computed": computed(lg, e)}

    @router.post("/api/explog/entries")
    def create_entry(payload: dict = JSON_BODY) -> dict[str, Any]:
        lg = log()
        new_id = str(payload.get("id") or "").strip()
        path = entry_path(lg, new_id)
        if path.exists():
            raise HTTPException(409, f"entry {new_id} already exists")
        like, method = payload.get("like"), payload.get("method")
        if like:
            if like not in lg.entries:
                raise HTTPException(404, f"no entry {like} to copy")
            entry = copy_recipe(lg.entries[like], new_id, like)
        elif method in METHODS:
            entry = {"id": new_id, "status": "planned", "method": method}
        else:
            raise HTTPException(400, f"give a method ({', '.join(METHODS)}) or an entry to copy")
        return {"id": new_id, "version": save(path, entry, None)}

    @router.put("/api/explog/entries/{entry_id}")
    def put_entry(entry_id: str, payload: dict = JSON_BODY) -> dict[str, Any]:
        lg = log()
        path = entry_path(lg, entry_id)
        entry = prune(payload.get("entry") or {})
        if entry.get("id") != entry_id:
            raise HTTPException(400, f"the entry's id {entry.get('id')!r} is not {entry_id!r}")
        version = save(path, entry, payload.get("version"))
        lg = log()
        return {"version": version, "findings": _findings(check_entry(lg, entry, entry_id)),
                "computed": computed(lg, entry)}

    @router.delete("/api/explog/entries/{entry_id}")
    def delete_entry(entry_id: str) -> dict[str, Any]:
        lg = log()
        entry_path(lg, entry_id)
        try:
            moved = remove_entry(lg, entry_id)
        except FileNotFoundError as exc:
            raise HTTPException(404, f"no entry {entry_id}") from exc
        refresh_index()
        return {"moved_to": str(moved.relative_to(root))}

    @router.post("/api/explog/check")
    def check(payload: dict = JSON_BODY) -> dict[str, Any]:
        """Findings and computed amounts for a draft, without saving it.  Not pruned, so an
        item just added and still blank keeps its place in the paths the page highlights."""
        lg = log()
        entry = payload.get("entry") or {}
        return {"findings": _findings(check_entry(lg, entry, payload.get("id"))),
                "computed": computed(lg, entry)}

    @router.get("/api/explog/entries/{entry_id}/card", response_class=PlainTextResponse)
    def entry_card(entry_id: str) -> str:
        lg = log()
        if entry_id not in lg.entries:
            raise HTTPException(404, f"no entry {entry_id}")
        return card(lg, entry_id)

    @router.get("/api/explog/diff")
    def entry_diff(a: str, b: str) -> dict[str, Any]:
        lg = log()
        missing = [i for i in (a, b) if i not in lg.entries]
        if missing:
            raise HTTPException(404, f"no entry {', '.join(missing)}")
        return {"a": a, "b": b, "rows": diff(lg, a, b)}

    def table_routes(name: str, checker) -> None:
        def path(lg: Log) -> Path:
            return getattr(lg, f"{name}_path")

        def derived(lg: Log) -> dict[str, Any]:
            out = {"findings": {k: _findings(f) for k, f in checker(lg).items()}}
            if name == "reagents":
                out["molar_mass"] = {code: reagent_mw(lg, code) for code in lg.reagents}
            return out

        @router.get(f"/api/explog/{name}", name=f"get_{name}")
        def get_table() -> dict[str, Any]:
            lg = log()
            return {"table": getattr(lg, name), "version": version_of(path(lg)), **derived(lg)}

        @router.put(f"/api/explog/{name}", name=f"put_{name}")
        def put_table(payload: dict = JSON_BODY) -> dict[str, Any]:
            lg = log()
            table = payload.get("table")
            if not isinstance(table, dict):
                raise HTTPException(400, "expected a table of code -> fields")
            bad = [k for k in table if not ID_PATTERN.match(str(k))]
            if bad:
                raise HTTPException(400, f"not a valid code: {', '.join(map(str, bad))}")
            table = {k: prune(v) if isinstance(v, dict) else {} for k, v in table.items()}
            version = save(path(lg), table, payload.get("version"))
            return {"version": version, **derived(log())}

    table_routes("reagents", check_reagents)
    table_routes("phases", check_phases)

    @router.get("/api/explog/index.md", response_class=PlainTextResponse)
    def index_md() -> str:
        return index(log())

    return router
