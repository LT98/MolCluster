"""The spec builder: a page that writes specs, so you never edit JSON by hand.

Boundary worth keeping straight.  The VIEWER (`ui/app.py`) is read-only over the
registry and stays that way — it opens `mode=ro` connections and cannot write.  The
BUILDER writes two things and nothing else: spec files and library entries on disk, and
rows in `runs`/`tasks` when you submit.  It never writes a structure, a geometry or a
label; those only ever arrive through `registry.api` from a worker.

That split is what lets a page submit hour-long jobs without becoming stateful: the run
state lives in the task table, so this page can be closed, reopened, or replaced by a
notebook without anything being lost.
"""
from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Body, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse

from mofsbu import library
from mofsbu.config import data_root
from mofsbu.energy.relax import ml_model_status, mode_status
from mofsbu.geometry.placer import GEOMETRIES
from mofsbu.spec import RUN_MODES, BuildSpec
from mofsbu._types import MofsbuError

STATIC = Path(__file__).parent / "static"
MAX_PREVIEW_ATOMS = 300


class _ReadOnlyRegistry:
    """Enough of `Registry` for the read-side job queries, over a read-only connection.

    `jobs.task_rows` and `jobs.outcome_summary` want a `Registry` because that is where
    the connection lives; they only ever read.  Handing them a real `Registry` here would
    open a WRITABLE connection from a page whose whole contract is that it does not write
    to the registry (see the module docstring).  This shim keeps that boundary intact and
    keeps the queries in one place rather than copied into the router.
    """

    __slots__ = ("conn",)

    def __init__(self, conn) -> None:
        self.conn = conn


def _depict(mol) -> str:
    """2D SVG depiction.  RDKit draws it; no chemistry library is needed in the browser."""
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    from rdkit import Chem

    flat = Chem.RemoveHs(Chem.Mol(mol))
    rdDepictor.Compute2DCoords(flat)
    drawer = rdMolDraw2D.MolDraw2DSVG(360, 260)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, flat)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def build_router(db_path: Path, store_path: Path, spec_dir: Path,
                 active: "ActiveDatabase | None" = None) -> APIRouter:
    """The write side of the app.

    `active` is the shared, switchable registry reference (`ui/active.py`).  Every
    handler below reads `current()` rather than closing over `db_path`, which is what
    lets the page point a run at a different database without restarting the server.
    The one exception is deliberate: `submit_run` snapshots the path BEFORE starting its
    thread, so switching databases mid-run cannot redirect a build already in flight.
    """
    from mofsbu.ui.active import ActiveDatabase, compute_state

    router = APIRouter()
    active = active or ActiveDatabase(db_path)
    running: dict[int, str] = {}

    def current() -> Path:
        return active.path

    @router.get("/builder", response_class=HTMLResponse)
    def builder_page() -> str:
        return (STATIC / "builder.html").read_text(encoding="utf-8")

    @router.get("/runs", response_class=HTMLResponse)
    def runs_page() -> str:
        """The run inspector.  Read-only over runs and tasks.

        This exists because the console was the only place a run's reasoning appeared,
        and a console that refreshes a progress line is not a place you can read a
        rejection.  Everything shown here is stored, so it is still readable tomorrow.
        """
        return (STATIC / "runs.html").read_text(encoding="utf-8")

    @router.get("/api/capabilities")
    def capabilities() -> dict[str, Any]:
        """What the pipeline can actually do.  The page renders from this, so an
        unimplemented option is labelled by the backend rather than by a hard-coded
        string in the markup that can drift out of date."""
        status = mode_status()
        models = ml_model_status()
        labels = {"construct": "Construct only",
                  "ml_go": "ML geometry optimisation (MACE-MP-0 / MACE-OMOL-0)",
                  "xtb_go": "GFN2-xTB geometry optimisation",
                  "dft_go": "DFT geometry optimisation"}
        return {
            "run_modes": [
                {"id": m, "available": status.get(m, {}).get("available", False),
                 "label": labels[m],
                 "backend": status.get(m, {}).get("backend"),
                 # The THEORY, not just the rung: `ml_go` is served by two models.
                 "method": status.get(m, {}).get("method"),
                 "note": status.get(m, {}).get("note", "")}
                for m in RUN_MODES],
            # Which theory the ML rung would be, not just whether it can run.  The page
            # renders the picker from this, so a model that is not installed here cannot
            # be selected and a model that is gets named rather than called "MACE".
            "ml_models": [{"id": key, **entry} for key, entry in models.items()],
            "max_degree_implemented": 1,
            "degree_note": ("degree 1 = complete one centre's coordination sphere. "
                            "degree 2+ = polynuclear / extended growth — not implemented, "
                            "needs M5 (assembly.join) and M6 (multi-centre placer)."),
            "geometries": sorted(GEOMETRIES),
            "metals_optional": True,
            "metals_note": ("metal centres are optional: a purely molecular construction "
                            "(COF, organic cage) is a first-class case. Joining molecule "
                            "to molecule is still M5, so metal-free runs currently produce "
                            "activation states and sites only."),
            "spec_dir": str(spec_dir),
            # Hardware, enumerated but never chosen for you.  See `/api/compute`.
            "compute": compute_state(),
            "database": str(current()),
        }

    # ── hardware: declared from the page, not from a shell ───────────────────

    @router.get("/api/compute")
    def get_compute() -> dict[str, Any]:
        """Which devices exist here, and which one this server is declared to use."""
        return compute_state()

    @router.post("/api/compute")
    def set_compute(payload: dict = Body(...)) -> dict[str, Any]:
        """Declare the device and worker count for subsequent runs.

        Possible at all because of one property of this architecture: `submit_run`
        executes on a background thread IN THIS PROCESS, and `compute_device()` reads
        the environment at execution time.  So setting it here reaches the work without
        any new process model.

        It remains a DECLARATION.  The list this page offers comes from asking torch
        what exists; picking from that list is a person choosing, which is the thing
        ground rule 9 requires.  Nothing here ever selects CUDA because a card is
        present.
        """
        from mofsbu.ui.active import declare_compute

        device = payload.get("device")
        workers = payload.get("workers")
        try:
            return declare_compute(device if device else None,
                                   int(workers) if workers else None)
        except (ValueError, TypeError) as exc:
            raise HTTPException(400, str(exc)) from exc

    # ── which database a run goes to ─────────────────────────────────────────

    def _selectable() -> list[Path]:
        """The databases this server will agree to point at: everything under the data
        root, plus the one it was LAUNCHED with, plus the one it is on now.

        A whitelist rather than a name, because the two are not interchangeable.  The
        launch database can live anywhere (`--db /somewhere/else/registry.db`), so it
        cannot be addressed as a name under the data root — and two databases in
        different folders can share a basename, which is exactly the case that would let
        a selection land on the wrong file.  Selection therefore quotes a full path, and
        the path must be one this function already listed: a value from the browser is
        never joined onto a root and opened.

        `db_path` is in the list permanently and that is the point.  Listing only
        `current()` made the set change as you moved: switching off a launch database
        that lived outside the data root dropped it from the whitelist, so the way back
        was a 404 and the only fix was restarting the server with the flag again.  A
        one-way switch is worse than no switch.
        """
        from mofsbu.ui.__main__ import discover_databases

        out: list[Path] = []
        seen: set[Path] = set()
        for path in [*discover_databases(), Path(db_path), current()]:
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                out.append(path)
        return out

    @router.get("/api/databases")
    def list_databases() -> dict[str, Any]:
        """Every registry under the data root, with a row count.

        The count is the point.  "Which one is my real one" is the question actually
        being asked, and a filename answers it only if you happened to name it well —
        `registry.db` and `registry_old.db` look identical and differ by 12,000 rows.
        Counted over a read-only connection, and a file that cannot be read is reported
        as such rather than omitted.
        """
        from mofsbu.ui.app import open_read_only

        here = current().resolve()
        out = []
        seen = set()
        for path in _selectable():
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            item: dict[str, Any] = {
                "name": path.name, "path": str(resolved),
                # Basenames collide across folders, so the page needs something it can
                # show that is actually unique.  The parent is enough and is readable;
                # the full path is in the title attribute.
                "folder": str(resolved.parent),
                "in_data_root": resolved.parent == data_root().resolve(),
                "active": resolved == here,
                "exists": path.exists(),
                "size_mb": round(path.stat().st_size / 1e6, 2) if path.exists() else 0.0,
            }
            if path.exists():
                try:
                    con = open_read_only(path)
                    try:
                        item["n_structures"] = con.execute(
                            "SELECT COUNT(*) FROM structures").fetchone()[0]
                        item["n_runs"] = con.execute(
                            "SELECT COUNT(*) FROM runs").fetchone()[0]
                    finally:
                        con.close()
                except Exception as exc:                           # noqa: BLE001
                    item["error"] = f"{type(exc).__name__}: {exc}"
            else:
                item["note"] = "not created yet — it appears on the first run"
            out.append(item)
        return {"databases": out, "active": str(here), "data_root": str(data_root())}

    @router.post("/api/databases")
    def create_database(payload: dict = Body(...)) -> dict[str, Any]:
        """Make a new, empty registry and switch to it.

        Cheap and complete: `Registry(path).migrate()` builds the whole schema from
        nothing, which is how every registry in this project has ever been made.  The
        name is validated as a path COMPONENT — the same rule `save_spec` applies, for
        the same reason.
        """
        from mofsbu.registry import BlobStore, Registry
        from mofsbu.ui.active import _validate_name

        try:
            name = _validate_name(payload.get("name") or "")
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc
        path = data_root() / name
        created = not path.exists()
        if created:
            with Registry(path, BlobStore(store_path)) as reg:
                reg.migrate("created from the builder page")
        active.switch(path)
        return {"database": str(path), "name": name, "created": created,
                "note": ("new, empty registry — it is now the active one" if created
                         else "that database already existed; switched to it")}

    @router.post("/api/database/select")
    def select_database(payload: dict = Body(...)) -> dict[str, Any]:
        """Point the whole app at an existing registry, named by its full path.

        The viewer, the run inspector and the builder all read the same reference, so
        this moves all three together — a page showing structures from one database and
        runs from another would be worse than no switch at all.  It does NOT widen
        anyone's access: the viewer still opens `mode=ro` connections on the new path.

        The path is matched against `_selectable()` rather than trusted.  That is what
        keeps an arbitrary string from the browser from being opened as a database,
        while still allowing the launch database, which may sit outside the data root
        and therefore cannot be addressed by name.
        """
        wanted = str(payload.get("path") or payload.get("name") or "").strip()
        if not wanted:
            raise HTTPException(400, "say which database")
        allowed = {str(p.resolve()): p for p in _selectable()}
        # Accept a bare name too, but only when it is unambiguous — which is the whole
        # reason the page sends a path.
        if wanted not in allowed:
            matches = [p for s, p in allowed.items() if p.name == wanted]
            if len(matches) == 1:
                wanted = str(matches[0].resolve())
            elif len(matches) > 1:
                raise HTTPException(409,
                    f"{len(matches)} databases are called {wanted!r} "
                    f"({', '.join(str(m.parent) for m in matches)}); say which by path")
            else:
                raise HTTPException(404, f"{wanted!r} is not one of the databases this "
                                         f"server offers; use 'new…' to create one")
        path = allowed[wanted]
        if not path.exists():
            raise HTTPException(404, f"{path} does not exist yet")
        active.switch(path)
        return {"database": str(path), "name": path.name, "switched": True}

    @router.post("/api/molecule/preview")
    def preview(payload: dict = Body(...)) -> dict[str, Any]:
        """Parse a SMILES and report what the pipeline sees, before anything is committed."""
        from mofsbu.geometry.embed import embed_molecule
        from mofsbu.graph.from_mol import mol_from_smiles
        from mofsbu.registry import formula as formula_of
        from mofsbu.graph.from_mol import from_rdkit
        from mofsbu.sites.model import chelate_pockets, perceive
        from mofsbu.sites.protomers import enumerate_protomers

        smiles = (payload.get("smiles") or "").strip()
        if not smiles:
            raise HTTPException(400, "no SMILES given")
        multiplicity = int(payload.get("multiplicity") or 1)
        try:
            mol = mol_from_smiles(smiles)
        except MofsbuError as exc:
            return {"ok": False, "error": str(exc)}
        if mol.GetNumAtoms() > MAX_PREVIEW_ATOMS:
            return {"ok": False, "error": f"{mol.GetNumAtoms()} atoms; preview caps at "
                                          f"{MAX_PREVIEW_ATOMS}"}
        try:
            graph = from_rdkit(mol, charge=None, multiplicity=multiplicity)
            embedded = embed_molecule(mol, seed=7)
            sites = perceive(embedded)
            pockets = chelate_pockets(embedded, sites)
            protomers = enumerate_protomers(mol, multiplicity=multiplicity)
            return {
                "ok": True, "smiles": smiles, "svg": _depict(mol),
                "formula": formula_of(graph), "n_atoms": len(graph),
                "charge": graph.charge,
                "donors": [{"idx": s.atom_idx, "type": s.donor_type,
                            "labile": s.labile, "live_dof": s.live_dof} for s in sites],
                "pockets": [{"donors": list(p.donors), "descriptor": p.descriptor,
                             "feasible_by": p.feasible_by} for p in pockets],
                "protomers": [{"charge": pr.charge, "label": pr.label,
                               "selections": len(pr.selections)} for pr in protomers],
            }
        except Exception as exc:                                   # noqa: BLE001
            return {"ok": False, "error": f"{type(exc).__name__}: {exc}"}

    # ── library ──────────────────────────────────────────────────────────────

    @router.get("/api/library")
    def get_library() -> dict[str, Any]:
        return {"molecules": [e.__dict__ for e in library.load()]}

    @router.post("/api/library")
    def add_to_library(payload: dict = Body(...)) -> dict[str, Any]:
        name = (payload.get("name") or "").strip()
        smiles = (payload.get("smiles") or "").strip()
        if not name or not smiles:
            raise HTTPException(400, "name and smiles are both required")
        try:
            entry = library.add(name, smiles,
                                multiplicity=int(payload.get("multiplicity") or 1),
                                note=payload.get("note") or "")
        except MofsbuError as exc:
            raise HTTPException(400, str(exc)) from exc
        return {"added": entry.__dict__, "path": str(library.LIBRARY_PATH)}

    @router.delete("/api/library/{name}")
    def remove_from_library(name: str) -> dict[str, Any]:
        return {"removed": library.remove(name)}

    # ── spec files ───────────────────────────────────────────────────────────

    def _spec_dirs() -> list[tuple[str, Path]]:
        """Where seed specs live: the tracked reference set, and what this UI has saved."""
        from mofsbu.config import REPO_ROOT

        return [("reference", REPO_ROOT / "data" / "reference"), ("saved", spec_dir)]

    @router.get("/api/specs")
    def list_specs() -> dict[str, Any]:
        out = []
        for source, directory in _spec_dirs():
            if not directory.exists():
                continue
            for path in sorted(directory.glob("*.json")):
                try:
                    raw = json.loads(path.read_text(encoding="utf-8"))
                    if not BuildSpec.looks_like_spec(raw):
                        continue          # the molecule library and friends are not specs
                    spec = BuildSpec.from_dict(raw)
                except Exception as exc:                           # noqa: BLE001
                    out.append({"source": source, "name": path.name, "ok": False,
                                "error": f"{type(exc).__name__}: {exc}"})
                    continue
                out.append({
                    "source": source, "name": path.name, "ok": True,
                    "digest": spec.digest[:12], "note": spec.note,
                    "molecules": [m.name for m in spec.molecules],
                    "metals": [m.symbol for m in spec.metals],
                    "degree": spec.degree, "run_mode": spec.run_mode,
                    "ml_model": spec.ml_model,
                })
        return {"specs": out}

    @router.get("/api/specs/{source}/{name}")
    def read_spec(source: str, name: str) -> dict[str, Any]:
        if "/" in name or "\\" in name or ".." in name:
            raise HTTPException(400, "bad filename")
        for label, directory in _spec_dirs():
            if label != source:
                continue
            path = directory / name
            if not path.exists():
                raise HTTPException(404, f"no spec {name} in {label}")
            return {"name": name, "source": label, "spec": BuildSpec.load(path).to_dict()}
        raise HTTPException(404, f"unknown source {source!r}")

    # ── specs and runs ───────────────────────────────────────────────────────

    def _parse(payload: dict) -> BuildSpec:
        try:
            return BuildSpec.from_dict(payload)
        except (TypeError, ValueError) as exc:
            raise HTTPException(400, f"bad spec: {exc}") from exc

    @router.post("/api/spec")
    def save_spec(payload: dict = Body(...)) -> dict[str, Any]:
        spec = _parse(payload.get("spec") or {})
        name = (payload.get("filename") or f"spec_{spec.digest[:8]}.json").strip()
        if "/" in name or "\\" in name:
            raise HTTPException(400, "filename must not contain a path")
        path = spec.save(spec_dir / name)
        return {"path": str(path), "digest": spec.digest}

    @router.post("/api/runs")
    def submit_run(payload: dict = Body(...)) -> JSONResponse:
        """Plan the run, then execute it on a background thread.

        Returns as soon as the tasks exist, so the page never waits on the work.  Progress
        is read back from the task table, which is where it lives.

        The optional `device` field is the page's half of ground rule 9.  The value is
        applied to this process (`compute_device()` reads the environment at execution
        time, and the work happens on a thread here) AND written to the run row, so a
        stored run still says which hardware produced it.  A page that only applied it
        would leave every result unattributable.
        """
        from mofsbu.registry import BlobStore, Registry
        from mofsbu.runner import plan, work

        spec = _parse(payload.get("spec") or {})
        device = payload.get("device") or None
        workers = payload.get("workers")
        if device is not None or workers is not None:
            from mofsbu.ui.active import declare_compute

            try:
                declare_compute(device, int(workers) if workers else None)
            except (ValueError, TypeError) as exc:
                raise HTTPException(400, str(exc)) from exc

        # Snapshot the database BEFORE the thread starts: switching databases while a
        # run is executing must not redirect it into a different file part-way through.
        run_db = current()
        with Registry(run_db, BlobStore(store_path)) as reg:
            reg.migrate("builder")
            try:
                run_id, n_tasks = plan(reg, spec)
                from mofsbu.registry.jobs import get_diagnostics

                diagnostics = get_diagnostics(reg, run_id)
            except NotImplementedError as exc:
                raise HTTPException(501, str(exc)) from exc
            except MofsbuError as exc:
                raise HTTPException(400, str(exc)) from exc

        def execute() -> None:
            from mofsbu.registry import relabel_all
            from mofsbu.registry.jobs import finish_run

            try:
                with Registry(run_db, BlobStore(store_path)) as r2:
                    work(r2, spec, run_id)
                    relabel_all(r2)
                    finish_run(r2, run_id)
                running.pop(run_id, None)
            except Exception as exc:                               # noqa: BLE001
                running[run_id] = f"{type(exc).__name__}: {exc}"

        running[run_id] = "running"
        threading.Thread(target=execute, daemon=True, name=f"mofsbu-run-{run_id}").start()
        from mofsbu.config import compute_device

        return JSONResponse({"run_id": run_id, "tasks": n_tasks, "digest": spec.digest,
                             "database": run_db.name, "device": compute_device(),
                             "diagnostics": diagnostics}, status_code=202)

    def _ro():
        from mofsbu.ui.app import open_read_only

        path = current()
        if not path.exists():
            raise HTTPException(503, f"registry not found: {path}")
        return open_read_only(path)

    @router.get("/api/runs/{run_id}")
    def get_run(run_id: int) -> dict[str, Any]:
        """One run: its spec, why candidates were not queued, and what happened to the
        ones that were.  The planner's skips and the workers' rejections are different
        facts about the same run and belong on one page."""
        from mofsbu.registry.jobs import get_diagnostics, outcome_summary

        con = _ro()
        try:
            shim = _ReadOnlyRegistry(con)
            row = con.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
            if row is None:
                raise HTTPException(404, f"no run {run_id}")
            item = dict(row)
            spec_json = item.pop("spec_json", "")
            try:
                spec = json.loads(spec_json)
            except json.JSONDecodeError:
                spec = None
            item["diagnostics"] = json.loads(item.pop("diagnostics_json", "[]") or "[]")
            item["thread"] = running.get(run_id, "")
            return {"run": item, "spec": spec,
                    "summary": outcome_summary(shim, run_id),
                    "planner_skips": get_diagnostics(shim, run_id)}
        finally:
            con.close()

    @router.post("/api/runs/{run_id}/cancel")
    def cancel_run(run_id: int) -> dict[str, Any]:
        """Ask a run to stop.  Cooperative, and it writes — so not the read-only path.

        Nothing is killed.  The run is flagged, workers finish the task in hand and stop
        claiming, and everything already computed stays in the registry.  A build stopped
        on purpose is recorded as `cancelled`, never as `failed`: a run list where every
        abandoned experiment reads as a crash is a run list nobody trusts.
        """
        from mofsbu.registry import BlobStore, Registry
        from mofsbu.registry.jobs import request_cancel, task_counts

        with Registry(current(), BlobStore(store_path)) as reg:
            row = reg.conn.execute("SELECT status FROM runs WHERE id=?",
                                   (run_id,)).fetchone()
            if row is None:
                raise HTTPException(404, f"no run {run_id}")
            changed = request_cancel(reg, run_id)
            counts = task_counts(reg, run_id)
        return {"run_id": run_id, "requested": changed, "was": row["status"],
                "pending": counts.get("pending", 0),
                "in_flight": counts.get("claimed", 0),
                "note": ("stopping — the task in hand will finish, then the workers stop "
                         "claiming" if changed else
                         f"nothing to stop: this run is {row['status']}")}

    @router.post("/api/runs/{run_id}/resume")
    def resume_run_route(run_id: int) -> dict[str, Any]:
        """Put a cancelled run's remaining tasks back in the queue."""
        from mofsbu.registry import BlobStore, Registry
        from mofsbu.registry.jobs import resume_run

        from mofsbu.registry.jobs import task_counts

        with Registry(current(), BlobStore(store_path)) as reg:
            was = reg.conn.execute("SELECT status FROM runs WHERE id=?",
                                   (run_id,)).fetchone()
            if was is None:
                raise HTTPException(404, f"no run {run_id}")
            revived = resume_run(reg, run_id)
            pending = task_counts(reg, run_id).get("pending", 0)
        note = (f"{pending} task(s) waiting — resubmit the spec to execute them"
                if pending else "nothing left to do in this run")
        return {"run_id": run_id, "revived": revived, "was": was["status"],
                "pending": pending, "note": note}

    @router.get("/api/runs/{run_id}/tasks")
    def get_run_tasks(run_id: int, status: str | None = None,
                      code: str | None = None,
                      limit: int = 300, offset: int = 0) -> dict[str, Any]:
        """The tasks themselves, with what was attempted and what came back.

        `code` filters by the grouping key, which is the whole reason it exists: after
        the summary says "×212 qc_clash", this is how you look at those 212.
        """
        from mofsbu.registry.jobs import task_rows

        con = _ro()
        try:
            shim = _ReadOnlyRegistry(con)
            rows = task_rows(shim, run_id, status=status,
                             limit=max(1, min(limit, 2000)), offset=max(0, offset))
            if code:
                rows = [r for r in rows if (r.get("error_code") or "(uncoded)") == code]
            total = con.execute("SELECT COUNT(*) FROM tasks WHERE run_id=?",
                                (run_id,)).fetchone()[0]
            return {"run_id": run_id, "total": total, "returned": len(rows), "tasks": rows}
        finally:
            con.close()

    @router.get("/api/runs")
    def list_runs(limit: int = 20) -> dict[str, Any]:
        from mofsbu.ui.app import open_read_only

        path = current()
        if not path.exists():
            return {"runs": []}
        con = open_read_only(path)
        try:
            # `device` arrived after some registries were written, and this connection is
            # read-only so it cannot migrate one.  Select it only where it exists.
            from mofsbu.ui.app import has_column

            device_col = "r.device" if has_column(con, "runs", "device") else "'' AS device"
            rows = con.execute(
                f"SELECT r.id, r.status, r.note, r.created_at, r.finished_at, {device_col},"
                "  (SELECT COUNT(*) FROM tasks t WHERE t.run_id = r.id) AS n_tasks,"
                "  (SELECT COUNT(*) FROM tasks t WHERE t.run_id = r.id AND t.status='done')"
                "    AS n_done,"
                "  (SELECT COUNT(*) FROM tasks t WHERE t.run_id = r.id"
                "    AND t.status='rejected') AS n_rejected,"
                "  (SELECT COUNT(*) FROM tasks t WHERE t.run_id = r.id AND t.status='failed')"
                "    AS n_failed, r.diagnostics_json"
                " FROM runs r ORDER BY r.id DESC LIMIT ?", (limit,)).fetchall()
            out = []
            import json as _json

            for row in rows:
                item = dict(row)
                item["thread"] = running.get(row["id"], "")
                item["diagnostics"] = _json.loads(item.pop("diagnostics_json", "[]") or "[]")
                out.append(item)
            return {"runs": out}
        finally:
            con.close()

    # ── one entry: re-run it, or hide it ─────────────────────────────────────
    # Two requests that look like one and are not.  Re-running is nearly free, because
    # `construct` is a deterministic function of the choice vector and the seed (D13) —
    # so "run this again" is replaying a task that is already stored.  Deleting is the
    # dangerous one, and the danger is not in the UI: see `registry.api.set_hidden`.

    @router.get("/api/structures/{structure_id}/origin")
    def structure_origin(structure_id: int) -> dict[str, Any]:
        """The task and run that produced this structure, and whether it can be re-run.

        Asked by the page before it offers the button, so a structure that arrived by
        some other route (ingested, seeded, built before tasks recorded `structure_id`)
        gets an explanation instead of a control that fails when pressed.
        """
        con = _ro()
        try:
            row = con.execute(
                "SELECT t.id AS task_id, t.run_id, t.kind, t.status, t.structure_created,"
                "       r.spec_json IS NOT NULL AS has_spec"
                " FROM tasks t JOIN runs r ON r.id = t.run_id"
                " WHERE t.structure_id = ?"
                " ORDER BY t.structure_created DESC, t.id ASC LIMIT 1", (structure_id,)
            ).fetchone()
            n_routes = con.execute(
                "SELECT COUNT(*) FROM reactions WHERE product_structure_id = ?",
                (structure_id,)).fetchone()[0]
            hidden = 0
            if any(c[1] == "hidden" for c in con.execute("PRAGMA table_info(structures)")):
                got = con.execute("SELECT hidden FROM structures WHERE id = ?",
                                  (structure_id,)).fetchone()
                hidden = int(got["hidden"]) if got else 0
            return {
                "structure_id": structure_id,
                "can_rerun": row is not None and bool(row["has_spec"]),
                "task_id": row["task_id"] if row else None,
                "run_id": row["run_id"] if row else None,
                "kind": row["kind"] if row else None,
                "why_not": None if row is not None else
                    "no task in this database records building this structure — it was "
                    "ingested or seeded rather than constructed, so there is nothing to "
                    "replay",
                "n_incoming_routes": n_routes,
                "hidden": bool(hidden),
                # The delete policy, decided here so the page does not have to encode it.
                "can_hide": n_routes <= 1,
                "hide_warning": None if n_routes <= 1 else
                    f"{n_routes} provenance routes arrive at this structure; hiding it "
                    f"removes it from listings but leaves every route intact",
            }
        finally:
            con.close()

    @router.post("/api/structures/{structure_id}/rerun")
    def rerun_structure(structure_id: int) -> JSONResponse:
        """Queue the one task that built this structure, again.

        Not a new enumeration: the originating task's payload and its run's spec are
        both stored, so this is a run of exactly one task with exactly the inputs that
        produced the row you are looking at.

        Under D2 the usual outcome is that the registry recognises the identity and
        writes nothing.  That is the CORRECT result and the page reports it as
        "already present" rather than as a no-op — a re-run that changes nothing is how
        you confirm determinism, and it should look like a confirmation.
        """
        from mofsbu.registry import BlobStore, Registry
        from mofsbu.registry.jobs import add_task, create_run
        from mofsbu.runner import work

        con = _ro()
        try:
            row = con.execute(
                "SELECT t.id AS task_id, t.run_id, t.kind, t.payload_json, r.spec_json"
                " FROM tasks t JOIN runs r ON r.id = t.run_id"
                " WHERE t.structure_id = ?"
                " ORDER BY t.structure_created DESC, t.id ASC LIMIT 1", (structure_id,)
            ).fetchone()
        finally:
            con.close()
        if row is None:
            raise HTTPException(404,
                f"nothing in this database records building structure {structure_id}, "
                f"so there is no task to replay")
        spec = _parse({**json.loads(row["spec_json"]),
                       "note": f"re-run of structure {structure_id} "
                               f"(task {row['task_id']}, run {row['run_id']})"})

        run_db = current()
        with Registry(run_db, BlobStore(store_path)) as reg:
            reg.migrate("builder")
            run_id = create_run(reg, spec, note=spec.note)
            task_id = add_task(reg, run_id, row["kind"], json.loads(row["payload_json"]))
            reg.conn.commit()

        def execute() -> None:
            from mofsbu.registry import relabel_all
            from mofsbu.registry.jobs import finish_run

            try:
                with Registry(run_db, BlobStore(store_path)) as r2:
                    work(r2, spec, run_id)
                    relabel_all(r2)
                    finish_run(r2, run_id)
                running.pop(run_id, None)
            except Exception as exc:                               # noqa: BLE001
                running[run_id] = f"{type(exc).__name__}: {exc}"

        running[run_id] = "running"
        threading.Thread(target=execute, daemon=True,
                         name=f"mofsbu-rerun-{run_id}").start()
        return JSONResponse({"run_id": run_id, "task_id": task_id, "tasks": 1,
                             "source_task": row["task_id"], "source_run": row["run_id"],
                             "structure_id": structure_id}, status_code=202)

    @router.post("/api/structures/{structure_id}/hide")
    def hide_structure(structure_id: int, payload: dict = Body(default={})) -> dict[str, Any]:
        """Hide a structure, or bring a hidden one back.  Nothing is ever deleted.

        Refused by default when more than one provenance route arrives at the structure,
        and refused WITH THE ROUTES rather than with a warning — the question "is this
        safe to remove" is answerable only by seeing what else reached it.  `force`
        overrides, and even then this is a flag on a row: the geometries, the blobs and
        every reaction edge survive, and unhiding is a second press of the same button.
        """
        from mofsbu.registry import BlobStore, Registry, incoming_routes, set_hidden
        from mofsbu.registry.api import RegistryError

        hidden = bool(payload.get("hidden", True))
        with Registry(current(), BlobStore(store_path)) as reg:
            reg.migrate("builder")
            routes = incoming_routes(reg, structure_id)
            if hidden and len(routes) > 1 and not payload.get("force"):
                raise HTTPException(409, {
                    "message": f"structure {structure_id} was reached by "
                               f"{len(routes)} different routes; hiding it removes it "
                               f"from listings only, but confirm you mean this one",
                    "routes": routes,
                })
            try:
                result = set_hidden(reg, structure_id, hidden,
                                    reason=str(payload.get("reason") or ""))
            except RegistryError as exc:
                raise HTTPException(404, str(exc)) from exc
            reg.conn.commit()
        result["routes_detail"] = routes
        result["note"] = ("hidden from listings — nothing was deleted, and the same "
                          "button brings it back" if hidden else "restored to the listing")
        return result

    return router
