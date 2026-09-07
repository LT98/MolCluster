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
from mofsbu.energy.relax import mode_status
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


def build_router(db_path: Path, store_path: Path, spec_dir: Path) -> APIRouter:
    router = APIRouter()
    running: dict[int, str] = {}

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
        labels = {"construct": "Construct only",
                  "ml_go": "ML geometry optimisation (MACE)",
                  "xtb_go": "GFN2-xTB geometry optimisation",
                  "dft_go": "DFT geometry optimisation"}
        return {
            "run_modes": [
                {"id": m, "available": status.get(m, {}).get("available", False),
                 "label": labels[m],
                 "backend": status.get(m, {}).get("backend"),
                 "note": status.get(m, {}).get("note", "")}
                for m in RUN_MODES],
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
        }

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
        """
        from mofsbu.registry import BlobStore, Registry
        from mofsbu.runner import plan, work

        spec = _parse(payload.get("spec") or {})
        with Registry(db_path, BlobStore(store_path)) as reg:
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
                with Registry(db_path, BlobStore(store_path)) as r2:
                    work(r2, spec, run_id)
                    relabel_all(r2)
                    finish_run(r2, run_id)
                running.pop(run_id, None)
            except Exception as exc:                               # noqa: BLE001
                running[run_id] = f"{type(exc).__name__}: {exc}"

        running[run_id] = "running"
        threading.Thread(target=execute, daemon=True, name=f"mofsbu-run-{run_id}").start()
        return JSONResponse({"run_id": run_id, "tasks": n_tasks, "digest": spec.digest,
                             "diagnostics": diagnostics}, status_code=202)

    def _ro():
        from mofsbu.ui.app import open_read_only

        if not db_path.exists():
            raise HTTPException(503, f"registry not found: {db_path}")
        return open_read_only(db_path)

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

        if not db_path.exists():
            return {"runs": []}
        con = open_read_only(db_path)
        try:
            rows = con.execute(
                "SELECT r.id, r.status, r.note, r.created_at, r.finished_at,"
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

    return router
