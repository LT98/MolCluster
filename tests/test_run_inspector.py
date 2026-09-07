"""A run has to be readable after the console has scrolled away (rev 18).

Three things were invisible and are now recorded on the task row: WHY a candidate was
refused (structured, groupable, with atoms and elements), whether a task WROTE anything
or recognised an identity the registry already had, and what the embed actually did.
"""
from __future__ import annotations

import json

import pytest

from mofsbu.registry import BlobStore, Registry
from mofsbu.registry.jobs import (
    add_task, claim_task, complete_task, create_run, fail_task, outcome_summary,
    task_rows,
)
from mofsbu.runner import run
from mofsbu.spec import BuildSpec

THQ = "OC1=C(O)C(=O)C(O)=C(O)C1=O"


def spec_with(co_ligand: str | None, note: str = "") -> BuildSpec:
    return BuildSpec.from_dict({
        "spec_version": 3,
        "molecules": [{"name": "THQ", "smiles": THQ, "max_deprotonations": 1}],
        "metals": [{"symbol": "Fe", "oxidation_state": 3, "spin_class": "hs",
                    "multiplicity": 6}],
        "coordination": [6], "co_ligand": co_ligand, "binding": ["chelate"],
        "note": note,
    })


@pytest.fixture
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "store")) as r:
        r.migrate("test")
        yield r


# ── the record itself ────────────────────────────────────────────────────────

def test_rejection_carries_a_code_and_a_detail(reg):
    """Free text cannot be grouped: two clashes 0.01 A apart are one finding and two
    strings.  The code is what makes "×212 qc_clash" a single line."""
    run_id = create_run(reg, spec_with("O"))
    tid = add_task(reg, run_id, "place", {"x": 1})
    claim_task(reg, run_id)
    fail_task(reg, tid, "QC FAILED: ...", rejected=True, code="qc_clash",
              detail={"qc": {"elements": ["Fe", "I"]}})
    row = task_rows(reg, run_id)[0]
    assert row["status"] == "rejected"
    assert row["error_code"] == "qc_clash"
    assert row["detail"]["qc"]["elements"] == ["Fe", "I"]


def test_completion_records_whether_anything_was_written(reg):
    """D2 idempotency means a task that writes nothing is the SUCCESS case.  Until it
    was recorded, a run of pure duplicates looked exactly like a run of discoveries."""
    run_id = create_run(reg, spec_with("O"))
    a = add_task(reg, run_id, "place", {})
    b = add_task(reg, run_id, "place", {})
    claim_task(reg, run_id); claim_task(reg, run_id)
    complete_task(reg, a, structure_created=True, geometry_created=True)
    complete_task(reg, b, structure_created=False, geometry_created=False)
    s = outcome_summary(reg, run_id)
    assert s["new_structures"] == 1
    assert s["reused_structures"] == 1


def test_summary_groups_causes_and_counts_retries(reg):
    run_id = create_run(reg, spec_with("O"))
    ids = [add_task(reg, run_id, "place", {}) for _ in range(4)]
    for _ in ids:
        claim_task(reg, run_id)
    fail_task(reg, ids[0], "a", rejected=True, code="qc_clash")
    fail_task(reg, ids[1], "b", rejected=True, code="qc_clash")
    fail_task(reg, ids[2], "c", rejected=True, code="co_ligand_no_donor")
    fail_task(reg, ids[3], "boom", code="IndexError")
    by = {(g["status"], g["code"]): g["count"] for g in outcome_summary(reg, run_id)["by_code"]}
    assert by[("rejected", "qc_clash")] == 2
    assert by[("rejected", "co_ligand_no_donor")] == 1
    assert by[("failed", "IndexError")] == 1


# ── through a real run ───────────────────────────────────────────────────────

def test_a_real_run_records_distances_and_their_provenance(reg):
    """The number that decided the geometry, and whether it was calibrated, are stored
    on the task — so a QC verdict can be read against the distance that produced it."""
    result = run(reg, spec_with("[I-]", "iodide co-ligand"))
    done = [t for t in task_rows(reg, result["run_id"])
            if t["status"] == "done" and t["kind"] == "place"]
    assert done, "the iodide sphere should build at all — that is the rev 18 fix"
    dists = {(d["metal"], d["donor"]): d for d in done[0]["detail"]["distances"]}
    assert dists[("Fe", "I")]["d"] == pytest.approx(2.70, abs=0.01)
    assert dists[("Fe", "O")]["d"] == pytest.approx(2.10, abs=0.01)
    assert all(not d["estimated"] for d in dists.values())


def test_a_co_ligand_with_no_donor_is_an_answer_not_a_crash(reg):
    """This was `perceive(co)[0]` and an IndexError: a traceback, reported as `failed`,
    for what is a plain statement about the co-ligand."""
    result = run(reg, spec_with("CC", "ethane as co-ligand"))
    rejected = [t for t in task_rows(reg, result["run_id"]) if t["status"] == "rejected"]
    assert rejected
    assert all(t["error_code"] == "co_ligand_no_donor" for t in rejected)
    assert result["counts"].get("failed") is None            # not a machinery failure
    assert "hint" in rejected[0]["detail"]


def test_embed_provenance_is_recorded(reg):
    """Which force field relaxed it, and whether ETKDG had to fall back to random
    coordinates.  Both change what the geometry is worth and neither was visible."""
    result = run(reg, spec_with("O"))
    row = task_rows(reg, result["run_id"])[0]
    assert row["detail"]["embed"]["ff"] in ("MMFF94", "UFF")
    assert row["detail"]["embed"]["retried"] in (True, False)


# ── the API the page reads ───────────────────────────────────────────────────

def test_inspector_endpoints(tmp_path):
    from fastapi.testclient import TestClient

    from mofsbu.ui.app import create_app

    db, store = tmp_path / "r.db", tmp_path / "store"
    with Registry(db, BlobStore(store)) as reg:
        reg.migrate("test")
        result = run(reg, spec_with("[I-]"))
    client = TestClient(create_app(db, store))

    assert client.get("/runs").status_code == 200

    body = client.get(f"/api/runs/{result['run_id']}").json()
    assert body["spec"]["co_ligand"] == "[I-]"
    assert "summary" in body and "planner_skips" in body

    tasks = client.get(f"/api/runs/{result['run_id']}/tasks").json()
    assert tasks["total"] == result["tasks"]
    assert {"payload", "detail", "error_code", "structure_created"} <= set(tasks["tasks"][0])

    assert client.get("/api/runs/999999").status_code == 404


def test_inspector_endpoints_do_not_write(tmp_path):
    """The builder writes specs and queues; the inspector reads.  It shares the router
    but must not share the write path — a page that can corrupt a running job's task
    table is not a diagnostic tool."""
    import sqlite3

    from fastapi.testclient import TestClient

    from mofsbu.ui.app import create_app

    db, store = tmp_path / "r.db", tmp_path / "store"
    with Registry(db, BlobStore(store)) as reg:
        reg.migrate("test")
        result = run(reg, spec_with("O"))
    client = TestClient(create_app(db, store))
    before = sqlite3.connect(db).execute("SELECT COUNT(*) FROM tasks").fetchone()[0]
    client.get(f"/api/runs/{result['run_id']}/tasks")
    client.get(f"/api/runs/{result['run_id']}")
    after = sqlite3.connect(db).execute("SELECT COUNT(*) FROM tasks").fetchone()[0]
    assert before == after
