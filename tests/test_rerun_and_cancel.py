"""Re-running an unchanged spec must not recompute, and a run must be stoppable.

The bug: `queue_relax` emitted a relax task for every completed build task without asking
whether that exact relaxation already existed.  Re-running an unchanged spec rebuilds the
same constructs, recognises them by identity (D2), hands back the geometry ids it already
had — and every one of them was queued for relaxation again.  In the user's registry that
showed as relax tasks pointing three deep at the same source geometry, and 124 completed
relaxations whose result `put_geometry` then discarded as a duplicate: the calculation had
already been paid for by the time anything noticed.

The distinction that matters: dedup is keyed on (source geometry, method), NOT on the
structure.  Twelve raw constructs of one identity are twelve different starting points and
may relax into different minima; collapsing those would throw away real chemistry.
"""
from __future__ import annotations

import pytest

from mofsbu._types import Fidelity
from mofsbu.energy import relax as relax_mod
from mofsbu.energy.backends import NullBackend
from mofsbu.registry import Registry
from mofsbu.registry.jobs import (
    cancel_requested, finalise_cancel, request_cancel, resume_run, task_counts,
)
from mofsbu.runner import existing_relaxation, plan, run, work
from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


@pytest.fixture()
def null_ml(monkeypatch):
    monkeypatch.setattr(relax_mod, "backend_for",
                        lambda fidelity, **kw: NullBackend(fidelity=fidelity))


def spec(**kw) -> BuildSpec:
    return BuildSpec(
        molecules=(MoleculeSpec(name="AcOH", smiles="CC(=O)O", max_deprotonations=1),),
        metals=(MetalSpec(symbol="Zn", oxidation_state=2),),
        coordination=(4,), ligands_per_metal=(1,), binding=("mono",), co_ligand="O", **kw)


def _counts(reg):
    return {
        "relax_tasks": reg.conn.execute(
            "SELECT COUNT(*) c FROM tasks WHERE kind='relax'").fetchone()["c"],
        "relaxed_geoms": reg.conn.execute(
            "SELECT COUNT(*) c FROM geometries WHERE relaxed_from IS NOT NULL"
        ).fetchone()["c"],
    }


# ── the fix ──────────────────────────────────────────────────────────────────

def test_rerunning_an_unchanged_spec_queues_no_new_relaxations(reg, null_ml):
    """The headline: the second run must cost nothing."""
    run(reg, spec(run_mode="ml_go"))
    first = _counts(reg)
    assert first["relax_tasks"] > 0, "the first run must actually relax something"

    run(reg, spec(run_mode="ml_go"))
    second = _counts(reg)
    assert second["relax_tasks"] == first["relax_tasks"], (
        f"the re-run queued {second['relax_tasks'] - first['relax_tasks']} more "
        f"relaxations for work already done")
    assert second["relaxed_geoms"] == first["relaxed_geoms"]


def test_a_third_run_is_free_too(reg, null_ml):
    run(reg, spec(run_mode="ml_go"))
    run(reg, spec(run_mode="ml_go"))
    before = _counts(reg)
    run(reg, spec(run_mode="ml_go"))
    assert _counts(reg) == before


def test_dedup_keys_on_the_source_geometry_not_the_structure(reg, null_ml):
    """Two constructs of one identity are two starting points, and both deserve a relax.

    This is the line the fix must not cross.  Collapsing per structure would turn the
    twelve raw constructs of one identity into one relaxation and silently discard the
    second minimum — which is the information D10/D11 exist to preserve.
    """
    run(reg, spec(run_mode="ml_go"))
    row = reg.conn.execute(
        "SELECT structure_id, COUNT(*) n FROM geometries WHERE relaxed_from IS NOT NULL "
        "GROUP BY structure_id ORDER BY n DESC LIMIT 1").fetchone()
    if row is None or row["n"] < 2:
        pytest.skip("this fixture produced no structure with two distinct constructs")
    parents = [r["relaxed_from"] for r in reg.conn.execute(
        "SELECT relaxed_from FROM geometries WHERE structure_id=? AND relaxed_from IS NOT NULL",
        (row["structure_id"],))]
    assert len(set(parents)) == len(parents), "two relaxations shared a starting geometry"


def test_existing_relaxation_finds_what_was_computed(reg, null_ml):
    run(reg, spec(run_mode="ml_go"))
    row = reg.conn.execute(
        "SELECT id, structure_id, relaxed_from FROM geometries "
        "WHERE relaxed_from IS NOT NULL LIMIT 1").fetchone()
    found = existing_relaxation(reg, row["relaxed_from"], Fidelity.ML, row["structure_id"])
    assert found == row["id"]


def test_a_different_level_of_theory_is_not_a_duplicate(reg, null_ml):
    """xTB after MACE is new work, not a repeat, and must not be skipped."""
    run(reg, spec(run_mode="ml_go"))
    row = reg.conn.execute(
        "SELECT structure_id, relaxed_from FROM geometries WHERE relaxed_from IS NOT NULL "
        "LIMIT 1").fetchone()
    assert existing_relaxation(reg, row["relaxed_from"], Fidelity.XTB,
                               row["structure_id"]) is None


def test_asking_whether_it_was_computed_does_not_create_a_method_row(reg, null_ml):
    """`find_method_id` must not register the row it was asked about."""
    plan(reg, spec(run_mode="ml_go"))
    before = reg.conn.execute("SELECT COUNT(*) c FROM methods").fetchone()["c"]
    existing_relaxation(reg, 1, Fidelity.XTB, 1)
    assert reg.conn.execute("SELECT COUNT(*) c FROM methods").fetchone()["c"] == before


# ── the kill switch ──────────────────────────────────────────────────────────

def test_cancelling_stops_the_worker_and_keeps_what_finished(reg, null_ml):
    run_id, n = plan(reg, spec(run_mode="ml_go"))
    reg.conn.commit()
    work(reg, spec(run_mode="ml_go"), run_id, limit=3)      # a few tasks in
    built = reg.conn.execute("SELECT COUNT(*) c FROM structures").fetchone()["c"]
    assert built > 0

    assert request_cancel(reg, run_id) is True
    assert cancel_requested(reg, run_id) is True

    did = work(reg, spec(run_mode="ml_go"), run_id)          # must do nothing more
    assert did == 0, f"{did} task(s) ran after the stop was requested"
    assert reg.conn.execute(
        "SELECT COUNT(*) c FROM structures").fetchone()["c"] == built, (
        "cancelling destroyed structures that were already built")


def test_a_cancelled_run_is_not_a_failed_run(reg, null_ml):
    """`rejected` vs `failed` again: a build you stopped is not a build that broke."""
    run_id, _ = plan(reg, spec(run_mode="ml_go"))
    reg.conn.commit()
    request_cancel(reg, run_id)
    finalise_cancel(reg, run_id)
    status = reg.conn.execute("SELECT status FROM runs WHERE id=?", (run_id,)).fetchone()
    assert status["status"] == "cancelled"
    counts = task_counts(reg, run_id)
    assert counts.get("failed", 0) == 0
    assert counts.get("cancelled", 0) > 0


def test_finish_run_cannot_relabel_a_stop_as_done(reg, null_ml):
    """The worker loop exits normally on cancel; that must not read as success."""
    from mofsbu.registry.jobs import finish_run

    run_id, _ = plan(reg, spec(run_mode="ml_go"))
    reg.conn.commit()
    request_cancel(reg, run_id)
    assert finish_run(reg, run_id) == "cancelled"


def test_claim_refuses_while_a_stop_is_in_flight(reg, null_ml):
    """Defence in depth for a separate worker process mid-loop."""
    from mofsbu.registry.jobs import claim_task

    run_id, _ = plan(reg, spec(run_mode="ml_go"))
    reg.conn.commit()
    request_cancel(reg, run_id)
    assert claim_task(reg, run_id) is None


def test_a_cancelled_run_can_be_resumed(reg, null_ml):
    run_id, n = plan(reg, spec(run_mode="ml_go"))
    reg.conn.commit()
    request_cancel(reg, run_id)
    finalise_cancel(reg, run_id)
    revived = resume_run(reg, run_id)
    assert revived > 0
    assert work(reg, spec(run_mode="ml_go"), run_id, limit=2) == 2


def test_cancelling_a_finished_run_says_so_rather_than_pretending(reg, null_ml):
    out = run(reg, spec(run_mode="ml_go"))
    assert request_cancel(reg, out["run_id"]) is False
