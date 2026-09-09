"""M7b: the runner actually relaxes, and the ladder it builds is inspectable.

The regression this file exists for: rev 16 shipped the planner's pre-flight without the
executor.  `plan()` validated `run_mode`, `execute()` ignored it, and a 500-structure
`ml_go` run produced 281 RAW geometries, no ML row and no `mace` entry in `methods` while
reporting `done`.  Nothing failed; the option simply was not connected to anything.

Every test here runs against the Null backend, so the gate is the same on a laptop with
no quantum chemistry stack as on the workstation.
"""
from __future__ import annotations

import pytest

from mofsbu._types import Fidelity
from mofsbu.config import compute_device
from mofsbu.energy import backends as backends_mod
from mofsbu.energy import relax as relax_mod
from mofsbu.energy.backends import NullBackend
from mofsbu.registry import Registry
from mofsbu.runner import queue_relax, read_xyz, run
from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


@pytest.fixture()
def null_ml(monkeypatch):
    """Serve every rung from the test double, so the wiring is what is under test."""
    monkeypatch.setattr(relax_mod, "backend_for",
                        lambda fidelity, **kw: NullBackend(fidelity=fidelity))
    return NullBackend(fidelity=Fidelity.ML)


def acetate_spec(**kw) -> BuildSpec:
    return BuildSpec(
        molecules=(MoleculeSpec(name="AcOH", smiles="CC(=O)O", max_deprotonations=1),),
        metals=(MetalSpec(symbol="Zn", oxidation_state=2, multiplicity=1),),
        coordination=(4,), ligands_per_metal=(1,), binding=("mono",), co_ligand="O",
        **kw)


# ── the gate ─────────────────────────────────────────────────────────────────

def test_ml_go_actually_produces_relaxed_geometries(reg, null_ml):
    result = run(reg, acetate_spec(run_mode="ml_go"))
    assert result["status"] in {"done", "partial"}

    kinds = dict(reg.conn.execute(
        "SELECT kind, COUNT(*) n FROM tasks WHERE run_id=? GROUP BY kind",
        (result["run_id"],)).fetchall())
    assert kinds.get("relax", 0) > 0, (
        "no relax tasks were queued — this is exactly the rev 16 failure: the mode was "
        "accepted and nothing downstream acted on it")

    by_fidelity = dict(reg.conn.execute(
        "SELECT fidelity, COUNT(*) n FROM geometries GROUP BY fidelity").fetchall())
    assert by_fidelity.get(int(Fidelity.ML), 0) > 0, (
        f"every geometry is still a construct: {by_fidelity}")


def test_the_relaxed_row_is_chained_to_the_construct_it_came_from(reg, null_ml):
    run(reg, acetate_spec(run_mode="ml_go"))
    rows = reg.conn.execute(
        "SELECT id, structure_id, fidelity, relaxed_from, energy FROM geometries "
        "WHERE fidelity=?", (int(Fidelity.ML),)).fetchall()
    assert rows
    for row in rows:
        assert row["relaxed_from"] is not None, "a relaxed geometry with no parent"
        parent = reg.conn.execute("SELECT structure_id, fidelity FROM geometries WHERE id=?",
                                  (row["relaxed_from"],)).fetchone()
        assert parent["structure_id"] == row["structure_id"]
        assert parent["fidelity"] < row["fidelity"], "the ladder must go upwards"
        assert row["energy"] is not None


def test_the_construct_is_never_replaced(reg, null_ml):
    """One identity, several realisations (§6.2) — the cheap one is kept."""
    run(reg, acetate_spec(run_mode="ml_go"))
    structure_id = reg.conn.execute(
        "SELECT structure_id FROM geometries WHERE fidelity=? LIMIT 1",
        (int(Fidelity.ML),)).fetchone()["structure_id"]
    rungs = [r["fidelity"] for r in reg.conn.execute(
        "SELECT fidelity FROM geometries WHERE structure_id=? ORDER BY fidelity",
        (structure_id,))]
    assert len(rungs) >= 2 and rungs[0] < int(Fidelity.ML)


def test_best_geometry_moves_up_to_the_relaxed_one(reg, null_ml):
    run(reg, acetate_spec(run_mode="ml_go"))
    row = reg.conn.execute(
        "SELECT best_fidelity FROM structures WHERE best_fidelity=? LIMIT 1",
        (int(Fidelity.ML),)).fetchone()
    assert row is not None, "best_geometry_id never followed the ladder up"


def test_every_relaxed_energy_carries_its_method(reg, null_ml):
    """Ground rule 3, at the point a number is first produced."""
    run(reg, acetate_spec(run_mode="ml_go"))
    orphans = reg.conn.execute(
        "SELECT COUNT(*) c FROM geometries WHERE energy IS NOT NULL AND method_id IS NULL"
    ).fetchone()["c"]
    assert orphans == 0
    codes = {r["code"] for r in reg.conn.execute("SELECT code FROM methods")}
    assert "null" in codes, f"the backend that ran left no method row: {codes}"


# ── failure isolation ────────────────────────────────────────────────────────

def test_a_failed_relax_costs_nothing_already_built(reg, monkeypatch):
    """The reason relaxation is its own task and not a step inside the build."""
    class Exploding(NullBackend):
        def relax(self, *a, **kw):
            raise RuntimeError("SCF did not converge")

    monkeypatch.setattr(relax_mod, "backend_for",
                        lambda fidelity, **kw: Exploding(fidelity=fidelity))
    result = run(reg, acetate_spec(run_mode="ml_go"))

    constructs = reg.conn.execute(
        "SELECT COUNT(*) c FROM geometries WHERE fidelity < ?", (int(Fidelity.ML),)
    ).fetchone()["c"]
    assert constructs > 0, "a failed relaxation destroyed the structures already built"
    failed = reg.conn.execute(
        "SELECT COUNT(*) c FROM tasks WHERE run_id=? AND kind='relax' AND status='failed'",
        (result["run_id"],)).fetchone()["c"]
    assert failed > 0
    assert reg.conn.execute(
        "SELECT COUNT(*) c FROM geometries WHERE fidelity=?", (int(Fidelity.ML),)
    ).fetchone()["c"] == 0


# ── construct mode stays untouched ───────────────────────────────────────────

def test_construct_mode_queues_no_relaxation(reg, null_ml):
    result = run(reg, acetate_spec(run_mode="construct"))
    assert reg.conn.execute(
        "SELECT COUNT(*) c FROM tasks WHERE run_id=? AND kind='relax'",
        (result["run_id"],)).fetchone()["c"] == 0


def test_queue_relax_refuses_to_chain_off_a_relax(reg):
    """Otherwise one relax queues another and the run never terminates."""
    class T:
        kind, run_id, payload = "relax", 1, {}

    class Out:
        structure_id, geometry_id = 1, 1

    assert queue_relax(reg, acetate_spec(run_mode="ml_go"), T(), Out()) is None


# ── plumbing ─────────────────────────────────────────────────────────────────

def test_xyz_round_trips_through_the_store():
    text = "2\nname\nO  0.00000000 0.0 0.0\nH  0.96000000 0.0 0.0\n"
    symbols, coords = read_xyz(text)
    assert symbols == ["O", "H"]
    assert coords[1][0] == pytest.approx(0.96)


def test_device_is_declared_not_detected(monkeypatch):
    monkeypatch.delenv("MOFSBU_DEVICE", raising=False)
    assert compute_device() == "cpu"
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda")
    assert compute_device() == "cuda"
    monkeypatch.setenv("MOFSBU_DEVICE", "gpu0")
    with pytest.raises(ValueError, match="not a device torch would recognise"):
        compute_device()


def test_the_device_reaches_the_stored_method(monkeypatch):
    """Two runs of one model on different hardware can differ; that must be visible."""
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda")
    spec = backends_mod.MACEBackend().method_spec(charge=0, multiplicity=1)
    assert spec.extras["device"] == "cuda"


def test_relaxation_is_marked_as_wired():
    """If this flag is ever flipped back, the modes must go unavailable with it."""
    assert relax_mod.RELAXATION_IS_EXECUTED is True
    assert relax_mod.mode_status()["dft_go"]["available"] is False
