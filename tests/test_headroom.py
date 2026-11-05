"""The headroom check: a worker short of memory leaves the task pending instead of failing it."""
from __future__ import annotations

import pytest

from mofsbu import config, headroom
from mofsbu._types import MofsbuError
from mofsbu.registry import BlobStore, Registry
from mofsbu.registry.jobs import add_task, create_run, task_counts
from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec, PocketPredicate


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as r:
        r.migrate()
        yield r


@pytest.fixture(autouse=True)
def no_floors(monkeypatch):
    monkeypatch.delenv("MOFSBU_MIN_FREE_RAM_GIB", raising=False)
    monkeypatch.delenv("MOFSBU_MIN_FREE_VRAM_GIB", raising=False)


def spec(**over) -> BuildSpec:
    base = dict(
        molecules=(MoleculeSpec("catechol", "Oc1ccccc1O", 1),),
        metals=(MetalSpec("Zn", 2, "ls"),), degree=1, coordination=(4,),
        binding=("chelate",), pocket=PocketPredicate(ring_size=5, n_anionic=2),
        ligands_per_metal=(1,), seed=7)
    base.update(over)
    return BuildSpec(**base)


def test_no_floor_declared_reads_nothing(monkeypatch):
    def boom(*_a):
        raise AssertionError("read memory with no floor declared")

    monkeypatch.setattr(headroom, "free_ram_gib", boom)
    monkeypatch.setattr(headroom, "free_vram_gib", boom)
    assert headroom.shortfall(gpu=True) is None


@pytest.mark.parametrize("value", ["lots", "-1", "nan"])
def test_a_malformed_floor_is_refused_where_it_is_read(monkeypatch, value):
    monkeypatch.setenv("MOFSBU_MIN_FREE_RAM_GIB", value)
    with pytest.raises(ValueError, match="MOFSBU_MIN_FREE_RAM_GIB"):
        config.min_free_ram_gib()


def test_ram_under_the_floor_is_a_reason(monkeypatch):
    monkeypatch.setenv("MOFSBU_MIN_FREE_RAM_GIB", "8")
    monkeypatch.setattr(headroom, "free_ram_gib", lambda: 3.0)
    assert "RAM" in headroom.shortfall(gpu=False)
    monkeypatch.setattr(headroom, "free_ram_gib", lambda: 12.0)
    assert headroom.shortfall(gpu=False) is None


def test_the_vram_floor_applies_only_to_a_worker_feeding_the_device(monkeypatch):
    monkeypatch.setenv("MOFSBU_MIN_FREE_VRAM_GIB", "6")
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:1")
    monkeypatch.setattr(headroom, "free_vram_gib", lambda device: 2.0)
    assert "cuda:1" in headroom.shortfall(gpu=True)
    assert headroom.shortfall(gpu=False) is None


def test_a_floor_that_cannot_be_measured_is_refused_not_skipped(monkeypatch):
    monkeypatch.setenv("MOFSBU_MIN_FREE_VRAM_GIB", "6")
    monkeypatch.setenv("MOFSBU_DEVICE", "cpu")
    with pytest.raises(MofsbuError, match="not a CUDA device"):
        headroom.shortfall(gpu=True)


def test_free_ram_is_read_on_linux():
    ram = headroom.free_ram_gib()
    if ram is None:
        pytest.skip("no /proc/meminfo here")
    assert ram > 0


def test_a_worker_short_of_memory_holds_the_task_then_runs_it(reg, monkeypatch):
    """Nothing is claimed while there is a shortfall, the heartbeat keeps the run live,
    and the task runs once memory frees, rather than being claimed and failing."""
    from mofsbu import runner

    run_id = create_run(reg, spec())
    tid = add_task(reg, run_id, "place", {})
    reg.conn.commit()
    reasons = iter(["3.0 GiB RAM free < MOFSBU_MIN_FREE_RAM_GIB=8"] * 3)
    seen = []

    def short(*, gpu):
        status = reg.conn.execute("SELECT status FROM tasks WHERE id=?", (tid,)).fetchone()[0]
        beat = reg.conn.execute("SELECT heartbeat_at FROM runs WHERE id=?",
                                (run_id,)).fetchone()[0]
        seen.append((status, beat))
        return next(reasons, None)

    monkeypatch.setattr(headroom, "shortfall", short)
    monkeypatch.setattr(runner, "execute", lambda reg, task, spec: runner.Outcome(None, None))
    assert runner.work(reg, spec(), run_id, poll=0.01) == 1
    assert [s for s, _ in seen[:4]] == ["pending"] * 4, "claimed during a shortfall"
    assert seen[1][1] is not None, "a held-back worker stopped stamping the run"
    row = reg.conn.execute("SELECT status, attempts FROM tasks WHERE id=?", (tid,)).fetchone()
    assert (row["status"], row["attempts"]) == ("done", 1)


def test_an_empty_queue_still_ends_the_worker_during_a_shortfall(reg, monkeypatch):
    from mofsbu import runner

    run_id = create_run(reg, spec())
    reg.conn.commit()
    monkeypatch.setattr(headroom, "shortfall", lambda *, gpu: "always short")
    assert runner.work(reg, spec(), run_id, poll=0.01) == 0
    assert task_counts(reg, run_id) == {}
