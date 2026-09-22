"""Specs, the task queue, and the two ground rules they exist to serve."""
from __future__ import annotations

import json
import os
import socket

import pytest

from mofsbu import config
from mofsbu.registry import BlobStore, Registry, find
from mofsbu.registry.jobs import (
    add_task, claim_task, complete_task, create_run, fail_task, finish_run,
    reset_stale_claims, task_counts,
)
from mofsbu.runner import plan, run
from mofsbu.spec import (
    SPEC_VERSION, BuildSpec, MetalSpec, MoleculeSpec, PocketPredicate)


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as r:
        r.migrate()
        yield r


def _dead_pid() -> int:
    """A pid that certainly is not running: one that ran, exited, and was reaped.

    Skips where the question cannot be asked.  `worker_alive` is POSIX-only on purpose —
    `os.kill` on Windows terminates rather than probes — so a test that depends on a
    definitive "that process is gone" has nothing to assert there.
    """
    import subprocess
    import sys

    if os.name != "posix":
        pytest.skip("process liveness is POSIX-only; the heartbeat covers the rest")
    proc = subprocess.Popen([sys.executable, "-c", ""])
    proc.wait()
    return proc.pid


def catechol_spec(**over) -> BuildSpec:
    base = dict(
        molecules=(MoleculeSpec("catechol", "Oc1ccccc1O", 1),),
        metals=(MetalSpec("Zn", 2, "ls"),), degree=1, coordination=(4,),
        binding=("chelate",), pocket=PocketPredicate(ring_size=5, n_anionic=2),
        ligands_per_metal=(1,), seed=7)
    base.update(over)
    return BuildSpec(**base)


# ── a build is data ──────────────────────────────────────────────────────────

def test_spec_round_trips_and_hashes(tmp_path):
    spec = catechol_spec(note="round trip")
    path = spec.save(tmp_path / "spec.json")
    assert BuildSpec.load(path) == spec
    assert BuildSpec.load(path).digest == spec.digest


def test_identical_specs_share_a_digest_and_different_ones_do_not():
    assert catechol_spec().digest == catechol_spec().digest
    assert catechol_spec().digest != catechol_spec(coordination=(6,)).digest


# ── the queue ────────────────────────────────────────────────────────────────

def test_a_task_is_claimed_exactly_once(reg):
    """Two workers pulling at once must not get the same row."""
    run_id = create_run(reg, catechol_spec())
    for i in range(6):
        add_task(reg, run_id, "place", {"i": i})
    claimed = []
    while (task := claim_task(reg, run_id)) is not None:
        claimed.append(task.id)
        complete_task(reg, task.id)
    assert len(claimed) == len(set(claimed)) == 6
    assert claim_task(reg, run_id) is None


def test_rejected_is_not_failed(reg):
    """A candidate that fails QC is an answer, not a crash.

    Conflating them makes 'the run produced nothing' and 'the run broke' look identical.
    """
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    add_task(reg, run_id, "place", {})
    a, b = claim_task(reg, run_id), claim_task(reg, run_id)
    fail_task(reg, a.id, "QC said no", rejected=True)
    complete_task(reg, b.id)
    assert task_counts(reg, run_id) == {"rejected": 1, "done": 1}
    assert finish_run(reg, run_id) == "done"          # rejections do not fail the run


def test_a_dead_worker_returns_its_task_to_the_queue(reg):
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    claim_task(reg, run_id)                            # claimed, then the worker dies
    assert claim_task(reg, run_id) is None
    assert reset_stale_claims(reg, run_id) == 1
    assert claim_task(reg, run_id) is not None


# ── a run whose process is gone ──────────────────────────────────────────────
# A run row is written BY a process, so the one state it cannot record is "the process
# stopped existing".  Killed with its shell, the row keeps saying `pending` and the page
# keeps calling it ongoing — forever, and with its claimed tasks stranded where no worker
# will take them.

def test_a_run_whose_worker_is_gone_is_interrupted_not_ongoing(reg):
    from mofsbu.registry.jobs import INTERRUPTED, run_liveness, sweep_interrupted

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {"i": 0})
    add_task(reg, run_id, "place", {"i": 1})
    # A worker on this host claims one task and is killed: the claim stays behind with
    # its pid on it, and that pid is the evidence.
    claim_task(reg, run_id, worker=f"{socket.gethostname()}:{_dead_pid()}")
    reg.conn.commit()

    assert run_liveness(reg, run_id)["verdict"] == "interrupted"
    swept = sweep_interrupted(reg)
    assert [s["run_id"] for s in swept] == [run_id]
    assert swept[0]["returned_claims"] == 1

    row = reg.conn.execute("SELECT status FROM runs WHERE id=?", (run_id,)).fetchone()
    assert row["status"] == INTERRUPTED
    # The stranded task is claimable again — it was never in flight, only held.
    assert task_counts(reg, run_id).get("claimed") is None
    assert claim_task(reg, run_id) is not None


def test_a_live_worker_is_never_swept(reg):
    """The dangerous direction.  Sweeping a run that IS being worked on hands its tasks
    to a second worker and has the work done twice."""
    from mofsbu.registry.jobs import run_liveness, sweep_interrupted

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    claim_task(reg, run_id)                   # claimed by THIS process, which is alive
    reg.conn.commit()

    assert run_liveness(reg, run_id)["verdict"] == "live"
    assert sweep_interrupted(reg) == []
    assert task_counts(reg, run_id).get("claimed") == 1


def test_a_worker_on_another_host_is_unknown_not_dead(reg):
    """`claimed_by` is host:pid, so the question is only answerable on that host.

    From anywhere else the honest answer is "cannot tell", and a sweep must treat it as
    live: the two-machine workflow has one registry reachable from both ends, and
    guessing `dead` there would reset tasks the workstation is running.
    """
    from mofsbu.registry.jobs import run_liveness, sweep_interrupted

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    claim_task(reg, run_id, worker=f"some-other-box:{_dead_pid()}")
    reg.conn.commit()

    verdict = run_liveness(reg, run_id)
    assert verdict["verdict"] == "unknown" and "some-other-box" in verdict["reason"]
    assert sweep_interrupted(reg) == []


def test_a_long_task_is_not_mistaken_for_a_dead_one(reg, monkeypatch):
    """The heartbeat goes stale during an hour of xTB, and that must decide nothing.

    A claimed task is only ever judged by its process, never by the clock — otherwise the
    slowest tasks in the queue, which are the expensive ones, are exactly the ones that
    get reset and recomputed.
    """
    from mofsbu.registry import jobs
    from mofsbu.registry.jobs import run_liveness

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "relax", {})
    claim_task(reg, run_id)                   # this process: alive, and busy
    reg.conn.execute("UPDATE runs SET heartbeat_at=? WHERE id=?",
                     ("2020-01-01T00:00:00+00:00", run_id))
    reg.conn.commit()
    monkeypatch.setattr(jobs, "STALE_AFTER", 1.0)

    assert run_liveness(reg, run_id)["verdict"] == "live"
    assert jobs.sweep_interrupted(reg) == []


def test_a_run_that_never_claimed_anything_is_judged_by_its_heartbeat(reg, monkeypatch):
    """Killed between planning and the first claim: no pid anywhere to ask about.

    Nothing is claimed, so nothing can be disturbed by the verdict — which is exactly
    when the clock is safe to use.
    """
    from mofsbu.registry import jobs
    from mofsbu.registry.jobs import INTERRUPTED, run_liveness

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    reg.conn.commit()
    assert run_liveness(reg, run_id)["verdict"] == "live"      # just planned

    reg.conn.execute("UPDATE runs SET heartbeat_at=? WHERE id=?",
                     ("2020-01-01T00:00:00+00:00", run_id))
    reg.conn.commit()
    monkeypatch.setattr(jobs, "STALE_AFTER", 1.0)
    assert run_liveness(reg, run_id)["verdict"] == "interrupted"
    assert jobs.sweep_interrupted(reg)[0]["now"] == INTERRUPTED


def test_a_run_that_finished_without_being_closed_out_is_done_not_interrupted(reg):
    """Killed after the last task: the work IS all there, only the closing write is not.

    Reporting that as interrupted would invite someone to re-run a complete run.
    """
    from mofsbu.registry.jobs import run_liveness, sweep_interrupted

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    complete_task(reg, claim_task(reg, run_id).id)
    reg.conn.commit()

    assert run_liveness(reg, run_id)["verdict"] == "unfinalised"
    assert sweep_interrupted(reg)[0]["now"] == "done"


def test_a_stop_that_was_asked_for_stays_cancelled_even_if_the_worker_dies(reg):
    """The worker dying on the way out does not turn a decision into an accident."""
    from mofsbu.registry.jobs import request_cancel, sweep_interrupted

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    claim_task(reg, run_id, worker=f"{socket.gethostname()}:{_dead_pid()}")
    request_cancel(reg, run_id)

    assert sweep_interrupted(reg)[0]["now"] == "cancelled"
    assert reg.conn.execute("SELECT status FROM runs WHERE id=?",
                            (run_id,)).fetchone()["status"] == "cancelled"


def test_an_interrupted_run_can_be_resumed(reg):
    """The state nobody chose is the one where continuing matters most."""
    from mofsbu.registry.jobs import resume_run, sweep_interrupted

    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {"i": 0})
    add_task(reg, run_id, "place", {"i": 1})
    claim_task(reg, run_id, worker=f"{socket.gethostname()}:{_dead_pid()}")
    reg.conn.commit()
    sweep_interrupted(reg)

    resume_run(reg, run_id)
    assert reg.conn.execute("SELECT status FROM runs WHERE id=?",
                            (run_id,)).fetchone()["status"] == "running"
    assert task_counts(reg, run_id).get("pending") == 2


def test_a_worker_stamps_the_run_it_is_working_on(reg):
    """Without a stamp per task, a worker part-way through a long queue is
    indistinguishable from one that died at the first."""
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    reg.conn.execute("UPDATE runs SET heartbeat_at=NULL WHERE id=?", (run_id,))
    claim_task(reg, run_id)
    assert reg.conn.execute("SELECT heartbeat_at FROM runs WHERE id=?",
                            (run_id,)).fetchone()["heartbeat_at"] is not None


# ── one queue, two kinds of hardware ─────────────────────────────────────────

def test_a_worker_claims_only_the_kinds_it_is_asked_for(reg):
    """A GPU worker and a CPU worker pull from one queue and must not take each other's
    work: sharing it is what put construction behind the optimiser."""
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {"i": 0})
    add_task(reg, run_id, "relax", {"i": 1})
    add_task(reg, run_id, "place", {"i": 2})

    relaxer = claim_task(reg, run_id, kinds=("relax",))
    assert relaxer.kind == "relax"
    assert claim_task(reg, run_id, kinds=("relax",)) is None      # only the one

    builder = claim_task(reg, run_id, exclude_kinds=("relax",))
    assert builder.kind == "place"
    assert claim_task(reg, run_id, exclude_kinds=("relax",)).kind == "place"
    assert claim_task(reg, run_id, exclude_kinds=("relax",)) is None


def test_an_empty_kind_list_claims_nothing(reg):
    """`kinds=()` is 'this worker takes nothing', which is not the same as no filter."""
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    assert claim_task(reg, run_id, kinds=()) is None
    assert claim_task(reg, run_id) is not None


# ── ground rule 8: the laptop is the default ─────────────────────────────────

def test_parallelism_is_opt_in(monkeypatch):
    monkeypatch.delenv("MOFSBU_PROFILE", raising=False)
    monkeypatch.delenv("MOFSBU_WORKERS", raising=False)
    assert config.machine_profile() == config.LAPTOP
    assert config.max_workers() == 1
    assert not config.parallel_enabled()

    monkeypatch.setenv("MOFSBU_WORKERS", "4")
    assert config.max_workers() == 4 and config.parallel_enabled()

    monkeypatch.delenv("MOFSBU_WORKERS")
    monkeypatch.setenv("MOFSBU_PROFILE", "workstation")
    assert config.machine_profile() == config.WORKSTATION
    assert config.max_workers() == max(1, (os.cpu_count() or 2) - 1)


def test_a_run_on_an_undeclared_machine_stays_in_process(reg, monkeypatch):
    monkeypatch.delenv("MOFSBU_PROFILE", raising=False)
    monkeypatch.delenv("MOFSBU_WORKERS", raising=False)
    import multiprocessing as mp

    def explode(*_a, **_k):                      # pragma: no cover - must never run
        raise AssertionError("spawned a process on an undeclared machine")

    # Both, because the pool is started from the spawn context rather than from the
    # module: a forked child cannot re-initialise CUDA, which is the one configuration
    # the split exists for.
    monkeypatch.setattr(mp, "Process", explode)
    monkeypatch.setattr(mp.get_context("spawn"), "Process", explode)
    summary = run(reg, catechol_spec())
    assert summary["workers"] == 1
    assert summary["counts"].get("done", 0) > 0


# ── the GPU is one device; the cores are not ─────────────────────────────────

def test_the_pool_dedicates_a_worker_to_the_relax_queue_on_a_gpu(monkeypatch):
    """The reported bug: a declared GPU relaxed while ONE core built for it.

    The declared total is divided, never exceeded — the workstation profile already
    leaves a core for the machine to stay responsive, and adding a GPU feeder on top
    would take it back.
    """
    from mofsbu.runner import plan_workers

    monkeypatch.delenv("MOFSBU_RELAX_WORKERS", raising=False)
    monkeypatch.setenv("MOFSBU_WORKERS", "8")

    monkeypatch.setenv("MOFSBU_DEVICE", "cpu")
    cpu = plan_workers()
    assert (cpu.total, cpu.build, cpu.relax) == (8, 8, 0), (
        "on CPU, relaxation is core work like everything else and splitting the pool "
        "would only reserve a core for it")

    monkeypatch.setenv("MOFSBU_DEVICE", "cuda")
    gpu = plan_workers()
    assert (gpu.total, gpu.build, gpu.relax) == (8, 7, 1)
    assert gpu.build > 1, "the whole point: construction is not down to one core"

    # One card stays one card however wide the declared pool is: a second process feeding
    # it divides its memory rather than multiplying its throughput.
    monkeypatch.setenv("MOFSBU_WORKERS", "64")
    wide = plan_workers()
    assert (wide.total, wide.build, wide.relax) == (64, 63, 1)

    monkeypatch.setenv("MOFSBU_WORKERS", "8")
    monkeypatch.setenv("MOFSBU_RELAX_WORKERS", "3")            # a multi-GPU box says so
    assert plan_workers().relax == 3
    monkeypatch.setenv("MOFSBU_RELAX_WORKERS", "0")            # or opts back out
    assert plan_workers().relax == 0


def test_a_construct_run_reserves_nothing_for_a_queue_it_will_not_have(monkeypatch):
    """`construct` queues no relaxations, so a worker held back for them idles all run."""
    from mofsbu.runner import plan_workers

    monkeypatch.setenv("MOFSBU_DEVICE", "cuda")
    monkeypatch.setenv("MOFSBU_WORKERS", "8")
    assert plan_workers(relaxes=False).build == 8
    assert plan_workers(relaxes=False).relax == 0


def test_an_undrained_queue_is_not_a_finished_run(reg):
    """A worker that dies records nothing, so its tasks are still pending — and a run
    summarised from only the survivors' successes must not read as `done`."""
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {"i": 0})
    add_task(reg, run_id, "place", {"i": 1})
    complete_task(reg, claim_task(reg, run_id).id)
    assert finish_run(reg, run_id) == "failed"

    complete_task(reg, claim_task(reg, run_id).id)
    assert finish_run(reg, run_id) == "done"


def test_the_relax_share_never_consumes_the_builders(monkeypatch):
    """All-relax would starve the queue that feeds it — one busy core again."""
    from mofsbu.runner import plan_workers

    monkeypatch.setenv("MOFSBU_DEVICE", "cuda")
    monkeypatch.setenv("MOFSBU_WORKERS", "2")
    monkeypatch.setenv("MOFSBU_RELAX_WORKERS", "9")
    pool = plan_workers()
    assert (pool.total, pool.build, pool.relax) == (2, 1, 1)


def test_an_explicit_worker_count_is_a_declaration(monkeypatch):
    """`--workers 4` on an unconfigured machine means four, not 'laptop, so one'.

    Ground rule 8 is that parallelism is never ASSUMED; a number typed on the command
    line or in the page is somebody assuming it on purpose.
    """
    from mofsbu.runner import plan_workers

    monkeypatch.delenv("MOFSBU_PROFILE", raising=False)
    monkeypatch.delenv("MOFSBU_WORKERS", raising=False)
    assert plan_workers().total == 1                  # nothing declared: still one
    assert plan_workers(4).total == 4


def test_a_relax_worker_waits_for_work_the_builders_have_not_queued_yet(reg, monkeypatch):
    """An empty relax queue early in a run means 'not yet', not 'never'.

    Relax tasks are queued BY the build tasks as they finish, so a relax worker that
    stopped at the first empty claim would exit seconds into a run and leave the
    accelerator idle for the rest of it.
    """
    import threading

    from mofsbu import runner

    monkeypatch.setattr(runner, "execute",
                        lambda reg, task, spec: runner.Outcome(None, None))
    run_id = create_run(reg, catechol_spec())
    reg.conn.commit()

    builders_running = threading.Event()
    builders_running.set()

    def queue_one_late() -> None:
        with Registry(reg.db_path, BlobStore(reg.store.root)) as other:
            add_task(other, run_id, "relax", {"late": True})
            other.conn.commit()
        builders_running.clear()

    threading.Timer(0.3, queue_one_late).start()
    done = runner.work(reg, catechol_spec(), run_id, kinds=("relax",),
                       wait_while=builders_running.is_set, poll=0.05)
    assert done == 1, "the worker gave up before the builders had queued anything"
    assert task_counts(reg, run_id).get("done") == 1


def test_a_declared_pool_really_drains_the_queue_in_parallel(reg, monkeypatch):
    """End to end through spawned processes: every task is executed exactly once.

    Cheap to run and worth the seconds — the two things that only break across a process
    boundary are the claim (two workers, one row) and the writes (several processes, one
    SQLite file), and neither is exercised by an in-process worker.
    """
    monkeypatch.setenv("MOFSBU_WORKERS", "3")
    monkeypatch.setenv("MOFSBU_PROFILE", "workstation")
    summary = run(reg, catechol_spec(coordination=(4, 6), run_mode="construct"))
    assert summary["workers"] == 3
    counts = summary["counts"]
    assert counts.get("pending", 0) == 0 and counts.get("claimed", 0) == 0
    assert counts.get("done", 0) > 0
    claimants = {r[0] for r in reg.conn.execute(
        "SELECT DISTINCT claimed_by FROM tasks WHERE run_id=?", (summary["run_id"],))}
    # Every task was executed by a worker process, not by the process that planned the
    # run.  Not "by more than one of them": spawn takes about a second to come up and a
    # queue this short can be finished by whichever worker is ready first, which is
    # parallelism working rather than parallelism failing.
    assert claimants and f":{os.getpid()}" not in "".join(claimants)
    # Attempts stay at one per task: a second claim of the same row is the failure this
    # queue's IMMEDIATE transaction exists to prevent.
    assert summary["summary"]["retried_tasks"] == 0


# ── ground rule 7: missing functions are stubs that shout ────────────────────

def test_degree_above_one_refuses_rather_than_doing_something_smaller(reg):
    from mofsbu.assembly.join import NotBuiltYet

    with pytest.raises(NotBuiltYet, match="M5"):
        plan(reg, catechol_spec(degree=2))


def test_the_assembly_interfaces_that_are_still_stubs_all_raise():
    """Signatures are settled so callers can be written; bodies are scheduled work.

    The list shortened as M5 landed — `compatible` left at S2, `join` at S3, `grow` at S4,
    each covered by its own file (`test_compatible.py`, `test_join.py`,
    `test_construct.py`). `place_multicentre` is what is left, and the point of the test is
    unchanged: a scheduled body raises rather than returning a plausible-looking value.
    """
    from mofsbu._types import NotBuiltYet
    from mofsbu.geometry.placer import place_multicentre

    with pytest.raises(NotBuiltYet):
        place_multicentre([1, 2], [], [])


# ── end to end ───────────────────────────────────────────────────────────────

def test_a_spec_produces_structures(reg):
    summary = run(reg, catechol_spec(ligands_per_metal=(1, 2)))
    assert summary["status"] == "done"
    labels = [r["display_label"] for r in find(reg, limit=100)]
    assert any(lbl.startswith("Zn[") for lbl in labels)
    assert any("catechol" in lbl for lbl in labels)     # fragment alias applied on relabel


# ── a run that builds nothing must say why ───────────────────────────────────

def anthrarufin_spec(**over) -> BuildSpec:
    """The user's actual failing case: a bidentate ligand, no co-ligand, CN 4 and 6."""
    base = dict(
        molecules=(MoleculeSpec("ATF", "C1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C(=CC=C3)O", 1),),
        metals=(MetalSpec("Cu", 2, "ls"),), degree=1, coordination=(4, 6),
        binding=("chelate", "mono"), co_ligand=None, ligands_per_metal=(1,))
    base.update(over)
    return BuildSpec(**base)


def test_a_run_that_queues_nothing_explains_itself(reg):
    """Silence was the bug: a chelate needs 2 more sites at CN 4 and there is no
    co-ligand, so every candidate was skipped and the run looked like a success that
    simply produced no complexes."""
    from mofsbu.registry.jobs import get_diagnostics
    from mofsbu.runner import plan

    run_id, n_tasks = plan(reg, anthrarufin_spec())
    diagnostics = get_diagnostics(reg, run_id)
    assert diagnostics, "a run that builds no complexes must record why"
    assert any("co-ligand" in d["reason"] for d in diagnostics)
    assert all(d["hint"] for d in diagnostics)
    # ligands are still built; it is the coordination step that had nothing to do
    kinds = {r["kind"] for r in reg.conn.execute(
        "SELECT DISTINCT kind FROM tasks WHERE run_id=?", (run_id,))}
    assert kinds == {"ligand"}


def test_setting_a_co_ligand_builds_the_complexes(reg):
    from mofsbu.registry import find

    summary = run(reg, anthrarufin_spec(co_ligand="O"))
    assert summary["counts"].get("done", 0) > 0
    assert find(reg, metal="Cu", limit=50), "Cu complexes should now exist"


def test_the_placer_refusing_is_a_rejection_not_a_failure(reg):
    """A bidentate ligand cannot span a linear two-coordinate centre.

    Saying so is correct behaviour; recording it as `failed` would make a run full of
    sound chemistry look like a run full of bugs.
    """
    summary = run(reg, anthrarufin_spec(allow_unsaturated=True))
    assert summary["counts"].get("failed", 0) == 0
    assert summary["counts"].get("rejected", 0) > 0


# ── specs are files, so the schema needs migrations ──────────────────────────

def test_a_version_1_spec_still_loads(tmp_path):
    """Changing the dataclass without a migration silently breaks every spec on disk."""
    import json

    v1 = {"spec_version": 1, "molecules": [{"name": "x", "smiles": "O", "multiplicity": 1,
                                            "max_deprotonations": None}],
          "metals": [], "degree": 1, "coordination": [4], "geometries": None,
          "ligands_per_metal": [1], "co_ligand": "O", "binding": ["chelate"],
          "pocket": {"ring_size": None, "n_anionic": None, "min_anionic": None,
                     "rigid": None},
          "seed": 0, "relax_to": None, "note": "old"}
    path = tmp_path / "old.json"
    path.write_text(json.dumps(v1))
    spec = BuildSpec.load(path)
    assert spec.run_mode == "construct" and spec.allow_unsaturated is False
    # v4 added `ml_model`.  A v1 spec never expressed a choice of ML potential, so it
    # comes forward as None ("use what this machine declares") rather than being
    # back-filled with "mace-mp-0", which would invent a decision nobody made.
    assert spec.ml_model is None


def test_a_version_3_spec_comes_forward_without_inventing_a_model(tmp_path):
    spec = BuildSpec.from_dict({
        "spec_version": 3, "run_mode": "ml_go",
        "molecules": [{"name": "x", "smiles": "O", "multiplicity": 1,
                       "max_deprotonations": None}]})
    assert spec.ml_model is None
    assert spec.to_dict()["spec_version"] == SPEC_VERSION
    # v6 adds mixed-ligand enumeration.  A v3 spec meant one molecule per coordination
    # sphere, so it comes forward homoleptic — defaulting it to anything else would
    # silently multiply the size of every run already on disk.
    assert spec.max_distinct_ligands == 1
    # v7 adds the pathway ladder, and the same rule applies: an old spec asked for the
    # products, not for the intermediates they are assembled from.
    assert spec.pathways is False


def test_a_version_4_spec_drops_the_stale_metal_multiplicity():
    """v5 removes `MetalSpec.multiplicity` — it was never kept in sync with
    `oxidation_state`/`spin_class`, which is exactly how a Cu(II) centre ended up
    stored as a singlet.  A v4 spec's stored value is dropped, not trusted."""
    spec = BuildSpec.from_dict({
        "spec_version": 4,
        "molecules": [{"name": "x", "smiles": "O", "multiplicity": 1,
                       "max_deprotonations": None}],
        "metals": [{"symbol": "Cu", "oxidation_state": 2, "spin_class": "ls",
                    "multiplicity": 1}]})
    assert not hasattr(spec.metals[0], "multiplicity")
    assert spec.metals[0].spin_class == "ls"


def test_a_spec_may_name_the_ml_potential_and_a_typo_is_refused():
    """`ml_go` is a rung, not a theory; the spec is where the theory is recorded."""
    mol = [{"name": "x", "smiles": "O", "multiplicity": 1, "max_deprotonations": None}]
    spec = BuildSpec.from_dict({"spec_version": 4, "molecules": mol,
                                "run_mode": "ml_go", "ml_model": "MACE-OMOL-0"})
    assert spec.ml_model == "MACE-OMOL-0"
    # Two specs differing only in the model are two different runs, not one.
    other = BuildSpec.from_dict({"spec_version": 4, "molecules": mol,
                                 "run_mode": "ml_go", "ml_model": "mace-mp-0"})
    assert spec.digest != other.digest
    with pytest.raises(ValueError, match="unknown ML model"):
        BuildSpec.from_dict({"spec_version": 4, "molecules": mol,
                             "run_mode": "ml_go", "ml_model": "mace-omol-1"})


def test_a_future_spec_version_is_refused_not_guessed():
    with pytest.raises(ValueError, match="newer than this build"):
        BuildSpec.from_dict({"spec_version": 999, "molecules": [
            {"name": "x", "smiles": "O", "multiplicity": 1, "max_deprotonations": None}]})


def test_an_unknown_field_is_refused():
    with pytest.raises(ValueError, match="unknown spec field"):
        BuildSpec.from_dict({"spec_version": 2, "molecules": [
            {"name": "x", "smiles": "O", "multiplicity": 1, "max_deprotonations": None}],
            "wishful_thinking": True})


def test_non_spec_json_is_not_mistaken_for_a_spec():
    assert not BuildSpec.looks_like_spec({"molecules": [{"name": "a", "smiles": "O"}]})
    assert BuildSpec.looks_like_spec({"spec_version": 2, "molecules": []})


# ── heteroleptic enumeration: more than one ligand kind on a centre ──────────

def _place_compositions(reg, spec):
    from mofsbu.runner import plan
    import json as _json

    run_id, _ = plan(reg, spec)
    out = []
    for r in reg.conn.execute(
            "SELECT payload_json FROM tasks WHERE run_id=? AND kind='place'", (run_id,)):
        p = _json.loads(r["payload_json"])
        out.append(tuple(sorted(
            (spec.molecules[c["molecule"]].name, bool(c["selection"]), c["count"])
            for c in p["components"])))
    return set(out)


def _two_ligand_spec(**over):
    base = dict(
        molecules=(MoleculeSpec("A", "CC(=O)C", 1), MoleculeSpec("B", "Cl", 1)),
        metals=(MetalSpec("Mg", 2, "ls"),), coordination=(4,),
        ligands_per_metal=(2,), binding=("mono",), co_ligand=None,
        allow_unsaturated=True, geometries=("tetrahedral",))
    base.update(over)
    return BuildSpec(**base)


def test_one_distinct_ligand_is_homoleptic_and_is_the_default(reg):
    """The default must reproduce what every spec written before this field meant.

    Two molecules in a spec used to be enumerated independently and could never share a
    coordination sphere; that is still what happens unless you ask for otherwise, so no
    stored run changes size or content underneath you.
    """
    spec = _two_ligand_spec()
    assert spec.max_distinct_ligands == 1
    combos = _place_compositions(reg, spec)
    assert all(len(c) == 1 for c in combos), f"a default spec mixed ligands: {combos}"


def test_raising_the_cap_mixes_two_molecules_on_one_centre(reg):
    combos = _place_compositions(reg, _two_ligand_spec(max_distinct_ligands=2))
    mixed = {c for c in combos if len({name for name, _, _ in c}) > 1}
    assert mixed, "asked for 2 distinct ligand kinds and got only homoleptic spheres"
    assert any(len(c) == 1 for c in combos), "the cap is a maximum, not a target"


def test_two_protomers_of_one_molecule_are_two_kinds(reg):
    """HCl and Cl- may share a metal: a partially deprotonated set is real chemistry
    for a polyprotic linker, and this is the decision that allows it."""
    combos = _place_compositions(reg, _two_ligand_spec(max_distinct_ligands=2))
    same_mol_mixed = {
        c for c in combos
        if len(c) == 2 and len({name for name, _, _ in c}) == 1
        and len({deprot for _, deprot, _ in c}) == 2}
    assert same_mol_mixed, "no complex mixed a protonated and deprotonated form"


def test_a_composition_never_exceeds_the_requested_ligand_count(reg):
    for combo in _place_compositions(reg, _two_ligand_spec(max_distinct_ligands=2)):
        assert sum(count for _, _, count in combo) == 2


def test_max_distinct_ligands_must_be_at_least_one():
    with pytest.raises(ValueError, match="at least 1"):
        _two_ligand_spec(max_distinct_ligands=0)


# ── counts written as ranges ─────────────────────────────────────────────────
# A sweep of one to three ligand copies is ONE experiment, and typing it out as a list is
# the notation getting in the way of the question. The expansion happens in `from_dict`,
# so the page, a hand-written file and the estimate all read it the same way.

def test_a_range_expands_and_the_spec_stores_the_integers():
    from mofsbu.spec import int_series

    assert int_series("1~3") == (1, 2, 3)
    assert int_series("4-6, 8") == (4, 5, 6, 8)
    assert int_series("1..2") == (1, 2)
    assert int_series([4, 6]) == (4, 6)
    assert int_series("2,2,3") == (2, 3)               # duplicates collapse
    spec = catechol_spec(coordination="4~6", ligands_per_metal="1~2")
    assert spec.coordination == (4, 5, 6) and spec.ligands_per_metal == (1, 2)
    # The reproducibility record carries what was enumerated, not how it was phrased.
    assert json.loads(spec.to_json())["coordination"] == [4, 5, 6]


def test_a_range_written_as_a_spec_file_round_trips(tmp_path):
    spec = BuildSpec.from_dict({**catechol_spec().to_dict(), "ligands_per_metal": "1~3"})
    assert BuildSpec.load(spec.save(tmp_path / "range.json")).ligands_per_metal == (1, 2, 3)


@pytest.mark.parametrize("bad, says", [
    ("3~1", "counts down"),
    ("1~500", "past the"),
    ("0", "start at 1"),
    ("four", "whole number"),
    ("", "empty"),
])
def test_a_range_that_cannot_mean_what_it_says_refuses(bad, says):
    """Refusing beats truncating: a typo that would queue 500 coordination numbers is
    the case the ceiling exists for, and it says so rather than quietly building some."""
    with pytest.raises(ValueError, match=says):
        catechol_spec(coordination=bad)
