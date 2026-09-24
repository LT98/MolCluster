"""Runs and tasks: the queue that separates asking for work from doing it.

The claim is the only interesting part.  `claim_task` takes an IMMEDIATE transaction and
flips exactly one pending row to `claimed`, so two workers pulling at the same moment
cannot get the same task.  That is what lets the same code run as one in-process worker
on a laptop and as N processes on the workstation, and later under a launcher that knows
about MPI ranks or GPUs, without any caller changing.

SQLite is honest for a workstation with local disk.  On a shared cluster filesystem its
locking is not dependable, so nothing here exposes SQLite in its signatures — moving the
queue to Postgres should not touch a single caller.
"""
from __future__ import annotations

import json
import os
import socket
from dataclasses import dataclass
from typing import Any, Iterator, Sequence

from datetime import datetime, timezone

from mofsbu.registry.db import Registry, utcnow
from mofsbu.spec import BuildSpec

PENDING, CLAIMED, DONE, FAILED, REJECTED = "pending", "claimed", "done", "failed", "rejected"
#: How many times one claim looks again after another worker took the row it chose.
CLAIM_ATTEMPTS = 8
# A build you stopped on purpose is not a build that failed.  Keeping them apart is the
# same rule as `rejected` vs `failed`: a run list where every abandoned experiment reads
# as a crash is a run list nobody trusts.
CANCELLING, CANCELLED = "cancelling", "cancelled"
# And a build that was killed is neither.  Nothing asked it to stop and nothing went
# wrong with the chemistry: the process it was running in stopped existing, which is a
# third thing and the only one that can simply be continued.
INTERRUPTED = "interrupted"

#: How long a run with nothing in flight may go unstamped before it is taken to be over.
#: Only ever applied when NO task is claimed — see `run_liveness`, which refuses to call a
#: run dead while a worker might be inside a long relaxation.
STALE_AFTER = 300.0


def worker_id() -> str:
    return f"{socket.gethostname()}:{os.getpid()}"


def worker_alive(claimed_by: str | None) -> bool | None:
    """Is the process that claimed this task still there?  None = cannot tell from here.

    `claimed_by` is `host:pid`, so the question is only answerable on the host that wrote
    it — from anywhere else the honest answer is "unknown", not "no".  Three-valued for
    exactly that reason: a `False` here is used as PROOF that a run is over, and a guess
    dressed as proof would hand another worker tasks that are still being worked on.

    POSIX only.  `os.kill(pid, 0)` is the standard liveness probe there; on Windows
    `os.kill` maps onto `TerminateProcess`, so asking the question would answer it, and
    this returns None instead and lets the heartbeat decide.
    """
    if not claimed_by or ":" not in claimed_by:
        return None
    host, _, pid_text = claimed_by.rpartition(":")
    if host != socket.gethostname() or os.name != "posix":
        return None
    try:
        pid = int(pid_text)
    except ValueError:
        return None
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True            # alive; it just is not ours to signal
    except OSError:
        return None
    return True


def _age_seconds(stamp: str | None) -> float | None:
    """Seconds since an ISO timestamp, or None if there is nothing to measure."""
    if not stamp:
        return None
    try:
        when = datetime.fromisoformat(stamp)
    except ValueError:
        return None
    if when.tzinfo is None:
        when = when.replace(tzinfo=timezone.utc)
    return (datetime.now(timezone.utc) - when).total_seconds()


def touch_run(reg: Registry, run_id: int) -> None:
    """Record that a worker is still here.  One indexed UPDATE, on a row already open."""
    reg.conn.execute("UPDATE runs SET heartbeat_at=? WHERE id=?", (utcnow(), run_id))


@dataclass(frozen=True)
class Task:
    id: int
    run_id: int
    kind: str
    payload: dict[str, Any]
    attempts: int


def create_run(reg: Registry, spec: BuildSpec, *, note: str = "",
               device: str | None = None, workers: int | None = None) -> int:
    """Open a run row.

    `device` and `workers` are recorded, not applied: what the machine will actually use
    is read from the environment at execution time (`config.compute_device`,
    `config.max_workers`), and this is the record of what that was when the run was
    submitted.  Left to the caller because the caller is the only one that knows whether
    a choice was expressed — `None` means "not stated", which is different from "cpu".
    """
    if device is None:
        from mofsbu.config import compute_device

        try:
            device = compute_device()
        except ValueError:
            device = ""                      # a malformed declaration is not a run error
    if workers is None:
        from mofsbu.config import max_workers

        workers = max_workers()
    cur = reg.conn.execute(
        "INSERT INTO runs (spec_digest, spec_json, status, note, host, device, workers,"
        " created_at) VALUES (?,?,?,?,?,?,?,?)",
        (spec.digest, spec.to_json(indent=None), PENDING, note or spec.note,
         socket.gethostname(), device, workers, utcnow()))
    return int(cur.lastrowid)


def set_diagnostics(reg: Registry, run_id: int, diagnostics: list[dict[str, Any]]) -> None:
    """Record why candidates were not queued.  Read back by the UI and the CLI."""
    reg.conn.execute("UPDATE runs SET diagnostics_json=? WHERE id=?",
                     (json.dumps(diagnostics), run_id))


def get_diagnostics(reg: Registry, run_id: int) -> list[dict[str, Any]]:
    row = reg.conn.execute("SELECT diagnostics_json FROM runs WHERE id=?", (run_id,)).fetchone()
    return json.loads(row["diagnostics_json"]) if row else []


def add_task(reg: Registry, run_id: int, kind: str, payload: dict[str, Any],
             *, priority: int = 0) -> int:
    cur = reg.conn.execute(
        "INSERT INTO tasks (run_id, kind, payload_json, status, priority, created_at) "
        "VALUES (?,?,?,?,?,?)",
        (run_id, kind, json.dumps(payload, sort_keys=True), PENDING, priority, utcnow()))
    return int(cur.lastrowid)


def claim_task(reg: Registry, run_id: int | None = None, *,
               worker: str | None = None,
               kinds: Sequence[str] | None = None,
               exclude_kinds: Sequence[str] | None = None) -> Task | None:
    """Atomically take one pending task, or None if there are none left.

    The row is taken by a compare-and-set (`UPDATE … WHERE status='pending'`) under
    IMMEDIATE, so two workers that chose the same row cannot both have it: the loser sees
    zero rows changed and looks again.  Without that, two workers claim the same row and
    do the same job twice.

    `kinds` / `exclude_kinds` narrow what this worker will take, which is what makes a
    heterogeneous pool possible: the tasks a GPU should run and the tasks that should
    keep the cores busy are different rows in one queue, and a worker that claims either
    puts them back in competition for the same machine.  Neither filter changes the
    claim's atomicity — they are predicates on the same single-row UPDATE.
    """
    if kinds is not None and not kinds:
        # An empty allow-list is "this worker takes nothing", not "take anything" — and
        # `kind IN ()` is not SQL.
        return None
    if run_id is not None and cancel_requested(reg, run_id):
        # Defence in depth.  `work()` checks before each claim, but a separate worker
        # PROCESS can be mid-loop when the stop arrives, and the claim is the one
        # chokepoint every worker passes through.
        return None

    who = worker or worker_id()
    conn = reg.conn
    where = ["status = ?"]
    args: list[Any] = [PENDING]
    if run_id:
        where.append("run_id = ?")
        args.append(run_id)
    if kinds is not None:
        where.append(f"kind IN ({','.join('?' * len(kinds))})")
        args += list(kinds)
    if exclude_kinds:
        where.append(f"kind NOT IN ({','.join('?' * len(exclude_kinds))})")
        args += list(exclude_kinds)
    # A step that joins onto another task's product waits for it (B25).  Only `grow`
    # names a parent, so only a `grow` row pays for reading its payload.
    where.append(
        "(kind != 'grow' OR NOT EXISTS (SELECT 1 FROM tasks dep WHERE dep.id IN ("
        "json_extract(tasks.payload_json, '$.parent_task'), "
        "json_extract(tasks.payload_json, '$.ligand_task')) AND dep.status IN (?, ?)))")
    args += [PENDING, CLAIMED]
    sql = ("SELECT id, run_id, kind, payload_json, attempts FROM tasks WHERE "
           + " AND ".join(where) + " ORDER BY priority DESC, id ASC LIMIT 1")

    # The search is a plain read — WAL lets it run beside a writer — and the write lock is
    # held only for the compare-and-set that takes the row.  Searching inside the lock held
    # it for up to ~11 ms per claim on a 13k-task queue, and 32 workers polling that way
    # starved every real write past BUSY_TIMEOUT.  A row another worker took first is
    # simply looked for again.
    #
    # sqlite3 opens an implicit transaction for DML, so an explicit BEGIN IMMEDIATE
    # collides with it.  The connection is switched to manual control for the duration of
    # the claim only: making it global would turn every multi-statement write (put_structure
    # inserts a structure, its metals and its fragments) into separate autocommits and
    # lose their atomicity.
    conn.commit()
    previous = conn.isolation_level
    conn.isolation_level = None
    try:
        for _ in range(CLAIM_ATTEMPTS):
            row = conn.execute(sql, tuple(args)).fetchone()
            if row is None:
                return None
            conn.execute("BEGIN IMMEDIATE")
            try:
                now = utcnow()
                taken = conn.execute(
                    "UPDATE tasks SET status=?, claimed_by=?, claimed_at=?, "
                    "attempts=attempts+1 WHERE id=? AND status=?",
                    (CLAIMED, who, now, row["id"], PENDING)).rowcount
                if taken:
                    # In the same transaction as the claim, because it is the same fact: a
                    # worker is here.  `run_liveness` reads it to tell a run nobody is
                    # working on from one whose worker has nothing to report yet.
                    conn.execute("UPDATE runs SET heartbeat_at=? WHERE id=?",
                                 (now, row["run_id"]))
                conn.execute("COMMIT")
            except Exception:
                conn.execute("ROLLBACK")
                raise
            if taken:
                return Task(int(row["id"]), int(row["run_id"]), row["kind"],
                            json.loads(row["payload_json"]), int(row["attempts"]) + 1)
        return None
    finally:
        conn.isolation_level = previous


def release_task(reg: Registry, task_id: int) -> None:
    """Put a claimed task back in the queue, as if it had never been taken.

    For a task that did not fail on its own account — the database was busy — so it is
    retried by whichever worker claims it next rather than recorded as a failure.
    `attempts` keeps counting, so a task that keeps being returned is visible.
    """
    reg.conn.execute(
        "UPDATE tasks SET status=?, claimed_by=NULL, claimed_at=NULL WHERE id=? AND status=?",
        (PENDING, task_id, CLAIMED))
    reg.conn.commit()


def waiting_on_live_work(reg: Registry, run_id: int, *,
                         kinds: Sequence[str] | None = None,
                         exclude_kinds: Sequence[str] | None = None) -> bool:
    """True when this worker's share of the queue is not empty, only blocked, and the
    work it waits for is held by a worker that is still alive.

    `claim_task` holds back a step until the task it joins onto has finished, so an empty
    claim no longer means "nothing left".  Waiting is right while someone is building the
    parent; if every claim in the run belongs to a dead worker, waiting would never end, so
    this says no and the worker exits as it did before.
    """
    where, args = ["run_id = ?", "status = ?"], [run_id, PENDING]
    if kinds is not None:
        where.append(f"kind IN ({','.join('?' * len(kinds))})")
        args += list(kinds)
    if exclude_kinds:
        where.append(f"kind NOT IN ({','.join('?' * len(exclude_kinds))})")
        args += list(exclude_kinds)
    if reg.conn.execute("SELECT 1 FROM tasks WHERE " + " AND ".join(where) + " LIMIT 1",
                        tuple(args)).fetchone() is None:
        return False
    return any(worker_alive(r["claimed_by"]) is not False for r in reg.conn.execute(
        "SELECT DISTINCT claimed_by FROM tasks WHERE run_id=? AND status=?",
        (run_id, CLAIMED)))


def complete_task(reg: Registry, task_id: int, *, structure_id: int | None = None,
                  geometry_id: int | None = None,
                  structure_created: bool | None = None,
                  geometry_created: bool | None = None,
                  detail: dict[str, Any] | None = None) -> None:
    """Finish a task, recording whether it WROTE anything.

    `structure_created=False` is the idempotency case (D2): the task ran, produced a
    structure the registry already had, and wrote no new row.  That is a success, and it
    used to be indistinguishable from a task that built something new — which is the
    difference between "the enumerator covered ground it had already covered" and "the
    enumerator produced 400 structures".
    """
    reg.conn.execute(
        "UPDATE tasks SET status=?, structure_id=?, geometry_id=?, error=NULL, "
        "error_code=NULL, structure_created=?, geometry_created=?, detail_json=?, "
        "finished_at=? WHERE id=?",
        (DONE, structure_id, geometry_id,
         None if structure_created is None else int(structure_created),
         None if geometry_created is None else int(geometry_created),
         json.dumps(detail) if detail else None, utcnow(), task_id))


def fail_task(reg: Registry, task_id: int, error: str, *, rejected: bool = False,
              code: str | None = None, detail: dict[str, Any] | None = None) -> None:
    """`rejected` means the task was well-formed and the answer was no.

    A candidate that fails QC is not a failure of the machinery, and conflating the two
    makes a run that produced nothing look identical to a run that crashed.

    `code` and `detail` are the second half of that distinction.  The message is for
    reading; the code is for GROUPING (200 rejections that are all one cause should say
    so in one line); the detail is for diagnosis without re-running anything — which
    atoms, which elements, which distance against which limit.
    """
    reg.conn.execute(
        "UPDATE tasks SET status=?, error=?, error_code=?, detail_json=?, "
        "finished_at=? WHERE id=?",
        (REJECTED if rejected else FAILED, error[:2000], code,
         json.dumps(detail) if detail else None, utcnow(), task_id))


def task_rows(reg: Registry, run_id: int, *, status: str | None = None,
              limit: int = 500, offset: int = 0) -> list[dict[str, Any]]:
    """Task records with their JSON columns decoded.  The inspector's data source."""
    sql = ("SELECT * FROM tasks WHERE run_id=?" + (" AND status=?" if status else "")
           + " ORDER BY id LIMIT ? OFFSET ?")
    args: tuple[Any, ...] = ((run_id, status, limit, offset) if status
                             else (run_id, limit, offset))
    out = []
    for row in reg.conn.execute(sql, args):
        item = dict(row)
        item["payload"] = json.loads(item.pop("payload_json") or "{}")
        item["detail"] = json.loads(item.pop("detail_json") or "null")
        out.append(item)
    return out


def outcome_summary(reg: Registry, run_id: int) -> dict[str, Any]:
    """A run at a glance: what happened, how often, and to what.

    Everything here is one GROUP BY over the task table — which is the point of putting
    the codes and flags in columns rather than in a message.  Reading a run no longer
    means scrolling a console that has already refreshed.
    """
    counts = task_counts(reg, run_id)
    by_code = [
        {"status": r["status"], "code": r["error_code"] or "(uncoded)", "count": r["n"],
         "example": r["example"], "example_task": r["example_task"]}
        for r in reg.conn.execute(
            "SELECT status, error_code, COUNT(*) AS n, MIN(id) AS example_task, "
            "       MIN(error) AS example "
            "FROM tasks WHERE run_id=? AND status IN ('failed','rejected') "
            "GROUP BY status, error_code ORDER BY n DESC", (run_id,))]
    reuse = reg.conn.execute(
        "SELECT SUM(structure_created=1) AS new_structures, "
        "       SUM(structure_created=0) AS reused_structures, "
        "       SUM(geometry_created=1) AS new_geometries, "
        "       SUM(geometry_created=0) AS reused_geometries, "
        "       SUM(attempts > 1)       AS retried "
        "FROM tasks WHERE run_id=? AND status='done'", (run_id,)).fetchone()
    # A retried embed means ETKDG failed and the random-coordinates fallback carried the
    # geometry.  Cheap to spot with LIKE; the alternative is decoding every detail blob.
    embeds = reg.conn.execute(
        "SELECT COUNT(*) AS n FROM tasks WHERE run_id=? AND detail_json LIKE ?",
        (run_id, '%"retried": true%')).fetchone()
    # The near-miss ledger.  `geometry.qc.MARGINAL_OVERLAP` is a compute budget set from
    # one measurement, and the only way to know whether it is set well is to count how
    # often the gamble paid: built_despite_qc is the bet, qc_failed_after_relax is the
    # loss.  Reported per run so the threshold can be widened on evidence rather than on
    # the feeling that it is probably fine.
    marginal = reg.conn.execute(
        "SELECT SUM(detail_json LIKE '%\"built_despite_qc\"%')   AS built, "
        "       SUM(detail_json LIKE '%\"qc_failed_after_relax\"%') AS unresolved "
        "FROM tasks WHERE run_id=?", (run_id,)).fetchone()
    built = marginal["built"] or 0
    unresolved = marginal["unresolved"] or 0
    relax_reused = reg.conn.execute(
        "SELECT COUNT(*) AS n FROM tasks WHERE run_id=? "
        "AND json_extract(detail_json, '$.relax_reused') IS NOT NULL", (run_id,)).fetchone()
    return {
        "counts": counts,
        "by_code": by_code,
        "new_structures": reuse["new_structures"] or 0,
        "reused_structures": reuse["reused_structures"] or 0,
        "new_geometries": reuse["new_geometries"] or 0,
        "reused_geometries": reuse["reused_geometries"] or 0,
        "retried_tasks": reuse["retried"] or 0,
        "embed_retries": embeds["n"] or 0,
        "marginal_built": built,
        "marginal_unresolved": unresolved,
        "marginal_rescued": max(0, built - unresolved),
        "relax_reused": relax_reused["n"] or 0,
        "timing": task_timing(reg, run_id),
    }


def task_timing(reg: Registry, run_id: int) -> dict[str, dict[str, Any]]:
    """Per task kind: how many were timed, total seconds, median and p90 in ms.

    From `detail.duration_ms` — the worker's own wall time for the task, including a
    failed or rejected one, since that time was spent too.  A task from before timing was
    recorded has none and is not counted, rather than counted as zero.
    """
    by_kind: dict[str, list[int]] = {}
    for r in reg.conn.execute(
            "SELECT kind, json_extract(detail_json, '$.duration_ms') AS ms FROM tasks "
            "WHERE run_id=? AND json_extract(detail_json, '$.duration_ms') IS NOT NULL",
            (run_id,)):
        by_kind.setdefault(r["kind"], []).append(int(r["ms"]))
    out: dict[str, dict[str, Any]] = {}
    for kind, ms in sorted(by_kind.items()):
        ms.sort()
        out[kind] = {"n": len(ms), "total_s": round(sum(ms) / 1000, 1),
                     "median_ms": ms[len(ms) // 2],
                     "p90_ms": ms[min(len(ms) - 1, int(0.9 * len(ms)))]}
    return out


TERMINAL_RUN_STATUSES = (DONE, FAILED, CANCELLED, INTERRUPTED)


def run_liveness(reg: Registry, run_id: int) -> dict[str, Any]:
    """Is anything actually working on this run?  Reads only; decides nothing.

    A run row records that work was STARTED and that it FINISHED.  Nothing writes the
    third outcome — the process stopped existing — so a killed worker leaves a row that
    is indistinguishable from one that is busy, and the page has no choice but to call
    it live.  This is the missing fact, and it is derived from two signals rather than
    one because neither is sufficient alone:

    * the **claimants** (`host:pid`), which on their own host answer the question exactly;
    * the **heartbeat**, which answers it approximately and is all there is from another
      host or on Windows.

    Three verdicts, and the middle one is the point.  `live` and `interrupted` are
    claims; `unknown` is what is returned whenever the evidence does not support either,
    and a caller that sweeps must treat it as live.  Reporting "cannot tell from here" is
    the only honest answer about a process on a machine this one cannot see.

    Note the asymmetry in what each signal is allowed to prove.  A dead pid is proof and
    is acted on immediately.  A stale heartbeat is NOT, while any task is claimed: a
    worker inside an hour-long relaxation stamps nothing for an hour, and sweeping it
    would hand its tasks to a second worker and have the work done twice.  So the
    heartbeat only ever decides a run with nothing in flight — where there is, by
    definition, nothing to disturb.
    """
    row = reg.conn.execute(
        "SELECT status, heartbeat_at, created_at FROM runs WHERE id=?",
        (run_id,)).fetchone()
    if row is None:
        return {"verdict": "unknown", "reason": f"no run {run_id}", "workers": []}
    if row["status"] in TERMINAL_RUN_STATUSES:
        return {"verdict": "finished", "reason": f"the run is {row['status']}",
                "workers": [], "status": row["status"]}

    counts = task_counts(reg, run_id)
    claimants = [r["claimed_by"] for r in reg.conn.execute(
        "SELECT DISTINCT claimed_by FROM tasks WHERE run_id=? AND status=?",
        (run_id, CLAIMED))]
    workers = [{"worker": c, "alive": worker_alive(c)} for c in claimants if c]
    age = _age_seconds(row["heartbeat_at"] or row["created_at"])

    if workers:
        if any(w["alive"] for w in workers):
            return {"verdict": "live", "reason": "a worker that claimed a task is running",
                    "workers": workers, "age": age}
        if all(w["alive"] is False for w in workers):
            return {"verdict": "interrupted",
                    "reason": (f"every process holding a task is gone "
                               f"({', '.join(w['worker'] for w in workers)})"),
                    "workers": workers, "age": age}
        return {"verdict": "unknown",
                "reason": ("tasks are held by a process this machine cannot ask about "
                           f"({', '.join(w['worker'] for w in workers)})"),
                "workers": workers, "age": age}

    if counts.get(PENDING):
        if age is not None and age > STALE_AFTER:
            return {"verdict": "interrupted",
                    "reason": (f"{counts[PENDING]} task(s) are waiting and nothing has "
                               f"claimed one for {int(age)}s"),
                    "workers": [], "age": age}
        return {"verdict": "live", "reason": "tasks are waiting to be claimed",
                "workers": [], "age": age}

    # Nothing pending, nothing claimed: the work IS done and only the closing write is
    # missing.  That is a finished run, not an interrupted one, and `sweep_interrupted`
    # closes it out with the status its tasks earned.
    return {"verdict": "unfinalised",
            "reason": "every task is accounted for but the run was never closed out",
            "workers": [], "age": age}


def sweep_interrupted(reg: Registry, *, run_id: int | None = None) -> list[dict[str, Any]]:
    """Close out runs nothing is working on any more.  Returns what it changed.

    Called when a process TAKES OVER a database — a server starting, a run starting, a
    resume — rather than on a timer.  That is the moment the question is both worth
    asking and safely answerable: something new is about to work here, so a run that no
    longer has a process is a run whose tasks should be claimable again.

    Two outcomes, and they are different facts.  A run whose work is all accounted for is
    `finish_run`ed with the status its tasks earned — it really did finish, and only the
    closing write was lost. A run with work left becomes `interrupted`, and its claimed
    tasks go back to `pending` so the next worker can take them: the process holding them
    is gone, so they are not in flight, they are stranded.
    """
    swept: list[dict[str, Any]] = []
    sql = ("SELECT id, status FROM runs WHERE status NOT IN "
           f"({','.join('?' * len(TERMINAL_RUN_STATUSES))})")
    args: list[Any] = list(TERMINAL_RUN_STATUSES)
    if run_id is not None:
        sql += " AND id=?"
        args.append(run_id)
    for row in reg.conn.execute(sql, tuple(args)).fetchall():
        rid = int(row["id"])
        verdict = run_liveness(reg, rid)
        if verdict["verdict"] == "unfinalised":
            swept.append({"run_id": rid, "was": "unfinalised", "now": finish_run(reg, rid),
                          "returned_claims": 0, "reason": verdict["reason"]})
            continue
        if verdict["verdict"] != "interrupted":
            continue
        if row["status"] == CANCELLING:
            # Asked to stop, and it stopped — the worker dying on the way out does not
            # turn a decision into an accident.  `cancelled`, as the person requested.
            counts = finalise_cancel(reg, rid)
            swept.append({"run_id": rid, "was": CANCELLING, "now": CANCELLED,
                          "returned_claims": counts["returned_claims"],
                          "reason": "the stop was requested; the worker is gone"})
            continue
        returned = reset_stale_claims(reg, rid)
        reg.conn.execute("UPDATE runs SET status=?, finished_at=? WHERE id=?",
                         (INTERRUPTED, utcnow(), rid))
        swept.append({"run_id": rid, "was": "running", "now": INTERRUPTED,
                      "returned_claims": returned, "reason": verdict["reason"]})
    if swept:
        reg.conn.commit()
    return swept


def request_cancel(reg: Registry, run_id: int) -> bool:
    """Ask a run to stop.  Cooperative: nothing is killed.

    Sets the run to `cancelling` and returns whether the request changed anything.  The
    workers notice at their next claim, finish the task in hand, and stop — so every
    structure already computed stays in the registry and no write is interrupted
    part-way.  That matters more here than a fast stop: a worker killed mid-write is
    exactly the concurrency case ground rule 1's single insert path exists to avoid, and
    an hour of xTB thrown away to save ten seconds is a bad trade.
    """
    cur = reg.conn.execute(
        "UPDATE runs SET status=? WHERE id=? AND status IN (?,?,?)",
        (CANCELLING, run_id, PENDING, "running", CANCELLING))
    reg.conn.commit()
    return cur.rowcount > 0


def cancel_requested(reg: Registry, run_id: int) -> bool:
    """Has a stop been asked for?  Read fresh every time — the flag arrives from another
    process (the UI thread, or another worker), so a cached answer is a worker that
    never stops."""
    row = reg.conn.execute("SELECT status FROM runs WHERE id=?", (run_id,)).fetchone()
    return bool(row) and row["status"] in (CANCELLING, CANCELLED)


def finalise_cancel(reg: Registry, run_id: int) -> dict[str, int]:
    """Close out a stopped run: pending work is marked cancelled, not failed.

    Claimed-but-unfinished tasks are returned to `pending` first, so a resumed run picks
    them up rather than leaving them stranded in a state no one will claim.
    """
    returned = reset_stale_claims(reg, run_id)
    cur = reg.conn.execute(
        "UPDATE tasks SET status=?, finished_at=? WHERE run_id=? AND status=?",
        (CANCELLED, utcnow(), run_id, PENDING))
    reg.conn.execute("UPDATE runs SET status=?, finished_at=? WHERE id=?",
                     (CANCELLED, utcnow(), run_id))
    reg.conn.commit()
    return {"cancelled_tasks": cur.rowcount, "returned_claims": returned}


RESUMABLE = (CANCELLING, CANCELLED, INTERRUPTED)


def resume_run(reg: Registry, run_id: int) -> int:
    """Undo a stop.  Returns how many tasks were revived.

    Three states arrive here and all are legitimate.  A run already `cancelled` has its
    tasks parked in `cancelled` and they are put back to `pending`.  A run still
    `cancelling` — stopped, but no worker has closed it out yet — has tasks that never
    left `pending`, so nothing is revived and clearing the flag is the whole job.  Zero
    revived is therefore a success, not a no-op, which is why the caller is told which
    case it was rather than just a count.

    An `interrupted` run is the third: its tasks were returned to `pending` by the sweep
    that noticed it, so resuming it is also just the flag — and it is the state where
    resuming matters most, because nobody chose it.

    A task that FAILED only because the database was busy is revived too: that says
    nothing about the task (B26), and a resume is exactly when it should be tried again.
    Any other failure stays failed — it is a finding, or a bug, not a queue position.

    A run left with nothing to do is closed out with the status its tasks earned rather
    than marked `running`: a running run with no work and nobody executing it is the
    state that left the page saying "stopping" for ever.
    """
    row = reg.conn.execute("SELECT status FROM runs WHERE id=?", (run_id,)).fetchone()
    if row is None or row["status"] not in RESUMABLE:
        return 0
    cur = reg.conn.execute(
        "UPDATE tasks SET status=?, finished_at=NULL WHERE run_id=? AND status=?",
        (PENDING, run_id, CANCELLED))
    busy = reg.conn.execute(
        "UPDATE tasks SET status=?, finished_at=NULL, error=NULL, error_code=NULL "
        "WHERE run_id=? AND status=? AND error LIKE '%database is locked%'",
        (PENDING, run_id, FAILED))
    reg.conn.execute("UPDATE runs SET status=?, finished_at=NULL, heartbeat_at=? "
                     "WHERE id=?", ("running", utcnow(), run_id))
    reg.conn.commit()
    if not task_counts(reg, run_id).get(PENDING):
        finish_run(reg, run_id)
    return cur.rowcount + busy.rowcount


def finish_run(reg: Registry, run_id: int) -> str:
    row = reg.conn.execute("SELECT status FROM runs WHERE id=?", (run_id,)).fetchone()
    if row and row["status"] in (CANCELLING, CANCELLED):
        # A stop that was asked for must not be relabelled `done` by the worker loop
        # exiting normally, which is exactly what it does on the way out.
        finalise_cancel(reg, run_id)
        return CANCELLED
    counts = task_counts(reg, run_id)
    # Work still in the queue is not a finished run.  A pool of processes can lose one,
    # and its share of the queue stays `pending` while the survivors' successes are all
    # that is left to summarise — `done` would be the report that hides it.
    unfinished = counts.get(PENDING, 0) + counts.get(CLAIMED, 0)
    status = FAILED if (counts.get(FAILED) or unfinished) else DONE
    reg.conn.execute("UPDATE runs SET status=?, finished_at=? WHERE id=?",
                     (status, utcnow(), run_id))
    return status


def task_counts(reg: Registry, run_id: int) -> dict[str, int]:
    return {r["status"]: r["n"] for r in reg.conn.execute(
        "SELECT status, COUNT(*) AS n FROM tasks WHERE run_id=? GROUP BY status",
        (run_id,))}


def iter_tasks(reg: Registry, run_id: int, status: str | None = None) -> Iterator[Any]:
    sql = "SELECT * FROM tasks WHERE run_id=?" + (" AND status=?" if status else "")
    yield from reg.conn.execute(sql, (run_id, status) if status else (run_id,))


def reset_stale_claims(reg: Registry, run_id: int) -> int:
    """Return claimed-but-unfinished tasks to the queue (a worker died mid-task)."""
    cur = reg.conn.execute(
        "UPDATE tasks SET status=?, claimed_by=NULL WHERE run_id=? AND status=?",
        (PENDING, run_id, CLAIMED))
    return cur.rowcount
