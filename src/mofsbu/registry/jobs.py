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
from typing import Any, Iterator

from mofsbu.registry.db import Registry, utcnow
from mofsbu.spec import BuildSpec

PENDING, CLAIMED, DONE, FAILED, REJECTED = "pending", "claimed", "done", "failed", "rejected"
# A build you stopped on purpose is not a build that failed.  Keeping them apart is the
# same rule as `rejected` vs `failed`: a run list where every abandoned experiment reads
# as a crash is a run list nobody trusts.
CANCELLING, CANCELLED = "cancelling", "cancelled"


def worker_id() -> str:
    return f"{socket.gethostname()}:{os.getpid()}"


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
               worker: str | None = None) -> Task | None:
    """Atomically take one pending task, or None if there are none left.

    IMMEDIATE acquires the write lock before reading, so the select-then-update cannot
    interleave with another worker doing the same thing.  Without it two workers happily
    claim the same row and do the same job twice.
    """
    if run_id is not None and cancel_requested(reg, run_id):
        # Defence in depth.  `work()` checks before each claim, but a separate worker
        # PROCESS can be mid-loop when the stop arrives, and the claim is the one
        # chokepoint every worker passes through.
        return None

    who = worker or worker_id()
    conn = reg.conn
    # sqlite3 opens an implicit transaction for DML, so an explicit BEGIN IMMEDIATE
    # collides with it.  The connection is switched to manual control for the duration of
    # the claim only: making it global would turn every multi-statement write (put_structure
    # inserts a structure, its metals and its fragments) into separate autocommits and
    # lose their atomicity.
    conn.commit()
    previous = conn.isolation_level
    conn.isolation_level = None
    try:
        conn.execute("BEGIN IMMEDIATE")
        sql = ("SELECT id, run_id, kind, payload_json, attempts FROM tasks "
               "WHERE status = ?" + (" AND run_id = ?" if run_id else "") +
               " ORDER BY priority DESC, id ASC LIMIT 1")
        args = (PENDING, run_id) if run_id else (PENDING,)
        row = conn.execute(sql, args).fetchone()
        if row is None:
            conn.execute("COMMIT")
            return None
        conn.execute(
            "UPDATE tasks SET status=?, claimed_by=?, claimed_at=?, attempts=attempts+1 "
            "WHERE id=?", (CLAIMED, who, utcnow(), row["id"]))
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    finally:
        conn.isolation_level = previous
    return Task(int(row["id"]), int(row["run_id"]), row["kind"],
                json.loads(row["payload_json"]), int(row["attempts"]) + 1)


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
    }


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


def resume_run(reg: Registry, run_id: int) -> int:
    """Undo a stop.  Returns how many tasks were revived.

    Two states arrive here and both are legitimate.  A run already `cancelled` has its
    tasks parked in `cancelled` and they are put back to `pending`.  A run still
    `cancelling` — stopped, but no worker has closed it out yet — has tasks that never
    left `pending`, so nothing is revived and clearing the flag is the whole job.  Zero
    revived is therefore a success, not a no-op, which is why the caller is told which
    case it was rather than just a count.
    """
    row = reg.conn.execute("SELECT status FROM runs WHERE id=?", (run_id,)).fetchone()
    if row is None or row["status"] not in (CANCELLING, CANCELLED):
        return 0
    cur = reg.conn.execute(
        "UPDATE tasks SET status=?, finished_at=NULL WHERE run_id=? AND status=?",
        (PENDING, run_id, CANCELLED))
    reg.conn.execute("UPDATE runs SET status=?, finished_at=NULL WHERE id=?",
                     ("running", run_id))
    reg.conn.commit()
    return cur.rowcount


def finish_run(reg: Registry, run_id: int) -> str:
    row = reg.conn.execute("SELECT status FROM runs WHERE id=?", (run_id,)).fetchone()
    if row and row["status"] in (CANCELLING, CANCELLED):
        # A stop that was asked for must not be relabelled `done` by the worker loop
        # exiting normally, which is exactly what it does on the way out.
        finalise_cancel(reg, run_id)
        return CANCELLED
    counts = task_counts(reg, run_id)
    status = FAILED if counts.get(FAILED) else DONE
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
