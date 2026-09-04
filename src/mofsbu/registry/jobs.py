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


def worker_id() -> str:
    return f"{socket.gethostname()}:{os.getpid()}"


@dataclass(frozen=True)
class Task:
    id: int
    run_id: int
    kind: str
    payload: dict[str, Any]
    attempts: int


def create_run(reg: Registry, spec: BuildSpec, *, note: str = "") -> int:
    cur = reg.conn.execute(
        "INSERT INTO runs (spec_digest, spec_json, status, note, host, created_at) "
        "VALUES (?,?,?,?,?,?)",
        (spec.digest, spec.to_json(indent=None), PENDING, note or spec.note,
         socket.gethostname(), utcnow()))
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
                  geometry_id: int | None = None) -> None:
    reg.conn.execute(
        "UPDATE tasks SET status=?, structure_id=?, geometry_id=?, error=NULL, "
        "finished_at=? WHERE id=?",
        (DONE, structure_id, geometry_id, utcnow(), task_id))


def fail_task(reg: Registry, task_id: int, error: str, *, rejected: bool = False) -> None:
    """`rejected` means the task was well-formed and the answer was no.

    A candidate that fails QC is not a failure of the machinery, and conflating the two
    makes a run that produced nothing look identical to a run that crashed.
    """
    reg.conn.execute(
        "UPDATE tasks SET status=?, error=?, finished_at=? WHERE id=?",
        (REJECTED if rejected else FAILED, error[:2000], utcnow(), task_id))


def finish_run(reg: Registry, run_id: int) -> str:
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
