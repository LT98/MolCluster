"""Schema changes must reach databases that already exist.

`CREATE TABLE IF NOT EXISTS` does nothing to a table that is already there, so every
column added after a database was first written was silently missing from it.  That
surfaced as `no such column: diagnostics_json` on a registry created one revision earlier
— a schema file is not a migration.
"""
from __future__ import annotations

import sqlite3

import pytest

from mofsbu.registry import BlobStore, Registry
from mofsbu.registry.db import SCHEMA_PATH


def _old_schema_without(column: str, table: str) -> str:
    """The current schema with one column removed — an 'older' database."""
    lines = SCHEMA_PATH.read_text(encoding="utf-8").splitlines()
    kept, in_table, dropped = [], False, False
    for line in lines:
        stripped = line.strip()
        if stripped.upper().startswith(f"CREATE TABLE IF NOT EXISTS {table.upper()}") or \
                stripped.startswith(f"CREATE TABLE IF NOT EXISTS {table}"):
            in_table = True
        if in_table and stripped.startswith(column):
            dropped = True
            continue
        if in_table and stripped.startswith(");"):
            in_table = False
        kept.append(line)
    assert dropped, f"test setup: {column} not found in {table}"
    return "\n".join(kept)


def test_a_column_added_later_reaches_an_existing_database(tmp_path):
    db = tmp_path / "old.db"
    con = sqlite3.connect(db)
    con.executescript(_old_schema_without("diagnostics_json", "runs"))
    con.commit()
    con.close()

    before = {r[1] for r in sqlite3.connect(db).execute("PRAGMA table_info(runs)")}
    assert "diagnostics_json" not in before

    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate("upgrade")

    after = {r[1] for r in sqlite3.connect(db).execute("PRAGMA table_info(runs)")}
    assert "diagnostics_json" in after


def test_the_upgrade_is_recorded(tmp_path):
    db = tmp_path / "old.db"
    con = sqlite3.connect(db)
    con.executescript(_old_schema_without("diagnostics_json", "runs"))
    con.commit()
    con.close()
    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
    notes = [r[0] for r in sqlite3.connect(db).execute("SELECT note FROM migrations")]
    assert any("runs.diagnostics_json" in n for n in notes)


def test_a_stale_view_is_replaced(tmp_path):
    """`CREATE VIEW IF NOT EXISTS` leaves an outdated definition in place."""
    db = tmp_path / "stale.db"
    con = sqlite3.connect(db)
    con.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
    con.execute("DROP VIEW v_structures")
    con.execute("CREATE VIEW v_structures AS SELECT id FROM structures")   # an old shape
    con.commit()
    con.close()

    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
        columns = {d[0] for d in reg.conn.execute("SELECT * FROM v_structures").description}
    assert "display_label" in columns and "min_depth" in columns


def test_migrating_twice_changes_nothing(tmp_path):
    db = tmp_path / "r.db"
    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
        first = sorted(r[1] for r in reg.conn.execute("PRAGMA table_info(runs)"))
        n_before = reg.conn.execute("SELECT COUNT(*) FROM migrations").fetchone()[0]
    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
        second = sorted(r[1] for r in reg.conn.execute("PRAGMA table_info(runs)"))
        n_after = reg.conn.execute("SELECT COUNT(*) FROM migrations").fetchone()[0]
    assert first == second and n_before == n_after


def test_an_unexpected_extra_column_is_reported_not_dropped(tmp_path):
    """Destructive changes are never applied automatically."""
    db = tmp_path / "extra.db"
    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
        reg.conn.execute("ALTER TABLE runs ADD COLUMN somebodys_experiment TEXT")
    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
        columns = {r[1] for r in reg.conn.execute("PRAGMA table_info(runs)")}
        notes = [r[0] for r in reg.conn.execute("SELECT note FROM migrations")]
    assert "somebodys_experiment" in columns
    assert any("not in the schema" in n for n in notes)


def test_an_extra_column_does_not_break_every_later_migration(tmp_path):
    """The note about an extra column must not be re-inserted on every migrate.

    An extra column stays extra, so the note is re-derived every time `migrate()` runs.
    It used to be written at a hard-coded `version = 0` with no conflict handling — the
    same defect `_next_migration_slot` fixes for added columns, left in the one place
    that did not get it.  So the SECOND migrate of any registry carrying an unexpected
    column raised `UNIQUE constraint failed: migrations.version` and rolled the whole
    migration back.

    That is not a cosmetic failure.  `ui.builder.submit_run` migrates before it plans,
    so it turned every run submission after the first into a 500 — on a registry that
    was itself perfectly healthy.
    """
    db = tmp_path / "extra.db"
    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
        reg.conn.execute("ALTER TABLE runs ADD COLUMN somebodys_experiment TEXT")
        # A second table with an extra column too: two notes in ONE migration collided
        # with each other, which is the same bug seen from the other side.
        reg.conn.execute("ALTER TABLE tasks ADD COLUMN another_experiment TEXT")

    for _ in range(3):                                 # the 2nd is what used to raise
        with Registry(db, BlobStore(tmp_path / "store")) as reg:
            reg.migrate()

    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        versions = [r[0] for r in reg.conn.execute("SELECT version FROM migrations")]
        notes = [r[0] for r in reg.conn.execute(
            "SELECT note FROM migrations WHERE note LIKE '%not in the schema%'")]
        columns = {r[1] for r in reg.conn.execute("PRAGMA table_info(runs)")}

    assert len(versions) == len(set(versions)), "migration versions must stay unique"
    # Both tables reported, each exactly once however many times we migrated.
    assert len(notes) == 2, notes
    assert any("runs" in n for n in notes) and any("tasks" in n for n in notes)
    assert "somebodys_experiment" in columns, "and nothing was dropped"


def test_several_columns_can_be_added_in_one_revision(tmp_path):
    """Every added column got the SAME migration version — `SCHEMA_VERSION * 1000 +
    len(added)` — so a revision adding more than one violated the UNIQUE constraint on
    `migrations.version` and the whole migration rolled back.  No revision had ever added
    two, so it sat undiscovered until rev 18 added four to `tasks`; the failure mode is
    that every existing registry refuses to open after a `git pull`.
    """
    db = tmp_path / "old.db"
    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()
        for column in ("error_code", "detail_json", "structure_created",
                       "geometry_created"):
            reg.conn.execute(f"ALTER TABLE tasks DROP COLUMN {column}")
        reg.conn.execute(
            "INSERT INTO tasks (run_id, kind, payload_json, status, created_at) "
            "SELECT 1, 'place', '{}', 'done', '2020-01-01' WHERE EXISTS "
            "(SELECT 1 FROM runs)")

    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()                                  # this is what used to raise
        columns = {r[1] for r in reg.conn.execute("PRAGMA table_info(tasks)")}
    assert {"error_code", "detail_json", "structure_created", "geometry_created"} <= columns

    with Registry(db, BlobStore(tmp_path / "store")) as reg:
        reg.migrate()                                  # and it stays idempotent
        versions = [r[0] for r in reg.conn.execute("SELECT version FROM migrations")]
    assert len(versions) == len(set(versions))
