"""The SQLite registry: connection, migration, and version bookkeeping."""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from mofsbu.registry.store import BlobStore
from mofsbu.versions import ALGO_VERSIONS

SCHEMA_VERSION = 1
SCHEMA_PATH = Path(__file__).with_name("schema.sql")


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class Registry:
    """Owns the connection and the blob store.  Use as a context manager."""

    def __init__(self, db_path: str | Path, store: BlobStore | None = None,
                 *, journal_mode: str = "auto") -> None:
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.store = store or BlobStore(self.db_path.parent / "store")
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA foreign_keys = ON")
        self.journal_mode = self._set_journal_mode(journal_mode)

    def _set_journal_mode(self, requested: str) -> str:
        """WAL where possible, degrading rather than failing.

        A rollback journal needs to be DELETED on commit.  In sandboxes and on some
        network/synced filesystems that delete is refused, and SQLite reports it as a
        bare `disk I/O error` at commit time — far from its cause.  Probe once, here,
        where the fallback is cheap and the diagnosis is obvious.
        """
        import sqlite3 as _sq

        for mode in ([requested] if requested != "auto" else ["WAL", "TRUNCATE", "MEMORY"]):
            try:
                self.conn.execute(f"PRAGMA journal_mode = {mode}")
                self.conn.execute("CREATE TABLE IF NOT EXISTS _probe (x INTEGER)")
                self.conn.execute("INSERT INTO _probe VALUES (1)")
                self.conn.commit()
                self.conn.execute("DROP TABLE _probe")
                self.conn.commit()
                return mode
            except _sq.Error:
                self.conn.rollback()
                continue
        raise RuntimeError(
            f"no usable SQLite journal mode at {self.db_path}. If a stale "
            f"'{self.db_path.name}-journal' file is present and cannot be deleted, remove it."
        )

    # -- lifecycle ---------------------------------------------------------

    def migrate(self, note: str = "") -> None:
        """Bring an existing database up to the current schema.

        `CREATE TABLE IF NOT EXISTS` creates what is missing and does NOTHING to a table
        that already exists, so every column added after a database was first written was
        silently absent from it — which surfaced as `no such column: diagnostics_json` on a
        registry created one revision earlier.  A schema file alone is not a migration.

        So: create what is missing, recreate the views (they are derived definitions and
        `CREATE VIEW IF NOT EXISTS` will not update a stale one), then reconcile columns
        against the target schema.  Additive changes are applied; anything destructive is
        reported rather than guessed at.
        """
        schema = SCHEMA_PATH.read_text(encoding="utf-8")
        self.conn.executescript(schema)
        self._recreate_views(schema)
        # Before the column reconciliation, not after: the rebuild recreates the table
        # from the target schema, so it delivers `role` and `slot` itself and there is
        # nothing left for `_reconcile_columns` to add.  Running it the other way round
        # would add the columns to a table whose UNIQUE key still could not hold them.
        rebuilt = self._rebuild_site_catalog(schema)
        added = self._reconcile_columns(schema)
        for i, (table, column) in enumerate(added, start=1):
            # One row per column, and therefore one VERSION per column.  This used to be
            # `SCHEMA_VERSION * 1000 + len(added)` — the same number for every column in
            # the batch — which is a UNIQUE violation the moment a revision adds more
            # than one.  It never had, so the bug sat here until rev 18 added four to
            # `tasks` at once and every existing registry refused to migrate.
            self.conn.execute(
                "INSERT INTO migrations (version, applied_at, note) VALUES (?,?,?) "
                "ON CONFLICT(version) DO NOTHING",
                (SCHEMA_VERSION * 1000 + self._next_migration_slot(i), utcnow(),
                 f"added {table}.{column}"))
        if rebuilt:
            self.conn.execute(
                "INSERT INTO migrations (version, applied_at, note) VALUES (?,?,?) "
                "ON CONFLICT(version) DO NOTHING",
                (SCHEMA_VERSION * 1000 + self._next_migration_slot(len(added) + 1),
                 utcnow(), "rebuilt site_catalog: UNIQUE now includes slot"))
        cur = self.conn.execute("SELECT 1 FROM migrations WHERE version = ?", (SCHEMA_VERSION,))
        if cur.fetchone() is None:
            self.conn.execute(
                "INSERT INTO migrations (version, applied_at, note) VALUES (?,?,?)",
                (SCHEMA_VERSION, utcnow(), note or f"schema v{SCHEMA_VERSION}"),
            )
        self.record_algo_versions()
        self.conn.commit()

    def _next_migration_slot(self, offset: int) -> int:
        """A free slot above whatever this database has already recorded.

        Migration versions are only ever bookkeeping — they say what was applied, in
        order — so the requirement is uniqueness, not a globally meaningful number.
        """
        top = self.conn.execute(
            "SELECT COALESCE(MAX(version), 0) FROM migrations WHERE version >= ?",
            (SCHEMA_VERSION * 1000,)).fetchone()[0]
        return max(int(top) - SCHEMA_VERSION * 1000, 0) + offset

    def _target_schema(self, schema: str) -> sqlite3.Connection:
        """The schema as it SHOULD be, built in memory so SQLite itself parses it."""
        target = sqlite3.connect(":memory:")
        target.executescript(schema)
        return target

    def _recreate_views(self, schema: str) -> None:
        """Views are pure projections, so replacing them is always safe — and necessary,
        because `CREATE VIEW IF NOT EXISTS` leaves a stale definition in place."""
        target = self._target_schema(schema)
        try:
            for row in target.execute(
                    "SELECT name, sql FROM sqlite_master WHERE type='view'"):
                self.conn.execute(f"DROP VIEW IF EXISTS {row[0]}")
                self.conn.execute(row[1])
        finally:
            target.close()

    def _rebuild_site_catalog(self, schema: str) -> bool:
        """Widen `site_catalog`'s UNIQUE key to include `slot`.  Returns whether it ran.

        The one migration in this file that `_reconcile_columns` cannot do.  Adding a
        column is `ALTER TABLE ADD COLUMN`; CHANGING A CONSTRAINT is not expressible in
        SQLite at all, so the table has to be rebuilt — and a rebuild is exactly the kind
        of destructive step the rest of `migrate` refuses to guess at, which is why this
        one is written out explicitly, guarded, and reported.

        Why the key has to widen: a metal carries several vacant coordination vertices on
        ONE atom, so `UNIQUE (structure_id, canonical_idx)` allows one site per atom and a
        four-coordinate metal with three vacancies needs three rows on the same
        `canonical_idx`.  Donors are unaffected — they are all slot 0, so their uniqueness
        is unchanged, and every existing row migrates to `role='donor', slot=0` which is
        what it already meant.

        `site_state.site_id` references `site_catalog(id)`, so the rebuild PRESERVES `id`
        rather than letting SQLite reassign it; otherwise every per-geometry state row
        would silently point at a different site.  Foreign keys are disabled for the
        swap — the reference is by name and would otherwise be enforced against the
        half-built table mid-rename — and re-enabled afterwards, which the caller's
        `PRAGMA foreign_keys = ON` in `__init__` does not do for us on this connection.
        """
        row = self.conn.execute(
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='site_catalog'"
        ).fetchone()
        if row is None or "slot" in (row[0] or ""):
            return False                        # fresh database, or already widened

        target = self._target_schema(schema)
        try:
            want = target.execute(
                "SELECT sql FROM sqlite_master WHERE type='table' AND name='site_catalog'"
            ).fetchone()[0]
        finally:
            target.close()

        have = [r[1] for r in self.conn.execute("PRAGMA table_info(site_catalog)")]
        carried = ", ".join(have)               # only the columns this database actually has
        self.conn.execute("PRAGMA foreign_keys = OFF")
        try:
            self.conn.execute("DROP TABLE IF EXISTS site_catalog__new")
            self.conn.execute(want.replace("site_catalog", "site_catalog__new", 1))
            self.conn.execute(
                f"INSERT INTO site_catalog__new ({carried}) SELECT {carried} FROM site_catalog")
            self.conn.execute("DROP TABLE site_catalog")
            self.conn.execute("ALTER TABLE site_catalog__new RENAME TO site_catalog")
        finally:
            self.conn.execute("PRAGMA foreign_keys = ON")
        return True

    def _reconcile_columns(self, schema: str) -> list[tuple[str, str]]:
        """Add columns present in the target schema and missing from this database."""
        target = self._target_schema(schema)
        added: list[tuple[str, str]] = []
        try:
            tables = [r[0] for r in target.execute(
                "SELECT name FROM sqlite_master WHERE type='table' "
                "AND name NOT LIKE 'sqlite_%'")]
            for table in tables:
                want = {r[1]: r for r in target.execute(f"PRAGMA table_info({table})")}
                have = {r[1] for r in self.conn.execute(f"PRAGMA table_info({table})")}
                if not have:
                    continue                       # the CREATE above will have made it
                for name, spec in want.items():
                    if name in have:
                        continue
                    _cid, _name, coltype, notnull, default, _pk = spec
                    clause = f"{name} {coltype}"
                    if default is not None:
                        clause += f" DEFAULT {default}"
                    if notnull and default is not None:
                        clause += " NOT NULL"
                    elif notnull:
                        # SQLite cannot add a NOT NULL column without a default; adding it
                        # nullable is the honest compromise and is reported.
                        pass
                    self.conn.execute(f"ALTER TABLE {table} ADD COLUMN {clause}")
                    added.append((table, name))
                extra = have - set(want)
                if extra:
                    self.conn.execute(
                        "INSERT INTO migrations (version, applied_at, note) VALUES (?,?,?)",
                        (0, utcnow(),
                         f"NOTE: {table} has column(s) {sorted(extra)} not in the schema; "
                         "left alone (destructive changes are never applied automatically)"))
        finally:
            target.close()
        return added

    def record_algo_versions(self) -> None:
        for name, version in ALGO_VERSIONS.items():
            self.conn.execute(
                "INSERT INTO algo_versions (name, version, recorded_at) VALUES (?,?,?) "
                "ON CONFLICT(name) DO UPDATE SET version=excluded.version, "
                "recorded_at=excluded.recorded_at WHERE algo_versions.version != excluded.version",
                (name, version, utcnow()),
            )

    def stale_algo_versions(self) -> dict[str, tuple[str, str]]:
        """{name: (stored, current)} for recipes that have changed under the data.

        Ground rule 6: a version bump does not silently re-label rows; it makes them
        detectably stale so they can be recomputed on purpose.
        """
        out: dict[str, tuple[str, str]] = {}
        for row in self.conn.execute("SELECT name, version FROM algo_versions"):
            current = ALGO_VERSIONS.get(row["name"])
            if current is not None and current != row["version"]:
                out[row["name"]] = (row["version"], current)
        return out

    def close(self) -> None:
        self.conn.commit()
        self.conn.close()

    def __enter__(self) -> Registry:
        return self

    def __exit__(self, *exc) -> None:
        if exc[0] is None:
            self.conn.commit()
        else:
            self.conn.rollback()
        self.conn.close()

    # -- small helpers -----------------------------------------------------

    def count(self, table: str) -> int:
        if not table.isidentifier():
            raise ValueError(table)
        return self.conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]

    def __repr__(self) -> str:      # pragma: no cover
        return f"Registry({self.db_path}, {self.count('structures')} structures)"
