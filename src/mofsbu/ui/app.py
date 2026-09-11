"""FastAPI application for the read-only registry viewer.

Design constraints, all load-bearing:

* **Read-only by construction.**  Every connection is opened with
  ``file:...?mode=ro`` plus ``PRAGMA query_only``.  There is no write path in this
  module at all, so a write attempt fails loudly at the SQLite layer.
* **Reads the ``v_structures`` view, never the base tables**, so the registry's
  internal layout can change underneath the viewer (schema.sql says so explicitly).
* **Parameterised SQL only.**  The only identifiers that ever reach a query string
  are drawn from module-level whitelists (:data:`SORTABLE`, :data:`BRIDGE_CLASSES`).
* **The database is a live reference, not a captured constant.**  Every endpoint asks
  :class:`~mofsbu.ui.active.ActiveDatabase` where the registry is at request time, so the
  builder's database switch moves the viewer with it.  It is still a *path* — the
  connections opened from it are as read-only as they ever were.
* **Hidden structures are filtered, not deleted.**  ``hidden`` is a soft-delete flag
  (see ``registry/api.set_hidden`` for why a real delete is refused).  Listings exclude
  it by default and ``include_hidden=1`` brings it back, so nothing becomes unreachable.
* **Reserved columns degrade, they do not crash.**  ``donor_types``,
  ``binding_modes``, ``n_open_sites`` (M4/M5) and ``reactions.depth`` are NULL in a
  current database.  Their filters keep strict SQL semantics (NULL never matches),
  and ``/api/filters`` reports which of them are unpopulated so the UI can disable
  the control rather than omit it.
"""
from __future__ import annotations

import json
import sqlite3
import urllib.parse
from pathlib import Path
from typing import Any, Callable, Iterable

from fastapi import Depends, FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse

from mofsbu.registry.store import BlobStore

STATIC = Path(__file__).parent / "static"

#: Columns a client may sort by.  Anything else is a 400 — never interpolated.
SORTABLE: frozenset[str] = frozenset({
    # display_label is absent on purpose: it is a projection, not a sort key.
    "id", "formula", "l0_composition", "block_id",
    "n_metals", "n_atoms", "net_charge", "multiplicity", "n_dative_bonds",
    "max_bridge_class", "has_metal_metal", "n_open_sites",
    "best_fidelity", "best_energy", "best_n_atoms", "best_converged",
    "n_geometries", "n_incoming_routes", "created_at",
})

BRIDGE_CLASSES: tuple[str, ...] = ("none", "terminal", "mu2", "mu3", "muN")

#: Reserved-until-a-later-milestone columns the UI must show but may not rely on.
RESERVED = {
    "donor_types": ("structures", "M4"),
    "binding_modes": ("structures", "M5"),
    "n_open_sites": ("structures", "M4"),
    "depth": ("reactions", "M8"),
    "choice_vector_json": ("geometries", "M5"),
}

MAX_LIMIT = 1000


# ── connection ───────────────────────────────────────────────────────────────

def read_only_uri(db_path: Path | str) -> str:
    """``file:`` URI that SQLite will only ever open for reading."""
    posix = Path(db_path).expanduser().resolve().as_posix()
    return "file:" + urllib.parse.quote(posix, safe="/:") + "?mode=ro"


def open_read_only(db_path: Path | str) -> sqlite3.Connection:
    """A connection that cannot write: read-only URI *and* ``query_only``."""
    # check_same_thread=False is required, not a shortcut: FastAPI runs a synchronous
    # `yield` dependency in one threadpool thread and the endpoint that consumes it in
    # another, so a per-request connection legitimately crosses threads.  It is never
    # used by two threads at once — one request owns it from open to close — and it is
    # read-only besides.  Without this every request fails with "SQLite objects created
    # in a thread can only be used in that same thread".
    con = sqlite3.connect(read_only_uri(db_path), uri=True, check_same_thread=False)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    return con


def _rows(cur: Iterable[sqlite3.Row]) -> list[dict[str, Any]]:
    return [dict(r) for r in cur]


def has_column(con: sqlite3.Connection, table: str, column: str) -> bool:
    """Does this database have that column yet?

    The viewer opens `mode=ro` and therefore CANNOT migrate: a registry written before a
    column existed stays that way until something writable touches it.  So every query
    over a column added after M3.5 has to ask first, or a viewer pointed at an older
    file dies on `no such column` — which is a packaging accident presented as a bug in
    the page.  `table` is never user input; it comes from a literal in this module.
    """
    try:
        return any(r[1] == column for r in con.execute(f"PRAGMA table_info({table})"))
    except sqlite3.Error:
        return False


def _like(term: str) -> str:
    """Escape LIKE wildcards so a user's ``%`` searches for a literal ``%``."""
    return "%" + term.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "%"


# ── filter → SQL ─────────────────────────────────────────────────────────────

def _where(params: dict[str, Any]) -> tuple[str, list[Any]]:
    """Build the WHERE clause for /api/structures.  Values are always bound."""
    sql: list[str] = []
    args: list[Any] = []

    def add(clause: str, *vals: Any) -> None:
        sql.append(clause)
        args.extend(vals)

    # ── axis 1: composition / metal / charge / spin ──
    if q := params.get("q"):
        # Free text is STRUCTURAL: whole-structure formula, the L0 key, or the formula of
        # any ligand fragment.  It deliberately does not search display_label — a label is
        # composed after retrieval and must never be an input to it, or a query starts
        # depending on whether someone had named the thing.
        add("(v.formula LIKE ? ESCAPE '\\' OR v.l0_composition LIKE ? ESCAPE '\\'"
            " OR v.id IN (SELECT structure_id FROM structure_fragments"
            "             WHERE formula LIKE ? ESCAPE '\\'))", *([_like(q)] * 3))
    if fragment := params.get("fragment"):
        add("v.id IN (SELECT structure_id FROM structure_fragments WHERE fragment_l1 = ?)",
            fragment)
    if metals := params.get("metal"):
        # `metals` is stored comma-delimited and comma-bracketed (',Cu,Cu,')
        # precisely so an element match is a LIKE and not a join.
        add("(" + " OR ".join(["v.metals LIKE ?"] * len(metals)) + ")",
            *[f"%,{m},%" for m in metals])
    if (x := params.get("n_metals_min")) is not None:
        add("v.n_metals >= ?", x)
    if (x := params.get("n_metals_max")) is not None:
        add("v.n_metals <= ?", x)
    if (x := params.get("charge_min")) is not None:
        add("v.net_charge >= ?", x)
    if (x := params.get("charge_max")) is not None:
        add("v.net_charge <= ?", x)
    if (x := params.get("multiplicity")) is not None:
        add("v.multiplicity = ?", x)

    # ── axis 2: connectivity / binding mode ──
    if bridges := params.get("max_bridge_class"):
        add("v.max_bridge_class IN (" + ",".join("?" * len(bridges)) + ")", *bridges)
    if (x := params.get("has_metal_metal")) is not None:
        add("v.has_metal_metal = ?", int(x))
    if donors := params.get("donor_type"):                       # RESERVED M4
        add("(" + " OR ".join(["v.donor_types LIKE ?"] * len(donors)) + ")",
            *[f"%,{d},%" for d in donors])
    if (x := params.get("n_open_sites_min")) is not None:        # RESERVED M4
        add("v.n_open_sites >= ?", x)

    # ── axis 3: energy / fidelity ──
    if (x := params.get("fidelity_min")) is not None:
        add("COALESCE(v.best_fidelity, v.best_geom_fidelity) >= ?", x)
    if (x := params.get("converged")) is not None:
        add("v.best_converged = ?", int(x))
    if (x := params.get("energy_min")) is not None:
        add("v.best_energy >= ?", x)
    if (x := params.get("energy_max")) is not None:
        add("v.best_energy <= ?", x)

    # ── axis 4: provenance ──
    if (x := params.get("has_route")) is not None:
        add("v.n_incoming_routes > 0" if x else "v.n_incoming_routes = 0")
    if (x := params.get("depth_max")) is not None:               # RESERVED M8
        add("(SELECT MIN(r.depth) FROM reactions r"
            " WHERE r.product_structure_id = v.id) <= ?", x)

    if tag := params.get("tag"):
        add("EXISTS (SELECT 1 FROM structure_tags t"
            " WHERE t.structure_id = v.id AND t.tag = ?)", tag)

    # The soft delete.  Hidden rows are excluded from every listing unless asked for,
    # which is the only thing "hidden" means — the row, its geometries and its
    # provenance are all still there, and `include_hidden` is how you get back to a
    # structure you hid by mistake.
    if params.get("hidden_supported") and not params.get("include_hidden"):
        add("v.hidden = 0")

    return (" WHERE " + " AND ".join(sql)) if sql else "", args


DEPTH_EXPR = ("(SELECT MIN(r.depth) FROM reactions r "
              "WHERE r.product_structure_id = v.id) AS min_depth")


# ── app ──────────────────────────────────────────────────────────────────────

def create_app(db_path: Path, store_root: Path,
               active: "ActiveDatabase | None" = None) -> FastAPI:
    """Build the viewer app.  No server is started; tests use ``TestClient``.

    `active` is the switchable registry reference shared with the builder.  It defaults
    to a fresh one wrapping `db_path`, so the two-argument call every test and script
    already makes keeps working and simply cannot switch.
    """
    from mofsbu.ui.active import ActiveDatabase

    db_path = Path(db_path)
    active = active or ActiveDatabase(db_path)
    store = BlobStore(store_root)

    app = FastAPI(
        title="mofsbu viewer",
        description="Read-only browser over the mofsbu structure registry (M3.5).",
        version="0.1",
    )
    # The builder writes specs, library entries and run/task rows.  It is mounted here for
    # convenience but shares nothing with the read-only registry connections above: see
    # ui/builder.py for why that boundary is kept.
    from mofsbu.ui.builder import build_router

    # The blob store is deliberately NOT switched alongside the database.  It is
    # content-addressed: a geometry's .xyz is keyed by the hash of its own text, so two
    # registries that both contain a structure point at one blob rather than at two
    # copies of it.  Giving each database its own store would duplicate every shared
    # geometry and gain nothing.
    app.include_router(build_router(db_path, store_root,
                                    Path(db_path).parent / "specs", active=active))
    app.state.db_path = db_path                  # where it STARTED; see active for now
    app.state.active = active
    app.state.store = store
    # Exposed so callers (and the read-only test) can get a connection the same way
    # the endpoints do.  There is deliberately no writable counterpart.
    app.state.connect = lambda: open_read_only(active.path)

    def db() -> Any:
        path = active.path
        if not path.exists():
            raise HTTPException(503, f"registry not found: {path}")
        con = open_read_only(path)
        try:
            yield con
        finally:
            con.close()

    Con = Depends(db)

    # ── listing ──────────────────────────────────────────────────────────────
    @app.get("/api/structures")
    def list_structures(                                          # noqa: PLR0913
        con: sqlite3.Connection = Con,
        q: str | None = Query(None, description="free text over formula / L0 / fragment formula"),
        fragment: str | None = Query(None, description="fragment L1 hash — structural search"),
        metal: list[str] | None = Query(None, description="element symbol, repeatable"),
        n_metals_min: int | None = None,
        n_metals_max: int | None = None,
        charge_min: int | None = None,
        charge_max: int | None = None,
        multiplicity: int | None = None,
        max_bridge_class: list[str] | None = Query(None),
        has_metal_metal: bool | None = None,
        donor_type: list[str] | None = Query(None, description="populated by put_sites"),
        n_open_sites_min: int | None = Query(None,
            description="open sites on the best geometry; NULL until state is computed"),
        fidelity_min: int | None = None,
        converged: bool | None = None,
        energy_min: float | None = None,
        energy_max: float | None = None,
        has_route: bool | None = None,
        depth_max: int | None = Query(None, description="RESERVED M8"),
        tag: str | None = None,
        include_hidden: bool = Query(False,
            description="include soft-deleted structures (nothing is ever really deleted)"),
        sort: str = "id",
        order: str = "asc",
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        if sort not in SORTABLE:
            raise HTTPException(400, f"cannot sort by {sort!r}; allowed: "
                                     f"{', '.join(sorted(SORTABLE))}")
        if order.lower() not in ("asc", "desc"):
            raise HTTPException(400, "order must be 'asc' or 'desc'")
        for b in max_bridge_class or []:
            if b not in BRIDGE_CLASSES:
                raise HTTPException(400, f"unknown bridge class {b!r}; allowed: "
                                         f"{', '.join(BRIDGE_CLASSES)}")
        limit = max(1, min(int(limit), MAX_LIMIT))
        offset = max(0, int(offset))

        # A registry written before the soft delete existed has no `hidden` column, and
        # the viewer cannot migrate it (mode=ro).  Nothing can be hidden in such a file,
        # so the filter is simply not applied rather than the query failing.
        hidden_supported = has_column(con, "structures", "hidden")
        clause, args = _where({**locals(), "hidden_supported": hidden_supported})
        total = con.execute(f"SELECT COUNT(*) FROM v_structures v{clause}", args).fetchone()[0]
        # `sort`/`order` are whitelisted above; every value is still bound.
        rows = _rows(con.execute(
            f"SELECT v.*, {DEPTH_EXPR} FROM v_structures v{clause}"
            f" ORDER BY v.{sort} IS NULL, v.{sort} {order.upper()}, v.id ASC"
            f" LIMIT ? OFFSET ?", [*args, limit, offset]))
        n_hidden = con.execute("SELECT COUNT(*) FROM structures WHERE hidden = 1"
                               ).fetchone()[0] if hidden_supported else 0
        return {"total": total, "limit": limit, "offset": offset,
                "sort": sort, "order": order.lower(), "rows": rows,
                "n_hidden": n_hidden, "include_hidden": bool(include_hidden)}

    # ── one structure ────────────────────────────────────────────────────────
    @app.get("/api/structures/{structure_id}")
    def get_structure(structure_id: int, con: sqlite3.Connection = Con) -> dict[str, Any]:
        row = con.execute(
            f"SELECT v.*, {DEPTH_EXPR} FROM v_structures v WHERE v.id = ?",
            (structure_id,)).fetchone()
        if row is None:
            raise HTTPException(404, f"no structure {structure_id}")
        row = dict(row)
        geoms = _rows(con.execute(
            """SELECT g.id, g.fidelity, g.energy, g.converged, g.n_atoms, g.coords_hash,
                      g.l3_conformer_id, g.created_at, g.relaxed_from,
                      g.choice_vector_json IS NOT NULL AS has_choice_vector,
                      m.code, m.code_version, m.method, m.solvent
               FROM geometries g LEFT JOIN methods m ON m.id = g.method_id
               WHERE g.structure_id = ?
               ORDER BY g.fidelity DESC, g.id ASC""", (structure_id,)))
        for g in geoms:
            g["is_best"] = g["id"] == row["best_geometry_id"]
            # cheap, and it stops the UI offering a geometry whose blob is gone
            g["blob_present"] = _blob_ok(g["coords_hash"])
        routes = _rows(con.execute(
            """SELECT id, kind, intermediate, depth, note, fidelity, created_at
               FROM reactions WHERE product_structure_id = ? ORDER BY id""",
            (structure_id,)))
        return {"structure": row, "geometries": geoms, "routes": routes}

    def _blob_ok(digest: str | None) -> bool:
        if not digest:
            return False
        try:
            return store.has(digest)
        except ValueError:                    # not a sha256 digest → treat as absent
            return False

    # ── coordinates ──────────────────────────────────────────────────────────
    @app.get("/api/geometries/{geometry_id}/xyz", response_class=PlainTextResponse)
    def geometry_xyz(geometry_id: int, con: sqlite3.Connection = Con) -> PlainTextResponse:
        row = con.execute("SELECT coords_hash FROM geometries WHERE id = ?",
                          (geometry_id,)).fetchone()
        if row is None:
            raise HTTPException(404, f"no geometry {geometry_id}")
        digest = row["coords_hash"]
        try:
            text = store.get_text(digest)
        except (KeyError, ValueError, OSError):
            # a database can outlive its blobs (store pruned, DB copied alone).
            # That is a 404 about a missing document, not a server fault.
            raise HTTPException(404, f"blob {digest} not in store") from None
        return PlainTextResponse(text, media_type="text/plain; charset=utf-8")

    # ── run spec (copy-paste) ────────────────────────────────────────────────
    @app.get("/api/structures/{structure_id}/spec")
    def get_spec(structure_id: int, geometry_id: int | None = None,
                 con: sqlite3.Connection = Con) -> JSONResponse:
        row = con.execute(
            "SELECT id, l0_composition, l1_graph_hash, l2_isomer_tag, display_label,"
            " formula, best_geometry_id FROM v_structures WHERE id = ?",
            (structure_id,)).fetchone()
        if row is None:
            raise HTTPException(404, f"no structure {structure_id}")
        if geometry_id is not None:
            g = con.execute("SELECT id, choice_vector_json, seed, choice_vector_digest"
                            " FROM geometries WHERE id = ? AND structure_id = ?",
                            (geometry_id, structure_id)).fetchone()
            if g is None:
                raise HTTPException(404, f"no geometry {geometry_id} on structure {structure_id}")
        else:
            g = con.execute(
                """SELECT id, choice_vector_json, seed, choice_vector_digest
                   FROM geometries WHERE structure_id = ? AND choice_vector_json IS NOT NULL
                   ORDER BY fidelity DESC, id ASC LIMIT 1""", (structure_id,)).fetchone()
        cv = g["choice_vector_json"] if g is not None else None
        if cv:
            try:
                spec: Any = json.loads(cv)
            except json.JSONDecodeError:
                spec = {"raw_choice_vector": cv}
            if isinstance(spec, dict):
                spec.setdefault("structure_id", row["id"])
                spec.setdefault("geometry_id", g["id"])
                spec.setdefault("seed", g["seed"])
        else:
            # Empty is the expected pre-M5 state, not an error: choice vectors are
            # written by `construct` from M5 onward.
            spec = {
                "note": "no choice vector recorded (pre-M5 structure)",
                "structure_id": row["id"],
                "l0": row["l0_composition"],
                "l1": row["l1_graph_hash"],
                "l2": row["l2_isomer_tag"],
                "display_label": row["display_label"] or row["formula"],
                "geometry_id": g["id"] if g is not None else row["best_geometry_id"],
            }
        return JSONResponse(spec)

    # ── sidebar vocabulary ───────────────────────────────────────────────────
    @app.get("/api/filters")
    def filters(con: sqlite3.Connection = Con) -> dict[str, Any]:
        def split_set(column: str) -> list[str]:
            out: set[str] = set()
            for (raw,) in con.execute(
                    f"SELECT {column} FROM v_structures WHERE {column} IS NOT NULL"):
                out.update(p for p in str(raw).split(",") if p)
            return sorted(out)

        agg = con.execute(
            """SELECT MIN(net_charge) cmin, MAX(net_charge) cmax,
                      MIN(best_energy) emin, MAX(best_energy) emax,
                      MIN(n_metals) mmin, MAX(n_metals) mmax, COUNT(*) n
               FROM v_structures""").fetchone()
        present = [r["max_bridge_class"] for r in con.execute(
            "SELECT DISTINCT max_bridge_class FROM v_structures")]
        inactive = []
        for col, (table, milestone) in RESERVED.items():
            n = con.execute(f"SELECT COUNT(*) FROM {table} WHERE {col} IS NOT NULL").fetchone()[0]
            if not n:
                inactive.append({"column": col, "table": table, "populated_by": milestone})
        return {
            "metals": split_set("metals"),
            "bridge_classes": [b for b in BRIDGE_CLASSES if b in present],
            "fidelities": [r[0] for r in con.execute(
                "SELECT DISTINCT fidelity FROM geometries ORDER BY fidelity")],
            "donor_types": split_set("donor_types"),          # RESERVED M4 → []
            "multiplicities": [r[0] for r in con.execute(
                "SELECT DISTINCT multiplicity FROM v_structures ORDER BY multiplicity")],
            "tags": [r[0] for r in con.execute(
                "SELECT DISTINCT tag FROM structure_tags ORDER BY tag")],
            "charge": {"min": agg["cmin"], "max": agg["cmax"]},
            "energy": {"min": agg["emin"], "max": agg["emax"]},
            "n_metals": {"min": agg["mmin"], "max": agg["mmax"]},
            "n_structures": agg["n"],
            "sortable": sorted(SORTABLE),
            "reserved_inactive": inactive,
        }

    @app.get("/api/meta")
    def meta(con: sqlite3.Connection = Con) -> dict[str, Any]:
        return {
            "db": str(active.path),
            "db_name": active.path.name,
            "store": str(store.root),
            "read_only": True,
            "algo_versions": {r["name"]: r["version"] for r in
                              con.execute("SELECT name, version FROM algo_versions")},
        }

    # ── the page ─────────────────────────────────────────────────────────────
    @app.get("/", response_class=HTMLResponse)
    def index() -> HTMLResponse:
        page = STATIC / "index.html"
        if not page.exists():                       # pragma: no cover - packaging error
            raise HTTPException(500, "index.html missing from mofsbu/ui/static")
        return HTMLResponse(page.read_text(encoding="utf-8"))

    return app
