"""The controls that moved out of the shell and onto the page (UI backlog items 5-7).

Three things used to require a terminal — choosing a device, choosing a database, and
undoing a mistake — which ruled out most of the people who use the viewer.  These tests
cover the parts of that move where being wrong would be expensive:

* a device or database choice must not widen anyone's access, and must not let a value
  from the browser name a file the server never offered;
* the soft delete must actually be soft: the row, its geometries and every provenance
  edge survive, and the structure stays reachable;
* re-running one entry must replay the task that built it, not compose a new spec.
"""
from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from mofsbu.registry import BlobStore, Registry, incoming_routes, set_hidden  # noqa: E402
from mofsbu.registry.api import RegistryError  # noqa: E402
from mofsbu.ui.active import ActiveDatabase, _validate_name, available_devices  # noqa: E402
from mofsbu.ui.app import create_app, has_column  # noqa: E402
from seed_demo_registry import seed  # noqa: E402


@pytest.fixture
def registry(tmp_path, monkeypatch):
    """A seeded registry in an isolated data root.

    `MOFSBU_DATA` is redirected so `discover_databases()` and "new database…" cannot
    reach the real `data/` folder from a test.
    """
    monkeypatch.setenv("MOFSBU_DATA", str(tmp_path))
    db, store = tmp_path / "demo.db", tmp_path / "store"
    seed(db, store, n=12)
    with Registry(db, BlobStore(store)) as reg:
        reg.migrate("test")           # bring the seeded file up to the current schema
    return db, store


@pytest.fixture
def active(registry):
    return ActiveDatabase(registry[0])


@pytest.fixture
def client(registry, active):
    return TestClient(create_app(registry[0], registry[1], active))


# ── the soft delete ──────────────────────────────────────────────────────────

def test_hiding_removes_it_from_listings_and_nothing_else(client, registry):
    """The whole claim of a soft delete, checked rather than asserted in a comment."""
    db, store = registry
    before = client.get("/api/structures", params={"limit": 100}).json()
    sid = before["rows"][0]["id"]
    detail_before = client.get(f"/api/structures/{sid}").json()

    r = client.post(f"/api/structures/{sid}/hide", json={"hidden": True})
    assert r.status_code == 200, r.text

    after = client.get("/api/structures", params={"limit": 100}).json()
    assert after["total"] == before["total"] - 1
    assert sid not in [row["id"] for row in after["rows"]]
    assert after["n_hidden"] == 1

    # ...and nothing was actually removed.
    detail_after = client.get(f"/api/structures/{sid}").json()
    assert len(detail_after["geometries"]) == len(detail_before["geometries"])
    assert len(detail_after["routes"]) == len(detail_before["routes"])
    with sqlite3.connect(db) as con:
        assert con.execute("SELECT COUNT(*) FROM structures WHERE id=?", (sid,)).fetchone()[0] == 1


def test_hidden_structures_come_back(client):
    sid = client.get("/api/structures", params={"limit": 1}).json()["rows"][0]["id"]
    client.post(f"/api/structures/{sid}/hide", json={"hidden": True})

    shown = client.get("/api/structures",
                       params={"limit": 100, "include_hidden": True}).json()
    assert sid in [row["id"] for row in shown["rows"]]

    client.post(f"/api/structures/{sid}/hide", json={"hidden": False})
    plain = client.get("/api/structures", params={"limit": 100}).json()
    assert sid in [row["id"] for row in plain["rows"]]
    assert plain["n_hidden"] == 0


def test_hiding_a_structure_with_several_routes_is_refused_with_the_routes(client, registry):
    """The policy that makes this safe: a multi-route structure is not hidden silently.

    77% of structures in the real registry have more than one incoming edge, so this is
    the normal case and not an edge case — and the refusal has to name the other routes,
    because "is this safe to remove" is not answerable without them.
    """
    db, store = registry
    sid = client.get("/api/structures", params={"limit": 1}).json()["rows"][0]["id"]
    with Registry(db, BlobStore(store)) as reg:
        from mofsbu.registry import Provenance, put_reaction

        # Two routes to the same product: exactly what an enumeration produces, because
        # re-deriving an identity you already have is most of what an enumeration does.
        while len(incoming_routes(reg, sid)) < 2:
            put_reaction(reg, sid, Provenance(kind="assembly", note="another way there"))
        reg.conn.commit()
        assert len(incoming_routes(reg, sid)) > 1

    r = client.post(f"/api/structures/{sid}/hide", json={"hidden": True})
    assert r.status_code == 409
    detail = r.json()["detail"]
    assert len(detail["routes"]) > 1
    assert "routes" in detail["message"]

    # Still listed: a refusal must not half-apply.
    rows = client.get("/api/structures", params={"limit": 100}).json()["rows"]
    assert sid in [row["id"] for row in rows]

    # And it is overridable, because the user may genuinely mean it.
    assert client.post(f"/api/structures/{sid}/hide",
                       json={"hidden": True, "force": True}).status_code == 200


def test_set_hidden_refuses_an_unknown_structure(registry):
    db, store = registry
    with Registry(db, BlobStore(store)) as reg:
        with pytest.raises(RegistryError):
            set_hidden(reg, 999_999, True)


def test_a_registry_without_the_hidden_column_still_lists(tmp_path, monkeypatch):
    """A viewer pointed at an older registry must not die on `no such column`.

    The viewer opens `mode=ro`, so it CANNOT migrate a file into having the column —
    a registry written before the soft delete existed stays that way until something
    writable touches it.  Renaming the column away reproduces exactly that file.
    """
    monkeypatch.setenv("MOFSBU_DATA", str(tmp_path))
    db, store = tmp_path / "old.db", tmp_path / "store"
    seed(db, store, n=5)

    with sqlite3.connect(db) as w:
        w.execute("ALTER TABLE structures RENAME COLUMN hidden TO hidden_x")
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        assert not has_column(con, "structures", "hidden")
    finally:
        con.close()

    client = TestClient(create_app(db, store))
    r = client.get("/api/structures", params={"limit": 10})
    assert r.status_code == 200, r.text
    assert r.json()["total"] == 5
    # Nothing can be hidden in such a file, so the filter is skipped rather than failing.
    assert r.json()["n_hidden"] == 0


# ── re-running one entry ─────────────────────────────────────────────────────

def test_origin_explains_why_a_seeded_structure_cannot_be_rerun(client):
    """A control that would fail when pressed should be disabled WITH a reason.

    Seeded structures were never built by a task, so there is nothing to replay — and
    that is a fact about provenance, not a bug.
    """
    sid = client.get("/api/structures", params={"limit": 1}).json()["rows"][0]["id"]
    body = client.get(f"/api/structures/{sid}/origin").json()
    assert body["can_rerun"] is False
    assert body["why_not"] and "ingested or seeded" in body["why_not"]

    r = client.post(f"/api/structures/{sid}/rerun")
    assert r.status_code == 404


def test_origin_reports_routes_and_hidden_state(client, registry):
    sid = client.get("/api/structures", params={"limit": 1}).json()["rows"][0]["id"]
    client.post(f"/api/structures/{sid}/hide", json={"hidden": True})
    body = client.get(f"/api/structures/{sid}/origin").json()
    assert body["hidden"] is True
    assert body["n_incoming_routes"] >= 0
    assert body["can_hide"] == (body["n_incoming_routes"] <= 1)


# ── choosing a database ──────────────────────────────────────────────────────

@pytest.mark.parametrize("bad", ["../escape", "a/b", "a\\b", "", "   ", ".hidden"])
def test_a_database_name_is_a_path_component(bad):
    """Same rule `save_spec` applies, for the same reason: it arrives from a browser."""
    with pytest.raises(ValueError):
        _validate_name(bad)


def test_a_database_name_gets_its_extension():
    assert _validate_name("scratch") == "scratch.db"
    assert _validate_name("scratch.db") == "scratch.db"


def test_listing_databases_reports_row_counts(client, registry):
    """The count is the feature: two plausible filenames do not say which is the real one."""
    body = client.get("/api/databases").json()
    assert body["databases"]
    active = [d for d in body["databases"] if d["active"]]
    assert len(active) == 1
    assert active[0]["n_structures"] == 12
    assert "folder" in active[0]


def test_switching_moves_the_viewer_too(client, registry, tmp_path):
    """A page showing structures from one database and runs from another would be worse
    than no switch at all, so the switch has to move everything."""
    r = client.post("/api/databases", json={"name": "scratch"})
    assert r.status_code == 200, r.text
    assert r.json()["created"] is True

    assert client.get("/api/meta").json()["db"].endswith("scratch.db")
    assert client.get("/api/structures", params={"limit": 5}).json()["total"] == 0

    # ...and back again.  A one-way switch is a trap.
    back = client.post("/api/database/select", json={"path": str(registry[0])})
    assert back.status_code == 200, back.text
    assert client.get("/api/structures", params={"limit": 5}).json()["total"] == 12


def test_the_launch_database_stays_selectable_after_switching_away(tmp_path, monkeypatch):
    """The bug this guards: listing only the CURRENT database made the whitelist change
    as you moved, so a launch database outside the data root dropped out of it and the
    way back was a 404."""
    monkeypatch.setenv("MOFSBU_DATA", str(tmp_path / "root"))
    (tmp_path / "root").mkdir()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    db, store = elsewhere / "launch.db", tmp_path / "store"
    seed(db, store, n=4)

    active = ActiveDatabase(db)
    client = TestClient(create_app(db, store, active))
    client.post("/api/databases", json={"name": "other"})
    assert client.get("/api/meta").json()["db"].endswith("other.db")

    back = client.post("/api/database/select", json={"path": str(db.resolve())})
    assert back.status_code == 200, back.text
    assert client.get("/api/structures", params={"limit": 5}).json()["total"] == 4


def test_a_path_the_server_never_offered_is_refused(client):
    """A value from the browser is never joined onto a root and opened."""
    r = client.post("/api/database/select", json={"path": "/etc/passwd"})
    assert r.status_code == 404


def test_creating_a_database_refuses_a_path(client):
    assert client.post("/api/databases", json={"name": "../escape"}).status_code == 400


def test_switching_does_not_make_the_viewer_writable(client, registry, tmp_path):
    """The read-only guarantee is per-connection and must survive a switch."""
    client.post("/api/databases", json={"name": "scratch"})
    client.post("/api/database/select", json={"path": str(registry[0])})
    con = client.app.state.connect()
    try:
        with pytest.raises(sqlite3.OperationalError):
            con.execute("UPDATE structures SET formula='nope'")
    finally:
        con.close()


# ── choosing a device ────────────────────────────────────────────────────────

def test_cpu_is_always_offered_and_always_first(client):
    devices = client.get("/api/compute").json()["devices"]
    assert devices[0]["id"] == "cpu"
    assert [d["id"] for d in devices] == [d["id"] for d in available_devices()]


def test_a_device_is_declared_not_detected(client, monkeypatch):
    """Ground rule 9.  Listing what exists and letting a person pick is a declaration;
    defaulting to CUDA because a card is present is not — so the default never moves."""
    monkeypatch.delenv("MOFSBU_DEVICE", raising=False)
    assert client.get("/api/compute").json()["device"] == "cpu"


def test_declaring_a_device_sticks(client, monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", "cpu")
    body = client.post("/api/compute", json={"device": "cuda:1", "workers": 4}).json()
    assert body["device"] == "cuda:1"
    assert body["workers"] == 4
    # Asking for more than one worker IS declaring a workstation; the two must not
    # disagree, or the CLI and the page report different things about one machine.
    assert body["profile"] == "workstation"


def test_a_typo_is_refused_and_leaves_the_declaration_alone(client, monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", "cpu")
    r = client.post("/api/compute", json={"device": "banana"})
    assert r.status_code == 400
    assert "not a device torch would recognise" in r.json()["detail"]
    # The failed attempt must not have left a broken value behind.
    assert client.get("/api/compute").json()["device"] == "cpu"


def test_capabilities_carries_the_hardware_and_the_database(client):
    """The page renders from this, so an option it cannot offer cannot be mislabelled."""
    caps = client.get("/api/capabilities").json()
    assert caps["compute"]["devices"][0]["id"] == "cpu"
    assert caps["database"].endswith(".db")
