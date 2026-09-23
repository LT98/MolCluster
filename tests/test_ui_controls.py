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

import json
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


# ── how big is this run, before it is one ────────────────────────────────────
# "How many structures will this build?" was answerable only by submitting and reading
# the task count back, which is a decision made after the thing it decides.

THQ = "OC1=C(O)C(=O)C(O)=C(O)C1=O"


def spec_payload(**over) -> dict:
    spec = {
        "spec_version": 7,
        "molecules": [{"name": "THQ", "smiles": THQ, "max_deprotonations": 1}],
        "metals": [{"symbol": "Ni", "oxidation_state": 2, "spin_class": "hs"}],
        "coordination": "4~6", "ligands_per_metal": "1~2", "binding": ["chelate"],
        "co_ligand": "O", "run_mode": "construct",
    }
    spec.update(over)
    return spec


def test_the_estimate_answers_and_writes_nothing(client):
    before = len(client.get("/api/runs").json()["runs"])
    body = client.post("/api/estimate", json={"spec": spec_payload()}).json()
    assert body["ok"] and body["tasks"] > 0
    assert body["by_kind"]["place"] > 0
    assert len(client.get("/api/runs").json()["runs"]) == before


def test_a_wider_range_is_visibly_a_bigger_run(client):
    """The number is what the control is for: `1~3` should look more expensive than `1`."""
    def tasks(**over):
        return client.post("/api/estimate", json={"spec": spec_payload(**over)}
                           ).json()["tasks"]

    assert tasks(ligands_per_metal="1~3") > tasks(ligands_per_metal="1")
    assert tasks(pathways=True) > tasks(pathways=False)


def test_a_range_typed_on_the_page_is_stored_expanded(client, tmp_path):
    """The page sends what was typed; the spec file carries the integers it means."""
    saved = client.post("/api/spec", json={"spec": spec_payload(), "filename": "range.json"})
    assert saved.status_code == 200
    spec = json.loads(Path(saved.json()["path"]).read_text())
    assert spec["coordination"] == [4, 5, 6] and spec["ligands_per_metal"] == [1, 2]


def test_an_unfinished_spec_gets_a_reason_rather_than_a_500(client):
    """The page asks at every edit, so a spec mid-typing is the normal input."""
    body = client.post("/api/estimate", json={"spec": spec_payload(molecules=[])}).json()
    assert body["ok"] is False and "molecule" in body["reason"]
    bad = client.post("/api/estimate", json={"spec": spec_payload(coordination="3~1")}).json()
    assert bad["ok"] is False and "counts down" in bad["reason"]


# ── a run whose process is gone must stop reading as ongoing ─────────────────

def _strand_a_run(db, store):
    """Leave behind exactly what a killed shell leaves: a `pending` run holding a task
    claimed by a pid that no longer exists, and no ending written anywhere."""
    import os
    import socket
    import subprocess
    import sys

    from mofsbu.registry.jobs import add_task, claim_task, create_run
    from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec

    if os.name != "posix":
        pytest.skip("process liveness is POSIX-only; the heartbeat covers the rest")
    gone = subprocess.Popen([sys.executable, "-c", ""])
    gone.wait()

    spec = BuildSpec(molecules=(MoleculeSpec("THQ", THQ, 1),),
                     metals=(MetalSpec("Ni", 2, "hs"),), coordination=(6,),
                     ligands_per_metal=(1,), binding=("chelate",))
    with Registry(db, BlobStore(store)) as reg:
        reg.migrate("test")
        run_id = create_run(reg, spec)
        add_task(reg, run_id, "place", {"i": 0})
        add_task(reg, run_id, "place", {"i": 1})
        claim_task(reg, run_id, worker=f"{socket.gethostname()}:{gone.pid}")
        reg.conn.commit()
    return run_id


def test_a_run_killed_with_its_shell_no_longer_reads_as_ongoing(registry, active):
    """The reported symptom.  The run row still says `pending` — it was written by a
    process that is no longer able to correct it — and the page believed it."""
    db, store = registry
    run_id = _strand_a_run(db, store)

    # Starting a server is the clearest statement that the previous one is not running.
    client = TestClient(create_app(db, store, active))
    row = next(r for r in client.get("/api/runs").json()["runs"] if r["id"] == run_id)
    assert row["status"] == "interrupted"
    assert row["liveness"] is None, "a closed-out run has nothing left to be live about"

    detail = client.get(f"/api/runs/{run_id}").json()["run"]
    assert detail["status"] == "interrupted"
    # And the stranded task went back to the queue rather than staying claimed forever.
    counts = {r["status"]: r["n"] for r in [
        {"status": t["status"], "n": 1} for t in
        client.get(f"/api/runs/{run_id}/tasks").json()["tasks"]]}
    assert "claimed" not in counts


def test_a_worker_that_dies_while_the_server_is_up_is_reported_on_the_read_path(
        registry, active):
    """No restart to wait for: asking a pid whether it exists writes nothing, so the
    listing can tell the truth about a run it has not swept."""
    db, store = registry
    client = TestClient(create_app(db, store, active))     # sweeps; nothing to sweep yet
    run_id = _strand_a_run(db, store)                      # then the worker dies

    row = next(r for r in client.get("/api/runs").json()["runs"] if r["id"] == run_id)
    assert row["status"] == "pending", "the row itself cannot know"
    assert row["liveness"]["verdict"] == "interrupted"
    assert "gone" in row["liveness"]["reason"]


def test_resuming_an_interrupted_run_returns_its_tasks_to_the_queue(registry, active):
    db, store = registry
    run_id = _strand_a_run(db, store)
    client = TestClient(create_app(db, store, active))

    body = client.post(f"/api/runs/{run_id}/resume").json()
    assert body["was"] == "interrupted"
    assert body["pending"] == 2, "the claimed task is stranded, not in flight"

    with Registry(db, BlobStore(store)) as reg:
        assert reg.conn.execute("SELECT status FROM runs WHERE id=?",
                                (run_id,)).fetchone()["status"] == "running"


# ── the worker count on the page has to reach the work ───────────────────────

def test_a_run_submitted_from_the_page_uses_the_declared_pool(client, monkeypatch):
    """The reported bug, at its source.

    The page's worker count was written to the run row and applied to nothing: `submit_run`
    called `work` directly, which is one worker whatever the machine had been declared to
    be. Every run from the browser — the only way most people start one — was serial.
    """
    import threading

    from mofsbu import runner

    seen: dict = {}
    called = threading.Event()

    def fake_execute_run(reg, spec, run_id, *, workers=None):
        seen.update(run_id=run_id, workers=workers,
                    pool=runner.plan_workers(workers,
                                             relaxes=spec.run_mode != "construct"))
        called.set()
        return seen["pool"]

    monkeypatch.setattr(runner, "execute_run", fake_execute_run)
    monkeypatch.setattr(runner, "work", lambda *a, **k: pytest.fail(
        "the page went back to draining the queue with a single in-process worker"))
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda")
    compute = client.post("/api/compute",
                          json={"device": "cuda", "workers": 4}).json()
    # The hardware panel states the division, so "4 workers" cannot read as four
    # relaxations at once on one card.
    assert (compute["build_workers"], compute["relax_workers"]) == (3, 1)
    assert "relaxing" in compute["worker_note"]

    # A `construct` run has no relax queue, so all four build — the page reports the
    # division this run will use, not the one the machine would use for another.
    body = client.post("/api/runs", json={"spec": spec_payload()}).json()
    assert called.wait(10), "the run thread never started the work"
    assert seen["run_id"] == body["run_id"]
    pool = seen["pool"]
    assert (pool.total, pool.build, pool.relax) == (4, 4, 0)
    assert body["workers"] == 4 and body["build_workers"] == 4
    assert body["relax_workers"] == 0


def test_the_page_says_what_the_notation_and_the_ladder_mean(client):
    """Both are rendered from capabilities, so the page cannot describe them wrongly."""
    caps = client.get("/api/capabilities").json()
    assert "1~3" in caps["range_note"]
    assert "join" in caps["pathways_note"]


def test_the_co_ligand_count_is_offered_and_an_old_spec_reads_as_fill(client):
    """D26: the page renders the choice and its default from capabilities; a spec posted
    without the field (the page labels its specs v1, fixtures here v7) plans `fill`."""
    caps = client.get("/api/capabilities").json()
    assert caps["co_ligand_counts"] == ["fill", "range"]
    assert caps["co_ligand_counts_default"] == "range"
    assert caps["co_ligand_window_default"] == 2
    assert "never the bare metal" in caps["co_ligand_note"]

    def estimate(**over):
        return client.post("/api/estimate", json={"spec": spec_payload(**over)}).json()

    old = estimate()
    assert old["co_ligand_counts"] == "fill" and "co_ligand" not in old["by_kind"]
    new = estimate(co_ligand_counts="range", co_ligand_window=2)
    assert new["by_kind"]["co_ligand"] == 1 and new["tasks"] > old["tasks"]
    assert new["capped"] == []
    # The page labels its specs v1 (B5): the field it sends must survive the migration.
    page = estimate(spec_version=1, co_ligand_counts="range", co_ligand_window=2)
    assert page["co_ligand_counts"] == "range" and page["tasks"] == new["tasks"]
