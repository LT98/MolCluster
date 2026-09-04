"""Tests for the M3.5 read-only viewer (``mofsbu.ui``).

Every test runs against a freshly seeded *demo* registry in a temp directory, so the
suite never touches ``data/`` and never needs the real registry writer.  The app is
exercised through ``TestClient`` — no server is started.
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

from mofsbu.registry.store import BlobStore  # noqa: E402
from mofsbu.ui.app import SORTABLE, create_app  # noqa: E402
from seed_demo_registry import seed  # noqa: E402


@pytest.fixture(scope="module")
def registry(tmp_path_factory) -> tuple[Path, Path]:
    d = tmp_path_factory.mktemp("demo_registry")
    db, store = d / "demo_registry.db", d / "store"
    seed(db, store, n=40)
    return db, store


@pytest.fixture(scope="module")
def client(registry) -> TestClient:
    return TestClient(create_app(*registry))


def total(client: TestClient, **params) -> int:
    r = client.get("/api/structures", params=params)
    assert r.status_code == 200, r.text
    return r.json()["total"]


# ── shape ────────────────────────────────────────────────────────────────────

def test_index_serves_the_page(client):
    r = client.get("/")
    assert r.status_code == 200
    assert "3Dmol" in r.text and "mofsbu" in r.text


def test_list_shape(client):
    body = client.get("/api/structures", params={"limit": 5}).json()
    assert body["total"] == 40
    assert len(body["rows"]) == 5
    row = body["rows"][0]
    for key in ("id", "display_label", "formula", "metals", "net_charge", "multiplicity",
                "max_bridge_class", "best_fidelity", "best_energy", "n_geometries",
                "n_incoming_routes", "l1_graph_hash"):
        assert key in row


def test_every_row_carries_a_derived_label(client):
    rows = client.get("/api/structures", params={"limit": 1000}).json()["rows"]
    assert all(r["display_label"] for r in rows)   # derived, so never empty
    assert all(r["formula"] for r in rows)          # the fallback is always available


def test_detail_shape(client):
    body = client.get("/api/structures/4").json()
    assert body["structure"]["id"] == 4
    assert isinstance(body["routes"], list)
    for g in body["geometries"]:
        for key in ("id", "fidelity", "energy", "converged", "method", "n_atoms",
                    "is_best", "blob_present"):
            assert key in g


def test_detail_404(client):
    assert client.get("/api/structures/999999").status_code == 404


def test_xyz_round_trip(client, registry):
    db, store_root = registry
    con = sqlite3.connect(db)
    store = BlobStore(store_root)
    gid, digest = next((g, h) for g, h in
                       con.execute("SELECT id, coords_hash FROM geometries") if store.has(h))
    r = client.get(f"/api/geometries/{gid}/xyz")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/plain")
    lines = r.text.splitlines()
    assert int(lines[0]) == len(lines) - 2          # a real, parseable .xyz
    assert r.text == store.get_text(digest)


def test_filters_endpoint(client):
    f = client.get("/api/filters").json()
    assert {"Cu", "Zn", "Fe"} <= set(f["metals"])
    assert set(f["bridge_classes"]) <= {"none", "terminal", "mu2", "mu3", "muN"}
    assert f["fidelities"] and f["n_structures"] == 40
    assert f["charge"]["min"] <= f["charge"]["max"]
    assert set(f["sortable"]) == set(SORTABLE)


def test_meta(client):
    assert client.get("/api/meta").json()["read_only"] is True


# ── each filter axis actually narrows ────────────────────────────────────────

def test_axis_composition_narrows(client):
    everything = total(client)
    cu = total(client, metal="Cu")
    zn = total(client, metal="Zn")
    assert 0 < cu < everything
    assert cu + zn == total(client, metal=["Cu", "Zn"])          # repeatable → OR
    # free text is structural: whole formula, L0, or a FRAGMENT formula.
    # "paddlewheel" exists only in labels, so it must find nothing.
    assert total(client, q="paddlewheel") == 0
    assert 0 < total(client, q="CHO2") < everything
    assert total(client, n_metals_min=2) < everything
    assert total(client, n_metals_min=2) + total(client, n_metals_max=1) == everything
    assert total(client, charge_min=1) < everything
    assert 0 < total(client, multiplicity=1) < everything


def test_axis_connectivity_narrows(client):
    everything = total(client)
    mu3 = total(client, max_bridge_class="mu3")
    assert 0 < mu3 < everything
    assert mu3 <= total(client, max_bridge_class=["mu3", "mu2"])
    yes, no = total(client, has_metal_metal=True), total(client, has_metal_metal=False)
    assert 0 < yes < everything and yes + no == everything


def test_axis_energy_fidelity_narrows(client):
    everything = total(client)
    assert total(client, fidelity_min=4) < total(client, fidelity_min=1) < everything
    assert 0 < total(client, converged=True) < everything
    f = client.get("/api/filters").json()["energy"]
    midpoint = (f["min"] + f["max"]) / 2
    assert 0 < total(client, energy_min=midpoint) < everything
    assert total(client, energy_max=f["min"] - 1) == 0


def test_axis_provenance_narrows(client):
    everything = total(client)
    with_route = total(client, has_route=True)
    assert 0 < with_route < everything
    assert with_route + total(client, has_route=False) == everything


def test_sort_changes_order(client):
    asc = [r["id"] for r in client.get(
        "/api/structures", params={"sort": "best_energy", "order": "asc", "limit": 5}
    ).json()["rows"]]
    desc = [r["id"] for r in client.get(
        "/api/structures", params={"sort": "best_energy", "order": "desc", "limit": 5}
    ).json()["rows"]]
    assert asc != desc


def test_pagination(client):
    page1 = client.get("/api/structures", params={"limit": 10}).json()
    page2 = client.get("/api/structures", params={"limit": 10, "offset": 10}).json()
    assert page1["total"] == page2["total"] == 40
    assert {r["id"] for r in page1["rows"]}.isdisjoint({r["id"] for r in page2["rows"]})
    assert len(client.get("/api/structures", params={"limit": 99999}).json()["rows"]) == 40


# ── injection / validation ───────────────────────────────────────────────────

@pytest.mark.parametrize("bad_sort", [
    "id; DROP TABLE structures",
    "(SELECT 1)",
    "rowid",
    "l1_graph_hash) --",
])
def test_sort_rejects_anything_off_the_whitelist(client, bad_sort):
    r = client.get("/api/structures", params={"sort": bad_sort})
    assert r.status_code == 400
    assert total(client) == 40          # and nothing was executed


def test_order_and_bridge_class_are_validated(client):
    assert client.get("/api/structures", params={"order": "asc; DELETE"}).status_code == 400
    assert client.get("/api/structures", params={"max_bridge_class": "mu9"}).status_code == 400


def test_wildcards_in_free_text_are_literal(client):
    assert total(client, q="%") == 0     # a bare % must not match everything


# ── read-only ────────────────────────────────────────────────────────────────

def test_connection_is_read_only(registry):
    app = create_app(*registry)
    con = app.state.connect()
    with pytest.raises(sqlite3.OperationalError, match="readonly"):
        con.execute("INSERT INTO structure_tags (structure_id, tag) VALUES (1, 'x')")
    with pytest.raises(sqlite3.OperationalError):
        con.execute("DELETE FROM structures")
    with pytest.raises(sqlite3.OperationalError):
        con.execute("CREATE TABLE oops (a)")
    assert con.execute("SELECT COUNT(*) FROM v_structures").fetchone()[0] == 40
    con.close()


def test_app_exposes_no_write_helper():
    import mofsbu.ui.app as mod
    src = Path(mod.__file__).read_text()
    for verb in ("INSERT INTO", "UPDATE ", "DELETE FROM", "DROP "):
        assert verb not in src.upper().replace("READONLY", "")


# ── graceful degradation on RESERVED columns ─────────────────────────────────

def test_reserved_columns_are_null_but_do_not_break_the_list(client):
    """Site columns land in M4/M5; until then they are NULL and must not break anything."""
    rows = client.get("/api/structures", params={"limit": 5}).json()["rows"]
    for r in rows:
        assert r["n_open_sites"] is None and r["donor_types"] is None
        assert r["binding_modes"] is None
    # the filters still answer, they just select nothing yet
    assert client.get("/api/structures", params={"n_open_sites_min": 1}).status_code == 200
    assert total(client, n_open_sites_min=1) == 0
    assert total(client, donor_type="carboxylate_O") == 0


def test_provenance_depth_is_live_not_reserved(client):
    """Reaction edges are written from M3 onward, so `depth` carries real data.

    Only the pathway SCORES are reserved for M8.  A structure with no incoming route
    has a NULL min_depth, which must not be confused with depth 0.
    """
    everything = total(client)
    with_route = total(client, has_route=True)
    assert 0 < with_route < everything
    assert 0 < total(client, depth_max=3) <= with_route


def test_filters_reports_which_reserved_columns_are_inactive(client):
    f = client.get("/api/filters").json()
    inactive = {r["column"] for r in f["reserved_inactive"]}
    assert {"donor_types", "binding_modes", "n_open_sites", "choice_vector_json"} <= inactive
    assert f["donor_types"] == []       # empty vocabulary, not a crash


def test_structure_without_geometries_still_renders(client):
    body = client.get("/api/structures", params={"limit": 1000}).json()
    empties = [r for r in body["rows"] if r["n_geometries"] == 0]
    assert empties, "seed should include structures with no geometry"
    d = client.get(f"/api/structures/{empties[0]['id']}").json()
    assert d["geometries"] == []
    assert d["structure"]["best_energy"] is None


# ── missing blob ─────────────────────────────────────────────────────────────

def test_missing_blob_is_404_not_500(client, registry):
    db, store_root = registry
    store = BlobStore(store_root)
    con = sqlite3.connect(db)
    orphans = [g for g, h in con.execute("SELECT id, coords_hash FROM geometries")
               if not store.has(h)]
    assert orphans, "seed should include a geometry whose blob is absent"
    assert client.get(f"/api/geometries/{orphans[0]}/xyz").status_code == 404
    assert client.get("/api/geometries/999999/xyz").status_code == 404
    # the structure detail still loads and flags the geometry
    sid = con.execute("SELECT structure_id FROM geometries WHERE id = ?",
                      (orphans[0],)).fetchone()[0]
    d = client.get(f"/api/structures/{sid}").json()
    assert any(g["blob_present"] is False for g in d["geometries"])


# ── spec ─────────────────────────────────────────────────────────────────────

def test_spec_placeholder_when_no_choice_vector(client):
    spec = client.get("/api/structures/4/spec").json()
    assert spec["structure_id"] == 4
    assert "pre-M5" in spec["note"]
    assert spec["l1"] and len(spec["l1"]) == 64


def test_spec_404_for_unknown_structure(client):
    assert client.get("/api/structures/999999/spec").status_code == 404


# ── threading ────────────────────────────────────────────────────────────────

def test_a_connection_can_be_used_from_another_thread(registry):
    """The property that actually broke, tested directly.

    FastAPI runs a synchronous `yield` dependency in one threadpool thread and the
    endpoint that consumes it in another, so a per-request connection is opened on one
    thread and executed on a different one.  sqlite3 refuses that unless the connection
    was created with check_same_thread=False, and the viewer 500'd on every single
    request while this suite stayed green.

    Note this is deliberately NOT written as "fire concurrent requests at TestClient":
    that reuses a single worker thread, so it passes with or without the fix and guards
    nothing.  Opening here and executing there is the real condition.
    """
    import threading

    from mofsbu.ui.app import open_read_only

    db_path, _store = registry
    con = open_read_only(db_path)
    outcome: dict[str, object] = {}

    def use_it() -> None:
        try:
            outcome["count"] = con.execute("SELECT COUNT(*) FROM v_structures").fetchone()[0]
        except Exception as exc:                                    # noqa: BLE001
            outcome["error"] = exc

    thread = threading.Thread(target=use_it)
    thread.start()
    thread.join()
    con.close()

    assert "error" not in outcome, f"connection cannot cross threads: {outcome.get('error')}"
    assert outcome["count"] >= 0


def test_endpoints_answer_under_concurrent_requests(client):
    """A smoke test only — see the note above about why this cannot catch the thread bug."""
    from concurrent.futures import ThreadPoolExecutor

    urls = ["/api/structures?limit=10", "/api/filters", "/api/structures/1", "/api/meta"]
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(client.get, urls[i % len(urls)]) for i in range(32)]
        assert {f.result().status_code for f in futures} == {200}
