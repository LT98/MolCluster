"""Experiment log: blank vs n/a vs unknown, problems, computed amounts, and the page's API."""
from __future__ import annotations

import json

import pytest
from fastapi.testclient import TestClient

from mofsbu.explog.chem import molar_mass
from mofsbu.explog.model import check_entry, computed, copy_recipe, load_log
from mofsbu.ui.app import create_app

REAGENTS = {
    "CuOAc": {"name": "copper(II) acetate monohydrate", "kind": "metal_salt",
              "formula": "Cu(C2H3O2)2·H2O"},
    "ATF": {"name": "anthrarufin", "kind": "ligand", "formula": "C14H8O4"},
    "EDTA": {"name": "EDTA", "kind": "additive", "formula": "C10H16N2O8"},
    "EDA": {"name": "ethylenediamine", "kind": "additive", "formula": "C2H8N2"},
    "NH4OH": {"name": "ammonia solution", "kind": "additive", "formula": "NH3",
              "density_g_ml": 0.9, "wt_pct": 28},
    "HCl": {"name": "hydrochloric acid", "kind": "additive", "formula": "HCl",
            "molarity_m": 1.0},
    "DMF": {"name": "DMF", "kind": "solvent", "formula": "C3H7NO"},
    "H2O": {"name": "water", "kind": "solvent", "formula": "H2O"},
}


def _milled(**over) -> dict:
    entry = {
        "id": "2-01", "status": "done", "method": "ball_mill",
        "preparations": [{
            "label": "jar charge", "form": "solid_mixture",
            "components": [
                {"code": "CuOAc", "role": "metal_source", "mass_mg": 199.65, "ratio": 2},
                {"code": "EDTA", "role": "additive", "mass_mg": 146.12, "ratio": 1},
                {"code": "ATF", "role": "ligand", "mass_mg": 120.11, "ratio": 1},
                {"code": "DMF", "role": "solvent", "volume_ml": 0.1},
            ]}],
        "steps": [{"action": "mill", "what": "jar sealed and milled", "duration_h": 48}],
        "apparatus": {"liquid": "DMF", "liquid_ul": 93.1},
        "characterization": {"pxrd": {"phase": "A"}},
    }
    entry.update(over)
    return entry


def _solution() -> dict:
    return {
        "id": "1-01", "status": "done", "method": "solution",
        "preparations": [
            {"label": "A", "form": "solution", "components": [
                {"code": "CuOAc", "role": "metal_source", "mass_mg": 99.83, "ratio": 1},
                {"code": "H2O", "role": "solvent", "volume_ml": 50}]},
            {"label": "B", "form": "suspension", "components": [
                {"code": "ATF", "role": "ligand", "mass_mg": 120.11, "ratio": 1},
                {"code": "NH4OH", "role": "additive", "volume_ml": 1.0, "ratio": 30},
                {"code": "HCl", "role": "additive", "volume_ml": 2.0, "ratio": 4},
                {"code": "DMF", "role": "solvent", "volume_ml": 50}]},
        ],
    }


@pytest.fixture
def root(tmp_path):
    folder = tmp_path / "experiments"
    (folder / "entries").mkdir(parents=True)
    (folder / "reagents.json").write_text(json.dumps(REAGENTS))
    (folder / "phases.json").write_text(json.dumps({"A": {"description": "x"},
                                                    "B": {"description": "y"}}))
    return folder


@pytest.fixture
def client(tmp_path, root):
    return TestClient(create_app(tmp_path / "r.db", tmp_path / "store", explog_root=root))


def _write(root, entry: dict) -> None:
    (root / "entries" / f"{entry['id']}.json").write_text(json.dumps(entry))


def test_molar_mass_reads_groups_and_hydrates():
    assert molar_mass("H2O") == pytest.approx(18.015, abs=1e-2)
    assert molar_mass("Cu(C2H3O2)2·H2O") == pytest.approx(199.65, abs=0.05)
    assert molar_mass("NiCl2·6H2O") == pytest.approx(237.69, abs=0.05)
    for bad in ("Xq2", "Cu(OH", "2", ""):
        with pytest.raises(ValueError):
            molar_mass(bad)


def test_blank_is_pending_and_sentinels_are_not(root):
    entry = _milled()
    entry["apparatus"].update(mill=None, pauses="n/a", jar_material="unknown")
    found = check_entry(load_log(root), entry, "2-01")
    assert "apparatus.mill" in found.pending
    assert not any(p.startswith(("apparatus.pauses", "apparatus.jar_material"))
                   for p in found.pending)
    assert "cells (nothing recorded)" in found.pending
    assert "steps[mill].temperature_c" in found.pending
    assert found.problems == []


def test_planned_entry_is_not_asked_for_outcomes(root):
    found = check_entry(load_log(root), _milled(status="planned"), "2-01")
    assert not any(p.startswith(("product", "characterization", "cells")) for p in found.pending)


def test_solvent_needs_no_ratio_but_every_component_needs_an_amount(root):
    entry = _milled()
    entry["preparations"][0]["components"][3] = {"code": "DMF", "role": "solvent",
                                                 "note": "few drops"}
    found = check_entry(load_log(root), entry, "2-01")
    base = "preparations[jar charge].components[DMF]"
    assert f"{base}.ratio" not in found.pending
    assert f"{base}.amount (mass_mg, volume_ml or mmol)" in found.pending


def test_problems_are_reported_never_filled(root):
    entry = _milled(id="2-02", surprise=1)
    comps = entry["preparations"][0]["components"]
    comps[1]["code"] = "EDTX"
    entry["steps"][0]["duration_h"] = "two days"
    entry["characterization"]["pxrd"] = {"phase": "C", "file": "raw/missing.xy"}
    entry["date"] = "2026-13-01"
    text = "\n".join(check_entry(load_log(root), entry, "2-01").problems)
    for expected in ("not match the file name", "surprise: not a field", "'EDTX' is not in",
                     "steps[mill].duration_h: expected a number", "phase 'C' is not in",
                     "raw/missing.xy not found", "date: expected a date"):
        assert expected in text


def test_amounts_concentrations_and_ratio_check(root):
    log = load_log(root)
    calc = computed(log, _solution())
    mmol = {r["code"]: r["mmol"] for r in calc["components"]}
    assert mmol["CuOAc"] == pytest.approx(0.5, rel=1e-3)
    assert mmol["NH4OH"] == pytest.approx(1.0 * 0.9 * 1000 * 0.28 / 17.031, rel=1e-3)
    assert mmol["HCl"] == pytest.approx(2.0)
    prep_a, prep_b = calc["preparations"]
    assert prep_a["conc_m"]["CuOAc"] == pytest.approx(0.01, rel=1e-3)       # 0.5 mmol / 50 mL
    assert prep_b["solvent_ml"] == 50
    assert calc["mixture"]["solutes"]["CuOAc"]["conc_m"] == pytest.approx(0.005, rel=1e-3)
    entry = _solution()
    assert check_entry(log, entry, "1-01").problems == []
    entry["preparations"][1]["components"][0]["mass_mg"] = 60.0           # half the ligand
    problems = check_entry(log, entry, "1-01").problems
    assert len(problems) == 1 and "components[ATF]: ratio says 1" in problems[0]


def test_eta_needs_every_solid_mass(root):
    log = load_log(root)
    assert computed(log, _milled())["eta_ul_per_mg"] == pytest.approx(
        93.1 / (199.65 + 146.12 + 120.11))
    entry = _milled()
    entry["preparations"][0]["components"][1]["mass_mg"] = None
    assert computed(log, entry)["eta_ul_per_mg"] is None


def test_copy_recipe_keeps_the_recipe_and_drops_the_outcome():
    new = copy_recipe(_milled(date="2026-10-01", product={"appearance": "dark"}), "2-04", "2-01")
    assert new["id"] == "2-04" and new["derived_from"] == "2-01" and new["status"] == "planned"
    assert new["preparations"][0]["components"][1]["code"] == "EDTA"
    assert "product" not in new and "characterization" not in new and "date" not in new


# ── the page's API ───────────────────────────────────────────────────────────

def test_page_is_served_and_in_the_tab_strip(client):
    assert client.get("/experiments").status_code == 200
    assert '"/experiments"' in client.get("/static/chrome.js").text
    spec = client.get("/api/explog/spec").json()
    assert "ball_mill" in spec["apparatus"] and spec["sentinels"] == ["n/a", "unknown"]


def test_create_save_and_conflict(client, root):
    assert client.post("/api/explog/entries", json={"id": "../x", "method": "solution"}
                       ).status_code == 400
    r = client.post("/api/explog/entries", json={"id": "4-01", "method": "solution"})
    assert r.status_code == 200, r.text
    assert client.post("/api/explog/entries", json={"id": "4-01", "method": "solution"}
                       ).status_code == 409
    got = client.get("/api/explog/entries/4-01").json()
    entry = dict(got["entry"], title="  first try ", operator="", workup={"drying": {}})
    r = client.put("/api/explog/entries/4-01", json={"entry": entry, "version": got["version"]})
    assert r.status_code == 200, r.text
    stored = json.loads((root / "entries" / "4-01.json").read_text())
    assert stored["title"] == "first try"
    assert "operator" not in stored and "workup" not in stored              # blanks pruned
    assert "| 4-01 |" in (root / "INDEX.md").read_text()
    stale = client.put("/api/explog/entries/4-01", json={"entry": entry,
                                                         "version": got["version"]})
    assert stale.status_code == 409


def test_check_does_not_save_and_diff_pairs_by_role(client, root):
    _write(root, _milled())
    other = _milled(id="2-02")
    other["preparations"][0]["components"][1] = {
        "code": "EDA", "role": "additive", "mass_mg": 30.05, "ratio": 1}
    other["characterization"]["pxrd"]["phase"] = "B"
    _write(root, other)
    r = client.post("/api/explog/check", json={"entry": _milled(), "id": "2-01"}).json()
    assert r["findings"]["problems"] == [] and r["computed"]["eta_ul_per_mg"]
    rows = {x["field"]: (x["a"], x["b"]) for x in
            client.get("/api/explog/diff", params={"a": "2-01", "b": "2-02"}).json()["rows"]}
    assert rows["preparations[1].components[additive].code"] == ("EDTA", "EDA")
    assert rows["characterization.pxrd.phase"] == ("A", "B")
    assert not any("components[ligand]" in f for f in rows)
    listing = client.get("/api/explog/entries").json()["entries"]
    assert {e["id"]: e["phase"] for e in listing} == {"2-01": "A", "2-02": "B"}
    assert "η" in client.get("/api/explog/entries/2-01/card").text


def test_delete_moves_aside(client, root):
    _write(root, _milled())
    r = client.delete("/api/explog/entries/2-01")
    assert r.status_code == 200 and r.json()["moved_to"].startswith("entries/.deleted/")
    assert not (root / "entries" / "2-01.json").exists()
    assert list((root / "entries" / ".deleted").glob("2-01.*.json"))
    assert client.delete("/api/explog/entries/2-01").status_code == 404


def test_reagent_table_round_trip(client):
    got = client.get("/api/explog/reagents").json()
    assert got["molar_mass"]["CuOAc"] == pytest.approx(199.65, abs=0.05)
    assert "cas" in got["findings"]["CuOAc"]["pending"]
    table = dict(got["table"], NiCl={"name": "nickel chloride", "formula": "NiCl2·6H2O",
                                     "kind": "metal_salt"})
    r = client.put("/api/explog/reagents", json={"table": table, "version": got["version"]})
    assert r.status_code == 200, r.text
    assert r.json()["molar_mass"]["NiCl"] == pytest.approx(237.69, abs=0.05)
    assert client.put("/api/explog/reagents", json={"table": {"a b": {}},
                                                    "version": r.json()["version"]}
                      ).status_code == 400
