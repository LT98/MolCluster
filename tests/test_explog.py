"""Experiment log: blank vs n/a vs unknown, problems, computed amounts, round trips."""
from __future__ import annotations

import pytest
import yaml

from mofsbu.explog.chem import molar_mass
from mofsbu.explog.cli import main
from mofsbu.explog.emit import render_entry
from mofsbu.explog.model import check_entry, eta_ul_per_mg, line_mmol, load_log
from mofsbu.explog.render import diff, index, show
from mofsbu.explog.schema import entry_schema

REAGENTS = {
    "CuOAc": {"name": "copper(II) acetate monohydrate", "kind": "metal_salt",
              "formula": "Cu(C2H3O2)2·H2O"},
    "ATF": {"name": "anthrarufin", "kind": "ligand", "formula": "C14H8O4"},
    "EDTA": {"name": "EDTA", "kind": "additive", "formula": "C10H16N2O8"},
    "EDA": {"name": "ethylenediamine", "kind": "additive", "formula": "C2H8N2"},
    "DMF": {"name": "DMF", "kind": "solvent", "formula": "C3H7NO"},
}


def _entry(**over) -> dict:
    entry = {
        "id": "2-01", "status": "done", "method": "ball_mill",
        "reagents": [
            {"code": "CuOAc", "role": "metal_source", "mass_mg": 199.65, "ratio": 2},
            {"code": "EDTA", "role": "additive", "mass_mg": 146.12, "ratio": 1},
            {"code": "ATF", "role": "ligand", "mass_mg": 120.11, "ratio": 1},
            {"code": "DMF", "role": "solvent", "volume_ml": 0.1},
        ],
        "conditions": {"time_h": 48, "liquid": "DMF", "liquid_ul": 93.1},
        "characterization": {"pxrd": {"phase": "A"}},
    }
    entry.update(over)
    return entry


@pytest.fixture
def root(tmp_path):
    (tmp_path / "entries").mkdir()
    (tmp_path / "raw").mkdir()
    (tmp_path / "reagents.yaml").write_text(yaml.safe_dump(REAGENTS))
    (tmp_path / "phases.yaml").write_text(yaml.safe_dump({"A": {"description": "x"},
                                                          "B": {"description": "y"}}))
    return tmp_path


def _write(root, entry: dict) -> None:
    (root / "entries" / f"{entry['id']}.yaml").write_text(render_entry(entry))


def test_molar_mass_reads_groups_and_hydrates():
    assert molar_mass("H2O") == pytest.approx(18.015, abs=1e-2)
    assert molar_mass("Cu(C2H3O2)2·H2O") == pytest.approx(199.65, abs=0.05)
    assert molar_mass("NiCl2·6H2O") == pytest.approx(237.69, abs=0.05)
    for bad in ("Xq2", "Cu(OH", "2", ""):
        with pytest.raises(ValueError):
            molar_mass(bad)


def test_blank_is_pending_and_sentinels_are_not(root):
    entry = _entry()
    entry["conditions"]["mill"] = None
    entry["conditions"]["pauses"] = "n/a"
    entry["conditions"]["jar_material"] = "unknown"
    _write(root, entry)
    found = check_entry(load_log(root), "2-01")
    assert "conditions.mill" in found.pending
    assert not any(p.startswith(("conditions.pauses", "conditions.jar_material"))
                   for p in found.pending)
    assert "cells (nothing recorded)" in found.pending
    assert found.problems == []


def test_planned_entry_is_not_asked_for_outcomes(root):
    _write(root, _entry(status="planned"))
    found = check_entry(load_log(root), "2-01")
    assert not any(p.startswith(("product", "characterization", "cells")) for p in found.pending)


def test_solvent_needs_no_ratio_but_every_line_needs_an_amount(root):
    entry = _entry()
    entry["reagents"][3] = {"code": "DMF", "role": "solvent", "note": "few drops"}
    _write(root, entry)
    found = check_entry(load_log(root), "2-01")
    assert "reagents[DMF].ratio" not in found.pending
    assert "reagents[DMF].amount (mass_mg, volume_ml or mmol)" in found.pending


def test_problems_are_reported_never_filled(root):
    entry = _entry(id="2-02", surprise=1)
    entry["reagents"][1]["code"] = "EDTX"
    entry["conditions"]["time_h"] = "two days"
    entry["characterization"]["pxrd"] = {"phase": "C", "file": "raw/missing.xy"}
    (root / "entries" / "2-01.yaml").write_text(render_entry(entry))
    found = check_entry(load_log(root), "2-01")
    text = "\n".join(found.problems)
    for expected in ("not match the file name", "surprise: not a field", "'EDTX' is not in",
                     "conditions.time_h: expected a number", "phase 'C' is not in",
                     "raw/missing.xy not found"):
        assert expected in text


def test_mmol_from_mass_and_ratio_check(root):
    _write(root, _entry())
    log = load_log(root)
    lines = log.entries["2-01"]["reagents"]
    assert line_mmol(log, lines[0]) == pytest.approx(1.0, rel=1e-3)
    assert check_entry(log, "2-01").problems == []
    lines[2]["mass_mg"] = 60.0                      # half the ATF the 2:1:1 ratio asks for
    problems = check_entry(log, "2-01").problems
    assert len(problems) == 1 and problems[0].startswith("reagents[ATF]: ratio says 0.5")


def test_eta_needs_every_solid_mass(root):
    entry = _entry()
    assert eta_ul_per_mg(entry) == pytest.approx(93.1 / (199.65 + 146.12 + 120.11))
    entry["reagents"][1]["mass_mg"] = None
    assert eta_ul_per_mg(entry) is None


def test_fmt_round_trips_values_and_keeps_unknown_keys(root):
    entry = _entry(surprise="keep me", notes="line one\nline two\n")
    _write(root, entry)
    reread = yaml.safe_load(render_entry(load_log(root).entries["2-01"]))
    assert reread["surprise"] == "keep me"
    assert reread["notes"] == "line one\nline two\n"
    assert reread["reagents"][0] == {**entry["reagents"][0], "volume_ml": None, "mmol": None,
                                     "note": None}
    assert reread["conditions"]["liquid_ul"] == 93.1


def test_new_like_copies_recipe_and_blanks_outcome(root):
    entry = _entry(date=None, product={"appearance": "dark"})
    _write(root, entry)
    assert main(["--root", str(root), "new", "2-04", "--like", "2-01"]) == 0
    assert main(["--root", str(root), "new", "2-04", "--like", "2-01"]) == 2   # no overwrite
    new = load_log(root).entries["2-04"]
    assert new["derived_from"] == "2-01" and new["status"] == "planned"
    assert new["reagents"][1]["code"] == "EDTA"
    assert new["product"]["appearance"] is None
    assert new["characterization"]["pxrd"]["phase"] is None
    assert (root / ".schema" / "entry.schema.json").exists()


def test_views_render(root):
    _write(root, _entry())
    other = _entry(id="2-02")
    other["reagents"][1] = {"code": "EDA", "role": "additive", "mass_mg": 30.05, "ratio": 1}
    other["characterization"]["pxrd"]["phase"] = "B"
    _write(root, other)
    log = load_log(root)
    changed = diff(log, "2-01", "2-02")
    assert "| reagents[additive].code | EDTA | EDA |" in changed
    assert "| characterization.pxrd.phase | A | B |" in changed
    assert "reagents[ligand]" not in changed
    card = show(log, "2-01")
    assert "η (computed)" in card and "(calc)" in card
    overview = index(log)
    assert "| A | x |  | 2-01 |" in overview and "| B | y |  | 2-02 |" in overview


def test_schema_accepts_an_entry_and_rejects_a_typo(root):
    jsonschema = pytest.importorskip("jsonschema")
    _write(root, _entry())
    log = load_log(root)
    schema = entry_schema(log)
    jsonschema.validate(log.entries["2-01"], schema)
    bad = dict(log.entries["2-01"], conditions={"time_hours": 48})
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)
