"""Explicit construction routes in a spec (v9): load, refuse, build, replay, relax.

The S5 spec in `data/reference/` is the fixture: it is the file the M6 route test is run
from, so a change that breaks it breaks here first.  It builds from the ground up — every
reagent placed and stored by the pipeline's own paths, every join on blocks read back from
the registry, every vertex a formate takes first freed by releasing its water.
"""
from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from mofsbu._types import Fidelity
from mofsbu.energy import relax as relax_mod
from mofsbu.energy.backends import NullBackend
from mofsbu.energy.reference import check_balance
from mofsbu.registry import BlobStore, Registry, geometry_xyz, get_graph
from mofsbu.runner import enumerate_plan, existing_relaxation, read_xyz, run
from mofsbu.spec import SPEC_VERSION, BuildSpec, RouteSpec

S5 = Path(__file__).resolve().parents[1] / "data" / "reference" / "spec_m6_s5_omol.json"


def s5(**kw) -> BuildSpec:
    return dataclasses.replace(BuildSpec.load(S5), **kw)


@pytest.fixture()
def null_ml(monkeypatch):
    monkeypatch.setattr(relax_mod, "backend_for",
                        lambda fidelity, **kw: NullBackend(fidelity=fidelity))


def _run(tmp: Path, spec: BuildSpec, name: str = "r") -> tuple[Registry, dict]:
    reg = Registry(tmp / f"{name}.db", BlobStore(tmp / f"{name}_store"))
    reg.migrate("test")
    return reg, run(reg, spec, workers=1)


def _species(reg: Registry) -> dict[str, dict]:
    """Species by `route/name` — route B re-stores R2 and its waters under the same names."""
    out = {}
    for row in reg.conn.execute("SELECT payload_json, detail_json FROM tasks WHERE "
                                "kind='route' AND status='done' ORDER BY id"):
        route = json.loads(row["payload_json"])["name"]
        for sp in json.loads(row["detail_json"])["species"]:
            out[f"{route}/{sp['name']}"] = sp
    return out


def _xyz(reg, gid):
    sym, xyz = read_xyz(geometry_xyz(reg, gid))
    return sym, np.asarray(xyz, dtype=float)


def _edge(reg, product_sid):
    rid = reg.conn.execute("SELECT id FROM reactions WHERE product_structure_id=? AND "
                           "EXISTS (SELECT 1 FROM reaction_reagents rr WHERE rr.reaction_id="
                           "reactions.id) ORDER BY id DESC LIMIT 1", (product_sid,)).fetchone()[0]
    rows = reg.conn.execute("SELECT structure_id, stoich, role FROM reaction_reagents WHERE "
                            "reaction_id=?", (rid,)).fetchall()
    return rid, {(r["structure_id"], r["role"]): r["stoich"] for r in rows}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    reg, summary = _run(tmp_path_factory.mktemp("s5"), s5(run_mode="construct"))
    yield reg, summary, _species(reg)
    reg.close()


# ── loading ──────────────────────────────────────────────────────────────────

def test_the_s5_spec_loads_at_the_current_version():
    spec = BuildSpec.load(S5)
    assert spec.to_dict()["spec_version"] == SPEC_VERSION == 9
    assert [r.name for r in spec.routes] == ["A_sequential", "B_two_on_one"]
    assert (spec.relax_fmax, spec.relax_steps) == (0.05, 2000)


def test_a_v8_spec_migrates_to_no_routes_and_the_relax_defaults():
    d = json.loads(S5.read_text())
    for k in ("routes", "relax_fmax", "relax_steps"):
        d.pop(k)
    spec = BuildSpec.from_dict({**d, "spec_version": 8})
    assert spec.routes == () and spec.relax_fmax is None and spec.relax_steps is None


@pytest.mark.parametrize("step, says", [
    ({"op": "teleport", "as": "x"}, "unknown op"),
    ({"op": "metal", "as": "cu1", "symbol": "Cu"}, r"missing \['cn'\]"),
    ({"op": "release", "as": "w"}, r"missing \['ligand'\]"),
])
def test_a_malformed_step_is_refused_at_load(step, says):
    with pytest.raises(ValueError, match=says):
        RouteSpec(name="r", steps=(step, {"op": "store", "as": "s", "block": "x"}))


def test_a_route_that_stores_nothing_is_refused():
    with pytest.raises(ValueError, match="stores nothing"):
        RouteSpec(name="r", steps=({"op": "load", "as": "x", "species": "y"},))


def test_a_sphere_naming_a_molecule_the_spec_lacks_is_refused():
    d = json.loads(S5.read_text())
    d["molecules"] = [{"name": "acetate", "smiles": "CC(=O)[O-]"}]
    with pytest.raises(ValueError, match="'formate' is not in this spec's molecules"):
        BuildSpec.from_dict(d)


def test_a_load_names_exactly_one_source():
    d = json.loads(S5.read_text())
    d["routes"][0]["steps"][3] = {"op": "load", "as": "r1", "species": "R1", "block_id": "x"}
    with pytest.raises(ValueError, match="exactly one of"):
        BuildSpec.from_dict(d)


def test_each_route_is_planned_as_one_task():
    kinds = [t.kind for t in enumerate_plan(BuildSpec.load(S5)).tasks]
    assert kinds.count("route") == 2


# ── building from the ground up ──────────────────────────────────────────────

def test_both_routes_build(built):
    _, summary, sp = built
    assert summary["counts"] == {"done": 3}             # free formate + two routes
    assert {"A_sequential/proto_dimer", "A_sequential/product",
            "B_two_on_one/two_on_one"} <= set(sp)


def test_every_step_consumes_stored_reagents_and_balances(built):
    """The edges name registry entries — the reagents loaded, the waters released."""
    reg, _, sp = built
    a, b = (lambda n: sp[f"A_sequential/{n}"]["structure_id"]), \
        (lambda n: sp[f"B_two_on_one/{n}"]["structure_id"])
    water = a("water_1")
    cases = {
        a("proto_dimer"): {(a("R1"), "reagent"): 1, (a("R2"), "reagent"): 1,
                           (water, "leaving"): 1},
        a("product"): {(a("proto_dimer"), "reagent"): 1, (a("formate"), "reagent"): 1,
                       (water, "leaving"): 2},
        b("two_on_one"): {(b("R3"), "reagent"): 1, (b("R2"), "reagent"): 1,
                          (water, "leaving"): 2},
    }
    for product, expected in cases.items():
        rid, rows = _edge(reg, product)
        assert rows == expected
        assert check_balance(reg, rid).balanced


def test_the_built_species_are_qc_clean_and_have_the_expected_bridges(built):
    reg, _, sp = built
    for name, charge, mu2 in (("A_sequential/proto_dimer", 3, 1),
                              ("A_sequential/product", 2, 2),
                              ("B_two_on_one/two_on_one", 2, 1)):
        g = get_graph(reg, sp[name]["structure_id"])
        assert g.net_charge() == charge, name
        assert sum(g.fragment_bridge_class(f).value == "mu2"
                   for f in g.ligand_fragments()) == mu2, name
        qc = json.loads(reg.conn.execute("SELECT qc_json FROM geometries WHERE id=?",
                                         (sp[name]["geometry_id"],)).fetchone()[0])
        assert qc["ok"], (name, qc)


def test_both_routes_end_at_one_composition(built):
    reg, _, sp = built
    l0 = [get_graph(reg, sp[n]["structure_id"]).element_counts()
          for n in ("A_sequential/product", "B_two_on_one/two_on_one")]
    assert l0[0] == l0[1]


def test_a_perturbed_start_moves_only_what_it_names_by_what_it_says(built):
    reg, _, sp = built
    base, pulled = sp["A_sequential/proto_dimer"], sp["A_sequential/proto_dimer_pull_1.0"]
    shift = _xyz(reg, pulled["geometry_id"])[1] - _xyz(reg, base["geometry_id"])[1]
    moved = sorted(set(base["rows"]["r2"]))
    assert np.allclose(np.linalg.norm(shift[moved], axis=1), 1.0, atol=1e-5)
    assert np.allclose(np.delete(shift, moved, axis=0), 0.0, atol=1e-6)
    assert pulled["structure_id"] == base["structure_id"]


def test_a_route_replays_bit_identically(built, tmp_path):
    reg_a, _, sa = built
    reg_b, _ = _run(tmp_path, s5(run_mode="construct"), "b")
    sb = _species(reg_b)
    assert sa.keys() == sb.keys()
    for name in sa:
        assert _xyz(reg_a, sa[name]["geometry_id"])[1].tolist() == \
            _xyz(reg_b, sb[name]["geometry_id"])[1].tolist(), name


# ── relaxing ─────────────────────────────────────────────────────────────────

def test_every_stored_geometry_is_relaxed_with_the_specs_convergence(tmp_path, null_ml):
    reg, _ = _run(tmp_path, s5())
    payloads = {json.loads(r["payload_json"])["geometry_id"]: json.loads(r["payload_json"])
                for r in reg.conn.execute("SELECT payload_json FROM tasks WHERE kind='relax'")}
    for sp in _species(reg).values():
        p = payloads[sp["geometry_id"]]
        assert (p["fmax"], p["steps"]) == (0.05, 2000)


def test_a_spec_without_relax_settings_queues_the_payload_it_always_did(tmp_path, null_ml):
    reg, _ = _run(tmp_path, s5(relax_fmax=None, relax_steps=None))
    for r in reg.conn.execute("SELECT payload_json FROM tasks WHERE kind='relax'"):
        assert not {"fmax", "steps"} & set(json.loads(r["payload_json"]))


def test_a_looser_relaxation_is_not_reused_for_a_tighter_request(tmp_path, null_ml):
    reg, _ = _run(tmp_path, s5())
    row = reg.conn.execute("SELECT id, structure_id, relaxed_from FROM geometries "
                           "WHERE relaxed_from IS NOT NULL LIMIT 1").fetchone()
    reg.conn.execute("UPDATE tasks SET detail_json=json_set(detail_json, '$.relax.fmax', 0.15)"
                     " WHERE kind='relax' AND geometry_id=?", (row["id"],))
    args = (reg, row["relaxed_from"], Fidelity.ML, row["structure_id"])
    assert existing_relaxation(*args, fmax=0.05) is None
    assert existing_relaxation(*args, fmax=0.20) == row["id"]
    assert existing_relaxation(*args) == row["id"]
