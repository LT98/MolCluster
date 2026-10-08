"""`RouteSpec.vary`: a construction choice that decides an outcome is enumerated, not pinned."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from mofsbu.registry import BlobStore, Registry
from mofsbu.runner import enumerate_plan, run
from mofsbu.spec import BuildSpec, RouteSpec

ENUM = Path(__file__).resolve().parents[1] / "data" / "reference" / "spec_m6_s5_enumerate.json"
STEPS = ({"op": "metal", "as": "cu", "symbol": "Cu", "cn": "$cn", "geometry": "$g"},
         {"op": "store", "as": "s", "block": "cu"})


def test_every_combination_becomes_one_named_route():
    r = RouteSpec(name="r", steps=STEPS, vary={"cn": [4, 6], "g": ["square_planar", "octahedral"]})
    names = [r.resolved(v).name for v in r.variants()]
    assert len(names) == 4 == len(set(names))
    assert names[0] == 'r[cn=4,g="square_planar"]'


def test_only_exact_placeholders_are_substituted():
    r = RouteSpec(name="r", steps=STEPS + ({"op": "store", "as": "t", "block": "$cn-ish"},),
                  vary={"cn": [4], "g": ["square_planar"]})
    steps = r.resolved({"cn": 4, "g": "square_planar"}).steps
    assert steps[0]["cn"] == 4 and steps[2]["block"] == "$cn-ish"


@pytest.mark.parametrize("vary, says", [
    ({"cn": []}, "lists no values"),
    ({"cn": [4], "g": ["x"], "unused": [1]}, "used by no step"),
])
def test_a_bad_vary_is_refused(vary, says):
    with pytest.raises(ValueError, match=says):
        RouteSpec(name="r", steps=STEPS, vary=vary)


def test_the_enumeration_spec_plans_every_variant_of_both_routes():
    tasks = [t for t in enumerate_plan(BuildSpec.load(ENUM)).tasks if t.kind == "route"]
    assert len(tasks) == 32
    assert len({t.payload["name"] for t in tasks}) == 32


def test_both_s5_routes_reach_one_node_under_the_same_choice(tmp_path):
    """The S5 finding as a test: syn, upright, same face, second Cu syn — both routes build
    the bis-bridged node QC-clean; anti binding builds neither."""
    d = json.loads(ENUM.read_text())
    d["routes"][0]["vary"] = {"lp": [0, 1], "face": [0], "az": [23], "j": [{"toward": "r1.m0"}]}
    d["routes"][1]["vary"] = {"lp": [0, 1], "face2": [0], "az": [23], "j": [{"toward": "r3.m0"}]}
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as reg:
        reg.migrate("test")
        run(reg, BuildSpec.from_dict(d), workers=1)
        out = {}
        for row in reg.conn.execute("SELECT status, structure_id, json_extract(payload_json, "
                                    "'$.name') AS name FROM tasks WHERE kind='route'"):
            out[row["name"]] = (row["status"], row["structure_id"])
    a_syn = next(v for k, v in out.items() if k.startswith("A[") and "lp=1" in k)
    a_anti = next(v for k, v in out.items() if k.startswith("A[") and "lp=0" in k)
    b_syn = next(v for k, v in out.items() if k.startswith("B[") and "lp=1" in k)
    assert a_syn[0] == "done" and b_syn[0] == "done"
    assert a_anti[0] == "failed"                     # anti binding fails Route A too
    assert a_syn[1] is not None and b_syn[1] is not None
