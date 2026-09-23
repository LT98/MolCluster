"""How many co-ligands a product carries, and the steps that add them (D26).

`co_ligand_counts="fill"` puts a co-ligand on every vertex the ligands leave; `"range"`, the
default for a v8 spec, also builds the counts down to `co_ligand_window` (default 2) fewer, in
the same polyhedron. With `pathways`, gaining one co-ligand is a ladder step, no rung leaves
more than `window` vertices empty, and the lowest co-ligand state inside the window is the
root — never the bare metal.

The claims under test: an older spec replays the run it planned before (pinned against the
planner as it stood before the field existed), the new default plans the lower co-ligand
rungs and the steps joining them into one ladder per metal/CN/geometry, and a ladder cut
short by the task cap says so, with the size it would have been.
"""
from __future__ import annotations

import collections
import hashlib
import json
from pathlib import Path

import pytest

from mofsbu import runner
from mofsbu.registry import BlobStore, Registry, find, incoming_routes
from mofsbu.runner import enumerate_plan, estimate, run
from mofsbu.spec import SPEC_VERSION, BuildSpec, MetalSpec, MoleculeSpec

SLICE = Path(__file__).resolve().parents[1] / "data" / "reference" / "spec_ni_thq_cl_slice.json"

#: A v7 spec with a co-ligand, both CNs and pathways on — the shape the default changes.
OLD_V7 = {
    "spec_version": 7,
    "molecules": [{"name": "catechol", "smiles": "Oc1ccccc1O", "multiplicity": 1,
                   "max_deprotonations": 0}],
    "metals": [{"symbol": "Ni", "oxidation_state": 2, "spin_class": "hs"}],
    "coordination": [4, 6], "ligands_per_metal": [1, 2], "binding": ["mono", "chelate"],
    "co_ligand": "O", "pathways": True, "seed": 7,
}

#: `(task count, sha256 of every task's kind/payload/priority/refs)` as planned by the
#: runner before `co_ligand_counts` existed (origin/main at caf861f).  Pinned rather than
#: recomputed: comparing the migrated spec against an explicit `fill` spec would only show
#: the migration agrees with itself.
BEFORE = {
    "old_v7": (64, "541cc284e1e04eeb77ff4ebd9f1d8ffbfabb3001ae86c72f0952add84ee11b6c"),
    "slice": (722, "4bf4dfe92f952ca3d148d479cc60b44a290254019162183c36e722ebcb8bdba9"),
}


def plan_digest(spec: BuildSpec) -> tuple[int, str]:
    tasks = enumerate_plan(spec).tasks
    body = [[t.kind, t.payload, t.priority, list(t.task_refs)] for t in tasks]
    return len(tasks), hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest()


def aqua_spec(**over) -> BuildSpec:
    """Ni(II) + one monodentate catechol at CN 4: three waters when full."""
    base = dict(
        molecules=(MoleculeSpec("catechol", "Oc1ccccc1O", 1, 0),),
        metals=(MetalSpec("Ni", 2, "hs"),), coordination="4", ligands_per_metal="1",
        binding=("mono",), co_ligand="O", geometries=("square_planar",), seed=7,
        pathways=True)
    base.update(over)
    return BuildSpec(**base)


def places(planned) -> list[dict]:
    return [t.payload for t in planned.tasks if t.kind == "place"]


def ladder(planned):
    """Place indices, grow edges child -> parents, and undirected neighbours."""
    parents: dict[int, set[int]] = collections.defaultdict(set)
    near: dict[int, set[int]] = collections.defaultdict(set)
    for t in planned.tasks:
        if t.kind == "grow":
            child, parent = t.payload["child_task"], t.payload["parent_task"]
            parents[child].add(parent)
            near[child].add(parent)
            near[parent].add(child)
    nodes = [i for i, t in enumerate(planned.tasks) if t.kind == "place"]
    return nodes, parents, near


# ── old specs replay ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("version", range(1, SPEC_VERSION))
def test_every_older_version_reads_as_fill(version):
    spec = BuildSpec.from_dict({**OLD_V7, "spec_version": version})
    assert spec.co_ligand_counts == "fill"


def test_an_old_spec_plans_exactly_the_tasks_it_planned_before():
    assert plan_digest(BuildSpec.from_dict(dict(OLD_V7))) == BEFORE["old_v7"]


def test_the_reference_slice_plans_exactly_the_tasks_it_planned_before():
    """The spec the MVP registry was built from, read from its file as saved (v7)."""
    assert plan_digest(BuildSpec.load(SLICE)) == BEFORE["slice"]


def test_a_new_spec_defaults_to_range_with_a_window_of_two_and_round_trips():
    spec = aqua_spec()
    assert (spec.co_ligand_counts, spec.co_ligand_window) == ("range", 2)
    again = BuildSpec.from_json(spec.to_json())
    assert again == spec and again.digest == spec.digest


@pytest.mark.parametrize("over", [{"co_ligand_counts": "some"}, {"co_ligand_window": -1},
                                  {"co_ligand_window": True}, {"co_ligand_window": None}])
def test_a_bad_value_is_refused_with_the_choices(over):
    with pytest.raises(ValueError, match="co_ligand"):
        aqua_spec(**over)


# ── the new default plans the lower states ───────────────────────────────────

def test_range_builds_the_counts_inside_the_window_in_the_same_polyhedron():
    """Full down to two fewer by default; the uncovered vertices are left empty, CN kept."""
    products = [p for p in places(enumerate_plan(aqua_spec(pathways=False)))
                if p["components"]]
    # Two kinds (either phenol O binds), each at 3, 2 and 1 waters.
    assert sorted(p["n_co"] for p in products) == [1, 1, 2, 2, 3, 3]
    assert all(p["cn"] == 4 and p["n_co"] + p["n_vacant"] == 3 for p in products)

    wide = [p for p in places(enumerate_plan(aqua_spec(pathways=False, co_ligand_window=3)))
            if p["components"]]
    assert sorted({p["n_co"] for p in wide}) == [0, 1, 2, 3]


def test_a_higher_count_comes_only_from_a_cn_the_spec_asked_for():
    """Above full needs more vertices: CN 6 in the list gives 3..5 waters; CN 4 alone never
    builds a fourth, and no CN outside the list is invented."""
    def counts(cns):
        out: dict[int, set[int]] = collections.defaultdict(set)
        for p in places(enumerate_plan(aqua_spec(pathways=False, coordination=cns,
                                                 geometries=None))):
            out[p["cn"]].add(p["n_co"])
        return dict(out)

    assert counts("4") == {4: {1, 2, 3}}
    assert counts("4,6") == {4: {1, 2, 3}, 6: {3, 4, 5}}


def test_range_does_not_need_allow_unsaturated():
    """The field is the request for the vacant states; `allow_unsaturated` governs only
    the case with no co-ligand to fill with, and that case is unchanged."""
    lean = places(enumerate_plan(aqua_spec(pathways=False, allow_unsaturated=False)))
    assert any(p["n_vacant"] for p in lean)
    none = aqua_spec(pathways=False, co_ligand=None)
    assert not places(enumerate_plan(none)), "no co-ligand, not allowed unsaturated: skipped"


def test_the_free_co_ligand_is_planned_as_a_reagent():
    planned = enumerate_plan(aqua_spec())
    assert [t.payload for t in planned.tasks if t.kind == "co_ligand"] == [{"smiles": "O"}]
    assert not [t for t in enumerate_plan(aqua_spec(co_ligand_counts="fill")).tasks
                if t.kind == "co_ligand"]


def test_any_named_co_ligand_is_varied_not_only_water():
    spec = aqua_spec(pathways=False, co_ligand="CC#N")      # acetonitrile as the solvent
    assert sorted({p["n_co"] for p in places(enumerate_plan(spec)) if p["components"]}) \
        == [1, 2, 3]


# ── and joins them into one ladder ───────────────────────────────────────────

def test_a_co_ligand_step_adds_one_onto_the_rung_below():
    planned = enumerate_plan(aqua_spec())
    co_ix = next(i for i, t in enumerate(planned.tasks) if t.kind == "co_ligand")
    steps = [t.payload for t in planned.tasks
             if t.kind == "grow" and "co_ligand" in t.payload["component"]]
    assert steps
    for step in steps:
        parent = planned.tasks[step["parent_task"]].payload
        child = planned.tasks[step["child_task"]].payload
        assert step["ligand_task"] == co_ix, "the free co-ligand is the step's reagent"
        assert parent["components"] == child["components"]
        assert parent["n_co"] == child["n_co"] - 1
        assert parent["n_vacant"] == child["n_vacant"] + 1
        assert (parent["cn"], parent["geometry"]) == (child["cn"], child["geometry"])


def test_no_rung_leaves_more_than_the_window_empty_and_no_bare_metal():
    for window in (1, 2, 3):
        rungs = places(enumerate_plan(aqua_spec(co_ligand_window=window)))
        assert max(p["n_vacant"] for p in rungs) == window
        assert not [p for p in rungs if not p["components"] and not p["n_co"]], window


def test_the_seeds_are_one_chain_rooted_at_the_lowest_state_in_the_window():
    """Ni(H2O)3 <- Ni(H2O)2 + H2O, and Ni(H2O)2 has nothing below it: at CN 4 with a
    window of 2 it is the lowest co-ligand state of the ligand-free centre."""
    planned = enumerate_plan(aqua_spec())
    nodes, parents, _ = ladder(planned)
    seeds = {planned.tasks[i].payload["n_co"]: i for i in nodes
             if not planned.tasks[i].payload["components"]}
    assert sorted(seeds) == [2, 3]
    assert parents[seeds[3]] == {seeds[2]}
    assert not parents[seeds[2]], "the root has no rung below it"


@pytest.mark.parametrize("over", [{}, {"coordination": "4,6", "geometries": None},
                                  {"ligands_per_metal": "1~2"}])
def test_every_metal_cn_geometry_is_one_connected_ladder(over):
    """Every seed and every product of one metal, CN and polyhedron is reachable from every
    other through steps — the property the MVP registry's two unrelated seeds lacked."""
    planned = enumerate_plan(aqua_spec(**over))
    nodes, _, near = ladder(planned)
    groups: dict[tuple, list[int]] = collections.defaultdict(list)
    for i in nodes:
        p = planned.tasks[i].payload
        groups[(p["metal"], p["cn"], p["geometry"])].append(i)
    for key, members in groups.items():
        seen, todo = set(), [members[0]]
        while todo:
            j = todo.pop()
            if j not in seen:
                seen.add(j)
                todo += list(near[j])
        assert seen >= set(members), key


def test_a_window_of_zero_plans_no_step_and_says_so():
    planned = enumerate_plan(aqua_spec(co_ligand_window=0))
    assert not [t for t in planned.tasks if t.kind == "grow"]
    assert any("window 0" in d["reason"] for d in planned.diagnostics)


def test_a_charged_co_ligand_is_not_stepped_and_says_why():
    """B20: `place` leaves a co-ligand's charge out, so a join-built step would fork the
    node. The planner says so instead of queueing steps that land somewhere else."""
    planned = enumerate_plan(aqua_spec(co_ligand="[Cl-]"))
    assert not [t for t in planned.tasks
                if t.kind == "grow" and "co_ligand" in t.payload["component"]]
    assert any("charged" in d["reason"] and "B20" in d["hint"]
               for d in planned.diagnostics), planned.diagnostics
    lower = {p["n_co"] for p in places(planned) if p["components"]}
    assert lower == {1, 2, 3}, "the lower counts are still built"


# ── the cap says it bit ──────────────────────────────────────────────────────

def test_a_ladder_cut_by_the_cap_says_so_and_how_big_it_was(monkeypatch):
    uncapped = len(enumerate_plan(aqua_spec()).tasks)
    monkeypatch.setattr(runner, "MAX_PATHWAY_TASKS", 8)
    planned = enumerate_plan(aqua_spec())
    cap = [d for d in planned.diagnostics if "cap" in d]
    assert len(cap) == 1 and cap[0]["cap"] == 8
    assert cap[0]["uncapped_tasks"] == uncapped > len(planned.tasks)
    assert "co_ligand_window" in cap[0]["hint"]
    guess = estimate(aqua_spec())
    assert guess["capped"] and guess["capped"][0]["uncapped_tasks"] == uncapped


def test_an_uncapped_estimate_says_nothing_was_cut():
    guess = estimate(aqua_spec())
    assert guess["ok"] and guess["capped"] == []
    assert (guess["co_ligand_counts"], guess["co_ligand_window"]) == ("range", 2)
    assert guess["by_kind"]["co_ligand"] == 1


# ── performed: the steps land on the placed nodes ────────────────────────────

@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as r:
        r.migrate()
        yield r


def test_the_co_ligand_ladder_is_built_and_connected(reg):
    summary = run(reg, aqua_spec())
    assert summary["status"] == "done"
    assert not summary["counts"].get("failed") and not summary["counts"].get("rejected"), (
        summary["counts"])

    rows = {r["id"]: dict(r) for r in find(reg, limit=200)}
    by_formula = {r["formula"]: r["id"] for r in rows.values()}
    assert "Ni" not in by_formula, "no bare-metal row"
    complexes = {i for i, r in rows.items() if r["n_metals"]}
    # Ni(H2O)2, Ni(H2O)3, Ni(cat)(H2O)1..3 — the two catechol kinds are one identity (D2).
    assert len(complexes) == 5, sorted(rows[i]["display_label"] for i in complexes)

    # The water step cites the rung below and the free water.
    steps = [r for r in incoming_routes(reg, by_formula["NiH6O3"]) if r["reagent_ids"]]
    assert [sorted(s["reagent_ids"]) for s in steps] == [
        sorted([by_formula["NiH4O2"], by_formula["H2O"]])]

    # One ladder: every complex reaches every other through reagent-bearing routes.
    near: dict[int, set[int]] = collections.defaultdict(set)
    for sid in complexes:
        for route in incoming_routes(reg, sid):
            for other in route["reagent_ids"]:
                if other in complexes:
                    near[sid].add(other)
                    near[other].add(sid)
    seen, todo = set(), [by_formula["NiH4O2"]]
    while todo:
        j = todo.pop()
        if j not in seen:
            seen.add(j)
            todo += list(near[j])
    assert seen == complexes
