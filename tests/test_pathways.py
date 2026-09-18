"""The ladder a ligand-count sweep builds, and the steps between its rungs.

A run that asks for one to three copies of a ligand builds three complexes, and until now
nothing in the registry said the second is the first plus a ligand. `pathways` plans the
intermediate below each product and performs the addition with `assembly.join`, so the
relationship is an edge that was carried out rather than an inference from two formulas.

The claims under test are the ones that make that edge worth having: the step lands on the
identity the `place` path built (D2 — one node, two routes), it carries its reagents and
its choice vector, and where the relationship is NOT an addition it is refused with the
reason rather than faked.

`estimate` is here too, because it is the same enumeration seen from the other side: the
page's number has to be the number submitting produces.
"""
from __future__ import annotations

import pytest

from mofsbu.registry import BlobStore, Registry, find, incoming_routes
from mofsbu.registry.jobs import task_rows
from mofsbu.runner import enumerate_plan, estimate, plan, run
from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec, PocketPredicate
from mofsbu.versions import ALGO_VERSIONS


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as r:
        r.migrate()
        yield r


def ladder_spec(**over) -> BuildSpec:
    """Zn + catechol, one to two copies — the smallest sweep with a rung to climb."""
    base = dict(
        molecules=(MoleculeSpec("catechol", "Oc1ccccc1O", 1),),
        metals=(MetalSpec("Zn", 2, "ls"),), coordination="4",
        ligands_per_metal="1~2", binding=("chelate",), co_ligand=None,
        pocket=PocketPredicate(ring_size=5, n_anionic=2), seed=7, pathways=True,
        geometries=("tetrahedral",))
    base.update(over)
    return BuildSpec(**base)


def labelled(reg, fragment: str) -> dict:
    """The one structure whose label carries this fragment count."""
    rows = [r for r in find(reg, limit=200) if fragment in r["display_label"]]
    assert len(rows) == 1, [r["display_label"] for r in rows]
    return rows[0]


# ── the ladder is planned ────────────────────────────────────────────────────

def test_a_sweep_plans_a_step_between_its_rungs():
    """M(L) and M(L)2 are both asked for; the step from one to the other is the addition."""
    steps = [t for t in enumerate_plan(ladder_spec()).tasks if t.kind == "grow"]
    assert steps, "a two-rung sweep has a step in it"
    assert all(t.priority < 0 for t in steps), "rungs are built before steps are taken"
    assert all(t.payload["component"]["count"] == 1 for t in steps), "one copy per step"


def test_without_pathways_nothing_changes():
    """The flag is off by default and a spec that does not ask for the ladder plans
    exactly the run it planned before the field existed."""
    plain = enumerate_plan(ladder_spec(pathways=False))
    assert plain.by_kind().get("grow", 0) == 0
    lean = [t.payload for t in plain.tasks if t.kind == "place"]
    rich = [t.payload for t in enumerate_plan(ladder_spec()).tasks if t.kind == "place"]
    # The ladder ADDS intermediates; it never changes or drops a product that was asked for.
    for payload in lean:
        assert payload in rich


def test_the_intermediate_below_each_product_is_the_unsaturated_one():
    """Not the co-ligand-filled one: filling those vertices would make the step a
    substitution, which a join cannot express."""
    planned = enumerate_plan(ladder_spec())
    steps = [t for t in planned.tasks if t.kind == "grow"]
    assert steps
    for task in steps:
        parent = planned.tasks[task.payload["parent_task"]].payload
        child = planned.tasks[task.payload["child_task"]].payload
        denticity = task.payload["component"]["denticity"]
        # The vertices the added copy will occupy are free on the parent and taken on the
        # child, and the co-ligands are untouched — which is what "addition" means here.
        assert parent["n_vacant"] - child["n_vacant"] == denticity
        assert parent["n_co"] == child["n_co"]
        assert parent["cn"] == child["cn"] and parent["geometry"] == child["geometry"]


# ── the step is performed, and it converges ──────────────────────────────────

def test_the_step_reaches_the_structure_the_place_task_built(reg):
    """D2 after a round trip: the join-built product and the placed one are ONE row,
    and the step is a second incoming route to it."""
    summary = run(reg, ladder_spec())
    assert summary["status"] == "done"
    assert not summary["counts"].get("failed"), summary["counts"]

    product = labelled(reg, "]2")                     # Zn[catechol…]2
    routes = incoming_routes(reg, product["id"])
    placed = [r for r in routes if not r["reagent_ids"]]
    steps = [r for r in routes if r["reagent_ids"]]
    assert placed and steps, "one node, two kinds of route"
    assert product["n_incoming_routes"] == len(routes)


def test_the_step_carries_its_reagents_and_its_choice_vector(reg):
    """What makes the edge a pathway rather than a note: the rung it came from, the
    ligand that was added, and the vector that regenerates the move."""
    run(reg, ladder_spec())
    product, rung = labelled(reg, "]2"), labelled(reg, "Zn[catechol-2Hcfg(0,1)sep(3)] ")
    ligand = labelled(reg, "[catechol-2Hcfg(0,1)sep(3)] q-2")
    step = next(r for r in incoming_routes(reg, product["id"]) if r["reagent_ids"])
    assert sorted(step["reagent_ids"]) == sorted([rung["id"], ligand["id"]])
    # The recipe version, read rather than spelled: what this asserts is that the digest
    # is a VERSIONED key — a bare hash could not say it was stale — and pinning the
    # literal would turn every honest bump into a test edit.
    assert step["choice_vector_digest"].startswith(f"{ALGO_VERSIONS['choice_vector']}:")
    assert step["depth"] == 2                          # two ligand copies, two additions


def test_the_ladder_reaches_down_to_the_bare_centre(reg):
    """Every rung is reached from the one below it, so the chain is walkable end to end."""
    run(reg, ladder_spec())
    complexes = {r["id"]: r["display_label"] for r in find(reg, metal="Zn", limit=200)}
    chain, row_id = [], labelled(reg, "]2")["id"]
    while True:
        chain.append(complexes[row_id])
        step = next((r for r in incoming_routes(reg, row_id) if r["reagent_ids"]), None)
        if step is None:
            break
        below = [i for i in step["reagent_ids"] if i in complexes and i != row_id]
        assert len(below) == 1, below
        row_id = below[0]
    assert len(chain) == 3, chain                      # Zn(L)2 <- Zn(L) <- Zn


def test_the_intermediates_keep_their_open_vertices(reg):
    """A rung is a coordinatively unsaturated species and is stored as one — the vacancies
    are what the next step joins onto, so losing them would break the chain."""
    run(reg, ladder_spec())
    rung = labelled(reg, "Zn[catechol-2Hcfg(0,1)sep(3)] ")
    assert rung["n_open_sites"] == 2


# ── where it does not make sense, it says so ─────────────────────────────────

def test_a_saturated_parent_is_refused_as_a_substitution(reg):
    """The honest half of "where it makes sense". With the vertices filled by a co-ligand,
    getting from one rung to the next EXCHANGES ligands, and a join cannot express that."""
    from mofsbu.registry.jobs import add_task, claim_task, complete_task, create_run
    from mofsbu.runner import _Rejected, execute

    spec = ladder_spec(co_ligand="O")
    run_id = create_run(reg, spec)
    # Built by hand, because the planner never plans this step: it is the case the
    # planner avoids by choosing the unsaturated rung as the parent.
    saturated = next(t.payload for t in enumerate_plan(spec).tasks
                     if t.kind == "place" and t.payload["n_co"] and not t.payload["n_vacant"])
    parent_id = add_task(reg, run_id, "place", saturated)
    parent = claim_task(reg, run_id)
    out = execute(reg, parent, spec)
    complete_task(reg, parent.id, structure_id=out.structure_id,
                  geometry_id=out.geometry_id)

    add_task(reg, run_id, "grow", {
        "parent_task": parent_id, "child_task": None, "ligand_task": None,
        "component": dict(saturated["components"][0], count=1),
        "metal": 0, "cn": saturated["cn"], "geometry": saturated["geometry"]})
    with pytest.raises(_Rejected) as exc:
        execute(reg, claim_task(reg, run_id), spec)
    assert exc.value.code == "pathway_no_open_vertex"
    assert "SUBSTITUTION" in str(exc.value)


def test_a_step_whose_rung_was_rejected_says_which_one(reg):
    """The step is not wrong — the rung below it is missing, and a reader needs to be
    sent to the task that explains why rather than to this one."""
    from mofsbu.registry.jobs import add_task, claim_task, create_run, fail_task
    from mofsbu.runner import _Rejected, execute

    spec = ladder_spec()
    run_id = create_run(reg, spec)
    rung = add_task(reg, run_id, "place", {"components": []})
    fail_task(reg, rung, "QC said no", rejected=True, code="qc_clash")
    add_task(reg, run_id, "grow", {"parent_task": rung, "component": {}})
    with pytest.raises(_Rejected) as exc:
        execute(reg, claim_task(reg, run_id), spec)
    assert exc.value.code == "pathway_parent_missing"
    assert "rejected" in str(exc.value)


# ── the estimate is the same enumeration ─────────────────────────────────────

def test_the_estimate_is_what_submitting_queues(reg):
    """If these two ever disagree the number on the page is decoration."""
    spec = ladder_spec()
    guess = estimate(spec)
    _, n_tasks = plan(reg, spec)
    assert guess["ok"] and guess["tasks"] == n_tasks
    assert guess["by_kind"]["grow"] >= 1


def test_the_estimate_grows_with_the_range():
    """The point of showing it: a wider sweep is visibly a bigger run, before submitting."""
    one = estimate(ladder_spec(ligands_per_metal="1"))["tasks"]
    two = estimate(ladder_spec(ligands_per_metal="1~2"))["tasks"]
    assert two > one


def test_the_estimate_counts_the_relaxations_as_a_ceiling():
    """`queue_relax` skips what a registry has already relaxed, and this has no registry
    to ask, so the number is stated as an upper bound rather than counted as certain."""
    plain = estimate(ladder_spec())
    assert plain["relaxations"] == 0                     # construct relaxes nothing


def test_an_unrunnable_spec_is_reported_rather_than_raised():
    """The page asks while the spec is being edited; half-finished is the normal input."""
    refused = estimate(ladder_spec(degree=2))
    assert refused["ok"] is False and "S4.1" in refused["reason"]


def test_the_estimate_writes_nothing(reg, tmp_path):
    before = reg.count("structures"), reg.count("runs")
    estimate(ladder_spec())
    assert (reg.count("structures"), reg.count("runs")) == before
