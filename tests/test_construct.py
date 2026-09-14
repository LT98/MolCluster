"""M5/S4: construction as a decision tree, and the replay that makes a leaf regenerable.

Three groups, matching §6.7's three kinds of choice: the Kind-C inferences that must
branch or refuse, the Kind-B tree itself, and replay — the property that turns a stored
conformer from frozen coordinates into an object that can be rebuilt at a higher fidelity.
"""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu._types import AmbiguousSpecError, NotBuiltYet
from mofsbu.assembly.construct import (
    GEOMETRIES_FOR_CN, construct, enumerate_constructions, geometry_branches,
    ligand_block, metal_block, metal_block_branches, replayable, resolve_geometry)
from mofsbu.assembly.join import grow
from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.identity.keys import l1_graph_hash


@pytest.fixture(scope="module")
def aqua():
    return ligand_block(embed_molecule(mol_from_smiles("O")), charge=0, name="aqua")


@pytest.fixture(scope="module")
def acetate():
    return ligand_block(embed_molecule(mol_from_smiles("CC(=O)[O-]")), charge=-1,
                        name="acetate")


@pytest.fixture()
def zn():
    return metal_block("Zn", 2, cn=4, geometry="tetrahedral")


# ── Kind C: inference branches or refuses, never defaults ────────────────────

def test_cn_four_does_not_determine_a_geometry():
    """The headline Kind-C case: tetrahedral and square planar are different structures."""
    with pytest.raises(AmbiguousSpecError) as exc:
        resolve_geometry(4)
    assert "square_planar" in str(exc.value) and "tetrahedral" in str(exc.value)


def test_cn_six_is_determined_and_says_so():
    assert resolve_geometry(6) == "octahedral"
    assert geometry_branches(6) == ("octahedral",)


def test_the_ambiguity_is_enumerable_rather_than_merely_refused(zn):
    """Ground rule 5 is "branch explicitly", not "give up"."""
    branches = metal_block_branches("Zn", 2, cn=4)
    assert len(branches) == 2
    names = {b.graph.name for b in branches}
    assert names == {"Zn(tetrahedral)", "Zn(square_planar)"}
    # ...and they really are different centres, not two labels on one.
    axes = [tuple(sorted(tuple(np.round(s.frame["axis"], 6)) for s in b.sites))
            for b in branches]
    assert axes[0] != axes[1]


def test_every_cn_in_the_placers_table_is_reachable():
    """Derived from the placer's own table, so a geometry added there cannot go missing."""
    assert set(GEOMETRIES_FOR_CN) == {2, 3, 4, 5, 6}
    for cn, options in GEOMETRIES_FOR_CN.items():
        for name in options:
            assert resolve_geometry(cn, name) == name


def test_a_geometry_that_is_not_that_cn_is_refused():
    with pytest.raises(ValueError, match="not a CN-6 geometry"):
        resolve_geometry(6, "tetrahedral")


def test_multiplicity_is_derived_from_the_spin_class_not_defaulted():
    """A high-spin Fe(III) is a sextet. A module that preaches "do not guess" may not."""
    assert metal_block("Fe", 3, cn=6, spin_class="hs").graph.multiplicity == 6
    assert metal_block("Fe", 3, cn=6, spin_class="ls").graph.multiplicity == 2


def test_a_ligand_with_no_conformer_is_refused():
    """No conformer means no frames, and a site without a frame cannot be joined."""
    with pytest.raises(AmbiguousSpecError, match="no conformer"):
        ligand_block(mol_from_smiles("O"), charge=0, name="aqua")


# ── Kind B: the tree ─────────────────────────────────────────────────────────

def test_every_open_vertex_is_its_own_branch(zn, aqua):
    """Four leaves for four vertices, and they share one L1 — which is left standing.

    Collapsing them here by hashing the product would also merge near-degenerate cis and
    trans, which D10 exists to keep apart. Telling a symmetry duplicate from a real branch
    needs L2 and θ_geom clustering (S5/S6).
    """
    tree = enumerate_constructions(zn, [aqua], degree=1)
    assert len(tree) == 4
    assert len({l1_graph_hash(c.block.graph) for c in tree}) == 1


def test_a_live_torsion_branches_and_a_free_one_does_not(zn, aqua, acetate):
    """`TORSION_FREE` is the declared guard against conformer explosion, measured here.

    Asserted on the recorded wells rather than on a leaf count, because the count also
    carries the donor factor — acetate perceives as TWO equivalent oxygens, so its tree is
    four times aqua's (2 donors x 2 wells) and a bare `2 *` would be testing the wrong
    multiplier while looking like it tested this one.
    """
    free = enumerate_constructions(zn, [aqua], degree=1)
    live = enumerate_constructions(zn, [acetate], degree=1)
    assert {s["torsion_well"] for c in free for s in c.steps} == {0}
    assert {s["torsion_well"] for c in live for s in c.steps} == {0, 1}

    n_free_donors, n_live_donors = len(aqua.open_donors()), len(acetate.open_donors())
    assert len(free) == len(zn.open_vacancies()) * n_free_donors * 1
    assert len(live) == len(zn.open_vacancies()) * n_live_donors * 2


def test_the_tree_is_deterministic(zn, aqua, acetate):
    """A capped tree that depended on iteration order would truncate differently each run."""
    first = enumerate_constructions(zn, [aqua, acetate], degree=1)
    second = enumerate_constructions(zn, [aqua, acetate], degree=1)
    assert [c.digest for c in first] == [c.digest for c in second]


def test_depth_two_contains_its_depth_one_parents(zn, aqua):
    tree = enumerate_constructions(zn, [aqua], degree=2, max_products=100)
    assert {len(c.steps) for c in tree} == {1, 2}


def test_the_cap_is_reported_and_not_silent(zn, aqua, acetate):
    """A silently truncated tree is worse than a big one."""
    tree = enumerate_constructions(zn, [aqua, acetate], degree=2, max_products=12)
    assert tree.capped and len(tree) <= 12
    assert "capped" in tree.explain()


def test_a_tree_with_nothing_to_join_explains_itself(aqua, acetate):
    """Silence was a real bug once; an enumerator is exactly where it would come back."""
    tree = enumerate_constructions(aqua, [acetate], degree=1)     # a ligand has no vertices
    assert len(tree) == 0
    assert "nothing was offered" in tree.explain()


def test_a_tree_whose_branches_were_all_refused_says_why(zn, aqua):
    tree = enumerate_constructions(zn, [aqua], degree=1, modes=("chelate",))
    assert len(tree) == 0
    assert tree.refusals
    assert "every branch was refused" in tree.explain()


def test_degree_beyond_one_needs_geometry(zn, aqua):
    """Openness is a per-geometry fact, so a graph-only product cannot be grown onto."""
    with pytest.raises(ValueError, match="degree > 1 needs geometry"):
        enumerate_constructions(zn, [aqua], degree=2, with_geometry=False)


def test_growth_onto_a_second_centre_still_waits_for_the_placer(zn, acetate):
    """The S3 seam is not softened by being reached through the enumerator."""
    second = metal_block("Cu", 2, cn=4, geometry="square_planar")
    bridged = enumerate_constructions(zn, [acetate], degree=1).leaves[0].block
    with pytest.raises(NotBuiltYet, match="place_multicentre"):
        enumerate_constructions(bridged, [second], degree=1)


# ── replay ───────────────────────────────────────────────────────────────────

def test_replaying_a_leaf_reproduces_it_exactly(zn, aqua, acetate):
    """Exit gate 3: a stored conformer is a regenerable object, not frozen coordinates."""
    tree = enumerate_constructions(zn, [aqua, acetate], degree=2, max_products=40)
    leaf = next(c for c in tree if len(c.steps) == 2)
    again = construct(zn, [aqua, acetate], leaf)
    assert again.digest == leaf.digest
    assert l1_graph_hash(again.block.graph) == l1_graph_hash(leaf.block.graph)
    assert np.array_equal(again.block.geometry, leaf.block.geometry)


def test_replay_works_from_the_stored_choice_vector(zn, aqua):
    """The path survives the round trip through the form a `geometries` row holds."""
    leaf = enumerate_constructions(zn, [aqua], degree=1).leaves[0]
    stored = replayable(leaf)
    rebuilt = construct(zn, [aqua], stored.data)
    assert np.array_equal(rebuilt.block.geometry, leaf.block.geometry)


def test_the_same_steps_in_either_order_reach_one_node(zn, aqua, acetate):
    """"One node, two routes" where the routes are generated, not merely recognised."""
    tree = enumerate_constructions(zn, [aqua, acetate], degree=2, max_products=60)
    leaf = next(c for c in tree if len({s["partner"] for s in c.steps}) == 2)
    forward = construct(zn, [aqua, acetate], leaf.steps)
    backward = construct(zn, [aqua, acetate], tuple(reversed(leaf.steps)))
    assert forward.key == backward.key
    assert l1_graph_hash(forward.block.graph) == l1_graph_hash(backward.block.graph)


def test_a_step_that_does_not_say_which_partner_refuses(zn, aqua, acetate):
    leaf = enumerate_constructions(zn, [aqua], degree=1).leaves[0]
    orphaned = [{k: v for k, v in leaf.steps[0].items() if k != "partner"}]
    with pytest.raises(AmbiguousSpecError, match="does not say which partner"):
        construct(zn, [aqua, acetate], orphaned)


def test_a_step_with_no_donor_pair_refuses(zn, aqua):
    with pytest.raises(AmbiguousSpecError, match="nothing to replay"):
        construct(zn, [aqua], [{"partner": 0, "mode": "mono"}])


def test_replay_against_different_partners_refuses_rather_than_improvising(zn, aqua, acetate):
    """The seed and the partners are part of the path; a path alone is not a structure."""
    leaf = enumerate_constructions(zn, [acetate], degree=1).leaves[0]
    with pytest.raises(AmbiguousSpecError, match="recorded against"):
        construct(zn, [aqua], leaf.steps)          # aqua has no atom at acetate's donor


def test_a_choice_vector_with_no_steps_is_not_a_path(zn, aqua):
    with pytest.raises(AmbiguousSpecError, match="must carry its `steps`"):
        construct(zn, [aqua], {"op": "construct", "n_steps": 2})


# ── grow, on top of it ───────────────────────────────────────────────────────

def test_grow_carries_the_path_not_the_last_join(zn, aqua):
    results = grow(zn, (aqua,), degree=2, max_products=50)
    assert results
    deepest = max(results, key=lambda r: r.choice_vector["n_steps"])
    assert deepest.choice_vector["op"] == "construct"
    assert deepest.choice_vector["n_steps"] == 2
    assert sorted(deepest.atom_map) == zn.graph.nodes()      # the seed's map, identity


def test_grow_strain_is_the_worst_step_on_the_path(zn, aqua):
    """A construction is as viable as its weakest point — the reading M8 gives barriers."""
    for result in grow(zn, (aqua,), degree=2, max_products=50):
        assert result.strain == pytest.approx(0.0)           # every mono join is strainless
