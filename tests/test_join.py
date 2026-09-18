"""M5/S3: the recursive assembly operation itself.

Two blocks in, one block out, plus the atom maps and the choice vector that explain it.
The tests are grouped by the four claims the milestone rests on — the bond is real, the
bookkeeping survives the join, the product is reproducible, and the M5/M6 seam holds.
"""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu._types import AmbiguousSpecError, NotBuiltYet
from mofsbu.assembly.choice import ChoiceVector
from mofsbu.assembly.join import (
    BuildingBlock, IncompatibleJoin, join, join_chelate)
from mofsbu.geometry.embed import coordinates, embed_molecule
from mofsbu.geometry.placer import site_vectors
from mofsbu.graph._types import EdgeType, TypedGraph
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles
from mofsbu.identity.keys import l1_graph_hash
from mofsbu.sites.model import perceive, vacancy_sites
from mofsbu.sites.state import SiteStatus, refresh_state


def _block(graph: TypedGraph, sites, coords) -> BuildingBlock:
    states = refresh_state(sites, coords, graph=graph,
                           symbols=[graph.label(i).element for i in graph.nodes()])
    return BuildingBlock(graph=graph, sites=tuple(sites), geometry=coords,
                         state={(s.atom_idx, s.slot): s for s in states})


def metal_block(symbol: str = "Zn", charge: int = 2, cn: int = 6,
                geometry: str = "octahedral") -> BuildingBlock:
    """One metal atom at the origin with every coordination vertex empty."""
    g = TypedGraph(charge=charge, multiplicity=1, name=symbol)
    g.add_atom(symbol, formal_charge=charge, oxidation_state=charge, spin_class="ls")
    directions = list(site_vectors(geometry, cn, 1.0))
    return _block(g, vacancy_sites(0, [0.0, 0.0, 0.0], directions), np.zeros((1, 3)))


def ligand_block(smiles: str, charge: int, name: str) -> BuildingBlock:
    mol = embed_molecule(mol_from_smiles(smiles))
    g = from_rdkit(mol, charge=charge, multiplicity=1, name=name)
    return _block(g, perceive(mol), coordinates(mol))


@pytest.fixture(scope="module")
def aqua():
    return ligand_block("O", 0, "aqua")


@pytest.fixture(scope="module")
def acetate():
    return ligand_block("CC(=O)[O-]", -1, "acetate")


def first_join(ligand: BuildingBlock, metal: BuildingBlock, **kw):
    return join(ligand, metal, ligand.open_donors()[0], metal.open_vacancies()[0], **kw)


# ── the bond is real ─────────────────────────────────────────────────────────

def test_the_join_makes_a_dative_bond_at_the_pair_distance(aqua):
    result = first_join(aqua, metal_block())
    product = result.block
    donor = result.atom_map[aqua.open_donors()[0].atom_idx]
    metal_idx = product.graph.metals()[0]
    assert product.graph.edge_type(donor, metal_idx) is EdgeType.DATIVE
    measured = float(np.linalg.norm(product.geometry[donor] - product.geometry[metal_idx]))
    assert measured == pytest.approx(result.choice_vector["d_ml"], abs=1e-6)


def test_the_metal_stays_put_and_the_ligand_moves(aqua):
    """The block with the vacancy is the anchor, so a growing cluster keeps one frame."""
    metal = metal_block()
    result = first_join(aqua, metal)
    metal_idx = result.block.graph.metals()[0]
    assert np.allclose(result.block.geometry[metal_idx], metal.geometry[0])
    ligand_rows = result.block.geometry[:len(aqua.graph)]
    assert not np.allclose(ligand_rows, aqua.geometry)


def test_charge_adds_and_multiplicity_is_derived(acetate):
    product = first_join(acetate, metal_block()).block.graph
    assert product.charge == acetate.graph.net_charge() + 2
    assert product.multiplicity == 1                 # two singlets, no unpaired electrons


def test_a_parent_with_no_multiplicity_refuses_to_guess(aqua):
    metal = metal_block()
    metal.graph.multiplicity = None
    with pytest.raises(AmbiguousSpecError, match="no multiplicity"):
        first_join(aqua, metal)


# ── the bookkeeping survives ─────────────────────────────────────────────────

def test_one_vertex_is_consumed_and_the_metals_others_survive(aqua):
    """A metal's vacancies share an atom, so consuming by atom index would take all six."""
    metal = metal_block(cn=6)
    result = first_join(aqua, metal)
    vacancies = [s for s in result.block.sites if s.is_vacancy]
    assert len(vacancies) == 5
    assert {s.slot for s in vacancies} == {1, 2, 3, 4, 5}      # slot 0 was filled


def test_the_bound_donor_is_kept_and_reads_as_occupied(aqua):
    """D5: `site_catalog` says a site EXISTS, `site_state` says what it is doing.

    A bound carboxylate oxygen is still a donor site. Dropping it would lose the answer to
    "what bound here, and how easily did it activate" — which is what a reaction edge needs
    to explain itself.
    """
    result = first_join(aqua, metal_block())
    donor = result.atom_map[aqua.open_donors()[0].atom_idx]
    state = result.block.state[(donor, 0)]
    assert state.status is SiteStatus.OCCUPIED
    assert state.ease is not None                 # still scored, because still a question


def test_the_open_count_falls_by_two(aqua):
    metal = metal_block()
    before = len(aqua.open_sites()) + len(metal.open_sites())
    assert len(first_join(aqua, metal).block.open_sites()) == before - 2


def test_frames_travel_with_the_block_that_moved(acetate):
    """A rigid move preserves everything a frame encodes, so the frames move with it.

    If they did not, the product's surviving donor would claim to sit where its parent
    sat before the join, and the next join would align against a stale position.
    """
    result = first_join(acetate, metal_block())
    for site in result.block.sites:
        if site.is_vacancy or site.frame is None:
            continue
        assert np.allclose(site.frame["origin"], result.block.geometry[site.atom_idx],
                           atol=1e-9)
        assert float(np.linalg.norm(site.frame["axis"])) == pytest.approx(1.0)


def test_both_atom_maps_are_injective_and_total(acetate):
    result = first_join(acetate, metal_block())
    for parent, mapping in ((acetate, result.atom_map), (metal_block(), result.partner_atom_map)):
        assert sorted(mapping) == parent.graph.nodes()
        assert len(set(mapping.values())) == len(mapping)
    assert not (set(result.atom_map.values()) & set(result.partner_atom_map.values()))


# ── the product is reproducible ──────────────────────────────────────────────

def test_the_same_product_from_two_orders_is_one_node(aqua):
    """The D2 claim at graph level: identity is on the node, sequence is on the edges."""
    metal = metal_block()
    ligand_first = first_join(aqua, metal).block.graph
    metal_first = join(metal, aqua, metal.open_vacancies()[0],
                       aqua.open_donors()[0]).block.graph
    assert [metal_first.label(i).element for i in metal_first.nodes()] != \
           [ligand_first.label(i).element for i in ligand_first.nodes()]
    assert l1_graph_hash(ligand_first) == l1_graph_hash(metal_first)


def test_replaying_the_same_join_reproduces_the_coordinates(acetate):
    metal = metal_block()
    first = first_join(acetate, metal).block.geometry
    second = first_join(acetate, metal).block.geometry
    assert np.array_equal(first, second)


def test_a_live_torsion_branches_and_a_free_one_does_not(acetate, aqua):
    """The well index is the conformer coordinate (D13) — and TORSION_FREE never branches."""
    metal = metal_block()
    syn = first_join(acetate, metal, torsion_well=0).block.geometry
    anti = first_join(acetate, metal, torsion_well=1).block.geometry
    assert not np.allclose(syn, anti)

    free_a = first_join(aqua, metal, torsion_well=0).block.geometry
    free_b = first_join(aqua, metal, torsion_well=1).block.geometry
    assert np.array_equal(free_a, free_b)


def test_the_choice_vector_keys_the_join(acetate):
    metal = metal_block()
    syn = ChoiceVector(first_join(acetate, metal, torsion_well=0).choice_vector)
    anti = ChoiceVector(first_join(acetate, metal, torsion_well=1).choice_vector)
    again = ChoiceVector(first_join(acetate, metal, torsion_well=0).choice_vector)
    assert syn.digest == again.digest and syn.digest != anti.digest
    assert ChoiceVector.from_json(syn.to_json()).digest == syn.digest


def test_an_out_of_range_well_wraps_rather_than_failing(acetate):
    """Wells are a cyclic set, so an index past the end is the same branch, not an error."""
    metal = metal_block()
    assert np.array_equal(first_join(acetate, metal, torsion_well=2).block.geometry,
                          first_join(acetate, metal, torsion_well=0).block.geometry)


# ── the seam, and the refusals ───────────────────────────────────────────────

def test_a_second_metal_completes_the_graph_AND_the_geometry(acetate):
    """D20: the bridged dimer EMERGES from two one-contact joins — nothing solved for it.

    This used to be the M5/M6 seam and it raised here, on the premise that a two-centre
    product needs a constrained placer. It does not. One donor onto one vertex is one
    contact, and one contact is satisfied by a rigid move of the donor's block whatever
    either block already contains, so the pose is determined and the M···M that falls out
    is a measurement rather than an input.

    `with_geometry=False` still works and still yields a graph-only product — that is the
    caller who wants an identity without paying for coordinates, not a fallback.
    """
    bridging = first_join(acetate, metal_block(cn=4, geometry="tetrahedral")).block
    free_donor = bridging.open_donors()[0]
    second = metal_block("Cu", 2, 4, "square_planar")

    result = join(bridging, second, free_donor, second.open_vacancies()[0])
    assert len(result.block.graph.metals()) == 2
    assert result.block.geometry is not None
    assert l1_graph_hash(result.block.graph)                 # an identity is derivable
    assert result.choice_vector["op"] == "join"

    graph_only = join(bridging, second, free_donor, second.open_vacancies()[0],
                      with_geometry=False)
    assert graph_only.block.geometry is None
    assert (l1_graph_hash(graph_only.block.graph)
            == l1_graph_hash(result.block.graph))


def test_a_graph_only_product_refuses_to_call_its_sites_open(acetate):
    """No state means `open` is unknown, and unknown is not the same as everything."""
    bridging = first_join(acetate, metal_block(cn=4, geometry="tetrahedral")).block
    second = metal_block("Cu", 2, 4, "square_planar")
    product = join(bridging, second, bridging.open_donors()[0],
                   second.open_vacancies()[0], with_geometry=False).block
    with pytest.raises(NotBuiltYet, match="no site state"):
        product.open_sites()


def test_two_donors_cannot_be_joined_and_the_verdict_says_why(aqua, acetate):
    with pytest.raises(IncompatibleJoin) as exc:
        join(aqua, acetate, aqua.open_donors()[0], acetate.open_donors()[0])
    assert "two donors" in str(exc.value)
    assert exc.value.verdict.feasible is False


def test_a_site_passed_with_the_wrong_block_is_refused(aqua):
    """Swapped arguments build an atom map onto the wrong parent, which nothing detects."""
    metal = metal_block()
    with pytest.raises(ValueError, match="not on the block it was passed with"):
        join(aqua, metal, metal.open_vacancies()[0], aqua.open_donors()[0])


def test_a_geometry_in_the_wrong_shape_is_refused(aqua):
    metal = metal_block()
    broken = BuildingBlock(graph=aqua.graph, sites=aqua.sites,
                           geometry=np.zeros((2, 3)), state=aqua.state)
    with pytest.raises(ValueError, match="graph has 3 atoms"):
        first_join(broken, metal)


# ── two points at once: the chelate join ─────────────────────────────────────
# `chelate_compatible` has judged pairs since S2 and `join` placed one donor, so every
# enumeration was monodentate and a bidentate ligand could be scored but never built.
# These cover the second contact: what it fixes, what it consumes, and what it refuses.

@pytest.fixture(scope="module")
def catecholate():
    """A chelator — both phenolate oxygens, converging on one centre."""
    from mofsbu.sites.perception import deprotonate
    from mofsbu.sites.protomers import labile_sites

    mol = mol_from_smiles("Oc1ccccc1O")
    base, _ = deprotonate(mol, labile_sites(mol))
    embedded = embed_molecule(base)
    graph = from_rdkit(embedded, charge=-2, multiplicity=1, name="catecholate")
    return _block(graph, perceive(embedded), coordinates(embedded))


def pocket(block: BuildingBlock):
    return [s for s in block.open_donors() if s.donor_type == "phenolate_O"][:2]


def vertex_pair(metal: BuildingBlock, *, trans: bool):
    """Two vertices that are opposite each other, or two that are not."""
    vacancies = metal.open_vacancies()
    first = vacancies[0]
    axis = np.asarray(first.frame["axis"], dtype=float)
    for other in vacancies[1:]:
        cosine = float(np.dot(axis, np.asarray(other.frame["axis"], dtype=float)))
        if (cosine < -0.95) is trans:
            return [first, other]
    raise AssertionError("no such pair on this polyhedron")


def test_a_chelate_puts_both_donors_at_their_own_pair_distance(catecholate):
    """The bite mismatch goes into the ANGLE, not into the bond lengths."""
    from mofsbu.geometry.distances import metal_donor_distance

    metal = metal_block(cn=6)
    result = join_chelate(metal, catecholate, vertex_pair(metal, trans=False),
                          pocket(catecholate))
    graph, coords = result.block.graph, np.asarray(result.block.geometry)
    centre = graph.metals()[0]
    bonded = [j for i, j, et in graph.edges() if i == centre and et is EdgeType.DATIVE]
    bonded += [i for i, j, et in graph.edges() if j == centre and et is EdgeType.DATIVE]
    assert len(bonded) == 2
    target = metal_donor_distance("Zn", "O").value
    for donor in bonded:
        assert float(np.linalg.norm(coords[donor] - coords[centre])) == pytest.approx(
            target, abs=1e-6)


def test_both_vertices_are_consumed_and_the_donors_read_as_occupied(catecholate):
    metal = metal_block(cn=6)
    product = join_chelate(metal, catecholate, vertex_pair(metal, trans=False),
                           pocket(catecholate)).block
    assert len(product.open_vacancies()) == len(metal.open_vacancies()) - 2
    bound = [s for s in product.sites if s.donor_type == "phenolate_O"]
    assert len(bound) == 2
    assert all(product.state[BuildingBlock.state_key(s)].status is SiteStatus.OCCUPIED
               for s in bound)


def test_a_chelate_that_cannot_reach_a_trans_pair_says_so(catecholate):
    """D10's case from the assembly side: a pocket that spans cis cannot span trans."""
    metal = metal_block(cn=6)
    with pytest.raises(IncompatibleJoin) as exc:
        join_chelate(metal, catecholate, vertex_pair(metal, trans=True),
                     pocket(catecholate))
    assert "cannot reach" in str(exc.value)
    assert exc.value.verdict.strain > 1.0


def test_the_choice_vector_records_both_ends_and_no_torsion_well(catecholate):
    """The roll is DETERMINED by the second contact, so there is no well to record."""
    metal = metal_block(cn=6)
    cv = join_chelate(metal, catecholate, vertex_pair(metal, trans=False),
                      pocket(catecholate)).choice_vector
    assert cv["mode"] == "chelate" and "torsion_well" not in cv
    assert [d["atom"] for d in cv["donors"]] == [s.atom_idx for s in pocket(catecholate)]
    assert len(cv["vacancies"]) == 2 and len(cv["d_ml"]) == 2


def test_a_chelate_join_needs_two_donors_and_two_vertices(catecholate):
    metal = metal_block(cn=6)
    with pytest.raises(IncompatibleJoin, match="two donors from one block"):
        join_chelate(metal, catecholate, vertex_pair(metal, trans=False),
                     [pocket(catecholate)[0], metal.open_vacancies()[0]])
    with pytest.raises(ValueError, match="two sites on each block"):
        join_chelate(metal, catecholate, metal.open_vacancies()[:3], pocket(catecholate))


def test_either_donor_order_reaches_one_node(catecholate):
    """Which donor took which vertex is a route, not an identity (D2)."""
    metal = metal_block(cn=6)
    pair = vertex_pair(metal, trans=False)
    donors = pocket(catecholate)
    forward = join_chelate(metal, catecholate, pair, donors).block
    backward = join_chelate(metal, catecholate, pair, donors[::-1]).block
    assert l1_graph_hash(forward.graph) == l1_graph_hash(backward.graph)
