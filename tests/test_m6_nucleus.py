"""M6/S4: the declared nucleus, and which vertices a centre leaves open.

Two halves of one idea — the caller says what the starting geometry leaves open — and
they are tested together because they fail together.  A nucleus declared without reserved
vertices is a dimer nothing can bridge, and reserved vertices with no way to bond two
metals are an arrangement for a partner that cannot exist.

The distance is the part worth stating plainly.  Every other join in this package reports
M...M as an OUTPUT of ligand geometry (D20), which is what makes a bridged dimer's
2.673 A a measurement rather than an assumption.  This is the one route with no ligand to
derive it from, so the number is DECLARED — and `metal_metal_distance` refuses to supply
a motif-specific one, because the element symbols do not distinguish a quadruply-bonded
Mo2 at 2.09 A from a Cu2 paddlewheel at 2.62 A.

The reserved-vertex case is the defect the pathway ladder hit: a CN-6 centre carrying
four co-ligands came back with its two vacancies *trans*, where a ~90 deg chelate cannot
reach them, and the step refused with `chelate_cannot_span`.  Both arrangements are
asserted here, so "cis" is a measured angle rather than a claim about index conventions.
"""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu._types import AmbiguousSpecError
from mofsbu.assembly.construct import (
    enumerate_constructions, ligand_block, metal_block)
from mofsbu.assembly.join import (
    IncompatibleJoin, compatible, join_metal_metal, metal_metal_compatible)
from mofsbu.geometry._linalg import angle_between
from mofsbu.geometry.distances import metal_metal_distance
from mofsbu.geometry.embed import embed_molecule
from mofsbu.geometry.placer import (
    LigandPlacement, cis_vertices, place_mononuclear, site_vectors)
from mofsbu.graph._types import EdgeType
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.state import SiteStatus

#: Literature-typical for the Cu2 paddlewheel, and DECLARED here rather than looked up:
#: that is the whole point of the route under test.
CU_CU = 2.62


def cu(cn: int = 6) -> object:
    return metal_block("Cu", 2, cn=cn, geometry="octahedral")


def vacancies_of(block) -> list:
    return [s for s in block.sites if s.is_vacancy]


def aqua() -> LigandPlacement:
    return LigandPlacement(mol=embed_molecule(mol_from_smiles("O"), seed=3),
                           donor_idxs=(0,), donor_types=("aqua_O",),
                           mode=BindingMode.MONODENTATE, name="H2O")


@pytest.fixture(scope="module")
def dimer():
    """A declared Cu2 nucleus: two octahedral Cu(II), bonded vertex to vertex."""
    a, b = cu(), cu()
    return join_metal_metal(a, b, vacancies_of(a)[0], vacancies_of(b)[0],
                            d_mm=CU_CU).block


# --- the distance is declared, and nothing invents one -----------------------------

def test_a_declared_distance_is_used_exactly_and_says_it_was_declared():
    d = metal_metal_distance("Cu", "Cu", override=CU_CU)
    assert d.value == CU_CU
    assert d.source == "override"
    assert not d.estimated


def test_without_a_declaration_the_distance_is_a_covalent_estimate_that_admits_it():
    d = metal_metal_distance("Cu", "Cu")
    assert d.value == pytest.approx(2.64)
    assert d.source == "covalent-radii"
    assert d.estimated, "an estimate that does not report itself as one is worse than none"


def test_a_metal_metal_bond_is_not_lengthened_like_a_dative_one():
    from mofsbu.geometry.distances import COVALENT_RADII, DATIVE_LENGTHENING

    d = metal_metal_distance("Cu", "Cu")
    assert d.value == pytest.approx(2 * COVALENT_RADII["Cu"])
    assert d.value != pytest.approx(2 * COVALENT_RADII["Cu"] + DATIVE_LENGTHENING)


def test_naming_a_motif_is_refused_because_there_is_no_table_to_look_it_up_in():
    with pytest.raises(AmbiguousSpecError) as exc:
        metal_metal_distance("Mo", "Mo", motif="paddlewheel")
    message = str(exc.value)
    assert "paddlewheel" in message
    # The refusal carries the two numbers that make the point, and says what to do.
    assert "2.09" in message and "2.62" in message
    assert "override" in message


def test_a_declaration_outranks_a_motif_rather_than_being_refused_alongside_it():
    d = metal_metal_distance("Mo", "Mo", motif="paddlewheel", override=2.09)
    assert (d.value, d.source) == (2.09, "override")


# --- vacancy to vacancy is a metal-metal join -------------------------------------

def test_the_monodentate_verdict_refuses_two_vertices_and_names_what_does_join_them():
    a, b = cu(), cu()
    verdict = compatible(vacancies_of(a)[0], vacancies_of(b)[0])
    assert not verdict.feasible
    assert "join_metal_metal" in verdict.reason


def test_two_vertices_on_two_blocks_are_a_feasible_metal_metal_pair():
    a, b = cu(), cu()
    verdict = metal_metal_compatible(vacancies_of(a)[0], vacancies_of(b)[0],
                                     elements=("Cu", "Cu"), d_mm=CU_CU)
    assert verdict.feasible
    assert verdict.strain == 0.0
    assert verdict.d_ml == CU_CU
    assert verdict.mode == EdgeType.METAL_METAL.value


def test_an_estimated_separation_is_feasible_but_says_so_in_its_verdict():
    a, b = cu(), cu()
    verdict = metal_metal_compatible(vacancies_of(a)[0], vacancies_of(b)[0],
                                     elements=("Cu", "Cu"))
    assert verdict.feasible
    assert "ESTIMATE" in verdict.reason


def test_a_donor_is_refused_by_the_metal_metal_verdict_and_pointed_back_at_join():
    a, b = cu(), cu()
    ligand = ligand_block(embed_molecule(mol_from_smiles("O"), seed=3), charge=0)
    donor = next(s for s in ligand.sites if not s.is_vacancy)
    verdict = metal_metal_compatible(donor, vacancies_of(b)[0], elements=("O", "Cu"))
    assert not verdict.feasible
    assert "`join`" in verdict.reason


def test_the_nucleus_is_two_metals_one_mm_edge_and_the_distance_that_was_asked_for(dimer):
    assert len(dimer.graph) == 2
    edges = [(u, v, e) for u, v, e in dimer.graph.edges()]
    assert len(edges) == 1
    u, v, etype = edges[0]
    assert etype is EdgeType.METAL_METAL, "identity has to SEE an M-M bond as one"
    xyz = np.asarray(dimer.geometry)
    assert float(np.linalg.norm(xyz[0] - xyz[1])) == pytest.approx(CU_CU)


def test_each_metal_spends_the_vertex_that_faces_the_other(dimer):
    open_vertices = vacancies_of(dimer)
    assert len(open_vertices) == 10, "a CN-6 pair has 12 vertices; the join consumes 2"
    per_atom = {i: sum(1 for s in open_vertices if s.atom_idx == i) for i in (0, 1)}
    assert per_atom == {0: 5, 1: 5}


def test_the_moved_block_reports_its_vertices_in_the_dimers_coordinate_system(dimer):
    """The anchor never moves, so the partner's frames must have travelled with it.

    A frame left in the block's own coordinates would point at where that metal used to
    be, and the next join would be judged against a stale direction.
    """
    xyz = np.asarray(dimer.geometry)
    for site in vacancies_of(dimer):
        origin = np.asarray(site.frame["origin"], dtype=float)
        assert origin == pytest.approx(xyz[site.atom_idx]), (
            "a vacancy's origin IS its metal's position")
    partner = [s for s in vacancies_of(dimer) if s.atom_idx == 1]
    assert all(np.asarray(s.frame["origin"])[0] == pytest.approx(CU_CU) for s in partner)


def test_no_vertex_of_the_dimer_is_occluded_by_the_other_metal(dimer):
    """B8's occlusion rule is the thing that could quietly make this route useless."""
    for site in vacancies_of(dimer):
        state = dimer.state[(site.atom_idx, site.slot)]
        assert state.status is SiteStatus.OPEN, (
            f"vertex {site.slot} on atom {site.atom_idx} came back {state.status}")


def test_the_choice_vector_records_the_distance_and_where_the_number_came_from():
    a, b = cu(), cu()
    declared = join_metal_metal(a, b, vacancies_of(a)[0], vacancies_of(b)[0], d_mm=CU_CU)
    assert declared.choice_vector["d_mm"] == CU_CU
    assert declared.choice_vector["d_mm_source"] == "override"
    c, d = cu(), cu()
    estimated = join_metal_metal(c, d, vacancies_of(c)[0], vacancies_of(d)[0])
    assert estimated.choice_vector["d_mm_source"] == "covalent-radii"


def test_a_block_cannot_bond_to_itself(dimer):
    a = cu()
    sites = vacancies_of(a)
    with pytest.raises(IncompatibleJoin) as exc:
        join_metal_metal(a, a, sites[0], sites[1])
    assert "same block" in str(exc.value)


def test_a_roll_about_the_new_axis_leaves_the_separation_alone():
    a, b = cu(), cu()
    straight = join_metal_metal(a, b, vacancies_of(a)[0], vacancies_of(b)[0], d_mm=CU_CU)
    c, d = cu(), cu()
    rolled = join_metal_metal(c, d, vacancies_of(c)[0], vacancies_of(d)[0],
                              d_mm=CU_CU, roll_deg=45.0)
    for result in (straight, rolled):
        xyz = np.asarray(result.block.geometry)
        assert float(np.linalg.norm(xyz[0] - xyz[1])) == pytest.approx(CU_CU)
    assert rolled.choice_vector["roll_deg"] == 45.0


# --- and the nucleus is something the assembly path can grow on --------------------

def test_enumerate_constructions_grows_onto_both_metals_of_the_dimer(dimer):
    water = ligand_block(embed_molecule(mol_from_smiles("O"), seed=3), charge=0,
                         name="H2O")
    tree = enumerate_constructions(dimer, [water], degree=1)
    assert not tree.capped
    assert tree.refusals == {}
    assert len(tree) == len(vacancies_of(dimer)), (
        "one product per open vertex, or the dimer's vertices are not all reachable")
    assert {c.steps[-1]["vacancy"]["atom"] for c in tree.leaves} == {0, 1}


def test_the_metal_metal_bond_survives_growing_a_ligand_onto_the_nucleus(dimer):
    water = ligand_block(embed_molecule(mol_from_smiles("O"), seed=3), charge=0,
                         name="H2O")
    tree = enumerate_constructions(dimer, [water], degree=2)
    assert len(tree) > 0
    for leaf in tree.leaves:
        xyz = np.asarray(leaf.block.geometry)
        assert float(np.linalg.norm(xyz[0] - xyz[1])) == pytest.approx(CU_CU)
        assert any(e is EdgeType.METAL_METAL for _u, _v, e in leaf.block.graph.edges())


# --- which vertices stay open ------------------------------------------------------

def test_cis_vertices_measures_the_angle_rather_than_trusting_the_table_order():
    for geometry, cn in (("octahedral", 6), ("square_planar", 4),
                         ("square_pyramidal", 5), ("trigonal_bipyramidal", 5)):
        picked = cis_vertices(geometry, cn, 2)
        targets = site_vectors(geometry, cn, 1.0)
        assert angle_between(targets[picked[0]], targets[picked[1]]) == pytest.approx(90.0)


def test_a_trigonal_bipyramid_picks_an_axial_equatorial_pair_not_two_equatorials():
    """120 deg apart is the closest pair a table-order answer would have returned."""
    assert cis_vertices("trigonal_bipyramidal", 5, 2) == (0, 3)


def test_a_facial_triple_is_mutually_cis():
    picked = cis_vertices("octahedral", 6, 3)
    targets = site_vectors("octahedral", 6, 1.0)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        assert angle_between(targets[picked[i]], targets[picked[j]]) == pytest.approx(90.0)


def test_four_co_ligands_leave_their_vacancies_trans_when_nobody_says_otherwise():
    result = place_mononuclear("Cu", [aqua() for _ in range(4)],
                               geometry="octahedral", cn=6)
    assert result.ok
    v = np.asarray(result.vacancies)
    assert len(v) == 2
    assert angle_between(v[0], v[1]) == pytest.approx(180.0), (
        "this trans pair is the measured defect `reserve` exists to fix")


def test_reserving_a_cis_pair_returns_one():
    reserved = cis_vertices("octahedral", 6, 2)
    result = place_mononuclear("Cu", [aqua() for _ in range(4)],
                               geometry="octahedral", cn=6, reserve=reserved)
    assert result.ok
    v = np.asarray(result.vacancies)
    assert angle_between(v[0], v[1]) == pytest.approx(90.0)


def test_the_reserved_vertices_are_the_ones_that_come_back_vacant():
    reserved = (2, 4)
    result = place_mononuclear("Cu", [aqua() for _ in range(4)],
                               geometry="octahedral", cn=6, reserve=reserved)
    targets = site_vectors("octahedral", 6, 1.0)
    got = {tuple(np.round(x, 6)) for x in result.vacancies}
    assert got == {tuple(np.round(targets[i], 6)) for i in reserved}


def test_reserving_nothing_leaves_the_choice_vector_exactly_as_it_was():
    """An argument that exists must not re-label structures built before it did."""
    plain = place_mononuclear("Cu", [aqua() for _ in range(4)],
                              geometry="octahedral", cn=6)
    explicit = place_mononuclear("Cu", [aqua() for _ in range(4)],
                                 geometry="octahedral", cn=6, reserve=())
    assert "reserve" not in plain.choice_vector
    assert "reserve" not in explicit.choice_vector
    assert plain.choice_vector == explicit.choice_vector


def test_a_reservation_travels_in_the_choice_vector_because_it_moves_the_ligands():
    result = place_mononuclear("Cu", [aqua() for _ in range(4)],
                               geometry="octahedral", cn=6, reserve=(4, 2))
    assert result.choice_vector["reserve"] == [2, 4]


def test_reserving_more_vertices_than_the_ligands_can_spare_is_refused_with_the_count():
    with pytest.raises(ValueError) as exc:
        place_mononuclear("Cu", [aqua() for _ in range(4)],
                          geometry="octahedral", cn=6, reserve=(0, 1, 2))
    assert "reserved" in str(exc.value)


def test_a_vertex_that_is_not_on_the_polyhedron_is_refused_rather_than_dropped():
    with pytest.raises(ValueError) as exc:
        place_mononuclear("Cu", [aqua()], geometry="octahedral", cn=6, reserve=(6,))
    assert "0..5" in str(exc.value)
