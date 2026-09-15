"""M5/S6: L3 is provenance-primary and geometry only verifies (D11).

The ordering is the design, so it is what the tests are mostly about: a choice vector
decides identity, and geometry is allowed to collapse duplicates and reconcile
divergence/convergence — never to name anything. A clustering pass that ran first could
merge two branches a builder deliberately made, which is the mistake D11 exists to prevent.

`theta_geom` is passed explicitly in every clustering test below, so none of them depends on
the default. The default's own calibration is guarded separately, at the bottom.
"""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu.assembly.construct import ligand_block, metal_block
from mofsbu.assembly.join import join
from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.identity.conformers import (
    CALIBRATION_KIND_A_MAX, CALIBRATION_KIND_B_MIN, DEFAULT_ENERGY_WINDOW,
    DEFAULT_THETA_GEOM, Conformer, cluster, core_atoms, core_rmsd, l3_conformer_id)
from mofsbu.identity.keys import l1_graph_hash
from mofsbu.identity.keys import l3_conformer_id as key_l3

CV_A = {"op": "join", "mode": "mono", "torsion_well": 0, "vertex": 0}
CV_B = {"op": "join", "mode": "mono", "torsion_well": 1, "vertex": 0}


@pytest.fixture(scope="module")
def built():
    """Four real products: two choice vectors, each sampled at two embedding seeds."""
    out = {}
    for label, slot in (("v0", 0), ("v1", 1)):
        for seed in (0xC0FFEE, 11):
            aqua = ligand_block(embed_molecule(mol_from_smiles("O"), seed=seed),
                                charge=0, name="aqua")
            pyridine = ligand_block(embed_molecule(mol_from_smiles("c1ccncc1"), seed=seed),
                                    charge=0, name="pyridine")
            block = metal_block("Fe", 3, cn=6, spin_class="hs")
            for ligand, vertex_slot in ((aqua, 5), (pyridine, slot)):
                vertex = next(s for s in block.open_vacancies() if s.slot == vertex_slot)
                block = join(block, ligand, vertex, ligand.open_donors()[0]).block
            out[(label, seed)] = block
    return out


def conformer(block, cv) -> Conformer:
    return Conformer(choice_vector=cv, graph=block.graph, geometry=block.geometry)


# ── the label ────────────────────────────────────────────────────────────────

def test_the_label_is_the_choice_vector_and_nothing_else():
    assert l3_conformer_id(CV_A) == l3_conformer_id(dict(reversed(list(CV_A.items()))))
    assert l3_conformer_id(CV_A) != l3_conformer_id(CV_B)


def test_the_label_does_not_move_when_the_coordinates_do():
    """Provenance-primary, stated as a property: geometry is not an input to the name."""
    assert l3_conformer_id(CV_A, geom=np.zeros((3, 3))) == l3_conformer_id(CV_A)


def test_no_provenance_means_no_l3_identity():
    """Correct rather than degraded: a provenance-based key needs a provenance."""
    assert l3_conformer_id(None) == ""
    assert l3_conformer_id({}) == ""


def test_the_key_module_delegates_here():
    assert key_l3(CV_A) == l3_conformer_id(CV_A)


# ── what gets compared ───────────────────────────────────────────────────────

def test_the_core_is_the_rigid_part_and_the_coordination_sphere(built):
    """Not the whole molecule: a floppy tail gives a big number and means nothing (§4.2)."""
    block = built[("v0", 0xC0FFEE)]
    core = core_atoms(block.graph)
    assert block.graph.metals()[0] in core
    assert len(core) < len(block.graph)                 # something was excluded
    donors = block.graph.neighbors(block.graph.metals()[0], __import__(
        "mofsbu.graph._types", fromlist=["EdgeType"]).EdgeType.DATIVE)
    assert set(donors) <= set(core)                     # every donor is in it


def test_rmsd_is_zero_for_a_structure_against_itself(built):
    block = built[("v0", 0xC0FFEE)]
    assert core_rmsd(block.graph, block.geometry,
                     block.graph, block.geometry) == pytest.approx(0.0, abs=1e-9)


def test_rmsd_ignores_where_the_structure_happens_to_sit(built):
    """Best-fit superposition, so a translation or rotation is not a difference."""
    block = built[("v0", 0xC0FFEE)]
    rotated = block.geometry @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]]).T + 17.0
    assert core_rmsd(block.graph, block.geometry,
                     block.graph, rotated) == pytest.approx(0.0, abs=1e-6)


# ── clustering: collapse, converge, diverge ──────────────────────────────────

def test_samples_of_one_choice_vector_are_one_conformer(built):
    """Kind-A collapse — and it needs no threshold at RAW, because the spread is zero."""
    group = [conformer(built[("v0", seed)], CV_A) for seed in (0xC0FFEE, 11)]
    clusters = cluster(group, theta_geom=0.1)
    assert len(clusters) == 1
    assert clusters[0].size == 2
    assert clusters[0].merged_by == "choice-vector"


def test_the_join_absorbs_the_embeddings_orientation(built):
    """What a frame-determined join buys, stated at the magnitude it actually delivers.

    Before `_place_donor_block` used the frame's `ref` it used only the axis, so the roll
    about the new bond was settled by `rotation_between`'s minimal rotation — which depends
    on how the ligand happened to be oriented in its own coordinates. Re-embedding rotates
    a ligand rigidly, so the SAME choice vector produced cores **1.56 A** apart.

    It is not zero now, and the residual is worth naming rather than tolerating: ~6e-06 A
    of the ligand's OWN internal geometry, which MMFF does not converge bit-identically
    from different starts. That is five orders of magnitude below the placement defect it
    replaced, and it is the only Kind-A noise the RAW construct path has.
    """
    a, b = built[("v0", 0xC0FFEE)], built[("v0", 11)]
    rmsd = core_rmsd(a.graph, a.geometry, b.graph, b.geometry)
    assert rmsd < 1e-4
    assert rmsd < 1.56 / 10_000                # against the defect this replaced


def test_two_choice_vectors_that_reach_one_geometry_are_reconciled(built):
    """Convergence: many choices, one well. The only merge geometry is allowed to make."""
    block = built[("v0", 0xC0FFEE)]
    group = [conformer(block, CV_A), conformer(block, CV_B)]
    clusters = cluster(group, theta_geom=0.1)
    assert len(clusters) == 1
    assert clusters[0].merged_by == "geometry"
    assert "converged onto this one" in clusters[0].reasons[0]


def test_two_choice_vectors_at_different_geometries_stay_apart(built):
    """The D10 case in L3's clothing: different branches are different conformers."""
    group = [conformer(built[("v0", 0xC0FFEE)], CV_A),
             conformer(built[("v1", 0xC0FFEE)], CV_B)]
    assert l1_graph_hash(group[0].graph) == l1_graph_hash(group[1].graph)
    assert len(cluster(group, theta_geom=0.1)) == 2


def test_a_near_degenerate_pair_is_not_merged_by_a_threshold(built):
    """Provenance beats geometry: distinct choice vectors survive any theta_geom.

    Turn the threshold up until it would swallow the whole structure and the two still do
    not merge, because the first pass never consults geometry at all... and the second pass
    only ever merges what the first pass left, which for identical coordinates is the
    convergence case tested above. Here the coordinates genuinely differ.
    """
    group = [conformer(built[("v0", 0xC0FFEE)], CV_A),
             conformer(built[("v1", 0xC0FFEE)], CV_A)]       # SAME choice vector
    assert len(cluster(group, theta_geom=1e-6)) == 2         # diverged: two wells
    assert len(cluster(group, theta_geom=99.0)) == 1         # one well at a huge threshold


def test_one_choice_vector_in_two_wells_gets_two_labels(built):
    """Divergence: same provenance, different minima. The suffix disambiguates them."""
    group = [conformer(built[("v0", 0xC0FFEE)], CV_A),
             conformer(built[("v1", 0xC0FFEE)], CV_A)]
    clusters = cluster(group, theta_geom=1e-6)
    assert {c.label for c in clusters} == {l3_conformer_id(CV_A) + "#0",
                                           l3_conformer_id(CV_A) + "#1"}
    assert "diverged into 2 wells" in clusters[0].reasons[-1]


def test_an_energy_window_stops_a_bad_geometry_hiding_in_a_good_ones_identity(built):
    """Similar shape, very different energy: not a duplicate, a worse structure."""
    block = built[("v0", 0xC0FFEE)]
    near = Conformer(CV_A, block.graph, block.geometry, energy=-100.0)
    far = Conformer(CV_B, block.graph, block.geometry, energy=-40.0)
    assert len(cluster([near, far], theta_geom=0.1, energy_window=5.0)) == 2
    assert len(cluster([near, far], theta_geom=0.1, energy_window=100.0)) == 1


def test_a_missing_energy_is_not_evidence(built):
    """D18's rule applied here: an absent number never gates anything out."""
    block = built[("v0", 0xC0FFEE)]
    scored = Conformer(CV_A, block.graph, block.geometry, energy=-100.0)
    unscored = Conformer(CV_B, block.graph, block.geometry, energy=None)
    assert len(cluster([scored, unscored], theta_geom=0.1, energy_window=1.0)) == 1


def test_candidates_without_geometry_are_grouped_by_provenance_alone(built):
    """No coordinates, no verifier — and that is the whole answer, not half of one."""
    group = [Conformer(CV_A), Conformer(CV_A), Conformer(CV_B)]
    clusters = cluster(group, theta_geom=0.1)
    assert sorted(c.size for c in clusters) == [1, 2]
    assert all(c.merged_by == "choice-vector" for c in clusters)


# ── C2: the threshold's calibration, guarded ─────────────────────────────────

def test_theta_geom_sits_in_the_valley_it_was_calibrated_from():
    """The C2 call, pinned so it cannot drift away from the data that set it.

    Not a re-measurement — that needs 12 xTB relaxations and lives in the work plan. This
    asserts the number stays consistent with the populations it was placed between, so
    editing it without revisiting the calibration fails here rather than silently changing
    what "the same conformer" means.
    """
    assert CALIBRATION_KIND_A_MAX < DEFAULT_THETA_GEOM < CALIBRATION_KIND_B_MIN
    above = DEFAULT_THETA_GEOM / CALIBRATION_KIND_A_MAX
    below = CALIBRATION_KIND_B_MIN / DEFAULT_THETA_GEOM
    assert above > 3.0 and below > 3.0            # margin on BOTH sides, not just one
    assert abs(above - below) < 1.0               # ...and roughly balanced: a geometric mean


def test_the_energy_window_is_deliberately_unset():
    """Half of C2 is called and half is not, and the code says which is which.

    A default of 0.0 would silently reject every candidate; a plausible-looking number would
    bake in the rigid-core defect that produced the 16.6 kcal/mol Kind-A spread it would
    have been calibrated from.
    """
    assert DEFAULT_ENERGY_WINDOW is None
