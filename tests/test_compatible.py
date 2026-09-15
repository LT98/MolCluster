"""M5/S2: can these two open sites be joined, and at what strain?

Written as a verdict table, the same way the M2 discrimination tests are: pairs that must
join, pairs that must not, and — because an opaque `False` is what makes an assembly step
undebuggable three milestones later — an assertion on WHY in every negative case.

The headline pair is the last section. A chelate that spans a cis vertex pair comfortably
cannot reach a trans one, which is the D10 anthrarufin distinction seen from the assembly
side: cis and trans differ in what can be built next, not merely in what they are called.
"""
from __future__ import annotations

import pytest

from mofsbu.assembly.join import (
    MAX_BITE_MISMATCH_DEG, chelate_compatible, compatible)
from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.frames import BindingMode, binding_modes, live_dof
from mofsbu.sites.model import Site, perceive, vacancy_sites


def donor(donor_type: str, atom_idx: int = 1, frame: dict | None = None) -> Site:
    """A donor site with the real tables behind it and a frame only where one matters."""
    return Site(atom_idx=atom_idx, donor_type=donor_type, labile=True, charge_after=0,
                live_dof=live_dof(donor_type).value,
                binding_modes=tuple(m.value for m in binding_modes(donor_type)),
                frame=frame)


def vertices(*directions, metal_idx: int = 0) -> list[Site]:
    return vacancy_sites(metal_idx, [0.0, 0.0, 0.0], list(directions))


CIS = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))        # 90 deg apart
TRANS = ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0))     # 180 deg apart


@pytest.fixture(scope="module")
def acetate_donors():
    """A real carboxylate's two oxygens, carrying the frames perception gave them."""
    mol = embed_molecule(mol_from_smiles("CC(=O)[O-]"))
    sites = [s for s in perceive(mol) if s.donor_type == "carboxylate_O"]
    assert len(sites) == 2, "acetate should perceive as one group of two equivalent O"
    return sites


# ── one donor, one vertex ────────────────────────────────────────────────────

def test_a_donor_goes_on_a_vacant_vertex():
    verdict = compatible(donor("carboxylate_O"), vertices(*CIS)[0])
    assert verdict.feasible and verdict.strain == 0.0
    assert verdict.mode == BindingMode.MONODENTATE.value


def test_a_single_point_join_reports_no_strain_rather_than_inventing_one():
    """The blocks are in different coordinate systems; the rigid move between them is free.

    Any number here would be about where the two blocks happen to sit, which is not a fact
    about whether they can be joined.
    """
    verdict = compatible(donor("pyridyl_N"), vertices(*CIS)[0])
    assert verdict.strain == 0.0
    assert "rigid move" in verdict.reason


def test_a_live_torsion_branches_and_a_free_one_does_not():
    """The well count IS the number of L3 siblings the bond will generate (D11)."""
    assert compatible(donor("carboxylate_O"), vertices(*CIS)[0]).wells == (0.0, 180.0)
    assert compatible(donor("aqua_O"), vertices(*CIS)[0]).wells == (0.0,)


def test_two_donors_do_not_bond_to_each_other():
    verdict = compatible(donor("carboxylate_O"), donor("pyridyl_N", atom_idx=4))
    assert not verdict.feasible
    assert "two donors" in verdict.reason and "vacancy" in verdict.reason


def test_two_vertices_are_a_metal_metal_bond_and_not_a_join():
    a, b = vertices(*CIS)
    verdict = compatible(a, b)
    assert not verdict.feasible
    assert "metal-metal" in verdict.reason


def test_a_mode_the_donor_does_not_offer_is_refused_by_name():
    verdict = compatible(donor("pyridyl_N"), vertices(*CIS)[0],
                         mode=BindingMode.CHELATE.value)
    assert not verdict.feasible
    assert "pyridyl_N" in verdict.reason and "mono" in verdict.reason


def test_one_vertex_cannot_host_a_chelate():
    """A chelate needs two vertices, so this is a different call — and it says so."""
    verdict = compatible(donor("carboxylate_O"), vertices(*CIS)[0],
                         mode=BindingMode.CHELATE.value)
    assert not verdict.feasible
    assert "vacancy does not bind" in verdict.reason


def test_the_distance_is_a_property_of_the_pair():
    """And the element comes from the caller, never from the donor type's name."""
    o = compatible(donor("carboxylate_O"), vertices(*CIS)[0],
                   partner="Cu", donor_element="O")
    n = compatible(donor("pyridyl_N"), vertices(*CIS)[0],
                   partner="Cu", donor_element="N")
    assert o.d_ml is not None and n.d_ml is not None and o.d_ml != n.d_ml
    # No partner named, no distance invented.
    assert compatible(donor("carboxylate_O"), vertices(*CIS)[0]).d_ml is None


def test_a_frameless_donor_still_answers_but_says_what_it_could_not_use():
    verdict = compatible(donor("carboxylate_O"), vertices(*CIS)[0])
    assert verdict.feasible
    assert "no frame" in verdict.reason


def test_a_fallback_frame_is_flagged_and_not_refused():
    """An inferred axis is a real answer from a donor with nothing to orient it."""
    guessed = donor("aqua_O", frame={"origin": [0.0, 0.0, 0.0], "axis": [0.0, 0.0, 1.0],
                                     "ref": [1.0, 0.0, 0.0], "mode": "fallback"})
    verdict = compatible(guessed, vertices(*CIS)[0])
    assert verdict.feasible
    assert "fallback" in verdict.reason


# ── two donors, two vertices: where the geometry starts to bite ──────────────

def test_a_chelate_spans_cis_and_cannot_reach_trans(acetate_donors):
    """The headline: cis and trans differ in what can be built next, not just in name."""
    cis = chelate_compatible(acetate_donors, vertices(*CIS))
    trans = chelate_compatible(acetate_donors, vertices(*TRANS))
    assert cis.feasible and not trans.feasible
    assert cis.strain < trans.strain
    assert "cannot reach" in trans.reason


def test_feasible_is_exactly_strain_within_one(acetate_donors):
    """Strain is the mismatch over the tolerance, so the policy lives in one constant."""
    for vac in (vertices(*CIS), vertices(*TRANS)):
        verdict = chelate_compatible(acetate_donors, vac)
        assert verdict.feasible == (verdict.strain <= 1.0)


def test_the_bite_angle_populations_stay_far_apart(acetate_donors):
    """Guards the calibration behind `MAX_BITE_MISMATCH_DEG` rather than the threshold.

    Acetate is the strained end of what really chelates. If its mismatch against a cis
    pair ever climbs to within reach of its mismatch against a trans pair, the constant
    is no longer separating two populations and needs re-measuring — which is a different
    and louder failure than one case flipping.
    """
    cis = chelate_compatible(acetate_donors, vertices(*CIS)).strain * MAX_BITE_MISMATCH_DEG
    trans = chelate_compatible(acetate_donors, vertices(*TRANS)).strain * MAX_BITE_MISMATCH_DEG
    assert cis < 35.0 < 85.0 < trans


def test_two_vertices_on_different_metals_are_a_bridge_not_a_chelate(acetate_donors):
    split = [vertices((1.0, 0.0, 0.0))[0], vertices((0.0, 1.0, 0.0), metal_idx=9)[0]]
    verdict = chelate_compatible(acetate_donors, split)
    assert not verdict.feasible
    assert "mu2 bridge" in verdict.reason


def test_a_donor_type_that_does_not_chelate_is_refused(acetate_donors):
    pyridyls = [donor("pyridyl_N", atom_idx=i, frame=acetate_donors[0].frame)
                for i in (1, 2)]
    verdict = chelate_compatible(pyridyls, vertices(*CIS))
    assert not verdict.feasible and "does not chelate" in verdict.reason


def test_one_atom_cannot_chelate_by_itself(acetate_donors):
    twice = [acetate_donors[0], acetate_donors[0]]
    verdict = chelate_compatible(twice, vertices(*CIS))
    assert not verdict.feasible and "by itself" in verdict.reason


def test_a_missing_frame_raises_instead_of_reading_as_impossible(acetate_donors):
    """`unknown` reported as `no` is the failure mode `open_sites` already refuses."""
    with pytest.raises(ValueError, match="unknown answer, not a negative"):
        chelate_compatible([acetate_donors[0], donor("carboxylate_O", atom_idx=7)],
                           vertices(*CIS))


def test_a_bridge_is_not_asked_of_this_function(acetate_donors):
    with pytest.raises(ValueError, match="mu2/mu3"):
        chelate_compatible(acetate_donors, vertices(*CIS)[:1])
