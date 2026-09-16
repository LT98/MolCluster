"""The lone-pair lobe as a branch — M6/S1, mechanism A's enabling step.

A bridge is the first operation that cares *which* lone pair of an sp2 donor binds.
`sites.frames.site_frame` has always known there are two and says so in its own comment;
`sites.model.perceive` kept the first and discarded the rest, so every stored site pointed
its metal at the syn lobe and no bridging geometry was reachable from the assembly path.

The numbers here are the calibration, not decoration: the four lobe combinations of a
formate's two oxygens are the three textbook carboxylate bridging modes, and they are
separated by more than a bond length.  A test that only asserted "there are two lobes"
would pass on two lobes that pointed the same way.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from mofsbu.assembly.join import _transform_frame
from mofsbu.geometry._linalg import axis_rotation
from mofsbu.geometry.distances import metal_donor_distance
from mofsbu.sites.frames import BindingMode, LiveDOF, live_dof, lone_pair_frames, torsion_wells
from mofsbu.sites.model import frame_lobes, perceive
from mofsbu.graph.from_mol import mol_from_smiles


# Literature, for the three modes the lobes select between.  A syn-syn carboxylate is what
# holds a paddlewheel together at ~2.6 A; anti-anti is the extended bridge of a chain.
PADDLEWHEEL_CU_CU = 2.62


@pytest.fixture(scope="module")
def formate():
    mol = mol_from_smiles("[O-]C=O", embed=True)
    oxygens = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O"]
    return mol, oxygens


def implied_separation(mol, oxygens, wells, metal="Cu") -> float:
    """How far apart the two metals a bridging ligand implies would sit.

    Each donor's frame says the metal lies at `origin + d * axis`.  For two donors that is
    two metal positions, and their separation is the M...M distance the bridge dictates —
    the quantity that distinguishes the bridging modes and the one a paddlewheel's identity
    turns on.
    """
    conf = mol.GetConformer()
    d = metal_donor_distance(metal, "O").value
    points = []
    for idx, well in zip(oxygens, wells):
        frame = lone_pair_frames(mol, idx, "carboxylate_O", conf)[well]
        points.append(np.array(frame.origin) + d * np.array(frame.axis_hat))
    return float(np.linalg.norm(points[0] - points[1]))


# ── the lobes exist, and they point somewhere different ──────────────────────

def test_an_sp2_donor_offers_two_lone_pairs_and_a_determined_one_offers_one(formate):
    """The tuple length is the answer to 'can this atom bridge on its own'."""
    mol, oxygens = formate
    conf = mol.GetConformer()
    assert len(lone_pair_frames(mol, oxygens[0], "carboxylate_O", conf)) == 2

    # Water's oxygen has two hydrogens, so its bonding determines one direction.  That is
    # why a mu2-aqua bridge is not expressible by this model and is M6's mechanism B.
    water = mol_from_smiles("O", embed=True)
    o = next(a.GetIdx() for a in water.GetAtoms() if a.GetSymbol() == "O")
    assert len(lone_pair_frames(water, o, "aqua_O", water.GetConformer())) == 1


def test_the_two_lobes_are_not_the_same_direction(formate):
    """Guards the collapse rule: two lobes that agreed would be one lobe reported twice."""
    mol, oxygens = formate
    lobes = lone_pair_frames(mol, oxygens[0], "carboxylate_O", mol.GetConformer())
    cos = float(np.dot(np.array(lobes[0].axis_hat), np.array(lobes[1].axis_hat)))
    assert cos < 0.9, f"the two lone pairs point {np.degrees(np.arccos(cos)):.0f} deg apart"


# ── the calibration: the lobes ARE the bridging modes ────────────────────────

@pytest.mark.parametrize("wells, expected, mode", [
    ((1, 1), 2.673, "syn-syn"),
    ((0, 1), 5.148, "syn-anti"),
    ((1, 0), 5.148, "anti-syn"),
    ((0, 0), 5.516, "anti-anti"),
])
def test_the_lobe_combinations_are_the_three_bridging_modes(formate, wells, expected, mode):
    mol, oxygens = formate
    assert implied_separation(mol, oxygens, wells) == pytest.approx(expected, abs=1e-3), mode


def test_one_combination_is_the_paddlewheel_and_the_others_are_nowhere_near_it(formate):
    """The separation that matters, and the margin that makes it unambiguous.

    Pinning the *gap* rather than only the winner: a change that moved every lobe would
    still satisfy "syn-syn is closest to 2.62", and would not satisfy this.
    """
    mol, oxygens = formate
    syn_syn = implied_separation(mol, oxygens, (1, 1))
    others = [implied_separation(mol, oxygens, w) for w in ((0, 0), (0, 1), (1, 0))]

    assert abs(syn_syn - PADDLEWHEEL_CU_CU) < 0.10, (
        f"syn-syn implies Cu...Cu {syn_syn:.3f} A against a literature {PADDLEWHEEL_CU_CU}")
    assert min(others) - syn_syn > 2.0, (
        "the non-bridging combinations must be far away, not merely worse; closest is "
        f"{min(others):.3f} A against syn-syn's {syn_syn:.3f}")


def test_the_span_is_a_property_of_the_ligand_not_of_the_mode():
    """Benzoate is a carboxylate too and does not offer the same bite.

    So a bridging tolerance cannot be one number per mode — it is measured per ligand,
    which is why `bridge_compatible` has to take the frames rather than a table.
    """
    mol = mol_from_smiles("[O-]C(=O)c1ccccc1", embed=True)
    oxygens = [s.atom_idx for s in perceive(mol) if s.donor_type == "carboxylate_O"]
    assert implied_separation(mol, oxygens, (1, 1)) == pytest.approx(2.461, abs=1e-3)


# ── perception records them, and they survive the trip ───────────────────────

def test_perceive_records_every_lobe_and_leaves_lobe_zero_where_it_was(formate):
    """Additive by construction: `frame` is unchanged, so nothing that reads it moves."""
    mol, oxygens = formate
    conf = mol.GetConformer()
    site = next(s for s in perceive(mol) if s.atom_idx == oxygens[0])

    assert len(frame_lobes(site.frame)) == 2
    assert site.frame["axis"] == list(lone_pair_frames(
        mol, oxygens[0], "carboxylate_O", conf)[0].axis_hat)
    assert frame_lobes(site.frame)[0]["axis"] == site.frame["axis"]


def test_a_frame_with_no_lobe_record_reads_as_having_exactly_one():
    """Every site stored before this existed still answers the question."""
    old = {"origin": [0.0, 0.0, 0.0], "axis": [1.0, 0.0, 0.0], "ref": [0.0, 1.0, 0.0],
           "mode": "tilted"}
    assert frame_lobes(old) == [old]
    assert frame_lobes(None) == []


def test_the_lobes_travel_through_the_registry_as_json(formate):
    """`frame_json` is an opaque blob, so this needed no schema change — prove it."""
    mol, oxygens = formate
    site = next(s for s in perceive(mol) if s.atom_idx == oxygens[0])
    assert len(frame_lobes(json.loads(json.dumps(site.frame)))) == 2


def test_a_rigid_move_carries_every_lobe_with_it(formate):
    """A lobe left behind would point at where the metal used to be.

    That is not cosmetic: the second bridge onto an already-moved block is judged against
    these directions, so a stale lobe is a wrong verdict rather than a wrong picture.
    """
    mol, oxygens = formate
    site = next(s for s in perceive(mol) if s.atom_idx == oxygens[0])
    rot = axis_rotation(np.array([0.3, 0.5, 0.8]), 1.1)
    moved = _transform_frame(site.frame, rot, np.zeros(3), np.array([4.0, -2.0, 1.0]))

    before, after = frame_lobes(site.frame), frame_lobes(moved)
    assert len(after) == len(before) == 2
    for old, new in zip(before, after):
        assert np.allclose(rot @ np.array(old["axis"]), np.array(new["axis"]))
    # The relationship between the lobes is what the move must preserve, and it is the
    # thing a per-lobe transform could get right individually and still break.
    def angle(frames):
        return float(np.dot(np.array(frames[0]["axis"]), np.array(frames[1]["axis"])))
    assert angle(after) == pytest.approx(angle(before), abs=1e-12)


# ── a two-point mode does not branch its torsion ─────────────────────────────

@pytest.mark.parametrize("mode", [BindingMode.CHELATE,
                                  BindingMode.BRIDGE_MU2, BindingMode.BRIDGE_MU3])
def test_a_two_point_mode_has_one_torsion_well(mode):
    """The fit sets the roll, so branching it would emit identical siblings.

    `join_chelate` spells the rule out: the donor-donor line onto the vertex-vertex line,
    then the roll fixed by requiring the donor axes to point back at the metal.  A bridge
    is that same fit with the two vertices on different metals, so it inherits the rule.
    """
    assert torsion_wells("carboxylate_O", mode) == (0.0,)


def test_monodentate_still_branches_and_torsion_free_still_never_does():
    """The two-point rule must not leak into the one-point case, or M5's tree collapses."""
    assert torsion_wells("carboxylate_O", BindingMode.MONODENTATE) == (0.0, 180.0)
    assert live_dof("aqua_O") is LiveDOF.TORSION_FREE
    assert torsion_wells("aqua_O", BindingMode.MONODENTATE) == (0.0,)
