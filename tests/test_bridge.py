"""The lone-pair lobe as a branch — M6/S1, mechanism A's enabling step.

A bridge is the first operation that cares *which* lone pair of an sp2 donor binds.
`sites.frames.site_frame` has always known there are two and says so in its own comment;
`sites.model.perceive` kept the first and discarded the rest, so every stored site pointed
its metal at the syn lobe and no bridging geometry was reachable from the assembly path.

The numbers here are the calibration, not decoration: the four lobe combinations of a
formate's two oxygens are the three textbook carboxylate bridging modes, and they are
separated by more than a bond length.  A test that only asserted "there are two lobes"
would pass on two lobes that pointed the same way.

Two layers, measured separately and required to agree.  What the FRAMES imply is the
first half; what two ordinary `join` calls actually BUILD is the second, and they land on
the same four numbers.  That agreement is D20 in its smallest form — a bridged dimer is
not placed by a solver, it is what a sequence of one-contact joins leaves behind, and its
M...M is an output to validate rather than an input to impose.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from mofsbu.assembly.construct import construct, ligand_block, metal_block
from mofsbu.assembly.join import (
    MAX_BRIDGE_SPAN_MISMATCH_A, IncompatibleJoin, _transform_frame, bridge_compatible,
    compatible, join, join_bridge)
from mofsbu.geometry._linalg import axis_rotation
from mofsbu.geometry.distances import metal_donor_distance
from mofsbu.geometry.embed import embed_molecule
from mofsbu.geometry.qc import qc
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


# ── the same three modes, reached through the assembly path ──────────────────
#
# The section above measures what the FRAMES imply.  This one measures what two ordinary
# `join` calls actually build, and the two agreeing is the point: a bridge is not a new
# operation, it is a sequence of one-contact joins whose M...M is an output (D20).


def bridged_dimer(lobe_a: int, lobe_b: int, smiles: str = "[O-]C=O", *,
                  cn: int = 4, geom: str = "square_planar"):
    """Formate onto one Cu, then its free oxygen onto a second.  Two ordinary joins."""
    ligand = ligand_block(embed_molecule(mol_from_smiles(smiles), seed=7),
                          charge=-1, name="bridge")
    donors = sorted((s for s in ligand.open_donors()
                     if s.donor_type == "carboxylate_O"), key=lambda s: s.atom_idx)
    first_metal = metal_block("Cu", 2, cn=cn, geometry=geom)
    first = join(ligand, first_metal, donors[0], first_metal.open_vacancies()[0],
                 lone_pair=lobe_a)

    free = next(s for s in first.block.open_donors()
                if s.atom_idx == first.atom_map[donors[1].atom_idx])
    second_metal = metal_block("Cu", 2, cn=cn, geometry=geom)
    return join(first.block, second_metal, free, second_metal.open_vacancies()[0],
                lone_pair=lobe_b)


def metal_separation(block) -> float:
    g, coords = block.graph, block.geometry
    rows = {node: k for k, node in enumerate(g.nodes())}
    m1, m2 = g.metals()
    return float(np.linalg.norm(coords[rows[m1]] - coords[rows[m2]]))


@pytest.mark.parametrize("lobes, expected, mode", [
    ((1, 1), 2.673, "syn-syn"),
    ((0, 1), 5.148, "syn-anti"),
    ((1, 0), 5.148, "anti-syn"),
    ((0, 0), 5.516, "anti-anti"),
])
def test_the_assembly_path_reaches_every_bridging_mode(lobes, expected, mode):
    """Separately enumerable, and each one lands where its frames said it would.

    Nothing here imposes a Cu···Cu. The two joins each satisfy one contact by a rigid
    move, and the separation is read off the product afterwards — which is what makes it
    a validation rather than a constraint (D20).
    """
    result = bridged_dimer(*lobes)
    assert len(result.block.graph.metals()) == 2
    assert metal_separation(result.block) == pytest.approx(expected, abs=1e-3), mode


def test_the_syn_syn_route_lands_on_the_paddlewheel_separation():
    """The one measurement M6 exists to reach, against the literature and with margin."""
    syn_syn = metal_separation(bridged_dimer(1, 1).block)
    others = [metal_separation(bridged_dimer(*w).block)
              for w in ((0, 0), (0, 1), (1, 0))]
    assert abs(syn_syn - PADDLEWHEEL_CU_CU) < 0.10, (
        f"syn-syn built Cu...Cu {syn_syn:.3f} A against a literature {PADDLEWHEEL_CU_CU}")
    assert min(others) - syn_syn > 2.0


def test_the_lobe_is_recorded_so_a_replay_reproduces_it_rather_than_the_default():
    """A choice that did not travel in the vector would silently replay as lobe 0 —
    which is a different structure wearing the original's provenance, and 2.8 A away."""
    result = bridged_dimer(1, 1)
    assert result.choice_vector["donor"]["lone_pair"] == 1

    ligand = ligand_block(embed_molecule(mol_from_smiles("[O-]C=O"), seed=7),
                          charge=-1, name="bridge")
    metal = metal_block("Cu", 2, cn=4, geometry="square_planar")
    donor = sorted((s for s in ligand.open_donors()
                    if s.donor_type == "carboxylate_O"), key=lambda s: s.atom_idx)[0]
    step = dict(join(ligand, metal, donor, metal.open_vacancies()[0],
                     lone_pair=1).choice_vector, partner=0)
    again = construct(ligand, [metal], [step])
    assert again.steps[0]["donor"]["lone_pair"] == 1
    assert np.array_equal(
        again.block.geometry,
        join(ligand, metal, donor, metal.open_vacancies()[0], lone_pair=1).block.geometry)


def test_a_donor_with_one_lobe_ignores_the_index_rather_than_failing():
    """An aqua has one direction, so `lone_pair=1` is not an error — it is the same bond.

    The index wraps over what the donor ACTUALLY offers, so a caller enumerating lobes
    uniformly across a mixed donor set does not have to special-case the determined ones.
    """
    water = ligand_block(embed_molecule(mol_from_smiles("O"), seed=7), charge=0,
                         name="aqua")
    donor = water.open_donors()[0]
    metal = metal_block("Zn", 2, cn=4, geometry="tetrahedral")
    a = join(water, metal, donor, metal.open_vacancies()[0], lone_pair=0)
    b = join(water, metal, donor, metal.open_vacancies()[0], lone_pair=1)
    assert np.array_equal(a.block.geometry, b.block.geometry)
    assert b.choice_vector["donor"]["lone_pair"] == 0
    assert a.compatibility.lone_pairs == 1


def test_the_verdict_says_a_choice_is_being_made():
    """Two real lobes and no default: the verdict reports the branch rather than hiding
    it behind whichever one `perceive` happened to store first."""
    ligand = ligand_block(embed_molecule(mol_from_smiles("[O-]C=O"), seed=7),
                          charge=-1, name="formate")
    metal = metal_block("Cu", 2, cn=4, geometry="square_planar")
    donor = next(s for s in ligand.open_donors() if s.donor_type == "carboxylate_O")
    verdict = compatible(donor, metal.open_vacancies()[0], partner="Cu", donor_element="O")
    assert verdict.lone_pairs == 2
    assert "lone_pair" in verdict.reason


# ── join_bridge: one ligand, two metals, one move ────────────────────────────


def cross_metal_vertices(dimer_block):
    """Every (vertex on metal A, vertex on metal B) pair a bridge could take."""
    m_a, m_b = dimer_block.graph.metals()
    per = {}
    for s in dimer_block.open_vacancies():
        per.setdefault(s.atom_idx, []).append(s)
    return [(va, vb) for va in per[m_a] for vb in per[m_b]]


def second_bridge_onto(dimer_block, slots=None):
    """A second formate, and the vertex pair a caller would actually pick for it.

    Ranked by `bridge_compatible`, because ranking on one criterion and placing under
    another is how a step chooses the pair it then cannot build — `chelate_reach` says so
    one layer down and it is no less true here.
    """
    lig = ligand_block(embed_molecule(mol_from_smiles("[O-]C=O"), seed=7),
                       charge=-1, name="bridge2")
    donors = sorted((s for s in lig.open_donors() if s.donor_type == "carboxylate_O"),
                    key=lambda s: s.atom_idx)
    pairs = cross_metal_vertices(dimer_block)
    if slots is not None:
        pair = next((a, b) for a, b in pairs if (a.slot, b.slot) == tuple(slots))
    else:
        pair = min(pairs, key=lambda vs: bridge_compatible(
            donors, list(vs), partner="Cu", donor_elements=["O", "O"]).strain)
    return lig, donors, list(pair)


def test_a_bridge_binds_each_donor_to_its_OWN_metal():
    """The line `join_chelate` cannot draw: its two bonds both go to `vacancies[0]`.

    µ-ness is then derived from the two dative edges rather than labelled (D14), so the
    product reporting two mu2 fragments is the graph agreeing with the operation.
    """
    dimer = bridged_dimer(1, 1).block
    lig, donors, pair = second_bridge_onto(dimer)
    result = join_bridge(lig, dimer, donors, pair, lone_pairs=(0, 0))

    g = result.block.graph
    assert len(g.metals()) == 2
    bridges = [f for f in g.ligand_fragments()
               if g.fragment_bridge_class(f).value == "mu2"]
    assert len(bridges) == 2, "the new ligand did not come out bridging"
    assert result.choice_vector["mode"] == "mu2"
    assert [v["metal"] for v in result.choice_vector["vacancies"]] == ["Cu", "Cu"]


def test_two_vertices_on_ONE_metal_are_refused_and_named_as_a_chelate():
    dimer = bridged_dimer(1, 1).block
    m_a = dimer.graph.metals()[0]
    same = [s for s in dimer.open_vacancies() if s.atom_idx == m_a][:2]
    lig, donors, _ = second_bridge_onto(dimer)
    verdict = bridge_compatible(donors, same, partner="Cu", donor_elements=["O", "O"])
    assert not verdict.feasible
    assert "chelate" in verdict.reason and "join_chelate" in verdict.reason


def test_one_atom_on_two_metals_is_refused_as_the_other_mechanism():
    """Mechanism B is not this operation, and saying so by name is the point of §4."""
    dimer = bridged_dimer(1, 1).block
    lig, donors, pair = second_bridge_onto(dimer)
    verdict = bridge_compatible([donors[0], donors[0]], pair, partner="Cu",
                                donor_elements=["O", "O"])
    assert not verdict.feasible
    assert "SINGLE-atom bridge" in verdict.reason


def test_a_donor_that_does_not_bridge_is_refused_by_name():
    """`pyridyl_N` declares `mono` only, and the refusal says which modes it does offer."""
    bipy = ligand_block(embed_molecule(mol_from_smiles("c1ccc(-c2ccccn2)nc1"), seed=7),
                        charge=0, name="bipy")
    nitrogens = sorted((s for s in bipy.open_donors() if s.donor_type == "pyridyl_N"),
                       key=lambda s: s.atom_idx)[:2]
    assert len(nitrogens) == 2
    dimer = bridged_dimer(1, 1).block
    _lig, _donors, pair = second_bridge_onto(dimer)
    verdict = bridge_compatible(nitrogens, pair, partner="Cu",
                                donor_elements=["N", "N"])
    assert not verdict.feasible
    assert "does not bridge" in verdict.reason and "mono" in verdict.reason


def test_the_verdict_is_necessary_and_not_sufficient_and_qc_is_the_arbiter():
    """Measured, and it is the reason the tolerance is a bound rather than a derivation.

    On a CN-6 dimer the *smallest* span mismatch in the whole candidate set — 0.183 Å,
    better than the square-pyramidal case that builds cleanly — places the ligand's carbon
    1.4 Å from the far metal. The span test cannot see that, because where a mismatch goes
    depends on the vertex axes and not on its size. So a caller ranks on this verdict and
    then *checks*, which is what `join_chelate`'s own docstring says one layer down.
    """
    dimer = bridged_dimer(1, 1, cn=6, geom="octahedral").block
    lig, donors, pair = second_bridge_onto(dimer)      # the BEST-ranked pair
    verdict = bridge_compatible(donors, pair, partner="Cu", donor_elements=["O", "O"])
    assert verdict.feasible, verdict.reason
    assert "Necessary and not sufficient" in verdict.reason

    result = join_bridge(lig, dimer, donors, pair)
    g, coords = result.block.graph, result.block.geometry
    rows = {node: k for k, node in enumerate(g.nodes())}
    symbols = [g.label(i).element for i in g.nodes()]
    bonded = {(rows[i], rows[j]) for i, j, _t in g.edges()}
    report = qc(symbols, coords, bonded)
    assert not report.ok, "qc must catch what the span test cannot see"


def test_a_span_nothing_can_absorb_is_refused_with_its_number():
    """A ligand far too short for the vertices it is offered, refused and quantified."""
    dimer = bridged_dimer(0, 0).block          # anti-anti: the metals are 5.5 Å apart
    lig, donors, _ = second_bridge_onto(dimer)

    pairs = cross_metal_vertices(dimer)
    verdicts = [(bridge_compatible(donors, list(vs), partner="Cu",
                                   donor_elements=["O", "O"]), vs) for vs in pairs]
    worst, vs = max(verdicts, key=lambda pv: pv[0].strain)
    assert not worst.feasible
    assert "cannot reach" in worst.reason
    # The number is in the verdict, not only in the exception's type.
    assert f"{worst.strain * MAX_BRIDGE_SPAN_MISMATCH_A:.3f} A of span mismatch" in worst.reason

    with pytest.raises(IncompatibleJoin) as exc:
        join_bridge(lig, dimer, donors, list(vs))
    assert "cannot reach" in str(exc.value)


def test_the_lobes_are_recorded_for_both_donors():
    """A bridge's discrete choice is which lobe each end used — there is no torsion well
    to record, because the second contact determines the roll."""
    dimer = bridged_dimer(1, 1).block
    lig, donors, pair = second_bridge_onto(dimer)
    cv = join_bridge(lig, dimer, donors, pair, lone_pairs=(1, 0)).choice_vector
    assert [d["lone_pair"] for d in cv["donors"]] == [1, 0]
    assert "torsion_well" not in cv


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
