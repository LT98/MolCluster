"""Typed-graph invariants (PLAN §0 ground rules 1-6)."""
from __future__ import annotations

import pytest

import fixtures as fx
from mofsbu.graph import BridgeClass, EdgeType, TypedGraph
from mofsbu._types import AmbiguousSpecError, GraphValidationError


def test_dative_must_run_donor_to_metal():
    g = TypedGraph(charge=0, multiplicity=1)
    o, c = g.add_atom("O"), g.add_atom("C")
    with pytest.raises(GraphValidationError):
        g.add_bond(o, c, EdgeType.DATIVE)


def test_two_metals_may_not_be_covalent():
    g = TypedGraph(charge=0, multiplicity=1)
    a = g.add_atom("Cu", oxidation_state=2, spin_class="hs")
    b = g.add_atom("Cu", oxidation_state=2, spin_class="hs")
    with pytest.raises(GraphValidationError):
        g.add_bond(a, b, EdgeType.COVALENT)
    g.add_bond(a, b, EdgeType.METAL_METAL)


def test_per_centre_labels_are_metal_only():
    g = TypedGraph(charge=0, multiplicity=1)
    with pytest.raises(GraphValidationError):
        g.add_atom("O", oxidation_state=2)


def test_charge_and_multiplicity_are_never_guessed():
    """Ground rule 5: ambiguity raises, it does not default."""
    from mofsbu.identity import l0_composition

    no_charge = TypedGraph(charge=None, multiplicity=1)
    no_charge.add_atom("O")
    with pytest.raises(AmbiguousSpecError):
        l0_composition(no_charge)

    no_spin = TypedGraph(charge=0, multiplicity=None)
    no_spin.add_atom("O")
    with pytest.raises(AmbiguousSpecError):
        l0_composition(no_spin)

    assert l0_composition(fx.water()) == "H2O_q0_s1"   # both given -> fine


def test_bridge_class_is_derived_not_stored():
    """mu2 is a consequence of the dative edges, never an input (D14)."""
    bridged = fx.zn2_bridged_formates()
    chelated = fx.zn2_chelated_formates()
    bridging_os = [i for i in bridged.nodes()
                   if bridged.bridge_class(i) is BridgeClass.TERMINAL]
    assert len(bridging_os) == 4                     # each O serves exactly one Zn
    pw = fx.cu_paddlewheel()
    # every OXYGEN is terminal — each donates to one Cu ...
    assert all(pw.bridge_class(i) is not BridgeClass.MU2 for i in pw.nodes())
    # ... but every carboxylate LIGAND spans two, which is what mu2 means.
    frags = pw.ligand_fragments()
    assert len(frags) == 4
    assert all(pw.fragment_bridge_class(f) is BridgeClass.MU2 for f in frags)
    assert pw.max_bridge_class() is BridgeClass.MU2

    trimer = fx.fe3_mu3_oxo()
    mu3 = [i for i in trimer.nodes() if trimer.bridge_class(i) is BridgeClass.MU3]
    assert len(mu3) == 1                             # the central oxo, and only it
    assert trimer.max_bridge_class() is BridgeClass.MU3
    assert fx.water().max_bridge_class() is BridgeClass.NONE


def test_json_round_trip_is_byte_identical(all_fixtures):
    for name, g in all_fixtures.items():
        assert TypedGraph.from_json(g.to_json()).to_json() == g.to_json(), name


def test_bond_order_is_not_in_the_hashed_view():
    """Resonance forms must not be different structures."""
    from mofsbu.identity import l1_graph_hash

    a = fx.formate()
    b = fx.formate()
    c_idx = next(i for i in b.nodes() if b.label(i).element == "C")
    o_idx = b.neighbors(c_idx)[0]
    b.add_bond(c_idx, o_idx, EdgeType.COVALENT, order=2.0)   # same bond, order bumped
    assert l1_graph_hash(a) == l1_graph_hash(b)
