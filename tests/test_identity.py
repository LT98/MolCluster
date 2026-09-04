"""The M2 exit gate: what L0/L1 must and must not tell apart."""
from __future__ import annotations

import json

import pytest

import fixtures as fx
from mofsbu.identity import block_id, identity, l0_composition, l1_graph_hash, l2_isomer_tag

PERMUTATION_SEEDS = range(8)


# -- gate 1: build-order invariance --------------------------------------------

def test_l1_is_build_order_invariant(all_fixtures):
    for name, g in all_fixtures.items():
        ref = l1_graph_hash(g)
        for seed in PERMUTATION_SEEDS:
            assert l1_graph_hash(g.permuted(seed)) == ref, f"{name} seed={seed}"


def test_l0_is_build_order_invariant(all_fixtures):
    for name, g in all_fixtures.items():
        ref = l0_composition(g)
        for seed in PERMUTATION_SEEDS:
            assert l0_composition(g.permuted(seed)) == ref, f"{name} seed={seed}"


# -- gate 2: connectivity discrimination ---------------------------------------

def test_bridging_and_chelating_are_different_structures():
    """Same atoms, different binding mode -> different L1."""
    a, b = fx.zn2_bridged_formates(), fx.zn2_chelated_formates()
    assert l0_composition(a) == l0_composition(b)     # identical composition
    assert l1_graph_hash(a) != l1_graph_hash(b)


def test_protomers_are_different_structures():
    """'Activate = deprotonate' is the core primitive, so protomers must be distinct."""
    assert l1_graph_hash(fx.formate()) != l1_graph_hash(fx.formic_acid())
    assert l1_graph_hash(fx.btc(True)) != l1_graph_hash(fx.btc(False))


# -- gate 3: L1 must NOT over-discriminate -------------------------------------

def test_l1_does_not_separate_cis_from_trans():
    """cis/trans-Pt(NH3)2Cl2 share a connectivity.

    This test must keep passing.  Configurational isomers are separated at L2 by
    downstream-assembly relevance (D10), NOT by making L1 finer — a finer L1 would
    also split things that are genuinely one structure.
    """
    a, b = fx.pt_ammine_dichloride("a"), fx.pt_ammine_dichloride("b")
    assert l1_graph_hash(a) == l1_graph_hash(b)
    assert l2_isomer_tag(a) == l2_isomer_tag(b) == ""      # stub until M5


# -- gate 4: per-centre labels (D12) -------------------------------------------

def test_mixed_valence_differs_at_l0():
    homo, mixed = fx.fe3_mu3_oxo((3, 3, 3)), fx.fe3_mu3_oxo((2, 3, 3))
    assert l0_composition(homo) != l0_composition(mixed)
    assert "Fe(+2,hs)" in l0_composition(mixed)
    assert l1_graph_hash(homo) != l1_graph_hash(mixed)     # labels are hashed too


def test_spin_state_is_part_of_identity():
    hs, ls = fx.fe_hexaaqua("hs"), fx.fe_hexaaqua("ls")
    assert l0_composition(hs) != l0_composition(ls)
    assert l1_graph_hash(hs) != l1_graph_hash(ls)


def test_paddlewheel_is_one_node_with_two_equivalent_centres():
    pw = fx.cu_paddlewheel()
    assert l0_composition(pw) == "Cu2_C4H4O8_q0_s1|Cu(+2,hs),Cu(+2,hs)"
    assert len(pw.metals()) == 2


# -- gate 5: two-machine parity ------------------------------------------------

def test_golden_hashes(all_fixtures, golden_path):
    """Laptop and workstation must agree, and hashes must not drift silently.

    Regenerate deliberately with `python scripts/regen_golden.py` and commit the
    diff alongside the ALGO_VERSIONS bump that justifies it.
    """
    golden = json.loads(golden_path.read_text())
    current = {name: identity(g)["l1"] for name, g in all_fixtures.items()}
    assert current == golden["l1"], "L1 hashes drifted — bump ALGO_VERSIONS if intended"
    assert {n: identity(g)["l0"] for n, g in all_fixtures.items()} == golden["l0"]


def test_block_id_shape():
    g = fx.formate()
    bid = block_id(l0_composition(g), l1_graph_hash(g))
    assert bid.count("/") == 3 and bid.startswith("CHO2_q-1_s1/")
