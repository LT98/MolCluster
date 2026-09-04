"""Canonicalisation: invariance, exactness, and the WL collisions it exists to fix."""
from __future__ import annotations

import pytest

import fixtures as fx
from mofsbu.graph import (
    canonical_certificate, canonical_index_map, canonical_order, is_isomorphic,
    vf2_isomorphic, wl_hash,
)

PERMUTATION_SEEDS = range(8)


def test_canonical_order_is_permutation_invariant(all_fixtures):
    """The property every site record depends on (D5)."""
    for name, g in all_fixtures.items():
        ref = canonical_certificate(g)
        for seed in PERMUTATION_SEEDS:
            assert canonical_certificate(g.permuted(seed)) == ref, f"{name} seed={seed}"


def test_canonical_index_map_covers_every_atom(all_fixtures):
    for name, g in all_fixtures.items():
        m = canonical_index_map(g)
        assert sorted(m) == g.nodes(), name
        assert sorted(m.values()) == list(range(len(g))), name


@pytest.mark.parametrize("build_a,build_b,label", [
    (fx.zn2_bridged_formates, fx.zn2_chelated_formates, "mu2-bridging vs chelating"),
    (fx.cyclohexane, fx.two_cyclopropanes, "one 6-ring vs two 3-rings"),
])
def test_wl_collides_where_the_certificate_does_not(build_a, build_b, label):
    """Why L1 is the certificate and not the WL hash (D16).

    Both pairs are non-isomorphic yet 1-WL-indistinguishable: every atom has the
    same local environment, so the colours are stable immediately and no number of
    iterations separates them.  The bridging-vs-chelating pair is the one that
    matters chemically — it is a binding mode this project must tell apart.
    """
    a, b = build_a(), build_b()
    assert wl_hash(a) == wl_hash(b), f"{label}: expected a WL collision to demonstrate"
    assert not is_isomorphic(a, b), label
    assert not vf2_isomorphic(a, b), label


def test_certificate_agrees_with_vf2_on_every_fixture_pair(all_fixtures):
    """Independent cross-check of the collision resolver."""
    items = sorted(all_fixtures.items())
    for i, (na, ga) in enumerate(items):
        for nb, gb in items[i:]:
            assert is_isomorphic(ga, gb) == vf2_isomorphic(ga, gb), f"{na} vs {nb}"


def test_symmetric_species_stay_tractable():
    """Automorphism pruning: [Fe(H2O)6]2+ has |Aut| = 46080."""
    import time

    g = fx.fe_hexaaqua()
    t = time.perf_counter()
    canonical_order(g)
    assert (time.perf_counter() - t) < 2.0
