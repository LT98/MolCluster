"""Composing steps into a route: the target at zero, and what the numbers are OF.

`energy/routes.py` prices one edge and `tests/test_reference.py` guards the rules it
applies. This file is about the arithmetic on top — a running sum down a chain — and
about the two ways a running sum can lie: by summing steps that do not actually join,
and by putting two routes on one axis when they are not measured against the same thing.

(`tests/test_pathways.py` is the runner's construction ladder, which is a different
pathway entirely.)
"""
from __future__ import annotations

import pytest

from mofsbu._types import Fidelity, MethodSpec, ReferenceSchemeError
from mofsbu.energy.reference import put_balanced_reaction
from mofsbu.pathways.route import price_path
from mofsbu.registry import Registry

from test_reference import (  # noqa: E402 - the species and the storer are already built
    _store, aqua_ion, hydronium, hydroxo_complex, water,
)

XTB = MethodSpec(code="tblite", code_version="0.4.0", method="GFN2-xTB")
FF = MethodSpec(code="test", code_version="0", method="toy-FF")


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


@pytest.fixture()
def chain(reg):
    """[Ni(H2O)]2+ -> [Ni(H2O)2]2+ -> [Ni(H2O)3]2+, each step adding one water.

    Energies are round numbers so every assertion below is exact rather than a
    tolerance: each step is -5 eV, so the target sits at 0 and its precursors at +5
    and +10.
    """
    one = _store(reg, aqua_ion(1), energy=-40.0)
    two = _store(reg, aqua_ion(2), energy=-60.0)
    three = _store(reg, aqua_ion(3), energy=-80.0)
    wat = _store(reg, water(), energy=-15.0)
    first = put_balanced_reaction(reg, two, reagents=[(one, 1), (wat, 1)])
    second = put_balanced_reaction(reg, three, reagents=[(two, 1), (wat, 1)])
    return {"one": one, "two": two, "three": three, "water": wat,
            "first": first, "second": second}


def test_the_target_sits_at_zero_and_every_node_is_the_running_sum(reg, chain):
    """The anchor is the one node every route to a product shares.

    Read leftwards, y is minus the energy still to be released, so precursors of an
    exothermic assembly sit ABOVE the thing they fall to.
    """
    p = price_path(reg, [chain["three"], chain["two"], chain["one"]],
                   [chain["second"], chain["first"]])
    assert p["can_price"]
    assert [n["y"] for n in p["nodes"]] == [0.0, 5.0, 10.0]
    assert p["total_dE"] == pytest.approx(-10.0)
    # The worst step is the most endothermic one, not the largest in magnitude.
    assert p["worst_step"]["dE"] == pytest.approx(-5.0)


def test_the_pieces_not_yet_used_are_carried_alongside(reg, chain):
    """A y value is the energy of the node PLUS its spectators, or it is not comparable
    to anything.  Each step back leaves one more water unconsumed."""
    p = price_path(reg, [chain["three"], chain["two"], chain["one"]],
                   [chain["second"], chain["first"]])
    free = [n["free_pieces"] for n in p["nodes"]]
    assert free[0] == []
    assert [t["stoich"] for t in free[1]] == [1]
    assert [t["stoich"] for t in free[2]] == [2]
    assert {t["structure_id"] for t in free[2]} == {chain["water"]}


def test_a_chain_whose_steps_do_not_join_is_refused_by_name(reg, chain):
    """Summing two edges that do not meet produces a number for a route that does not
    exist.  It is the one failure a running sum cannot show on its face, so it is
    checked rather than assumed."""
    with pytest.raises(ReferenceSchemeError, match="does not consume structure"):
        price_path(reg, [chain["three"], chain["one"]], [chain["second"]])

    with pytest.raises(ReferenceSchemeError, match="does not produce structure"):
        price_path(reg, [chain["two"], chain["one"]], [chain["second"]])


def test_a_walk_and_its_edges_have_to_be_the_same_length(reg, chain):
    with pytest.raises(ReferenceSchemeError, match="need 1 edges"):
        price_path(reg, [chain["three"], chain["two"]], [])


def test_a_step_with_no_number_leaves_a_gap_and_never_a_zero(reg):
    """D18: absent is not zero.  A missing dE must not read as "this step is free" —
    everything to the left of it has no y at all, and the path says how many."""
    one = _store(reg, aqua_ion(1), energy=-40.0)
    two = _store(reg, aqua_ion(2), energy=-60.0)
    three = _store(reg, aqua_ion(3), energy=-80.0)
    # Only an FF energy for the water in the first step, so that step cannot clear the
    # ML floor while the second one can.
    wat = _store(reg, water(), energy=-15.0, method=FF, fidelity=Fidelity.FF)
    first = put_balanced_reaction(reg, two, reagents=[(one, 1), (wat, 1)])
    second = put_balanced_reaction(reg, three, reagents=[(two, 1), (wat, 1)])

    p = price_path(reg, [three, two, one], [second, first])
    assert not p["can_price"]
    assert p["n_unpriced"] == 2
    assert p["total_dE"] is None
    assert [n["y"] for n in p["nodes"]] == [0.0, None, None]
    # The reason survives onto the path rather than being flattened to a bare refusal.
    assert "has no energy" in p["why_not"] and str(wat) in p["why_not"]


def test_a_path_reports_the_rung_its_steps_were_priced_at(reg, chain):
    """A route is only as good as its weakest step, so the rung it advertises is the
    lowest one across them — never the best, which would describe a step and not a path.

    (A chain that mixes rungs cannot be built here on purpose: `energy_for_terms` pins
    an equation to the THEORY of its first term, so a step is single-method by
    construction and the mixing this guards against can only happen BETWEEN steps.)
    """
    p = price_path(reg, [chain["three"], chain["two"], chain["one"]],
                   [chain["second"], chain["first"]])
    rungs = {s["fidelity"] for s in p["steps"]}
    assert rungs == {int(Fidelity.XTB)}
    assert p["fidelity"] == min(rungs)
    assert p["fidelity_name"] == "XTB"


def test_shedding_a_proton_changes_what_the_numbers_are_measured_against(reg):
    """Two routes may be compared only where the same things have left along the way.

    An assembly step sheds nothing; a deprotonation sheds a hydronium.  After one of
    those the node's y is the energy of a system the other route never has, so the
    basis differs and the page has something to say about it.
    """
    one = _store(reg, aqua_ion(1), energy=-40.0)
    two = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-18.0)
    hydroxo = _store(reg, hydroxo_complex(), energy=-55.0)

    assembly = put_balanced_reaction(reg, two, reagents=[(one, 1), (wat, 1)])
    deprot = put_balanced_reaction(reg, hydroxo, reagents=[(two, 1)],
                                   solvents=[(wat, 1)], leaving=[(h3o, 1)],
                                   kind="deprotonation")

    plain = price_path(reg, [two, one], [assembly])
    shedding = price_path(reg, [hydroxo, two, one], [deprot, assembly])

    assert plain["basis"] == "—", "an assembly route sheds nothing"
    assert shedding["basis"] == f"{h3o}x1"
    assert shedding["basis"] != plain["basis"]
    # The node BEFORE the deprotonation is still on the plain basis: routes diverge at
    # the step that sheds, not at the route level.
    assert [n["basis"] for n in shedding["nodes"]] == ["—", f"{h3o}x1", f"{h3o}x1"]
    # Water carried the proton; it is not a piece that was built in.
    assert shedding["steps"][0]["added"] == []
    assert [t["structure_id"] for t in shedding["steps"][0]["solvent"]] == [wat]
    assert [t["structure_id"] for t in shedding["steps"][0]["shed"]] == [h3o]


def test_a_derived_leg_is_marked_inferred_and_carries_no_reaction(reg):
    """Option 2 on C15: a structure the runner built whole has no recorded precursor,
    so the only answer to what it is made of is derivation.  It is priced and walkable
    — and never labelled as something that was recorded."""
    one = _store(reg, aqua_ion(1), energy=-40.0)
    two = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)

    p = price_path(reg, [two, one], [f"d{wat}"])
    assert p["can_price"]
    assert p["steps"][0]["origin"] == "inferred"
    assert p["steps"][0]["reaction_id"] is None
    assert p["nodes"][1]["y"] == pytest.approx(5.0)
    assert [t["structure_id"] for t in p["steps"][0]["added"]] == [wat]


def test_a_derived_leg_that_never_split_that_way_is_refused(reg):
    one = _store(reg, aqua_ion(1), energy=-40.0)
    two = _store(reg, aqua_ion(2), energy=-60.0)
    _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-18.0)

    with pytest.raises(ReferenceSchemeError, match="does not split into"):
        price_path(reg, [two, one], [f"d{h3o}"])
