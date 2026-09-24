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
from mofsbu.energy.reference import CHARGE_SEPARATION, COORDINATION_CHANGE, put_balanced_reaction
from mofsbu.graph import EdgeType, TypedGraph
from mofsbu.pathways.route import price_path
from mofsbu.registry import Registry, outgoing_routes

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


# ── walking both ways: an exchange composed at a shared parent ───────────────
#
# The shape of the Ni/tHQ/Cl case on data/mvp_ni_thq_cl.db (NiCl2(H2O)2 back to
# Ni(H2O)2 2+, forward to Ni(tHQ-)2(H2O)2), rebuilt with hydroxide and chloride and round
# energies so every number below is exact.

def chloride() -> TypedGraph:
    g = TypedGraph(charge=-1, multiplicity=1, name="Cl-")
    g.add_atom("Cl", formal_charge=-1)
    return g


def hydroxide() -> TypedGraph:
    g = TypedGraph(charge=-1, multiplicity=1, name="OH-")
    o = g.add_atom("O", formal_charge=-1)
    g.add_bond(o, g.add_atom("H"), EdgeType.COVALENT)
    return g


def diaqua(*anions: str) -> TypedGraph:
    """[Ni(H2O)2(X)n] with each X an OH- or Cl- bound through its anionic atom."""
    g = TypedGraph(charge=2 - len(anions), multiplicity=3,
                   name="Ni(H2O)2" + "".join(f"({x})" for x in anions))
    m = g.add_atom("Ni", oxidation_state=2, spin_class="hs")
    for _ in range(2):
        o = g.add_atom("O")
        g.add_bond(o, m, EdgeType.DATIVE)
        for _ in range(2):
            g.add_bond(o, g.add_atom("H"), EdgeType.COVALENT)
    for x in anions:
        d = g.add_atom("O" if x == "OH" else "Cl", formal_charge=-1)
        g.add_bond(d, m, EdgeType.DATIVE)
        if x == "OH":
            g.add_bond(d, g.add_atom("H"), EdgeType.COVALENT)
    return g


@pytest.fixture()
def exchange(reg):
    """P = [Ni(H2O)2]2+, the shared parent.  OH- builds A then B off it; Cl- builds C
    then D.  Every edge is an assembly, so every edge changes the coordination count.

        rA  P + OH -> A   -5        rC  P + Cl -> C   -6
        rB  A + OH -> B   -7        rD  C + Cl -> D   -4
        rE  A + Cl -> E   -2        rE2 C + OH -> E  -1   (E = [Ni(H2O)2(OH)Cl])
    """
    s = {"P": _store(reg, aqua_ion(2), energy=-60.0),
         "OH": _store(reg, hydroxide(), energy=-10.0),
         "Cl": _store(reg, chloride(), energy=-20.0),
         "A": _store(reg, diaqua("OH"), energy=-75.0),
         "B": _store(reg, diaqua("OH", "OH"), energy=-92.0),
         "C": _store(reg, diaqua("Cl"), energy=-86.0),
         "D": _store(reg, diaqua("Cl", "Cl"), energy=-110.0),
         "E": _store(reg, diaqua("OH", "Cl"), energy=-97.0)}
    for name, product, parts in (("rA", "A", ("P", "OH")), ("rB", "B", ("A", "OH")),
                                 ("rC", "C", ("P", "Cl")), ("rD", "D", ("C", "Cl")),
                                 ("rE", "E", ("A", "Cl")), ("rE2", "E", ("C", "OH"))):
        s[name] = put_balanced_reaction(reg, s[product], kind="assembly",
                                        reagents=[(s[p], 1) for p in parts])
    return s


def _walk(s):
    """D <- C <- P -> A -> B: back to the shared parent, then forward."""
    return ([s["D"], s["C"], s["P"], s["A"], s["B"]],
            [s["rD"], s["rC"], f"c{s['rA']}", f"c{s['rB']}"])


def test_a_consumed_by_leg_enters_the_route_at_minus_its_dE(reg, exchange):
    """Walked with its arrow, an edge runs in reverse along the route: the step's dE is
    the edge's with the sign flipped, and the edge's own number is kept beside it."""
    p = price_path(reg, *_walk(exchange))
    assert p["can_price"]
    assert [s["direction"] for s in p["steps"]] == ["made_from", "made_from",
                                                     "consumed_by", "consumed_by"]
    assert [s["reaction_dE"] for s in p["steps"]] == pytest.approx([-4, -6, -5, -7])
    assert [s["dE"] for s in p["steps"]] == pytest.approx([-4, -6, 5, 7])
    assert [n["y"] for n in p["nodes"]] == pytest.approx([0, 4, 10, 5, -2])
    assert p["total_dE"] == pytest.approx(2.0)
    assert p["steps"][2]["via"] == f"c{exchange['rA']}"


def test_the_signed_sum_is_the_net_equation_priced_whole(reg, exchange):
    """B + 2 Cl- -> D + 2 OH-.  Summing the legs with backward ones negated and pricing
    that one equation directly are the same number, or one of them is wrong."""
    p = price_path(reg, *_walk(exchange))
    eq = p["net_equation"]
    s = exchange
    assert {(t["structure_id"], t["side"], t["stoich"]) for t in eq["terms"]} == {
        (s["D"], "product", 1), (s["OH"], "product", 2),
        (s["B"], "reagent", 1), (s["Cl"], "reagent", 2)}
    assert eq["terms"][0]["structure_id"] == s["D"], "the target pins the theory"
    assert eq["balanced"] is True
    assert eq["dE"] == pytest.approx(p["total_dE"])
    assert eq["agrees_with_steps"] is True


def test_a_route_is_judged_on_its_net_equation_not_its_legs(reg, exchange):
    """Every leg forms or breaks a metal-donor bond, so every leg is non-isodesmic.  The
    exchange as a whole keeps the count, so the route is — and the legs keep saying what
    they are, because a per-step caveat is still true of the step."""
    p = price_path(reg, *_walk(exchange))
    assert all(s["isodesmic"] is False for s in p["steps"])
    assert all(COORDINATION_CHANGE in s["caveats"] for s in p["steps"])
    assert p["isodesmic"] is True
    assert p["caveats"] == []
    assert p["net_equation"]["dative_delta"] == 0


def test_a_route_that_only_builds_up_is_still_not_isodesmic(reg, exchange):
    """The net equation is not a loophole: P + 2 Cl -> D forms two bonds from nothing."""
    s = exchange
    p = price_path(reg, [s["D"], s["C"], s["P"]], [s["rD"], s["rC"]])
    assert p["isodesmic"] is False
    # Three ions become one neutral complex, so the gas-phase charge rule (C18) applies too.
    assert sorted(p["caveats"]) == sorted([COORDINATION_CHANGE, CHARGE_SEPARATION])
    assert p["net_equation"]["dative_delta"] == 2


def test_what_a_backward_leg_added_is_released_and_nets_against_what_was_taken(
        reg, exchange):
    """Run in reverse, an edge gives back the piece it built in.  The tally is one
    signed multiset, so a chloride taken up on the way out and given back on the way in
    cancels rather than being listed on both sides."""
    s = exchange
    p = price_path(reg, *_walk(s))
    last = p["nodes"][-1]
    assert [(t["structure_id"], t["stoich"]) for t in last["free_pieces"]] == [(s["Cl"], 2)]
    assert [(t["structure_id"], t["stoich"]) for t in last["shed"]] == [(s["OH"], 2)]
    assert last["basis"] == f"{s['OH']}x2"
    assert [t["structure_id"] for t in p["steps"][2]["shed"]] == [s["OH"]]
    assert p["steps"][2]["added"] == []

    # D <- C <- P -> A -> E: two chlorides taken up, one given back by rE run in reverse.
    q = price_path(reg, [s["D"], s["C"], s["P"], s["A"], s["E"]],
                   [s["rD"], s["rC"], f"c{s['rA']}", f"c{s['rE']}"])
    end = q["nodes"][-1]
    assert [(t["structure_id"], t["stoich"]) for t in end["free_pieces"]] == [(s["Cl"], 1)]
    assert [(t["structure_id"], t["stoich"]) for t in end["shed"]] == [(s["OH"], 1)]
    assert {(t["structure_id"], t["side"], t["stoich"])
            for t in q["net_equation"]["terms"]} == {
        (s["D"], "product", 1), (s["OH"], "product", 1),
        (s["E"], "reagent", 1), (s["Cl"], "reagent", 1)}
    assert q["total_dE"] == pytest.approx(-3.0)
    assert q["net_equation"]["dE"] == pytest.approx(-3.0)


def test_a_pivot_is_flagged_and_no_step_is_singled_out_across_it(reg, exchange):
    """The shared parent is where the walk turned.  Its y is the energy of the bare
    parent plus every ligand, which is bookkeeping for the exchange and not a height the
    route climbs — so it is marked, and no leg next to it is reported as the worst."""
    p = price_path(reg, *_walk(exchange))
    assert p["pivots"] == [2]
    assert [n["pivot"] for n in p["nodes"]] == [False, False, True, False, False]
    assert p["nodes"][2]["pivot_kind"] == "parent"
    assert p["worst_step"] is None
    assert "pivot" in p["worst_step_why"]

    # Turning at a shared PRODUCT — forward to E, then back down its other edge — is a
    # pivot too, of the other kind: A -> E <- C, net C + OH -> A + Cl.
    s = exchange
    q = price_path(reg, [s["A"], s["E"], s["C"]], [f"c{s['rE']}", s["rE2"]])
    assert q["pivots"] == [1]
    assert q["nodes"][1]["pivot_kind"] == "product"
    assert q["total_dE"] == pytest.approx(1.0)
    assert q["net_equation"]["dE"] == pytest.approx(1.0)
    assert q["isodesmic"] is True


def test_a_composed_route_is_neither_recorded_nor_inferred(reg, exchange):
    """Its net equation is not a row, and it was not derived from composition either:
    it was composed from rows, which are named as the witness."""
    s = exchange
    p = price_path(reg, *_walk(s))
    assert p["origin"] == "composed"
    assert p["witness"] == [
        {"reaction_id": s["rD"], "direction": "made_from"},
        {"reaction_id": s["rC"], "direction": "made_from"},
        {"reaction_id": s["rA"], "direction": "consumed_by"},
        {"reaction_id": s["rB"], "direction": "consumed_by"}]
    assert all(st["origin"] == "recorded" for st in p["steps"])

    one = price_path(reg, [s["C"], s["P"]], [s["rC"]])
    assert one["origin"] == "recorded", "one edge read as written IS the row"
    backwards = price_path(reg, [s["P"], s["C"]], [f"c{s['rC']}"])
    assert backwards["origin"] == "composed", "a row read backwards is not the row"
    assert backwards["total_dE"] == pytest.approx(6.0)


def test_a_consumed_by_leg_has_to_join_the_nodes_it_is_given_for(reg, exchange):
    s = exchange
    with pytest.raises(ReferenceSchemeError, match="does not produce structure"):
        price_path(reg, [s["P"], s["D"]], [f"c{s['rC']}"])
    with pytest.raises(ReferenceSchemeError, match="does not consume structure"):
        price_path(reg, [s["A"], s["C"]], [f"c{s['rC']}"])


def test_walking_an_edge_straight_back_is_refused(reg, exchange):
    """C <- P -> C along the same edge is not a route; it is the trail's back button."""
    s = exchange
    with pytest.raises(ReferenceSchemeError, match="straight back"):
        price_path(reg, [s["C"], s["P"], s["C"]], [s["rC"], f"c{s['rC']}"])


def test_outgoing_routes_are_the_edges_that_take_a_structure_as_a_reagent(reg, exchange):
    s = exchange
    got = outgoing_routes(reg, s["P"])
    assert sorted(r["id"] for r in got) == sorted([s["rA"], s["rC"]])
    assert {r["product_structure_id"] for r in got} == {s["A"], s["C"]}
    assert sorted(r["id"] for r in outgoing_routes(reg, s["Cl"])) == sorted(
        [s["rC"], s["rD"], s["rE"]])


def test_the_proton_couple_is_not_a_hub(reg):
    """Water in and H3O+ out of a deprotonation carry a proton; they are not built on.
    So neither lists the edge as consuming it, and a consumed-by leg through the carrier
    is refused by name — while the protonated parent still walks forward."""
    two = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-18.0)
    hydroxo = _store(reg, hydroxo_complex(), energy=-55.0)
    deprot = put_balanced_reaction(reg, hydroxo, reagents=[(two, 1)],
                                   solvents=[(wat, 1)], leaving=[(h3o, 1)],
                                   kind="deprotonation")
    assert outgoing_routes(reg, wat) == []
    assert outgoing_routes(reg, h3o) == []
    assert [r["id"] for r in outgoing_routes(reg, two)] == [deprot]

    with pytest.raises(ReferenceSchemeError, match="proton carrier"):
        price_path(reg, [wat, hydroxo], [f"c{deprot}"])

    # Protonation, walked: the proton is taken up and its carrier released.
    p = price_path(reg, [two, hydroxo], [f"c{deprot}"])
    step = p["steps"][0]
    assert [t["structure_id"] for t in step["added"]] == [h3o]
    assert [t["structure_id"] for t in step["shed"]] == [wat]
    assert p["nodes"][1]["basis"] == f"{wat}x1"
    assert step["dE"] == pytest.approx(-step["reaction_dE"])


# ── where a released proton goes ─────────────────────────────────────────────

def _deprotonation_with_acetate(reg):
    from mofsbu.energy.protons import link_protomers

    from test_reference import _from_smiles

    s = {"acid": _store(reg, aqua_ion(2), energy=-100.0),
         "base": _store(reg, hydroxo_complex(), energy=-95.0),
         "wat": _store(reg, water(), energy=-15.0),
         "h3o": _store(reg, hydronium(), energy=-14.0),
         "hoac": _store(reg, _from_smiles("CC(=O)O", "acetic acid"), energy=-50.0),
         "oac": _store(reg, _from_smiles("CC(=O)[O-]", "acetate"), energy=-49.0)}
    link_protomers(reg, water_id=s["wat"], hydronium_id=s["h3o"])
    s["rid"] = reg.conn.execute(
        "SELECT r.id FROM reactions r JOIN reaction_reagents rr ON rr.reaction_id = r.id "
        "WHERE r.product_structure_id = ? AND rr.structure_id = ?",
        (s["base"], s["acid"])).fetchone()[0]
    return s


def test_a_proton_handed_to_acetate_instead_of_water_is_exact(reg):
    """With acetate as the sink, a deprotonation step becomes AH + OAc- -> A- + HOAc:
    the step's dE is the direct equation's, the basis names HOAc, and water is gone."""
    from mofsbu.pathways.route import proton_sinks

    s = _deprotonation_with_acetate(reg)
    assert [(k["base"], k["acid"]) for k in proton_sinks(reg)] == [(s["oac"], s["hoac"])]
    walk = ([s["base"], s["acid"]], [s["rid"]])
    to_water = price_path(reg, *walk, min_fidelity=Fidelity.XTB)
    to_acetate = price_path(reg, *walk, min_fidelity=Fidelity.XTB, proton_sink=s["oac"])
    assert to_water["total_dE"] == pytest.approx(6.0)
    assert to_acetate["total_dE"] == pytest.approx((-95.0 + -50.0) - (-100.0 + -49.0))
    assert to_acetate["proton_sink"]["dE_per_proton"] == pytest.approx(-2.0)
    assert to_acetate["nodes"][1]["basis"] == f"{s['hoac']}x1"
    net = {t["structure_id"]: (t["side"], t["stoich"])
           for t in to_acetate["net_equation"]["terms"]}
    assert net == {s["base"]: ("product", 1), s["hoac"]: ("product", 1),
                   s["acid"]: ("reagent", 1), s["oac"]: ("reagent", 1)}
    assert to_acetate["net_equation"]["agrees_with_steps"] is True


def test_a_sink_the_registry_does_not_hold_is_refused_by_name(reg):
    s = _deprotonation_with_acetate(reg)
    with pytest.raises(ReferenceSchemeError, match="not a proton sink"):
        price_path(reg, [s["base"], s["acid"]], [s["rid"]], min_fidelity=Fidelity.XTB,
                   proton_sink=s["wat"])
