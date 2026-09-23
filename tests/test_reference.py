"""M7's real gate: the reference scheme refuses subtractions it cannot justify.

The case under test is the one that actually went wrong.  `legacy/xtb_energy_metal.py`
computed E(complex) - E(M^q+) - E(ligand anion), an equation that balances in atoms and
NOT in charge, and the resulting Ni/EDTA sequestration verdict reversed once a medium was
added (`docs/reports/solvation_report.md`).  Every test here is a way of asking: would
this module have let that number be produced silently?
"""
from __future__ import annotations

import sqlite3
from dataclasses import replace

import pytest

from mofsbu._types import Fidelity, MethodSpec, ReferenceSchemeError, split_medium
from mofsbu.energy.backends import NullBackend
from mofsbu.energy.reference import (
    BARE_ION, CHARGE_SEPARATION, COORDINATION_CHANGE, NAKED_POLYANION, check_balance,
    check_reference_quality, put_balanced_reaction, reaction_balanced_energy,
    Term, quality_for_terms, reaction_terms, store_reaction_energy,
)
from mofsbu.energy.protons import deprotonation_pairs, link_protomers
from mofsbu.energy.routes import decompositions, price_incoming_routes, price_reaction
from mofsbu.graph import EdgeType, TypedGraph
from mofsbu.registry import (
    Provenance, ReadOnlyRegistry, Registry, RegistryError, get_graph, put_geometry,
    put_reaction, put_solvation_correction, put_structure, solvation_corrections,
)

XTB = MethodSpec(code="tblite", code_version="0.4.0", method="GFN2-xTB")
XTB_WATER = MethodSpec(code="tblite", code_version="0.4.0", method="GFN2-xTB",
                       solvent="water")
MACE = MethodSpec(code="mace", code_version="0.3.6/medium", method="MACE-MP-0",
                  extras={"charge_blind": True, "spin_blind": True})
# Same package, same rung of the ladder, different theory — and, unlike MP-0, given the
# total charge and spin multiplicity, so the reference scheme has no reason to refuse it.
OMOL = MethodSpec(code="mace", code_version="0.3.14/extra_large", method="MACE-OMOL-0")


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


# ── small hand-built species: a metal, a monodentate donor, and their complex ──

def water(charge: int = 0) -> TypedGraph:
    g = TypedGraph(charge=charge, multiplicity=1, name="H2O")
    o = g.add_atom("O")
    for _ in range(2):
        g.add_bond(o, g.add_atom("H"), EdgeType.COVALENT)
    return g


def hydroxide() -> TypedGraph:
    g = TypedGraph(charge=-1, multiplicity=1, name="OH-")
    o = g.add_atom("O", formal_charge=-1)
    g.add_bond(o, g.add_atom("H"), EdgeType.COVALENT)
    return g


def hydronium() -> TypedGraph:
    g = TypedGraph(charge=1, multiplicity=1, name="H3O+")
    o = g.add_atom("O", formal_charge=1)
    for _ in range(3):
        g.add_bond(o, g.add_atom("H"), EdgeType.COVALENT)
    return g


def bare_ion() -> TypedGraph:
    g = TypedGraph(charge=2, multiplicity=3, name="Ni2+")
    g.add_atom("Ni", oxidation_state=2, spin_class="hs")
    return g


def aqua_ion(n: int = 2) -> TypedGraph:
    """[Ni(H2O)n]2+ — the reference species that replaces the bare ion."""
    g = TypedGraph(charge=2, multiplicity=3, name=f"Ni(H2O){n}")
    m = g.add_atom("Ni", oxidation_state=2, spin_class="hs")
    for _ in range(n):
        o = g.add_atom("O")
        g.add_bond(o, m, EdgeType.DATIVE)
        for _ in range(2):
            g.add_bond(o, g.add_atom("H"), EdgeType.COVALENT)
    return g


def hydroxo_complex() -> TypedGraph:
    """[Ni(OH)(H2O)]+ — one water deprotonated on the metal."""
    g = TypedGraph(charge=1, multiplicity=3, name="Ni(OH)(H2O)+")
    m = g.add_atom("Ni", oxidation_state=2, spin_class="hs")
    oh = g.add_atom("O", formal_charge=-1)
    g.add_bond(oh, m, EdgeType.DATIVE)
    g.add_bond(oh, g.add_atom("H"), EdgeType.COVALENT)
    o = g.add_atom("O")
    g.add_bond(o, m, EdgeType.DATIVE)
    for _ in range(2):
        g.add_bond(o, g.add_atom("H"), EdgeType.COVALENT)
    return g


def _store(reg, g, *, energy, method=XTB, fidelity=Fidelity.XTB, converged=True):
    put = put_structure(reg, g)
    lines = [str(len(g)), g.name]
    for k, i in enumerate(g.nodes()):
        lines.append(f"{g.label(i).element} {k:.4f} 0.0000 0.0000")
    put_geometry(reg, put.id, "\n".join(lines) + "\n", fidelity=fidelity, method=method,
                 energy=energy, converged=converged)
    return put.id


# ── balance ──────────────────────────────────────────────────────────────────

def test_the_legacy_equation_balances_and_is_still_refused(reg):
    """E(complex) - E(bare ion) - E(free anion) BALANCES.  That is exactly the trap.

    The archived report calls this scheme "charge-conserving" and it is: Ni(2+) + OH(-)
    + H2O has the same atoms and the same total charge as [Ni(OH)(H2O)](+).  A balance
    check alone passes it.  What must catch it is the second question — are the two sides
    comparable chemistry? — because the reference side has a naked cation, a free anion,
    and not one metal-donor bond where the product has two.
    """
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    ion = _store(reg, bare_ion(), energy=-40.0)
    anion = _store(reg, hydroxide(), energy=-20.0)
    wat = _store(reg, water(), energy=-15.0)

    rid = put_balanced_reaction(reg, product,
                                reagents=[(ion, 1), (anion, 1), (wat, 1)])
    assert check_balance(reg, rid).balanced, "the legacy equation does balance"

    quality = check_reference_quality(reg, rid)
    assert not quality.isodesmic
    # Ni2+ + OH- pair up: the ion count changes too, which gas phase gets wrong by eV.
    assert {i.code for i in quality.issues} == {BARE_ION, COORDINATION_CHANGE,
                                                CHARGE_SEPARATION}
    assert quality.dative_delta == 2          # two bonds appear out of nothing

    with pytest.raises(ReferenceSchemeError, match="not isodesmic"):
        reaction_balanced_energy(reg, rid)

    # Reproducing the legacy number on purpose is allowed, and labelled.
    deliberate = reaction_balanced_energy(reg, rid, strict=False)
    assert deliberate.dE == pytest.approx(-100.0 - (-40.0 + -20.0 + -15.0))
    assert deliberate.quality.isodesmic is False
    assert "NOT isodesmic" in deliberate.describe()


def test_a_naked_polyanion_is_flagged(reg):
    """The E(BTC(3-)) term: charge with no metal under it."""
    g = TypedGraph(charge=-3, multiplicity=1, name="L3-")
    c = g.add_atom("C")
    for _ in range(3):
        g.add_bond(c, g.add_atom("O", formal_charge=-1), EdgeType.COVALENT)
    trianion = _store(reg, g, energy=-50.0)

    prod = TypedGraph(charge=-1, multiplicity=3, name="NiL-")
    m = prod.add_atom("Ni", oxidation_state=2, spin_class="hs")
    pc = prod.add_atom("C")
    for k in range(3):
        o = prod.add_atom("O", formal_charge=-1)
        prod.add_bond(pc, o, EdgeType.COVALENT)
        if k == 0:
            prod.add_bond(o, m, EdgeType.DATIVE)
    complexed = _store(reg, prod, energy=-160.0)
    ion = _store(reg, bare_ion(), energy=-40.0)

    rid = put_balanced_reaction(reg, complexed, reagents=[(ion, 1), (trianion, 1)])
    codes = {i.code for i in check_reference_quality(reg, rid).issues}
    assert NAKED_POLYANION in codes and BARE_ION in codes


def test_ligand_exchange_conserves_coordination(reg):
    """The shape the scheme wants: the metal is coordinated on both sides."""
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    quality = check_reference_quality(reg, rid)
    assert quality.isodesmic, quality.describe()
    assert quality.dative_delta == 0


def test_a_proton_balanced_exchange_is_accepted(reg):
    """[Ni(H2O)2]2+ -> [Ni(OH)(H2O)]+ + H3O+ ... balanced only once the proton has a carrier.

    Written with a bare H+ it does not balance; written with the proton on a water it
    does.  That is the whole discipline: no species is left to fend for itself.
    """
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)

    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    report = check_balance(reg, rid)
    assert report.balanced, report.describe()
    assert report.charge_delta == 0 and report.element_delta == {}


def test_balance_report_names_what_is_missing(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    with pytest.raises(ReferenceSchemeError) as exc:
        put_balanced_reaction(reg, product, reagents=[(aqua, 1)])
    message = str(exc.value)
    assert "H" in message and "charge" in message


def test_an_unbalanced_reaction_is_not_left_behind(reg):
    """A rejected write must not leave a half-inserted row for the next reader to find."""
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    before = reg.conn.execute("SELECT COUNT(*) c FROM reactions").fetchone()["c"]
    with pytest.raises(ReferenceSchemeError):
        put_balanced_reaction(reg, product, reagents=[(aqua, 1)])
    after = reg.conn.execute("SELECT COUNT(*) c FROM reactions").fetchone()["c"]
    assert after == before
    assert reg.conn.execute("SELECT COUNT(*) c FROM reaction_reagents").fetchone()["c"] == 0


def test_leaving_groups_sit_on_the_product_side(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    sides = {t.role: t.side for t in reaction_terms(reg, rid)}
    assert sides["leaving"] == "product"
    assert sides["reagent"] == "reagent"


# ── the energy ───────────────────────────────────────────────────────────────

@pytest.fixture()
def balanced(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    return rid, {"product": product, "aqua": aqua, "water": wat, "h3o": h3o}


def test_energy_is_products_minus_reagents(reg, balanced):
    rid, _ = balanced
    result = reaction_balanced_energy(reg, rid, strict=False)
    assert result.dE == pytest.approx((-100.0 + -14.0) - (-60.0 + -15.0))
    assert result.method.method == "GFN2-xTB"
    assert result.fidelity is Fidelity.XTB
    assert result.balance.balanced


def test_stoichiometry_is_honoured(reg):
    """2 H2O -> H3O+ + OH-: the coefficients have to reach the arithmetic."""
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    oh = _store(reg, hydroxide(), energy=-13.0)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    energy = reaction_balanced_energy(reg, rid, strict=False)
    assert energy.dE == pytest.approx((-14.0 + -13.0) - 2 * -15.0)


def test_mixed_levels_of_theory_are_refused(reg, balanced):
    """One species relaxed in water and the rest in gas phase is not a reaction energy."""
    rid, ids = balanced
    put_geometry(reg, ids["water"], "3\nH2O-aq\nO 9.0 0.0 0.0\nH 9.9 0.0 0.0\nH 8.7 0.9 0.0\n",
                 fidelity=Fidelity.XTB, method=XTB_WATER, energy=-16.0, converged=True)
    # Asking for the aqueous medium now finds water solvated and everything else absent.
    with pytest.raises(ReferenceSchemeError, match="no energy"):
        reaction_balanced_energy(reg, rid, strict=False, solvent="water")


def test_a_missing_energy_is_named_not_skipped(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    bare = put_structure(reg, water()).id          # stored, but never given an energy
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (bare, 1)],
                                leaving=[(h3o, 1)])
    with pytest.raises(ReferenceSchemeError, match="has no energy"):
        reaction_balanced_energy(reg, rid, strict=False)


def test_the_null_backend_cannot_score_a_reaction_by_accident(reg):
    null = NullBackend().method_spec(charge=0, multiplicity=1)
    wat = _store(reg, water(), energy=-15.0, method=null, fidelity=Fidelity.RAW)
    h3o = _store(reg, hydronium(), energy=-14.0, method=null, fidelity=Fidelity.RAW)
    oh = _store(reg, hydroxide(), energy=-13.0, method=null, fidelity=Fidelity.RAW)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    with pytest.raises(ReferenceSchemeError, match="non-physical"):
        reaction_balanced_energy(reg, rid, strict=False)
    energy = reaction_balanced_energy(reg, rid, strict=False, allow_null=True)
    assert energy.dE == pytest.approx(-27.0 + 30.0)


def test_a_charge_blind_backend_is_refused_on_a_charged_equation(reg):
    """MACE cannot see the difference between H3O+ and OH-; it must not be asked to."""
    wat = _store(reg, water(), energy=-15.0, method=MACE, fidelity=Fidelity.ML)
    h3o = _store(reg, hydronium(), energy=-14.0, method=MACE, fidelity=Fidelity.ML)
    oh = _store(reg, hydroxide(), energy=-13.0, method=MACE, fidelity=Fidelity.ML)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    with pytest.raises(ReferenceSchemeError, match="charge-blind"):
        reaction_balanced_energy(reg, rid, strict=False)


def test_a_charge_blind_backend_is_fine_on_a_neutral_equation(reg):
    """2 H2O -> H2O + H2O is trivial, neutral, and legal for an MLIP."""
    wat = _store(reg, water(), energy=-15.0, method=MACE, fidelity=Fidelity.ML)
    rid = put_balanced_reaction(reg, wat, reagents=[(wat, 2)], leaving=[(wat, 1)])
    assert reaction_balanced_energy(reg, rid).dE == pytest.approx(0.0)


def test_a_charge_aware_mlip_is_accepted_on_a_charged_equation(reg):
    """MACE-OMOL-0 is handed the charge, so the charge-blind refusal must not fire.

    This is the point of the whole exercise: the refusal keys on what the STORED method
    says it could see, not on the string "mace".  Getting that wrong in either direction
    is expensive — refusing OMOL-0 leaves the charged species at xTB cost forever, and
    accepting MP-0 produces a number that looks the same and is not.
    """
    wat = _store(reg, water(), energy=-15.0, method=OMOL, fidelity=Fidelity.ML)
    h3o = _store(reg, hydronium(), energy=-14.0, method=OMOL, fidelity=Fidelity.ML)
    oh = _store(reg, hydroxide(), energy=-13.0, method=OMOL, fidelity=Fidelity.ML)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    energy = reaction_balanced_energy(reg, rid, strict=False)
    assert energy.dE == pytest.approx(-14.0 - 13.0 - 2 * -15.0)
    assert energy.method.method == "MACE-OMOL-0"


def test_two_mlips_on_one_rung_are_never_mixed_into_one_equation(reg):
    """The failure this guards is silent: both are fidelity=ML, so nothing sorts them.

    Water has BOTH an MP-0 and an OMOL-0 energy — which is the normal state of a
    registry once a second model is in use.  The old rule ("highest rung, then lowest
    energy") would have picked whichever model had the deeper reference, for every
    species, and produced an equation assembled out of two theories with no complaint.
    The equation is now pinned to the theory of its first term.
    """
    wat = _store(reg, water(), energy=-15.0, method=OMOL, fidelity=Fidelity.ML)
    _store(reg, water(), energy=-3000.0, method=MACE, fidelity=Fidelity.ML)
    h3o = _store(reg, hydronium(), energy=-14.0, method=OMOL, fidelity=Fidelity.ML)
    oh = _store(reg, hydroxide(), energy=-13.0, method=OMOL, fidelity=Fidelity.ML)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    energy = reaction_balanced_energy(reg, rid, strict=False)
    # -3000 is the numerically lowest energy in the table and must not appear.
    assert energy.dE == pytest.approx(-14.0 - 13.0 - 2 * -15.0)
    assert energy.method.method == "MACE-OMOL-0"


def test_a_species_missing_from_the_pinned_theory_is_named_not_substituted(reg):
    """Having SOME ML energy is not having one in the equation's theory."""
    wat = _store(reg, water(), energy=-15.0, method=OMOL, fidelity=Fidelity.ML)
    h3o = _store(reg, hydronium(), energy=-14.0, method=OMOL, fidelity=Fidelity.ML)
    oh = _store(reg, hydroxide(), energy=-13.0, method=MACE, fidelity=Fidelity.ML)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    with pytest.raises(ReferenceSchemeError, match="MACE-OMOL-0"):
        reaction_balanced_energy(reg, rid, strict=False)


def test_the_weakest_rung_sets_the_fidelity(reg):
    """An xTB product minus a raw-construct reagent is a raw-construct number."""
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0, fidelity=Fidelity.RAW)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    assert reaction_balanced_energy(reg, rid, strict=False).fidelity is Fidelity.RAW


def test_unconverged_geometries_are_flagged_and_can_be_refused(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0, converged=False)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    assert reaction_balanced_energy(reg, rid, strict=False).all_converged is False
    with pytest.raises(ReferenceSchemeError, match="upper bound"):
        reaction_balanced_energy(reg, rid, strict=False, allow_unconverged=False)


# ── writing it back ──────────────────────────────────────────────────────────

def test_a_non_isodesmic_energy_cannot_be_stored_silently(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    ion = _store(reg, bare_ion(), energy=-40.0)
    anion = _store(reg, hydroxide(), energy=-20.0)
    wat = _store(reg, water(), energy=-15.0)
    rid = put_balanced_reaction(reg, product, reagents=[(ion, 1), (anion, 1), (wat, 1)])
    result = reaction_balanced_energy(reg, rid, strict=False)
    with pytest.raises(ReferenceSchemeError, match="refusing to store"):
        store_reaction_energy(reg, rid, result)
    store_reaction_energy(reg, rid, result, force=True)
    note = reg.conn.execute("SELECT note FROM reactions WHERE id=?", (rid,)).fetchone()["note"]
    assert "NOT-ISODESMIC" in note and BARE_ION in note


def test_storing_an_energy_records_its_method(reg, balanced):
    """Ground rule 3, at the reaction level: no bare float in the dG column."""
    rid, _ = balanced
    result = reaction_balanced_energy(reg, rid, strict=False)
    store_reaction_energy(reg, rid, result, force=True)     # gas-phase charge separation
    row = reg.conn.execute("SELECT dG, method_id, fidelity FROM reactions WHERE id=?",
                           (rid,)).fetchone()
    assert row["dG"] == pytest.approx(result.dE)
    assert row["method_id"] is not None
    assert row["fidelity"] == int(Fidelity.XTB)
    method = reg.conn.execute("SELECT * FROM methods WHERE id=?",
                              (row["method_id"],)).fetchone()
    assert method["method"] == "GFN2-xTB"


# ── pricing an edge for a reader ─────────────────────────────────────────────
#
# `reaction_balanced_energy` answers by raising.  A panel listing a dozen routes needs
# the same answer as data, one entry per edge, and needs it without a writable registry.

def test_a_priced_route_carries_the_caveat_with_the_number(reg, balanced):
    rid, _ = balanced
    route = price_reaction(reg, rid, min_fidelity=Fidelity.XTB)
    assert route["can_price"] is True
    assert route["total_dE"] == pytest.approx(route["steps"][0]["dE"])
    step = route["steps"][0]
    assert step["fidelity_name"] == "XTB"
    assert "GFN2-xTB" in step["method"]
    assert step["isodesmic"] is True          # this fixture is a proton-balanced exchange


def test_the_floor_admits_a_higher_rung(reg, balanced):
    """The regression that matters: a floor, not a rung.

    `_energy_row` matches a requested fidelity EXACTLY, so pricing must resolve the
    equation at whatever rung it has and compare afterwards.  Asking the reference
    scheme for `ML` would refuse this xTB equation, which is the opposite of a floor.
    """
    rid, _ = balanced
    route = price_reaction(reg, rid, min_fidelity=Fidelity.ML)
    assert route["can_price"] is True
    assert route["steps"][0]["fidelity"] == int(Fidelity.XTB)


def test_a_rung_below_the_floor_is_refused_and_the_rung_is_named(reg):
    MMFF = MethodSpec(code="rdkit", code_version="2026.03", method="ETKDGv3+MMFF")
    product = _store(reg, hydroxo_complex(), energy=-100.0, method=MMFF,
                     fidelity=Fidelity.FF)
    aqua = _store(reg, aqua_ion(2), energy=-60.0, method=MMFF, fidelity=Fidelity.FF)
    wat = _store(reg, water(), energy=-15.0, method=MMFF, fidelity=Fidelity.FF)
    h3o = _store(reg, hydronium(), energy=-14.0, method=MMFF, fidelity=Fidelity.FF)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    route = price_reaction(reg, rid, min_fidelity=Fidelity.ML)
    assert route["can_price"] is False
    assert route["total_dE"] is None
    # Naming the rung that IS there is the difference between "no" and "not yet".
    assert "FF" in route["why_not"] and "ML" in route["why_not"]


def test_an_unpriceable_route_still_reports_its_diagnosis(reg):
    """A refusal is a result: the caveats survive the number failing to appear.

    Built through `put_reaction`, which records what was constructed and does not check
    balance — the path every assembly edge in a real run comes from, and the reason a
    provenance panel meets unbalanced equations at all.
    """
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    ion = _store(reg, bare_ion(), energy=-40.0)
    rid = put_reaction(reg, product, Provenance(kind="assembly", reagent_ids=(ion,)))
    route = price_reaction(reg, rid)
    assert route["can_price"] is False
    assert "unbalanced" in route["why_not"]
    step = route["steps"][0]
    assert step["isodesmic"] is False
    assert COORDINATION_CHANGE in step["caveats"]
    # The unpriced shape carries every key the priced one does, so no caller branches.
    assert set(step) == {"reaction_id", "dE", "fidelity", "fidelity_name", "method",
                         "isodesmic", "caveats", "all_converged", "terms",
                         "can_price", "why_not"}


def test_pricing_needs_no_writable_registry(reg, balanced, tmp_path):
    """The viewer holds `PRAGMA query_only = ON`; `Registry` would write to probe it."""
    rid, _ = balanced
    reg.conn.commit()
    con = sqlite3.connect(f"file:{reg.db_path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA query_only = ON")
    try:
        reader = ReadOnlyRegistry(conn=con, store=reg.store)
        route = price_reaction(reader, rid, min_fidelity=Fidelity.XTB)
        assert route["can_price"] is True
        assert con.execute("PRAGMA query_only").fetchone()[0] == 1
    finally:
        con.close()


def test_every_incoming_edge_gets_an_entry(reg, balanced):
    rid, ids = balanced
    routes = price_incoming_routes(reg, ids["product"], min_fidelity=Fidelity.XTB)
    assert [r["reaction_id"] for r in routes] == [rid]
    assert routes[0]["kind"] == "reaction"


# ── what this could have been made from ──────────────────────────────────────

def test_decompositions_finds_a_split_that_was_never_recorded(reg):
    """A `place` edge records a construction and cites nothing; the pieces still exist.

    `[M(H2O)2]2+` is `[M(H2O)]2+` plus a water, and both are in the registry — so the
    question "what adds up to this" has an answer even though no edge asserts it.
    """
    product = _store(reg, aqua_ion(2), energy=-100.0)
    rung = _store(reg, aqua_ion(1), energy=-80.0)
    wat = _store(reg, water(), energy=-15.0)

    found = decompositions(reg, product, min_fidelity=Fidelity.XTB)
    split = next(d for d in found
                 if sorted(p["structure_id"] for p in d["parts"]) == sorted([rung, wat]))
    assert split["origin"] == "inferred"
    assert split["recorded_as"] is None           # nothing wrote this down
    assert split["can_price"] is True
    assert split["total_dE"] == pytest.approx(-100.0 - (-80.0 + -15.0))


def test_a_derived_split_says_when_it_is_already_recorded(reg):
    product = _store(reg, aqua_ion(2), energy=-100.0)
    rung = _store(reg, aqua_ion(1), energy=-80.0)
    wat = _store(reg, water(), energy=-15.0)
    rid = put_reaction(reg, product, Provenance(reagent_ids=(rung, wat)))

    split = next(d for d in decompositions(reg, product, min_fidelity=Fidelity.XTB)
                 if sorted(p["structure_id"] for p in d["parts"]) == sorted([rung, wat]))
    assert split["recorded_as"] == rid


def test_protomers_of_one_molecule_are_found_by_arithmetic(reg):
    """`[M(H2O)2]2+` and `[M(OH)(H2O)]+` differ by one proton and one charge."""
    acid = _store(reg, aqua_ion(2), energy=-100.0)
    base = _store(reg, hydroxo_complex(), energy=-95.0)
    assert (acid, base, 1) in deprotonation_pairs(reg)


def test_a_deprotonation_edge_is_isodesmic_and_charge_separating(reg):
    """No metal-donor bond changes across a proton transfer, so the bond-type rules have
    nothing to object to — which is the whole reason the couple exists rather than a bare
    H+.  It still separates one ion into two, so in gas phase `strict` refuses it (C18)."""
    acid = _store(reg, aqua_ion(2), energy=-100.0)
    base = _store(reg, hydroxo_complex(), energy=-95.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)

    written = link_protomers(reg, water_id=wat, hydronium_id=h3o)
    edge = next(w for w in written
                if (w["protonated"], w["deprotonated"]) == (acid, base))
    assert edge["why_not"] is None

    with pytest.raises(ReferenceSchemeError, match=CHARGE_SEPARATION):
        reaction_balanced_energy(reg, edge["reaction_id"], strict=True)
    energy = reaction_balanced_energy(reg, edge["reaction_id"], strict=False)
    assert energy.quality.isodesmic is True
    assert [i.code for i in energy.quality.issues] == [CHARGE_SEPARATION]
    assert energy.dE == pytest.approx((-95.0 + -14.0) - (-100.0 + -15.0))


def test_linking_twice_writes_nothing_the_second_time(reg):
    _store(reg, aqua_ion(2), energy=-100.0)
    _store(reg, hydroxo_complex(), energy=-95.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)

    first = link_protomers(reg, water_id=wat, hydronium_id=h3o)
    assert any(w["reaction_id"] is not None for w in first)
    assert link_protomers(reg, water_id=wat, hydronium_id=h3o) == []


def test_the_couple_is_not_linked_through_itself(reg):
    """`H3O+ + H2O -> H2O + H3O+` is balanced, isodesmic and exactly zero."""
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    assert (h3o, wat, 1) in deprotonation_pairs(reg)          # the arithmetic sees it
    assert link_protomers(reg, water_id=wat, hydronium_id=h3o) == []   # and declines it


def test_a_split_below_the_floor_is_refused_not_hidden(reg):
    """Absent is absent: a candidate with no ML energy is listed with its reason."""
    product = _store(reg, aqua_ion(2), energy=-100.0)
    _store(reg, aqua_ion(1), energy=-80.0)
    _store(reg, water(), energy=-15.0)
    found = decompositions(reg, product, min_fidelity=Fidelity.DFT)
    assert found                                   # still offered, not silently dropped
    assert all(not d["can_price"] for d in found)
    assert all("XTB" in d["why_not"] and "DFT" in d["why_not"] for d in found)


# ── a medium on a stored energy (C16–C18, docs/WORKPLAN_solvation.md) ─────────

ALPB = "alpb:water"


def _xtb_alpb(charge, multiplicity, version="0.7.0"):
    return MethodSpec(code="tblite", code_version=version, method="GFN2-xTB", solvent=ALPB,
                      charge=charge, multiplicity=multiplicity)


def _geometry_of(reg, sid):
    return int(reg.conn.execute(
        "SELECT best_geometry_id FROM structures WHERE id=?", (sid,)).fetchone()[0])


def _correct(reg, sid, dG, *, version="0.7.0"):
    g = get_graph(reg, sid)
    return put_solvation_correction(
        reg, _geometry_of(reg, sid), method=_xtb_alpb(g.charge, g.multiplicity, version),
        e_gas=-10.0, e_solv=-10.0 + dG)


@pytest.fixture()
def deprotonation(reg):
    """[Ni(H2O)2]2+ + H2O -> [Ni(OH)(H2O)]+ + H3O+ on OMOL-0: one ion becomes two."""
    ids = {"product": _store(reg, hydroxo_complex(), energy=-100.0, method=OMOL,
                             fidelity=Fidelity.ML),
           "aqua": _store(reg, aqua_ion(2), energy=-60.0, method=OMOL, fidelity=Fidelity.ML),
           "water": _store(reg, water(), energy=-15.0, method=OMOL, fidelity=Fidelity.ML),
           "h3o": _store(reg, hydronium(), energy=-14.0, method=OMOL, fidelity=Fidelity.ML)}
    rid = put_balanced_reaction(reg, ids["product"], reagents=[(ids["aqua"], 1),
                                                               (ids["water"], 1)],
                                leaving=[(ids["h3o"], 1)])
    return rid, ids


DG = {"product": -2.0, "aqua": -6.0, "water": -0.4, "h3o": -4.4}


def test_charge_separation_blocks_in_gas_and_travels_as_a_caveat_in_a_continuum(
        reg, deprotonation):
    rid, _ = deprotonation
    gas = check_reference_quality(reg, rid)
    assert gas.isodesmic                              # no bond-type rule is broken
    assert gas.ion_delta == 1 and not gas.acceptable
    assert [(i.code, i.blocking) for i in gas.issues] == [(CHARGE_SEPARATION, True)]
    wet = check_reference_quality(reg, rid, medium=ALPB)
    assert wet.acceptable
    assert [(i.code, i.blocking) for i in wet.issues] == [(CHARGE_SEPARATION, False)]
    assert "caveat" in wet.describe()


def test_the_rule_reads_only_the_terms_so_a_net_equation_is_weighed_the_same(reg):
    """A composed route's net equation: the carrier on both sides cancels, the ions count."""
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    oh = _store(reg, hydroxide(), energy=-13.0)
    through = (Term(h3o, 1, "reagent", "reagent"), Term(h3o, 1, "leaving", "product"))
    assert quality_for_terms(reg, through).ion_delta == 0
    net = (Term(wat, 2, "reagent", "reagent"), Term(h3o, 1, "leaving", "product"),
           Term(oh, 1, "leaving", "product"))
    assert quality_for_terms(reg, net).ion_delta == 2
    assert quality_for_terms(reg, net[::-1]).ion_delta == 2      # order is irrelevant


def test_a_gas_phase_charge_separation_needs_force_to_be_stored(reg, deprotonation):
    rid, _ = deprotonation
    with pytest.raises(ReferenceSchemeError, match="not usable in gas phase"):
        reaction_balanced_energy(reg, rid)
    result = reaction_balanced_energy(reg, rid, strict=False)
    with pytest.raises(ReferenceSchemeError, match="refusing to store"):
        store_reaction_energy(reg, rid, result)
    store_reaction_energy(reg, rid, result, force=True)
    note = reg.conn.execute("SELECT note FROM reactions WHERE id=?", (rid,)).fetchone()["note"]
    assert f"CAVEAT:{CHARGE_SEPARATION}" in note and "NOT-ISODESMIC" not in note


def test_a_corrected_medium_adds_each_species_correction_on_its_own_geometry(
        reg, deprotonation):
    rid, ids = deprotonation
    for name, sid in ids.items():
        assert _correct(reg, sid, DG[name]).created
    energy = reaction_balanced_energy(reg, rid, solvent=ALPB)       # strict: a caveat now
    gas = (-100.0 + -14.0) - (-60.0 + -15.0)
    assert energy.dE == pytest.approx(gas + DG["product"] + DG["h3o"] - DG["aqua"]
                                      - DG["water"])
    assert energy.method.method == "MACE-OMOL-0" and energy.method.solvent == ALPB
    assert energy.method.extras["solvation_correction"]["method"] == "GFN2-xTB"
    assert "(caveat)" in energy.describe()


def test_a_species_without_a_correction_is_named_not_left_in_gas(reg, deprotonation):
    rid, ids = deprotonation
    for name in ("product", "aqua", "h3o"):
        _correct(reg, ids[name], DG[name])
    with pytest.raises(ReferenceSchemeError, match=f"{ids['water']} .*no energy.*{ALPB}"):
        reaction_balanced_energy(reg, rid, solvent=ALPB)


def test_a_direct_and_a_corrected_energy_are_never_mixed(reg, deprotonation):
    rid, ids = deprotonation
    for name in ("product", "aqua", "h3o"):
        _correct(reg, ids[name], DG[name])
    put_geometry(reg, ids["water"], "3\nH2O-aq\nO 9.0 0.0 0.0\nH 9.9 0.0 0.0\nH 8.7 0.9 0.0\n",
                 fidelity=Fidelity.XTB, method=_xtb_alpb(0, 1), energy=-16.0, converged=True)
    with pytest.raises(ReferenceSchemeError, match="by one route"):
        reaction_balanced_energy(reg, rid, solvent=ALPB)


def test_the_correction_theory_is_pinned_across_the_equation(reg, deprotonation):
    rid, ids = deprotonation
    for name in ("product", "aqua", "h3o"):
        _correct(reg, ids[name], DG[name])
    _correct(reg, ids["water"], DG["water"], version="0.6.0")
    with pytest.raises(ReferenceSchemeError, match="pinned to the first term's correction"):
        reaction_balanced_energy(reg, rid, solvent=ALPB)


def test_a_correction_is_never_borrowed_from_a_sibling_geometry(reg):
    """The energy and its correction come from ONE geometry, whichever that is."""
    wat = _store(reg, water(), energy=-15.0, method=OMOL, fidelity=Fidelity.ML)
    first = _geometry_of(reg, wat)
    second = put_geometry(reg, wat, "3\nH2O-b\nO 5.0 0.0 0.0\nH 5.9 0.0 0.0\nH 4.7 0.9 0.0\n",
                          fidelity=Fidelity.ML, method=OMOL, energy=-14.5,
                          converged=True).id
    put_solvation_correction(reg, second, method=_xtb_alpb(0, 1), e_gas=-10.0, e_solv=-10.3)
    assert _geometry_of(reg, wat) == first != second    # the correction is not on the best
    rid = put_balanced_reaction(reg, wat, reagents=[(wat, 2)], leaving=[(wat, 1)])
    energy = reaction_balanced_energy(reg, rid, solvent=ALPB)
    assert {round(e, 9) for _, _, e in energy.contributions} == {round(-14.5 - 0.3, 9)}


def test_a_correction_states_its_model_charge_and_starting_point(reg):
    wat = _store(reg, water(), energy=-15.0, method=OMOL, fidelity=Fidelity.ML)
    gid = _geometry_of(reg, wat)
    bare = MethodSpec(code="tblite", code_version="0.7.0", method="GFN2-xTB",
                      solvent="water", charge=0, multiplicity=1)
    with pytest.raises(RegistryError, match="model:solvent"):
        put_solvation_correction(reg, gid, method=bare, e_gas=-1.0, e_solv=-1.4)
    with pytest.raises(RegistryError, match="gas phase"):
        put_solvation_correction(reg, gid, method=replace(bare, solvent=None),
                                 e_gas=-1.0, e_solv=-1.4)
    with pytest.raises(RegistryError, match="charge"):
        put_solvation_correction(reg, gid, method=_xtb_alpb(1, 1), e_gas=-1.0, e_solv=-1.4)
    with pytest.raises(RegistryError, match="must be stated"):
        put_solvation_correction(reg, gid, method=replace(_xtb_alpb(0, 1), charge=None),
                                 e_gas=-1.0, e_solv=-1.4)
    solvated = put_geometry(reg, wat, "3\nH2O-aq\nO 9.0 0.0 0.0\nH 9.9 0.0 0.0\nH 8.7 0.9 0.0\n",
                            fidelity=Fidelity.XTB, method=_xtb_alpb(0, 1), energy=-16.0).id
    with pytest.raises(RegistryError, match="counts solvation twice"):
        put_solvation_correction(reg, solvated, method=_xtb_alpb(0, 1), e_gas=-1.0,
                                 e_solv=-1.4)


def test_a_repeated_correction_is_one_row_and_a_disagreeing_one_is_refused(reg):
    wat = _store(reg, water(), energy=-15.0, method=OMOL, fidelity=Fidelity.ML)
    gid = _geometry_of(reg, wat)
    first = put_solvation_correction(reg, gid, method=_xtb_alpb(0, 1), e_gas=-1.0,
                                     e_solv=-1.4)
    again = put_solvation_correction(reg, gid, method=_xtb_alpb(0, 1), e_gas=-1.0,
                                     e_solv=-1.4)
    assert first.created and not again.created and first.id == again.id
    [row] = solvation_corrections(reg, gid, ALPB)
    assert row["dG_solv"] == pytest.approx(-0.4) and row["medium"] == ALPB
    assert solvation_corrections(reg, gid, "gbsa:water") == []          # absent, not zero
    with pytest.raises(RegistryError, match="disagrees"):
        put_solvation_correction(reg, gid, method=_xtb_alpb(0, 1), e_gas=-1.0, e_solv=-1.5)


def test_a_medium_must_name_its_model():
    assert split_medium("ALPB:Water") == ("alpb", "water")
    for token in ("water", ":water", "alpb:"):
        with pytest.raises(ValueError, match="model:solvent"):
            split_medium(token)
