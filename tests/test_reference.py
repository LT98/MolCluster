"""M7's real gate: the reference scheme refuses subtractions it cannot justify.

The case under test is the one that actually went wrong.  `legacy/xtb_energy_metal.py`
computed E(complex) - E(M^q+) - E(ligand anion), an equation that balances in atoms and
NOT in charge, and the resulting Ni/EDTA sequestration verdict reversed once a medium was
added (`docs/reports/solvation_report.md`).  Every test here is a way of asking: would
this module have let that number be produced silently?
"""
from __future__ import annotations

import pytest

from mofsbu._types import Fidelity, MethodSpec, ReferenceSchemeError
from mofsbu.energy.backends import NullBackend
from mofsbu.energy.reference import (
    BARE_ION, COORDINATION_CHANGE, NAKED_POLYANION, check_balance,
    check_reference_quality, put_balanced_reaction, reaction_balanced_energy,
    reaction_terms, store_reaction_energy,
)
from mofsbu.graph import EdgeType, TypedGraph
from mofsbu.registry import Registry, put_geometry, put_structure

XTB = MethodSpec(code="tblite", code_version="0.4.0", method="GFN2-xTB")
XTB_WATER = MethodSpec(code="tblite", code_version="0.4.0", method="GFN2-xTB",
                       solvent="water")
MACE = MethodSpec(code="mace", code_version="0.3.6/medium", method="MACE-MP-0",
                  extras={"charge_blind": True})


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
    assert {i.code for i in quality.issues} == {BARE_ION, COORDINATION_CHANGE}
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
    result = reaction_balanced_energy(reg, rid)
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
    assert reaction_balanced_energy(reg, rid).dE == pytest.approx((-14.0 + -13.0) - 2 * -15.0)


def test_mixed_levels_of_theory_are_refused(reg, balanced):
    """One species relaxed in water and the rest in gas phase is not a reaction energy."""
    rid, ids = balanced
    put_geometry(reg, ids["water"], "3\nH2O-aq\nO 9.0 0.0 0.0\nH 9.9 0.0 0.0\nH 8.7 0.9 0.0\n",
                 fidelity=Fidelity.XTB, method=XTB_WATER, energy=-16.0, converged=True)
    # Asking for the aqueous medium now finds water solvated and everything else absent.
    with pytest.raises(ReferenceSchemeError, match="no energy"):
        reaction_balanced_energy(reg, rid, solvent="water")


def test_a_missing_energy_is_named_not_skipped(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    bare = put_structure(reg, water()).id          # stored, but never given an energy
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (bare, 1)],
                                leaving=[(h3o, 1)])
    with pytest.raises(ReferenceSchemeError, match="has no energy"):
        reaction_balanced_energy(reg, rid)


def test_the_null_backend_cannot_score_a_reaction_by_accident(reg):
    null = NullBackend().method_spec(charge=0, multiplicity=1)
    wat = _store(reg, water(), energy=-15.0, method=null, fidelity=Fidelity.RAW)
    h3o = _store(reg, hydronium(), energy=-14.0, method=null, fidelity=Fidelity.RAW)
    oh = _store(reg, hydroxide(), energy=-13.0, method=null, fidelity=Fidelity.RAW)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    with pytest.raises(ReferenceSchemeError, match="non-physical"):
        reaction_balanced_energy(reg, rid)
    assert reaction_balanced_energy(reg, rid, allow_null=True).dE == pytest.approx(-27.0 + 30.0)


def test_a_charge_blind_backend_is_refused_on_a_charged_equation(reg):
    """MACE cannot see the difference between H3O+ and OH-; it must not be asked to."""
    wat = _store(reg, water(), energy=-15.0, method=MACE, fidelity=Fidelity.ML)
    h3o = _store(reg, hydronium(), energy=-14.0, method=MACE, fidelity=Fidelity.ML)
    oh = _store(reg, hydroxide(), energy=-13.0, method=MACE, fidelity=Fidelity.ML)
    rid = put_balanced_reaction(reg, h3o, reagents=[(wat, 2)], leaving=[(oh, 1)])
    with pytest.raises(ReferenceSchemeError, match="charge-blind"):
        reaction_balanced_energy(reg, rid)


def test_a_charge_blind_backend_is_fine_on_a_neutral_equation(reg):
    """2 H2O -> H2O + H2O is trivial, neutral, and legal for an MLIP."""
    wat = _store(reg, water(), energy=-15.0, method=MACE, fidelity=Fidelity.ML)
    rid = put_balanced_reaction(reg, wat, reagents=[(wat, 2)], leaving=[(wat, 1)])
    assert reaction_balanced_energy(reg, rid).dE == pytest.approx(0.0)


def test_the_weakest_rung_sets_the_fidelity(reg):
    """An xTB product minus a raw-construct reagent is a raw-construct number."""
    product = _store(reg, hydroxo_complex(), energy=-100.0)
    aqua = _store(reg, aqua_ion(2), energy=-60.0, fidelity=Fidelity.RAW)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    assert reaction_balanced_energy(reg, rid).fidelity is Fidelity.RAW


def test_unconverged_geometries_are_flagged_and_can_be_refused(reg):
    product = _store(reg, hydroxo_complex(), energy=-100.0, converged=False)
    aqua = _store(reg, aqua_ion(2), energy=-60.0)
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    rid = put_balanced_reaction(reg, product, reagents=[(aqua, 1), (wat, 1)],
                                leaving=[(h3o, 1)])
    assert reaction_balanced_energy(reg, rid).all_converged is False
    with pytest.raises(ReferenceSchemeError, match="upper bound"):
        reaction_balanced_energy(reg, rid, allow_unconverged=False)


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
    result = reaction_balanced_energy(reg, rid)
    store_reaction_energy(reg, rid, result)
    row = reg.conn.execute("SELECT dG, method_id, fidelity FROM reactions WHERE id=?",
                           (rid,)).fetchone()
    assert row["dG"] == pytest.approx(result.dE)
    assert row["method_id"] is not None
    assert row["fidelity"] == int(Fidelity.XTB)
    method = reg.conn.execute("SELECT * FROM methods WHERE id=?",
                              (row["method_id"],)).fetchone()
    assert method["method"] == "GFN2-xTB"
