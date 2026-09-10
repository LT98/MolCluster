"""M7: the backend protocol, the spin convention, and what a MethodSpec must carry.

Nothing here needs xTB or MACE installed — that is the point of the protocol.  The two
tests that do need them are marked and skip cleanly, so this file is the same gate on
the laptop and on the workstation.
"""
from __future__ import annotations

import pytest

from mofsbu._types import EnergyBackendUnavailable, Fidelity, MethodSpec
from mofsbu.assembly.join import NotBuiltYet
from mofsbu.energy.backends import (
    ML_BACKENDS, MACEBackend, MACEOmolBackend, NullBackend, XTBBackend, available_backends,
    backend_for, check_spin, d_electrons, electron_count, high_spin_multiplicity,
    minimal_multiplicity, ml_backend_key,
)
from mofsbu.energy.relax import (
    MODE_FIDELITY, available_modes, ml_model_status, mode_status, relax_geometry,
)

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])


# ── spin ─────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("symbol,charge,d,mult", [
    ("Fe", 3, 5, 6),        # the Fe(III) sextet the archived run used, as a rule not an env var
    ("Fe", 2, 6, 5),        # high-spin d6 quintet
    ("Ni", 2, 8, 3),        # d8 triplet — NOT the minimal-spin singlet the Ni run assumed
    ("Cu", 2, 9, 2),
    ("Zn", 2, 10, 1),       # d10 closed shell: high spin and low spin coincide
    ("Mn", 2, 5, 6),
    ("Cr", 3, 3, 4),
])
def test_high_spin_multiplicity(symbol, charge, d, mult):
    assert d_electrons(symbol, charge) == d
    assert high_spin_multiplicity(symbol, charge) == mult


def test_high_spin_and_minimal_spin_disagree_for_nickel():
    """The two conventions are different answers, and the code says which one it used.

    The archived Ni(II) numbers are minimal-spin (singlet); the Fe(III) ones are high-spin
    (sextet).  That is a real difference in what was computed, and it now lives in the
    MethodSpec stored with each energy instead of in an environment variable that was set
    in a shell six months ago.
    """
    symbols = ["Ni"]
    assert minimal_multiplicity(symbols, 2) == 1
    assert high_spin_multiplicity("Ni", 2) == 3


def test_impossible_ion_is_refused():
    with pytest.raises(ValueError, match="not a real ion"):
        d_electrons("Ni", 12)
    with pytest.raises(ValueError, match="not a transition metal"):
        d_electrons("C", 0)


def test_multiplicity_must_match_the_electron_count():
    """An even electron count cannot be a doublet, and xTB will not tell you so."""
    assert electron_count(["O", "H", "H"], 0) == 10
    check_spin(["O", "H", "H"], 0, 1)
    with pytest.raises(ValueError, match="cannot give multiplicity"):
        check_spin(["O", "H", "H"], 0, 2)
    check_spin(["O", "H", "H"], 1, 2)          # 9 electrons: doublet is right


# ── the protocol ─────────────────────────────────────────────────────────────

def test_null_backend_is_deterministic_and_declares_itself():
    b = NullBackend()
    a = b.single_point(*WATER, charge=0, multiplicity=1)
    again = b.single_point(*WATER, charge=0, multiplicity=1)
    assert a.energy == again.energy
    # It must be impossible to mistake one of these for a real number in the registry.
    assert a.method.code == "null"
    assert "not-physical" in a.method.method


def test_every_backend_reports_a_method_spec_with_charge_and_spin():
    for backend in (NullBackend(), XTBBackend(), MACEBackend()):
        spec = backend.method_spec(charge=-2, multiplicity=3)
        assert isinstance(spec, MethodSpec)
        assert spec.charge == -2 and spec.multiplicity == 3
        assert spec.code_version, f"{backend.name} must pin a code version (ground rule 6)"


def test_charge_blind_backends_say_so_in_the_stored_method():
    """MACE-MP-0 cannot see formal charge; the method row carries that, not a docstring."""
    assert MACEBackend().method_spec(charge=-3, multiplicity=1).extras["charge_blind"] is True
    assert MACEBackend().method_spec(charge=-3, multiplicity=1).extras["spin_blind"] is True
    assert "charge_blind" not in XTBBackend().method_spec(charge=-3, multiplicity=1).extras


def test_mace_refuses_a_solvent_rather_than_ignoring_it():
    with pytest.raises(ValueError, match="no solvation model"):
        MACEBackend().method_spec(charge=0, multiplicity=1, solvent="water")
    with pytest.raises(ValueError, match="no solvation model"):
        MACEOmolBackend().method_spec(charge=0, multiplicity=1, solvent="water")


# ── two MACE models on one rung ──────────────────────────────────────────────
# MACE-OMOL-0 is not "a newer MACE".  It is trained on OMol25 and takes total charge and
# spin multiplicity as inputs, which moves it across the line the reference scheme draws.
# These tests exist because the thing that would hurt is not the model failing to load —
# it is the two models being treated as one theory by anything downstream.

def test_omol_is_not_charge_blind_and_mp_is():
    """The flag is the whole point: same package, same rung, different theory."""
    assert MACEBackend().charge_aware is False
    assert MACEOmolBackend().charge_aware is True
    assert MACEOmolBackend().spin_aware is True
    omol = MACEOmolBackend().method_spec(charge=-3, multiplicity=6)
    assert "charge_blind" not in omol.extras
    assert "spin_blind" not in omol.extras


def test_the_two_mace_models_are_never_the_same_theory():
    """`same_theory` gates every subtraction; it must separate these two."""
    mp = MACEBackend().method_spec(charge=0, multiplicity=1)
    omol = MACEOmolBackend().method_spec(charge=0, multiplicity=1)
    assert mp.code == omol.code == "mace"          # same package ...
    assert mp.method != omol.method                # ... different theory
    assert not mp.same_theory(omol)


def test_a_stored_method_names_the_training_set():
    """`code='mace'` alone cannot say which reference the number is against."""
    assert "Materials Project" in MACEBackend().method_spec(
        charge=0, multiplicity=1).extras["training_set"]
    assert "OMol25" in MACEOmolBackend().method_spec(
        charge=0, multiplicity=1).extras["training_set"]


def test_omol_checks_the_multiplicity_it_is_handed():
    """A model that is GIVEN a spin state must be given a reachable one."""
    with pytest.raises(ValueError, match="cannot give multiplicity"):
        MACEOmolBackend().single_point(*WATER, charge=0, multiplicity=2)


def test_the_ml_rung_is_a_declared_choice_not_a_lookup(monkeypatch):
    monkeypatch.delenv("MOFSBU_ML_MODEL", raising=False)
    assert backend_for(Fidelity.ML).method == "MACE-MP-0"
    assert backend_for(Fidelity.ML, ml_model="mace-omol-0").method == "MACE-OMOL-0"
    monkeypatch.setenv("MOFSBU_ML_MODEL", "mace-omol-0")
    assert backend_for(Fidelity.ML).method == "MACE-OMOL-0"
    # the spec still wins over the machine
    assert backend_for(Fidelity.ML, ml_model="mace-mp-0").method == "MACE-MP-0"


def test_an_unknown_ml_model_refuses_rather_than_defaulting():
    """A typo that silently selected the charge-blind model would be invisible."""
    with pytest.raises(ValueError, match="unknown ML model"):
        ml_backend_key("mace-omol-1")        # close, and not a model that exists
    with pytest.raises(ValueError, match="unknown ML model"):
        backend_for(Fidelity.ML, ml_model="omol25")


def test_omol_reports_an_old_mace_as_not_installed_not_as_a_crash(monkeypatch):
    """mace-torch 0.3.6 imports fine and has no `mace_omol`; that is a missing INSTALL.

    Guarded on `mace.calculators`, NOT on `mace`.  The top-level package imports without
    torch and the calculators subpackage does not, so a machine with a half-installed ML
    stack (mace present, torch absent) passed `importorskip("mace")` and then died on the
    next line — a collection error where a skip was meant.  Which is the same mistake the
    code under test exists to prevent, made in the test that tests it: "importable" and
    "usable" are different claims.
    """
    pytest.importorskip("torch", reason="the ML stack is not installed on this machine")
    calculators = pytest.importorskip(
        "mace.calculators", reason="mace-torch not installed on this machine")

    monkeypatch.delattr(calculators, "mace_omol", raising=False)
    assert MACEOmolBackend().available() is False
    assert "0.3.14" in MACEOmolBackend().install_hint()


def test_ml_model_status_names_both_models(monkeypatch):
    monkeypatch.delenv("MOFSBU_ML_MODEL", raising=False)
    status = ml_model_status()
    assert set(status) == set(ML_BACKENDS)
    assert status["mace"]["declared"] is True and status["mace_omol"]["declared"] is False
    assert status["mace_omol"]["charge_aware"] is True
    for entry in status.values():
        if not entry["available"]:
            assert entry["note"], "an unavailable model must say what would fix it"


def test_xtb_refuses_an_unknown_solvent():
    with pytest.raises(ValueError, match="unknown solvent"):
        XTBBackend().method_spec(charge=0, multiplicity=1, solvent="unobtanium")


def test_same_theory_ignores_charge_and_spin_but_not_the_medium():
    """Two species in one equation differ in charge; that must not block the comparison."""
    gas_a = XTBBackend().method_spec(charge=0, multiplicity=1)
    gas_b = XTBBackend().method_spec(charge=-3, multiplicity=6)
    water = XTBBackend().method_spec(charge=0, multiplicity=1, solvent="water")
    assert gas_a.same_theory(gas_b)
    assert not gas_a.same_theory(water)


# ── availability is a property of the machine, not of the code ───────────────

def test_dft_is_a_missing_body_not_a_missing_install():
    """Ground rule 8: the two are different failures and get different exceptions."""
    with pytest.raises(NotBuiltYet, match="DFT"):
        backend_for(Fidelity.DFT)


def test_raw_and_ff_have_no_backend_because_they_are_not_energies():
    for rung in (Fidelity.RAW, Fidelity.FF):
        with pytest.raises(NotBuiltYet):
            backend_for(rung)


def test_available_modes_covers_every_declared_run_mode():
    modes = available_modes()
    assert set(modes) == set(MODE_FIDELITY)
    assert modes["construct"] is True
    assert modes["dft_go"] is False            # no external code is wired up


def test_mode_status_explains_a_disabled_mode():
    """A disabled option must carry its reason, so the UI cannot invent one."""
    for mode, entry in mode_status().items():
        if not entry["available"]:
            assert entry["note"], f"{mode} is disabled with no reason given"


def test_relax_refuses_an_unavailable_backend_with_an_install_hint():
    if XTBBackend().available():
        pytest.skip("xtb is installed here; the refusal path needs it absent")
    with pytest.raises(EnergyBackendUnavailable, match="conda install"):
        relax_geometry(WATER[1], WATER[0], charge=0, multiplicity=1, target=Fidelity.XTB)


def test_available_backends_answers_for_all_four():
    got = available_backends()
    assert set(got) == {"xtb", "mace", "mace_omol", "null"}
    assert got["null"] is True


# ── the real thing, when the stack is there ──────────────────────────────────

def test_xtb_relaxes_water():
    pytest.importorskip("tblite.ase", reason="tblite not installed on this machine")
    result = XTBBackend().relax(*WATER, charge=0, multiplicity=1, fmax=0.05, steps=100)
    assert result.converged
    assert result.energy <= (result.initial_energy or 0.0) + 1e-9
    assert result.method.method == "GFN2-xTB"
    assert result.to_xyz("water").splitlines()[0] == "3"


def test_xtb_pins_a_real_version_not_the_word_unknown():
    """A methods row saying `tblite-unknown` cannot distinguish two builds (ground rule 6)."""
    pytest.importorskip("tblite.ase", reason="tblite not installed on this machine")
    assert XTBBackend().code_version() != "unknown"


def test_xtb_sees_the_solvent():
    pytest.importorskip("tblite.ase", reason="tblite not installed on this machine")
    b = XTBBackend()
    gas = b.single_point(*WATER, charge=0, multiplicity=1)
    aq = b.single_point(*WATER, charge=0, multiplicity=1, solvent="water")
    assert gas.energy != aq.energy
    assert aq.method.solvent == "water"
