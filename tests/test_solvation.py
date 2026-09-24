"""The continuum-correction producer (WORKPLAN_solvation S2): computed, stored, resumable."""
from __future__ import annotations

import pytest

from mofsbu.energy.backends import get_backend
from mofsbu.energy.reference import Term, energy_for_terms
from mofsbu.energy.solvation import correct_registry
from mofsbu.registry import Registry

from test_reference import _store, hydronium, hydroxide, water  # noqa: E402

pytestmark = pytest.mark.skipif(not get_backend("xtb").available(),
                                reason="tblite not installed")


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


def test_every_energy_gets_a_correction_once_and_prices_in_the_medium(reg):
    wat = _store(reg, water(), energy=-15.0)
    h3o = _store(reg, hydronium(), energy=-14.0)
    oh = _store(reg, hydroxide(), energy=-10.0)

    first = correct_registry(reg, "alpb:water")
    assert (first["candidates"], first["written"], first["refused"]) == (3, 3, {})
    assert correct_registry(reg, "alpb:water")["candidates"] == 0      # resumable

    # Autoionisation, 2 H2O -> H3O+ + OH-: two ions from two neutrals.
    terms = (Term(h3o, 1, "product", "product"), Term(oh, 1, "leaving", "product"),
             Term(wat, 2, "reagent", "reagent"))
    gas = energy_for_terms(reg, terms, fidelity=None, strict=False, subject="x")
    solv = energy_for_terms(reg, terms, fidelity=None, strict=False, subject="x",
                            solvent="alpb:water")
    # Ions are stabilised by the continuum far more than neutrals: the gap closes.
    assert solv.dE < gas.dE


def test_a_medium_with_no_producer_is_refused_by_name(reg):
    with pytest.raises(Exception, match="no producer"):
        correct_registry(reg, "gbsa:water")
