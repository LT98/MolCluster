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


def test_dmf_is_a_medium_and_differs_from_water(reg):
    """The solvothermal syntheses run in DMF/H2O; a water continuum alone misstates them."""
    _store(reg, water(), energy=-15.0)
    _store(reg, hydroxide(), energy=-10.0)
    assert correct_registry(reg, "alpb:dmf")["written"] == 2
    assert correct_registry(reg, "alpb:water")["written"] == 2       # its own set
    by_medium = dict(reg.conn.execute(
        "SELECT m.solvent, SUM(sc.dG_solv) FROM solvation_corrections sc JOIN methods m "
        "ON m.id = sc.method_id GROUP BY m.solvent").fetchall())
    assert set(by_medium) == {"alpb:dmf", "alpb:water"}
    assert by_medium["alpb:dmf"] != pytest.approx(by_medium["alpb:water"])


def test_the_script_splits_the_geometries_over_workers(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    db = tmp_path / "r.db"
    with Registry(db) as r:
        r.migrate()
        for g, e in ((water(), -15.0), (hydronium(), -14.0), (hydroxide(), -10.0)):
            _store(r, g, energy=e)
        r.conn.commit()
    script = Path(__file__).resolve().parents[1] / "scripts" / "solvate.py"
    out = subprocess.run([sys.executable, str(script), "--db", str(db), "--store",
                          str(tmp_path / "store"), "--medium", "alpb:dmf", "--workers", "2"],
                         capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
    assert "3 of 3 geometries corrected in alpb:dmf" in out.stdout
    with Registry(db) as r:
        assert correct_registry(r, "alpb:dmf")["candidates"] == 0


def test_a_late_worker_still_gets_its_share(tmp_path):
    """#98: the parent slices one list.  A worker started after another has committed must
    not re-query a shorter list and skip what fell between the slices."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import solvate

    db = tmp_path / "r.db"
    with Registry(db) as r:
        r.migrate()
        for g, e in ((water(), -15.0), (hydronium(), -14.0), (hydroxide(), -10.0)):
            _store(r, g, energy=e)
        r.conn.commit()
        from mofsbu.energy.solvation import uncorrected
        todo = uncorrected(r, "alpb:dmf")
    env = {"OMP_NUM_THREADS": "1"}
    first = solvate._shard(str(db), str(tmp_path / "store"), "alpb:dmf", todo[0::2], env)
    late = solvate._shard(str(db), str(tmp_path / "store"), "alpb:dmf", todo[1::2], env)
    assert first["written"] + late["written"] == 3
    with Registry(db) as r:
        assert not uncorrected(r, "alpb:dmf")
