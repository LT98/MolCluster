"""M1: the descriptor tables, and the coverage rule that keeps them honest.

The gate here is `test_every_perceivable_donor_type_has_a_row`.  `sites.perception` can
emit a donor type from three places — the labile SMARTS list, `_classify_anionic` and
`_classify_neutral` — and a type with no descriptor row is a site the ease model will
later refuse to score, discovered at that point rather than this one.  Adding a donor
type without adding its descriptors now fails here, which is the same trick
`test_donor_patterns.py` uses on the SMARTS.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from mofsbu.assembly.join import NotBuiltYet
from mofsbu.descriptors import (
    UnknownDescriptor, activation_ease, donor, donor_table, hsab_match,
    known_donor_types, metal, metal_table, sync_to_registry,
)
from mofsbu._types import Fidelity
from mofsbu.registry import Registry
from mofsbu.sites import perception
from mofsbu.sites.frames import live_dof
from mofsbu.sites.model import Site

SRC = Path(perception.__file__)


def perceivable_donor_types() -> set[str]:
    """Every donor type the perception code can produce, read out of the source.

    Reading the source rather than importing a list is deliberate: the three producers
    are a data table, a chain of `return "..."` in `_classify_anionic`, and an f-string
    for halides.  A hand-maintained union would be a fourth place to forget.
    """
    text = SRC.read_text(encoding="utf-8")
    types = {name for name, _, _ in perception.LABILE_DONOR_PATTERNS}
    types |= set(re.findall(r'return "(\w+_\w+)"', text))
    if 'return f"halide_{sym}"' in text:
        types |= {f"halide_{x}" for x in ("F", "Cl", "Br", "I")}
    return types


# ── the coverage gate ────────────────────────────────────────────────────────

def test_every_perceivable_donor_type_has_a_row():
    missing = sorted(perceivable_donor_types() - known_donor_types())
    assert not missing, (
        f"donor types perception can emit with no descriptor row: {missing}. "
        f"Add them to data/reference/donor_descriptors.tsv with a source, or the ease "
        f"model will refuse to score those sites later instead of now.")


def test_no_orphan_descriptor_rows():
    """A row for a type nothing can perceive is dead weight that reads as coverage."""
    orphans = sorted(known_donor_types() - perceivable_donor_types())
    assert not orphans, f"descriptor rows for donor types nothing emits: {orphans}"


# ── provenance ───────────────────────────────────────────────────────────────

def test_every_row_carries_a_source():
    for d in donor_table().values():
        assert d.source, f"{d.donor_type} has no source"
        assert d.source_version, f"{d.donor_type} has no source_version"
    for m in metal_table().values():
        assert m.source and m.source_version, f"{m.symbol}{m.charge:+d} lacks provenance"


def test_a_pka_always_carries_its_spread():
    """A curated estimate and a tabulated value must not look alike."""
    for d in donor_table().values():
        if d.pka is not None:
            assert d.pka_sigma is not None and d.pka_sigma > 0, d.donor_type


def test_estimates_are_labelled_as_estimates():
    """`estimate` in the source column is the honesty mechanism; it must be in use."""
    sources = {d.source for d in donor_table().values()}
    assert "estimate" in sources, (
        "no row is marked `estimate` — either every value is tabulated, which is not "
        "credible for this donor set, or the labelling has been dropped")


# ── the meaning of a missing pKa ─────────────────────────────────────────────

def test_neutral_donors_have_no_pka_and_say_so():
    """Empty pKa means 'no activation step', not 'we forgot'."""
    for donor_type in ("amine_N", "pyridyl_N", "ether_O", "carbonyl_O", "aqua_O",
                       "nitrile_N", "thioether_S", "phosphine_P", "imine_N"):
        d = donor(donor_type)
        assert d.pka is None, f"{donor_type} has a pKa; it has no proton to lose"
        assert d.needs_activation is False


def test_anionic_donors_carry_their_conjugate_acid_pka():
    assert donor("carboxylate_O").pka == pytest.approx(4.76, abs=0.5)
    assert donor("phenolate_O").pka == pytest.approx(10.0, abs=0.5)
    assert donor("thiolate_S").pka == pytest.approx(10.6, abs=0.5)
    for donor_type in ("carboxylate_O", "phenolate_O", "thiolate_S", "alkoxide_O"):
        assert donor(donor_type).needs_activation is True


def test_the_pka_column_runs_in_one_direction_only():
    """Sulfonate is a strong acid, alkoxide a very weak one; the order proves the sign.

    If the column ever quietly acquired conjugate-acid basicities for neutral donors,
    this ordering would stop meaning anything.
    """
    assert donor("sulfonate_O").pka < donor("carboxylate_O").pka
    assert donor("carboxylate_O").pka < donor("phenolate_O").pka
    assert donor("phenolate_O").pka < donor("alkoxide_O").pka


# ── one home for the torsion model ───────────────────────────────────────────

def test_live_dof_comes_from_frames_not_from_the_file():
    for d in donor_table().values():
        assert d.live_dof is live_dof(d.donor_type)
    assert "live_dof" not in Path(
        "data/reference/donor_descriptors.tsv").read_text().split("\n")[-40].split("\t")


# ── lookups refuse rather than default ───────────────────────────────────────

def test_unknown_donor_raises_and_names_the_file():
    with pytest.raises(UnknownDescriptor, match="donor_descriptors.tsv"):
        donor("unobtanium_Q")


def test_unknown_ion_points_at_the_generator_not_the_output():
    with pytest.raises(UnknownDescriptor, match="build_metal_descriptors"):
        metal("Fe", 7)


# ── the metal half ───────────────────────────────────────────────────────────

def test_metal_rows_are_chemically_sane():
    fe3 = metal("Fe", 3)
    assert fe3.d_electrons == 5 and fe3.hsab == "hard"
    assert fe3.ionic_radius == pytest.approx(0.645, abs=0.01)      # Shannon HS CN6
    cu2 = metal("Cu", 2)
    assert cu2.exchange_lability == "fast"                          # Jahn-Teller
    assert metal("Cr", 3).exchange_lability == "slow"               # the inert benchmark
    assert metal("Zn", 2).ionic_radius < metal("Cd", 2).ionic_radius


def test_main_group_ions_have_no_d_count():
    for symbol, charge in (("Al", 3), ("Mg", 2), ("Ca", 2)):
        assert metal(symbol, charge).d_electrons is None


def test_preferred_cn_is_parsed_as_numbers():
    assert metal("Cu", 2).preferred_cn == (4, 5, 6)
    assert metal("Fe", 3).preferred_cn == (6,)


# ── into the registry ────────────────────────────────────────────────────────

def test_sync_writes_both_tables_and_is_idempotent(tmp_path):
    with Registry(tmp_path / "r.db") as reg:
        reg.migrate()
        n_d, n_m = sync_to_registry(reg)
        assert n_d == len(donor_table()) and n_m == len(metal_table())
        sync_to_registry(reg)                       # a correction must replace, not stack
        rows = reg.conn.execute("SELECT COUNT(*) c FROM donor_descriptors").fetchone()["c"]
        assert rows == n_d
        stored = reg.conn.execute(
            "SELECT live_dof FROM donor_descriptors WHERE donor_type='carboxylate_O'"
        ).fetchone()["live_dof"]
        assert stored == live_dof("carboxylate_O").value


# ── C5 is resolved (D18); C7 is not, and the difference is visible ──────────

def test_the_partner_free_floor_runs_and_the_partner_term_still_refuses():
    """D18 resolved C5 only.  A caller who asks for a partner must be told no.

    Returning the intrinsic number under a partner-shaped call is the failure this
    guards: the answer would be indistinguishable from a partner-aware one, and the
    whole point of C7's factorization is that those are different quantities.
    """
    site = Site(atom_idx=0, donor_type="carboxylate_O", labile=True, charge_after=-1,
                live_dof="live", binding_modes=("mono",))
    record = activation_ease(site)
    assert record.fidelity is Fidelity.HEURISTIC and record.scored

    with pytest.raises(NotBuiltYet, match="C7"):
        activation_ease(site, partner=metal("Cu", 2))
    with pytest.raises(NotBuiltYet, match="C7"):
        hsab_match(donor("carboxylate_O"), metal("Cu", 2))
