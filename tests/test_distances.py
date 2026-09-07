"""Metal-donor distances are a property of the PAIR (rev 18).

The bug these tests exist for: `placer.D_ML` gave one distance per metal, so a Zn(II)
centre placed an iodide at 2.00 A — oxygen's distance, 0.6 A too short — and the clash
check then refused it.  Every halide co-ligand was rejected by construction, with a
message ("QC FAILED: 2 clash(es), closest 1.40 A") that named neither iodine nor the
distance.
"""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu.geometry.distances import (
    BASE_MO, DELTA_BY_ELEMENT, metal_donor_distance,
)
from mofsbu.geometry.embed import embed_molecule
from mofsbu.geometry.placer import LigandPlacement, place_mononuclear, site_vectors
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.model import perceive


def ligand(smiles: str, name: str | None = None) -> LigandPlacement:
    mol = embed_molecule(mol_from_smiles(smiles), seed=3)
    site = perceive(mol)[0]
    return LigandPlacement(mol=mol, donor_idxs=(site.atom_idx,),
                           donor_types=(site.donor_type,), name=name or smiles)


def bonds(result) -> list[float]:
    return [float(np.linalg.norm(result.coords[i] - result.coords[result.metal_idx]))
            for i in result.donor_idxs]


# ── the model ────────────────────────────────────────────────────────────────

def test_oxygen_is_the_base_and_nothing_moved():
    """Every geometry built before this change used an O or N donor.  Those distances
    are the base of the new model, so they must be bit-identical, not merely close."""
    assert DELTA_BY_ELEMENT["O"] == 0.0
    for metal, d in BASE_MO.items():
        assert metal_donor_distance(metal, "O").value == d


def test_heavier_donors_are_further_out():
    """The whole point: the ordering M-O < M-Cl < M-Br < M-I is chemistry, and the old
    table asserted they were all equal."""
    d = [metal_donor_distance("Zn", el).value for el in ("O", "N", "Cl", "Br", "I")]
    assert d == sorted(d)
    assert d[0] < d[-1] - 0.5                     # Zn-O 2.00 vs Zn-I 2.60


@pytest.mark.parametrize("metal,element,expected", [
    ("Zn", "O", 2.00), ("Zn", "N", 2.05), ("Zn", "S", 2.32),
    ("Zn", "Cl", 2.25), ("Zn", "I", 2.60),
    ("Cu", "N", 2.03), ("Fe", "Cl", 2.35),
])
def test_calibrated_pairs(metal, element, expected):
    """Spot values against CSD-typical divalent bond lengths.  If someone retunes the
    offsets, this is where it has to be argued rather than noticed later as drift."""
    assert metal_donor_distance(metal, element).value == pytest.approx(expected, abs=0.01)


def test_unknown_pair_is_estimated_and_says_so():
    """An estimate that does not announce itself is worse than no number: a QC verdict
    on it reads as a statement about the chemistry."""
    d = metal_donor_distance("Zn", "Xe")
    assert d.estimated and d.source == "covalent-radii"
    assert metal_donor_distance("Zn", "I").estimated is False


def test_override_wins_and_is_labelled():
    d = metal_donor_distance("Zn", "I", override=2.0)
    assert d.value == 2.0 and d.source == "override"


# ── the placer uses it ───────────────────────────────────────────────────────

def test_each_donor_is_placed_at_its_own_distance():
    """A mixed sphere has more than one M-L distance.  This is the assertion the old
    code could not have made, because there was only ever one number."""
    result = place_mononuclear(
        "Zn", [ligand("O", "H2O"), ligand("O", "H2O"),
               ligand("[I-]", "I"), ligand("[I-]", "I")],
        geometry="tetrahedral")
    d = sorted(round(x, 3) for x in bonds(result))
    assert d == [2.0, 2.0, 2.6, 2.6]


def test_halide_co_ligands_are_no_longer_refused_by_construction():
    """The regression that started this: an all-iodide sphere.  It is not a statement
    about whether Zn/I is a good idea — it is that the refusal must come from the
    chemistry and not from a table that never heard of iodine."""
    for smiles in ("[I-]", "[Cl-]", "[Br-]"):
        result = place_mononuclear("Zn", [ligand(smiles) for _ in range(4)],
                                   geometry="tetrahedral")
        assert result.ok, f"{smiles}: {result.report}"


def test_d_ml_override_still_pins_every_bond():
    """`d_ml` is what the older tests use to pin an exact length; it must still win."""
    result = place_mononuclear("Zn", [ligand("[I-]") for _ in range(4)],
                               geometry="tetrahedral", d_ml=2.0)
    assert np.allclose(bonds(result), 2.0, atol=1e-6)


def test_site_vectors_accepts_one_distance_per_vertex():
    vecs = site_vectors("tetrahedral", 4, [2.0, 2.0, 2.6, 2.6])
    assert [round(float(np.linalg.norm(v)), 3) for v in vecs] == [2.0, 2.0, 2.6, 2.6]
    with pytest.raises(ValueError):
        site_vectors("tetrahedral", 4, [2.0, 2.0])


def test_a_calibrated_sphere_is_not_flagged_as_a_guess():
    """The caveat travels with the geometry — and must not be attached to pairs that
    ARE calibrated, or it stops meaning anything.  Every donor perception can currently
    produce (O, N, S, P, Se, halide) is in the table, so a built sphere carries no note;
    the estimate path is exercised directly in `test_unknown_pair_is_estimated_and_says_so`.
    """
    result = place_mononuclear("Zn", [ligand("O") for _ in range(3)] + [ligand("[I-]")],
                               geometry="tetrahedral")
    assert not [n for n in result.report.notes if "estimated" in n]
