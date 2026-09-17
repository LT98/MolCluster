"""The choice vector's key: stable, seed-free, and actually written (M5/S1).

Two layers, and the second is the one that caught a real defect: the digest recipe in
`assembly.choice`, and the registry's promise that a stored geometry carries the digest of
the vector it was built from.  Before S1 the writer read a `"digest"` key that no producer
has ever set, so every stored geometry had a NULL digest and `ix_geometries_choice`
indexed nothing.  `test_put_geometry_stores_a_digest` is the gate on that.
"""
from __future__ import annotations

import json
import math

import pytest

import fixtures as fx
from mofsbu._types import Fidelity, MethodSpec
from mofsbu.assembly.choice import (
    FLOAT_DP, ChoiceVector, ChoiceVectorError, canonical, canonical_json, digest_of)
from mofsbu.registry import (
    Provenance, Registry, backfill_choice_digests, geometries_from_choice, put_geometry,
    put_structure)
from mofsbu.versions import ALGO_VERSIONS

BUILD = MethodSpec(code="construct", code_version="0", method="raw-construct")


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r

#: The shape `place_mononuclear` actually emits, trimmed to what matters here.
PLACER_CV = {
    "metal": "Cu", "geometry": "square_planar", "cn": 4, "d_ml": 1.98, "seed": 7,
    "ligands": [
        {"ligand": "acetate", "mode": "monodentate", "donors": [1],
         "torsion_well": 0, "torsion_deg": 0.0, "azimuth_step": 2, "oop_step": 1,
         "azimuth_deg": 30.0, "oop_deg": 7.5,
         "distances": [{"d": 1.98, "source": "element-offset", "metal": "Cu",
                        "donor": "O", "estimated": False}]},
        {"ligand": "aqua", "mode": "monodentate", "donors": [0],
         "torsion_well": 0, "azimuth_step": 0, "oop_step": 0,
         "distances": [{"d": 2.0, "source": "covalent-radii", "metal": "Cu",
                        "donor": "O", "estimated": True}]},
    ],
}


# ── the recipe ───────────────────────────────────────────────────────────────

def test_key_order_is_not_a_choice():
    """Two dicts that differ only in insertion order are one choice."""
    a = {"metal": "Cu", "cn": 4, "geometry": "sq"}
    b = {"geometry": "sq", "metal": "Cu", "cn": 4}
    assert digest_of(a) == digest_of(b)
    assert ChoiceVector(a) == ChoiceVector(b)


def test_float_noise_below_the_quantum_is_not_a_new_conformer():
    """Ground rule 9: the laptop and the workstation must agree on the key."""
    noise = 10 ** -(FLOAT_DP + 6)
    assert digest_of({"d_ml": 1.98}) == digest_of({"d_ml": 1.98 + noise})
    # ...but a difference the recipe claims to keep is kept.
    assert digest_of({"d_ml": 1.98}) != digest_of({"d_ml": 1.98 + 10 ** -FLOAT_DP})


def test_negative_zero_is_zero():
    assert digest_of({"oop_deg": -0.0}) == digest_of({"oop_deg": 0.0})


def test_the_seed_is_not_part_of_the_key():
    """A seed is Kind-A noise (§6.7), not a branch — so it must not split the family.

    This is the whole reason `seed` is its own column: 'every geometry from this choice'
    has to return both samples for clustering to have anything to collapse.
    """
    one = dict(PLACER_CV, seed=1)
    two = dict(PLACER_CV, seed=99)
    assert digest_of(one) == digest_of(two)
    assert "seed" not in canonical(one)
    # Excluded from the key, not discarded: it is still recoverable for replay.
    assert ChoiceVector(one).seed == 1
    assert ChoiceVector(PLACER_CV, seed=3).seed == 3


def test_a_vector_carrying_its_own_digest_still_digests_to_that_digest():
    """Otherwise the stored-JSON round trip is not a fixed point."""
    d = digest_of(PLACER_CV)
    assert digest_of(dict(PLACER_CV, digest=d)) == d


def test_ligand_order_is_a_choice():
    """Lists are not sorted: ligand 0 and ligand 1 sit on different vertices."""
    swapped = dict(PLACER_CV, ligands=list(reversed(PLACER_CV["ligands"])))
    assert digest_of(swapped) != digest_of(PLACER_CV)


def test_tuples_and_lists_are_one_thing():
    """JSON has one sequence type, so a stored vector must survive the trip."""
    assert digest_of({"donors": (1, 2)}) == digest_of({"donors": [1, 2]})


def test_true_is_not_one():
    """`bool` is an `Integral`; coercing it would merge two values JSON keeps apart."""
    assert digest_of({"estimated": True}) != digest_of({"estimated": 1})
    assert json.loads(canonical_json({"estimated": True}))["estimated"] is True


def test_numpy_scalars_are_accepted_as_the_numbers_they_are():
    np = pytest.importorskip("numpy")
    assert digest_of({"cn": np.int64(4)}) == digest_of({"cn": 4})
    assert digest_of({"d_ml": np.float64(1.98)}) == digest_of({"d_ml": 1.98})


def test_a_failed_search_cannot_be_stored_as_a_reproducible_choice():
    for bad in (float("nan"), math.inf):
        with pytest.raises(ChoiceVectorError):
            digest_of({"azimuth_deg": bad})


def test_a_value_with_no_json_form_is_refused_by_name():
    with pytest.raises(ChoiceVectorError) as exc:
        digest_of({"ligands": {1, 2}})
    assert "ligands" in str(exc.value)          # the message names WHERE, not just what


def test_the_digest_carries_its_recipe_version():
    assert digest_of(PLACER_CV).startswith(ALGO_VERSIONS["choice_vector"] + ":")


def test_round_trip_through_storage_preserves_the_key():
    cv = ChoiceVector(PLACER_CV)
    assert ChoiceVector.from_json(cv.to_json()).digest == cv.digest


def test_a_choice_vector_is_usable_as_a_dict_key():
    """L3 dedup groups by choice before it ever measures a geometry (D11)."""
    assert len({ChoiceVector(PLACER_CV), ChoiceVector(dict(PLACER_CV, seed=42))}) == 1


def test_coerce_takes_whatever_a_producer_emits():
    assert ChoiceVector.coerce(None) is None
    assert ChoiceVector.coerce(PLACER_CV).digest == ChoiceVector(PLACER_CV).digest
    cv = ChoiceVector(PLACER_CV)
    assert ChoiceVector.coerce(cv) is cv


def test_mutating_the_caller_dict_afterwards_does_not_move_the_key():
    source = dict(PLACER_CV)
    cv = ChoiceVector(source)
    before = cv.digest
    source["metal"] = "Zn"
    assert cv.digest == before


# ── the registry, which is where the seam was ────────────────────────────────

def _xyz(jitter: float = 0.0) -> str:
    """Water, displaced — a different geometry of the SAME structure each time."""
    rows = [("O", 0.0), ("H", 0.96), ("H", 1.92)]
    return "3\nwater\n" + "".join(
        f"{el} {x + jitter:.4f} 0.0000 0.0000\n" for el, x in rows)


def test_put_geometry_stores_a_digest(reg):
    """The regression gate: the writer computes the key, the producer cannot forget it."""
    sid = put_structure(reg, fx.water(), provenance=Provenance(kind="assembly")).id
    put = put_geometry(reg, sid, _xyz(), fidelity=Fidelity.RAW, method=BUILD,
                       choice_vector=PLACER_CV)
    row = reg.conn.execute(
        "SELECT choice_vector_digest, choice_vector_json, seed FROM geometries WHERE id=?",
        (put.id,)).fetchone()
    assert row["choice_vector_digest"] == digest_of(PLACER_CV)
    assert row["seed"] == 7                     # lifted from the vector, not lost
    assert ChoiceVector.from_json(row["choice_vector_json"]).digest == digest_of(PLACER_CV)


def test_the_same_choice_at_several_seeds_is_one_findable_family(reg):
    """What `ix_geometries_choice` is for, and what clustering will consume."""
    sid = put_structure(reg, fx.water(), provenance=Provenance(kind="assembly")).id
    for seed in (1, 2, 3):
        put_geometry(reg, sid, _xyz(0.1 * seed), fidelity=Fidelity.RAW, method=BUILD,
                     choice_vector=dict(PLACER_CV, seed=seed))
    found = geometries_from_choice(reg, PLACER_CV)
    assert len(found) == 3
    assert sorted(r["seed"] for r in found) == [1, 2, 3]


def test_backfill_keys_old_rows_and_leaves_keyless_ones_alone(reg):
    """The 99-rows-with-a-vector-and-no-key case, reproduced."""
    sid = put_structure(reg, fx.water(), provenance=Provenance(kind="assembly")).id
    with_cv = put_geometry(reg, sid, _xyz(), fidelity=Fidelity.RAW, method=BUILD,
                           choice_vector=PLACER_CV).id
    without = put_geometry(reg, sid, _xyz(0.5), fidelity=Fidelity.RAW, method=BUILD).id
    # Put the database back into its pre-S1 state: a vector, no key, no seed.
    reg.conn.execute(
        "UPDATE geometries SET choice_vector_digest=NULL, seed=NULL, choice_vector_json=? "
        "WHERE id=?", (json.dumps(PLACER_CV), with_cv))
    reg.conn.commit()

    assert backfill_choice_digests(reg, dry_run=True) == 1
    assert reg.conn.execute(
        "SELECT choice_vector_digest FROM geometries WHERE id=?",
        (with_cv,)).fetchone()[0] is None       # dry run wrote nothing

    assert backfill_choice_digests(reg) == 1
    row = reg.conn.execute(
        "SELECT choice_vector_digest, seed FROM geometries WHERE id=?",
        (with_cv,)).fetchone()
    assert row["choice_vector_digest"] == digest_of(PLACER_CV)
    assert row["seed"] == 7
    # A geometry with no vector has nothing to derive a key from and keeps its NULL.
    assert reg.conn.execute(
        "SELECT choice_vector_digest FROM geometries WHERE id=?",
        (without,)).fetchone()[0] is None
    assert backfill_choice_digests(reg) == 0        # idempotent
