"""#15: donor sets of three or more, built from pockets.

The claims: a set is connected through chelate rings the pocket layer accepted, holds one
donor per group, never pairs two donors that a rigid path keeps apart, and is reported
once per symmetry class of the ligand.
"""
from __future__ import annotations

import collections

import pytest

from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.model import donor_sets, perceive

EDTA = "[O-]C(=O)CN(CC(=O)[O-])CCN(CC(=O)[O-])CC(=O)[O-]"
THQ = "O=C(C(O)=C1O)C(O)=C(O)C1=O"


def sizes(smiles: str, k: int = 6) -> dict[int, int]:
    mol = embed_molecule(mol_from_smiles(smiles), seed=7)
    return dict(collections.Counter(d.denticity for d in
                                    donor_sets(mol, perceive(mol), max_denticity=k)))


@pytest.mark.parametrize("smiles,expect", [
    ("NCCNCCN", {3: 1}),                                   # dien
    ("NCCN(CCN)CCN", {3: 1, 4: 1}),                        # tren
    ("c1ccnc(c1)-c1cccc(n1)-c1ccccn1", {3: 1}),            # terpy
    (EDTA, {3: 2, 4: 2, 5: 1, 6: 1}),
])
def test_a_polydentate_ligand_offers_one_set_per_symmetry_class(smiles, expect):
    assert sizes(smiles) == expect


@pytest.mark.parametrize("smiles", ["NCCN", "Oc1ccccc1O"])
def test_a_bidentate_ligand_offers_no_set(smiles):
    """Pairs stay `chelate_pockets`' job."""
    assert sizes(smiles) == {}


def test_a_rigid_ring_of_donors_pointing_outward_wraps_nothing():
    """tHQ's oxygens chelate in pairs; two meta oxygens on one rigid ring cannot meet."""
    assert sizes(THQ) == {}


def test_a_carboxylate_contributes_one_oxygen():
    mol = embed_molecule(mol_from_smiles(EDTA), seed=7)
    for s in donor_sets(mol, perceive(mol), max_denticity=6):
        carbons = [next(n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors())
                   for i in s.donors if s.donor_types[s.donors.index(i)] == "carboxylate_O"]
        assert len(carbons) == len(set(carbons)), s


def test_the_cap_caps():
    assert sizes(EDTA, k=4) == {3: 2, 4: 2}
    assert sizes(EDTA, k=2) == {}


# ── planned and built ────────────────────────────────────────────────────────

def _dien_spec(**over):
    from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec

    base = dict(molecules=(MoleculeSpec("dien", "NCCNCCN", 1, 0),),
                metals=(MetalSpec("Ni", 2, "hs"),), coordination="6", ligands_per_metal="1",
                binding=("chelate",), co_ligand="O", geometries=("octahedral",), seed=7,
                co_ligand_counts="fill")
    base.update(over)
    return BuildSpec(**base)


def _arrangements_planned(spec) -> list:
    from mofsbu.runner import enumerate_plan

    return sorted(t.payload.get("arrangement", -1) for t in enumerate_plan(spec).tasks
                  if t.kind == "place"
                  and any(c["denticity"] == 3 for c in t.payload["components"]))


def test_fac_and_mer_are_planned_as_two_placements():
    assert _arrangements_planned(_dien_spec()) == [0, 1]


def test_a_spec_from_before_v10_plans_pairs_only():
    from mofsbu.spec import BuildSpec

    old_dict = _dien_spec().to_dict()
    old_dict.pop("max_denticity")
    old = BuildSpec.from_dict({**old_dict, "spec_version": 9})
    assert old.max_denticity == 2
    assert _arrangements_planned(old) == []
    assert _arrangements_planned(_dien_spec(max_denticity=2)) == []


def test_both_dien_isomers_are_built(tmp_path):
    from mofsbu.registry import BlobStore, Registry
    from mofsbu.runner import run

    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as reg:
        reg.migrate()
        run(reg, _dien_spec(run_mode="construct"))
        rows = reg.conn.execute(
            "SELECT status, payload_json FROM tasks WHERE kind = 'place'").fetchall()
    import json
    built = sorted(json.loads(r["payload_json"]).get("arrangement", -1)
                   for r in rows if r["status"] == "done")
    assert built[-2:] == [0, 1], built
