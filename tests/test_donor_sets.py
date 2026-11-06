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
