"""Donor perception and the frame model, on the Zn/THQ system."""
from __future__ import annotations

from collections import Counter

import numpy as np
import pytest

from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import from_smiles, mol_from_smiles
from mofsbu.sites.frames import LiveDOF, live_dof, torsion_wells
from mofsbu.sites.model import chelate_pockets, perceive
from mofsbu.sites.perception import deprotonate, find_donor_sites

THQ = "OC1=C(O)C(=O)C(O)=C(O)C1=O"


def donors(smiles: str) -> dict[str, int]:
    return dict(Counter(s.donor_type for s in find_donor_sites(mol_from_smiles(smiles))))


def test_thq_donors_need_the_two_added_types():
    """The legacy patterns saw nothing on THQ: its ring is not aromatic and C=O was
    never a donor.  Both halves of its chelate pocket depend on the additions."""
    assert donors(THQ) == {"enol_O": 4, "carbonyl_O": 2}


def test_thq_ring_is_not_aromatic_which_is_why_phenolate_misses_it():
    mol = mol_from_smiles(THQ)
    assert not any(a.GetIsAromatic() for a in mol.GetAtoms())


@pytest.mark.parametrize("smiles,expected", [
    ("O", {"aqua_O": 1}),
    ("N", {"amine_N": 1}),
    ("[NH4+]", {}),                                   # ammonium donates nothing
    ("c1ccncc1", {"pyridyl_N": 1}),
    ("c1cc[nH]c1", {"azolate_N": 1}),
    ("CC(C)=O", {"carbonyl_O": 1}),
    ("OC(=O)c1cc(C(O)=O)cc(C(O)=O)c1", {"carboxylate_O": 6}),
])
def test_perception_across_donor_chemistry(smiles, expected):
    assert donors(smiles) == expected


def test_hydrogen_counting_survives_explicit_hydrogens():
    """`GetTotalNumHs()` counts only implicit Hs, so on these molecules it is always 0.

    Reading it directly turned 'is this still protonated?' into 'no, never' — water
    stopped being a donor and an N-H heteroaromatic would have passed as pyridyl.
    """
    assert donors("O") == {"aqua_O": 1}
    assert "pyridyl_N" not in donors("c1cc[nH]c1")


def test_activation_does_not_destroy_the_site_it_creates():
    """A recalled structure is usually already deprotonated; its anionic donor must
    still be perceived, or 'activate = deprotonate' loses the site."""
    mol = embed_molecule(mol_from_smiles(THQ))
    enols = [s for s in find_donor_sites(mol) if s.donor_type == "enol_O"]
    anion, _ = deprotonate(mol, enols[:1])
    found = Counter(s.donor_type for s in find_donor_sites(anion))
    assert found["enolate_O"] == 1
    assert found["enol_O"] == 3


def test_chelate_pockets_are_described_structurally_not_named():
    """Ring size, how many donors are anionic, and which types — no chemical labels."""
    mol = embed_molecule(mol_from_smiles(THQ))
    pockets = chelate_pockets(mol, perceive(mol))
    assert pockets
    for pocket in pockets:
        assert 4 <= pocket.ring_size <= 7
        assert set(pocket.donor_types) <= {"enol_O", "enolate_O", "carbonyl_O"}
        assert pocket.descriptor.startswith(f"ring{pocket.ring_size}/")


def test_frames_are_orthonormal_and_carry_a_live_dof():
    mol = embed_molecule(mol_from_smiles(THQ))
    for site in perceive(mol):
        axis = np.array(site.frame["axis"])
        ref = np.array(site.frame["ref"])
        assert abs(np.linalg.norm(axis) - 1) < 1e-6
        assert abs(float(np.dot(axis, ref))) < 1e-6
        assert site.live_dof in (LiveDOF.TORSION_LIVE.value, LiveDOF.TORSION_FREE.value)


def test_only_live_torsions_branch():
    """The live-DOF tag is the guard against combinatorial explosion."""
    assert len(torsion_wells("aqua_O")) == 1
    assert live_dof("aqua_O") is LiveDOF.TORSION_FREE
    assert len(torsion_wells("carboxylate_O")) == 2
    assert live_dof("carboxylate_O") is LiveDOF.TORSION_LIVE


def test_ammonium_keeps_its_charge_on_the_nitrogen():
    """Localised charge stays on the atom; delocalised charge goes to the graph (D15)."""
    nh4 = from_smiles("[NH4+]", multiplicity=1, name="ammonium")
    assert nh4.localised_charge() == 1 and nh4.net_charge() == 1
    thq = from_smiles(THQ, multiplicity=1, name="thq")
    assert thq.localised_charge() == 0
