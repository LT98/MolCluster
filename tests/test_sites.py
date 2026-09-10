"""Donor perception and the frame model, on the Zn/THQ system."""
from __future__ import annotations

from collections import Counter

import numpy as np
import pytest

from build_routes import ROUTE_CASES, build, donor_set, routes
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


# ── perception does not depend on the route that reached the structure ───────

@pytest.mark.parametrize("name,smiles,donor_type,denticity", ROUTE_CASES,
                         ids=[c[0] for c in ROUTE_CASES])
def test_every_build_route_to_one_identity_perceives_the_same_donors(
        name, smiles, donor_type, denticity):
    """The seam `registry.api.catalog_drift` was reporting, closed at its source.

    Binding a metal through one oxygen of a delocalised group forces `to_rdkit` to render
    the OTHER oxygen as the anion, so the two routes to one identity arrive as two
    resonance forms.  D15 hashes them the same on purpose; perception used to read the
    bond orders and disagree, so a structure's donor set depended on which build got
    there first.  Atom indices are directly comparable here — see `build_routes`.
    """
    mol, all_routes = routes(smiles, donor_type, denticity)
    assert len(all_routes) > 1, f"{name}: needs at least two routes to compare"
    perceived = {r: donor_set(build(mol, r, donor_type)) for r in all_routes}
    first = perceived[all_routes[0]]
    for route, found in perceived.items():
        assert found == first, (
            f"{name}: binding through {route} perceives {sorted(found)}, binding through "
            f"{all_routes[0]} perceives {sorted(first)}. Same identity, different catalog.")


@pytest.mark.parametrize("name,smiles,donor_type,denticity", ROUTE_CASES,
                         ids=[c[0] for c in ROUTE_CASES])
def test_a_coordinated_group_keeps_every_oxygen_it_had_when_free(
        name, smiles, donor_type, denticity):
    """Binding one oxygen must not delete the others from the catalog.

    The original symptom: an acetate bound through its C=O oxygen perceived ONE
    carboxylate donor instead of two, because the bound oxygen no longer looked like
    anything the per-atom valence rules recognised.
    """
    mol, all_routes = routes(smiles, donor_type, denticity)
    free = sum(1 for s in perceive(mol) if s.donor_type == donor_type)
    for route in all_routes:
        bound = sum(1 for s in perceive(build(mol, route, donor_type))
                    if s.donor_type == donor_type)
        assert bound == free, (
            f"{name}: {free} {donor_type} free, {bound} after binding through {route}")


def test_ammonium_keeps_its_charge_on_the_nitrogen():
    """Localised charge stays on the atom; delocalised charge goes to the graph (D15)."""
    nh4 = from_smiles("[NH4+]", multiplicity=1, name="ammonium")
    assert nh4.localised_charge() == 1 and nh4.net_charge() == 1
    thq = from_smiles(THQ, multiplicity=1, name="thq")
    assert thq.localised_charge() == 0
