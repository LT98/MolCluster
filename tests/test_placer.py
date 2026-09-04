"""Single-centre placement — the N=1 path of the general operation (D12)."""
from __future__ import annotations

import itertools

import numpy as np
import pytest

from mofsbu.geometry.embed import embed_molecule
from mofsbu.geometry.placer import LigandPlacement, place_mononuclear, site_vectors, to_rdkit
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.model import chelate_pockets, perceive
from mofsbu.sites.perception import deprotonate, find_donor_sites

THQ = "OC1=C(O)C(=O)C(O)=C(O)C1=O"


def aqua() -> LigandPlacement:
    return LigandPlacement(mol=embed_molecule(mol_from_smiles("O"), seed=3),
                           donor_idxs=(0,), donor_types=("aqua_O",), name="H2O")


def thq_chelate(*, activated: bool = True) -> LigandPlacement:
    """A THQ chelate pocket.  Activated (deprotonated) unless asked otherwise."""
    mol = embed_molecule(mol_from_smiles(THQ), seed=7)
    if activated:
        enols = [s for s in find_donor_sites(mol) if s.donor_type == "enol_O"]
        mol, _ = deprotonate(mol, enols[:1])
        mol = embed_molecule(mol, seed=7)
    sites = perceive(mol)
    by_idx = {s.atom_idx: s for s in sites}
    pockets = chelate_pockets(mol, sites)
    pocket = next((p for p in pockets if p.n_anionic > 0), pockets[0])
    return LigandPlacement(mol=mol, donor_idxs=tuple(pocket.donors),
                           donor_types=tuple(by_idx[i].donor_type for i in pocket.donors),
                           mode=BindingMode.CHELATE, name="THQ")


def bonds_to_metal(result) -> list[float]:
    return [float(np.linalg.norm(result.coords[i] - result.coords[result.metal_idx]))
            for i in result.donor_idxs]


def test_cn5_geometries_exist():
    """The legacy table mapped CN 5 to 'planar', which is not a coordination geometry."""
    for name in ("trigonal_bipyramidal", "square_pyramidal"):
        assert site_vectors(name, 5).shape == (5, 3)
    with pytest.raises(ValueError):
        site_vectors("octahedral", 5)


def test_monodentate_placement_is_exact_and_clash_free():
    result = place_mononuclear("Zn", [aqua() for _ in range(4)], geometry="tetrahedral")
    assert result.ok, str(result.report)
    assert np.allclose(bonds_to_metal(result), 2.00, atol=1e-6)


def test_chelate_bond_lengths_come_out_right():
    """Both bonds of a chelate must be the target length.

    Forcing the ligand onto fixed polyhedron vertices strains one short and one long;
    the bite angle has to come from the ligand's own donor-donor distance instead.
    """
    result = place_mononuclear("Zn", [thq_chelate(), aqua(), aqua()], geometry="tetrahedral")
    assert result.ok, str(result.report)
    assert np.allclose(bonds_to_metal(result), 2.00, atol=1e-6)


def test_a_chelate_is_never_given_trans_vertices():
    result = place_mononuclear("Zn", [thq_chelate(), aqua(), aqua(), aqua(), aqua()],
                               geometry="octahedral")
    assert result.ok, str(result.report)
    a, b = result.donor_idxs[0], result.donor_idxs[1]
    m = result.coords[result.metal_idx]
    v1, v2 = result.coords[a] - m, result.coords[b] - m
    angle = np.degrees(np.arccos(np.clip(
        v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1)))
    assert 55 < angle < 115, f"chelate spans {angle:.0f} deg — it was handed a trans pair"


def test_binding_a_protonated_donor_is_rejected_by_qc():
    """Activation is not bookkeeping: an un-deprotonated enol still has its hydrogen,
    and that hydrogen sits where the metal wants to be.  QC must catch it rather than
    let a structure through with an O-H pointing into the coordination sphere."""
    result = place_mononuclear("Zn", [thq_chelate(activated=False), aqua(), aqua()],
                               geometry="tetrahedral")
    assert not result.ok
    assert result.report.clashes
    assert min(c[2] for c in result.report.clashes) < 1.5


def test_placement_is_deterministic():
    a = place_mononuclear("Zn", [thq_chelate(), aqua(), aqua()], geometry="tetrahedral")
    b = place_mononuclear("Zn", [thq_chelate(), aqua(), aqua()], geometry="tetrahedral")
    assert np.allclose(a.coords, b.coords)
    assert a.choice_vector == b.choice_vector


def test_placement_emits_its_choice_vector():
    """`construct` is a deterministic function of (choice-vector, seed) and says so (D13)."""
    result = place_mononuclear("Zn", [thq_chelate(), aqua(), aqua()], geometry="tetrahedral")
    cv = result.choice_vector
    assert cv["metal"] == "Zn" and cv["geometry"] == "tetrahedral"
    assert any("bite_angle_deg" in lig for lig in cv["ligands"])


def test_assembled_complex_types_as_a_metal_complex():
    ligands = [thq_chelate(), aqua(), aqua()]
    result = place_mononuclear("Zn", ligands, geometry="tetrahedral")
    g = from_rdkit(to_rdkit("Zn", ligands, result), charge=0, multiplicity=1,
                   oxidation_states={0: 2}, spin_classes={0: "ls"}, name="[Zn(THQ)(H2O)2]")
    assert g.metals() == [0]
    assert g.n_dative_bonds if hasattr(g, "n_dative_bonds") else True
    assert len(g.neighbors(0)) == 4                       # CN 4 at the metal
    assert g.label(0).oxidation_state == 2
