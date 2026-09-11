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
#: Di-tert-butyl ketone.  One small donor carrying two large rotors, which is what makes
#: it the case that found the unconstrained rotations: the tert-butyls reach far enough
#: to hit a neighbouring vertex, and they flank the carbonyl closely enough to sit on the
#: metal itself.
DTBK = "C(C)(C)(C)C(C(C)(C)(C))=O"


def aqua() -> LigandPlacement:
    return LigandPlacement(mol=embed_molecule(mol_from_smiles("O"), seed=3),
                           donor_idxs=(0,), donor_types=("aqua_O",), name="H2O")


def monodentate(smiles: str, donor_type: str, name: str, seed: int = 7,
                **kw) -> LigandPlacement:
    mol = embed_molecule(mol_from_smiles(smiles), seed=seed)
    idx = next(s.idx for s in find_donor_sites(mol) if s.donor_type == donor_type)
    return LigandPlacement(mol=mol, donor_idxs=(idx,), donor_types=(donor_type,),
                           mode=BindingMode.MONODENTATE, name=name, **kw)


def dtbk(name: str = "dtBK", **kw) -> LigandPlacement:
    return monodentate(DTBK, "carbonyl_O", name, **kw)


def chloride(name: str = "Cl") -> LigandPlacement:
    return monodentate("[Cl-]", "halide_Cl", name, seed=3)


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


def test_a_bulky_ligand_does_not_sit_on_its_neighbour():
    """The case that named the bug: dtBK and a chloride on adjacent tetrahedral vertices.

    Pinning the donor leaves the ligand free to spin about the M-L axis, and left
    arbitrary that spin put a tert-butyl group on top of the chloride —
    `C[dtBK]...Cl 1.43 A` against a 2.14 A limit, over by 0.71.
    """
    result = place_mononuclear("Mg", [dtbk(), chloride(), aqua(), aqua()],
                               geometry="tetrahedral")
    assert result.ok, str(result.report)


def test_two_bulky_ligands_do_not_rotate_into_each_other():
    """Legacy's own symptom, on separate vertices: 'bipy ligands rotating into each
    other around Fe3+ despite using separate sites'.  Two dtBK on one centre reproduced
    it as `C[dtBK]...C[dtBK] 0.79 A` — two ligands sharing the same space entirely."""
    result = place_mononuclear("Mg", [dtbk("dtBK_a"), dtbk("dtBK_b"), aqua(), aqua()],
                               geometry="tetrahedral")
    assert result.ok, str(result.report)


def test_the_metal_is_swung_out_of_a_hindered_donors_plane():
    """A spin about the M-L axis cannot move the metal, because the metal is ON that
    axis.  dtBK puts a methyl hydrogen 1.41 A from the metal (limit 1.46) at every
    azimuth and in both of the frame's in-plane lone pairs, so clearing it needs the
    other rotation: swinging the metal off the donor's sp2 plane, which leaves the
    M-O-C angle exactly where the frame set it."""
    result = place_mononuclear("Mg", [dtbk(), chloride(), aqua(), aqua()],
                               geometry="tetrahedral")
    assert result.choice_vector["ligands"][0]["oop_step"] != 0
    m = result.coords[result.metal_idx]
    closest = min(float(np.linalg.norm(result.coords[i] - m))
                  for i, s in enumerate(result.symbols) if s == "H")
    assert closest > 1.465, f"a hydrogen sits {closest:.3f} A from the metal"


def test_replaying_the_emitted_indices_reproduces_the_geometry():
    """D13: construct is a deterministic function of its choice-vector.  The search has
    to leave the vector sufficient — handing the emitted indices back must give the same
    coordinates without searching again, or the stored vector is not the whole recipe."""
    ligands = [dtbk(), chloride(), aqua(), aqua()]
    first = place_mononuclear("Mg", ligands, geometry="tetrahedral")
    chosen = first.choice_vector["ligands"][0]

    replay = [dtbk(azimuth_step=chosen["azimuth_step"], oop_step=chosen["oop_step"]),
              chloride(), aqua(), aqua()]
    second = place_mononuclear("Mg", replay, geometry="tetrahedral")
    n = first.atom_offsets[1] - first.atom_offsets[0]
    sl = slice(first.atom_offsets[0], first.atom_offsets[0] + n)
    assert np.allclose(first.coords[sl], second.coords[sl])
    assert second.choice_vector["ligands"][0]["azimuth_step"] == chosen["azimuth_step"]


def test_the_search_never_merges_two_torsion_wells():
    """The wells are the stored conformer coordinate (D13), and D11 leans on two
    structures that differ only in their choice-vector being different structures.  A
    search free to sweep the whole circle would erase that — both wells would find the
    same optimum and the index would stop meaning anything.  Offsets are therefore capped
    at half a well spacing, so a well can never be carried onto its neighbour."""
    a = place_mononuclear("Mg", [dtbk(), chloride(), aqua(), aqua()],
                          geometry="tetrahedral")
    flipped = dtbk()
    flipped.torsion_well = 1
    b = place_mononuclear("Mg", [flipped, chloride(), aqua(), aqua()],
                          geometry="tetrahedral")

    n = a.atom_offsets[1] - a.atom_offsets[0]
    sl = slice(a.atom_offsets[0], a.atom_offsets[0] + n)
    assert not np.allclose(a.coords[sl], b.coords[sl])
    apart = abs(a.choice_vector["ligands"][0]["azimuth_deg"]
                - b.choice_vector["ligands"][0]["azimuth_deg"])
    assert 90.0 < apart < 270.0, f"wells ended up {apart:.0f} deg apart"


def test_placement_is_deterministic():
    a = place_mononuclear("Zn", [thq_chelate(), aqua(), aqua()], geometry="tetrahedral")
    b = place_mononuclear("Zn", [thq_chelate(), aqua(), aqua()], geometry="tetrahedral")
    assert np.allclose(a.coords, b.coords)
    assert a.choice_vector == b.choice_vector


def test_the_orientation_search_is_deterministic():
    """The reinstated search is the part most at risk of being stochastic, which is what
    D13 replaced legacy's version to avoid.  It samples a fixed grid in a fixed order and
    consults no seed, so repeating it is not merely likely to agree — it must."""
    def build():
        return [dtbk(), chloride(), aqua(), aqua()]

    a = place_mononuclear("Mg", build(), geometry="tetrahedral")
    b = place_mononuclear("Mg", build(), geometry="tetrahedral")
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


# ── coordinative unsaturation: vacant vertices are a RESULT, not a smaller centre ──

def test_an_unsaturated_centre_keeps_the_polyhedron_it_was_asked_for():
    """Two donors on a tetrahedron is not a linear two-coordinate centre.

    The placer used to be handed `cn = n_donors`, so "CN 4, two ligands" could only be
    built by rebuilding it as CN 2 — which is a different molecule with different bond
    angles, and the one thing it is NOT is the species the next assembly step attaches
    to.  `cn` is now independent of the donor count and the surplus vertices come back.
    """
    result = place_mononuclear("Zn", [aqua(), aqua()], geometry="tetrahedral", cn=4)
    assert len(result.donor_idxs) == 2
    assert len(result.vacancies) == 2
    assert result.cn == 4

    # the vacant directions really are the unused tetrahedral vertices
    import numpy as np

    verts = site_vectors("tetrahedral", 4, 1.0)
    for vac in result.vacancies:
        assert min(float(np.linalg.norm(np.array(vac) - v)) for v in verts) < 1e-9
    # and they are tetrahedral to each other, not collinear
    a, b = (np.array(v) for v in result.vacancies)
    angle = np.degrees(np.arccos(np.clip(a @ b, -1, 1)))
    assert 100.0 < angle < 120.0, "vacancies should sit at the tetrahedral angle"


def test_a_saturated_call_is_unchanged_and_has_no_vacancies():
    result = place_mononuclear("Zn", [aqua() for _ in range(4)], geometry="tetrahedral")
    assert result.vacancies == () and result.cn == 4


def test_more_donors_than_vertices_is_refused():
    with pytest.raises(ValueError, match="cannot fit"):
        place_mononuclear("Zn", [aqua() for _ in range(4)], geometry="linear", cn=2)


def test_a_chelate_still_gets_cis_vertices_when_the_rest_are_vacant():
    """The bite-angle assignment must not be disturbed by leaving vertices empty."""
    result = place_mononuclear("Zn", [thq_chelate()], geometry="octahedral", cn=6)
    assert len(result.donor_idxs) == 2 and len(result.vacancies) == 4
    import numpy as np

    d1, d2 = (result.coords[i] - result.coords[0] for i in result.donor_idxs)
    cos = float(d1 @ d2 / (np.linalg.norm(d1) * np.linalg.norm(d2)))
    assert np.degrees(np.arccos(np.clip(cos, -1, 1))) < 120.0, "chelate went trans"
