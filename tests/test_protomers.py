"""Activation states are an enumerated branch, and identity decides how many there are."""
from __future__ import annotations

import pytest

from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.model import chelate_pockets, find_pocket, perceive
from mofsbu.sites.protomers import enumerate_protomers, labile_sites

THQ = "OC1=C(O)C(=O)C(O)=C(O)C1=O"


@pytest.fixture(scope="module")
def protomers():
    return enumerate_protomers(mol_from_smiles(THQ))


def test_thq_has_four_acidic_hydroxyls():
    assert len(labile_sites(mol_from_smiles(THQ), "enol_O")) == 4


def test_sixteen_subsets_collapse_to_seven_protomers(protomers):
    """Nothing tells it the molecule's symmetry; the L1 hash discovers it.

    Four hydroxyls give 16 subsets.  Symmetry-equivalent selections land on one node, so
    what survives is 7 distinct ligands.
    """
    assert sum(len(p.selections) for p in protomers) == 16
    assert len(protomers) == 7


def test_every_charge_state_is_present(protomers):
    assert sorted({p.charge for p in protomers}) == [-4, -3, -2, -1, 0]


def test_the_dianion_splits_into_three_configurations(protomers):
    """Three distinct arrangements, discovered by Aut(parent) rather than named.

    The descriptors are generated orbit representatives; nothing asserts "ortho" here,
    because a taxonomy would have to be extended for every scaffold that follows.
    """
    dianions = [p for p in protomers if p.charge == -2]
    assert len(dianions) == 3
    assert len({p.configuration for p in dianions}) == 3
    assert all(len(p.selections) == 2 for p in dianions)
    assert len({p.separation for p in dianions}) == 3      # separations distinguish them too


def test_equivalence_is_decided_by_the_automorphism_group(protomers):
    """The collapse agrees with L1 — two independent routes to the same answer.

    Selections in one orbit of Aut(parent) must land on one L1 hash, and vice versa.  If
    these ever disagree, one of the two is wrong and the descriptor is not trustworthy.
    """
    by_config: dict[str, set[str]] = {}
    for p in protomers:
        by_config.setdefault(f"{p.n_deprotonated}:{p.configuration}", set()).add(p.l1)
    assert all(len(hashes) == 1 for hashes in by_config.values())
    assert len(by_config) == len(protomers)


def test_descriptors_generalise_to_other_scaffolds():
    """A taxonomy would have nothing to say about these; the orbit key does."""
    from mofsbu.graph.from_mol import mol_from_smiles as smiles

    # BTC's three carboxylates are mutually equivalent, so every doubly-deprotonated
    # selection is ONE configuration — a ring-distance table would have called them meta.
    btc = enumerate_protomers(smiles("OC(=O)c1cc(C(O)=O)cc(C(O)=O)c1"))
    dianions = [p for p in btc if p.charge == -2]
    assert len(dianions) == 1 and len(dianions[0].selections) == 3

    # A four-membered ring has no ortho/meta/para at all; it still gets descriptors.
    four = enumerate_protomers(smiles("OC1=C(O)C(O)=C1O"))
    assert len({p.configuration for p in four if p.charge == -2}) == 2


def test_singly_and_triply_deprotonated_are_each_one_structure(protomers):
    """All four single deprotonations are the same molecule; so are all four triples."""
    for charge in (-1, -3):
        matches = [p for p in protomers if p.charge == charge]
        assert len(matches) == 1 and len(matches[0].selections) == 4


def test_exactly_one_dianion_configuration_can_chelate(protomers):
    """Only one of the three arrangements puts its two anionic oxygens on adjacent carbons.

    Asserted structurally: the chelating one is identified by HAVING a dioxolene pocket,
    and independently by having the smallest separation profile.  Neither test knows the
    word "ortho", which is why this keeps working on scaffolds nobody has classified.
    """
    dianions = [p for p in protomers if p.charge == -2]
    chelating = [p for p in dianions
                 if find_pocket(embed_molecule(p.mol, seed=7),
                                perceive(embed_molecule(p.mol, seed=7)),
                                ring_size=5, n_anionic=2) is not None]
    assert len(chelating) == 1
    assert chelating[0].separation == min(p.separation for p in dianions)


def test_neutral_thq_has_no_bindable_pocket():
    """Its pockets are all still protonated: activation is a prerequisite, not a label."""
    mol = embed_molecule(mol_from_smiles(THQ), seed=7)
    sites = perceive(mol)
    assert all(p.n_anionic == 0 for p in chelate_pockets(mol, sites))
    assert find_pocket(mol, sites, ring_size=5, n_anionic=2) is None


def test_protomer_identity_is_stable_across_selections(protomers):
    """Two equivalent selections give one L1, which is what made the collapse happen."""
    assert len({p.l1 for p in protomers}) == len(protomers)
