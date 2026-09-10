"""The build routes that reach one identity through different resonance forms.

Shared by the perception tests (`test_donor_patterns`, `test_sites`) and the registry
gate (`test_sites_state`), because they are assertions about the same structures at two
different layers and a second copy of this builder would let them drift apart.

A delocalised group's oxygens are equivalent, so binding a metal through any ONE of them
gives the same identity — D15 excludes bond order from the L1 hash exactly so that C=O
and C-O(-) do not split it.  `to_rdkit` still has to render some resonance form, and it
renders whichever one puts the anion where the metal is not.  Each route below is one of
those renderings; every assertion built on them is that the layer under test cannot tell
them apart.

Atom indices are comparable ACROSS routes without canonicalising: the ligand mol is one
object placed several ways, and `to_rdkit` numbers the metal first and then the ligands in
placement order, so only which atom the dative bond lands on changes.
"""
from __future__ import annotations

import itertools

from rdkit import Chem

from mofsbu.geometry.embed import embed_molecule
from mofsbu.geometry.placer import LigandPlacement, place_mononuclear, to_rdkit
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.model import perceive

#: (name, ligand SMILES, the group's donor type, denticity).  Every case is a group whose
#: oxygens are equivalent under D15 and inequivalent under RDKit's bond orders.
ROUTE_CASES: list[tuple[str, str, str, int]] = [
    ("acetate",           "CC(=O)[O-]",         "carboxylate_O", 1),
    ("benzoate",          "[O-]C(=O)c1ccccc1",  "carboxylate_O", 1),
    ("methanesulfonate",  "CS(=O)(=O)[O-]",     "sulfonate_O",   1),
    ("methylphosphonate", "CP(=O)([O-])[O-]",   "phosphonate_O", 1),
    # The bidentate case: the same two oxygens, handed to the placer in both orders.  A
    # chelate has no "which oxygen got the metal" asymmetry to hide behind, so if the
    # donor set still moves here it moved for some other reason.
    ("acetate-chelate",   "CC(=O)[O-]",         "carboxylate_O", 2),
]

METAL = "Zn"
#: Filled out to four with aqua co-ligands, so every case is one tetrahedral centre and
#: the only thing varying between routes is which oxygen the metal is bonded to.
COORDINATION = 4


def ligand(smiles: str) -> Chem.Mol:
    return embed_molecule(mol_from_smiles(smiles), seed=7)


def donor_idxs(mol: Chem.Mol, donor_type: str) -> tuple[int, ...]:
    """The group's oxygens, as the free ligand perceives them."""
    return tuple(s.atom_idx for s in perceive(mol) if s.donor_type == donor_type)


def routes(smiles: str, donor_type: str, denticity: int) -> tuple[Chem.Mol, list[tuple[int, ...]]]:
    """The ligand, and every way it can present `denticity` of its equivalent donors."""
    mol = ligand(smiles)
    return mol, list(itertools.permutations(donor_idxs(mol, donor_type), denticity))


def build(mol: Chem.Mol, donors: tuple[int, ...], donor_type: str) -> Chem.Mol:
    """One assembled centre, binding `mol` through exactly the donors given."""
    water = ligand("O")
    mode = BindingMode.CHELATE if len(donors) > 1 else BindingMode.MONODENTATE
    placements = [LigandPlacement(mol=mol, donor_idxs=tuple(donors),
                                  donor_types=(donor_type,) * len(donors),
                                  mode=mode, name="L")]
    placements += [LigandPlacement(mol=water, donor_idxs=(0,), donor_types=("aqua_O",),
                                   name="O")
                   for _ in range(COORDINATION - len(donors))]
    result = place_mononuclear(METAL, placements, geometry="tetrahedral")
    return to_rdkit(METAL, placements, result)


def donor_set(mol: Chem.Mol) -> set[tuple[int, str, int, bool]]:
    """What a structure is recorded as having: which atom, typed how, at what charge."""
    return {(s.atom_idx, s.donor_type, s.charge_after, s.labile) for s in perceive(mol)}
