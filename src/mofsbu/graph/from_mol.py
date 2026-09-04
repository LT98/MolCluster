"""RDKit molecules -> typed graphs.  The molecular input path.

Two things need care crossing this boundary, and both are identity decisions:

**Charge placement (D15).**  RDKit puts formal charge on specific atoms.  For an
ammonium nitrogen that is correct — the charge really is on that atom.  For a
deprotonated phenol or carboxylate it is a bookkeeping convention: the charge is
delocalised over a pi system, and pinning it to one oxygen makes that oxygen
inequivalent to its partner, which would give one molecule several hashes depending on
which resonance form the SMILES happened to be written in.  So charge on a
pi-conjugated atom is lifted to the graph level and only genuinely localised charge
stays on the atom.

**Bond typing.**  Metal-ligand bonds become DATIVE, metal-metal bonds METAL_METAL, and
everything else COVALENT.  Bond ORDER is carried through as an annotation but never
enters identity, so Kekule forms and C=O/C-O resonance do not split a structure.
"""
from __future__ import annotations

from rdkit import Chem

from mofsbu.graph._types import METALS, EdgeType, NodeLabel, TypedGraph
from mofsbu._types import AmbiguousSpecError, GraphValidationError


def _is_delocalised(atom: Chem.Atom) -> bool:
    """Is this atom's formal charge spread over a pi system rather than sitting on it?

    True for aromatic atoms and for any atom adjacent to a multiple bond — a phenolate
    oxygen, a carboxylate oxygen, an enolate on a quinone.  False for a saturated centre
    with no neighbouring pi system, such as an ammonium nitrogen.
    """
    if atom.GetIsAromatic():
        return True
    for bond in atom.GetBonds():
        if bond.GetBondTypeAsDouble() > 1.0:
            return True
        other = bond.GetOtherAtom(atom)
        if other.GetIsAromatic():
            return True
        if any(b.GetBondTypeAsDouble() > 1.0 for b in other.GetBonds()):
            return True
    return False


def from_rdkit(
    mol: Chem.Mol,
    *,
    charge: int | None = None,
    multiplicity: int | None = None,
    name: str = "",
    oxidation_states: dict[int, int] | None = None,
    spin_classes: dict[int, str] | None = None,
) -> TypedGraph:
    """Convert a hydrogen-explicit RDKit mol into a typed graph.

    `charge` defaults to the mol's total formal charge.  `multiplicity` is never
    guessed (ground rule 5) — for a metal-free closed-shell molecule pass 1 explicitly.
    """
    if any(a.GetAtomicNum() > 1 and a.GetTotalNumHs() > 0 for a in mol.GetAtoms()):
        raise GraphValidationError(
            "hydrogens must be explicit before conversion; call Chem.AddHs first"
        )
    total = Chem.GetFormalCharge(mol) if charge is None else charge
    g = TypedGraph(charge=total, multiplicity=multiplicity, name=name)
    oxidation_states = oxidation_states or {}
    spin_classes = spin_classes or {}

    remap: dict[int, int] = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        is_metal = atom.GetSymbol() in METALS
        local = 0 if _is_delocalised(atom) and not is_metal else atom.GetFormalCharge()
        remap[idx] = g.add_atom(
            atom.GetSymbol(),
            formal_charge=0 if is_metal else local,
            oxidation_state=oxidation_states.get(idx) if is_metal else None,
            spin_class=spin_classes.get(idx) if is_metal else None,
        )

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtom(), bond.GetEndAtom()
        mi, mj = i.GetSymbol() in METALS, j.GetSymbol() in METALS
        if mi and mj:
            etype = EdgeType.METAL_METAL
        elif mi or mj:
            etype = EdgeType.DATIVE
        else:
            etype = EdgeType.COVALENT
        g.add_bond(remap[i.GetIdx()], remap[j.GetIdx()], etype,
                   order=bond.GetBondTypeAsDouble() or 1.0)
    return g


def mol_from_smiles(smiles: str, *, embed: bool = False, seed: int = 0xC0FFEE) -> Chem.Mol:
    """Parse SMILES to a hydrogen-explicit mol, optionally with 3D coordinates."""
    raw = Chem.MolFromSmiles(smiles)
    if raw is None:
        raise GraphValidationError(f"invalid SMILES: {smiles!r}")
    mol = Chem.AddHs(raw)
    if embed:
        from mofsbu.geometry.embed import embed_molecule

        mol = embed_molecule(mol, seed=seed)
    return mol


def from_smiles(
    smiles: str,
    *,
    charge: int | None = None,
    multiplicity: int | None = None,
    name: str = "",
) -> TypedGraph:
    if multiplicity is None:
        raise AmbiguousSpecError(
            f"multiplicity not given for {smiles!r}; state it explicitly rather than "
            "assuming a closed shell (ground rule 5)"
        )
    return from_rdkit(mol_from_smiles(smiles), charge=charge,
                      multiplicity=multiplicity, name=name or smiles)
