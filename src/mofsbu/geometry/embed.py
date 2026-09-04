"""Real 3D coordinates for molecules, via RDKit distance geometry.

This is genuine molecular geometry — ETKDG followed by an MMFF/UFF relaxation — unlike
`geometry.layout`, which only draws a graph.  It applies to organic fragments; metal
complexes are assembled from these by the placer.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from mofsbu._types import MofsbuError


class EmbeddingError(MofsbuError):
    pass


def embed_molecule(mol: Chem.Mol, *, seed: int = 0xC0FFEE, relax: bool = True) -> Chem.Mol:
    """Return a copy carrying one 3D conformer."""
    mol = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) < 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) < 0:
            raise EmbeddingError(f"could not embed {Chem.MolToSmiles(mol)!r}")
    if relax:
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:                                   # noqa: BLE001
            pass                                            # a rough geometry still beats none
    return mol


def coordinates(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])


def set_coordinates(mol: Chem.Mol, coords: np.ndarray) -> Chem.Mol:
    from rdkit.Geometry import Point3D

    mol = Chem.Mol(mol)
    conf = mol.GetConformer()
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    return mol


def to_xyz(mol: Chem.Mol, comment: str = "") -> str:
    coords = coordinates(mol)
    lines = [str(mol.GetNumAtoms()), comment]
    for atom, (x, y, z) in zip(mol.GetAtoms(), coords):
        lines.append(f"{atom.GetSymbol():<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"
