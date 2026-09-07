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


def embed_with_report(mol: Chem.Mol, *, seed: int = 0xC0FFEE,
                      relax: bool = True) -> tuple[Chem.Mol, dict]:
    """Embed, and say what it took.

    Three things used to happen silently here, and each of them changes how much the
    resulting geometry is worth:

    * **the retry** — ETKDG failing and the random-coordinates fallback succeeding gives
      a markedly worse starting geometry, so a QC clash afterwards is as likely to be
      about the embed as about the chemistry;
    * **the force field** — MMFF where parameters exist, UFF where they do not; a
      molecule relaxed under UFF is not comparable to one relaxed under MMFF;
    * **the swallowed exception** — the relaxation raising and being ignored left an
      unrelaxed conformer that looks exactly like a relaxed one.

    The report is a plain dict so it can go into the task record and the run inspector
    without inventing a schema for it.
    """
    mol = Chem.Mol(mol)
    report: dict = {"seed": seed, "retried": False, "ff": None, "relaxed": False}
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) < 0:
        report["retried"] = True
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) < 0:
            raise EmbeddingError(f"could not embed {Chem.MolToSmiles(mol)!r}")
    if relax:
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                report["ff"] = "MMFF94"
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                report["ff"] = "UFF"
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            report["relaxed"] = True
        except Exception as exc:                            # noqa: BLE001
            # a rough geometry still beats none — but it is recorded as rough
            report["relax_error"] = f"{type(exc).__name__}: {exc}"
    return mol, report


def embed_molecule(mol: Chem.Mol, *, seed: int = 0xC0FFEE, relax: bool = True) -> Chem.Mol:
    """Return a copy carrying one 3D conformer.  See `embed_with_report` for the how."""
    return embed_with_report(mol, seed=seed, relax=relax)[0]


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
