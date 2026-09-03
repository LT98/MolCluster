"""
EBU visualization utilities — notebook-only, not required for the DFT pipeline.

Provides:
    EBUVisualizer.draw_3d(mol, metal_symbol, ...)   — interactive 3D via py3Dmol (mol input)
    EBUVisualizer.summary(mol, metal_symbol, label) — atom/bond count summary
    view_xyz(source, ...)                           — interactive 3D via py3Dmol (XYZ input)

Requires: pip install py3Dmol
"""

import os as _os
from collections import Counter
from copy import deepcopy
from typing import Optional

import numpy as np
from rdkit import Chem

PT = Chem.GetPeriodicTable()

_SANITIZE_NO_PROPS = (
    Chem.SanitizeFlags.SANITIZE_ALL ^
    Chem.SanitizeFlags.SANITIZE_PROPERTIES
)


def _prep_xyz(mol: Chem.Mol) -> str:
    """Strip Hs, return XYZ string. Safe for mols with DATIVE bonds or aromatic flags."""
    rw   = Chem.RemoveAllHs(deepcopy(mol), sanitize=False)
    conf = rw.GetConformer()
    lines = [str(rw.GetNumAtoms()), 'EBU']
    for atom in rw.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f'{atom.GetSymbol():<4s}  {p.x:14.8f}  {p.y:14.8f}  {p.z:14.8f}')
    return '\n'.join(lines) + '\n'


class EBUVisualizer:

    @staticmethod
    def draw_3d(
        mol          : Chem.Mol,
        metal_symbol : str,
        width        : int  = 620,
        height       : int  = 460,
        spin         : bool = False,
    ):
        """
        Interactive 3D view via py3Dmol from a Chem.Mol with a 3D conformer.
        Metal atoms → large orange sphere; donor atoms → blue sphere.
        """
        try:
            import py3Dmol
        except ImportError:
            raise ImportError('pip install py3Dmol')
        if mol.GetNumConformers() == 0:
            raise RuntimeError('No 3D conformer — run GeometryPlacer.place() first.')

        mol_noH = Chem.RemoveAllHs(deepcopy(mol), sanitize=False)
        xyz_str = _prep_xyz(mol)

        an = PT.GetAtomicNumber(metal_symbol)
        metal_s, donor_s = [], []
        for atom in mol_noH.GetAtoms():
            s = atom.GetIdx()
            if atom.GetAtomicNum() == an:
                metal_s.append(s)
            elif any(nb.GetAtomicNum() == an for nb in atom.GetNeighbors()):
                donor_s.append(s)

        v = py3Dmol.view(width=width, height=height)
        v.addModel(xyz_str, 'xyz')
        v.setStyle({'stick': {'radius': 0.15, 'colorscheme': 'grayCarbon'}})
        for s in metal_s:
            v.setStyle({'serial': s}, {'sphere': {'radius': 0.55, 'color': '0xE8671A'},
                        'stick':  {'radius': 0.15, 'colorscheme': 'grayCarbon'}})
        for s in donor_s:
            v.setStyle({'serial': s},
                       {'sphere': {'radius': 0.28, 'color': '0x2E8FE8'},
                        'stick':  {'radius': 0.15, 'colorscheme': 'grayCarbon'}})
        v.zoomTo()
        if spin:
            v.spin('y', 0.5)
        v.show()
        return v

    @staticmethod
    def summary(mol: Chem.Mol, metal_symbol: str, label: str = ''):
        """Print atom count, formula, and M-L bond count for an EBU mol."""
        comp    = Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
        formula = ''.join(f"{e}{n if n > 1 else ''}" for e, n in sorted(comp.items()))
        an      = PT.GetAtomicNumber(metal_symbol)
        n_m     = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == an)
        n_ml    = sum(1 for b in mol.GetBonds()
                      if (b.GetBeginAtom().GetAtomicNum() == an) ^
                         (b.GetEndAtom().GetAtomicNum() == an))
        print(f"\n{'─'*50}\n  EBU: {label}\n  Formula: {formula}\n  "
              f"Atoms: {mol.GetNumAtoms()}  |  {n_m}x{metal_symbol}  |  {n_ml} M-L bonds\n{'─'*50}")


def view_xyz(
    source,
    width         : int            = 540,
    height        : int            = 380,
    spin          : bool           = False,
    metal_symbols                  = None,
    style         : str            = 'stick',
):
    """
    Visualise an XYZ structure inside the notebook using py3Dmol.

    source        : XYZ string OR path to an .xyz file
    metal_symbols : str or list — rendered as large orange spheres
    style         : 'stick' | 'sphere' | 'line'

    Requires: pip install py3Dmol
    """
    try:
        import py3Dmol
    except ImportError:
        raise ImportError('pip install py3Dmol')

    if isinstance(source, str) and '\n' not in source and _os.path.isfile(source):
        with open(source) as fh:
            xyz = fh.read()
    else:
        xyz = source

    v = py3Dmol.view(width=width, height=height)
    v.addModel(xyz, 'xyz')
    v.setStyle({'stick': {'radius': 0.14, 'colorscheme': 'grayCarbon'}})

    if metal_symbols:
        syms = [metal_symbols] if isinstance(metal_symbols, str) else list(metal_symbols)
        for sym in syms:
            v.setStyle({'elem': sym},
                       {'sphere': {'radius': 0.55, 'color': '0xE8671A'},
                        'stick':  {'radius': 0.14, 'colorscheme': 'grayCarbon'}})
    v.zoomTo().spin('y', 0.5).show()
    return v
