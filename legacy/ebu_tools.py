import warnings, numpy as np
import os as _os
import copy
from itertools import product as _iproduct

from copy import deepcopy
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, RWMol

PT = Chem.GetPeriodicTable()

# SMARTS: (group_smarts, labile_H_smarts, donor_atom_pos, attachment_atom_pos)
#
# donor_atom_pos may be a single int or a list[int].
#   Single int  → one donor per functional group (existing behaviour).
#   List[int]   → multiple donors per group; index 0 is the PRIMARY donor.
#                 Donors are ordered in the flat list as: all primaries first,
#                 then all secondaries — so with_sites(n ≤ n_groups) is backward
#                 compatible (picks one donor per group, as before).
#
# carboxylate example  [CX3](=O)[OH]:
#   pos 0 = carboxylate C, pos 1 = keto O (=O), pos 2 = hydroxyl O (-OH)
#   Primary donor (pos 2): the deprotonated -O — same as the old single-donor behaviour.
#   Secondary donor (pos 1): the keto =O — equally valid coordination site.
FG_DEF: Dict[str, tuple] = {
    #                   group SMARTS         H-removal  donor      attachment
    'carboxylate'  : ('[CX3](=O)[OH]',    '[OH]',     [2, 1],    0),
    'carboxylate_C': ('[CX3](=O)[OH]',    '[OH]',     0,         0),
    'hydroxy'      : ('[OH]',             '[OH]',     0,         0),
    'pyridyl'      : ('n',                 None,       0,         0),
    'imidazolate'  : ('[nH]',              '[nH]',     0,         0),
    'amine'        : ('[NX3;H2]',          '[NX3;H2]', 0,         0),
    'catecholate'  : ('[c][OH]',           '[OH]',     1,         0),
    'thiol'        : ('[SX2H]',            '[SX2H]',   0,         0),
}

def _site_vectors(geometry: str, n: int, d: float = 2.05) -> np.ndarray:
    g = geometry.lower()
    if g == 'planar':
        angle = np.radians([(360 / n) * i for i in range(n)])
        return d * np.column_stack([np.cos(angle), np.sin(angle), np.zeros(n)])
    elif g == 'tetrahedral':
        return d * np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float) / np.sqrt(3)
    elif g == 'octahedral':
        return d * np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], float)
    else:
        raise ValueError(f'Unknown geometry {geometry!r}. Options: planar, tetrahedral, octahedral')


class LigandBuilder:
    """
    Load a ligand from SMILES, find coordination sites, remove labile H.

    Parameters
    ----------
    smiles    : ligand SMILES
    fg_type   : coordination group type — key in FG_DEF
    max_sites : activate at most this many sites (None = all)
    name      : display label
    """

    def __init__(
        self,
        smiles    : str,
        fg_type   : str           = 'carboxylate',
        max_sites : Optional[int] = None,
        name      : str           = 'ligand',
    ):
        if fg_type not in FG_DEF:
            raise ValueError(f'Unknown fg_type {fg_type!r}. Options: {list(FG_DEF)}')
        self.name      = name
        self.fg_type   = fg_type
        self.max_sites = max_sites

        raw = Chem.MolFromSmiles(smiles)
        if raw is None:
            raise ValueError(f'Invalid SMILES: {smiles!r}')

        (self.mol,
         self.donor_indices, self.attachment_indices,
         self.unique_donor_indices, self.unique_attachment_indices,
         self._canon_ranks) = self._activate(raw)

        self.n_sites  = len(self.donor_indices)
        self.n_unique = len(self.unique_donor_indices)
        print(f'[LigandBuilder] {name!r} | {fg_type} | '
              f'{self.n_sites} site(s) ({self.n_unique} unique) | '
              f'donors: {self.donor_indices} | '
              f'unique donors: {self.unique_donor_indices}')

    def _activate(self, raw):
        smarts, h_smarts, donor_pos, attach_pos = FG_DEF[self.fg_type]
        donor_pos_list = [donor_pos] if isinstance(donor_pos, int) else list(donor_pos)
        n_layers = len(donor_pos_list)

        mol = Chem.RWMol(Chem.AddHs(deepcopy(raw)))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)

        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        if self.max_sites is not None:
            matches = matches[:self.max_sites]

        # Per-group storage: group_donors[g][layer], group_attach[g]
        group_donors = []
        group_attach = []
        h_remove     = []

        h_q = Chem.MolFromSmarts(h_smarts) if h_smarts is not None else None

        for match in matches:
            group_donors.append([match[dp] for dp in donor_pos_list])
            group_attach.append(match[attach_pos])
            if h_q is not None:
                match_set = set(match)
                for hm in mol.GetSubstructMatches(h_q):
                    h_bearer = hm[-1]
                    if h_bearer in match_set:
                        for nb in mol.GetAtomWithIdx(h_bearer).GetNeighbors():
                            if nb.GetAtomicNum() == 1:
                                h_remove.append(nb.GetIdx())
                                break
                        break

        for hidx in sorted(set(h_remove), reverse=True):
            mol.RemoveAtom(hidx)
            group_donors = [
                [d - (1 if d > hidx else 0) for d in dlist]
                for dlist in group_donors
            ]
            group_attach = [a - (1 if a > hidx else 0) for a in group_attach]

        try:
            Chem.SanitizeMol(mol,
                Chem.SanitizeFlags.SANITIZE_ALL ^
                Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass

        final_mol = mol.GetMol()

        # Flat lists: all primaries first, then layer 1, layer 2, …
        donor_idxs  = []
        attach_idxs = []
        for layer in range(n_layers):
            for g in range(len(group_donors)):
                donor_idxs.append(group_donors[g][layer])
                attach_idxs.append(group_attach[g])

        # Deduplicate symmetry-equivalent sites via canonical atom ranks
        canon_ranks = list(Chem.CanonicalRankAtoms(final_mol, breakTies=False))
        seen_keys          = set()
        unique_donor_idxs  = []
        unique_attach_idxs = []
        for d, a in zip(donor_idxs, attach_idxs):
            key = (canon_ranks[d], canon_ranks[a])
            if key not in seen_keys:
                seen_keys.add(key)
                unique_donor_idxs.append(d)
                unique_attach_idxs.append(a)

        return final_mol, donor_idxs, attach_idxs, unique_donor_idxs, unique_attach_idxs, canon_ranks

    def get_outward_vectors(self) -> List[np.ndarray]:
        """
        For each active donor, return a unit vector pointing from the donor
        toward its heavy-atom bonding environment (i.e. anti-lone-pair direction).

        GeometryPlacer._rot aligns this vector with the site target (M→donor),
        which places the metal on the lone-pair side of the donor — the correct
        approach geometry for every supported functional group:

          sp-O / sp-S  (carboxylate, hydroxy, catecholate, thiol)
              one C neighbour → outward = C − donor  (along C−X bond)
          aromatic-N  (pyridyl, imidazolate)
              two C neighbours → outward = mean(C₁,C₂) − N  (bisects C−N−C inward)
          amine-N
              one or more C neighbours → outward = mean(Cₙ) − N
        """
        conf = self.mol.GetConformer()
        vectors = []
        for donor_idx in self.donor_indices:
            d_pos = np.array(conf.GetAtomPosition(donor_idx))
            atom  = self.mol.GetAtomWithIdx(donor_idx)
            nb_heavy = [
                np.array(conf.GetAtomPosition(nb.GetIdx()))
                for nb in atom.GetNeighbors()
                if nb.GetAtomicNum() > 1
            ]
            if nb_heavy:
                outward = np.mean(nb_heavy, axis=0) - d_pos
            else:
                # Isolated donor — fall back to centroid direction
                heavy_pos = np.array([
                    list(conf.GetAtomPosition(a.GetIdx()))
                    for a in self.mol.GetAtoms() if a.GetAtomicNum() > 1
                ])
                outward = heavy_pos.mean(axis=0) - d_pos
            norm = np.linalg.norm(outward)
            vectors.append(outward / max(norm, 1e-9))
        return vectors

    def with_sites(self, n: int) -> 'LigandBuilder':
        """
        Return a copy that uses only the first n donor sites.

        Donors are ordered primaries-first, so with_sites(n_groups) selects
        exactly one primary donor per group — backward-compatible behaviour.
        """
        c = copy.copy(self)
        c.donor_indices      = self.donor_indices[:n]
        c.attachment_indices = self.attachment_indices[:n]
        c.n_sites            = n
        seen_keys = set()
        c.unique_donor_indices      = []
        c.unique_attachment_indices = []
        for d, a in zip(c.donor_indices, c.attachment_indices):
            key = (self._canon_ranks[d], self._canon_ranks[a])
            if key not in seen_keys:
                seen_keys.add(key)
                c.unique_donor_indices.append(d)
                c.unique_attachment_indices.append(a)
        c.n_unique = len(c.unique_donor_indices)
        base   = self.name.split('[')[0]
        c.name = f'{base}[s{n}]'
        return c

    def with_donor(self, slot: int) -> 'LigandBuilder':
        """
        Return a copy that uses only the unique donor at *slot* index.

        Indexes into unique_donor_indices so symmetry-equivalent duplicates
        are skipped; slot 0 = primary unique donor, slot 1 = secondary, etc.
        """
        if slot >= len(self.unique_donor_indices):
            raise IndexError(
                f'slot {slot} out of range '
                f'(0..{len(self.unique_donor_indices) - 1})'
            )
        c = copy.copy(self)
        d = self.unique_donor_indices[slot]
        a = self.unique_attachment_indices[slot]
        c.donor_indices             = [d]
        c.attachment_indices        = [a]
        c.unique_donor_indices      = [d]
        c.unique_attachment_indices = [a]
        c.n_sites  = 1
        c.n_unique = 1
        base   = self.name.split('[')[0]
        c.name = f'{base}[d{slot}]'
        return c


METAL_CHARGE: Dict[str, int] = {
    'Zn':2,'Cu':2,'Fe':3,'Co':2,'Ni':2,'Mn':2,'Cd':2,
    'Zr':4,'Hf':4,'Al':3,'Cr':3,'In':3,'Eu':3,'Y':3,'Ti':4,'V':3,
}

_SANITIZE_NO_PROPS = (
    Chem.SanitizeFlags.SANITIZE_ALL ^
    Chem.SanitizeFlags.SANITIZE_PROPERTIES
)


class EBUBuilder:
    """
    Assemble EBU: metal(s) + full ligand molecules joined by DATIVE bonds.

    Parameters
    ----------
    metal_symbol    : element symbol
    oxidation_state : formal charge; default from METAL_CHARGE table.
    """

    def __init__(self, metal_symbol: str, oxidation_state: Optional[int] = None):
        self.symbol     = metal_symbol
        self.charge     = (oxidation_state if oxidation_state is not None
                           else METAL_CHARGE.get(metal_symbol, 2))
        self.atomic_num = PT.GetAtomicNumber(metal_symbol)

    def _new_metal(self) -> Chem.Atom:
        a = Chem.Atom(self.atomic_num)
        a.SetFormalCharge(self.charge)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        return a

    def build(
        self,
        ligands          : List[LigandBuilder],
        sites_per_ligand : Optional[List[int]] = None,
        n_ligands        : Optional[int]        = None,
    ) -> Chem.Mol:
        """
        Mononuclear EBU.

        Each ligand is appended exactly once; multiple dative bonds for
        multi-site ligands all point to that single metal atom.
        atom count = 1 metal + sum(atoms per ligand).
        """
        if n_ligands is not None:
            ligands = ligands[:n_ligands]
        if sites_per_ligand is None:
            sites_per_ligand = [l.n_sites for l in ligands]

        combo     = RWMol()
        metal_idx = combo.AddAtom(self._new_metal())

        for lig, n_use in zip(ligands, sites_per_ligand):
            offset = combo.GetNumAtoms()
            combo  = RWMol(Chem.CombineMols(combo, lig.mol))
            for slot in range(min(n_use, lig.n_sites)):
                combo.AddBond(metal_idx,
                              offset + lig.donor_indices[slot],
                              Chem.BondType.DATIVE)

        try:
            Chem.SanitizeMol(combo, _SANITIZE_NO_PROPS)
        except Exception:
            pass
        combo.GetAtomWithIdx(metal_idx).SetFormalCharge(self.charge)

        mol        = combo.GetMol()
        n_ml       = sum(min(n, l.n_sites) for l, n in zip(ligands, sites_per_ligand))
        n_expected = 1 + sum(l.mol.GetNumAtoms() for l in ligands)
        status     = '✓' if mol.GetNumAtoms() == n_expected else f'✗ expected {n_expected}'
        print(f'[EBUBuilder] {self.symbol}({self.charge}+) '
              f'+ {len(ligands)} ligand(s) | {n_ml} M-L bonds | '
              f'{mol.GetNumAtoms()} atoms {status}')
        return mol


class GeometryPlacer:

    @staticmethod
    def _rot(v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
        """Rodrigues rotation matrix: maps unit vector v_from → v_to."""
        f = v_from / np.linalg.norm(v_from)
        t = v_to   / np.linalg.norm(v_to)
        c = np.cross(f, t)
        d = float(np.dot(f, t))
        s = np.linalg.norm(c)
        if s < 1e-8:
            if d > 0:
                return np.eye(3)
            perp = np.array([1, 0, 0]) if abs(f[0]) < 0.9 else np.array([0, 1, 0])
            ax = np.cross(f, perp); ax /= np.linalg.norm(ax)
            return 2 * np.outer(ax, ax) - np.eye(3)
        K = np.array([[0, -c[2], c[1]], [c[2], 0, -c[0]], [-c[1], c[0], 0]])
        return np.eye(3) + K + K @ K * (1 - d) / s**2

    @staticmethod
    def place(
        ebu_mol      : Chem.Mol,
        metal_symbol : str,
        geometry     : str,
        ligands      : List[LigandBuilder],
        d_ml         : float = 2.05,
    ) -> Chem.Mol:
        """
        Assign 3D coordinates to a built EBU.

        Each ligand occupies one arm of the coordination geometry.
        The ligand is rotated so its primary outward vector aligns with the
        arm's site vector, then translated so the donor lands at d_ml from
        the metal.
        """
        n_sites = len(ligands)
        sites   = _site_vectors(geometry, n_sites, d=d_ml)
        an      = PT.GetAtomicNumber(metal_symbol)
        m_idxs  = [a.GetIdx() for a in ebu_mol.GetAtoms() if a.GetAtomicNum() == an]

        conf = Chem.Conformer(ebu_mol.GetNumAtoms())
        conf.SetAtomPosition(m_idxs[0], (0., 0., 0.))
        if len(m_idxs) > 1:
            conf.SetAtomPosition(m_idxs[1], (0., 0., 2.6))

        offset = len(m_idxs)
        for lig_i, lig in enumerate(ligands):
            lig_conf    = lig.mol.GetConformer()
            n_lig_atoms = lig.mol.GetNumAtoms()
            target      = sites[lig_i]

            outward = lig.get_outward_vectors()[0]
            R       = GeometryPlacer._rot(outward, target)

            donor_local = np.array(lig_conf.GetAtomPosition(lig.donor_indices[0]))
            translation = target - R @ donor_local

            for j in range(n_lig_atoms):
                p_local = np.array(lig_conf.GetAtomPosition(j))
                conf.SetAtomPosition(offset + j, (R @ p_local + translation).tolist())

            offset += n_lig_atoms

        rw = RWMol(ebu_mol)
        rw.RemoveAllConformers()
        rw.AddConformer(conf, assignId=True)
        return rw.GetMol()

    @staticmethod
    def to_ase(mol_3d: Chem.Mol):
        """Convert a placed EBU mol to an ASE Atoms object (vacuum box, no PBC)."""
        import ase
        conf  = mol_3d.GetConformer()
        syms  = [a.GetSymbol() for a in mol_3d.GetAtoms()]
        pos   = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_3d.GetNumAtoms())])
        atoms = ase.Atoms(symbols=syms, positions=pos)
        atoms.center(vacuum=12.)
        atoms.pbc = False
        return atoms


class EBUExporter:
    """Export a placed EBU mol (with 3D conformer) to structure formats."""

    @staticmethod
    def _require_conformer(mol: Chem.Mol):
        if mol.GetNumConformers() == 0:
            raise RuntimeError('No 3D conformer — run GeometryPlacer.place() first.')

    @staticmethod
    def to_xyz(mol_3d: Chem.Mol, path: Optional[str] = None, comment: str = 'EBU') -> str:
        """
        Write XYZ format. Returns the XYZ string; also writes to *path* if given.
        No external dependencies beyond RDKit.
        """
        EBUExporter._require_conformer(mol_3d)
        conf  = mol_3d.GetConformer()
        lines = [str(mol_3d.GetNumAtoms()), comment]
        for atom in mol_3d.GetAtoms():
            p = conf.GetAtomPosition(atom.GetIdx())
            lines.append(f'{atom.GetSymbol():<4s}  {p.x:14.8f}  {p.y:14.8f}  {p.z:14.8f}')
        xyz = '\n'.join(lines) + '\n'
        if path is not None:
            with open(path, 'w') as fh:
                fh.write(xyz)
            print(f'XYZ → {path}')
        return xyz


def enumerate_sbus(
    builder        : 'EBUBuilder',
    ligand_options : List[List['LigandBuilder']],
    geometry       : str,
    out_dir        : str,
    label_prefix   : str,
    constraint     : Optional[callable]  = None,
    metal_symbol   : Optional[str]       = None,
    d_ml           : float               = 2.05,
    display_fn     : Optional[callable]  = None,
) -> Dict[str, Tuple['Chem.Mol', 'Chem.Mol']]:
    """
    Enumerate SBU configurations by Cartesian product over per-ligand variants.

    Parameters
    ----------
    builder        : EBUBuilder for the target metal.
    ligand_options : Per-position lists of LigandBuilder variants to try.
                     e.g. [[atf.with_sites(1), atf.with_sites(2)],
                            [edta.with_sites(n) for n in range(1, 5)]]
    geometry       : Coordination geometry ('planar', 'tetrahedral', 'octahedral').
    out_dir        : Directory for XYZ output files (created if absent).
    label_prefix   : String prepended to every output label and filename.
    constraint     : Optional callable(combo: list[LigandBuilder]) -> bool.
                     Return False to skip a combination.
    metal_symbol   : Override metal symbol for GeometryPlacer (default: builder.symbol).
    d_ml           : Metal-donor bond length in Å.
    display_fn     : Optional callable(mol, mol3d, label) invoked after each EBU is built.

    Returns
    -------
    dict mapping label -> (mol_graph, mol_3d)
    """
    _os.makedirs(out_dir, exist_ok=True)
    sym     = metal_symbol or builder.symbol
    results: Dict[str, Tuple] = {}

    for combo in _iproduct(*ligand_options):
        combo = list(combo)
        if constraint is not None and not constraint(combo):
            continue

        label    = label_prefix + '_' + '_'.join(lig.name for lig in combo)
        xyz_path = _os.path.join(out_dir, f'{label}.xyz')

        sites = [lig.n_sites for lig in combo]
        mol   = builder.build(combo, sites_per_ligand=sites)
        mol3d = GeometryPlacer.place(mol, sym, geometry, combo, d_ml=d_ml)
        EBUExporter.to_xyz(mol3d, path=xyz_path, comment=label)

        if display_fn is not None:
            display_fn(mol, mol3d, label)

        results[label] = (mol, mol3d)

    print(f'{len(results)} EBU(s) written → {out_dir}')
    return results


# ── Extended builders ────────────────────────────────────────────────────────
# Not part of the core mononuclear pipeline; kept here for future MOF nodes.

class PaddlewheelBuilder(EBUBuilder):
    """Binuclear M₂ paddlewheel builder."""

    def build_paddlewheel(
        self,
        bridging_ligand : LigandBuilder,
        n_bridging      : int                           = 4,
        axial_ligands   : Optional[List[LigandBuilder]] = None,
    ) -> Chem.Mol:
        """
        Binuclear M₂ paddlewheel.
        Bridging ligand needs >= 2 donor sites: site[0]→M1, site[1]→M2.
        axial_ligands: exactly 2 LigandBuilders (one per metal) or None.
        """
        if bridging_ligand.n_sites < 2:
            raise ValueError(
                f'Bridging ligand needs >= 2 donor sites (got {bridging_ligand.n_sites}).')

        combo = RWMol()
        m1 = combo.AddAtom(self._new_metal())
        m2 = combo.AddAtom(self._new_metal())

        for _ in range(n_bridging):
            off   = combo.GetNumAtoms()
            combo = RWMol(Chem.CombineMols(combo, bridging_ligand.mol))
            combo.AddBond(m1, off + bridging_ligand.donor_indices[0], Chem.BondType.DATIVE)
            combo.AddBond(m2, off + bridging_ligand.donor_indices[1], Chem.BondType.DATIVE)

        if axial_ligands is not None:
            if len(axial_ligands) != 2:
                raise ValueError('Pass exactly 2 axial ligands.')
            for midx, ax in zip([m1, m2], axial_ligands):
                off   = combo.GetNumAtoms()
                combo = RWMol(Chem.CombineMols(combo, ax.mol))
                combo.AddBond(midx, off + ax.donor_indices[0], Chem.BondType.DATIVE)

        combo.AddBond(m1, m2, Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(combo, _SANITIZE_NO_PROPS)
        except Exception:
            pass
        for midx in [m1, m2]:
            combo.GetAtomWithIdx(midx).SetFormalCharge(self.charge)

        mol = combo.GetMol()
        print(f'[PaddlewheelBuilder] {self.symbol}₂ | '
              f'{n_bridging} bridging + {len(axial_ligands) if axial_ligands else 0} axial '
              f'| {mol.GetNumAtoms()} atoms')
        return mol


# ── Materials Project integration ────────────────────────────────────────────

def ligand_from_mp(
    mp_id   : str,
    api_key : Optional[str] = None,
) -> Tuple[Optional[str], object]:
    """
    Fetch a ligand from the Materials Project by ID.

    Parameters
    ----------
    mp_id   : MP ID — molecule DB ('mpcule-...') or bulk material ('mp-XXXX')
    api_key : MP API key; falls back to the MP_API_KEY environment variable.

    Returns
    -------
    (smiles, pmg_object)
        smiles     — SMILES string if available, else None
        pmg_object — pymatgen Molecule (molecule DB) or Structure (bulk)

    Usage
    -----
        smiles, _ = ligand_from_mp('mpcule-66383f8c-8dc8-416f-b0de-b8cbdea66b66')
        lig = LigandBuilder(smiles, fg_type='carboxylate', name='BDC-MP')

    Requires: pip install mp-api pymatgen
    """
    try:
        from mp_api.client import MPRester
    except ImportError:
        raise ImportError('pip install mp-api pymatgen')

    key = api_key or _os.environ.get('MP_API_KEY')
    if not key:
        raise ValueError(
            'No MP API key. Pass api_key= or set MP_API_KEY env var.\n'
            'Get a free key at https://next-gen.materialsproject.org/api'
        )

    with MPRester(key) as mpr:
        if not mp_id.startswith('mp-'):
            try:
                docs = mpr.molecules.search(molecule_ids=[mp_id], fields=['molecule', 'smiles'])
                if docs:
                    doc    = docs[0]
                    smiles = getattr(doc, 'smiles', None)
                    pmg    = getattr(doc, 'molecule', None)
                    print(f'[ligand_from_mp] Molecule {mp_id}: SMILES={smiles}')
                    return smiles, pmg
            except Exception as e:
                print(f'[ligand_from_mp] molecules endpoint failed ({e}), trying materials...')

        doc = mpr.materials.get_data_by_id(mp_id, fields=['structure', 'formula_pretty'])
        if doc is None:
            raise ValueError(f'No record found for {mp_id!r} on Materials Project.')
        print(f'[ligand_from_mp] Material {mp_id}: {doc.formula_pretty}')
        return None, doc.structure
