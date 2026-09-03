import warnings, numpy as np
import os as _os
import copy
from itertools import product as _iproduct

from copy import deepcopy
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, RWMol
from rdkit.Geometry import Point3D

PT = Chem.GetPeriodicTable()

# SMARTS: (group_smarts, labile_H_smarts, donor_atom_pos, attachment_atom_pos)
FG_DEF: Dict[str, Tuple[str, Optional[str], int, int]] = {
    #                   group SMARTS         H-removal  donor  attachment
    'carboxylate'  : ('[CX3](=O)[OH]',    '[OH]',     2,     0),
    'carboxylate_C': ('[CX3](=O)[OH]',    '[OH]',     0,     0),
    'hydroxy'      : ('[OH]',             '[OH]',     0,     0),
    'pyridyl'      : ('n',                 None,       0,     0),
    'imidazolate'  : ('[nH]',              '[nH]',     0,     0),
    'amine'        : ('[NX3;H2]',          '[NX3;H2]', 0,     0),
    'catecholate'  : ('[c][OH]',           '[OH]',     1,     0),
    'thiol'        : ('[SX2H]',            '[SX2H]',   0,     0),
}

# Default per-CN geometry options used by enumerate_sbus.
# CN=4 gets both square-planar and tetrahedral; CN=6 gets octahedral only.
CN_GEOMETRIES: Dict[int, List[str]] = {
    2: ['planar'],
    3: ['planar'],
    4: ['planar', 'tetrahedral'],
    5: ['planar'],
    6: ['octahedral'],
}


def _site_vectors(geometry: str, n: int, d: float = 2.05) -> np.ndarray:
    """
    Return n coordination-site vectors at distance d from the metal origin.
    tetrahedral and octahedral require n=4 and n=6 respectively.
    """
    g = geometry.lower()
    if g == 'planar':
        angle = np.radians([(360 / n) * i for i in range(n)])
        return d * np.column_stack([np.cos(angle), np.sin(angle), np.zeros(n)])
    elif g == 'tetrahedral':
        if n != 4:
            raise ValueError(f'tetrahedral requires CN=4, got CN={n}')
        return d * np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float) / np.sqrt(3)
    elif g == 'octahedral':
        if n != 6:
            raise ValueError(f'octahedral requires CN=6, got CN={n}')
        return d * np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], float)
    else:
        raise ValueError(f'Unknown geometry {geometry!r}. Options: planar, tetrahedral, octahedral')


def _kabsch_align(
    p_local : np.ndarray,   # (k, 3) donor positions in ligand's MMFF frame
    targets : np.ndarray,   # (k, 3) target donor positions in metal frame
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Kabsch algorithm: optimal rigid-body alignment of k source points onto k targets.
    Returns (R, t) such that R @ p_local[i] + t ≈ targets[i] for all i.
    """
    c_p = p_local.mean(axis=0)
    c_t = targets.mean(axis=0)
    P   = p_local - c_p
    T   = targets  - c_t
    H   = P.T @ T
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1., 1., d]) @ U.T
    t = c_t - R @ c_p
    return R, t


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

        self.mol, self.donor_indices, self.attachment_indices = self._activate(raw)
        self.n_sites = len(self.donor_indices)
        print(f'[LigandBuilder] {name!r} | {fg_type} | '
              f'{self.n_sites} site(s) | donor idx: {self.donor_indices} | '
              f'attach idx: {self.attachment_indices}')

    def _activate(self, raw):
        smarts, h_smarts, donor_pos, attach_pos = FG_DEF[self.fg_type]

        mol = Chem.RWMol(Chem.AddHs(deepcopy(raw)))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)

        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        if self.max_sites is not None:
            matches = matches[: self.max_sites]

        donor_idxs  = []
        attach_idxs = []
        h_remove    = []

        for match in matches:
            donor_idxs.append(match[donor_pos])
            attach_idxs.append(match[attach_pos])
            if h_smarts is not None:
                q2 = Chem.MolFromSmarts(h_smarts)
                for m2 in mol.GetSubstructMatches(q2):
                    candidate = m2[-1]
                    atom = mol.GetAtomWithIdx(candidate)
                    for nb in atom.GetNeighbors():
                        if nb.GetAtomicNum() == 1:
                            h_remove.append(nb.GetIdx())
                            break
                    if len(h_remove) >= len(donor_idxs):
                        break

        for hidx in sorted(set(h_remove), reverse=True):
            mol.RemoveAtom(hidx)
            donor_idxs  = [d - 1 if d > hidx else d for d in donor_idxs]
            attach_idxs = [a - 1 if a > hidx else a for a in attach_idxs]

        try:
            Chem.SanitizeMol(mol,
                Chem.SanitizeFlags.SANITIZE_ALL ^
                Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            pass

        return mol.GetMol(), donor_idxs, attach_idxs

    def get_outward_vectors(self) -> List[np.ndarray]:
        conf = self.mol.GetConformer()
        heavy_pos = np.array([
            list(conf.GetAtomPosition(a.GetIdx()))
            for a in self.mol.GetAtoms()
            if a.GetAtomicNum() > 1
        ])
        centroid = heavy_pos.mean(axis=0)

        vectors = []
        for donor_idx in self.donor_indices:
            donor_pos = np.array(conf.GetAtomPosition(donor_idx))
            outward   = centroid - donor_pos
            norm      = np.linalg.norm(outward)

            if norm < 1e-6:
                atom   = self.mol.GetAtomWithIdx(donor_idx)
                nb_pos = [
                    np.array(conf.GetAtomPosition(b.GetOtherAtomIdx(donor_idx)))
                    for b in atom.GetBonds()
                    if self.mol.GetAtomWithIdx(
                        b.GetOtherAtomIdx(donor_idx)).GetAtomicNum() > 1
                ]
                outward = (np.mean(nb_pos, axis=0) - donor_pos) if nb_pos \
                          else np.array([1., 0., 0.])
                norm    = np.linalg.norm(outward)

            vectors.append(outward / max(norm, 1e-9))

        return vectors

    def with_sites(self, n: int) -> 'LigandBuilder':
        """Return a copy that uses only the first n donor sites."""
        c = copy.copy(self)
        c.donor_indices      = self.donor_indices[:n]
        c.attachment_indices = self.attachment_indices[:n]
        c.n_sites            = n
        base   = self.name.split('[')[0]
        c.name = f'{base}[s{n}]'
        return c

    def with_donor(self, slot: int) -> 'LigandBuilder':
        """Return a copy that uses only the donor at *slot* index."""
        c = copy.copy(self)
        c.donor_indices      = [self.donor_indices[slot]]
        c.attachment_indices = [self.attachment_indices[slot]]
        c.n_sites            = 1
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
    def _build_coord_map(
        ligands   : List[LigandBuilder],
        sites     : np.ndarray,
        metal_idx : int,
    ) -> Tuple[Dict[int, Point3D], Dict[int, np.ndarray]]:
        """
        Build the coordMap for constrained embedding and a parallel numpy dict.

        Returns
        -------
        coord_map_pt3d : {atom_idx: Point3D}  for ETKDGv3
        coord_map_np   : {atom_idx: np.array} for convenience
        """
        coord_map_pt3d = {metal_idx: Point3D(0.0, 0.0, 0.0)}
        coord_map_np   = {metal_idx: np.zeros(3)}

        cursor = metal_idx + 1   # atoms start right after the (single) metal
        site_i = 0
        for lig in ligands:
            for i in range(lig.n_sites):
                global_idx = cursor + lig.donor_indices[i]
                target     = sites[site_i]
                coord_map_pt3d[global_idx] = Point3D(*target.tolist())
                coord_map_np[global_idx]   = target
                site_i += 1
            cursor += lig.mol.GetNumAtoms()

        return coord_map_pt3d, coord_map_np

    @staticmethod
    def _rigid_place(
        ebu_mol      : Chem.Mol,
        metal_symbol : str,
        geometry     : str,
        ligands      : List[LigandBuilder],
        d_ml         : float,
    ) -> Chem.Mol:
        """
        Fallback placement: per-ligand rigid-body transform.
        Single-donor ligands use Rodrigues; multi-donor use Kabsch.
        No conformational flexibility — backbone atoms are not relaxed.
        """
        n_total = sum(lig.n_sites for lig in ligands)
        sites   = _site_vectors(geometry, n_total, d=d_ml)

        an     = PT.GetAtomicNumber(metal_symbol)
        m_idxs = [a.GetIdx() for a in ebu_mol.GetAtoms() if a.GetAtomicNum() == an]

        conf = Chem.Conformer(ebu_mol.GetNumAtoms())
        conf.SetAtomPosition(m_idxs[0], (0., 0., 0.))
        if len(m_idxs) > 1:
            conf.SetAtomPosition(m_idxs[1], (0., 0., 2.6))

        cursor = 0
        offset = len(m_idxs)
        for lig in ligands:
            lig_conf    = lig.mol.GetConformer()
            n_lig_atoms = lig.mol.GetNumAtoms()
            k           = lig.n_sites

            lig_targets  = sites[cursor : cursor + k]
            cursor      += k

            donors_local = np.array([
                list(lig_conf.GetAtomPosition(lig.donor_indices[i]))
                for i in range(k)
            ])

            if k == 1:
                outward = lig.get_outward_vectors()[0]
                target  = lig_targets[0]
                R       = GeometryPlacer._rot(outward, target / np.linalg.norm(target))
                t       = target - R @ donors_local[0]
            else:
                R, t = _kabsch_align(donors_local, lig_targets)

            for j in range(n_lig_atoms):
                p = np.array(lig_conf.GetAtomPosition(j))
                conf.SetAtomPosition(offset + j, (R @ p + t).tolist())

            offset += n_lig_atoms

        rw = RWMol(ebu_mol)
        rw.RemoveAllConformers()
        rw.AddConformer(conf, assignId=True)
        return rw.GetMol()

    @staticmethod
    def _place_ligand_dg(
        lig     : 'LigandBuilder',
        targets : np.ndarray,
        d_ml    : float,
    ) -> Optional[Chem.Mol]:
        """
        Embed a single ligand with donor atoms near target positions using distance
        geometry (RDKit rdDistGeom / bounds-matrix API).

        A phantom He atom is added at the metal site, connected to every donor, and
        the bounds matrix is overridden with tight metal–donor distance constraints
        plus donor–donor distance constraints derived from the target geometry.  DG
        then finds a backbone conformation that simultaneously satisfies all of these
        constraints — naturally "wrapping" flexible ligands (e.g. EDTA) around the
        metal without needing a manually defined bending axis.

        The returned mol has the phantom atom removed; all coordinates are expressed
        in the frame where the metal lies at the origin.  Absolute orientation is
        arbitrary — the caller must Kabsch/Rodrigues-align donor positions onto the
        exact target vectors before placing atoms into the EBU conformer.

        Returns None if embedding fails after three random seeds.
        """
        from rdkit.Chem import rdDistGeom

        k = len(lig.donor_indices)

        # ── Add phantom He atom bonded to every donor ──────────────────────
        rwmol  = Chem.RWMol(deepcopy(lig.mol))
        ph_idx = rwmol.AddAtom(Chem.Atom(2))     # He: no valence rules
        for d_i in lig.donor_indices:
            rwmol.AddBond(ph_idx, d_i, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(rwmol, _SANITIZE_NO_PROPS)
        except Exception:
            pass
        mol = rwmol.GetMol()

        # ── Override bounds matrix ─────────────────────────────────────────
        # Convention: bm[lo][hi]  (lo < hi) = upper bound
        #             bm[hi][lo]  (hi > lo) = lower bound
        bm     = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        tol_ml = 0.05   # Å — tight metal–donor window

        for d_i in lig.donor_indices:
            lo, hi = sorted([ph_idx, d_i])
            bm[lo][hi] = d_ml + tol_ml
            bm[hi][lo] = d_ml - tol_ml

        # Donor–donor bounds from the target geometry distances
        tol_dd = 0.3    # Å — looser to allow backbone flexibility
        for i in range(k):
            for j in range(i + 1, k):
                d_i = lig.donor_indices[i]
                d_j = lig.donor_indices[j]
                dd  = float(np.linalg.norm(targets[i] - targets[j]))
                lo, hi = sorted([d_i, d_j])
                new_ub = dd + tol_dd
                new_lb = max(0.1, dd - tol_dd)
                bm[lo][hi] = min(bm[lo][hi], new_ub)
                bm[hi][lo] = max(bm[hi][lo], new_lb)
                if bm[hi][lo] > bm[lo][hi]:   # guard: lb must not exceed ub
                    mid       = 0.5 * (bm[lo][hi] + bm[hi][lo])
                    bm[lo][hi] = mid + 0.1
                    bm[hi][lo] = mid - 0.1

        rdDistGeom.DoTriangleSmoothing(bm)

        # ── Try a few random seeds ─────────────────────────────────────────
        cid = -1
        for seed in [42, 123, 7]:
            cid = rdDistGeom.EmbedMolecule(mol, bm, randomSeed=seed)
            if cid >= 0:
                break
        if cid < 0:
            return None

        # ── Translate so phantom (metal) sits at origin ────────────────────
        emb_conf = mol.GetConformer()
        ph_pos   = np.array(emb_conf.GetAtomPosition(ph_idx))

        out      = Chem.RWMol(deepcopy(lig.mol))
        out_conf = Chem.Conformer(lig.mol.GetNumAtoms())
        for j in range(lig.mol.GetNumAtoms()):
            p = np.array(emb_conf.GetAtomPosition(j)) - ph_pos
            out_conf.SetAtomPosition(j, p.tolist())
        out.RemoveAllConformers()
        out.AddConformer(out_conf, assignId=True)
        return out.GetMol()

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

        Strategy
        --------
        Per ligand:
          1. Distance-geometry embedding with phantom metal (see _place_ligand_dg):
             produces a backbone conformation that wraps around the metal site.
          2. Kabsch / Rodrigues alignment of the DG-embedded donor positions onto
             the exact target coordination-site vectors.
          3. Fallback to rigid-body placement (old Kabsch on free-ligand conformer)
             if DG fails for a particular ligand.

        Whole-assembly post-processing:
          4. Build a UFF-compatible scratch copy (metal → C, dative → single),
             fix carboxylate resonance so donor-O has SINGLE bond to carboxylate C,
             then UFF-minimise with metal + all donors frozen.
          5. Copy relaxed coordinates back into the original ebu_mol topology.
        """
        n_total   = sum(lig.n_sites for lig in ligands)
        sites     = _site_vectors(geometry, n_total, d=d_ml)
        an        = PT.GetAtomicNumber(metal_symbol)
        metal_idx = next(a.GetIdx() for a in ebu_mol.GetAtoms() if a.GetAtomicNum() == an)

        conf = Chem.Conformer(ebu_mol.GetNumAtoms())
        conf.SetAtomPosition(metal_idx, (0., 0., 0.))

        cursor = 0
        offset = metal_idx + 1   # ligand atoms follow the metal in EBUBuilder order

        for lig in ligands:
            k       = lig.n_sites
            targets = sites[cursor : cursor + k]   # (k, 3)
            cursor += k

            placed = GeometryPlacer._place_ligand_dg(lig, targets, d_ml)

            if placed is not None:
                # Kabsch / Rodrigues: orient DG conformation onto exact target sites
                src_conf    = placed.GetConformer()
                donors_dg   = np.array([list(src_conf.GetAtomPosition(lig.donor_indices[i]))
                                        for i in range(k)])
                if k == 1:
                    d0 = donors_dg[0]
                    R  = GeometryPlacer._rot(d0, targets[0])   # normalises internally
                    t  = targets[0] - R @ d0
                else:
                    R, t = _kabsch_align(donors_dg, targets)
                positions = [R @ np.array(src_conf.GetAtomPosition(j)) + t
                             for j in range(lig.mol.GetNumAtoms())]
            else:
                warnings.warn(f'[GeometryPlacer] DG failed for {lig.name!r}; '
                              'using rigid-body fallback.')
                src_conf     = lig.mol.GetConformer()
                donors_local = np.array([list(src_conf.GetAtomPosition(lig.donor_indices[i]))
                                         for i in range(k)])
                if k == 1:
                    outward = lig.get_outward_vectors()[0]
                    R = GeometryPlacer._rot(outward, targets[0] / np.linalg.norm(targets[0]))
                    t = targets[0] - R @ donors_local[0]
                else:
                    R, t = _kabsch_align(donors_local, targets)
                positions = [R @ np.array(src_conf.GetAtomPosition(j)) + t
                             for j in range(lig.mol.GetNumAtoms())]

            for j, pos in enumerate(positions):
                conf.SetAtomPosition(offset + j, pos.tolist())
            offset += lig.mol.GetNumAtoms()

        # ── Write initial geometry ─────────────────────────────────────────
        rw = RWMol(ebu_mol)
        rw.RemoveAllConformers()
        rw.AddConformer(conf, assignId=True)
        mol3d = rw.GetMol()

        # ── Collect pinned indices (metal + all donors in assembled EBU) ───
        coord_map_pt3d, _ = GeometryPlacer._build_coord_map(ligands, sites, metal_idx)
        pinned = list(coord_map_pt3d.keys())

        # ── Build UFF-compatible scratch copy ──────────────────────────────
        rw2    = RWMol(deepcopy(mol3d))
        m_atom = rw2.GetAtomWithIdx(metal_idx)
        m_atom.SetAtomicNum(6)
        m_atom.SetFormalCharge(0)
        m_atom.SetNoImplicit(False)   # implicit Hs satisfy C valence

        for bond in rw2.GetBonds():
            if bond.GetBondType() == Chem.BondType.DATIVE:
                bond.SetBondType(Chem.BondType.SINGLE)

        # Fix carboxylate resonance: donor-O must carry SINGLE bond to the
        # carboxylate C (RDKit may place C=O on the donor O; after dative→single
        # that gives O explicit valence 3, which sanitization rejects).
        donor_idxs = [idx for idx in pinned if idx != metal_idx]
        for d_idx in donor_idxs:
            d_atom = rw2.GetAtomWithIdx(d_idx)
            if d_atom.GetAtomicNum() != 8:
                continue
            for bond in d_atom.GetBonds():
                nb_idx = bond.GetOtherAtomIdx(d_idx)
                if nb_idx == metal_idx:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    bond.SetBondType(Chem.BondType.SINGLE)
                    nb_atom = rw2.GetAtomWithIdx(nb_idx)
                    for b2 in nb_atom.GetBonds():
                        alt_idx = b2.GetOtherAtomIdx(nb_idx)
                        if alt_idx != d_idx and rw2.GetAtomWithIdx(alt_idx).GetAtomicNum() == 8:
                            b2.SetBondType(Chem.BondType.DOUBLE)
                            break
                    break

        try:
            Chem.SanitizeMol(rw2, _SANITIZE_NO_PROPS)
        except Exception:
            return mol3d   # keep DG geometry if sanitization fails

        # ── UFF backbone relaxation with pinned atoms frozen ───────────────
        try:
            ff = AllChem.UFFGetMoleculeForceField(rw2)
            if ff is not None:
                for idx in pinned:
                    ff.AddFixedPoint(idx)
                ff.Minimize(maxIts=500)
                opt_conf = rw2.GetConformer()
                mol_conf = mol3d.GetConformer()
                for i in range(mol3d.GetNumAtoms()):
                    p = opt_conf.GetAtomPosition(i)
                    mol_conf.SetAtomPosition(i, (p.x, p.y, p.z))
        except Exception:
            pass   # keep DG geometry if UFF fails

        return mol3d

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
    out_dir        : str,
    label_prefix   : str,
    geometries                   = CN_GEOMETRIES,
    constraint     : Optional[callable]  = None,
    metal_symbol   : Optional[str]       = None,
    d_ml           : float               = 2.05,
    display_fn     : Optional[callable]  = None,
) -> Dict[str, Tuple['Chem.Mol', 'Chem.Mol']]:
    """
    Enumerate SBU configurations by Cartesian product over per-ligand variants
    AND over coordination geometries for each total donor count (CN).

    Parameters
    ----------
    builder        : EBUBuilder for the target metal.
    ligand_options : Per-position lists of LigandBuilder variants to try.
                     e.g. [[atf.with_sites(1), atf.with_sites(2)],
                            [edta.with_sites(n) for n in range(1, 5)]]
    out_dir        : Directory for XYZ output files (created if absent).
    label_prefix   : Prepended to every output label and filename.
    geometries     : Controls which geometries to try for each donor-count combo.
                     str                 → same geometry for every combo
                     Dict[int,List[str]] → maps total CN to geometry list;
                                           defaults to CN_GEOMETRIES which tries both
                                           planar and tetrahedral at CN=4, octahedral at CN=6
                     callable(combo)     → returns List[str] for that combo
    constraint     : Optional callable(combo: list[LigandBuilder]) -> bool.
                     Return False to skip a combination entirely.
    metal_symbol   : Override metal symbol for GeometryPlacer (default: builder.symbol).
    d_ml           : Metal–donor bond length in Å.
    display_fn     : Optional callable(mol, mol3d, label) invoked after each EBU is built.

    Returns
    -------
    dict mapping label -> (mol_graph, mol_3d)
    """
    _os.makedirs(out_dir, exist_ok=True)
    sym     = metal_symbol or builder.symbol
    results : Dict[str, Tuple] = {}

    for combo in _iproduct(*ligand_options):
        combo = list(combo)
        if constraint is not None and not constraint(combo):
            continue

        cn = sum(lig.n_sites for lig in combo)

        if callable(geometries):
            geom_list = geometries(combo)
        elif isinstance(geometries, dict):
            geom_list = geometries.get(cn, ['planar'])
        else:
            geom_list = [str(geometries)]

        for geom in geom_list:
            label    = (label_prefix + '_'
                        + '_'.join(lig.name for lig in combo)
                        + f'_{geom}')
            xyz_path = _os.path.join(out_dir, f'{label}.xyz')

            sites = [lig.n_sites for lig in combo]
            mol   = builder.build(combo, sites_per_ligand=sites)
            mol3d = GeometryPlacer.place(mol, sym, geom, combo, d_ml=d_ml)
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
