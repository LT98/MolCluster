"""
ebu_core.py -- rebuilt EBU (extended building unit) pipeline core.

Supersedes ebu_tools.py / ebu_tools_v2.py. Same overall shape (LigandBuilder,
EBUBuilder, GeometryPlacer, EBUExporter, enumerate_sbus) so existing notebook
cells migrate with small edits, but donor perception is now generic
(see donor_perception.py) instead of a hardcoded single-fg_type SMARTS
whitelist, and several concrete correctness bugs found in the old code are
fixed here:

  1. H-removal for multi-site same-type ligands used to remove the SAME
     first hydroxyl H over and over. Fixed by scoping H-lookup to each
     individual SMARTS match (donor_perception.py).

  2. Deprotonated donors never got a formal charge, so RDKit's sanitizer
     silently regenerated a phantom implicit H. Fixed by setting
     FormalCharge/NoImplicit explicitly on activation.

  3. The distance-geometry multidentate wrap called RDKit APIs that have
     moved in current RDKit, and was unguarded -- hard-crashed instead of
     running. Fixed with a version-robust shim below.

  4. The phantom-metal atom used by that DG wrap was bonded to donor atoms
     with a real SINGLE bond, over-valencing neutral donors (tertiary amine
     N, aromatic pyridyl N). Fixed with a BondType.ZERO connector.

  5. Multi-site ligands were ALWAYS forced through the same-metal chelate
     DG wrap regardless of whether their donor sites are geometrically
     compatible with chelating one metal. A rigid, widely-separated donor
     pair (e.g. ATF's two catecholate-type O ~7.6 Ang apart across a fused
     ring system, or BDC's/BTC's carboxylates across a benzene ring) is not
     a chelate -- it's a bridging linker meant to reach two different metal
     centers -- and forcing it into one metal's coordination sphere produced
     unphysical torsion of what should be a rigid backbone. Fixed by
     checking chelate-feasibility per donor pair (rigidity + distance)
     before attempting the joint wrap; incompatible ligands fall back to
     single-anchor rigid placement instead of being distorted.

  6. Placing several independent (monodentate) ligands around one metal left
     one rotational degree of freedom per ligand (rotation around its own
     metal-donor axis) completely unconstrained, so e.g. three bipy ligands
     around Fe3+ could rotate straight into each other despite occupying
     distinct coordination sites. Fixed with a cheap azimuthal rotation
     search that picks the least-clashing orientation against
     already-placed atoms.

  7. Nothing checked whether a requested combination was geometrically
     sane at all (e.g. two bulky ligands forced onto adjacent sites of a
     tight 5-coordinate center). geometry_qc.py adds a fast clash/bond/
     ring-planarity check; enumerate_sbus() now discards combinations that
     fail it by default (see `discard_invalid`).
"""
import warnings, numpy as np
import os as _os
import copy
from itertools import product as _iproduct, combinations as _icombinations

from copy import deepcopy
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, RWMol
from rdkit.Geometry import Point3D

from donor_perception import find_donor_sites, activate_sites, DonorSite
from geometry_qc import qc_ebu, QCResult

PT = Chem.GetPeriodicTable()

_SANITIZE_NO_PROPS = (
    Chem.SanitizeFlags.SANITIZE_ALL ^
    Chem.SanitizeFlags.SANITIZE_PROPERTIES
)
_SANITIZE_NO_PROPS_NO_KEKULIZE = (
    _SANITIZE_NO_PROPS ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
)


# ── Distance-geometry API shim ──────────────────────────────────────────────

def _triangle_smooth(bm: np.ndarray) -> None:
    try:
        from rdkit.Chem import rdDistGeom
        rdDistGeom.DoTriangleSmoothing(bm)
        return
    except AttributeError:
        pass
    import rdkit.DistanceGeometry as _DG
    _DG.DoTriangleSmoothing(bm)


def _embed_with_bounds(mol: Chem.Mol, bm: np.ndarray, seed: int) -> int:
    from rdkit.Chem import rdDistGeom
    if hasattr(rdDistGeom, 'EmbedParameters'):
        params = rdDistGeom.EmbedParameters()
        params.SetBoundsMat(bm)
        params.randomSeed = seed
        params.useRandomCoords = True
        return rdDistGeom.EmbedMolecule(mol, params)
    return rdDistGeom.EmbedMolecule(mol, bm, randomSeed=seed)


# Ideal metal-donor-substituent (M-D-A) angle for donors that have exactly
# ONE real bonded neighbor. For these (a deprotonated carboxylate/phenolate/
# sulfonate/etc. O, or a thiolate S) the coordination direction is NOT
# determined by existing bond geometry alone -- one bond fixes only an axis,
# not a specific direction, so an assumed idealized angle is required. Most
# of these are resonance-delocalized, sp2-ish anionic O donors where ~120 deg
# is the standard textbook value (the same angle real monodentate carboxylate
# M-O-C contacts cluster around in crystal structures); alkoxide/thiolate on
# an sp3 carbon skew closer to tetrahedral.
IDEAL_MDA_ANGLE_DEG: Dict[str, float] = {
    'carboxylate_O':  120.0,
    'phenolate_O':    120.0,
    'sulfonate_O':    120.0,
    'sulfinate_O':    115.0,
    'phosphonate_O':  120.0,
    'boronate_O':     120.0,
    'sulfonamide_N':  109.5,
    'azolate_N':      126.0,
    'alkoxide_O':     109.5,
    'thiolate_S':     100.0,
}
DEFAULT_MDA_ANGLE_DEG = 120.0


def _donor_placement_frame(mol: Chem.Mol, conf, donor_idx: int, donor_type: str):
    """
    Determine how to orient a ligand's rigid body so its donor atom's local
    bonding geometry makes chemical sense relative to the metal, instead of
    forcing the donor's one real bond to point straight through the metal
    (a 180 deg M-donor-substituent angle, which is what naively aligning the
    donor->neighbor bond vector to the metal-donor axis produces).

    Returns (mode, ref_hat, tilt_deg):
      mode == 'aligned' -- donor has >=2 real bonded neighbors (heavy atoms
        or H). Their bond geometry already fully determines the coordination
        direction: the "missing" position completing the local trigonal/
        tetrahedral arrangement is -sum(unit bond vectors). Point that
        exactly at the metal (tilt 0) -- this reproduces correct angles to
        ALL existing bonds simultaneously (verified numerically: a pyridyl N
        placed this way lands both M-N-C angles at ~121.7 deg, not the ~180
        deg you'd get by (mis)using the raw, unnormalized mean-of-neighbor-
        positions vector some earlier code used).
      mode == 'tilted' -- donor has exactly 1 real bonded neighbor. Direction
        is ambiguous from bonding alone, so tilt away from that single bond
        by the donor type's ideal M-D-A angle (see IDEAL_MDA_ANGLE_DEG),
        leaving azimuth around that angle's cone free -- that remaining
        rotational DOF is exactly what _best_azimuthal_rotation searches
        over for clash avoidance.
      mode == 'fallback' -- donor has no bonded neighbors at all (shouldn't
        normally occur); use the whole-ligand centroid as a stand-in.
    """
    atom = mol.GetAtomWithIdx(donor_idx)
    d_pos = np.array(conf.GetAtomPosition(donor_idx))
    nbrs = list(atom.GetNeighbors())

    if len(nbrs) >= 2:
        bond_vecs = [np.array(conf.GetAtomPosition(nb.GetIdx())) - d_pos for nb in nbrs]
        unit_vecs = [v / max(np.linalg.norm(v), 1e-9) for v in bond_vecs]
        ref = -sum(unit_vecs)
        norm = np.linalg.norm(ref)
        if norm < 1e-6:
            # degenerate (near-linear substituent arrangement) -- fall back
            # to the single strongest bond direction, tilted mode.
            ref = -unit_vecs[0]
            norm = np.linalg.norm(ref)
            theta = IDEAL_MDA_ANGLE_DEG.get(donor_type, DEFAULT_MDA_ANGLE_DEG)
            return 'tilted', ref / norm, theta
        return 'aligned', ref / norm, 0.0

    if len(nbrs) == 1:
        bond_vec = np.array(conf.GetAtomPosition(nbrs[0].GetIdx())) - d_pos
        ref_hat = bond_vec / max(np.linalg.norm(bond_vec), 1e-9)
        theta = IDEAL_MDA_ANGLE_DEG.get(donor_type, DEFAULT_MDA_ANGLE_DEG)
        return 'tilted', ref_hat, theta

    heavy_pos = np.array([list(conf.GetAtomPosition(a.GetIdx()))
                           for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
    outward = heavy_pos.mean(axis=0) - d_pos
    norm = np.linalg.norm(outward)
    return 'fallback', outward / max(norm, 1e-9), 0.0


def _outward_vector(mol: Chem.Mol, conf, donor_idx: int) -> np.ndarray:
    """Legacy helper, still used by LigandBuilder.get_outward_vectors() for
    diagnostics; NOT used for placement anymore (see _donor_placement_frame)."""
    atom = mol.GetAtomWithIdx(donor_idx)
    d_pos = np.array(conf.GetAtomPosition(donor_idx))
    nb_heavy = [np.array(conf.GetAtomPosition(nb.GetIdx()))
                for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
    if nb_heavy:
        outward = np.mean(nb_heavy, axis=0) - d_pos
    else:
        heavy_pos = np.array([list(conf.GetAtomPosition(a.GetIdx()))
                               for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
        outward = heavy_pos.mean(axis=0) - d_pos
    norm = np.linalg.norm(outward)
    return outward / max(norm, 1e-9)


CN_GEOMETRIES: Dict[int, List[str]] = {
    2: ['planar'],
    3: ['planar'],
    4: ['planar', 'tetrahedral'],
    5: ['planar'],
    6: ['octahedral'],
}


def _site_vectors(geometry: str, n: int, d: float = 2.05) -> np.ndarray:
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


def _kabsch_align(p_local: np.ndarray, targets: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    c_p = p_local.mean(axis=0)
    c_t = targets.mean(axis=0)
    P = p_local - c_p
    T = targets - c_t
    H = P.T @ T
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1., 1., d]) @ U.T
    t = c_t - R @ c_p
    return R, t


def _rot(v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
    f = v_from / np.linalg.norm(v_from)
    t = v_to / np.linalg.norm(v_to)
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


def _rotate_around_axis(points: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    """Rodrigues rotation of an (N,3) point array around a unit `axis` through the origin."""
    k = axis / np.linalg.norm(axis)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    return (points * cos_t
            + np.cross(k, points) * sin_t
            + np.outer(points @ k, k) * (1 - cos_t))


def _rotation_matrix_around_axis(axis: np.ndarray, theta: float) -> np.ndarray:
    """3x3 Rodrigues rotation matrix for angle `theta` (radians) around a unit `axis`."""
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def _arbitrary_perpendicular(v_hat: np.ndarray) -> np.ndarray:
    """Any unit vector perpendicular to `v_hat` (choice is arbitrary -- the
    azimuthal rotation search later sweeps through all such choices anyway)."""
    ref = np.array([1., 0., 0.]) if abs(v_hat[0]) < 0.9 else np.array([0., 1., 0.])
    perp = np.cross(v_hat, ref)
    return perp / np.linalg.norm(perp)


def _best_azimuthal_rotation(
    positions    : np.ndarray,
    axis         : np.ndarray,
    existing_pos : np.ndarray,
    n_samples    : int = 24,
) -> np.ndarray:
    """
    Placing a ligand with a single pinned donor atom leaves one rotational
    degree of freedom (spin around the metal-donor axis) completely
    unconstrained by the alignment math. Left arbitrary, that DOF regularly
    swings a rigid ligand straight into a neighboring one that is occupying
    a geometrically distinct coordination site (observed: bipy ligands
    rotating into each other around Fe3+ despite using separate sites).

    Sample n_samples angles around `axis`, keep whichever orientation
    maximizes the minimum pairwise distance to atoms already placed for
    this EBU (metal + previously placed ligands). Greedy/sequential, not a
    global optimum, but cheap and removes the worst clashes in practice.
    """
    if existing_pos.shape[0] == 0:
        return positions

    best_positions = positions
    best_score = -np.inf
    for theta in np.linspace(0, 2 * np.pi, n_samples, endpoint=False):
        cand = _rotate_around_axis(positions, axis, theta)
        # cheap min-distance via broadcasting (fine at EBU scale: tens of atoms)
        diff = cand[:, None, :] - existing_pos[None, :, :]
        d = np.sqrt((diff ** 2).sum(axis=-1))
        score = d.min()
        if score > best_score:
            best_score = score
            best_positions = cand
    return best_positions


# ── Chelate vs. bridging feasibility ────────────────────────────────────────

def _interior_rotatable_bonds(mol: Chem.Mol, path: List[int]) -> int:
    """
    Count rotatable single bonds strictly BETWEEN two donor atoms' anchor
    points -- excluding the bond directly touching either donor, since
    spinning a donor about its own attachment bond re-orients the donor
    locally but does not change the span between the two donors (this is
    what an anthraquinone-type ligand's exocyclic C-O bonds do: technically
    a rotatable single bond, but rotating it cannot shorten a rigid
    fused-ring backbone's ~7.6 Ang end-to-end distance).
    """
    if len(path) <= 3:
        return 0
    count = 0
    for a1, a2 in zip(path[1:-2], path[2:-1]):
        b = mol.GetBondBetweenAtoms(a1, a2)
        if (b.GetBondType() == Chem.BondType.SINGLE and not b.IsInRing()
                and mol.GetAtomWithIdx(a1).GetDegree() > 1
                and mol.GetAtomWithIdx(a2).GetDegree() > 1):
            count += 1
    return count


def _pair_chelate_compatible(
    mol            : Chem.Mol,
    conf,
    d_i, d_j       : int,
    far_threshold  : float = 3.5,   # a single metal at ~2.05 Ang bonds can span <=~4.1 Ang (trans); beyond this a rigid pair is BRIDGING, not chelating (fixes rigid BTC/BDC ring distortion)
    min_rotatable  : int   = 3,
) -> Tuple[bool, str]:
    """
    Decide whether donor atoms d_i, d_j can plausibly both chelate the SAME
    metal, using the ligand's own free-conformer geometry:

      - if they are already close (<= far_threshold, roughly the diameter of
        a typical single-metal coordination sphere) -- compatible regardless
        of flexibility (covers rigid-but-close cases like true catechol).
      - if they are far apart but the path between them has several
        (>= min_rotatable) interior rotatable bonds -- compatible; a
        flexible chain can fold to reach even if its current/extended
        conformer happens to be spread out (this is EDTA: its default
        embedded conformer is fairly extended, distances up to ~9.9 Ang
        between arms, but 4-7 rotatable bonds mean it can and does fold).
      - otherwise (far AND rigid) -- NOT compatible; this is a bridging
        pair, not a chelating one.
    """
    path = list(Chem.GetShortestPath(mol, d_i, d_j))
    dist = float(np.linalg.norm(
        np.array(conf.GetAtomPosition(d_i)) - np.array(conf.GetAtomPosition(d_j))))
    if dist <= far_threshold:
        return True, f'close ({dist:.1f} Ang)'
    n_rot = _interior_rotatable_bonds(mol, path)
    if n_rot >= min_rotatable:
        return True, f'flexible ({n_rot} interior rotatable bonds, {dist:.1f} Ang extended)'
    return False, f'rigid and far apart ({dist:.1f} Ang, only {n_rot} interior rotatable bond(s))'


def _ligand_chelate_compatible(lig: 'LigandBuilder') -> Tuple[bool, List[str]]:
    """Check every donor pair; ligand is chelate-compatible only if ALL pairs are."""
    if lig.n_sites <= 1:
        return True, []
    conf = lig.mol.GetConformer()
    notes = []
    all_ok = True
    for d_i, d_j in _icombinations(lig.donor_indices, 2):
        ok, reason = _pair_chelate_compatible(lig.mol, conf, d_i, d_j)
        notes.append(f'{d_i}-{d_j}: {reason}')
        all_ok = all_ok and ok
    return all_ok, notes


# ── LigandBuilder ────────────────────────────────────────────────────────────

class LigandBuilder:
    """
    Load a ligand from SMILES and auto-detect every donor site (generic
    perception -- see donor_perception.py). No fg_type is required.

    Parameters
    ----------
    smiles      : ligand SMILES
    donor_types : optional filter -- restrict to these donor_type tags.
    max_sites   : activate at most this many sites (None = all)
    name        : display label
    """

    def __init__(
        self,
        smiles     : str,
        donor_types: Optional[List[str]] = None,
        max_sites  : Optional[int]       = None,
        name       : str                 = 'ligand',
    ):
        self.name      = name
        self.donor_types_filter = donor_types
        self.max_sites = max_sites

        raw = Chem.MolFromSmiles(smiles)
        if raw is None:
            raise ValueError(f'Invalid SMILES: {smiles!r}')

        (self.mol, self.donor_indices, self.attachment_indices,
         self.site_types, self.unique_donor_indices,
         self.unique_attachment_indices, self._canon_ranks) = self._activate(raw)

        self.n_sites  = len(self.donor_indices)
        self.n_unique = len(self.unique_donor_indices)
        type_counts = {}
        for t in self.site_types:
            type_counts[t] = type_counts.get(t, 0) + 1
        types_str = ', '.join(f'{v}x{k}' for k, v in type_counts.items())
        print(f'[LigandBuilder] {name!r} | {self.n_sites} site(s) ({self.n_unique} unique) | '
              f'{types_str} | donors: {self.donor_indices}')

    def _activate(self, raw):
        mol = Chem.AddHs(deepcopy(raw))
        cid = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if cid < 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            pass

        sites = find_donor_sites(mol)
        if self.donor_types_filter is not None:
            sites = [s for s in sites if s.donor_type in self.donor_types_filter]
        if self.max_sites is not None:
            sites = sites[: self.max_sites]

        final_mol, sites = activate_sites(mol, sites)

        try:
            Chem.SanitizeMol(final_mol, _SANITIZE_NO_PROPS_NO_KEKULIZE)
        except Exception:
            pass

        donor_idxs  = [s.donor_idx for s in sites]
        attach_idxs = list(donor_idxs)
        site_types  = [s.donor_type for s in sites]

        canon_ranks = list(Chem.CanonicalRankAtoms(final_mol, breakTies=False))
        seen_keys = set()
        unique_donor_idxs, unique_attach_idxs = [], []
        for d, a in zip(donor_idxs, attach_idxs):
            key = (canon_ranks[d], canon_ranks[a])
            if key not in seen_keys:
                seen_keys.add(key)
                unique_donor_idxs.append(d)
                unique_attach_idxs.append(a)

        return (final_mol, donor_idxs, attach_idxs, site_types,
                unique_donor_idxs, unique_attach_idxs, canon_ranks)

    def get_outward_vectors(self) -> List[np.ndarray]:
        conf = self.mol.GetConformer()
        return [_outward_vector(self.mol, conf, d) for d in self.donor_indices]

    def with_sites(self, n: int) -> 'LigandBuilder':
        c = copy.copy(self)
        c.donor_indices      = self.donor_indices[:n]
        c.attachment_indices = self.attachment_indices[:n]
        c.site_types         = self.site_types[:n]
        c.n_sites            = n
        seen_keys = set()
        c.unique_donor_indices, c.unique_attachment_indices = [], []
        for d, a in zip(c.donor_indices, c.attachment_indices):
            key = (self._canon_ranks[d], self._canon_ranks[a])
            if key not in seen_keys:
                seen_keys.add(key)
                c.unique_donor_indices.append(d)
                c.unique_attachment_indices.append(a)
        c.n_unique = len(c.unique_donor_indices)
        base = self.name.split('[')[0]
        c.name = f'{base}[s{n}]'
        return c

    def with_donor(self, slot: int) -> 'LigandBuilder':
        if slot >= len(self.unique_donor_indices):
            raise IndexError(f'slot {slot} out of range (0..{len(self.unique_donor_indices)-1})')
        c = copy.copy(self)
        d = self.unique_donor_indices[slot]
        a = self.unique_attachment_indices[slot]
        idx_in_full = self.donor_indices.index(d)
        c.donor_indices = [d]
        c.attachment_indices = [a]
        c.site_types = [self.site_types[idx_in_full]]
        c.unique_donor_indices = [d]
        c.unique_attachment_indices = [a]
        c.n_sites = 1
        c.n_unique = 1
        base = self.name.split('[')[0]
        c.name = f'{base}[d{slot}]'
        return c

    def with_types(self, donor_types: List[str]) -> 'LigandBuilder':
        c = copy.copy(self)
        keep = [i for i, t in enumerate(self.site_types) if t in donor_types]
        c.donor_indices      = [self.donor_indices[i] for i in keep]
        c.attachment_indices = [self.attachment_indices[i] for i in keep]
        c.site_types         = [self.site_types[i] for i in keep]
        c.n_sites            = len(c.donor_indices)
        seen_keys = set()
        c.unique_donor_indices, c.unique_attachment_indices = [], []
        for d, a in zip(c.donor_indices, c.attachment_indices):
            key = (self._canon_ranks[d], self._canon_ranks[a])
            if key not in seen_keys:
                seen_keys.add(key)
                c.unique_donor_indices.append(d)
                c.unique_attachment_indices.append(a)
        c.n_unique = len(c.unique_donor_indices)
        base = self.name.split('[')[0]
        c.name = f"{base}[{'+'.join(donor_types)}]"
        return c


METAL_CHARGE: Dict[str, int] = {
    'Zn':2,'Cu':2,'Fe':3,'Co':2,'Ni':2,'Mn':2,'Cd':2,'Ca':2,'Mg':2,
    'Zr':4,'Hf':4,'Al':3,'Cr':3,'In':3,'Eu':3,'Y':3,'Ti':4,'V':3,
}


# ── EBUBuilder ───────────────────────────────────────────────────────────────

class EBUBuilder:
    """Assemble EBU: metal(s) + full ligand molecules joined by DATIVE bonds."""

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
                combo.AddBond(metal_idx, offset + lig.donor_indices[slot], Chem.BondType.DATIVE)

        try:
            Chem.SanitizeMol(combo, _SANITIZE_NO_PROPS_NO_KEKULIZE)
        except Exception:
            pass
        combo.GetAtomWithIdx(metal_idx).SetFormalCharge(self.charge)

        mol        = combo.GetMol()
        n_ml       = sum(min(n, l.n_sites) for l, n in zip(ligands, sites_per_ligand))
        n_expected = 1 + sum(l.mol.GetNumAtoms() for l in ligands)
        status     = 'OK' if mol.GetNumAtoms() == n_expected else f'expected {n_expected}'
        print(f'[EBUBuilder] {self.symbol}({self.charge}+) '
              f'+ {len(ligands)} ligand(s) | {n_ml} M-L bonds | '
              f'{mol.GetNumAtoms()} atoms [{status}]')
        return mol


# ── GeometryPlacer ───────────────────────────────────────────────────────────

class GeometryPlacer:

    @staticmethod
    def _place_ligand_dg(lig: 'LigandBuilder', targets: np.ndarray, d_ml: float) -> Optional[Chem.Mol]:
        k = len(lig.donor_indices)
        rwmol = Chem.RWMol(deepcopy(lig.mol))
        ph_idx = rwmol.AddAtom(Chem.Atom(2))
        for d_i in lig.donor_indices:
            rwmol.AddBond(ph_idx, d_i, Chem.BondType.ZERO)
        try:
            Chem.SanitizeMol(rwmol, _SANITIZE_NO_PROPS_NO_KEKULIZE)
        except Exception:
            pass
        mol = rwmol.GetMol()

        try:
            from rdkit.Chem import rdDistGeom
            bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        except Exception:
            return None

        tol_ml = 0.05
        for d_i in lig.donor_indices:
            lo, hi = sorted([ph_idx, d_i])
            bm[lo][hi] = d_ml + tol_ml
            bm[hi][lo] = d_ml - tol_ml

        tol_dd = 0.3
        for i in range(k):
            for j in range(i + 1, k):
                d_i, d_j = lig.donor_indices[i], lig.donor_indices[j]
                dd = float(np.linalg.norm(targets[i] - targets[j]))
                lo, hi = sorted([d_i, d_j])
                new_ub = dd + tol_dd
                new_lb = max(0.1, dd - tol_dd)
                bm[lo][hi] = min(bm[lo][hi], new_ub)
                bm[hi][lo] = max(bm[hi][lo], new_lb)
                if bm[hi][lo] > bm[lo][hi]:
                    mid = 0.5 * (bm[lo][hi] + bm[hi][lo])
                    bm[lo][hi] = mid + 0.1
                    bm[hi][lo] = mid - 0.1

        try:
            _triangle_smooth(bm)
        except Exception as e:
            warnings.warn(f'[GeometryPlacer] triangle smoothing failed for {lig.name!r}: {e}')
            return None

        cid = -1
        for seed in [42, 123, 7]:
            try:
                cid = _embed_with_bounds(mol, bm, seed)
            except Exception as e:
                warnings.warn(f'[GeometryPlacer] DG embed raised for {lig.name!r} (seed {seed}): {e}')
                cid = -1
            if cid >= 0:
                break
        if cid < 0:
            return None

        emb_conf = mol.GetConformer()
        ph_pos = np.array(emb_conf.GetAtomPosition(ph_idx))

        out = Chem.RWMol(deepcopy(lig.mol))
        out_conf = Chem.Conformer(lig.mol.GetNumAtoms())
        for j in range(lig.mol.GetNumAtoms()):
            p = np.array(emb_conf.GetAtomPosition(j)) - ph_pos
            out_conf.SetAtomPosition(j, p.tolist())
        out.RemoveAllConformers()
        out.AddConformer(out_conf, assignId=True)
        return out.GetMol()

    @staticmethod
    def _rigid_place_single_anchor(lig: 'LigandBuilder', anchor_target: np.ndarray) -> np.ndarray:
        """
        Rigid-body placement pinning ONLY the ligand's first donor to
        `anchor_target`; every other declared donor (if any) is left
        wherever the ligand's own undistorted geometry naturally puts it.

        This is the "bridging" placement mode: used when a multi-site
        ligand's sites are not geometrically compatible with chelating one
        metal (see _ligand_chelate_compatible), and also the normal path for
        any plain monodentate ligand (k==1). The dative bonds to any other
        declared sites still exist in the molecular graph -- this only
        affects 3D geometry, not connectivity/formula -- but their actual
        position will generally land far from d_ml, which the QC step
        (geometry_qc) will correctly flag if that combination is then
        queried as if it were a real multi-point chelate.

        The rotation is built from the donor's own local bonding geometry
        (see _donor_placement_frame) rather than forcing the donor's
        neighbor bond to point straight through the metal -- that naive
        approach produces an exactly-180-degree M-donor-neighbor angle for
        any donor with a single real substituent (e.g. a monodentate
        carboxylate M-O-C contact), which is not physically reasonable
        (verified: such donors should sit around ~120 deg, not linear).
        """
        src_conf = lig.mol.GetConformer()
        donor_idx = lig.donor_indices[0]
        donor_type = lig.site_types[0]
        v_hat = anchor_target / np.linalg.norm(anchor_target)

        mode, ref_hat, tilt_deg = _donor_placement_frame(lig.mol, src_conf, donor_idx, donor_type)
        # Align ref_hat (the donor's own coordination/lone-pair direction) to
        # -v_hat, i.e. the donor->metal direction (metal sits at +v_hat*d_ml
        # from the donor's perspective at the origin-relative target, so the
        # direction FROM the donor back TO the metal is -v_hat). Aligning to
        # +v_hat here would point the donor's bonding geometry away from the
        # metal instead of at it.
        R0 = _rot(ref_hat, -v_hat)
        donor0_local = np.array(src_conf.GetAtomPosition(donor_idx))
        local = np.array([list(src_conf.GetAtomPosition(j))
                          for j in range(lig.mol.GetNumAtoms())])
        heavy = np.array([a.GetAtomicNum() > 1 and a.GetIdx() != donor_idx
                          for a in lig.mol.GetAtoms()])

        def _emit(R):
            t = anchor_target - R @ donor0_local
            return (R @ local.T).T + t

        if mode == 'tilted' and tilt_deg > 1e-6:
            # Setting the M-donor-substituent angle by tilting `tilt_deg` about a
            # perpendicular to the metal-donor axis fixes the ANGLE but leaves the
            # tilt-axis azimuth free -- and that azimuth decides which way the rest
            # of the ligand (e.g. a carboxylate's benzene ring) swings. An arbitrary
            # choice can fold the ligand body straight onto the metal (ring C/H
            # ~1-2 Ang from M). Crucially the later _best_azimuthal_rotation spins
            # about the metal-donor axis, which PRESERVES every atom's distance to
            # the metal, so it cannot undo an inward fold. So choose the tilt-axis
            # azimuth here that pushes the ligand body FARTHEST from the metal.
            base_perp = _arbitrary_perpendicular(v_hat)
            best_pos, best_score = None, -np.inf
            for phi in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                axis = _rotation_matrix_around_axis(v_hat, phi) @ base_perp
                R = _rotation_matrix_around_axis(axis, np.radians(tilt_deg)) @ R0
                pos = _emit(R)
                score = np.linalg.norm(pos[heavy], axis=1).min() if heavy.any() else 0.0
                if score > best_score:
                    best_score, best_pos = score, pos
            return best_pos
        return _emit(R0)

    @staticmethod
    def place(
        ebu_mol      : Chem.Mol,
        metal_symbol : str,
        geometry     : str,
        ligands      : List['LigandBuilder'],
        d_ml         : float = 2.05,
        avoid_clashes: bool = True,
    ) -> Chem.Mol:
        n_total = sum(lig.n_sites for lig in ligands)
        sites   = _site_vectors(geometry, n_total, d=d_ml)
        an      = PT.GetAtomicNumber(metal_symbol)
        metal_idx = next(a.GetIdx() for a in ebu_mol.GetAtoms() if a.GetAtomicNum() == an)

        conf = Chem.Conformer(ebu_mol.GetNumAtoms())
        conf.SetAtomPosition(metal_idx, (0., 0., 0.))

        cursor = 0
        offset = metal_idx + 1
        pinned = [metal_idx]
        n_dg_success = 0
        n_dg_attempted = 0
        n_bridging = 0

        # Heavy-atom positions placed so far, for the azimuthal clash search.
        placed_heavy = [np.zeros((1, 3))]   # metal itself
        lig_records = []   # (start_offset, lig, single_anchor_axis) for refinement

        for lig in ligands:
            k = lig.n_sites
            targets = sites[cursor: cursor + k]
            cursor += k
            start = offset

            single_anchor_axis = None   # set when there's exactly one true rotational DOF
            placed = None

            if k > 1:
                compatible, notes = _ligand_chelate_compatible(lig)
                if not compatible:
                    n_bridging += 1
                    warnings.warn(
                        f"[GeometryPlacer] {lig.name!r}: sites {lig.donor_indices} are not "
                        f"geometrically compatible with chelating one metal ({'; '.join(notes)}). "
                        f"Treating as a bridging ligand -- only donor {lig.donor_indices[0]} is "
                        f"pinned to its target site; the ligand's own rigid geometry determines "
                        f"where its other declared donor(s) end up (they will likely NOT land at "
                        f"bonding distance -- run geometry_qc.qc_ebu to catch this).")
                    positions = GeometryPlacer._rigid_place_single_anchor(lig, targets[0])
                    single_anchor_axis = targets[0] / np.linalg.norm(targets[0])
                else:
                    n_dg_attempted += 1
                    try:
                        placed = GeometryPlacer._place_ligand_dg(lig, targets, d_ml)
                    except Exception as e:
                        warnings.warn(f'[GeometryPlacer] DG placement raised for {lig.name!r}: {e}')
                    if placed is not None:
                        n_dg_success += 1
                        src_conf = placed.GetConformer()
                        donors_dg = np.array([list(src_conf.GetAtomPosition(lig.donor_indices[i]))
                                               for i in range(k)])
                        R, t = _kabsch_align(donors_dg, targets)
                        positions = np.array([R @ np.array(src_conf.GetAtomPosition(j)) + t
                                              for j in range(lig.mol.GetNumAtoms())])
                    else:
                        warnings.warn(f'[GeometryPlacer] DG failed for {lig.name!r}; using rigid '
                                      f'multi-point fallback (Kabsch on the free conformer).')
                        src_conf = lig.mol.GetConformer()
                        donors_local = np.array([list(src_conf.GetAtomPosition(lig.donor_indices[i]))
                                                 for i in range(k)])
                        R, t = _kabsch_align(donors_local, targets)
                        positions = np.array([R @ np.array(src_conf.GetAtomPosition(j)) + t
                                              for j in range(lig.mol.GetNumAtoms())])
            else:
                positions = GeometryPlacer._rigid_place_single_anchor(lig, targets[0])
                single_anchor_axis = targets[0] / np.linalg.norm(targets[0])

            if avoid_clashes and single_anchor_axis is not None:
                existing = np.concatenate(placed_heavy, axis=0)
                heavy_mask = np.array([a.GetAtomicNum() > 1 for a in lig.mol.GetAtoms()])
                if heavy_mask.any():
                    positions = _best_azimuthal_rotation(positions, single_anchor_axis, existing)

            for j, pos in enumerate(positions):
                conf.SetAtomPosition(offset + j, pos.tolist())
            for i in range(k):
                pinned.append(offset + lig.donor_indices[i])

            heavy_pos = positions[np.array([a.GetAtomicNum() > 1 for a in lig.mol.GetAtoms()])]
            placed_heavy.append(heavy_pos)

            offset += lig.mol.GetNumAtoms()
            lig_records.append((start, lig, single_anchor_axis))

        # ── Iterative all-against-all orientation refinement ────────────────
        # The per-ligand placement above is greedy/sequential: each ligand only
        # avoids the atoms of ligands placed BEFORE it, so with several bulky
        # monodentate linkers around one metal the last-placed ones (and any
        # mis-oriented ring) get stuck in avoidable clashes. Re-optimize every
        # single-anchor ligand's spin about its own metal-donor axis against
        # ALL other atoms, a few passes (coordinate descent), so the ligands
        # splay apart cooperatively instead of one-at-a-time.
        if avoid_clashes and sum(1 for _, _, ax in lig_records if ax is not None) > 1:
            cur = {st: np.array([list(conf.GetAtomPosition(st + j))
                                 for j in range(lig.mol.GetNumAtoms())])
                   for st, lig, ax in lig_records}
            for _pass in range(4):
                for st, lig, ax in lig_records:
                    if ax is None:
                        continue
                    others = [np.zeros((1, 3))]
                    for st2, lig2, _ in lig_records:
                        if st2 == st:
                            continue
                        hm2 = np.array([a.GetAtomicNum() > 1 for a in lig2.mol.GetAtoms()])
                        if hm2.any():
                            others.append(cur[st2][hm2])
                    existing = np.concatenate(others, axis=0)
                    cur[st] = _best_azimuthal_rotation(cur[st], ax, existing, n_samples=72)
            for st, lig, ax in lig_records:
                if ax is None:
                    continue
                for j in range(lig.mol.GetNumAtoms()):
                    conf.SetAtomPosition(st + j, cur[st][j].tolist())

        rw = RWMol(ebu_mol)
        rw.RemoveAllConformers()
        rw.AddConformer(conf, assignId=True)
        mol3d = rw.GetMol()

        # NOTE: no force-field relaxation pass here by design. UFF was
        # previously used to locally polish the analytic placement, but it
        # required disguising the metal atom as a neutral carbon so RDKit's
        # organic-only UFF parameter table would accept it -- carbon's
        # standard valence caps at 4, so this silently threw
        # AtomValenceException (caught by a blanket except) for ANY EBU with
        # coordination number >= 5, and separately for any anionic
        # (deprotonated) donor once its dative bond was converted to a real
        # bond (formal charge -1 donors are only allowed 1 real bond, and
        # the conversion gave them 2). In practice that meant relaxation was
        # silently skipped for nearly every candidate this pipeline cares
        # about, most importantly every multidentate chelate. Rather than
        # keep patching a force field that fundamentally can't represent a
        # metal's coordination sphere, we rely solely on the analytic/DG
        # placement above (now with corrected donor bond-angle geometry --
        # see _donor_placement_frame) plus geometry_qc's post-hoc validation.
        # A real ML interatomic potential (see energy_model.py) is the
        # intended path for actual geometry refinement + energy estimation
        # going forward, used as an explicit opt-in step, not silently
        # inside this function.

        msgs = []
        if n_dg_attempted:
            msgs.append(f'{n_dg_success}/{n_dg_attempted} chelate DG wrap(s) succeeded')
        if n_bridging:
            msgs.append(f'{n_bridging} ligand(s) placed as bridging (single-anchor)')
        if msgs:
            print(f'[GeometryPlacer] ' + '; '.join(msgs))

        return mol3d

    @staticmethod
    def to_ase(mol_3d: Chem.Mol):
        import ase
        conf = mol_3d.GetConformer()
        syms = [a.GetSymbol() for a in mol_3d.GetAtoms()]
        pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_3d.GetNumAtoms())])
        atoms = ase.Atoms(symbols=syms, positions=pos)
        atoms.center(vacuum=12.)
        atoms.pbc = False
        return atoms


# ── EBUExporter ──────────────────────────────────────────────────────────────

class EBUExporter:
    @staticmethod
    def _require_conformer(mol: Chem.Mol):
        if mol.GetNumConformers() == 0:
            raise RuntimeError('No 3D conformer -- run GeometryPlacer.place() first.')

    @staticmethod
    def to_xyz(mol_3d: Chem.Mol, path: Optional[str] = None, comment: str = 'EBU') -> str:
        EBUExporter._require_conformer(mol_3d)
        conf = mol_3d.GetConformer()
        lines = [str(mol_3d.GetNumAtoms()), comment]
        for atom in mol_3d.GetAtoms():
            p = conf.GetAtomPosition(atom.GetIdx())
            lines.append(f'{atom.GetSymbol():<4s}  {p.x:14.8f}  {p.y:14.8f}  {p.z:14.8f}')
        xyz = '\n'.join(lines) + '\n'
        if path is not None:
            with open(path, 'w') as fh:
                fh.write(xyz)
            print(f'XYZ -> {path}')
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
    validate       : bool                = True,
    discard_invalid: bool                = True,
    clash_scale    : float               = 0.65,
    ring_max_deviation: float            = 0.25,
) -> Dict[str, Tuple['Chem.Mol', 'Chem.Mol', Optional['QCResult']]]:
    """
    Enumerate SBU configurations. Returns dict label -> (mol_graph, mol_3d, qc_result).

    validate        : run geometry_qc.qc_ebu on each candidate (clash check,
                       M-donor bond sanity, rigid-ring-planarity check).
    discard_invalid : if True (default), combinations that fail QC are NOT
                       written to disk and NOT included in the returned
                       dict -- only their reason is printed. Set False to
                       keep everything (with qc_result attached) for
                       inspection.
    """
    _os.makedirs(out_dir, exist_ok=True)
    sym     = metal_symbol or builder.symbol
    results : Dict[str, Tuple] = {}
    n_discarded = 0

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
            label    = (label_prefix + '_' + '_'.join(lig.name for lig in combo) + f'_{geom}')
            xyz_path = _os.path.join(out_dir, f'{label}.xyz')

            sites = [lig.n_sites for lig in combo]
            mol   = builder.build(combo, sites_per_ligand=sites)
            mol3d = GeometryPlacer.place(mol, sym, geom, combo, d_ml=d_ml)

            qc_result = None
            if validate:
                qc_result = qc_ebu(mol3d, sym, combo, d_ml,
                                    clash_scale=clash_scale, ring_max_deviation=ring_max_deviation)
                if not qc_result.is_valid:
                    print(f'[QC] {label}: {qc_result}')
                    if discard_invalid:
                        n_discarded += 1
                        continue

            EBUExporter.to_xyz(mol3d, path=xyz_path, comment=label)

            if display_fn is not None:
                display_fn(mol, mol3d, label)

            results[label] = (mol, mol3d, qc_result)

    msg = f'{len(results)} EBU(s) written -> {out_dir}'
    if n_discarded:
        msg += f'  ({n_discarded} discarded on QC failure)'
    print(msg)
    return results


# ── Autonomous ligand-count / dentation search ──────────────────────────────
#
# enumerate_sbus() above requires the caller to hand-declare a fixed number
# of "slots" and, for each slot, which with_sites(k)/with_donor(k) variants
# to try -- e.g. Test 3's original setup hard-coded "2 BTC slots + 2 bipy
# slots" and enumerated combinations within that fixed shape. That doesn't
# generalize to "try every sensible number of BTC copies (1, 2, 3, ...) on
# this metal, at every dentation that fits its coordination number" without
# the caller doing that combinatorics by hand ahead of time.
#
# autonomous_enumerate_sbus() does that search itself: given one or more
# ligand SPECIES (full LigandBuilder objects, not pre-sliced with_sites()
# variants) and a metal, it searches over (a) how many copies of each
# species to include, and (b) what dentation each individual copy uses,
# such that the total coordination number lands on one of the metal's valid
# CN values (from `geometries`). This is the direct generalization of
# "1 BTC : 1 metal, 2 BTC : 1 metal, 3 BTC : 1 metal, etc." to any ligand
# set.
#
# A note on charge: a ligand's donor sites are deprotonated once, at
# LigandBuilder construction time (donor_perception.activate_sites) --
# with_sites(k) only changes how many of those already-deprotonated donors
# form a dative bond to THIS metal, it does not reprotonate the rest. So
# e.g. BTC is a fixed -3 anion regardless of whether with_sites(1), (2), or
# (3) is used for binding: the other, non-coordinating carboxylates are
# still deprotonated (in a real crystal they'd typically bind a DIFFERENT
# neighboring metal -- this is exactly the bridging-ligand chemistry
# _ligand_chelate_compatible already flags). That means "make this single,
# mononuclear SBU exactly charge-neutral" is not generally an achievable, or
# even physically meaningful, target for a multiprotic bridging linker like
# BTC -- true charge balance for those happens across the whole periodic
# framework, not one isolated metal center. Given that, charge is reported
# for every candidate (net formal charge = metal charge + sum of each
# included ligand's own fixed charge) but is NOT filtered on by default;
# pass charge_balance=True only for ligand sets where single-SBU neutrality
# is actually the right criterion (e.g. a simple mono-anionic monodentate
# ligand whose entire charge is "used up" by the one site that binds).

def _ligand_formal_charge(lig: 'LigandBuilder') -> int:
    """Net formal charge of this ligand AS BUILT (all its identified donor
    sites deprotonated at construction time) -- independent of with_sites(k),
    since with_sites() only restricts which donors form a dative bond."""
    return Chem.GetFormalCharge(lig.mol)


def _dentation_assignments(max_dentations: List[int], target_total: int):
    """
    Yield every tuple of per-copy dentation (each in [1, max_dentations[i]])
    that sums to target_total. One entry per ligand COPY (not per species) --
    e.g. max_dentations=[3, 3] for two BTC copies, target_total=4 yields
    (1,3), (2,2), (3,1).
    """
    if not max_dentations:
        if target_total == 0:
            yield ()
        return
    first, rest = max_dentations[0], max_dentations[1:]
    rest_min = len(rest)          # each remaining copy needs >=1 site
    rest_max = sum(rest)
    for d in range(1, first + 1):
        remaining = target_total - d
        if remaining < rest_min or remaining > rest_max:
            continue
        for tail in _dentation_assignments(rest, remaining):
            yield (d,) + tail


def _copy_count_vectors(n_species: int, max_copies: int):
    """Yield every (k_0, ..., k_{n-1}) with 0 <= k_i <= max_copies, excluding
    the all-zero vector (need at least one ligand copy somewhere)."""
    ranges = [range(0, max_copies + 1) for _ in range(n_species)]
    for combo in _iproduct(*ranges):
        if any(combo):
            yield combo


def autonomous_enumerate_sbus(
    builder                : 'EBUBuilder',
    ligand_species          : List['LigandBuilder'],
    out_dir                 : str,
    label_prefix            : str,
    geometries                          = CN_GEOMETRIES,
    max_copies_per_species  : int       = 3,
    charge_balance          : bool      = False,
    target_charge           : int       = 0,
    d_ml                    : float     = 2.05,
    validate                : bool      = True,
    discard_invalid         : bool      = True,
    clash_scale             : float     = 0.65,
    ring_max_deviation      : float     = 0.25,
    display_fn              : Optional[callable] = None,
    energy_model                        = None,   # Optional[energy_model.MLEnergyModel]
    relax_with_energy_model : bool      = True,
) -> Dict[str, tuple]:
    """
    Autonomously search over how many copies of each ligand species to use,
    and at what dentation each copy binds, so the total coordination number
    lands on one of `geometries`'s valid CN values -- generalizing e.g.
    "1 BTC : 1 metal, 2 BTC : 1 metal, 3 BTC : 1 metal" to any ligand set
    without hand-declaring with_sites() slots per combination.

    ligand_species: full LigandBuilder objects (NOT pre-sliced with with_sites
    or with_donor) -- one entry per chemically distinct ligand you want the
    search to consider using any number of copies of.

    Returns label -> (mol_graph, mol3d, qc_result, formation_energy_result_or_None),
    same shape as enumerate_sbus but with an extra formation-energy slot.
    """
    _os.makedirs(out_dir, exist_ok=True)
    charges  = [_ligand_formal_charge(l) for l in ligand_species]
    max_dent = [l.n_sites for l in ligand_species]
    n_species = len(ligand_species)

    if isinstance(geometries, dict):
        cn_targets = sorted(geometries.keys())
    else:
        cn_targets = list(range(2, 7))

    results: Dict[str, tuple] = {}
    n_attempted = 0
    n_discarded = 0
    seen_recipes = set()   # avoid duplicate (species,dentation)-multiset / geom combos

    for copy_counts in _copy_count_vectors(n_species, max_copies_per_species):
        net_charge_fixed = builder.charge + sum(k * c for k, c in zip(copy_counts, charges))
        if charge_balance and net_charge_fixed != target_charge:
            continue

        per_copy_species = []
        per_copy_max_dent = []
        for i, k in enumerate(copy_counts):
            per_copy_species += [i] * k
            per_copy_max_dent += [max_dent[i]] * k
        if not per_copy_max_dent:
            continue

        min_possible_cn = len(per_copy_max_dent)
        max_possible_cn = sum(per_copy_max_dent)

        for cn in cn_targets:
            if cn < min_possible_cn or cn > max_possible_cn:
                continue
            for dent_assignment in _dentation_assignments(per_copy_max_dent, cn):
                recipe_key = (tuple(sorted(zip(per_copy_species, dent_assignment))), cn)
                if recipe_key in seen_recipes:
                    continue
                seen_recipes.add(recipe_key)
                n_attempted += 1

                combo = [ligand_species[sp_idx].with_sites(d)
                         for sp_idx, d in zip(per_copy_species, dent_assignment)]

                geom_list = geometries.get(cn, ['planar']) if isinstance(geometries, dict) else [str(geometries)]
                for geom in geom_list:
                    label = (label_prefix + '_' + '_'.join(l.name for l in combo)
                             + f'_{geom}_q{net_charge_fixed:+d}')
                    sites = [l.n_sites for l in combo]
                    mol = builder.build(combo, sites_per_ligand=sites)
                    mol3d = GeometryPlacer.place(mol, builder.symbol, geom, combo, d_ml=d_ml)

                    qc_result = None
                    if validate:
                        qc_result = qc_ebu(mol3d, builder.symbol, combo, d_ml,
                                            clash_scale=clash_scale, ring_max_deviation=ring_max_deviation)
                        if not qc_result.is_valid:
                            print(f'[QC] {label}: {qc_result}')
                            if discard_invalid:
                                n_discarded += 1
                                continue

                    xyz_path = _os.path.join(out_dir, f'{label}.xyz')
                    EBUExporter.to_xyz(mol3d, path=xyz_path, comment=label)
                    if display_fn is not None:
                        display_fn(mol, mol3d, label)

                    formation = None
                    if energy_model is not None:
                        try:
                            formation = energy_model.estimate_formation_energy(
                                mol3d, builder.symbol, combo,
                                relax_ebu=relax_with_energy_model,
                                relax_refs=relax_with_energy_model,
                            )
                            print(f'[Energy] {label}: formation energy ~ '
                                  f'{formation.formation_energy_eV:.3f} eV')
                        except ImportError as e:
                            print(f'[Energy] {label}: skipped ({e})')

                    results[label] = (mol, mol3d, qc_result, formation)

    msg = (f'{len(results)}/{n_attempted} EBU(s) kept -> {out_dir}'
           f'  ({n_discarded} discarded on QC failure)')
    print(msg)

    if energy_model is not None:
        ranked = sorted(
            ((lbl, r[3].formation_energy_eV) for lbl, r in results.items() if r[3] is not None),
            key=lambda x: x[1],
        )
        if ranked:
            print('Ranked by approximate formation energy (lowest = most stable, screening-level only):')
            for lbl, e in ranked:
                print(f'  {e:10.3f} eV  {lbl}')

    return results
