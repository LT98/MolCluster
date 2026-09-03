"""
Generic donor / labile-group perception.

Replaces the old single-`fg_type`-per-molecule SMARTS whitelist. Instead of
requiring every functional group to be pre-registered with a fixed
(group SMARTS, H-removal SMARTS, donor-position-in-match) tuple, this module
scans a molecule in ONE pass and finds every plausible metal-coordinating
site, tagged by chemical type. A molecule with two different donor
chemistries at once (e.g. one end carboxylate, other end pyridyl; or EDTA's
carboxylate arms + backbone amines) is handled automatically -- no need to
pick a single fg_type up front.

Two categories of donor:

  LABILE   -- an acidic proton must be removed to expose an anionic donor
             (carboxylate, phenolate/catecholate, sulfonate, phosphonate,
             boronate, thiolate, azolate ring N-H, ...).  Curated as a flat,
             appendable SMARTS list (pattern convention: the LAST atom in
             the match is always the donor heavy atom, directly bonded to
             the removable H -- this keeps H-lookup unambiguous and scoped
             to that specific match, unlike the old code which re-scanned
             from the top of the global match list every time).

  NEUTRAL  -- an atom that already has an available lone pair and needs no
             modification: aromatic sp2 N with no H (pyridyl-type), sp2
             imine N, nitrile N, sp3 amine N (primary/secondary/tertiary --
             this is what makes EDTA's backbone N's visible), ether O,
             thioether S, phosphine P. Classified by direct valence/
             hybridization/formal-charge rules rather than one SMARTS per
             group, so it generalizes to donor atoms nobody hand-registered
             in advance (e.g. it does NOT require a special-cased pattern
             to find a tertiary amine -- any N with a free lone pair and no
             formal charge qualifies, the same logic a chemist would use).

             Delocalized/unavailable lone pairs are excluded on purpose:
             amide N (conjugated into C=O), protonated N (formal charge),
             aromatic N-H whose lone pair is in the ring pi system (that
             tautomer is instead reachable via the LABILE azolate_N entry,
             i.e. deprotonate to expose the anionic donor).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from rdkit import Chem

# -- Category 1: labile acidic protons -> anionic donor on removal ---------
# (donor_type, SMARTS, formal_charge_after_deprotonation)
# Convention: the donor heavy atom (bonded to the removable H) is the LAST
# atom in the SMARTS match.
LABILE_DONOR_PATTERNS: List[Tuple[str, str, int]] = [
    ('carboxylate_O',  '[CX3](=O)[OX2H1]',              -1),
    ('phenolate_O',    '[c][OX2H1]',                     -1),   # incl. catecholate (2 matches/ring)
    ('sulfonate_O',    '[SX4](=[OX1])(=[OX1])[OX2H1]',   -1),
    ('sulfinate_O',    '[SX3](=[OX1])[OX2H1]',           -1),
    ('phosphonate_O',  '[PX4](=[OX1])([#6,#8])[OX2H1]',  -1),
    ('boronate_O',     '[BX3]([#6,#8])[OX2H1]',          -1),
    ('thiolate_S',     '[#6][SX2H1]',                    -1),
    ('sulfonamide_N',  '[SX4](=[OX1])(=[OX1])[NX3H1]',   -1),
    ('azolate_N',      '[nX3H1]',                        -1),   # pyrazole/imidazole/triazole/tetrazole NH
    ('alkoxide_O',     '[CX4][OX2H1]',                   -1),   # plain aliphatic alcohol (lowest priority)
]

_LABILE_PATTERNS_COMPILED = [
    (name, Chem.MolFromSmarts(smarts), charge) for name, smarts, charge in LABILE_DONOR_PATTERNS
]


@dataclass
class DonorSite:
    donor_idx : int             # atom index of the coordinating atom (current mol indexing)
    donor_type: str             # e.g. 'carboxylate_O', 'pyridyl_N', 'amine_N'
    labile    : bool            # True if this site required H removal to activate
    h_idx     : Optional[int]   # the specific H atom to remove (labile sites only)
    charge_after: int           # formal charge the donor atom should carry once activated


def _find_labile_sites(mol: Chem.Mol) -> List[DonorSite]:
    claimed = set()
    sites: List[DonorSite] = []
    for name, patt, charge in _LABILE_PATTERNS_COMPILED:
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            donor = match[-1]
            if donor in claimed:
                continue
            donor_atom = mol.GetAtomWithIdx(donor)
            h_idx = next((nb.GetIdx() for nb in donor_atom.GetNeighbors()
                          if nb.GetAtomicNum() == 1), None)
            if h_idx is None:
                continue   # pattern matched but no explicit H present (shouldn't happen post-AddHs)
            claimed.add(donor)
            sites.append(DonorSite(donor, name, True, h_idx, charge))
    return sites


def _classify_neutral_atom(atom: Chem.Atom) -> Optional[str]:
    if atom.GetFormalCharge() != 0:
        return None
    sym = atom.GetSymbol()

    if sym == 'N':
        if atom.GetIsAromatic():
            return 'pyridyl_N' if atom.GetTotalNumHs() == 0 else None
        # amide N: lone pair delocalized into an adjacent C=O -> not a good donor
        is_amide = any(
            nb.GetSymbol() == 'C' and any(
                b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(nb).GetSymbol() == 'O'
                for b in nb.GetBonds())
            for nb in atom.GetNeighbors())
        if is_amide:
            return None
        if any(b.GetBondTypeAsDouble() == 3 and b.GetOtherAtom(atom).GetSymbol() == 'C'
               for b in atom.GetBonds()):
            return 'nitrile_N'
        if any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds()):
            return 'imine_N'
        return 'amine_N'   # primary / secondary / tertiary sp3 amine, incl. EDTA-type backbone N

    if sym == 'O':
        if atom.GetTotalNumHs() > 0:
            return None   # still protonated; handled (if acidic) by the labile category
        heavy_nb = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
        if len(heavy_nb) == 2 and all(b.GetBondTypeAsDouble() == 1 for b in atom.GetBonds()):
            return 'ether_O'
        return None

    if sym == 'S':
        if atom.GetTotalNumHs() == 0:
            heavy_nb = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
            if len(heavy_nb) == 2 and all(b.GetBondTypeAsDouble() == 1 for b in atom.GetBonds()):
                return 'thioether_S'
        return None

    if sym == 'P':
        if atom.GetDegree() == 3 and atom.GetTotalNumHs() == 0:
            has_po = any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds())
            if not has_po:
                return 'phosphine_P'
        return None

    return None


def _find_neutral_sites(mol: Chem.Mol, exclude: set) -> List[DonorSite]:
    sites = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in exclude:
            continue
        dtype = _classify_neutral_atom(atom)
        if dtype is not None:
            sites.append(DonorSite(idx, dtype, False, None, atom.GetFormalCharge()))
    return sites


def find_donor_sites(mol: Chem.Mol) -> List[DonorSite]:
    """
    Scan a fully-Hs-explicit RDKit mol and return every donor site found,
    labile and neutral together, sorted by atom index for determinism.
    """
    labile = _find_labile_sites(mol)
    neutral = _find_neutral_sites(mol, exclude={s.donor_idx for s in labile})
    return sorted(labile + neutral, key=lambda s: s.donor_idx)


def activate_sites(mol: Chem.Mol, sites: List[DonorSite]) -> Tuple[Chem.Mol, List[DonorSite]]:
    """
    Remove the labile H for every labile site, set correct formal charges,
    and return a new mol together with the sites re-indexed to match it.

    Uses atom map numbers for safe re-indexing instead of manual
    index-shift arithmetic (the old code's `d - 1 if d > hidx else d`
    approach, done once per removed atom, is exactly the kind of
    bookkeeping that silently drifts wrong on molecules with several
    labile sites -- see EDTA, which has 4).
    """
    rw = Chem.RWMol(mol)

    # Tag every atom with a stable map number = current index + 1 (0 is "unset").
    for atom in rw.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    site_map_nums = [(s, s.donor_idx + 1) for s in sites]
    h_to_remove = sorted({s.h_idx for s in sites if s.labile and s.h_idx is not None}, reverse=True)
    for hidx in h_to_remove:
        rw.RemoveAtom(hidx)

    # Rebuild index lookup from map number -> new index.
    remap = {atom.GetAtomMapNum(): atom.GetIdx() for atom in rw.GetAtoms()}

    new_sites: List[DonorSite] = []
    for s, map_num in site_map_nums:
        new_idx = remap[map_num]
        if s.labile:
            donor_atom = rw.GetAtomWithIdx(new_idx)
            donor_atom.SetFormalCharge(s.charge_after)
            donor_atom.SetNoImplicit(True)
            donor_atom.SetNumExplicitHs(0)
        new_sites.append(DonorSite(new_idx, s.donor_type, s.labile, None, s.charge_after))

    for atom in rw.GetAtoms():
        atom.SetAtomMapNum(0)

    return rw.GetMol(), new_sites
