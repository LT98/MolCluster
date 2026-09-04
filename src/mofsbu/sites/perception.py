"""Donor perception: which atoms of a molecule can bind a metal, and how.

Ported from `legacy/donor_perception.py`, which is the best-conditioned code in the old
tree: one pass over the molecule, every plausible donor tagged by chemical type, no
per-molecule functional-group whitelist.  Two categories:

  LABILE   an acidic proton must go before an anionic donor is exposed.  Curated as an
           appendable SMARTS list; by convention the donor heavy atom is the LAST atom
           of the match, which keeps the H lookup scoped to that match.
  NEUTRAL  an available lone pair already.  Classified by valence and hybridisation
           rather than one pattern per group, so it generalises to donors nobody
           registered in advance.

Two donor types are ADDED here that the legacy list does not cover, both of which the
Zn/THQ test system needs and neither of which is exotic:

* **`enol_O`** — a hydroxyl on an sp2 carbon that is part of a C=C.  The legacy labile
  list only had `phenolate_O` (aromatic carbon) and `alkoxide_O` (sp3 carbon).
  Tetrahydroxy-1,4-benzoquinone's ring is a cyclohexadiene-dione, not aromatic, so all
  four of its hydroxyls — the acidic ones, the whole point of the molecule — were
  invisible to the old perception.
* **`carbonyl_O`** — a ketone or quinone oxygen.  The legacy neutral classifier returned
  None for any O with no hydrogens that is not an ether, so C=O was never a donor.  It
  is one, and in a hydroxyquinone it is half of the chelate pocket: the enolate and the
  adjacent carbonyl are what a metal actually binds, exactly as in the anthrarufin-Cu
  case the design doc uses to justify L2.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

from rdkit import Chem

# (donor_type, SMARTS, formal charge on the donor after deprotonation)
# Order matters: earlier patterns claim an atom first.
LABILE_DONOR_PATTERNS: list[tuple[str, str, int]] = [
    ("carboxylate_O", "[CX3](=O)[OX2H1]", -1),
    ("phenolate_O",   "[c][OX2H1]", -1),
    ("enol_O",        "[#6X3]=[#6X3][OX2H1]", -1),      # ADDED — enol/hydroxyquinone OH
    ("sulfonate_O",   "[SX4](=[OX1])(=[OX1])[OX2H1]", -1),
    ("sulfinate_O",   "[SX3](=[OX1])[OX2H1]", -1),
    ("phosphonate_O", "[PX4](=[OX1])([#6,#8])[OX2H1]", -1),
    ("boronate_O",    "[BX3]([#6,#8])[OX2H1]", -1),
    ("thiolate_S",    "[#6][SX2H1]", -1),
    ("sulfonamide_N", "[SX4](=[OX1])(=[OX1])[NX3H1]", -1),
    ("azolate_N",     "[nX3H1]", -1),
    ("alkoxide_O",    "[CX4][OX2H1]", -1),
]

_LABILE = [(name, Chem.MolFromSmarts(smarts), q) for name, smarts, q in LABILE_DONOR_PATTERNS]
_CARBOXYLATE_C = Chem.MolFromSmarts("[CX3](=[OX1])[OX1,OX2]")


@dataclass(frozen=True)
class DonorSite:
    idx: int                  # atom index in the mol this was perceived from
    donor_type: str
    labile: bool
    h_idx: int | None         # the specific H to remove (labile only)
    charge_after: int


def _find_labile(mol: Chem.Mol) -> list[DonorSite]:
    claimed: set[int] = set()
    sites: list[DonorSite] = []
    for name, patt, charge in _LABILE:
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            donor = match[-1]
            if donor in claimed:
                continue
            atom = mol.GetAtomWithIdx(donor)
            h_idx = next((nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1),
                         None)
            if h_idx is None:
                continue
            claimed.add(donor)
            sites.append(DonorSite(donor, name, True, h_idx, charge))
    return sites


def _n_hydrogens(atom: Chem.Atom) -> int:
    """Hydrogens on this atom, explicit or implicit.

    `GetTotalNumHs()` alone counts only IMPLICIT hydrogens, so on the
    hydrogens-explicit molecules this package works with it returns 0 for every atom —
    which silently turns "is this oxygen still protonated?" into "no, never".  The
    SMARTS patterns above are unaffected (SMARTS `H1` counts total hydrogens), which is
    why the bug shows up only in the valence-rule classifier.
    """
    return atom.GetTotalNumHs() + sum(
        1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1)


def _classify_anionic(atom: Chem.Atom) -> str | None:
    """An already-deprotonated donor.

    "Activate = deprotonate" is the core primitive, so a structure recalled from the
    registry is usually ALREADY activated: its donor carries a negative charge and no
    hydrogen.  Classifying only neutral atoms means re-perceiving such a structure finds
    nothing — activation would destroy the very site it creates.  The type is recovered
    from the environment, mirroring the labile patterns that would have produced it.
    """
    if atom.GetFormalCharge() >= 0 or _n_hydrogens(atom) > 0:
        return None
    sym = atom.GetSymbol()
    heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
    if sym == "O":
        if not heavy:
            return "hydroxide_O"
        carbon = heavy[0]
        if carbon.GetSymbol() != "C":
            return {"S": "sulfonate_O", "P": "phosphonate_O",
                    "B": "boronate_O"}.get(carbon.GetSymbol(), "alkoxide_O")
        oxygens = [nb for nb in carbon.GetNeighbors() if nb.GetSymbol() == "O"]
        if len(oxygens) >= 2:
            return "carboxylate_O"
        if carbon.GetIsAromatic():
            return "phenolate_O"
        if any(b.GetBondTypeAsDouble() == 2 for b in carbon.GetBonds()):
            return "enolate_O"
        return "alkoxide_O"
    if sym == "S":
        return "thiolate_S"
    if sym == "N":
        return "amide_N"
    return None


def _classify_neutral(atom: Chem.Atom) -> str | None:
    if atom.GetFormalCharge() != 0:
        return _classify_anionic(atom)
    sym = atom.GetSymbol()

    if sym == "N":
        if atom.GetIsAromatic():
            return "pyridyl_N" if _n_hydrogens(atom) == 0 else None
        is_amide = any(
            nb.GetSymbol() == "C" and any(
                b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(nb).GetSymbol() == "O"
                for b in nb.GetBonds())
            for nb in atom.GetNeighbors())
        if is_amide:
            return None
        if any(b.GetBondTypeAsDouble() == 3 and b.GetOtherAtom(atom).GetSymbol() == "C"
               for b in atom.GetBonds()):
            return "nitrile_N"
        if any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds()):
            return "imine_N"
        return "amine_N"

    if sym == "O":
        n_h = _n_hydrogens(atom)
        if n_h > 0:
            heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
            if not heavy:
                return "aqua_O"                 # free or coordinated water
            if n_h == 1 and len(heavy) == 1:
                return None                     # a hydroxyl: the labile list's business
            return None
        heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
        bonds = atom.GetBonds()
        if len(heavy) == 1 and any(b.GetBondTypeAsDouble() == 2 for b in bonds):
            # ADDED: a ketone / quinone oxygen is a donor, and in a hydroxyquinone it is
            # half of the chelate pocket.  Carboxylate oxygens are re-labelled by the
            # caller, which knows the whole carboxylate group rather than one atom.
            return "carbonyl_O"
        if len(heavy) == 2 and all(b.GetBondTypeAsDouble() == 1 for b in bonds):
            return "ether_O"
        return None

    if sym == "S":
        if _n_hydrogens(atom) == 0:
            heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
            if len(heavy) == 2 and all(b.GetBondTypeAsDouble() == 1 for b in atom.GetBonds()):
                return "thioether_S"
        return None

    if sym == "P":
        if atom.GetDegree() == 3 and _n_hydrogens(atom) == 0:
            if not any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds()):
                return "phosphine_P"
    return None


def _carboxylate_oxygens(mol: Chem.Mol) -> set[int]:
    out: set[int] = set()
    if _CARBOXYLATE_C is None:
        return out
    for match in mol.GetSubstructMatches(_CARBOXYLATE_C, uniquify=True):
        out.update(a for a in match if mol.GetAtomWithIdx(a).GetSymbol() == "O")
    return out


def find_donor_sites(mol: Chem.Mol) -> list[DonorSite]:
    """Every plausible donor in one pass, sorted by atom index for determinism."""
    labile = _find_labile(mol)
    claimed = {s.idx for s in labile}
    carboxylate_o = _carboxylate_oxygens(mol)
    neutral: list[DonorSite] = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in claimed:
            continue
        dtype = _classify_neutral(atom)
        if dtype == "carbonyl_O" and idx in carboxylate_o:
            dtype = "carboxylate_O"
        if dtype is not None:
            neutral.append(DonorSite(idx, dtype, False, None, atom.GetFormalCharge()))
    return sorted(labile + neutral, key=lambda s: s.idx)


def deprotonate(mol: Chem.Mol, sites: list[DonorSite]) -> tuple[Chem.Mol, list[DonorSite]]:
    """Remove the labile H of each given site and set the resulting formal charges.

    Re-indexing goes through atom map numbers rather than index arithmetic: the old
    `d - 1 if d > hidx else d` shuffle, applied once per removed H, drifts wrong on any
    molecule with several labile sites — and THQ has four.
    """
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    tagged = [(s, s.idx + 1) for s in sites]
    for h_idx in sorted({s.h_idx for s in sites if s.labile and s.h_idx is not None},
                        reverse=True):
        rw.RemoveAtom(h_idx)

    remap = {atom.GetAtomMapNum(): atom.GetIdx() for atom in rw.GetAtoms()}
    out: list[DonorSite] = []
    for site, tag in tagged:
        new_idx = remap[tag]
        if site.labile:
            atom = rw.GetAtomWithIdx(new_idx)
            atom.SetFormalCharge(site.charge_after)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        out.append(replace(site, idx=new_idx, h_idx=None))

    for atom in rw.GetAtoms():
        atom.SetAtomMapNum(0)
    mol_out = rw.GetMol()
    Chem.SanitizeMol(mol_out, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    return mol_out, out
