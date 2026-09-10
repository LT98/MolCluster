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

  GROUP    a delocalised oxo-acid whose oxygens are one donor set rather than several
           atoms that happen to be adjacent.  Matched as a whole group, because
           classifying its oxygens one at a time reads bond order and bond order is
           exactly what resonance moves around.

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

from mofsbu.graph._types import METALS

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
    # acyl sulfonamide before plain sulfonamide: both end on the same N, and the first
    # pattern to claim an atom wins, so the more specific one has to come first or it
    # can never fire.
    ("acyl_sulfonamide_N",
     "[NX3H1;$([NX3H1]([CX3]=[OX1])[SX4](=[OX1])(=[OX1]))]", -1),   # ~4, COOH bioisostere
    ("sulfonamide_N", "[SX4](=[OX1])(=[OX1])[NX3H1]", -1),
    ("azolate_N",     "[nX3H1]", -1),
    ("hydroxamate_O", "[CX3](=[OX1])[NX3;H1][OX2H1]", -1),          # pKa ~9, siderophore motif
    ("oxime_O",       "[#6X3]=[NX2][OX2H1]", -1),                   # ~11
    ("nhydroxy_O",    "[NX3;H0,H1]([#6])[OX2H1]", -1),              # N-hydroxy heterocycles
    # Written as a single-atom recursive match on purpose.  The donor is the LAST atom of
    # the match (see the module docstring), so the obvious spelling
    # `[CX3](=[OX1])[NX3H1][CX3]=[OX1]` ends on a carbonyl oxygen that carries no H: the
    # pattern matches, the h_idx check then discards it, and it looks like coverage while
    # being dead.  `$(...)` puts the environment inside the donor atom instead.
    ("imide_N",       "[NX3H1;$([NX3H1]([CX3]=[OX1])[CX3]=[OX1])]", -1),  # succinimide, ~9
    ("thiophosphate_S", "[PX4](=[OX1,SX1])[SX2H1]", -1),
    ("selenol_Se",    "[#6][SeX2H1]", -1),                          # ~5, beats thiol
    ("peroxy_O",      "[#6][OX2][OX2H1]", -1),                      # hydroperoxide / peracid
    ("thiolate_S",    "[#6][SX2H1]", -1),
    ("alkoxide_O",    "[CX4][OX2H1]", -1),
    ("hydrohalide_X", "[F,Cl,Br,I;H1;X1]", -1),
]

_LABILE = [(name, Chem.MolFromSmarts(smarts), q) for name, smarts, q in LABILE_DONOR_PATTERNS]

# A carboxylate's two oxygens are ONE donor set, not two atoms that happen to sit next to
# each other.  Classified atom by atom they come out differently depending on which
# resonance form RDKit is holding: `[O-]C(=O)C` yields two carboxylate oxygens and
# `[O]=C([O-])C` yields one, and which of those a stored structure gets depends on which
# build route reached it first.  D15 hashes the two forms the SAME on purpose — bond order
# is excluded from the L1 hash precisely so resonance does not split an identity — so a
# site catalog that tells them apart is not a function of the identity it hangs off.  That
# is the seam `registry.api.catalog_drift` was reporting, and this table is what closes it.
#
# The GROUP is matched, never the individual oxygen: the central atom by SMARTS, then every
# terminal oxygen on it, and all of them get one donor type and one charge.  Nothing here
# reads a bond order, which is what makes it resonance-invariant, and metal neighbours are
# ignored so an oxygen that is already coordinated is still part of its own group.
#
# The taxonomy is deliberately COARSE — carbonate's oxygens come out `carboxylate_O`,
# sulfate's `sulfonate_O`, a phosphate diester's `phosphonate_O`.  Each of those shares a
# donor element, a frame model and a set of binding modes with the group it is named for,
# and the alternative is a name per oxo-acid, which is a promise to have anticipated every
# one of them — the mistake `Pocket` refuses to make.  Split a name out when something
# downstream needs to act on the difference, not before.
#
# (donor_type, SMARTS for the group's CENTRAL atom, terminal oxygens the group takes)
DELOCALISED_GROUPS: list[tuple[str, str, int]] = [
    ("carboxylate_O", "[CX3]", 2),      # carboxylic acid, carboxylate, carbonate
    ("sulfonate_O",   "[SX4]", 3),      # sulfonate, sulfate
    ("sulfinate_O",   "[SX3]", 2),      # sulfinate
    ("phosphonate_O", "[PX4]", 2),      # phosphinate, phosphonate, phosphate
    ("nitro_O",       "[NX3]", 2),      # nitro, nitrate
    ("boronate_O",    "[BX3]", 2),      # boronate, borate
]

_GROUPS = [(name, Chem.MolFromSmarts(smarts), n) for name, smarts, n in DELOCALISED_GROUPS]


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


def _terminal_oxygens(mol: Chem.Mol, centre: int) -> list[Chem.Atom]:
    """Oxygens on `centre` with no heavy neighbour of their own but a metal.

    An ester's -OR oxygen has a second heavy neighbour and is not part of the delocalised
    set; a coordinated oxygen's only extra neighbour is the metal it donates to, and that
    is coordination rather than constitution, so it still is.  Protonated oxygens count
    here — a carboxylic ACID is the same group as its carboxylate, which is what lets the
    C=O of `CC(=O)O` be typed by the group while the O-H stays the labile list's business.
    """
    return [
        nb for nb in mol.GetAtomWithIdx(centre).GetNeighbors()
        if nb.GetSymbol() == "O"
        and all(far.GetIdx() == centre or far.GetSymbol() in METALS
                for far in nb.GetNeighbors() if far.GetAtomicNum() > 1)
    ]


def _find_delocalised(mol: Chem.Mol, claimed: set[int]) -> list[DonorSite]:
    """Whole-group donors: every terminal oxygen of an oxo-acid, typed identically."""
    sites: list[DonorSite] = []
    centres: set[int] = set()
    for name, patt, n_oxygens in _GROUPS:
        if patt is None:
            continue
        for (centre,) in mol.GetSubstructMatches(patt, uniquify=True):
            if centre in centres:
                continue
            oxygens = _terminal_oxygens(mol, centre)
            if len(oxygens) < n_oxygens:
                continue
            centres.add(centre)
            # ONE charge for the whole group, read off the group and not off whichever
            # oxygen the resonance form parked the minus sign on.  Anionic or not is the
            # only distinction a donor makes (`LABILE_DONOR_PATTERNS` says -1 per site for
            # a diprotic acid too), so a doubly-deprotonated phosphonate is two -1 donors
            # rather than one atom carrying -2.
            group_charge = (mol.GetAtomWithIdx(centre).GetFormalCharge()
                            + sum(o.GetFormalCharge() for o in oxygens))
            charge_after = -1 if group_charge < 0 else 0
            for oxygen in oxygens:
                if oxygen.GetIdx() in claimed or _n_hydrogens(oxygen) > 0:
                    continue        # still protonated: the labile list activates it
                sites.append(DonorSite(oxygen.GetIdx(), name, False, None, charge_after))
    return sites


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
    if sym in ("F", "Cl", "Br", "I"):
        # A halide arrives already anionic -- as a co-ligand `[Cl-]`, or as what
        # `hydrohalide_X` leaves behind.  Without this branch the site created by
        # deprotonating HX is invisible the moment the structure is read back, which
        # breaks perceive-once (D5): activation would destroy the site it creates.
        return f"halide_{sym}"
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
            # half of the chelate pocket.  An oxo-acid's oxygens never reach this line —
            # `_find_delocalised` has already claimed them as a group, which is what keeps
            # a carboxylate from being read as "one carbonyl and one anion".
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


def find_donor_sites(mol: Chem.Mol) -> list[DonorSite]:
    """Every plausible donor in one pass, sorted by atom index for determinism.

    Three passes, and the order is the point.  The labile acids claim their protonated
    donors first; the delocalised groups then claim whole oxo-acid oxygen sets; the
    per-atom valence rules take what is left.  Groups run BEFORE the per-atom rules
    because the per-atom rules are the ones that read bond order, and an oxygen typed
    from its own bond order is typed from whichever resonance form happened to arrive.
    """
    labile = _find_labile(mol)
    claimed = {s.idx for s in labile}
    grouped = _find_delocalised(mol, claimed)
    claimed |= {s.idx for s in grouped}
    neutral: list[DonorSite] = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in claimed:
            continue
        dtype = _classify_neutral(atom)
        if dtype is not None:
            neutral.append(DonorSite(idx, dtype, False, None, atom.GetFormalCharge()))
    return sorted(labile + grouped + neutral, key=lambda s: s.idx)


def deprotonate(mol: Chem.Mol, sites: list[DonorSite]) -> tuple[Chem.Mol, list[DonorSite]]:
    """Remove the labile H of each given site and set the resulting formal charges.

    Re-indexing goes through atom map numbers
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
