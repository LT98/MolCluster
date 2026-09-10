"""The `Site` record: what perception produces and the registry stores."""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from rdkit import Chem

from mofsbu.sites.frames import BindingMode, LiveDOF, binding_modes, live_dof, site_frame
from mofsbu.sites.perception import DonorSite, find_donor_sites


@dataclass(frozen=True)
class Site:
    atom_idx: int                       # index in the molecule / graph this came from
    donor_type: str
    labile: bool
    charge_after: int
    live_dof: str
    binding_modes: tuple[str, ...]
    frame: dict | None = None


def perceive(mol: Chem.Mol, *, with_frames: bool = True) -> list[Site]:
    """Donor perception plus the frame model (D13), for one molecule."""
    conf = mol.GetConformer() if (with_frames and mol.GetNumConformers()) else None
    out: list[Site] = []
    for donor in find_donor_sites(mol):
        frame = (site_frame(mol, donor.idx, donor.donor_type, conf).to_dict()
                 if conf is not None else None)
        out.append(Site(
            atom_idx=donor.idx,
            donor_type=donor.donor_type,
            labile=donor.labile,
            charge_after=donor.charge_after,
            live_dof=live_dof(donor.donor_type).value,
            binding_modes=tuple(m.value for m in binding_modes(donor.donor_type)),
            frame=frame,
        ))
    return out


MIN_CHELATE_RING = 4          # a carboxylate biting through both its own oxygens
MAX_CHELATE_RING = 7          # beyond this the ring is too floppy to call a chelate
# A pocket is feasible if a metal placed where the two donors' lone pairs meet is at a
# sane distance from both and subtends a sane bite angle.  This replaced a raw
# "the two predictions must land within X angstrom of each other" test, which is not
# scale-aware and which rejected the enolate/enolate pocket of doubly-deprotonated THQ at
# miss = 1.61 against a threshold of 1.60 — a hundredth of an angstrom deciding whether
# the project's central chelate exists.  The free dianion's two O(-) repel and splay,
# which inflates the raw miss without making the pocket unreachable.
MIN_BITE_DEG, MAX_BITE_DEG = 45.0, 125.0
# Calibrated, not guessed.  Measured over rigid pockets that must be found (catechol,
# doubly-deprotonated THQ, 1,10-phenanthroline) and rigid pairs that must not be
# (resorcinol, benzene-1,4-diol): the positives sit at 1.08-1.14 and the negatives at
# 1.55-1.63, so 1.35 separates them with margin on both sides.  Re-measure if this ever
# has to move; a threshold with no measurement behind it is how the enolate/enolate
# pocket got rejected by a hundredth of an angstrom.
MIN_MD_RATIO, MAX_MD_RATIO = 0.55, 1.35      # of the nominal metal-donor distance
# ...and the metal must not be standing inside the ligand.  Without this a para- or
# meta-diol "chelates" through a 6- or 7-ring whose metal sits in the middle of the
# benzene it is supposedly binding to.  Same reasoning as the clash check in QC, applied
# before a pocket is proposed rather than after a structure is built.
MIN_METAL_CLEARANCE = 1.75                   # angstrom to any non-donor heavy atom


@dataclass(frozen=True)
class Pocket:
    """Two donors that can close a chelate ring, described by what they ARE.

    No pocket is given a chemical name.  "Dioxolene", "acac-like" and the rest are a
    taxonomy, and a taxonomy is a promise to have anticipated every motif — the same
    mistake as classifying substitution patterns as ortho/meta/para.  A pocket is
    therefore reported as the ring it would close, how many of its donors are anionic,
    and which donor types are involved.  Selection is by predicate on those facts.
    """

    donors: tuple[int, int]
    donor_types: tuple[str, ...]
    ring_size: int                # atoms in the chelate ring, metal included
    n_anionic: int
    convergence: float            # how far the two donors' implied metal positions miss, A
    n_rotatable: int              # rotatable bonds on the donor-to-donor path
    wells: tuple[int, int] = (0, 0)   # torsion wells at which the donors converge

    @property
    def rigid(self) -> bool:
        return self.n_rotatable == 0

    @property
    def feasible_by(self) -> str:
        """Why this pocket can close: it already points the right way, or it can turn."""
        return "rigid-converged" if self.rigid else f"flexible({self.n_rotatable})"

    @property
    def descriptor(self) -> str:
        return (f"ring{self.ring_size}/{self.n_anionic}anionic/"
                + "+".join(sorted(self.donor_types)))


def _rotatable_on_path(mol: Chem.Mol, path: tuple[int, ...]) -> int:
    """Single, acyclic bonds strictly BETWEEN the two donors — interior bonds only.

    These are what let a ligand turn into its chelating conformation.  A free ligand is
    almost never embedded in that pose — bipyridine's ETKDG conformer is anti and rotates
    to syn only when there is something to bind — so judging a pocket by the free
    conformer's geometry rejects most real chelators.

    The bonds at each END of the path are excluded, and that is not a detail: a donor
    atom lies ON its own bond axis, so rotating that bond moves its substituents and not
    the donor.  Counting them made every hydroxyl look flexible and let a seven-membered
    "chelate" across a benzene ring pass as reachable.
    """
    interior = list(zip(path, path[1:]))[1:-1]
    n = 0
    for u, v in interior:
        bond = mol.GetBondBetweenAtoms(u, v)
        if bond is None:
            continue
        if (bond.GetBondTypeAsDouble() == 1.0 and not bond.IsInRing()
                and not bond.GetIsAromatic()):
            n += 1
    return n


def _convergence(mol: Chem.Mol, conf, a: int, b: int, d_ml: float) -> float:
    """Closest the two donors can agree on a metal position, over their torsion wells.

    Each donor's frame says the metal lies at `origin + d_ml * axis`, and an sp2 donor
    has two such axes (its two lone pairs).  The question is whether SOME combination of
    wells converges, not whether an arbitrarily chosen pair does — so all four are tried
    and the winning wells are returned for the placer to reuse.

    Each donor's frame says the metal lies at `origin + d_ml * axis`.  If the two
    predictions coincide the pair can chelate; if they point apart, no ring closes no
    matter what the path length says.  This is the frame model (D13) earning its keep:
    a lone outward vector per donor would give the same answer only by luck, and the
    path-length test alone happily proposes a seven-membered chelate across a benzene
    ring whose two oxygens face in opposite directions.
    """
    if conf is None:
        return float("inf"), (0, 0), 0.0, 0.0, 0.0
    best, wells, bite, ratio, clear = float("inf"), (0, 0), 0.0, 0.0, 0.0
    for wa, wb in itertools.product((0, 1), repeat=2):
        fa = site_frame(mol, a, mol.GetAtomWithIdx(a).GetSymbol(), conf,
                        heavy_only=True, well=wa)
        fb = site_frame(mol, b, mol.GetAtomWithIdx(b).GetSymbol(), conf,
                        heavy_only=True, well=wb)
        oa, ob = np.array(fa.origin), np.array(fb.origin)
        pa = oa + d_ml * np.array(fa.axis_hat)
        pb = ob + d_ml * np.array(fb.axis_hat)
        miss = float(np.linalg.norm(pa - pb))
        if miss >= best:
            continue
        metal = (pa + pb) / 2.0                     # where a metal would actually sit
        clearance = float("inf")
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            if i in (a, b) or atom.GetAtomicNum() <= 1:
                continue
            clearance = min(clearance,
                            float(np.linalg.norm(metal - np.array(conf.GetAtomPosition(i)))))
        d1, d2 = float(np.linalg.norm(metal - oa)), float(np.linalg.norm(metal - ob))
        v1, v2 = oa - metal, ob - metal
        cos = float(np.dot(v1, v2) / max(np.linalg.norm(v1) * np.linalg.norm(v2), 1e-9))
        best, wells = miss, (wa, wb)
        bite = float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))
        ratio = max(d1, d2) / d_ml if d_ml else 0.0
        clear = clearance
    return best, wells, bite, ratio, clear


def chelate_pockets(
    mol: Chem.Mol,
    sites: list[Site],
    *,
    min_ring: int = MIN_CHELATE_RING,
    max_ring: int = MAX_CHELATE_RING,
    d_ml: float = 2.0,
    check_geometry: bool = True,
) -> list[Pocket]:
    """Donor pairs that could close a chelate ring, found by shortest path.

    Ring size is the through-ligand path plus the metal that closes it, so a carboxylate
    biting through its own two oxygens is a 4-ring, a catecholate or enolate/carbonyl pair
    is a 5-ring, and a salicylate (phenolate + carboxylate across the arene) is a 6-ring.

    The previous version only accepted donors whose attachment atoms were bonded or
    identical, which finds 4- and 5-rings and silently misses every 6-ring — including
    salicylate, whose 6-membered chelate is the reason anyone uses it.
    """
    by_idx = {s.atom_idx: s for s in sites}
    conf = mol.GetConformer() if mol.GetNumConformers() else None
    out: list[Pocket] = []
    for a, b in itertools.combinations(sorted(by_idx), 2):
        # NB: no check that either donor's `binding_modes` lists "chelate".  That table is
        # a per-donor-type taxonomy, and whether two donors can chelate is a property of
        # the PAIR — ring size, rigidity, whether the frames meet — not of either donor
        # alone.  Gating on it silently denied bipyridine and ethylenediamine any pocket
        # at all, which is a strange thing for a chelate finder to do.
        path = Chem.GetShortestPath(mol, a, b)
        if not path:
            continue
        ring_size = len(path) + 1                       # + the metal closing it
        if not (min_ring <= ring_size <= max_ring):
            continue
        convergence, wells, bite, md_ratio, clearance = _convergence(mol, conf, a, b, d_ml)
        n_rotatable = _rotatable_on_path(mol, path)
        # A RIGID path cannot reorganise, so its donors must already face a common point;
        # a path with rotatable bonds can turn into the chelating conformation, so the
        # free conformer's geometry says nothing and only the ring size constrains it.
        # Without this split the test is either useless (accepting a 7-ring across a
        # benzene whose oxygens point apart) or ruinous (rejecting bipy and en).
        if check_geometry and n_rotatable == 0:
            if not (MIN_BITE_DEG <= bite <= MAX_BITE_DEG):
                continue
            if not (MIN_MD_RATIO <= md_ratio <= MAX_MD_RATIO):
                continue
            if clearance < MIN_METAL_CLEARANCE:
                continue                        # the metal would be inside the ligand
        n_anionic = sum(1 for i in (a, b)
                        if mol.GetAtomWithIdx(i).GetFormalCharge() < 0)
        out.append(Pocket((a, b), (by_idx[a].donor_type, by_idx[b].donor_type),
                          ring_size, n_anionic, convergence, n_rotatable, wells))
    return sorted(out, key=lambda p: (p.ring_size, p.donors))


def find_pockets(
    mol: Chem.Mol,
    sites: list[Site],
    *,
    ring_size: int | None = None,
    n_anionic: int | None = None,
    min_anionic: int | None = None,
    donor_types: set[str] | None = None,
    rigid: bool | None = None,
) -> list[Pocket]:
    """Pockets matching a structural predicate.  Ask for properties, not for a name."""
    out = []
    for pocket in chelate_pockets(mol, sites):
        if rigid is not None and pocket.rigid is not rigid:
            continue
        if ring_size is not None and pocket.ring_size != ring_size:
            continue
        if n_anionic is not None and pocket.n_anionic != n_anionic:
            continue
        if min_anionic is not None and pocket.n_anionic < min_anionic:
            continue
        if donor_types is not None and not set(pocket.donor_types) <= donor_types:
            continue
        out.append(pocket)
    return out


def find_pocket(mol: Chem.Mol, sites: list[Site], **predicate) -> Pocket | None:
    matches = find_pockets(mol, sites, **predicate)
    return matches[0] if matches else None


#: A pocket closing a ring this size or larger spans two DIFFERENT functional groups.
#: Below it, the "pocket" is a single group biting through its own two donors — a
#: carboxylate's O,O 4-ring — which is not a neighbour effect at all.
INTERGROUP_RING = 5


def shifting_pocket_donors(mol: Chem.Mol, sites: list[Site]) -> frozenset[int]:
    """Donors whose pKa a NEIGHBOURING group plausibly shifts.  D18's provisional rule.

    Not simply "is in a pocket".  A carboxylate closes a 4-membered ring through its own
    two oxygens, so every carboxylate donor is in a pocket, and flagging all of them says
    the table is wrong about the one value it is most confident in — benzoate's pKa IS
    4.2, and the 4.76 in the table is that number.  Flagging everything is the same as
    flagging nothing.

    What D18 actually means by provisional is the salicylate/anthrarufin case: a donor
    whose acidity is moved by a DIFFERENT group next to it — the peri-OH the quinone
    H-bonds, the phenol ortho to a carboxylate.  Those close 5-, 6- and 7-membered rings,
    so the ring size is the discriminator, and it is a structural fact rather than a
    chemical name (the same reason `Pocket` has a descriptor and not a taxonomy).
    """
    return frozenset(
        idx
        for pocket in chelate_pockets(mol, sites)
        if pocket.ring_size >= INTERGROUP_RING
        for idx in pocket.donors)
