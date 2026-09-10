"""A site is a FRAME, not a vector (D13).

A single outward vector cannot specify a join: it fixes where the metal sits but says
nothing about the torsion (roll) about the new bond, so the same "site" reproduces
different structures run to run.  A frame — origin, outward axis, and a reference
direction that defines the torsion zero — fixes the relative pose completely, and a
`live_dof` tag says whether that torsion is a real degree of freedom or a don't-care.

`legacy/ebu_core._donor_placement_frame` computed the axis and then threw it away, and
`_best_azimuthal_rotation` searched the torsion stochastically for clash avoidance.
Both are promoted here into a stored, reproducible site model.  One substantive change:
where the legacy code left the azimuth of a single-neighbour donor free for a random
search, this places it in the plane the neighbour's own substituents define — which is
where an sp2 lone pair actually points, and which makes a hydroxyquinone's chelate
pocket come out right instead of by luck.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np
from rdkit import Chem

EPS = 1e-9


class LiveDOF(str, Enum):
    """Whether the torsion about the forming bond is a real coordinate."""

    TORSION_LIVE = "live"      # donor has a plane: syn/anti, cis/trans are distinguishable
    TORSION_FREE = "free"      # near axially symmetric: torsion is a don't-care


class BindingMode(str, Enum):
    MONODENTATE = "mono"
    CHELATE = "chelate"
    BRIDGE_MU2 = "mu2"
    BRIDGE_MU3 = "mu3"


# Ideal metal-donor-attachment angle, degrees.  Derived from donor hybridisation; the
# frame uses it only where the bonding geometry leaves the direction ambiguous.
IDEAL_MDA_ANGLE: dict[str, float] = {
    "carboxylate_O": 120.0, "carbonyl_O": 120.0, "enol_O": 120.0, "enolate_O": 120.0,
    "hydroxide_O": 109.5, "amide_N": 120.0,
    "phenolate_O": 120.0, "alkoxide_O": 109.5, "aqua_O": 104.5, "ether_O": 110.0,
    "amine_N": 109.5, "pyridyl_N": 120.0, "imine_N": 120.0, "nitrile_N": 180.0,
    "azolate_N": 120.0, "sulfonamide_N": 120.0,
    "thiolate_S": 100.0, "thioether_S": 100.0, "phosphine_P": 109.5,
    "sulfonate_O": 120.0, "sulfinate_O": 120.0, "phosphonate_O": 120.0, "boronate_O": 120.0,
    "nitro_O": 120.0,
}
DEFAULT_MDA_ANGLE = 120.0

# A donor with a defined plane has a live torsion; a near-axially-symmetric one does not.
# This tag is the guard against combinatorial explosion: only live torsions branch.
_TORSION_FREE_DONORS = frozenset({"aqua_O", "amine_N", "alkoxide_O", "nitrile_N", "thiolate_S"})

_BINDING_MODES: dict[str, tuple[BindingMode, ...]] = {
    "carboxylate_O": (BindingMode.MONODENTATE, BindingMode.CHELATE, BindingMode.BRIDGE_MU2),
    "carbonyl_O":    (BindingMode.MONODENTATE, BindingMode.CHELATE),
    "enol_O":        (BindingMode.MONODENTATE, BindingMode.CHELATE, BindingMode.BRIDGE_MU2),
    "enolate_O":     (BindingMode.MONODENTATE, BindingMode.CHELATE, BindingMode.BRIDGE_MU2),
    "hydroxide_O":   (BindingMode.MONODENTATE, BindingMode.BRIDGE_MU2, BindingMode.BRIDGE_MU3),
    "phenolate_O":   (BindingMode.MONODENTATE, BindingMode.CHELATE, BindingMode.BRIDGE_MU2),
    "alkoxide_O":    (BindingMode.MONODENTATE, BindingMode.BRIDGE_MU2, BindingMode.BRIDGE_MU3),
    "aqua_O":        (BindingMode.MONODENTATE, BindingMode.BRIDGE_MU2),
    "pyridyl_N":     (BindingMode.MONODENTATE,),
    "imine_N":       (BindingMode.MONODENTATE,),
    "amine_N":       (BindingMode.MONODENTATE,),
    "azolate_N":     (BindingMode.MONODENTATE, BindingMode.BRIDGE_MU2),
    # The other delocalised oxo-acids, listed for the same reason carboxylate is: their
    # oxygens are equivalent, so O,O-chelation and mu2 bridging are available to them in
    # exactly the way they are to a carboxylate.  Before perception typed these groups as
    # groups, a sulfonate's S=O oxygens were reaching this table as `carbonyl_O` and
    # picking up "chelate" from that row by accident; a phosphonate mu3-bridges layered
    # frameworks and is the one entry here that goes further than carboxylate.
    "sulfonate_O":   (BindingMode.MONODENTATE, BindingMode.CHELATE, BindingMode.BRIDGE_MU2),
    "sulfinate_O":   (BindingMode.MONODENTATE, BindingMode.CHELATE, BindingMode.BRIDGE_MU2),
    "phosphonate_O": (BindingMode.MONODENTATE, BindingMode.CHELATE, BindingMode.BRIDGE_MU2,
                      BindingMode.BRIDGE_MU3),
    "nitro_O":       (BindingMode.MONODENTATE, BindingMode.CHELATE),
}


def live_dof(donor_type: str) -> LiveDOF:
    return LiveDOF.TORSION_FREE if donor_type in _TORSION_FREE_DONORS else LiveDOF.TORSION_LIVE


def binding_modes(donor_type: str) -> tuple[BindingMode, ...]:
    return _BINDING_MODES.get(donor_type, (BindingMode.MONODENTATE,))


def torsion_wells(donor_type: str, mode: BindingMode = BindingMode.MONODENTATE) -> tuple[float, ...]:
    """Discrete torsion minima, in degrees.  The conformer-generating coordinate (D11).

    A live torsion on a planar donor has two wells 180 degrees apart (syn / anti); a
    don't-care torsion has one, so it never branches.
    """
    if live_dof(donor_type) is LiveDOF.TORSION_FREE:
        return (0.0,)
    if mode is BindingMode.CHELATE:
        return (0.0,)                       # the chelate ring fixes it
    return (0.0, 180.0)


@dataclass(frozen=True)
class SiteFrame:
    """Position plus two orthogonal directions: enough to fix a join completely."""

    origin: tuple[float, float, float]
    axis_hat: tuple[float, float, float]     # donor -> where the metal goes
    ref_hat: tuple[float, float, float]      # perpendicular to axis; torsion zero
    mode: str                                # aligned | tilted | fallback

    def as_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (np.array(self.origin), np.array(self.axis_hat), np.array(self.ref_hat))

    def to_dict(self) -> dict:
        return {"origin": list(self.origin), "axis": list(self.axis_hat),
                "ref": list(self.ref_hat), "mode": self.mode}

    @classmethod
    def from_dict(cls, d: dict) -> SiteFrame:
        return cls(tuple(d["origin"]), tuple(d["axis"]), tuple(d["ref"]), d["mode"])


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > EPS else np.array([0.0, 0.0, 1.0])


def _perpendicular_to(v: np.ndarray) -> np.ndarray:
    trial = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    return _unit(np.cross(v, trial))


def _rotate(v: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    axis = _unit(axis)
    return (v * np.cos(theta)
            + np.cross(axis, v) * np.sin(theta)
            + axis * np.dot(axis, v) * (1 - np.cos(theta)))


def site_frame(mol: Chem.Mol, donor_idx: int, donor_type: str, conf=None, *,
               heavy_only: bool = False, well: int = 0) -> SiteFrame:
    """The frame for one donor, from the molecule's own 3D geometry.

    `well` selects which torsion well the outward axis sits in.  For a donor whose
    direction is fully determined by its bonding (two or more neighbours) there is one
    well and this is ignored; for an sp2 donor with a single neighbour there are two
    in-plane lone pairs and `well` chooses between them.

    `heavy_only` computes the frame as if the donor were deprotonated: an acidic O-H
    hydrogen sits in the lone-pair direction and drags the outward axis with it, so a
    protonated catechol's two oxygens appear to face away from each other even though the
    catecholate they become is the textbook chelator.  Asking "could this pocket bind
    once activated" therefore has to ignore the proton that activation removes.
    """
    conf = conf if conf is not None else mol.GetConformer()
    atom = mol.GetAtomWithIdx(donor_idx)
    origin = np.array(conf.GetAtomPosition(donor_idx))
    neighbours = list(atom.GetNeighbors())
    if heavy_only:
        heavy = [nb for nb in neighbours if nb.GetAtomicNum() > 1]
        if heavy:
            neighbours = heavy

    if len(neighbours) >= 2:
        # The bonding geometry already determines where the missing position is: the
        # direction that completes the local trigonal/tetrahedral arrangement.
        units = [_unit(np.array(conf.GetAtomPosition(nb.GetIdx())) - origin) for nb in neighbours]
        axis = -sum(units)
        if float(np.linalg.norm(axis)) > 1e-6:
            ref = units[0] - np.dot(units[0], _unit(axis)) * _unit(axis)
            return SiteFrame(tuple(origin), tuple(_unit(axis)), tuple(_unit(ref)), "aligned")
        axis_seed, mode = -units[0], "tilted"           # near-linear: fall through to tilted
    elif len(neighbours) == 1:
        axis_seed, mode = _unit(np.array(conf.GetAtomPosition(neighbours[0].GetIdx())) - origin), \
            "tilted"
    else:
        heavy = np.array([list(conf.GetAtomPosition(a.GetIdx()))
                          for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
        axis = origin - heavy.mean(axis=0) if len(heavy) else np.array([0.0, 0.0, 1.0])
        axis = _unit(axis)
        return SiteFrame(tuple(origin), tuple(axis), tuple(_perpendicular_to(axis)), "fallback")

    # One neighbour: the direction is ambiguous from bonding alone, so tilt away from the
    # bond by the donor's ideal angle.  Choose the rotation axis from the neighbour's own
    # substituents where they exist, which puts the metal in the sp2 plane — where the
    # lone pair is — instead of at an arbitrary azimuth.
    neighbour = neighbours[0]
    others = [nb for nb in neighbour.GetNeighbors() if nb.GetIdx() != donor_idx]
    plane_normal = None
    if others:
        a = _unit(np.array(conf.GetAtomPosition(others[0].GetIdx()))
                  - np.array(conf.GetAtomPosition(neighbour.GetIdx())))
        candidate = np.cross(axis_seed, a)
        if float(np.linalg.norm(candidate)) > 1e-6:
            plane_normal = _unit(candidate)
    if plane_normal is None:
        plane_normal = _perpendicular_to(axis_seed)

    # Rotating by the ideal angle about the plane normal gives the right ANGLE, but the
    # SIDE is a genuine degree of freedom: an sp2 donor has two in-plane lone-pair lobes,
    # syn and anti to its neighbour's substituents, and BOTH are real.  That is exactly
    # what `live_dof = TORSION_LIVE` means, so the side is the torsion well index rather
    # than something to be resolved by a rule.  Picking one arbitrarily made a chelate
    # pocket appear or vanish depending on which substituent RDKit happened to list first;
    # picking the "clearest" one always chose anti and broke every chelate.
    theta = np.radians(IDEAL_MDA_ANGLE.get(donor_type, DEFAULT_MDA_ANGLE))
    sign = 1.0 if well % 2 == 0 else -1.0
    axis = _unit(_rotate(-axis_seed, sign * plane_normal, np.pi - theta))
    ref = _unit(np.cross(plane_normal, axis))
    return SiteFrame(tuple(origin), tuple(axis), tuple(ref), mode)
