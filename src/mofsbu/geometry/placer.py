"""Placing ligands around a coordination centre.

This is the **N = 1 path of the general operation** (D12): one centre is not a special
case, it is the smallest case.  The multi-centre placer of M6 adds inter-centre
constraints on top of exactly this machinery; nothing here branches on the number of
metals.

Placement is deterministic, driven by the stored site frames rather than a stochastic
search: a ligand is rotated so its donor's outward axis points at the metal, then rolled
to a chosen torsion well.  That is what makes a built structure reproducible from its
choice-vector instead of only from its coordinates.

Pinning the donor and its outward axis does not fix the pose.  Two rotations survive, and
until `placement 2` both were left wherever the alignment arithmetic happened to drop
them:

* the **azimuth**, a spin of the ligand about the M-L axis.  `torsion_wells` names the
  chemically distinguished values of it, but between the wells the ligand is free, and
  that freedom regularly swung a bulky ligand into a neighbour on a different vertex
  (`legacy/ebu_core._best_azimuthal_rotation` was written for exactly this).
* the **out-of-plane angle**, a swing of the metal about the donor-neighbour bond.
  `sites.frames` puts it at zero — the metal in the donor's sp2 plane, where the lone
  pair is — and for a hindered donor that plane is precisely where the substituents are.
  A di-tert-butyl ketone puts a methyl hydrogen 1.41 A from the metal at every azimuth,
  because the metal lies ON the azimuthal axis and no spin can move it.

Both are now searched over a FIXED grid and the winning indices are recorded in the
choice vector, so the search is a deterministic function of its inputs and replay is
exact.  Neither is allowed to wander off its chemistry: the azimuth stays within half a
well spacing of the well it was given, so wells never merge into each other, and the
out-of-plane swing is capped at `MAX_OOP_DEG`.  A donor whose torsion the well table
calls a don't-care has one well, hence a full turn of azimuthal freedom — which is what
"don't-care" means.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from mofsbu._types import AmbiguousSpecError
from mofsbu.geometry.distances import (
    BASE_MO, DEFAULT_BASE_MO, Distance, donor_elements, metal_donor_distance,
)
# Aliased to the spellings this module has always used, so the arithmetic below reads
# exactly as it did when it was written and reviewed.  The definitions moved to
# `_linalg`; nothing about what they compute did.
from mofsbu.geometry._linalg import (
    angle_between as _angle_between, axis_rotation as _axis_rotation, kabsch as _kabsch,
    perpendicular as _perpendicular, rotation_between as _rotation_between, unit as _unit,
)
from mofsbu.geometry.embed import coordinates
from mofsbu.geometry.qc import QCReport, clash_limit, qc
from mofsbu.sites.frames import BindingMode, SiteFrame, site_frame, torsion_wells

# Ideal metal-donor distances, angstrom.  This table is now the M-O BASE of the model in
# `geometry.distances`, not the whole of it: the distance a donor is actually placed at
# depends on the donor's element as well as the metal.  Kept under its old name because
# it is what several tests and the choice vectors of every structure built so far refer
# to, and because as the O-donor case it is still exactly right.
D_ML: dict[str, float] = BASE_MO
DEFAULT_D_ML = DEFAULT_BASE_MO

GEOMETRIES: dict[str, dict[int, str]] = {
    "linear": {2: "linear"}, "bent": {2: "bent"}, "trigonal": {3: "trigonal"},
    "tetrahedral": {4: "tetrahedral"}, "square_planar": {4: "square_planar"},
    "trigonal_bipyramidal": {5: "trigonal_bipyramidal"},
    "square_pyramidal": {5: "square_pyramidal"}, "octahedral": {6: "octahedral"},
}


def site_vectors(geometry: str, n: int, d: float = 2.05, *,
                 angle_deg: float | None = None) -> np.ndarray:
    """Unit coordination directions for a coordination number, scaled to `d`.

    The legacy table mapped CN 5 to 'planar', which is not a coordination geometry;
    both real CN-5 polyhedra are here instead.

    **`bent` is the one entry whose name does not determine its vertices**, so it takes
    `angle_deg` and refuses without it.  Every other polyhedron here fixes its own angles
    — a tetrahedron is 109.47 and nothing may ask it for 104 — which is why the argument
    is rejected for them rather than ignored.  A bent bridge is a real span: a µ2-hydroxide
    sits near 100 deg and a bent µ2-oxo well above it, and a single baked-in number would
    be the kind of tolerance-turned-folklore the M6 risk table names.
    """
    g = geometry.lower()
    if angle_deg is not None and g != "bent":
        raise ValueError(
            f"{geometry} fixes its own angles, so angle_deg={angle_deg} has nothing to "
            f"set. Only 'bent' takes one.")
    s3 = 1 / np.sqrt(3)
    table = {
        "linear": np.array([[1, 0, 0], [-1, 0, 0]], float),
        "bent": None,                      # built below, from the angle the caller states
        "trigonal": np.array([[1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [-0.5, -np.sqrt(3) / 2, 0]]),
        "square_planar": np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], float),
        "tetrahedral": np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) * s3,
        "trigonal_bipyramidal": np.array(
            [[1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [-0.5, -np.sqrt(3) / 2, 0],
             [0, 0, 1], [0, 0, -1]], float),
        "square_pyramidal": np.array(
            [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1]], float),
        # ordered so that consecutive vertices are CIS, not trans — see _assign_targets
        "octahedral": np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0],
                                [0, 0, 1], [0, 0, -1]], float),
    }
    if g not in table:
        raise ValueError(f"unknown coordination geometry {geometry!r}; have {sorted(table)}")
    if g == "bent":
        if angle_deg is None:
            raise AmbiguousSpecError(
                "a bent centre is not determined by its coordination number — the angle "
                "IS the geometry, and it is what separates a mu2-hydroxide near 100 deg "
                "from a bent mu2-oxo far above it. Pass angle_deg; there is no default "
                "that is not a guess about the chemistry (ground rule 5).")
        half = np.radians(float(angle_deg)) / 2.0
        vecs = np.array([[np.cos(half), np.sin(half), 0.0],
                         [np.cos(half), -np.sin(half), 0.0]], float)
    else:
        vecs = table[g]
    if len(vecs) != n:
        raise ValueError(f"{geometry} has {len(vecs)} sites, asked for {n}")
    # `d` may be one distance or one per vertex.  Unit vectors are the honest internal
    # form — the polyhedron is a set of DIRECTIONS, and only the donor decides how far
    # along one it sits — but the scalar signature is kept because the vertices are also
    # used as a length-scaled point set by callers and tests.
    scale = np.asarray(d, dtype=float)
    if scale.ndim == 0:
        return float(scale) * vecs
    if scale.shape != (n,):
        raise ValueError(f"{scale.shape[0]} distances for {n} sites")
    return scale[:, None] * vecs


def _check_reserved(reserve: Sequence[int] | None, cn: int) -> set[int]:
    """`reserve` as a validated set of vertex indices.  An index out of range RAISES.

    Silently dropping one would reserve a different arrangement than the caller asked
    for and report success, which is the failure this whole argument exists to prevent.
    """
    if reserve is None:
        return set()
    out: set[int] = set()
    for v in reserve:
        i = int(v)
        if not 0 <= i < cn:
            raise ValueError(
                f"vertex {i} is not on a CN-{cn} polyhedron, whose vertices are "
                f"0..{cn - 1}. `reserve` indexes site_vectors(geometry, cn).")
        out.add(i)
    return out


def cis_vertices(geometry: str, cn: int, k: int = 2, *,
                 angle_deg: float | None = None) -> tuple[int, ...]:
    """`k` mutually-CIS vertex indices of this polyhedron — the set to `reserve`.

    Cis is what a chelate can span and what a bridge needs (~90 deg on an octahedron);
    trans is what `_assign_targets` leaves behind when nobody says otherwise.  The set
    returned is the one whose widest internal angle is smallest, so "cis" is MEASURED off
    the polyhedron rather than read from a table of index conventions that a reordering of
    `site_vectors` would silently invalidate.

    Ties are broken by the lowest indices, so the answer is deterministic and a build that
    reserves a cis pair replays to the same pair.
    """
    from itertools import combinations

    if k < 2:
        raise ValueError(f"cis is a relation between at least two vertices; got k={k}")
    targets = site_vectors(geometry, cn, 1.0, angle_deg=angle_deg)
    if k > cn:
        raise ValueError(f"{geometry} has {cn} vertices; cannot pick {k} of them")

    def widest(combo: tuple[int, ...]) -> float:
        return max(_angle_between(targets[i], targets[j])
                   for i, j in combinations(combo, 2))

    return min(combinations(range(cn), k), key=lambda c: (widest(c), c))


def reserve_for(geometry: str, cn: int, k: int, ligands: Sequence[LigandPlacement], *,
                angle_deg: float | None = None) -> tuple[int, ...]:
    """The `k` vertices to leave empty so that the ligands can still be placed (B21).

    `cis_vertices` alone chooses the empty set before any chelate is served, and on an
    octahedron that can strand a chelate on a trans pair: four vertices out of six always
    contain a trans pair, so every set ties at 180 deg and the lowest indices leave the
    other trans pair — the only one a chelate cannot span — to the ligand.  Candidates are
    tried in `cis_vertices`'s own order (most cis first, lowest indices on a tie) and the
    first that `_assign_targets` accepts wins, so a set that already worked is unchanged.
    None accepted: `cis_vertices`'s answer, and the placer refuses it with its reason.
    """
    from itertools import combinations

    targets = site_vectors(geometry, cn, 1.0, angle_deg=angle_deg)

    def widest(combo: tuple[int, ...]) -> float:
        return max(_angle_between(targets[i], targets[j])
                   for i, j in combinations(combo, 2))

    for combo in sorted(combinations(range(cn), k), key=lambda c: (widest(c), c)):
        try:
            _assign_targets(targets, list(ligands), reserved=combo)
        except ValueError:
            continue
        return combo
    return cis_vertices(geometry, cn, k, angle_deg=angle_deg)


def bridging_metal_positions(mol: Chem.Mol, atom_idx: int, conf=None, *, geometry: str,
                             n_metals: int, d_m: float,
                             angle_deg: float | None = None) -> np.ndarray:
    """Where the metals go around a SINGLE-atom bridge — mechanism B (WORKPLAN_M6 §4).

    The reframing this function is: **a bridging atom is a centre whose vertices are metal
    positions.**  `site_vectors` never cared which element sits at the middle, so what was
    missing was not machinery but a caller willing to point it at a ligand atom.

    **It reads the bonding rather than a `SiteFrame`, and that is the whole difference
    between the two questions.**  A frame answers "where does ONE metal go", and for a
    single-neighbour donor its axis is a lone-pair lobe tilted ~120 deg off the bond.
    Using that as the bisector of a µ2 bridge swings one of the two metals back toward the
    substituent — measured on hydroxide, 1.68 A from its own proton.  A bridge's bisector
    is the direction AWAY from everything the atom is already bonded to, which is a
    property of its bonding and not of any one lobe.

    Three cases, by how much the bonding constrains the result:

    * **Two or more neighbours** (a bridging aqua): the bisector points away from their
      sum and the metals open PERPENDICULAR to the substituent plane, where water's lone
      pairs are.  Opening in-plane would put a metal on top of an O-H.
    * **One neighbour** (a hydroxide): the bisector is the bond reversed.  The azimuth
      about it is a genuine free rotation — nothing else on the atom orients it — so it is
      chosen deterministically and is arbitrary, which is worth knowing rather than
      hiding.
    * **None** (a bare µ3/µ4 oxo): every direction is free.  Honest rather than sloppy:
      such a bridge's skeleton is fixed by its ANGLES, so every orientation yields the
      same set of M···M distances, which is exactly why µ3 and µ4 land on the literature
      values without anything to orient against.
    """
    conf = conf if conf is not None else mol.GetConformer()
    origin = np.array(conf.GetAtomPosition(atom_idx))
    units = [_unit(np.array(conf.GetAtomPosition(nb.GetIdx())) - origin)
             for nb in mol.GetAtomWithIdx(atom_idx).GetNeighbors()]

    if not units:
        bisector = np.array([1.0, 0.0, 0.0])
        spread = np.array([0.0, 1.0, 0.0])
    else:
        total = sum(units)
        bisector = _unit(-total) if float(np.linalg.norm(total)) > 1e-6 else _unit(-units[0])
        # Perpendicular to the substituent PLANE where there is one; otherwise any
        # perpendicular, deterministically.
        normal = np.cross(units[0], units[1]) if len(units) >= 2 else np.zeros(3)
        spread = (_unit(normal) if float(np.linalg.norm(normal)) > 1e-6
                  else _perpendicular(bisector))
    spread = spread - float(np.dot(spread, bisector)) * bisector
    spread = (_unit(spread) if float(np.linalg.norm(spread)) > 1e-9
              else _perpendicular(bisector))
    # The polyhedron below is symmetric about +x and opens in the xy-plane, so the
    # orientation wanted is (+x -> bisector, +y -> the direction the metals open into).
    rot = np.column_stack([bisector, spread, np.cross(bisector, spread)])
    local = site_vectors(geometry, n_metals, d_m, angle_deg=angle_deg)
    return origin + local @ rot.T


# A five-membered chelate ring (M-D-C-C-D) subtends roughly this at the metal.  A
# polydentate ligand handed two trans vertices cannot close its ring and folds into
# itself, which is a placement bug that shows up only as a pile of clashes downstream.
IDEAL_BITE_ANGLE = 82.0
BITE_ANGLE_RANGE = (55.0, 115.0)


# The orientation grids.  These numbers are part of the stored recipe: an index only
# means something against a grid of a known size, so changing either is an
# `ALGO_VERSIONS['placement']` bump, not a tweak.
N_AZIMUTH = 24          # as many samples as legacy took; the step depends on the wells
N_OOP = 7               # 0, +-15, +-30, +-45

#: How far the metal may be swung out of the donor's sp2 plane.  Not a free parameter:
#: the lone pair really is in that plane, so this is a hindrance allowance, not a claim
#: that the metal binds the pi face.  45 deg is enough for the cases that motivated it
#: (a di-tert-butyl ketone clears its own tert-butyl by 30) and small enough that the
#: donor's own geometry still recognisably decides where the metal goes.
MAX_OOP_DEG = 45.0


def _centred_grid(step: float, n: int) -> list[float]:
    """`n` multiples of `step` ordered 0, +1, -1, +2, -2, ...

    Index 0 is the undisplaced pose, so a stored index of 0 means exactly what placement
    meant before there was a search, and the argmax — which keeps the first of equal
    scores — breaks ties towards the least displaced orientation rather than towards
    whichever angle the grid happened to list first.
    """
    out = [0.0]
    k = 1
    while len(out) < n:
        out.append(+k * step)
        if len(out) < n:
            out.append(-k * step)
        k += 1
    return out[:n]


def _oop_angles() -> list[float]:
    """The out-of-plane grid: in-plane first, then symmetric swings up to MAX_OOP_DEG."""
    return _centred_grid(MAX_OOP_DEG / (N_OOP // 2), N_OOP)


def _azimuth_offsets(n_wells: int) -> list[float]:
    """Azimuthal offsets from a well, spanning half a well spacing either side.

    The bound is what keeps the well meaningful.  Wells 180 deg apart move by just under
    +-90, so no offset can carry one well's orientation onto the other's and two
    structures that differ only in `torsion_well` still differ in their coordinates —
    which D11 relies on when it says near-degenerate isomers survive because their
    choice-vectors differ.  A single-well donor is one the table calls a don't-care, and
    gets the whole turn.
    """
    half = 180.0 / max(n_wells, 1)
    # `N_AZIMUTH + 1` so the extreme offset falls strictly INSIDE the half spacing.  With
    # `2*half/N_AZIMUTH` the last sample lands exactly on the boundary, and that boundary
    # belongs to both wells at once — well 0 reaching +90 and well 1 reaching -90 name the
    # same orientation — which is the one way this grid could still let two wells coincide.
    return _centred_grid(2 * half / (N_AZIMUTH + 1), N_AZIMUTH)


def _clearance_scorer(symbols: list[str], other_symbols: list[str],
                      ignore: np.ndarray):
    """Score a candidate pose by its smallest distance/limit ratio.  Below 1.0 clashes.

    The RATIO is what gets maximised, not the raw separation, because the ratio is what
    QC actually tests (`qc.clash_limit`): C...Cl is allowed 2.14 A and H...Cl only 1.48,
    so the furthest-apart pose and the legal one are regularly different poses.

    `ignore` masks pairs that are not contacts at all — the donor and the metal it is
    bonded to.  Masking them is not tidiness: that pair sits at a fixed distance by
    construction, so leaving it in pins the minimum to a constant, and a constant
    objective makes every candidate score identically and the search silently return
    whatever it looked at first.  Any rotation-invariant pair does that, which is the
    flaw in scoring by raw distance to everything the way `legacy` did.
    """
    limits = np.array([[clash_limit(a, b) for b in other_symbols] for a in symbols])

    def score(placed: np.ndarray, other: np.ndarray) -> float:
        if placed.size == 0 or other.size == 0:
            return float(np.inf)
        d = np.linalg.norm(placed[:, None, :] - other[None, :, :], axis=-1)
        ratio = d / limits
        ratio[ignore] = np.inf
        return float(ratio.min())

    return score


def _assign_targets(targets: np.ndarray, ligands: list, *, reserved: Sequence[int] = (),
                    ) -> tuple[list[np.ndarray], list[int]]:
    """Give each ligand coordination vertices its denticity can actually reach.

    Polydentate ligands are served first and take the vertex set whose angles sit
    closest to a chelate bite angle; monodentate ligands take whatever is left.  Without
    this the vertices are handed out in table order, and on an octahedron that puts a
    bidentate ligand across a trans pair.

    Returns the per-ligand vertex sets AND the indices of the vertices nothing claimed.
    That second value used to be dropped on the floor, which is why an unsaturated centre
    had to be rebuilt at a smaller coordination number instead of simply having empty
    vertices: the information that they existed was computed and then discarded here.

    `reserved` vertices are withheld before any ligand is served and come back among the
    vacancies: this is *which* vertices stay open, as against how many (S4).
    """
    from itertools import combinations

    held = set(reserved)
    remaining = [i for i in range(len(targets)) if i not in held]
    order = sorted(range(len(ligands)), key=lambda i: -ligands[i].denticity)
    assigned: dict[int, np.ndarray] = {}
    for i in order:
        k = ligands[i].denticity
        if k == 1:
            pick = [remaining[0]]
        else:
            def cost(combo: tuple[int, ...]) -> float:
                angles = [_angle_between(targets[a], targets[b])
                          for a, b in combinations(combo, 2)]
                return sum(abs(x - IDEAL_BITE_ANGLE) for x in angles)

            options = list(combinations(remaining, k))
            pick = list(min(options, key=cost))
            worst = max((_angle_between(targets[a], targets[b])
                         for a, b in combinations(pick, 2)), default=0.0)
            if not (BITE_ANGLE_RANGE[0] <= worst <= BITE_ANGLE_RANGE[1]):
                raise ValueError(
                    f"no vertex set on this geometry can host a {k}-dentate ligand: "
                    f"best bite angle {worst:.0f} deg, need "
                    f"{BITE_ANGLE_RANGE[0]:.0f}-{BITE_ANGLE_RANGE[1]:.0f}")
        for v in pick:
            remaining.remove(v)
        assigned[i] = targets[pick]
    # Ascending, so a vacancy's slot number is a property of the polyhedron rather than
    # of the order this function happened to hand vertices out in.
    return [assigned[i] for i in range(len(ligands))], sorted(set(remaining) | held)


def _oop_axis(mol: Chem.Mol, donor_idx: int, origin: np.ndarray, conf) -> np.ndarray | None:
    """The donor-neighbour bond, about which the metal swings out of the donor's plane.

    Rotating about this bond leaves the metal-donor-neighbour angle exactly where the
    frame put it — `IDEAL_MDA_ANGLE` survives the search untouched, only its plane moves.
    A donor with no heavy neighbour (a bare halide) has no such bond and no plane to be
    out of, so it gets None and the swing collapses to a single candidate.
    """
    heavy = [nb for nb in mol.GetAtomWithIdx(donor_idx).GetNeighbors()
             if nb.GetAtomicNum() > 1]
    if not heavy:
        return None
    return _unit(origin - np.array(conf.GetAtomPosition(heavy[0].GetIdx())))


def _monodentate_pose(xyz: np.ndarray, origin: np.ndarray, axis: np.ndarray,
                      target: np.ndarray, oop_axis: np.ndarray | None,
                      oop_deg: float, azimuth_deg: float) -> np.ndarray:
    """Place a single-donor ligand: swing the metal off-plane, aim it, then spin."""
    if oop_axis is not None and oop_deg:
        axis = _axis_rotation(oop_axis, np.radians(oop_deg)) @ axis
    rot = _rotation_between(axis, _unit(-target))
    rot = _axis_rotation(_unit(-target), np.radians(azimuth_deg)) @ rot
    return (rot @ (xyz - origin).T).T + target


@dataclass
class LigandPlacement:
    """One ligand, which of its donors bind, and how."""

    mol: Chem.Mol
    donor_idxs: tuple[int, ...]
    donor_types: tuple[str, ...]
    mode: BindingMode = BindingMode.MONODENTATE
    torsion_well: int = 0
    #: Index into the azimuth / out-of-plane grids, or None to let the placer search.
    #:
    #: These are the two rotations the alignment arithmetic leaves undetermined (see the
    #: module docstring).  `None` means "choose one and tell me which", and is the
    #: default because the alternative — leaving them wherever the arithmetic lands — is
    #: the bug this pair of fields exists to fix.  Passing the integers back reproduces a
    #: stored geometry exactly and skips the search, which is what makes a replay a
    #: replay rather than a re-derivation.
    azimuth_step: int | None = None
    oop_step: int | None = None
    name: str = ""

    @property
    def denticity(self) -> int:
        return len(self.donor_idxs)


#: How many times the refinement sweeps the ligands.  Part of the stored recipe.
#:
#: One sweep is not enough on its own and zero is what the bug was: a ligand can only be
#: scored against ligands that already exist, so in a single forward pass the first
#: ligand placed is optimised against nothing at all.  That is exactly the losing case —
#: a bulky ketone assigned the first vertex would keep its arbitrary azimuth and the
#: co-ligands would then be dropped on top of it.  A second sweep gives every ligand a
#: fully-placed neighbourhood to avoid.  Sweeps stop early once nothing moves, so the
#: usual cost is two.
N_REFINE_SWEEPS = 2


@dataclass
class _Monodentate:
    """Everything needed to re-pose one single-donor ligand during refinement."""

    base: int
    n_atoms: int
    xyz: np.ndarray
    origin: np.ndarray
    axis: np.ndarray
    target: np.ndarray
    oop_axis: np.ndarray | None
    well_deg: float
    az_grid: list[float]
    oop_grid: list[float]
    az_fixed: int | None
    oop_fixed: int | None
    donor_atom: int
    choice: dict

    @property
    def searching(self) -> bool:
        return self.az_fixed is None or self.oop_fixed is None


@dataclass
class PlacementResult:
    symbols: list[str]
    coords: np.ndarray
    metal_idx: int
    donor_idxs: list[int]
    bonded: set[tuple[int, int]]
    report: QCReport
    choice_vector: dict = field(default_factory=dict)
    atom_offsets: list[int] = field(default_factory=list)
    #: Which ligand each atom came from, parallel to `symbols`.  A clash is only
    #: interpretable if you know whether it is intra-ligand or between two ligands.
    owners: list[str] = field(default_factory=list)
    #: The M-L distance each donor was placed at, and where that number came from.
    donor_distances: list[Distance] = field(default_factory=list)
    #: Unit directions of the coordination vertices NOTHING was placed on.
    #:
    #: A coordinatively unsaturated centre is not a smaller centre.  Two ketones on a
    #: tetrahedral Mg leave two vacant vertices at tetrahedral angles; calling that a
    #: 2-coordinate linear complex describes a different molecule, and it is the one the
    #: next assembly step would build on.  These directions are what makes the vacancy
    #: addressable: origin at the metal, axis along the empty vertex — a frame, the same
    #: shape as any other site (D13).
    vacancies: tuple[tuple[float, float, float], ...] = ()

    @property
    def ok(self) -> bool:
        return self.report.ok

    @property
    def cn(self) -> int:
        """Coordination number of the POLYHEDRON, occupied vertices plus vacant ones."""
        return len(self.donor_idxs) + len(self.vacancies)

    def to_xyz(self, comment: str = "") -> str:
        lines = [str(len(self.symbols)), comment]
        for sym, (x, y, z) in zip(self.symbols, self.coords):
            lines.append(f"{sym:<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(lines) + "\n"


def _refine_azimuths(arr: np.ndarray, symbols: list[str],
                     refinable: list[_Monodentate]) -> None:
    """Choose each single-donor ligand's two free rotations.  Edits `arr` in place.

    Deterministic throughout: fixed grids, fixed sweep order, and `>` rather than `>=`
    so an unbeaten score keeps the earliest — hence least displaced — candidate.  Nothing
    here consults the seed, so the result is a function of the choice vector alone, and
    feeding the emitted indices back in reproduces the coordinates without searching.
    """
    if not any(r.searching for r in refinable):
        return
    for _sweep in range(N_REFINE_SWEEPS):
        moved = False
        for r in refinable:
            if not r.searching:
                continue
            mine = np.arange(r.base, r.base + r.n_atoms)
            others = np.setdiff1d(np.arange(len(symbols)), mine)
            other_symbols = [symbols[i] for i in others]
            # The one pair here that is a BOND, not a contact.  See `_clearance_scorer`:
            # left in, its fixed length pins the minimum and flattens the whole search.
            ignore = np.zeros((r.n_atoms, len(others)), dtype=bool)
            ignore[r.donor_atom - r.base, np.searchsorted(others, 0)] = True
            score_fn = _clearance_scorer([symbols[i] for i in mine], other_symbols, ignore)

            az_range = (range(len(r.az_grid)) if r.az_fixed is None
                        else [r.az_fixed % len(r.az_grid)])
            oop_range = (range(len(r.oop_grid)) if r.oop_fixed is None
                         else [r.oop_fixed % len(r.oop_grid)])
            best: tuple[float, int, int, np.ndarray] | None = None
            for oi in oop_range:
                for ai in az_range:
                    placed = _monodentate_pose(
                        r.xyz, r.origin, r.axis, r.target, r.oop_axis,
                        r.oop_grid[oi], r.well_deg + r.az_grid[ai])
                    score = score_fn(placed, arr[others])
                    if best is None or score > best[0]:
                        best = (score, ai, oi, placed)
            assert best is not None
            _score, ai, oi, placed = best
            if (ai, oi) != (r.choice["azimuth_step"], r.choice["oop_step"]):
                moved = True
            arr[r.base:r.base + r.n_atoms] = placed
            r.choice["azimuth_step"] = ai
            r.choice["oop_step"] = oi
            r.choice["azimuth_deg"] = round(r.well_deg + r.az_grid[ai], 3)
            r.choice["oop_deg"] = round(r.oop_grid[oi], 3)
        if not moved:
            break


def place_mononuclear(
    metal: str,
    ligands: list[LigandPlacement],
    *,
    geometry: str,
    cn: int | None = None,
    d_ml: float | None = None,
    reserve: Sequence[int] | None = None,
    seed: int = 0,
) -> PlacementResult:
    """Assemble one coordination centre from ligands whose frames are already known.

    Each donor is placed at ITS OWN distance from the metal (`geometry.distances`), read
    from the donor atom's element.  A centre carrying a water and an iodide has two
    different M-L distances and always did; giving both 2.00 A put the iodide 0.6 A too
    close and every halide co-ligand was then refused by the clash check.  `d_ml`
    overrides all of them, which is what the tests that pin an exact bond length use.

    `cn` is the coordination number of the POLYHEDRON, which may exceed the number of
    donors supplied: the surplus vertices come back as `vacancies`.  It defaults to the
    donor count, so a saturated call is unchanged.  Passing it is how you ask for a
    coordinatively unsaturated centre and get one — a tetrahedral Mg with two ketones and
    two empty vertices, rather than the 2-coordinate linear complex that "just build it
    at the CN you can fill" silently produces instead.

    Single-donor ligands then have their two undetermined rotations chosen by
    `_refine_azimuths` and reported in the choice vector as `azimuth_step` / `oop_step`.
    Handing those integers back on `LigandPlacement` replays the geometry exactly and
    skips the search.  A vacancy is empty space, so a ligand is free to rotate INTO one —
    the search sees only the atoms that are actually there.

    `reserve` names vertices of `site_vectors(geometry, cn)` that no ligand may take, so
    the caller decides WHICH vertices stay open rather than only how many.  Without it the
    assignment fills in its own order and a CN-6 centre carrying four co-ligands comes
    back with its two vacancies *trans*, where a ~90 deg chelate cannot reach them and the
    next step refuses with `chelate_cannot_span`.  A bridge wants the same control for the
    same reason: the vertex facing the partner metal is not available to anything else.
    `cis_vertices` names a mutually-cis set to hand in here.  It is part of the choice
    vector, because it changes which vertex each ligand occupies.
    """
    n_sites = sum(lig.denticity for lig in ligands)
    cn = n_sites if cn is None else int(cn)
    reserved = _check_reserved(reserve, cn)
    if cn < n_sites:
        raise ValueError(
            f"{n_sites} donor site(s) cannot fit a CN-{cn} polyhedron; the ligands ask "
            f"for more vertices than the geometry has")
    if cn - len(reserved) < n_sites:
        raise ValueError(
            f"{n_sites} donor site(s) cannot fit the {cn - len(reserved)} vertex/vertices "
            f"a CN-{cn} polyhedron has left once {sorted(reserved)} are reserved. Reserve "
            f"fewer, or build at a higher coordination number.")
    # The polyhedron chooses DIRECTIONS; the donor chooses how far along one it sits.
    # So assignment happens on unit vectors — it scores angles, which are scale-free —
    # and the scaling is applied per donor afterwards.
    unit_targets = site_vectors(geometry, cn, 1.0)
    per_ligand_units, vacant = _assign_targets(unit_targets, ligands, reserved=reserved)

    per_ligand_d: list[list[Distance]] = []
    for lig in ligands:
        elements = donor_elements(lig.mol, lig.donor_idxs)
        per_ligand_d.append([metal_donor_distance(metal, el, override=d_ml)
                             for el in elements])
    per_ligand_targets = [units * np.array([x.value for x in dists])[:, None]
                          for units, dists in zip(per_ligand_units, per_ligand_d)]
    # Reported as the centre's nominal distance where one number is still wanted (the
    # choice vector's `d_ml`, the M-O base).  It is the base, explicitly, not a mean of
    # whatever happened to be coordinated.
    d = float(d_ml) if d_ml is not None else D_ML.get(metal, DEFAULT_D_ML)

    symbols: list[str] = [metal]
    coords: list[np.ndarray] = [np.zeros(3)]
    bonded: set[tuple[int, int]] = set()
    donor_idxs: list[int] = []
    donor_targets: list[float] = []
    donor_sources: list[str] = []
    owners: list[str] = [f"{metal}(centre)"]
    offsets: list[int] = []
    choices: list[dict] = []
    refinable: list[_Monodentate] = []

    for lig, assigned, dists in zip(ligands, per_ligand_targets, per_ligand_d):
        base = len(symbols)
        offsets.append(base)
        mol, conf = lig.mol, lig.mol.GetConformer()
        xyz = coordinates(mol)
        frames = [site_frame(mol, idx, dtype, conf)
                  for idx, dtype in zip(lig.donor_idxs, lig.donor_types)]
        if lig.denticity == 1:
            frame = frames[0]
            origin, axis, _ref = frame.as_arrays()
            target = assigned[0]
            wells = torsion_wells(lig.donor_types[0], lig.mode)
            well_deg = wells[lig.torsion_well % len(wells)]
            oop_axis = _oop_axis(mol, lig.donor_idxs[0], origin, conf)
            az_grid = _azimuth_offsets(len(wells))
            oop_grid = _oop_angles() if oop_axis is not None else [0.0]
            # Index 0 of both grids is the undisplaced pose, so this first placement is
            # exactly what the placer did before there was a search.  The search itself
            # runs afterwards, once every ligand has somewhere to be.
            ai = 0 if lig.azimuth_step is None else lig.azimuth_step % len(az_grid)
            oi = 0 if lig.oop_step is None else lig.oop_step % len(oop_grid)
            placed = _monodentate_pose(xyz, origin, axis, target, oop_axis,
                                       oop_grid[oi], well_deg + az_grid[ai])
            choice = {"ligand": lig.name, "mode": lig.mode.value,
                      "donors": list(lig.donor_idxs),
                      "torsion_well": lig.torsion_well,
                      "torsion_deg": well_deg,
                      "azimuth_step": ai, "oop_step": oi,
                      "azimuth_deg": round(well_deg + az_grid[ai], 3),
                      "oop_deg": round(oop_grid[oi], 3),
                      "distances": [x.to_dict() for x in dists]}
            choices.append(choice)
            refinable.append(_Monodentate(
                base=base, n_atoms=mol.GetNumAtoms(), xyz=xyz, origin=origin, axis=axis,
                target=target, oop_axis=oop_axis, well_deg=well_deg, az_grid=az_grid,
                oop_grid=oop_grid, az_fixed=lig.azimuth_step, oop_fixed=lig.oop_step,
                donor_atom=base + lig.donor_idxs[0], choice=choice))
        elif lig.denticity == 2:
            # A chelate's bite angle is a property of the LIGAND, not of the idealised
            # polyhedron.  Forcing its two donors onto fixed vertices strains one bond
            # short and the other long — here that showed up as Zn-O spanning 1.85 to
            # 2.39 A when both should be 2.00.  So the vertices only choose the sector;
            # the actual separation comes from the ligand's own donor-donor distance,
            # which makes the two-point alignment exact and both bonds correct.
            p_donors = np.array([f.origin for f in frames])
            d_oo = float(np.linalg.norm(p_donors[1] - p_donors[0]))
            d1, d2 = dists[0].value, dists[1].value
            bisector = _unit(assigned.mean(axis=0))
            normal = np.cross(_unit(assigned[0]), _unit(assigned[1]))
            if float(np.linalg.norm(normal)) < 1e-6:
                normal = _perpendicular(bisector)
            normal = _unit(normal)
            # The bite angle now comes from the law of cosines rather than a symmetric
            # arcsine, because the two donors may sit at DIFFERENT distances (an N,O
            # chelate; an S,O chelate).  With d1 == d2 this reduces exactly to the old
            # `2*asin(d_oo/2d)`, so nothing already built moves.
            cos_bite = (d1 * d1 + d2 * d2 - d_oo * d_oo) / (2.0 * d1 * d2)
            gamma = float(np.arccos(max(-1.0, min(1.0, cos_bite))))
            t1 = d1 * _unit(_axis_rotation(normal, +gamma / 2) @ bisector)
            t2 = d2 * _unit(_axis_rotation(normal, -gamma / 2) @ bisector)
            fitted = np.array([t1, t2])

            # exact two-point alignment (anchored on donor 0, so donor 1 lands on t2 by
            # construction), then the remaining roll about the donor-donor axis is set
            # by pointing the ligand's bulk away from the metal
            rot1 = _rotation_between(p_donors[1] - p_donors[0], fitted[1] - fitted[0])
            moved = (rot1 @ (xyz - p_donors[0]).T).T
            axis = _unit(fitted[1] - fitted[0])
            mid = fitted.mean(axis=0)
            best, best_score = None, -np.inf
            for k in range(72):
                theta = 2 * np.pi * k / 72
                cand = (_axis_rotation(axis, theta) @ moved.T).T + t1
                score = float(np.dot(_unit(cand.mean(axis=0) - mid), _unit(mid)))
                if score > best_score:
                    best, best_score = cand, score
            placed = best
            bite = np.degrees(gamma)
            choices.append({"ligand": lig.name, "mode": lig.mode.value,
                            "donors": list(lig.donor_idxs), "bite_angle_deg": round(bite, 1),
                            "donor_donor_A": round(d_oo, 3),
                            "distances": [x.to_dict() for x in dists]})
        else:
            p_donors = np.array([f.origin for f in frames])
            centroid = xyz.mean(axis=0)
            outward = _unit(assigned.mean(axis=0))
            span = float(np.linalg.norm(centroid - p_donors.mean(axis=0)))
            rot, trans = _kabsch(np.vstack([p_donors, centroid]),
                                 np.vstack([assigned, assigned.mean(axis=0) + outward * span]))
            placed = (rot @ xyz.T).T + trans
            choices.append({"ligand": lig.name, "mode": lig.mode.value,
                            "donors": list(lig.donor_idxs),
                            "distances": [x.to_dict() for x in dists]})

        owner = lig.name or f"ligand{len(offsets)}"
        for atom, position in zip(mol.GetAtoms(), placed):
            symbols.append(atom.GetSymbol())
            coords.append(position)
            owners.append(owner)
        for bond in mol.GetBonds():
            bonded.add((base + bond.GetBeginAtomIdx(), base + bond.GetEndAtomIdx()))
        for idx, dist in zip(lig.donor_idxs, dists):
            donor_idxs.append(base + idx)
            donor_targets.append(dist.value)
            donor_sources.append(dist.source)
            bonded.add((0, base + idx))

    arr = np.array(coords)
    _refine_azimuths(arr, symbols, refinable)
    report = qc(symbols, arr, bonded, metal_idx=0, donor_idxs=donor_idxs,
                d_ml=donor_targets, owners=owners, sources=donor_sources)
    estimated = sorted({x.donor_element for dists in per_ligand_d for x in dists
                        if x.estimated})
    if estimated:
        # Not a failure — a caveat that has to travel with the geometry, so a QC verdict
        # on a pair nobody calibrated is not read as a statement about the chemistry.
        report.notes.append("M-L distance estimated from covalent radii for: "
                            + ", ".join(estimated))
    cv = {"metal": metal, "geometry": geometry, "cn": cn,
          "d_ml": d, "seed": seed, "ligands": choices,
          "donor_distances": [x.to_dict() for dists in per_ligand_d for x in dists]}
    if reserved:
        # Only when something was reserved, so a saturated build's choice vector — and
        # therefore its stored identity — is unchanged by this argument existing.
        cv["reserve"] = sorted(reserved)
    return PlacementResult(symbols, arr, 0, donor_idxs, bonded, report,
                           choice_vector=cv,
                           atom_offsets=offsets, owners=owners,
                           donor_distances=[x for dists in per_ligand_d for x in dists],
                           vacancies=tuple(tuple(float(c) for c in unit_targets[v])
                                           for v in vacant))


def to_rdkit(metal: str, ligands: list[LigandPlacement], result: PlacementResult,
             *, metal_charge: int = 0) -> Chem.Mol:
    """The placed assembly as one RDKit mol, with DATIVE metal-ligand bonds.

    Going back through RDKit rather than building a typed graph by hand means the
    assembled structure is typed by exactly the same code path as a molecule read from
    SMILES — there is one conversion, so there is one place for it to be wrong.
    """
    from rdkit.Geometry import Point3D

    rw = Chem.RWMol()
    m = Chem.Atom(metal)
    m.SetFormalCharge(metal_charge)
    m.SetNoImplicit(True)
    rw.AddAtom(m)

    for lig, base in zip(ligands, result.atom_offsets):
        remap = {}
        for atom in lig.mol.GetAtoms():
            new = Chem.Atom(atom.GetAtomicNum())
            new.SetFormalCharge(atom.GetFormalCharge())
            new.SetNoImplicit(True)
            new.SetNumExplicitHs(0)
            new.SetIsAromatic(atom.GetIsAromatic())
            remap[atom.GetIdx()] = rw.AddAtom(new)
        for bond in lig.mol.GetBonds():
            rw.AddBond(remap[bond.GetBeginAtomIdx()], remap[bond.GetEndAtomIdx()],
                       bond.GetBondType())
        for idx in lig.donor_idxs:
            rw.AddBond(remap[idx], 0, Chem.BondType.DATIVE)

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(result.coords):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                     ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    return mol


@dataclass
class Center:
    """One centre of a polynuclear node — and `element` is not always a metal.

    A BRIDGING ATOM IS A CENTRE whose vertices are metal positions (WORKPLAN_M6 §4).  That
    is the whole of mechanism B: a mu3-oxo has no covalent neighbour and so no lone-pair
    axis to point at a second metal, but `trigonal` vertices at 1.90 A put three metals
    3.291 A apart, which is the literature number.  Nothing here is metal-specific, which
    is why this dataclass can carry both halves of a node.

    `cn` and `local_geometry` have no defaults on purpose (ground rule 5): which
    polyhedron a centre sits in is precisely the ambiguity that has to branch rather than
    be guessed, and a default here would be the guess.  The rest are the per-centre labels
    D12 insists on — a mixed-valence Fe3 node is three centres with different
    `oxidation_state`, not one node with an average.
    """

    element: str
    cn: int
    local_geometry: str
    charge: int = 0
    oxidation_state: int | None = None
    spin_class: str | None = None


def place_multicentre(centers: list[Center], joins: list,
                      constraints: list["InterCentreConstraint"], *, seed: int = 0):
    """Reconcile the centres of a polynuclear node.  NOT IMPLEMENTED.

    Ground rule 7: settled signature, scheduled body (M6, slice S3).

    **This is not where a polynuclear node comes from (D20).**  A node is what a sequence
    of joins produces, and M-M distance is an output to validate, not an input to impose:
    a Cu paddlewheel builds QC-clean from well-1 site frames and `_linalg.kabsch` with no
    solver anywhere.  What this function is for is the case where TWO determinants fix the
    same M...M and disagree — an oxo-centred cluster, where the bridging centre asks for
    3.291 A and the syn-syn carboxylate offers 2.696.  That is one scalar per edge, not a
    general constrained optimisation, and the reconciliation is explicit policy with the
    residual reported as strain, never a hidden average.

    A skeleton that is DETERMINED is built, not searched.  An under-determined set must
    raise `AmbiguousSpecError` naming what is missing rather than falling into a minimiser.

    `joins` is deliberately untyped: what the placer receives from the join path is S3's
    call, and WORKPLAN_M6 §7 lists the `Join` name in this signature as undefined.  It is
    still undefined rather than guessed.
    """
    from mofsbu.assembly.join import NotBuiltYet

    raise NotBuiltYet(
        f"geometry.placer.place_multicentre ({len(centers)} centres) — M6. "
        "Use place_mononuclear for a single centre."
    )


@dataclass
class InterCentreConstraint:
    """Distance / angle relationship between two coordination centres (M6)."""

    centres: tuple[int, int]
    mm_distance: float | None = None
    #: The acceptance window, when one has been curated for this node
    #: (`data/reference/node_cases.tsv`).  Absent means nobody has measured it, which is
    #: not the same as a zero-width window — `window()` falls back to a tolerance rather
    #: than refusing everything.
    mm_lo: float | None = None
    mm_hi: float | None = None
    bridge_bite_deg: float | None = None
    #: Does this pair carry a METAL_METAL edge?  **C10 is open** and this field is the
    #: "declare" half of its leaning (branch or declare, never default).  Not inferable
    #: from `mm_distance`: a paddlewheel's Cu-Cu at 2.62 A is an edge and an Fe3-oxo
    #: trimer's Fe...Fe at 3.29 A is not, and the two hand-written fixtures differ at L1 by
    #: exactly that edge.  `EdgeType.METAL_METAL` is in the certificate, so a placer that
    #: guessed it from a distance threshold would be guessing the identity of its product.
    metal_metal_bond: bool = False

    def window(self, *, tol: float = 0.15) -> tuple[float, float] | None:
        """The range this pair may sit in, or None if the constraint sets no distance."""
        if self.mm_lo is not None and self.mm_hi is not None:
            return (self.mm_lo, self.mm_hi)
        if self.mm_distance is None:
            return None
        return (self.mm_distance - tol, self.mm_distance + tol)
