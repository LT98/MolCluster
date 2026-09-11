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

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from mofsbu.geometry.distances import (
    BASE_MO, DEFAULT_BASE_MO, Distance, donor_elements, metal_donor_distance,
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
    "linear": {2: "linear"}, "trigonal": {3: "trigonal"},
    "tetrahedral": {4: "tetrahedral"}, "square_planar": {4: "square_planar"},
    "trigonal_bipyramidal": {5: "trigonal_bipyramidal"},
    "square_pyramidal": {5: "square_pyramidal"}, "octahedral": {6: "octahedral"},
}


def site_vectors(geometry: str, n: int, d: float = 2.05) -> np.ndarray:
    """Unit coordination directions for a coordination number, scaled to `d`.

    The legacy table mapped CN 5 to 'planar', which is not a coordination geometry;
    both real CN-5 polyhedra are here instead.
    """
    g = geometry.lower()
    s3 = 1 / np.sqrt(3)
    table = {
        "linear": np.array([[1, 0, 0], [-1, 0, 0]], float),
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


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])


def _perpendicular(v: np.ndarray) -> np.ndarray:
    trial = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    return _unit(np.cross(v, trial))


def _rotation_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector `a` onto unit vector `b`."""
    a, b = _unit(a), _unit(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if float(np.linalg.norm(v)) < 1e-9:
        if c > 0:
            return np.eye(3)
        trial = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = _unit(np.cross(a, trial))
        return _axis_rotation(axis, np.pi)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1 / (1 + c))


def _axis_rotation(axis: np.ndarray, theta: float) -> np.ndarray:
    axis = _unit(axis)
    x, y, z = axis
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ])


def _kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation + translation least-squares fitting point set `p` onto `q`."""
    cp, cq = p.mean(axis=0), q.mean(axis=0)
    h = (p - cp).T @ (q - cq)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, cq - rot @ cp


# A five-membered chelate ring (M-D-C-C-D) subtends roughly this at the metal.  A
# polydentate ligand handed two trans vertices cannot close its ring and folds into
# itself, which is a placement bug that shows up only as a pile of clashes downstream.
IDEAL_BITE_ANGLE = 82.0
BITE_ANGLE_RANGE = (55.0, 115.0)


def _angle_between(a: np.ndarray, b: np.ndarray) -> float:
    cos = float(np.dot(_unit(a), _unit(b)))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))


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


def _assign_targets(targets: np.ndarray, ligands: list) -> tuple[list[np.ndarray], list[int]]:
    """Give each ligand coordination vertices its denticity can actually reach.

    Polydentate ligands are served first and take the vertex set whose angles sit
    closest to a chelate bite angle; monodentate ligands take whatever is left.  Without
    this the vertices are handed out in table order, and on an octahedron that puts a
    bidentate ligand across a trans pair.

    Returns the per-ligand vertex sets AND the indices of the vertices nothing claimed.
    That second value used to be dropped on the floor, which is why an unsaturated centre
    had to be rebuilt at a smaller coordination number instead of simply having empty
    vertices: the information that they existed was computed and then discarded here.
    """
    from itertools import combinations

    remaining = list(range(len(targets)))
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
    return [assigned[i] for i in range(len(ligands))], remaining


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
    """
    n_sites = sum(lig.denticity for lig in ligands)
    cn = n_sites if cn is None else int(cn)
    if cn < n_sites:
        raise ValueError(
            f"{n_sites} donor site(s) cannot fit a CN-{cn} polyhedron; the ligands ask "
            f"for more vertices than the geometry has")
    # The polyhedron chooses DIRECTIONS; the donor chooses how far along one it sits.
    # So assignment happens on unit vectors — it scores angles, which are scale-free —
    # and the scaling is applied per donor afterwards.
    unit_targets = site_vectors(geometry, cn, 1.0)
    per_ligand_units, vacant = _assign_targets(unit_targets, ligands)

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
    return PlacementResult(symbols, arr, 0, donor_idxs, bonded, report,
                           choice_vector={"metal": metal, "geometry": geometry,
                                          "cn": cn,
                                          "d_ml": d, "seed": seed, "ligands": choices,
                                          "donor_distances": [
                                              x.to_dict() for dists in per_ligand_d
                                              for x in dists]},
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


def place_multicentre(centers: list[Center], joins: list[Join],
                      constraints: list["InterCentreConstraint"], *, seed: int = 0):
    """Place several coordination centres with inter-centre constraints.  NOT IMPLEMENTED.

    Ground rule 7: settled signature, scheduled body (M6, the plan's declared headline
    cost).  `place_mononuclear` above is the N=1 path of exactly this operation — there is
    no separate mononuclear code path to reconcile, which is the point of D12.

    Needs: M-M distance and bridge bite-angle constraints solved jointly with each centre's
    local geometry, plus `geometry.qc.check_intercentre`.
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
    bridge_bite_deg: float | None = None
