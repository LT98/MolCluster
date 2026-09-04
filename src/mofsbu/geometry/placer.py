"""Placing ligands around a coordination centre.

This is the **N = 1 path of the general operation** (D12): one centre is not a special
case, it is the smallest case.  The multi-centre placer of M6 adds inter-centre
constraints on top of exactly this machinery; nothing here branches on the number of
metals.

Placement is deterministic, driven by the stored site frames rather than a stochastic
search: a ligand is rotated so its donor's outward axis points at the metal, then rolled
to a chosen torsion well.  That is what makes a built structure reproducible from its
choice-vector instead of only from its coordinates.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem

from mofsbu.geometry.embed import coordinates
from mofsbu.geometry.qc import QCReport, qc
from mofsbu.sites.frames import BindingMode, SiteFrame, site_frame, torsion_wells

# Ideal metal-donor distances, angstrom.  Coarse, per metal; refined by relaxation.
D_ML: dict[str, float] = {"Zn": 2.00, "Cu": 1.98, "Ni": 2.06, "Co": 2.08, "Fe": 2.10,
                          "Mn": 2.18, "Cr": 2.00, "Mg": 2.10, "Ca": 2.40, "Cd": 2.28}
DEFAULT_D_ML = 2.05

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
    return d * vecs


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


def _assign_targets(targets: np.ndarray, ligands: list) -> list[np.ndarray]:
    """Give each ligand coordination vertices its denticity can actually reach.

    Polydentate ligands are served first and take the vertex set whose angles sit
    closest to a chelate bite angle; monodentate ligands take whatever is left.  Without
    this the vertices are handed out in table order, and on an octahedron that puts a
    bidentate ligand across a trans pair.
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
    return [assigned[i] for i in range(len(ligands))]


@dataclass
class LigandPlacement:
    """One ligand, which of its donors bind, and how."""

    mol: Chem.Mol
    donor_idxs: tuple[int, ...]
    donor_types: tuple[str, ...]
    mode: BindingMode = BindingMode.MONODENTATE
    torsion_well: int = 0
    name: str = ""

    @property
    def denticity(self) -> int:
        return len(self.donor_idxs)


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

    @property
    def ok(self) -> bool:
        return self.report.ok

    def to_xyz(self, comment: str = "") -> str:
        lines = [str(len(self.symbols)), comment]
        for sym, (x, y, z) in zip(self.symbols, self.coords):
            lines.append(f"{sym:<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(lines) + "\n"


def place_mononuclear(
    metal: str,
    ligands: list[LigandPlacement],
    *,
    geometry: str,
    d_ml: float | None = None,
    seed: int = 0,
) -> PlacementResult:
    """Assemble one coordination centre from ligands whose frames are already known."""
    d = d_ml if d_ml is not None else D_ML.get(metal, DEFAULT_D_ML)
    n_sites = sum(lig.denticity for lig in ligands)
    targets = site_vectors(geometry, n_sites, d)

    per_ligand_targets = _assign_targets(targets, ligands)

    symbols: list[str] = [metal]
    coords: list[np.ndarray] = [np.zeros(3)]
    bonded: set[tuple[int, int]] = set()
    donor_idxs: list[int] = []
    offsets: list[int] = []
    choices: list[dict] = []

    for lig, assigned in zip(ligands, per_ligand_targets):
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
            # the donor's outward axis must point at the metal, which sits at the origin
            rot = _rotation_between(axis, _unit(-target))
            wells = torsion_wells(lig.donor_types[0], lig.mode)
            theta = np.radians(wells[lig.torsion_well % len(wells)])
            rot = _axis_rotation(_unit(-target), theta) @ rot
            placed = (rot @ (xyz - origin).T).T + target
            choices.append({"ligand": lig.name, "mode": lig.mode.value,
                            "donors": list(lig.donor_idxs),
                            "torsion_well": lig.torsion_well,
                            "torsion_deg": wells[lig.torsion_well % len(wells)]})
        elif lig.denticity == 2:
            # A chelate's bite angle is a property of the LIGAND, not of the idealised
            # polyhedron.  Forcing its two donors onto fixed vertices strains one bond
            # short and the other long — here that showed up as Zn-O spanning 1.85 to
            # 2.39 A when both should be 2.00.  So the vertices only choose the sector;
            # the actual separation comes from the ligand's own donor-donor distance,
            # which makes the two-point alignment exact and both bonds correct.
            p_donors = np.array([f.origin for f in frames])
            d_oo = float(np.linalg.norm(p_donors[1] - p_donors[0]))
            bisector = _unit(assigned.mean(axis=0))
            normal = np.cross(_unit(assigned[0]), _unit(assigned[1]))
            if float(np.linalg.norm(normal)) < 1e-6:
                normal = _perpendicular(bisector)
            normal = _unit(normal)
            half = float(np.arcsin(min(1.0, d_oo / (2.0 * d))))
            t1 = d * _unit(_axis_rotation(normal, +half) @ bisector)
            t2 = d * _unit(_axis_rotation(normal, -half) @ bisector)
            fitted = np.array([t1, t2])

            # exact two-point alignment, then the remaining roll about the donor-donor
            # axis is set by pointing the ligand's bulk away from the metal
            rot1 = _rotation_between(p_donors[1] - p_donors[0], fitted[1] - fitted[0])
            moved = (rot1 @ (xyz - p_donors.mean(axis=0)).T).T
            axis = _unit(fitted[1] - fitted[0])
            mid = fitted.mean(axis=0)
            best, best_score = None, -np.inf
            for k in range(72):
                theta = 2 * np.pi * k / 72
                cand = (_axis_rotation(axis, theta) @ moved.T).T + mid
                score = float(np.dot(_unit(cand.mean(axis=0) - mid), _unit(mid)))
                if score > best_score:
                    best, best_score = cand, score
            placed = best
            bite = np.degrees(2 * half)
            choices.append({"ligand": lig.name, "mode": lig.mode.value,
                            "donors": list(lig.donor_idxs), "bite_angle_deg": round(bite, 1),
                            "donor_donor_A": round(d_oo, 3)})
        else:
            p_donors = np.array([f.origin for f in frames])
            centroid = xyz.mean(axis=0)
            outward = _unit(assigned.mean(axis=0))
            span = float(np.linalg.norm(centroid - p_donors.mean(axis=0)))
            rot, trans = _kabsch(np.vstack([p_donors, centroid]),
                                 np.vstack([assigned, assigned.mean(axis=0) + outward * span]))
            placed = (rot @ xyz.T).T + trans
            choices.append({"ligand": lig.name, "mode": lig.mode.value,
                            "donors": list(lig.donor_idxs)})

        for atom, position in zip(mol.GetAtoms(), placed):
            symbols.append(atom.GetSymbol())
            coords.append(position)
        for bond in mol.GetBonds():
            bonded.add((base + bond.GetBeginAtomIdx(), base + bond.GetEndAtomIdx()))
        for idx in lig.donor_idxs:
            donor_idxs.append(base + idx)
            bonded.add((0, base + idx))

    arr = np.array(coords)
    report = qc(symbols, arr, bonded, metal_idx=0, donor_idxs=donor_idxs, d_ml=d)
    return PlacementResult(symbols, arr, 0, donor_idxs, bonded, report,
                           choice_vector={"metal": metal, "geometry": geometry,
                                          "d_ml": d, "seed": seed, "ligands": choices},
                           atom_offsets=offsets)


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
