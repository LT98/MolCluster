"""L2: the configurational isomer tag — cis/trans, fac/mer, Delta/Lambda (D10, M5/S5).

L1 asks what is bonded to what; L2 asks how those bonds are *arranged in space*, and that
is a question the graph cannot answer.  cis- and trans-Pt(NH3)2Cl2 have one connectivity
and one L1 hash, which `test_l1_does_not_separate_cis_from_trans` exists to keep true — a
finer L1 would also split things that are genuinely one structure.  So the split happens
here, and it needs coordinates.

**No geometry, no tag.**  `l2_isomer_tag(g)` with no coordinates returns `""` and that is
the correct answer, not a degraded one: nothing about the graph distinguishes the isomers,
so claiming a tag would be inventing one.  It has a consequence worth stating plainly —
**L2 is not a pure function of a structure row.**  A structure ingested without coordinates
legitimately carries `""` and one built from a placement carries a real tag, and those are
two rows under `UNIQUE (l0, l1, l2)`.  That is the same shape as the `site_catalog` seam
(§6 of the architecture map) and it is why `put_structure` takes `l2=` as a parameter: the
builder has coordinates in hand at insert time and the registry does not.

**What gets compared.**  Donors are grouped into *classes* by the certificate of the ligand
fragment they belong to, so "the two ammines" is a fact about the ligands rather than about
the element.  Each fragment INSTANCE in a class is then reduced to one direction — the mean
of its donors' unit vectors from the metal — and the arrangement of those instance
directions is what gets classified.  Reducing to the instance is what makes the same rule
serve a bis-chelate as well as MA2B2: two monodentate ammines and two bidentate diamines
are both "two instances of one class", and asking whether they sit cis or trans is the same
question in both cases.  A class with only one instance is skipped, because a chelate's own
bite angle is not isomerism — it is a ligand binding the way that ligand binds.
"""
from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np

from mofsbu.graph._types import EdgeType, TypedGraph
from mofsbu.graph.canon import canonical_index_map, certificate_digest

#: Above this angle two donors are trans, below it they are cis.  The midpoint of the two
#: ideal values (90 and 180), which is the only defensible place for it while the
#: distinction is barrier-separated (D10): real cis and trans geometries sit 45 degrees
#: either side, so nothing hinges on the exact number and everything would hinge on it if
#: the two populations ever met in the middle.  A distorted structure that lands near 135
#: is a structure whose isomer nobody should be reading off a threshold, which is why
#: `AMBIGUOUS` exists rather than a silent rounding.
CIS_TRANS_SPLIT_DEG = 135.0

#: How far from the split an angle must sit before the call is made at all.  Inside this
#: band the arrangement is reported as ambiguous and the centre goes untagged, so a bad
#: geometry produces a missing answer instead of a confident wrong one.
AMBIGUITY_BAND_DEG = 15.0

AMBIGUOUS = "?"


def _coords(geom: Any, n_atoms: int) -> np.ndarray | None:
    """Coordinates as (n_atoms, 3) in graph node order, or None if there are none."""
    if geom is None:
        return None
    geom = getattr(geom, "coords", geom)
    arr = np.asarray(geom, dtype=float)
    if arr.shape != (n_atoms, 3):
        raise ValueError(
            f"geometry is {arr.shape} but the graph has {n_atoms} atoms. Row k must be "
            f"the atom at graph.nodes()[k]; a geometry in another order would read the "
            f"isomer off the wrong atoms and label it with confidence.")
    return arr


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])


def _angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    cos = float(np.dot(_unit(a), _unit(b)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def _relation(angle: float) -> str:
    """cis, trans, or an admission that this geometry does not say."""
    if abs(angle - CIS_TRANS_SPLIT_DEG) < AMBIGUITY_BAND_DEG:
        return AMBIGUOUS
    return "trans" if angle > CIS_TRANS_SPLIT_DEG else "cis"


def _fragment_of(g: TypedGraph) -> dict[int, tuple[int, ...]]:
    """{atom: the ligand fragment it belongs to}."""
    out: dict[int, tuple[int, ...]] = {}
    for fragment in g.ligand_fragments():
        for atom in fragment:
            out[atom] = fragment
    return out


def _class_label(fragment: tuple[int, ...], donor_element: str) -> str:
    """A readable, order-invariant name for "this kind of ligand arm".

    Element plus fragment size, which reads as `N4` for an ammine and `Cl1` for a chloride.
    Two different ligands can collide on that label, so the CLASS itself is keyed on the
    fragment's certificate and the label is only how the class is spelled; collisions are
    broken by a rank suffix at the point the tag is assembled.
    """
    return f"{donor_element}{len(fragment)}"


def _chirality(axis: np.ndarray, pairs: Sequence[tuple[np.ndarray, np.ndarray]]) -> str:
    """Handedness of a tris-chelate propeller: `Delta` (right-handed) or `Lambda`.

    IUPAC defines Delta as the RIGHT-handed helix.  A right-handed helix — take
    `(cos t, sin t, ct)` with `c > 0` — advances along its axis while turning
    anticlockwise about that same axis by the right-hand rule.  So orient the axis from
    the lower donor triangle to the upper one, and for each chelate ask which way it turns
    going from its lower donor to its upper one: anticlockwise about the axis is
    right-handed, which is Delta.

    The orientation of the axis is arbitrary and the answer does not depend on it, which is
    the check that this is measuring chirality rather than a point of view: flipping the
    axis swaps which donor of each pair is "upper" AND negates the axis, and the two sign
    changes cancel.  `test_the_handedness_does_not_depend_on_which_end_you_look_from` pins
    that.

    Returns `AMBIGUOUS` when the three chelates disagree, which means the centre is not
    the propeller this function assumes and no label should be attached to it.
    """
    signs = []
    for lower, upper in pairs:
        lo = lower - np.dot(lower, axis) * axis        # project out the axial component
        up = upper - np.dot(upper, axis) * axis
        if float(np.linalg.norm(lo)) < 1e-6 or float(np.linalg.norm(up)) < 1e-6:
            return AMBIGUOUS                            # a donor sits on the axis itself
        signs.append(float(np.dot(axis, np.cross(_unit(lo), _unit(up)))))
    if all(s > 1e-6 for s in signs):
        return "Delta"
    if all(s < -1e-6 for s in signs):
        return "Lambda"
    return AMBIGUOUS


def _centre_tag(g: TypedGraph, coords: np.ndarray, metal: int) -> str:
    """The arrangement tag for one coordination centre, or `''` if it has no distinction."""
    donors = g.neighbors(metal, EdgeType.DATIVE)
    if len(donors) < 2:
        return ""
    origin = coords[metal]
    fragments = _fragment_of(g)

    # class key -> {fragment instance -> [donor unit vectors]}
    classes: dict[tuple[str, str], dict[tuple[int, ...], list[np.ndarray]]] = {}
    for donor in donors:
        fragment = fragments.get(donor, (donor,))
        cert = certificate_digest(g.subgraph(fragment))
        key = (cert, _class_label(fragment, g.label(donor).element))
        classes.setdefault(key, {}).setdefault(fragment, []).append(
            _unit(coords[donor] - origin))

    labels: dict[str, list[tuple[str, str]]] = {}
    for (cert, label), instances in classes.items():
        if len(instances) < 2:
            # One instance: a chelate biting the way it bites, not an isomeric choice.
            continue
        arrangement = _arrangement(instances)
        if arrangement:
            labels.setdefault(label, []).append((cert, arrangement))

    parts: list[str] = []
    for label in sorted(labels):
        entries = sorted(labels[label])                  # by certificate: order-invariant
        for rank, (_cert, arrangement) in enumerate(entries):
            # A rank suffix only when two different ligands spell their label the same.
            spelled = label if len(entries) == 1 else f"{label}#{rank}"
            parts.append(f"{spelled}={arrangement}")
    return ",".join(parts)


def _arrangement(instances: dict[tuple[int, ...], list[np.ndarray]]) -> str:
    """cis/trans, fac/mer, or Delta/Lambda for the instances of one donor class."""
    ordered = sorted(instances.items())
    directions = [_unit(np.sum(vectors, axis=0)) for _frag, vectors in ordered]
    chelating = [len(vectors) >= 2 for _frag, vectors in ordered]

    if len(directions) == 2:
        return _relation(_angle_deg(directions[0], directions[1]))

    if len(directions) == 3:
        if all(chelating):
            # Three bidentate arms on one centre: the distinction is handedness, not
            # fac/mer.
            #
            # The C3 axis is the NORMAL TO THE PLANE the three arm directions lie in, and
            # emphatically not their sum: the six donors of a fully coordinated octahedron
            # sum to zero, so the arm directions do too, and "the sum of the arms" hands
            # back a zero vector for every tris-chelate there is.  The three arms of a
            # propeller sit at 120 degrees in the plane perpendicular to the axis, which is
            # what makes the normal the right construction and the sum the wrong one.
            axis = np.cross(directions[1] - directions[0], directions[2] - directions[0])
            if float(np.linalg.norm(axis)) < 1e-6:
                return AMBIGUOUS              # the arms are collinear: not a propeller
            axis = _unit(axis)
            pairs = []
            for _frag, vectors in ordered:
                by_axis = sorted(vectors, key=lambda v: float(np.dot(v, axis)))
                pairs.append((by_axis[0], by_axis[-1]))
            return _chirality(axis, pairs)
        relations = [_relation(_angle_deg(a, b))
                     for i, a in enumerate(directions) for b in directions[i + 1:]]
        if AMBIGUOUS in relations:
            return AMBIGUOUS
        n_trans = relations.count("trans")
        if n_trans == 0:
            return "fac"                      # all mutually cis: one triangular face
        if n_trans == 1:
            return "mer"                      # one trans pair: a meridian
        return AMBIGUOUS                      # not a recognisable CN-6 arrangement
    return ""                                 # 1, or 4+: no two-way distinction to draw


def isomer_tag(g: TypedGraph, geom: Any | None = None) -> str:
    """The L2 tag for a whole structure: every centre that has an arrangement to report.

    Centres are keyed by CANONICAL index, not build index, so shuffling the input atom
    order cannot move the tag — the same invariance L1 is already held to, and the reason
    a tag can be part of an identity key at all.
    """
    coords = _coords(geom, len(g))
    if coords is None:
        return ""
    canon = canonical_index_map(g)
    parts = []
    for metal in g.metals():
        tag = _centre_tag(g, coords, metal)
        if tag and AMBIGUOUS not in tag:
            parts.append((canon[metal], tag))
    return ";".join(f"{c}:{t}" for c, t in sorted(parts))
