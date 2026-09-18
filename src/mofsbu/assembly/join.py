"""Joining building blocks — the recursive assembly operation (D1, D13).  PARTLY BUILT.

Ground rule 7: the signatures are the design decision and they are settled here, so
callers can be written against them now.  The remaining bodies are scheduled work (M5/S3
onward).  Every one of those raises; none returns a plausible-looking value, because a
stub that quietly answers is indistinguishable from a working function until the answers
are wrong.

Where the milestone stands:

* `sites.frames` — DONE.  A join needs both partners' frames and a torsion well, and
  those are stored rather than re-derived.
* `assembly.choice` — DONE (S1).  The vector a join has to emit, and its key.
* `compatible()` / `chelate_compatible()` — DONE (S2).  Frame-alignment feasibility
  between two open sites, and the same question `sites.model.chelate_pockets` answers
  within one molecule, generalised to two blocks.
* `join()` — align, bond, recompute the product's open sites by inheriting through the
  atom map, and emit the choice-vector that regenerates it.  DONE (S3).
* `join_chelate()` / `chelate_reach()` — the two-point version: one ligand across two
  vertices of one centre.  DONE.  The verdict it places under is distance geometry (can a
  rigid ligand reach both vertices at their own bond lengths) rather than the pocket's
  frame-implied bite, and the difference matters — see `chelate_reach`.
* `geometry.placer.place_multicentre` — RECONCILIATION, for the edges where two
  determinants fix one M...M and disagree (D20).  Not a prerequisite for a polynuclear
  product: a one-contact join is a rigid move however many metals either block carries,
  so a bridge EMERGES from a sequence of them and its M...M is an output to validate.
* `join_bridge()` — one donor across vertices of DIFFERENT metals in one move.  M6/S1,
  not built; `chelate_compatible` refuses that case by name and that refusal is where it
  begins.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Sequence

import numpy as np

from mofsbu._types import AmbiguousSpecError, AtomMap, MofsbuError, NotBuiltYet
from mofsbu.energy.backends import combined_multiplicity
from mofsbu.graph._types import EdgeType, TypedGraph
from mofsbu.geometry.distances import metal_donor_distance
from mofsbu.geometry._linalg import (
    angle_between as _angle_deg, axis_rotation as _axis_rotation,
    rotation_between as _rotation_between, unit as _unit,
)
from mofsbu.sites.frames import BindingMode, torsion_wells
from mofsbu.sites.inherit import inherit_sites, merge_inherited
from mofsbu.sites.model import Site, frame_lobes
from mofsbu.sites.state import refresh_state


#: Re-exported: it moved to `_types` (where the rest of the hierarchy lives) when this
#: module grew real dependencies on `sites` and `descriptors`, which had been importing
#: it from here.  Every existing `from mofsbu.assembly.join import NotBuiltYet` still
#: works, and none of them has to care that the definition moved.
__all__ = ["BuildingBlock", "Compatibility", "IncompatibleJoin", "JoinResult",
           "NotBuiltYet", "chelate_compatible", "chelate_reach", "compatible", "grow",
           "join", "join_chelate"]


@dataclass(frozen=True)
class BuildingBlock:
    """A molecule, a metal, or an assembled fragment — one interface for all three (D1)."""

    graph: TypedGraph
    sites: tuple[Site, ...]
    geometry: Any | None = None
    structure_id: int | None = None
    #: Per-geometry state from `sites.state.refresh_state`, keyed by `(atom_idx, slot)`.
    #: The slot is in the key because a metal carries several vacancies on ONE atom, so
    #: an atom-keyed dict silently keeps whichever of them was inserted last.  Absent
    #: means "not computed for this block", which is not the same as "nothing is open".
    state: dict[tuple[int, int], Any] | None = None

    @staticmethod
    def state_key(site: Any) -> tuple[int, int]:
        """The key a site's state is stored under.  One definition, used by both sides."""
        return (site.atom_idx, getattr(site, "slot", 0))

    def open_sites(self) -> tuple[Site, ...]:
        """The sites a join may actually use — donors AND vacant metal vertices.

        Open means: not already dative-bonded to a metal, and not sterically walled off
        (`sites.state.SiteStatus`).  Both halves matter — "both ends are unoccupied" is
        the test §6.3 explicitly says is not sufficient, and a site the metal cannot
        reach is not a site a join can use however free its valence looks.

        The two kinds come back in one tuple on purpose.  A join needs a donor on one
        block and somewhere on the other block to put it, and on a metal that somewhere is
        a vacancy — so `compatible(a, b)` can take two `Site`s and ask one frame-alignment
        question, rather than needing a separate accessor and a separate predicate for the
        metal's side of every bond.

        A block with no state raises rather than returning every site.  Treating unknown
        as open is how an assembly step would confidently join onto a buried donor.
        """
        if self.state is None:
            raise NotBuiltYet(
                "this BuildingBlock carries no site state, so 'open' is unknown. Run "
                "sites.state.refresh_state on its geometry (or load it with "
                "registry.get_site_state) — returning every site would silently treat "
                "occupied and buried donors as available.")
        from mofsbu.sites.state import SiteStatus

        return tuple(s for s in self.sites
                     if getattr(self.state.get(self.state_key(s)), "status", None)
                     is SiteStatus.OPEN)

    def open_vacancies(self) -> tuple[Site, ...]:
        """Just the metal's usable empty vertices — where an incoming ligand can go."""
        return tuple(s for s in self.open_sites() if s.is_vacancy)

    def open_donors(self) -> tuple[Site, ...]:
        """Just the usable donor atoms — what this block can offer another one."""
        return tuple(s for s in self.open_sites() if not s.is_vacancy)


@dataclass(frozen=True)
class Compatibility:
    """The verdict on one proposed bond.

    `feasible` / `strain` / `reason` are the settled three.  The rest are additive and
    defaulted: a caller written against the original three is unaffected, and `join` needs
    them — it has to know which mode was agreed, how many torsion wells the bond branches
    over (that count IS the number of L3 siblings it will generate, D11), and at what
    distance to place the donor.

    `strain` is unitless and monotone by construction: it is the mismatch divided by the
    tolerance, so `feasible` is exactly `strain <= 1` for every geometric verdict.  That
    keeps the one policy number in a named constant rather than spread across comparisons.
    """

    feasible: bool
    strain: float
    reason: str
    mode: str = BindingMode.MONODENTATE.value
    wells: tuple[float, ...] = ()
    d_ml: float | None = None
    #: How many in-plane lobes the donor offers — the size of the OTHER discrete branch a
    #: join chooses from, alongside the torsion well.  One for a donor whose direction its
    #: bonding fixes; two for an sp2 donor, and for those two the choice is worth 2.8 A of
    #: M...M (the difference between a syn-syn bridge and an anti-anti one).
    lone_pairs: int = 1


#: How far a ligand's bite angle may sit from the separation of the two vertices it is
#: asked to span.  The two populations it separates are not close together.  Against an
#: octahedral CIS pair (90 deg): acetate measures 59.7 deg here, a 30.3 deg mismatch, and
#: it is the STRAINED end of what really forms — the common chelators (acac ~92, en ~85,
#: bipy ~78) sit within 12 deg.  Against a TRANS pair (180 deg) every one of them misses
#: by 88 or more; acetate by 120.  So anything from ~35 to ~60 draws the same line, and 40
#: is taken from the low half so that the strained-but-real case passes with margin while
#: nothing comes near trans.  `test_the_bite_angle_populations_stay_far_apart` pins the
#: separation, so a donor type whose geometry moves cannot quietly cross it.
MAX_BITE_MISMATCH_DEG = 40.0

#: Nominal metal-donor distance used to locate the implied metal when the caller has not
#: named a partner.  Same default, for the same reason, as `sites.model.chelate_pockets`:
#: the bite ANGLE is what the verdict turns on, and it moves very little with this number.
NOMINAL_D_ML = 2.0


@dataclass(frozen=True)
class JoinResult:
    """The product, plus everything needed to explain and reproduce it.

    `atom_map` is the FIRST argument's map — `a`'s `{parent_idx: child_idx}` — because `a`
    is the block being grown and provenance follows it. `partner_atom_map` is `b`'s, and it
    is a separate field rather than a merged one because the two index spaces are
    different and merging them is precisely the mistake `sites.inherit` is written to
    prevent. Additive and defaulted: a caller written against the settled four is
    unaffected.
    """

    block: BuildingBlock
    atom_map: AtomMap
    choice_vector: dict[str, Any]
    strain: float
    partner_atom_map: AtomMap = field(default_factory=dict)
    compatibility: "Compatibility | None" = None


class IncompatibleJoin(MofsbuError):
    """A join the site pair does not permit.  Carries the verdict that refused it."""

    def __init__(self, verdict: "Compatibility") -> None:
        super().__init__(verdict.reason)
        self.verdict = verdict


def _lone_pairs(site: Site) -> int:
    """How many directions this donor offers.  One for anything its bonding determines.

    The count IS the answer to "can this atom bridge on its own" (`lone_pair_frames`),
    and it is the size of the Kind-B branch a join chooses from.
    """
    return max(len(frame_lobes(site.frame)), 1)


def _frame_of(site: Site, lone_pair: int = 0
              ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """`(origin, axis, ref)` from a site's stored frame, or None if it has none.

    `lone_pair` selects among an sp2 donor's in-plane lobes; 0 is the stored `frame`
    itself, so every existing caller reads exactly what it always read.  Which lobe binds
    is a discrete choice with no default answer in the chemistry — it is the difference
    between a syn-syn bridge and an anti-anti one — so it is an index into a set, the same
    shape as a torsion well (D13).
    """
    lobes = frame_lobes(site.frame)
    frame = lobes[lone_pair % len(lobes)] if lobes else None
    if not frame or "axis" not in frame or "origin" not in frame:
        return None
    return (np.asarray(frame["origin"], dtype=float),
            _unit(np.asarray(frame["axis"], dtype=float)),
            _unit(np.asarray(frame.get("ref", [0.0, 0.0, 1.0]), dtype=float)))


def _triad(axis: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """An orthonormal frame as a matrix whose COLUMNS are its axes.

    `ref` is re-orthogonalised against `axis` rather than trusted: a stored frame is
    perpendicular by construction, but it has been through a transform and a JSON round
    trip, and a `ref` that has drifted a degree off would tilt every join by a degree.
    """
    a = _unit(axis)
    r = ref - np.dot(ref, a) * a
    if float(np.linalg.norm(r)) < 1e-9:
        r = np.cross(a, [1.0, 0.0, 0.0] if abs(a[0]) < 0.9 else [0.0, 1.0, 0.0])
    r = _unit(r)
    return np.column_stack([a, r, np.cross(a, r)])


def _rotate_about(v: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    return _axis_rotation(axis, theta) @ v


def _metal_symbol(partner: Any | None) -> str | None:
    """The partner a join bonds to, as an element symbol, however it was handed over."""
    if partner is None:
        return None
    if isinstance(partner, str):
        return partner
    for attr in ("element", "symbol"):
        value = getattr(partner, attr, None)
        if isinstance(value, str):
            return value
    return None


def _roles(a: Site, b: Site) -> tuple[Site, Site] | None:
    """`(donor, vacancy)` in that order, or None if this pair is not one of each."""
    if a.is_vacancy and not b.is_vacancy:
        return b, a
    if b.is_vacancy and not a.is_vacancy:
        return a, b
    return None


def compatible(a: Site, b: Site, *, partner: Any | None = None,
               mode: str = BindingMode.MONODENTATE.value,
               donor_element: str | None = None) -> Compatibility:
    """Can these two open sites be joined at acceptable strain?

    Frame-alignment feasibility, not merely "both are open" — and *openness is not
    rechecked here*.  `BuildingBlock.open_sites` already refuses to guess at it, and a
    predicate that silently re-answered a question its caller has already answered is how
    two definitions of "open" start to drift.  Pass sites that are open.

    **What a single-point join can and cannot be strained by.**  The two sites live in two
    different blocks, so the rigid-body transform that brings them together is free: for
    one donor onto one vacancy there is always a placement realising the alignment exactly,
    and any "strain" reported for it would be a number about the blocks' current positions
    in their own coordinate systems, which mean nothing to each other.  So this returns
    `strain = 0.0` for a monodentate join and says why, rather than manufacturing a score.
    The geometric content of compatibility appears the moment a join has to satisfy TWO
    constraints at once — that is `chelate_compatible` below, and it is where the same
    convergence question `chelate_pockets` answers within one molecule becomes an
    inter-block one.

    What can still refuse a single-point join, and does:

    * **Roles.**  A join puts a donor on a vacant coordination vertex.  Donor-to-donor and
      vacancy-to-vacancy are not bonds this model forms, and both are reported by name
      rather than as a bare `False`.
    * **Mode.**  The requested mode must be offered by both sites.  A vacancy offers
      `mono` only — a chelate needs two vertices, which is a different call.

    `partner` is the metal being bonded to (a symbol, or anything carrying `.element` /
    `.symbol`).  With `donor_element` it yields the pair's M-D distance, which `join` needs
    to place the donor.  The element is taken from the caller, never parsed out of
    `donor_type`: `geometry.distances` is explicit that `halide_I` and `hydrohalide_X` name
    the same iodine and that inferring an element from a type name is how a new donor type
    silently inherits oxygen's bond length.
    """
    pair = _roles(a, b)
    if pair is None:
        kind = "two vacancies" if a.is_vacancy else "two donors"
        what = ("a bond between two coordination vertices is a metal-metal bond, which is "
                "not formed by joining sites (D12 puts it in the graph, and M6 places it)"
                if a.is_vacancy else
                "two donors do not bond to each other; a join needs somewhere to put one, "
                "and on a metal that somewhere is a vacancy")
        return Compatibility(False, float("inf"), f"{kind}: {what}", mode=mode)

    donor, vacancy = pair
    for site, label in ((donor, donor.donor_type or "donor"), (vacancy, "vacancy")):
        if mode not in site.binding_modes:
            offered = ", ".join(site.binding_modes) or "nothing"
            return Compatibility(
                False, float("inf"),
                f"{label} does not bind {mode!r}; it offers {offered}", mode=mode)

    d_ml = None
    metal = _metal_symbol(partner)
    if metal is not None and donor_element is not None:
        d_ml = metal_donor_distance(metal, donor_element).value

    wells = torsion_wells(donor.donor_type, BindingMode(mode))
    lobes = _lone_pairs(donor)
    notes = []
    frame = _frame_of(donor)
    if frame is None:
        notes.append("donor has no frame (perceived without a geometry), so the verdict "
                     "is from roles and modes alone")
    elif lobes > 1:
        # Not a refusal and not a preference: both lobes are real, and this says so
        # rather than letting the caller assume the stored frame is the only direction.
        notes.append(f"this donor offers {lobes} in-plane lone pairs, so `lone_pair=` is "
                     "a choice with no default answer in the chemistry — the metal goes "
                     "where the caller puts it")
    if frame is not None and (donor.frame or {}).get("mode") == "fallback":
        # Not a refusal.  A fallback frame is a real answer from a donor with no
        # neighbours to orient it, and the caller should know the axis was inferred from
        # the molecule's centroid rather than from bonding.
        notes.append("donor frame is a fallback (axis inferred from the molecule's "
                     "centroid, not from its bonding)")
    reason = (f"{donor.donor_type or 'donor'} onto a vacant vertex: a rigid move realises "
              f"this alignment exactly, so there is no residual strain to report")
    if notes:
        reason += " — " + "; ".join(notes)
    return Compatibility(True, 0.0, reason, mode=mode, wells=wells, d_ml=d_ml,
                         lone_pairs=lobes)


def chelate_compatible(donors: Sequence[Site], vacancies: Sequence[Site], *,
                       partner: Any | None = None,
                       donor_elements: Sequence[str] | None = None,
                       tolerance_deg: float = MAX_BITE_MISMATCH_DEG) -> Compatibility:
    """Can this donor pair span these two vertices — the inter-block convergence question.

    This is where a join acquires geometry it cannot rigid-move away.  Two quantities,
    each an INTERNAL property of its own block and therefore invariant under the free
    transform that brings the blocks together:

    * the **bite angle** the ligand offers — the angle its two donors subtend at the metal
      position their own frames imply (D13: the frames are stored, so this is read, not
      re-searched);
    * the **vertex separation** the metal requires — the angle between the two vacancy
      axes, ~90 deg for a cis octahedral pair and 180 deg for a trans one.

    Their mismatch is the strain, and it is what makes cis/trans a feasibility question
    rather than a labelling one: a chelate that spans a cis pair comfortably cannot reach
    a trans pair at all, which is the D10 anthrarufin case seen from the assembly side.

    Frames are required here and their absence RAISES rather than returning infeasible —
    the whole verdict is derived from them, and "unknown" reported as "impossible" is the
    failure mode `BuildingBlock.open_sites` already refuses.
    """
    if len(donors) != 2 or len(vacancies) != 2:
        raise ValueError(
            f"a chelate spans exactly two donors and two vertices; got {len(donors)} and "
            f"{len(vacancies)}. Bridging modes (mu2/mu3) put ONE donor on vertices of "
            f"DIFFERENT metals and are a different question (M6).")
    mode = BindingMode.CHELATE.value
    if donors[0].atom_idx == donors[1].atom_idx:
        return Compatibility(False, float("inf"),
                             f"both donors are atom {donors[0].atom_idx}; one atom cannot "
                             "chelate a metal by itself", mode=mode)
    if vacancies[0].atom_idx != vacancies[1].atom_idx:
        return Compatibility(
            False, float("inf"),
            f"the two vertices are on different metals ({vacancies[0].atom_idx} and "
            f"{vacancies[1].atom_idx}); one ligand across two centres is a mu2 bridge, "
            "not a chelate", mode=mode)
    for donor in donors:
        if mode not in donor.binding_modes:
            return Compatibility(
                False, float("inf"),
                f"{donor.donor_type} does not chelate; it offers "
                f"{', '.join(donor.binding_modes)}", mode=mode)

    frames = [_frame_of(s) for s in (*donors, *vacancies)]
    if any(f is None for f in frames):
        raise ValueError(
            "chelate compatibility is computed FROM the frames — bite angle on the ligand "
            "side, vertex separation on the metal side — and at least one site has none. "
            "Perceive with a geometry (or load the stored frames) before asking; a "
            "missing frame is an unknown answer, not a negative one.")
    (o_a, ax_a, _), (o_b, ax_b, _), (_, v_a, _), (_, v_b, _) = frames

    metal = _metal_symbol(partner)
    if metal is not None and donor_elements is not None and len(donor_elements) == 2:
        d = [metal_donor_distance(metal, el).value for el in donor_elements]
    else:
        d = [NOMINAL_D_ML, NOMINAL_D_ML]
    # Where each donor's own frame says the metal goes, and therefore where it actually
    # would go: the midpoint, the same construction `sites.model._convergence` uses.
    p_a, p_b = o_a + d[0] * ax_a, o_b + d[1] * ax_b
    implied_metal = (p_a + p_b) / 2.0
    miss = float(np.linalg.norm(p_a - p_b))
    bite = _angle_deg(o_a - implied_metal, o_b - implied_metal)
    separation = _angle_deg(v_a, v_b)

    mismatch = abs(bite - separation)
    strain = mismatch / tolerance_deg
    feasible = mismatch <= tolerance_deg
    verdict = "spans" if feasible else "cannot reach"
    reason = (f"the pocket offers a {bite:.1f} deg bite and the vertices sit {separation:.1f} "
              f"deg apart, so it {verdict} them (mismatch {mismatch:.1f} deg against a "
              f"{tolerance_deg:.0f} deg tolerance); the two stored frames put the metal "
              f"{miss:.2f} A apart, which is reported and not judged — whether a pocket "
              f"converges at all is `sites.model.chelate_pockets`' question and it searches "
              f"all four torsion wells to answer it, while these frames are the stored ones")
    return Compatibility(feasible, strain, reason, mode=mode,
                         wells=torsion_wells(donors[0].donor_type, BindingMode.CHELATE),
                         d_ml=d[0] if metal is not None else None)


def _coords_of(block: BuildingBlock, label: str) -> np.ndarray | None:
    """A block's coordinates as an (n_atoms, 3) array, in its own atom order.

    That IS the convention — `BuildingBlock.geometry` is typed `Any` because the object
    that carries coordinates differs by producer, so the one thing fixed here is the
    shape and the ordering: row `k` is the atom at `graph.nodes()[k]`.  Anything that
    `numpy` can read as (n, 3) works; a `.coords` attribute is unwrapped first so a
    placer result can be handed over whole.
    """
    geometry = block.geometry
    if geometry is None:
        return None
    geometry = getattr(geometry, "coords", geometry)
    coords = np.asarray(geometry, dtype=float)
    if coords.shape != (len(block.graph), 3):
        raise ValueError(
            f"{label} geometry is {coords.shape}, but its graph has {len(block.graph)} "
            f"atoms. Row k must be the atom at graph.nodes()[k] — a geometry in a "
            f"different order would put the new bond on the wrong atom.")
    return coords


def _transform_frame(frame: dict | None, rot: np.ndarray, origin: np.ndarray,
                     target: np.ndarray) -> dict | None:
    """Carry a frame through the rigid move its block just made.

    `sites.inherit` says frames survive a join untouched, and its reasoning holds for
    every atom the join did not move.  A join DOES move one whole block, so that block's
    frames move with it — which is not an exception to the rule but the same rule applied
    to the block that travelled.  The move is rigid, so every relationship a frame encodes
    (the outward direction, the torsion zero, the angle between them) is preserved
    exactly; nothing is re-derived and no search re-runs.
    """
    if not frame:
        return frame
    out = dict(frame)
    if "origin" in frame:
        out["origin"] = [float(x) for x in
                         rot @ (np.asarray(frame["origin"], dtype=float) - origin) + target]
    for key in ("axis", "ref"):
        if key in frame:
            out[key] = [float(x) for x in rot @ np.asarray(frame[key], dtype=float)]
    # The other lone-pair lobes are frames too, and a lobe left behind by the move would
    # point at where the metal used to be — so a second bridge onto an already-moved block
    # would be judged against a stale direction.  Recursing is the whole implementation:
    # a lobe has no lobes of its own, so this bottoms out immediately.
    if frame.get("lone_pairs"):
        out["lone_pairs"] = [_transform_frame(lobe, rot, origin, target)
                             for lobe in frame["lone_pairs"]]
    return out


def _merged_graph(a: TypedGraph, b: TypedGraph, name: str,
                  ) -> tuple[TypedGraph, AtomMap, AtomMap]:
    """`a`'s atoms then `b`'s, with both parents' `{parent: child}` maps.

    Built through `to_dict`/`from_dict` so bond order travels too.  The ordering is
    deterministic but otherwise arbitrary, and that is safe rather than lucky: L1 is a
    canonical certificate, so the same product built in either order hashes the same —
    which is the "one node, two routes" gate, and the reason a join may pick whatever
    ordering is simplest to reason about.
    """
    da, db = a.to_dict(), b.to_dict()
    offset = len(da["atoms"])
    map_a = {atom["i"]: k for k, atom in enumerate(da["atoms"])}
    map_b = {atom["i"]: offset + k for k, atom in enumerate(db["atoms"])}
    merged = {
        "schema": da["schema"], "charge": None, "multiplicity": None, "name": name,
        "atoms": ([dict(atom, i=map_a[atom["i"]]) for atom in da["atoms"]]
                  + [dict(atom, i=map_b[atom["i"]]) for atom in db["atoms"]]),
        "bonds": ([dict(bond, i=map_a[bond["i"]], j=map_a[bond["j"]]) for bond in da["bonds"]]
                  + [dict(bond, i=map_b[bond["i"]], j=map_b[bond["j"]]) for bond in db["bonds"]]),
    }
    return TypedGraph.from_dict(merged), map_a, map_b


def join(a: BuildingBlock, b: BuildingBlock, site_a: Site, site_b: Site, *,
         mode: str = "mono", torsion_well: int = 0, lone_pair: int = 0, seed: int = 0,
         with_geometry: bool = True) -> JoinResult:
    """Join two blocks at a pair of sites and return the product.

    Aligns by frame, adds the typed bond, carries both atom maps, inherits the parents'
    sites minus those consumed, and emits the choice vector that regenerates the result.

    **The ligand moves and the metal stays put.**  The block carrying the vacancy is the
    anchor: its coordinates, and therefore every frame on it, are untouched, and the block
    carrying the donor is rigidly transformed onto it.  That is what makes repeated joins
    cheap and stable — a growing cluster keeps one coordinate system across every addition
    instead of being re-expressed after each one.

    **What is consumed, and what merely changes status.**  The VACANCY is consumed: a
    filled vertex is not a vertex, and it is consumed by `(atom_idx, slot)` so the metal's
    other vertices survive.  The DONOR is not.  §6.3 says "minus sites consumed by the new
    bond", and the vacancy is unambiguously that; a bound carboxylate oxygen, however, is
    still a donor site that exists — which is exactly the distinction D5 draws between
    `site_catalog` (which sites exist, geometry-free) and `site_state` (what they are
    doing).  So the donor is inherited and `refresh_state` marks it `OCCUPIED` off the
    dative bond in the product graph.  The open-site count falls by two either way; the
    difference is that the product can still answer "what bound here, and how easily did
    it activate", which is what a reaction edge needs to explain itself.

    **Which lone pair the metal binds is a choice, and it is this caller's.**  An sp2
    donor has two in-plane lobes pointing in genuinely different directions, and for a
    carboxylate bridging two metals the choice is worth 2.8 A of M...M — anti-anti at
    5.516, syn-anti at 5.148, syn-syn at 2.673, which is the paddlewheel.  `lone_pair` is
    an index into the set the donor actually offers, wrapped modulo its size, and it
    travels in the choice vector so a replay reproduces the lobe rather than the default.
    A donor whose direction its bonding determines has one lobe and ignores it.

    **A multi-centre product is no longer a stopping point (D20).**  One donor onto one
    vertex is one contact, and one contact is satisfied by a rigid move of the donor's
    block whatever either block already contains — so the geometry is determined and the
    M...M that results is an OUTPUT, which is the whole of the emergence claim.  What
    still needs `place_multicentre` is the case this function does not do: two contacts
    that must be satisfied at once across centres whose separation two determinants
    disagree about.
    """
    pair = _roles(site_a, site_b)
    if pair is None:
        raise IncompatibleJoin(compatible(site_a, site_b, mode=mode))
    donor_site, vacancy_site = pair
    donor_block, metal_block = (a, b) if donor_site is site_a else (b, a)
    for site, block, which in ((donor_site, donor_block, "donor"),
                               (vacancy_site, metal_block, "vacancy")):
        if site not in block.sites:
            raise ValueError(
                f"the {which} site is not on the block it was passed with. `join(a, b, "
                f"site_a, site_b)` reads site_a off a and site_b off b; swapping them "
                f"builds an atom map onto the wrong parent, which nothing downstream can "
                f"detect.")

    metal_element = metal_block.graph.label(vacancy_site.atom_idx).element
    donor_element = donor_block.graph.label(donor_site.atom_idx).element
    verdict = compatible(site_a, site_b, partner=metal_element, mode=mode,
                         donor_element=donor_element)
    if not verdict.feasible:
        raise IncompatibleJoin(verdict)

    product, map_a, map_b = _merged_graph(
        a.graph, b.graph, name=f"{a.graph.name or 'a'}+{b.graph.name or 'b'}")
    map_donor = map_a if donor_block is a else map_b
    map_metal = map_b if donor_block is a else map_a
    product.add_bond(map_donor[donor_site.atom_idx], map_metal[vacancy_site.atom_idx],
                     EdgeType.DATIVE)

    # Charge and multiplicity are ASKED for, never assumed: `net_charge` raises when a
    # parent never declared one, which is ground rule 5 doing its job one layer up.
    product.charge = a.graph.net_charge() + b.graph.net_charge()
    if a.graph.multiplicity is None or b.graph.multiplicity is None:
        raise AmbiguousSpecError(
            "a parent has no multiplicity, so the product's is not derivable. Unpaired "
            "electrons add (energy.backends.combined_multiplicity) and there is no "
            "default that is not a guess about spin state.")
    product.multiplicity = combined_multiplicity(a.graph.multiplicity, b.graph.multiplicity)

    wells = verdict.wells or (0.0,)
    well_deg = wells[torsion_well % len(wells)]
    lobe = int(lone_pair) % max(verdict.lone_pairs, 1)
    d_ml = verdict.d_ml if verdict.d_ml is not None else NOMINAL_D_ML
    choice_vector = {
        "op": "join", "mode": mode,
        "torsion_well": int(torsion_well) % len(wells), "torsion_deg": well_deg,
        "donor": {"atom": donor_site.atom_idx, "type": donor_site.donor_type,
                  "element": donor_element, "lone_pair": lobe,
                  "block": donor_block.structure_id},
        "vacancy": {"atom": vacancy_site.atom_idx, "slot": vacancy_site.slot,
                    "metal": metal_element, "block": metal_block.structure_id},
        "d_ml": round(float(d_ml), 6),
        "order": "donor-block-first" if donor_block is a else "metal-block-first",
        "seed": int(seed),
    }

    donor_sites = donor_block.sites
    coords = None
    if with_geometry:
        coords, donor_sites = _place_donor_block(
            donor_block, donor_site, vacancy_site, d_ml=d_ml, well_deg=well_deg,
            lone_pair=lobe)

    sites = merge_inherited(
        inherit_sites(donor_sites, map_donor),
        inherit_sites(metal_block.sites, map_metal,
                      consumed=[(vacancy_site.atom_idx, vacancy_site.slot)]),
    )

    geometry, state = None, None
    if coords is not None:
        metal_coords = _coords_of(metal_block, "the metal block")
        if metal_coords is not None:
            geometry = (np.vstack([coords, metal_coords]) if donor_block is a
                        else np.vstack([metal_coords, coords]))
            symbols = [product.label(i).element for i in product.nodes()]
            states = refresh_state(sites, geometry, graph=product, symbols=symbols)
            state = {(s.atom_idx, s.slot): s for s in states}

    return JoinResult(
        block=BuildingBlock(graph=product, sites=tuple(sites), geometry=geometry,
                            state=state),
        atom_map=map_a, choice_vector=choice_vector, strain=verdict.strain,
        partner_atom_map=map_b, compatibility=verdict)


def _place_donor_block(donor_block: BuildingBlock, donor_site: Site, vacancy_site: Site,
                       *, d_ml: float, well_deg: float, lone_pair: int = 0,
                       ) -> tuple[np.ndarray | None, tuple[Site, ...]]:
    """Rigidly move the donor's block onto the vacant vertex.  The N=1 alignment.

    Three constraints, which between them fix the pose completely — which is D13's whole
    argument for a frame over a vector:

    1. the donor atom lands at `metal + d_ml * vertex_axis`;
    2. its chosen lone-pair axis points back at the metal (antiparallel to the vertex
       axis).  WHICH lobe that is comes from the caller, because an sp2 donor's two are
       genuinely different directions and nothing in the donor decides between them;
    3. the roll about the new bond is the chosen torsion WELL, not whatever the arithmetic
       happened to leave — an index into a discrete set, so it replays exactly.
    """
    coords = _coords_of(donor_block, "the donor block")
    donor_frame = _frame_of(donor_site, lone_pair)
    vacancy_frame = _frame_of(vacancy_site)
    if coords is None or donor_frame is None or vacancy_frame is None:
        missing = ("coordinates" if coords is None else
                   "a donor frame" if donor_frame is None else "a vacancy frame")
        raise ValueError(
            f"cannot align this join: {missing} is absent. Alignment is computed from the "
            f"stored frames (D13) and applied to real coordinates; pass "
            f"with_geometry=False to build the graph-level product instead.")

    o_d, ax_d, ref_d = donor_frame
    o_v, ax_v, ref_v = vacancy_frame
    # FRAME onto frame, not axis onto axis.  Matching axes alone leaves the roll about the
    # new bond undetermined, and `rotation_between` then settles it with its minimal
    # rotation — which depends on how the ligand happened to be oriented in its own
    # coordinates.  Re-embedding a ligand rotates it rigidly, so the SAME choice vector
    # produced products whose rings sat at different azimuths: measured at 2.3 A per ring
    # atom, 1.56 A core RMSD, for a ligand whose own geometry was bit-identical between
    # the two runs.  That is precisely what D13 says a lone outward vector cannot do and a
    # frame can, and this function was using the vector half of the frame it was handed.
    rot = _triad(-ax_v, _rotate_about(ref_v, ax_v, math.radians(well_deg))) @ _triad(
        ax_d, ref_d).T
    target = o_v + d_ml * ax_v
    moved = (rot @ (coords - o_d).T).T + target
    sites = tuple(replace(s, frame=_transform_frame(s.frame, rot, o_d, target))
                  for s in donor_block.sites)
    return moved, sites


def join_chelate(a: BuildingBlock, b: BuildingBlock, sites_a: Sequence[Site],
                 sites_b: Sequence[Site], *, seed: int = 0,
                 with_geometry: bool = True) -> JoinResult:
    """Join two blocks at TWO points at once — one ligand across two vertices of one metal.

    The same operation as `join` and a different geometry problem.  A single-point join has
    a free rigid move and therefore no residual strain (`compatible` says so and why); two
    points at once must satisfy both constraints with one rigid move, so the bite angle the
    ligand offers and the separation of the vertices it is asked to span have to agree.
    `chelate_compatible` is the verdict and it is asked first; this places what it passed.

    **How the mismatch is absorbed.**  The pose is fixed by three things, in this order:
    the donor-donor line onto the vertex-vertex line, the roll about that line set by
    requiring the ligand's own donor axes to point back at the metal, and then a push
    along the bisector so that BOTH donors sit at their metal-donor distance.  The residual
    goes into the bite angle, which is the quantity `chelate_compatible` already bounded,
    rather than into the M-D bond lengths, which `geometry.qc` checks and which an
    optimiser cannot be asked to repair (`QCReport.marginal` excludes bad bonds for exactly
    this reason).

    Both vertices are consumed; both donors are inherited and become `OCCUPIED` off the
    product's dative bonds — the same split `join` documents, for the same reason.

    Both vertices are on ONE metal — `chelate_compatible` refuses the other case by name,
    and that refusal is where a µ2 bridge begins.
    """
    if len(sites_a) != 2 or len(sites_b) != 2:
        raise ValueError(
            f"a chelate join takes two sites on each block; got {len(sites_a)} and "
            f"{len(sites_b)}. One donor onto one vertex is `join`; one donor across two "
            f"METALS is a mu2 bridge and is M6.")
    # Each side has to be entirely one role: two donors here, two vertices there.  A
    # mixed pair is refused as a pairing rather than caught later as an atom-map error,
    # because "one of these four is on the wrong side" is the thing the caller got wrong.
    a_vacant, b_vacant = {s.is_vacancy for s in sites_a}, {s.is_vacancy for s in sites_b}
    if len(a_vacant) != 1 or len(b_vacant) != 1 or a_vacant == b_vacant:
        return _refuse_chelate(sites_a, sites_b, "this pairing is not two donors on one "
                                                 "block and two vertices on the other")
    donors, vacancies = (sites_b, sites_a) if True in a_vacant else (sites_a, sites_b)
    donor_block, metal_block = (b, a) if True in a_vacant else (a, b)
    for pair, block, which in ((donors, donor_block, "donor"),
                               (vacancies, metal_block, "vacancy")):
        for site in pair:
            if site not in block.sites:
                raise ValueError(
                    f"a {which} site is not on the block it was passed with. "
                    f"`join_chelate(a, b, sites_a, sites_b)` reads sites_a off a and "
                    f"sites_b off b; swapping them builds an atom map onto the wrong "
                    f"parent, which nothing downstream can detect.")

    metal_element = metal_block.graph.label(vacancies[0].atom_idx).element
    donor_elements = [donor_block.graph.label(s.atom_idx).element for s in donors]
    verdict = chelate_reach(donors, vacancies, partner=metal_element,
                            donor_elements=donor_elements)
    if not verdict.feasible:
        raise IncompatibleJoin(verdict)

    product, map_a, map_b = _merged_graph(
        a.graph, b.graph, name=f"{a.graph.name or 'a'}+{b.graph.name or 'b'}")
    map_donor = map_a if donor_block is a else map_b
    map_metal = map_b if donor_block is a else map_a
    for donor in donors:
        product.add_bond(map_donor[donor.atom_idx], map_metal[vacancies[0].atom_idx],
                         EdgeType.DATIVE)

    product.charge = a.graph.net_charge() + b.graph.net_charge()
    if a.graph.multiplicity is None or b.graph.multiplicity is None:
        raise AmbiguousSpecError(
            "a parent has no multiplicity, so the product's is not derivable. Unpaired "
            "electrons add (energy.backends.combined_multiplicity) and there is no "
            "default that is not a guess about spin state.")
    product.multiplicity = combined_multiplicity(a.graph.multiplicity, b.graph.multiplicity)

    distances = [metal_donor_distance(metal_element, el).value for el in donor_elements]
    choice_vector = {
        "op": "join", "mode": BindingMode.CHELATE.value,
        # No torsion index: the roll about the new bonds is DETERMINED by the second
        # contact, so there is no well to choose.  What replaces it as the discrete
        # choice is which donor went to which vertex, and that is recorded in order.
        "donors": [{"atom": s.atom_idx, "type": s.donor_type, "element": el,
                    "block": donor_block.structure_id}
                   for s, el in zip(donors, donor_elements)],
        "vacancies": [{"atom": s.atom_idx, "slot": s.slot, "metal": metal_element,
                       "block": metal_block.structure_id} for s in vacancies],
        "d_ml": [round(float(x), 6) for x in distances],
        "order": "donor-block-first" if donor_block is a else "metal-block-first",
        "seed": int(seed),
    }

    donor_sites = donor_block.sites
    coords = None
    if with_geometry:
        coords, donor_sites = _place_chelating_block(donor_block, donors, vacancies,
                                                     distances=distances)

    sites = merge_inherited(
        inherit_sites(donor_sites, map_donor),
        inherit_sites(metal_block.sites, map_metal,
                      consumed=[(v.atom_idx, v.slot) for v in vacancies]),
    )

    geometry, state = None, None
    if coords is not None:
        metal_coords = _coords_of(metal_block, "the metal block")
        if metal_coords is not None:
            geometry = (np.vstack([coords, metal_coords]) if donor_block is a
                        else np.vstack([metal_coords, coords]))
            symbols = [product.label(i).element for i in product.nodes()]
            states = refresh_state(sites, geometry, graph=product, symbols=symbols)
            state = {(s.atom_idx, s.slot): s for s in states}

    return JoinResult(
        block=BuildingBlock(graph=product, sites=tuple(sites), geometry=geometry,
                            state=state),
        atom_map=map_a, choice_vector=choice_vector, strain=verdict.strain,
        partner_atom_map=map_b, compatibility=verdict)


def chelate_reach(donors: Sequence[Site], vacancies: Sequence[Site], *,
                  partner: Any | None = None,
                  donor_elements: Sequence[str] | None = None,
                  tolerance_deg: float = MAX_BITE_MISMATCH_DEG) -> Compatibility:
    """Can a RIGID ligand reach both vertices with both bonds at their own length?

    A different question from `chelate_compatible`'s, and the one a PLACEMENT has to
    answer.  That function compares the direction the pocket's frames converge on with
    the separation of the two vertices — the pocket's intent — and reports, without
    judging, how far apart the two frames put the metal.  When those frames diverge (a
    real case: `sites.model.chelate_pockets` establishes convergence by searching torsion
    wells, and a stored frame is one well, not the converging one) that angle says
    "trans" about a pair of donors sitting 3 A apart, which no rigid move can stretch
    across an octahedron.

    So the verdict here is distance geometry, computed from the same three numbers the
    placement uses: the donors are a fixed distance apart, each has to sit at its own
    metal-donor distance, and those two facts fix the angle the pair actually subtends at
    the centre.  Comparing THAT with the vertex separation is what decides cis from trans
    for a rigid placement — and it agrees with `chelate_compatible` for every pocket
    whose frames do converge.  The pocket verdict travels in the reason, because "these
    frames say something different" is worth seeing.

    This is the gate `join_chelate` uses, and the one a caller should rank candidate
    vertex pairs by — ranking on one criterion and placing under another is how a step
    picks the pair it then cannot build.
    """
    metal = _metal_symbol(partner)
    pocket = chelate_compatible(donors, vacancies, partner=partner,
                                donor_elements=donor_elements)
    if not pocket.feasible and math.isinf(pocket.strain):
        # Structural, not geometric: one atom chelating itself, two vertices on two
        # different metals, a donor that does not chelate at all.  None of those is a
        # placement that could be improved by moving anything.
        return pocket
    if metal is None or not donor_elements or len(donor_elements) != 2:
        # Without the pair's elements there is no bond length to place the donors at, so
        # the pocket's own verdict is the only one available.
        return pocket
    frames = [_frame_of(s) for s in (*donors, *vacancies)]
    if any(f is None for f in frames):
        raise ValueError(
            "a chelate placement is computed FROM the frames and at least one site has "
            "none. Perceive with a geometry before joining; a missing frame is an "
            "unknown answer, not a negative one.")
    (p_1, _, _), (p_2, _, _), (_, v_1, _), (_, v_2, _) = frames
    d_mean = sum(metal_donor_distance(metal, el).value for el in donor_elements) / 2.0
    half = float(np.linalg.norm(p_2 - p_1)) / 2.0
    subtended = 2.0 * math.degrees(math.asin(min(half / d_mean, 1.0)))
    separation = _angle_deg(v_1, v_2)
    mismatch = abs(subtended - separation)
    feasible = mismatch <= tolerance_deg
    reason = (
        f"the two donors sit {2 * half:.2f} A apart, so at {d_mean:.2f} A bonds they "
        f"subtend {subtended:.1f} deg at the centre, and the vertices sit "
        f"{separation:.1f} deg apart: the ligand {'reaches' if feasible else 'cannot reach'} "
        f"both (mismatch {mismatch:.1f} deg against a {tolerance_deg:.0f} deg tolerance)")
    if pocket.feasible != feasible:
        reason += f" — note the pocket's own frames disagree: {pocket.reason}"
    return Compatibility(feasible, mismatch / tolerance_deg, reason,
                         mode=BindingMode.CHELATE.value, wells=pocket.wells,
                         d_ml=pocket.d_ml)


def _refuse_chelate(sites_a: Sequence[Site], sites_b: Sequence[Site],
                    what: str) -> JoinResult:
    roles = [("vacancy" if s.is_vacancy else "donor") for s in (*sites_a, *sites_b)]
    raise IncompatibleJoin(Compatibility(
        False, float("inf"),
        f"{what} (got {', '.join(roles)}): a chelate join puts two donors from one block "
        f"on two vertices of one metal on the other",
        mode=BindingMode.CHELATE.value))


def _place_chelating_block(donor_block: BuildingBlock, donors: Sequence[Site],
                           vacancies: Sequence[Site], *, distances: Sequence[float],
                           ) -> tuple[np.ndarray | None, tuple[Site, ...]]:
    """Rigidly move a chelating block onto two vacant vertices.  The N=2 alignment."""
    coords = _coords_of(donor_block, "the donor block")
    frames = [_frame_of(s) for s in (*donors, *vacancies)]
    if coords is None or any(f is None for f in frames):
        raise ValueError(
            "cannot align this chelate join: coordinates or a frame are absent. The pose "
            "is computed from the stored frames (D13) and applied to real coordinates; "
            "pass with_geometry=False to build the graph-level product instead.")
    (p_1, ax_1, _), (p_2, ax_2, _), (metal, v_1, _), (_, v_2, _) = frames
    pocket_out = _pocket_outward(donor_block, coords, donors, (ax_1, ax_2))

    t_1 = metal + distances[0] * v_1
    t_2 = metal + distances[1] * v_2
    span, target_span = _unit(p_2 - p_1), _unit(t_2 - t_1)
    rot = _rotation_between(span, target_span)

    # The roll about the new span: the ligand's own donor axes point at where it expects
    # the metal, so turn them onto where the metal actually is.  Everything here is
    # measured PERPENDICULAR to the span, which is the only direction the roll can move.
    ligand_view = _perpendicular(rot @ pocket_out, target_span)
    metal_view = _perpendicular(metal - (t_1 + t_2) / 2.0, target_span)
    if ligand_view is not None and metal_view is not None:
        theta = math.atan2(float(np.dot(np.cross(ligand_view, metal_view), target_span)),
                           float(np.dot(ligand_view, metal_view)))
        rot = _axis_rotation(target_span, theta) @ rot

    # Then slide along the perpendicular until both donors sit at their M-D distance.
    # The bite mismatch has to go somewhere and this is the choice: it goes into the
    # ANGLE, which `chelate_compatible` has already bounded, rather than into the bond
    # lengths, which QC checks and which no optimiser can be asked to undo.
    #
    # `metal_view` is the direction from the donors' midpoint to the metal, and it is
    # undefined for a TRANS pair — the two vertices are collinear through the centre, so
    # their midpoint IS the centre and points nowhere.  That is not a degenerate request:
    # a ligand with a wide enough bite really does span trans (D10's cis/trans question
    # from the other side), and the midpoint of its donors then sits on the metal.  The
    # ligand's own view of where the metal goes is what orients it in that case, and it
    # is the same direction in every case the two are both defined.
    half = float(np.linalg.norm(p_2 - p_1)) / 2.0
    d_mean = (distances[0] + distances[1]) / 2.0
    toward = metal_view if metal_view is not None else ligand_view
    if toward is None:
        toward = _perpendicular(np.array([1.0, 0.0, 0.0]), target_span)
        if toward is None:
            toward = _perpendicular(np.array([0.0, 1.0, 0.0]), target_span)
    # Offset zero means the donors straddle the centre, which is what a bite wider than
    # the bond length asks for; it is reported by QC rather than fudged to a minimum.
    offset = math.sqrt(max(d_mean ** 2 - half ** 2, 0.0))
    target = metal - offset * toward
    moved = (rot @ (coords - (p_1 + p_2) / 2.0).T).T + target
    sites = tuple(replace(s, frame=_transform_frame(s.frame, rot, (p_1 + p_2) / 2.0, target))
                  for s in donor_block.sites)
    return moved, sites


def _pocket_outward(block: BuildingBlock, coords: np.ndarray, donors: Sequence[Site],
                    axes: Sequence[np.ndarray]) -> np.ndarray:
    """Which way out of the ligand the metal lies — the pocket's own convergence direction.

    Taken from each donor's BONDING (the direction away from what it is attached to),
    summed over the pair, rather than from the stored frame axes.  The frames are one
    torsion well of a possible several, and for an sp2 donor with a single neighbour the
    two in-plane lobes point to opposite sides: an ortho diolate's two stored axes can
    come out nearly antiparallel, whose sum then points INTO the ring.  Rolling the
    ligand onto that vector puts the metal underneath the ring rather than in the pocket,
    which is a ligand wrapped around the centre and a wall of clashes.

    `sites.model.chelate_pockets` answers the same question by searching all four well
    combinations; this is the cheap local form of that answer, and it needs no search
    because the direction away from a donor's substituents does not depend on a well.
    A donor with no bonded neighbour (a bare halide) has no such direction and keeps its
    frame axis.
    """
    rows = {node: i for i, node in enumerate(block.graph.nodes())}
    out = np.zeros(3)
    for site, axis in zip(donors, axes):
        neighbours = [rows[n] for n in block.graph.neighbors(site.atom_idx) if n in rows]
        direction = axis
        if neighbours:
            away = coords[rows[site.atom_idx]] - coords[neighbours].mean(axis=0)
            if float(np.linalg.norm(away)) > 1e-9:
                direction = _unit(away)
        out = out + direction
    return _unit(out) if float(np.linalg.norm(out)) > 1e-9 else _unit(axes[0])


def _perpendicular(v: np.ndarray, axis: np.ndarray) -> np.ndarray | None:
    """`v` with its component along `axis` removed, or None if nothing is left of it."""
    out = np.asarray(v, dtype=float) - float(np.dot(v, axis)) * axis
    return None if float(np.linalg.norm(out)) < 1e-9 else _unit(out)


def grow(seed_block: BuildingBlock, partners: tuple[BuildingBlock, ...], *,
         degree: int, max_products: int = 1000, seed: int = 0,
         modes: tuple[str, ...] = (BindingMode.MONODENTATE.value,)) -> list[JoinResult]:
    """Repeatedly join partners onto a block, `degree` additions deep.

    This is what "iterate cluster formation to a certain degree" means, and it is a thin
    front on `assembly.construct.enumerate_constructions` — the branch tree is the real
    object and this is the shape the settled signature promised.  Each result's choice
    vector is the whole PATH (not the last join), and its `strain` is the worst step on
    that path: a construction is as viable as its weakest point, which is the same reading
    M8 will give a barrier.

    `BranchTree` is what to call when the refusals or the `capped` flag matter — a list
    cannot say "there were more" or "every branch was refused, here is why", and both of
    those are things a caller usually needs to know before believing a short answer.

    Growth ONTO a second metal centre still stops where S3 draws the line: the graph is
    built, and positioning it waits for `geometry.placer.place_multicentre` (M6).
    """
    # Imported here, not at module scope: `construct` is built ON this module, so the
    # dependency only runs in this direction at call time.
    from mofsbu.assembly.construct import constructions_as_joins, enumerate_constructions

    tree = enumerate_constructions(seed_block, tuple(partners), degree=degree,
                                   modes=tuple(modes), max_products=max_products)
    return constructions_as_joins(tree)
