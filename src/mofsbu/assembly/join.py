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
  atom map, and emit the choice-vector that regenerates it.  **S3.**
* `geometry.placer.place_multicentre` — inter-centre constraints (M-M distance, bridge
  bite angle).  The plan's declared headline cost (M6).
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
from mofsbu.sites.model import Site
from mofsbu.sites.state import refresh_state


#: Re-exported: it moved to `_types` (where the rest of the hierarchy lives) when this
#: module grew real dependencies on `sites` and `descriptors`, which had been importing
#: it from here.  Every existing `from mofsbu.assembly.join import NotBuiltYet` still
#: works, and none of them has to care that the definition moved.
__all__ = ["BuildingBlock", "Compatibility", "IncompatibleJoin", "JoinResult",
           "NotBuiltYet", "chelate_compatible", "compatible", "grow", "join"]


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


def _frame_of(site: Site) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """`(origin, axis, ref)` from a site's stored frame, or None if it has none."""
    frame = site.frame
    if not frame or "axis" not in frame or "origin" not in frame:
        return None
    return (np.asarray(frame["origin"], dtype=float),
            _unit(np.asarray(frame["axis"], dtype=float)),
            _unit(np.asarray(frame.get("ref", [0.0, 0.0, 1.0]), dtype=float)))


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
    notes = []
    frame = _frame_of(donor)
    if frame is None:
        notes.append("donor has no frame (perceived without a geometry), so the verdict "
                     "is from roles and modes alone")
    elif (donor.frame or {}).get("mode") == "fallback":
        # Not a refusal.  A fallback frame is a real answer from a donor with no
        # neighbours to orient it, and the caller should know the axis was inferred from
        # the molecule's centroid rather than from bonding.
        notes.append("donor frame is a fallback (axis inferred from the molecule's "
                     "centroid, not from its bonding)")
    reason = (f"{donor.donor_type or 'donor'} onto a vacant vertex: a rigid move realises "
              f"this alignment exactly, so there is no residual strain to report")
    if notes:
        reason += " — " + "; ".join(notes)
    return Compatibility(True, 0.0, reason, mode=mode, wells=wells, d_ml=d_ml)


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
         mode: str = "mono", torsion_well: int = 0, seed: int = 0,
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

    **The M5/M6 seam.**  A product with two or more metal centres gets its graph, its atom
    maps, its inherited sites and its choice vector — everything an identity, a registry
    row and a provenance edge need — and then raises `NotBuiltYet` at the geometry step,
    naming `place_multicentre`.  Pass `with_geometry=False` to take the graph-level product
    and stop there, which is how a paddlewheel's GRAPH is reachable from the assembly path
    before the placer that can position one exists.
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
    d_ml = verdict.d_ml if verdict.d_ml is not None else NOMINAL_D_ML
    choice_vector = {
        "op": "join", "mode": mode,
        "torsion_well": int(torsion_well) % len(wells), "torsion_deg": well_deg,
        "donor": {"atom": donor_site.atom_idx, "type": donor_site.donor_type,
                  "element": donor_element,
                  "block": donor_block.structure_id},
        "vacancy": {"atom": vacancy_site.atom_idx, "slot": vacancy_site.slot,
                    "metal": metal_element, "block": metal_block.structure_id},
        "d_ml": round(float(d_ml), 6),
        "order": "donor-block-first" if donor_block is a else "metal-block-first",
        "seed": int(seed),
    }

    donor_sites = donor_block.sites
    coords = None
    n_metals = len(product.metals())
    if with_geometry:
        if n_metals > 1:
            raise NotBuiltYet(
                f"this join makes a {n_metals}-centre product, and positioning one needs "
                f"geometry.placer.place_multicentre (M6): the M-M distance and the bridge "
                f"bite angle have to be solved together with each centre's local "
                f"geometry, which single-centre alignment cannot do. The graph, the atom "
                f"maps, the inherited sites and the choice vector are all derivable — "
                f"pass with_geometry=False to take them.")
        coords, donor_sites = _place_donor_block(
            donor_block, donor_site, vacancy_site, d_ml=d_ml, well_deg=well_deg)

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
                       *, d_ml: float, well_deg: float,
                       ) -> tuple[np.ndarray | None, tuple[Site, ...]]:
    """Rigidly move the donor's block onto the vacant vertex.  The N=1 alignment.

    Three constraints, which between them fix the pose completely — which is D13's whole
    argument for a frame over a vector:

    1. the donor atom lands at `metal + d_ml * vertex_axis`;
    2. its outward axis points back at the metal (antiparallel to the vertex axis);
    3. the roll about the new bond is the chosen torsion WELL, not whatever the arithmetic
       happened to leave — an index into a discrete set, so it replays exactly.
    """
    coords = _coords_of(donor_block, "the donor block")
    donor_frame, vacancy_frame = _frame_of(donor_site), _frame_of(vacancy_site)
    if coords is None or donor_frame is None or vacancy_frame is None:
        missing = ("coordinates" if coords is None else
                   "a donor frame" if donor_frame is None else "a vacancy frame")
        raise ValueError(
            f"cannot align this join: {missing} is absent. Alignment is computed from the "
            f"stored frames (D13) and applied to real coordinates; pass "
            f"with_geometry=False to build the graph-level product instead.")

    o_d, ax_d, _ = donor_frame
    o_v, ax_v, _ = vacancy_frame
    rot = _axis_rotation(ax_v, math.radians(well_deg)) @ _rotation_between(ax_d, -ax_v)
    target = o_v + d_ml * ax_v
    moved = (rot @ (coords - o_d).T).T + target
    sites = tuple(replace(s, frame=_transform_frame(s.frame, rot, o_d, target))
                  for s in donor_block.sites)
    return moved, sites


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
