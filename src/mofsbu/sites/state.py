"""Per-geometry site state — the second tier of D5's split.

`site_catalog` says which donor sites EXIST and is perceived once per structure.  This
module says what they are DOING in one particular geometry, and is recomputed whenever a
better geometry arrives.  §6.3: "accessibility genuinely changes when a buried site opens
after a higher-fidelity relax", which is the whole reason the two tiers are separate
tables instead of one.

**`refresh_state` never perceives.**  It takes the sites it is given.  Re-running
perception here would re-derive the catalog from a moved geometry, and the moment two
geometries of one structure disagreed about which atoms are donors, the canonical indices
would stop meaning the same thing across the fidelity ladder — which is precisely what
storing sites by canonical index (D5) exists to prevent.  There is a test that counts
`find_donor_sites` calls and fails if this function makes any.

**Status is three-valued, and the third value is not a synonym for the second.**
`OCCUPIED` is a graph fact: this donor has a dative bond to a metal.  `BLOCKED` is a
geometry fact: nothing is bonded here, but nothing can reach it either.  Collapsing them
into "not open" would lose the distinction the assembly layer needs — an occupied site is
finished, a blocked one might open after a relaxation, and D11 makes a flip of exactly
that flag the primary trigger for a new conformer.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Sequence

from mofsbu._types import Fidelity, MethodSpec
from mofsbu.descriptors.ease import EaseRecord, activation_ease, occlusion
from mofsbu.graph._types import EdgeType, TypedGraph
from mofsbu.versions import ALGO_VERSIONS


class SiteStatus(str, Enum):
    OPEN = "open"
    OCCUPIED = "occupied"
    BLOCKED = "blocked"


#: Above this fraction of the approach cone blocked, a donor with no metal on it is
#: reported as sterically BLOCKED rather than open.  Chosen so that a donor whose entire
#: hemisphere of approach is walled off reads as blocked while ordinary substituent
#: crowding does not: a bare aqua O sits near 0.0, a carboxylate O in an open carboxylate
#: near 0.2-0.4, and a donor pointing into its own ring system above 0.8.
BLOCKED_OCCLUSION = 0.80

STATE_METHOD = MethodSpec(
    code="heuristic", code_version=ALGO_VERSIONS["site_state"],
    method="cone-occlusion+graph-status",
    extras={"algo": ALGO_VERSIONS["site_state"]})


@dataclass(frozen=True)
class SiteState:
    """One site, in one geometry.  Maps 1:1 onto a `site_state` row."""

    atom_idx: int
    status: SiteStatus
    buried_vol: float | None
    in_pocket: bool
    ease: EaseRecord | None
    pka: float | None = None
    fukui: float | None = None
    marginal_de: float | None = None

    @property
    def is_open(self) -> bool:
        return self.status is SiteStatus.OPEN


def buried_volume(coords: Any, site: Any, *, radius: float = 3.5,
                  symbols: Sequence[str] | None = None) -> float:
    """Blocked fraction of this site's approach cone, in [0, 1].

    Named `buried_volume` because that is the name §6.4 and the `site_state.buried_vol`
    column use, but it is a cone-occlusion fraction, not a %V_bur volume integral — see
    `descriptors.ease.occlusion` for why the cone is cast along the site's own axis
    instead of over a sphere.  The two agree in ordering and not in units, so the column
    is a ranking coordinate and nothing here converts it to a percentage it is not.
    """
    frame = getattr(site, "frame", None)
    if not frame or "axis" not in frame:
        # No frame means no direction to ask about.  A sphere-average would be a
        # different quantity wearing this one's name, so the answer is "unknown".
        return None                                              # type: ignore[return-value]
    return occlusion(coords, site.atom_idx, frame["axis"], radius=radius,
                     elements=symbols)


def _occupied(graph: TypedGraph | None, atom_idx: int) -> bool:
    """Is this donor already dative-bonded to a metal?"""
    if graph is None or atom_idx >= len(graph):
        return False
    return any(graph.is_metal(j)
               for j in graph.neighbors(atom_idx, EdgeType.DATIVE))


def refresh_state(
    sites: Sequence[Any],
    coords: Any,
    *,
    graph: TypedGraph | None = None,
    symbols: Sequence[str] | None = None,
    pocket_donors: frozenset[int] = frozenset(),
    fidelity: Fidelity = Fidelity.RAW,
) -> list[SiteState]:
    """Recompute per-geometry state for an already-perceived site list.

    `pocket_donors` is the set of atom indices that sit in a chelate pocket, from
    `sites.model.chelate_pockets`.  It is passed IN rather than computed here so that this
    function stays a pure function of (sites, coords) and the pocket search — which needs
    the molecule, not just coordinates — happens once at the call site.

    `fidelity` is the fidelity of the GEOMETRY the state was computed on, and it travels
    onto the row so a state derived from a raw construct is never mistaken for one derived
    from an xTB minimum.  It is not the ease record's own fidelity, which is HEURISTIC
    until something computes a component (D18).
    """
    out: list[SiteState] = []
    for site in sites:
        occluded = buried_volume(coords, site, symbols=symbols)
        in_pocket = site.atom_idx in pocket_donors

        if _occupied(graph, site.atom_idx):
            status = SiteStatus.OCCUPIED
        elif occluded is not None and occluded >= BLOCKED_OCCLUSION:
            status = SiteStatus.BLOCKED
        else:
            status = SiteStatus.OPEN

        state = SiteState(atom_idx=site.atom_idx, status=status,
                          buried_vol=occluded, in_pocket=in_pocket, ease=None)
        # The ease record reads `buried_vol` and `in_pocket` off the state, so the state
        # is built first and the ease attached second.  An occupied site is still scored:
        # "how easily would this have activated" stays a meaningful question about a bond
        # that already formed, and it is what a reaction edge needs to explain itself.
        ease = _ease_or_none(state, site)
        out.append(SiteState(atom_idx=site.atom_idx, status=status,
                             buried_vol=occluded, in_pocket=in_pocket, ease=ease,
                             pka=_pka_of(site)))
    return out


def _ease_or_none(state: SiteState, site: Any) -> EaseRecord | None:
    """Score the site, or record that the floor could not.

    A donor type with no row in the M1 table is a real gap (a perception rule exists that
    the descriptor table has not caught up with).  It is left unscored rather than given a
    guessed pKa — `EaseRecord.scored` is False and the caller can see why.
    """
    from mofsbu.descriptors.ease import EaseModelError

    try:
        return activation_ease(site, state=state)
    except EaseModelError:
        return None


def _pka_of(site: Any) -> float | None:
    from mofsbu.descriptors import tables

    try:
        return tables.donor(site.donor_type).pka
    except tables.UnknownDescriptor:
        return None
