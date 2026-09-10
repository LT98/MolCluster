"""Activation ease — the C5 model, resolved as **D18**.

"Activation" = deprotonating a donor to expose it, and that donor's readiness to form the
next bond (§6.4).  The honest quantity is a set of NAMED COMPONENTS at whatever fidelity
is available; the scalar is a sort key layered on top and never a substitute (D6).

**What D18 settles, and why each part is the way it is.**

* **The floor is zero-QM and unconditional.**  Every perceived site gets a record from the
  M1 tables alone, so no structure carries an empty ease field.  A site the model cannot
  score is reported as *unscored*, never as *hard* — a missing number that renders as a
  low one is the failure this whole layer exists to avoid.
* **Absent components stay absent.**  The scalar is a weighted mean over the components
  actually PRESENT, renormalised by their weights.  Nothing is defaulted to 0, or to 0.5,
  or to "neutral": a defaulted component is a number the model did not produce, and it
  would be indistinguishable from one it did.
* **`confidence` answers "how much of the model ran", not "how right is it".**  It is
  `coverage x sharpness`: how much of the weighted model was available, times how sharp
  the table row behind it is.  A pKa carrying a 4-unit spread (`hydrohalide_X`) cannot
  produce a confident anything, and the table's own `pka_sigma` column is what says so.
* **`provisional` marks sites whose TABLE VALUE IS THE WRONG QUESTION**, not sites that
  merely scored badly: a donor in a chelate pocket has a pKa shifted by its neighbour
  (the anthrarufin peri-OH, whose value the quinone H-bond moves), and a wide-sigma row is
  the table saying it is a class estimate.  Those are the sites worth spending QM on, and
  flagging them is the floor's most useful output — it is a work list, not a verdict.

**The rung above the floor is MACE-OMOL-0, and that is new.**  Deprotonation is a
charge-changing quantity: `A-H -> A(-) + H(+)`.  A charge-blind potential cannot see it
even in principle, so before OMOL-0 the ladder for the primary component ran
`table -> (nothing) -> xTB` — the ML rung, the one that makes screening affordable, could
not serve the component the model is mostly made of.  That gap is why the floor had to
carry the whole model and why C5 stayed open this long.  `deprotonation_energy()` below is
that rung, and it takes any charge-aware backend rather than naming one.

**C7 is still open and is not resolved here.**  `hsab_match` still raises, and
`activation_ease(partner=...)` refuses rather than quietly returning the partner-free
number under a partner-shaped signature.  §6.4 is explicit that the intrinsic scalar is
the partner-free default with the per-partner correction *layered on*, so the floor is
well-defined without C7 — but a caller who asked for a partner and silently got no partner
term would have no way to tell.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from mofsbu._types import Fidelity, MethodSpec, MofsbuError
from mofsbu.assembly.join import NotBuiltYet
from mofsbu.descriptors import tables
from mofsbu.descriptors.tables import DonorDescriptor, MetalDescriptor
from mofsbu.versions import ALGO_VERSIONS


class EaseModelError(MofsbuError):
    """The ease model was asked for something it does not define."""


# ── the constants D18 fixes ──────────────────────────────────────────────────

#: Weights over the four named components of §6.4.  Deprotonation dominates because §6.4
#: calls it the primary component and because it is the only one the floor can always
#: produce.  These are the DENOMINATOR too: a record missing `marginal_de` is renormalised
#: over the weights it does have, so weights are ratios between components and never a
#: hidden penalty for a component that was not computed.
WEIGHTS: dict[str, float] = {
    "deprotonation": 0.55,
    "steric":        0.25,
    "electronic":    0.10,
    "marginal_de":   0.10,
}

#: pKa at which the deprotonation sub-score is 0.5, and the width of the logistic.
#: Anchored, not tuned: these two numbers place the donor classes this project actually
#: builds with where chemistry says they belong under coordination-synthesis conditions
#: (DMF/water, modulator or mild base).  Carboxylate (4.76) -> 0.85 "goes readily";
#: phenolate (9.99) -> 0.50 "needs help"; alkoxide (15.9) -> 0.12 "not without strong
#: base".  If the conditions model ever becomes real (see `Conditions`), CENTRE is the
#: parameter it moves, which is why it is one named constant and not a buried literal.
PKA_CENTRE = 10.0
PKA_WIDTH = 3.0

#: `confidence` decays with the table's own honest spread.  sigma=0.3 (halide_F, a
#: measured parent acid) -> 0.90; sigma=0.8 (carboxylate, a real class) -> 0.77;
#: sigma=2.5 (azolate, a wide class) -> 0.43; sigma=4.0 (hydrohalide_X, HF to HI) -> 0.26.
SIGMA_SCALE = 3.0

#: At or above this spread the row is a class estimate rather than a value, and the site
#: self-flags for promotion regardless of what it scored.
PROVISIONAL_SIGMA = 2.0

HEURISTIC_METHOD = MethodSpec(
    code="heuristic",
    code_version=ALGO_VERSIONS["ease_model"],
    method="pka-table-floor",
    extras={"algo": ALGO_VERSIONS["ease_model"],
            "descriptor_tables": ALGO_VERSIONS["descriptor_tables"]},
)


@dataclass(frozen=True)
class Conditions:
    """The medium an activation is asked about.  A placeholder with a settled name.

    Ground rule 7: the parameter exists in the signature because §6.4 defines ease as
    `activation_ease(site, partner, conditions)` and callers are written against that.
    What it does NOT have is a model — a solvent-dependent pKa shift is its own piece of
    chemistry and inventing one here would put a number in the record that nothing
    measured.  Passing a non-default `Conditions` therefore raises rather than being
    accepted and ignored.
    """

    solvent: str | None = None
    temperature_k: float = 298.15

    @property
    def is_default(self) -> bool:
        return self.solvent is None and self.temperature_k == 298.15


@dataclass(frozen=True)
class EaseRecord:
    """Named components first, scalar second.

    §6.6: "activation ease as a single scalar is lossy — the named-component design
    mitigates".  The components are the record; `scalar` is a sort key.
    """

    components: dict[str, float]      # deprotonation / electronic / steric / marginal_dE
    scalar: float
    fidelity: Fidelity
    confidence: float
    provisional: bool
    method: MethodSpec
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def coverage(self) -> float:
        """Fraction of the weighted model that actually produced a number."""
        got = sum(WEIGHTS[k] for k in self.components if k in WEIGHTS)
        return got / sum(WEIGHTS.values())

    @property
    def scored(self) -> bool:
        """False when nothing could be computed — distinct from 'scored badly'."""
        return bool(self.components)


# ── the components ───────────────────────────────────────────────────────────


def deprotonation_score(donor: DonorDescriptor) -> tuple[float, str]:
    """pKa -> ease of removing the proton, in [0, 1].  Returns (score, note).

    A neutral donor is the interesting case.  Its `pka` is EMPTY in the table and the
    table's header says why: there is no acid whose deprotonation exposes it, so there is
    no activation step to be easy or hard.  That scores 1.0 — already available — and it
    is emphatically not a missing value.  Scoring it 0, or dropping the component, would
    rank a pyridyl N below a carboxylic acid at the one thing pyridine is good at.
    """
    if donor.pka is None:
        return 1.0, f"{donor.donor_type}: neutral donor, no activation step"
    score = 1.0 / (1.0 + math.exp((donor.pka - PKA_CENTRE) / PKA_WIDTH))
    return score, f"{donor.donor_type}: pKa {donor.pka:.2f}"


def sharpness(donor: DonorDescriptor) -> float:
    """How much the table trusts its own pKa, in (0, 1].

    A neutral donor has no pKa to be uncertain about, so its sharpness is 1.0 — the
    absence is a fact about the chemistry, not a gap in the table.
    """
    if donor.pka is None or donor.pka_sigma is None:
        return 1.0
    return math.exp(-abs(donor.pka_sigma) / SIGMA_SCALE)


def occlusion(coords: Any, site_idx: int, axis: Any, *,
              radius: float = 3.5, n_rays: int = 64, half_angle_deg: float = 75.0,
              elements: Any = None) -> float:
    """Fraction of the approach cone that is blocked, in [0, 1].

    §6.4 lists steric accessibility as "cone/solid angle or buried volume %V_bur —
    geometric, cheap".  This is the cone form, and it is cast along the site's own OUTWARD
    AXIS rather than over a whole sphere, which is the frame model (D13) earning its keep
    a second time: a buried-volume sphere centred on the donor counts the ligand the donor
    is attached to as blocking, and every donor then looks equally crowded.  What matters
    is whether a metal can get in ALONG THE DIRECTION THE LONE PAIR POINTS.
    """
    import numpy as np

    coords = np.asarray(coords, dtype=float)
    origin = coords[site_idx]
    axis = np.asarray(axis, dtype=float)
    axis = axis / max(float(np.linalg.norm(axis)), 1e-9)

    # Everything except the donor itself.  Hydrogens included: an H sitting in the
    # approach path really does block it, and this is the check that notices the proton
    # a deprotonation would remove.
    others = [i for i in range(len(coords)) if i != site_idx]
    if not others:
        return 0.0
    blockers = coords[others]
    radii = np.full(len(others), 1.5)
    if elements is not None:
        radii = np.array([_VDW.get(str(elements[i]), 1.5) for i in others])

    # A deterministic Fibonacci cap around the axis — same rays every call, so the number
    # is reproducible and two runs of the same geometry cannot disagree.
    blocked = 0
    cos_max = math.cos(math.radians(half_angle_deg))
    perp = _perp(axis)
    perp2 = np.cross(axis, perp)
    for k in range(n_rays):
        # cos(theta) uniform in [cos_max, 1] gives rays uniform over the cap's solid angle
        cos_t = cos_max + (1.0 - cos_max) * (k + 0.5) / n_rays
        sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
        phi = 2.399963229728653 * k                      # golden angle, in radians
        ray = cos_t * axis + sin_t * (math.cos(phi) * perp + math.sin(phi) * perp2)
        rel = blockers - origin
        along = rel @ ray
        near = (along > 0) & (along < radius)
        if not near.any():
            continue
        perp_dist = np.linalg.norm(rel[near] - np.outer(along[near], ray), axis=1)
        if (perp_dist < radii[near]).any():
            blocked += 1
    return blocked / n_rays


_VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
        "Cl": 1.75, "Br": 1.85, "I": 1.98, "Se": 1.90, "B": 1.92, "Si": 2.10}


def _perp(v: Any) -> Any:
    import numpy as np

    trial = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    out = trial - v * float(v @ trial)
    return out / max(float(np.linalg.norm(out)), 1e-9)


# ── the record ───────────────────────────────────────────────────────────────


def activation_ease(site: Any, *, partner: MetalDescriptor | None = None,
                    conditions: Conditions | None = None,
                    state: Any | None = None) -> EaseRecord:
    """How easily does this site activate?  The D18 floor.

    `site` is a `sites.model.Site`.  `state` may carry a per-geometry steric number
    (`buried_vol`, from `sites.state.refresh_state`); without one the record is
    deprotonation-only and says so through `coverage`, rather than pretending a
    geometry-free answer covered the geometric component.
    """
    if partner is not None:
        raise NotBuiltYet(
            "activation_ease(partner=...) needs C7, which is still an open checkpoint. "
            "D18 resolved C5 — the partner-FREE intrinsic floor — only. The per-partner "
            "correction layers on top of this record (§6.4) and is not folded into it; "
            "returning the intrinsic number under a partner-shaped call would be "
            "indistinguishable from a partner-aware one.")
    if conditions is not None and not conditions.is_default:
        raise NotBuiltYet(
            "activation_ease(conditions=...) has a settled signature and no model. A "
            "solvent-shifted pKa is its own chemistry; the table's values are aqueous, "
            "25 C, and saying otherwise would put an unmeasured number in the record.")

    donor_type = getattr(site, "donor_type", None)
    if donor_type is None:
        raise EaseModelError(f"{site!r} is not a Site: no donor_type")
    try:
        donor = tables.donor(donor_type)
    except tables.UnknownDescriptor as exc:
        # A perceived donor with no table row is a real gap and is reported as one. The
        # alternative — an ease record built on a guessed pKa — is a number that looks
        # exactly like the ones that came from the curated table.
        raise EaseModelError(
            f"no descriptor row for donor_type {donor_type!r}; the floor cannot score a "
            f"donor the M1 table has never heard of") from exc

    components: dict[str, float] = {}
    notes: list[str] = []

    score, note = deprotonation_score(donor)
    components["deprotonation"] = score
    notes.append(note)

    buried = getattr(state, "buried_vol", None) if state is not None else None
    if buried is not None:
        components["steric"] = max(0.0, min(1.0, 1.0 - float(buried)))
        notes.append(f"steric from occlusion {float(buried):.2f}")

    # `electronic` needs a partner (C7) and `marginal_de` needs energies; neither is
    # available at the floor, and both are therefore absent rather than zero.

    total = sum(WEIGHTS[k] for k in components)
    scalar = sum(WEIGHTS[k] * v for k, v in components.items()) / total if total else 0.0

    # A donor with no pKa has nothing for the table to be wrong ABOUT, so it never
    # self-flags however it is placed: promoting a pyridyl N to xTB to refine a
    # deprotonation it does not undergo is spending QM to learn nothing.
    activatable = donor.pka is not None
    in_pocket = (bool(getattr(state, "in_pocket", False)) and activatable
                 if state is not None else False)
    wide = (activatable and donor.pka_sigma is not None
            and abs(donor.pka_sigma) >= PROVISIONAL_SIGMA)
    if in_pocket:
        notes.append("neighbouring group shifts this pKa: table value is the wrong question")
    if wide:
        notes.append(f"pKa spread {donor.pka_sigma} is a class estimate")

    coverage = total / sum(WEIGHTS.values())
    return EaseRecord(
        components=components,
        scalar=scalar,
        fidelity=Fidelity.HEURISTIC,
        confidence=coverage * sharpness(donor),
        provisional=bool(in_pocket or wide),
        method=HEURISTIC_METHOD,
        notes=tuple(notes),
    )


# ── the rung above the floor ─────────────────────────────────────────────────


def deprotonation_energy(backend: Any, symbols: Any, positions: Any, *,
                         charge: int, multiplicity: int,
                         deprotonated_positions: Any, deprotonated_symbols: Any,
                         deprotonated_multiplicity: int) -> float:
    """ΔE of `A-H -> A(-) + H(+)`, in eV, from a charge-aware backend.

    The reason this can exist at the ML rung at all: a deprotonation CHANGES THE CHARGE,
    so a charge-blind potential is not merely inaccurate about it, it is blind to it —
    both sides look like the same collection of elements minus an H.  MACE-OMOL-0 takes
    total charge and multiplicity as inputs, so it can be asked; MACE-MP-0 cannot, and is
    refused here rather than returning the difference of two numbers that never saw the
    charge.

    The bare proton contributes no electronic energy (it has no electrons), so this is
    the raw electronic ΔE and NOT a ΔG: no ZPE, no thermal correction, no solvation of
    H(+).  It is therefore comparable BETWEEN SITES computed the same way and is not an
    absolute acidity — which is exactly the use the ease model puts it to, and the reason
    it is returned as a component rather than converted into a pKa it cannot support.
    """
    if not getattr(backend, "charge_aware", False):
        raise EaseModelError(
            f"{getattr(backend, 'name', backend)!r} is charge-blind; a deprotonation "
            f"energy is a charge CHANGE and this backend cannot see one. Use a "
            f"charge-aware backend (MACE-OMOL-0, xTB).")
    protonated = backend.single_point(symbols, positions, charge=charge,
                                      multiplicity=multiplicity)
    anion = backend.single_point(deprotonated_symbols, deprotonated_positions,
                                 charge=charge - 1,
                                 multiplicity=deprotonated_multiplicity)
    return float(anion.energy - protonated.energy)


def hsab_match(donor: DonorDescriptor, metal: MetalDescriptor) -> float:
    """Hard/soft compatibility of a donor and an ion, in [0, 1].  C7's factorized core."""
    raise NotBuiltYet(
        "descriptors.ease.hsab_match — C7 is an open checkpoint. The factorization is the "
        "leaning (compute at query time, never store a site x metal matrix), but the "
        "match function itself is unratified. D18 resolved C5 and deliberately did not "
        "reach for this one: the floor is partner-free by design (§6.4).")
