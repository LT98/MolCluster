"""Activation ease — the C5 model.  NOT IMPLEMENTED.

Ground rule 8: the signatures are settled here so `sites`, the runner and the UI can be
written against them, and the bodies raise.  C5 (the heuristic-tier floor) and C7 (the
partner-dependence factorization) are both still OPEN CHECKPOINTS — deliberately left
open when the descriptor tables landed, because the tables are facts and this is a
policy, and shipping the policy before it is ratified is how a leaning becomes a
decision by accident.

What already exists, so this is a small job when the call is made:

* `descriptors.tables.donor()` gives pKa + spread + HSAB class + denticity, and says via
  `needs_activation` whether there is a proton to remove at all.
* `descriptors.tables.metal()` gives the partner side: HSAB class, ionic radius,
  preferred CN, exchange lability.
* `sites.model.Site` carries the perceived donor and its frame; `site_state` is where an
  ease value would be written, and is still unpopulated (M4, second half).

What C5 still has to decide, and what a reader should NOT assume from this file:

* how the named components combine into `scalar` — the design doc is explicit that a
  single scalar is lossy, so the components dict is the primary record and the scalar is
  a convenience;
* what `confidence` means numerically when a pKa carries a 2.5-unit spread;
* which donors self-flag `provisional` — the leaning is "in a chelate or H-bond pocket",
  which `sites.model.chelate_pockets` can already answer.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from mofsbu._types import Fidelity, MethodSpec
from mofsbu.assembly.join import NotBuiltYet
from mofsbu.descriptors.tables import DonorDescriptor, MetalDescriptor


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


def activation_ease(site: Any, *, partner: MetalDescriptor | None = None,
                    conditions: Any | None = None, state: Any | None = None) -> EaseRecord:
    """How easily does this site activate, optionally against a specific metal?"""
    raise NotBuiltYet(
        "descriptors.ease.activation_ease — C5 is an open checkpoint. The descriptor "
        "tables it needs are built (M1); what is missing is the decision about how the "
        "components combine and what confidence means, which is a call to make, not code "
        "to write.")


def hsab_match(donor: DonorDescriptor, metal: MetalDescriptor) -> float:
    """Hard/soft compatibility of a donor and an ion, in [0, 1].  C7's factorized core."""
    raise NotBuiltYet(
        "descriptors.ease.hsab_match — C7 is an open checkpoint. The factorization is the "
        "leaning (compute at query time, never store a site x metal matrix), but the "
        "match function itself is unratified.")
