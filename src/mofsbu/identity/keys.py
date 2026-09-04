"""The composite L0-L3 identity key (design doc §4).

    BLOCK_ID = L0_composition / L1_graph_hash / L2_isomer_tag / L3_conformer_id

L0 and L1 are live.  L2 and L3 ship with their final signatures and stub bodies —
they land in the schema now (M3) and are filled in at M5, so the column and the call
site never have to change.
"""
from __future__ import annotations

from typing import Any

from mofsbu.graph.canon import certificate_digest, wl_hash
from mofsbu.graph._types import TypedGraph
from mofsbu._types import AmbiguousSpecError
from mofsbu.versions import ALGO_VERSIONS

L2_UNSET = ""
L3_UNSET = ""


def hill_formula(counts: dict[str, int]) -> str:
    """Hill order: C, then H, then everything else alphabetically."""
    parts: list[str] = []
    for el in ("C", "H"):
        if counts.get(el):
            parts.append(f"{el}{counts[el]}" if counts[el] > 1 else el)
    for el in sorted(k for k in counts if k not in ("C", "H")):
        parts.append(f"{el}{counts[el]}" if counts[el] > 1 else el)
    return "".join(parts)


def l0_composition(g: TypedGraph) -> str:
    """Composition + charge + spin + per-centre labels.

    Format:  ``Cu2_C4H4O8_q0_s3|Cu(+2,hs),Cu(+2,hs)``
             ``C9H3O6_q-3_s1``            (metal-free: no centre block)

    Per-centre oxidation state and spin are part of L0, not a single global pair
    (D12): a mixed-valence Fe(II)/Fe(III)/Fe(III) trimer is a different structure
    from its all-Fe(III) analogue even though formula, charge and multiplicity may
    coincide.

    Raises `AmbiguousSpecError` if the multiplicity was never specified — an
    identity key may not contain a guessed spin state (ground rule 5).
    """
    if g.multiplicity is None:
        raise AmbiguousSpecError(
            f"{g!r} has no multiplicity; set it explicitly (TypedGraph(multiplicity=...) "
            "or .with_multiplicity(...)) before asking for an identity key"
        )
    charge = g.net_charge()      # raises AmbiguousSpecError if unset
    counts = g.element_counts()
    metal_idxs = g.metals()
    metal_counts: dict[str, int] = {}
    for i in metal_idxs:
        el = g.label(i).element
        metal_counts[el] = metal_counts.get(el, 0) + 1
    organic = {k: v for k, v in counts.items() if k not in metal_counts}

    chunks: list[str] = []
    if metal_counts:
        chunks.append("".join(f"{el}{metal_counts[el]}" for el in sorted(metal_counts)))
    if organic:
        chunks.append(hill_formula(organic))
    chunks.append(f"q{charge}")
    chunks.append(f"s{g.multiplicity}")
    base = "_".join(chunks)

    if not metal_idxs:
        return base
    centres = sorted(
        "{}({},{})".format(
            g.label(i).element,
            "?" if g.label(i).oxidation_state is None else f"{g.label(i).oxidation_state:+d}",
            g.label(i).spin_class or "?",
        )
        for i in metal_idxs
    )
    return f"{base}|{','.join(centres)}"


def l1_graph_hash(g: TypedGraph) -> str:
    """Canonical hash of the typed graph.  Build-order invariant (D3).

    Exact: computed from the canonical certificate, not from Weisfeiler-Lehman.
    WL merges mu2-bridging with chelating carboxylate (see `wl_hash`), which is a
    distinction this project is built to draw, so WL is kept only as the fast bucket
    index below.  Proposed as D16.

    Returns the bare hex digest; the recipe version is pinned in
    `ALGO_VERSIONS['l1_certificate']` and stored in its own column, never
    concatenated into the key.
    """
    return certificate_digest(g)


def wl_index(g: TypedGraph) -> str:
    """Cheap bucket key for narrowing candidates before the exact comparison."""
    return wl_hash(g)


def l2_isomer_tag(g: TypedGraph, geom: Any | None = None) -> str:
    """Configurational isomer tag: cis/trans, fac/mer, Delta/Lambda (D10).

    STUB until M5.  Signature is final: it takes the graph plus an optional geometry,
    because the distinction L2 draws is spatial and invisible to the graph alone.
    """
    return L2_UNSET


def l3_conformer_id(choice_vector: Any | None = None, geom: Any | None = None) -> str:
    """Conformer id — provenance-primary, geometry-verifier (D11).

    STUB until M5.  The label is the construction choice-vector; geometric clustering
    only collapses stochastic duplicates and reconciles divergence/convergence.
    """
    return L3_UNSET


def block_id(l0: str, l1: str, l2: str = L2_UNSET, l3: str = L3_UNSET) -> str:
    """`L0/L1/L2/L3`, trailing empty levels kept so the arity is always readable."""
    return "/".join((l0, l1, l2, l3))


def identity(g: TypedGraph, *, geom: Any | None = None) -> dict[str, str]:
    """Everything the registry needs to key a structure row."""
    l0 = l0_composition(g)
    l1 = l1_graph_hash(g)
    l2 = l2_isomer_tag(g, geom)
    return {
        "l0": l0,
        "l1": l1,
        "l2": l2,
        "wl_index": wl_index(g),
        "block_id": block_id(l0, l1, l2),
        "algo_l0": ALGO_VERSIONS["l0_composition"],
        "algo_l1": ALGO_VERSIONS["l1_certificate"],
        "algo_l2": ALGO_VERSIONS["l2_isomer_tag"],
        "algo_canon": ALGO_VERSIONS["canonical_order"],
    }
