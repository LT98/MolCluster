"""Pinned algorithm versions.

Ground rule 6 (docs/PLAN_implementation.md §0): every stored value carries the version of
the recipe that produced it. Bump a version here when its algorithm changes; never
silently re-label existing rows.
"""
from __future__ import annotations

ALGO_VERSIONS: dict[str, str] = {
    "graph_schema":    "1",       # NodeLabel/EdgeType layout + what enters a hash
    "l0_composition":  "1",       # L0 string format
    "l1_certificate":  "cert1",   # sha256 over the IR canonical certificate — the L1 key
    "wl_index":        "wl3-nx",  # networkx WL, 3 iterations — fast bucket index only
    "canonical_order": "ir1",     # individualisation-refinement canonical labelling
    "l2_isomer_tag":   "0-stub",  # M5
    "l3_conformer_id": "0-stub",  # M5
    "energy_backends":  "1",      # backend protocol + how a MethodSpec is filled in
    "reference_scheme": "balanced1",   # reaction-balanced energies (M7)
    "spin_convention":  "hs1",    # high_spin_multiplicity: d-count table + charge
}


def version_block() -> str:
    return "\n".join(f"{k:18s} {v}" for k, v in sorted(ALGO_VERSIONS.items()))
