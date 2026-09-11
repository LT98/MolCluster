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
    # 2: the MACE backends became a family (MP-0 / OMOL-0).  A method row now records
    #    `training_set`, and `spin_blind` alongside `charge_blind`, so a stored ML number
    #    says which foundation model made it instead of only "mace".  Old rows keep
    #    algo=1 and are never re-labelled (ground rule 6); they are simply not the same
    #    method as anything computed from here on, which is correct — they are not.
    "energy_backends":  "2",      # backend protocol + how a MethodSpec is filled in
    "reference_scheme": "balanced1",   # reaction-balanced energies (M7)
    "spin_convention":  "hs1",    # high_spin_multiplicity: d-count table + charge
    "descriptor_tables": "1",     # curated donor table + generated metal table (M1, C8)
    # The C5 floor as ratified in D18: which components exist, how they are weighted,
    # what confidence and provisional mean.  Bumping this is what marks stored ease rows
    # as products of a different policy — the numbers are only comparable within a
    # version, because the weights ARE the model.
    "ease_model":       "floor1",  # pKa-table floor + cone occlusion (M4, C5/D18)
    "site_state":       "1",       # status rule + buried-volume convention (M4)
    # 2: pinning a donor and its outward axis leaves two rotations undetermined — the
    #    ligand's spin about the M-L axis, and the metal's swing out of the donor's
    #    plane — and `1` left both wherever the alignment arithmetic dropped them.  They
    #    are now searched over fixed grids and the winning indices are recorded, so the
    #    coordinates of every monodentate placement differ from a `1` one.  Old rows keep
    #    algo=1 and are never re-labelled (ground rule 6): they are not a worse version
    #    of this recipe, they are a different one, and only comparable among themselves.
    #    The grid sizes and the out-of-plane cap are part of the recipe, because an index
    #    only means anything against a grid of a known size.
    "placement":        "2",       # how a ligand is oriented on a coordination vertex
}


def version_block() -> str:
    return "\n".join(f"{k:18s} {v}" for k, v in sorted(ALGO_VERSIONS.items()))
