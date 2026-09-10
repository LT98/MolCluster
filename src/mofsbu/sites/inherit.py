"""Carry sites across a join, instead of perceiving the product from scratch.

§6.3, last paragraph: "when a new structure is built by joining stored blocks, inherit
the parents' sites through the assembly's **atom map** — minus sites consumed by the new
bond, plus any newly exposed ones.  Perception is incremental, never from scratch."

Three reasons that is not merely an optimisation:

* **Provenance.**  The atom map is carried on every `reactions` edge (D8), so an inherited
  site knows which parent it descends from — "this axial site descends from the Cu that
  entered at step 1" is a query, not a reconstruction.
* **Stability.**  Re-perceiving a product can disagree with the parent it came from, which
  is why the registry watches for it at all (`registry.api.catalog_drift`).  Inheritance
  cannot drift, because it never asks the question twice.
* **Frames survive.**  A site's frame is the expensive, geometry-derived part (D13).  The
  parent's frame is still correct for every atom the join did not move, and re-deriving it
  from the product's coordinates would put a stochastic search back in the path the frame
  model exists to take it out of.

What this module deliberately does NOT do is decide which sites a join CONSUMES.  That is
`assembly.join`'s business — it is the thing that knows a bond was formed — and passing
the consumed set in keeps this function a pure remap.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Iterable, Sequence

from mofsbu._types import AtomMap
from mofsbu.sites.model import Site


def inherit_sites(parent_sites: Sequence[Site], atom_map: AtomMap,
                  consumed: Iterable[int] = ()) -> list[Site]:
    """Remap a parent's sites onto a product, dropping the ones the join used up.

    `atom_map` is the parent's convention throughout: `{parent_idx: child_idx}` (see
    `_types.AtomMap`), and `consumed` is given in PARENT indices too.  Mixing the two
    index spaces is the obvious way to get this wrong, so there is exactly one rule —
    everything you pass in is about the parent, everything you get back is about the
    child.

    A parent site whose atom has no entry in the map is dropped: the atom is not in the
    product, so neither is its site.  That covers a leaving group without needing a
    special case for it.
    """
    consumed_set = set(consumed)
    out: list[Site] = []
    for site in parent_sites:
        if site.atom_idx in consumed_set:
            continue                      # the join used this donor to make its bond
        child_idx = atom_map.get(site.atom_idx)
        if child_idx is None:
            continue                      # atom is not in the product at all
        out.append(replace(site, atom_idx=child_idx))
    return sorted(out, key=lambda s: s.atom_idx)


def merge_inherited(*groups: Sequence[Site]) -> list[Site]:
    """Combine the inherited site lists of several parents into one product list.

    Two parents cannot legitimately contribute a site on the same product atom — that
    would mean the join mapped two different parent atoms onto one child — so a collision
    is a broken atom map and is reported rather than silently resolved by ordering.
    """
    seen: dict[int, Site] = {}
    for group in groups:
        for site in group:
            if site.atom_idx in seen and seen[site.atom_idx] != site:
                raise ValueError(
                    f"two parents claim product atom {site.atom_idx}: "
                    f"{seen[site.atom_idx].donor_type} vs {site.donor_type}. The atom map "
                    f"is not injective, which makes the product's site model ambiguous.")
            seen[site.atom_idx] = site
    return sorted(seen.values(), key=lambda s: s.atom_idx)
