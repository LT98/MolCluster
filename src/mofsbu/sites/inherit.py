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
                  consumed: Iterable[int | tuple[int, int]] = ()) -> list[Site]:
    """Remap a parent's sites onto a product, dropping the ones the join used up.

    `atom_map` is the parent's convention throughout: `{parent_idx: child_idx}` (see
    `_types.AtomMap`), and `consumed` is given in PARENT indices too.  Mixing the two
    index spaces is the obvious way to get this wrong, so there is exactly one rule —
    everything you pass in is about the parent, everything you get back is about the
    child.

    **`consumed` takes an atom index OR an `(atom_idx, slot)` pair, and on a metal it has
    to be the pair.**  A metal carries several vacant coordination vertices on ONE atom —
    that is why `Site.slot` exists and why `site_catalog`'s UNIQUE key was widened to
    include it — so consuming a vacancy by atom index alone would drop every OTHER
    vacancy on the same metal at the same time.  A four-coordinate metal would lose three
    open vertices to its first join and read as coordinatively saturated after one ligand.
    An int still means "this atom, whatever slot", which is right for a donor: a donor is
    always slot 0 and a consumed donor is consumed entire.

    A parent site whose atom has no entry in the map is dropped: the atom is not in the
    product, so neither is its site.  That covers a leaving group without needing a
    special case for it.
    """
    consumed_atoms = {c for c in consumed if isinstance(c, int)}
    consumed_slots = {tuple(c) for c in consumed if not isinstance(c, int)}
    out: list[Site] = []
    for site in parent_sites:
        if site.atom_idx in consumed_atoms:
            continue                      # the join used this donor to make its bond
        if (site.atom_idx, getattr(site, "slot", 0)) in consumed_slots:
            continue                      # ...or this ONE of the metal's vertices
        child_idx = atom_map.get(site.atom_idx)
        if child_idx is None:
            continue                      # atom is not in the product at all
        out.append(replace(site, atom_idx=child_idx))
    return sorted(out, key=lambda s: (s.atom_idx, getattr(s, "slot", 0)))


def merge_inherited(*groups: Sequence[Site]) -> list[Site]:
    """Combine the inherited site lists of several parents into one product list.

    Two parents cannot legitimately contribute a site at the same product **(atom, slot)**
    — that would mean the join mapped two different parent atoms onto one child — so a
    collision is a broken atom map and is reported rather than silently resolved by
    ordering.

    Keyed on the pair and not on the atom, for the same reason `consumed` accepts a pair:
    a metal's several vacancies all sit on one atom, and an atom-keyed check reads them as
    one parent contradicting itself.  It would have refused to merge any metal block
    carrying more than one open vertex — which is every metal a join is interesting on.
    """
    seen: dict[tuple[int, int], Site] = {}
    for group in groups:
        for site in group:
            key = (site.atom_idx, getattr(site, "slot", 0))
            if key in seen and seen[key] != site:
                raise ValueError(
                    f"two parents claim product atom {site.atom_idx} slot {key[1]}: "
                    f"{seen[key].donor_type or seen[key].role} vs "
                    f"{site.donor_type or site.role}. The atom map is not injective, "
                    f"which makes the product's site model ambiguous.")
            seen[key] = site
    return sorted(seen.values(), key=lambda s: (s.atom_idx, getattr(s, "slot", 0)))
