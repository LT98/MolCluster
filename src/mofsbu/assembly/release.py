"""Release one bound ligand from a block — the step `join` runs backwards.

Bridging chemistry displaces co-ligands: a carboxylate cannot take a vertex a water still
holds.  `release` removes one ligand fragment, and the vertex its donor occupied becomes an
open vacancy pointing where that donor was, so the next join can use it.  The released
ligand comes back as its own block, so it can be stored and named on the edge as `leaving`
(`energy.reference.PRODUCT_SIDE_ROLES`): complex → complex′ + L, balanced by construction.

Rigid, like every assembly step: nothing else moves.  A donor bound to more than one metal
is a bridge, and dissolving a bridge is a different question — refused by name.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Sequence

import numpy as np

from mofsbu._types import MofsbuError
from mofsbu.assembly.join import BuildingBlock
from mofsbu.graph._types import EdgeType

__all__ = ["ReleaseRefused", "Released", "release"]


class ReleaseRefused(MofsbuError):
    """A release the bonding does not permit.  The message says which and why."""


@dataclass(frozen=True)
class Released:
    block: BuildingBlock               # what remains, with the freed vertex open
    ligand: BuildingBlock              # the released fragment, on its own
    atom_map: dict[int, int]           # parent atom -> remaining-block atom
    ligand_map: dict[int, int]         # parent atom -> released-block atom
    choice_vector: dict[str, Any]


def _remap_sites(sites, keep_map: dict[int, int]):
    return tuple(replace(s, atom_idx=keep_map[s.atom_idx]) for s in sites
                 if s.atom_idx in keep_map)


def release(block: BuildingBlock, atoms: Sequence[int], *,
            multiplicity: int = 1) -> Released:
    """Remove the ligand fragment containing `atoms` from `block`; open its vertex.

    `atoms` may name the whole fragment or any atom of it.  `multiplicity` is the released
    ligand's own (1 for water); the remaining block keeps the parent's, which is right
    when the ligand leaves closed-shell and is the caller's to state otherwise.
    """
    from mofsbu.sites.model import vacancy_sites
    from mofsbu.sites.state import refresh_state

    g = block.graph
    if block.geometry is None:
        raise ReleaseRefused("release needs coordinates: the freed vertex points where "
                               "the donor was, and a graph-only block has no 'where'")
    frag = next((tuple(f) for f in g.ligand_fragments() if set(atoms) & set(f)), None)
    if frag is None:
        raise ReleaseRefused(f"atoms {sorted(atoms)} are not in any ligand fragment; a "
                               f"metal is not released")
    bonds = [(a, m) for a in frag for m in g.neighbors(a, EdgeType.DATIVE) if g.is_metal(m)]
    if not bonds:
        raise ReleaseRefused(f"fragment {sorted(frag)} is not bound to a metal; there is "
                               f"nothing to release it from")
    if len(bonds) > 1:
        raise ReleaseRefused(
            f"fragment {sorted(frag)} holds {len(bonds)} metal contacts {bonds}: it is a "
            f"chelate or a bridge, and releasing it is not one vertex opening — refused")
    donor, metal = bonds[0]

    nodes = list(g.nodes())
    row = {n: k for k, n in enumerate(nodes)}
    xyz = np.asarray(block.geometry, dtype=float)
    keep = [n for n in nodes if n not in frag]
    lig = sorted(frag)
    frag_charge = sum(g.label(a).formal_charge for a in lig)

    rest = g.subgraph(keep).with_charge(g.net_charge() - frag_charge)
    rest = rest.with_multiplicity(g.multiplicity) if g.multiplicity is not None else rest
    rest.name = g.name
    keep_map = {old: new for new, old in enumerate(sorted(keep))}
    piece = g.subgraph(lig).with_charge(frag_charge).with_multiplicity(multiplicity)
    lig_map = {old: new for new, old in enumerate(lig)}

    # The freed vertex: on the metal, pointing at where the donor was.
    direction = xyz[row[donor]] - xyz[row[metal]]
    used = [s.slot for s in block.sites if s.is_vacancy and s.atom_idx == metal]
    vac = vacancy_sites(keep_map[metal], xyz[row[metal]], [direction])[0]
    vac = replace(vac, slot=(max(used) + 1) if used else 0)

    rest_xyz = xyz[[row[n] for n in sorted(keep)]]
    rest_sites = _remap_sites(block.sites, keep_map) + (vac,)
    rest_sym = [rest.label(n).element for n in rest.nodes()]
    states = refresh_state(rest_sites, rest_xyz, graph=rest, symbols=rest_sym)
    remaining = BuildingBlock(graph=rest, sites=rest_sites, geometry=rest_xyz,
                              state={(s.atom_idx, getattr(s, "slot", 0)): s for s in states})

    lig_xyz = xyz[[row[n] for n in lig]]
    lig_sites = _remap_sites(block.sites, lig_map)
    lig_sym = [piece.label(n).element for n in piece.nodes()]
    lstates = refresh_state(lig_sites, lig_xyz, graph=piece, symbols=lig_sym)
    released = BuildingBlock(graph=piece, sites=lig_sites, geometry=lig_xyz,
                             state={(s.atom_idx, getattr(s, "slot", 0)): s for s in lstates})

    cv = {"op": "release", "metal": g.label(metal).element,
          "donor_element": g.label(donor).element, "n_atoms_released": len(lig),
          "vertex_slot": vac.slot, "from_structure": block.structure_id}
    return Released(remaining, released, keep_map, lig_map, cv)
