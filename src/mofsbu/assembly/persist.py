"""Storing what assembly builds — the seam between M5 and the registry (M5/S7).

Nothing here writes to the database.  Every insert goes through `registry.api`, which is
the only write surface (ground rule 1); this module is the orchestration that knows the
ORDER those calls have to happen in and what assembly knows that the registry does not.

Three things it supplies that the registry cannot derive for itself:

* **The L2 tag.**  `put_structure` calls `l2_isomer_tag(g)` with no coordinates, because a
  structure row is created *before* its geometry exists.  Only the builder is holding the
  coordinates at insert time, which is why `put_structure` has always taken `l2=` — the
  parameter was the right interface before there was anything to put in it.
* **The sites, by inheritance.**  A product's donors are carried through the atom map, not
  re-perceived: `sites.perception` treats a metal as an ordinary heavy neighbour, so
  re-perceiving an assembled complex returns NO donors at all (`docs/ISSUES.md` 4).  The
  runner's `_record_sites` perceives because it builds from a molecule; this path must not.
* **The provenance.**  Which blocks went in, the atom map, the choice-vector digest, and
  the depth — D8's edge, carrying what D2 keeps off the node.

**Idempotent by construction, not by checking.**  Storing the same construction twice
yields one `structures` row and two `reactions` edges, because that is what
`put_structure` does with a repeated identity.  That is the D2 claim and S7's first exit
gate, and this module gets it for free by not trying to be clever about it.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np

from mofsbu._types import Fidelity, MethodSpec
from mofsbu.assembly.choice import ChoiceVector
from mofsbu.identity.keys import l2_isomer_tag
from mofsbu.registry.api import (
    Provenance, Put, canonical_map, put_geometry, put_site_state, put_sites, put_structure)
from mofsbu.registry.db import Registry
from mofsbu.sites.state import refresh_state
from mofsbu.versions import ALGO_VERSIONS

#: What produced these coordinates.  A construct is not a relaxation and says so: RAW is
#: the bottom geometry rung, and an assembled product arrives there until something spends
#: compute on it.
CONSTRUCT = MethodSpec(code="construct", code_version=ALGO_VERSIONS["placement"],
                       method="frame-aligned-join",
                       extras={"choice_vector": ALGO_VERSIONS["choice_vector"]})


@dataclass(frozen=True)
class Stored:
    """What came back from storing one construction."""

    structure_id: int
    geometry_id: int
    reaction_id: int | None
    structure_created: bool
    geometry_created: bool
    l2: str
    choice_digest: str
    n_sites: int
    n_states: int


def to_xyz(graph: Any, coords: Any, comment: str = "") -> str:
    """`.xyz` text in graph node order — the order everything else keys against."""
    coords = np.asarray(coords, dtype=float)
    nodes = graph.nodes()
    lines = [str(len(nodes)), comment.replace("\n", " ")]
    for row, i in zip(coords, nodes):
        lines.append(f"{graph.label(i).element} "
                     f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}")
    return "\n".join(lines) + "\n"


def store_block(reg: Registry, block: Any, *, choice_vector: dict | ChoiceVector | None = None,
                reagent_ids: Iterable[int] = (), depth: int | None = None,
                note: str = "", tags: Sequence[str] = (), seed: int | None = None,
                fidelity: Fidelity = Fidelity.RAW, method: MethodSpec = CONSTRUCT,
                energy: float | None = None) -> Stored:
    """Persist one assembled `BuildingBlock`: identity, geometry, sites, provenance.

    The order is not arbitrary.  The structure row has to exist before a geometry can hang
    off it; `put_sites` keys against the canonical order frozen at that insert, so it
    cannot run earlier; and `site_state` references both a site row and a geometry row, so
    it runs last.  Getting this wrong does not fail loudly — it writes state against the
    wrong geometry — which is why it is written down once here instead of at each call
    site.
    """
    graph, coords = block.graph, block.geometry
    cv = ChoiceVector.coerce(choice_vector)
    tag = l2_isomer_tag(graph, coords) if coords is not None else ""

    prov = Provenance(kind="assembly", reagent_ids=tuple(reagent_ids), note=note,
                      choice_vector_digest=cv.digest if cv else None, depth=depth)
    put: Put = put_structure(reg, graph, l2=tag, provenance=prov, tags=tuple(tags))

    if coords is None:
        raise ValueError(
            "this block has no geometry, so there is nothing to store as one. A "
            "graph-only product is a real thing (a two-centre join stops there until M6) "
            "but it is a structure row without a geometry, and storing it as though it "
            "had coordinates would invent them. Use put_structure directly if that is "
            "what you want.")

    geom = put_geometry(reg, put.id, to_xyz(graph, coords, note or graph.name),
                        fidelity=fidelity, method=method, energy=energy,
                        choice_vector=cv.data if cv else None,
                        seed=seed if seed is not None else (cv.seed if cv else None))

    # Sites are INHERITED, never re-perceived here — see the module docstring.
    sites = list(block.sites)
    n_sites = put_sites(reg, put.id, sites, algo=f"inherited/{ALGO_VERSIONS['graph_schema']}")
    symbols = [graph.label(i).element for i in graph.nodes()]
    states = block.state
    if states is None:
        states = {(s.atom_idx, getattr(s, "slot", 0)): s
                  for s in refresh_state(sites, coords, graph=graph, symbols=symbols,
                                         fidelity=fidelity)}
    n_states = put_site_state(reg, put.id, geom.id, list(states.values()), fidelity=fidelity)

    reaction_id = _edge_for(reg, put.id, prov)
    reg.conn.commit()
    return Stored(put.id, geom.id, reaction_id, put.created, geom.created, tag,
                  cv.digest if cv else "", n_sites, n_states)


def _edge_for(reg: Registry, structure_id: int, prov: Provenance) -> int | None:
    """The provenance edge `put_structure` just wrote, so the caller can point at it.

    `put_structure` records the edge itself; this reads back the newest one rather than
    writing a second, because two edges for one act of building would make
    `n_incoming_routes` count the storing rather than the routes.
    """
    row = reg.conn.execute(
        "SELECT id FROM reactions WHERE product_structure_id=? ORDER BY id DESC LIMIT 1",
        (structure_id,)).fetchone()
    return int(row["id"]) if row else None


def store_construction(reg: Registry, construction: Any, *, reagent_ids: Iterable[int] = (),
                       tags: Sequence[str] = (), **kw: Any) -> Stored:
    """Persist a `construct.Construction`, path and all.

    The choice vector stored on the geometry is the WHOLE construction's, not the last
    join's, and `depth` is the number of steps — so a query can ask "what was built in
    three moves" and "rebuild this exactly" of the same row.
    """
    return store_block(reg, construction.block,
                       choice_vector=construction.choice_vector,
                       reagent_ids=reagent_ids, depth=len(construction.steps),
                       tags=tags, **kw)


def canonical_sites(reg: Registry, structure_id: int, sites: Sequence[Any]) -> dict[int, Any]:
    """`{canonical index: site}` for an assembled block — the registry's own coordinates.

    Useful when comparing what was stored against what was built: the catalog is keyed by
    canonical index and a `BuildingBlock` is keyed by build index, and conflating the two
    is the mistake `sites.inherit` is written to prevent one layer down.
    """
    cmap = canonical_map(reg, structure_id)
    return {cmap[s.atom_idx]: s for s in sites}
