"""The only write surface (ground rule 1) and the query surface over it.

`put_structure` is the point where a typed graph and its canonical map become a
persisted record.  It is idempotent on the identity key: inserting the same structure
twice yields ONE row and TWO provenance edges, which is the D2 claim made real —
identity on the node, sequence on the edges.
"""
from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, NamedTuple

from mofsbu.graph import BridgeClass, EdgeType, TypedGraph, canonical_order, certificate_digest
from mofsbu.identity import block_id, hill_formula, l0_composition, l2_isomer_tag, wl_index
from mofsbu.naming import compose_label, decompose
from mofsbu.registry.db import Registry, utcnow
# MethodSpec lives in `_types` beside `Fidelity` so that `energy` (which produces
# numbers) and this module (which stores them) share one definition without either
# importing the other.  Re-exported here because callers have always imported it
# from `mofsbu.registry`.
from mofsbu._types import Fidelity, MethodSpec, MofsbuError
from mofsbu.versions import ALGO_VERSIONS

_BRIDGE_RANK = {BridgeClass.NONE: 0, BridgeClass.TERMINAL: 1, BridgeClass.MU2: 2,
                BridgeClass.MU3: 3, BridgeClass.MU_N: 4}


class Put(NamedTuple):
    """`created` distinguishes a fresh row from a recognised duplicate."""

    id: int
    created: bool


class RegistryError(MofsbuError):
    pass


@dataclass(frozen=True)
class Provenance:
    """A DAG edge: how this structure came to be.  Never part of its identity."""

    kind: str = "assembly"          # ingest | assembly | reaction
    reagent_ids: tuple[int, ...] = ()
    note: str = ""
    atom_map: dict[int, int] | None = None
    choice_vector_digest: str | None = None
    depth: int | None = None
    intermediate: bool = False


# ── derivation ───────────────────────────────────────────────────────────────

def _set_column(values: Iterable[str]) -> str:
    """Bracketed-comma form, so LIKE '%,Cu,%' cannot match a prefix of another value."""
    vals = sorted(values)
    return "," + ",".join(vals) + "," if vals else ""


def formula(g: TypedGraph) -> str:
    counts = g.element_counts()
    metal_counts: dict[str, int] = {}
    for i in g.metals():
        el = g.label(i).element
        metal_counts[el] = metal_counts.get(el, 0) + 1
    organic = {k: v for k, v in counts.items() if k not in metal_counts}
    head = "".join(f"{el}{metal_counts[el]}" if metal_counts[el] > 1 else el
                   for el in sorted(metal_counts))
    return head + hill_formula(organic)


def display_label(g: TypedGraph, aliases: dict[str, str] | None = None) -> str:
    """Composed from the GRAPH, never from an author-supplied name.

    `TypedGraph.name` is deliberately not consulted.  It was, and the result looked
    convincing exactly where it was least trustworthy: a hand-typed name made
    author-written examples read perfectly while anything the enumerator generated
    degraded to a bare formula that collides across isomers.  A label has to be a
    function of the structure or it is not a label, it is a memory of what someone
    called it once.
    """
    return compose_label(g, aliases=aliases)


def _derive_columns(g: TypedGraph) -> dict[str, Any]:
    metal_nodes = g.metals()
    edges = g.edges()
    return {
        "formula": formula(g),
        "metals": _set_column(g.label(i).element for i in metal_nodes),
        "n_metals": len(metal_nodes),
        "n_atoms": len(g),
        "net_charge": g.net_charge(),
        "multiplicity": g.multiplicity,
        "ox_states": _set_column(
            f"{g.label(i).element}{g.label(i).oxidation_state:+d}"
            for i in metal_nodes if g.label(i).oxidation_state is not None
        ),
        "max_bridge_class": g.max_bridge_class().value,
        "n_dative_bonds": sum(1 for _u, _v, et in edges if et is EdgeType.DATIVE),
        "has_metal_metal": int(any(et is EdgeType.METAL_METAL for _u, _v, et in edges)),
        "display_label": display_label(g),
    }


# ── methods ──────────────────────────────────────────────────────────────────

def find_method_id(reg: Registry, spec: MethodSpec) -> int | None:
    """The id of an existing method row, or None.  Never inserts.

    `method_id` registers on miss, which is right when you are about to store a number
    and wrong when you are asking "has this already been computed?" — that question is
    asked before deciding whether to spend an hour of xTB, and it must not have the side
    effect of creating the row it was looking for.
    """
    row = reg.conn.execute(
        "SELECT id FROM methods WHERE code=? AND code_version=? AND method=? "
        "AND solvent IS ? AND charge IS ? AND multiplicity IS ? AND extras_json=?",
        (spec.code, spec.code_version, spec.method, spec.solvent, spec.charge,
         spec.multiplicity, json.dumps(spec.extras, sort_keys=True, separators=(",", ":")))
    ).fetchone()
    return None if row is None else int(row["id"])


def method_id(reg: Registry, spec: MethodSpec) -> int:
    extras = json.dumps(spec.extras, sort_keys=True, separators=(",", ":"))
    key = (spec.code, spec.code_version, spec.method, spec.solvent,
           spec.charge, spec.multiplicity, extras)
    row = reg.conn.execute(
        "SELECT id FROM methods WHERE code=? AND code_version=? AND method=? "
        "AND solvent IS ? AND charge IS ? AND multiplicity IS ? AND extras_json=?", key
    ).fetchone()
    if row:
        return row["id"]
    try:
        cur = reg.conn.execute(
            "INSERT INTO methods (code, code_version, method, solvent, charge, multiplicity, "
            "extras_json) VALUES (?,?,?,?,?,?,?)", key)
    except sqlite3.IntegrityError:                      # another worker registered it first
        row = reg.conn.execute(
            "SELECT id FROM methods WHERE code=? AND code_version=? AND method=? "
            "AND solvent IS ? AND charge IS ? AND multiplicity IS ? AND extras_json=?",
            key).fetchone()
        if row is None:
            raise
        return int(row["id"])
    return int(cur.lastrowid)


# ── structures ───────────────────────────────────────────────────────────────

def put_structure(
    reg: Registry,
    g: TypedGraph,
    *,
    l2: str | None = None,
    provenance: Provenance | None = None,
    tags: Iterable[str] = (),
) -> Put:
    """Persist a typed graph as a structure record.  Idempotent on (L0, L1, L2).

    Stores the canonical order alongside the hash rather than recomputing it later:
    it is the coordinate system every site annotation keys against, and it is only
    deterministic relative to a pinned algorithm version (ground rule 6).
    """
    order = canonical_order(g)                       # computed ONCE, reused below
    l0 = l0_composition(g)                           # raises if charge/spin unset
    l1 = certificate_digest(g, order)
    tag2 = l2_isomer_tag(g) if l2 is None else l2

    row = reg.conn.execute(
        "SELECT id FROM structures WHERE l0_composition=? AND l1_graph_hash=? "
        "AND l2_isomer_tag=?", (l0, l1, tag2)
    ).fetchone()

    if row is not None:
        sid, created = int(row["id"]), False
    else:
        cols = _derive_columns(g)
        try:
            cur = reg.conn.execute(
            "INSERT INTO structures (l0_composition, l1_graph_hash, l2_isomer_tag, block_id, "
            " wl_index, typed_graph_hash, canonical_order_json, algo_l0, algo_l1, algo_l2, "
            " algo_canon, formula, metals, n_metals, n_atoms, net_charge, multiplicity, "
            " ox_states, max_bridge_class, n_dative_bonds, has_metal_metal, display_label, "
            " created_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (l0, l1, tag2, block_id(l0, l1, tag2), wl_index(g),
             reg.store.put(g.to_json()), json.dumps(order),
             ALGO_VERSIONS["l0_composition"], ALGO_VERSIONS["l1_certificate"],
             ALGO_VERSIONS["l2_isomer_tag"], ALGO_VERSIONS["canonical_order"],
             cols["formula"], cols["metals"], cols["n_metals"], cols["n_atoms"],
             cols["net_charge"], cols["multiplicity"], cols["ox_states"],
             cols["max_bridge_class"], cols["n_dative_bonds"], cols["has_metal_metal"],
                 cols["display_label"], utcnow()),
            )
        except sqlite3.IntegrityError:
            # Lost a race: another worker inserted this identity between our SELECT and
            # our INSERT.  The UNIQUE constraint is what makes that safe rather than
            # silently duplicating, and re-reading is the correct resolution — the other
            # worker's row IS this structure.  Idempotency on its own is not concurrency
            # safety; check-then-act needs the constraint to arbitrate.
            existing = reg.conn.execute(
                "SELECT id FROM structures WHERE l0_composition=? AND l1_graph_hash=? "
                "AND l2_isomer_tag=?", (l0, l1, tag2)).fetchone()
            if existing is None:
                raise
            sid, created = int(existing["id"]), False
        else:
            sid, created = int(cur.lastrowid), True
            _put_structure_metals(reg, sid, g)
            _put_structure_fragments(reg, sid, g)

    for tag in tags:
        reg.conn.execute(
            "INSERT OR IGNORE INTO structure_tags (structure_id, tag) VALUES (?,?)", (sid, tag)
        )
    if provenance is not None:
        put_reaction(reg, sid, provenance)
    return Put(sid, created)


def _put_structure_metals(reg: Registry, sid: int, g: TypedGraph) -> None:
    tally: dict[tuple[str, int | None, str | None], int] = {}
    for i in g.metals():
        lab = g.label(i)
        key = (lab.element, lab.oxidation_state, lab.spin_class)
        tally[key] = tally.get(key, 0) + 1
    for (sym, ox, spin), n in tally.items():
        reg.conn.execute(
            "INSERT OR REPLACE INTO structure_metals (structure_id, symbol, oxidation_state, "
            "spin_class, count) VALUES (?,?,?,?,?)", (sid, sym, ox, spin, n)
        )


def _put_structure_fragments(reg: Registry, sid: int, g: TypedGraph) -> None:
    """The searchable decomposition.  Queries go through here; labels come later."""
    for frag in decompose(g):
        reg.conn.execute(
            "INSERT OR REPLACE INTO structure_fragments (structure_id, fragment_l1, formula, "
            "count, bridge_class, n_metals_bound) VALUES (?,?,?,?,?,?)",
            (sid, frag.l1, frag.formula, frag.count, frag.bridge_class, frag.n_metals_bound))


def alias_fragment(reg: Registry, fragment_l1: str, name: str, *, source: str = "") -> None:
    """Give a fragment a human name.  Display only — no query depends on it."""
    reg.conn.execute(
        "INSERT INTO fragment_aliases (fragment_l1, name, source, created_at) VALUES (?,?,?,?) "
        "ON CONFLICT(fragment_l1) DO UPDATE SET name=excluded.name, source=excluded.source",
        (fragment_l1, name, source, utcnow()))


def fragment_aliases(reg: Registry) -> dict[str, str]:
    return {r["fragment_l1"]: r["name"]
            for r in reg.conn.execute("SELECT fragment_l1, name FROM fragment_aliases")}


def relabel_all(reg: Registry) -> int:
    """Recompute every cached label from its graph and the current aliases.

    The cache exists so a listing need not re-parse thousands of graphs; this is how it
    is refreshed after aliases change.  Labels are still never searched.
    """
    aliases = fragment_aliases(reg)
    n = 0
    for row in reg.conn.execute("SELECT id, typed_graph_hash FROM structures"):
        g = TypedGraph.from_json(reg.store.get(row["typed_graph_hash"]))
        reg.conn.execute("UPDATE structures SET display_label=? WHERE id=?",
                         (compose_label(g, aliases=aliases), row["id"]))
        n += 1
    return n


def put_reaction(reg: Registry, product_id: int, prov: Provenance) -> int:
    cur = reg.conn.execute(
        "INSERT INTO reactions (product_structure_id, kind, intermediate, atom_map_json, "
        " choice_vector_digest, depth, note, created_at) VALUES (?,?,?,?,?,?,?,?)",
        (product_id, prov.kind, int(prov.intermediate),
         json.dumps(prov.atom_map) if prov.atom_map else None,
         prov.choice_vector_digest, prov.depth, prov.note, utcnow()),
    )
    rid = int(cur.lastrowid)
    for reagent in prov.reagent_ids:
        reg.conn.execute(
            "INSERT OR IGNORE INTO reaction_reagents (reaction_id, structure_id) VALUES (?,?)",
            (rid, reagent),
        )
    return rid


# ── geometries ───────────────────────────────────────────────────────────────

def put_geometry(
    reg: Registry,
    structure_id: int,
    xyz_text: str,
    *,
    fidelity: Fidelity,
    method: MethodSpec,
    energy: float | None = None,
    converged: bool | None = None,
    relaxed_from: int | None = None,
    choice_vector: dict | None = None,
    seed: int | None = None,
    qc: dict | None = None,
) -> Put:
    """Store one geometry realisation.  `method` is required, not optional.

    Fidelity is a property of the geometry, never of the structure (D4), so one
    identity can carry a raw construct, an xTB relaxation and a DFT relaxation at once.
    """
    struct = reg.conn.execute(
        "SELECT n_atoms FROM structures WHERE id=?", (structure_id,)
    ).fetchone()
    if struct is None:
        raise RegistryError(f"no structure {structure_id}")

    n_atoms = int(xyz_text.strip().splitlines()[0])
    if n_atoms != struct["n_atoms"]:
        raise RegistryError(
            f"geometry has {n_atoms} atoms, structure {structure_id} has {struct['n_atoms']}"
        )

    mid = method_id(reg, method)
    coords_hash = reg.store.put_text(xyz_text)
    row = reg.conn.execute(
        "SELECT id FROM geometries WHERE structure_id=? AND coords_hash=? AND method_id=?",
        (structure_id, coords_hash, mid),
    ).fetchone()
    if row is not None:
        return Put(int(row["id"]), False)

    try:
        cur = reg.conn.execute(
            "INSERT INTO geometries (structure_id, coords_hash, n_atoms, fidelity, method_id, "
            " energy, converged, relaxed_from, choice_vector_digest, choice_vector_json, seed, "
            " qc_json, created_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (structure_id, coords_hash, n_atoms, int(fidelity), mid, energy,
         None if converged is None else int(converged), relaxed_from,
         (choice_vector or {}).get("digest"),
         json.dumps(choice_vector) if choice_vector else None,
             seed, json.dumps(qc) if qc else None, utcnow()),
        )
    except sqlite3.IntegrityError:                      # concurrent identical geometry
        existing = reg.conn.execute(
            "SELECT id FROM geometries WHERE structure_id=? AND coords_hash=? AND method_id=?",
            (structure_id, coords_hash, mid)).fetchone()
        if existing is None:
            raise
        return Put(int(existing["id"]), False)
    gid = int(cur.lastrowid)
    _refresh_best_geometry(reg, structure_id)
    return Put(gid, True)


def _refresh_best_geometry(reg: Registry, structure_id: int) -> None:
    """Maintain the best-available pointer.  The ONLY writer of these two columns.

    Highest fidelity wins; among equals, the lowest energy; a converged geometry beats
    an unconverged one at the same fidelity.
    """
    row = reg.conn.execute(
        "SELECT id, fidelity FROM geometries WHERE structure_id=? "
        "ORDER BY fidelity DESC, (converged IS 1) DESC, "
        "         CASE WHEN energy IS NULL THEN 1 ELSE 0 END, energy ASC, id ASC LIMIT 1",
        (structure_id,),
    ).fetchone()
    if row is None:
        reg.conn.execute(
            "UPDATE structures SET best_geometry_id=NULL, best_fidelity=NULL WHERE id=?",
            (structure_id,))
    else:
        reg.conn.execute(
            "UPDATE structures SET best_geometry_id=?, best_fidelity=? WHERE id=?",
            (row["id"], row["fidelity"], structure_id))


# ── sites (M4) ───────────────────────────────────────────────────────────────

def put_sites(reg: Registry, structure_id: int, sites: list, *, algo: str = "perception/1") -> int:
    """Write the site catalog for a structure.  Perceive once (D5).

    Sites are keyed by CANONICAL atom index, taken from the order frozen at insert, so
    they survive recall and re-ordering.  Re-running perception on the same structure
    replaces the catalog rather than duplicating it.
    """
    reg.conn.execute("DELETE FROM site_catalog WHERE structure_id=?", (structure_id,))
    cmap = canonical_map(reg, structure_id)
    n = 0
    for site in sites:
        canonical_idx = cmap[site.atom_idx]
        reg.conn.execute(
            "INSERT INTO site_catalog (structure_id, canonical_idx, donor_type, labile, "
            " charge_after, live_dof, binding_modes, frame_json, algo_perception) "
            "VALUES (?,?,?,?,?,?,?,?,?)",
            (structure_id, canonical_idx, site.donor_type, int(site.labile),
             site.charge_after, site.live_dof, _set_column(site.binding_modes),
             json.dumps(site.frame) if site.frame else None, algo))
        n += 1
    # n_perceived_donors is what was actually counted.  n_open_sites needs site_state,
    # which is not populated yet, so it stays NULL rather than being filled with a
    # different quantity that happens to be an integer.
    reg.conn.execute(
        "UPDATE structures SET donor_types=?, n_perceived_donors=?, n_open_sites=NULL "
        "WHERE id=?", (_set_column({s.donor_type for s in sites}), n, structure_id))
    return n


def get_sites(reg: Registry, structure_id: int) -> list:
    return list(reg.conn.execute(
        "SELECT * FROM site_catalog WHERE structure_id=? ORDER BY canonical_idx",
        (structure_id,)))


# ── reads ────────────────────────────────────────────────────────────────────

def get_structure(reg: Registry, ref: int | str) -> sqlite3.Row:
    col = "id" if isinstance(ref, int) else "block_id"
    row = reg.conn.execute(f"SELECT * FROM v_structures WHERE {col}=?", (ref,)).fetchone()
    if row is None:
        raise RegistryError(f"no structure with {col}={ref!r}")
    return row


def get_graph(reg: Registry, structure_id: int) -> TypedGraph:
    """Round-trip the stored graph back into an object."""
    row = reg.conn.execute(
        "SELECT typed_graph_hash FROM structures WHERE id=?", (structure_id,)).fetchone()
    if row is None:
        raise RegistryError(f"no structure {structure_id}")
    return TypedGraph.from_json(reg.store.get(row["typed_graph_hash"]))


def canonical_map(reg: Registry, structure_id: int) -> dict[int, int]:
    """{atom index -> canonical index}, as frozen at insert.  Never recomputed."""
    row = reg.conn.execute(
        "SELECT canonical_order_json FROM structures WHERE id=?", (structure_id,)).fetchone()
    if row is None:
        raise RegistryError(f"no structure {structure_id}")
    return {v: k for k, v in enumerate(json.loads(row["canonical_order_json"]))}


def best_geometry(reg: Registry, structure_id: int,
                  min_fidelity: Fidelity | None = None) -> sqlite3.Row | None:
    if min_fidelity is None:
        row = reg.conn.execute(
            "SELECT g.* FROM geometries g JOIN structures s ON s.best_geometry_id = g.id "
            "WHERE s.id=?", (structure_id,)).fetchone()
        return row
    return reg.conn.execute(
        "SELECT * FROM geometries WHERE structure_id=? AND fidelity>=? "
        "ORDER BY fidelity DESC, energy ASC LIMIT 1", (structure_id, int(min_fidelity))
    ).fetchone()


def geometry_xyz(reg: Registry, geometry_id: int) -> str:
    row = reg.conn.execute(
        "SELECT coords_hash FROM geometries WHERE id=?", (geometry_id,)).fetchone()
    if row is None:
        raise RegistryError(f"no geometry {geometry_id}")
    return reg.store.get_text(row["coords_hash"])


def export_xyz(reg: Registry, geometry_id: int, path: str | Path) -> Path:
    """Files are a VIEW of the registry, never a record in it."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(geometry_xyz(reg, geometry_id), encoding="utf-8")
    return p


# display_label is deliberately absent: it is a projection, not a query surface.
_SORTABLE = {"id", "formula", "net_charge", "n_metals", "n_atoms", "multiplicity",
             "best_fidelity", "best_energy", "created_at"}


def find(
    reg: Registry,
    *,
    l0: str | None = None,
    l1: str | None = None,
    formula_like: str | None = None,
    metal: str | None = None,
    n_metals: int | None = None,
    charge: int | None = None,
    max_bridge_class: str | None = None,
    has_metal_metal: bool | None = None,
    min_fidelity: Fidelity | None = None,
    has_route: bool | None = None,
    tag: str | None = None,
    fragment_l1: str | None = None,
    fragment_formula: str | None = None,
    sort: str = "id",
    order: str = "asc",
    limit: int = 100,
    offset: int = 0,
) -> list[sqlite3.Row]:
    """Query across the four filter axes.  Parameterised throughout."""
    if sort not in _SORTABLE:
        raise RegistryError(f"cannot sort on {sort!r}; allowed: {sorted(_SORTABLE)}")
    if order.lower() not in ("asc", "desc"):
        raise RegistryError(f"bad order {order!r}")

    where: list[str] = []
    args: list[Any] = []
    for clause, value in (
        ("l0_composition = ?", l0),
        ("l1_graph_hash = ?", l1),
        ("formula LIKE ?", f"%{formula_like}%" if formula_like else None),
        ("metals LIKE ?", f"%,{metal}," if metal else None),
        ("n_metals = ?", n_metals),
        ("net_charge = ?", charge),
        ("max_bridge_class = ?", max_bridge_class),
        ("best_fidelity >= ?", None if min_fidelity is None else int(min_fidelity)),
    ):
        if value is not None:
            where.append(clause)
            args.append(value)
    if has_metal_metal is not None:
        where.append("has_metal_metal = ?")
        args.append(int(has_metal_metal))
    if has_route is not None:
        where.append("n_incoming_routes > 0" if has_route else "n_incoming_routes = 0")
    if tag is not None:
        where.append("id IN (SELECT structure_id FROM structure_tags WHERE tag = ?)")
        args.append(tag)
    # Structural search: "which structures contain this ligand" is an indexed join on the
    # fragment's own identity, never a string match on a name.
    if fragment_l1 is not None:
        where.append("id IN (SELECT structure_id FROM structure_fragments "
                     "WHERE fragment_l1 = ?)")
        args.append(fragment_l1)
    if fragment_formula is not None:
        where.append("id IN (SELECT structure_id FROM structure_fragments WHERE formula = ?)")
        args.append(fragment_formula)

    sql = "SELECT * FROM v_structures"
    if where:
        sql += " WHERE " + " AND ".join(where)
    sql += f" ORDER BY {sort} {order.upper()} LIMIT ? OFFSET ?"
    return list(reg.conn.execute(sql, (*args, limit, offset)))
