"""Integrity checking: does every stored row still agree with its own graph?

Ground rule 1 says the registry is only ever written through `registry.api`.  That is a
rule, and rules are not enforcement — the first consumer of this schema wrote rows with
raw INSERTs and produced a database describing molecules that cannot exist (a
paddlewheel with no metal-metal bond, `muN` bridging on a two-metal node, eight dative
bonds recorded as three).  Every column was decorative because nothing ever checked a
column against the graph it claimed to summarise.

This module is that check.  It re-derives everything derivable from the stored typed
graph and reports disagreements.  It is cheap enough to run in CI and after any bulk
operation, and it is the difference between "the writer is correct" being a belief and
being a measurement.

Two severities, deliberately distinguished:

* **error** — the row contradicts its graph.  Something wrote around the API, or the
  API has a bug.  Either way the data is wrong.
* **stale** — the row is self-consistent but was produced by an older algorithm version.
  Expected after a deliberate version bump; means "recompute", not "corrupt".
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterator, Literal

from mofsbu.graph import TypedGraph, canonical_order, certificate_digest
from mofsbu.identity import l0_composition, wl_index
from mofsbu.registry.db import Registry
from mofsbu.versions import ALGO_VERSIONS

Severity = Literal["error", "stale"]


@dataclass(frozen=True)
class Problem:
    severity: Severity
    check: str
    detail: str
    structure_id: int | None = None
    geometry_id: int | None = None

    def __str__(self) -> str:
        where = ""
        if self.structure_id is not None:
            where = f" structure {self.structure_id}"
        if self.geometry_id is not None:
            where += f" geometry {self.geometry_id}"
        return f"[{self.severity}] {self.check}{where}: {self.detail}"


def _cmp(problems: list[Problem], sid: int, check: str, stored, recomputed) -> None:
    if stored != recomputed:
        problems.append(Problem("error", check, f"stored {stored!r} != graph {recomputed!r}",
                                structure_id=sid))


def verify_structure(reg: Registry, row) -> list[Problem]:
    """Re-derive every derivable column for one structure and compare."""
    from mofsbu.registry.api import _derive_columns

    sid = int(row["id"])
    problems: list[Problem] = []

    # 1. the graph must be there and must parse
    try:
        blob = reg.store.get(row["typed_graph_hash"])
    except KeyError:
        return [Problem("error", "graph_blob_missing",
                        f"typed_graph_hash {row['typed_graph_hash'][:16]}… not in the store",
                        structure_id=sid)]
    try:
        g = TypedGraph.from_json(blob)
    except Exception as exc:                                   # noqa: BLE001
        return [Problem("error", "graph_blob_unparseable",
                        f"{type(exc).__name__}: {exc}", structure_id=sid)]

    # 2. version drift is reported, never silently absorbed (ground rule 6)
    versions = {"algo_l0": "l0_composition", "algo_l1": "l1_certificate",
                "algo_l2": "l2_isomer_tag", "algo_canon": "canonical_order"}
    current_algo = True
    for column, name in versions.items():
        if row[column] != ALGO_VERSIONS[name]:
            current_algo = False
            problems.append(Problem(
                "stale", "algo_version",
                f"{name} row={row[column]!r} current={ALGO_VERSIONS[name]!r}",
                structure_id=sid))

    # 3. identity must reproduce from the graph — only meaningful at the current version
    if current_algo:
        order = canonical_order(g)
        _cmp(problems, sid, "l1_graph_hash", row["l1_graph_hash"], certificate_digest(g, order))
        _cmp(problems, sid, "wl_index", row["wl_index"], wl_index(g))
        try:
            _cmp(problems, sid, "l0_composition", row["l0_composition"], l0_composition(g))
        except Exception as exc:                               # noqa: BLE001
            problems.append(Problem("error", "l0_composition",
                                    f"graph cannot produce an L0: {exc}", structure_id=sid))
        stored_order = json.loads(row["canonical_order_json"])
        if sorted(stored_order) != g.nodes():
            problems.append(Problem("error", "canonical_order",
                                    "stored order is not a permutation of the graph's atoms",
                                    structure_id=sid))
        elif stored_order != order:
            problems.append(Problem("error", "canonical_order",
                                    "stored order disagrees with the canonical one",
                                    structure_id=sid))

    # 4. every denormalised column must be what the graph says
    derived = _derive_columns(g)
    for column in ("formula", "metals", "n_metals", "n_atoms", "net_charge", "multiplicity",
                   "ox_states", "max_bridge_class", "n_dative_bonds", "has_metal_metal"):
        _cmp(problems, sid, column, row[column], derived[column])

    # 5. the cached label must be what the graph derives, and the searchable fragment
    #    decomposition must match it.  A label is a projection; if the cache has drifted,
    #    something wrote a label by hand, which is exactly what must not happen.
    from mofsbu.naming import compose_label, decompose

    aliases = {r["fragment_l1"]: r["name"]
               for r in reg.conn.execute("SELECT fragment_l1, name FROM fragment_aliases")}
    expected_label = compose_label(g, aliases=aliases)
    if row["display_label"] != expected_label:
        problems.append(Problem("error", "display_label",
                                f"stored {row['display_label']!r} != derived "
                                f"{expected_label!r}", structure_id=sid))
    expected_fragments = {f.l1: (f.formula, f.count) for f in decompose(g)}
    stored_fragments = {r["fragment_l1"]: (r["formula"], r["count"])
                        for r in reg.conn.execute(
                            "SELECT fragment_l1, formula, count FROM structure_fragments "
                            "WHERE structure_id=?", (sid,))}
    if stored_fragments != expected_fragments:
        problems.append(Problem("error", "structure_fragments",
                                f"stored {stored_fragments} != graph {expected_fragments}",
                                structure_id=sid))

    # 6. the best-geometry cache must match the geometries it caches
    best = reg.conn.execute(
        "SELECT id, fidelity FROM geometries WHERE structure_id=? "
        "ORDER BY fidelity DESC, (converged IS 1) DESC, "
        "         CASE WHEN energy IS NULL THEN 1 ELSE 0 END, energy ASC, id ASC LIMIT 1",
        (sid,)).fetchone()
    expect_id = best["id"] if best else None
    expect_fid = best["fidelity"] if best else None
    if row["best_geometry_id"] != expect_id or row["best_fidelity"] != expect_fid:
        problems.append(Problem(
            "error", "best_geometry",
            f"row=({row['best_geometry_id']}, {row['best_fidelity']}) "
            f"geometries say ({expect_id}, {expect_fid})", structure_id=sid))

    # 7. geometries must describe this structure and their coordinates must exist
    for geom in reg.conn.execute(
            "SELECT id, coords_hash, n_atoms, method_id FROM geometries WHERE structure_id=?",
            (sid,)):
        gid = int(geom["id"])
        if geom["n_atoms"] != derived["n_atoms"]:
            problems.append(Problem("error", "geometry_atom_count",
                                    f"{geom['n_atoms']} atoms vs structure's "
                                    f"{derived['n_atoms']}", structure_id=sid, geometry_id=gid))
        if geom["method_id"] is None:
            problems.append(Problem("error", "geometry_method",
                                    "no method recorded (ground rule 3)",
                                    structure_id=sid, geometry_id=gid))
        if not reg.store.has(geom["coords_hash"]):
            problems.append(Problem("error", "coords_blob_missing",
                                    f"{geom['coords_hash'][:16]}… not in the store",
                                    structure_id=sid, geometry_id=gid))
    return problems


def iter_problems(reg: Registry) -> Iterator[Problem]:
    for row in reg.conn.execute("SELECT * FROM structures ORDER BY id"):
        yield from verify_structure(reg, row)


def verify(reg: Registry) -> list[Problem]:
    """Every disagreement between the stored rows and the graphs they summarise."""
    return list(iter_problems(reg))


def summarise(problems: list[Problem]) -> str:
    if not problems:
        return "registry is consistent: every column agrees with its graph"
    by_check: dict[tuple[str, str], int] = {}
    for p in problems:
        by_check[(p.severity, p.check)] = by_check.get((p.severity, p.check), 0) + 1
    lines = [f"{len(problems)} problem(s):"]
    for (sev, check), n in sorted(by_check.items()):
        lines.append(f"  {sev:6s} {check:24s} x{n}")
    return "\n".join(lines)
