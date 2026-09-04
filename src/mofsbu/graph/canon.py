"""Canonicalisation: the WL hash (identity key) and canonical atom order (site keys).

Three layers, cheapest first (design doc §4.1 / D3):

1. `wl_hash`      — networkx Weisfeiler-Lehman.  Fast, fixed length, pure Python, so
                    it produces the SAME key on the laptop and the workstation.  This
                    is the PRIMARY stored key (PLAN §0 rule 9).
2. `canonical_order` — individualisation-refinement canonical labelling.  A pure-Python
                    stand-in for nauty that runs on both machines.  It fixes the atom
                    ordering used as the coordinate system for site annotations (D5),
                    and its certificate is exact, so it also settles WL collisions.
3. `is_isomorphic` — exact check via canonical certificate, with `vf2_isomorphic`
                    available as an independent cross-check in tests.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict

import networkx as nx
from networkx.algorithms.isomorphism import categorical_edge_match, categorical_node_match

from mofsbu.graph._types import TypedGraph
from mofsbu._types import CanonicalisationError
from mofsbu.versions import ALGO_VERSIONS

DEFAULT_WL_ITERATIONS = 3
DEFAULT_MAX_LEAVES = 200_000

Cells = list[list[int]]
Certificate = tuple[tuple[str, ...], tuple[tuple[int, int, str], ...]]


# -- layer 1: WL hash ----------------------------------------------------------

def wl_hash(g: TypedGraph, iterations: int = DEFAULT_WL_ITERATIONS) -> str:
    """Fast bucket INDEX — deliberately not the identity key.

    1-WL cannot separate a mu2-bridging carboxylate from a chelating one (an
    8-membered M-O-C-O-M-O-C-O ring versus two 4-membered chelate rings): every
    atom has the same local environment in both, so the colours are stable from the
    first iteration and more iterations do not help.  That is a distinction this
    project exists to make, so WL is demoted to an index and `l1_graph_hash` uses the
    exact canonical certificate instead.  See `tests/test_canon.py` for the
    demonstration.
    """
    return nx.weisfeiler_lehman_graph_hash(
        g.to_networkx(), node_attr="lbl", edge_attr="elbl", iterations=iterations
    )


# -- layer 2: individualisation-refinement canonical labelling -----------------

def _initial_cells(g: TypedGraph) -> Cells:
    buckets: dict[str, list[int]] = defaultdict(list)
    for v in g.nodes():
        buckets[g.label(v).key()].append(v)
    return [sorted(buckets[k]) for k in sorted(buckets)]


def _refine(adj: dict[int, list[tuple[str, int]]], cells: Cells) -> Cells:
    """1-WL colour refinement of an ordered partition, to a stable coloring."""
    cells = [list(c) for c in cells]
    while True:
        colour = {v: ci for ci, cell in enumerate(cells) for v in cell}
        out: Cells = []
        changed = False
        for cell in cells:
            if len(cell) == 1:
                out.append(cell)
                continue
            sigs: dict[tuple, list[int]] = defaultdict(list)
            for v in cell:
                sigs[tuple(sorted((et, colour[u]) for et, u in adj[v]))].append(v)
            if len(sigs) == 1:
                out.append(cell)
                continue
            changed = True
            for sig in sorted(sigs):
                out.append(sorted(sigs[sig]))
        cells = out
        if not changed:
            return cells


def _adj_pairs(g: TypedGraph, v: int) -> list[tuple[str, int]]:
    return sorted((g.edge_type(v, u).value, u) for u in g.neighbors(v))


def certificate(g: TypedGraph, order: list[int]) -> Certificate:
    """Isomorphism-invariant fingerprint of `g` under one atom ordering.

    Built position by position so that fixing the first k atoms determines exactly
    the first k entries.  That prefix property is what makes branch-and-bound in
    `canonical_order` sound.
    """
    return tuple(_cert_entries(g, order))


def _cert_entries(g: TypedGraph, order: list[int]) -> list[tuple]:
    pos = {v: k for k, v in enumerate(order)}
    out = []
    for k, v in enumerate(order):
        back = tuple(sorted(
            (pos[u], g.edge_type(v, u).value)
            for u in g.neighbors(v) if u in pos and pos[u] < k
        ))
        out.append((g.label(v).key(), back))
    return out


class _UnionFind:
    __slots__ = ("parent",)

    def __init__(self, items) -> None:
        self.parent = {i: i for i in items}

    def find(self, a: int) -> int:
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def _known_prefix(cells: Cells) -> list[int]:
    """The leading run of atoms whose canonical positions are already fixed."""
    out: list[int] = []
    for c in cells:
        if len(c) != 1:
            break
        out.append(c[0])
    return out


def canonical_order(g: TypedGraph, *, max_leaves: int = DEFAULT_MAX_LEAVES) -> list[int]:
    return _canonicalise(g, max_leaves=max_leaves)[0]


def automorphism_generators(g: TypedGraph, *,
                            max_leaves: int = DEFAULT_MAX_LEAVES) -> list[dict[int, int]]:
    """A generating set for Aut(g), discovered during canonical labelling.

    The search already finds these — every leaf that ties the best certificate IS an
    automorphism — and already uses them to prune.  Promoting them out of the search is
    what lets symmetry be *asked about* rather than only exploited internally: it is how
    "which of these substitution patterns are the same" gets answered without anyone
    writing down a taxonomy.
    """
    return _canonicalise(g, max_leaves=max_leaves)[1]


def automorphism_group(g: TypedGraph, *, max_size: int = 100_000) -> list[dict[int, int]]:
    """Aut(g) in full, by closure over the generators.

    Bounded: a graph whose symmetry group exceeds `max_size` raises rather than hanging.
    """
    nodes = g.nodes()
    identity = tuple(range(len(nodes)))
    index = {v: k for k, v in enumerate(nodes)}
    gens = []
    for perm in automorphism_generators(g):
        gens.append(tuple(index[perm.get(v, v)] for v in nodes))

    seen = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for gen in gens:
            composed = tuple(gen[i] for i in current)
            if composed not in seen:
                if len(seen) >= max_size:
                    raise CanonicalisationError(
                        f"automorphism group of {g!r} exceeds {max_size} elements")
                seen.add(composed)
                frontier.append(composed)
    return [{nodes[i]: nodes[p[i]] for i in range(len(nodes))} for p in sorted(seen)]


def atom_orbits(g: TypedGraph) -> list[tuple[int, ...]]:
    """Atoms grouped into symmetry-equivalent classes under Aut(g)."""
    parent = {v: v for v in g.nodes()}

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for perm in automorphism_generators(g):
        for v, w in perm.items():
            ra, rb = find(v), find(w)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    groups: dict[int, list[int]] = {}
    for v in g.nodes():
        groups.setdefault(find(v), []).append(v)
    return sorted(tuple(sorted(vs)) for vs in groups.values())


def _canonicalise(g: TypedGraph, *,
                  max_leaves: int = DEFAULT_MAX_LEAVES) -> tuple[list[int], list[dict[int, int]]]:
    """Canonical atom ordering: `order[k]` is the atom at canonical position k.

    Individualisation-refinement, minimising the certificate over the search tree.
    Two prunings keep it tractable on the symmetric species that MOF nodes are made
    of (six equivalent aqua ligands, four equivalent paddlewheel bridges):

    * **branch and bound** on the certificate prefix — a partial ordering already
      worse than the best complete one cannot be completed into a better one;
    * **automorphism pruning** — when a leaf ties the best certificate, the
      permutation between them is an automorphism of the graph.  At any node, the
      automorphisms that fix every atom individualised so far generate orbits, and
      only one atom per orbit needs to be tried.  Restricting to path-fixing
      automorphisms is what keeps this sound (they are stabiliser elements, so the
      subtrees really are equivalent).

    The result defines the canonical atom indices every site record keys against
    (D5), so its recipe is pinned in `ALGO_VERSIONS['canonical_order']`: changing it
    changes every stored key.
    """
    if len(g) == 0:
        return [], []
    adj = {v: _adj_pairs(g, v) for v in g.nodes()}
    best: list[tuple] | None = None
    best_order: list[int] | None = None
    autos: list[dict[int, int]] = []
    leaves = 0

    def search(cells: Cells, path: list[int]) -> None:
        nonlocal best, best_order, leaves
        cells = _refine(adj, cells)
        prefix = _known_prefix(cells)
        entries = _cert_entries(g, prefix)
        if best is not None:
            head = best[:len(entries)]
            if entries > head:
                return                      # bound: cannot beat the incumbent
        target = next((i for i, c in enumerate(cells) if len(c) > 1), None)
        if target is None:
            leaves += 1
            if leaves > max_leaves:
                raise CanonicalisationError(
                    f"canonical labelling exceeded {max_leaves} leaves on {g!r}"
                )
            order = [c[0] for c in cells]
            cert = list(_cert_entries(g, order))
            if best is None or cert < best:
                best, best_order = cert, order
            elif cert == best and best_order is not None:
                autos.append({best_order[i]: order[i] for i in range(len(order))})
            return

        cell = cells[target]
        usable = [p for p in autos if all(p.get(v, v) == v for v in path)]
        uf = _UnionFind(cell)
        for p in usable:
            for v in cell:
                w = p.get(v, v)
                if w in uf.parent:
                    uf.union(v, w)
        seen: set[int] = set()
        for v in sorted(cell):
            root = uf.find(v)
            if root in seen:
                continue                    # same orbit as a sibling already tried
            seen.add(root)
            rest = [u for u in cell if u != v]
            search(cells[:target] + [[v], rest] + cells[target + 1:], path + [v])

    search(_initial_cells(g), [])
    assert best_order is not None
    return best_order, autos


def canonical_certificate(g: TypedGraph) -> Certificate:
    """The exact canonical form.  Equal certificates <=> isomorphic graphs."""
    return certificate(g, canonical_order(g))


def certificate_digest(g: TypedGraph, order: list[int] | None = None) -> str:
    """sha256 over the canonical certificate — the L1 key producer.

    Exact (no collisions) and fixed length.  Pure Python, so the laptop and the
    workstation agree (ground rule 9); pynauty stays a verifier and must never
    produce a stored key, because its canonical form differs from this one.

    Pass `order` when the caller has already canonicalised, to avoid paying for it
    twice — the registry needs both the digest and the order on every insert.
    """
    cert = certificate(g, order) if order is not None else canonical_certificate(g)
    payload = json.dumps(cert, sort_keys=True,
                         separators=(",", ":"), default=list).encode()
    return hashlib.sha256(payload).hexdigest()


def canonical_index_map(g: TypedGraph) -> dict[int, int]:
    """{build-order atom index -> canonical index}.  What site records key against."""
    return {v: k for k, v in enumerate(canonical_order(g))}


# -- layer 3: isomorphism ------------------------------------------------------

def is_isomorphic(a: TypedGraph, b: TypedGraph) -> bool:
    """Exact test via canonical certificate — the WL collision resolver."""
    return canonical_certificate(a) == canonical_certificate(b)


def vf2_isomorphic(a: TypedGraph, b: TypedGraph) -> bool:
    """Independent VF2 check.  Used to cross-validate `is_isomorphic` in tests."""
    return nx.is_isomorphic(
        a.to_networkx(), b.to_networkx(),
        node_match=categorical_node_match("lbl", None),
        edge_match=categorical_edge_match("elbl", None),
    )


def algo_versions() -> dict[str, str]:
    return {k: ALGO_VERSIONS[k]
            for k in ("l1_certificate", "wl_index", "canonical_order", "graph_schema")}
