"""Record which `mofsbu` function calls which, while a real operation runs.

A DYNAMIC call graph, not a static one, and the difference is the point: a static reader
lists every call that appears in the source, including the branches this path never takes
and the dispatch it cannot resolve.  This records what actually ran, so the graph is the
process flow for one operation rather than a map of the package.

    conda run -n ebu python scripts/callgraph.py bridge --format mermaid

Scenarios are registered in `SCENARIOS`; add one by writing a function that performs the
operation and returns a short label.  Output is Mermaid (renders anywhere Markdown does)
or JSON (edges and call counts, for anything else).

Nothing here is imported by `src/` — it is a reading tool, and it is in `scripts/` for the
same reason `regen_golden.py` is.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


# ── the recorder ─────────────────────────────────────────────────────────────

class CallRecorder:
    """Caller -> callee edges for functions defined under `src/mofsbu`.

    `sys.setprofile` fires on every call including C functions; we keep Python frames
    whose code object lives under the package, and read the caller off `frame.f_back`.
    A frame whose caller is outside the package is an ENTRY POINT, which is worth showing
    — it is where the operation crosses into the code being read.
    """

    def __init__(self, root: Path, package: str = "mofsbu") -> None:
        self.root = str(root)
        self.package = package
        self.calls: Counter[str] = Counter()
        self.depth: dict[str, int] = {}
        #: The full in-package stack at each call, not just the immediate caller.  Edges
        #: are derived from these at render time, so hiding a function RE-PARENTS its
        #: callees onto the nearest visible ancestor instead of orphaning them — which
        #: matters because the interesting callers (`_place_bridging_block`) are exactly
        #: the private ones the default filter removes.
        self.paths: Counter[tuple[str, ...]] = Counter()
        self._stack: list[str] = []

    def _name(self, frame) -> str | None:
        code = frame.f_code
        if not code.co_filename.startswith(self.root):
            return None
        module = frame.f_globals.get("__name__", "")
        if not module.startswith(self.package):
            return None
        qual = getattr(code, "co_qualname", code.co_name)
        return f"{module}.{qual}"

    def profile(self, frame, event, _arg):
        if event == "call":
            name = self._name(frame)
            if name is None:
                return
            self.calls[name] += 1
            self.depth[name] = min(self.depth.get(name, 99), len(self._stack))
            self._stack.append(name)
            self.paths[tuple(self._stack)] += 1
        elif event == "return" and self._stack:
            name = self._name(frame)
            if name is not None and self._stack[-1] == name:
                self._stack.pop()

    def __enter__(self):
        sys.setprofile(self.profile)
        return self

    def __exit__(self, *_exc):
        sys.setprofile(None)
        self._stack.clear()
        return False


# ── scenarios ────────────────────────────────────────────────────────────────

def scenario_bridge() -> str:
    """Two ordinary joins into a bridged dimer, then a mu2 bridge across it."""
    from mofsbu.assembly.construct import ligand_block, metal_block
    from mofsbu.assembly.join import bridge_compatible, join, join_bridge
    from mofsbu.geometry.embed import embed_molecule
    from mofsbu.graph.from_mol import mol_from_smiles

    def formate(name):
        return ligand_block(embed_molecule(mol_from_smiles("[O-]C=O"), seed=7),
                            charge=-1, name=name)

    def carbox(block):
        return sorted((s for s in block.open_donors()
                       if s.donor_type == "carboxylate_O"), key=lambda s: s.atom_idx)

    lig = formate("bridge1")
    donors = carbox(lig)
    m1 = metal_block("Cu", 2, cn=4, geometry="square_planar")
    first = join(lig, m1, donors[0], m1.open_vacancies()[0], lone_pair=1)
    free = next(s for s in first.block.open_donors()
                if s.atom_idx == first.atom_map[donors[1].atom_idx])
    m2 = metal_block("Cu", 2, cn=4, geometry="square_planar")
    dimer = join(first.block, m2, free, m2.open_vacancies()[0], lone_pair=1).block

    lig2 = formate("bridge2")
    d2 = carbox(lig2)
    m_a, m_b = dimer.graph.metals()
    per: dict[int, list] = {}
    for s in dimer.open_vacancies():
        per.setdefault(s.atom_idx, []).append(s)
    pairs = [(va, vb) for va in per[m_a] for vb in per[m_b]]
    best = min(pairs, key=lambda vs: bridge_compatible(
        d2, list(vs), partner="Cu", donor_elements=["O", "O"]).strain)
    join_bridge(lig2, dimer, d2, list(best), lone_pairs=(0, 0))
    return "join x2 -> bridged dimer, then join_bridge"


def scenario_identity() -> str:
    """A fixture graph through the identity layers."""
    from mofsbu import examples
    from mofsbu.identity.keys import l0_composition, l1_graph_hash

    g = examples.cu_paddlewheel()
    l0_composition(g)
    l1_graph_hash(g)
    return "cu_paddlewheel -> L0 + L1"


def scenario_perceive() -> str:
    """Donor perception and the frame model for one molecule."""
    from mofsbu.geometry.embed import embed_molecule
    from mofsbu.graph.from_mol import mol_from_smiles
    from mofsbu.sites.model import perceive

    perceive(embed_molecule(mol_from_smiles("[O-]C(=O)c1ccccc1"), seed=7))
    return "benzoate -> perceive"


SCENARIOS = {
    "bridge": scenario_bridge,
    "identity": scenario_identity,
    "perceive": scenario_perceive,
}


# ── rendering ────────────────────────────────────────────────────────────────

def _short(name: str) -> str:
    """`mofsbu.assembly.join.join_bridge` -> `join.join_bridge`."""
    if name == "<entry>":
        return name
    parts = name.split(".")
    return ".".join(parts[-2:]) if len(parts) > 2 else name


def _module(name: str) -> str:
    if name == "<entry>":
        return "<entry>"
    parts = name.split(".")
    return ".".join(parts[1:-1]) if len(parts) > 2 else name


def _visible(rec: CallRecorder, *, min_calls: int, private: bool,
             max_depth: int | None, modules: tuple[str, ...] = ()) -> set[str]:
    """Which functions survive onto the diagram.

    Three filters, and the defaults are chosen so the result reads as a process flow
    rather than as a census.  Comprehensions and lambdas go unconditionally: `<locals>`
    frames are the single largest group by count and none of them is a step in the
    story.  `_private` helpers go by default for the same reason at one remove.  And
    `max_depth` cuts by position in the call tree, which is what "show me the top level"
    actually means — a leaf called from everywhere is not top level however often it runs.
    """
    keep = {n for n, c in rec.calls.items()
            if c >= min_calls and "<locals>" not in n and not n.endswith(".<module>")}
    if not private:
        keep = {n for n in keep if not n.split(".")[-1].startswith("_")}
    if max_depth is not None:
        keep = {n for n in keep if rec.depth.get(n, 99) <= max_depth}
    if modules:
        # Narrowing to a subsystem is the common reading question ("what does assembly
        # do here"), and it works because `collapsed_edges` re-parents through whatever
        # is hidden — so the calls that leave and re-enter the subsystem still connect.
        keep = {n for n in keep
                if any(n.startswith(f"mofsbu.{m}.") or n.startswith(f"{m}.")
                       for m in modules)}
    return keep


def to_mermaid(rec: CallRecorder, *, min_calls: int = 1, private: bool = False,
               max_depth: int | None = None, modules: tuple[str, ...] = ()) -> str:
    """A flowchart, one subgraph per module, edges labelled with their call counts."""
    keep = _visible(rec, min_calls=min_calls, private=private, max_depth=max_depth,
                    modules=modules)
    ids: dict[str, str] = {}
    by_module: dict[str, list[str]] = {}
    for name in sorted(keep):
        ids[name] = f"n{len(ids)}"
        by_module.setdefault(_module(name), []).append(name)

    lines = ["flowchart TD"]
    for module, names in sorted(by_module.items()):
        safe = module.replace(".", "_")
        lines.append(f'  subgraph {safe}["{module}"]')
        for name in names:
            calls = rec.calls[name]
            label = _short(name) + (f" ×{calls}" if calls > 1 else "")
            lines.append(f'    {ids[name]}["{label}"]')
        lines.append("  end")

    edges = collapsed_edges(rec, keep)
    roots = {b for (a, b) in edges if a == "<entry>"}
    for callee in sorted(roots):
        lines.append(f"  START(( )) --> {ids[callee]}")
    for (caller, callee), n in sorted(edges.items(),
                                      key=lambda kv: (-kv[1], kv[0])):
        if caller == "<entry>":
            continue
        arrow = f"-- {n} -->" if n > 1 else "-->"
        lines.append(f"  {ids[caller]} {arrow} {ids[callee]}")
    return "\n".join(lines)


def collapsed_edges(rec: CallRecorder, keep: set[str]) -> Counter[tuple[str, str]]:
    """Edges between VISIBLE functions, with hidden frames collapsed through.

    Each recorded path is filtered down to the functions on the diagram and then read as
    consecutive pairs, so a call that passed through two hidden helpers becomes one edge
    rather than two dangling ones.  A visible function with no visible ancestor is an
    entry point and is attached to `START`.
    """
    edges: Counter[tuple[str, str]] = Counter()
    for path, n in rec.paths.items():
        if path[-1] not in keep:
            continue
        visible = [name for name in path if name in keep]
        edges[(visible[-2] if len(visible) > 1 else "<entry>", visible[-1])] += n
    return edges


def to_json(rec: CallRecorder, keep: set[str]) -> str:
    edges = collapsed_edges(rec, keep)
    return json.dumps({
        "calls": {n: c for n, c in rec.calls.most_common() if n in keep},
        "edges": [{"from": a, "to": b, "n": n} for (a, b), n in edges.most_common()],
    }, indent=2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("scenario", choices=sorted(SCENARIOS))
    ap.add_argument("--format", choices=("mermaid", "json"), default="mermaid")
    ap.add_argument("--min-calls", type=int, default=1,
                    help="drop functions called fewer times than this")
    ap.add_argument("--private", action="store_true",
                    help="include _private helpers (off by default: they dominate)")
    ap.add_argument("--max-depth", type=int, default=None,
                    help="keep only functions first reached at this call depth or less")
    ap.add_argument("--modules", default="",
                    help="comma-separated module prefixes to keep, e.g. assembly,sites")
    args = ap.parse_args()
    modules = tuple(m.strip() for m in args.modules.split(",") if m.strip())

    run = SCENARIOS[args.scenario]
    # Run once untraced first.  Imports execute module and class bodies, and those arrive
    # as `call` events indistinguishable from real work — a first-run graph is mostly the
    # import graph.  The warm pass also fills the `lru_cache`d reference tables, so what
    # gets recorded is the steady-state flow rather than one-off table loading.
    run()
    rec = CallRecorder(str(SRC / "mofsbu"))
    with rec:
        label = run()

    if args.format == "json":
        keep = _visible(rec, min_calls=args.min_calls, private=args.private,
                        max_depth=args.max_depth, modules=modules)
        print(to_json(rec, keep))
        return
    print(f"%% {args.scenario}: {label}")
    print(f"%% {len(rec.calls)} functions, {sum(rec.calls.values())} calls")
    print(to_mermaid(rec, min_calls=args.min_calls, private=args.private,
                     max_depth=args.max_depth, modules=modules))


if __name__ == "__main__":
    main()
