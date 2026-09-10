"""Labels are DERIVED FROM RETRIEVED ROWS.  They are never an input to retrieval.

The rule, and the reason for it: you search by structure — metals, fragments,
connectivity, charge — get rows back, and only then compose something readable.  Going
the other way means a query depends on someone having typed the right string, which
fails silently the moment structures are generated rather than authored, and fails
catastrophically once two different species compose the same name.  Nothing in
`registry.find` or the viewer's filters touches a label.

Everything here is a pure function of the typed graph plus an optional alias table.
An alias is display-only: it can be wrong, missing, or changed, and no query moves.

**Placeholders are deliberate and loud.**  Where the pipeline does not yet know
something, the label says so — `?L2` means configurational isomers are not being
distinguished yet (M5), so two rows carrying the same label may be different species.
A label that quietly omits what it cannot compute is a debugging trap.
"""
from __future__ import annotations

from dataclasses import dataclass

from mofsbu.graph._types import BridgeClass, TypedGraph
from mofsbu.identity import hill_formula, l1_graph_hash

UNNAMED = ""                       # no alias known: the fragment's formula is used
L2_PLACEHOLDER = "?L2"             # isomer tag not computed yet (M5)
L3_PLACEHOLDER = "?L3"             # conformer id not computed yet (M5)
# `SITES_PLACEHOLDER = "?sites"` lived here while `site_state` was unpopulated. M4's
# second half populates it (D18), so the placeholder is gone rather than left behind
# unreferenced — a fence around a hole that has been filled is just a wrong comment.

_MU = {BridgeClass.MU2: "mu2", BridgeClass.MU3: "mu3", BridgeClass.MU_N: "muN"}


@dataclass(frozen=True)
class Fragment:
    """One ligand fragment of a structure, identified by its own graph hash."""

    l1: str
    formula: str
    count: int
    bridge_class: str
    n_metals_bound: int

    def name(self, aliases: dict[str, str] | None = None) -> str:
        alias = (aliases or {}).get(self.l1)
        return alias if alias else self.formula


def fragment_formula(g: TypedGraph, atoms: tuple[int, ...]) -> str:
    counts: dict[str, int] = {}
    for i in atoms:
        el = g.label(i).element
        counts[el] = counts.get(el, 0) + 1
    return hill_formula(counts)


def fragment_hash(g: TypedGraph, atoms: tuple[int, ...]) -> str:
    """Identity of the fragment ALONE, so the same ligand hashes the same everywhere.

    The fragment is lifted out with its own charge and multiplicity left unset — they
    belong to the parent structure, not to a piece of it — so the hash is over
    connectivity and atom labels only.
    """
    sub = g.subgraph(atoms)
    sub.charge, sub.multiplicity = 0, 1
    return l1_graph_hash(sub)


def decompose(g: TypedGraph) -> list[Fragment]:
    """Group a structure's ligand fragments by identity.  Pure function of the graph."""
    tally: dict[str, dict] = {}
    for atoms in g.ligand_fragments():
        key = fragment_hash(g, atoms)
        metals = set()
        for i in atoms:
            metals.update(g.bridging_metals(i))
        entry = tally.get(key)
        if entry is None:
            tally[key] = {
                "l1": key, "formula": fragment_formula(g, atoms), "count": 1,
                "bridge_class": g.fragment_bridge_class(atoms).value,
                "n_metals_bound": len(metals),
            }
        else:
            entry["count"] += 1
    return [Fragment(**e) for e in
            sorted(tally.values(), key=lambda e: (-e["count"], e["formula"], e["l1"]))]


def metal_part(g: TypedGraph) -> str:
    counts: dict[str, int] = {}
    for i in g.metals():
        el = g.label(i).element
        counts[el] = counts.get(el, 0) + 1
    return "".join(f"{el}{counts[el]}" if counts[el] > 1 else el for el in sorted(counts))


def compose_label(
    g: TypedGraph,
    *,
    aliases: dict[str, str] | None = None,
    l2_tag: str = "",
    include_placeholders: bool = True,
) -> str:
    """Compose a readable label from the graph alone.

    No author-supplied name is consulted.  Anything the pipeline cannot yet compute is
    marked rather than omitted, so an incomplete label is visibly incomplete.
    """
    parts: list[str] = []
    metals = metal_part(g)
    if metals:
        parts.append(metals)

    for frag in decompose(g):
        prefix = ""
        bridge = frag.bridge_class
        if bridge in ("mu2", "mu3", "muN"):
            prefix = f"{bridge}-"
        suffix = f"{frag.count}" if frag.count > 1 else ""
        parts.append(f"[{prefix}{frag.name(aliases)}]{suffix}")

    if not parts:
        parts.append(hill_formula(g.element_counts()))

    charge = g.charge if g.charge is not None else "?"
    label = "".join(parts) + f" q{charge}"
    if g.multiplicity and g.multiplicity > 1:
        label += f" s{g.multiplicity}"
    if include_placeholders:
        label += f" {l2_tag}" if l2_tag else f" {L2_PLACEHOLDER}"
    return label


def unnamed_fragments(g: TypedGraph, aliases: dict[str, str] | None = None) -> list[str]:
    """Fragment hashes with no alias — what a naming pass would need to be given."""
    aliases = aliases or {}
    return [f.l1 for f in decompose(g) if f.l1 not in aliases]
