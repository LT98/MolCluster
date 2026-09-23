"""Which structures are the same molecule with a proton removed, and what that costs.

The registry stores every protomer of a ligand as a first-class row and records nothing
connecting them: `runner`'s `ligand` branch calls `put_structure` with tags and no
provenance, so tHQ and tHQ(-1H) are two unrelated structures that happen to share a name.
The relationship is recoverable — one is the other minus an H and minus a charge — and
once it is written down as an edge the existing reference scheme prices it.

Why a couple rather than a bare proton: `AH -> A- + H+` cannot be weighed. A proton has
no electrons, so no gas-phase method has an energy for it and `check_balance` reads charge
off a graph with no atoms. `AH + H2O -> A- + H3O+` balances, and because no metal-donor
bond changes across the arrow it is **isodesmic** — which every assembly edge in this
project fails to be. These are the first numbers here that `strict=True` accepts.
"""
from __future__ import annotations

from typing import Any

from mofsbu._types import MofsbuError
from mofsbu.energy.reference import put_balanced_reaction
from mofsbu.graph._types import EdgeType, TypedGraph
from mofsbu.identity import l1_graph_hash
from mofsbu.registry.api import get_graph

__all__ = ["deprotonation_pairs", "labile_hydrogens", "link_protomers",
           "PROTON_DONOR", "PROTON_ACCEPTOR"]

#: The couple, by the `TypedGraph.name` `examples` gives them.
PROTON_DONOR = "water"          # takes the proton away: H2O -> H3O+
PROTON_ACCEPTOR = "hydronium"


def labile_hydrogens(g: TypedGraph) -> list[int]:
    """Every H bonded to one non-carbon, non-metal atom: the protons a graph can lose."""
    out = []
    for i in g.nodes():
        if g.label(i).element != "H":
            continue
        nbrs = g.neighbors(i)
        if (len(nbrs) == 1 and not g.is_metal(nbrs[0])
                and g.label(nbrs[0]).element != "C"
                and g.edge_type(i, nbrs[0]) is EdgeType.COVALENT):
            out.append(i)
    return out


def deprotonation_pairs(reg: Any) -> list[tuple[int, int, int]]:
    """`(protonated_id, deprotonated_id, 1)` for every pair exactly one proton apart.

    Exact, not by formula: each labile H of the protonated graph is removed
    (`TypedGraph.without_proton`) and the result is looked up by L1 hash, charge and
    multiplicity. So an edge never also moves which atom binds the metal, and the atoms
    that differ are exactly one H — which is what makes its dE a deprotonation energy.
    Two protons apart is two edges through the intermediate, never one.
    """
    rows = [r for r in reg.conn.execute(
        "SELECT id, l1_graph_hash, net_charge, multiplicity FROM structures WHERE hidden = 0")]
    index: dict[tuple[str, int, int], list[int]] = {}
    for r in rows:
        index.setdefault((r["l1_graph_hash"], int(r["net_charge"] or 0),
                          int(r["multiplicity"] or 1)), []).append(int(r["id"]))

    out: list[tuple[int, int, int]] = []
    for r in rows:
        g = get_graph(reg, int(r["id"]))
        seen: set[str] = set()
        for h in labile_hydrogens(g):
            # The charge left behind is recorded on the atom or delocalised, depending on
            # the molecule (`from_rdkit`); the rest of the graph must match exactly either way.
            for localised in (True, False):
                key = l1_graph_hash(g.without_proton(h, localised=localised))
                if key in seen:       # a symmetry-equivalent proton: same product
                    continue
                seen.add(key)
                for partner in index.get((key, int(r["net_charge"] or 0) - 1,
                                          int(r["multiplicity"] or 1)), ()):
                    out.append((int(r["id"]), partner, 1))
    return out


def link_protomers(reg: Any, *, water_id: int, hydronium_id: int,
                   limit: int | None = None) -> list[dict[str, Any]]:
    """Write one balanced deprotonation edge per pair, skipping what already exists.

    Returns a row per pair with the reaction id, or the reason it was refused. A pair
    whose equation does not balance on the graphs is reported, not raised: the point is
    to link what can be linked, and to say what could not.
    """
    existing = _already_linked(reg)
    # The couple is itself a protomer pair, and linking it through itself states
    # `H3O+ + H2O -> H2O + H3O+` — balanced, isodesmic, and exactly zero. True and empty.
    couple = {int(water_id), int(hydronium_id)}
    written: list[dict[str, Any]] = []
    for protonated, deprotonated, n in deprotonation_pairs(reg):
        if (protonated, deprotonated) in existing or {protonated, deprotonated} & couple:
            continue
        if limit is not None and len(written) >= limit:
            break
        entry: dict[str, Any] = {"protonated": protonated,
                                 "deprotonated": deprotonated, "n": n,
                                 "reaction_id": None, "why_not": None}
        try:
            entry["reaction_id"] = put_balanced_reaction(
                reg, deprotonated,
                reagents=[(protonated, 1)],
                # The couple is the proton's CARRIER, not a piece being built in: water
                # takes the H+ away and leaves as H3O+.  Recording it as a reagent makes
                # a reader see "added water", which is true of the arithmetic and wrong
                # about the chemistry.  `solvent` is reagent-side, so balance is
                # unchanged (REAGENT_SIDE_ROLES).
                solvents=[(water_id, n)],
                leaving=[(hydronium_id, n)],
                kind="deprotonation",
                note=f"deprotonation x{n} (H2O/H3O+ couple)")
            existing.add((protonated, deprotonated))
        except MofsbuError as exc:
            entry["why_not"] = str(exc).splitlines()[0]
        written.append(entry)
    return written


def _already_linked(reg: Any) -> set[tuple[int, int]]:
    """Pairs a deprotonation edge already connects, so re-running writes nothing new."""
    out: set[tuple[int, int]] = set()
    for row in reg.conn.execute(
            "SELECT id, product_structure_id FROM reactions WHERE kind = 'deprotonation'"):
        # Only the acid is a `reagent`; the couple is solvent + leaving.
        for r in reg.conn.execute(
                "SELECT structure_id FROM reaction_reagents "
                "WHERE reaction_id = ? AND role = 'reagent'", (row["id"],)):
            out.add((int(r["structure_id"]), int(row["product_structure_id"])))
    return out
