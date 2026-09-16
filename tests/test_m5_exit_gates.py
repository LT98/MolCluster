"""M5's exit gates, at the registry — the three claims the milestone is judged on.

`PLAN_implementation.md` §M5 lists them: one node reached by two routes, cis/trans surviving
as distinct records, and a conformer regenerable from its stored provenance. Earlier slices
showed each of these in memory; this file shows them after a round trip through SQLite,
which is where they actually have to hold.

**One substitution, stated rather than hidden.** The plan names the anthrarufin–Cu pair for
gate 2, and this file meets it with Pt(OH₂)₂Cl₂ — monodentate, buildable, and a genuine
cis/trans pair. The claim under test is unchanged: near-degenerate configurational isomers
persist as separate records carrying different choice vectors, and geometric clustering does
not merge them. The substitution is no longer *forced* — `join_chelate` places a two-point
ligand now (`tests/test_join.py`) — so re-deriving this gate on anthrarufin is available work
rather than blocked work.
"""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu._types import AmbiguousSpecError
from mofsbu.assembly.choice import ChoiceVector
from mofsbu.assembly.construct import (
    construct, enumerate_constructions, ligand_block, metal_block, resolve_geometry)
from mofsbu.assembly.join import join
from mofsbu.assembly.persist import store_block, store_construction
from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.identity.conformers import (
    DEFAULT_THETA_GEOM, Conformer, cluster, core_rmsd)
from mofsbu.registry import (
    Registry, geometries_from_choice, geometry_xyz, get_structure, incoming_routes)


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "gates.db") as r:
        r.migrate()
        yield r


@pytest.fixture(scope="module")
def parts():
    return {
        "aqua": ligand_block(embed_molecule(mol_from_smiles("O")), charge=0, name="aqua"),
        "chloride": ligand_block(embed_molecule(mol_from_smiles("[Cl-]")), charge=-1,
                                 name="chloride"),
    }


def platinum():
    return metal_block("Pt", 2, cn=4, geometry="square_planar", spin_class="ls")


def occupy(metal, pairs):
    block = metal
    for slot, ligand in pairs:
        vertex = next(s for s in block.open_vacancies() if s.slot == slot)
        block = join(block, ligand, vertex, ligand.open_donors()[0]).block
    return block


def xyz_coords(text: str) -> np.ndarray:
    return np.array([[float(x) for x in line.split()[1:4]]
                     for line in text.strip().splitlines()[2:]])


# ── gate 1: one node, two routes ─────────────────────────────────────────────

def test_one_structure_two_routes(reg, parts):
    """The D2 claim after a round trip: identity on the node, sequence on the edges.

    Two build orders, two provenance edges, ONE structures row — and the denormalised
    `n_incoming_routes` agrees with `incoming_routes`, which is the column a reader will
    actually reach for.
    """
    aqua, chloride = parts["aqua"], parts["chloride"]
    forward = occupy(platinum(), [(0, aqua), (2, chloride)])
    backward = occupy(platinum(), [(2, chloride), (0, aqua)])

    first = store_block(reg, forward, choice_vector={"order": "aqua-first"}, note="route A")
    second = store_block(reg, backward, choice_vector={"order": "chloride-first"},
                         note="route B")

    assert first.structure_id == second.structure_id
    assert first.structure_created and not second.structure_created

    routes = incoming_routes(reg, first.structure_id)
    assert len(routes) == 2
    assert get_structure(reg, first.structure_id)["n_incoming_routes"] == 2
    # ...and the two edges are distinguishable, which is what makes them worth keeping.
    assert len({r["choice_vector_digest"] for r in routes}) == 2


def test_the_two_routes_stored_two_geometries_under_one_identity(reg, parts):
    """One structure, many geometries (D4) — the fidelity ladder's shape, reached by routes."""
    aqua, chloride = parts["aqua"], parts["chloride"]
    a = store_block(reg, occupy(platinum(), [(0, aqua), (2, chloride)]),
                    choice_vector={"order": "aqua-first"})
    b = store_block(reg, occupy(platinum(), [(2, chloride), (0, aqua)]),
                    choice_vector={"order": "chloride-first"})
    n = reg.conn.execute("SELECT COUNT(*) FROM geometries WHERE structure_id=?",
                         (a.structure_id,)).fetchone()[0]
    assert a.structure_id == b.structure_id and n == 2


# ── gate 2: cis and trans survive ────────────────────────────────────────────

def test_cis_and_trans_are_two_records(reg, parts):
    """One L1, two rows — the whole point of L2 being in the identity key (D10)."""
    aqua, chloride = parts["aqua"], parts["chloride"]
    cis = store_block(reg, occupy(platinum(), [(0, aqua), (1, aqua),
                                               (2, chloride), (3, chloride)]),
                      choice_vector={"isomer": "cis"}, note="cis")
    trans = store_block(reg, occupy(platinum(), [(0, aqua), (2, aqua),
                                                 (1, chloride), (3, chloride)]),
                        choice_vector={"isomer": "trans"}, note="trans")

    assert cis.structure_id != trans.structure_id
    rows = [get_structure(reg, s.structure_id) for s in (cis, trans)]
    assert rows[0]["l1_graph_hash"] == rows[1]["l1_graph_hash"]      # L1 does NOT split them
    assert rows[0]["l2_isomer_tag"] != rows[1]["l2_isomer_tag"]      # L2 does
    assert "cis" in rows[0]["l2_isomer_tag"] and "trans" in rows[1]["l2_isomer_tag"]


def test_clustering_does_not_merge_them_however_near_degenerate(reg, parts):
    """D10's sharpest claim: the discriminator is relevance, NOT the energy gap.

    So the pair is handed identical energies and a wide-open energy window — if ΔE were
    doing the work, that would merge them. It does not, because ΔE never merges anything:
    the window only ever *prevents* a merge. What keeps them apart is that they are 1.9 Å
    apart over the rigid core, an order of magnitude past the calibrated θ_geom.

    (An earlier version of this test asserted they survive θ_geom = 99 Å as well. They do
    not, and should not: two choice vectors whose geometries really do coincide are the
    convergence case D11 asks clustering to reconcile. 99 Å is not "near-degenerate", it is
    "every structure is one structure".)
    """
    aqua, chloride = parts["aqua"], parts["chloride"]
    cis = occupy(platinum(), [(0, aqua), (1, aqua), (2, chloride), (3, chloride)])
    trans = occupy(platinum(), [(0, aqua), (2, aqua), (1, chloride), (3, chloride)])
    candidates = [
        Conformer({"isomer": "cis"}, cis.graph, cis.geometry, energy=-100.0),
        Conformer({"isomer": "trans"}, trans.graph, trans.geometry, energy=-100.0),
    ]
    assert len(cluster(candidates, theta_geom=DEFAULT_THETA_GEOM, energy_window=99.0)) == 2

    separation = core_rmsd(cis.graph, cis.geometry, trans.graph, trans.geometry)
    assert separation > 10 * DEFAULT_THETA_GEOM        # not a near miss on the threshold


# ── gate 3: replay ───────────────────────────────────────────────────────────

def test_a_stored_conformer_is_regenerable_from_its_row(reg, parts):
    """If this fails, stored conformers are frozen coordinates rather than objects.

    The round trip is the real one: build, store, read the row back, rebuild from the
    `choice_vector_json` and `seed` the row carries, and compare against the stored `.xyz`.
    """
    aqua = parts["aqua"]
    seed_block = platinum()
    tree = enumerate_constructions(seed_block, [aqua], degree=2, max_products=20)
    leaf = next(c for c in tree if len(c.steps) == 2)
    stored = store_construction(reg, leaf, reagent_ids=(), tags=("gate3",))

    row = reg.conn.execute(
        "SELECT choice_vector_json, seed, coords_hash FROM geometries WHERE id=?",
        (stored.geometry_id,)).fetchone()
    replayed = construct(seed_block, [aqua], ChoiceVector.from_json(
        row["choice_vector_json"]).data)

    original = xyz_coords(geometry_xyz(reg, stored.geometry_id))
    assert replayed.block.geometry == pytest.approx(original, abs=1e-6)
    assert replayed.digest == leaf.digest


def test_the_stored_digest_finds_the_geometry_it_built(reg, parts):
    """`ix_geometries_choice` doing the job it was created for three milestones ago."""
    aqua = parts["aqua"]
    leaf = enumerate_constructions(platinum(), [aqua], degree=1, max_products=5).leaves[0]
    stored = store_construction(reg, leaf)
    found = geometries_from_choice(reg, leaf.choice_vector)
    assert [r["id"] for r in found] == [stored.geometry_id]


# ── the negative the gate list ends on ───────────────────────────────────────

def test_an_ambiguous_spec_refuses_rather_than_defaulting():
    """Ground rule 5 at the gate: CN 4 does not determine a coordination geometry."""
    with pytest.raises(AmbiguousSpecError) as exc:
        resolve_geometry(4)
    assert "square_planar" in str(exc.value) and "tetrahedral" in str(exc.value)


def test_a_block_without_coordinates_is_not_stored_as_though_it_had_them(reg, parts):
    """A two-centre product is real and has no geometry until M6; storing one would invent it."""
    aqua = parts["aqua"]
    mono = occupy(platinum(), [(0, aqua)])
    second = metal_block("Cu", 2, cn=4, geometry="square_planar", spin_class="ls")
    graph_only = join(mono, second, mono.open_donors()[0] if mono.open_donors() else
                      next(s for s in mono.sites if not s.is_vacancy),
                      second.open_vacancies()[0], with_geometry=False).block
    with pytest.raises(ValueError, match="would invent them"):
        store_block(reg, graph_only, choice_vector={"op": "join"})
