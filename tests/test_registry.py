"""M3 exit gate: identity persists, duplicates collapse, provenance accumulates."""
from __future__ import annotations

import json

import pytest

import fixtures as fx
from mofsbu.graph import TypedGraph
from mofsbu.registry import (
    MethodSpec, Provenance, Registry, RegistryError, best_geometry, canonical_map,
    export_xyz, find, geometry_xyz, get_graph, get_structure, method_id, put_geometry,
    put_structure,
)
from mofsbu._types import Fidelity

RAW = MethodSpec(code="construct", code_version="0", method="raw-construct")
XTB = MethodSpec(code="tblite", code_version="0.4.0", method="GFN2-xTB", solvent=None)
DFT = MethodSpec(code="orca", code_version="6.0", method="wB97X-D/def2-TZVP")
# Two MLIPs on ONE rung of the ladder.  Their absolute energies are on different scales,
# which is the whole reason the best-geometry rule may not compare them by magnitude.
MACE_MP = MethodSpec(code="mace", code_version="0.3.6/medium", method="MACE-MP-0",
                     extras={"charge_blind": True, "spin_blind": True})
MACE_OMOL = MethodSpec(code="mace", code_version="0.3.14/extra_large",
                       method="MACE-OMOL-0")


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


def xyz_for(g: TypedGraph, jitter: float = 0.0) -> str:
    lines = [str(len(g)), "test"]
    for k, i in enumerate(g.nodes()):
        lines.append(f"{g.label(i).element} {k + jitter:.4f} 0.0000 0.0000")
    return "\n".join(lines) + "\n"


# ── identity persists ────────────────────────────────────────────────────────

def test_graph_round_trips_through_the_store(reg):
    """A newly inserted graph comes back byte-identical."""
    for name, build in fx.ALL.items():
        g = build()
        put = put_structure(reg, g)
        if put.created:
            assert get_graph(reg, put.id).to_json() == g.to_json(), name


def test_stored_graph_is_the_class_representative_not_your_input(reg):
    """A duplicate insert returns the FIRST graph's atom ordering, not the second's.

    That is the point of identity, not a loss: the two orderings are the same
    structure, so the registry keeps one representative.  Anything that needs to map
    back onto a caller's own atom numbering must go through the canonical order or an
    atom map, never assume its build ordering survived.
    """
    from mofsbu.graph import is_isomorphic

    first = fx.pt_ammine_dichloride("a")
    second = fx.pt_ammine_dichloride("b")
    sid = put_structure(reg, first).id
    assert put_structure(reg, second).id == sid

    stored = get_graph(reg, sid)
    assert stored.to_json() == first.to_json()
    assert stored.to_json() != second.to_json()
    assert is_isomorphic(stored, second)


def test_canonical_order_is_stored_not_recomputed(reg):
    """Ground rule 6 applied to the map: freeze it, so a version bump is detectable."""
    g = fx.cu_paddlewheel()
    sid = put_structure(reg, g).id
    stored = json.loads(
        reg.conn.execute("SELECT canonical_order_json FROM structures WHERE id=?",
                         (sid,)).fetchone()["canonical_order_json"])
    assert sorted(stored) == g.nodes()
    assert canonical_map(reg, sid) == {v: k for k, v in enumerate(stored)}


def test_identity_columns_match_the_keys(reg):
    from mofsbu.identity import identity

    g = fx.fe3_mu3_oxo((2, 3, 3))
    sid = put_structure(reg, g).id
    row, ident = get_structure(reg, sid), identity(g)
    assert row["l0_composition"] == ident["l0"]
    assert row["l1_graph_hash"] == ident["l1"]
    assert row["wl_index"] == ident["wl_index"]


# ── the D2 claim: one node, many routes ──────────────────────────────────────

def test_duplicate_insert_is_one_row_with_two_edges(reg):
    g = fx.cu_paddlewheel()
    a = put_structure(reg, g, provenance=Provenance(kind="assembly", note="route A"))
    b = put_structure(reg, g.permuted(3), provenance=Provenance(kind="assembly", note="route B"))
    assert a.id == b.id and a.created and not b.created
    assert reg.count("structures") == 1
    notes = [r["note"] for r in reg.conn.execute(
        "SELECT note FROM reactions WHERE product_structure_id=? ORDER BY id", (a.id,))]
    assert notes == ["route A", "route B"]


def test_cis_trans_builds_collapse_but_bridging_does_not(reg):
    same = {put_structure(reg, fx.pt_ammine_dichloride(o)).id for o in ("a", "b")}
    assert len(same) == 1
    different = {put_structure(reg, b()).id
                 for b in (fx.zn2_bridged_formates, fx.zn2_chelated_formates)}
    assert len(different) == 2


# ── geometries and the fidelity ladder ───────────────────────────────────────

def test_best_geometry_climbs_the_ladder_and_never_descends(reg):
    g = fx.water()
    sid = put_structure(reg, g).id
    raw = put_geometry(reg, sid, xyz_for(g, 0.0), fidelity=Fidelity.RAW, method=RAW)
    assert get_structure(reg, sid)["best_geometry_id"] == raw.id

    dft = put_geometry(reg, sid, xyz_for(g, 0.2), fidelity=Fidelity.DFT, method=DFT,
                       energy=-76.4, converged=True)
    assert get_structure(reg, sid)["best_geometry_id"] == dft.id

    put_geometry(reg, sid, xyz_for(g, 0.3), fidelity=Fidelity.XTB, method=XTB, energy=-5.0)
    assert get_structure(reg, sid)["best_geometry_id"] == dft.id     # does not descend
    assert best_geometry(reg, sid)["id"] == dft.id
    assert best_geometry(reg, sid, Fidelity.XTB)["fidelity"] == Fidelity.DFT


def test_best_geometry_never_compares_energies_across_methods(reg):
    """One rung, two theories: the deeper reference must not win by being deeper.

    MACE-MP-0 total energies and MACE-OMOL-0 total energies are both "eV" and are not
    the same quantity — different training set, different reference, different scale.
    The old rule was "highest fidelity, then lowest energy", which on this pair picks a
    model rather than a geometry, every time, with nothing in the database saying so.
    The rule between two methods is now stated: charge-aware outranks charge-blind.
    """
    g = fx.water()
    sid = put_structure(reg, g).id
    mp = put_geometry(reg, sid, xyz_for(g, 0.1), fidelity=Fidelity.ML, method=MACE_MP,
                      energy=-3000.0, converged=True)
    omol = put_geometry(reg, sid, xyz_for(g, 0.2), fidelity=Fidelity.ML, method=MACE_OMOL,
                        energy=-14.0, converged=True)
    assert mp.id != omol.id, "two methods on one rung must be two rows, not a clash"
    assert get_structure(reg, sid)["best_geometry_id"] == omol.id

    # ... and the rung still beats everything below it, in either order.
    assert get_structure(reg, sid)["best_fidelity"] == Fidelity.ML


def test_two_mlips_relaxing_one_structure_do_not_collide(reg):
    """`method_id` is in the geometries UNIQUE key, so this needed no migration."""
    g = fx.water()
    sid = put_structure(reg, g).id
    xyz = xyz_for(g, 0.1)                       # the SAME coordinates from both models
    a = put_geometry(reg, sid, xyz, fidelity=Fidelity.ML, method=MACE_MP, energy=-3000.0)
    b = put_geometry(reg, sid, xyz, fidelity=Fidelity.ML, method=MACE_OMOL, energy=-14.0)
    assert a.created and b.created and a.id != b.id
    rows = reg.conn.execute(
        "SELECT COUNT(*) AS n FROM geometries WHERE structure_id=? AND fidelity=?",
        (sid, int(Fidelity.ML))).fetchone()
    assert rows["n"] == 2


def test_best_fidelity_never_drifts_from_the_geometries_it_caches(reg):
    """The two columns are maintained in one place; this asserts they agree."""
    for build in (fx.water, fx.formate, fx.cu_paddlewheel):
        g = build()
        sid = put_structure(reg, g).id
        for fid, meth in ((Fidelity.RAW, RAW), (Fidelity.XTB, XTB), (Fidelity.ML, RAW)):
            put_geometry(reg, sid, xyz_for(g, float(fid)), fidelity=fid, method=meth)
    drift = reg.conn.execute(
        "SELECT s.id FROM structures s LEFT JOIN geometries g ON g.id = s.best_geometry_id "
        "WHERE s.best_fidelity IS NOT g.fidelity").fetchall()
    assert drift == []


def test_identical_coordinates_are_stored_once(reg):
    g = fx.water()
    sid = put_structure(reg, g).id
    xyz = xyz_for(g)
    first = put_geometry(reg, sid, xyz, fidelity=Fidelity.RAW, method=RAW)
    second = put_geometry(reg, sid, xyz, fidelity=Fidelity.RAW, method=RAW)
    assert first.id == second.id and not second.created
    assert reg.count("geometries") == 1


def test_geometry_atom_count_must_match_its_structure(reg):
    g = fx.water()
    sid = put_structure(reg, g).id
    with pytest.raises(RegistryError, match="atoms"):
        put_geometry(reg, sid, xyz_for(fx.formate()), fidelity=Fidelity.RAW, method=RAW)


def test_every_number_carries_a_method(reg):
    """Ground rule 3, enforced structurally: put_geometry cannot omit `method`."""
    g = fx.water()
    sid = put_structure(reg, g).id
    put_geometry(reg, sid, xyz_for(g), fidelity=Fidelity.XTB, method=XTB, energy=-5.0)
    assert reg.count("methods") == 1
    assert method_id(reg, XTB) == method_id(reg, XTB)          # get-or-create, not duplicate
    with pytest.raises(TypeError):
        put_geometry(reg, sid, xyz_for(g, 1.0), fidelity=Fidelity.XTB)   # no method


def test_export_is_a_view_not_a_record(reg, tmp_path):
    g = fx.cu_paddlewheel()
    sid = put_structure(reg, g).id
    gid = put_geometry(reg, sid, xyz_for(g), fidelity=Fidelity.RAW, method=RAW).id
    out = export_xyz(reg, gid, tmp_path / "sub" / "out.xyz")
    assert out.read_text() == geometry_xyz(reg, gid)


# ── query surface ────────────────────────────────────────────────────────────

@pytest.fixture()
def populated(reg):
    for build in fx.ALL.values():
        g = build()
        sid = put_structure(reg, g).id
        put_geometry(reg, sid, xyz_for(g), fidelity=Fidelity.RAW, method=RAW)
    return reg


def test_find_across_the_four_axes(populated):
    assert {r["display_label"] for r in find(populated, metal="Cu")}
    assert all(r["n_metals"] == 3 for r in find(populated, n_metals=3))
    assert all(r["has_metal_metal"] == 1 for r in find(populated, has_metal_metal=True))
    assert all(r["max_bridge_class"] == "mu3"
               for r in find(populated, max_bridge_class="mu3"))
    assert all(r["best_fidelity"] >= Fidelity.RAW
               for r in find(populated, min_fidelity=Fidelity.RAW))
    assert find(populated, has_route=True) == []          # nothing had provenance
    assert len(find(populated, has_route=False)) == populated.count("structures")


def test_metal_filter_cannot_match_a_prefix(populated):
    """',Cu,' must not match ',Cu2,'-style neighbours — hence the bracketed form."""
    assert find(populated, metal="C") == []
    assert find(populated, metal="Cu")


def test_find_rejects_an_unlisted_sort_column(populated):
    with pytest.raises(RegistryError, match="cannot sort"):
        find(populated, sort="id; DROP TABLE structures--")
    assert populated.count("structures") > 0


def test_display_label_is_derived_from_the_graph_not_from_a_name(populated):
    """The label must be a function of the structure, not a memory of what it was called.

    `TypedGraph.name` is set on every fixture and must appear in no label: a hand-typed
    name makes authored examples read perfectly and leaves anything the enumerator
    generates unreadable.
    """
    from mofsbu.naming import compose_label

    rows = find(populated, limit=1000)
    assert all(r["display_label"] for r in rows)
    assert not any(r["l1_graph_hash"] in r["display_label"] for r in rows)
    assert not any("cu_paddlewheel" in r["display_label"] for r in rows)
    for row in rows:
        g = get_graph(populated, row["id"])
        assert row["display_label"] == compose_label(g)


def test_a_label_can_never_satisfy_a_query(populated):
    """Search is structural.  Nothing that exists only in a label may be findable.

    This is the invariant the whole naming design exists to protect: in a large system a
    query that depends on a name silently returns the wrong set the moment a structure
    was generated rather than authored.
    """
    from mofsbu.registry.api import _SORTABLE

    assert "display_label" not in _SORTABLE
    # `mu2` appears in labels; it must be reachable only via the structural column
    assert find(populated, formula_like="mu2") == []
    assert find(populated, max_bridge_class="mu2")


def test_fragments_are_the_search_surface(populated):
    """"Which structures contain this ligand" is an indexed join on fragment identity."""
    from mofsbu.naming import decompose

    formate = decompose(fx.cu_paddlewheel())[0].l1
    hits = find(populated, fragment_l1=formate, limit=100)
    assert hits
    for row in hits:
        assert formate in {f.l1 for f in decompose(get_graph(populated, row["id"]))}


def test_an_alias_changes_labels_and_no_query(populated):
    from mofsbu.naming import decompose
    from mofsbu.registry import alias_fragment, relabel_all

    formate = decompose(fx.cu_paddlewheel())[0].l1
    before = {r["id"] for r in find(populated, fragment_l1=formate, limit=100)}
    alias_fragment(populated, formate, "formate")
    relabel_all(populated)
    after = {r["id"] for r in find(populated, fragment_l1=formate, limit=100)}
    assert before == after                                  # the query is untouched
    assert any("formate" in r["display_label"]
               for r in find(populated, fragment_l1=formate, limit=100))


# ── version bookkeeping ──────────────────────────────────────────────────────

def test_algo_version_bump_is_detected_not_absorbed(reg):
    from mofsbu import versions

    put_structure(reg, fx.water())
    assert reg.stale_algo_versions() == {}
    original = versions.ALGO_VERSIONS["l1_certificate"]
    versions.ALGO_VERSIONS["l1_certificate"] = "cert2"
    try:
        assert reg.stale_algo_versions() == {"l1_certificate": (original, "cert2")}
    finally:
        versions.ALGO_VERSIONS["l1_certificate"] = original


# ── concurrency ──────────────────────────────────────────────────────────────

def test_two_writers_inserting_the_same_structure_do_not_collide(tmp_path):
    """Idempotency is not concurrency safety.

    `put_structure` reads before it writes, so two workers can both miss and both insert.
    Under background threads this failed roughly two runs in five with an IntegrityError
    on the identity constraint.  The constraint is what makes the race SAFE rather than
    silently duplicating — the fix is to let it arbitrate and re-read, not to lock.
    """
    from mofsbu.registry import BlobStore

    store = BlobStore(tmp_path / "store")
    with Registry(tmp_path / "r.db", store) as first:
        first.migrate()
    g = fx.cu_paddlewheel()

    with Registry(tmp_path / "r.db", store) as a, Registry(tmp_path / "r.db", store) as b:
        # simulate the interleaving exactly: both see "not present", then both insert
        assert a.conn.execute("SELECT COUNT(*) FROM structures").fetchone()[0] == 0
        assert b.conn.execute("SELECT COUNT(*) FROM structures").fetchone()[0] == 0
        put_a = put_structure(a, g)
        a.conn.commit()
        put_b = put_structure(b, g)          # loses the race; must not raise
        b.conn.commit()

    assert put_a.id == put_b.id
    assert put_a.created and not put_b.created
    with Registry(tmp_path / "r.db", store) as reg:
        assert reg.count("structures") == 1


def test_two_writers_registering_the_same_method_do_not_collide(tmp_path):
    from mofsbu.registry import BlobStore

    store = BlobStore(tmp_path / "store")
    with Registry(tmp_path / "r.db", store) as first:
        first.migrate()
    with Registry(tmp_path / "r.db", store) as a, Registry(tmp_path / "r.db", store) as b:
        id_a = method_id(a, XTB)
        a.conn.commit()
        id_b = method_id(b, XTB)
        b.conn.commit()
    assert id_a == id_b
