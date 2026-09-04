"""The integrity checker must catch the class of bug that produced it.

A verifier that only ever says "fine" is worse than none, so every check here is
exercised by breaking something on purpose and asserting it is reported.
"""
from __future__ import annotations

import pytest

import fixtures as fx
from mofsbu.registry import MethodSpec, Registry, put_geometry, put_structure
from mofsbu.registry.verify import verify
from mofsbu._types import Fidelity

RAW = MethodSpec(code="construct", code_version="0", method="raw-construct")


def xyz_for(g, jitter: float = 0.0) -> str:
    lines = [str(len(g)), "test"]
    for k, i in enumerate(g.nodes()):
        lines.append(f"{g.label(i).element} {k + jitter:.4f} 0.0000 0.0000")
    return "\n".join(lines) + "\n"


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        for build in (fx.cu_paddlewheel, fx.fe3_mu3_oxo, fx.btc, fx.water):
            g = build()
            sid = put_structure(r, g).id
            put_geometry(r, sid, xyz_for(g), fidelity=Fidelity.RAW, method=RAW)
        yield r


def checks(problems) -> set[str]:
    return {p.check for p in problems}


def test_a_registry_written_through_the_api_is_consistent(reg):
    assert verify(reg) == []


def test_a_hand_poked_column_is_caught(reg):
    """The actual bug: a row whose columns describe chemistry its graph does not."""
    reg.conn.execute("UPDATE structures SET has_metal_metal=0, n_dative_bonds=3, "
                     "max_bridge_class='muN' WHERE has_metal_metal=1")
    found = checks(verify(reg))
    assert {"has_metal_metal", "n_dative_bonds", "max_bridge_class"} <= found


def test_a_hand_written_label_is_caught(reg):
    """Labels are derived; a stored one that disagrees means something wrote it by hand."""
    reg.conn.execute("UPDATE structures SET display_label='the good one' WHERE id=1")
    assert "display_label" in checks(verify(reg))


def test_a_wrong_fragment_decomposition_is_caught(reg):
    reg.conn.execute("UPDATE structure_fragments SET count=count+1 WHERE structure_id=1")
    assert "structure_fragments" in checks(verify(reg))


def test_a_forged_identity_is_caught(reg):
    reg.conn.execute("UPDATE structures SET l1_graph_hash='0'*64 WHERE id=1")
    assert "l1_graph_hash" in checks(verify(reg))


def test_a_hand_written_l0_format_is_caught(reg):
    reg.conn.execute("UPDATE structures SET l0_composition='C4H6Cu2O9|q-2|m6' WHERE id=1")
    assert "l0_composition" in checks(verify(reg))


def test_a_missing_graph_blob_is_caught(reg):
    row = reg.conn.execute("SELECT typed_graph_hash FROM structures WHERE id=1").fetchone()
    reg.store.path(row["typed_graph_hash"]).unlink()
    assert "graph_blob_missing" in checks(verify(reg))


def test_a_missing_coordinate_blob_is_caught(reg):
    row = reg.conn.execute("SELECT coords_hash FROM geometries WHERE id=1").fetchone()
    reg.store.path(row["coords_hash"]).unlink()
    assert "coords_blob_missing" in checks(verify(reg))


def test_a_stale_best_geometry_pointer_is_caught(reg):
    reg.conn.execute("UPDATE structures SET best_geometry_id=NULL, best_fidelity=NULL WHERE id=1")
    assert "best_geometry" in checks(verify(reg))


def test_a_geometry_from_another_molecule_is_caught(reg):
    """put_geometry refuses this, so it can only arrive by writing around the API."""
    reg.conn.execute("UPDATE geometries SET n_atoms = n_atoms + 1 WHERE id=1")
    assert "geometry_atom_count" in checks(verify(reg))


def test_a_reordered_canonical_map_is_caught(reg):
    reg.conn.execute("UPDATE structures SET canonical_order_json='[0,1,2]' WHERE id=1")
    assert "canonical_order" in checks(verify(reg))


def test_a_version_bump_is_stale_not_corrupt(reg):
    """A deliberate algorithm change means 'recompute', not 'the data is wrong'."""
    reg.conn.execute("UPDATE structures SET algo_l1='cert0' WHERE id=1")
    problems = verify(reg)
    stale = [p for p in problems if p.check == "algo_version"]
    assert stale and all(p.severity == "stale" for p in stale)
    # and identity is NOT re-checked against the current recipe for that row
    assert not [p for p in problems if p.structure_id == 1 and p.check == "l1_graph_hash"]
