"""D30 / B28: a geometry keeps its rows as built and records their link to the graph.

Two builds of one species number its atoms differently (two SMILES spellings, a join versus a
placement).  Each geometry's coordinates are stored exactly as its build wrote them, with the
link from canonical positions to its rows; sites are keyed by canonical position through the
build's own map; every reader goes through the link.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import fixtures as fx
from mofsbu._types import Fidelity
from mofsbu.assembly.persist import load_block
from mofsbu.geometry.connectivity import connectivity_check
from mofsbu.geometry.embed import embed_molecule, to_xyz
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles
from mofsbu.registry import (MethodSpec, Registry, RegistryError, geometry_graph,
                             geometry_link, geometry_xyz, method_id, put_geometry,
                             put_site_state, put_structure)
from mofsbu.registry.recover import recover_atom_links
from mofsbu.registry.verify import verify
from mofsbu.runner import _perceive_sites, _record_sites

FF = MethodSpec(code="rdkit", code_version="test", method="ETKDGv3+MMFF")
OMOL = MethodSpec(code="mace", code_version="test", method="MACE-OMOL-0")
RAW = MethodSpec(code="construct", code_version="0", method="raw-construct")


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


def parse(xyz: str) -> tuple[list[str], np.ndarray]:
    lines = xyz.strip().splitlines()
    rows = [ln.split() for ln in lines[2:2 + int(lines[0])]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def build(reg, smiles: str, name: str):
    """What the runner's ligand task does: embed, perceive, store with the build's link."""
    mol = embed_molecule(mol_from_smiles(smiles), seed=7)
    g = from_rdkit(mol, charge=0, multiplicity=1, name=name)
    perceived = _perceive_sites(mol, g, Fidelity.FF)
    put = put_structure(reg, g)
    text = to_xyz(mol, name)
    geom = put_geometry(reg, put.id, text, fidelity=Fidelity.FF, method=FF,
                        canonical_order=put.canonical_order)
    _record_sites(reg, put.id, geom.id, perceived, Fidelity.FF,
                  canonical_order=put.canonical_order)
    return put, geom, text


@pytest.mark.parametrize("first,second", [("OC1CCCCC1", "C1CC(O)CCC1"),
                                          ("CC(=O)O", "OC(C)=O")])
def test_a_second_spelling_is_stored_as_built_and_linked(reg, first, second):
    a, _ga, _ = build(reg, first, "L")
    b, gb, text_b = build(reg, second, "L")
    assert b.id == a.id and not b.created
    assert geometry_xyz(reg, gb.id) == text_b                 # coordinates exactly as built
    sym, xyz = parse(text_b)
    graph = geometry_graph(reg, gb.id)                        # numbered like gb's rows
    assert [graph.label(k).element for k in graph.nodes()] == sym
    assert not connectivity_check(sym, xyz, graph).broken


def test_site_states_of_a_second_build_sit_on_its_own_rows(reg):
    build(reg, "CC(=O)O", "AcOH")
    b, gb, text_b = build(reg, "OC(C)=O", "AcOH")
    sym, xyz = parse(text_b)
    link = geometry_link(reg, gb.id)
    rows = reg.conn.execute(
        "SELECT sc.canonical_idx, ss.frame_json FROM site_state ss "
        "JOIN site_catalog sc ON sc.id = ss.site_id WHERE ss.geometry_id=? "
        "AND ss.frame_json IS NOT NULL", (gb.id,)).fetchall()
    assert rows
    for r in rows:
        row = link[r["canonical_idx"]]
        assert sym[row] == "O"
        assert np.allclose(json.loads(r["frame_json"])["origin"], xyz[row], atol=1e-3)


def test_e1_judges_a_second_build_and_a_relaxation_inherits_its_link(reg):
    build(reg, "OC1CCCCC1", "CyOH")
    b, gb, text_b = build(reg, "C1CC(O)CCC1", "CyOH")
    verdict = reg.conn.execute("SELECT connectivity_changed FROM geometries WHERE id=?",
                               (gb.id,)).fetchone()[0]
    assert verdict == 0
    sym, xyz = parse(text_b)
    relaxed_text = "\n".join([str(len(sym)), "relaxed"] + [
        f"{s} {x + 0.01:.5f} {y:.5f} {z:.5f}" for s, (x, y, z) in zip(sym, xyz)]) + "\n"
    relaxed = put_geometry(reg, b.id, relaxed_text, fidelity=Fidelity.ML, method=OMOL,
                           energy=-1.0, converged=True, relaxed_from=gb.id)
    assert geometry_link(reg, relaxed.id) == geometry_link(reg, gb.id)
    assert reg.conn.execute("SELECT connectivity_changed FROM geometries WHERE id=?",
                            (relaxed.id,)).fetchone()[0] == 0


def test_load_block_of_a_second_build_puts_its_donor_sites_on_oxygen(reg):
    build(reg, "CC(=O)O", "AcOH")
    b, gb, _ = build(reg, "OC(C)=O", "AcOH")
    block = load_block(reg, b.id, gb.id)
    donors = [s for s in block.sites if s.role == "donor"]
    assert donors and all(block.graph.label(s.atom_idx).element == "O" for s in donors)
    for s in donors:
        if s.frame:
            assert np.allclose(s.frame["origin"], block.geometry[s.atom_idx], atol=1e-3)


# ── rows stored before the link existed ──────────────────────────────────────

def _legacy_second_build(reg):
    """What the pre-D30 writer left: a second build's rows and states, no link."""
    first, _g1, _ = build(reg, "CC(=O)O", "AcOH")
    mol = embed_molecule(mol_from_smiles("OC(C)=O"), seed=7)
    g = from_rdkit(mol, charge=0, multiplicity=1, name="AcOH")
    _perceived, sites, states = _perceive_sites(mol, g, Fidelity.FF)
    sid = put_structure(reg, g).id
    text = to_xyz(mol, "AcOH")

    def insert(xyz_text, fidelity, method, relaxed_from=None, cv=None, energy=None):
        cur = reg.conn.execute(
            "INSERT INTO geometries (structure_id, coords_hash, n_atoms, fidelity, method_id, "
            "energy, converged, relaxed_from, choice_vector_json, created_at) "
            "VALUES (?,?,?,?,?,?,?,?,?,?)",
            (sid, reg.store.put_text(xyz_text), len(g), int(fidelity), method_id(reg, method),
             energy, 1, relaxed_from, None if cv is None else json.dumps(cv), "2026-01-01"))
        return int(cur.lastrowid)

    gid = insert(text, Fidelity.FF, FF)
    put_site_state(reg, sid, gid, states, fidelity=Fidelity.FF,      # the B28 write: no link
                   frames={(x.atom_idx, x.slot): x.frame for x in sites if x.frame})
    sym, xyz = parse(text)
    moved = "\n".join([str(len(sym)), "r"] + [f"{s} {x + 0.01:.5f} {y:.5f} {z:.5f}"
                                              for s, (x, y, z) in zip(sym, xyz)]) + "\n"
    child = insert(moved, Fidelity.ML, OMOL, relaxed_from=gid, energy=-1.0)
    reg.conn.execute("UPDATE geometries SET canonical_order_json=NULL, "
                     "connectivity_changed=NULL, connectivity_json=NULL")
    reg.conn.commit()
    return sid, gid, child


def test_recovery_dry_run_counts_and_writes_nothing(reg):
    _sid, gid, child = _legacy_second_build(reg)
    out = recover_atom_links(reg)
    assert (out["stored_numbering"], out["other_numbering"], out["unmatched"]) == (1, 2, [])
    assert reg.conn.execute("SELECT COUNT(*) FROM geometries "
                            "WHERE canonical_order_json IS NOT NULL").fetchone()[0] == 0
    with pytest.raises(RegistryError, match="recover_atom_links"):
        geometry_graph(reg, gid)


def test_recovery_links_without_touching_coordinates_and_drops_miskeyed_states(reg):
    sid, gid, child = _legacy_second_build(reg)
    hashes = dict(reg.conn.execute("SELECT id, coords_hash FROM geometries").fetchall())
    out = recover_atom_links(reg, dry_run=False)
    assert out["other_numbering"] == 2 and out["states_dropped"] > 0
    assert dict(reg.conn.execute("SELECT id, coords_hash FROM geometries").fetchall()) == hashes
    for g_id in (gid, child):
        sym, xyz = parse(geometry_xyz(reg, g_id))
        assert not connectivity_check(sym, xyz, geometry_graph(reg, g_id)).broken
    assert geometry_link(reg, child) == geometry_link(reg, gid)
    assert reg.conn.execute("SELECT COUNT(*) FROM site_state WHERE geometry_id=?",
                            (gid,)).fetchone()[0] == 0
    assert recover_atom_links(reg, dry_run=False)["geometries"] == 0


def test_a_pulled_back_start_is_linked_only_in_the_stored_numbering(reg):
    g = fx.hexaaqua()
    put = put_structure(reg, g)
    coords = parse(fx.realised_xyz(g))[1]
    axis = coords[1] / np.linalg.norm(coords[1])
    pulled = coords.copy()
    pulled[[1, 2, 3]] += 1.5 * axis

    def xyz_of(c):
        return "\n".join([str(len(g)), "p"] + [f"{g.label(k).element} {x:.5f} {y:.5f} {z:.5f}"
                                               for k, (x, y, z) in enumerate(c)]) + "\n"

    def legacy(c):
        cur = reg.conn.execute(
            "INSERT INTO geometries (structure_id, coords_hash, n_atoms, fidelity, method_id, "
            "choice_vector_json, created_at) VALUES (?,?,?,?,?,?,?)",
            (put.id, reg.store.put_text(xyz_of(c)), len(g), int(Fidelity.RAW),
             method_id(reg, RAW), json.dumps({"perturb": {"by": 1.5}}), "2026-01-01"))
        return int(cur.lastrowid)

    ok = legacy(pulled)
    scrambled = pulled.copy()
    scrambled[[1, 4]] = scrambled[[4, 1]]                 # two water oxygens' rows swapped
    bad = legacy(scrambled)
    reg.conn.commit()
    out = recover_atom_links(reg)
    assert out["stored_numbering"] == 1 and out["unmatched"] == [bad] and ok != bad


def test_verify_reports_an_unlinked_geometry_and_a_wrong_link(reg):
    _sid, gid, _child = _legacy_second_build(reg)
    stale = [p for p in verify(reg) if p.check == "geometry_link"]
    assert stale and all(p.severity == "stale" for p in stale)
    recover_atom_links(reg, dry_run=False)
    assert not [p for p in verify(reg) if p.check == "geometry_link"]
    link = geometry_link(reg, gid)
    sym, _ = parse(geometry_xyz(reg, gid))
    a = 0
    b = next(c for c in range(len(link)) if sym[link[c]] != sym[link[a]])
    link[a], link[b] = link[b], link[a]                   # an O row where a C belongs
    reg.conn.execute("UPDATE geometries SET canonical_order_json=? WHERE id=?",
                     (json.dumps(link), gid))
    errors = [p for p in verify(reg) if p.check == "geometry_link"]
    assert errors and errors[0].severity == "error"
