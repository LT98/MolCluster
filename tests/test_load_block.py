"""`persist.load_block`: a stored structure read back as a joinable block (`store_block` reversed)."""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu._types import Fidelity, NotBuiltYet
from mofsbu.assembly.construct import ligand_block, metal_block
from mofsbu.assembly.join import bridge_compatible, join, join_bridge
from mofsbu.assembly.persist import CONSTRUCT, load_block, store_block, to_xyz
from mofsbu.geometry._linalg import axis_rotation
from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.registry import BlobStore, Registry, put_geometry


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "store")) as r:
        r.migrate("test")
        yield r


def _formate():
    lig = ligand_block(embed_molecule(mol_from_smiles("[O-]C=O"), seed=7), charge=-1, name="f")
    return lig, sorted(lig.open_donors(), key=lambda s: s.atom_idx)


def proto_dimer():
    lig, d = _formate()
    m1 = metal_block("Cu", 2, cn=4, geometry="square_planar")
    first = join(lig, m1, d[0], m1.open_vacancies()[0], lone_pair=1)
    free = next(s for s in first.block.open_donors() if s.atom_idx == first.atom_map[d[1].atom_idx])
    m2 = metal_block("Cu", 2, cn=4, geometry="square_planar")
    return join(first.block, m2, free, m2.open_vacancies()[0], lone_pair=1).block


def _second_bridge(block):
    lig, d = _formate()
    m = block.graph.metals()
    vac = block.open_vacancies()
    pairs = [[p, q] for p in vac for q in vac if p.atom_idx == m[0] and q.atom_idx == m[1]]
    pair = min(pairs, key=lambda pq: (bridge_compatible(d, pq, partner="Cu",
                                                        donor_elements=["O", "O"]).strain,
                                      [(v.atom_idx, v.slot) for v in pq]))
    return join_bridge(lig, block, d, pair, lone_pairs=(0, 0)).block


def test_a_stored_block_loads_back_as_itself(reg):
    built = proto_dimer()
    st = store_block(reg, built)
    loaded = load_block(reg, st.structure_id, st.geometry_id)
    assert np.allclose(loaded.geometry, built.geometry, atol=1e-6)
    key = lambda b: sorted((s.role, s.slot, s.donor_type) for s in b.sites)  # noqa: E731
    assert key(loaded) == key(built)
    assert ({k: v.status for k, v in loaded.state.items()}
            == {k: v.status for k, v in built.state.items()})


def test_joining_onto_a_loaded_block_is_joining_onto_the_built_one(reg):
    built = proto_dimer()
    st = store_block(reg, built)
    a = _second_bridge(built)
    b = _second_bridge(load_block(reg, st.structure_id, st.geometry_id))
    assert np.allclose(np.asarray(a.geometry), np.asarray(b.geometry), atol=1e-5)


def test_a_rigid_copy_carries_its_frames(reg):
    built = proto_dimer()
    st = store_block(reg, built)
    rot = axis_rotation(np.array([0.3, 0.4, 0.866]), 0.7)
    moved = np.asarray(built.geometry) @ rot.T + np.array([1.0, -2.0, 0.5])
    gid = put_geometry(reg, st.structure_id, to_xyz(built.graph, moved), fidelity=Fidelity.RAW,
                       method=CONSTRUCT).id
    loaded = load_block(reg, st.structure_id, gid)
    vac = next(s for s in loaded.sites if s.is_vacancy)
    orig = next(s for s in built.sites if s.is_vacancy and s.slot == vac.slot
                and s.atom_idx == vac.atom_idx)
    assert np.allclose(vac.frame["axis"], rot @ np.asarray(orig.frame["axis"]), atol=1e-5)


def test_a_second_geometry_of_one_identity_loads_with_its_own_frames(reg):
    """Two poses of one node: the catalog keeps the first's frames, the state rows keep each's."""
    first = proto_dimer()
    store_block(reg, first)
    rot = axis_rotation(np.array([0.0, 0.0, 1.0]), 0.3)
    bent = np.asarray(first.geometry).copy()
    bent[: len(bent) // 2] = bent[: len(bent) // 2] @ rot.T     # not a rigid copy of the first
    from dataclasses import replace

    from mofsbu.assembly.join import _transform_frame
    moved = set(list(first.graph.nodes())[: len(bent) // 2])
    second = replace(first, geometry=bent, state=None, sites=tuple(
        replace(s, frame=_transform_frame(s.frame, rot, np.zeros(3), np.zeros(3)))
        if s.atom_idx in moved else s for s in first.sites))
    st = store_block(reg, second)
    loaded = load_block(reg, st.structure_id, st.geometry_id)
    want = {(s.atom_idx, s.slot): s.frame["axis"] for s in second.sites if s.frame}
    got = {(s.atom_idx, s.slot): s.frame["axis"] for s in loaded.sites if s.frame}
    assert all(np.allclose(got[k], want[k]) for k in want)


def test_a_geometry_its_frames_do_not_describe_is_refused(reg):
    built = proto_dimer()
    st = store_block(reg, built)
    bent = np.asarray(built.geometry).copy()
    bent[0] += [0.3, 0.0, 0.0]                      # not a rigid copy any more
    gid = put_geometry(reg, st.structure_id, to_xyz(built.graph, bent), fidelity=Fidelity.RAW,
                       method=CONSTRUCT).id
    with pytest.raises(NotBuiltYet, match="M6/S7"):
        load_block(reg, st.structure_id, gid)
