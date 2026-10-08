"""Releasing a bound ligand (`assembly.release`) — the water-loss step `join` runs backwards."""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu.assembly.construct import ligand_block, metal_block
from mofsbu.assembly.join import join
from mofsbu.assembly.release import ReleaseRefused, release
from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.identity import identity


def _water(seed: int):
    return ligand_block(embed_molecule(mol_from_smiles("O"), seed=seed), charge=0, name="aqua")


def tetraaqua():
    """[Cu(H2O)4]2+ by four joins, and the atom index of each water's oxygen."""
    block = metal_block("Cu", 2, cn=4, geometry="square_planar")
    oxygens = []
    for k in range(4):
        lig = _water(11 + k)
        res = join(lig, block, lig.open_donors()[0], block.open_vacancies()[0])
        oxygens = [res.partner_atom_map[o] for o in oxygens] + [res.atom_map[0]]
        block = res.block
    return block, oxygens


def _pos(block, atom):
    order = {n: k for k, n in enumerate(block.graph.nodes())}
    return np.asarray(block.geometry)[order[atom]]


def test_a_released_water_opens_its_vertex_pointing_where_it_was():
    block, oxygens = tetraaqua()
    metal = block.graph.metals()[0]
    want = _pos(block, oxygens[2]) - _pos(block, metal)
    out = release(block, [oxygens[2]])
    assert out.block.graph.element_counts() == {"Cu": 1, "O": 3, "H": 6}
    assert out.ligand.graph.element_counts() == {"O": 1, "H": 2}
    assert (out.block.graph.net_charge(), out.ligand.graph.net_charge()) == (2, 0)
    (vac,) = out.block.open_vacancies()
    assert np.allclose(vac.frame["axis"], want / np.linalg.norm(want))
    # nothing that stayed has moved
    for old, new in out.atom_map.items():
        assert np.allclose(_pos(block, old), _pos(out.block, new))


def test_release_then_rejoin_is_the_same_node():
    block, oxygens = tetraaqua()
    out = release(block, [oxygens[0]])
    lig = _water(99)
    back = join(lig, out.block, lig.open_donors()[0], out.block.open_vacancies()[0]).block
    assert identity(back.graph)["l1"] == identity(block.graph)["l1"]


def test_a_bridge_is_refused_by_name():
    lig = ligand_block(embed_molecule(mol_from_smiles("[O-]C=O"), seed=7), charge=-1, name="f")
    d = sorted(lig.open_donors(), key=lambda s: s.atom_idx)
    m1 = metal_block("Cu", 2, cn=4, geometry="square_planar")
    first = join(lig, m1, d[0], m1.open_vacancies()[0], lone_pair=1)
    free = next(s for s in first.block.open_donors() if s.atom_idx == first.atom_map[d[1].atom_idx])
    m2 = metal_block("Cu", 2, cn=4, geometry="square_planar")
    dimer = join(first.block, m2, free, m2.open_vacancies()[0], lone_pair=1).block
    with pytest.raises(ReleaseRefused, match="bridge"):
        release(dimer, [first.atom_map[d[0].atom_idx]])


def test_an_unbound_fragment_and_a_metal_are_refused():
    block, _ = tetraaqua()
    with pytest.raises(ReleaseRefused, match="metal is not released"):
        release(block, [block.graph.metals()[0]])
    with pytest.raises(ReleaseRefused, match="not bound to a metal"):
        release(_water(1), [0])
