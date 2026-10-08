"""A route join states the face it means (`torsion: {bring, toward}`), because a well index
does not name a face (B27): cis vertices' torsion references point opposite ways."""
from __future__ import annotations

import numpy as np

from mofsbu.assembly.route_steps import execute_route
from mofsbu.spec import BuildSpec, MoleculeSpec, RouteSpec


def _two_formates_on_one_cu(torsion=None):
    second = {"op": "join", "donor": "f2.d0", "onto": "cu1", "slot": 1, "lone_pair": 1}
    if torsion is not None:
        second["torsion"] = torsion
    route = RouteSpec(name="r", steps=(
        {"op": "metal", "as": "cu1", "symbol": "Cu", "oxidation_state": 2, "cn": 4,
         "geometry": "square_planar"},
        {"op": "ligand", "as": "f1", "molecule": "formate"},
        {"op": "join", "donor": "f1.d0", "onto": "cu1", "slot": 0, "lone_pair": 1},
        {"op": "ligand", "as": "f2", "molecule": "formate"},
        second,
        {"op": "store", "as": "out", "block": "cu1"},
    ))
    spec = BuildSpec(molecules=(MoleculeSpec(name="formate", smiles="[O-]C=O"),))
    got = {}

    def keep(ws, step, s, bid):
        got["block"], got["rows"] = ws.blocks[bid], ws.rows(bid)
        return type("Stored", (), {"structure_id": 0, "geometry_id": 0, "kind": "built",
                                   "name": s["as"]})()

    execute_route(spec, route, keep)
    xyz = np.asarray(got["block"].geometry)
    free1, free2 = (xyz[got["rows"][k][0]] for k in ("f1.d1", "f2.d1"))
    return float(np.linalg.norm(free1 - free2))


def test_well_zero_on_cis_vertices_puts_two_formates_on_opposite_faces():
    """B27, measured: the default well is not one face across vertices."""
    assert _two_formates_on_one_cu() > 4.0


def test_torsion_toward_puts_the_second_formate_on_the_first_ones_face():
    same = _two_formates_on_one_cu({"bring": "f2.d1", "toward": "f1.d1"})
    assert same < 3.5 < _two_formates_on_one_cu()


def test_a_torsion_index_is_still_accepted():
    assert _two_formates_on_one_cu(0) == _two_formates_on_one_cu()
