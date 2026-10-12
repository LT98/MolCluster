"""Gate E1: a geometry that broke or formed a bond is flagged and never selected for energy.

`docs/WORKPLAN_energy.md` §3 stage 2.  The check is `geometry.connectivity`; the registry
records its verdict on every geometry above RAW, and energy selection reads it.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import fixtures as fx
from mofsbu._types import Fidelity, MofsbuError, ReferenceSchemeError
from mofsbu.energy.reference import Term, _energy_row, energy_for_terms
from mofsbu.geometry.connectivity import BROKEN, FORMED, connectivity_check
from mofsbu.graph import EdgeType, TypedGraph
from mofsbu.registry import (MethodSpec, Registry, backfill_connectivity, best_geometry,
                             get_structure, put_geometry, put_structure)
from mofsbu.registry.verify import verify

RAW = MethodSpec(code="construct", code_version="0", method="raw-construct")
OMOL = MethodSpec(code="mace", code_version="0.3.14/extra_large", method="MACE-OMOL-0")

AXES = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def hexaaqua_coords(d_mo: float = 2.10) -> np.ndarray:
    """[M(H2O)6] in graph order: M, then (O, H, H) per water, each water on one axis."""
    xyz = [np.zeros(3)]
    for a in AXES:
        perp = np.cross(a, [0.3, 0.5, 0.8])
        perp /= np.linalg.norm(perp)
        o = d_mo * a
        xyz += [o, o + 0.60 * a + 0.76 * perp, o + 0.60 * a - 0.76 * perp]
    return np.array(xyz)


def to_xyz(graph, coords) -> str:
    lines = [str(len(graph)), graph.name]
    lines += [f"{graph.label(k).element} {x:.5f} {y:.5f} {z:.5f}"
              for k, (x, y, z) in enumerate(coords)]
    return "\n".join(lines) + "\n"


def symbols(graph) -> list[str]:
    return [graph.label(k).element for k in graph.nodes()]


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


# ── the check ────────────────────────────────────────────────────────────────

def test_an_intact_geometry_reports_no_change():
    g = fx.hexaaqua()
    report = connectivity_check(symbols(g), hexaaqua_coords(), g)
    assert not report.changed, report.describe()
    assert report.to_dict()["thresholds"] == {"broken": BROKEN, "formed": FORMED}


def test_a_water_drifting_off_the_metal_is_a_broken_dative_bond():
    g = fx.hexaaqua()
    xyz = hexaaqua_coords()
    xyz[1:4] += 2.0 * AXES[0]                     # first water out to ~4.1 A
    report = connectivity_check(symbols(g), xyz, g)
    assert report.changed and not report.formed
    [bond] = report.broken
    assert (bond.i, bond.j, bond.etype) == (0, 1, EdgeType.DATIVE)
    assert bond.distance > bond.limit == pytest.approx(BROKEN * bond.ideal)


def test_a_proton_crossing_to_another_oxygen_breaks_one_bond_and_forms_another():
    g = fx.hexaaqua()
    xyz = hexaaqua_coords()
    o_a, h_a, o_b = 1, 2, 7                       # H of the +x water onto the +y water's O
    xyz[h_a] = xyz[o_b] + 0.98 * (xyz[o_a] - xyz[o_b]) / np.linalg.norm(xyz[o_a] - xyz[o_b])
    report = connectivity_check(symbols(g), xyz, g)
    assert [(c.i, c.j) for c in report.broken] == [(o_a, h_a)]
    assert [(c.i, c.j, c.etype) for c in report.formed] == [(h_a, o_b, EdgeType.COVALENT)]


def test_a_donor_swinging_onto_the_metal_forms_a_dative_bond():
    """kappa1-formate in the graph; the free O ends up at Zn-O distance: kappa2 formed."""
    g = TypedGraph(charge=1, multiplicity=1, name="[Zn(k1-O2CH)]+")
    zn = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    c, o1, o2, h = (g.add_atom(e) for e in ("C", "O", "O", "H"))
    for x in (o1, o2, h):
        g.add_bond(c, x, EdgeType.COVALENT)
    g.add_bond(o1, zn, EdgeType.DATIVE)
    xyz = np.array([[0.0, 0.0, 0.0], [2.45, 0.0, 0.0], [1.80, 1.10, 0.0],
                    [1.80, -1.10, 0.0], [3.55, 0.0, 0.0]])
    report = connectivity_check(symbols(g), xyz, g)
    assert not report.broken
    assert [(c.i, c.j, c.etype) for c in report.formed] == [(zn, o2, EdgeType.DATIVE)]


def test_two_metals_coming_close_form_nothing_but_a_declared_edge_can_break():
    g = fx.zn2_bridged_formates()                 # two Zn, no declared M-M edge (D25)
    assert not any(et is EdgeType.METAL_METAL for _u, _v, et in g.edges())
    close = connectivity_check(symbols(g), _spread(g, g.metals(), 2.4), g)
    assert not any(c.etype is EdgeType.METAL_METAL for c in close.formed)

    pw = fx.cu_paddlewheel()                      # Cu-Cu declared
    far = connectivity_check(symbols(pw), _spread(pw, pw.metals(), 6.0), pw)
    assert [(c.i, c.j) for c in far.broken if c.etype is EdgeType.METAL_METAL] == \
        [tuple(pw.metals())]


def _spread(g, metals, d_mm) -> np.ndarray:
    """Every atom far from every other, except the two metals at `d_mm`."""
    xyz = np.array([[10.0 * k, 0.0, 0.0] for k in range(len(g))])
    xyz[metals[1]] = xyz[metals[0]] + [0.0, d_mm, 0.0]
    return xyz


def test_atoms_in_another_order_are_refused_not_judged():
    g = fx.hexaaqua()
    sym = symbols(g)
    sym[0], sym[1] = sym[1], sym[0]
    with pytest.raises(MofsbuError, match="one atom order"):
        connectivity_check(sym, hexaaqua_coords(), g)


# ── the registry records it, and energy selection reads it ───────────────────

def _stored_hexaaqua(reg):
    g = fx.hexaaqua()
    sid = put_structure(reg, g).id
    raw = put_geometry(reg, sid, to_xyz(g, hexaaqua_coords()), fidelity=Fidelity.RAW,
                       method=RAW)
    return g, sid, raw.id


def _torn(g):
    xyz = hexaaqua_coords(2.05)
    xyz[1:4] += 2.0 * AXES[0]
    return to_xyz(g, xyz)


def _verdict(reg, gid):
    row = reg.conn.execute("SELECT connectivity_changed, connectivity_json FROM geometries "
                           "WHERE id=?", (gid,)).fetchone()
    return row["connectivity_changed"], json.loads(row["connectivity_json"] or "null")


def test_a_raw_construct_is_not_checked_and_a_relaxed_one_is(reg):
    g, sid, raw = _stored_hexaaqua(reg)
    assert _verdict(reg, raw) == (None, None)
    ok = put_geometry(reg, sid, to_xyz(g, hexaaqua_coords(2.05)), fidelity=Fidelity.ML,
                      method=OMOL, energy=-10.0, converged=True, relaxed_from=raw)
    changed, detail = _verdict(reg, ok.id)
    assert changed == 0 and detail["version"] == "conn1" and not detail["changed"]


def test_a_changed_geometry_is_flagged_and_never_best(reg):
    g, sid, raw = _stored_hexaaqua(reg)
    torn = put_geometry(reg, sid, _torn(g), fidelity=Fidelity.ML, method=OMOL,
                        energy=-11.0, converged=True, relaxed_from=raw)
    changed, detail = _verdict(reg, torn.id)
    assert changed == 1 and detail["broken"][0]["etype"] == "dat"
    assert get_structure(reg, sid)["best_geometry_id"] == raw
    assert best_geometry(reg, sid, Fidelity.ML) is None
    assert not [p for p in verify(reg) if p.check == "best_geometry"]


def test_an_equation_refuses_a_species_whose_only_energy_changed_connectivity(reg):
    g, sid, raw = _stored_hexaaqua(reg)
    put_geometry(reg, sid, _torn(g), fidelity=Fidelity.ML, method=OMOL, energy=-11.0,
                 converged=True, relaxed_from=raw)
    with pytest.raises(ReferenceSchemeError, match="connectivity_changed"):
        energy_for_terms(reg, (Term(sid, 1, "reagent", "reagent"),
                               Term(sid, 1, "product", "product")), fidelity=Fidelity.ML)


def test_an_intact_sibling_is_used_when_the_lower_energy_one_changed(reg):
    g, sid, raw = _stored_hexaaqua(reg)
    put_geometry(reg, sid, _torn(g), fidelity=Fidelity.ML, method=OMOL, energy=-11.0,
                 converged=True, relaxed_from=raw)
    ok = put_geometry(reg, sid, to_xyz(g, hexaaqua_coords(2.05)), fidelity=Fidelity.ML,
                      method=OMOL, energy=-10.0, converged=True, relaxed_from=raw)
    assert get_structure(reg, sid)["best_geometry_id"] == ok.id
    assert _energy_row(reg, sid, fidelity=Fidelity.ML, solvent=None)["id"] == ok.id


def test_backfill_flags_rows_written_before_e1_and_moves_the_pointer(reg):
    g, sid, raw = _stored_hexaaqua(reg)
    torn = put_geometry(reg, sid, _torn(g), fidelity=Fidelity.ML, method=OMOL,
                        energy=-11.0, converged=True, relaxed_from=raw)
    reg.conn.execute("UPDATE geometries SET connectivity_changed=NULL, connectivity_json=NULL")
    reg.conn.execute("UPDATE structures SET best_geometry_id=?, best_fidelity=? WHERE id=?",
                     (torn.id, int(Fidelity.ML), sid))
    expect = {"checked": 1, "changed": 1, "unchecked": 0, "refused": {}}
    assert backfill_connectivity(reg, dry_run=True) == expect
    assert _verdict(reg, torn.id) == (None, None)
    assert backfill_connectivity(reg) == expect
    assert _verdict(reg, torn.id)[0] == 1
    assert get_structure(reg, sid)["best_geometry_id"] == raw
    assert backfill_connectivity(reg)["checked"] == 0


def test_a_relaxation_of_a_build_in_another_numbering_is_unchecked_not_excluded(reg):
    """Until B28 is fixed, a verdict read under the wrong numbering would exclude a good energy."""
    g, sid, _raw = _stored_hexaaqua(reg)
    xyz = hexaaqua_coords()
    o1, o2 = 1, 4                                  # two water oxygens: same elements, other numbering
    xyz[[o1, o2]] = xyz[[o2, o1]]
    raw2 = put_geometry(reg, sid, to_xyz(g, xyz), fidelity=Fidelity.RAW, method=RAW)
    relaxed = put_geometry(reg, sid, to_xyz(g, xyz + 0.01), fidelity=Fidelity.ML, method=OMOL,
                           energy=-10.0, converged=True, relaxed_from=raw2.id)
    changed, detail = _verdict(reg, relaxed.id)
    assert changed is None and "B28" in detail["unchecked"]
    assert _energy_row(reg, sid, fidelity=Fidelity.ML, solvent=None)["id"] == relaxed.id
