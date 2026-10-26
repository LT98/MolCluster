"""D31: a connectivity-changed relaxation is filed under the species it became.

The flagged geometry stays as history; the species it became is its own structure, joined to
the one it was built as by a `relaxation` edge that balances.  Single-piece metal-donor bond
changes and proton transfers are filed; a separation is held back pending the second-shell
probe; a skeleton change, a charged or open-shell leaver, or metals ending up apart are
refused with the reason.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import fixtures as fx
from mofsbu._types import Fidelity
from mofsbu.assembly.derive import DeriveRefused, derive_species, file_derived
from mofsbu.energy import backends as backends_mod
from mofsbu.energy import relax as relax_mod
from mofsbu.energy.backends import NullBackend
from mofsbu.energy.reference import check_balance
from mofsbu.graph import EdgeType, TypedGraph
from mofsbu.graph.canon import canonical_certificate
from mofsbu.registry import (
    MethodSpec,
    Registry,
    get_graph,
    get_structure,
    put_geometry,
    put_structure,
)
from mofsbu.runner import run
from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec

RAW = MethodSpec(code="construct", code_version="0", method="raw-construct")
OMOL = MethodSpec(code="mace", code_version="test", method="MACE-OMOL-0")
COV, DAT = EdgeType.COVALENT, EdgeType.DATIVE


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db") as r:
        r.migrate()
        yield r


def xyz_of(g, coords) -> str:
    return "\n".join([str(len(g)), g.name] + [
        f"{g.label(k).element} {x:.5f} {y:.5f} {z:.5f}"
        for k, (x, y, z) in enumerate(coords)]) + "\n"


def verdict_row(reg, gid):
    row = reg.conn.execute("SELECT connectivity_changed, connectivity_json FROM geometries "
                           "WHERE id=?", (gid,)).fetchone()
    return row["connectivity_changed"], json.loads(row["connectivity_json"] or "null")


# ── kappa1 -> kappa2: a donor closes onto the metal ──────────────────────────

def zn_formate_diaqua() -> TypedGraph:
    """[Zn(k1-O2CH)(H2O)2]+ : Zn 0, C 1, O 2 (bound), O 3 (free, C=O), H 4, waters 5-10."""
    g = TypedGraph(charge=1, multiplicity=1, name="[Zn(k1-HCOO)(H2O)2]+")
    zn = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    c, o1, o2, h = (g.add_atom(e) for e in ("C", "O", "O", "H"))
    g.add_bond(c, o1, COV)
    g.add_bond(c, o2, COV, order=2.0)
    g.add_bond(c, h, COV)
    g.add_bond(o1, zn, DAT)
    for _ in range(2):
        o, h1, h2 = g.add_atom("O"), g.add_atom("H"), g.add_atom("H")
        g.add_bond(o, h1, COV)
        g.add_bond(o, h2, COV)
        g.add_bond(o, zn, DAT)
    return g


def zn_coords(kappa2: bool) -> np.ndarray:
    o2 = [1.65, -1.10, 0.0] if kappa2 else [2.90, -1.10, 0.0]
    return np.array([[0, 0, 0], [2.45, 0, 0], [1.65, 1.10, 0], o2, [3.55, 0, 0],
                     [-2.05, 0, 0], [-2.65, 0.76, 0], [-2.65, -0.76, 0],
                     [0, 0, 2.05], [0, 0.76, 2.65], [0, -0.76, 2.65]], float)


def test_a_donor_closing_onto_the_metal_is_filed_as_the_chelate(reg):
    g = zn_formate_diaqua()
    sid = put_structure(reg, g).id
    raw = put_geometry(reg, sid, xyz_of(g, zn_coords(False)), fidelity=Fidelity.RAW,
                       method=RAW)
    relaxed = put_geometry(reg, sid, xyz_of(g, zn_coords(True)), fidelity=Fidelity.ML,
                           method=OMOL, energy=-42.0, converged=True, relaxed_from=raw.id)
    assert verdict_row(reg, relaxed.id)[0] == 1

    out = file_derived(reg, relaxed.id)
    assert len(out["filed"]) == 1
    new_sid, new_gid = out["filed"][0], out["geometries"][0]
    assert new_sid != sid
    new_graph = get_graph(reg, new_sid)
    assert len(new_graph.neighbors(next(iter(new_graph.metals())), DAT)) == 4   # k2 + 2 aqua
    same_blob = reg.conn.execute("SELECT a.coords_hash = b.coords_hash FROM geometries a, "
                                 "geometries b WHERE a.id=? AND b.id=?",
                                 (relaxed.id, new_gid)).fetchone()[0]
    assert same_blob
    row = reg.conn.execute("SELECT energy, fidelity, derived_from FROM geometries WHERE id=?",
                           (new_gid,)).fetchone()
    assert (row["energy"], row["fidelity"], row["derived_from"]) == (-42.0, int(Fidelity.ML),
                                                                     relaxed.id)
    assert verdict_row(reg, new_gid)[0] == 0                   # intact against its own graph
    assert verdict_row(reg, relaxed.id)[0] == 1                # the history stays flagged
    assert get_structure(reg, new_sid)["best_geometry_id"] == new_gid
    edge = reg.conn.execute("SELECT kind FROM reactions WHERE id=?", (out["edge"],)).fetchone()
    assert edge["kind"] == "relaxation" and check_balance(reg, out["edge"]).balanced
    assert file_derived(reg, relaxed.id)["skipped"] == "already filed"


@pytest.mark.parametrize("pending", [True, False], ids=["inside-a-transaction", "none-open"])
def test_an_edge_that_fails_to_balance_leaves_nothing_filed(reg, monkeypatch, pending):
    """All or none: the structure, geometry and edge are rolled back together, so a retry
    files the species instead of skipping it as already filed, and nothing is committed."""
    from mofsbu._types import MofsbuError
    from mofsbu.energy import reference

    g = zn_formate_diaqua()
    sid = put_structure(reg, g).id
    raw = put_geometry(reg, sid, xyz_of(g, zn_coords(False)), fidelity=Fidelity.RAW,
                       method=RAW)
    relaxed = put_geometry(reg, sid, xyz_of(g, zn_coords(True)), fidelity=Fidelity.ML,
                           method=OMOL, energy=-42.0, converged=True, relaxed_from=raw.id)
    if not pending:
        reg.conn.commit()

    def counts():
        return [reg.conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
                for t in ("structures", "geometries", "reactions")]

    before = counts()

    class Unbalanced:
        balanced = False

        def describe(self):
            return "unbalanced (test)"

    monkeypatch.setattr(reference, "check_balance", lambda reg_, rid: Unbalanced())
    with pytest.raises(MofsbuError, match="must conserve atoms"):
        file_derived(reg, relaxed.id)
    assert counts() == before
    assert reg.conn.in_transaction == pending, "file_derived committed or opened a transaction"

    monkeypatch.undo()
    assert len(file_derived(reg, relaxed.id)["filed"]) == 1


# ── a water leaves: the separation is held back ──────────────────────────────

AXES = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def hexaaqua_coords(pull: float = 0.0) -> np.ndarray:
    xyz = [np.zeros(3)]
    for k, a in enumerate(AXES):
        perp = np.cross(a, [0.3, 0.5, 0.8])
        perp /= np.linalg.norm(perp)
        o = (2.10 + (pull if k == 0 else 0.0)) * a
        xyz += [o, o + 0.60 * a + 0.76 * perp, o + 0.60 * a - 0.76 * perp]
    return np.array(xyz)


def test_a_water_leaving_is_held_back_and_nothing_is_filed(reg):
    g = fx.hexaaqua()
    sid = put_structure(reg, g).id
    raw = put_geometry(reg, sid, xyz_of(g, hexaaqua_coords()), fidelity=Fidelity.RAW,
                       method=RAW)
    relaxed = put_geometry(reg, sid, xyz_of(g, hexaaqua_coords(2.0)), fidelity=Fidelity.ML,
                           method=OMOL, energy=-80.0, converged=True, relaxed_from=raw.id)
    assert verdict_row(reg, relaxed.id)[0] == 1
    out = file_derived(reg, relaxed.id)
    assert "separation" in out["refused"]
    assert reg.conn.execute("SELECT COUNT(*) FROM geometries WHERE derived_from IS NOT NULL"
                            ).fetchone()[0] == 0
    assert reg.conn.execute("SELECT COUNT(*) FROM reactions WHERE kind='relaxation'"
                            ).fetchone()[0] == 0


# ── a proton crosses: formal charges follow the localisation rule ────────────

def test_a_proton_transfer_from_aqua_to_formate_gives_hydroxo_and_formic_acid():
    g = zn_formate_diaqua()
    o_water, h_water, o_free = 5, 6, 3
    verdict = {"changed": True,
               "broken": [{"i": o_water, "j": h_water, "sym_i": "O", "sym_j": "H",
                           "etype": "cov"}],
               "formed": [{"i": h_water, "j": o_free, "sym_i": "H", "sym_j": "O",
                           "etype": "cov"}]}
    [piece] = derive_species(g, verdict).pieces
    new = piece.graph
    assert new.label(o_water).formal_charge == -1              # localised: hydroxide
    assert new.label(o_free).formal_charge == 0                # delocalised: stays on graph
    assert new.net_charge() == 1

    want = TypedGraph(charge=1, multiplicity=1, name="[Zn(k1-HCOOH)(OH)(H2O)]+")
    zn = want.add_atom("Zn", oxidation_state=2, spin_class="ls")
    c, o1, o2, h = (want.add_atom(e) for e in ("C", "O", "O", "H"))
    for x in (o1, o2, h):
        want.add_bond(c, x, COV)
    want.add_bond(o1, zn, DAT)
    hacid = want.add_atom("H")
    want.add_bond(o2, hacid, COV)
    oh, hh = want.add_atom("O", formal_charge=-1), want.add_atom("H")
    want.add_bond(oh, hh, COV)
    want.add_bond(oh, zn, DAT)
    o, h1, h2 = want.add_atom("O"), want.add_atom("H"), want.add_atom("H")
    want.add_bond(o, h1, COV)
    want.add_bond(o, h2, COV)
    want.add_bond(o, zn, DAT)
    assert canonical_certificate(new) == canonical_certificate(want)


# ── refused, with the reason ─────────────────────────────────────────────────

def test_a_broken_skeleton_is_refused():
    g = zn_formate_diaqua()
    verdict = {"changed": True, "formed": [],
               "broken": [{"i": 1, "j": 4, "sym_i": "C", "sym_j": "H", "etype": "cov"}]}
    with pytest.raises(DeriveRefused, match="skeleton"):
        derive_species(g, verdict)


def test_a_charged_leaver_is_refused():
    g = TypedGraph(charge=1, multiplicity=1, name="[ZnCl(H2O)]+")
    zn = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    cl = g.add_atom("Cl", formal_charge=-1)
    g.add_bond(cl, zn, DAT)
    o, h1, h2 = g.add_atom("O"), g.add_atom("H"), g.add_atom("H")
    g.add_bond(o, h1, COV)
    g.add_bond(o, h2, COV)
    g.add_bond(o, zn, DAT)
    verdict = {"changed": True, "formed": [],
               "broken": [{"i": cl, "j": zn, "sym_i": "Cl", "sym_j": "Zn", "etype": "dat"}]}
    with pytest.raises(DeriveRefused, match="formal charge"):
        derive_species(g, verdict)


def test_metals_ending_up_apart_are_refused():
    g = fx.zn2_bridged_formates()
    zb = g.metals()[1]
    broken = [{"i": o, "j": zb, "sym_i": "O", "sym_j": "Zn", "etype": "dat"}
              for o in g.neighbors(zb, DAT)]
    with pytest.raises(DeriveRefused, match="metals"):
        derive_species(g, {"changed": True, "broken": broken, "formed": []})


# ── wired into the relax task ────────────────────────────────────────────────

class PullingBackend(NullBackend):
    """Relaxes by pulling the water nearest atom 0 (the metal) 2.5 Å further out."""

    def relax(self, symbols, positions, *, charge, multiplicity, solvent=None,
              fmax=0.20, steps=250):
        xyz = np.array(positions, dtype=float)
        sym = list(symbols)
        waters = []
        for o, s in enumerate(sym):
            hs = [h for h, t in enumerate(sym) if t == "H"
                  and np.linalg.norm(xyz[h] - xyz[o]) < 1.15]
            if s == "O" and len(hs) == 2 and sym[0] not in ("O", "H", "C"):
                waters.append((float(np.linalg.norm(xyz[o] - xyz[0])), o, hs))
        if waters:
            _d, o, hs = min(waters)
            out = (xyz[o] - xyz[0]) / np.linalg.norm(xyz[o] - xyz[0])
            for a in (o, *hs):
                xyz[a] = xyz[a] + 2.5 * out
        return super().relax(sym, xyz, charge=charge, multiplicity=multiplicity,
                             solvent=solvent, fmax=fmax, steps=steps)


def test_a_relax_task_records_a_separation_as_held_back(reg, monkeypatch):
    monkeypatch.setattr(relax_mod, "backend_for",
                        lambda fidelity, **kw: PullingBackend(fidelity=fidelity))
    spec = BuildSpec(
        molecules=(MoleculeSpec(name="AcOH", smiles="CC(=O)O", max_deprotonations=1),),
        metals=(MetalSpec(symbol="Zn", oxidation_state=2),),
        coordination=(4,), ligands_per_metal=(1,), binding=("mono",), co_ligand="O",
        run_mode="ml_go")
    result = run(reg, spec)
    assert result["status"] in {"done", "partial"}
    details = [json.loads(r["detail_json"]) for r in reg.conn.execute(
        "SELECT detail_json FROM tasks WHERE kind='relax' AND status='done'")]
    derived = [d["derived"] for d in details if "derived" in d]
    assert derived and any("separation" in d.get("refused", "") for d in derived), derived
    assert reg.conn.execute("SELECT COUNT(*) FROM reactions WHERE kind='relaxation'"
                            ).fetchone()[0] == 0
    assert backends_mod is not None
