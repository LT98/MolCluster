"""M4's exit gate: the three properties the site model has to have to be trusted.

These were the gate items with no tests behind them for three milestones.  Each one
guards a claim the rest of the project reads as settled:

* **canonical-index stability** — D5 stores sites by canonical index so they survive
  recall and reordering.  If that fails, every `site_state` row points at the wrong atom
  after a permutation and nothing downstream would notice.
* **perceive-once** — §6.3 splits the two tiers precisely so perception is not repeated.
  `refresh_state` re-perceiving would let two geometries of one structure disagree about
  which atoms are donors.
* **frame reproducibility** — the whole point of promoting `_donor_placement_frame` out of
  the placer (D13) is that the frame is stored rather than re-found by stochastic search.
  A frame that does not survive the database round-trip has not been promoted anywhere.

A fourth property joined them once perception was made resonance-invariant: **the catalog
is a function of the identity it hangs off**, and `catalog_drift` returning empty for the
build routes in `build_routes.py` is the gate on it.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from build_routes import ROUTE_CASES, build, routes
from mofsbu._types import Fidelity
from mofsbu.descriptors.ease import PKA_CENTRE, activation_ease
from mofsbu.geometry.embed import embed_molecule, to_xyz
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles
from mofsbu.registry import (
    BlobStore, MethodSpec, Registry, RegistryError, canonical_map, catalog_drift,
    get_site_state, get_sites, put_geometry, put_site_state, put_sites, put_structure,
)
from mofsbu.sites import perception
from mofsbu.sites.frames import site_frame
from mofsbu.sites.inherit import inherit_sites, merge_inherited
from mofsbu.sites.model import chelate_pockets, perceive
from mofsbu.sites.state import SiteStatus, refresh_state

FF = MethodSpec(code="rdkit", code_version="2026.03", method="ETKDGv3+MMFF")

SALICYLIC = "OC(=O)c1ccccc1O"


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as r:
        r.migrate("test")
        yield r


def _store(reg, smiles: str, name: str = "L"):
    mol = embed_molecule(mol_from_smiles(smiles), seed=7)
    graph = from_rdkit(mol, charge=None, multiplicity=1, name=name)
    put = put_structure(reg, graph, tags=[name])
    geom = put_geometry(reg, put.id, to_xyz(mol, name), fidelity=Fidelity.FF, method=FF)
    sites = perceive(mol)
    put_sites(reg, put.id, sites)
    return mol, graph, put.id, geom.id, sites


# ── gate 1: canonical-index stability ────────────────────────────────────────

def test_shuffling_the_input_order_gives_the_same_site_catalog(reg):
    """Same molecule, atoms written in a different order -> identical catalog rows.

    Not "the same sites in a different order" — identical, because the canonical index
    is what is stored and it is a property of the graph, not of the input ordering.
    """
    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    sites = perceive(mol)
    graph = from_rdkit(mol, charge=None, multiplicity=1, name="sal")

    put_a = put_structure(reg, graph, tags=["a"])
    put_sites(reg, put_a.id, sites)
    rows_a = [(r["canonical_idx"], r["donor_type"], r["live_dof"], r["binding_modes"])
              for r in get_sites(reg, put_a.id)]

    # A permutation of the SAME graph is the same identity (D2), so it lands on the same
    # structure row; the catalog must already agree with what the shuffled order perceives.
    shuffled = graph.permuted(5)
    put_b = put_structure(reg, shuffled, tags=["b"])
    assert put_b.id == put_a.id, "a permutation is not a new identity"

    cmap_a = canonical_map(reg, put_a.id)
    assert sorted(cmap_a[s.atom_idx] for s in sites) == sorted(r[0] for r in rows_a)
    assert rows_a == sorted(rows_a), "catalog is returned in canonical order"
    assert len({r[0] for r in rows_a}) == len(rows_a), "one row per canonical index"


def test_a_reordered_molecule_perceives_the_same_donor_types(reg):
    """Perception itself must not depend on input atom order."""
    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    base = sorted(s.donor_type for s in perceive(mol))
    from rdkit import Chem

    order = list(range(mol.GetNumAtoms()))[::-1]
    reordered = Chem.RenumberAtoms(mol, order)
    assert sorted(s.donor_type for s in perceive(reordered)) == base


# ── gate 2: perceive once ────────────────────────────────────────────────────

def test_refresh_state_never_perceives(monkeypatch):
    """The counter §6.3's two-tier split is worth.

    `refresh_state` is handed sites and must use them.  If it ever calls perception, two
    geometries of one structure can disagree about which atoms are donors, and the
    canonical indices in `site_state` stop meaning the same thing across the ladder.
    """
    calls = []
    real = perception.find_donor_sites
    counted = lambda *a, **k: (calls.append(1), real(*a, **k))[1]      # noqa: E731
    # Patched in BOTH namespaces on purpose: `sites.model` binds the name at import, so
    # patching only the defining module would leave the call this test is counting
    # completely unobserved — and the test would pass by never seeing anything.
    monkeypatch.setattr(perception, "find_donor_sites", counted)
    monkeypatch.setattr("mofsbu.sites.model.find_donor_sites", counted)

    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    sites = perceive(mol)                       # the one legitimate perception
    assert len(calls) == 1

    conf = mol.GetConformer()
    coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
               conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    refresh_state(sites, coords, symbols=[a.GetSymbol() for a in mol.GetAtoms()])
    refresh_state(sites, coords, symbols=[a.GetSymbol() for a in mol.GetAtoms()])
    assert len(calls) == 1, "refresh_state re-perceived"


def test_the_catalog_is_written_once_per_structure(reg):
    """A second build of one identity keeps the first catalog — and its state rows.

    This is the bug that put state on the wrong geometry: `put_sites` used to DELETE the
    catalog first, and the cascade took `site_state` with it.  Under D2, re-deriving an
    identity the registry already has is the EXPECTED outcome for most of an enumeration,
    so that fired constantly.
    """
    mol, graph, sid, gid_a, sites = _store(reg, SALICYLIC, "sal")
    conf = mol.GetConformer()
    coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
               conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    states = refresh_state(sites, coords, graph=graph)
    put_site_state(reg, sid, gid_a, states, fidelity=Fidelity.FF)
    assert len(get_site_state(reg, sid, gid_a)) == len(sites)

    # A second geometry of the same structure, and a second call to put_sites.
    gid_b = put_geometry(reg, sid, to_xyz(mol, "sal2"), fidelity=Fidelity.RAW,
                         method=FF).id
    put_sites(reg, sid, sites)
    assert len(get_site_state(reg, sid, gid_a)) == len(sites), (
        "writing the catalog again destroyed the first geometry's state")
    put_site_state(reg, sid, gid_b, states, fidelity=Fidelity.RAW)
    assert len(get_site_state(reg, sid, gid_b)) == len(sites)


def _store_route(reg, name: str, mol, route: tuple[int, ...], donor_type: str):
    """One build route, all the way through the registry: identity, geometry, sites."""
    complex_mol = build(mol, route, donor_type)
    graph = from_rdkit(complex_mol, charge=None, multiplicity=1, name=name)
    put = put_structure(reg, graph, tags=[name])
    put_geometry(reg, put.id, to_xyz(complex_mol, name), fidelity=Fidelity.FF, method=FF)
    return put.id, perceive(complex_mol)


@pytest.mark.parametrize("name,smiles,donor_type,denticity", ROUTE_CASES,
                         ids=[c[0] for c in ROUTE_CASES])
def test_no_build_route_drifts_from_the_stored_catalog(
        reg, name, smiles, donor_type, denticity):
    """The gate on `site_catalog` being a function of the identity it hangs off.

    `catalog_drift` exists because it was not.  D15 excludes bond order from the L1 hash
    so that C=O and C-O(-) hash identically, perception read bond order, and a monodentate
    acetate therefore recorded one carboxylate donor or two depending on which oxygen the
    first build route bound through.  `put_sites` keeps the first catalog, so a drift here
    is not cosmetic: it is the second route being told it was wrong about its own donors.

    Empty is the whole assertion.  If perception ever becomes bond-order-sensitive again
    this is where it shows up, and it shows up as a list naming the atoms that moved.
    """
    mol, all_routes = routes(smiles, donor_type, denticity)
    assert len(all_routes) > 1, f"{name}: needs at least two routes to compare"

    structure_ids, drifts = [], []
    for route in all_routes:
        sid, sites = _store_route(reg, name, mol, route, donor_type)
        drifts.append((route, catalog_drift(reg, sid, sites)))
        put_sites(reg, sid, sites)
        structure_ids.append(sid)

    assert len(set(structure_ids)) == 1, (
        f"{name}: the routes {all_routes} landed on structures {structure_ids} — they are "
        f"supposed to be ONE identity, so this test is not testing what it claims to. "
        f"Fix the fixture before reading anything into the drift below.")
    assert all(not drift for _, drift in drifts), (
        f"{name}: routes disagree with the stored catalog: "
        + "; ".join(f"{route} -> {drift}" for route, drift in drifts if drift))


# ── gate 3: a frame survives the database ────────────────────────────────────

def test_a_stored_frame_reproduces_the_re_derived_frame(reg):
    """1e-6, which is the tolerance the plan names.

    If this drifts, a stored geometry is not regenerable from its site model and D13's
    promotion of the frame out of the stochastic search has bought nothing.
    """
    mol, _graph, sid, _gid, sites = _store(reg, SALICYLIC, "sal")
    conf = mol.GetConformer()
    by_canonical = {r["canonical_idx"]: r for r in get_sites(reg, sid)}
    cmap = canonical_map(reg, sid)

    assert sites, "nothing to check"
    for site in sites:
        row = by_canonical[cmap[site.atom_idx]]
        stored = json.loads(row["frame_json"])
        fresh = site_frame(mol, site.atom_idx, site.donor_type, conf).to_dict()
        for key in ("origin", "axis", "ref"):
            assert np.allclose(np.array(stored[key]), np.array(fresh[key]), atol=1e-6), (
                f"{site.donor_type} frame {key} did not survive the round trip")


def test_frames_are_orthonormal_after_the_round_trip(reg):
    _mol, _graph, sid, _gid, _sites = _store(reg, SALICYLIC, "sal")
    for row in get_sites(reg, sid):
        frame = json.loads(row["frame_json"])
        axis, ref = np.array(frame["axis"]), np.array(frame["ref"])
        assert np.isclose(np.linalg.norm(axis), 1.0, atol=1e-6)
        assert np.isclose(np.linalg.norm(ref), 1.0, atol=1e-6)
        assert abs(float(axis @ ref)) < 1e-6


# ── gate 4: known-molecule coverage, as the plan words it ────────────────────

@pytest.mark.parametrize("name,smiles,expect", [
    ("bipy", "c1ccc(-c2ccccn2)nc1", {"pyridyl_N": 2}),
    ("aqua", "O", {"aqua_O": 1}),
    ("methylamine", "CN", {"amine_N": 1}),
    ("EDTA", "OC(=O)CN(CC(O)=O)CCN(CC(O)=O)CC(O)=O",
     {"carboxylate_O": 8, "amine_N": 2}),
])
def test_known_molecule_donor_coverage(name, smiles, expect):
    from collections import Counter

    mol = embed_molecule(mol_from_smiles(smiles), seed=7)
    assert dict(Counter(s.donor_type for s in perceive(mol))) == expect, name


def test_btc_has_three_carboxylate_pockets_each_offering_mono_chelate_bridge():
    """The plan words this per GROUP; perception is per ATOM, and both are true.

    BTC has six carboxylate oxygens and three carboxylate GROUPS.  The group-level view
    the exit gate describes is the pocket layer's, which is where it belongs: whether two
    donors can chelate is a property of the PAIR, not of either donor alone.
    """
    mol = embed_molecule(mol_from_smiles("OC(=O)c1cc(C(O)=O)cc(C(O)=O)c1"), seed=7)
    sites = perceive(mol)
    assert sum(1 for s in sites if s.donor_type == "carboxylate_O") == 6

    pockets = chelate_pockets(mol, sites)
    assert len(pockets) == 3, "one pocket per carboxylate group"
    for pocket in pockets:
        assert pocket.ring_size == 4 and pocket.rigid
    modes = {m for s in sites for m in s.binding_modes}
    assert {"mono", "chelate", "mu2"} <= modes


def test_torsion_live_and_free_split_the_way_the_gate_says():
    live = perceive(embed_molecule(mol_from_smiles("c1ccc(-c2ccccn2)nc1"), seed=7))
    assert {s.live_dof for s in live} == {"live"}, "pyridyl N branches torsion"
    for smiles in ("O", "CN"):
        free = perceive(embed_molecule(mol_from_smiles(smiles), seed=7))
        assert {s.live_dof for s in free} == {"free"}, f"{smiles} must not branch"


# ── the state layer itself ───────────────────────────────────────────────────

def test_a_donor_dative_bonded_to_a_metal_reads_as_occupied():
    """The status rule itself, on a graph built to exercise it.

    Deliberately NOT routed through `perceive` on an assembled complex.  The reason used
    to be the resonance seam, which is now closed for delocalised groups; what is left is
    narrower and still real — the per-atom valence rules count a metal as an ordinary
    heavy neighbour, so a coordinated aqua oxygen is not perceived as a donor at all.
    Either way a test that went that way would be testing perception's treatment of bound
    atoms rather than the status rule, and would pass or fail for the wrong reason.
    """
    from mofsbu.graph._types import EdgeType, TypedGraph
    from mofsbu.sites.model import Site

    g = TypedGraph(charge=2, multiplicity=1, name="probe")
    metal = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    bound = g.add_atom("O")
    free = g.add_atom("O")
    g.add_bond(bound, metal, EdgeType.DATIVE)

    def site(idx: int) -> Site:
        return Site(atom_idx=idx, donor_type="aqua_O", labile=False, charge_after=0,
                    live_dof="free", binding_modes=("mono",))

    coords = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [-8.0, 0.0, 0.0]]
    states = {s.atom_idx: s for s in
              refresh_state([site(bound), site(free)], coords, graph=g)}
    assert states[bound].status is SiteStatus.OCCUPIED
    assert states[free].status is SiteStatus.OPEN
    assert not states[bound].is_open


def test_occupied_and_blocked_are_different_facts():
    """D11 makes a flip of the open-site flag the primary conformer trigger, so a buried
    site and a bonded one must not collapse into one 'not open'."""
    assert SiteStatus.OCCUPIED is not SiteStatus.BLOCKED
    assert {s.value for s in SiteStatus} == {"open", "occupied", "blocked"}


def test_buried_volume_stays_in_range_on_a_real_complex():
    from mofsbu.geometry.placer import LigandPlacement, place_mononuclear, to_rdkit
    from mofsbu.sites.frames import BindingMode

    lig = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    lig_sites = perceive(lig)
    pocket = chelate_pockets(lig, lig_sites)[0]
    by_idx = {s.atom_idx: s for s in lig_sites}
    placements = [LigandPlacement(
        mol=lig, donor_idxs=tuple(pocket.donors),
        donor_types=tuple(by_idx[i].donor_type for i in pocket.donors),
        mode=BindingMode.CHELATE, name="sal")]
    water = embed_molecule(mol_from_smiles("O"), seed=3)
    placements += [LigandPlacement(mol=water, donor_idxs=(0,), donor_types=("aqua_O",),
                                   name="O") for _ in range(2)]
    result = place_mononuclear("Cu", placements, geometry="tetrahedral")
    complex_mol = to_rdkit("Cu", placements, result)
    graph = from_rdkit(complex_mol, charge=2, multiplicity=2)

    sites = perceive(complex_mol)
    conf = complex_mol.GetConformer()
    coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
               conf.GetAtomPosition(i).z] for i in range(complex_mol.GetNumAtoms())]
    states = refresh_state(sites, coords, graph=graph,
                           symbols=[a.GetSymbol() for a in complex_mol.GetAtoms()])
    assert states, "the complex perceived no sites at all"
    assert all(s.buried_vol is None or 0.0 <= s.buried_vol <= 1.0 for s in states)


def _provisional_fraction(smiles: str) -> tuple[int, int]:
    from mofsbu.sites.model import shifting_pocket_donors

    mol = embed_molecule(mol_from_smiles(smiles), seed=7)
    sites = perceive(mol)
    conf = mol.GetConformer()
    coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
               conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    states = refresh_state(sites, coords,
                           pocket_donors=shifting_pocket_donors(mol, sites),
                           symbols=[a.GetSymbol() for a in mol.GetAtoms()])
    return (sum(1 for s in states if s.ease is not None and s.ease.provisional),
            len(states))


@pytest.mark.parametrize("name,smiles", [
    ("salicylic acid", SALICYLIC),          # pKa1 2.97 vs benzoic 4.20 — real H-bond shift
    ("catechol", "Oc1ccccc1O"),             # pKa1 9.25 vs phenol 9.99
])
def test_a_neighbour_shifted_donor_asks_to_be_promoted(name, smiles):
    """§6.6's named case, and the two molecules where it is measurably true."""
    flagged, total = _provisional_fraction(smiles)
    assert flagged == total > 0, f"{name}: every donor here has a shifted pKa"


@pytest.mark.parametrize("name,smiles", [
    ("benzoic acid", "OC(=O)c1ccccc1"),     # the table's 4.76 IS this molecule's value
    ("BTC", "OC(=O)c1cc(C(O)=O)cc(C(O)=O)c1"),
    ("4-hydroxybenzoic", "OC(=O)c1ccc(O)cc1"),   # para: too far to interact
    ("bipy", "c1ccc(-c2ccccn2)nc1"),        # neutral donors: no pKa to be wrong about
    ("phenol", "Oc1ccccc1"),
])
def test_an_unshifted_donor_does_not_waste_qm_on_itself(name, smiles):
    """The half that makes the flag worth having.

    A carboxylate closes a 4-ring through its own two oxygens, so a naive "is in a
    pocket" rule flags every carboxylate in existence — and a work list containing
    everything is the same as no work list.  `shifting_pocket_donors` requires an
    INTER-GROUP ring (>= 5), which is the structural signature of the neighbour effect
    D18 actually means.
    """
    flagged, total = _provisional_fraction(smiles)
    assert total > 0 and flagged == 0, f"{name}: nothing here needs a computed pKa"


def test_ease_absent_components_are_absent_not_zero():
    """D18's central rule.  A geometry-free record must not claim a steric number."""
    from mofsbu.sites.model import Site

    site = Site(atom_idx=0, donor_type="carboxylate_O", labile=True, charge_after=-1,
                live_dof="live", binding_modes=("mono",))
    bare = activation_ease(site)
    assert "steric" not in bare.components
    assert bare.coverage < 1.0 and bare.fidelity is Fidelity.HEURISTIC

    class _State:
        buried_vol = 0.25
        in_pocket = False

    withgeom = activation_ease(site, state=_State())
    assert "steric" in withgeom.components
    assert withgeom.coverage > bare.coverage
    assert withgeom.confidence > bare.confidence, "more of the model ran"


def test_a_neutral_donor_is_easy_not_missing():
    """`pka` empty means "no activation step", and that is an answer, not a gap."""
    from mofsbu.sites.model import Site

    pyridyl = Site(atom_idx=0, donor_type="pyridyl_N", labile=False, charge_after=0,
                   live_dof="live", binding_modes=("mono",))
    record = activation_ease(pyridyl)
    assert record.components["deprotonation"] == 1.0
    assert "no activation step" in " ".join(record.notes)


def test_the_pka_anchor_orders_the_donor_classes_the_way_chemistry_does():
    """The one thing the scalar convention has to get right to be worth having."""
    from mofsbu.sites.model import Site

    def scalar(donor_type: str) -> float:
        return activation_ease(Site(atom_idx=0, donor_type=donor_type, labile=True,
                                    charge_after=-1, live_dof="live",
                                    binding_modes=("mono",))).scalar

    assert scalar("carboxylate_O") > scalar("phenolate_O") > scalar("alkoxide_O")
    assert PKA_CENTRE == 10.0
    assert 0.45 < scalar("phenolate_O") < 0.55, "phenol sits at the anchor point"


def test_a_geometry_may_not_be_stored_at_the_heuristic_rung(reg):
    """HEURISTIC is an ease rung, not a geometry rung (D18)."""
    mol = embed_molecule(mol_from_smiles("O"), seed=7)
    graph = from_rdkit(mol, charge=0, multiplicity=1, name="water")
    put = put_structure(reg, graph, tags=["w"])
    with pytest.raises(RegistryError, match="HEURISTIC"):
        put_geometry(reg, put.id, to_xyz(mol, "water"), fidelity=Fidelity.HEURISTIC,
                     method=FF)


# ── inheritance, which is what M5 consumes ───────────────────────────────────

def test_inherited_sites_are_remapped_and_consumed_ones_dropped():
    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    parent = perceive(mol)
    assert len(parent) >= 2
    consumed = {parent[0].atom_idx}
    atom_map = {s.atom_idx: s.atom_idx + 100 for s in parent}

    child = inherit_sites(parent, atom_map, consumed)
    assert len(child) == len(parent) - 1
    assert all(s.atom_idx >= 100 for s in child)
    # the frame came across untouched — that is the point (D13)
    survivor = next(s for s in parent if s.atom_idx != parent[0].atom_idx)
    inherited = next(s for s in child if s.atom_idx == survivor.atom_idx + 100)
    assert inherited.frame == survivor.frame
    assert inherited.donor_type == survivor.donor_type


def test_an_atom_not_in_the_product_loses_its_site():
    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    parent = perceive(mol)
    partial = {s.atom_idx: s.atom_idx for s in parent[1:]}   # first atom left behind
    assert len(inherit_sites(parent, partial)) == len(parent) - 1


def test_two_parents_claiming_one_product_atom_is_a_broken_atom_map():
    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    sites = perceive(mol)
    a = inherit_sites(sites[:1], {sites[0].atom_idx: 0})
    b = inherit_sites(sites[1:2], {sites[1].atom_idx: 0})
    with pytest.raises(ValueError, match="not injective"):
        merge_inherited(a, b)


def test_open_sites_refuses_a_block_with_no_state():
    from mofsbu.assembly.join import BuildingBlock, NotBuiltYet

    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    graph = from_rdkit(mol, charge=None, multiplicity=1, name="sal")
    block = BuildingBlock(graph=graph, sites=tuple(perceive(mol)))
    with pytest.raises(NotBuiltYet, match="no site state"):
        block.open_sites()


def test_open_sites_uses_the_state_when_it_has_one():
    from mofsbu.assembly.join import BuildingBlock

    mol = embed_molecule(mol_from_smiles(SALICYLIC), seed=7)
    graph = from_rdkit(mol, charge=None, multiplicity=1, name="sal")
    sites = perceive(mol)
    conf = mol.GetConformer()
    coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
               conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    states = refresh_state(sites, coords, graph=graph)
    by_idx = {s.atom_idx: s for s in states}
    block = BuildingBlock(graph=graph, sites=tuple(sites), state=by_idx)
    opened = block.open_sites()
    assert {o.atom_idx for o in opened} <= {s.atom_idx for s in sites}
    assert all(by_idx[o.atom_idx].status is SiteStatus.OPEN for o in opened)
    assert len(opened) == sum(1 for s in states if s.is_open)
