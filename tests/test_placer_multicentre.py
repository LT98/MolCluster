"""M6 — the multi-centre placer.  Its exit gate, written before the placer exists.

PLAN §7 item 3 says fixtures before code, because that is what made M2 cheap: it forces
the commitments at the point where they are still free.  This module is that move for M6.
Three kinds of test live here and they are not interchangeable:

* **Ground truth that holds today.**  What the target graphs are, and whether the curated
  numbers in `data/reference/node_cases.tsv` are mutually consistent.  A typo in a
  literature value would otherwise be discovered as a pile of clashes out of a placer that
  is working correctly, which is the most expensive possible way to find it.
* **The exit gate** (PLAN §2, M6): a placer-built paddlewheel and Fe3-oxo land on the
  **same L1 hash as the hand-written fixtures**, and QC passes.  Marked `@m6`.
* **Invariants and negatives** — D12's "no mononuclear branch", D13's determinism and
  replay, and the refusals.  Also `@m6`.

`@m6` is a STRICT xfail on `NotBuiltYet`.  That is deliberate and it has three effects:
the suite stays green today; a stub that starts returning something plausible instead of
raising turns the test red rather than passing quietly; and the moment M6 actually works,
every one of these XPASSes and the suite fails until the marker is removed.  Nobody can
land M6 and leave the gate unexercised.

Where a test needs a seam that is not built, it calls the seam by its settled name and
lets `NotBuiltYet` out.  Those names are in `geometry/placer.py` and `geometry/qc.py`
alongside their raising bodies (ground rule 8).
"""
from __future__ import annotations

import numpy as np
import pytest

import fixtures as fx
from mofsbu.assembly.join import NotBuiltYet
from mofsbu.descriptors.tables import known_donor_types
from mofsbu.geometry.distances import BASE_MO, base_distance
from mofsbu.geometry.embed import embed_molecule
from mofsbu.geometry.placer import (
    Center, InterCentreConstraint, Join, LigandPlacement, place_mononuclear,
    place_multicentre, site_vectors, to_rdkit_multicentre,
)
from mofsbu.geometry.qc import check_intercentre, clash_limit, qc
from mofsbu.graph import EdgeType
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles
from mofsbu.identity import l1_graph_hash
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.model import perceive
from node_cases import NodeCase, carboxylate_o_o, case, load_node_cases

CASES = load_node_cases()
GATE_CASES = [c for c in CASES if c.fixture]

m6 = pytest.mark.xfail(raises=NotBuiltYet, strict=True,
                       reason="M6: the multi-centre placer is a settled interface with a "
                              "scheduled body.  Remove this marker when it has one.")


# ── the ligands every node in the set is built from ──────────────────────────

def formate(name: str) -> LigandPlacement:
    """HCOO-, presenting both oxygens.  The bridge in every row of the table.

    Its donors are given as a pair and the JOINS decide what that pair means: two joins on
    one centre is a chelate, two joins on different centres is a mu2 bridge.  Nothing on
    the ligand says which, which is D14 held one layer below the graph.
    """
    mol = embed_molecule(mol_from_smiles("[O-]C=O"), seed=7)
    idxs = tuple(s.atom_idx for s in perceive(mol) if s.donor_type == "carboxylate_O")
    assert len(idxs) == 2, f"formate perceived {len(idxs)} carboxylate oxygens"
    return LigandPlacement(mol=mol, donor_idxs=idxs,
                           donor_types=("carboxylate_O",) * 2,
                           mode=BindingMode.BRIDGE_MU2, name=name)


def central_oxo(name: str = "mu-O") -> LigandPlacement:
    """The bare oxide at the middle of the Fe3 and Zn4 nodes.

    Typed `hydroxide_O` because that is what perception calls a lone oxygen today and
    there is no `oxo_O` row in `donor_descriptors.tsv` — see
    `test_the_central_oxo_has_no_donor_type_of_its_own_yet`, which is here so that gap is
    not discovered by the placer instead.
    """
    mol = embed_molecule(mol_from_smiles("[O-2]"), seed=7)
    return LigandPlacement(mol=mol, donor_idxs=(0,), donor_types=("hydroxide_O",),
                           mode=BindingMode.BRIDGE_MU3, name=name)


def _centres(c: NodeCase, charge: int, oxidation_state: int, spin_class: str,
             n: int | None = None) -> list[Center]:
    return [Center(c.metal, cn=c.cn, local_geometry=c.local_geometry, charge=charge,
                   oxidation_state=oxidation_state, spin_class=spin_class)
            for _ in range(n if n is not None else c.n_metals)]


def _bridges(pairs: list[tuple[int, int]]) -> list[Join]:
    """One carboxylate per centre pair: two joins naming the same block."""
    joins: list[Join] = []
    for k, (a, b) in enumerate(pairs):
        lig = formate(f"HCOO{k}")
        joins.append(Join(center=a, block=lig, site=0, mode=BindingMode.BRIDGE_MU2))
        joins.append(Join(center=b, block=lig, site=1, mode=BindingMode.BRIDGE_MU2))
    return joins


def paddlewheel_inputs(c: NodeCase | None = None):
    """Cu2(mu-O2CH)4 as (centres, joins, constraints) — the first ground-truth target."""
    c = c or case("cu_paddlewheel")
    centres = _centres(c, charge=2, oxidation_state=2, spin_class="hs")
    joins = _bridges([(0, 1)] * c.n_bridges)
    constraints = [InterCentreConstraint((0, 1), mm_distance=c.d_mm,
                                         mm_lo=c.d_mm_lo, mm_hi=c.d_mm_hi,
                                         metal_metal_bond=c.mm_bond)]
    return centres, joins, constraints


def fe3_oxo_inputs(c: NodeCase | None = None):
    """Fe3(mu3-O)(mu-O2CH)6 — two carboxylates over each edge, one oxo over all three."""
    c = c or case("fe3_mu3_oxo_333")
    centres = _centres(c, charge=3, oxidation_state=3, spin_class="hs")
    edges = [(0, 1), (1, 2), (2, 0)]
    joins = _bridges([e for e in edges for _ in range(2)])
    oxo = central_oxo()
    joins += [Join(center=k, block=oxo, site=0, mode=BindingMode.BRIDGE_MU3)
              for k in range(3)]
    constraints = [InterCentreConstraint(e, mm_distance=c.d_mm, mm_lo=c.d_mm_lo,
                                         mm_hi=c.d_mm_hi, metal_metal_bond=c.mm_bond)
                   for e in edges]
    return centres, joins, constraints


def zn4o_inputs(c: NodeCase | None = None):
    """Zn4O(mu-O2CH)6 — one carboxylate over each of the six tetrahedron edges."""
    c = c or case("zn4o")
    centres = _centres(c, charge=2, oxidation_state=2, spin_class="ls")
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    joins = _bridges(edges)
    oxo = central_oxo()
    joins += [Join(center=k, block=oxo, site=0, mode=BindingMode.BRIDGE_MU3)
              for k in range(4)]
    constraints = [InterCentreConstraint(e, mm_distance=c.d_mm, mm_lo=c.d_mm_lo,
                                         mm_hi=c.d_mm_hi, metal_metal_bond=c.mm_bond)
                   for e in edges]
    return centres, joins, constraints


NODE_BUILDERS = {"cu_paddlewheel": paddlewheel_inputs,
                 "fe3_mu3_oxo_333": fe3_oxo_inputs,
                 "zn4o": zn4o_inputs}


def built_graph(name: str):
    """Place a node and type it — the whole path the exit gate measures.

    Placer, graph extraction and identity are all in here on purpose.  The gate is that
    the three of them AGREE, and a test of the placer alone could not say that.
    """
    c = case(name)
    centres, joins, constraints = NODE_BUILDERS[name](c)
    result = place_multicentre(centres, joins, constraints, seed=0)
    mol = to_rdkit_multicentre(centres, joins, constraints, result)
    fixture = fx.ALL[c.fixture]()
    return from_rdkit(mol, charge=fixture.charge, multiplicity=fixture.multiplicity,
                      name=name), result


# ══ 1. ground truth: what the curated table says, checked against itself ══════

def test_the_table_covers_the_nodes_the_exit_gate_names():
    """PLAN §2 names the paddlewheel and the Fe3-oxo trimer as M6's ground truth."""
    names = {c.name for c in CASES}
    assert {"cu_paddlewheel", "fe3_mu3_oxo_333"} <= names
    assert len(GATE_CASES) >= 3, "a hash gate needs a fixture to hash against"


def test_every_case_naming_a_fixture_has_one(all_fixtures):
    for c in GATE_CASES:
        assert c.fixture in all_fixtures, f"{c.name} names a fixture that does not exist"


def test_the_gate_targets_are_pinned_in_the_golden_hashes(golden_path):
    """The hash a built node must reproduce has to be a FIXED number, not one recomputed
    from the same fixture in the same process.  Otherwise a drift in the certificate
    recipe moves the target and the gate silently follows it."""
    import json

    golden = json.loads(golden_path.read_text())
    for c in GATE_CASES:
        assert c.fixture in golden["l1"], f"{c.fixture} has no golden L1 hash"


def test_the_mu_oxo_targets_are_geometrically_consistent():
    """For a symmetric mu-oxo node, M...M = 2 d sin(theta/2) is a definition, not a fit.

    So this is a check on the TABLE: if the M-oxo distance, the M-O-M angle and the M...M
    distance in a row do not satisfy it, one of the three was typed wrong, and a placer
    handed all three would be given a triangle that does not close.
    """
    checked = 0
    for c in CASES:
        if not c.has_central_oxo:
            continue
        assert abs(c.mm_from_oxo() - c.d_mm) < 0.05, (
            f"{c.name}: oxo geometry implies M...M {c.mm_from_oxo():.3f} A but the row "
            f"says {c.d_mm:.2f} ({c.source})")
        checked += 1
    assert checked >= 2


def test_a_syn_syn_carboxylate_can_span_every_target():
    """The bridge is the same 2.22 A bite on every row; only its tilt changes.

    M...M runs from 2.09 A (Mo-Mo, quadruple bond) to 3.34 A (mixed-valence Fe3) across
    this table — a 1.25 A range spanned by one rigid group, which it manages by leaning
    the M-M-O angle from just past perpendicular down to about 75 deg.  A row the bridge
    cannot reach at all makes `bridge_tilt_deg` raise, and that is the failure mode worth
    catching here rather than in a clash report.
    """
    assert 2.20 < carboxylate_o_o() < 2.25
    for c in CASES:
        tilt = c.bridge_tilt_deg()
        assert 70.0 <= tilt <= 95.0, (
            f"{c.name}: a syn-syn carboxylate would need M-M-O {tilt:.0f} deg to span "
            f"M...M {c.d_mm} A at M-O {c.d_m_o} A")


def test_the_window_brackets_the_ideal_on_every_row():
    for c in CASES:
        assert c.d_mm_lo < c.d_mm < c.d_mm_hi, f"{c.name}: ideal outside its own window"
        assert c.d_mm_hi - c.d_mm_lo <= 0.40, (
            f"{c.name}: a {c.d_mm_hi - c.d_mm_lo:.2f} A window is wide enough to accept "
            "a node that is simply wrong")


def test_every_case_asks_for_a_local_geometry_that_exists_at_its_cn():
    """`cn` is the polyhedron's, vacancies included — the same meaning `place_mononuclear`
    gave it, so a node that leaves axial sites open still names a real polyhedron."""
    for c in CASES:
        assert site_vectors(c.local_geometry, c.cn).shape == (c.cn, 3)


def test_the_M_O_targets_agree_with_the_distance_model():
    """The curated bridging M-O distance and `geometry.distances` must be the same claim.

    They are derived independently — one from node crystallography, one from a calibrated
    element-offset model — so agreeing is evidence and disagreeing is a question worth
    answering before the placer uses either.
    """
    for c in CASES:
        if c.metal not in BASE_MO:
            continue
        assert abs(c.d_m_o - base_distance(c.metal)) <= 0.15, (
            f"{c.name}: table says M-O {c.d_m_o}, distances model says "
            f"{base_distance(c.metal):.2f}")


def test_the_distance_model_does_not_know_every_metal_in_the_table():
    """A named gap, so it is a decision rather than a surprise.

    Rh and Mo have no `BASE_MO` row, so a paddlewheel built for either would silently get
    the 2.05 A default. Neither is on M6's critical path — the exit gate is Cu and Fe —
    but a family row exists for each so the gap is visible. Add the rows and update this
    test in the same commit.
    """
    missing = sorted({c.metal for c in CASES if c.metal not in BASE_MO})
    assert missing == ["Mo", "Rh"], (
        f"the set of metals `geometry.distances` cannot place has changed: {missing}")


# ══ 2. ground truth: what the target GRAPHS are ══════════════════════════════

def test_the_paddlewheel_target_is_one_node_with_four_mu2_bridges():
    g = fx.cu_paddlewheel()
    assert len(g.metals()) == 2
    assert any(t is EdgeType.METAL_METAL for _i, _j, t in g.edges())
    frags = g.ligand_fragments()
    assert len(frags) == 4
    assert all(g.fragment_bridge_class(f).value == "mu2" for f in frags)


def test_the_fe3_target_is_held_together_without_a_metal_metal_edge():
    """The discrimination that makes `metal_metal_bond` a declared input rather than
    something a placer could read off a distance: this node's Fe...Fe is 3.29 A and
    carries no edge, the paddlewheel's Cu-Cu is 2.62 A and carries one."""
    g = fx.fe3_mu3_oxo((3, 3, 3))
    assert len(g.metals()) == 3
    assert not any(t is EdgeType.METAL_METAL for _i, _j, t in g.edges())
    assert g.max_bridge_class().value == "mu3"


def test_the_zn4o_target_reaches_four_metals_through_one_oxygen():
    """mu-ness is COUNTED, not tagged (D14).  One graph carrying both MU_N and MU2 is the
    cheapest demonstration that nothing anywhere wrote the bridge order down."""
    g = fx.zn4o()
    assert len(g.metals()) == 4
    assert g.max_bridge_class().value == "muN"
    frags = g.ligand_fragments()
    oxo = [f for f in frags if len(f) == 1]
    assert len(oxo) == 1 and len(g.bridging_metals(oxo[0][0])) == 4
    assert sum(1 for f in frags if g.fragment_bridge_class(f).value == "mu2") == 6


@pytest.mark.parametrize("c", GATE_CASES, ids=[c.name for c in GATE_CASES])
def test_the_table_and_the_fixture_graph_agree(c: NodeCase):
    """Every structural column is checkable against the hand-written graph.  They are two
    statements of the same node, so a disagreement means one of them is a typo — and the
    table is what the placer will be driven by."""
    g = fx.ALL[c.fixture]()
    assert len(g.metals()) == c.n_metals
    assert {g.label(i).element for i in g.metals()} == {c.metal}
    has_mm = any(t is EdgeType.METAL_METAL for _i, _j, t in g.edges())
    assert has_mm is c.mm_bond, f"{c.name}: mm_bond column disagrees with the fixture"
    bridges = sum(1 for f in g.ligand_fragments()
                  if g.fragment_bridge_class(f).value == "mu2")
    assert bridges == c.n_bridges


@pytest.mark.parametrize("c", GATE_CASES, ids=[c.name for c in GATE_CASES])
def test_the_vacancy_count_follows_from_cn_and_the_fixture(c: NodeCase):
    """`n_vacancies` is not free bookkeeping: CN minus the dative bonds the fixture
    actually draws at a metal is how many vertices the placer must leave empty.  A
    coordinatively unsaturated centre is not a smaller centre."""
    g = fx.ALL[c.fixture]()
    for metal in g.metals():
        occupied = len(g.neighbors(metal, EdgeType.DATIVE))
        assert c.cn - occupied == c.n_vacancies, (
            f"{c.name}: metal {metal} has {occupied} donors at CN {c.cn}, so "
            f"{c.cn - occupied} vertices are vacant, not {c.n_vacancies}")


# ══ 3. prerequisites M6 has to land before the gate can pass ═════════════════

def test_the_central_oxo_has_no_donor_type_of_its_own_yet():
    """A bare oxide perceives as `hydroxide_O`, which is a different ligand.

    The mu3/mu4 oxo is a donor type with no row in `donor_descriptors.tsv`, and a type
    with no row fails a descriptor test (AGENTS.md §7), so M6 cannot simply start emitting
    one.  **Delete this test in the commit that adds the `oxo_O` row**, and update
    `central_oxo()` above with it.
    """
    assert "oxo_O" not in known_donor_types(), (
        "oxo_O now has a descriptor row — point central_oxo() at it and delete this test")


def test_the_multicentre_interfaces_exist_and_all_raise():
    """Ground rule 8, with REAL arguments.

    `test_jobs.py` already checks these raise when handed junk. This checks they raise
    when handed a fully-formed paddlewheel — which is the case where a stub that had
    quietly learned to return something plausible would not be caught.
    """
    centres, joins, constraints = paddlewheel_inputs()
    for call in (lambda: place_multicentre(centres, joins, constraints),
                 lambda: to_rdkit_multicentre(centres, joins, constraints, None),
                 lambda: check_intercentre(np.zeros((2, 3)), constraints,
                                           metal_idxs=[0, 1])):
        with pytest.raises(NotBuiltYet, match="M6"):
            call()


def test_qc_today_cannot_see_a_node_whose_centres_are_far_too_close():
    """Why `check_intercentre` has to exist at all, demonstrated rather than asserted.

    Two Cu at 1.90 A is roughly 0.7 A inside any real paddlewheel and would tear apart
    under a relax. `check_clashes` allows the pair 1.74 A and `check_metal_bonds` only
    looks at metal-DONOR pairs, so the existing report comes back clean. A node built at
    the wrong M-M distance is therefore indistinguishable from a good one today.
    """
    symbols = ["Cu", "Cu"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.90]])
    assert clash_limit("Cu", "Cu") < 1.90            # below the clash floor, so invisible
    report = qc(symbols, coords, bonded=set(), metal_idx=0, donor_idxs=[])
    assert report.ok, "qc already catches this — rewrite the M6 intercentre tests"


def test_a_qc_report_can_carry_centre_centre_findings():
    """The M6 verdict travels in the same `QCReport` as everything else, so a caller does
    not have to know to ask a second object whether the node is sane."""
    from mofsbu.geometry.qc import BadIntercentre, QCReport

    bad = BadIntercentre((0, 1), 1.90, 2.50, 2.75, ("Cu", "Cu"), "cu_paddlewheel")
    report = QCReport(ok=False, bad_intercentres=[bad])
    assert report.code == "qc_intercentre"
    assert not report.marginal, "a node built at the wrong M-M distance is not a near miss"
    assert "Cu" in report.elements_involved()
    assert "2.50-2.75" in str(report)


# ══ 4. the exit gate ═════════════════════════════════════════════════════════

@m6
@pytest.mark.parametrize("name", sorted(NODE_BUILDERS), ids=sorted(NODE_BUILDERS))
def test_a_placed_node_hashes_to_the_hand_written_fixture(name, golden_path):
    """**The strong exit gate** (PLAN §2, M6).

    One assertion proves three things at once: the placer built the right connectivity,
    the graph extraction read it back faithfully, and identity agrees with both. Any one
    of them wrong and the hash moves. The target is the GOLDEN hash, not a freshly
    recomputed one, so a drift in the certificate recipe cannot move the goalposts.
    """
    import json

    golden = json.loads(golden_path.read_text())
    g, _result = built_graph(name)
    assert l1_graph_hash(g) == golden["l1"][case(name).fixture]


@m6
@pytest.mark.parametrize("name", sorted(NODE_BUILDERS), ids=sorted(NODE_BUILDERS))
def test_a_placed_node_passes_qc(name):
    """The other half of the gate: no clashes, and the centres where the table says.

    A node can hash correctly and still be geometric nonsense — the hash sees connectivity
    and nothing else. This is the test that the coordinates mean something.
    """
    c = case(name)
    _g, result = built_graph(name)
    assert result.report.ok, str(result.report)
    assert not result.report.bad_intercentres
    metals = [i for i, s in enumerate(result.symbols) if s == c.metal]
    assert len(metals) == c.n_metals
    for a, b in [(metals[i], metals[j])
                 for i in range(len(metals)) for j in range(i + 1, len(metals))]:
        d = float(np.linalg.norm(result.coords[a] - result.coords[b]))
        assert c.d_mm_lo <= d <= c.d_mm_hi, (
            f"{name}: {c.metal}...{c.metal} {d:.2f} A, window "
            f"{c.d_mm_lo}-{c.d_mm_hi} ({c.source})")


@m6
def test_the_paddlewheel_keeps_its_two_axial_vacancies():
    """The fixture is the anhydrous node: each Cu's apical vertex is open.

    Losing them would not change the hash — a vacancy is not a graph edge — so nothing in
    the gate above would notice, and the next assembly step has nowhere to attach. Note
    this needs `vacancies` to say WHICH centre each belongs to; today it is a flat tuple
    of directions, which is sufficient for one centre and not for two.
    """
    centres, joins, constraints = paddlewheel_inputs()
    result = place_multicentre(centres, joins, constraints)
    assert len(result.vacancies) == 2 * case("cu_paddlewheel").n_vacancies
    a, b = (np.array(v) for v in result.vacancies[:2])
    cos = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
    assert cos < -0.9, "the two axial sites should point in opposite directions"


# ══ 5. invariants and refusals ═══════════════════════════════════════════════

@m6
def test_one_centre_through_the_multicentre_path_matches_place_mononuclear():
    """D12 and ground rule 4: mononuclear is N=1, not a special case.

    The placer docstring already claims there is no separate mononuclear code path to
    reconcile. This is the test that would fail if M6 quietly grew one — and it is the
    cheapest guard there is against the `if n_metals == 1` branch the ground rules forbid.
    """
    ligands = [formate("HCOO0")]
    mono = place_mononuclear("Zn", ligands, geometry="tetrahedral", cn=4)
    centres = [Center("Zn", cn=4, local_geometry="tetrahedral", charge=2,
                      oxidation_state=2, spin_class="ls")]
    joins = [Join(center=0, block=ligands[0], site=s, mode=BindingMode.CHELATE)
             for s in (0, 1)]
    multi = place_multicentre(centres, joins, [])
    assert np.allclose(mono.coords, multi.coords, atol=1e-9)


@m6
def test_a_node_is_deterministic_and_emits_its_choice_vector():
    """D13 — `construct` is a deterministic function of (choice-vector, seed) and says so."""
    a = place_multicentre(*paddlewheel_inputs(), seed=0)
    b = place_multicentre(*paddlewheel_inputs(), seed=0)
    assert np.allclose(a.coords, b.coords)
    assert a.choice_vector == b.choice_vector
    cv = a.choice_vector
    assert cv["centres"] and len(cv["centres"]) == 2
    assert any("mm_distance" in str(k) or "mm_distance" in str(v)
               for k, v in cv.items()), "the M-M target has to be in the replay recipe"


@m6
def test_replaying_a_nodes_choice_vector_reproduces_its_coordinates():
    """If this fails, a stored polynuclear node is frozen coordinates rather than a
    regenerable object — which is M5's exit gate 3, one level up."""
    centres, joins, constraints = paddlewheel_inputs()
    first = place_multicentre(centres, joins, constraints, seed=0)
    for join, emitted in zip(joins, first.choice_vector["ligands"]):
        join.block.azimuth_step = emitted.get("azimuth_step")
        join.block.oop_step = emitted.get("oop_step")
    second = place_multicentre(centres, joins, constraints, seed=0)
    assert np.allclose(first.coords, second.coords)


@m6
def test_the_metal_metal_edge_is_declared_and_never_inferred_from_distance():
    """Two paddlewheels that differ only in `metal_metal_bond` must be different
    structures, and only the declared one may match the fixture.

    If a placer decided the edge from how close the metals came out, this would collapse
    into one hash and identity would be a function of a bond-length table.
    """
    c = case("cu_paddlewheel")
    hashes = []
    for bonded in (True, False):
        centres, joins, _ = paddlewheel_inputs()
        constraints = [InterCentreConstraint((0, 1), mm_distance=c.d_mm,
                                             mm_lo=c.d_mm_lo, mm_hi=c.d_mm_hi,
                                             metal_metal_bond=bonded)]
        result = place_multicentre(centres, joins, constraints)
        mol = to_rdkit_multicentre(centres, joins, constraints, result)
        hashes.append(l1_graph_hash(from_rdkit(mol, charge=0, multiplicity=1)))
    assert hashes[0] != hashes[1]
    assert hashes[0] == l1_graph_hash(fx.cu_paddlewheel())


@m6
def test_a_centre_with_no_local_geometry_raises_rather_than_defaulting():
    """Ground rule 5 / invariant 4: ambiguity branches, it does not default.

    There is no answer to "what shape is a CN-5 Cu" that a placer is entitled to pick on
    its own — both CN-5 polyhedra are real and they build different molecules.
    """
    from mofsbu._types import AmbiguousSpecError

    _centres_ok, joins, constraints = paddlewheel_inputs()
    centres = [Center("Cu", cn=5, local_geometry="", charge=2, oxidation_state=2,
                      spin_class="hs") for _ in range(2)]
    with pytest.raises(AmbiguousSpecError):
        place_multicentre(centres, joins, constraints)


@m6
def test_a_constraint_naming_a_centre_that_does_not_exist_is_refused():
    centres, joins, _ = paddlewheel_inputs()
    with pytest.raises((ValueError, IndexError)):
        place_multicentre(centres, joins, [InterCentreConstraint((0, 7), mm_distance=2.6)])


# ══ 6. check_intercentre on its own, so M6 can be built in two steps ═════════
#
# These call nothing but the QC seam, so they go green as soon as `check_intercentre`
# has a body — before the placer does.  Coordinates are written by hand here for exactly
# that reason: a QC test that needs a working placer to produce its input cannot fail
# independently of it, and then it is not really testing QC.

def _two_centres(d: float) -> np.ndarray:
    return np.array([[0.0, 0.0, 0.0], [0.0, 0.0, d]])


@m6
def test_check_intercentre_accepts_a_node_at_its_literature_distance():
    c = case("cu_paddlewheel")
    found = check_intercentre(
        _two_centres(c.d_mm),
        [InterCentreConstraint((0, 1), mm_distance=c.d_mm, mm_lo=c.d_mm_lo,
                               mm_hi=c.d_mm_hi, metal_metal_bond=True)],
        metal_idxs=[0, 1], symbols=["Cu", "Cu"])
    assert found == []


@m6
@pytest.mark.parametrize("d", [1.90, 3.60], ids=["squeezed", "stretched"])
def test_check_intercentre_refuses_a_node_outside_its_window(d):
    """Both directions, because only one of them is ever a clash.

    A squeezed node might eventually trip the clash check with a big enough metal; a
    stretched one never will, and it is the one that quietly produces two mononuclear
    complexes sharing a record.
    """
    c = case("cu_paddlewheel")
    found = check_intercentre(
        _two_centres(d),
        [InterCentreConstraint((0, 1), mm_distance=c.d_mm, mm_lo=c.d_mm_lo,
                               mm_hi=c.d_mm_hi, metal_metal_bond=True)],
        metal_idxs=[0, 1], symbols=["Cu", "Cu"])
    assert len(found) == 1
    assert found[0].centres == (0, 1)
    assert abs(found[0].distance - d) < 1e-6
    assert (found[0].lo, found[0].hi) == (c.d_mm_lo, c.d_mm_hi)


@m6
def test_check_intercentre_falls_back_to_a_tolerance_when_no_window_is_given():
    """A constraint may carry only an ideal — `mm_lo`/`mm_hi` are the curated refinement
    of it, not a requirement.  Absent is not zero: a missing window must not read as a
    zero-width one that refuses everything."""
    found = check_intercentre(
        _two_centres(2.62), [InterCentreConstraint((0, 1), mm_distance=2.62)],
        metal_idxs=[0, 1], symbols=["Cu", "Cu"])
    assert found == []
