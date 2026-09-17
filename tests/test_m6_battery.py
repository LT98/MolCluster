"""The M6 battery — its targets, checked against each other before anything builds them.

This is [`WORKPLAN_M6.md`](../docs/WORKPLAN_M6.md) **S0(b)**: "the battery, and the
perception fixes it depends on, as failing tests".  It is not the exit gates — those are
S7's `test_m6_exit_gates.py`, across the whole battery, and they do not exist yet.

`data/reference/node_cases.tsv` is the battery as data, and adding a row needs no change
here.  What this module asks of those rows:

* **Are they self-consistent?**  A bridging atom treated as a centre (§4) fixes M···M by
  `2 d sin(θ/2)`; that is a definition, so checking it checks the *table*.  A typo in a
  target would otherwise be discovered as a pile of clashes out of a placer that is working
  correctly, which is the most expensive possible way to find it.
* **Do they agree with the fixture graphs?**  `mm_bond`, nuclearity, bridge count and
  vacancy count are all recoverable from `mofsbu.examples`, and a disagreement means one of
  the two is a typo.
* **Do they agree with the code?**  §4's µ3 = 3.291 Å and µ4 = 3.168 Å are reproduced here
  from `site_vectors` rather than copied, so a change to a local geometry moves the test.
* **Are the stated blockers real?**  Every `blocked_on` row is checked against the defect
  it names, so the battery cannot quietly claim to be waiting on something already fixed.

**D20 governs what is an input.**  M···M is an output to validate, not a constraint to
impose: the paddlewheel builds from well-1 frames and `kabsch` with no solver, and its
`collision()` is `None` because only one determinant fixes its distance.  The `@m6` tests
here are confined to the case that *does* need `place_multicentre` — an oxo-centred cluster
where two determinants disagree by about half an ångström — and to `check_intercentre`,
which validates afterwards.

`@m6` is a STRICT xfail on `NotBuiltYet`: the suite stays green today, a stub that starts
answering plausibly instead of raising turns red rather than passing quietly, and the moment
S3/S6 land these XPASS and the suite fails until the marker comes off.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

import fixtures as fx
from mofsbu.assembly.join import NotBuiltYet
from mofsbu.geometry.distances import BASE_MO, base_distance
from mofsbu.geometry.embed import embed_molecule
from mofsbu.geometry.placer import Center, InterCentreConstraint, place_multicentre, site_vectors
from mofsbu.geometry.qc import check_intercentre, clash_limit, qc
from mofsbu.graph import EdgeType
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.model import perceive
from node_cases import NodeCase, carboxylate_o_o, case, load_node_cases

CASES = load_node_cases()
GATE_CASES = [c for c in CASES if c.fixture]
DETERMINED = [c for c in CASES if c.determined]

m6 = pytest.mark.xfail(raises=NotBuiltYet, strict=True,
                       reason="M6 S3/S6: settled interface, scheduled body.  Remove this "
                              "marker when it has one.")


def skeleton_mm(c: NodeCase) -> list[float]:
    """The M...M distances a bridging centre's own polyhedron produces.

    Computed from `site_vectors` rather than from the arithmetic in `node_cases`, so this
    is the code's answer and not the table's — which is the point of comparing them.
    """
    v = site_vectors(c.mu_geometry, c.n_metals, c.d_m_mu)
    return [float(np.linalg.norm(v[i] - v[j]))
            for i in range(len(v)) for j in range(i + 1, len(v))]


# ══ 1. the battery is the one WORKPLAN_M6 §6 describes ═══════════════════════

def test_the_table_covers_both_mechanisms_and_three_nuclearities():
    """§1 draws the A/B line first, because a milestone that treats them as one thing will
    build A, declare victory on the paddlewheel, and discover B at the trimer."""
    assert {c.mechanism for c in CASES} >= {"A", "B", "A+B"}
    assert {c.n_metals for c in CASES} >= {2, 3, 4}
    assert {c.mm_bond for c in CASES} == {True, False}
    assert len([c for c in CASES if c.blocked_on]) >= 2, "the known gaps should be rows"


def test_every_case_naming_a_fixture_has_one(all_fixtures):
    for c in GATE_CASES:
        assert c.fixture in all_fixtures, f"{c.name} names a fixture that does not exist"


def test_the_gate_targets_are_pinned_in_the_golden_hashes(golden_path):
    """Gate 2 compares a built node against a FIXED hash, not one recomputed from the same
    fixture in the same process — otherwise a drift in the certificate recipe moves the
    target and the gate silently follows it."""
    golden = json.loads(golden_path.read_text())
    for c in GATE_CASES:
        assert c.fixture in golden["l1"], f"{c.fixture} has no golden L1 hash"
    assert golden["l1"]["cu_paddlewheel"].startswith("aed8d6a834418f38")
    assert golden["l1"]["fe3_mu3_oxo_333"].startswith("9a3022deefc43dae")


# ══ 2. mechanism B: the bridging atom as a centre (§4) ═══════════════════════

@pytest.mark.parametrize("c", DETERMINED, ids=[c.name for c in DETERMINED])
def test_a_bridging_centres_geometry_determines_the_skeleton(c: NodeCase):
    """`M...M = 2 d sin(theta/2)` is a definition for a symmetric bridge, so this checks the
    TABLE: if the M-bridge distance, the angle and the M...M in a row disagree, one of the
    three was typed wrong and the skeleton it describes does not close."""
    assert abs(c.mm_from_bridging_centre() - c.d_mm) < 0.06, (
        f"{c.name}: the bridging centre implies M...M "
        f"{c.mm_from_bridging_centre():.3f} A but the row says {c.d_mm:.2f}")


@pytest.mark.parametrize("c", DETERMINED, ids=[c.name for c in DETERMINED])
def test_the_skeleton_the_code_produces_is_the_one_the_table_claims(c: NodeCase):
    """§4's reframing, checked against `site_vectors` instead of quoted.

    `site_vectors` and `vacancy_sites` are already generic over which atom sits at the
    centre, so mechanism B is mostly a change of viewpoint — and the viewpoint lands the
    numbers exactly: µ3-O trigonal at 1.90 A gives 3.291, µ4-O tetrahedral at 1.94 gives
    3.168. A change to either local geometry moves this test, which is what makes it a
    check on the code rather than a restatement of the workplan.
    """
    if c.mu_geometry not in ("trigonal", "tetrahedral"):
        pytest.skip(f"{c.mu_geometry} is not expressible in GEOMETRIES yet")
    edges = skeleton_mm(c)
    assert max(edges) - min(edges) < 1e-9, "a symmetric bridge gives one edge length"
    # 5e-4, not zero: the table stores the angle to two decimals (109.47 for a tetrahedron)
    # while `site_vectors` builds from the exact vertices.  Tight enough that a wrong angle
    # fails, loose enough that a rounded one does not.
    assert abs(edges[0] - c.mm_from_bridging_centre()) < 5e-4


def test_the_mu3_and_mu4_skeletons_land_on_the_literature_values():
    """The two measurements WORKPLAN_M6 §4 reports, reproduced from the code."""
    assert abs(skeleton_mm(case("fe3_mu3_oxo_333"))[0] - 3.291) < 0.001
    assert abs(skeleton_mm(case("zn4o"))[0] - 3.168) < 0.001


def test_a_bent_bridging_centre_cannot_be_expressed_at_all():
    """§4's µ2 row fails for a stated reason: `GEOMETRIES` has only `linear` at CN 2.

    Same class of defect as the legacy table mapping CN 5 to 'planar', which this package
    already fixed once. Delete this test in the commit that adds the bent geometry.
    """
    from mofsbu.geometry.placer import GEOMETRIES

    at_cn2 = {name for name, by_cn in GEOMETRIES.items() if 2 in by_cn}
    assert at_cn2 == {"linear"}, f"CN 2 now offers {sorted(at_cn2)}"
    with pytest.raises(ValueError, match="unknown coordination geometry"):
        site_vectors("bent", 2)
    hydroxo = case("cu2_mu2_hydroxide")
    assert hydroxo.mu_geometry == "bent" and hydroxo.blocked_on
    assert "bent" in hydroxo.blocked_on


# ══ 3. mechanism A, and where the two collide (§5) ═══════════════════════════

def test_the_paddlewheel_has_only_one_determinant_and_so_never_collides():
    """Why the paddlewheel never needs `place_multicentre`, stated as a property of the
    table rather than as a claim in a docstring: nothing but the bridges fixes its M···M,
    so there is no second answer to reconcile."""
    pw = case("cu_paddlewheel")
    assert pw.mechanism == "A" and not pw.has_bridging_centre
    assert pw.collision() is None


def test_the_two_determinants_disagree_by_the_measured_amount():
    """§5, as data. Both numbers are measured, and their gap is M6's real work — 0.595 Å on
    Fe₃ and 0.491 Å on Zn₄O. A rigid ligand cannot open its O–C–O to close it, so this is a
    known limitation to report as strain, never to average away."""
    assert abs(case("fe3_mu3_oxo_333").collision() - 0.595) < 0.005
    assert abs(case("zn4o").collision() - 0.491) < 0.005


def test_the_span_is_a_property_of_the_bridging_ligand():
    """Formate offers 2.673 Å where benzoate offers 2.461 (§3). The same motif with a
    different bridge misses the literature Cu–Cu from the other side, which is why the
    benzoate analogue is in the battery rather than being assumed to follow."""
    formate = case("cu_paddlewheel").d_mm_bridge
    benzoate = case("cu2_benzoate_paddlewheel").d_mm_bridge
    assert abs(formate - 2.673) < 0.001 and abs(benzoate - 2.461) < 0.001
    assert formate > benzoate


def test_a_syn_syn_carboxylate_can_reach_every_target_that_has_one():
    """The bridge is the same 2.22 Å bite on every row; only its tilt against the M–M axis
    changes. A row it cannot reach at all makes `bridge_reach_deg` raise, and that is the
    failure mode worth catching here rather than in a clash report."""
    assert 2.20 < carboxylate_o_o() < 2.25
    checked = 0
    for c in CASES:
        if c.d_mm is None or c.d_m_o is None:
            continue
        reach = c.bridge_reach_deg()
        assert 70.0 <= reach <= 95.0, (
            f"{c.name}: a syn-syn carboxylate would need M-M-O {reach:.0f} deg to span "
            f"M...M {c.d_mm} A at M-O {c.d_m_o} A")
        checked += 1
    assert checked >= 5


# ══ 4. the targets against the code and the fixtures ═════════════════════════

def test_the_window_brackets_the_ideal_wherever_one_is_given():
    for c in CASES:
        if c.d_mm is None:
            assert c.d_mm_lo is None and c.d_mm_hi is None, (
                f"{c.name}: a window with no ideal — absent is not zero")
            continue
        assert c.d_mm_lo < c.d_mm < c.d_mm_hi, f"{c.name}: ideal outside its own window"
        assert round(c.d_mm_hi - c.d_mm_lo, 6) <= 0.40, (
            f"{c.name}: a {c.d_mm_hi - c.d_mm_lo:.2f} A window is wide enough to accept "
            "a node that is simply wrong")


def test_every_case_asks_for_a_local_geometry_that_exists_at_its_cn():
    """`cn` is the polyhedron's, vacancies included — the same meaning `place_mononuclear`
    gave it, so a node that leaves axial sites open still names a real polyhedron."""
    for c in CASES:
        assert site_vectors(c.local_geometry, c.cn).shape == (c.cn, 3)


def test_the_M_O_targets_agree_with_the_distance_model():
    """Curated bridging M–O and `geometry.distances` are derived independently — one from
    node crystallography, one from a calibrated element-offset model — so agreeing is
    evidence and disagreeing is a question worth answering before either is used."""
    for c in CASES:
        if c.d_m_o is None or c.metal not in BASE_MO:
            continue
        assert abs(c.d_m_o - base_distance(c.metal)) <= 0.15, (
            f"{c.name}: table says M-O {c.d_m_o}, distances model says "
            f"{base_distance(c.metal):.2f}")


def test_the_distance_model_does_not_know_every_metal_in_the_table():
    """A named gap, so it is a decision rather than a surprise. S4 adds
    `metal_metal_distance` following `metal_donor_distance` exactly; Rh and Mo are in the
    table so that function has more than one metal to be right about. Update this test in
    the commit that adds their rows."""
    missing = sorted({c.metal for c in CASES if c.metal not in BASE_MO})
    assert missing == ["Mo", "Rh"], (
        f"the set of metals `geometry.distances` cannot place has changed: {missing}")


@pytest.mark.parametrize("c", GATE_CASES, ids=[c.name for c in GATE_CASES])
def test_the_table_and_the_fixture_graph_agree(c: NodeCase):
    """Every structural column is recoverable from the hand-written graph. They are two
    statements of the same node, so a disagreement means one of them is a typo."""
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
    """`n_vacancies` is not free bookkeeping: CN minus the dative bonds the fixture draws
    at a metal is how many vertices a build must leave empty, and a coordinatively
    unsaturated centre is not a smaller centre."""
    g = fx.ALL[c.fixture]()
    for metal in g.metals():
        occupied = len(g.neighbors(metal, EdgeType.DATIVE))
        assert c.cn - occupied == c.n_vacancies, (
            f"{c.name}: metal {metal} has {occupied} donors at CN {c.cn}, so "
            f"{c.cn - occupied} vertices are vacant, not {c.n_vacancies}")


def test_the_paddlewheel_target_is_one_node_with_four_mu2_bridges():
    g = fx.cu_paddlewheel()
    assert len(g.metals()) == 2
    assert any(t is EdgeType.METAL_METAL for _i, _j, t in g.edges())
    frags = g.ligand_fragments()
    assert len(frags) == 4
    assert all(g.fragment_bridge_class(f).value == "mu2" for f in frags)


def test_the_fe3_target_is_held_together_without_a_metal_metal_edge():
    """The discrimination C10 has to settle: this node's Fe···Fe is 3.29 Å and carries no
    edge, the paddlewheel's Cu–Cu is 2.62 Å and carries one. Both are in the certificate,
    so a distance threshold applied silently would be a wrong identity."""
    g = fx.fe3_mu3_oxo((3, 3, 3))
    assert len(g.metals()) == 3
    assert not any(t is EdgeType.METAL_METAL for _i, _j, t in g.edges())
    assert g.max_bridge_class().value == "mu3"


def test_the_zn4o_target_reaches_four_metals_through_one_oxygen():
    """µ-ness is COUNTED, not tagged (D14). One graph carrying both MU_N and MU2 is the
    cheapest demonstration that nothing anywhere wrote the bridge order down."""
    g = fx.zn4o()
    assert len(g.metals()) == 4
    assert g.max_bridge_class().value == "muN"
    frags = g.ligand_fragments()
    oxo = [f for f in frags if len(f) == 1]
    assert len(oxo) == 1 and len(g.bridging_metals(oxo[0][0])) == 4
    assert sum(1 for f in frags if g.fragment_bridge_class(f).value == "mu2") == 6


# ══ 5. the stated blockers are real ══════════════════════════════════════════

def test_the_hydroxide_bridge_is_blocked_by_b13():
    """A bare hydroxide perceives zero donors, so the µ2-OH row cannot build. Checked
    rather than cited, so the battery cannot claim to be waiting on a fixed defect."""
    assert perceive(embed_molecule(mol_from_smiles("[OH-]"), seed=7)) == []
    assert "B13" in case("cu2_mu2_hydroxide").blocked_on


def test_the_pyrazolate_bridge_is_blocked_by_b14():
    """Pyrazolate's two nitrogens are equivalent by resonance and type differently, and
    `pyridyl_N` declares `mono` only — so the bridge is refused on one end whatever the
    other offers, and mechanism A has no non-carboxylate test until it is fixed."""
    types = {s.donor_type for s in perceive(embed_molecule(mol_from_smiles("c1cc[n-]n1"),
                                                           seed=7))}
    assert len(types) > 1, "pyrazolate now types its two N alike — B14 is fixed"
    assert "B14" in case("cu2_pyrazolate").blocked_on


def test_a_blocked_row_carries_no_measured_numbers_it_cannot_have():
    """Nothing has been built for a blocked row, so nothing can have been measured off one.
    A `d_mm_bridge` on a row that cannot bridge would be a number someone made up."""
    for c in CASES:
        if c.blocked_on:
            assert c.d_mm_bridge is None, f"{c.name} is blocked but carries a measurement"


# ══ 6. what qc() cannot see today, and the shape of the verdict ══════════════

def test_qc_today_cannot_see_a_node_whose_centres_are_far_too_close():
    """Why S6 exists, demonstrated rather than asserted.

    Two Cu at 1.90 Å is roughly 0.7 Å inside any real paddlewheel and would tear apart
    under a relax. `check_clashes` allows the pair 1.74 Å and `check_metal_bonds` only looks
    at metal–DONOR pairs, so the existing report comes back clean.
    """
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.90]])
    assert clash_limit("Cu", "Cu") < 1.90            # below the clash floor, so invisible
    report = qc(["Cu", "Cu"], coords, bonded=set(), metal_idx=0, donor_idxs=[])
    assert report.ok, "qc already catches this — rewrite the S6 tests"


def test_a_qc_report_carries_centre_centre_findings_in_the_bad_bond_shape():
    """S6 asks for `BadBond`'s shape one level up, and for a wrong M···M to disqualify a
    report from `marginal` for `BadBond`'s reason: the node was BUILT wrong, and relaxing
    does not recover the geometry that was asked for."""
    from mofsbu.geometry.qc import BadIntercentre, QCReport

    bad = BadIntercentre((0, 1), 1.90, 2.62, 0.13, ("Cu", "Cu"), "cu_paddlewheel")
    report = QCReport(ok=False, bad_intercentres=[bad])
    assert report.code == "qc_intercentre"
    assert not report.marginal
    assert "Cu" in report.elements_involved()
    assert "target 2.62" in str(report) and "cu_paddlewheel" in str(report)


def test_the_multicentre_interfaces_exist_and_all_raise():
    """Ground rule 8, with real arguments rather than the junk `test_jobs.py` passes —
    which is the case where a stub that had quietly learned to answer would not be caught."""
    centres, constraints = fe3_oxo_skeleton()
    for call in (lambda: place_multicentre(centres, [], constraints),
                 lambda: check_intercentre(np.zeros((3, 3)), constraints,
                                           metal_idxs=[0, 1, 2])):
        with pytest.raises(NotBuiltYet, match="M6"):
            call()


# ══ 7. S3 — reconciliation, the one thing place_multicentre is for ═══════════

def oxo_cluster_skeleton(name: str):
    """An oxo-centred cluster as centres + constraints.

    The bridging atom is a `Center` like any other (§4) — that is the whole reframing, and
    it is why `Center.element` is not called `metal`.
    """
    c = case(name)
    centres = [Center(c.metal, cn=c.cn, local_geometry=c.local_geometry,
                      oxidation_state=3 if c.metal == "Fe" else 2,
                      spin_class="hs" if c.metal == "Fe" else "ls")
               for _ in range(c.n_metals)]
    centres.append(Center("O", cn=c.n_metals, local_geometry=c.mu_geometry, charge=-2))
    pairs = [(i, j) for i in range(c.n_metals) for j in range(i + 1, c.n_metals)]
    constraints = [InterCentreConstraint(p, mm_distance=c.d_mm, mm_lo=c.d_mm_lo,
                                         mm_hi=c.d_mm_hi, metal_metal_bond=c.mm_bond)
                   for p in pairs]
    return centres, constraints


def fe3_oxo_skeleton():
    return oxo_cluster_skeleton("fe3_mu3_oxo_333")


@m6
@pytest.mark.parametrize("name", ["fe3_mu3_oxo_333", "zn4o"])
def test_a_determined_skeleton_is_built_rather_than_searched(name):
    """S3's first half. Two determinants, one answer each, and the reconciliation is
    explicit policy — so the metals come out on the bridging centre's own polyhedron and
    the coordinates are reproducible, not the output of a minimiser."""
    c = case(name)
    centres, constraints = oxo_cluster_skeleton(name)
    result = place_multicentre(centres, [], constraints)
    metals = [i for i, s in enumerate(result.symbols) if s == c.metal]
    assert len(metals) == c.n_metals
    edges = [float(np.linalg.norm(result.coords[a] - result.coords[b]))
             for i, a in enumerate(metals) for b in metals[i + 1:]]
    assert max(edges) - min(edges) < 1e-6, "a symmetric skeleton has one edge length"
    assert abs(edges[0] - c.mm_from_bridging_centre()) < 0.01


@m6
@pytest.mark.parametrize("name", ["fe3_mu3_oxo_333", "zn4o"])
def test_the_residual_is_reported_as_strain_and_never_averaged(name):
    """§5 and the risk register's third row. The rigid model is wrong here by a measured
    ~0.5 Å; the honest response is to report it and let relaxation close it, never to
    split the difference and store a number neither determinant asked for."""
    c = case(name)
    result = place_multicentre(*[x for x in oxo_cluster_skeleton(name)][:1],
                              [], oxo_cluster_skeleton(name)[1])
    strain = result.choice_vector.get("strain") or getattr(result, "strain", None)
    assert strain is not None, "the disagreement has to be reported, not absorbed"
    assert abs(max(strain.values()) - c.collision()) < 0.05


@m6
def test_an_under_determined_skeleton_raises_rather_than_minimising():
    """S3's exit condition, and ground rule 5. Two centres and no distance between them is
    not a node with a default separation — it is a question nobody answered."""
    from mofsbu._types import AmbiguousSpecError

    centres = [Center("Cu", cn=5, local_geometry="square_pyramidal", oxidation_state=2,
                      spin_class="hs") for _ in range(2)]
    with pytest.raises(AmbiguousSpecError):
        place_multicentre(centres, [], [InterCentreConstraint((0, 1))])


@m6
def test_a_centre_with_no_local_geometry_raises_rather_than_defaulting():
    """There is no answer to "what shape is a CN-5 Cu" a placer may pick on its own: both
    CN-5 polyhedra are real and they build different molecules."""
    from mofsbu._types import AmbiguousSpecError

    centres = [Center("Cu", cn=5, local_geometry="", oxidation_state=2, spin_class="hs")
               for _ in range(2)]
    with pytest.raises(AmbiguousSpecError):
        place_multicentre(centres, [], [InterCentreConstraint((0, 1), mm_distance=2.62)])


@m6
def test_a_constraint_naming_a_centre_that_does_not_exist_is_refused():
    centres, _ = fe3_oxo_skeleton()
    with pytest.raises((ValueError, IndexError)):
        place_multicentre(centres, [], [InterCentreConstraint((0, 99), mm_distance=3.29)])


# ══ 8. S6 — check_intercentre on its own, so it can land before S3 ═══════════
#
# These call nothing but the QC seam, so they go green as soon as `check_intercentre` has
# a body — before any placer does.  Coordinates are written by hand for exactly that
# reason: a QC test that needs a working placer to produce its input cannot fail
# independently of it, and then it is not really testing QC.

def _two_centres(d: float) -> np.ndarray:
    return np.array([[0.0, 0.0, 0.0], [0.0, 0.0, d]])


def _paddlewheel_constraint() -> InterCentreConstraint:
    c = case("cu_paddlewheel")
    return InterCentreConstraint((0, 1), mm_distance=c.d_mm, mm_lo=c.d_mm_lo,
                                 mm_hi=c.d_mm_hi, metal_metal_bond=True)


@m6
def test_check_intercentre_accepts_a_node_at_its_literature_distance():
    found = check_intercentre(_two_centres(case("cu_paddlewheel").d_mm),
                              [_paddlewheel_constraint()],
                              metal_idxs=[0, 1], symbols=["Cu", "Cu"])
    assert found == []


@m6
@pytest.mark.parametrize("d", [1.90, 3.60], ids=["squeezed", "stretched"])
def test_check_intercentre_refuses_a_node_outside_its_window(d):
    """Both directions, because only one of them is ever a clash. A squeezed node might
    eventually trip the clash check with a big enough metal; a stretched one never will,
    and it is the one that quietly produces two mononuclear complexes sharing a record."""
    c = case("cu_paddlewheel")
    found = check_intercentre(_two_centres(d), [_paddlewheel_constraint()],
                              metal_idxs=[0, 1], symbols=["Cu", "Cu"])
    assert len(found) == 1
    assert found[0].centres == (0, 1)
    assert abs(found[0].distance - d) < 1e-6
    lo, hi = found[0].window
    assert lo <= c.d_mm <= hi and found[0].source


@m6
def test_check_intercentre_falls_back_to_a_tolerance_when_no_window_is_given():
    """A constraint may carry only an ideal — `mm_lo`/`mm_hi` are the curated refinement of
    it, not a requirement. Absent is not zero: a missing window must not read as a
    zero-width one that refuses everything."""
    found = check_intercentre(_two_centres(2.62),
                              [InterCentreConstraint((0, 1), mm_distance=2.62)],
                              metal_idxs=[0, 1], symbols=["Cu", "Cu"])
    assert found == []
