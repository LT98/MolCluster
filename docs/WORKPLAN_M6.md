# M6 work plan — polynuclear nodes that emerge from joins

`PLAN_implementation.md` §M6 says *what* the milestone is and what its exit gates are. This
file is the *order of work*: which slice lands first, what each one can be tested against on
its own, and the decisions that have to be made before any code is written. It is scratch —
when M6 lands, its conclusions go into the design doc's Decision Ledger and the plan's
changelog, and this file is deleted.

**The one-sentence scope:** after M6 a polynuclear node is something the assembly path
*reaches* — sequentially, or from a declared nucleus — across both bridge mechanisms and a
range of nuclearities, and a route that cannot form one is **refused with the number that
refused it**.

**What changed about the milestone, before any of it is built.** DESIGN §6.1 describes M6 as
generalising the placer to "a set of coordination centers with inter-center constraints" and
calls it "the headline engineering cost of the whole plan". Measured against a battery of real
motifs, that is not where the cost is. A Cu paddlewheel builds QC-clean from machinery that
already exists, and hashes to its M2 fixture at both L0 and L1, with no solver anywhere. What
the battery *does* turn up is that there are **two bridge mechanisms and only one of them is
supported**, and that for oxo-centred clusters two independent determinants fix the same M···M
distance and disagree by about half an angstrom. That disagreement is M6's real problem, and it
is one scalar per edge rather than a general constrained optimisation.

---

## 1. The two mechanisms

| | **A — multi-atom bridge** | **B — single-atom bridge** |
|---|---|---|
| shape | two donor atoms, one per metal | one donor atom, 2–4 metals |
| examples | carboxylate O,O · pyrazolate N,N | µ2-OH · µ2-O · µ3-O · µ4-O |
| geometry comes from | the two donors' lone-pair axes | the bridging atom's own local geometry |
| needed by | paddlewheel, Cu₂(µ-O₂CR)₂ | Fe₃-µ₃-oxo, Zn₄O, hydroxo dimers |
| state today | reachable, once the well branch exists (§3) | **not expressible at all** (§4) |

Drawing this line first is the point of the file. The two mechanisms share a verdict shape and
a QC path and share almost no geometry, and a milestone that treats them as one thing will
build A, declare victory on the paddlewheel, and discover B when it reaches the trimer.

---

## 2. What M6 inherits (do not rebuild any of this)

| Piece | Where | What it gives M6 |
|---|---|---|
| `site_frame(..., well=)` | `sites/frames.py` | an sp2 donor's **two** in-plane lone pairs — which is the whole of mechanism A's geometry (§3) |
| `_convergence`'s `product((0,1), repeat=2)` | `sites/model.py` | the precedent for branching over both donors' wells; it already does this for chelate pockets |
| `kabsch` (with its reflection guard) | `geometry/_linalg.py` | the two-point fit that places a bridge across two centres |
| `site_vectors` / `GEOMETRIES` | `geometry/placer.py` | local geometries — and they are generic over *which atom* sits at the centre, which is mechanism B (§4) |
| `vacancy_sites(idx, origin, directions)` | `sites/model.py` | already takes a bare index, so a **ligand** atom can carry vertices |
| the denticity-2 branch of `place_mononuclear` | `geometry/placer.py:526` | law-of-cosines two-point alignment, for donors at two different distances |
| **`join_chelate` / `chelate_reach`** | `assembly/join.py` | **the two-point join already exists** — one ligand across two vertices of *one* metal, with the mismatch absorbed into the bite angle rather than the M–D bonds. M6's µ2 bridge is this operation with the two vertices on *different* metals, and the function says so itself: *"one donor across two METALS is a mu2 bridge and is M6."* |
| `chelate_compatible`'s verdict shape | `assembly/join.py` | strain = mismatch / tolerance, `feasible ⇔ strain ≤ 1`, one policy number in one named constant |
| `Distance` with `source` / `estimated` | `geometry/distances.py` | a number that carries where it came from — the pattern M–M distances must follow |
| `BranchTree.refusals` | `assembly/construct.py` | how a failed route says why, instead of being silently absent |
| `inherit_sites` / `merge_inherited` keyed on `(atom_idx, slot)` | `sites/inherit.py` | site inheritance that already survives a metal carrying several vertices |

---

## 3. Mechanism A is one branch away

`site_frame` already knows a carboxylate oxygen has two in-plane lone pairs. Measured over
formate, the four well combinations are the three textbook bridging modes:

| wells | implied M···M | axis angle | what it is |
|---|---|---|---|
| (0,0) | 5.516 Å | 109° | anti-anti |
| (0,1) / (1,0) | 5.148 Å | 131° | syn-anti |
| **(1,1)** | **2.673 Å** | 11° | **syn-syn** — literature paddlewheel Cu···Cu is 2.62 |

**`sites.model.perceive` calls `site_frame` without `well=`**, so every stored `Site` carries
well 0 and only the anti-anti bridge is reachable. That single omission is why no polynuclear
node has ever emerged from the assembly path.

The span is a property of the ligand, not of the mode: benzoate offers **2.461 Å** where
formate offers 2.673.

**A whole paddlewheel builds from this and nothing else** — well-1 frames, `kabsch` for the
two-point fit, the existing distance table, the existing `qc`:

```
Cu...Cu 2.673 A  (lit. 2.62)
Cu-O    min 1.980  max 1.980   (target 1.98)
O-Cu-O basal angles  [89.5, 89.5, 89.5, 89.5]   (square-pyramidal wants 90)
QC: ok
L1 match vs examples.cu_paddlewheel(): True   |   L0 match: True
```

So the M6 exit gate's hard half is already answered: **the graph extraction and identity layers
agree with the fixtures**, and a red gate is the placer or its wiring, never identity.

---

## 4. Mechanism B: the bridging atom is a centre

A µ3-oxo has no covalent neighbours at all, and a bridging aqua has two — so neither has a
lone-pair axis to point at a second metal. `site_frame` returns `mode="aligned"` for water (two
neighbours determine one direction), the `well` argument is ignored, and both wells give the
**same** axis: implied M···M of **0.000 Å**. A µ2-aqua bridge is geometrically unbuildable
today even though `aqua_O` declares `mu2`.

The reframing that fixes it: **a bridging atom is a centre whose vertices are metal positions.**
`site_vectors` and `vacancy_sites` are already generic over the index, so this is mostly a
change of viewpoint — and it lands the skeletons exactly:

| bridge | local geometry | M–O–M | M···M | literature |
|---|---|---|---|---|
| µ3-O | trigonal | 120.0° | **3.291 Å** | 3.29 ✅ |
| µ4-O | tetrahedral | 109.5° | **3.168 Å** | 3.17 ✅ |
| µ2-OH | — | 180.0° | 3.900 Å | ~3.0 at ~100° ❌ |

The µ2 row fails for a reason worth stating plainly: **`GEOMETRIES` has only `linear` at CN 2**,
so a bent bridging oxygen cannot be expressed. Same class of defect as the legacy table mapping
CN 5 to `'planar'`, which this package already fixed once.

---

## 5. Where the two mechanisms collide — the real M6 problem

In a paddlewheel the bridges alone fix M···M and nothing contradicts them. In an oxo-centred
cluster the central atom fixes it too, and the two answers disagree:

| cluster | edge required by the oxo | offered by syn-syn formate | mismatch |
|---|---|---|---|
| Fe₃(µ₃-O)(µ-O₂CH)₆ | 3.291 Å | 2.696 Å | **0.595 Å** |
| Zn₄O(O₂CH)₆ | 3.168 Å | 2.677 Å | **0.491 Å** |

Real structures resolve this in the ligand: the O–C–O angle opens and the M–O–C angles adjust
until the carboxylate spans the edge the oxo dictates. A **rigid** ligand cannot do either, so
the model is quantifiably wrong here by about half an angstrom, and the honest response is to
report it as strain and let relaxation close it — not to average the two and store a number
neither determinant asked for.

This is what `place_multicentre` is for, and it is the measured reason to keep it: the
paddlewheel never needs it, and Fe₃-oxo and Zn₄O cannot be built without it. Note what it is
*not* — a general solver. Each disagreement is one scalar on one edge.

---

## 6. The battery — partly built

M5's cheapest move was fixtures before code, and a parallel session did that half already:
**`data/reference/node_cases.tsv` is the interface for adding an M6 case**, read by
`tests/node_cases.py` and checked by `tests/test_m6_battery.py`, with `examples.zn4o` and its
golden hashes landed alongside. Add a row and the whole set is re-checked; nothing in the test
module needs changing. Same arrangement as `ligand_cases.tsv` one layer down.

Two things about that table to keep in view. Its numbers are **typical of the named compound
class, not a refinement of one deposited structure** — no CIF was consulted — which is why
every row carries a window rather than only an ideal, and why the tests treat it as a shape
check. Before any of them backs a quantitative claim, check the ideal against the CSD and
narrow the window in the same commit. And three rows (Cr/Rh/Mo paddlewheels) carry
`fixture = -` deliberately, so the M–M target is not a two-point table; Rh and Mo are there
precisely because `geometry.distances.BASE_MO` does not know them, which is the gap S4 fills.

The exit-gate tests for the nodes are already written and **strict-xfailed on `NotBuiltYet`**,
so they flip to passing the moment a body lands and the strict marker forces the marker's
removal.

**Every row is now in the table** — S0 closed the two that were owed. The two that cannot carry
numbers carry `blocked_on` instead, and the blockers are *asserted*, not cited, so a fixed
defect cannot leave a row waiting for it.

| cluster | n | bridges | what it stresses | today |
|---|---|---|---|---|
| Cu₂(µ-O₂CH)₄ paddlewheel | 2 | A ×4 + M–M | the baseline | fixture + row ✅, builds QC-clean, L0+L1 ✅ |
| Fe₃(µ₃-O)(µ-O₂CH)₆ | 3 | B + A ×6 | both mechanisms; 0.60 Å collision | fixture + row ✅, skeleton exact |
| Fe₃ mixed-valence (2,3,3) | 3 | as above | per-centre labels in L0 and L1 | fixture + row ✅ |
| Zn₄O(O₂CH)₆ | 4 | B(µ4) + A ×6 | nuclearity 4; 0.49 Å collision; no vacancy at all | fixture + row ✅ |
| Cr₂ / Rh₂ / Mo₂ paddlewheels | 2 | A ×4 + M–M | M–M is not a two-point table; Rh and Mo are absent from `BASE_MO` | rows ✅, no fixture by design |
| Zn₂(µ-O₂CH)₂ | 2 | A ×2 | under-bridged dimer; route discrimination | fixture + row ✅ (`zn2_bridged_formate`); `d_mm` deliberately absent — nobody has measured it, and absent is not zero |
| Cu₂(µ-O₂C–Ph)₄ | 2 | A ×4 | ligand-dependence of the span (2.46 vs 2.67 Å) | row ✅, no fixture; benzoate perceives ✅ |
| Cu₂(µ-pyrazolate)₂ | 2 | A, N,N | mechanism A beyond carboxylate | row ✅, `blocked_on` [B14](BUGS.md#b14) |
| Cu₂(µ-OH)₂ | 2 | B ×2 | bent bridging centre | row ✅, `blocked_on` [B13](BUGS.md#b13) + no bent CN-2 |

Between them: both mechanisms, nuclearity 2/3/4, with and without an M–M bond, carboxylate and
non-carboxylate, symmetric and mixed-valence.

---

## 7. What is missing, and what is wrong

| File | State |
|---|---|
| ~~`sites/model.py::perceive`~~ | ✅ **landed (S1)** — every lone pair is stored, `frame` is still lobe 0 byte for byte, and the rest travel under `frame["lone_pairs"]` with no schema change. `perception` bumped to `3` |
| ~~`sites/frames.py::torsion_wells`~~ | ✅ **landed (S1)** — `_TWO_POINT_MODES` covers chelate and both bridges: a mode whose roll is an *output* of the fit must not branch over it, or it emits siblings with identical coordinates |
| `sites/model.py::vacancy_sites` | a vacancy offers `mono` only. **Not blocking** — `bridge_compatible` checks the donors' modes, as `chelate_compatible` does — but `compatible` still refuses `mu2` by blaming the vacancy instead of naming `bridge_compatible` |
| ~~`assembly/join.py::chelate_compatible`~~ | ✅ **landed (S1)** — `bridge_compatible` exists and refuses the same-metal case back, by name, pointing at `join_chelate`. Its verdict is a **distance**, and it is necessary-not-sufficient: qc stays the arbiter |
| ~~`assembly/join.py::join_chelate`~~ | ✅ **landed (S1)** — `join_bridge` is that body with the two vertices on different metals: one DATIVE per donor to its OWN metal, no bisector slide, both lobes recorded |
| ~~`assembly/join.py::join`~~ | ✅ **landed (S1)** — the `n_metals > 1` guard is gone under D20, and `lone_pair=` is a recorded coordinate. A bridged dimer now builds from two one-contact joins at 2.673 Å. `choice_vector` → `cv2` |
| `geometry/placer.py::place_mononuclear` | fills vertices in its own order, so the caller cannot say which to leave open — a CN-6 centre with four co-ligands comes back with its two vacancies **trans**, and a ~90° chelate cannot reach them (`chelate_cannot_span`). Declared placer work by the pathway ladder that hit it |
| `assembly/join.py::compatible` | refuses vacancy↔vacancy, naming M6 as what will place it |
| `geometry/placer.py::place_multicentre` | **signature settled, body raises.** `Center` and `InterCentreConstraint` (with `metal_metal_bond` and `window()`) exist; `Center.element` is deliberately not always a metal, because a bridging atom is a centre (§4). **`Join` stays named-but-undefined on purpose** — what the placer receives from the join path is S1's output and therefore S3's call, so guessing it now would settle the wrong end first |
| ~~`geometry/placer.py::GEOMETRIES`~~ | ✅ **landed (S2)** — `bent` at CN 2, taking `angle_deg` and refusing without it. CN 2 is now a Kind-C branch |
| `geometry/distances.py` | no M–M distances in `src/`; the *targets* are curated in `data/reference/node_cases.tsv`, and Rh/Mo are deliberately absent from `BASE_MO` |
| `geometry/qc.py::check_intercentre` | **signature settled, body raises** — `BadIntercentre` exists |
| `geometry/placer.py::to_rdkit` | hardcodes one metal at index 0 and never writes an M–M bond. The fix is a **widening of this function, not a second one** — one conversion means one place to be wrong, which is its own docstring's argument (S7) |
| ~~`sites/perception.py`~~ | ✅ **fixed (S2)** — `[OH-]` perceives `hydroxide_O`; [B13 archived](archive/BUGS_resolved.md). Turned up [B15](BUGS.md#b15): `[O-2]` types as a hydroxide |
| `sites/perception.py` | pyrazolate's two equivalent N type differently — filed as [B14](BUGS.md#b14) |
| ~~`assembly/construct.py`~~ | ✅ **fixed (S0)** — `metal_block` writes `formal_charge=0`; [B12 archived](archive/BUGS_resolved.md) |
| ~~`examples.py`~~ | ✅ Zn₄O landed; `examples.ALL` holds 17 and `PLAN_implementation.md` §6 now says so |

The two perception entries are defects rather than unbuilt scope, so they are in `BUGS.md` and
outlive this file. B14 is the same class that `cff47a8` fixed for oxo-acids — perception is
resonance-invariant, an oxo-acid is one donor group — applied to azolates, which that change did
not cover.

---

## 8. Slice 0 — closeout, the battery, and five calls · **M** · ✅ **LANDED**

**(a) M5's bookkeeping, which never happened.** ✅ Ledger entries **D21** (θ_geom = 0.15 Å,
carrying the calibration table so the measurement outlives the workplan that made it), **D22**
(the L2-wiring call, and the digest backfill folded in beside it — re-derive annotations,
version addresses) and **D23** (`MAX_BITE_MISMATCH_DEG`), each with a changelog line.
`ALGO_VERSIONS["l3_conformer_id"]` is `"conf1"`. M5 moved to `archive/PLAN_completed.md`;
`WORKPLAN_M5.md` deleted; `CODE_ARCHITECTURE.md` §6 now says what B2 actually is, which is not
what it said.

**B2 is sharper than "the stub is not filled in", and the sharpening is the useful part.** The
classifier is built and `store_block` derives a tag whenever it has coordinates — but `runner`
passes `l2=""` deliberately, because it is the only path that does not tag, and a path that
tagged alone would file its products away from the node every other route reaches. So B2 is a
*wiring* seam, not an unbuilt one, and closing it moves both paths together.

**(b) The battery** (§6). ✅ Complete, and wider than this slice asked for: all nine rows are in
`data/reference/node_cases.tsv`, including the two mechanism-A rows §6 listed as owed
(`cu2_benzoate_paddlewheel`, and `zn2_bridged_formate` as the under-bridged dimer) and the two
blocked ones. The blocked pair is checked rather than cited —
`test_the_hydroxide_bridge_is_blocked_by_b13` and its pyrazolate twin assert the *defect*, so
the battery cannot go on claiming to wait for something already fixed.

**(c) D20 — emergence primary, constraints for reconciliation.** ✅ Called. Supersedes §6.1's
framing, which is struck there and preserved in `archive/DESIGN_history.md`. M···M is an
**output to validate** wherever one mechanism determines it (§3), and `place_multicentre` is
retained for the clusters where two determinants collide (§5) — a measured reason rather than a
blanket "headline cost".

**(d) C9 — what multiplicity does a polynuclear node carry?** ✅ Called as **D24: stated, never
combined.** The two ground-truth fixtures disagree — `cu_paddlewheel` declares 1 (AF-coupled
d⁹–d⁹) where `combined_multiplicity` gives 3; `fe3_mu3_oxo` declares 16 and the additive rule
agrees — and the resolution is that *both are right*, because the additive rule has a
precondition the paddlewheel does not meet. Coupling is not recoverable from the centres, so
the multicentre path requires a stated multiplicity and raises otherwise, which is already what
`from_rdkit` does. Pinned by
`test_m6_battery.py::test_multiplicity_is_declared_because_coupling_is_not_derivable`, against
the golden L0 strings so a derived answer fails at the identity rather than at the arithmetic.

**(e) C10 — when is there an M–M edge?** ✅ Called as **D25: declared, never inferred.** At
2.673 Å two Cu are bonded; at 5.516 Å they are not. `EdgeType.METAL_METAL` is in the
certificate, so a wrong answer is a wrong identity, and a distance threshold applied silently is
exactly the kind of inference ground rule 5 forbids.

**Independently reached, which is the strongest evidence a gate is real.** The parallel
session's reference table arrived at the same place from the data side — `node_cases.tsv`
carries `mm_bond` as a column and says it is "a chemical decision, NOT derivable from d_mm",
and `InterCentreConstraint.metal_metal_bond` says a placer that guessed it from distance would
be guessing the identity of its product. The Fe₃ trimer at 3.29 Å has no edge and the Cu₂
paddlewheel at 2.62 Å has one; the two fixtures differ at L1 by exactly that, and
`test_the_table_and_the_fixture_graph_agree` holds the column and the graph to each other.

**(f) [B12](archive/BUGS_resolved.md) — `metal_block` labels metals differently from every
other producer.** ✅ Fixed: it writes `formal_charge=0` like `from_rdkit` and `examples.py`
(D15, charge is graph-level). **What the fix turned up is worth more than the fix.** 683 tests
passed before and after and no golden hash moved — every stored hash descends from
`examples.py`, and nothing had ever asserted on an identity produced by the enumerator. D2's
claim was being tested along one route. So it ships with a producer-agreement test pinning the
label string itself.

**(g) The plan-B trigger,** which `PLAN_implementation.md` §4 and §7 require be written down
when M6 starts: ✅ **2026-10-14**, recorded in §4's risk row. Kept because the plan demands a
date rather than a mood, while recording that §3 and §4 have largely retired the risk it guards
— the fallback was "templates instead of a constrained placer", and the constrained placer
turned out not to be on the critical path.

---

## 9. The slices

### S1 — mechanism A: the well branch and the two-point bridge join · **L** · *mostly landed*

**The lone-pair well, as perception (landed first).** `perceive` records every lobe rather
than only the first, so the anti lobe exists at all; `lone_pair_frames` answers "how many
directions does this donor offer" in one place, and its tuple length *is* the answer to "can
this atom bridge on its own" — a determined donor like aqua collapses to one lobe, which is
§4's problem stated as data. `torsion_wells(_, BRIDGE_MU2)` has its own entry, alongside
chelate, under a rule worth keeping: **a mode whose roll is an output of the fit must not
branch over it**, or the tree emits siblings whose coordinates are identical. `perception`
bumped to `3`, and a `2` catalog is not readable as "this donor has one lobe" — which is why
it is a bump and not a backfill.

**The lone-pair well, as a join coordinate (landed second) — and the measurement that
reframes the rest of the slice.** `join(..., lone_pair=k)` selects the lobe, `compatible`
reports how many are on offer, and the index travels in the choice vector so a replay
reproduces the lobe rather than the default. With it, **the bridged dimer is reachable through
two ordinary joins** — formate onto one Cu, its free oxygen onto a second — and all four lobe
combinations land exactly on the frame-implied numbers:

| lobes | built Cu···Cu | mode |
|---|---|---|
| (1,1) | **2.673 Å** | syn-syn — the paddlewheel, 0.05 Å from literature |
| (0,1) / (1,0) | 5.148 Å | syn-anti |
| (0,0) | 5.516 Å | anti-anti |

So `join`'s `n_metals > 1` guard is gone, which was **S5's job and turns out to belong here**:
the guard was written on the premise D20 replaced. One donor onto one vertex is *one contact*,
and one contact is satisfied by a rigid move of the donor's block whatever either block already
carries — so nothing about a second centre makes the pose undetermined, and the M···M is an
output. `choice_vector` bumped to `cv2`: a key over a larger set of coordinates is a different
key, and that is true even though lobe 0 replays every `cv1` path byte-identically.

**`bridge_compatible` and `join_bridge` (landed third).** One ligand across two vertices of
*different* metals in one move. The verdict keeps `chelate_compatible`'s shape but compares a
**distance** — the donors' own separation against the separation the two target points require
— because two vertices on two centres have no common origin to subtend an angle at.

**Three things the building turned up, none of them in the plan:**

* **The mismatch does not land where the chelate's does, and the tolerance had to be
  re-derived because of it.** A chelate slides along its bisector so the residual goes into the
  bite angle and the M–D bonds keep their length — deliberate, because bond lengths are what QC
  checks. A bridge has no bisector, so where the residual goes falls out of the two vertex
  axes. Measured: a **0.380 Å** span mismatch on a square-planar dimer produced Cu–O bonds of
  **1.989 Å against a 1.980 target** — 0.009 Å, not the 0.190 a per-bond split predicts. It
  went into the angles. So `MAX_BRIDGE_SPAN_MISMATCH_A` is QC's own bond tolerance used as a
  **bound** ("a mismatch bigger than the slack one bond gets has nowhere to go"), not as a
  derivation, and the constant says so.
* **The verdict is necessary and not sufficient.** On a CN-6 dimer the *smallest* mismatch in
  the whole candidate set — 0.183 Å, better than the square-pyramidal pair that builds cleanly
  — puts the ligand's carbon 1.41 Å from the far metal and its far oxygen 0.76 Å from it. A
  span test cannot see that, so a caller ranks on the verdict and then runs `qc`, which is
  exactly what `chelate_reach` says one layer down. Pinned as a test rather than a caveat.
* **The second bridge is not blocked after all — the earlier reading of the vertex-orientation
  measurement was wrong.** A second formate places across a sequentially-built dimer and the
  product is **QC-clean** at square-planar and square-pyramidal, two µ2 bridges on the graph.
  What the arbitrary vertex orientation costs is not feasibility but *quality and choice*: the
  usable pairs are accidents of where `site_vectors` happened to point, most cross-metal pairs
  are refused, and nothing arranges four bridges at 90° around the M···M axis. That is still
  S4's job; it is a worse paddlewheel rather than no second bridge.

**Still owed:** `vacancy_sites` offers `mono` only. `bridge_compatible` does not consult it —
it checks the donors' modes, as `chelate_compatible` does — so this is not blocking, but
`compatible` still refuses `mu2` with "the vacancy does not bind it" where it should say that a
bridge is two contacts on two metals and name `bridge_compatible`.

**The two-point join is not invented here, it is generalised.** `join_chelate` already places
one ligand across two vertices with one rigid move, absorbing the residual into the bite angle
rather than the M–D bonds; it merely binds both donors to `vacancies[0].atom_idx`. A µ2 bridge
is the same body with the two vertices on different metals. So the slice is: relax that
assumption, and decide what the verdict compares when the two vertices no longer share an
origin — for a chelate it is the bite against the vertex separation *angle*; for a bridge the
vertices have no common centre, so it is the donor–donor **distance** against the
vertex-to-vertex distance.

*Exit:* ✅ the three bridging modes are separately enumerable and reproduce 2.67 / 5.15 /
5.52 Å **through the assembly path**, not only from the frames; ✅ replay from the emitted
vector is bit-identical; ✅ `join_chelate`'s own tests pass unchanged, because a chelate is the
case where the two vertices happen to share a metal; ✅ a span that cannot reach the vertices
refuses **by name and with its number**. Still owed: the **whole** paddlewheel and its benzoate
analogue QC-clean — two bridges build, four need S4's vertex arrangement.

### S2 — mechanism B: the bridging atom as a centre · **M** · ✅ **LANDED**

Three pieces, and the third one is where the change of viewpoint actually bit.

**[B13](archive/BUGS_resolved.md) closed.** `_classify_anionic` opened with "anionic means
deprotonated", which is true of every donor it recovers except the one this slice needs:
hydroxide is anionic **and** still carries its proton, so the guard dropped it. The giveaway
was that the function already had an unreachable `hydroxide_O` branch. It also turned up a
separate defect — `[O-2]` types as `hydroxide_O` too, so an oxo reports a hydroxide's pKa —
filed as [B15](BUGS.md#b15) rather than folded in, because giving oxo its own type is a change
to `donor_descriptors.tsv`.

**`bent` added to `GEOMETRIES`, and it is the one entry whose name does not fix its vertices**,
so `site_vectors` takes `angle_deg` and refuses without it (ground rule 5). A single baked-in
number would have been exactly the tolerance-turned-folklore §10 warns about: a µ2-hydroxide
sits near 100° and a bent µ2-oxo well above it. CN 2 therefore becomes a Kind-C branch, as CN 4
already is — the difference between linear and bent is 0.9 Å of M···M.

**`bridging_metal_positions` reads the BONDING, not a `SiteFrame`, and that is the whole
distinction between mechanism A and mechanism B.** A frame answers "where does *one* metal go",
and for a single-neighbour donor its axis is a lone-pair lobe tilted ~120° off the bond. Using
that as a bridge's bisector swung one metal of a µ2-hydroxide to **1.68 Å from its own proton**.
A bridge's bisector is the direction away from everything the atom is already bonded to. With
that, all three exit numbers land:

| bridge | geometry | M–O–M | M···M | clearance to substituents |
|---|---|---|---|---|
| µ2-OH | bent | 100.0° | **2.988 Å** | 2.69 Å |
| µ2-aqua | bent | 104.5° | **3.084 Å** (was 0.000) | 2.48 Å |
| µ3-oxo | trigonal | 120.0° | **3.291 Å** | — |
| µ4-oxo | tetrahedral | 109.47° | **3.168 Å** | — |

A bare oxo has nothing to orient against and every azimuth is free — honest rather than sloppy,
since such a skeleton is fixed by its *angles* and every orientation gives the same M···M set.

**And a data defect the slice tripped over:** `cu2_mu2_hydroxide` had been one field SHORT in
`node_cases.tsv` since the day it was written, so `csv.DictReader` shifted every value after the
gap one column left — its note was being read as its `source`, and `note` was `None`. Nothing
caught it because nothing asserted on those two columns. `tests/node_cases.py` now refuses a row
whose width disagrees with the header.

*Exit:* ✅ µ3 and µ4 skeletons reproduce 3.291 and 3.168 Å; ✅ a µ2-OH bridge builds at 100° and
2.988 Å; ✅ a bridging aqua reports 3.084 Å instead of zero.

### S3 — reconciliation, and `place_multicentre` · **L**

`Center` and `InterCentreConstraint` — named in the stub's signature and never defined.
Skeletons that are *determined* are **built, not searched**: two centres and a distance; a µ₃
bridge as an equilateral triangle at M···M = √3·d. An under-determined set raises
`AmbiguousSpecError` naming what is missing, rather than falling into a minimiser. The
reconciliation is explicit policy — which determinant wins, and the residual reported as strain
— never a hidden average.

*Exit:* Fe₃-oxo and Zn₄O build with the residual reported at the measured ~0.5 Å; the
paddlewheel path never calls this; an under-determined constraint set refuses.

### S4 — the declared nucleus, and reserved vertices · **M**

Two halves of one idea: **the caller says what the starting geometry leaves open.**

Vacancy↔vacancy becomes a `METAL_METAL` join instead of a refusal. Its placement needs an M–M
distance, which is a legitimate **input** here because a nucleus is being declared:
`metal_metal_distance(m1, m2, *, motif) -> Distance` in `geometry/distances.py`, following
`metal_donor_distance`'s shape.

**Called, 2026-09-18: there is no curated M–M table, and the caller states the number.** The
two tempting sources both fail on inspection. `node_cases.tsv`'s Cr/Rh/Mo rows say in their own
notes that they exist for this table — but that file's header says its `d_mm` values are
*literature-typical for a compound class*, no CIF consulted, which is why every row carries a
window; they are what a built node is measured **against**, and driving a placer from them
would close the loop and make gate 4 compare a number with itself. A second curated TSV
duplicating them would only move the problem and add two files to keep in step.

So `metal_metal_distance` returns a covalent-radii estimate that marks itself `estimated`, and
**raises rather than inventing a motif-specific number** — a quadruply-bonded Mo₂ at 2.09 Å and
a Cu₂ paddlewheel at 2.62 Å are not the same question, and nothing in the elements distinguishes
them. That is consistent rather than restrictive: D20 says a distance is an input only where it
is *declared*, and a declaration comes from the caller, not from a table the caller did not
write.

And `place_mononuclear` learns to take **which vertices to reserve**. It currently fills them in
its own order, so a CN-6 centre carrying four co-ligands comes back with its two vacancies
*trans* and a ~90° chelate cannot reach them — the pathway ladder hit exactly this and refused
the step with `chelate_cannot_span`, correctly naming it placer work. A bridge needs the same
thing for the same reason: the next bridge wants a *cis* pair, and the vertex that faces the
partner metal is not available to anything else.

*Exit:* a Cu₂ block reports its remaining vacancies with frames in the dimer's own coordinate
system, and `enumerate_constructions` grows ligands onto it — prove this rather than assume it,
since B8's occlusion rule is exactly what could mark a dimer's vertices `BLOCKED`. A centre
asked to reserve a cis pair returns one, and the ladder's co-ligand-saturated series connects.

### S5 — ~~lift the `n_metals` guard~~, and discriminate the routes · **S**

**The guard is already gone** — it came out with S1, because that is the slice whose
measurement proved it was refusing a determined placement. What is left here is the second
half: with both routes running, the model says which works:

| route | result |
|---|---|
| `Cu₂(µ-HCOO) + HCOO` | closes at **2.673 Å**, 0.05 Å from literature |
| `Cu(HCOO)₂ + Cu` | bottoms out at **1.464 Å** — refused |

The second is a genuinely failed route, and the measurement says why: pinning both formates to
one Cu *before* the second exists forces their free oxygens to agree about a metal they were
never placed for. Best over a continuous scan of both rolls **and** the out-of-plane swing —
1.819 Å roll-only, 1.464 Å with the swing, on square-planar and octahedral alike; 3.273 Å
tetrahedral, which is right, since 109° vertices splay the formates further than a 90° pair.
Irreducible in a rigid model.

*Exit:* the failed route is refused **with its number**, through `BranchTree.refusals`, not
silently absent.

### S6 — `qc.check_intercentre` · **S**

Validate the M···M the joins produced against the table, plus bridge angles, in the `BadBond`
shape — measured, target, tol, source, `describe()`, `to_dict()`. A wrong M···M disqualifies a
report from `marginal` for `BadBond`'s reason: the node was *built* wrong, and relaxing does
not recover the geometry that was asked for.

### S7 — extraction, and the exit gates · **M**

Multi-metal `to_rdkit` — one DATIVE per join, SINGLE per M–M edge — then the existing
`from_rdkit`, so an assembled node is typed by the same path as a molecule read from SMILES.
`tests/test_m6_exit_gates.py`, written first and failing, across the battery rather than one
motif:

| # | Gate |
|---|---|
| 1 | sequential and nucleus-first routes reach **one** node with **two** provenance edges |
| 2 | paddlewheel L1 == `aed8d6a834418f38…`; Fe₃-oxo L1 == `9a3022deefc43dae…` |
| 3 | every battery row builds QC-clean, or is refused with a stated number |
| 4 | each cluster's M···M is within literature range — reported as a measurement, not asserted |
| 5 | an xTB relax keeps each node intact: re-extract afterwards, same L1 (skipped without tblite, as `scripts/regress_m7.py` does) |
| 6 | `Cu(HCOO)₂ + Cu` is refused and says 1.46 Å |

Gate 1 is the strong one. It is M5's "one node, two routes" generalised to polynuclear, and it
is M8's Path A vs Path B made buildable at M6.

### S8 — persist and the round trip · **S**

`assembly.persist.store_construction` for polynuclear nodes, proved **after** a round trip
through SQLite rather than in memory, the way `tests/test_m5_exit_gates.py` does.

---

## 10. Risks, with the escape hatch named in advance

| Risk | Early signal | Hatch |
|---|---|---|
| **L3 explosion from the well branch** | leaf counts far above S4's measured table (405 leaves over 5 distinct L1 at CN 6 degree 2) | the declared hatches already exist and S1 must show them still holding: choice-vector dedup before geometric clustering, the per-`(L1,L2)` cap, `TORSION_FREE` never branching |
| **A tolerance becomes folklore** | a bridging threshold picked from the figures in this file | those are *points*, not populations. Calibrate from the battery's distribution the way θ_geom was, and pin the separation rather than the threshold |
| **The rigid ligand model is wrong by ~0.5 Å on oxo-centred clusters** | Fe₃-oxo and Zn₄O QC-fail on bond lengths | it is a known, measured limitation (§5): report it as strain, close it by relaxation, and never average two determinants into a number neither asked for |
| **Identity retrofit pressure** | a proposal to add a flag to L1 so a bridge works | the M2 discrimination table is the contract. §3 shows the fixtures already match, so any pressure here is a bug elsewhere |
| **M6 overruns** | the battery still will not build | **2026-10-14** → `geometry/templates.py` (§8g) |

---

## 11. Bookkeeping when M6 lands

- ✅ `ALGO_VERSIONS`: `perception` → `3` (the lobes are part of what a site records) and
  `choice_vector` → `cv2` (a join records which lobe it bound), both in S1. `placement` does
  **not** move: lobe 0 reproduces every earlier geometry byte for byte, so the orientation
  model is unchanged and only the set of recorded coordinates grew. Old rows keep their own
  version's answer and are never re-labelled (ground rule 6 / D19).
- ✅ Design doc: **D20**, **D24** (C9) and **D25** (C10) written, with a changelog entry; §6.1's
  constraint framing struck there and preserved in `archive/DESIGN_history.md`.
- `CODE_ARCHITECTURE.md`: flip `place_multicentre`'s 🔴 row; add the two bridge mechanisms to §4
  if either becomes an invariant.
- ✅ `BUGS.md`: B12 closed and archived. Still owed: the hydroxide and pyrazolate perception
  entries (B13, B14) move when they are fixed.
- `PLAN_implementation.md`: M6 moves to `archive/PLAN_completed.md`. ✅ §6's fixture list is
  corrected to what `examples.ALL` actually holds.
- Delete this file.

---

## 12. Order, at a glance

```
S0 battery + decisions ─┬─► S1 mechanism A ──┬─► S3 reconciliation ─► S7 gates ─► S8 persist
         ✅             ├─► S2 mechanism B ──┘         ▲        ▲
                        │        ✅                    │        │
                        └─► S4 nucleus ────────────────┘        │
                                                   S5 routes · S6 qc
```

S0 and S2 are closed. S1 built the lone-pair branch, the join coordinate and `join_bridge`; what
it has not reached is the whole four-bridge paddlewheel, which needs S4's vertex arrangement —
measured, not assumed, so **S4 comes before the rest of S1**.

S1 and S2 are independently testable and can be done in either order; S3 needs both, because
reconciliation is by definition what happens where they meet. S5 is small and can land as soon
as S1 does — it is what makes the failed route *say so*, which is worth having early.
