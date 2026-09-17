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

## 6. The battery, written before the code

M5's cheapest move was fixtures before code. M6's battery does not exist yet: `examples.ALL`
holds 16 entries and **Zn₄O is not among them**, nor are mononuclear Zn/BDC, BDC, bipy, EDTA,
anthrarufin, or the cis/trans anthrarufin–Cu pair that `PLAN_implementation.md` §6 lists as the
shared fixture set. Building the battery is the first half of S0, and each row is a failing
target first.

| cluster | n | bridges | what it stresses | today |
|---|---|---|---|---|
| Cu₂(µ-O₂CH)₄ paddlewheel | 2 | A ×4 + M–M | the baseline | builds, QC ok, L0+L1 ✅ |
| Cu₂(µ-O₂CH)₂ | 2 | A ×2 | under-bridged dimer; route discrimination | `zn2_bridged` is the Zn analogue |
| Cu₂(µ-O₂C–Ph)₄ | 2 | A ×4 | ligand-dependence of the span (2.46 vs 2.67 Å) | benzoate perceives ✅ |
| Cu₂(µ-pyrazolate)₂ | 2 | A, N,N | mechanism A beyond carboxylate | N typing wrong (§7) |
| Cu₂(µ-OH)₂ | 2 | B ×2 | bent bridging centre | OH⁻ unperceived; no bent CN-2 |
| Fe₃(µ₃-O)(µ-O₂CH)₆ | 3 | B + A ×6 | both mechanisms; 0.60 Å collision | fixture ✅, skeleton exact |
| Fe₃ mixed-valence (2,3,3) | 3 | as above | per-centre labels in L0 and L1 | fixture ✅ |
| Zn₄O(O₂CH)₆ | 4 | B(µ4) + A ×6 | nuclearity 4; 0.49 Å collision | fixture missing |

Between them: both mechanisms, nuclearity 2/3/4, with and without an M–M bond, carboxylate and
non-carboxylate, symmetric and mixed-valence.

---

## 7. What is missing, and what is wrong

| File | State |
|---|---|
| `sites/model.py::perceive` | stores well 0 only — mechanism A is unreachable (§3) |
| `sites/frames.py::torsion_wells` | no `BRIDGE_MU2` branch; falls through to the monodentate `(0.0, 180.0)` |
| `sites/model.py::vacancy_sites` | a vacancy offers `mono` only, so `compatible` refuses `mu2` on the metal side |
| `assembly/join.py::chelate_compatible` | refuses the two-different-metals case **by name** — that refusal is the new `bridge_compatible` |
| `assembly/join.py::join_chelate` | hardcodes `vacancies[0].atom_idx` for both bonds, i.e. one metal. The µ2 bridge is the same body with two |
| `assembly/join.py::join` | raises for `n_metals > 1` |
| `geometry/placer.py::place_mononuclear` | fills vertices in its own order, so the caller cannot say which to leave open — a CN-6 centre with four co-ligands comes back with its two vacancies **trans**, and a ~90° chelate cannot reach them (`chelate_cannot_span`). Declared placer work by the pathway ladder that hit it |
| `assembly/join.py::compatible` | refuses vacancy↔vacancy, naming M6 as what will place it |
| `geometry/placer.py::place_multicentre` | stub; `Center` and `Join` are named in its signature and **do not exist** |
| `geometry/placer.py::GEOMETRIES` | no bent CN-2 |
| `geometry/distances.py` | no M–M distances anywhere in `src/` or `data/reference/` |
| `geometry/qc.py` | no `check_intercentre` |
| `geometry/placer.py::to_rdkit` | hardcodes one metal at index 0; never writes an M–M bond |
| `sites/perception.py` | `[OH-]` perceives **zero** donors — filed as [B13](BUGS.md#b13) |
| `sites/perception.py` | pyrazolate's two equivalent N type differently — filed as [B14](BUGS.md#b14) |
| `assembly/construct.py` | `metal_block`'s metal charge gives one species two L1 — filed as [B12](BUGS.md#b12) |
| `examples.py` | Zn₄O absent; `PLAN_implementation.md` §6's fixture list was substantially aspirational (now corrected) |

Three of those are defects rather than unbuilt scope, so they are in `BUGS.md` and outlive this
file. B14 is the same class that `cff47a8` fixed for oxo-acids — perception is
resonance-invariant, an oxo-acid is one donor group — applied to azolates, which that change did
not cover.

---

## 8. Slice 0 — closeout, the battery, and five calls · **M**

**(a) M5's bookkeeping, which never happened.** `WORKPLAN_M5.md` §6 lists it: ledger D-numbers
for **C2** (θ_geom = 0.15 Å), for the L2-wiring/backfill call, and for `MAX_BITE_MISMATCH_DEG`,
each with a changelog line; `ALGO_VERSIONS["l3_conformer_id"]` off `"0-stub"`; the remaining
`CODE_ARCHITECTURE.md` rows flipped; `WORKPLAN_M5.md` deleted. M6 adds its own ledger entries,
so the ledger has to be current before it does.

**(b) The battery** (§6) and the perception fixes it depends on — [B13](BUGS.md#b13) for the
hydroxide bridge, [B14](BUGS.md#b14) for the pyrazolate one — as failing tests.

**(c) D20 — emergence primary, constraints for reconciliation.** Supersedes §6.1's framing.
M···M is an **output to validate** wherever one mechanism determines it (§3), and
`place_multicentre` is retained for the clusters where two determinants collide (§5) — a
measured reason rather than a blanket "headline cost".

**(d) C9 — what multiplicity does a polynuclear node carry?** The two ground-truth fixtures
disagree: `cu_paddlewheel` declares 1 (AF-coupled d⁹–d⁹) where `combined_multiplicity` — which
`join` uses — gives 3; `fe3_mu3_oxo` declares 16 and the additive rule agrees. Multiplicity is
in L0, so this decides identity. *Leaning:* coupling is a Kind-C branch, not a derivation —
whether two d⁹ centres give a singlet or a triplet is not recoverable from the centres, so the
multicentre path should **require** a stated multiplicity and raise otherwise. That is already
what `from_rdkit` does, and it makes both fixtures correct instead of one of them wrong.

**(e) C10 — when is there an M–M edge?** At 2.673 Å two Cu are bonded; at 5.516 Å they are not.
`EdgeType.METAL_METAL` is in the certificate, so a wrong answer is a wrong identity, and a
distance threshold applied silently is exactly the kind of inference ground rule 5 forbids.
Branch or declare; never default.

**(f) [B12](BUGS.md#b12) — `metal_block` labels metals differently from every other producer.**
`from_rdkit` and `examples.py` set a metal's `formal_charge` to 0 (D15, charge is graph-level);
`metal_block` writes `formal_charge=oxidation_state`, and that string feeds `NodeLabel.key()`.
Same species, two L1 hashes. Settle it before the gate, because the gate compares against
fixtures written the other way.

**(g) The plan-B trigger,** which `PLAN_implementation.md` §4 and §7 require be written down
when M6 starts: **2026-10-14**. Kept because the plan demands a date rather than a mood, while
recording that §3 and §4 have largely retired the risk it guards — the fallback was "templates
instead of a constrained placer", and the constrained placer turned out not to be on the
critical path.

---

## 9. The slices

### S1 — mechanism A: the well branch and the two-point bridge join · **L**

The lone-pair well becomes a Kind-B branch carried on the site or the verdict, for a
`TORSION_LIVE` donor in a bridging mode, following `_convergence`'s existing 4-way precedent.
`torsion_wells(_, BRIDGE_MU2)` gets its own entry instead of falling through to the monodentate
pair. `chelate_compatible`'s named refusal becomes `bridge_compatible`, keeping the verdict
shape.

**The two-point join is not invented here, it is generalised.** `join_chelate` already places
one ligand across two vertices with one rigid move, absorbing the residual into the bite angle
rather than the M–D bonds; it merely binds both donors to `vacancies[0].atom_idx`. A µ2 bridge
is the same body with the two vertices on different metals. So the slice is: relax that
assumption, and decide what the verdict compares when the two vertices no longer share an
origin — for a chelate it is the bite against the vertex separation *angle*; for a bridge the
vertices have no common centre, so it is the donor–donor **distance** against the
vertex-to-vertex distance.

*Exit:* the three bridging modes are separately enumerable and reproduce 2.67 / 5.15 / 5.52 Å;
the paddlewheel and its benzoate analogue build QC-clean; a bite that cannot span the vertices
refuses **by name**; replay from the emitted vector is bit-identical; `join_chelate`'s own
tests still pass unchanged, because a chelate is the case where the two vertices happen to
share a metal.

### S2 — mechanism B: the bridging atom as a centre · **M**

A bridging atom carries its own local geometry and vertices that are metal positions. Add the
bent CN-2 geometry `GEOMETRIES` lacks. Perception must first return a bare hydroxide at all.

*Exit:* µ3 and µ4 skeletons reproduce **3.291** and **3.168 Å**; a µ2-OH dimer builds at ~100°
and ~3.0 Å; a bridging aqua stops reporting an implied M···M of zero.

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
`metal_donor_distance` exactly — curated table, covalent-radii fallback that marks itself
`estimated` and says so in `source`.

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

### S5 — lift the `n_metals` guard, and discriminate the routes · **S**

`join` raises for `n_metals > 1` because the author assumed multi-centre needed a constraint
solve. Under D20 it does not. With the guard gone both routes run, and the model says which
works:

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

- `ALGO_VERSIONS`: `perception` bumps (the well is now part of what a site records);
  `placement` bumps if `join`'s orientation model changes. Old rows keep their own version's
  answer and are never re-labelled (ground rule 6 / D19).
- Design doc: D-numbers for **D20**, **C9** and **C10**, with a changelog line each. §6.1's
  constraint framing is superseded, not deleted — move it to `archive/DESIGN_history.md`.
- `CODE_ARCHITECTURE.md`: flip `place_multicentre`'s 🔴 row; add the two bridge mechanisms to §4
  if either becomes an invariant.
- `BUGS.md`: close B12; move the hydroxide and pyrazolate perception entries if they were filed
  there.
- `PLAN_implementation.md`: M6 moves to `archive/PLAN_completed.md`; §6's fixture list is
  corrected to what `examples.ALL` actually holds.
- Delete this file.

---

## 12. Order, at a glance

```
S0 battery + decisions ─┬─► S1 mechanism A ──┬─► S3 reconciliation ─► S7 gates ─► S8 persist
                        ├─► S2 mechanism B ──┘         ▲        ▲
                        └─► S4 nucleus ────────────────┘        │
                                                   S5 routes · S6 qc
```

S1 and S2 are independently testable and can be done in either order; S3 needs both, because
reconciliation is by definition what happens where they meet. S5 is small and can land as soon
as S1 does — it is what makes the failed route *say so*, which is worth having early.
