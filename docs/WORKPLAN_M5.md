# M5 work plan — recursive assembly (the N=1 path of the general operation)

`PLAN_implementation.md` §M5 says *what* the milestone is and what its exit gates are.
This file is the *order of work*: which slice lands first, what each one can be tested
against on its own, and the two decisions that have to be made before any code is written.
It is scratch — when M5 lands, its conclusions go into the design doc's Decision Ledger and
the plan's changelog, and this file is deleted.

**The one-sentence scope:** after M5 you can take stored blocks, join them, and rebuild the
product exactly from its provenance — for products with **one metal centre**. Two centres is
where M6 begins, and that seam is drawn explicitly in slice 3 rather than discovered halfway
through.

---

## 1. What M5 inherits (do not rebuild any of this)

| Piece | Where | What it gives a join |
|---|---|---|
| `SiteFrame`, `live_dof`, `binding_modes`, `torsion_wells` | `sites/frames.py` | both partners' frames and the discrete well set — D13's whole point is that a join reads these rather than re-searching |
| `Site` with `role = donor \| vacancy`, `slot` | `sites/model.py` | one type for both ends of a bond, so `compatible(a, b)` is one question (landed in `bf9efe9`) |
| `BuildingBlock.open_sites / open_donors / open_vacancies` | `assembly/join.py` | the input set for any join; raises rather than guessing when state is absent |
| `inherit_sites`, `merge_inherited` | `sites/inherit.py` | the product's site list, without re-perceiving |
| `place_mononuclear(..., cn=...)` returning `vacancies` + `choice_vector` | `geometry/placer.py` | the N=1 geometry path **and** a working precedent for what a choice vector looks like |
| `choice_vector_digest` / `choice_vector_json` / `seed` columns; `reactions.choice_vector_digest`; `Provenance(atom_map=..., depth=...)` | `registry/schema.sql`, `registry/api.py` | storage is already shaped for this; M5 adds no migration for the core path |

## 2. What is missing

| File | State |
|---|---|
| ~~`assembly/choice.py`~~ | ✅ **landed (S1)** — `ChoiceVector`, canonical form, `cv1:` digest, JSON round trip |
| ~~`assembly/construct.py`~~ | ✅ **landed (S4)** — branch tree, replay, Kind-C refusals |
| `assembly/join.py::compatible` | ✅ **landed (S2)** — plus `chelate_compatible`, the two-constraint case |
| `assembly/join.py::join` | ✅ **landed (S3)** — frame alignment, atom maps, site inheritance, choice vector |
| `assembly/join.py::grow` | ✅ **landed (S4)** — a thin front on `enumerate_constructions` |
| `identity/keys.py::l2_isomer_tag` | ✅ **landed (S5)** — `identity/isomers.py`, version `"iso1"`; builder wiring deferred |
| `identity/keys.py::l3_conformer_id` | ✅ **landed (S6)** — `identity/conformers.py`; provenance-primary + clustering |
| C2 — θ_geom | ✅ **called (S6)** — 0.15 Å, calibrated on xTB-relaxed geometries |
| C2 — the energy window | open, deliberately: the measurement is contaminated by ISSUES 6c |
| S7 registry wiring | ✅ **landed** — `assembly/persist.py`; all three exit gates pass |

---

## 3. Slice 0 — two decisions, before any code · **S**

**(a) L2 backfill vs. version bump.** ~~The moment `l2_isomer_tag` returns a real tag,
structures that were one row become two.~~ **Reframed by S5, and no longer blocking.** The
premise was wrong: `put_structure` calls `l2_isomer_tag(g)` with no coordinates — it cannot
have any, since the structure row is created *before* its geometry — and nothing in `src/`
passes `l2=`. So building the classifier split nothing, and the corpus is untouched.

The decision survives in a smaller and later form: **when the builder is wired to pass its
coordinates, what happens to the rows that predate it.** Measured rather than guessed — of
the 39 stored structures (all `l2=""` at `algo_l2="0-stub"`), **18 would gain a real tag**
and 21 would still tag `""`. Ground rule 6 says old rows are never re-labelled, so the
default is that they stay as un-tagged ancestors distinguishable by `algo_l2`; re-tagging
them is possible (they all have geometries) but has to be an explicit migration with a
D-number. Decide it next to S4.1, where the wiring actually happens.

**(b) The choice-vector digest seam — a live bug, found while planning.** *(Writer fixed in
S1; the backfill decision below is still open.)*
`registry.api.put_geometry` writes `(choice_vector or {}).get("digest")`, but nothing in the
tree ever puts a `"digest"` key in a choice vector. The measured consequence:

```
geometries: 191 rows | 99 carry choice_vector_json | 0 carry choice_vector_digest
```

So `ix_geometries_choice` indexes nothing, and "find every geometry built from this choice"
— the query replay depends on — cannot run today. M5 owns this either way; decide now
whether the 99 existing JSON rows get their digests computed retroactively (they are
regenerable, so this is cheap and safe) or are left NULL as pre-M5 artifacts.

---

## 4. The slices

### S1 — `assembly/choice.py`: the vocabulary · **M** · ✅ **LANDED**

Everything downstream emits one of these, so it landed first.

**What shipped** — `ChoiceVector` (canonicalising, serialisable, hashable), `canonical`,
`canonical_json`, `digest_of`, `ChoiceVector.from_json` / `.coerce`. The key set was not
invented: `place_mononuclear` already emitted `{metal, geometry, cn, d_ml, seed,
ligands[{ligand, mode, torsion_well, azimuth_step, oop_step, …}], donor_distances}`, and
`azimuth_step`/`oop_step` are exactly the "stochastic search reduced to a recorded index"
pattern D13 asks for — so `ChoiceVector` wraps a mapping and every producer stays a
producer of dicts.

**Three things the recipe decides, each with a test:**

* **The seed is out of the key.** It is a Kind-A coordinate (§6.7) — noise to dedup, not a
  branch — and `ix_geometries_choice` exists to return the whole sampled family so
  clustering has something to collapse. Putting the seed in would make that family
  unfindable and every duplicate look like a distinct conformer. It is not discarded:
  `ChoiceVector.seed` lifts it, and `put_geometry` falls back to it for the column.
* **A self-referential `digest` key is out**, or the stored-JSON round trip is not a fixed
  point.
* **Floats quantise to 6 dp**, and the version prefixes the digest (`cv1:…`) rather than
  living in a column — a choice digest is not part of an address, so a stale one should be
  self-identifying. `ALGO_VERSIONS["choice_vector"]` is the pin.

**Correction to this plan as first written:** the exit test was stated as
`replay(digest_of(cv)) == cv`, which is not a thing — a digest is one-way. The round trip
that matters is through *storage*: `ChoiceVector.from_json(cv.to_json()).digest == cv.digest`.
Turning a vector back into coordinates is `construct`'s job (S4) and needs the row's seed as
well, which is the second reason the seed is a column.

**Slice 0(b) closed on the writer side.** `put_geometry` now digests the vector instead of
reading a `"digest"` key nobody wrote. Tests: `tests/test_choice.py` (19), with
`test_put_geometry_stores_a_digest` as the regression gate and
`test_the_same_choice_at_several_seeds_is_one_findable_family` covering the index.
`registry.geometries_from_choice` is the query the index was always for.

**The data half is CALLED: no backfill.** The pre-S1 rows stay as they are — un-keyed
artifacts of the runs that produced them. `registry.backfill_choice_digests` stays in the
tree (written, tested, idempotent, never called by `migrate`) because the same situation
recurs on any registry restored from an older copy, and a dry run reports what it would
touch without touching it: **99 of 191 rows** in `data/registry.db` today. Nothing was
modified.

---

### S2 — `compatible(a, b)`: frame-alignment feasibility · **M** · ✅ **LANDED**

**What shipped** — `compatible(a, b, *, partner, mode, donor_element)` and
`chelate_compatible(donors, vacancies, …)` in `assembly/join.py`, both returning
`Compatibility`. Tests: `tests/test_compatible.py` (18), written as a verdict table with an
assertion on the `reason` string in every negative case.

**The finding that shaped the slice: a single-point join has no strain to report.** The two
sites live in different blocks, so the rigid-body transform that brings them together is
free — for one donor onto one vacancy there is *always* a placement realising the alignment
exactly, and any number reported for it would describe where the two blocks happen to sit in
their own coordinate systems, which is not a fact about whether they can be joined. So
`compatible` returns `strain = 0.0` and says so, and what can still refuse a mono join is
roles (donor→vacancy; donor↔donor and vacancy↔vacancy each refused by name) and modes.

**The geometry appears when a join has to satisfy two constraints at once**, which is
`chelate_compatible` — and the plan's own sentence about "the convergence question
`chelate_pockets` answers within one molecule, generalised to two blocks" turns out to mean
this precisely. Two quantities, each internal to its own block and therefore invariant under
that free transform: the **bite angle** the ligand offers (from its two stored frames — D13
earning its keep, the frames are read, not re-searched) and the **vertex separation** the
metal requires (the angle between the two vacancy axes). Their mismatch is the strain.

That makes cis/trans a *feasibility* question rather than a labelling one, which is D10's
anthrarufin case seen from the assembly side. Measured, not asserted: acetate offers a
**59.7° bite**, misses a **90° cis** pair by 30.3° (feasible) and a **180° trans** pair by
120.3° (refused). Acetate is the strained end of what really chelates — acac/en/bipy sit
within 12° of a cis pair — so `MAX_BITE_MISMATCH_DEG = 40` separates two populations that are
50° apart, and a test pins the *separation* rather than the threshold, so the louder failure
is the one that fires first.

`strain` is unitless and monotone by construction — mismatch over tolerance — so
`feasible ⇔ strain ≤ 1` and the one policy number stays in one named constant.

**Two boundaries drawn on purpose:** openness is not rechecked (`open_sites` already refuses
to guess at it, and a second definition of "open" is how two definitions drift), and a
missing frame **raises** in `chelate_compatible` rather than returning infeasible — the
verdict is derived from the frames, and "unknown" reported as "impossible" is exactly what
`open_sites` refuses to do.

`tests/test_jobs.py::test_the_assembly_interfaces_exist_and_all_raise` was the stub gate and
is now `…_that_are_still_stubs_all_raise`: `compatible` left the list, `join` / `grow` /
`place_multicentre` stay on it.

---

### S3 — `join()`: the operation · **L** · ✅ **LANDED** — *the milestone's core*

**What shipped** — `join(a, b, site_a, site_b, *, mode, torsion_well, seed, with_geometry)`
returning `JoinResult`. Aligns by frame, adds the `DATIVE` bond, carries both parents' atom
maps, inherits their sites, and emits the choice vector. Tests: `tests/test_join.py` (19).

**The alignment is the join** — not a call back into the placer. Three constraints fix the
pose completely, which is D13's argument for a frame over a vector: the donor lands at
`metal + d_ml·vertex_axis`, its outward axis points back at the metal, and the roll about the
new bond is the chosen torsion **well index**, so it replays exactly. `place_mononuclear`
stays the thing that fills a whole sphere in one shot; `join` adds one ligand to a block that
already exists.

**The ligand moves, the metal stays put.** The block carrying the vacancy is the anchor, so a
growing cluster keeps one coordinate system across every addition instead of being
re-expressed after each one — which is what makes repeated joins cheap in S4.

**A call §6.3's wording does not settle: the vacancy is consumed, the donor is not.** A filled
vertex is not a vertex. A bound carboxylate oxygen, however, is still a donor site that
*exists* — and that is exactly the line D5 draws between `site_catalog` (which sites exist)
and `site_state` (what they are doing), so the donor is inherited and `refresh_state` marks it
`OCCUPIED` off the dative bond. The open count falls by two under either reading; the
difference is that the product can still answer "what bound here, and how easily did it
activate", which is what a reaction edge needs to explain itself.

**Two defects in `sites/inherit.py`, both latent until a vacancy reached it** (it predates
vacancy sites by a commit). Consumption was by atom index, and a metal's vertices all share
one atom — so the first join would have taken *all six* of an octahedral metal's vacancies and
left it reading as saturated. `merge_inherited` keyed collisions on the atom too, so it would
have refused to merge any metal carrying more than one open vertex, i.e. every metal worth
joining to. Both now key on `(atom_idx, slot)`, the pair `BuildingBlock.state_key` already
uses.

**Frames travel with the block that moved.** `inherit_sites` says frames survive untouched and
its reasoning holds for every atom the join did not move — but a join moves one whole block, so
that block's frames move with it. The transform is rigid, so every relationship a frame encodes
survives exactly and nothing is re-derived.

**`NotBuiltYet` moved to `_types`.** `descriptors.ease` imported it from `assembly.join` at
module scope, so the cycle became real the moment `join` needed `sites` and `descriptors`;
everything else was already dodging it with function-local imports. `assembly.join` re-exports
the name, so no existing import changed.

**The M5/M6 seam holds, and is tested:** a two-centre product gets its graph, both atom maps,
its inherited sites and its choice vector, then raises `NotBuiltYet` naming `place_multicentre`.
`with_geometry=False` takes the graph-level product — so a µ2-bridged bimetallic **is already
reachable from the assembly path and has a derivable L1**, with no placer in sight.

**Exit-gate evidence already in hand:** the same product built ligand-first and metal-first has
different atom orderings and **one L1 hash** (gate 1, at graph level); replaying a join
reproduces coordinates bit-identically (the foundation of gate 3).

---

### S4 — `construct.py`: deterministic build + branch tree · **L** · ✅ **LANDED**

**What shipped** — `assembly/construct.py`: `construct` (replay a recorded path),
`enumerate_constructions` → `BranchTree`, the Kind-C helpers (`geometry_branches`,
`resolve_geometry`, `metal_block`, `metal_block_branches`, `ligand_block`), and `grow` off
its stub as a thin front on the tree. Tests: `tests/test_construct.py` (25).

**Kind C got teeth, and the sharp case is CN 4.** It does not determine a polyhedron —
tetrahedral and square planar are both CN 4 and are different structures with different
chemistry, different frames and (once S5 lands) different L2 tags. `resolve_geometry(4)`
raises and names both; `geometry_branches(4)` returns both; `metal_block_branches` builds
both. CN 6 *is* determined and says so. The CN→geometry table is derived from the placer's
own `GEOMETRIES` rather than retyped, so a geometry added there cannot go missing here.
`metal_block` also derives multiplicity from the d-count and spin class rather than
defaulting to 1 — a module whose argument is "do not guess" handing back a singlet for
high-spin Fe(III) would be the loudest possible contradiction.

**Replay is exact, including through storage.** A leaf rebuilt from its recorded path
reproduces coordinates **bit-identically** and keys to the same digest — and so does one
rebuilt from the stored choice-vector JSON. That is exit gate 3, met. A step that
under-determines its join raises instead of filling the gap: a replay that substituted a
default would return *a different structure wearing the original's provenance*, which is the
worst failure available at this layer.

**What the enumerator deliberately does NOT do.** Enumerating the six vertices of an
octahedron gives six leaves that share one L1. Collapsing them here by hashing the product
would be easy and wrong — near-degenerate cis and trans *also* share an L1, and keeping them
apart is what D10 exists for. Distinguishing a symmetry duplicate from a real branch needs L2
and θ_geom clustering, i.e. S5 and S6. So the output is every branch, a cap, and a `capped`
flag; a silently truncated tree is worse than a big one.

**Measured, and it is calibration data S6 will want:**

| centre | degree | leaves | distinct L1 |
|---|---|---|---|
| CN 4 tetrahedral | 1 | 20 | 2 |
| CN 4 tetrahedral | 2 | 170 | 5 |
| CN 6 octahedral | 1 | 30 | 2 |
| CN 6 octahedral | 2 | **405** | **5** |

81 leaves per distinct L1 at CN 6 degree 2, with two partners and monodentate only. That
ratio is the size of the job S6's clustering has to do, and it is why the cap is a reported
flag rather than a quiet slice.

**Two gaps the tests found, both fixed:** a mode a donor does not offer was filtered out
*silently*, so `modes=("chelate",)` against an aqua donor produced an empty tree explaining
only that "nothing was offered" — skipped candidates are now counted into `refusals` with
their reason. And `_candidate_steps` only offered partner-donor→block-vertex, which made
polynuclear growth **unreachable through the enumerator while looking like it had merely
found nothing**; it now offers both directions, so adding a metal reaches the S3 seam and
gets M6's refusal instead of silence.

**`runner.plan` refused `degree > 1` by *calling* `grow(None, (), degree=…)` to borrow its
exception.** That became a lie the moment S4 built grow — the tripwire would have thrown an
`AttributeError` on the `None` seed instead of explaining anything. It now refuses directly,
and refuses something different: not "growth is not written" but "the run pipeline does not
drive it yet", which is S4.1.

---

### S4.1 — the enumerator behind the GUI · **M** · *deliberately deferred*

**Not part of this run, and split out on purpose so S4 stays about the algebra.** Once
`construct`/`enumerate_constructions` exist, the builder UI should drive them rather than
the single-shot placement path it drives today — and in the longer run the enumerator
becomes the builder's **default** tool, with direct placement demoted to the special case
it actually is (degree 1, one centre, every choice already made).

Scope when it is picked up, sketched only so nothing is rediscovered:

* `ui/builder.py` submits a spec; the branch tree is what a spec *means*, so the builder
  needs to show leaves before committing to a run — a preview of "this spec is 24
  structures, here is why" rather than a count that appears after planning.
* The UI's existing shape is a thin editor over `BuildSpec` (§2.5), which is the seam
  this plugs into. Do not let the enumerator become a second source of truth about a run.
* `docs/UI_BACKLOG.md` is where the interaction annoyances are already diagnosed; check it
  before designing new controls.
* The cap and the dedup rule are policy the UI must *display*, not silently apply — a
  builder that quietly drops branches is worse than one that refuses to.

---

### S5 — `l2_isomer_tag` for real · **M** · ✅ **LANDED (classifier); wiring deferred**

**What shipped** — `identity/isomers.py`: cis/trans, fac/mer and Δ/Λ from graph + geometry,
with `keys.l2_isomer_tag` delegating. `ALGO_VERSIONS["l2_isomer_tag"]` is off `"0-stub"` at
`"iso1"`. Tests: `tests/test_isomers.py` (11), and the structures under test are **built**
by S4 rather than hand-written, so "cis" means two ligands that really are 90° apart.

**No geometry, no tag — and that reframes S0(a) completely.** `put_structure` calls
`l2_isomer_tag(g)` with no coordinates (it cannot have any: the structure row is created
*before* its geometry), and nothing in `src/` passes `l2=` or calls `identity(g, geom=…)`.
So building the classifier splits **nothing**: all 39 stored structures still tag `""`, and
so does every new insert. The decision S0(a) was reserved for is not "does filling in L2
re-hash the corpus" but "when we wire the builder to pass its coordinates, what happens to
the rows that predate it" — and it can be made later, with the numbers below in hand.

**Measured against `data/registry.db`** by recomputing the tag from each structure's best
geometry:

| | count |
|---|---|
| structures, all currently `l2=""` at `algo_l2="0-stub"` | 39 |
| would gain a real tag once a builder passes geometry | **18** |
| would still tag `""` (no arrangement to report) | 21 |
| changed by this slice | **0** |

Sample of what they would become: `Mg[Cl]4[H2O]2 → O3=trans`,
`Mg[Cl]3[H2O]2[Cl⁻] → Cl2=mer,O3=trans`, `Mg[dtBK]2[Cl⁻]2[H2O]2 → Cl1=cis,O28=cis,O3=trans`.

**The consequence worth stating: L2 is not a pure function of a structure row.** A structure
ingested without coordinates legitimately carries `""`; one built from a placement carries a
tag; under `UNIQUE (l0, l1, l2)` those are two rows. Same shape as the `site_catalog` seam,
and it is why `put_structure` already takes `l2=` — the builder has coordinates at insert
time and the registry does not. That parameter was the right interface before there was
anything to put in it.

**How donors are compared.** Grouped into classes by the certificate of the ligand fragment
they belong to — so "the two ammines" is a fact about the ligands, not the element — then
each fragment *instance* is reduced to one direction and the arrangement of instances is
classified. That is what makes one rule serve MA2B2 and a bis-chelate alike. A class with
one instance is skipped: a chelate's own bite angle is not isomerism.

**A bug worth recording, because the obvious construction is wrong.** The pseudo-C3 axis of
a tris-chelate cannot be the sum of the arm directions — the six donors of a complete
octahedron sum to zero, so the arms do too, and every tris-chelate gets a zero vector. The
axis is the *normal to the plane the three arms lie in*. Δ/Λ came out as `""` until that was
fixed, which is the good failure: an unrecognisable propeller reports `AMBIGUOUS` and drops
out of the key rather than guessing a handedness.

**Δ/Λ rests on a stated convention**: IUPAC's Δ is the right-handed helix, and a right-handed
helix turns anticlockwise about its axis as it advances. The implementation is
orientation-independent (flipping the axis swaps which donor is "upper" *and* negates the
axis; the sign changes cancel) and that invariance is tested, as is mirror-antisymmetry.
Worth one check against a known crystal structure before a report leans on the absolute
label.

**Not done, deliberately:** wiring the builder to pass `l2=`. That is the step that splits
identities, it needs the S0(a) call, and it belongs next to S4.1 where the run pipeline
learns to drive the enumerator.

**Blocked on something else entirely:** the canonical test pair is Pt(NH₃)₂Cl₂ and it
**cannot be assembled** — `occlusion` reports every sp3 amine as sterically blocked, so
`open_donors()` excludes ammonia, methylamine and ethylenediamine outright. Aqua and chloride
stand in. See the note below.

---

### Found while building S5 — sp3 amine donors are unusable · *needs its own fix*

`descriptors.ease.occlusion` counts a donor's own covalently bonded neighbours as walls. A
bonded atom sits ~1.0 Å away, so for a ray leaving the donor at angle θ its perpendicular
distance is `d·sin(|θ_nb − θ|)` — below vdW(H) = 1.20 Å for separations up to ~73°. A bonded
neighbour therefore occludes a ±73° wedge *wherever it actually points*. Ammonia's hydrogens
sit **112.8° off the lone-pair axis** and still block **86%** of the 75° cone.

| donor | open | buried_vol |
|---|---|---|
| ammonia, methylamine | **0/1** | 0.86, 0.84 |
| ethylenediamine | **0/2** | 0.88, 0.86 |
| water | 1/1 | 0.41 |
| pyridine, bipy | ok | 0.56, 0.58–0.61 |
| acetate | 2/2 | 0.28–0.30 |

With `BLOCKED_OCCLUSION = 0.80` that marks every amine `BLOCKED`, so no amine can be joined
or enumerated — **ethylenediamine, the textbook chelator, has zero open donors.** Note also
that the constant's own docstring claims "a bare aqua O sits near 0.0"; it measures 0.41, so
the documented calibration no longer matches the code. Not fixed here: `buried_vol` is
stored, feeds the ease model's steric term and decides `n_open_sites`, so changing it moves
stored values and wants `ALGO_VERSIONS["site_state"]` bumped with a re-measured calibration —
its own change, not a drive-by inside S5.

---

### S6 — C2, then `l3_conformer_id` · **L** · ✅ **LANDED — C2 called on θ_geom, energy window deliberately open**

**What shipped** — `identity/conformers.py`: `l3_conformer_id` (provenance-primary),
`core_atoms` / `core_rmsd` (rigid core + coordination sphere), and `cluster`, which does the
three jobs D11 assigns geometry — collapse Kind-A duplicates, reconcile divergence, reconcile
convergence — and never names anything. Tests: `tests/test_conformers.py` (16).

**RAW constructs could not set θ_geom, and that is itself the first result.** The plan said
to put the threshold in the valley between the stochastic-duplicate peak and the real-branch
peak, and to say so if there is no valley. At RAW there is no valley because **there is no
stochastic-duplicate peak at all**:

| population (RAW) | n | min | median | max |
|---|---|---|---|---|
| Kind A — same choice vector, different embedding seed | 105 | 0.0000 | **0.0000** | **0.0000** |
| Kind B — different choice vector | 756 | 0.0000 | 0.7952 | 1.1646 |

A join is a deterministic function of its choice vector and absorbs the ligand's embedding
noise completely, so Kind-A spread is identically zero (6e-06 Å of the ligand's own MMFF
convergence). A threshold calibrated against that is calibrated against no noise. Kind B is
not a peak either — four discrete spikes, because a RAW construct only ever places ligands
on idealised polyhedron vertices.

**Relaxation re-introduces the noise, and then the valley is enormous.** 12 structures over
4 choice vectors × 3 embedding seeds, GFN2-xTB, all converged:

| population (xTB-relaxed) | n | min | median | max |
|---|---|---|---|---|
| Kind A | 12 | 0.0000 | 0.0000 | **0.0316** |
| Kind B | 54 | **0.6435** | 0.8320 | 1.1871 |

**Gap: 0.61 Å**, with nothing whatsoever in between — 9 of the 12 Kind-A pairs are at exactly
zero, and the Kind-B floor is 20× the Kind-A ceiling.

**C2, half called: θ_geom = 0.15 Å.** Placed at the geometric mean of the two bounds, which
puts the same multiplicative margin on each side — **4.7× above** the widest duplicate,
**4.3× below** the closest real branch. Geometric rather than arithmetic because these are
ratios of distances; an arithmetic midpoint would sit 20× above one population and 1.5× below
the other. Both bounds ship as named constants beside it (`CALIBRATION_KIND_A_MAX`,
`CALIBRATION_KIND_B_MIN`) and a test asserts the threshold stays between them with margin on
both sides, so the number cannot drift away from the data that set it.

**The energy window is deliberately NOT set, and the reason is a second finding.** Over the
same relaxed set the Kind-A *energy* spread reached **16.6 kcal/mol** between samples whose
cores agreed to 0.03 Å. All of that motion is outside the core — because the rigid-core rule
cuts a delocalised carboxylate C–O as if it were rotatable (ISSUES 6c), letting the whole
carboxylate swing on a charged complex. A window set from that data would bake the
core-definition defect into a stored threshold. `DEFAULT_ENERGY_WINDOW is None`, the
mechanism is built and tested (a missing energy never gates anything out — D18's rule), and
the number waits for the core fix.

---

**The calibration paid for itself by finding a defect in S3 instead.** The first run put
Kind A at a **median of 1.20 Å and a max of 1.57 Å** — wider than most of Kind B, which would
have made the whole exercise meaningless. Diagnosis: for one product pair the pyridine ring
carbons sat **2.3 Å apart** while its coordinating N moved 0.15 Å, and the free ligand's own
geometry was bit-identical between the two runs. The ligand was bound in the right place and
rotated about the M–N axis.

Cause: `_place_donor_block` computed `rotation_between(donor_axis, -vacancy_axis)` and
**discarded the frames' `ref` vectors**. Matching axes leaves the roll undetermined, and
`rotation_between` settles it with its minimal rotation — which depends on how the ligand
happened to be oriented in its own coordinate file. Re-embedding rotates a ligand rigidly, so
the same choice vector gave different rolls. That is exactly what `sites/frames.py` opens by
saying a lone outward vector cannot do, and the join was using the vector half of the frame
it was handed. Now it maps frame onto frame through an orthonormal triad, with the torsion
well applied to the vacancy's own `ref`. Kind-A spread: **1.56 Å → 0.0000 Å**.

The original procedure, for the record:

1. Build the M5 fixture set (S3/S4 produce it as a by-product).
2. Plot the pairwise RMSD distribution over **rigid core + coordination sphere only** — not
   whole-molecule RMSD, which §4.2 says is a bad sole trigger.
3. Put θ_geom in the valley between the stochastic-duplicate peak (Kind A) and the
   real-branch peak (Kind B). If there is no valley, say so in the write-up — that is
   information about the fixture set, and picking a number anyway is how a threshold becomes
   folklore.
4. Set the energy window against the same set. Note it gates *clustering*, so it wants a
   backend; M7's `reaction_balanced_energy` is already there.

Then `l3_conformer_id`: **provenance-primary** (the label is the choice vector),
**geometry-verifier** (clustering only collapses Kind-A duplicates and reconciles
divergence/convergence). Choice-vector dedup runs **before** geometric clustering — that
ordering is the declared escape hatch for the conformer-explosion risk, along with the
per-`(L1,L2)` cap and the rule that `TORSION_FREE` sites never branch.

**Deliverable beyond code** — the distribution plot and the chosen numbers, written up with a
D-number. A gate resolved only in someone's head is how the two documents drift apart.

---

### S7 — registry wiring + the three headline validations · **M** · ✅ **LANDED**

**What shipped** — `assembly/persist.py` (`store_block`, `store_construction`), and
`tests/test_m5_exit_gates.py` (8) proving the gates *after a round trip through SQLite*
rather than in memory. Nothing in it writes to the database: every insert goes through
`registry.api`, and this module is the orchestration that knows the order and supplies the
three things assembly knows that the registry cannot derive.

| # | Gate | Result |
|---|---|---|
| 1 | **One node, two routes** | ✅ two build orders → **one** `structures` row, **two** `reactions` edges, `n_incoming_routes == 2`, two geometries under one identity |
| 2 | **cis/trans survive** | ✅ one L1, two rows, different L2 tags; clustering keeps them apart at the calibrated θ_geom with identical energies and a wide-open window |
| 3 | **Replay** | ✅ build → store → read the row back → rebuild from its `choice_vector_json` → coordinates match the stored `.xyz` to 1e-6, same digest |
| — | Ambiguous spec refuses | ✅ `resolve_geometry(4)` raises and names both polyhedra |

**Three things `persist` supplies that the registry cannot derive.** The **L2 tag**, because
`put_structure` has no coordinates at insert time (the structure row precedes its geometry) —
which is why it has always taken `l2=`. The **sites, by inheritance**, because re-perceiving
an assembled complex returns no donors at all (ISSUES 4 / #16) — the runner's `_record_sites`
perceives because it builds from a molecule; this path must not. And the **provenance**: the
blocks that went in, the atom map, the choice-vector digest, the depth.

**A gap `incoming_routes` had.** It selected `id, kind, intermediate, depth, note,
created_at` — so two routes to one node came back looking identical apart from their id and
note, and "the same product reached two ways" was a claim a reader had to take on trust. It
now projects `choice_vector_digest` and `atom_map_json` too, which is what makes the edges
distinguishable and is the half of D2 that lives off the node.

**A correction to this plan's own gate-2 wording.** The first version of the test asserted
cis and trans survive θ_geom = 99 Å. They do not, and should not — two choice vectors whose
geometries genuinely coincide are exactly the convergence case D11 asks clustering to
reconcile, and 99 Å does not mean "near-degenerate", it means "every structure is one
structure". D10's claim is about the **energy** gap: the discriminator is relevance, not ΔE.
So the test now gives the pair identical energies and a wide-open window — if ΔE were doing
the work that would merge them — and they stay apart because their cores are **1.9 Å**
apart, 13× θ_geom.

**One substitution, stated rather than hidden.** The plan names anthrarufin–Cu for gate 2;
anthrarufin is bidentate and a chelate cannot be joined yet (ISSUES 3 / #15). The gate is met
with Pt(OH₂)₂Cl₂ — monodentate, buildable, a genuine cis/trans pair — and the test says so in
its docstring, with a note to re-run it against anthrarufin when #15 lands.

---

## 5. Risks, with the escape hatch named in advance

| Risk | Early signal | Hatch |
|---|---|---|
| **Identity retrofit pressure** — "just add a flag to L1 so the join works" | a proposed change to what enters the certificate | the M2 discrimination table is the contract. Changing it costs a version bump and a corpus re-hash, and that cost is the point. |
| **L3 explosion** | thousands of near-identical rows per `(L1,L2)` | choice-vector dedup before geometric clustering; cap per `(L1,L2)`; `TORSION_FREE` never branches |
| **C2 has no valley** | the RMSD histogram is unimodal | report it and set θ_geom from the open-site-flag trigger instead — D11 makes the flag flip the *primary* trigger and RMSD the secondary one, so this is a degradation, not a blocker |
| **S3 drifts into M6** | joins start needing inter-centre distances to decide feasibility | the seam in S3 is a hard stop: graph + identity + provenance complete, geometry raises. Two-centre geometry is not M5 scope under any schedule pressure. |

---

## 6. Bookkeeping when M5 lands

- `ALGO_VERSIONS`: `l2_isomer_tag` and `l3_conformer_id` off `"0-stub"`; ~~a new
  `choice_vector` entry~~ ✅ `"cv1"`, added in S1.
- Design doc: D-numbers for **C2**, for the S0(a) L2 backfill-vs-version call, and for
  `MAX_BITE_MISMATCH_DEG` (S2's one policy constant). Changelog line for each. S0(b) is
  already decided — no digest backfill — and needs only the changelog line.
- `CODE_ARCHITECTURE.md`: flip the four 🔴 rows (`l2_isomer_tag`, `l3_conformer_id`,
  `assembly/join.py`, `choice.py`/`construct.py`); close the "L2 is `''` everywhere" seam in §6.
- Delete this file.

## 7. Order, at a glance

```
S0 decisions ─► S1 choice.py ─► S2 compatible ─► S3 join ─► S4 construct ─► S5 L2
                                                    │                         │
                                                    └──────► S6 C2 + L3 ◄─────┘
                                                                 │
                                                                 ▼
                                                          S7 exit gates
```

S1 and S2 are independently testable and can be done in either order if S2 goes first without
emitting a vector. Everything after S3 depends on the atom map being right, so S3's tests are
the ones worth over-writing.
