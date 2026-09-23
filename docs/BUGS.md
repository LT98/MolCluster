# Bug tracker — open only

Everything here is **unresolved**. When something is fixed and has a test, move the entry to
[`archive/BUGS_resolved.md`](archive/BUGS_resolved.md) — do not leave it here struck through.
The point of this file is that its length is the size of the problem.

Scope: things that behave wrongly or report wrongly. *Unbuilt* functionality is not a bug —
that is a milestone, and it lives in [`PLAN_implementation.md`](PLAN_implementation.md).
A stub that raises loudly is working as designed (ground rule 8).

**Severity** — `correctness` a stored value or a claim is wrong · `honesty` the code is right
but what it reports is misleading · `cosmetic` it looks wrong and misleads nobody ·
`undecided` known behaviour that nobody has explicitly chosen.

| # | Severity | Where | One line |
|---|---|---|---|
| [B2](#b2) | honesty | `runner.py` | Every stored row's `l2_isomer_tag` is `''`, so cis and trans are one row. The classifier exists; the run pipeline declines to use it (**D19**, **D22**) |
| [B3](#b3) | honesty | registry data | 16 structures have no `site_state`; `n_open_sites` is NULL, not 0. **Self-healing (D19)** |
| [B4](#b4) | cosmetic | `ui/static/runs.html` | "What was attempted" renders `molecule undefined · undefined×0-dentate` |
| [B5](#b5) | undecided | `ui/static/builder.html` | The builder posts `spec_version: 1` and no test covers the migration it relies on |
| [B6](#b6) | undecided | `ui/static/index.html` | The registry page never refreshes, so its counts go stale silently |
| [B7](#b7) | blocked | `scripts/ingest.py` | The legacy `.xyz` corpus cannot be ingested: no charge, no multiplicity |
| [B8](#b8) | correctness | `descriptors/ease.py` | Every sp3 amine donor reads as sterically blocked, so no amine can be joined |
| [B9](#b9) | correctness | `identity/conformers.py` | The rigid-core rule cuts a delocalised carboxylate C–O as if it were rotatable |
| [B10](#b10) | undecided | `identity/isomers.py` | Δ/Λ is self-consistent but its absolute assignment is unverified |
| [B11](#b11) | cosmetic | `sites/model.py` | `vacancy_sites` normalises by hand, differently from `_linalg.unit` |
| [B14](#b14) | correctness | `sites/perception.py` | Pyrazolate's two equivalent N type differently, and one of them cannot bridge |
| [B15](#b15) | correctness | `sites/perception.py` | A bare oxide `[O-2]` types as `hydroxide_O`, so an oxo reports a hydroxide's pKa |
| [B16](#b16) | correctness | `scripts/run_spec.py` | `--store` is ignored: spawned workers re-resolve `data_root()`, so blobs land in the default store |
| [B17](#b17) | honesty | `ui/static/index.html` | Changing the geometry does not re-render the detail panel, so `method` and `converged` go stale |
| [B19](#b19) | honesty | `ui/static/graph.html` | One structure drawn at two heights by two routes, with the reason — a different basis — only in the legend |
| [B20](#b20) | correctness | `runner._build_sphere` | A charged co-ligand's charge is left out of the complex's net charge: Ni(II) + 3 `[Cl-]` is stored q+2, not q−1 |
| [B21](#b21) | correctness | `geometry/placer.py` | CN 6 with a chelate and ≥2 reserved-empty vertices is refused (`best bite angle 180`), though cis-M(L–L)2(□)2 exists |

---

## B2

**`l2_isomer_tag` is `''` everywhere, so cis and trans collapse into one row.** `honesty` ·
`identity/keys.py`

Measured 2026-09-14: 39 of 39 structures carry an empty L2. The stub is by design (M2 fixed
the signature, M5 fills the body), but the *consequence* is a live misreport — the registry
currently claims two isomers are one structure.

**Called: version bump, not backfill — D19.** `ALGO_VERSIONS["l2_isomer_tag"]` moves off
`0-stub` when M5 fills the body; rows written under the stub keep `l2_isomer_tag = ''` and are
never re-derived. Backfilling was refused because it would rewrite stored identities and every
`reactions` edge pointing at them — the one thing here that is not regenerable.

`structures.algo_l2` already records the generation per row (all 39 currently read `0-stub`),
so a `''` is readable as *"this predates L2"* rather than *"this has no isomerism"*.

**M5 landed and this is still open, for a reason worth stating precisely — the classifier is
built and the pipeline deliberately does not use it.** `identity/isomers.py` returns real tags,
`ALGO_VERSIONS["l2_isomer_tag"]` is `iso1`, and `assembly.persist.store_block` derives a tag
whenever it has coordinates. But `runner` passes `l2=""` on purpose (**D22**): the run pipeline
does not tag isomers anywhere else, so a route that tagged its own product would file it apart
from the node every other route reached — a worse failure than the one this entry describes,
because it breaks D2 rather than merely blurring it.

So what closes B2 is **wiring both paths together**, not filling in a stub. Accepted cost of
D19 when that happens: the same species built before and after can occupy two rows,
distinguishable by `structures.algo_l2`.

---

## B3

**16 structures have no `site_state`.** `honesty` · registry data

Measured 2026-09-14: 16 of 39 structures were built before M4's second half, so they have a
catalog and no state, and `n_open_sites` is NULL for them.

NULL is the correct value — absent is not zero (ground rule 9 / invariant 9) — and the bug is
that nothing stops a consumer from treating it as zero and reporting a fully-occupied structure.

**Now self-healing, and the cause turned out to be shared.** These 16 are the same 16 whose
`site_catalog` was *empty*, because perception dropped every donor the moment it coordinated
(archived). Under D19 a `perception/1` catalog is rewritten the next time anything touches its
structure, and `_record_sites` writes fresh state straight after — so each one is repaired as
it is next built or relaxed.

**What is still open** is the read side: nothing yet makes a NULL `n_open_sites` loud to a
consumer that treats it as 0. Until every row has been touched, the registry is a mix.

*(An earlier revision of this note said 569 of 594. That was a different registry; the numbers
above are a fresh query against `data/registry.db`.)*

---

## B4

**The run inspector's most useful column renders `undefined`.** `cosmetic` ·
`ui/static/runs.html:433`

For `place` tasks the "what was attempted" column reads
`molecule undefined · undefined×0-dentate`.

**Cause, confirmed:** the payload shape changed and the page did not follow it. A `place`
payload used to name one molecule and a count; it now carries a `components` list
(`runner._components` at [runner.py:246](../src/mofsbu/runner.py:246) still accepts both
shapes). The JS reads the flat `p.molecule` / `p.donors` / `p.n_ligands` keys, which no new
payload has.

Fix is to read `components` with the same both-shapes fallback the runner already implements,
rather than to teach the planner to write the old keys back.

---

## B5

**The builder posts `spec_version: 1` and relies on the migration chain.** `undecided` ·
`ui/static/builder.html:521`

It works, and it is arguably the right design — the page never has to know the current
version. But it is load-bearing behaviour that no test covers, and `SPEC_VERSION` is at 5.
A migration that silently stopped handling the v1→v2 hop would surface as wrong builds, not
as an error.

Either add the test that a v1 spec from the builder migrates to the current version with the
fields the page intended, or make the page post the current version. The first is better;
either is better than the present state.

**v8 made the hazard concrete.** The v7→v8 migration sets `co_ligand_counts: fill` on any
spec that omits it, so a page that forgot to send the field would silently plan `fill` while
showing `range`. The page sends it explicitly, and
`test_the_co_ligand_count_is_offered_and_an_old_spec_reads_as_fill` covers that one field
through a v1-labelled post; the other fields are still uncovered.

---

## B6

**The registry page never refreshes.** `undecided` · `ui/static/index.html`

`index.html` has no polling, so a count shown there can be stale after a run finishes in
another tab. Arguably correct — a registry view that moves under you while you are reading it
is worse — but nobody decided it, which is why it is here rather than in the design doc.

Note that `/runs` deliberately went the other way (rev 24, item 2 & 3 in the archive): it
polls, cheaply, and backs off. Whatever is decided here should be decided against that.

---

## B7

**The legacy `.xyz` corpus cannot be ingested.** `blocked, needs a call` ·
`scripts/ingest.py`

M3's real-data test is stuck on something structural rather than incidental: **an `.xyz` file
carries no charge and no multiplicity**, and `l0_composition` refuses to guess either (ground
rule 5). Bond perception from coordinates is straightforward; charge and spin are not
recoverable from geometry.

Options, in increasing order of effort:

* **(a)** a manifest CSV mapping file → (charge, multiplicity), part auto-filled from the
  legacy filenames that encode it (`_q-4`, `Zn2+`) and part filled by hand;
* **(b)** infer charge from perceived ligand protonation plus metal oxidation state — a guess
  wearing a rule;
* **(c)** skip the corpus and let M6 regenerate it.

**(a) is the honest one** and it keeps the filename as a migration *input* that is then
discarded, so meaning still does not live in filenames going forward.

Until this is called, M3's real-data test is the fixture round-trip instead, and the
duplicate-collapse count that ingesting the corpus would have produced — itself a result
worth having — does not exist.

---

## B8

**Every sp3 amine donor reads as sterically blocked, so no amine can be joined.**
`correctness` · `descriptors/ease.py`

`occlusion` counts a donor's own covalently bonded neighbours as walls, so `refresh_state`
marks amines `BLOCKED` and `BuildingBlock.open_donors()` excludes them entirely.

Measured 2026-09-15:

| donor | open donors | `buried_vol` |
|---|---|---|
| ammonia | **0 / 1** | 0.86 |
| methylamine | **0 / 1** | 0.84 |
| ethylenediamine | **0 / 2** | 0.88, 0.86 |
| water | 1 / 1 | 0.41 |
| pyridine | 1 / 1 | 0.56 |
| 2,2'-bipyridine | 2 / 2 | 0.58, 0.61 |
| acetate | 2 / 2 | 0.28, 0.30 |

**Structural, not a threshold to retune.** A ray leaving the donor at angle θ passes within
`d·sin(|θ_neighbour − θ|)` of a neighbour at distance `d`. With an N–H bond of 1.02 Å and
vdW(H) = 1.20 Å that is under the radius for angular separations up to **~73°**, so a bonded
neighbour occludes a ±73° wedge *regardless of where it points*. Ammonia's hydrogens sit
**112.8° off the lone-pair axis** and still block **86%** of the 75° cone. For a terminal
donor the number largely measures "does this donor have neighbours".

**Costs now:** ethylenediamine, the textbook chelator, has zero open donors; Pt(NH₃)₂Cl₂
cannot be assembled, so `tests/test_isomers.py` and `tests/test_m5_exit_gates.py` stand in
aqua and chloride.

**`BLOCKED_OCCLUSION`'s calibration note is part of the same bug.** It documents "a bare aqua
O sits near 0.0"; aqua measures **0.41**. The carboxylate figure (0.2–0.4) still holds. A
threshold whose stated calibration is wrong cannot be safely moved, so re-measure the whole
panel and rewrite the note from the new numbers as part of the fix.

**Fixing it moves stored values:** `buried_vol` lives in `site_state`, feeds the ease model's
steric term and decides `n_open_sites`. Bump `ALGO_VERSIONS["site_state"]` (now `"1"`) and
possibly `"ease_model"`. Under D19 the existing rows keep their own version's answer.

**Likely fix:** exclude atoms covalently bonded to the donor from the ray test — they are what
*defines* the axis and cannot be in the way of it — or damp their radius by distance.

Tracked as [#13](https://github.com/LT98/MolCluster/issues/13) and
[#14](https://github.com/LT98/MolCluster/issues/14).

---

## B9

**The rigid-core rule cuts a delocalised carboxylate C–O as if it were rotatable.**
`correctness` · `identity/conformers.py`

`rotatable()` reuses `sites.model`'s rule — single, acyclic, non-aromatic — so that "rigid"
means one thing across the codebase. That rule is order-blind by design (D15 keeps bond order
out of identity), so it cannot tell a carboxylate's delocalised C–O from an ether's, and cuts
both.

So the rigid core of an acetate complex stops at the coordinating oxygen and excludes the
carboxylate carbon, and an acetate torsion registers as **exactly zero** core RMSD. That
agrees with design §4.2's worked examples — a non-coordinating carboxyl torsion is explicitly
a *trivial* conformer — but it is right by accident, and the same cut would hide real motion
in an amide or an ester.

**What it currently blocks:** the L3 energy window. Over 12 xTB-relaxed structures the Kind-A
*energy* spread reached **16.6 kcal/mol** between samples whose cores agreed to 0.03 Å, all of
it motion outside the core. A window calibrated from those numbers would bake this defect into
a stored threshold, so `DEFAULT_ENERGY_WINDOW is None` until this is fixed. θ_geom beside it
**is** calibrated (0.15 Å) and unaffected.

**Error direction is safe:** a smaller core makes the comparison more permissive, never less —
it can fail to split two conformers and can never merge two that differ elsewhere.

Tracked as [#20](https://github.com/LT98/MolCluster/issues/20) and
[#19](https://github.com/LT98/MolCluster/issues/19).

---

## B10

**Δ/Λ is self-consistent but its absolute assignment is unverified.** `undecided` ·
`identity/isomers.py`

`_chirality` assigns Δ to the right-handed propeller, from IUPAC's definition of Δ as the
right-handed helix plus the geometric fact that a right-handed helix turns anticlockwise about
its axis as it advances.

**Tested:** mirror images get opposite labels; rotations do not move the label; the answer does
not depend on which end of the C3 axis you look from (flipping it swaps which donor counts as
"upper" *and* negates the axis — the sign changes cancel).

**Not tested:** that the label on a real Δ complex is `Delta` and not `Lambda`. A systematic
inversion would be invisible to every test here and visible in every report. One check against
a known crystal structure settles it. Until then the discrimination is trustworthy and the
absolute assignment is provisional.

Tracked as [#17](https://github.com/LT98/MolCluster/issues/17).

---

## B11

**`vacancy_sites` normalises by hand, differently from `_linalg.unit`.** `cosmetic` ·
`sites/model.py`

It divides by `max(norm, 1e-9)`; `geometry._linalg.unit` falls back to the z axis. Identical
for any real input, different for a zero-length direction vector — one returns a near-infinite
vector, the other a valid axis.

Left out of the `_linalg` consolidation deliberately: swapping it is a behaviour change in the
degenerate case, not a de-duplication. Worth deciding which answer is wanted rather than
leaving two.

Tracked as [#22](https://github.com/LT98/MolCluster/issues/22).

---

## B14

**Pyrazolate's two equivalent nitrogens type differently, and one of them cannot bridge.**
`correctness` · `sites/perception.py`

Measured 2026-09-16:

| molecule | donors |
|---|---|
| pyrazole `c1cc[nH]n1` | `azolate_N`, `pyridyl_N` |
| pyrazolate `c1cc[n-]n1` | **`amide_N`**, **`pyridyl_N`** |

The pyrazolate anion's two nitrogens are equivalent by resonance — the charge is delocalised
around the ring — so typing them as two different donor types is the same defect `cff47a8`
fixed for oxo-acids under the heading *perception is resonance-invariant, an oxo-acid is one
donor group*. Azolates were not covered by that change.

**What it costs:** `pyridyl_N` declares `mono` only, so a pyrazolate N,N bridge is refused on
one end whatever the other end offers. Pyrazolate-bridged dimers are the standard
non-carboxylate test of a bridging model, so M6 has no way to check that mechanism A generalises
beyond carboxylates until this is fixed.

---

## B15

**A bare oxide types as a hydroxide, so an oxo reports a hydroxide's pKa.** `correctness` ·
`sites/perception.py`

Found while fixing [B13](archive/BUGS_resolved.md). `_classify_anionic` returns
`"hydroxide_O"` for any anionic oxygen with no heavy neighbour, and that branch does not look
at the hydrogen count:

| molecule | perceived as | should be |
|---|---|---|
| `[OH-]` | `hydroxide_O` | `hydroxide_O` ✅ |
| `[O-2]` | `hydroxide_O` | an oxo — its own type |

The consequence is not cosmetic. `descriptors.tables.donor("hydroxide_O")` carries pKa 15.7,
which is the number for water losing its second proton; an oxo has already lost both and does
not deprotonate at all. So `activation_ease` scores a µ-oxo as though it had an activation step
it cannot have, and D18's rule that "a donor with no pKa never self-flags" is bypassed by
giving it someone else's.

**The fix is a new donor type, which makes it a reference-data change rather than a one-line
one.** `oxo_O` needs a row in `data/reference/donor_descriptors.tsv` with its own pKa (absent,
not large), HSAB class and denticity, plus an entry in `sites.frames._BINDING_MODES` offering
`mu2`/`mu3`/`mu4`. `donor()` raises on an unknown key by design, so adding the type without the
row fails loudly rather than silently — which is the right order to do it in.

**Not blocking M6/S2.** µ3 and µ4 oxo bridges are placed as **centres**
(`geometry.placer.bridging_metal_positions`), a path that never asks perception for a donor
type. It bites whenever an oxo-centred node is stored and its sites are scored.

---

## B16

**`run_spec.py --store` is ignored, and the run still succeeds.** `correctness` ·
`scripts/run_spec.py`

A run launched with `--db data/proof.db --store data/proof_store` wrote its database where it
was told and its blobs somewhere else: `data/store`, the default. `data/proof_store` was
created and left empty.

The parent honours the flag; the spawned workers do not. `execute_run` starts workers with the
**spawn** context, so a worker begins from nothing and re-reads `config.data_root()` — which
knows `MOFSBU_DATA` and the repo root, and nothing about a `--store` argument parsed in another
process. The database path survives because it is passed explicitly; the store path is not.

Nothing detects it. The blobs are content-addressed, so they are written and read back under
the same digest within a run, and the mismatch only shows up later when someone moves the
database and finds its geometries are not beside it.

**Fix direction:** either pass the store path into each worker the way the database path already
is, or refuse a `--store` that differs from `store_root()` rather than accepting a flag that
does not take effect. The second is smaller and honest; the first is what the flag promises.

---

## B17

**Changing the geometry leaves the detail panel showing the previous one's method.**
`honesty` · `ui/static/index.html`

The geometry `<select>`'s `onchange` calls `showGeometry()`, which repaints the 3D viewer and
sets `state.geomId`. It does not call `renderDetail()`. But `renderDetail` reads the selected
geometry once, at render time, and builds the `method` and `converged` rows from it.

So selecting a second geometry moves the structure on screen while the rows beside it go on
describing the first — most visibly when the two differ in exactly the way that matters, a
raw construct and its ML relaxation.

Found while adding per-route energies to the provenance panel; those are scoped to the
equation rather than to the selected geometry, so they do not inherit it. Any future
geometry-scoped row would.

**Fix direction:** have `onchange` re-render the panel rather than only the viewer, or move the
geometry-scoped rows into `showGeometry` so there is one writer for them.

---

## B21

**An octahedral centre carrying chelates and two or more empty vertices is refused.**
`correctness` · `geometry/placer.py` (`place_mononuclear` with `reserve=`)

Zn(II) + catechol (5-ring, dianionic chelate), CN 6 octahedral, construct: `place` refuses
Zn(cat)2 with two empty vertices and Zn(cat) with four, both `placer_refused` — *"no vertex
set on this geometry can host a 2-dentate ligand: best bite angle 180 deg, need 55-115"*.
Reproduced on main with `co_ligand: null, allow_unsaturated: true`, so it predates D26. The
geometry exists: reserve one cis pair and the four vertices left still hold two cis pairs.
The suspicion is the order — `cis_vertices` reserves the lowest-index cis pair and
`_assign_targets` then fills greedily — but that is unconfirmed.

Why it matters now: under `co_ligand_counts: range` (D26) these rungs are asked for by
default, so the refusal shows up as rejected `place` tasks and, above them, `grow` steps
rejected with `pathway_parent_missing`. The ladder is incomplete there, loudly.
`test_the_co_ligand_saturated_series_connects` pins `fill` for that reason.

---

## B20

**A charged co-ligand's charge is not counted in the complex's net charge.** `correctness` ·
`runner._build_sphere`

`charge = metal.oxidation_state + ligand_charge`, where `ligand_charge` sums the spec's
molecule components only; the `n_co` co-ligand copies contribute nothing. For water that is
invisible. For `[Cl-]` it is wrong: measured on a construct run (Ni(II) hs, catechol, CN 4,
`co_ligand: "[Cl-]"`), `place` stored Ni[Cl]3 and Ni[Cl]3[catechol] as **q+2**; the true net
charges are −1 and −1. The co-ligand's own graph (`[Cl] q-1`) is right.

It surfaced through D26: a co-ligand `grow` joins the free `[Cl-]` block onto the rung below,
and `join` sums block charges, so the step's product came out q+1 — a different identity from
the q+2 node `place` built — and the ladder forked (12 structures where 8 were expected). The
planner now **does not queue co-ligand steps for a charged co-ligand** and records a diagnostic
naming this entry; `_execute_grow` refuses a hand-built one with `co_ligand_charged`.

**Fix direction** — add `formal_charge(co_ligand) × n_co` to the charge in `_build_sphere`.
That changes the L0 identity of every stored complex built with a charged co-ligand, so it is a
registry decision (a versioned re-derivation, D19) and not a one-line patch; once it lands, the
planner guard and the executor refusal come out together.

---

## B19

**Two routes draw the same structure at two different heights, and the chart does not say
why.** `honesty` · `ui/static/graph.html`

Reported as a suspected bad reference on a deprotonation step. It is not one — the numbers are
right and self-consistent — but the chart presents them in a way that makes a reader conclude
otherwise, which is the defect.

Measured on `data/mvp_ni_thq_cl.db`, walking back from structure 17:

| route | walk | y(#81) | basis |
|---|---|---|---|
| A | `17 ← 33 ← 74 ← 81` | **+9.92468** | `83x1` |
| B | `17 ← 38 ← 77 ← 72 ← 81` | **+2.56853** | `83x2` |

The gap is **7.35615 eV**. The free-ligand deprotonation `tHQ + H2O → tHQ⁻ + H3O⁺` (reaction
633) is **7.35616 eV**. They agree to four decimal places, and that is the whole explanation:
route A sheds one proton on the way to the target and route B sheds two, so the two walks
arrive at structure 81 carrying different spectators and measured against different references.
Both are balanced. Neither number is wrong. They are energies of different systems, and
`pathways.route` already says so — `basis` is `83x1` against `83x2`, where 83 is H3O⁺.

What fails is where that is said. `basis` is reported per node by the API and rendered **only
in the route legend**, while the thing a reader actually looks at is two bars at two heights
under one name. The page has the information that would stop the misreading and does not put it
where the misreading happens.

This is the same class as the free-vs-bound dE caveat: the number travels with its
qualification everywhere, or it gets read as something it is not.

**Fix direction** — three parts, in order of how much they buy:

1. **Mark the node, not just the route.** When one `structure_id` appears at more than one
   height in a drawn chart, draw both with a shared marker and say on each what its basis is.
   The data is already there: `nodes[i].basis` per route.
2. **Name the difference in the units a chemist reads.** `83x1` vs `83x2` is a structure id and
   a count. It should render as `−1 H⁺` vs `−2 H⁺`, and the gap between two such nodes should
   be offered as *"these differ by one proton transfer"* rather than left as a subtraction the
   reader performs and then distrusts.
3. **Decide what the chart does about it.** Two options, and this is a real call rather than an
   oversight: either keep drawing both heights and label them, or offer a *common-basis* view
   that re-references every route to the most-shed basis, so the curves become directly
   comparable at the cost of no longer being the raw stored subtraction. The second is more
   useful and more dangerous; it must never be the default and must be badged when on.

**Not** in scope of the fix: the stored energies, the deprotonation edges, or
`pathways.route`'s arithmetic. All three were checked against this case and are correct.

**What consumed-by hops change here (still open).** The graph now walks edges forward as well as
back (`c<id>` legs, docs/gui/05_graph.md), and that moves this entry in two ways:

- **The basis can now differ by more than protons.** A route that turns at a shared parent
  sheds whatever the backward legs released — Cl⁻ in the Ni/tHQ/Cl exchange, basis `7x2` — so
  fix direction 2 cannot assume a basis token is `n H3O⁺`; it has to render any species.
  (The spectator tally is now netted, so a piece taken up and given back no longer appears on
  both sides; that makes a basis shorter, not different in kind.)
- **A proton transfer can now be walked as well as subtracted.** A deprotonation edge is offered
  as a consumed-by hop from its protonated parent, so a route can take up or give back a
  proton explicitly and the page prices it with the route's net equation. That makes the kind
  of gap in the table above something a reader can compose and check by hand; it is not fix
  direction 3's common-basis view, and it changes nothing about the defect itself — the chart
  still draws #81 at two heights with the reason only in the legend. Not re-measured on the
  table's routes.
