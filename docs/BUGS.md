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
| [B2](#b2) | honesty | `identity/keys.py` | `l2_isomer_tag` is `''` for all 39 structures, so cis and trans are one row. **Resolution decided (D19)** |
| [B3](#b3) | honesty | registry data | 16 structures have no `site_state`; `n_open_sites` is NULL, not 0. **Self-healing (D19)** |
| [B4](#b4) | cosmetic | `ui/static/runs.html` | "What was attempted" renders `molecule undefined · undefined×0-dentate` |
| [B5](#b5) | undecided | `ui/static/builder.html` | The builder posts `spec_version: 1` and no test covers the migration it relies on |
| [B6](#b6) | undecided | `ui/static/index.html` | The registry page never refreshes, so its counts go stale silently |
| [B7](#b7) | blocked | `scripts/ingest.py` | The legacy `.xyz` corpus cannot be ingested: no charge, no multiplicity |
| [B8](#b8) | correctness | `descriptors/ease.py` | Every sp3 amine donor reads as sterically blocked, so no amine can be joined |
| [B9](#b9) | correctness | `identity/conformers.py` | The rigid-core rule cuts a delocalised carboxylate C–O as if it were rotatable |
| [B10](#b10) | undecided | `identity/isomers.py` | Δ/Λ is self-consistent but its absolute assignment is unverified |
| [B11](#b11) | cosmetic | `sites/model.py` | `vacancy_sites` normalises by hand, differently from `_linalg.unit` |

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

**Still open until M5 lands**, because the misreport is live until then: the registry presently
claims cis and trans are one structure. Accepted cost of D19: the same species built before and
after M5 can occupy two rows.

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

