# Known issues

Defects and live risks, with what was measured rather than what was suspected. Scheduled
work lives in `PLAN_implementation.md` and `WORKPLAN_M5.md`; a milestone that has not been
built yet is not an issue and does not belong here.

Each entry says what is wrong, how it was measured, what it costs today, and what fixing it
would move — that last part matters because most of these touch a **pinned recipe version**,
and changing a stored quantity is never a drive-by.

Ordered by what is blocking work now.

---

## 1. Every sp3 amine donor reads as sterically blocked · **blocking**

`descriptors.ease.occlusion` counts a donor's own covalently bonded neighbours as walls, so
`sites.state.refresh_state` marks amines `BLOCKED` and `BuildingBlock.open_donors()` excludes
them. **No amine ligand can be joined or enumerated.**

Measured:

| donor | open donors | `buried_vol` |
|---|---|---|
| ammonia | **0 / 1** | 0.86 |
| methylamine | **0 / 1** | 0.84 |
| ethylenediamine | **0 / 2** | 0.88, 0.86 |
| water | 1 / 1 | 0.41 |
| pyridine | 1 / 1 | 0.56 |
| 2,2'-bipyridine | 2 / 2 | 0.58, 0.61 |
| acetate | 2 / 2 | 0.28, 0.30 |
| chloride, acetonitrile | ok | 0.00 |

The cause is structural, not a threshold to retune. A ray leaving the donor at angle θ passes
within `d·sin(|θ_neighbour − θ|)` of a neighbour at distance `d`. With an N–H bond of 1.02 Å
and vdW(H) = 1.20 Å that is under the radius for angular separations up to **~73°**, so a
bonded neighbour occludes a ±73° wedge *regardless of where it actually points*. Ammonia's
hydrogens sit **112.8° off the lone-pair axis** — nowhere near the approach path — and still
block **86%** of the 75° cone. For a terminal donor the number largely measures "does this
donor have neighbours", not "can a metal get in".

**Costs today:** ethylenediamine, the textbook chelator, has zero open donors; the canonical
cis/trans fixture Pt(NH₃)₂Cl₂ cannot be assembled at all (`tests/test_isomers.py` uses aqua
and chloride instead).

**Fixing it moves:** `buried_vol` is stored in `site_state`, feeds the ease model's steric
component and decides `n_open_sites`. Needs `ALGO_VERSIONS["site_state"]` bumped (currently
`"1"`), possibly `"ease_model"` (`"floor1"`), and a re-measured calibration — see issue 2.
Likely shape of the fix: exclude atoms covalently bonded to the donor from the ray test (they
are what *defines* the axis and cannot be in the way of it), or damp their radius by distance.

---

## 2. `BLOCKED_OCCLUSION`'s calibration note no longer matches the code

`sites/state.py` documents the threshold as: "a bare aqua O sits near 0.0, a carboxylate O in
an open carboxylate near 0.2-0.4, and a donor pointing into its own ring system above 0.8".

Measured today: **aqua O is 0.41**, carboxylate is 0.28–0.30, and *ammonia* — which is not
pointing into anything — is 0.86. The carboxylate figure still holds; the aqua figure is off
by the width of the band it was meant to anchor.

So the numbers the constant was chosen against are not the numbers the code produces. Whether
that is drift from a later frame change or a note written against `heavy_only` frames is not
established. Re-measure the whole donor panel and rewrite the note from it as part of issue 1;
a threshold whose stated calibration is wrong is a threshold nobody can safely move.

---

## 3. A chelate cannot be joined — only judged

`assembly.chelate_compatible` answers whether a donor pair can span two vertices, and
`assembly.join` places **one** donor. There is no two-point join, so a bidentate ligand cannot
actually be attached by the assembly path.

The failure mode is at least honest: `join(..., mode="chelate")` is refused by `compatible`
before any coordinates are produced, because a vacancy offers `mono` only.

```
vacancy does not bind 'chelate'; it offers mono
```

**Costs today:** `tests/test_isomers.py` builds its tris-chelate by hand rather than by
joining, because there is no other way to get one. Every enumeration is monodentate.

**Not the same as M6.** M6 is inter-*centre* geometry; this is two points on one centre, which
`chelate_compatible` already has the arithmetic for. It is a missing S3 case, not a scheduled
milestone.

---

## 4. Re-perceiving an assembled complex loses every donor

`sites.perception` treats a metal as an ordinary heavy neighbour, so a donor that is already
coordinated is not perceived as a donor.

Measured on Zn(H₂O)₂ built by the placer:

```
free water           -> [(0, 'aqua_O')]
the complex          -> []          # not "occupied": absent
```

Every donor in the complex vanishes, not merely the bound ones. `SiteStatus.OCCUPIED` can
therefore never be written by a fresh perception — only by inheritance, which is one of the
reasons M5/S3 inherits sites through the atom map instead of re-perceiving, and why
`registry.catalog_drift` exists to watch for the disagreement.

**Costs today:** nothing in the assembly path, which never re-perceives. It would bite any
future code that perceives a stored complex directly, and it makes `n_perceived_donors`
meaningless for assembled structures.

**Fixing it moves** `n_perceived_donors` for every assembled structure and would change
`site_catalog` contents, so it needs its own fixture set and a `perception` version bump.

---

## 5. Δ/Λ handedness rests on an unverified absolute convention

`identity.isomers._chirality` assigns Δ to the right-handed propeller, from IUPAC's definition
of Δ as the right-handed helix plus the geometric fact that a right-handed helix turns
anticlockwise about its axis as it advances.

What is **tested**: mirror images get opposite labels, rotations do not move the label, and
the result does not depend on which end of the C3 axis you look from (flipping it swaps which
donor counts as "upper" *and* negates the axis; the sign changes cancel).

What is **not tested**: that the label matching a real Δ complex is `Delta` and not `Lambda`.
A systematic inversion would be invisible to every test here and visible in every report. One
check against a known crystal structure settles it; until then treat the absolute assignment
as provisional, while trusting the discrimination.

---

## 6. 16 of 39 stored structures have no site records

`data/registry.db`: 16 structures have no `site_catalog` rows and therefore `n_open_sites IS
NULL`. They predate M4's second half.

NULL is the correct value — nothing was computed — but a query that reads NULL as 0 will
report them as coordinatively saturated, which is the opposite of what an absent measurement
means. Not a defect in the writer; a trap for readers.

---

## 6b. θ_geom is a placeholder, not a calibration

`identity.conformers.DEFAULT_THETA_GEOM = 0.25 Å` is an unset number wearing a default. C2
could not be called against the M5 fixture set because the Kind-A (stochastic duplicate)
population has **identically zero** spread there — a join is a deterministic function of its
choice vector and absorbs the ligand's embedding noise entirely, so there is no peak for a
threshold to sit above. Measured: Kind A n=105, max 0.0000 Å; Kind B n=756, median 0.7952 Å,
in four discrete spikes rather than a distribution.

Mitigated rather than hidden: `cluster` takes `theta_geom` as an argument, no test depends on
the default, and the docstring says what it is. The fixture set that can set it is M7's
relaxed one, where two samples of one choice vector walk to *almost* the same minimum.

The energy window is in the same position — it gates "a bad geometry mistaken for a real
minimum" and every RAW construct is a non-minimum. Mechanism built and tested; number absent.

---

## 6c. The rigid-core definition treats a delocalised C–O as rotatable

`identity.conformers.rotatable` reuses `sites.model`'s rule — single, acyclic, non-aromatic —
so that "rigid" means one thing across the codebase. That rule is order-blind by design (D15
keeps bond order out of identity), so it cannot tell a carboxylate's delocalised C–O from an
ether's genuinely rotatable one, and cuts both.

Consequence: the rigid core of an acetate complex stops at the coordinating oxygen and
excludes the carboxylate carbon, so an acetate torsion registers as **exactly zero** core
RMSD. That happens to agree with §4.2's worked examples (a non-coordinating carboxyl torsion
is explicitly a *trivial* conformer that should not spawn L3), so the answer is right — but
it is right by accident, and the same cut would hide a real motion in an amide or an ester.

The error direction is safe: a smaller core makes the comparison more permissive, never less,
so it can fail to split two conformers and can never merge two that differ elsewhere.

---

## 7. The branch tree has more leaves than structures

The enumerator emits one leaf per branch and deliberately does not merge leaves that share an
L1 — cis and trans share one too, and telling a symmetry duplicate from a real isomer needs L2
and θ_geom clustering.

Measured with two partners, monodentate only:

| centre | degree | leaves | distinct L1 |
|---|---|---|---|
| CN 4 tetrahedral | 1 | 20 | 2 |
| CN 4 tetrahedral | 2 | 170 | 5 |
| CN 6 octahedral | 1 | 30 | 2 |
| CN 6 octahedral | 2 | **405** | **5** |

81 leaves per distinct L1 at CN 6 degree 2. `max_products` caps it and `BranchTree.capped`
reports the cap, which is the declared guard — but the ratio is the size of the job S6's
clustering has to do, and it grows with degree.

---

## 8. `vacancy_sites` normalises by hand, differently

`sites/model.py::vacancy_sites` divides by `max(norm, 1e-9)` where `geometry._linalg.unit`
falls back to the z axis. Identical for any real input; different for a zero-length direction
vector, where one returns a near-infinite vector and the other returns a valid axis.

Left out of the `_linalg` consolidation on purpose: swapping it is a behaviour change in the
degenerate case, not a de-duplication. Worth deciding which answer is wanted rather than
leaving two.

---

## Known states that are not issues

Recorded so they are not rediscovered as bugs:

- **99 of 191 geometries carry a choice vector with no digest.** Decided, not overlooked: the
  pre-S1 rows stay un-keyed as artifacts of the runs that made them.
  `registry.backfill_choice_digests` exists and is deliberately not called by `migrate`.
- **All 39 structures carry `l2=""` at `algo_l2="0-stub"`.** The classifier is built but no
  writer passes it a geometry, so nothing has been tagged. 18 of the 39 would gain a tag once
  the builder is wired — that wiring, and what happens to the rows predating it, is scheduled
  work, not a defect.
- **`runner.plan` refuses `degree > 1`.** The branch tree and `grow` are built; the run
  pipeline does not drive them yet (S4.1).
