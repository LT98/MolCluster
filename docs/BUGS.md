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
| [B1](#b1) | honesty | `sites/perception.py` | A coordinated donor is perceived as no donor at all |
| [B2](#b2) | honesty | `identity/keys.py` | `l2_isomer_tag` is `''` for all 39 structures, so cis and trans are one row |
| [B3](#b3) | honesty | registry data | 16 structures have no `site_state`; `n_open_sites` is NULL, not 0 |
| [B4](#b4) | cosmetic | `ui/static/runs.html` | "What was attempted" renders `molecule undefined · undefined×0-dentate` |
| [B5](#b5) | undecided | `ui/static/builder.html` | The builder posts `spec_version: 1` and no test covers the migration it relies on |
| [B6](#b6) | undecided | `ui/static/index.html` | The registry page never refreshes, so its counts go stale silently |
| [B7](#b7) | blocked | `scripts/ingest.py` | The legacy `.xyz` corpus cannot be ingested: no charge, no multiplicity |

---

## B1

**A coordinated donor is perceived as no donor at all.** `honesty` ·
`sites/perception.py`

Perception counts a metal as an ordinary heavy neighbour. A coordinated aqua oxygen therefore
fails the donor test, no `site_catalog` row is written for it, and the
`SiteStatus.OCCUPIED` row that should exist on every geometry never appears. The sites that
are *occupied* are exactly the ones that vanish.

Not the same seam as the resonance one (closed — see the archive): that was perception reading
bond order, this is perception reading coordination. They were filed together and only one of
them was fixed.

**Why it is still open:** closing it moves `n_perceived_donors` for every assembled structure,
so it wants its own fixture set and a decision about whether stored catalogs are re-derived or
versioned.

---

## B2

**`l2_isomer_tag` is `''` everywhere, so cis and trans collapse into one row.** `honesty` ·
`identity/keys.py`

Measured 2026-09-14: 39 of 39 structures carry an empty L2. The stub is by design (M2 fixed
the signature, M5 fills the body), but the *consequence* is a live misreport — the registry
currently claims two isomers are one structure.

**The decision this blocks:** filling L2 in retroactively **splits identities**. Backfill
versus version bump has to be called before M5 touches `l2_isomer_tag`, because after that
point the split is happening either way and the only question is whether it was chosen.

---

## B3

**16 structures have no `site_state`.** `honesty` · registry data

Measured 2026-09-14: 16 of 39 structures were built before M4's second half, so they have a
catalog and no state, and `n_open_sites` is NULL for them.

NULL is the correct value — absent is not zero (ground rule 9 / invariant 9) — and the bug is
that nothing stops a consumer from treating it as zero and reporting a fully-occupied
structure. Either backfill the state for the old rows or make NULL loud at the read side.

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
