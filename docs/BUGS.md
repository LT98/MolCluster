# Bug checklist

One line per bug. **New bugs are filed as GitHub issues** (`gh issue create --label bug`) and
get a `[ ]` line here linking the issue: an issue is visible from every branch, a section in a
file only from the branches that have it. Put the severity and the module in the issue body,
using the vocabulary below. The `B…` entries were filed before 2026-10-19; their diagnosis and
fix direction are in [`BUGS_details.md`](BUGS_details.md).

**Fixing one** is a one-line change here — tick `[ ]` to `[x]` — plus, for a `B…` entry, a
one-line `**Fixed (rev N, PR #M).** …` note under its heading in `BUGS_details.md`. Nothing is
moved. Entries fixed before this layout are in [`archive/BUGS_resolved.md`](archive/BUGS_resolved.md).

Scope: things that behave wrongly or report wrongly. *Unbuilt* functionality is not a bug —
that is a milestone, and it lives in [`PLAN_implementation.md`](PLAN_implementation.md).
A stub that raises loudly is working as designed (ground rule 8).

**Severity** — `correctness` a stored value or a claim is wrong · `honesty` the code is right
but what it reports is misleading · `cosmetic` it looks wrong and misleads nobody ·
`undecided` known behaviour that nobody has explicitly chosen.

## Open

- [ ] **B2** · honesty · `runner.py` — Every stored row's `l2_isomer_tag` is `''`, so cis and trans are one row. The classifier exists; the run pipeline declines to use it (**D19**, **D22**) ([details](BUGS_details.md#b2))
- [ ] **B3** · honesty · registry data — 16 structures have no `site_state`; `n_open_sites` is NULL, not 0. **Self-healing (D19)** ([details](BUGS_details.md#b3))
- [ ] **B5** · undecided · `ui/static/builder.html` — The builder posts `spec_version: 1` and no test covers the migration it relies on ([details](BUGS_details.md#b5))
- [ ] **B6** · undecided · `ui/static/index.html` — The registry page never refreshes, so its counts go stale silently ([details](BUGS_details.md#b6))
- [ ] **B7** · blocked · `scripts/ingest.py` — The legacy `.xyz` corpus cannot be ingested: no charge, no multiplicity ([details](BUGS_details.md#b7))
- [x] **B8** · correctness · `descriptors/ease.py` — Every sp3 amine donor reads as sterically blocked, so no amine can be joined ([details](BUGS_details.md#b8), #13)
- [ ] **B9** · correctness · `identity/conformers.py` — The rigid-core rule cuts a delocalised carboxylate C–O as if it were rotatable ([details](BUGS_details.md#b9), #20)
- [ ] **B10** · undecided · `identity/isomers.py` — Δ/Λ is self-consistent but its absolute assignment is unverified ([details](BUGS_details.md#b10), #17)
- [ ] **B11** · cosmetic · `sites/model.py` — `vacancy_sites` normalises by hand, differently from `_linalg.unit` ([details](BUGS_details.md#b11), #22)
- [ ] **B14** · correctness · `sites/perception.py` — Pyrazolate's two equivalent N type differently, and one of them cannot bridge ([details](BUGS_details.md#b14))
- [ ] **B15** · correctness · `sites/perception.py` — A bare oxide `[O-2]` types as `hydroxide_O`, so an oxo reports a hydroxide's pKa ([details](BUGS_details.md#b15))
- [ ] **B16** · correctness · `scripts/run_spec.py` — `--store` is ignored: spawned workers re-resolve `data_root()`, so blobs land in the default store ([details](BUGS_details.md#b16))
- [ ] **B17** · honesty · `ui/static/index.html` — Changing the geometry does not re-render the detail panel, so `method` and `converged` go stale ([details](BUGS_details.md#b17))
- [ ] **B19** · honesty · `ui/static/graph.html` — One structure drawn at two heights by two routes, with the reason — a different basis — only in the legend ([details](BUGS_details.md#b19))
- [ ] **B27** · correctness · `sites/model.py` — A vacancy's torsion reference flips sign between cis vertices, so `join(torsion_well=0)` puts cis ligands on opposite faces ([details](BUGS_details.md#b27))
- [ ] **B29** · correctness · `identity/isomers.py` — `isomer_tag` keys each metal centre by canonical index, which is not stable for symmetry-equivalent atoms — one polynuclear species could get two L2 tags (latent until B2 is closed) ([details](BUGS_details.md#b29))
- [x] **#70** · correctness · `energy/dft.py` — `DFTBackend` ignores `MOFSBU_DEVICE`: GPU4PySCF spreads over every visible GPU ([issue](https://github.com/LT98/MolCluster/issues/70))
- [ ] **#53** · undecided · `ui/` — Searching the registry for a displayed label (`tHQ`, chloride) returns no entries ([issue](https://github.com/LT98/MolCluster/issues/53))

## Resolved

- [x] **B20** — a charged co-ligand's charge was left out of the complex's net charge ([details](archive/BUGS_resolved.md#b20--a-charged-co-ligands-charge-was-left-out-of-the-complexs-net-charge-))
- [x] **B28** — a later build's rows read through the first build's numbering ([details](archive/BUGS_resolved.md#b28--a-later-builds-rows-read-through-the-first-builds-numbering-))
- [x] **B26** — `database is locked` with many build workers ([details](archive/BUGS_resolved.md#b26--database-is-locked-with-many-build-workers-))
- [x] **B4** — the run inspector's "what was attempted" read `undefined` ([details](archive/BUGS_resolved.md#b4--the-run-inspectors-what-was-attempted-read-undefined-))
- [x] **B21** — an octahedral centre with chelates and empty vertices was refused ([details](archive/BUGS_resolved.md#b21--an-octahedral-centre-with-chelates-and-empty-vertices-was-refused-))
- [x] **B25** — a step was claimed while the rung it joins onto was still being built ([details](archive/BUGS_resolved.md#b25--a-step-was-claimed-while-the-rung-it-joins-onto-was-still-being-built-))
- [x] **B24** — a grow tried one orientation and rejected what `place` built easily ([details](archive/BUGS_resolved.md#b24--a-grow-tried-one-orientation-and-rejected-what-place-built-easily-))
- [x] **B23** — a "deprotonation" edge could also change which atom binds the metal ([details](archive/BUGS_resolved.md#b23--a-deprotonation-edge-could-also-change-which-atom-binds-the-metal-))
- [x] **B22** — an anion entered as one was registered as neutral ([details](archive/BUGS_resolved.md#b22--an-anion-entered-as-one-was-registered-as-neutral-))
- [x] **B18** — a `grow` step wrote an edge without the ligand it added ([details](archive/BUGS_resolved.md#b18--a-grow-step-wrote-an-edge-without-the-ligand-it-added-))
- [x] **B13** — a bare hydroxide perceived zero donors ([details](archive/BUGS_resolved.md#b13--a-bare-hydroxide-perceived-zero-donors-))
- [x] **B12** — `metal_block` gave a metal a formal charge nothing else does ([details](archive/BUGS_resolved.md#b12--metal_block-gave-a-metal-a-formal-charge-nothing-else-does-))
- [x] **Runner, issue #32** — a GPU relaxed while one core built for it ([details](archive/BUGS_resolved.md#runner-issue-32--a-gpu-relaxed-while-one-core-built-for-it-))
- [x] **Runner, issue #32 follow-up** — a run killed with its shell was ongoing for ever ([details](archive/BUGS_resolved.md#runner-issue-32-follow-up--a-run-killed-with-its-shell-was-ongoing-for-ever-))
- [x] **UI, rev 24** — seven usability defects, all closed ([details](archive/BUGS_resolved.md#ui-rev-24--seven-usability-defects-all-closed))
- [x] **Identity, M4 follow-up** — `site_catalog` was not a pure function of identity ([details](archive/BUGS_resolved.md#identity-m4-follow-up--site_catalog-was-not-a-pure-function-of-identity-))
- [x] **Registry, M4 follow-up** — `put_sites` deleted the state it had just written ([details](archive/BUGS_resolved.md#registry-m4-follow-up--put_sites-deleted-the-state-it-had-just-written-))
- [x] **Perception, pre-M5** — a donor was deleted at the moment it became occupied ([details](archive/BUGS_resolved.md#perception-pre-m5--a-donor-was-deleted-at-the-moment-it-became-occupied-))
- [x] **#69** (was B30) · correctness — constructs put an O–H hydrogen 1.3–2.0 Å from the metal, and QC passed them; fixed in rev 60 ([issue](https://github.com/LT98/MolCluster/issues/69))
