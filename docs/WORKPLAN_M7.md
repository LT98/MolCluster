# M7 work plan — the gate, and what the archive can actually be gated against

`PLAN_implementation.md` §M7 says *what* the milestone is: the backends, `MethodSpec`, the
reference scheme and the relax runner are built, and what is left is an exit gate that has a
harness and **no recorded result from either of its stages**. This file is the *order of work*:
which check lands first, what each one can be run against on its own, and the calls that have
to be made before any code is written. It is scratch — when M7 lands, its conclusions go into
the design doc's Decision Ledger and the plan's changelog, and this file is deleted.

**The one-sentence scope:** after M7 every archived number this project intends to build on has
been re-derived through the current stack at a stated tolerance, the re-derivation is **recorded
where the next milestone can read it**, and the numbers the archive cannot support are named
rather than inherited.

**What changed about the milestone, before any of it is built.** The plan describes the gate as
"re-run the Ni/BTC and Fe/BTC cases and reproduce the rankings within noise", and
`scripts/regress_m7.py` splits that into a cheap stage 1 and a stage 2 it declines to build for
want of the archived candidate set. Measured against the files actually in `legacy/`, that
framing is wrong in both directions. The archived **rankings are in the repo** — per candidate,
not just as the five-row table in the report — so the gate has 92 comparisons available to it
and not 5. And the archive **disagrees with itself** before any of it is recomputed: the same
candidate appears in several scenarios with different energies, by up to 1.29 eV. The expensive
half of M7 is not running xTB. It is deciding what "reproduce" means against a target that is
not self-consistent, and that decision costs no compute at all.

---

## 1. What the gate is against — what exists, and what does not

| Artifact | Named by | In the repo? |
|---|---|---|
| `fe_btc_refs.json`, `ni_btc_refs.json` | `regress_m7.py` stage 1 | ✅ `legacy/` — four reference energies each |
| `fe_btc_xtb_results.csv` | `fe_btc_report.md` §Files | ✅ `legacy/` — **49 rows**, per-candidate `E_form_eV` |
| `ni_btc_xtb_results.csv` | `ni_btc_report.md` §Files | ✅ `legacy/` — **43 rows**, per-candidate `E_form_eV` |
| `solvation_raw.csv`, `solvation_summary.csv` | `solvation_report.md` §Files | ✅ `legacy/` — complete, refs included |
| `ni_btc_outputs/`, `fe_btc_outputs/` (`.xyz`) | both reports §Files | ❌ absent |
| `*_landscape.png` | all three reports §Files | ❌ absent |
| `ni_btc_results.csv` | — | ✅ present, `formation_eV` column **empty** |

**`stage_rankings`'s docstring is wrong about one of these and right about the other two.** It
says of `ni_btc_outputs/` and `ni_btc_xtb_results.csv` that "neither is in the repo". The CSV is
in the repo, at `legacy/ni_btc_xtb_results.csv`, with a populated `E_form_eV` for all 43
candidates; the same holds for Fe. What is genuinely absent is the **geometries** and the
figures, and the empty `formation_eV` in `ni_btc_results.csv` is real. Correcting that docstring
is S0's first item, because the whole shape of stage 2 follows from it: the missing input is
geometry, not ranking, so stage 2 is a re-derivation *with a per-candidate target*, not a
re-derivation into the dark.

The label carries the specification: `Ni2+_A_BTC_only_BTC[s1]_BTC[s1]_planar_q-4` is metal,
oxidation state, scenario, the ligand multiset with each ligand's donor-site pattern, the local
geometry, and the total charge. That is a `BuildSpec` in a string. What it does not carry is
which `sites` pattern legacy's `[sN]` corresponds to — the one mapping S4 has to establish.

---

## 2. What M7 inherits (do not rebuild any of this)

| Piece | Where | What it gives M7 |
|---|---|---|
| `XTBBackend` with ALPB | `energy/backends.py` | the solvation sweep's continuum, already keyed by the same solvent names the archive used (`water` / `ethanol` / `dmf`) |
| `high_spin_multiplicity` | `energy/backends.py` | the Fe(III) sextet convention as a rule rather than an env var |
| `minimal_multiplicity` | `energy/backends.py` | the Ni(II) convention — its own docstring says it exists *because the M7 regression has to reproduce those numbers* |
| `check_spin` | `energy/backends.py` | catches an unreachable multiplicity before xTB returns a plausible number for it |
| `MethodSpec` on every stored energy | `_types.py`, `energy/backends.py` | the provenance a recorded gate result needs to be re-checkable |
| `check_reference_quality` | `energy/reference.py` | the three rules (bare ion, naked polyanion, coordination change) that classify the archived equation |
| `reaction_balanced_energy(..., strict=False)` | `energy/reference.py` | evaluating the archived equation without endorsing it; the result carries `isodesmic=False` onward |
| `data/reference/*.tsv` + `scripts/check_cases.py` + `tests/test_reference_cases.py` | repo-wide | **the established pattern for a curated table the suite re-checks**: add a row, no code changes. M7's zero-QM stage is this pattern applied to energies |
| `RELAXATION_IS_EXECUTED` and `mode_status()` | `energy/relax.py` | the precedent for reporting "the backend imports" and "the backend will run" as two different claims |
| `tests/test_energy.py`'s `importorskip("tblite.ase")` | tests | how a tblite-dependent check skips rather than fails on a machine without it |

---

## 3. Five things the archive says, measured with no QM at all

Every number below came out of `legacy/` and the Python standard library. None of it needed
tblite, and none of it needs a workstation.

**(a) The archived arithmetic is exact.** All 43 Ni and all 49 Fe `E_form_eV` re-derive from
that row's `E_relax_eV` and the stored `*_refs.json`, with the ligand multiset parsed out of the
label: worst disagreement **0.001 eV** on both sets. So the archived formation equation is
`E_relax − E(M^q+) − Σ E(anion)` exactly as the reports state it, with no undocumented term.

**(b) The solvation sweep re-derives completely.** All 14 rows of `solvation_summary.csv` —
`Ef_node`, `Ef_ion`, `Ef_edta`, `nucleation_dG`, `seq_margin` — reproduce from `solvation_raw.csv`
under the report's cluster-continuum equation and its equal-weight mixture rule, worst
disagreement **0.005 eV**, which is the summary's own rounding. The solvation half of M7's exit
gate needs **zero** xTB: its inputs are archived single points, and everything the report
concludes from them is arithmetic over those.

**(c) The two archived runs disagree with themselves on identical species.** The Ni run and the
Fe run each computed the same three free ligand anions at the same level of theory:

| species | Ni run | Fe run | Δ |
|---|---|---|---|
| BTC³⁻ | −1241.6447 | −1241.6396 | 0.005 eV |
| EDDA²⁻ | −1100.2638 | −1100.2531 | 0.011 eV |
| **EDTA⁴⁻** | −1805.8501 | −1805.7784 | **0.072 eV** |
| EDTA⁴⁻ (ALPB water) | −1829.5229 | −1829.4355 | **0.087 eV** |

That is the archive's own run-to-run spread on identical species at identical settings, largest
on the floppiest ligand, and it is the only empirical anchor M7 has for the word "noise". Stage 1's tolerance is **1.0 eV**, roughly ten
times it. The tolerance is defensible as what it says it is — a wrong-charge/wrong-multiplicity
detector — and it is not evidence of agreement. C11 below is what to do about that.

**(d) The archived candidate set is not self-consistent.** Scenario is a filter over one
candidate set, so a composition appearing in two scenarios should carry one energy. It does not:

| | rows | distinct compositions | appear >once | **carry different energies** | worst |
|---|---|---|---|---|---|
| Ni | 43 | 21 | 18 | **9** | **1.287 eV** |
| Fe | 49 | 24 | 19 | **11** | 0.615 eV |

The cause is in the reports' own caveats — "the Ni run used a 90→50-step relaxation cap and Fe a
50→35-step cap for throughput" — and the step counts in the CSV confirm it: `BTC[s1]_BTC[s1]_planar_q-4`
is 64 steps in scenarios A/B and 44 in C/D. The consequence is not a caveat. **The new stack is
identity-keyed**: one composition is one L1, one row, one energy. It cannot reproduce both
copies, so the gate has to be told which copy is the target. That is C12.

**(e) One reported conclusion is a step-cap artifact, and one reported column does not
reproduce.**

* `ni_btc_report.md` reports the target node as "additive-independent (−36.4 eV in BTC-only;
  −36.1 eV in every additive scenario)". Those two numbers are **the same composition**,
  `BTC[s1]_BTC[s1]_planar_q-4`, relaxed for 64 steps and for 44. The node's additive-independence
  is a real and stronger result than the report claims — the structure is literally identical, so
  the residual 0.26 eV measures the step cap and nothing chemical.
* `fe_btc_report.md`'s "states **below** node" column reads 0 / 4 / 0 / 4. Counting strictly
  below the target node (−58.97) over the archived CSV gives 0 / **5** / 0 / **6**. The same
  counting rule reproduces the Ni report's column exactly (0 / 3 / 0 / 4), and the Fe report's
  spans, candidate counts and "gap of deepest state" all reproduce (−3.836 → −3.8, −3.372 → −3.4).
  Rounding the target to −59.0 fixes B and not D; using scenario D's own node copy (−58.75, the
  under-relaxed one) gives 6. **No rule tried reproduces both Fe cells**, and the ambiguity is
  (d)'s duplicate energies showing up in a conclusion.

Neither of these changes a headline verdict. Both are exactly the class of thing the gate exists
to catch, and both were caught for the price of reading a CSV.

---

## 4. What that makes the gate — three stages, three costs

The plan and the harness describe two stages. There are three, and the cheapest one is not built:

| | what it checks | needs | cost | state |
|---|---|---|---|---|
| **Stage 0** — re-derivation | the archived numbers follow from the archived inputs, and the reports follow from the numbers | nothing but Python | seconds, in the suite | **not built**; §3 is its result, run by hand |
| **Stage 1** — references | the current `XTBBackend` is the calculator that produced the archived reference energies | tblite | minutes | **built, Fe only, never recorded** |
| **Stage 2** — rankings | the current construction path reaches the archived candidates and orders them the same way | tblite + the enumeration | hours | **not built; raises and says why** |

Stage 0 comes first because it is free, because it is the half that can be a permanent
regression guard rather than something run once on a workstation, and because — as §3 shows — it
finds things. It also fixes the gate's target before any compute is spent reproducing a number
that two archived rows disagree about.

**The recording surface is the actual blocker, not the compute.** "No result from either is
recorded anywhere" is what keeps M7 open; a stage-1 pass that scrolls past in a terminal leaves
the milestone exactly where it is. S3 is small and is the thing that closes the milestone.

---

## 5. The reference scheme is not cosmetic — and it moves exactly the comparison the claim rests on

The archived equation trips all three of `check_reference_quality`'s rules: a bare `Fe³⁺` / `Ni²⁺`
ion, `BTC³⁻` and `EDTA⁴⁻` as naked polyanions, and a metal–donor bond count that goes from zero to
six across the arrow. `strict=False` is needed to evaluate it, and the value carries
`isodesmic=False`. D17 already says so.

What follows from the arithmetic, and is worth stating before stage 2 is designed:

* Within a **fixed ligand multiset**, the reference terms are a constant. Two `EDTA[s4]` candidates,
  or the two BTC nodes, differ only in `E_relax`, so **their ordering is scheme-invariant** — any
  balanced scheme ranks them identically.
* Across **different multisets** it is not, and the reports' headline is a cross-multiset
  comparison: the Ni–EDTA chelate 8.5 eV below the Ni–BTC node. That gap is
  `[E(NiEDTA) − E(EDTA⁴⁻)] − [E(Ni(BTC)₂) − 2E(BTC³⁻)]`, and it is carried by the free-anion
  reference terms that D17 refuses.
* The reports already know this and hedge it the long way round. `E_form_per_bond` exists
  precisely because "total E_form rewards making more bonds", and per bond the ordering flips
  (node −18.2, chelate −7.5 eV/bond); the Ni report then defends the sink conclusion by appeal to
  Ni–EDTA being "a genuinely benchmark-stable complex", which is evidence from outside the
  calculation. An isodesmic ligand-exchange equation conserves the dative count by construction,
  which is the same hedge done once, in the arithmetic, instead of in prose.
* And the solvation run already reversed this verdict once: gas-phase −8.5 eV became −0.3 to
  +0.5 eV in every solvent. D17's whole justification is that this reversal happened.

So stage 2 compares **orderings within a multiset** as a reproduction claim, and reports the
cross-multiset ordering **twice** — once under the archived scheme at `strict=False`, once under
an isodesmic exchange equation — as a divergence to be read, never averaged. The value of the
gate is in the divergence.

---

## 6. The calls to make before writing code

**C11 — what does "reproduce within noise" mean?** Stage 1's 1.0 eV is ten times the archive's own
0.09 eV worst same-species spread (§3c), so a pass is nearly uninformative about agreement. Two
choices, and they are not the same test: a **per-species absolute tolerance** (is this the same
calculator?) or an **ordering statistic** (does the landscape rank the same way?).
*Leaning:* both, named separately and never conflated. Keep the 1.0 eV absolute check as the
calculator-identity test it already is and say so in its output; add an ordering criterion for
stage 2, because ordering is what M8 consumes and values are what the reference scheme moves.
Calibrate the absolute tolerance from §3c's distribution rather than from a round number — the
same move θ_geom got, and the M6 plan's "a tolerance becomes folklore" risk applies here word for
word.

**C12 — which copy of a duplicated archived candidate is the target?** Forced by §3d, and forced
*now* because it decides what stage 2 is scored against. Four candidate rules: lowest energy; the
`converged=True` copy; the copy from the highest step cap; or refuse and re-run. *Leaning:* the
highest step cap, because the disagreement is a convergence artifact with a known direction — a
capped relaxation is an upper bound, so the deepest copy is the most relaxed and the least
arbitrary. Whatever is chosen, the dedup must be **recorded as a column in the curated table**,
not applied in code, so a reader can see which 21 (Ni) and 24 (Fe) rows are the target and why.

**C13 — does M7 close on stages 0 and 1, with stage 2 as M8's entry cost?** `PLAN_implementation.md`
§7 asks exactly this and leaves it open. The plan sizes M7's remainder as **S**; stages 0, 1 and
the recording are S, and stage 2 is L. *Leaning:* close M7 on stages 0 + 1 + a recorded result,
and move stage 2 into M8 **as a written entry cost with its own exit line**, not as a deletion —
M8 takes these numbers as given, so the re-derivation has to happen before a claim rests on them,
and M8 is where the claim is. The cost of this choice: M7's headline becomes "the archive is
re-derivable and the calculator matches" rather than "the rankings reproduce", and
`PLAN_implementation.md` §5's promise — "trust the relative energies attached to any of it" —
has to be narrowed to say so honestly.

Rule for all three, from the plan §3: when a gate is called, **write the resolution into the
Decision Ledger with a new D-number and a changelog line**, then code against it.

---

## 7. The slices

### S0 — correct the record, and make the three calls · **S**

`stage_rankings`'s docstring (§1) — the Ni and Fe result CSVs are in `legacy/`; what is missing is
the geometries. The rest of that docstring stands, including the empty `formation_eV` in
`ni_btc_results.csv` and the re-derivation argument. C11, C12 and C13 get their D-numbers and
changelog lines. Nothing downstream is worth building against a docstring that says the target
does not exist.

*Exit:* the docstring matches the filesystem; three ledger entries; this file's §3 numbers cited
from the archive rather than from here.

### S1 — stage 0: the zero-QM re-derivation, in the suite · **S**

`data/reference/energy_cases.tsv` as **the interface for adding a case**, in the shape
`ligand_cases.tsv` and `node_cases.tsv` already established: one row per archived claim, with its
source file, the equation it asserts, the expected value and the tolerance. `scripts/check_energy_cases.py`
evaluates it and `tests/test_energy_cases.py` runs it, so a case added here is a permanent guard
rather than something that scrolled past once. Three families of row:

1. per-candidate `E_form` from `E_relax` + refs (92 rows, §3a);
2. the solvation summary from the raw single points (14 rows × 5 columns, §3b);
3. the reports' own tables — spans, candidate counts, below-node counts, deepest-gap — against the
   CSVs (§3e).

Family 3 is the one that fails today, on two Fe cells. Land it **failing**, with the measured
counts in the row, then resolve it under C12 and fix whichever of the report or the count is
wrong in the same commit.

*Exit:* `python -m pytest tests/test_energy_cases.py` passes on a machine with no scientific stack
at all; the Ni node's additive-independence is recorded as the step-cap identity it is; the Fe
below-node column reproduces or the report is corrected with the archived numbers beside it.

### S2 — stage 1, both metals, both conventions · **S**

`regress_m7.py` hardcodes `METAL, METAL_CHARGE = "Fe", 3` and `legacy/fe_btc_refs.json`.
`legacy/ni_btc_refs.json` is there with the same shape (its metal key is `Ni`, not `M`) and the
Ni run's minimal-spin convention is already `minimal_multiplicity`, which says in its own
docstring that it exists for this. So this is a parametrisation, not a build: `--metal ni|fe`,
the spin convention carried **per archived run** rather than per repo, and the reference-file key
handled rather than assumed.

*Exit:* both reference sets recompute; the per-species tolerance is C11's calibrated number with
the 1.0 eV calculator-identity check reported separately and labelled as such; the output states
which spin convention each run used.

### S3 — record the result where the next milestone reads it · **S**

The milestone's actual blocker. A gate result needs: the `MethodSpec` (`tblite` version, method,
solvent, charge, multiplicity), the machine, the date, the tolerance and the per-species
deltas — `--json` already emits most of it. Two surfaces, and they are not interchangeable:
`docs/reports/m7_regression_report.md` for the reading (what was run, what it means, what it does
not cover), and a committed JSON for the re-checking. Both, or the milestone reopens the first
time someone asks what the numbers were.

*Exit:* `PLAN_implementation.md`'s "no recorded result" sentence is false, and the report says in
its own words which of stage 0 / 1 / 2 it covers.

### S4 — stage 2: the re-derivation · **L** *(C13 decides whether here or at M8)*

The candidate set comes out of the labels (§1): metal, oxidation state, ligand multiset, local
geometry, charge. The work is in three places, none of them the xTB:

* **the `[sN]` mapping** — legacy's donor-site pattern indices onto this project's `sites`
  perception. Establish it against `ligand_cases.tsv`, which already pins what each ligand
  perceives, and refuse a label whose pattern has no counterpart rather than guessing one.
* **the enumeration** — a `BuildSpec` per scenario through `construct` and the relax runner. It is
  hours of compute and it is resumable work, so it belongs behind its own flag and in the jobs
  table, not inside a script function. `RELAXATION_IS_EXECUTED`'s cautionary note applies: a
  pre-flight that passes without an executor is how rev 16 wrote 281 unrelaxed structures under a
  label saying otherwise.
* **the scoring, twice** — §5. Archived equation at `strict=False`, and an isodesmic exchange
  equation, reported side by side. Within-multiset orderings are the reproduction claim;
  the cross-multiset ordering is the divergence to read.

*Exit:* every archived composition is either reached by the construction path or **refused with a
reason**, the way M6's failed route is refused with its number; within-multiset orderings match;
the two scoring schemes are reported side by side and the places they disagree are named.

### S5 — the solvation sweep · **M**

The ALPB path exists in `XTBBackend` and the archive's solvent keys are already its keys. S1
covers the *analysis*; this covers the *energies* — 39 species × 3 solvents = **117 ALPB single
points** reproducing `solvation_raw.csv`, at which point the two decision metrics follow from
S1's arithmetic for free.
Worth doing because the solvation run is the evidence for D17 and the one that reversed a
published verdict, so it is the number M8 leans on hardest.

*Exit:* the raw ALPB single points reproduce within C11's tolerance; `nucleation_dG` and
`seq_margin` follow from S1's already-tested arithmetic rather than from new code.

---

## 8. Risks, with the escape hatch named in advance

| Risk | Early signal | Hatch |
|---|---|---|
| **The gate passes and means nothing** | stage 1 green at 1.0 eV against a 0.09 eV archive spread | C11 — two named criteria, calibrated from §3c's distribution; never report "reproduces" off the absolute check alone |
| **Stage 2 expands into M8** | the `[sN]` mapping turns into perception work, or the enumeration needs construction features that do not exist | C13 is this hatch, taken in advance: close M7 on stages 0+1, carry stage 2 as M8's written entry cost with its own exit line |
| **A tolerance becomes folklore** | a threshold picked from a figure in this file | §3's numbers are points, not populations — the same warning M6's plan carries, and the same answer: calibrate from the distribution and pin the separation |
| **The archive is treated as ground truth** | a stage-2 mismatch is fixed by adjusting the new stack | the archive disagrees with itself by up to 1.29 eV (§3d) and reversed its own verdict once solvation was added. It is a **comparison**, not an oracle; a divergence is a finding, and which side is wrong is the question the gate asks |
| **No tblite anywhere** | neither machine can run stages 1, 2 or 5 | stage 0 is the whole of S1 and needs no scientific stack at all — the milestone's cheapest half stays runnable, and stages 1/2/5 skip the way `tests/test_energy.py` already skips |
| **Recording is deferred to "when it all passes"** | stage 1 run, result in a terminal | S3 is sized S and ordered before S4 for this reason. A stage-1 result recorded alone still closes something |

---

## 9. Bookkeeping when M7 lands

- Design doc: D-numbers for **C11**, **C12** and **C13** with a changelog line each — the next free
  numbers after M6's (D23 onward if M6 lands D20–D22). D17 is unchanged and is what §5 codes
  against.
- `PLAN_implementation.md`: the "Partly done" row and the §M7 section move to
  `archive/PLAN_completed.md` if C13 closes the milestone; §5's "trust the relative energies
  attached to any of it" is narrowed to what was actually gated; §7 item 1 is struck and replaced
  by whatever C13 decides about stage 2; §6's golden/regression bullet points at
  `tests/test_energy_cases.py`.
- `CODE_ARCHITECTURE.md`: an `energy_cases.tsv` row beside the other curated tables; the §"where
  to read about energies" pointer gains the recorded report.
- `BUGS.md`: nothing of M7's is a defect today. If S1 family 3 resolves as *the report is wrong*,
  the correction goes into the report with the archived numbers beside it, not into `BUGS.md` —
  a report is a record of a run, and correcting it is an erratum, not a bug fix.
- `data/reference/NOTES.md`: `energy_cases.tsv` gets its provenance paragraph, saying plainly that
  its targets are **archived results, not literature**, and that two of them disagree with
  themselves.
- Delete this file.

---

## 10. Order, at a glance

```
S0 record + C11/C12/C13 ─► S1 stage 0 (zero QM, in the suite) ─┬─► S2 stage 1 ─► S3 record ─► [C13]
                                                               └─► S5 solvation energies      │
                                                                                              ▼
                                                                              S4 stage 2 (here, or M8)
```

S1 needs nothing but S0's calls and can be done on any machine, including one with no scientific
stack — which is the argument for doing it first. S2 and S5 both need tblite and are independent
of each other. S3 can record a stage-1-only result and should, rather than waiting for a complete
one. S4 is gated on C13 and is the only slice that is not S or M.
