# M7 work plan — closing the energy gate

M7's code is built. Its exit gate has a harness and no recorded result, and that is all that
keeps the milestone open. This file lists the work that closes it. Delete the file when M7
lands; its conclusions go to the Decision Ledger and the plan's changelog.

`PLAN_implementation.md` §M7 says what the milestone is. This file says what to do.

---

## 1. The gate today

| Stage | Checks | Needs | State |
|---|---|---|---|
| 0 — re-derivation | the archived numbers follow from the archived inputs | Python only | **does not exist** |
| 1 — `--refs` | `XTBBackend` reproduces the archived reference energies | tblite, minutes | built for Fe only; never recorded |
| 2 — `--rankings` | the construction path reaches the archived candidates and orders them the same | tblite + enumeration, hours | raises `NotBuiltYet` |

Stage 0 is new in this plan. §3 and §4 are its output, run by hand. It needs no scientific stack,
so it goes in the test suite. Run it first: it fixes what stages 1 and 2 are measured against.

---

## 2. Archived artifacts: what is in `legacy/`, what is not

| File | In repo | Contents |
|---|---|---|
| `fe_btc_refs.json` | yes | 4 gas-phase reference energies (Fe³⁺, BTC³⁻, EDTA⁴⁻, EDDA²⁻) |
| `ni_btc_refs.json` | yes | same 4 for the Ni run; metal key is `Ni`, not `M` |
| `fe_btc_xtb_results.csv` | yes | 49 candidates: `E_sp_eV`, `E_relax_eV`, `converged`, `steps`, `E_form_eV` |
| `ni_btc_xtb_results.csv` | yes | 43 candidates, same columns |
| `solvation_raw.csv` | yes | 39 species × 3 solvents of ALPB single points, references included |
| `solvation_summary.csv` | yes | 14 rows, 5 derived columns |
| `ni_btc_results.csv` | yes | `formation_eV` column is empty |
| `ni_btc_outputs/`, `fe_btc_outputs/` (`.xyz`) | **no** | the candidate geometries |
| `*_landscape.png` | **no** | the figures |

The candidate label is a build specification. `Ni2+_A_BTC_only_BTC[s1]_BTC[s1]_planar_q-4` gives
metal, oxidation state, scenario, ligand multiset with each ligand's legacy donor-site pattern,
local geometry, total charge. It does not say which `sites` pattern legacy's `[sN]` means; S4
establishes that mapping.

---

## 3. Measurements taken from `legacy/`

All of these were produced with the Python standard library against tracked files. No xTB.

**(a) `E_form_eV` follows from the same row's `E_relax_eV`.** For all 43 Ni and 49 Fe rows,
`E_relax_eV − E(metal) − Σ E(ligand anion)`, with the ligand multiset read off the label and the
reference energies from `*_refs.json`, reproduces the stored `E_form_eV` to within 0.001 eV. The
archived equation is what the reports say it is, with no undocumented term.

**(b) The solvation summary is derivable from the raw single points.** All 14 rows and all 5
derived columns of `solvation_summary.csv` reproduce from `solvation_raw.csv` under the report's
formation equation and its equal-weight mixture rule. Worst disagreement 0.005 eV, which is the
summary's rounding.

**(c) The two archived runs disagree on identical species.**

| species | Ni run | Fe run | difference |
|---|---|---|---|
| BTC³⁻ gas | −1241.6447 | −1241.6396 | 0.005 eV |
| EDDA²⁻ gas | −1100.2638 | −1100.2531 | 0.011 eV |
| EDTA⁴⁻ gas | −1805.8501 | −1805.7784 | 0.072 eV |
| EDTA⁴⁻ ALPB | −1829.5229 (water) | −1829.4355 (water) | 0.087 / 0.093 / 0.096 eV (water / ethanol / DMF) |

0.096 eV is the largest same-species disagreement in the archive. Stage 1's tolerance is 1.0 eV.

**(d) The archived candidate set contains duplicates that disagree.** A composition appearing in
more than one scenario should carry one energy.

| | rows | distinct compositions | appear more than once | energies differ | largest difference |
|---|---|---|---|---|---|
| Ni | 43 | 21 | 18 | 9 | 1.287 eV |
| Fe | 49 | 24 | 19 | 11 | 0.615 eV |

Cause: the relaxation step cap differed between scenario groups, which both reports state as a
caveat. `BTC[s1]_BTC[s1]_planar_q-4` ran 64 steps in scenarios A/B and 44 in C/D.

Consequence: one composition is one L1 and one row in this registry, so a re-run produces one
energy per composition and cannot match both copies. C12 decides which copy is the target.

---

## 4. Errors found, and what to do about each

| # | Where | Error | Fix |
|---|---|---|---|
| 1 | `scripts/regress_m7.py:125` | the docstring says of `ni_btc_outputs/` and `ni_btc_xtb_results.csv` that "neither is in the repo". The CSV is in `legacy/`. | S0: correct the sentence to say the geometries are missing and the results CSV is not. The rest of that docstring is accurate, including the empty `formation_eV` in `ni_btc_results.csv`. |
| 2 | `docs/reports/fe_btc_report.md`, "states below node" column | reads 0 / 4 / 0 / 4. Counting strictly below −58.97 over `fe_btc_xtb_results.csv` gives 0 / 5 / 0 / 6. The same rule reproduces the Ni report's column exactly (0 / 3 / 0 / 4), and Fe's spans, candidate counts and deepest-gap figures (−3.836 → −3.8, −3.372 → −3.4) all reproduce. Rounding the target to −59.0 fixes scenario B and not D; using scenario D's own node copy (−58.75) gives 6. | S1: encode the counting rule and the measured counts as a case, land it failing, then resolve under C12 and correct either the report or the count in the same commit. |
| 3 | `docs/reports/ni_btc_report.md`, "additive-independent (−36.4 eV in BTC-only; −36.1 eV in every additive scenario)" | the two numbers are the same composition, `BTC[s1]_BTC[s1]_planar_q-4`, at 64 and 44 relaxation steps. The 0.26 eV gap measures the step cap. | S1: record the node's additive-independence as composition identity, and drop the 0.26 eV from any statement about additive sensitivity. |

Neither report error changes a headline verdict.

---

## 5. Decisions to make before writing code

Each needs a Decision Ledger entry with a D-number and a changelog line before it is coded
against. The next free numbers are after M6's (D23 onward, if M6 lands D20–D22).

**C11 — what "reproduce within noise" means.**
The 1.0 eV tolerance is 10× the archive's largest same-species disagreement (§3c), so passing it
shows the calculator is the same one, not that the numbers agree.
*Recommendation:* keep two criteria and report them separately. (i) An absolute per-species
tolerance, labelled in the output as a calculator-identity check, with its value calibrated from
§3c's distribution rather than left at 1.0. (ii) An ordering criterion for stage 2, since
orderings are what M8 consumes.
*Cost:* stage 1 can no longer be reported as "the numbers reproduce". It reports "same
calculator".

**C12 — which copy of a duplicated archived candidate is the target (§3d).**
*Options:* lowest energy · the `converged=True` copy · the copy from the highest step cap ·
refuse and re-run everything.
*Recommendation:* highest step cap. A capped relaxation is an upper bound, so the deepest copy is
the most relaxed, and the rule is stated rather than chosen per row.
*Required either way:* the choice is a column in `energy_cases.tsv`, not logic in a script, so a
reader can see which 21 Ni and 24 Fe rows are the target.

**C13 — does stage 2 close M7, or become M8's entry cost?**
`PLAN_implementation.md` §7 asks this and leaves it open. Stages 0, 1 and the recording are each
S. Stage 2 is L, so the plan's "S remaining" only holds if stage 2 moves.
*Recommendation:* close M7 on stages 0 + 1 + a recorded result. Move stage 2 into M8 with its own
exit line — M8 is where these numbers first carry a claim.
*Cost:* `PLAN_implementation.md` §5 currently promises that after M7 you can "trust the relative
energies attached to any of it". That must be narrowed to what was gated.

---

## 6. What the reference scheme does to stage 2

The archived equation trips all three rules in `check_reference_quality`: a bare metal ion, two
naked polyanions, and a metal–donor bond count going from 0 to 6 across the arrow. It needs
`strict=False` and carries `isodesmic=False`. That is D17, already decided.

What follows for stage 2:

- Within one ligand multiset the reference terms are a constant, so the ordering there is the same
  under any balanced scheme. Use those orderings as the reproduction test.
- Across different multisets they are not. The reports' headline — the Ni–EDTA chelate 8.5 eV
  below the Ni–BTC node — is `[E(NiEDTA) − E(EDTA⁴⁻)] − [E(Ni(BTC)₂) − 2E(BTC³⁻)]`, carried by the
  free-anion terms D17 refuses. The reports' own `E_form_per_bond` column reverses this ordering
  (node −18.2, chelate −7.5 eV/bond), and the solvation run moved the same gap from −8.5 eV to
  −0.3…+0.5 eV.

So stage 2 scores each candidate twice — archived equation at `strict=False`, and an isodesmic
ligand-exchange equation — and reports both. Where they disagree, report the disagreement. Do not
average them and do not pick one.

---

## 7. Slices

### S0 — correct the record and make the calls · **S**

1. Fix `scripts/regress_m7.py:125` per §4 row 1.
2. Write Decision Ledger entries for C11, C12 and C13, with a changelog line each.

**Done when:** the docstring matches the filesystem and three D-numbers exist.

### S1 — stage 0: the zero-QM re-derivation, in the suite · **S**

1. Add `data/reference/energy_cases.tsv`. Columns: `name`, `kind`, `source_file`, `subject`,
   `metric`, `expected`, `tol`, `note` — where `subject` identifies the row being checked (a
   candidate label, a metal+solvent system, or a report cell) and `metric` names the quantity.
   No column holds an expression to evaluate; the checker dispatches on `kind`. Header comment in
   the style of `node_cases.tsv`, stating that these targets are archived results rather than
   literature, and that some of them disagree with each other (§3d).
2. Add `scripts/check_energy_cases.py` with `load_cases()` and `evaluate(case)`, matching
   `scripts/check_cases.py`'s interface.
3. Add `tests/test_energy_cases.py`, parameterised over the rows, matching
   `tests/test_reference_cases.py`.
4. Populate three kinds of row:
   - `formation` — 92 rows, one per archived candidate: `E_relax_eV` minus the refs reproduces
     `E_form_eV` (§3a).
   - `solvation` — 14 rows: `Ef_node`, `Ef_ion`, `Ef_edta`, `nucleation_dG`, `seq_margin`
     reproduce from `solvation_raw.csv` (§3b).
   - `report` — the tables in `ni_btc_report.md` and `fe_btc_report.md`: per-scenario candidate
     count, span, states below node, deepest gap.
5. Land the two failing Fe `report` rows with their measured counts, then resolve under C12 and
   correct the report or the count (§4 row 2).
6. Correct §4 row 3 in `ni_btc_report.md`.

**Done when:** `python -m pytest tests/test_energy_cases.py` passes on a machine with no
scientific stack installed.

### S2 — stage 1 for both metals · **S**

`scripts/regress_m7.py` hardcodes `METAL, METAL_CHARGE = "Fe", 3` (line 53) and
`ARCHIVED = legacy/fe_btc_refs.json` (line 49), and reads the metal energy as `archived["M"]`
(line 80). `ni_btc_refs.json` keys it `Ni`.

1. Add `--metal ni|fe`; select the reference file and the metal key from it.
2. Carry the spin convention per archived run, not per repo: Fe used a sextet
   (`high_spin_multiplicity`), Ni used minimal spin (`minimal_multiplicity`, which says in its own
   docstring that it exists for this regression). Print which was used.
3. Apply C11: report the calibrated per-species tolerance and the 1.0 eV calculator-identity check
   as two separate lines.

**Done when:** `--refs --metal ni` and `--refs --metal fe` both run, and the output states each
run's spin convention.

### S3 — record the result · **S**

The milestone is blocked on this, not on compute.

1. Commit the `--json` output for each stage that has been run.
2. Write `docs/reports/m7_regression_report.md`: which stages ran, on which machine, on what date,
   the `MethodSpec` (tblite version, method, solvent, charge, multiplicity), the per-species
   deltas, and which stages are **not** covered.

**Done when:** `PLAN_implementation.md`'s "no recorded result" sentence is false.

### S4 — stage 2: the re-derivation · **L** · *gated on C13*

1. Map legacy's `[sN]` donor-site pattern indices onto this project's `sites` perception, checked
   against `ligand_cases.tsv`. A label whose pattern has no counterpart is refused, not guessed.
2. Build one `BuildSpec` per scenario and run it through `construct` and the relax runner. This is
   hours of compute and is resumable, so put it in the jobs table behind its own flag, not inside a
   script function. Do not ship the pre-flight without the executor; `RELAXATION_IS_EXECUTED` in
   `energy/relax.py` records what that cost last time.
3. Score every candidate twice per §6.
4. Compare against the deduplicated archived set from C12.

**Done when:** every archived composition is either reached or refused with a reason;
within-multiset orderings match; the two scoring schemes are reported side by side.

### S5 — the solvation energies · **M**

S1 covers the arithmetic; this covers the single points. `XTBBackend.SOLVENTS` already accepts
`water`, `ethanol` and `dmf`, which are the archive's three.

1. Recompute the 117 ALPB single points (39 species × 3 solvents) of `solvation_raw.csv`.
2. Compare at C11's tolerance. The two decision metrics then follow from S1's tested arithmetic.

**Done when:** the raw single points reproduce and no new analysis code was written to check them.

---

## 8. Risks

| Risk | Signal | What to do |
|---|---|---|
| Stage 1 passes and means nothing | green at 1.0 eV against a 0.096 eV archive spread | C11: two criteria, calibrated, reported separately |
| Stage 2 expands into M8 | the `[sN]` mapping turns into perception work, or the enumeration needs unbuilt features | C13, taken in advance: close M7 on stages 0 + 1 |
| A tolerance becomes folklore | a threshold taken from a figure in this file | §3's numbers are points, not distributions. Calibrate from the archive's spread, as θ_geom was |
| The archive is treated as ground truth | a stage-2 mismatch is fixed by changing the new stack | it disagrees with itself by up to 1.287 eV (§3d) and reversed its own verdict once solvation was added. A divergence is a finding; which side is wrong is the question |
| No machine has tblite | stages 1, 2 and 5 cannot run | stage 0 needs no scientific stack. Stages 1/2/5 skip the way `tests/test_energy.py` already skips |
| Recording waits for a complete pass | stage 1 run, result in a terminal | S3 is sized S and ordered before S4. A stage-1-only result still closes something |

---

## 9. Bookkeeping when M7 lands

- Decision Ledger: D-numbers for C11, C12 and C13, one changelog line each. D17 is unchanged.
- `PLAN_implementation.md`: narrow §5's "trust the relative energies attached to any of it" to
  what was gated; replace §7 item 1 with C13's outcome; point §6's golden/regression bullet at
  `tests/test_energy_cases.py`; move §M7 to `archive/PLAN_completed.md` if C13 closes it.
- `CODE_ARCHITECTURE.md`: add `energy_cases.tsv` beside the other curated tables; point the
  "where to read about energies" row at the recorded report.
- `data/reference/NOTES.md`: provenance paragraph for `energy_cases.tsv`.
- `BUGS.md`: nothing. A wrong number in a report is an erratum, corrected in the report with the
  archived figures beside it.
- Delete this file.

---

## 10. Order

```
S0 ─► S1 ─┬─► S2 ─► S3 ─► [C13] ─► S4
          └─► S5
```

S1 runs on any machine. S2 and S5 need tblite and are independent of each other. S3 should record
a stage-1-only result rather than wait. S4 is the only slice larger than M.
