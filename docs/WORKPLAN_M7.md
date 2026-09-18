# M7 work plan — energy over any route, and what to spend it on

M7's primitives are built: backends behind one protocol, `MethodSpec` on every number, the
fidelity ladder, and a reference scheme that refuses unsound subtractions. What is missing is not
a re-run of two archived datasets. It is two capabilities:

1. **Evaluate energies over any route in the reaction DAG**, not a curated candidate list.
2. **Decide which candidates get quantum mechanics**, so that (1) is affordable.

Today the code does neither: `reaction_balanced_energy` scores one reaction at a time and nothing
walks a route, and `queue_relax` relaxes everything that was built. The archived Ni/Fe BTC work
stops being the milestone and becomes its test set — the only data in this repo with known
answers to measure a selection policy against.

Delete this file when M7 lands. `PLAN_implementation.md` §M7 says what the milestone is; this
file says what to do. **It also resizes it: M7 is not S.**

---

## 1. What is built, what is missing

| Piece | General over any structure/route? | State |
|---|---|---|
| `energy/backends.py` | yes | built |
| `energy/relax.py` | yes | built |
| `energy/reference.py` — `reaction_balanced_energy(reg, reaction_id)` | one reaction, any reaction | built |
| walking the DAG to get a **route** | — | `registry.incoming_routes` returns one step; nothing composes steps |
| per-step energies along a route | — | does not exist |
| deciding which candidates get QM | — | does not exist (§2) |
| `scripts/regress_m7.py` | no — hardcoded to Fe/BTC | built, Fe only, never recorded |

`pathways/__init__.py` is a docstring with no module behind it.

---

## 2. Why it is narrow today

**Everything built gets relaxed.** `runner.queue_relax` (line 1153) emits one `relax` task per
completed build. Its only filter is `existing_relaxation` — a dedup against work already done at
that fidelity. Nothing asks whether the candidate is worth the cost. `estimate` (line 372) says so
outright: `relaxations = builds`.

**The one existing selection policy is the precedent to follow.** `runner.py:846` decides whether
a marginal-clash geometry is built *because a relaxation is coming*: a named threshold
(`qc.MARGINAL_OVERLAP`), a structured record carried with the geometry, and a re-check afterwards.
That is the right shape. It is currently the only decision of its kind in the codebase.

**The cheap screening score exists and nothing consumes it.** `descriptors/ease.py:245`
`activation_ease` is D18's zero-QM floor. It is written to produce "a QM work list rather than a
verdict" — its `provisional` flag marks sites where the table value is the wrong question. Only
`sites/state.py:181` calls it, to fill the site catalog. No scheduling decision reads it.

**Other zero-cost signals already produced and discarded at scheduling time:** the QC report on
every geometry, the strain/refusal verdicts from `join` and the placer, and L3's conformer
clustering at θ_geom = 0.15 Å (near-identical conformers do not each need their own xTB).

---

## 3. Boundary with M8

M8 is "reactions + pathways" and owns the claims. M7 owns the substrate. Keeping the line explicit
is what stops M7 from absorbing C6 and C7.

| M7 lands | M8 lands |
|---|---|
| `pathways/route.py` — enumerate routes over the reactions DAG | `pathways/proxy.py` — the barrier proxy (**C6**) |
| `pathways/evaluate.py` — per-step energies for a given route, with provenance | `pathways/score.py` — max barrier, cumulative ΔG, sink detection |
| the selection policy: which candidates get which fidelity | `pathways/search.py` and the A-vs-B driver |
| "here are the numbers for this route" | "this route is better than that one, and why" |

The two scoring questions are different and must not merge. M7's policy scores **is this candidate
worth computing** — a cost decision from descriptors that already exist. M8's proxy scores **how
hard is this step** — a chemistry claim. `activation_ease(partner=...)` already raises on C7, and
nothing in M7 may route around that.

**Route evaluation's first real test is M6's exit gate**, not the archive: M6 gate 1 produces one
paddlewheel node reached by a sequential route and a nucleus-first route. Two routes, one product,
both in the DAG. M7 puts per-step numbers on them; M8 says which is better.

---

## 4. The archive as a labelled test set

`legacy/` holds 43 Ni and 49 Fe candidates with energies, orderings and scenario groupings, plus a
complete solvation sweep. It is the only labelled data in the repo, which makes it the only place
a selection policy's **recall** can be measured: does the cheap screen keep the candidates that
turn out to matter?

Three defects have to be fixed before it can serve that purpose. All were measured from `legacy/`
with the Python standard library; no xTB.

**(a) The set contains duplicates that disagree.** A composition appearing in more than one
scenario should carry one energy.

| | rows | distinct compositions | appear more than once | energies differ | largest |
|---|---|---|---|---|---|
| Ni | 43 | 21 | 18 | 9 | 1.287 eV |
| Fe | 49 | 24 | 19 | 11 | 0.615 eV |

Cause: the relaxation step cap differed between scenario groups, which both reports state as a
caveat. `BTC[s1]_BTC[s1]_planar_q-4` ran 64 steps in scenarios A/B and 44 in C/D. A labelled set
whose labels disagree cannot measure recall. C12 decides which copy is the label.

**(b) Two reported figures do not follow from the data.**

| Where | Problem |
|---|---|
| `fe_btc_report.md`, "states below node" | reads 0 / 4 / 0 / 4. Counting strictly below −58.97 over `fe_btc_xtb_results.csv` gives 0 / 5 / 0 / 6. The same rule reproduces the Ni report's column exactly (0 / 3 / 0 / 4), and Fe's spans, candidate counts and deepest-gap figures (−3.836 → −3.8, −3.372 → −3.4) all reproduce. Rounding the target to −59.0 fixes scenario B and not D; using scenario D's own node copy (−58.75) gives 6. |
| `ni_btc_report.md`, "additive-independent (−36.4 eV in BTC-only; −36.1 eV in every additive scenario)" | the two numbers are the same composition at 64 and 44 relaxation steps. The 0.26 eV gap measures the step cap, not additive sensitivity. |

**(c) `scripts/regress_m7.py:125` says the Ni results CSV "is not in the repo".** It is, at
`legacy/ni_btc_xtb_results.csv`, with a populated `E_form_eV` for all 43 candidates; the same holds
for Fe. Absent are the geometries (`*_outputs/`) and the figures. The rest of that docstring is
accurate, including the empty `formation_eV` in `ni_btc_results.csv`.

**What does hold, and is worth keeping:** all 92 archived `E_form_eV` re-derive from that row's
`E_relax_eV` and the stored `*_refs.json` to within 0.001 eV; all 14 rows of
`solvation_summary.csv` re-derive from `solvation_raw.csv` to within 0.005 eV. The two runs differ
by up to 0.096 eV on identical free ligands (EDTA in DMF; BTC 0.005 eV) — the archive's own
same-species spread, and the only empirical anchor for the word "noise".

---

## 5. What the reference scheme does to any cross-candidate comparison

The archived equation trips all three rules in `check_reference_quality`: a bare metal ion, two
naked polyanions, and a metal–donor bond count going from 0 to 6 across the arrow. It needs
`strict=False` and carries `isodesmic=False`. That is D17, already decided. Two consequences that
apply to **any** route evaluation, not just the archive:

- Within one ligand multiset the reference terms are a constant, so the ordering is the same under
  any balanced scheme.
- Across different multisets it is not. The archived headline — the Ni–EDTA chelate 8.5 eV below
  the Ni–BTC node — is `[E(NiEDTA) − E(EDTA⁴⁻)] − [E(Ni(BTC)₂) − 2E(BTC³⁻)]`, carried by the
  free-anion terms D17 refuses. The reports' own `E_form_per_bond` column reverses that ordering
  (node −18.2, chelate −7.5 eV/bond), and the solvation run moved the same gap from −8.5 eV to
  −0.3…+0.5 eV.

So `evaluate` reports which scheme produced each number and refuses to compare across multisets
under a non-isodesmic one. Where two schemes disagree, report the disagreement. Do not average
them and do not pick one.

---

## 6. Decisions to make before writing code

Each needs a Decision Ledger entry with a D-number and a changelog line. Next free numbers are
after M6's (D23 onward, if M6 lands D20–D22).

**C11 — what the selection policy reads, and whether it is declared or inferred.**
Available at zero cost: QC status, `activation_ease` (including `provisional`), join/placer strain
verdicts, L3 conformer clustering, identity dedup.
*Recommendation:* declared in `BuildSpec`, never inferred. Ground rule 5 forbids a silent default,
and a policy that changes what gets computed is exactly the kind of inference that must branch or
be stated. Add a `selection` field with an explicit "everything" value so today's behaviour stays
expressible and stays the default until the policy is calibrated.
*Cost:* one more spec field and a migration.

**C12 — how a candidate that was not computed is recorded.**
*Recommendation:* as pruned, with the reason and the cheap score that pruned it. Never silently
absent and never scored zero. This is D18's absent-is-absent rule one level up: "we did not compute
this" and "this scored badly" must be distinguishable in the registry, or the first rendering of a
landscape becomes a lie.
*Forces:* a `pruned` record with `MethodSpec`-equivalent provenance for the screen itself.

**C13 — which copy of a duplicated archived candidate is the label (§4a).**
*Options:* lowest energy · the `converged=True` copy · the copy from the highest step cap · refuse
and re-run.
*Recommendation:* highest step cap. A capped relaxation is an upper bound, so the deepest copy is
the most relaxed, and the rule is stated rather than chosen per row. Record the choice as a column
in the curated table, not as logic in a script.

**C14 — the acceptance criterion for the policy.**
*Recommendation:* recall on the labelled set. The policy is acceptable when, over the archived
candidates, it retains every candidate in the top-N of the known ordering at a stated N and a
stated budget fraction. State both numbers; a policy with no measured recall is a guess about
which structures matter.
*Depends on:* C13, because recall against labels that disagree is not measurable.

**C15 — what "reproduce" means for the calculator check.**
The 1.0 eV tolerance in `regress_m7.py` is 10× the archive's largest same-species disagreement
(§4), so passing it shows the calculator is the same one, not that the numbers agree.
*Recommendation:* keep it, label it in the output as a calculator-identity check, and calibrate a
second, tighter per-species tolerance from the archive's own spread.

---

## 7. Slices

### S0 — correct the record, make the calls · **S**

1. Fix `scripts/regress_m7.py:125` per §4c.
2. Ledger entries for C11–C15, one changelog line each.
3. Correct §4b's two report figures, with the archived numbers beside them.

**Done when:** the docstring matches the filesystem and five D-numbers exist.

### S1 — the labelled set, checked in the suite · **S**

Zero-QM, so it runs anywhere and becomes a permanent guard.

1. Add `data/reference/energy_cases.tsv`. Columns: `name`, `kind`, `source_file`, `subject`,
   `metric`, `expected`, `tol`, `note` — `subject` identifies the row being checked (candidate
   label, metal+solvent system, or report cell), `metric` names the quantity. No column holds an
   expression; the checker dispatches on `kind`. Header comment in the style of `node_cases.tsv`.
2. Add `scripts/check_energy_cases.py` (`load_cases`, `evaluate`) matching
   `scripts/check_cases.py`, and `tests/test_energy_cases.py` matching
   `tests/test_reference_cases.py`.
3. Rows: `formation` (92, §4), `solvation` (14 × 5 columns), `report` (the per-scenario tables).
4. Apply C13's dedup as a column, so the 21 Ni and 24 Fe labelled candidates are readable.

**Done when:** `python -m pytest tests/test_energy_cases.py` passes on a machine with no
scientific stack, and the labelled set is 21 + 24 rows with one energy each.

### S2 — routes over the DAG · **M**

1. `pathways/route.py`: `routes_to(reg, structure_id, *, max_depth)` — compose `reactions` edges
   into ordered routes from starting materials to a product. A cycle or a missing reagent is
   refused by name, not dropped.
2. Route identity: two routes to one product are distinguishable by their step sequence and
   choice-vector digests, which `incoming_routes` already projects for this reason.

**Done when:** M6's paddlewheel product returns exactly two routes — sequential and
nucleus-first — and `examples`' mononuclear two-route structure returns two.

### S3 — energies along a route · **M**

1. `pathways/evaluate.py`: `route_energies(reg, route, *, fidelity)` → per-step `ReactionEnergy`
   via the existing `reaction_balanced_energy`, each carrying its `MethodSpec` and `isodesmic`
   flag.
2. A step whose energy is not stored at that fidelity is reported **absent**, never zero and never
   silently dropped (D18).
3. A step whose equation `check_reference_quality` refuses is reported refused, with the issue
   codes.
4. No scoring. No max barrier, no sink detection, no ranking — those are M8 (§3).

**Done when:** M6's two paddlewheel routes return per-step energies or a stated reason per step,
and a route with one unstored step returns a profile with a hole in it rather than a number.

### S4 — the selection policy · **L**

The slice that makes the rest affordable.

1. `BuildSpec.selection` per C11, defaulting to today's "everything" so no existing spec changes
   behaviour. Migration per ground rule 7.
2. A screen between build and relax, reading the §2 signals. `runner.queue_relax` (line 1153)
   consults it instead of queueing unconditionally.
3. Promote up the ladder rather than running one rung on everything: cheap fidelity across the
   set, expensive fidelity on survivors. `relaxed_from` already chains the rungs.
4. Pruned candidates recorded per C12.
5. `estimate` (line 372) reports a budget instead of `relaxations = builds`, so the page shows
   what the run will actually cost.

**Done when:** a spec with a policy computes strictly fewer relaxations than the same spec without
one; every skipped candidate has a record naming the reason and the score; and the default spec's
behaviour is byte-identical to today's.

### S5 — measure the policy · **M**

1. Replay the archived candidate set through the screen alone, no QM.
2. Report recall against C13's labels at the C14 criterion: which candidates the screen would have
   dropped, and where they sat in the known ordering.

**Done when:** the recall number exists and is recorded. A policy shipped without it is C14 unmet.

### S6 — the calculator check, both metals, recorded · **S**

1. `regress_m7.py` hardcodes `METAL, METAL_CHARGE = "Fe", 3` (line 53), `ARCHIVED` (line 49) and
   `archived["M"]` (line 80); `ni_btc_refs.json` keys it `Ni`. Add `--metal ni|fe`.
2. Carry the spin convention per archived run: Fe used a sextet (`high_spin_multiplicity`), Ni
   used minimal spin (`minimal_multiplicity`, whose docstring says it exists for this). Print
   which.
3. Apply C15: two tolerances, reported separately.
4. Record the result — `docs/reports/m7_regression_report.md` plus the committed `--json`: which
   stages ran, on what machine, on what date, the `MethodSpec`, the per-species deltas, and what
   is **not** covered.

**Done when:** `PLAN_implementation.md`'s "no recorded result" sentence is false.

### S7 — the archived rankings, re-derived · **L** · *defer to M8 unless S5 needs it*

Only the geometries are missing (§4c), so this is a re-derivation with a per-candidate target.
Needs: legacy's `[sN]` donor-site patterns mapped onto this project's `sites` perception (checked
against `ligand_cases.tsv`, refused where there is no counterpart); one `BuildSpec` per scenario
through `construct` and the relax runner, in the jobs table behind its own flag because it is
hours of resumable compute; scoring under both schemes per §5.

Do not ship a pre-flight without the executor — `RELAXATION_IS_EXECUTED` in `energy/relax.py`
records what that cost last time.

---

## 8. Risks

| Risk | Signal | What to do |
|---|---|---|
| M7 swallows M8 | a barrier proxy or a route ranking appears in `evaluate.py` | §3's table. M7 produces numbers; M8 produces verdicts |
| The policy hides real candidates | recall below C14's criterion on the labelled set | S5 exists to measure this before the policy ships. Default stays "everything" until it passes |
| A pruned candidate looks like a bad one | a landscape plot with holes in it that read as data | C12: pruned is a record, not an absence |
| A threshold becomes folklore | a screen cutoff taken from a figure in this file | §4's numbers are points, not distributions. Calibrate from the archive's spread, as θ_geom was |
| The archive is treated as ground truth | a mismatch is fixed by changing the new stack | it disagrees with itself by up to 1.287 eV (§4a) and reversed its own verdict once solvation was added. It is a labelled test set, not an oracle |
| No machine has tblite | S6, S7 and any QM cannot run | S1–S5 need no scientific stack except S3's stored energies. QM slices skip the way `tests/test_energy.py` already skips |
| Sizing | M7 planned as S and carrying S4 | it is not S. S4 is L and S2/S3/S5 are M. Say so in the plan |

---

## 9. Bookkeeping when M7 lands

- Decision Ledger: D-numbers for C11–C15, one changelog line each. D17 unchanged.
- `PLAN_implementation.md`: resize M7 from **S remaining**; move `pathways/reaction.py`'s DAG walk
  out of M8's module list into M7's; narrow §5's "trust the relative energies attached to any of
  it" to what was gated; point §6's golden/regression bullet at `tests/test_energy_cases.py`.
- `CODE_ARCHITECTURE.md`: rows for `pathways/route.py`, `pathways/evaluate.py`, the selection
  policy, and `energy_cases.tsv`.
- `data/reference/NOTES.md`: provenance paragraph for `energy_cases.tsv`, stating that its targets
  are archived results rather than literature and that some of them disagreed with each other.
- `BUGS.md`: nothing. A wrong number in a report is an erratum, corrected in the report.
- Delete this file.

---

## 10. Order

```
S0 ─► S1 ─┬─► S2 ─► S3 ─────────────► S6
          └─► S4 ─► S5 ─► [C14] ─► S7 (or M8)
```

S1 runs anywhere and fixes the labels everything else is measured against. S2 and S3 are the
general capability; S4 and S5 are what make it affordable and prove it is safe. S6 is small and
independent. S7 is the only slice that can be deferred wholesale.
