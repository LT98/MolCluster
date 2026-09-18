# M7 work plan — evaluate everything, but leave the decision points open

M7's primitives are built: backends behind one protocol, `MethodSpec` on every number, the
fidelity ladder, and a reference scheme that refuses unsound subtractions. Two capabilities are
missing: **energies over a whole route** rather than one reaction at a time, and a **place to
decide what is worth computing**.

**This plan builds the first and only the seam for the second.** M7 evaluates everything —
every candidate relaxed, every route enumerated, every step scored, nothing pruned and nothing
ranked. That is what the code does today; the work is to make it do so through named, recorded
decision points, so that a later policy plugs into a call site that already exists rather than
restructuring the pipeline.

Anything that decides *which* candidates matter is out of scope here. §4 says what is deferred
and to where.

Delete this file when M7 lands. `PLAN_implementation.md` §M7 says what the milestone is; this
file says what to do. **It resizes it from S to M.**

---

## 1. What is built, what is missing

| Piece | State |
|---|---|
| `energy/backends.py`, `energy/relax.py` | built, general |
| `energy/reference.py` — `reaction_balanced_energy(reg, reaction_id)` | built; one reaction, any reaction |
| walking the DAG to get a **route** | `registry.incoming_routes` returns one step; nothing composes steps |
| per-step energies along a route | does not exist |
| any decision point between building and computing | does not exist (§2) |
| `scripts/regress_m7.py` | built, Fe only, never recorded |

`pathways/__init__.py` is a docstring with no module behind it.

---

## 2. Where the decisions would go, and what is there now

`runner.queue_relax` (line 1153) emits one `relax` task per completed build. Its only filter is
`existing_relaxation`, a dedup against work already done at that fidelity. `estimate` (line 372)
states the consequence: `relaxations = builds`. Nothing else in the pipeline chooses.

The one decision of this kind that does exist is the precedent to copy: `runner.py:846` decides
whether a marginal-clash geometry is built *because a relaxation is coming* — a named threshold
(`qc.MARGINAL_OVERLAP`), a structured record carried with the geometry, and a re-check afterwards.

Signals already computed and discarded at scheduling time, which a later policy would read:
the QC report on every geometry, `activation_ease` (`descriptors/ease.py:245`, D18's zero-QM floor,
written to produce "a QM work list" and read only by the site catalog at `sites/state.py:181`),
join and placer strain verdicts, and L3 conformer clustering at θ_geom = 0.15 Å.

M7 does not read any of them. It puts the call site where they would be read.

---

## 3. The four seams

Each seam is a named policy resolved from a lookup table, exactly as `get_backend(name)` resolves
a backend from `_BACKENDS`. Each has a naive default that reproduces today's behaviour.

| # | Seam | Called by | Receives | Returns | Naive default | Plugs in later |
|---|---|---|---|---|---|---|
| 1 | `SelectionPolicy` | `runner.queue_relax:1153` | `reg`, `spec`, `task`, `Outcome` | compute at fidelity F, or skip with a reason | `"all"` — every build at `MODE_FIDELITY[run_mode]` | a screen reading the §2 signals; ladder promotion |
| 2 | `RouteFilter` | `pathways.route.routes_to` | the enumerated routes to a product | the subset to evaluate, plus a reason per route dropped | `"all"` — every route within `max_depth` | route pruning |
| 3 | `StepFilter` | `pathways.evaluate.route_energies` | the steps of one route | the steps to price | `"all"` — every step | skipping steps whose energy will not change an answer |
| 4 | `Ranking` | the caller of `route_energies` | evaluated routes | an order | `"none"` — DAG order, no ordering | M8's C6 barrier proxy, sink detection, A-vs-B |

**Seam 3 is not where D17 lives.** Refusing an unsound equation is a correctness rule and runs
whatever the policy says. A seam that could switch it off would be a way to get a number D17
exists to withhold.

**Seam 4 is defined and left empty on purpose.** It is M8's, and naming it now is what stops M8
from having to restructure `evaluate` to insert it.

### What makes a seam accessible rather than a comment

1. **Named, not callable.** A policy is selected by name and the name is recorded. Ground rule 5
   and replay: a spec records what produced a run, and a callable injected at runtime cannot be
   recorded, so a run using one could not be replayed from its record. A caller wanting custom
   behaviour registers it under a name, as an ML backend is registered.
2. **A verdict is always recorded, including the naive one.** The `"all"` policy emits a verdict
   per candidate — policy name, version, decision, reason — even though the decision is always
   *compute*. If the naive path records nothing, the first real policy has to invent the record,
   the schema and the reading surface at the same time, which is how a pruned candidate ends up
   indistinguishable from an absent one. This is D18's absent-is-absent rule one level up.
3. **Versioned.** `ALGO_VERSIONS` gains an entry per seam, so a stored verdict says which policy
   produced it (ground rule 6).
4. **Tested through the seam.** Each seam gets a test that injects a deliberately restrictive
   policy and asserts the pipeline honours it. Without one, a seam that is silently bypassed looks
   identical to a seam that works, until the day it matters.
5. **The default is provably today's behaviour.** A test asserts the naive path queues exactly
   what the current code queues.

### Where each name is declared

Seam 1 is a build-time choice, so it is a `BuildSpec` field with a migration, defaulting to the
name that reproduces earlier runs — the idiom `pathways` and `max_distinct_ligands` already
follow (`spec.py:_migrate`). Seams 2–4 are query-time: routes are evaluated when asked, not during
a build, so their names are arguments to the query functions and are recorded in whatever the
evaluation writes. Do not put a query-time policy in a build spec.

### What is not being built

No plugin discovery, no entry points, no config-file-loaded callables, no user-supplied Python.
The registry is a module-level dict, like `_BACKENDS`.

---

## 4. Deferred, and to where

| Deferred | To |
|---|---|
| any real selection policy — what it reads, how it scores | a later slice, once there is something to measure it against |
| the acceptance criterion for a policy (recall on a labelled set) | proposed with the policy, not before |
| barrier proxy, sink detection, route ranking, A-vs-B | **M8** — C6 and C7 |

M7 produces numbers. M8 produces verdicts. `activation_ease(partner=...)` already raises on C7 and
nothing in M7 may route around it.

Route evaluation's first test is **M6's exit gate**, not the archive: one paddlewheel node reached
by a sequential route and a nucleus-first route, both in the DAG.

---

## 5. The archived data, and the three defects in it

`legacy/` holds 43 Ni and 49 Fe candidates with energies and orderings, plus a complete solvation
sweep. M7 uses it as a fixed regression target, and it is the labelled set a future policy would
be measured against. All of the following was measured from `legacy/` with the Python standard
library; no xTB.

**(a) Duplicate candidates disagree.** A composition appearing in more than one scenario should
carry one energy.

| | rows | distinct compositions | appear more than once | energies differ | largest |
|---|---|---|---|---|---|
| Ni | 43 | 21 | 18 | 9 | 1.287 eV |
| Fe | 49 | 24 | 19 | 11 | 0.615 eV |

Cause: the relaxation step cap differed between scenario groups, which both reports state as a
caveat. `BTC[s1]_BTC[s1]_planar_q-4` ran 64 steps in scenarios A/B and 44 in C/D. C13 decides
which copy is the reference value.

**(b) Two reported figures do not follow from the data.**

| Where | Problem |
|---|---|
| `fe_btc_report.md`, "states below node" | reads 0 / 4 / 0 / 4. Counting strictly below −58.97 over `fe_btc_xtb_results.csv` gives 0 / 5 / 0 / 6. The same rule reproduces the Ni report's column exactly (0 / 3 / 0 / 4), and Fe's spans, candidate counts and deepest-gap figures (−3.836 → −3.8, −3.372 → −3.4) all reproduce. Rounding the target to −59.0 fixes scenario B and not D; using scenario D's own node copy (−58.75) gives 6. |
| `ni_btc_report.md`, "additive-independent (−36.4 eV in BTC-only; −36.1 eV in every additive scenario)" | the two numbers are the same composition at 64 and 44 relaxation steps. The 0.26 eV gap measures the step cap, not additive sensitivity. |

**(c) `scripts/regress_m7.py:125` says the Ni results CSV "is not in the repo".** It is, at
`legacy/ni_btc_xtb_results.csv`, with a populated `E_form_eV` for all 43 candidates; the same holds
for Fe. Absent are the geometries (`*_outputs/`) and the figures. The rest of that docstring is
accurate, including the empty `formation_eV` in `ni_btc_results.csv`.

**What holds:** all 92 archived `E_form_eV` re-derive from that row's `E_relax_eV` and the stored
`*_refs.json` to within 0.001 eV; all 14 rows of `solvation_summary.csv` re-derive from
`solvation_raw.csv` to within 0.005 eV. The two runs differ by up to 0.096 eV on identical free
ligands (EDTA in DMF; BTC 0.005 eV) — the archive's own same-species spread, and the only empirical
anchor for the word "noise".

---

## 6. What the reference scheme does to any comparison

The archived equation trips all three rules in `check_reference_quality`: a bare metal ion, two
naked polyanions, and a metal–donor bond count going from 0 to 6 across the arrow. It needs
`strict=False` and carries `isodesmic=False` (D17). Two consequences that apply to any evaluation,
not only the archive:

- Within one ligand multiset the reference terms are a constant, so the ordering is the same under
  any balanced scheme.
- Across different multisets it is not. The archived headline — the Ni–EDTA chelate 8.5 eV below
  the Ni–BTC node — is `[E(NiEDTA) − E(EDTA⁴⁻)] − [E(Ni(BTC)₂) − 2E(BTC³⁻)]`, carried by the
  free-anion terms D17 refuses. The reports' own `E_form_per_bond` column reverses that ordering
  (node −18.2, chelate −7.5 eV/bond), and the solvation run moved the same gap from −8.5 eV to
  −0.3…+0.5 eV.

So every stored number says which scheme produced it, and a cross-multiset comparison under a
non-isodesmic scheme is refused rather than returned with a caveat.

---

## 7. Decisions to make before writing code

Each needs a Decision Ledger entry with a D-number and a changelog line. Next free numbers are
after M6's (D23 onward, if M6 lands D20–D22).

**C11 — how a policy is selected.**
*Options:* a callable passed at the call site · a name resolved from a lookup table.
*Recommendation:* a name. A run must be replayable from its record (ground rule 5), and a callable
cannot be recorded. Reuses `get_backend`'s mechanism, including its error: unknown name, here are
the known ones.

**C12 — what a verdict record holds and where it lives.**
*Recommendation:* one row per candidate carrying policy name, policy version, decision, and reason,
written by every policy including the naive one. Not-computed and computed-and-poor must be
distinguishable in the registry.
*Cost:* a table and a migration for records that, under the naive policy, all say the same thing.
That cost is the point: it is paid once now, not at the same time as the first real policy.

**C13 — which copy of a duplicated archived candidate is the reference value (§5a).**
*Options:* lowest energy · the `converged=True` copy · the copy from the highest step cap · refuse
and re-run.
*Recommendation:* highest step cap. A capped relaxation is an upper bound, so the deepest copy is
the most relaxed, and the rule is stated rather than chosen per row. Record the choice as a column
in the curated table, not as logic in a script.

**C14 — what the calculator check means.**
The 1.0 eV tolerance in `regress_m7.py` is 10× the archive's largest same-species disagreement
(§5), so passing it shows the calculator is the same one, not that the numbers agree.
*Recommendation:* keep it, label it in the output as a calculator-identity check, and calibrate a
second, tighter per-species tolerance from the archive's own spread.

---

## 8. Slices

### S0 — correct the record, make the calls · **S**

1. Fix `scripts/regress_m7.py:125` per §5c.
2. Ledger entries for C11–C14, one changelog line each.
3. Correct §5b's two report figures, with the archived numbers beside them.

**Done when:** the docstring matches the filesystem and four D-numbers exist.

### S1 — the archived data, checked in the suite · **S**

Zero-QM, so it runs anywhere and becomes a permanent guard.

1. Add `data/reference/energy_cases.tsv`. Columns: `name`, `kind`, `source_file`, `subject`,
   `metric`, `expected`, `tol`, `note` — `subject` identifies the row being checked (candidate
   label, metal+solvent system, or report cell), `metric` names the quantity. No column holds an
   expression; the checker dispatches on `kind`. Header comment in the style of `node_cases.tsv`.
2. Add `scripts/check_energy_cases.py` (`load_cases`, `evaluate`) matching
   `scripts/check_cases.py`, and `tests/test_energy_cases.py` matching
   `tests/test_reference_cases.py`.
3. Rows: `formation` (92), `solvation` (14 × 5 columns), `report` (the per-scenario tables).
4. Apply C13's choice as a column, so the 21 Ni and 24 Fe reference values are readable.

**Done when:** `python -m pytest tests/test_energy_cases.py` passes on a machine with no
scientific stack installed.

### S2 — the policy mechanism · **S**

Build the lookup before the first seam uses it, so all four share one.

1. `pathways/policy.py`: a `Policy` protocol per seam kind, a module-level name→class dict, and
   `policy_for(kind, name)` raising on an unknown name with the known ones listed — `get_backend`'s
   shape.
2. A `Verdict` record per C12, and the table and migration to store it.
3. `ALGO_VERSIONS` entries, one per seam.

**Done when:** `policy_for` resolves the four naive policies and refuses a fifth name by listing
what exists.

### S3 — routes over the DAG, with seam 2 · **M**

1. `pathways/route.py`: `routes_to(reg, structure_id, *, max_depth, route_filter="all")` — compose
   `reactions` edges into ordered routes from starting materials to a product. A cycle or a missing
   reagent is refused by name, not dropped.
2. Two routes to one product are distinguishable by step sequence and choice-vector digests, which
   `incoming_routes` already projects for this reason.
3. Seam 2 with the `"all"` policy; a test injecting a restrictive one proves it is wired.

**Done when:** M6's paddlewheel product returns exactly two routes — sequential and nucleus-first —
and a restrictive filter injected in a test returns one, with a reason for the other.

### S4 — energies along a route, with seams 3 and 4 · **M**

1. `pathways/evaluate.py`: `route_energies(reg, route, *, fidelity, step_filter="all",
   ranking="none")` → per-step `ReactionEnergy` via the existing `reaction_balanced_energy`, each
   carrying its `MethodSpec` and `isodesmic` flag.
2. A step whose energy is not stored at that fidelity is reported **absent**, never zero and never
   silently dropped (D18). A step whose equation `check_reference_quality` refuses is reported
   refused, with the issue codes — independent of any policy (§3).
3. Seam 4 defined with the `"none"` ranking. No scoring, no max barrier, no sink detection.

**Done when:** M6's two paddlewheel routes return per-step energies or a stated reason per step; a
route with one unstored step returns a profile with a hole in it rather than a number; and an
injected step filter is honoured.

### S5 — seam 1 in the queue · **S**

1. `BuildSpec.selection: str = "all"`, with a migration defaulting older specs to `"all"` so an old
   spec re-run queues exactly what it queued before.
2. `runner.queue_relax` (line 1153) consults `policy_for("selection", spec.selection)` instead of
   queueing unconditionally. `existing_relaxation`'s dedup stays where it is — it is a correctness
   check, not a policy.
3. Every decision written as a verdict per C12.
4. `estimate` (line 372) reports the naive ceiling as it does now, and reports it as coming from
   the named policy.

**Done when:** a run under the default spec queues the same relax tasks as before the slice, each
with a verdict row; and a test injecting a restrictive policy queues fewer.

### S6 — the calculator check, both metals, recorded · **S**

1. `regress_m7.py` hardcodes `METAL, METAL_CHARGE = "Fe", 3` (line 53), `ARCHIVED` (line 49) and
   `archived["M"]` (line 80); `ni_btc_refs.json` keys it `Ni`. Add `--metal ni|fe`.
2. Carry the spin convention per archived run: Fe used a sextet (`high_spin_multiplicity`), Ni used
   minimal spin (`minimal_multiplicity`, whose docstring says it exists for this). Print which.
3. Apply C14: two tolerances, reported separately.
4. Record the result — `docs/reports/m7_regression_report.md` plus the committed `--json`: which
   stages ran, on what machine, on what date, the `MethodSpec`, the per-species deltas, and what is
   **not** covered.

**Done when:** `PLAN_implementation.md`'s "no recorded result" sentence is false.

### S7 — the archived rankings, re-derived · **L** · *defer to M8 unless something needs it sooner*

Only the geometries are missing (§5c), so this is a re-derivation with a per-candidate target.
Needs: legacy's `[sN]` donor-site patterns mapped onto this project's `sites` perception (checked
against `ligand_cases.tsv`, refused where there is no counterpart); one `BuildSpec` per scenario
through `construct` and the relax runner, in the jobs table behind its own flag because it is hours
of resumable compute; scoring under both schemes per §6.

Do not ship a pre-flight without the executor — `RELAXATION_IS_EXECUTED` in `energy/relax.py`
records what that cost last time.

---

## 9. Risks

| Risk | Signal | What to do |
|---|---|---|
| A seam becomes a plugin framework | discovery, entry points, or a callable loaded from config | §3's "what is not being built". A module dict and one lookup function |
| A seam is wired but bypassed | the naive default passes and nothing else was tried | every seam ships with a test injecting a restrictive policy |
| The naive policy changes behaviour | a run queues different tasks after S5 than before | S5's done-when is byte-equality against the current path, and the spec migration defaults old specs to `"all"` |
| M7 grows a policy anyway | a score appears in `policy.py` that ranks candidates | §4. The seam is the deliverable; the policy is not |
| M7 swallows M8 | a barrier proxy or a route ranking appears in `evaluate.py` | seam 4 exists so that M8 has somewhere to land. It stays `"none"` |
| A verdict table nobody reads | S5 lands and nothing surfaces the rows | acceptable for now and stated as such: the format is the deliverable, the reader arrives with the first real policy |
| The archive is treated as ground truth | a mismatch is fixed by changing the new stack | it disagrees with itself by up to 1.287 eV (§5a) and reversed its own verdict once solvation was added |
| No machine has tblite | S6, S7 and any QM cannot run | S0–S5 need no scientific stack beyond stored energies. QM slices skip as `tests/test_energy.py` already skips |

---

## 10. Bookkeeping when M7 lands

- Decision Ledger: D-numbers for C11–C14, one changelog line each. D17 unchanged.
- `ALGO_VERSIONS`: one entry per seam.
- `PLAN_implementation.md`: resize M7 from **S remaining** to **M**; move the DAG walk out of M8's
  module list into M7's; keep C6 and C7 with M8; point §6's golden/regression bullet at
  `tests/test_energy_cases.py`.
- `CODE_ARCHITECTURE.md`: rows for `pathways/route.py`, `pathways/evaluate.py`, `pathways/policy.py`
  and `energy_cases.tsv`; a §4 line for the seam rule if it becomes an invariant.
- `data/reference/NOTES.md`: provenance paragraph for `energy_cases.tsv`, stating that its targets
  are archived results rather than literature and that some of them disagreed with each other.
- `BUGS.md`: nothing. A wrong number in a report is an erratum, corrected in the report.
- Delete this file.

---

## 11. Order

```
S0 ─► S1
 └──► S2 ─┬─► S3 ─► S4
          └─► S5
S6 independent · S7 deferrable
```

S1 and S2 are independent of each other. S2 comes before any seam so all four share one mechanism.
S3 and S4 are the general capability. S5 is small because it is naive. S6 needs tblite and nothing
else. S7 is the only slice that can be deferred wholesale.
