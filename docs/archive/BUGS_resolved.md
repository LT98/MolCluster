# Resolved issues

**Archive.** `../BUGS.md` is the active tracker and holds only what is still open. An entry
moves here when it is fixed and covered by a test.

Each entry keeps its diagnosis as well as its fix. The diagnosis is the expensive part and it
explains the shape of the fix; several of these also record what building the fix turned up
that the diagnosis had not, which is the part worth re-reading before touching the same code.

Per-milestone engineering history lives in [`PLAN_completed.md`](PLAN_completed.md) §4; this
file is specifically *things that behaved wrongly*.

---

## B25 — a step was claimed while the rung it joins onto was still being built ✅

**`correctness` · `registry/jobs.claim_task` · found in the salt-study smoke run.**

The queue orders by priority and nothing else. With eight workers, a `grow` step could be
claimed while another worker still held its `parent_task`, and `_rung` then rejected it as
`pathway_parent_missing` with the parent "(claimed)" — a hole in the ladder that the chemistry
did not cause (13 steps in the NiCl₂ smoke run).

**Fix:** `claim_task` skips a task whose `parent_task` or `ligand_task` is pending or
claimed. An empty claim then no longer means an empty queue, so `work()` waits while
`waiting_on_live_work` says the blocking work belongs to a live worker, and exits as before
once every claim is held by a dead one. Re-measured: every remaining `pathway_parent_missing`
has a parent that was genuinely rejected. Test: `tests/test_jobs.py`.

---

## B24 — a grow tried one orientation and rejected what `place` built easily ✅

**`correctness` · `runner._execute_grow` · found in the salt-study smoke run.**

`grow` called `join` once — first open vertex, lone pair 0, torsion well 0 — and rejected the
step on any clash. `place` builds the same compositions from a whole-sphere layout and passed
them, so the ladder's recorded edges went missing exactly where a crowded sphere needed a
different pose. Measured on the NiCl₂/tHQ fallback spec (construct): **336 of 336** tHQ steps
onto an octahedral rung and 178 of 333 onto a tetrahedral one rejected, O···O down to 0.73 Å —
an sp² O has few wells, and each swung the ring's other oxygens into a cis water.

**Fix:** `_first_clear_join` tries every open vertex, lobe and well, first at the well's own
angle and then at roll offsets of ±30° steps (`ROLL_OFFSETS_DEG`), and takes the first product
that passes clash QC, else the least bad. `join(roll_deg=)` records the offset in the choice
vector only when non-zero, so no existing choice vector changes. Same spec: tetrahedral 336/336
built, octahedral 168 clean + 63 marginal (built in `ml_go`) of 336; the 105 left are the most
crowded rungs. Overall rejections 688 → 322 of 1778. Tests: `tests/test_join.py` (roll recorded
only when used), the pathway suites unchanged.

---

## B23 — a "deprotonation" edge could also change which atom binds the metal ✅

**`correctness` · `energy/protons.py` · found reviewing the MVP energies for the interim
presentation.**

`deprotonation_pairs` keyed structures on (heavy-atom formula, charge − H count), so it paired
every protonated row with every row of the same formula one proton lighter. On
`mvp_ni_thq_cl.db` that was 458 edges, a median of 8 and up to 21 partners per acid, and one
edge could also move which oxygen bound Ni or which geometry frame was compared. Same-charge
tHQ deprotonations spanned −3.65 to +5.24 eV: the number was a deprotonation plus an
isomerisation.

**Fix:** a pair is exact. `TypedGraph.without_proton(h)` removes one labile H (bonded to one
non-carbon, non-metal atom) and lowers the charge; the result is looked up by L1 hash, charge
and multiplicity. D15 records a delocalised charge as graph-level only, so both conventions
(charge on the atom, charge delocalised) are tried; everything else must match exactly. Two
protons apart is two edges through the intermediate. On the same registry: 87 pairs, at most
one partner per distinct proton, and Ni²⁺ tHQ deprotonations now span −0.61 to +1.92 eV.
The Ni–ClH species remain the outliers, which is B22's input rule to prevent, not this rule's.
Tests: `tests/test_reference.py` (a formula look-alike is refused; equivalent protons give one
edge).

---

## B22 — an anion entered as one was registered as neutral ✅

**`correctness` · `sites/protomers.py` · found planning the NiCl₂ vs Ni(OAc)₂ salt study.**

`enumerate_protomers` gave every protomer `charge = −k` (k = protons removed) and built the
parent at charge 0, so the input's own formal charge was discarded. `[Cl-]` and `CC(=O)[O-]`
both came out **q0**, and every complex built from them carried the error into its net charge
(`runner._build_sphere` sums component charges). Nothing refused it.

It had stayed hidden because the only chloride dataset (`mvp_ni_thq_cl.db`) entered chloride as
`Cl` (HCl) with one deprotonation — which charges correctly, but also creates Ni–ClH species
that do not exist in water, and those produced the most negative "deprotonation" energies in
that registry (−3.65 eV).

**Fix:** the charge is read off the molecule — `Chem.GetFormalCharge` of the parent and of each
deprotonated variant (`deprotonate` already sets the formal charges it creates). A neutral input
is unchanged, so no stored identity moves; an anion's zero-deprotonation state is labelled
`as_given` rather than `neutral`. Tests: `tests/test_protomers.py` (anion charges, acetic acid
0/−1, and Ni(Cl)(H₂O)₃ +1 / NiCl₂(H₂O)₂ 0 from `_build_sphere`).

---

## B18 — a `grow` step wrote an edge without the ligand it added ✅

**`correctness` · `runner.py` · found while pricing provenance edges in the viewer, fixed in
the same change.**

`grow` recorded two reagents — the rung below and the free ligand — and the ladder is built so
that those balance: `_parent_payloads` keeps `n_co` fixed, so the co-ligand waters cancel
across the arrow and `parent + ligand = product` in atoms and charge.

The ligand lookup did not hold up its end:

```python
ligand_sid = None
if payload.get("ligand_task") is not None:
    row = reg.conn.execute("SELECT structure_id FROM tasks WHERE id=?", …).fetchone()
    ligand_sid = None if row is None else row["structure_id"]
```

If the ligand task had no structure yet — not run, failed, or rejected — `ligand_sid` stayed
`None` and the list comprehension that builds `reagent_ids` filtered it straight back out. The
join still happened and the edge was still written, **with only the parent as reagent**. The
result is an edge unbalanced by exactly one ligand formula, with no metal in the delta — a
signature that distinguishes it from a `place` edge, which is unbalanced by its whole
composition.

The asymmetry was the tell: the parent rung goes through `_rung`, which raises `_Rejected`
with a reason naming the missing task. The ligand went through a bare `SELECT`. It now raises
the same way, with code `pathway_ligand_missing`.

**What the diagnosis turned up:** nothing counted these. A silently short edge is
indistinguishable from a correct one until something tries to *price* it, which nothing did
until the provenance panel. That is the argument for pricing being a first-class reader rather
than a report — it exercises invariants that only writers had been trusted with.

---

## B13 — a bare hydroxide perceived zero donors ✅

**`correctness` · `sites/perception.py` · found while planning M6, fixed in S2.**

| molecule | donors perceived, before |
|---|---|
| `[OH-]` | **none** |
| `C[O-]` methoxide | `alkoxide_O` |
| `CO` methanol | `alkoxide_O` |

Methoxide and methanol perceived normally, so the rule was dropping `[OH-]` specifically. A
hydroxide bridge is among the commonest motifs in polynuclear chemistry and the M6 battery
needs one, so this was not a curiosity.

**Cause: one assumption in `_classify_anionic`, true of every donor it recovers except this
one.** That function exists because "activate = deprotonate", so a structure recalled from the
registry is usually already activated — anionic and unprotonated. It opened with
`if atom.GetFormalCharge() >= 0 or _n_hydrogens(atom) > 0: return None`. Hydroxide is anionic
**and** still carries its proton, so the second clause dropped it. The dead giveaway was that
the function already had a `hydroxide_O` branch for an oxygen with no heavy neighbour — it was
simply unreachable.

**Fixed** by splitting the guard: an anion that is still protonated is a hydroxide when it has
no heavy neighbour, and otherwise still gets no donor rather than a guessed one. `[OH-]` now
perceives `hydroxide_O`, which already declared `mono`/`mu2`/`mu3` in `sites.frames`.

**What the fix turned up, and it is a separate defect: `[O-2]` types as `hydroxide_O` too.**
The `not heavy` branch does not distinguish a bare oxide from a hydroxide, so an oxo reports a
hydroxide's type and pKa. It is filed as [B15](../BUGS.md#b15) rather than folded in here,
because giving oxo its own type means a new curated row in `donor_descriptors.tsv` and that is
a change to the reference dataset, not to this function. It does not block S2: µ3/µ4 oxo
bridges are handled as **centres**, which never asks perception for a donor type.

**The neighbouring limitation recorded in the original entry is now also closed**, by a
different route: a donor with two neighbours has one determined axis, so a bridging aqua
implied an M···M of **0.000 Å**. Treating the bridging atom as a centre (`WORKPLAN_M6` §4,
`geometry.placer.bridging_metal_positions`) gives it 3.084 Å at 104.5°.

---

## B12 — `metal_block` gave a metal a formal charge nothing else does ✅

**`correctness` · `assembly/construct.py` · found while planning M6, fixed in S0.**

`graph.from_mol.from_rdkit` sets a metal's `formal_charge` to **0** — D15 puts charge on the
graph, not on the atom — and `examples.py` writes its metals the same way.
`assembly.construct.metal_block` wrote `formal_charge=oxidation_state`.

That string is not decoration: `NodeLabel.key()` emits `element/formal_charge/oxidation_state/
spin_class` for a metal, and `graph.canon._cert_entries` hashes it. Measured 2026-09-16 on a
two-Cu fragment identical but for that field:

| metal label | L1 |
|---|---|
| `Cu/+0/+2/hs` | `ef2229b65dcd3f97…` |
| `Cu/+2/+2/hs` | `b63727cfed06a87a…` |

So the same species reached through the runner (placer → `to_rdkit` → `from_rdkit`) and through
the enumerator (`metal_block`) landed on two different nodes. **That is D2's central claim — one
node, two routes — broken at the producer**, and it had not bitten only because nothing had yet
built the same species both ways.

**Fixed** by making `metal_block` write 0. The oxidation state is already carried in its own
label field, so nothing is lost and the charge stops being counted twice (atom, *and*
`TypedGraph.charge`). Gate:
`tests/test_construct.py::test_a_metal_is_labelled_the_same_way_by_every_producer` pins the
label string itself — `Cu/+0/+2/hs` — against `from_rdkit`'s, so the two producers cannot
drift apart again silently.

**What the fix turned up: the whole suite was blind to it.** 683 tests passed before and
after, and no golden hash moved, because every stored hash descends from `examples.py` and
nothing asserted on a `metal_block` identity. A correctness bug that no test can see is a
correctness bug that will be reintroduced, which is why the fix ships with a producer-agreement
test rather than with a hash update.

**It moves identities, and the cost is the one already accepted for [B2](../BUGS.md#b2).** Rows
written under the old behaviour keep their stored hash and are never re-derived (ground rule 6
/ D19), so the same species can occupy two rows across the change. In practice the cost was not
paid: `metal_block` has no caller in `src/` or `scripts/` — the enumerator is reached from tests
only — so no stored row descends from it.

---

## Runner, issue #32 — a GPU relaxed while one core built for it ✅

**Reported as** "process utilization unoptimized": with a GPU declared for the ML rung, the
CPU side of a run — enumeration, embedding, perception, placement, hashing — used one core,
so the accelerator spent the run waiting for something to relax.

**Two independent causes, and the second was the larger one.**

*One queue, two kinds of hardware.* `work()` claimed any pending task, so a pool of N workers
was N processes competing for the same rows. The optimiser and the constructor do not want
the same hardware: construction scales with cores, and an ML relaxation on one card is one
device's worth of work however many processes ask for it. With one pool the two halves of a
run interleave rather than overlap — and on a GPU, several workers on one card divide its
memory rather than multiplying its throughput.

*The page never spawned a worker at all.* `ui/builder.submit_run` called `work()` directly on
its background thread. The worker count typed on `/builder` was written to the run row,
reported back to the page, and applied to nothing. Every run started from the browser — which
is how most of them are started, and the only way to select the GPU without a shell — was one
process, whatever the machine had been declared to be. That is where the "only one core" came
from; the queue design was merely what kept it from being fixable by declaring more workers.

**Built** `claim_task` takes `kinds` / `exclude_kinds`, so one queue can feed an
inhomogeneous pool without losing the single-row claim. `runner.plan_workers` divides the
**declared** total — never exceeds it, the workstation profile's spare core stays spare — into
builders and a relax share from `config.relax_workers`: 1 on a declared accelerator, 0 on CPU,
where relaxation is core work like everything else. `runner.execute_run` is now the single
place the pool is built, and both the CLI and the page go through it.

Relax tasks are queued *by* build tasks as they finish, so a relax worker that stopped at its
first empty claim would exit seconds into a run. `work(wait_while=…)` lets it wait; the parent
sets an event when the builders are done, which is the only thing that ends the wait.

**What building it turned up.**

* **Fork cannot serve the case this is for.** Children were forked. A forked child cannot
  re-initialise CUDA, so the moment a GPU worker was a real GPU worker it would have died at
  its first relaxation. The pool uses the **spawn** context — also the right answer for the
  web server, which is threaded, and a fork from a threaded process is its own hazard.
* **A dead worker read as a finished run.** `finish_run` called a run `done` whenever no task
  was marked `failed`. With one in-process worker that was unreachable; with a pool it is not —
  a worker that dies records nothing, so its share of the queue stays `pending` while the
  survivors' successes are all there is to summarise. Leftover pending or claimed work is now
  `failed`, and `execute_run` names the process and its exit code.
* **Two small things that only bite with N writers**, both found by starting workers
  simultaneously rather than by reasoning: the journal-mode probe used a fixed table name, so
  two connections opening at once could drop each other's probe table and conclude the
  filesystem could not hold a WAL; and connections took sqlite3's default 5 s lock timeout,
  which a dozen processes committing structures exceed as a matter of throughput and report as
  `database is locked`. Per-connection probe name; `BUSY_TIMEOUT = 30 s`.
* **Fork was hiding a test-suite fault.** A forked worker inherits the parent's memory, so a
  test that substitutes a fake energy backend in this process had it apply inside the workers
  too. A spawned one starts from nothing: database, spec JSON and environment reach it, and
  in-process state does not. That turned a pre-existing leak into a hang — `POST /api/compute`
  writes `MOFSBU_WORKERS` into the process on purpose (it is how the page reaches the work),
  no test undid it, and a later test calling `run()` therefore spawned real workers that went
  looking for a real optimiser. It passed only because `test_rerun_and_cancel` sorts before
  `test_ui_controls`. An autouse fixture in `tests/conftest.py` now undoes every `MOFSBU_*`
  declaration at the end of the test that made it.

Ground rule 8 is unchanged: an undeclared machine still runs one in-process worker and spawns
nothing. What did change is that an explicit `workers=` — a `--workers` flag, the number typed
on the page — is now honoured as the declaration it is, instead of being overridden back to
one by the laptop default.

Gates: `tests/test_jobs.py` (kind filters, the division, the relax worker's wait, and a real
three-process run), `tests/test_ui_controls.py::test_a_run_submitted_from_the_page_uses_the_declared_pool`.

---

## Runner, issue #32 follow-up — a run killed with its shell was ongoing for ever ✅

**Was** A run row records that work STARTED and that it FINISHED. Nothing writes the third
outcome — *the process stopped existing* — so a run killed with its terminal kept the last
thing it managed to say, `pending`, and `runs.html` had no choice but to believe it: `const
live = ['pending','running','cancelling'].includes(run.status)`. The "running" label beside it
came from an in-memory dict in the server process, so a restart cleared that and left the row
alone. Worse than the display: the tasks that worker held stayed `claimed`, and nothing outside
the cancel path calls `reset_stale_claims`, so they were unclaimable — the run could not be
finished by anything, ever.

Not fixable by a `finally`. SIGHUP from a closed terminal, SIGKILL, a crash and a flat battery
all skip it, and the web server's run lives on a `daemon=True` thread, which is not reliably
given the chance to run one either. The missing thing was not a cleanup path; it was a FACT —
nothing in the schema said whether a process was still there.

**Built** Two signals, because neither is sufficient alone, and they are trusted asymmetrically:

* **the claimants.** `claimed_by` is `host:pid`, so on that host the question is answerable
  exactly. A dead pid is *proof* and is acted on at once. From anywhere else — and on Windows,
  where `os.kill` terminates rather than probes — the answer is `None`, and the verdict is
  `unknown` rather than `dead`. The two-machine workflow makes that distinction load-bearing:
  guessing "dead" about the workstation would reset tasks it is running.
* **`runs.heartbeat_at`**, stamped in the same transaction as each claim and after each task.
  It only ever decides a run with NOTHING claimed. A worker inside an hour of xTB stamps
  nothing for an hour, and sweeping it would hand its tasks to a second worker and pay for the
  work twice — so a claimed task is judged by its process or not at all.

`run_liveness` returns `live` / `interrupted` / `unknown` / `unfinalised` and writes nothing,
which is what lets the read-only listing tell the truth about a worker that died while the
server stayed up. `sweep_interrupted` acts on it when a process TAKES OVER a database — a
server starting, a run starting, a resume — rather than on a timer, because that is the moment
the question is both worth asking and safely answerable.

`interrupted` is a new run status and deliberately not one of the existing ones. It is not
`cancelled` (nobody chose it) and not `failed` (the work was fine); it is the one state that
simply resumes, which is what `resume_run` now accepts. Its stranded tasks go back to
`pending`: the process holding them is gone, so they were never in flight.

**What building it turned up.** A run killed *after its last task* is not interrupted at all —
the work is all there and only the closing write was lost. That is `unfinalised`, and the sweep
closes it out with `finish_run` so it reports the status its tasks earned. Reporting it as
interrupted would have invited someone to re-run a complete run.

**Not done here:** actually draining a resumed run's queue. Resume returns the tasks and says
`resubmit the spec to execute them`, exactly as it did for a cancelled run — a run whose queue
is refilled and then left alone will read `interrupted` again once the heartbeat ages out,
which is accurate. Executing it from the page is issue #33.

Gates: `tests/test_jobs.py` (the verdicts, including the two that must NOT sweep — a live
worker and a claimant on another host), `tests/test_ui_controls.py` (a run stranded with a real
dead pid, through a server start, the read path, and resume).

---

## UI, rev 24 — seven usability defects, all closed

None of these was a correctness bug: nothing here produced a wrong structure or a wrong
number. They were all "the page does not behave the way a person expects".

Covered by `tests/test_ui_controls.py` (25 tests) plus the existing viewer suite.

### 1. Arrow keys scroll the page instead of moving through results ✅

**Was** No keyboard handling on `index.html` at all — selection was mouse-only, so the
browser's default scroll was the only thing that happened. The feature was never written.

**Built** A `keydown` handler: `ArrowUp`/`ArrowDown` by one row, `PageUp`/`PageDown` by a
screenful, `Home`/`End` to the ends. (`PageUp`/`PageDown` were later rebound to flip
between pages of results — the prev/next buttons from the keyboard — because a screenful
of a page you can already see was the less useful of the two readings.) The table has `tabindex="0"` and a focus ring, and a
click on a row focuses it so the arrows work immediately afterwards rather than needing a
Tab press nobody would guess at.

Both details the diagnosis flagged were kept: the handler ignores events from any form
control (typing "d" in the search box does not move the selection), and `preventDefault`
is applied *only* to the keys it handles, so typing is untouched.

**Also fixed, because keyboard movement made it obvious:** `select()` used to call
`load()` purely to repaint the selection highlight — a full `/api/structures` round trip
per selection. Tolerable for a click; absurd for holding down an arrow key. The highlight
is now repainted locally.

**Deliberately not done** Arrowing off the end does not fetch the next page. Pagination is
a different action with its own control, and a key that sometimes issues a query and
sometimes does not is not predictable.

---

### 2 & 3. An expanded error collapses; the server log scrolls on its own ✅

These were one problem seen from two sides, and they got one fix.

**Was** `setInterval(… 4000)` with `$('#detail').innerHTML = …`, so every `<details>` was
destroyed and recreated every 4 s — a fresh `<details>` is closed, and reading a long
traceback was a race against the timer. From the server side the same loop issued three
requests a tick (`/api/runs`, `/api/runs/{id}`, `/api/runs/{id}/tasks`), which uvicorn
logged at ~45 lines a minute from an idle open tab, forever.

**Built**, and it went further than the diagnosis proposed:

1. **Skip the DOM write when nothing changed.** The payload is hashed and an identical
   tick writes nothing at all.
2. **Preserve open state across a real re-render.** Each `<details>` carries
   `data-key="task-{id}"`; the open set is read out immediately before the write and
   reapplied immediately after. This is the one that survives a genuinely live run.
3. **The steady state is now ONE request, not three.** Rather than merely stopping the
   poll on a terminal run, the tick asks `/api/runs` only — which already carries each
   run's status and done/rejected/failed counts — and pays for the two detail requests
   only when the selected run has actually moved.

   *This is the part the diagnosis got slightly wrong.* Simply stopping the poll on a
   finished run would have meant a run submitted from another tab never appeared in the
   list. Cheapening the tick achieves the same reduction without going blind.
4. **Back off.** 4 s while moving, 15 s when finished or unchanged. A failed fetch retries
   on the slow cadence rather than hammering a server that is already unwell.

`--log-level` now defaults to `warning`, so access lines are opt-in (`--log-level info`).
The builder page's own run table went from polling every 2.5 s to every 10 s: it was never
the place you watch a run move — `/runs` is.

**Measured:** ~45 requests/minute → ~4, and the expanded traceback stays open across both
a quiet tick and a forced rebuild.

---

### 4. "22 done / 36, 14 rejected" makes the reader do arithmetic ✅

**Was** Four independent numbers, with nothing saying that **a rejection is a settled
outcome** — the task ran, the chemistry answered no, nothing is pending. A finished run
read as 22-of-36 when it was 36-of-36 settled.

**Built** The headline from the diagnosis, verbatim in spirit:

```
settled 48 / 48   ·   42 built · 6 rejected · 0 failed
```

plus the three-segment progress bar (green/amber/red), which was nearly free once the
ratio existed. The per-outcome cards stay underneath. A finished run also says *"every
task has an answer — a rejection is one of them"*, and notes when it has stopped polling.

The distinction that mattered was kept: rejected is **not** merged with failed into one
"unsuccessful" bucket. It is merged with them into one **settled** bucket, which is a
different claim — and the one the rest of the project already makes (`finish_run` returns
`done` for a run full of rejections).

---

### 5. Choosing hardware requires a shell ✅

**Was** The startup prompt made the device *visible* but not *reachable*: answering it
needs a terminal, which most of the people who open the page do not have.

**Built** `GET/POST /api/compute` and a selector on `/builder`. Possible without a new
process model for the reason the diagnosis identified: `submit_run` executes on a thread
**in this process** and `compute_device()` reads the environment at execution time.

* the device list is enumerated honestly — `cpu` always, plus each visible CUDA device
  *by its own name* (`NVIDIA RTX 4000 Ada Generation (cuda:0)`, not a bare `cuda:0`), plus
  `mps` where applicable. A machine without torch reports cpu alone rather than offering a
  device that would fail at the first structure of an hour-long run;
* the choice is recorded in the **run row** (`runs.device`, `runs.workers`), not merely
  applied — ground rule 6 for hardware. Pre-existing runs report `''`, which is honest:
  they never expressed a choice, and back-filling `cpu` would invent provenance;
* a typo is refused with the same message the CLI gives, and leaves the previous
  declaration intact rather than a broken one.

**Declared, never detected** survives: listing the GPUs and letting a person choose is a
declaration. Nothing defaults to CUDA because a card exists — the default is still `cpu`,
and there is a test asserting it does not move.

As the diagnosis anticipated, the startup prompt is now opt-in behind `--prompt-device`
(`--no-prompt` still accepted so old command lines keep working).

---

### 6. A run always goes to whichever database the server was started with ✅

**Was** `create_app(db, store)` closed over one path; `submit_run` opened `Registry(db_path)`
on it. Separating exploratory runs from real ones meant restarting with a different flag —
item 5's problem again.

**Built** `ui/active.py::ActiveDatabase`: one mutable path that the viewer, the builder and
the run inspector all read at request time, so a switch moves all three together. A page
showing structures from one database and runs from another would be worse than no switch.

All three things the diagnosis said to get right were got right:

* **the viewer is still read-only.** `ActiveDatabase` holds a *path*, never a connection.
  The viewer keeps opening `mode=ro`, the builder keeps opening its own writable ones.
  Tested explicitly after a switch.
* **name validation** — `_validate_name` applies `save_spec`'s rule (a name is a path
  *component*) plus a character whitelist.
* **row counts in the dropdown** — `registry.db — 594 structures · 38 runs`.

**Two things the diagnosis had not anticipated, both found by testing:**

* **Basenames collide.** The launch database may live outside the data root
  (`--db /elsewhere/registry.db`) and share a name with one inside it. Selecting by *name*
  was therefore ambiguous and silently resolved to the wrong file. Selection now quotes a
  full path, matched against a whitelist the server itself produced; a bare name is
  accepted only when unambiguous and returns 409 naming both folders when it is not.
* **The switch was one-way.** The whitelist was built from `current()`, so switching off a
  launch database outside the data root dropped it out of the list and the way back was a
  404 — fixable only by restarting with the flag. The launch path is now in the whitelist
  permanently. There is a regression test for exactly this.

A run **snapshots the database before its thread starts**, so switching mid-run cannot
redirect a build already in flight. The blob store is deliberately *not* switched
alongside: it is content-addressed, so two registries share one copy of a geometry.

---

### 7. No way to re-run or delete a single entry ✅

Split, as the diagnosis said to.

#### Re-run — small, and it worked out as predicted

`POST /api/structures/{id}/rerun` finds the task that built the structure, takes its
payload and its run's spec, and queues a new run of exactly one task.
`GET .../origin` answers "can this be re-run at all" *before* the button is offered, so a
seeded or ingested structure gets an explanation rather than a control that fails when
pressed.

The D2 outcome is reported as the confirmation it is: *"already present — re-running
produced the same identity, so nothing was written."*

**What testing on the real registry turned up.** Re-running structure 120 wrote a *new*
structure (595) rather than recognising the identity. Not non-determinism: same L1 graph
hash, same formula, same 49 atoms — different **multiplicity, 1 vs 4**. That is the spec
v4→v5 migration doing exactly what its comment says, `MetalSpec.multiplicity` having been
dropped because it was never kept in sync with `spin_class`, so an old spec's Co(+2,hs)
carried a stale singlet and re-deriving it now yields the reachable quartet.

So the page distinguishes three outcomes rather than two: same identity (reproducible),
*same graph but a different derived identity* (a migration showing through — stated, with
the difference named), and a genuinely different graph (flagged as worth looking at).
A flat "wrote structure 595" would have hidden the only fact that explains it.

**Consequence worth knowing:** re-running anything built before spec v5 will produce a
second row carrying the corrected spin. Both rows are real. The registry very likely
already contains such pairs.

#### Delete — became hide, as recommended

The measurement in the diagnosis held up: **457 of 594 structures — 77% — have more than
one incoming `reactions` edge**, so deleting "one entry" usually severs some other route's
history. Everything in the registry is regenerable *except* provenance.

So there is no `delete_structure` in `registry/api.py`, and its absence is the decision:

* `structures.hidden` / `hidden_at` / `hidden_reason`, filtered out of listings;
* a hidden structure keeps its row, geometries, blobs, every `reactions` edge, **and its
  L0/L1/L2 identity** — so a later run that re-derives it is still recognised under D2
  rather than inserting a duplicate;
* hiding a structure with more than one incoming route is refused (409) **with the routes
  listed**, because "is this safe to remove" is not answerable without seeing them;
  `force` overrides after confirmation;
* `include_hidden=1` and a sidebar switch bring them back. The undo is part of the
  feature, not an advanced option.

A blanket "delete any row" button was not built, as recommended.

**Compatibility:** the viewer opens `mode=ro` and so cannot migrate a registry into having
the column. Every query over it checks `has_column` first, and a registry written before
the soft delete existed simply lists normally. Tested.

---

### Packaging — a double-clickable launcher (shipped alongside, not a bug)

`launch/mofsbu.sh`, `launch/mofsbu.bat`, `launch/install-desktop-entry.sh` — one
double-clickable thing, so starting the app is not three commands.

The load-bearing decision: **the interpreter is used by full path and never activated.**
`conda activate` needs a shell that has been `conda init`-ed, which a double-click does not
provide, and it fails differently depending on how conda was installed — none of the
failures saying so clearly.

The environment is **probed, not assumed**: `$MOFSBU_PYTHON`, then `$MOFSBU_ENV`, then the
name in `environment.yml`, then *any* conda env that can `import fastapi, uvicorn, rdkit`.
That last fallback is not hypothetical — `environment.yml` says the env is called `mofsbu`
and on the development machine it is not, so a launcher hard-coding the name would have
failed on the machine it was written on.

---

---

## Identity, M4 follow-up — `site_catalog` was not a pure function of identity ✅

**Was** D15 excludes bond order from the hash, so C=O / C–O⁻ resonance forms hash identically —
but perception *read* bond order. One identity reached by two routes therefore perceived two
different donor sets, and `put_sites` keeps the first catalog, so the disagreement was silent.

**Built** `sites.perception.DELOCALISED_GROUPS` types an oxo-acid **as a whole group**: match
the central atom, take every terminal oxygen on it (ignoring metal neighbours — coordination is
not constitution), give all of them one donor type and one charge. No bond order is read
anywhere in that path, which is what makes it invariant rather than patched.

**Wider than the acetate case that surfaced it.** A sulfonate's two S=O oxygens were typed
`carbonyl_O` and only its anionic one `sulfonate_O`; nitro came out as one `carbonyl_O` and one
`alkoxide_O`. Not near-misses, and route-dependent for the same reason acetate was.

The taxonomy is deliberately coarse — carbonate is `carboxylate_O`, sulfate is `sulfonate_O`,
a phosphate diester is `phosphonate_O`. A name per oxo-acid is a promise to have anticipated
every one of them. `nitro_O` is the one genuinely new type.

`registry.catalog_drift` stays as a guard rather than a known finding;
`tests/test_sites_state.py::test_no_build_route_drifts_from_the_stored_catalog` is the gate.

**Not closed by this:** perception still counts a metal as an ordinary heavy neighbour — filed
together, fixed separately, and still open as B1.

---

## Registry, M4 follow-up — `put_sites` deleted the state it had just written ✅

`put_sites` deleted before inserting and the FK cascade took `site_state` with it. Under D2,
re-deriving an identity the registry already has is the *expected* outcome for most of an
enumeration, so this fired constantly: a structure built twice kept state only on its second
geometry, and `n_open_sites` counted against a best geometry that no longer had any.

The catalog is geometry-independent; it is now written once per structure and kept.

---

## Perception, pre-M5 — a donor was deleted at the moment it became occupied ✅

**Was** `_classify_neutral` counts heavy neighbours to tell an ether from an alcohol and a
ketone from a carboxylate — and **a metal is a heavy neighbour**. So three donor types were
perceived while free and not perceived once bound:

| ligand | free | coordinated |
|---|---|---|
| water | `aqua_O` | — |
| THF | `ether_O` | — |
| acetone | `carbonyl_O` | — |

The site left `site_catalog` at exactly the moment it became occupied, so `refresh_state` had
no row to mark `OCCUPIED` and the bond that had just formed was recorded nowhere.

**Uneven, which is why it survived this long.** N-donors and anionic donors were never
affected — they are classified by aromaticity, bond order or formal charge, none of which a
metal neighbour perturbs. A pyridine complex looked perfect; an aqua complex was empty.

**Measured on the registry before the fix:** 34 of 35 metal-bearing structures had a catalog
smaller than their own dative-bond count, and **16 had entirely empty catalogs** —
`Mg[THF]4[H2O]2` recorded 6 dative bonds and 0 perceived donors.

**Why it was caught before M5 and not during.** M5's exit gate is the anthrarufin–Cu cis/trans
pair, and an anthrarufin peri-pocket is phenolate + quinone C=O. The C=O vanished the moment
the pocket closed:

```
anthrarufin free:       phenolate_O x2, carbonyl_O x2
anthrarufin-Cu bound:   phenolate_O x2, carbonyl_O x1
```

An assembled block's `open_sites()` was therefore wrong, and `join()` inheriting sites through
the atom map would have disagreed with what `runner._record_sites` re-perceives — for exactly
the donor types assembly creates. The gate could have passed on the surviving second pocket
while the bookkeeping underneath it was wrong.

**Built** `_constitutional_heavy()` excludes metals from the heavy-neighbour count, and the
valence-rule classifier uses it. This is not a new rule: `_terminal_oxygens` already drew the
same line for the delocalised path, in those words — *coordination rather than constitution*.
The fix makes the two paths agree.

**The part worth keeping.** The fix passed all 431 pre-existing tests **unchanged**, which is
itself the finding: nothing covered the coordinated case at all. Every perception test used a
free ligand. `test_a_donor_survives_being_coordinated` is the gate — 8 pairs, each asserting
the donor multiset is identical free and bound; 4 of them fail without the fix.

**Consequence, handled by D19:** `ALGO_VERSIONS["perception"]` 1 → 2, and a `perception/1`
catalog is rewritten the next time anything touches its structure rather than being kept under
D5's write-once rule. An old catalog is not a differently-worded answer to the same question;
it is the answer to a question the current recipe no longer asks. This also closes the
"569 structures have no `site_state`" complaint from the other direction: the 16 empty catalogs
and the 16 missing-state structures were the same 16, one cause.
