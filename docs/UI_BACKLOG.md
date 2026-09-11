# UI backlog — all seven done

Every item diagnosed in the previous revision of this file has been implemented. This
version keeps the diagnoses, because they are the expensive part and they explain the
shape of each fix, and records what was actually built underneath each one — including
the three places where building it turned up something the diagnosis had not.

None of these was a correctness bug: nothing here produced a wrong structure or a wrong
number. They were all "the page does not behave the way a person expects".

Covered by `tests/test_ui_controls.py` (25 tests) plus the existing viewer suite.

---

## 1. Arrow keys scroll the page instead of moving through results ✅

**Was** No keyboard handling on `index.html` at all — selection was mouse-only, so the
browser's default scroll was the only thing that happened. The feature was never written.

**Built** A `keydown` handler: `ArrowUp`/`ArrowDown` by one row, `PageUp`/`PageDown` by a
screenful, `Home`/`End` to the ends. The table has `tabindex="0"` and a focus ring, and a
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

## 2 & 3. An expanded error collapses; the server log scrolls on its own ✅

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

## 4. "22 done / 36, 14 rejected" makes the reader do arithmetic ✅

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

## 5. Choosing hardware requires a shell ✅

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

## 6. A run always goes to whichever database the server was started with ✅

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

## 7. No way to re-run or delete a single entry ✅

Split, as the diagnosis said to.

### Re-run — small, and it worked out as predicted

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

### Delete — became hide, as recommended

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

## Packaging (not previously in this list)

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

## Still not in this list

Unchanged from the previous revision, and still deliberately excluded:

* `index.html` has no polling, so the registry page can show a stale count after a run
  finishes in another tab. Arguably correct — a registry view that moves under you while
  you read it is worse — but it is a decision nobody made explicitly.
* The builder posts `spec_version: 1` and relies on the migration chain to bring it
  forward. It works, and it means the page never has to know the current version, but it
  is load-bearing behaviour that no test covers.

Noticed while implementing, and *not* fixed because it is outside this work:

* The run inspector's "what was attempted" column renders `molecule undefined ·
  undefined×0-dentate` for `place` tasks. The payload keys the JS reads (`p.molecule`,
  `p.donors`) are not the keys the planner writes. Cosmetic, on a column that is otherwise
  the most useful one on the page.
