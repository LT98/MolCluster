# UI backlog — known annoyances, diagnosed

Low priority by the user's own assessment, and none of them is a correctness bug: nothing
here produces a wrong structure or a wrong number. They are all "the page does not behave
the way a person expects", which is why they are written down rather than fixed in the
middle of chemistry work — and why each entry carries a **diagnosis**, not just a
complaint. The expensive part of a small fix is working out what causes it; that part is
done here.

Ordered by (annoyance ÷ effort), worst ratio first.

---

## 1. Arrow keys scroll the page instead of moving through results

**Where** `ui/static/index.html`

**Diagnosis** There is no keyboard handling on that page at all — `grep -n "keydown"` returns
nothing. Selection is mouse-only (`tr.onclick = () => select(row.id)`, index.html:286), so
the browser's default scroll behaviour is all that happens. Nothing is broken; the feature
was never written.

**Fix** A `keydown` listener on the results table: `ArrowDown`/`ArrowUp` move the selection
by one row, `preventDefault()` to stop the scroll, and `scrollIntoView({block:"nearest"})`
so the selected row stays visible without the page jumping. `Home`/`End` for first/last are
nearly free once the handler exists. `PageDown`/`PageUp` should move a screenful rather
than paginating, because pagination is a different action with a different control.

Two details worth getting right rather than discovering later:
* the handler must ignore key events originating in the filter inputs, or typing "d" in a
  search box moves the selection;
* the table needs `tabindex="0"` and a visible focus ring, or keyboard users cannot reach
  it in the first place.

**Effort** Small. One handler, ~25 lines.

---

## 2. An expanded error detail collapses every few seconds

**Where** `ui/static/runs.html:334`, `:124`, `:317`

**Diagnosis** The run inspector polls on a fixed timer —
`setInterval(() => { loadRuns(); if(state.runId !== null) render(); }, 4000)` — and
`render()` replaces the whole panel with `$('#detail').innerHTML = ...`. The
`<details><summary>full record</summary>` element (runs.html:317) is therefore *destroyed
and recreated* every 4 s, and a fresh `<details>` defaults to closed. Reading a long
traceback is a race against the timer.

**Fix**, cheapest correct version first:

1. **Do not re-render when nothing changed.** Hash the payload (`run` + `tasks`) and skip
   the DOM write when the hash matches the last one. This alone fixes the collapse for a
   finished run, which is the common reading case, and it also fixes item 3 below.
2. **Preserve open state across a real re-render.** Give each `<details>` a stable key
   (the task id), collect `[...$('#detail').querySelectorAll('details[open]')].map(d => d.dataset.key)`
   before the write and re-apply `open` after it.
3. **Stop polling a finished run.** When `run.status` is terminal and `thread !== 'running'`
   there is nothing left to poll for; clear the interval and leave the manual `refresh`
   button (runs.html:194) as the way back.

(1) and (3) are each a few lines. (2) is the one that survives a genuinely live run.

**Effort** Small–medium. Do (1) and (3) together; add (2) if it still annoys.

---

## 3. The server log scrolls constantly with no user action

**Where** same `setInterval`, `ui/static/runs.html:334`

**Diagnosis** Not a backend behaviour at all — it is the same 4-second poll seen from the
other side. Each tick issues `/api/runs`, `/api/runs/{id}` and `/api/runs/{id}/tasks`, and
uvicorn logs one line per request, so an idle open tab produces ~45 log lines a minute
forever. It maps to "a run inspector is open somewhere", which is why it looks unrelated
to anything you are doing.

**Fix** Item 2's (3) removes the steady state entirely. Beyond that, back the interval off
while the run is active but unchanged (4 s → 15 s after a few identical polls), and
consider `uvicorn --log-level warning` for the default launch so the access log is opt-in.
A single `/api/runs/{id}/summary` returning all three payloads would cut the line count
threefold, but that is a bigger change than the annoyance justifies.

**Effort** Small, and mostly the same edit as item 2.

---

## 4. "22 done / 36, 14 rejected" makes the reader do arithmetic

**Where** `ui/static/runs.html:140-148` (the `cards` array)

**Diagnosis** The cards report `done`, `rejected`, `failed` and `pending` as four
independent numbers. Nothing says that **a rejection is a settled outcome** — the task ran,
the chemistry answered "no", and there is nothing further to wait for. So a finished run
reads as 22-of-36 complete when it is in fact 36-of-36 settled, and the reader has to add
up the other cards to discover that. The project already draws this distinction correctly
everywhere else (`rejected` is deliberately not `failed`, and `finish_run` returns `done`
for a run full of rejections) — the UI just does not show it.

**Fix** A headline above the existing cards:

```
settled 36 / 36   ·   22 built · 14 rejected · 0 failed
```

where `settled = done + rejected + failed` and the total is `settled + pending + claimed`.
Keep the per-outcome cards underneath — they are the useful breakdown once you know the
run is finished. The progress bar the user asked about is then almost free, since the ratio
is already computed: one `<div>` with a percentage width, coloured by whether `failed > 0`.

Worth keeping: a rejection is not a failure, and the summary must not merge them into one
"unsuccessful" bucket. It should merge them into one **settled** bucket, which is a
different claim.

**Effort** Small. The numbers are all already in `outcome_summary`.

---

## 5. Choosing hardware requires a shell, which rules out most of the group

**Where** `ui/__main__.py::confirm_compute_settings`, `config.compute_device`

**Diagnosis** The startup prompt made the device *visible*, which was the previous problem,
but it did not make it *reachable*: it still assumes the person launching knows what a
terminal is. For a labmate who opens the app and uses the browser, there is currently no
path to the GPU at all.

**Fix** Put the control in the page, which is possible for a reason specific to this
architecture: `submit_run` executes the run on a background thread **in the same process**,
and `compute_device()` reads the environment at execution time. So a device selector in
`/builder` that sets the value for subsequent runs works without any new process model.

Shape of it:
* `GET /api/capabilities` already reports installed backends; extend it with the available
  devices (`cpu`, plus `cuda:N` for each visible GPU, `mps` where applicable) and which one
  is currently declared;
* a `<select>` next to the run-mode radio buttons, defaulting to the declared value;
* the choice is recorded **in the run's spec/metadata**, not just applied, so a stored run
  still says which device produced it. That is ground rule 6 applied to hardware, and it is
  the part that must not be skipped for convenience.

The one thing to preserve: **declared, never detected.** Listing the GPUs and letting a
person choose is still a declaration. Silently defaulting to CUDA because a card exists is
not, and would undo the reason the prompt was added.

Once the page can do it, the startup prompt should probably become opt-in
(`--prompt-device`) rather than the default, since it will then be the second-best way to
answer the same question.

**Effort** Medium — small on the page, slightly more to enumerate devices honestly and to
thread the value into the run record.

---

## 6. A run always goes to whichever database the server was started with

**Where** `ui/__main__.py::resolve_database`, `ui/builder.py::submit_run`

**Diagnosis** The database is chosen once, at launch, by `--db` or by discovery, and then
baked into the app: `create_app(db, store)` closes over it and `submit_run` opens
`Registry(db_path, ...)` on that single path. Nothing in the page can point a run
somewhere else, so separating exploratory runs from real ones means restarting the server
with a different flag — i.e. a shell, which is item 5's problem again.

The machinery to *find* the alternatives already exists and is already exposed:
`discover_databases()` lists every `.db` under the data root, newest first, and
`--list-db` prints them. It is only the selection that is missing.

**Fix** A database `<select>` on `/builder`, populated from `discover_databases()`, plus a
"new database…" option that takes a name and creates it. Creation is genuinely cheap —
`Registry(path).migrate()` builds a complete schema from nothing, which is how every
registry in the project has ever been made.

Three things to get right:
* **the viewer is read-only and must stay that way.** The registry page opens its
  connection with `open_read_only()`; adding a writable path for runs must not widen that.
  The write already happens on a different connection inside `submit_run`, so this is a
  matter of not "simplifying" the two into one later.
* **name validation**: `save_spec` already refuses a filename containing `/` or `\`, and
  the same check belongs here — a database name is a path component, not a path.
* **show the row count** next to each database in the dropdown. "Which one is my real
  one" is the actual question being asked, and a name alone frequently does not answer it.

**Effort** Small–medium. The listing exists; it is a select, a POST, and being careful
about the read-only connection.

---

## 7. No way to re-run or delete a single entry

**Where** `registry/api.py`, `ui/app.py`

**Diagnosis** Two different requests that look like one.

**Re-run** is nearly free and mostly exists: a geometry stores its `choice_vector` and
`seed`, and the whole point of D13 is that `construct` is a deterministic function of the
two. `/api/structures/{id}/spec` already returns the spec that produced a structure, for
copy-paste. Re-running one entry is therefore "take that spec, narrow it to this one
structure, queue it" — which is a button over machinery that is already there. Under D2 it
will usually recognise the identity and write nothing, which is the correct outcome and
should be *reported* as such rather than looking like a no-op.

**Delete** is the one to be careful about, and it is not a UI problem. The schema cascades
hard: deleting a structure takes its geometries, its `site_catalog`, its `site_state`, and
its `reactions` edges with it. That is right for a mistake and wrong for a structure some
other route also reached — `reactions` is a DAG, and **457 of 594 structures in the working
registry — 77% — have more than one incoming edge**, so deleting "one entry" usually
severs some other route's history. That is not an edge case to guard against; under D2 it
is the normal state of the registry, because re-deriving an identity you already have is
what most of an enumeration does.

**Fix** Split them:
* **Re-run**: a button on the structure page. Reuses the stored spec; shows the outcome
  including "already present, nothing written".
* **Delete**: offer it only for a structure with a single incoming provenance edge, and
  refuse (with the reason, and a list of the other routes) otherwise. Prefer a soft delete
  — a `hidden` flag filtered out of `v_structures` — over a real one: the registry is
  meant to be regenerable, but provenance that has been cascaded away is not recoverable,
  and "I deleted the wrong row" is a much worse afternoon than "the list has a filter on
  it".

A blanket "delete any row" button is the one thing here I would not build.

**Effort** Re-run: small. Delete: medium, and mostly deciding the policy rather than
writing it.

---

## Not in this list

Things noticed while diagnosing the above, deliberately excluded because they are not UI:

* `index.html` has no polling at all, so the registry page can show a stale count after a
  run finishes in another tab. Arguably correct — a registry view that moves under you
  while you read it is worse — but it is a decision nobody made explicitly.
* The builder posts `spec_version: 1` and relies on the migration chain to bring it
  forward. It works, and it means the page never has to know the current version, but it
  is load-bearing behaviour that no test covers.
