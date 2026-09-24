# The Registry viewer — `/` — and the Run inspector — `/runs`

## Registry viewer

Read-only over `v_structures`. Parameterised SQL only; the sole identifiers that reach a
query string come from module-level whitelists (`SORTABLE`, `BRIDGE_CLASSES`).

### Browse and filter — `GET /api/structures`, `GET /api/filters`

`/api/filters` reports the *live* domain of every control, computed from the database in
front of you rather than assumed:

| Filter | Source |
|---|---|
| metals | distinct values present |
| bridge classes | only those actually present, ordered by the whitelist |
| fidelities | distinct `geometries.fidelity` |
| multiplicities, tags | distinct values |
| charge, energy, `n_metals` | min/max ranges |
| donor types | populated at M4 |
| sortable columns | the `SORTABLE` whitelist |

**Reserved columns degrade, they do not crash.** `donor_types`, `binding_modes`,
`n_open_sites` and `reactions.depth` are NULL in a current database. Their filters keep
strict SQL semantics (NULL never matches), and `/api/filters` returns
`reserved_inactive` naming each unpopulated column and the milestone that will fill it —
so the UI disables the control and tells you why, instead of hiding it.

**Hidden structures are filtered, not deleted.** `hidden` is a soft-delete flag; listings
exclude it by default and `include_hidden=1` brings it back. Nothing becomes unreachable.

### Inspect one entry

- `GET /api/structures/{id}` — the full record.
- `GET /api/geometries/{id}/xyz` — coordinates as plain text, rendered in-page by **3Dmol**.
- `GET /api/structures/{id}/spec` — the spec that produced it.
- `GET /api/structures/{id}/origin` — its provenance.
- `GET /api/structures/{id}/xyz` — the best geometry's coordinates, without the record.
  It exists because `GET /api/structures/{id}` prices every incoming route and is far too
  dear to call at hover rate.

### The provenance panel

Incoming edges are a table — source / added / dE / rung / caveat / count. The source cell is
the edge as an equation, reagents on the left of the arrow and the product on the right, and
every species is a link that opens that structure whether or not it is on the current page.
Hovering a link previews its geometry.

Edges with the same reagents and the same choice vector are one edge reached by different
orderings of the same choices, so they collapse to a single row carrying a count; opening the
row lists the edges behind it. **The collapsing is in the browser** — the API still returns
every edge, and the route count on a listing is the uncollapsed one.

A dE always carries its reference-scheme caveat. Every assembly edge in this project is
non-isodesmic, and a number shown without that caveat would assert what the scheme exists to
refuse. An edge that cannot be priced shows `— (reason)`, never a blank.
- `GET /api/meta` — which database and blob store are being served, that it is read-only,
  and **the `algo_versions` table as stored**. That last one is how you tell a row built
  under `placement 2` from one built under `3`.

### Two write actions, and they live on the builder router

`POST /api/structures/{id}/rerun` and `POST /api/structures/{id}/hide`. Both are on the
*builder* router, not the viewer — the viewer genuinely cannot write. `hide` is the
soft-delete; a real delete is refused (see `registry/api.set_hidden`).

## Run inspector

`GET /api/runs`, `GET /api/runs/{id}`, `GET /api/runs/{id}/tasks`,
`POST /api/runs/{id}/cancel`, `POST /api/runs/{id}/resume`.

**Stop and resume.** *Stop* is cooperative: workers finish the task in hand and stop claiming,
and the run becomes `cancelled`. If nothing is executing the run — no executor in this server
and no live worker holding its tasks — there is nobody to notice the flag, so the stop closes
the run out at once ("stopped — nothing was executing this run") instead of saying *stopping*
for ever. *Resume* puts the cancelled tasks back in the queue — plus any that failed only
because the database was busy — and **executes them** with the spec stored on the run, the same
way a submit does. It will not start a second executor on a run this server is already running,
or one whose tasks are held by live workers elsewhere. A run from an earlier session (killed
shell, closed laptop) is resumable the same way: the start-up sweep marks it `interrupted`.

What makes this page the useful one:

- **Refusals are grouped by code, then drillable.** The summary says `×212 qc_clash`; the
  tasks endpoint takes `code=` so you can look at exactly those 212. Same for `status=`.
- **A refusal is an answer, not a breakage.** The pipeline distinguishes *rejected* from
  *failed*: a bidentate ligand that cannot span a linear centre is sound chemistry
  correctly refused, and a run full of that does not read as a run full of bugs. Codes you
  will see include `qc_clash`, `placer_refused`, `join_refused`,
  `pathway_no_open_vertex`, `pathway_parent_missing`, `co_ligand_no_donor`,
  `donor_index_stale` and `chelate_cannot_span`.
- **"never queued (planner)"** — a section for what the planner decided not to emit, so
  the absence of a product is explained rather than silent.
- **What a run did not recompute.** "reused (skipped)" counts builds that landed on an
  identity the registry already had; "relaxes skipped (already computed)" counts builds whose
  relaxation already existed at this level of theory (`detail.relax_reused` names the
  geometry). A second run into the same database — the NiCl₂ then Ni(OAc)₂ salt study — shows
  its saving here rather than as a smaller task count.
- **Time per task.** Each task records its own wall time (`detail.duration_ms`); the table
  gives n, total seconds, median and p90 per kind. Rejected tasks count, because their time
  was spent. Tasks from before this was recorded are left out rather than counted as zero.
- **"after the run"** — what `runner.finalise_run` did once the queue drained: for a
  `pathways` run, how many one-proton deprotonation edges it wrote and why any were refused.
- **Liveness is computed, never stored.** The page asks whether a pid still exists rather
  than believing the row — because a row is written by a process that, if it died, is by
  definition no longer able to correct it. Finished runs show no liveness verdict at all,
  so the ones that matter stay visible.

### What M6/S4 changed here

`chelate_cannot_span` is the code S4 was aimed at. A CN-6 octahedral ladder run with a
co-ligand and `pathways` on used to show **two rejected grow steps** under that code, with
a hole in the chain. It now shows **none**, and all three steps done. Nothing about the
page changed — the refusals simply stopped happening.
