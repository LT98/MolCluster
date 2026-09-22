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
- **Liveness is computed, never stored.** The page asks whether a pid still exists rather
  than believing the row — because a row is written by a process that, if it died, is by
  definition no longer able to correct it. Finished runs show no liveness verdict at all,
  so the ones that matter stay visible.

### What M6/S4 changed here

`chelate_cannot_span` is the code S4 was aimed at. A CN-6 octahedral ladder run with a
co-ligand and `pathways` on used to show **two rejected grow steps** under that code, with
a hole in the chain. It now shows **none**, and all three steps done. Nothing about the
page changed — the refusals simply stopped happening.
