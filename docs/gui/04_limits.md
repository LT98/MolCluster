# What the GUI cannot do yet

Written as limits rather than omissions, because the page states most of them itself —
`/api/capabilities` carries the notes the UI renders, so the boundary is backend-declared
rather than folklore.

## Degree is capped at 1

`BuildSpec.degree` defaults to `1` and `capabilities.max_degree_implemented` is `1`.
Degree 1 means *complete one centre's coordination sphere*. Polynuclear and extended
growth are not submittable from the page.

**The note the page shows is now partly stale, and it is worth correcting.** It reads:

> degree 2+ = polynuclear / extended growth — not implemented, needs M5 (assembly.join)
> and M6 (multi-centre placer).

M5 landed, and M6 is partway: `join`, `join_bridge` and now `join_metal_metal` build
genuine polynuclear products, and `enumerate_constructions` grows onto them. What is still
missing is `place_multicentre` (M6/S3 — reconciliation where two determinants fix one
M···M and disagree) and the S7 wiring that makes a polynuclear node extractable and
storable. So the cap is real, but its stated reason now under-reports what exists.

## The S4 work is not directly exposed

M6/S4 added a **declared nucleus** (`join_metal_metal`, an M–M distance the caller states)
and **reserved vertices** (`place_mononuclear(reserve=)`). Neither has a control on the
builder page:

- There is no field for a metal–metal distance, which is consistent — the function refuses
  to invent one and the page has nowhere to declare it.
- There is no vertex-reservation control, and there does not need to be: the runner now
  reserves a mutually-cis set automatically whenever a rung leaves two or more vertices
  open. That is the S4 capability reaching the GUI **as behaviour rather than as a
  control**, and it is why co-ligand ladder runs stopped producing `chelate_cannot_span`.

## The pathway ladder stops at 12000 tasks

`runner.MAX_PATHWAY_TASKS` caps the whole plan once the ladder walk starts. When it bites,
the rungs below some products are **not planned** — their steps and intermediates are
missing, so part of the graph is not one connected ladder. It is never silent:
the diagnostic carries `cap` and `uncapped_tasks` (the ladder's full size, from the same
walk run uncapped on a copy), `estimate` returns it under `capped`, and the builder page
leads its estimate with it.

The co-ligand range (D26) is the setting most likely to reach it. Measured on
`data/reference/spec_ni_thq_cl_slice.json` re-read as a v8 spec:

| `co_ligand_counts` / window | tasks | capped? |
|---|---|---|
| `fill` (what the v7 file plans) | 722 | no |
| `range`, window 0 (no steps) | 243 | no |
| `range`, window 1 | 1141 | no |
| `range`, window 2 (**default**) | 2041 | no |
| `range`, window 3 | 2521 | no |
| `range`, window 2, CN `4,6` | 3059 | no |

The cap is 12000: the salt study's acetate spec (`spec_salt_nioac2_thq.json`, tHQ up to
two protons, CN `4,6`, window 2) plans 11,064 tasks and fits uncut. Past it: narrow `max_distinct_ligands` or the CN list,
lower the window, or use `fill`.

## Placer refusals leave holes in a co-ligand ladder

A rung the placer refuses is a rejected `place`, and every step onto it is rejected with
`pathway_parent_missing`. [B21](../BUGS.md#b21): an octahedral centre with chelates and two or
more empty vertices is refused, and `range` asks for exactly those rungs.

## Columns that exist but are empty

`/api/filters` returns `reserved_inactive`, naming each column and the milestone that will
populate it. Currently NULL in a fresh database:

| Column | Table | Filled by |
|---|---|---|
| `donor_types` | structures | M4 |
| `n_open_sites` | structures | M4 |
| `binding_modes` | structures | M5 |
| `choice_vector_json` | geometries | M5 |
| `depth` | reactions | M8 |

Their filters keep strict SQL semantics — NULL never matches — so the control is disabled
with a reason rather than quietly returning nothing.

## Metal-free runs stop early

Metals are optional and a purely molecular construction is a first-class case, but joining
molecule to molecule is M5 work that is not on this path yet. A metal-free run currently
produces **activation states and sites only**.

## Run modes depend on what is installed

`dft_go` has no backend wired. `ml_go` and `xtb_go` are offered only when their backend
imports here — `mode_status()` decides, and the page reports `available: false` with a
note rather than failing at submit time. The `ml_go` model picker lists only installed
models.

## The energy graph does not search, score or rank

The graph page composes routes **you** walk. It will not find them for you, and it puts no
route above another.

- **No path search.** There is no "show me the cheapest route to this structure", and no
  "find me an exchange". Every step is chosen, one hop at a time, from what the registry
  records (plus derived splits when asked). Walking back to a shared parent and forward again
  composes an exchange only because you walked it.
- **No scoring, no ranking, no barrier.** The legend reports each route's total, its worst
  single step and its rung. It does not order routes, call one better, or estimate a barrier.
- **A pivot is not a barrier.** Where a composed route turns, the node's `y` is how the
  composition is booked — the bare parent plus every ligand on both sides — not a height the
  exchange passes over. Nothing on the page says whether the exchange is dissociative. A pivot's
  `y` must never feed a barrier proxy (C6), and the legend gives no worst step for such a route.

That is a boundary rather than an omission. Ranking needs a barrier proxy and a decision about
partner dependence — **C6** and **C7** in `DESIGN_registry_assembly.md` — and both are
explicitly scheduled for M8, at the point where the code first forces the call. A drawing tool
pre-empting them would be answering a question nobody has decided how to ask.

What the page will say instead is when two routes *cannot* be compared: it tracks what each has
shed and badges a route whose accounting differs from the others'. [B19](../BUGS.md#b19) records
where that warning is still in the wrong place.

**Derived splits are two-part only.** `decompositions` finds pairs of existing structures whose
atoms and charge sum to the target; a three-way split is the ladder M8 builds. And they are a
full scan of the structures table, so they are asked for rather than always on — except at a
node whose recorded edges cite nothing, where derivation is the only answer available.

## Two things that are deliberate, not missing

- **The viewer cannot write.** Not an oversight — `mode=ro` plus `PRAGMA query_only`, with
  no write path in the module.
- **A real delete does not exist.** `hide` is a soft-delete flag; `registry/api.set_hidden`
  refuses a true delete, and listings can always bring hidden rows back with
  `include_hidden=1`.
