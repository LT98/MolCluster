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

## Two things that are deliberate, not missing

- **The viewer cannot write.** Not an oversight — `mode=ro` plus `PRAGMA query_only`, with
  no write path in the module.
- **A real delete does not exist.** `hide` is a soft-delete flag; `registry/api.set_hidden`
  refuses a true delete, and listings can always bring hidden rows back with
  `include_hidden=1`.
