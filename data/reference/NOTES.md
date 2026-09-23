# Reference data — selective, git-tracked

`data/` is git-ignored by default (regenerable registry DB, blob store, compute outputs).
This `data/reference/` subfolder is the ONE place under `data/` that IS tracked, for data you
deliberately want on both machines.

## How to promote a file into git
    git add data/reference/<file>        # normal text/data (.xyz, .csv, .json)
    git add -f data/reference/<file>     # if it's caught by a global ignore (e.g. a .png)

## Log each promotion here
For every dataset added, note what it is and anything the OTHER machine must change to use it
(paths, env vars, a package version, a manual step). Example:

- 2026-09-03  fe_btc_node_reference.xyz  — DFT-relaxed Fe3-oxo node.
  Laptop access: none needed. Workstation: produced by scripts/run_node.py at commit <sha>.

- 2026-09-04  ligand_cases.tsv  — reference cases for donor perception, protomer
  enumeration and chelate detection.  **This is the interface for adding test cases:**
  add a row, then `python scripts/check_cases.py` (or `python -m pytest`, which runs the
  same rows as parametrised tests).  No code changes needed for a new case.
  Both machines: nothing to configure — it is tracked, unlike the rest of `data/`.
  Two rows are marked EXPECTATION CORRECTED: they failed on first run because the
  expectation was wrong, not the code (oxalic acid has 3 protomers not 4; a carboxylate
  does chelate through its own two oxygens).  Corrections are annotated in the row so the
  reasoning is not lost.

- 2026-09-17  node_cases.tsv  — the **M6 battery** (`docs/WORKPLAN_M6.md` §6) as data: the
  polynuclear nodes the milestone must reach, which bridge mechanism each one exercises, and
  the numbers a built one is measured against.  **This is the interface for adding an M6 test
  case:** add a row and `python -m pytest tests/test_m6_battery.py` re-checks the whole set.
  Read by `tests/node_cases.py`.
  Both machines: nothing to configure — tracked, like the rest of this folder.
  **Two kinds of number live in it and they are not interchangeable**, which the file header
  says at length and which is repeated here because it decides what the values may be used for:
  * **measured** (`source` starts with `WORKPLAN_M6`) — from the battery run recorded in that
    file, produced by this package's own geometry.  `d_mm_bridge` and the µ-oxo skeleton
    columns are of this kind.
  * **literature-typical** (`source` names a compound) — values typical of the named compound
    class.  No CIF was consulted and no single refinement is reproduced.  That is why every row
    carries a window rather than only an ideal.  Before any of them backs a quantitative claim,
    check against the CSD and narrow the window in the same commit.
  Under **D20** the M···M distance is an output to *validate*, not an input to impose, so these
  are what gate 4 measures a built node against — not what drives a placer.  The exception is a
  declared nucleus (S4), where the caller states the distance on purpose.
  Rows with `blocked_on` set cannot build yet ([B13](../../docs/BUGS.md#b13) for the hydroxide
  bridge, [B14](../../docs/BUGS.md#b14) for the pyrazolate one, and the missing bent CN-2
  geometry); a test checks each stated blocker is still real, so the table cannot claim to be
  waiting on something already fixed.  Three rows (Cr/Rh/Mo paddlewheels) are not battery rows
  at all — they are there for S4's `metal_metal_distance` table, and Rh and Mo are deliberately
  metals `geometry.distances.BASE_MO` does not know.

- **`spec_salt_nicl2_thq.json`, `spec_salt_nioac2_thq.json`** (and `_k1` fallbacks) — the
  NiCl₂ vs Ni(OAc)₂ salt study: which nickel salt forms THQ complexes more readily.
  Both machines: nothing to configure. Run both into one database, chloride first, so the
  second run reuses every THQ/water-only structure the first one built.
  **Input rule — every species is entered in the form it exists in water under the study
  conditions, and `max_deprotonations` says how far below that form it may go.** So chloride
  is `[Cl-]` with 0 (HCl does not exist in water, pKa ≈ −6 — entering it as `Cl` with one
  deprotonation is what put Ni–ClH species into `mvp_ni_thq_cl.db`); acetic acid is `CC(=O)O`
  with 1 (both HOAc and OAc⁻ are real at pKa ≈ 4.8); tHQ is neutral with 2 (1 in `_k1`).
  THQ is listed first in both specs because placement payloads name molecules by spec index,
  and reuse between the two runs depends on those payloads matching. Square-planar is left
  out: d8 Ni(II) there is low-spin, and these specs model Ni(II) high-spin only.
