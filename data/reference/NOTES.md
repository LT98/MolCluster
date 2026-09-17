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

- 2026-09-14  node_cases.tsv  — the polynuclear nodes the M6 multi-centre placer must build,
  and the M–M / M–O / µ-oxo numbers a built one is measured against.  **This is the interface
  for adding an M6 test case:** add a row and `python -m pytest tests/test_placer_multicentre.py`
  re-checks the whole set.  Read by `tests/node_cases.py`.
  Both machines: nothing to configure — tracked, like the rest of this folder.
  **Provenance, stated in the file's own header and repeated here because it matters:** these
  are values TYPICAL of the named compound class, carried from the structural-chemistry
  literature.  No CIF was consulted and no single refinement is being reproduced.  That is why
  every row carries a window rather than only an ideal, and why the tests treat them as a shape
  check.  Before any of these numbers backs a quantitative claim, check the ideal against the
  CSD and narrow the window in the same commit.
  Three rows (Cr/Rh/Mo paddlewheels) have `fixture = -`: they exist so the M–M target is not a
  two-point table, and Rh and Mo are deliberately there as metals `geometry.distances.BASE_MO`
  does not know.
