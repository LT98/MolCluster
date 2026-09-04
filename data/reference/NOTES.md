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
