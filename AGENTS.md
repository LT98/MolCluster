# Agent workflow

This repository is a Python/Conda project. Start from the repository root and read this file
before searching broadly.

## Reading order

1. `docs/CODE_ARCHITECTURE.md` — compact module map, invariants, and known seams.
2. `README.md` — setup, command entry points, data-root rules, and two-machine workflow.
3. The task-specific module listed below.
4. The focused tests before changing behavior.
5. `docs/PLAN_implementation.md` and `docs/DESIGN_registry_assembly.md` only when the task
   changes an invariant, milestone boundary, or data model decision.

## Document layout

Four active documents, each with one job. Nothing finished stays in an active document.

| | |
|---|---|
| `docs/CODE_ARCHITECTURE.md` | the map — read first |
| `docs/BUGS.md` | **open** defects only, IDs `B1…`. Fixing one means *moving* its entry to `docs/archive/BUGS_resolved.md`, not striking it through |
| `docs/PLAN_implementation.md` | **remaining** milestones and open decision gates. A milestone that meets its exit gate moves to `docs/archive/PLAN_completed.md` |
| `docs/DESIGN_registry_assembly.md` | current reasoning — decision ledger D1…, open checkpoints C…. Revisions move to `docs/archive/DESIGN_history.md` |

`docs/archive/` is reference, not scheduled work. Read it to find out why something is the way
it is; never to find out what to do next.

**Rationale goes in docs, not in the code.** A bug fix does not get its history written into
the file it touched. Source comments say what the code does and which invariant it upholds
(one or two lines; a `D…` / `B…` / ground-rule reference is ideal); the diagnosis, the
measurements, and what the fix turned up go in the archive entry. A module docstring may state
the module's contract — it is not a changelog.

Specifically, do not write: what the code *used to* be, which revision changed it, what the
old behaviour got wrong, or a re-argument of a decision the ledger already records.

Existing comments predate this rule and several modules still carry the long form. **Fix them
opportunistically** — when a task already has you editing a file, cut its bug narratives down
to the one-line form and move anything worth keeping to `docs/archive/BUGS_resolved.md`. Do not
open a file solely to trim it, and do not bundle a trim into a change that is otherwise about
behaviour without saying so in the summary.

Do not infer current branch state from old session notes. Run `git status --short --branch`,
`git fetch origin --prune`, and compare the checked-out branch with its tracking branch before
branch or merge work. A local branch may be stale even when the completed work exists in another
local ref or on `origin/main`.

## Task routing

| Task | Start here | Typical tests |
|---|---|---|
| Build planning/execution | `src/mofsbu/spec.py`, `runner.py` | `test_relax_pipeline.py`, `test_jobs.py`, `test_rerun_and_cancel.py` |
| Graph identity/canonicalization | `src/mofsbu/graph/`, `identity/` | `test_graph.py`, `test_canon.py`, `test_identity.py` |
| Donor/site behavior | `src/mofsbu/sites/` | `test_sites.py`, `test_donor_patterns.py`, `test_sites_state.py` |
| Geometry/QC | `src/mofsbu/geometry/` | `test_placer.py`, `test_distances.py` |
| Registry/schema/migrations | `src/mofsbu/registry/` | `test_registry.py`, `test_migration.py`, `test_verify.py` |
| Energy backends/references | `src/mofsbu/energy/` | `test_energy.py`, `test_reference.py` |
| Browser/API behavior | `src/mofsbu/ui/` | `test_ui.py`, `test_run_inspector.py` |

`src/mofsbu/registry/api.py` is the only write surface. Preserve graph-derived identity,
method/fidelity provenance, versioned migrations, and the distinction between absent and zero.
Missing milestone functionality must fail loudly with the repository's explicit exception
patterns rather than returning plausible placeholder data.

## Environment and validation

The preferred local environment is `ebu`:

```bash
conda run -n ebu bash scripts/check.sh
```

The gate runs the full pytest suite, registry verification, and ruff when installed. For an
iterative change, run the smallest focused test set first, then run the full gate before handoff:

```bash
conda run -n ebu python -m pytest tests/test_registry.py tests/test_migration.py -q
conda run -n ebu bash scripts/check.sh
```

If a dependency is missing, report the environment and install only the missing dependency;
do not silently switch interpreters. Warnings from optional scientific backends are acceptable
only when the relevant tests still pass.

## Search and data hygiene

Prefer searches scoped to `src`, `tests`, `scripts`, and `docs`. `legacy/` is reference material,
not the active implementation. `data/registry.db`, `data/store/`, archive directories, caches,
bytecode, and generated outputs are not source context. Use `git grep`, or `rg` when installed,
with explicit paths and globs instead of searching the repository root. The root `.ignore` keeps
common generated surfaces out of default ripgrep searches.

The registry and blob store are local/regenerable. Do not commit, delete, reset, or copy live
SQLite files or blob trees as part of a code change. Treat edits to curated files under
`data/reference/` as user data: inspect and preserve them unless the task explicitly changes
the reference dataset.

## Git safety

Do not reset, stash, rebase, force-push, or switch branches unless explicitly requested.
Do not commit unless asked. Before merging, fetch first and verify the target branch and
working-tree state. Preserve unrelated local edits and stop if a merge would overwrite them.
