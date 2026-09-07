# MOF SBU — assembly + structure registry

A method to predict and guide MOF synthesis design by building structures from molecules:
activate molecules at viable sites, expand structures indefinitely by joining more building
blocks, evaluate synthesis-route viability (reagents -> product, with intermediates and
conditions), and store everything for recall and further construction.

**Architecture & decisions:** [`docs/DESIGN_registry_assembly.md`](docs/DESIGN_registry_assembly.md)
(living design doc — decision ledger D1-D13, open checkpoints, roadmap).

**Implementation plan:** [`docs/PLAN_implementation.md`](docs/PLAN_implementation.md)
(module map + API signatures, milestone spine M0-M9 with exit gates, salvage ledger, decision gates).

## Layout

    src/mofsbu/        the package
      graph/           typed molecular graph + canonicalization (L0-L1)
      identity/        composite L0-L3 key + hashing
      sites/           donor perception, frames + live-DOF tags (perceive-once)
      geometry/        multi-center placer, fidelity-laddered geometries
      registry/        SQLite registry + content-addressed blob store
      energy/          xTB / MACE backends
      descriptors/     per-donor / per-metal tables
      assembly/        BuildingBlock, frame-alignment join, choice-vector
      pathways/        reactions + pathway scoring
    docs/              design doc + archived study reports
    tests/  scripts/   tests and drivers
    data/              (git-ignored) registry.db, blob store, outputs
    legacy/            archived pre-rewrite code (reference / salvage)

## Setup (both machines)

    conda env create -f environment.yml
    conda activate mofsbu
    pip install -e .        # editable install of the package

## Running things

Always run from the **repo root**, never from inside `src/` or a package folder: Python puts the
current directory first on `sys.path`, and a package folder on the path shadows stdlib modules.
(Our modules are named `_types.py` rather than `types.py` for exactly this reason.)

    conda activate mofsbu
    python -m pytest                              # the gate — must pass on both machines
    python scripts/seed_demo_registry.py          # demo data for the viewer
    python scripts/viewer.py                      # the viewer — no install, no arguments

The `scripts/` launchers put `src` on the path themselves, so a fresh clone works immediately.
For the shorter commands, install the package once per machine:

    pip install -e .
    python -m mofsbu.ui
    mofsbu-ui

Both find the registry themselves: the real one if it exists, otherwise the demo. A relative
`--db` is resolved against the working directory, then `MOFSBU_DATA`, then the repo root, so the
command behaves the same wherever it is run from.

## Reading a run (`/runs`)

The console shows a progress line, not a diagnosis. `/runs` is the inspector: what the run
planned, what it never queued and why, what it built, what it *reused* rather than built, and
what it refused — grouped by cause, with the evidence attached.

* **Causes are grouped.** Every rejection carries a machine-readable `error_code`, so 200
  refusals with one cause are one line. Click it to filter the task list.
* **A refusal names what it refused.** A QC clash records both atoms, their elements, which
  ligand each came from, the measured distance and the limit — not just "2 clash(es)".
* **New vs reused.** `reused structure` means the task ran and wrote nothing because the
  registry already had that identity. Under D2 that is a success; it used to be
  indistinguishable from building something new.
* **Re-attempted / embed retried / UFF fallback** are badges, because each changes what a
  geometry is worth and none of them were visible before.

Metal–donor distances are per *pair*, not per metal (`geometry/distances.py`): a centre carrying
a water and an iodide has two different M–L distances. A pair the table has not calibrated is
still placed, but its distance is labelled an estimate everywhere it appears.

## Energies (M7)

`python -m pytest` passes with no quantum-chemistry stack installed: the backends are behind
a protocol, and the tests that need one skip themselves. What a machine can actually run is
reported by the code rather than assumed — `mofsbu.energy.available_backends()` and
`mode_status()` are what the builder page's disabled options are drawn from.

    conda install -c conda-forge tblite-python ase     # GFN2-xTB (also: pip install tblite ase)
    pip install torch --index-url https://download.pytorch.org/whl/cpu && pip install mace-torch

`tblite` publishes working PyPI wheels, so the pip line is a real alternative to conda —
worth knowing, because `legacy/energy_model.py` records the opposite for the sandbox it was
written in. DFT has no backend wired up and asking for one raises rather than substituting
something cheaper.

Energies are only ever subtracted through `mofsbu.energy.reference`, which refuses equations
it cannot justify — see `docs/PLAN_implementation.md` rev 16 for why a *balanced* equation is
not enough. The archived Fe(III) references can be re-derived with:

    python scripts/regress_m7.py --refs

If that install fails with *"build backend is missing the 'build_editable' hook"*, the environment's
setuptools predates PEP 660 — `pip install -U setuptools` and retry.

The registry and blob store default to `<repo>/data/`. If the repo lives in OneDrive, point them
elsewhere — a sync client copying a live SQLite file corrupts it, and none of this needs syncing:

    setx MOFSBU_DATA %LOCALAPPDATA%\mofsbu        # windows
    export MOFSBU_DATA=~/.local/share/mofsbu      # linux

## Two-machine workflow (git is the bridge)

Laptop = editing / light dev; Linux workstation = heavy xTB/DFT/MACE runs.
Source + docs live in git; the `data/` store and compute outputs are regenerable and stay
out of git (sync them separately or regenerate on each machine).

    # edit anywhere ->
    git add -A && git commit -m "..." && git push
    # on the other machine ->
    git pull

Notebooks are kept out of git by default (they merge badly across machines); prefer `.py`
modules + `scripts/`. Use `git add -f` for a notebook you deliberately want tracked.
