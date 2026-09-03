# MOF SBU — assembly + structure registry

A method to predict and guide MOF synthesis design by building structures from molecules:
activate molecules at viable sites, expand structures indefinitely by joining more building
blocks, evaluate synthesis-route viability (reagents -> product, with intermediates and
conditions), and store everything for recall and further construction.

**Architecture & decisions:** [`docs/DESIGN_registry_assembly.md`](docs/DESIGN_registry_assembly.md)
(living design doc — decision ledger D1-D13, open checkpoints, roadmap).

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
