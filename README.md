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
      energy/          xTB / MACE-MP-0 / MACE-OMOL-0 backends
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

## Starting the viewer (no terminal needed)

For anyone who just wants to use the app, there is one file to run:

    ./launch/mofsbu.sh              # Linux / macOS — double-click, or run it
    launch\mofsbu.bat               # Windows — double-click

It finds the conda environment itself, starts the server and opens a browser. To get a
menu entry (or a desktop icon) so it never has to be found in a file manager:

    ./launch/install-desktop-entry.sh
    ./launch/install-desktop-entry.sh --gpu       # a second shortcut that uses the GPU

The launcher takes the two options worth having up front, and passes anything else
through to `scripts/viewer.py`:

    ./launch/mofsbu.sh --gpu --workers 8
    ./launch/mofsbu.sh --db scratch.db --port 8001

It does **not** run `conda activate`. The environment's interpreter is used by full path,
because `activate` needs a shell that has been `conda init`-ed and a double-click does not
give you one. The environment is found by trying the name in `environment.yml`, then
`$MOFSBU_ENV`, then **any** conda env that can import `fastapi`, `uvicorn` and `rdkit` —
so an environment created under a different name still works. Override it directly if you
prefer: `MOFSBU_PYTHON=/path/to/envs/NAME/bin/python ./launch/mofsbu.sh`.

Device and database are also changeable **inside the page** (`/builder`), so neither one
needs a restart or a flag. See "Choosing hardware and database" below.

## Running things (development)

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

**Settled, not done.** The headline reads `settled 36 / 36 · 22 built · 14 rejected · 0 failed`.
A rejection is a *settled* outcome — the task ran and the chemistry answered no — so a run full
of rejections is finished, not stuck at 22-of-36. It is deliberately not merged with `failed`:
those are different claims, and `finish_run` returns `done` for a run full of rejections.

**The page stops asking when there is nothing to ask about.** A finished run is a document; it is
not re-fetched or repainted, so an expanded traceback stays expanded. While a run is live the
panel is rebuilt, and open sections are restored by task id. The steady state is one request
every 15 s instead of three every 4 s — which is why the console no longer scrolls on its own.
The server also defaults to `--log-level warning`; pass `--log-level info` to get access lines back.

## Choosing hardware and database (`/builder`)

Both used to be decided once on the command line, which put them behind a terminal.

* **Device.** A selector lists what this machine actually has (`cpu`, plus each visible CUDA
  device by name, `mps` where applicable) and what is currently declared. It stays a
  *declaration*: enumerating the GPUs and letting a person choose is not detection, and nothing
  ever selects CUDA because a card is present. The chosen device is written into the run row,
  so a stored run still says which hardware produced it. It applies process-wide, so two runs
  going at once share one declaration.
* **Database.** A dropdown lists every registry under the data root **with its row count** —
  "which one is my real one" is the actual question, and two plausible filenames do not answer
  it. "new…" creates an empty one. Switching moves the viewer, the builder and the run inspector
  together. The viewer stays read-only across a switch; only paths the server already listed can
  be selected.

## Re-running and hiding one entry

* **Re-run** replays the task that built a structure — its payload and its run's spec are both
  stored, and `construct` is a deterministic function of the choice vector and seed (D13). The
  usual result is "already present, nothing written", which is reported as the confirmation it
  is. If the stored spec predates a migration the replay can land on a *different* identity
  (same graph, different derived multiplicity, say); the page says so rather than reporting a
  bare new row.
* **Hide** is this project's delete, and it is soft on purpose. The schema cascades hard —
  deleting a structure takes its geometries, `site_catalog`, `site_state` and its `reactions`
  edges — and 457 of 594 structures in the working registry (77%) have more than one incoming
  edge, so deleting "one entry" usually severs some other route's history. Everything here is
  regenerable except provenance. So a hidden structure keeps its row, its blobs, its edges and
  its L0/L1/L2 identity (it is still recognised under D2) and simply leaves the listing.
  Hiding a multi-route structure is refused until confirmed, and the refusal names the routes.
  "show hidden structures" in the sidebar is the way back.

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

Three more things are **declared, never detected**, so a build cannot decide on its own to
seize hardware or to switch theories:

    MOFSBU_PROFILE   laptop | workstation        # how many workers may run
    MOFSBU_DEVICE    cpu | cuda | cuda:N | mps   # where the MLIP runs
    MOFSBU_ML_MODEL  mace-mp-0 | mace-omol-0     # WHICH MLIP; default mace-mp-0

`ml_go` names a rung of the fidelity ladder, not a theory. **MACE-MP-0** (Materials
Project) is blind to formal charge and spin, so the reference scheme refuses it on any
charged equation; **MACE-OMOL-0** (OMol25, wB97M-V/def2-TZVPD) is given the total charge
and spin multiplicity and is accepted. Their energies are on different scales and are
never subtracted from one another — the `methods` row records which model produced each
number, and nothing compares energies across method rows. A spec may pin the model
(`ml_model`), and that beats the environment. MACE-OMOL-0 needs `mace-torch>=0.3.14`.

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
