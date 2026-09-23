# The Builder page — `/builder`

*"Writes a spec and queues a run; it never writes structures."* Four numbered sections on
one page, plus two side panels.

The page renders its options from `GET /api/capabilities` rather than from hard-coded
markup, so an unimplemented choice is labelled by the backend and cannot drift out of date
with the code.

## 1 · Molecule

- Enter a SMILES and **preview it before committing**: `POST /api/molecule/preview`
  embeds it, perceives its donor sites, and reports what it found. Capped at 300 atoms.
- Pick the **chelate pocket predicate** — ring size and number of anionic donors — which is
  what selects *which* pocket of a polydentate ligand is used.
- Protomer handling is part of the spec, so deprotonation is an enumerated branch rather
  than an assumption.

## 2 · In this run

- **Metals** — symbol, oxidation state, spin class. Metals are **optional**: a purely
  molecular construction (COF, organic cage) is a first-class case. Metal-free runs
  currently produce activation states and sites only.
- **Coordination number** and **ligands per metal**, both of which accept **ranges**:
  `1~3` means 1, 2 and 3; `1~3, 6` adds 6. A spec stores the expanded list. One range may
  span `MAX_RANGE_SPAN` values, and the page says so instead of leaving it to be found.
- **Geometries** — chosen from the real `GEOMETRIES` table, so `bent` and both CN-5
  polyhedra are offered because the placer has them.
- **Binding modes**, and an optional **co-ligand** given as the species that actually binds
  (`O` for water, `[OH-]` for hydroxide, `[I-]` rather than `I2`). A co-ligand with no
  perceivable donor is refused by name with that hint, not silently dropped.
- **Co-ligand count** (`co_ligand_counts`, D26) — applies to whatever co-ligand or solvent
  is named (water is only the default). `fill` puts one on every vertex the ligands leave;
  `range` (the default for a new spec) also builds the counts down to the **window** below
  full (default 2), the uncovered vertices left **empty in the same polyhedron** — a CN-4 Ni
  with one ligand and two waters is still a square plane with one open vertex, never a
  3-coordinate complex. A count above full needs a higher CN, so it appears only when the CN
  list has one; none is invented. With a co-ligand set, `allow_unsaturated` does not gate
  these — choosing `range` is the request for them. A spec saved before v8 loads as `fill`,
  which is what it planned.
- **Pathways** toggle — record how the rungs of a ligand-count sweep reach each other. The
  intermediate below each product is built and the step between them performed with
  `assembly.join`, so the registry holds the *route* and not just the endpoints. It adds
  the coordinatively unsaturated intermediates to the run. Under `range` a co-ligand is a
  step too — `Ni(H2O)2 + H2O → Ni(H2O)3`, the free co-ligand built and relaxed as a
  reagent — and no rung leaves more than the window's number of vertices empty, so the
  seeds of one metal and CN are one ladder on the graph page, rooted at the lowest
  co-ligand state inside the window. No bare-metal row is made. A window of 0 plans no step.

## 3 · Construction

- **Run mode**, one of four rungs, each reported with whether it is actually available
  here and which *theory* would serve it: `construct` (geometry only), `ml_go`
  (MACE-MP-0 / MACE-OMOL-0), `xtb_go` (GFN2-xTB), `dft_go`.
- For `ml_go`, the **model picker is built from what is installed** — a model that is not
  present cannot be selected, and one that is gets named rather than called "MACE".
- **Estimate before submitting**: `POST /api/estimate` runs the same enumeration
  submitting would, so the number on the page is the number you get. There is a test
  pinning those two together, because otherwise the number is decoration. When the
  pathway ladder hits its task cap (`MAX_PATHWAY_TASKS`, 4000) the estimate leads with a
  red line giving the cap and the size the uncapped ladder would have been — the run is
  then **incomplete**, not merely large (see [limits](04_limits.md)).

## Side panels

**Library** (`GET`/`POST /api/library`) — save molecules you reuse.
**Spec files** (`GET /api/specs`, `GET /api/specs/{source}/{name}`, `POST /api/spec`) —
write, list and reload named specs from `spec_dir`, so a run is reproducible from a file.

## Hardware and database, declared from the page

- `GET`/`POST /api/compute` — enumerate the devices that exist here and **declare** which
  one this server uses. Never chosen for you.
- `GET`/`POST /api/databases`, `POST /api/database/select` — create a registry, list them,
  and switch the *whole app* to one. The viewer follows the switch, because every endpoint
  asks `ActiveDatabase` where the registry is at request time rather than capturing it at
  start-up.

## Submitting

`POST /api/runs` queues the run. From there the Run inspector owns it — including
`cancel` and `resume`.
