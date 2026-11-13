# energy — overview

Read this, `docs/CODE_ARCHITECTURE.md` §4, and the module you are touching. The full plans
(`WORKPLAN_energy.md`, `WORKPLAN_solvation.md`) are only needed when changing *what* gets
computed, not *how* this package serves it.

## What it owns

Two halves, kept apart on purpose:

- **Producing numbers** (`backends`, `dft`, `relax`): run a backend, return an
  `EnergyResult`/`RelaxResult` carrying the `MethodSpec` that made it.
- **Deciding which numbers may be subtracted** (`reference`, `routes`): reaction-balanced
  energies that refuse rather than approximate.
- Also: `protons` (protomer couples as priced edges), `solvation` (continuum corrections
  stored beside an energy, C17).

## Connection points

Changing any name with a *Used by* entry below is a cross-package change: check those
callers first. Everything without one is internal or used only by tests.

## Rules this package must keep

1. **Every number carries its `MethodSpec` + `Fidelity`** — charge and multiplicity are
   part of the spec, never an environment default (invariant 2).
2. **`reference` refuses, with a reason, any equation that is:** unbalanced; balanced but
   not isodesmic; gas-phase when not usable in gas phase (unless `strict=False`); missing an
   energy for any term in the requested medium; scored by the null backend (unless
   `allow_null`); charged species scored by a charge-blind backend; or mixing levels of
   theory. Never drop a term or coerce to make a number appear.
3. **Refusals are data in `routes`.** `price_reaction` returns the diagnosis per edge; it
   does not raise on errors in `REPORTABLE`.
4. **Declared, never detected** (invariant 3). No GPU → DFT is unavailable, not a CPU
   fallback. Device and ML model come from `MOFSBU_DEVICE` / `MOFSBU_ML_MODEL`.
5. **`NotBuiltYet` ≠ `EnergyBackendUnavailable`** (invariant 5): missing body vs. missing
   install.
6. **Reported energies are the stage-5 recipe; MACE is a screen/accelerator.** Don't let a
   change make an ML number reportable where the plan says it screens.
7. **Bump `ALGO_VERSIONS`** for any recipe change; never relabel stored rows (invariant 6).

## Dependencies

- Imports: `_types`, `versions`, `config`, `graph`, `identity`, `geometry` (lazy, in
  `protons`), and `registry`, which it uses for reads and for writes through `registry/api.py` only
  (invariant 1).
- **Must not import `assembly`.** `assembly` imports `energy` for spin helpers, so the reverse
  edge makes a cycle. Shared exceptions live in `_types`.

## Tests

`tests/test_energy.py`, `test_reference.py`, `test_reference_cases.py`, `test_solvation.py`,
`test_relax_pipeline.py`, `test_path_energy.py`. Run `test_pathways.py` and `test_ui.py` as
well when the public surface changes.

<!-- BEGIN GENERATED: scripts/gen_overview.py — do not edit by hand -->

## API (generated)

`*` = re-exported from the package root. *Used by* lists importers in `src/` and `scripts/` outside this package; blank = internal or tests only.

### `backends` (811 lines)

Energy backends behind one protocol (M7).

| Name | Does | Used by |
|---|---|---|
| `class EnergyResult`* | One energy, with the thing that produced it attached. Ground rule 3. |  |
| `class RelaxResult`* | A relaxation: the moved coordinates AND the energy they belong to. |  |
| `class EnergyBackend`* | One interface for xTB, an MLIP, and the test double. |  |
| `d_electrons(symbol, charge) -> int`* | d-count for a transition-metal ion. Raises for anything else. | `scripts/build_metal_descriptors.py` |
| `high_spin_multiplicity(symbol, charge) -> int`* | Maximum-multiplicity ground state for a d-block ion: 2S+1 with S = unpaired/2. | `scripts/regress_m7.py` |
| `minimal_multiplicity(symbols, charge) -> int`* | Singlet if the electron count is even, doublet if odd. |  |
| `spin_class_multiplicity(symbol, charge, spin_class) -> int`* | Route a stated `spin_class` to the ion's own multiplicity. | `assembly.construct`, `runner`, `scripts/build.py` |
| `combined_multiplicity(*multiplicities) -> int`* | Total multiplicity of independent centres: unpaired electrons add, multiplicities don't. | `assembly.join`, `runner`, `scripts/build.py` |
| `electron_count(symbols, charge) -> int` | Total electrons: the sum of atomic numbers minus the net charge. |  |
| `check_spin(symbols, charge, multiplicity) -> None`* | A multiplicity has to be reachable from the electron count. |  |
| `class NullBackend`* | A deterministic, deliberately non-physical backend for tests. |  |
| `class XTBBackend`* | GFN2-xTB through tblite's ASE calculator. | `scripts/regress_m7.py` |
| `class MACEBackend`* | MACE-MP-0, a universal MLIP — fast, and blind to charge and spin. |  |
| `class MACEOmolBackend`* | MACE-OMOL-0 — the same architecture, trained on OMol25, and NOT charge-blind. |  |
| `class MACEMHBackend` | MACE-MH-1, a multi-head foundation model; one head per level of theory. |  |
| `class MACEPolarBackend` | MACE-POLAR-1: explicit long-range electrostatics, given total charge and spin. |  |
| `ml_backend_key(ml_model) -> str`* | Which ML backend a request means: the one named, else the one declared. |  |
| `get_backend(name, **kwargs) -> EnergyBackend`* | A backend by its registry key (`xtb`, `mace_omol`, `dft`, ...); unknown keys raise. | `runner`, `scripts/compare_relaxers.py`, `scripts/dft_arbitrate.py`, `scripts/dft_reopt.py`, `scripts/funnel_eval.py` |
| `backend_for(fidelity, *, ml_model, **kwargs) -> EnergyBackend`* | The backend that produces geometries at this rung, or a reason why not. | `runner` |
| `available_backends() -> dict[str, bool]`* | What this machine can actually run. The UI reads this rather than guessing. |  |

### `dft` (363 lines)

DFT through GPU4PySCF — the `Fidelity.DFT` rung (WORKPLAN_energy E4, stage 5b/5c).

| Name | Does | Used by |
|---|---|---|
| `class OxidationStateMismatch` | The SCF put a different number of unpaired electrons on the metal than L0 claims. |  |
| `contraction_engine() -> str` | Which engine GPU4PySCF contracts with here: 'cutensor', or the fallback it warned of. | `scripts/dft_reopt.py` |
| `class DFTBackend` | Kohn–Sham DFT on GPU4PySCF behind the `EnergyBackend` protocol. | `scripts/dft_arbitrate.py`, `scripts/dft_reopt.py` |
| `expected_unpaired(symbol, oxidation_state, spin_class) -> int` | Unpaired d electrons the graph claims for one metal centre. |  |
| `check_oxidation_state(result, metal_index, expected, *, tol) -> float` | Return the metal's spin population, or raise `OxidationStateMismatch`. |  |

### `protons` (182 lines)

Which structures are the same molecule with a proton removed, and what that costs.

| Name | Does | Used by |
|---|---|---|
| `ensure_species(reg, name, smiles, charge, *, target, ml_model) -> int` | Register a small closed-shell species if absent, with an energy at `target`. |  |
| `ensure_couple(reg, *, target, ml_model) -> tuple[int, int]` | `(water_id, hydronium_id)`, each with an energy at `target`. | `runner`, `scripts/link_protomers.py` |
| `labile_hydrogens(g) -> list[int]` | Every H bonded to one non-carbon, non-metal atom: the protons a graph can lose. |  |
| `deprotonation_pairs(reg) -> list[tuple[int, int, int]]` | `(protonated_id, deprotonated_id, 1)` for every pair exactly one proton apart. | `scripts/link_protomers.py` |
| `link_protomers(reg, *, water_id, hydronium_id, limit) -> list[dict[str, Any]]` | Write one balanced deprotonation edge per pair, skipping what already exists. | `runner`, `scripts/link_protomers.py` |

### `reference` (717 lines)

Reaction-balanced energies — the fix for the charged-ion reference problem (M7).

| Name | Does | Used by |
|---|---|---|
| `class Term`* | One species in the equation, with the sign the arrow gives it. | `pathways.route` |
| `class BalanceReport`* | Why an equation is or is not usable. Returned, not raised, so a UI can show it. |  |
| `class ReactionEnergy`* | A reaction energy that knows what it is allowed to claim. |  |
| `reaction_terms(reg, reaction_id) -> tuple[Term, ...]`* | Every species in the reaction, each on the side its role puts it. |  |
| `check_balance(reg, reaction_id) -> BalanceReport`* | Does this equation conserve atoms and charge? | `assembly.derive` |
| `balance_for_terms(reg, terms) -> BalanceReport` | `check_balance` for an equation that is not (yet) a row. | `pathways.route` |
| `class QualityIssue`* |  |  |
| `class ReferenceQuality`* | Whether the equation's two sides are comparable enough to subtract. |  |
| `check_reference_quality(reg, reaction_id, *, max_naked_charge, medium) -> ReferenceQuality`* | Are the reference species comparable to the product, or merely balanced? |  |
| `quality_for_terms(reg, terms, *, max_naked_charge, medium) -> ReferenceQuality` | `check_reference_quality` for an equation that is not (yet) a row. | `pathways.route` |
| `reaction_balanced_energy(reg, reaction, *, fidelity, solvent, allow_null, allow_unconverged, strict) -> ReactionEnergy`* | dE for a balanced reaction, or a refusal that says which rule was broken. |  |
| `energy_for_terms(reg, terms, *, fidelity, solvent, allow_null, allow_unconverged, strict, subject) -> ReactionEnergy` | `reaction_balanced_energy` for an equation that is not (yet) a row. | `pathways.route` |
| `store_reaction_energy(reg, reaction_id, energy, *, force) -> None`* | Record dE on the reaction, with the method and fidelity that produced it. |  |
| `put_balanced_reaction(reg, product_structure_id, *, reagents, leaving, chelators, solvents, kind, note, depth) -> int`* | Store a reaction and REFUSE it if it does not balance. |  |

### `relax` (147 lines)

Geometry optimisation — the caller-facing half of the energy layer (M7).

| Name | Does | Used by |
|---|---|---|
| `relax_geometry(coords, symbols, *, charge, multiplicity, target, solvent, max_steps, fmax, ml_model) -> RelaxResult`* | Optimise a geometry to the requested fidelity, returning coords + energy. | `runner` |
| `single_point(coords, symbols, *, charge, multiplicity, target, solvent, ml_model)`* | One energy, no relaxation. The other half of what the reference scheme needs. |  |
| `available_modes(ml_model) -> dict[str, bool]`* | Which run modes can actually execute, asked of the backends, not asserted. |  |
| `mode_status(ml_model) -> dict[str, dict[str, Any]]`* | `available_modes` plus the reason, for a UI that has to explain a disabled option. | `runner`, `ui.builder` |
| `ml_model_status() -> dict[str, dict[str, Any]]`* | Both ML foundation models, side by side, for the UI's model picker. | `ui.builder` |

### `routes` (252 lines)

What one provenance edge costs, and why it sometimes cannot be said.

| Name | Does | Used by |
|---|---|---|
| `price_reaction(reg, reaction_id, *, min_fidelity, solvent) -> dict[str, Any]`* | dE for one edge, or the stated reason there is none. | `pathways.route`, `ui.app` |
| `price_incoming_routes(reg, structure_id, *, min_fidelity, solvent) -> list[dict[str, Any]]`* | Every edge arriving at this structure, priced or refused, in `reactions` order. |  |
| `decompositions(reg, structure_id, *, min_fidelity, solvent, limit) -> list[dict[str, Any]]` | Pairs of registry structures whose atoms and charge add up to this one. | `pathways.route`, `ui.app` |

### `solvation` (90 lines)

Continuum corrections on stored energies — the producer (WORKPLAN_solvation S2, C17).

| Name | Does | Used by |
|---|---|---|
| `correct(reg, geometry_id, medium) -> Any` | Compute and store the continuum correction for one geometry. Returns the `Put`. |  |
| `correct_registry(reg, medium, *, min_fidelity, commit_every, limit) -> dict[str, Any]` | Correct every converged gas-phase energy at or above `min_fidelity` that lacks one. | `runner`, `scripts/solvate.py` |

<!-- END GENERATED -->
