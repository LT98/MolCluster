# Code architecture — the quick reference

`PLAN_implementation.md` is the detail (milestones, changelog, rationale) and
`DESIGN_registry_assembly.md` is the reasoning (decision ledger, subsystems). **This file is
the map.** Read it first; go to the other two when you need to know *why*.

Rule for keeping it useful: if something here needs a paragraph, it belongs in one of the
other two documents with a pointer from here.

---

## 1. The shape of the system, in one paragraph

A **build spec** (serialisable) is turned into **tasks** by `runner.plan`; a worker executes
them via `runner.work`. Executing a task means: enumerate protomers → perceive donor sites →
place ligands around a metal → extract a **typed graph** → hash it into a **composite
identity** → store the structure once, its geometries many times. Energies attach to
geometries, never to structures. Sites attach to structures once (`site_catalog`) and to each
geometry separately (`site_state`). Everything a number came from — code, version, method,
charge, spin — travels with it in a `methods` row.

---

## 2. Module map

| Path | Purpose | State |
|---|---|---|
| `_types.py` | `Fidelity` ladder, `MethodSpec`, exception hierarchy | ✅ |
| `spec.py` | `BuildSpec` — a run is DATA. Versioned + migrated (v5) | ✅ |
| `runner.py` | `plan` writes tasks, `work` executes them. The two never assume one process | ✅ |
| `config.py` | Data root, machine profile, device, ML model — all **declared, never detected** | ✅ |
| `versions.py` | `ALGO_VERSIONS` — pinned recipe versions (ground rule 6) | ✅ |
| **graph/** | | |
| `graph/_types.py` | `TypedGraph`, `NodeLabel`, `EdgeType{COVALENT,DATIVE,METAL_METAL}` | ✅ |
| `graph/canon.py` | Canonical certificate (individualisation-refinement) + WL bucket index | ✅ |
| `graph/from_mol.py` | RDKit mol → typed graph. The molecular input path | ✅ |
| **identity/** | | |
| `identity/keys.py` | L0 composition, L1 certificate hash, `block_id` | ✅ |
| ↳ `l2_isomer_tag` | cis/trans, fac/mer, Δ/Λ | 🔴 **stub → M5** |
| ↳ `l3_conformer_id` | choice-vector label + geometric verifier | 🔴 **stub → M5** |
| **sites/** | | |
| `sites/perception.py` | Which atoms can bind a metal, and what type of donor they are | ✅ |
| `sites/frames.py` | **A site is a FRAME, not a vector** (D13). `live_dof`, `binding_modes`, `torsion_wells` | ✅ |
| `sites/model.py` | `Site`, `perceive`, `chelate_pockets`, `shifting_pocket_donors` | ✅ |
| `sites/state.py` | Per-geometry `SiteStatus` + occlusion. **Never perceives** | ✅ |
| `sites/inherit.py` | Remap sites through an atom map; frames survive untouched | ✅ *(M5 consumes)* |
| `sites/protomers.py` | Protonation as an enumerated branch, not a default | ✅ |
| **descriptors/** | | |
| `descriptors/tables.py` | Curated donor TSV + generated metal TSV (C8) | ✅ |
| `descriptors/ease.py` | Activation-ease floor (**D18**). Components ▸ scalar | ✅ |
| ↳ `hsab_match` | partner term | 🔴 **stub → C7 open** |
| **geometry/** | | |
| `geometry/placer.py` | `place_mononuclear` — fills ONE coordination sphere in one shot | ✅ |
| ↳ `place_multicentre` | inter-centre constraints | 🔴 **stub → M6** |
| `geometry/qc.py` | Clash + distance checks; structured report | ✅ |
| `geometry/distances.py` | M–L target distance as a property of the *pair* | ✅ |
| `geometry/embed.py` | ETKDG + MMFF | ✅ |
| **energy/** | | |
| `energy/backends.py` | xTB / MACE-MP-0 / MACE-OMOL-0 / Null behind one protocol | ✅ |
| `energy/relax.py` | `relax_geometry`, `single_point`, `mode_status` | ✅ |
| `energy/reference.py` | **Refuses bad subtractions** (D17). Balance + isodesmic quality | ✅ |
| **assembly/** | | |
| `assembly/join.py` | `BuildingBlock` + `open_sites` ✅; `compatible`/`join`/`grow` | 🔴 **→ M5/M6** |
| `assembly/choice.py` | `ChoiceVector`, digest, replay | 🔴 **not written → M5** |
| `assembly/construct.py` | deterministic construct + branch-tree enumerator | 🔴 **not written → M5** |
| **registry/** | | |
| `registry/api.py` | **The only write surface** (ground rule 1) | ✅ |
| `registry/db.py` | Connection + additive migration (schema file is not a migration) | ✅ |
| `registry/store.py` | Content-addressed blob store for .xyz | ✅ |
| `registry/jobs.py` | Runs/tasks queue, claim/complete/cancel/resume | ✅ |
| `registry/verify.py` | Does every stored row still agree with its recipe? | ✅ |
| `naming.py` | Labels are **derived from retrieved rows**, never an input to retrieval | ✅ |
| `ui/` | Read-only viewer + spec builder (FastAPI + 3Dmol.js) | ✅ |
| `pathways/` | reaction DAG, path scoring | 🔴 **empty → M8** |

---

## 3. Data model, condensed

```
structures      one row per IDENTITY (L0/L1/L2).  UNIQUE(l0, l1, l2_isomer_tag)
  geometries    many per structure — the fidelity ladder.  best_geometry_id pointer
  site_catalog  one per (structure, canonical_idx).  GEOMETRY-INDEPENDENT, written once
  site_state    one per (site, geometry).  status + ease.  Recomputed per geometry
  reactions     provenance edges INTO a product.  atom_map_json, choice_vector_digest
methods         what produced a number.  Referenced by every stored value
runs / tasks    the queue.  plan() writes, work() drains
```

**Identity on the node, sequence on the edges** (D2). One structure reachable two ways is
one row with two incoming `reactions` edges — 441 such products exist today.

---

## 4. Invariants that actually bite

1. **One write surface.** All inserts go through `registry/api.py`. Nothing else opens a
   cursor to write.
2. **No bare floats.** Every stored number carries a `MethodSpec` + `Fidelity`.
3. **Declared, never detected.** Device (`MOFSBU_DEVICE`), workers (`MOFSBU_PROFILE`/
   `MOFSBU_WORKERS`), ML model (`MOFSBU_ML_MODEL`). A build never seizes hardware on its own.
4. **Ambiguity branches, it does not default.** CN, protonation, spin — enumerate the
   alternatives as separate calls.
5. **A stub never returns a plausible value.** `NotBuiltYet` = missing body;
   `EnergyBackendUnavailable` = missing install. Never collapse the two.
6. **Pinned recipe versions.** Bump `ALGO_VERSIONS`; never silently re-label existing rows.
7. **Fidelity is a property of a GEOMETRY**, never of a structure. Exception:
   `Fidelity.HEURISTIC = -1` is an *ease* rung and `put_geometry` refuses it.
8. **Perceive once per structure.** `site_catalog` is written once and kept;
   `refresh_state` never perceives (there is a test counting calls).
9. **Absent ≠ zero.** An uncomputed ease component is omitted, `n_open_sites` is NULL not 0,
   and a missing number never renders as a low one.

---

## 5. Decision ledger — one line each

Full text in `DESIGN_registry_assembly.md` §7.

| | |
|---|---|
| **D1** | Central object = recursive `BuildingBlock`; assembly re-exposes open sites |
| **D2** | Identity ≠ address ≠ provenance. Identity on node, sequence on edges |
| **D3/D16** | L1 = sha256 of a canonical **certificate**, not the WL hash (1-WL cannot separate µ2-bridging from chelating) |
| **D4** | Fidelity/coords live in child `geometries`; structure identity is fidelity-invariant |
| **D5** | Sites by **canonical index**; `site_catalog` (geometry-free) + `site_state` (per-geometry); inherit via atom map |
| **D6** | Ease = **named components** + derived scalar, tagged by fidelity |
| **D7/D12** | Polynuclear-native from v1 via the typed graph; mononuclear = N=1 |
| **D8** | Pathways = reaction DAG; viability = path-level score |
| **D9** | SQLite registry + content-addressed blob store |
| **D10** | L2 is in scope; discriminator is downstream **relevance**, not ΔE |
| **D11** | L3 = provenance-primary (choice-vector), geometry-verifier |
| **D13** | A site is **frame + live-DOF tag + binding-mode set**; torsion stored as a discrete well index |
| **D14** | Bridging is **derived**, not an edge type |
| **D15** | Net charge is **graph-level**; bond order excluded from the hash (resonance/Kekulé) |
| **D17** | An energy difference needs an **isodesmic** equation, not merely a balanced one |
| **D18** | Ease floor is zero-QM; **absent components stay absent**; `provisional` = "the table value is the wrong question" |

**Open checkpoints:** C2 (θ_geom + energy window — M5), C6 (barrier proxy — M8),
C7 (partner dependence — M8). *Resolved: C1→D10, C4→D12, C5→D18, C8→curated tables.*

---

## 6. Known seams (real tensions, not bugs)

- ~~**`site_catalog` is not a pure function of identity.**~~ **Closed.** D15 excludes bond
  order from the hash and perception used to read it, so two routes to one acetate complex
  perceived different donor sets. `sites.perception.DELOCALISED_GROUPS` now types an
  oxo-acid as a whole group — all its oxygens, one donor type, one charge, no bond order
  read. `registry.catalog_drift` stays as the guard (`put_sites` keeps the FIRST catalog,
  so a regression here is silent otherwise); `test_no_build_route_drifts_from_the_stored_catalog`
  is the gate.
- **Perception counts a metal as an ordinary heavy neighbour.** A coordinated aqua oxygen
  is therefore perceived as no donor at all, so its `SiteStatus.OCCUPIED` row is never
  written. Separate from the resonance seam above and not fixed with it: closing it moves
  `n_perceived_donors` for every assembled structure.
- **L2 is `''` everywhere**, so cis and trans currently collapse into one `structures` row.
  Filling it in retroactively **splits identities** — decide backfill-vs-version before M5
  touches `l2_isomer_tag`.
- **569 existing structures have no `site_state`** (built before M4's second half).
  `n_open_sites` is NULL for them — correct, but a query that treats NULL as 0 will lie.

---

## 7. Where to start reading, by task

| If you're touching… | Read first |
|---|---|
| identity / hashing | `graph/canon.py`, then DESIGN §4 |
| a new donor type | `sites/perception.py` + `data/reference/donor_descriptors.tsv` (a type with no row fails a test) |
| energies | `energy/reference.py` docstring — it explains what it refuses and why |
| the build pipeline | `runner.execute`, top to bottom |
| anything stored | `registry/api.py` — it is the only writer |
| the web UI | `docs/UI_BACKLOG.md` — the known annoyances are already diagnosed |
