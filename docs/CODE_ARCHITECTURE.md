# Code architecture — the quick reference

**This file is the map.** Read it first; go elsewhere when you need to know *why*.

| Document | Holds |
|---|---|
| `PLAN_implementation.md` | what is **left to build** — remaining milestones, exit gates, open decision gates |
| `DESIGN_registry_assembly.md` | the **reasoning** — decision ledger, subsystems, open checkpoints |
| `BUGS.md` | what is **built and behaving wrongly**. Open only |
| `archive/` | finished milestones, resolved bugs, and the full changelogs |

Rule for keeping this file useful: if something here needs a paragraph, it belongs in one of
the others with a pointer from here. Rule for keeping the set useful: an item leaves the active
document the moment it is done — it moves to `archive/`, it is not struck through in place.

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
| `spec.py` | `BuildSpec` — a run is DATA. Versioned + migrated (v7). Counts accept ranges (`"1~3"`), expanded at construction | ✅ |
| `runner.py` | `enumerate_plan` enumerates, `plan` writes tasks, `work` executes them. The two never assume one process | ✅ |
| ↳ `execute_run` | the only place the worker pool is built — CLI and page share it. `plan_workers` divides the declared workers between the build queue and the relax queue (spawn, never fork: CUDA) | ✅ |
| ↳ `estimate` | the same enumeration with no registry — what the builder page shows before you submit | ✅ |
| ↳ `grow` tasks | `spec.pathways`: the rung below each product, and the `join` between them | ✅ |
| `config.py` | Data root, machine profile, device, ML model — all **declared, never detected** | ✅ |
| `versions.py` | `ALGO_VERSIONS` — pinned recipe versions (ground rule 6) | ✅ |
| **graph/** | | |
| `graph/_types.py` | `TypedGraph`, `NodeLabel`, `EdgeType{COVALENT,DATIVE,METAL_METAL}` | ✅ |
| `graph/canon.py` | Canonical certificate (individualisation-refinement) + WL bucket index | ✅ |
| `graph/from_mol.py` | RDKit mol → typed graph. The molecular input path | ✅ |
| **identity/** | | |
| `identity/keys.py` | L0 composition, L1 certificate hash, `block_id` | ✅ |
| ↳ `l2_isomer_tag` | cis/trans, fac/mer, Δ/Λ | ✅ *(`identity/isomers.py`; needs a geometry, so `""` without one)* |
| ↳ `l3_conformer_id` | choice-vector label + geometric verifier | ✅ *(`identity/conformers.py`; θ_geom = 0.15 Å calibrated, energy window open)* |
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
| ↳ `place_multicentre` | **reconciliation only** (D20) — where two determinants fix one M···M and disagree. `Center` and `InterCentreConstraint` are defined; `Center.element` is not always a metal (a bridging atom is a centre) | 🔴 **stub → M6 S3** |
| `geometry/qc.py` | Clash + distance checks; structured report | ✅ |
| ↳ `check_intercentre` | validates the M···M the joins produced. `qc()` alone cannot see a squeezed node — 1.90 Å Cu···Cu is above the clash floor and is not a metal–donor pair | 🔴 **stub → M6 S6** |
| `geometry/distances.py` | M–L target distance as a property of the *pair* | ✅ |
| `geometry/embed.py` | ETKDG + MMFF | ✅ |
| `geometry/_linalg.py` | Pure rotation/alignment math, one copy. **Two rotation forms on purpose** — matrix and Rodrigues are not bit-identical and frames were built with the latter | ✅ |
| **energy/** | | |
| `energy/backends.py` | xTB / MACE-MP-0 / MACE-OMOL-0 / Null behind one protocol | ✅ |
| `energy/relax.py` | `relax_geometry`, `single_point`, `mode_status` | ✅ |
| `energy/reference.py` | **Refuses bad subtractions** (D17). Balance + isodesmic quality | ✅ |
| **assembly/** | | |
| `assembly/join.py` | `BuildingBlock`, `open_sites`, `compatible`/`chelate_compatible`, `join`, `join_chelate`/`chelate_reach`, `grow`. A join takes `lone_pair=` — which lobe of an sp2 donor binds is worth 2.8 Å of M···M and has no default in the chemistry | ✅ |
| ↳ `join_bridge` | one ligand across vertices of **different** metals in one move — the nucleus-first route. The sequential route needs no such thing (D20) | 🔴 **stub → M6 S1** |
| `assembly/choice.py` | `ChoiceVector`, canonical form, version-prefixed digest, JSON round trip | ✅ |
| `assembly/construct.py` | deterministic construct + branch-tree enumerator + Kind-C refusals | ✅ |
| `assembly/persist.py` | Stores an assembled block: L2 from its geometry, sites by inheritance, provenance. **Writes nothing itself** — orchestrates `registry.api` | ✅ |
| **registry/** | | |
| `registry/api.py` | **The only write surface** (ground rule 1) | ✅ |
| `registry/db.py` | Connection + additive migration (schema file is not a migration) | ✅ |
| `registry/store.py` | Content-addressed blob store for .xyz | ✅ |
| `registry/jobs.py` | Runs/tasks queue, claim/complete/cancel/resume | ✅ |
| ↳ `run_liveness` | is a process still behind this run? pid (proof, same host) + heartbeat (only when nothing is claimed). `sweep_interrupted` acts on it when a process takes the database over | ✅ |
| `registry/verify.py` | Does every stored row still agree with its recipe? | ✅ |
| `naming.py` | Labels are **derived from retrieved rows**, never an input to retrieval | ✅ |
| `ui/` | Read-only viewer + spec builder + run inspector (FastAPI + 3Dmol.js) | ✅ |
| `ui/static/chrome.{js,css}` | The three pages' shared tab strip and registry picker. `MOFSBU_PAGES` is the one list of pages; a new page is an entry there and a route | ✅ |
| `pathways/` | reaction DAG, path scoring | 🔴 **empty → M8** *(the DAG's edges are written now — `spec.pathways` — but nothing scores them)* |

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
                a run has FOUR endings: done, failed, cancelled (asked for) and
                interrupted (the process stopped existing — derived, never written
                by the process it happened to)
```

**Identity on the node, sequence on the edges** (D2). One structure reachable two ways is
one row with two incoming `reactions` edges — 441 such products exist today.

---

## 4. Invariants that actually bite

1. **One write surface.** All inserts go through `registry/api.py`. Nothing else opens a
   cursor to write.
2. **No bare floats.** Every stored number carries a `MethodSpec` + `Fidelity`.
3. **Declared, never detected.** Device (`MOFSBU_DEVICE`), workers (`MOFSBU_PROFILE`/
   `MOFSBU_WORKERS`, split by `MOFSBU_RELAX_WORKERS`), ML model (`MOFSBU_ML_MODEL`). A build
   never seizes hardware on its own — but a declared count is *used*, including one passed
   as `workers=`, and it is divided between the two queues rather than exceeded.
4. **Ambiguity branches, it does not default.** CN, protonation, spin — enumerate the
   alternatives as separate calls.
5. **A stub never returns a plausible value.** `NotBuiltYet` = missing body;
   `EnergyBackendUnavailable` = missing install. Never collapse the two.
6. **Pinned recipe versions.** Bump `ALGO_VERSIONS`; never silently re-label existing rows.
7. **Fidelity is a property of a GEOMETRY**, never of a structure. Exception:
   `Fidelity.HEURISTIC = -1` is an *ease* rung and `put_geometry` refuses it.
8. **Perceive once per structure.** `site_catalog` is written once and kept *within a
   perception version*; a catalog older than `ALGO_VERSIONS["perception"]` is rewritten the
   next time anything touches the structure, and its `site_state` goes with it (D19).
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
| **D3** | L1 = canonical hash of an explicit typed graph; never SMILES/InChI. *Producer superseded by D16* |
| **D4** | Fidelity/coords live in child `geometries`; structure identity is fidelity-invariant |
| **D5** | Sites by **canonical index**; `site_catalog` (geometry-free) + `site_state` (per-geometry); inherit via atom map |
| **D6** | Ease = **named components** + derived scalar, tagged by fidelity |
| **D7/D12** | Polynuclear-native from v1 via the typed graph; mononuclear = N=1 |
| **D8** | Pathways = reaction DAG; viability = path-level score |
| **D9** | SQLite registry + content-addressed blob store |
| **D10** | L2 is in scope; discriminator is downstream **relevance**, not ΔE |
| **D11** | L3 = provenance-primary (choice-vector), geometry-verifier |
| **D13** | A site is **frame + live-DOF tag + binding-mode set**; torsion stored as a discrete well index |
| **D14** | Bridging is **derived**, not an edge type (`EdgeType` is `{COVALENT, DATIVE, METAL_METAL}`) |
| **D15** | Net charge is **graph-level**; bond order excluded from the hash (resonance/Kekulé) |
| **D16** | L1 = sha256 of a canonical **certificate**, not the WL hash (1-WL cannot separate µ2-bridging from chelating). WL is a bucket index; nauty is unused |
| **D17** | An energy difference needs an **isodesmic** equation, not merely a balanced one |
| **D18** | Ease floor is zero-QM; **absent components stay absent**; `provisional` = "the table value is the wrong question" |
| **D19** | Re-derive what is only an **annotation** (`site_catalog`); **version** what is an address (identity). A stored identity keeps the answer its own recipe version gave |

**Open checkpoints:** C2 **half-resolved** (θ_geom = 0.15 Å, calibrated on xTB-relaxed
geometries; the energy window stays open — see ISSUES 6b), C6 (barrier proxy — M8),
C7 (partner dependence — M8). *Resolved: C1→D10, C4→D12, C5→D18, C8→curated tables.*

---

## 6. Known seams

Two places where the code is internally consistent and still reports something misleading.
Full diagnosis and current measurements in `BUGS.md`; one line each here.

| | Seam | Bites when |
|---|---|---|
| [B2](BUGS.md#b2) | `l2_isomer_tag` is `''` on every stored row, so cis and trans are one `structures` row | the classifier exists and `store_block` derives a tag where it has coordinates, but the run pipeline passes `l2=""` **on purpose** (D22) so its products land on the node every other route reaches. Both paths tag together or neither |
| [B3](BUGS.md#b3) | Structures built before M4's second half have no `site_state`, so `n_open_sites` is NULL | a query treats NULL as 0 and reports a fully-occupied structure. Self-heals as `perception/2` re-derivation reaches each one (D19) |

*Closed, and both worth reading before touching perception: `site_catalog` was not a pure
function of identity, and a coordinated donor was perceived as no donor at all —
`archive/BUGS_resolved.md`.*

---

## 7. Where to start reading, by task

| If you're touching… | Read first |
|---|---|
| identity / hashing | `graph/canon.py`, then DESIGN §4 |
| a new donor type | `sites/perception.py` + `data/reference/donor_descriptors.tsv` (a type with no row fails a test) |
| energies | `energy/reference.py` docstring — it explains what it refuses and why |
| the build pipeline | `runner.execute`, top to bottom |
| anything stored | `registry/api.py` — it is the only writer |
| the web UI | `docs/BUGS.md` — B4/B5/B6 are the open ones, already diagnosed |
