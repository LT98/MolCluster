# Implementation plan — what is left to build

> **Looking for the map, not the detail?** `CODE_ARCHITECTURE.md` is the quick reference:
> module table, invariants, decision ledger one-liners. This document is the detail behind
> it — remaining milestones, exit gates, open decision gates.

**This file holds only unfinished work.** M0–M4 are done; they and the full changelog live in
[`archive/PLAN_completed.md`](archive/PLAN_completed.md). Things that are *built and behaving
wrongly* are not milestones — they are in [`BUGS.md`](BUGS.md).

**Companion to** [`DESIGN_registry_assembly.md`](DESIGN_registry_assembly.md). That doc says
*what* is being built and *why* (ledger D1–D18). This doc says *what code exists*, *in what
order*, and *what has to be true before the next thing starts*.

**Form:** dependency-ordered milestones, no calendar dates. Each has the modules it lands, an
exit gate (a test that must pass, not a feeling), and an effort size (S ≈ a sitting or two,
M ≈ a week of focused evenings, L ≈ multi-week / headline cost).

## Where the work stands

| | |
|---|---|
| **Done** | M0 rails · M1 descriptors · M2 graph + identity · M3 registry · M3.5 viewer · M4 sites · **M5 assembly** |
| **Partly done** | **M7** — backends, `MethodSpec`, reference scheme and the relax runner all shipped; the regression exit gate has a harness and no recorded result |
| **Next** | **M6** polynuclear nodes — the first milestone that builds a real SBU. `WORKPLAN_M6.md` carries the slice-level detail |
| **Then** | M8 pathways (the actual contribution) · M9 folded in opportunistically |

```
[M0 · M1 · M2 · M3 · M3.5 · M4 · M5 — done] ─► M6 polynuclear ─┐
                                                               ├─► M8 pathways
                    M7 energy — built, exit gate never run ────┘
```

**Open decision gates:** C2's energy window (blocked on [B9](BUGS.md#b9)), C9 and C10 (forced by
M6), C6 and C7 (forced by M8). See §3. Nothing is pre-resolved here; the design doc's leanings
stand.

---

## 0. Ground rules (invariants — true from the first commit, never retrofitted)

These are cheap to hold now and expensive to add later. Every one of them is a test in M0–M3.

1. **One insert path.** Nothing writes to the registry except `registry.api.put_structure` /
   `put_geometry` / `put_reaction`. No script writes a `.xyz` as a primary record. Files are a
   *view*, exported from the registry.
2. **Identity comes only from the typed graph.** No filename, no build order, no SMILES, no
   dict key ever carries meaning. If you can't recover it from `structures` + `geometries`, it
   doesn't exist.
3. **Every number carries `(fidelity, method, algo_version)`.** No bare float is stored anywhere.
   A value without provenance is a bug, not a shortcut.
4. **No `if n_metals == 1` branch, ever** (D12). Mononuclear is N=1 through the general code path.
   A reviewer grepping for `n_metal`/`single_metal` should find nothing but tests.
5. **`construct` is deterministic** given `(choice_vector, seed)` and **emits its choice-vector**
   (D13). Any ambiguity (CN, geometry, protonation, spin) *branches explicitly* — a silent default
   is a bug. There is a test that asserts an ambiguous spec raises rather than defaults.
6. **Algorithm versions live in one file** (`mofsbu/versions.py`). Changing a hash recipe bumps a
   version and marks affected rows stale; it never silently re-labels old records.
7. **A schema file is not a migration.**  `CREATE TABLE IF NOT EXISTS` creates what is
   absent and does NOTHING to a table that already exists, so a column added later is
   silently missing from every database written before it.  `Registry.migrate` therefore
   reconciles: it builds the target schema in memory, adds columns this database lacks,
   and recreates views (a stale `CREATE VIEW IF NOT EXISTS` is never updated).  Additive
   changes are applied; anything destructive is recorded and left alone.  The same rule
   applies one level up to spec files, which are versioned with a migration.
8. **Missing functions get a STUB, never a downgrade.**  When a step needs something that
   does not exist yet, define the call with its real signature and have it raise — do not
   reshape the pipeline around what happens to be implementable today.  The interface is
   the design decision; the body is schedulable work.  Reshaping instead is how the
   protomer enumeration became "deprotonate one, then two" and how pocket selection became
   "take the first one".
   **A stub must fail loudly.**  `raise NotImplementedError("... — needs M5")`, or return a
   value that is visibly a placeholder (`?L2`).  It must never return a plausible default:
   a stub that quietly answers is indistinguishable from a working function until the
   answers are wrong, which is exactly how a registry full of impossible molecules got
   written and believed.
9. **Parallelism is opt-in, and the laptop is the default.**  Nothing spawns worker
   processes unless the machine is explicitly declared a workstation
   (`MOFSBU_PROFILE=workstation`, or `MOFSBU_WORKERS=N`).  The default is one in-process
   worker, so running a build on the laptop cannot swamp it.  Scaling up is a
   configuration change, never a code change.
10. **Two-machine parity.** The L1 key is computed by a path available on *both* machines — the
   canonical certificate is pure Python and depends on no optional native package (D16).
   Nothing that is installed on only one machine may ever produce a key, or the laptop and the
   workstation disagree about identity. The golden-hash test in M2 is what enforces this.

---

## 1. Package map and API surface

This section was originally a contract to code against, written before the code existed. Most
of it now exists, and where the built thing differs from the original sketch, **the code is the
contract** — the sketch is not a target to migrate toward. So §1 is split:

* **§1a** is a map of what is built, verified against the modules. Signatures are exact as of
  the last audit; treat the module as authoritative if they ever disagree again.
* **§1b** is design intent for what is *not* built. Those signatures are still proposals.

Types are indicative throughout (dataclasses, `from __future__ import annotations` everywhere).

---

### 1a. Built — the seams M5–M8 will call

#### `mofsbu/graph/` — the typed molecular graph (identity substrate)

```python
# _types.py
class EdgeType(str, Enum):
    COVALENT = "cov"; DATIVE = "dat"; METAL_METAL = "mm"
class BridgeClass(str, Enum):                 # DERIVED, never stored as an edge (D14)
    NONE; TERMINAL; MU2; MU3; MU_N

@dataclass(frozen=True)
class NodeLabel:
    element: str; formal_charge: int = 0
    oxidation_state: int | None = None        # per-center, D12 — set on metals
    spin_class: str | None = None             # per-center, D12 — e.g. "hs", "ls"
    role: str = "atom"

class TypedGraph:
    def add_atom(self, ...) -> int
    def add_bond(self, i: int, j: int, etype: EdgeType, *, order: float = 1.0) -> None
    def neighbors(self, i, etype: EdgeType | None = None) -> list[int]
    def metals(self) -> list[int]
    def bridging_metals(self, i: int) -> tuple[int, ...]   # µ-ness is COUNTED, not tagged
    def bridge_class(self, i: int) -> BridgeClass
    def ligand_fragments(self) -> list[tuple[int, ...]]
    def net_charge(self) -> int                            # graph-level (D15)
    def localised_charge(self) -> int                      # the per-atom sum, separately
    def subgraph(self, idxs) -> TypedGraph
    def relabel(self, mapping: dict[int, int]) -> TypedGraph
    def permuted(self, seed: int) -> TypedGraph            # the build-order-invariance test
    def to_networkx(self) -> nx.Graph
    def to_json(self) -> bytes / from_json(blob) -> TypedGraph
```

A µ2-carboxylate is **not** a duplicated fragment: its atoms appear once, and its two O atoms
each carry a `DATIVE` edge to a different metal. Nothing tags the pair — `bridge_class` counts
the distinct metals reached (D14).

```python
# canon.py — D16.  Pure Python: no optional native package may produce a key (ground rule 10)
def certificate_digest(g, order=None) -> str      # PRIMARY key producer
def canonical_order(g, *, max_leaves=200_000) -> list[int]   # individualisation-refinement
def canonical_certificate(g) -> Certificate
def canonical_index_map(g) -> dict[int, int]
def wl_hash(g, iterations=3) -> str               # fast BUCKET INDEX, not a key
def automorphism_generators(g, ...) / automorphism_group(g, ...) / atom_orbits(g)
def is_isomorphic(a, b) -> bool                   # certificate compare
def vf2_isomorphic(a, b) -> bool                  # collision resolver ONLY

# from_mol.py — the molecular input path
def from_rdkit(mol, ...) -> TypedGraph
def from_smiles(...) -> TypedGraph
def mol_from_smiles(smiles, *, embed=False, seed=0xC0FFEE) -> Chem.Mol
```

There is **no `nauty_canon`**, and `pynauty` is imported nowhere — dropped with D16.

#### `mofsbu/identity/` — the composite L0–L3 key

```python
def l0_composition(g) -> str            # "Cu2_C8H20O8_q0_s3|Cu:+2/hs,Cu:+2/hs"
def l1_graph_hash(g) -> str             # = certificate_digest(g); version ALGO_VERSIONS['l1_certificate']
def wl_index(g) -> str                  # the bucket, stored in its own column
def block_id(l0, l1, l2=L2_UNSET, l3=L3_UNSET) -> str
def identity(g, *, geom=None) -> dict[str, str]
def hill_formula(counts: dict[str, int]) -> str
```

`l2_isomer_tag` / `l3_conformer_id` exist with their final signatures and return the unset
placeholder — see §1b.

#### `mofsbu/registry/` — SQLite + content-addressed blobs

```python
class BlobStore:
    def __init__(self, root: str | Path)
    def put(self, data: bytes) -> str            # sha256 hex; store/ab/cdef….bin
    def put_text(self, text: str) -> str
    def get(self, digest) -> bytes / get_text(digest) -> str
    def path(self, digest) -> Path / has(digest) -> bool

class Registry:                                   # context manager; owns conn + migrations
    def __init__(self, db_path, store: BlobStore | None = None, ...)
    def migrate(self, note: str = "") -> None     # additive reconcile; a schema file is not a migration
    def record_algo_versions(self) / stale_algo_versions() -> dict[str, tuple[str, str]]
    def count(self, table: str) -> int
```

`api.py` is the **only write surface** (ground rule 1). It returns a `Put` NamedTuple, not a
bare id, because "was this newly inserted or recognised under D2" is the answer most callers
need:

```python
def put_structure(reg, g, *, l2=None, provenance: Provenance | None = None, tags=()) -> Put
def put_geometry(reg, structure_id, xyz_text: str, *, fidelity: Fidelity, method: MethodSpec,
                 energy=None, converged=None, relaxed_from=None,
                 choice_vector=None, seed=None, qc=None) -> Put
def put_reaction(reg, product_id: int, prov: Provenance) -> int   # ONE provenance object
def put_sites(reg, structure_id, sites, *, algo="perception/1", ...) # written ONCE per structure
def put_site_state(reg, structure_id, geometry_id, states, ...)     # per geometry
def get_structure(reg, ref: int | str) -> sqlite3.Row               # id or block_id
def get_graph(reg, structure_id) -> TypedGraph
def get_sites(reg, sid) / get_site_state(reg, sid, ...) / canonical_map(reg, sid)
def best_geometry(reg, structure_id, min_fidelity=None) -> sqlite3.Row | None
def catalog_drift(reg, structure_id, sites) -> list[str]            # the D15 guard
def set_hidden(reg, structure_id, hidden=True, ...)                 # hide IS the delete
def incoming_routes(reg, structure_id) -> list[dict]
def display_label(g, aliases=None) -> str                           # derived, never a lookup key
def export_xyz(reg, geometry_id, path) -> Path                      # files are a VIEW
def find(reg, *, l0=None, l1=None, formula_like=None, metal=None, n_metals=None, charge=None,
         max_bridge_class=None, has_metal_metal=None, min_fidelity=None, has_route=None,
         tag=None, fragment_l1=None, fragment_formula=None,
         sort="id", order="asc", limit=100, offset=0) -> list[sqlite3.Row]
```

`put_structure` is **idempotent** on `(L0, L1, L2)`: the same key returns the existing id and
adds a provenance edge. That single property is what makes "one node, two build routes" real.

Rows come back as `sqlite3.Row`, not as typed record objects — there is no `StructureRecord` or
`GeometryRecord`.

#### `mofsbu/sites/` — perceive once, refresh accessibility

```python
# perception.py — takes an RDKit mol, NOT a graph
@dataclass(frozen=True)
class DonorSite: ...
def find_donor_sites(mol) -> list[DonorSite]
def deprotonate(mol, sites) -> tuple[Chem.Mol, list[DonorSite]]     # (was sketched as `activate`)
LABILE_DONOR_PATTERNS / DELOCALISED_GROUPS                          # group-level typing (D15)

# frames.py — D13
class LiveDOF(str, Enum): TORSION_LIVE; TORSION_FREE
class BindingMode(str, Enum): MONODENTATE; CHELATE; BRIDGE_*; ...
@dataclass(frozen=True)
class SiteFrame: ...
def site_frame(mol, donor_idx, donor_type, conf=None, ...) -> SiteFrame
def live_dof(donor_type) -> LiveDOF
def binding_modes(donor_type) -> tuple[BindingMode, ...]
def torsion_wells(donor_type, mode=BindingMode.MONODENTATE) -> tuple[float, ...]

# model.py — the pocket layer M5 needs
@dataclass
class Site: ...                                   # DONOR or VACANCY
def perceive(mol, *, with_frames=True) -> list[Site]
def vacancy_sites(metal_idx, origin, directions) -> list[Site]
@dataclass
class Pocket: ...
def chelate_pockets(...) / find_pockets(...) / find_pocket(mol, sites, **predicate)
def shifting_pocket_donors(mol, sites) -> frozenset[int]

# state.py — NEVER perceives (there is a test counting calls)
class SiteStatus(str, Enum): OPEN; OCCUPIED; BLOCKED
@dataclass
class SiteState: ...
def refresh_state(...) -> list[SiteState]
def buried_volume(coords, site, *, radius=3.5, ...) -> float

# inherit.py — M5 consumes this
def inherit_sites(parent_sites, atom_map: AtomMap, ...) -> list[Site]
def merge_inherited(*groups) -> list[Site]

# protomers.py — protonation as an enumerated branch, not a default
@dataclass
class Protomer: ...
def enumerate_protomers(...) / labile_sites(mol, donor_type=None) / configuration_key(...)
```

#### `mofsbu/descriptors/` — the shared lookup substrate

TSV-backed (C8 resolved), loaded through functions rather than exposed as module-level dicts,
so a missing row raises `UnknownDescriptor` at the point of use:

```python
@dataclass(frozen=True)
class DonorDescriptor: donor_type; pka; pka_sigma; hsab; ... ; source; source_version
@dataclass(frozen=True)
class MetalDescriptor: symbol; charge; ionic_radius; hsab; preferred_cn; d_electrons; ...
def donor_table(path=None) -> dict[str, DonorDescriptor]        # data/reference/donor_descriptors.tsv
def metal_table(path=None) -> dict[tuple[str, int], MetalDescriptor]
def donor(donor_type) -> DonorDescriptor                        # raises UnknownDescriptor
def metal(symbol, charge) -> MetalDescriptor
def known_donor_types() / known_ions() / sync_to_registry(reg, ...)

# ease.py — D18
@dataclass(frozen=True)
class EaseRecord: components: dict[str, float]; scalar; fidelity; confidence; provisional; method
def activation_ease(site, *, partner=None, conditions=None, state=None) -> EaseRecord
def deprotonation_energy(backend, symbols, positions, ...) -> float    # the charge-aware rung
def hsab_match(donor, metal) -> float                                  # RAISES — C7 open
```

#### `mofsbu/geometry/` — local geometries, the placer, QC

There is **no `local.py`** and no `templates.py`; `site_vectors` lives in `placer.py`.

```python
# placer.py
def site_vectors(geometry: str, n: int, d: float = 2.05) -> np.ndarray
GEOMETRIES: dict[str, dict[int, str]]             # CN -> local geometry; CN=5 fixed
@dataclass
class LigandPlacement: ...                        # incl. azimuth_step / oop_step = the replay coords
@dataclass
class PlacementResult: ...
def place_mononuclear(...) -> PlacementResult     # fills ONE coordination sphere in one shot
def to_rdkit(metal, ligands, result, ...)

# distances.py — M–L target distance is a property of the PAIR
def metal_donor_distance(metal, donor_element, ...) -> Distance
def base_distance(metal) -> float

# qc.py
def check_clashes(symbols, coords, bonded, ...) -> list[Clash]
def check_metal_bonds(coords, metal_idx, donor_idxs, ...) -> list[BadBond]
def qc(symbols, coords, bonded, ...) -> QCReport

# embed.py
def embed_molecule(mol, *, seed=0xC0FFEE, relax=True) -> Chem.Mol
def embed_with_report(...) / coordinates(mol) / set_coordinates(mol, coords) / to_xyz(mol, comment="")
```

#### `mofsbu/energy/` — backends behind one protocol

```python
class Fidelity(IntEnum):
    HEURISTIC = -1      # ease only; put_geometry REFUSES it (D18)
    RAW = 0; FF = 1; ML = 2; XTB = 3; DFT = 4

@dataclass(frozen=True)
class MethodSpec:
    code; version; method; solvent; charge; multiplicity; extras

class EnergyBackend(Protocol):
    def single_point(...) -> EnergyResult
    def relax(...) -> RelaxResult
class NullBackend / XTBBackend / MACEBackend (MP-0) / MACEOmolBackend (OMOL-0)
#   charge_aware / spin_aware are CLASS ATTRIBUTES — "does this model see charge?" in one place

def backend_for(fidelity, *, ml_model=None, ...) -> EnergyBackend
def get_backend(name, **kw) / available_backends() -> dict[str, bool]
def high_spin_multiplicity(symbol, charge) -> int
def spin_class_multiplicity(...) / combined_multiplicity(...) / check_spin(...)

# relax.py
def relax_geometry(coords, symbols, *, charge, multiplicity, ...) -> ...
def single_point(...) / available_modes(ml_model=None) / mode_status(ml_model=None)

# reference.py — REFUSES bad subtractions (D17)
def reaction_terms(reg, reaction_id) -> tuple[Term, ...]
def check_balance(reg, reaction_id) -> BalanceReport
def check_reference_quality(reg, reaction_id, ...) -> ReferenceQuality
def reaction_balanced_energy(reg, ...) -> ReactionEnergy      # strict=False to reproduce legacy
def store_reaction_energy(...) / put_balanced_reaction(...)
```

#### The run layer, and `scripts/`

`spec.py` (`BuildSpec`, versioned + migrated, currently v5) is the input; `runner.plan` writes
tasks and `runner.work` drains them; `registry/jobs.py` is the queue. `naming.py` derives labels
from retrieved rows and is never an input to retrieval.

There is **no unified `mofsbu` CLI with subcommands.** `pyproject.toml` declares one console
script, `mofsbu-ui`. Everything else is a standalone file under `scripts/`: `build.py`,
`run_spec.py`, `verify_registry.py`, `regress_m7.py`, `seed_demo_registry.py`,
`build_metal_descriptors.py`, `check_cases.py`, `regen_golden.py`, `viewer.py`, `check.sh`.
Those scripts contain no logic — they parse args and call the package. Consolidating them behind
one entry point is M9, not a thing that already happened.

---

### 1b. Design intent — not built, signatures still proposals

Do not code against these as though they were real; each currently raises. `NotBuiltYet`
(a `NotImplementedError` subclass) means *missing body*; `EnergyBackendUnavailable` means
*missing install*. Never collapse the two (ground rule 8).

```python
# identity/keys.py — signatures are final and IN THE SCHEMA; bodies return the unset placeholder
def l2_isomer_tag(g, geom=None) -> str            # cis/trans, fac/mer, Δ/Λ        -> M5
def l3_conformer_id(choice_vector=None, geom=None) -> str                        # -> M5

# assembly/join.py — HALF built.  BuildingBlock and open_sites work today:
@dataclass
class BuildingBlock:
    structure_id; graph; geometry; sites; net_charge; provenance
    def open_sites(self, ...) -> list[Site]       # RAISES on a block with no state, by design
def compatible(a, b, *, partner=None) -> Compatibility     # raises NotBuiltYet -> M5
def join(a, b, site_a, site_b, *, ...) -> JoinResult       # raises NotBuiltYet -> M5
def grow(seed_block, partners, *, ...)                     # raises NotBuiltYet -> M5/M6

# assembly/choice.py — NOT WRITTEN -> M5
@dataclass(frozen=True)
class Choice: kind: Literal["A","B","C"]; name: str; value: Any
@dataclass(frozen=True)
class ChoiceVector:
    choices: tuple[Choice, ...]
    def digest(self) -> str
    def replay(self) -> ConstructSpec

# assembly/construct.py — NOT WRITTEN -> M5
def construct(spec, *, seed) -> ConstructResult            # deterministic; emits choice_vector
def enumerate_constructs(spec) -> Iterator[ConstructSpec]  # Kind-B/C branch tree, live-DOF gated

# geometry/placer.py — the multicentre half -> M6.  SIGNATURES NOW IN THE CODE (bodies raise);
# these are exact as of the M6 test-case prep, which is what forced each of them.
@dataclass                    # `cn`/`local_geometry` have NO defaults: that is the Kind-C branch
class Center: element; cn; local_geometry; charge=0; oxidation_state=None; spin_class=None
@dataclass                    # a BRIDGE is two Joins naming the same block — D14, nothing tags it
class Join: center; block: LigandPlacement; site: int; mode; torsion_well
@dataclass
class InterCentreConstraint:                                   # note the British spelling
    centres; mm_distance; mm_lo; mm_hi; bridge_bite_deg; metal_metal_bond
    def window(self, *, tol=0.15) -> tuple[float, float] | None
def place_multicentre(centers, joins, constraints, *, seed=0) -> PlacementResult
def to_rdkit_multicentre(centers, joins, constraints, result) -> Chem.Mol
# geometry/qc.py
class BadIntercentre(NamedTuple): centres; distance; lo; hi; symbols; source
def check_intercentre(coords, constraints, *, metal_idxs, symbols=None) -> list[BadIntercentre]

# geometry/templates.py — NOT WRITTEN.  M6's declared plan B (see §4)
def node_template(name) -> TemplateNode          # "cu_paddlewheel", "fe3_mu3_oxo", "zn4o"
def graft(template, joins, *, seed) -> PlacementResult

# pathways/ — EMPTY package -> M8
@dataclass
class ProxyRecord: concurrent_bond_changes; exchange_lability; coulomb_penalty; bep_estimate; ...
def barrier_proxy(reg, reaction) -> ProxyRecord
def score_path(reg, reaction_ids) -> PathScore   # max_barrier, cumulative_dG, sink_flag, sink_node
def compare_paths(reg, paths) -> PathComparison
def enumerate_paths(reg, target_id, reagent_pool, *, max_steps=6) -> list[list[int]]
```

---

## 2. Remaining milestones

### M5 — Recursive assembly (N=1 path of the general operation) · **L**

**Work** — `assembly/block.py`, `join.py` (frame-alignment compatibility + join), `choice.py`
(ChoiceVector, digest, replay), `construct.py` (deterministic construct + branch-tree enumerator),
atom-map site inheritance, and `identity/l2_isomer_tag` filled in for real.

**Decision gate: C2** — θ_geom (RMSD on rigid core + coordination sphere) and the energy window.
This is the first milestone that actually *manufactures* conformers, so it is the first point the
numbers are forced. Suggested way to set them rather than guess: build the fixture set, plot the
pairwise core-RMSD distribution, and put θ_geom in the valley between the stochastic-duplicate peak
(Kind A) and the real-branch peak (Kind B) — i.e. calibrate from data you now have.

**Exit gate — the three headline validations:**
1. **One node, two routes:** a structure built in two different orders lands on **one** L1 node with
   two incoming provenance edges. (This is the D2 claim, tested.)
2. **cis/trans survive:** anthrarufin–Cu cis and trans persist as distinct records carrying
   different choice-vectors, and are **not** merged by geometric clustering even when near-degenerate
   (this is D10/D11 tested — the discriminator is relevance, not ΔE).
3. **Replay:** reconstructing from a stored `(choice_vector, seed)` reproduces coordinates within
   tolerance. If this fails, stored conformers are frozen coordinates, not regenerable objects.

Plus: `construct` on an ambiguous spec (CN not determined) **raises/branches** rather than defaulting.

**Status: met.** S1–S7 built, all three exit gates pass (`tests/test_m5_exit_gates.py`), C2
called on θ_geom (0.15 Å, calibrated on xTB-relaxed geometries). `WORKPLAN_M5.md` carries the
slice-level detail.

**The milestone is done; its bookkeeping is not, and that is scheduled.** The ledger D-numbers
for C2, the L2-wiring call and `MAX_BITE_MISMATCH_DEG` were never written;
`ALGO_VERSIONS["l3_conformer_id"]` still reads `"0-stub"`. Until that lands, M5 stays here
rather than in `archive/PLAN_completed.md` — see [`WORKPLAN_M6.md`](WORKPLAN_M6.md) §8(a), which
owns it. What is genuinely still unfinished is the carry-over list below.

**Carried out of M5 — unbuilt, so not bugs (see `BUGS.md`'s scope note):**

* **The L3 energy window is unset.** Half of C2. Blocked on `BUGS.md` B9 rather than on missing
  data — the measurement available today is contaminated by the rigid-core defect.
* **Symmetry collapse in the enumerator.** 405 leaves over 5 distinct L1 at CN 6 degree 2. The
  enumerator deliberately does not merge leaves sharing an L1 (cis and trans share one too); now
  that L2 and θ_geom both exist, the collapse they were waiting for is buildable.
* **Driving the enumerator from the run pipeline** for `degree > 1` (S4.1). `runner.plan` still
  refuses it, now because the pipeline does not drive the *tree* rather than because the tree
  does not exist. The **ladder** is driven: `spec.pathways` plans the rung below each product
  and joins the step between them, which is S4.1 for the one shape a mononuclear sweep makes.
* **Which vertices an intermediate leaves empty is not chosen.** `place_mononuclear` fills
  vertices in its own order, so a CN-6 centre carrying four co-ligands comes back with its two
  empty vertices **trans** — and a chelating ligand that subtends ~90° cannot reach them, so
  that rung's step is refused (`chelate_cannot_span`) with the measurement in the task row. The
  ladder therefore connects the unsaturated series but not the co-ligand-saturated one. The fix
  is for the placer to take the vertices to reserve, which is placer work, not pathway work.

---

### M6 — Polynuclear nodes, emerging from joins · **L**

**The framing changed before any of it was built, and the measurements are in
[`WORKPLAN_M6.md`](WORKPLAN_M6.md).** This section described M6 as generalising the placer to
"coordination centers + inter-center constraints" and called it the headline engineering cost.
It is not where the cost is. A Cu paddlewheel builds QC-clean from machinery that already
exists — well-1 site frames plus `_linalg.kabsch` — and hashes to its M2 fixture at **both L0
and L1**, with no solver anywhere. A polynuclear node is what a sequence of joins *produces*;
M–M distance is an **output to validate**, not an input to impose (**D20**).

**What the battery turned up instead — there are two bridge mechanisms, and only one is
supported:**

* **A — multi-atom bridge** (carboxylate O,O; pyrazolate N,N): two donors, one per metal.
  Geometry comes from the donors' lone-pair axes, and `site_frame(..., well=)` already knows
  both — but `sites.model.perceive` stores only well 0, so the syn-syn bridge that makes
  paddlewheels is unreachable. One branch away.
* **B — single-atom bridge** (µ2-OH, µ3-O, µ4-O): one donor, two to four metals. Not
  expressible at all: a donor with two neighbours has one axis, a bare oxo has none, and
  `GEOMETRIES` has no bent CN-2. The fix is a reframing — **a bridging atom is a centre whose
  vertices are metal positions** — and it lands µ3 and µ4 skeletons on the literature values
  exactly.

**Where they collide is the real work.** In a paddlewheel the bridges alone fix M···M. In an
oxo-centred cluster the central atom fixes it too, and the two disagree — by 0.595 Å for
Fe₃(µ₃-O) and 0.491 Å for Zn₄O. A rigid ligand cannot open its O–C–O the way a real one does,
so this is a measured limitation, reported as strain and closed by relaxation. **This is what
`place_multicentre` is for**, and the reason to keep it: the paddlewheel never needs it, and the
oxo clusters cannot be built without it. It is one scalar per edge, not a general solver.

**Work** — the lone-pair well as a Kind-B branch; `bridge_compatible` and a two-point form of
`join`; bridging atoms as centres, plus the bent CN-2 geometry; reconciliation via
`place_multicentre`; vacancy↔vacancy joins for a declared nucleus, with
`metal_metal_distance` in `geometry/distances.py`; lifting `join`'s `n_metals > 1` guard;
`geometry/qc.py` extended with `check_intercentre`; multi-metal `to_rdkit`. *(The CN=5 item this
section used to carry is done — both CN-5 polyhedra are in `site_vectors` and tested.)*

**Ground-truth targets:** a battery, not one motif — paddlewheel, under-bridged Cu₂(µ-O₂CH)₂,
its benzoate analogue, a pyrazolate dimer, Cu₂(µ-OH)₂, Fe₃(µ₃-O) in both valence patterns, and
Zn₄O. Between them: both mechanisms, nuclearity 2/3/4, with and without an M–M bond,
carboxylate and not, symmetric and mixed-valence. **Most of it has to be built** —
`examples.ALL` holds 16 entries and Zn₄O is not among them (see §6).

**Exit gate — the strong one:** the **sequential** route (Cu₂(µ-HCOO) + HCOO) and the
**nucleus-first** route (declare the dimer, add bridges) reach **one** node with **two**
incoming provenance edges. That is M5's "one node, two routes" generalised to polynuclear, and
it is M8's Path A vs Path B made buildable here. Plus: the battery builds QC-clean or is refused
with a stated number; each cluster's M···M is within literature range; an xTB relax does not
tear a node apart; and the failed route `Cu(HCOO)₂ + Cu` is **refused with its measurement**
(1.46 Å), not silently absent.

**Plan B, declared in advance (see §4):** `geometry/templates.py` — place from a stored
reference node geometry and graft ligands onto it (`ebu_tools_v2.PaddlewheelBuilder` is the
seed). Templates are a *geometry source*, not a second identity path — they enter through the
same API. Trigger date **2026-10-14**, kept because the plan asks for a date rather than a mood,
while recording that the risk it guards has largely retired.

---

### M7 — Energy backends + a defensible reference scheme · **S remaining** · *the gate, not the build*

**Built already** (revs 16, 19, 20, 22) — `energy/backends.py` (xTB via tblite, MACE-MP-0,
MACE-OMOL-0, Null) behind one protocol; a `MethodSpec` on every stored number;
`high_spin_multiplicity`; the `relaxed_from` fidelity ladder exercised by a real relax runner;
and `energy/reference.py`, which **refuses** unbalanced and non-isodesmic equations (D17) —
the fix for the design doc's first known risk, gas-phase xTB on isolated highly-charged anions.

**What is left is the exit gate.** `scripts/regress_m7.py` is the harness and it is written in
two stages; **no result from either is recorded anywhere.**

**Exit gate — regression against your own archived work:** re-run the Ni/BTC and Fe/BTC cases
through the new stack and reproduce the rankings in `reports/ni_btc_report.md` and
`reports/fe_btc_report.md` within noise. Same for the solvation sweep. If the new stack can't
reproduce the old results, one of them is wrong and you want to know *now*, not in M8 — which
takes these numbers as given.

* **Stage 1 — `--refs`, minutes, needs tblite.** Recomputes the four reference energies the
  archived Fe(III) run subtracted (sextet Fe³⁺, and the BTC / EDTA / EDDA anions) against
  `legacy/fe_btc_refs.json`. Same method, same species, so this is a direct check that
  `XTBBackend` is the same calculator that produced the archived numbers. It is **not** a check
  of the ranking and does not pretend to be. Tolerance is 1 eV deliberately: the archived ligand
  geometries came from `ebu_core.LigandBuilder` and these come from `geometry.embed`, so a few
  tenths is two ETKDG conformers, not a disagreement about energy.
* **Stage 2 — `--rankings`, not built.** It raises and says why: it needs the archived candidate
  set, which means it is not a re-run but a re-derivation under a different construction path.
  That is the substantive half of this gate and the reason M7 is still open.

Note the archived equation needs `strict=False` to evaluate at all, and the value carries
`isodesmic=False`. That is the point: the comparison is between the old number and the new one,
not an endorsement of the old scheme.

---

### M8 — Reactions + pathways (the actual contribution) · **L**

**Work** — `pathways/reaction.py`, `proxy.py`, `score.py` (max barrier, cumulative ΔG,
rate-limiting step, **sink detection**), `search.py`, and the A-vs-B driver script.

**Decision gates: C6** (barrier-proxy form) and **C7** (partner dependence — factorized vs. stored
matrix). Both are forced here and not before.

**Exit gate — the first real result:** the Cu paddlewheel **Path A (nucleus-first) vs. Path B
(sequential)** comparison produces a ranked answer with per-step numbers, intermediates stored as
first-class registry nodes, and a written report in `reports/`. Secondary: sink detection
reproduces the EDTA sequestration finding from `solvation_report.md`, now generalized from
endpoints to intermediates.

---

### M9 — Usability skin · **S** · *fold in opportunistically*

CLI subcommands, a `report` command that renders a structure/pathway summary, optional
`viz.py` (port `ebu_viz.py`), and a notebook-free query helper. Do this in slivers as each
milestone lands, not as a block at the end.

## 3. Open decision gates (placed where code forces them)

| Gate | Milestone | Forced by | What you need in hand to decide |
|---|---|---|---|
| **C2** — L3 thresholds (θ_geom **called**, energy window open) | **M5**, at the first conformer generation | clustering can't run without numbers | θ_geom = 0.15 Å, calibrated. The window waits on [B9](BUGS.md#b9) — the measurement available today is contaminated by the rigid-core defect |
| **C9** — polynuclear multiplicity | **M6**, at the first manufactured node | multiplicity is in L0, so it decides identity | the two ground-truth fixtures disagree: `cu_paddlewheel` declares 1 (AF-coupled d⁹–d⁹) where `combined_multiplicity` gives 3; `fe3_mu3_oxo` declares 16 and the additive rule agrees. *Leaning: coupling is a Kind-C branch — require it stated, raise otherwise* |
| **C10** — when is there an M–M edge? | **M6**, alongside C9 | `METAL_METAL` is in the certificate, so a wrong answer is a wrong identity | at 2.673 Å two Cu are bonded, at 5.516 Å they are not. A silent distance threshold is exactly what ground rule 5 forbids — branch or declare |
| **C6** — barrier proxy | **M8**, before the first path score | `PathScore.max_barrier` needs a definition | whether the paddlewheel A-vs-B ordering is stable under the cheap proxies alone |
| **C7** — partner dependence | **M8**, alongside C6 | `ease(site, partner)` is called at query time | how many (site, partner) pairs you actually intend to screen — the factorization only pays off if that number is large |

*Resolved: C1 → D10 · C4 → D12 · C5 → D18 · C8 → curated tables.*

One further call is open and is **not** a milestone gate — the L2 backfill-versus-version
decision ([B2](BUGS.md#b2)). M5 forces it, because M5 is what fills `l2_isomer_tag` in.

Rule for all of them: when a gate is called, **write the resolution into the design doc's
Decision Ledger with a new D-number and a changelog line**, then code against it. A gate
resolved only in your head is how the two documents drift apart.

---

## 4. Risk register and the escape hatches

| Risk | Signal it's happening | Escape hatch |
|---|---|---|
| **M6 overruns** | the battery still will not build | switch to `geometry/templates.py` (grafting onto stored reference nodes). **Trigger: 2026-10-14** — written down at M6's start, as this row asks. Largely retired in advance: the paddlewheel builds QC-clean at 2.673 Å from existing machinery, so the constrained placer this hatch was written against is not on the critical path. What remains is reconciling oxo-centred clusters, measured at ~0.5 Å. |
| **Two machines disagree about identity** | a fixture hash differs laptop vs. workstation | already mitigated by ground rule 10 (D16: the certificate is pure Python; no optional native package can produce a key). The golden-hash test in M2 is what catches it — do not skip it. |
| **L3 conformer explosion** | thousands of near-identical rows per (L1, L2) | choice-vector dedup runs *before* geometric clustering; cap conformers per (L1,L2); `TORSION_FREE` sites never branch (that tag is the guard) |
| **xTB numbers can't carry route claims** | M7 regression reproduces rankings but absolute ΔG look implausible | keep M8 claims *relative and within-metal*; the reaction-balanced reference scheme is the gate on any quantitative statement |
| **Identity retrofit pressure** | a temptation in M5/M6 to "just add a flag" to L1 | the M2 discrimination table is the contract; changing it means a version bump and a re-hash of the corpus, which is exactly the cost that should make you think twice |
| **Scope creep into periodic frameworks** | topology/net questions start appearing in tickets | out of scope for v1 (design §11) — record them as future work; MOFid/MOFkey is the reference oracle when you get there |

---

## 5. What the tool can do after each remaining milestone

| After | You can… |
|---|---|
| *(today)* | store, dedupe, query and **see** every structure built; ask a stored structure what sites it has, which are open, and which are worth spending QM on; build a mononuclear structure from stored blocks and rebuild it exactly from its provenance (M5) |
| M6 | build real SBUs — paddlewheels, µ₃-oxo trimers, Zn₄O — not just mononuclear nodes, and be told which formation route cannot reach one |
| M7 *(gate)* | trust the relative energies attached to any of it, because they reproduce the archived results |
| M8 | compare two synthesis routes to the same product and say which is more viable, and why |

---

## 6. Test strategy (what "done" means, per layer)

**Fixture set.** Hand-written typed graphs live in `mofsbu/examples.py`, not `tests/fixtures/`,
so scripts, notebooks and the demo registry all draw on the same ones. `examples.ALL` currently
holds **16**: water, formate, formic acid, BTC (both protonation states), Cu paddlewheel,
Fe₃-µ₃-oxo in two valence patterns, the Zn₂ bridged/chelated discrimination pair, the
Pt(NH₃)₂Cl₂ build-order pair, Fe hexaaqua high/low spin, and the cyclohexane / two-cyclopropanes
1-WL collision. Parameterised families (`paddlewheel`, `hexaaqua`, `aqua_carboxylate`,
`bipyridine`, `metal_bipy`) sit alongside them.

**What this list used to claim and does not have:** Zn₄O, mononuclear Zn/BDC, BDC, bipy as a
`TypedGraph`, EDTA, anthrarufin, and the cis/trans anthrarufin–Cu pair. M6 needs several of
them and builds them as its first slice — see [`WORKPLAN_M6.md`](WORKPLAN_M6.md) §6.

Four kinds of test, each with a job:

- **Invariance/property tests** — atom-order shuffles, insert idempotency, `construct` determinism.
  These guard the ground rules in §0 and are the ones that catch silent identity corruption.
- **Discrimination tests** — pairs that *must* differ, and pairs that *must not*, at each of L0/L1/L2/L3.
  Write them as one table so the identity boundaries are readable in one place.
- **Golden/regression tests** — the WL hashes of the fixture set (two-machine parity), and the
  Ni/Fe/BTC + solvation numbers from `reports/` (M7).
- **Negative tests** — ambiguous spec raises instead of defaulting; a bad geometry is rejected by QC
  rather than silently stored; a version bump marks rows stale rather than re-labelling them.

Not aiming for coverage percentages. Aiming for: **every claim in the design doc's ledger has a test
that would fail if the claim stopped being true.**

---

## 7. What to do next, concretely

1. **Run M7's stage 1 and record the result** — `python scripts/regress_m7.py --refs --json …`
   on a machine with tblite. Minutes of compute, and it is a *check* rather than a build. Then
   decide whether stage 2's re-derivation is worth building now or at M8, when the numbers
   actually carry a claim.
2. **Call the L2 question** ([B2](BUGS.md#b2)) before writing M5 code, not during.
3. **M5 fixtures before M5 code** — the cis/trans anthrarufin–Cu pair and the two-route
   structure as failing tests first. This is the same move that made M2 cheap: it forces the
   D10/D11 commitments at the point where they are still free.
4. **Write the M6 fallback date down when M6 starts** (§4, first row). A plan-B with no trigger
   is not a plan-B.
