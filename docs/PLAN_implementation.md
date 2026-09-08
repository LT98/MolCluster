# Implementation Plan — code, functionalities, and milestone order

**Companion to** [`DESIGN_registry_assembly.md`](DESIGN_registry_assembly.md). That doc says *what*
is being built and *why* (ledger D1–D13). This doc says *what code exists*, *in what order*, and
*what has to be true before the next thing starts*.

**Form:** dependency-ordered milestones, no calendar dates. Each milestone has entry conditions,
the modules it lands, an exit gate (a test that must pass, not a feeling), and an effort size
(S ≈ a sitting or two, M ≈ a week of focused evenings, L ≈ multi-week / headline cost).

**Open checkpoints C2/C5/C6/C7/C8 are scheduled as decision gates** at the milestone where code
first forces the call — see §6. Nothing is pre-resolved here; the design doc's leanings stand.

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
10. **Two-machine parity.** The primary L1 hash is computed by a path available on *both* machines
   (networkx WL). `pynauty` is a Linux-side *verifier*, never the primary key producer — otherwise
   the laptop and the workstation disagree about identity.

---

## 1. Package map and API surface

Signatures below are the contract to code against. Types are indicative (dataclasses,
`from __future__ import annotations` everywhere).

### `mofsbu/graph/` — the typed molecular graph (identity substrate)

```python
# types.py
class EdgeType(StrEnum):
    COVALENT = "cov"; DATIVE = "dat"; MU2 = "mu2"; MU3 = "mu3"; METAL_METAL = "mm"

@dataclass(frozen=True)
class NodeLabel:
    element: str
    formal_charge: int = 0
    oxidation_state: int | None = None     # per-center, D12 — set on metals
    spin_class: str | None = None          # per-center, D12 — e.g. "hs", "ls"
    role: str = "atom"                     # "atom" | "metal" | "bridge_anchor"

class TypedGraph:
    def add_atom(self, label: NodeLabel) -> int
    def add_bond(self, i: int, j: int, etype: EdgeType, order: float = 1.0) -> None
    def neighbors(self, i: int, etype: EdgeType | None = None) -> list[int]
    def metals(self) -> list[int]
    def subgraph(self, idxs: Iterable[int]) -> TypedGraph
    def relabel(self, mapping: dict[int, int]) -> TypedGraph
    def to_networkx(self) -> nx.Graph
    @classmethod
    def from_rdkit(cls, mol: Chem.Mol, *, metal_idxs=(), dative_pairs=()) -> TypedGraph
    @classmethod
    def from_json(cls, blob: bytes) -> TypedGraph
    def to_json(self) -> bytes                       # the stored typed_graph_blob
```

```python
# canon.py
def wl_hash(g: TypedGraph, iterations: int = 3) -> str          # PRIMARY key producer (both machines)
def nauty_canon(g: TypedGraph) -> list[int]                     # optional; Linux verifier + canonical order
def canonical_order(g: TypedGraph) -> list[int]                 # nauty if available else WL-refine + deterministic tiebreak
def is_isomorphic(a: TypedGraph, b: TypedGraph) -> bool         # VF2 — collision resolver ONLY
```

Note: a µ2-carboxylate is **not** a duplicated fragment — its atoms appear once, and its two O
atoms each carry a `DATIVE` edge to a different metal, with the pair additionally tagged `MU2`.

### `mofsbu/identity/` — the composite L0–L3 key

```python
def l0_composition(g: TypedGraph) -> str        # "Cu2_C8H20O8_q0_s3|Cu:+2/hs,Cu:+2/hs"
def l1_graph_hash(g: TypedGraph) -> str         # wl_hash + algo_version tag
def l2_isomer_tag(g: TypedGraph, geom: Geometry | None) -> str   # cis/trans, fac/mer, Δ/Λ  (M5)
def l3_conformer_id(choice_vector: ChoiceVector, geom: Geometry) -> str   # provenance-primary (D11)
def block_id(l0: str, l1: str, l2: str = "", l3: str = "") -> str
```

`l2_isomer_tag` ships in M2 with a fixed signature returning `""`; M5 fills it in. The *signature*
must be right in M2 because it lands in the schema.

### `mofsbu/registry/` — SQLite + content-addressed blobs

```python
class BlobStore:
    def __init__(self, root: Path)
    def put(self, data: bytes) -> str            # sha256 hex; writes store/ab/cdef….bin
    def get(self, digest: str) -> bytes
    def path(self, digest: str) -> Path

class Registry:                                   # context manager; owns the sqlite conn + migrations
    def __init__(self, db_path: Path, store: BlobStore)
    def migrate(self) -> None

# api.py — the ONLY write surface
def put_structure(reg, g: TypedGraph, *, l2: str = "", provenance: Provenance | None = None) -> int
def put_geometry(reg, structure_id: int, coords: Coords, *, fidelity: Fidelity, method: MethodSpec,
                 energy: float | None = None, converged: bool | None = None,
                 relaxed_from: int | None = None, choice_vector: ChoiceVector | None = None) -> int
def put_reaction(reg, product_id: int, reagent_ids: list[int], *, atom_map: AtomMap,
                 conditions: Conditions, dG: float | None = None,
                 barrier_proxy: ProxyRecord | None = None) -> int
def get_structure(reg, ref: int | str) -> StructureRecord            # id or block_id
def best_geometry(reg, structure_id: int, min_fidelity: Fidelity | None = None) -> GeometryRecord
def find(reg, *, l0=None, l1=None, formula=None, has_metal=None, open_site_type=None) -> list[StructureRecord]
def export_xyz(reg, geometry_id: int, path: Path) -> Path            # files are a VIEW, not a record
```

`put_structure` is **idempotent** on `(L0, L1, L2)`: same key returns the existing id and adds a
provenance edge. That single property is what makes "one node, two build routes" real.

### `mofsbu/sites/` — perceive once, refresh accessibility

```python
# perception.py  (port of legacy/donor_perception.py, ~verbatim, re-keyed to canonical indices)
@dataclass(frozen=True)
class DonorSite:
    canonical_idx: int; donor_type: str; labile: bool
    h_canonical_idx: int | None; charge_after: int
def find_donor_sites(g: TypedGraph, mol: Chem.Mol) -> list[DonorSite]
def activate(mol: Chem.Mol, sites: list[DonorSite]) -> tuple[Chem.Mol, list[DonorSite]]

# frames.py  (D13 — promoted out of legacy GeometryPlacer)
@dataclass(frozen=True)
class SiteFrame:
    origin: Vec3; axis_hat: Vec3; ref_hat: Vec3; mode: Literal["aligned","tilted","fallback"]
class LiveDOF(StrEnum): TORSION_LIVE = "live"; TORSION_FREE = "free"
class BindingMode(StrEnum): MONODENTATE="mono"; CHELATE="chel"; BRIDGE_SYN_SYN="mu2ss"; ...
def site_frame(mol, conf, donor_idx: int, donor_type: str) -> SiteFrame       # <- _donor_placement_frame
def live_dof(donor_type: str) -> LiveDOF
def binding_modes(donor_type: str) -> tuple[BindingMode, ...]
def torsion_wells(donor_type: str, mode: BindingMode) -> tuple[float, ...]    # discrete wells (D13)

@dataclass(frozen=True)
class Site:                       # what site_catalog stores
    canonical_idx: int; donor_type: str; frame: SiteFrame
    live_dof: LiveDOF; binding_modes: tuple[BindingMode, ...]
    labile: bool; charge_after: int

# state.py
class SiteStatus(StrEnum): OPEN="open"; OCCUPIED="occ"; BLOCKED="blk"
def refresh_state(structure: StructureRecord, geom: GeometryRecord) -> list[SiteState]
def buried_volume(coords, site: Site, radius: float = 3.5) -> float

# inherit.py
def inherit_sites(parent_sites: list[Site], atom_map: AtomMap, consumed: set[int]) -> list[Site]
```

### `mofsbu/descriptors/` — the shared lookup substrate (§6.6)

```python
@dataclass(frozen=True)
class DonorDescriptor:
    donor_type: str; pka: float; pka_sigma: float; hsab: Literal["hard","borderline","soft"]
    default_denticity: int; live_dof: LiveDOF; source: str; source_version: str

@dataclass(frozen=True)
class MetalDescriptor:
    symbol: str; charge: int; ionic_radius: float; hsab: str
    preferred_cn: tuple[int, ...]; d_electrons: int
    exchange_lability: Literal["fast","moderate","slow"]; source: str; source_version: str

DONOR_TABLE: dict[str, DonorDescriptor]
METAL_TABLE: dict[tuple[str, int], MetalDescriptor]

# ease.py
@dataclass(frozen=True)
class EaseRecord:
    components: dict[str, float]      # deprotonation / electronic / steric / marginal_dE
    scalar: float; fidelity: Fidelity; confidence: float; provisional: bool; method: MethodSpec
def activation_ease(site: Site, *, partner: MetalDescriptor | None = None,
                    conditions: Conditions | None = None, state: SiteState | None = None) -> EaseRecord
def hsab_match(donor: DonorDescriptor, metal: MetalDescriptor) -> float     # C7 factorization
```

### `mofsbu/geometry/` — local geometries, the multi-center placer, QC

```python
# local.py
def site_vectors(geometry: str, n: int, d: float = 2.05) -> np.ndarray
#   planar | tetrahedral | octahedral | trigonal_bipyramidal | square_pyramidal | trigonal_planar
#   (legacy mapped CN=5 to "planar" — that is wrong and is fixed here)

# placer.py  — THE headline cost (D12)
@dataclass
class Center:
    element: str; charge: int; oxidation_state: int; spin_class: str
    cn: int; local_geometry: str
@dataclass
class Join:
    center: int; block: BuildingBlock; site: Site; mode: BindingMode; torsion_well: int
@dataclass
class InterCenterConstraint:
    centers: tuple[int, int]; mm_distance: float | None; bridge: BridgeSpec | None

class MultiCenterPlacer:
    def place(self, centers: list[Center], joins: list[Join],
              constraints: list[InterCenterConstraint], *, seed: int) -> PlacementResult
    # PlacementResult(coords, strain, per_join_torsion, qc: QCReport, converged: bool)

# templates.py  — plan-B / accelerator for M6 (see §7 risk)
def node_template(name: str) -> TemplateNode      # "cu_paddlewheel", "fe3_mu3_oxo", "zn4o"
def graft(template: TemplateNode, joins: list[Join], *, seed: int) -> PlacementResult

# qc.py  (port of legacy/geometry_qc.py + new polynuclear checks)
def check_clashes(coords, ...) -> QCResult
def check_ml_bonds(coords, center_idx, donor_idxs, ...) -> QCResult
def check_ring_planarity(free, placed, ...) -> QCResult
def check_intercenter(coords, constraints) -> QCResult        # NEW: M–M distance, bridge bite angle
def qc_report(...) -> QCReport
```

### `mofsbu/assembly/` — the recursive BuildingBlock

```python
@dataclass
class BuildingBlock:
    structure_id: int | None; graph: TypedGraph; geometry: Coords | None
    sites: list[Site]; net_charge: int; provenance: Provenance
    def open_sites(self, state: list[SiteState] | None = None) -> list[Site]

# choice.py  (D13)
@dataclass(frozen=True)
class Choice:                       # one branch point
    kind: Literal["A","B","C"]; name: str; value: Any
@dataclass(frozen=True)
class ChoiceVector:
    choices: tuple[Choice, ...]
    def digest(self) -> str
    def replay(self) -> ConstructSpec

# join.py
def compatible(a: Site, b: Site, *, partner: MetalDescriptor | None = None) -> Compatibility
    # frame-alignment feasibility + strain estimate, NOT just "both open"
def join(a: BuildingBlock, b: BuildingBlock, sa: Site, sb: Site, *,
         mode: BindingMode, torsion_well: int, seed: int) -> JoinResult
    # JoinResult(block, atom_map, choice_vector, strain, qc)

# construct.py
def construct(spec: ConstructSpec, *, seed: int) -> ConstructResult     # deterministic; emits choice_vector
def enumerate_constructs(spec: EnumSpec) -> Iterator[ConstructSpec]     # Kind-B/C branch tree, live-DOF gated
```

### `mofsbu/energy/` — backends behind one protocol

```python
class Fidelity(IntEnum): RAW = 0; FF = 1; ML = 2; XTB = 3; DFT = 4

@dataclass(frozen=True)
class MethodSpec:
    code: str; version: str; method: str; solvent: str | None
    charge: int; multiplicity: int; extras: dict

class EnergyBackend(Protocol):
    def single_point(self, atoms, spec: MethodSpec) -> float
    def relax(self, atoms, spec: MethodSpec, *, fmax: float, steps: int) -> RelaxResult

class XTBBackend(EnergyBackend): ...      # tblite  (port of legacy xtb_energy*.py)
class MACEBackend(EnergyBackend): ...     # port of legacy energy_model.py / dft_predict.py
class NullBackend(EnergyBackend): ...     # tests

def high_spin_multiplicity(symbol: str, charge: int) -> int
def reaction_balanced_energy(reg, reaction) -> float    # M7 — fixes the charged-ion reference problem
```

### `mofsbu/pathways/` — the route-design layer

```python
@dataclass
class ProxyRecord:
    concurrent_bond_changes: int; exchange_lability: str
    coulomb_penalty: float | None; bep_estimate: float | None
    scalar: float; fidelity: Fidelity

def barrier_proxy(reg, reaction) -> ProxyRecord
def score_path(reg, reaction_ids: list[int]) -> PathScore
    # PathScore(max_barrier, cumulative_dG, rate_limiting_step, sink_flag, sink_node)
def compare_paths(reg, paths: dict[str, list[int]]) -> PathComparison
def enumerate_paths(reg, target_id: int, reagent_pool: list[int], *, max_steps: int = 6) -> list[list[int]]
```

### `scripts/` — the thin CLI (one entry point, subcommands)

`mofsbu ingest` (legacy .xyz → registry) · `build` · `enumerate` · `relax` · `sites` · `route` ·
`report` · `export`. Scripts contain no logic — they parse args and call the package.

---

## 2. Milestone spine

Critical path: **M0 → M2 → M3 → M4 → M5 → M6 → M8**. M1 is off the critical path (do it whenever
you want a low-energy session). M7 can run in parallel with M5/M6 once M3 exists.

```
M0 rails ─► M2 graph+identity ─► M3 registry ─► M3.5 viewer ─► M4 sites ─► M5 assembly ─► M6 placer ─► M8 pathways
    │                                 │                                          │             │
    └──► M1 descriptors ──────────────┴──────────────────────────────────────────┘             │
                                      └──► M7 energy + reference ────────────────────────────--┘
```

---

### M0 — Consolidate to one lineage + build the rails · **S**

**Why first:** the design doc calls this out (§2 "code hygiene"). Two git commits and five parallel
script versions is a correctness hazard, not just untidiness.

**Work**
- Declare **`legacy/ebu_core.py` the single salvage source**. `ebu_tools.py`, `ebu_tools_v2.py`,
  `ebu_toolkit.py` and all notebooks are archive-only — read for ideas, never ported wholesale.
  (Exception: `ebu_tools_v2.PaddlewheelBuilder.build_paddlewheel` is the seed for M6 templates.)
- `mofsbu/versions.py`, `mofsbu/config.py` (paths to `data/registry.db`, `data/store/`),
  structured logging, `mofsbu/types.py` for the shared `Vec3/Coords/AtomMap/Provenance` aliases.
- `pytest` + `ruff` + a `make check` (or `scripts/check.sh`) that runs both. `tests/conftest.py`
  with a `tmp_registry` fixture.
- Delete the `.ipynb_checkpoints` tree from `legacy/`; add a `legacy/README.md` recording the
  salvage verdicts from §5 below.

**Exit gate** — `pip install -e .` + `pytest` green on **both** machines; `mofsbu --version` prints
the pinned algo-version block; one smoke test imports every subpackage.

---

### M1 — Descriptor layer · **S** · *(off critical path — do early anyway, it's cheap)*

**Why early:** it is a lookup table, it unblocks C5/C6/C7, and it means no structure ever has an
empty ease field.

**Work** — `descriptors/donors.py`, `descriptors/metals.py`, `descriptors/ease.py` (heuristic tier
only — no QM in this milestone). Coverage target: every `donor_type` that `perception.py` can emit,
and the metals in play (Zn²⁺, Cu²⁺, Ni²⁺, Fe²⁺/³⁺, Zr⁴⁺, Cr³⁺, Al³⁺). Every row carries
`source` + `source_version`.

**Decision gate: C8** (where the numbers come from, how provenance is pinned) — must be answered
before the first row is typed, or you'll be re-sourcing the table later.

**Exit gate** — a test that iterates every donor type in `LABILE_DONOR_PATTERNS` + the neutral
classifier and asserts a `DONOR_TABLE` row exists (**no silent gaps**); `activation_ease` returns a
`fidelity="heuristic"` record with `provisional=True` for in-pocket / H-bonded donors (the
anthrarufin peri-OH case).

---

### M2 — Typed graph + identity · **M** · *the foundation; get it wrong and everything downstream is wrong*

**Key insight that de-risks this milestone: identity needs no 3D.** You can hand-write the Cu
paddlewheel and Fe₃-µ₃-oxo graphs as fixtures and fully test polynuclear identity **before the
placer exists**. Do exactly that — it front-loads the D12 commitments where they're cheap.

**Work** — `graph/` (types, TypedGraph, from_rdkit with explicit dative pairs), `graph/canon.py`
(WL primary, nauty optional, VF2 collision), `identity/` (L0 with per-center ox-state/spin, L1,
L2 stub with final signature, L3 stub, `block_id`).

**Exit gate — all of these are tests:**
1. **Build-order invariance:** the same molecule assembled in 3 different atom orders → identical L1.
2. **Coordination-isomer discrimination:** µ2-bridging vs. two terminal carboxylates → **different L1**.
3. **L1 correctly does *not* separate cis/trans** of the same connectivity (that's L2's job) — this
   test documents the boundary and fails loudly if someone "fixes" L1 to over-discriminate.
4. **Mixed valence:** Fe₃(µ₃-O) with FeII/FeIII/FeIII vs. all-FeIII → **different L0** (D12).
5. **Two-machine parity:** the WL hash of the fixture set is identical on laptop and workstation
   (commit the expected hashes as a golden file).
6. **Collision path exercised:** a synthetic WL collision resolves correctly through VF2.

---

### M3 — Registry v1 + the walking skeleton · **M**

**Work** — `registry/schema.sql` (the §5 tables, with `algo_versions`, a `migrations` table, and
indices on `l1_graph_hash`, `(structure_id, fidelity)`), `BlobStore`, `api.py`,
`scripts/ingest.py`.

**Denormalised query columns on `structures`** (§2.5 — cheap now, a migration later):

| Axis | Columns |
|---|---|
| composition / metal / charge / spin | `formula`, `metals` (sorted symbols), `n_metals`, `net_charge`, `multiplicity`, `ox_states` |
| connectivity / binding mode | `max_bridge_class`, `binding_modes` (set), `donor_types` (set), `n_open_sites`, `cn_by_centre` |
| energy / fidelity | on `geometries`: `fidelity`, `converged`, `energy`, `method_id`; plus `best_fidelity` on `structures` |
| provenance / route | on `reactions`: `product_id`, `reagent_ids`, `choice_vector_digest`, `depth`; index both directions |

Plus **`display_label`** on `structures`: a derived, human-readable string
(`"Cu2(mu-O2CH)4 · paddlewheel · q0"`). Non-unique, regenerable, **never used for lookup** —
identity stays the hash.  This is the filename problem solved properly rather than reintroduced:
a UI that lists 64-character digests is unusable, but meaning must not creep back into the key.

**Walking skeleton (do this before polishing anything):** SMILES → TypedGraph → L0/L1 →
`put_structure` → `put_geometry(raw)` → `export_xyz` round-trip. A thin end-to-end path this early
is what surfaces integration mistakes while they're still cheap.

**Real-data test — BLOCKED, needs a call.** Ingesting `legacy/ebu_outputs/*.xyz` (90 files) is
stuck on something structural, not incidental: **an `.xyz` file carries no charge and no
multiplicity**, and `l0_composition` refuses to guess either (ground rule 5).  Bond perception
from coordinates is straightforward; charge and spin are not recoverable from geometry.  Options,
in increasing order of effort: (a) a manifest CSV mapping file -> (charge, multiplicity), part
auto-filled from the legacy filenames that encode it (`_q-4`, `Zn2+`) and part filled by hand;
(b) infer charge from perceived ligand protonation plus metal oxidation state, which is a guess
wearing a rule; (c) skip the corpus and let M6 regenerate it.  (a) is the honest one and keeps the
filename as a migration *input* that is then discarded — meaning still does not live in filenames
going forward.  Until this is decided, M3's real-data test is the fixture round-trip instead.

**Original text:** ingest the ~200 `legacy/ebu_outputs/*.xyz` files. Their *filenames* carry the
old meaning (`T3b_Zn2+_BTC[s2]_bipy[s1]_..._q-1`) — parse them into records, then **throw the
filename away**. This both migrates the corpus and proves the "meaning lives in rows" claim on
messy real data. Expect duplicates to collapse; that collapse count is itself a result worth noting.

**Exit gate** — inserting the same structure twice yields one row and two provenance edges;
blob round-trip is byte-identical; the best-geometry pointer moves when a higher-fidelity geometry
lands and does *not* move for an equal/lower one; the full legacy corpus ingests with zero
unparseable records (or an explicit, listed skip set).

---

### M3.5 — Output viewer · **S/M** · *sequenced before M4 on purpose*

**Why here and not at the end:** M4 is the site-frame work — outward axes, torsion references,
binding modes.  Debugging that from printed floats is miserable; you want arrows on atoms.  The
legacy ingest puts ~90 real geometries in the registry at M3, so the viewer has content the day
it exists.

**Work** — `mofsbu/ui/`: a FastAPI app over the registry (read-only) plus one self-contained HTML
page.  Endpoints: `GET /structures` (filter + sort + paginate over the columns above),
`GET /structures/{id}`, `GET /geometries/{id}/xyz` (from the blob store), `GET /spec/{id}`
(the run spec that produced it, for copy-paste).  Front end: filter sidebar, results table,
3Dmol.js pane.  `python -m mofsbu.ui` opens it.

**Read-only by decision.** The viewer cannot write to the registry, so it needs no job state,
no queue, no failure handling — and it stays disposable.  To act on something you select, copy
its spec into a notebook or the CLI.

**Exit gate** — filter the ingested corpus by metal and by binding mode and get the right counts;
select a row and see its geometry render; `display_label` is legible for every ingested record;
the app opens read-only (a write attempt fails loudly, and there is a test asserting that).

---

### M4 — Sites v1 · **M**

**Work** — `sites/perception.py` (port `donor_perception.py`, re-keyed to **canonical** indices
from `canonical_order`), `sites/frames.py` (**promote `_donor_placement_frame` out of the placer**
— D13's central move), `live_dof` / `binding_modes` / `torsion_wells` tables, `site_catalog` +
`site_state` writes, heuristic ease from M1, `sites/inherit.py`.

**Decision gate: C5** — confirm the zero-QM heuristic floor (pKa table + HSAB tag +
provisional-flag) is the accepted v1 floor before `site_state` rows start accumulating.

**Exit gate**
- **Canonical-index stability:** shuffle input atom order → identical `site_catalog` rows.
- **Perceive-once:** an assertion/counter proving `find_donor_sites` runs once per structure and
  `refresh_state` never re-perceives.
- **Frame reproducibility:** a frame read back from the DB reproduces the re-derived frame to 1e-6
  (this is the whole point of promoting it out of the stochastic search).
- Known-molecule coverage: BTC → 3 carboxylate sites each with mono/chelate/bridge modes;
  EDTA → 4 carboxylate O + 2 backbone N; bipy → 2 pyridyl N with `TORSION_LIVE`;
  aqua/primary-amine → `TORSION_FREE`.

---

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

---

### M6 — Multi-center geometry placer · **L** · **the headline cost and the headline risk**

**Work** — `geometry/placer.py`: coordination centers + inter-center constraints (M–M distance,
bridge bite angle) + per-center local geometry; `geometry/qc.py` extended with
`check_intercenter`; fix the CN=5 local geometries.

**Ground-truth targets:** Cu₂(µ-O₂CR)₄ paddlewheel and Fe₃(µ₃-O) trimer, built **from scratch**.

**Exit gate — the strong one:** the placer-built paddlewheel and Fe₃-oxo graphs land on the
**same L1 hash as the hand-written M2 fixtures**. That single test proves the placer, the graph
extraction, and identity all agree. Plus: QC passes (no clashes, M–M within literature range,
bridge angles sane) and an xTB relax doesn't tear the node apart.

**Plan B, declared in advance (see §7):** if constrained multi-center placement stalls, fall back to
`geometry/templates.py` — place from a stored reference node geometry and graft ligands onto it
(`ebu_tools_v2.PaddlewheelBuilder` is the seed). Honest, much cheaper, and it keeps M8 reachable.
Templates are a *geometry source*, not a second identity path — they enter through the same API.

---

### M7 — Energy backends + a defensible reference scheme · **M** · *(parallelizable with M5/M6)*

**Work** — `energy/backends.py` (XTB via tblite, MACE, Null) behind one protocol; `MethodSpec`
recorded on every geometry; `high_spin_multiplicity` ported; the fidelity ladder
(`relaxed_from` chains) exercised; **`reaction_balanced_energy`** — the fix for §11's first risk
(gas-phase xTB on isolated highly-charged anions).

**Exit gate — regression against your own archived work:** re-run the Ni/BTC and Fe/BTC cases
through the new stack and reproduce the rankings in `docs/reports/ni_btc_report.md` and
`fe_btc_report.md` within noise. Same for the solvation sweep. If the new stack can't reproduce the
old results, one of them is wrong and you want to know *now*, not in M8.

---

### M8 — Reactions + pathways (the actual contribution) · **L**

**Work** — `pathways/reaction.py`, `proxy.py`, `score.py` (max barrier, cumulative ΔG,
rate-limiting step, **sink detection**), `search.py`, and the A-vs-B driver script.

**Decision gates: C6** (barrier-proxy form) and **C7** (partner dependence — factorized vs. stored
matrix). Both are forced here and not before.

**Exit gate — the first real result:** the Cu paddlewheel **Path A (nucleus-first) vs. Path B
(sequential)** comparison produces a ranked answer with per-step numbers, intermediates stored as
first-class registry nodes, and a written report in `docs/reports/`. Secondary: sink detection
reproduces the EDTA sequestration finding from `solvation_report.md`, now generalized from
endpoints to intermediates.

---

### M9 — Usability skin · **S** · *fold in opportunistically*

CLI subcommands, a `report` command that renders a structure/pathway summary, optional
`viz.py` (port `ebu_viz.py`), and a notebook-free query helper. Do this in slivers as each
milestone lands, not as a block at the end.

---

## 2.5 Interface architecture — inputs and outputs are two problems

Notebooks were the surrogate because there was nowhere else to put an interface.  The durable
answer is not one app; it is a split, because the two halves have opposite risk profiles.

**Inputs — a declarative spec file, and every UI is a thin editor of it.**
D13 already says `construct` is a deterministic function of `(choice_vector, seed)`.  Make that
spec a YAML file and the consequences follow: a run is reproducible by handing someone the file,
a notebook cell is five lines, a CLI invocation is `mofsbu build spec.yaml`, and a future form is
a spec editor rather than a new source of truth.  Nothing about the interface becomes
load-bearing, which is what makes it safe to keep changing it while the science is still moving.
Notebooks stay the working input surface for now — they are fine at this job once the API behind
them is thin.

**Outputs — a read-only viewer over the registry.**
Read-only is the whole risk story: it cannot corrupt anything, it needs no job state, and it can
be thrown away and rebuilt in a day.  That is what makes a browser UI affordable now when a real
GUI is not.  FastAPI + a single self-contained HTML page + 3Dmol.js; `python -m mofsbu.ui`.

**3D carries over either way.**  Legacy `ebu_viz.py` used py3Dmol, which is the same 3Dmol.js the
browser page uses — the notebook path and the viewer path share a renderer, so reinstating
visualisation is not two pieces of work.

**What this costs M3:** the denormalised query columns and `display_label` above.  All four filter
axes are in scope, including provenance — those columns are nearly free to add now and require a
migration plus a re-walk of the reaction DAG later.

**What is deliberately not being built yet:** a desktop GUI (too much surface while iterating),
and any interface that launches computations (job state is where these tools rot).  Both stay
reachable: the spec file is the seam a richer input UI plugs into, and a runner can be added
behind the same API once the pipeline stops changing shape.

---

## 3. What the tool can *do* after each milestone (functionality ledger)

| After | You can… |
|---|---|
| M0 | install and test the package identically on both machines |
| M1 | ask "how easy is this donor to activate, and with which metal" from tables alone, zero compute |
| M2 | answer "are these two structures the same?" for arbitrary polynuclear species |
| M3 | store, recall, and deduplicate every structure you've ever built; the legacy corpus is queryable |
| M3.5 | **see it** — filter and search the registry in a browser, render any geometry in 3D |
| M4 | ask a stored structure "what sites do you have, which are open, how accessible, how easy" |
| M5 | build a structure from stored blocks, and rebuild it exactly from its provenance |
| M6 | build real SBUs — paddlewheels, µ₃-oxo trimers — not just mononuclear nodes |
| M7 | attach trustworthy relative energies to any of it |
| M8 | compare two synthesis routes to the same product and say which is more viable, and why |

---

## 4. Test strategy (what "done" means, per layer)

**Fixture set — build it once in M2, reuse everywhere.** `tests/fixtures/`:
hand-written typed graphs for Cu paddlewheel, Fe₃-µ₃-oxo (both valence patterns), Zn₄O,
mononuclear Zn/BDC; molecules BTC, BDC, bipy, EDTA, anthrarufin; the cis/trans anthrarufin–Cu pair.

Four kinds of test, each with a job:

- **Invariance/property tests** — atom-order shuffles, insert idempotency, `construct` determinism.
  These guard the ground rules in §0 and are the ones that catch silent identity corruption.
- **Discrimination tests** — pairs that *must* differ, and pairs that *must not*, at each of L0/L1/L2/L3.
  Write them as one table so the identity boundaries are readable in one place.
- **Golden/regression tests** — the WL hashes of the fixture set (two-machine parity), and the
  Ni/Fe/BTC + solvation numbers from `docs/reports/` (M7).
- **Negative tests** — ambiguous spec raises instead of defaulting; a bad geometry is rejected by QC
  rather than silently stored; a version bump marks rows stale rather than re-labelling them.

Not aiming for coverage percentages. Aiming for: **every claim in the design doc's ledger has a test
that would fail if the claim stopped being true.**

---

## 5. Salvage ledger — what happens to `legacy/`

| Legacy | Verdict | Destination |
|---|---|---|
| `donor_perception.py` | **Port near-verbatim** — best-conditioned module you have | `sites/perception.py` (re-key to canonical idx) |
| `ebu_core._donor_placement_frame` | **Promote** — this is D13's frame | `sites/frames.py` |
| `ebu_core._best_azimuthal_rotation` | Replace with discrete torsion wells; keep as a clash-scoring fallback | `assembly/join.py` |
| `ebu_core._kabsch_align`, `_rot`, `_rotate_around_axis`, `_arbitrary_perpendicular` | Port as-is (pure math) | `geometry/_linalg.py` |
| `ebu_core._place_ligand_dg`, `_rigid_place_single_anchor` | Port as the **N=1 path** inside the general placer | `geometry/placer.py` |
| `ebu_core._site_vectors`, `CN_GEOMETRIES` | Port **and fix CN=5** (currently mapped to "planar") | `geometry/local.py` |
| `ebu_core._pair_chelate_compatible`, `_ligand_chelate_compatible`, `_interior_rotatable_bonds` | Port into the compatibility test | `assembly/join.py` |
| `geometry_qc.py` | Port + extend with `check_intercenter` | `geometry/qc.py` |
| `ebu_toolkit.formation_energies`, `solvation_sweep` | **Rewrite** against the new reference scheme; keep as numeric oracle | `energy/` + M7 regression tests |
| `xtb_energy.py`, `xtb_energy_metal.py`, `solvation_model.py` | Superseded; retain as M7 oracles only | — |
| `energy_model.py`, `dft_predict.py` (MACE) | Port behind the backend protocol | `energy/backends.py` |
| `ebu_tools_v2.PaddlewheelBuilder` | **Seed for M6 templates** | `geometry/templates.py` |
| `ebu_viz.py` | Optional, M9 | `viz.py` |
| `ebu_tools.py`, `ebu_tools_v2.py`, `ebu_toolkit.py`, all `.ipynb` | **Archive only — do not port** | stay in `legacy/` |
| `autonomous_enumerate_sbus` | **Do not port** — superseded by the choice-vector branch tree | `assembly/construct.py` |

---

## 6. Decision gates (the open checkpoints, placed where code forces them)

| Gate | Milestone | Forced by | What you need in hand to decide |
|---|---|---|---|
| **C8** — descriptor sourcing/provenance | **M1**, before the first table row | you can't type a pKa without deciding where it came from | a shortlist of sources (curated vs. dataset-pulled) and how the version gets pinned |
| **C5** — activation-ease floor | **M4**, before `site_state` accumulates | `site_state` needs an ease field for `raw_construct` structures | the M1 table's coverage; the list of in-pocket donors that self-flag provisional |
| **C2** — L3 thresholds (θ_geom, energy window) | **M5**, at the first conformer generation | clustering can't run without numbers | the pairwise core-RMSD distribution over the M5 fixture set — calibrate, don't guess |
| **C6** — barrier proxy | **M8**, before the first path score | `PathScore.max_barrier` needs a definition | whether the paddlewheel A-vs-B ordering is stable under the cheap proxies alone |
| **C7** — partner dependence | **M8**, alongside C6 | `ease(site, partner)` is called at query time | how many (site, partner) pairs you actually intend to screen — the factorization only pays off if that number is large |

Rule for all five: when a gate is called, **write the resolution into the design doc's Decision
Ledger with a new D-number and a changelog line**, then code against it. A gate resolved only in
your head is how the two documents drift apart.

---

## 7. Risk register and the escape hatches

| Risk | Signal it's happening | Escape hatch |
|---|---|---|
| **M6 placer overruns** (the known headline cost) | two weeks in and the paddlewheel still won't converge to sane M–M distances | switch to `geometry/templates.py` (grafting onto stored reference nodes). Decide by a pre-set date, not by mood — write the date down when M6 starts. |
| **Two machines disagree about identity** | a fixture hash differs laptop vs. workstation | already mitigated by §0.7 (WL primary, nauty as verifier). The golden-hash test in M2 is what catches it — do not skip it. |
| **L3 conformer explosion** | thousands of near-identical rows per (L1, L2) | choice-vector dedup runs *before* geometric clustering; cap conformers per (L1,L2); `TORSION_FREE` sites never branch (that tag is the guard) |
| **xTB numbers can't carry route claims** | M7 regression reproduces rankings but absolute ΔG look implausible | keep M8 claims *relative and within-metal*; the reaction-balanced reference scheme is the gate on any quantitative statement |
| **Identity retrofit pressure** | a temptation in M5/M6 to "just add a flag" to L1 | the M2 discrimination table is the contract; changing it means a version bump and a re-hash of the corpus, which is exactly the cost that should make you think twice |
| **Scope creep into periodic frameworks** | topology/net questions start appearing in tickets | out of scope for v1 (design §11) — record them as future work; MOFid/MOFkey is the reference oracle when you get there |

---

## 8. What to do first, concretely

1. **M0**, in one sitting: rails, `versions.py`, `make check`, `legacy/README.md` salvage verdicts.
2. **C8 call**, then **M1** — a low-effort session that unblocks three later gates.
3. **M2 fixtures before M2 code**: hand-write the paddlewheel and Fe₃-oxo typed graphs and the
   discrimination table as failing tests. Then make them pass. This is the single highest-leverage
   move in the plan — it forces the D12 polynuclear commitments into the schema at the point where
   they're free, and it gives M6 a ground-truth target that already exists.

---

## 9. Changelog

- *(rev 19)* **M1 descriptor layer — and no run mode had ever relaxed anything.**

  **The bug first, because it invalidates a run's label.**  A 500-structure `ml_go` build
  left the GPU idle.  Not a device problem: `relax_geometry` and `single_point` are called
  from **nowhere outside `mofsbu/energy/`**.  `plan()` checked whether the backend was
  importable and then queued ordinary `place` tasks; `execute()` never branches on
  `run_mode` and only ever writes `fidelity=RAW`.  The registry says it plainly — 281 RAW
  geometries, 17 FF, **none at ML**, no `mace` row in `methods` — while `runs` records the
  spec as `ml_go` and the status as `done`.  Rev 16 shipped the pre-flight without the
  executor, which is worse than shipping neither: the check made the option look wired.
  * `energy.relax.RELAXATION_IS_EXECUTED = False` is the one flag, and `mode_status()`
    now reports **two independent facts** — `backend_available` (installed on this
    machine) and `wired` (the pipeline calls it).  Collapsing them is exactly how "MACE
    imports" reached the user as "MACE will run".  `available` is the AND of the two.
  * `plan()` refuses an unwired mode with `NotBuiltYet`, not `EnergyBackendUnavailable`:
    a missing body and a missing install are different failures, and installing MACE does
    not fix this one.
  * Latent behind it: `MACEBackend` defaults to `device="cpu"`, so even once wired the GPU
    stays idle until a device is chosen.  Both are M7b, sequenced next as a separate
    `relax` task kind so a failed relax retries without rebuilding.
  * The 281 stored geometries are **not corrupt** — they are correctly labelled RAW with
    the placer as their method.  What was wrong was the label on the run, not the data.

  **M1, with C8 resolved: the two halves are sourced differently on purpose.**
  * `data/reference/donor_descriptors.tsv` — 36 rows, **hand-curated**, one per donor type
    `sites.perception` can emit.  Aqueous pKa is context-dependent and no package covers
    these classes, so each row names its class of reference and rows that are judgement
    say `estimate` rather than borrowing a citation they do not have.
  * `data/reference/metal_descriptors.tsv` — 35 ions, **generated** by
    `scripts/build_metal_descriptors.py` from `mendeleev` (Shannon radii by charge, CN and
    spin state), with HSAB class, preferred CN and water-exchange lability curated inside
    that generator because no package carries them.  The package version lands in every
    row, so regenerating against a different `mendeleev` is a visible diff (ground rule 6).
    `mendeleev` is a dev dependency; the TSV is what ships.
  * **`pka` runs in one direction only**: the pKa of the acid whose deprotonation exposes
    this donor.  A neutral donor needing no activation (amine, pyridyl, ether) has an
    EMPTY pka — a meaningful answer, not a missing value.  Storing basicity in the same
    column would have put two quantities under one name; a test pins the ordering
    sulfonate < carboxylate < phenolate < alkoxide, which only holds if the sign is
    consistent throughout.
  * `live_dof` is **not** in the TSV despite being a schema column: it is defined by
    `sites.frames.live_dof()` and filled in at load, so the torsion model keeps one home.
  * `donor()` and `metal()` **raise on an unknown key**.  A default row would make "never
    characterised" indistinguishable from a measurement, inside a screening loop.
  * `descriptors/ease.py` is a declared stub naming the checkpoint that blocks it — C5 for
    `activation_ease`, C7 for `hsab_match`.  Both stay open: the tables are facts, the
    ease model is a policy, and shipping the policy unratified is how a leaning becomes a
    decision by accident.
  * `tests/test_descriptors.py` gates coverage **both ways** — every perceivable donor
    type has a row, and no row exists for a type nothing emits.  The perceivable set is
    read out of `perception.py` rather than hand-listed, so it cannot drift.

  268 passed, 1 skipped.

- *(rev 18)* **Every halide co-ligand was rejected by construction, and the message said
  "2 clash(es)".**  Trying iodide in place of the water that fills leftover coordination
  sites produced a run of rejections with no way to tell why.  The cause was one line:
  `placer.D_ML` held a distance per METAL and applied it to every donor, so a Zn(II)
  centre placed an iodide at 2.00 A.  Zn-I is 2.60.  The clash check then compared the
  short bond against van der Waals radii — iodine's 1.98 A against oxygen's 1.52 — and
  the two errors compounded.  Water passed the same test by 0.04 A, which is how a
  systematically wrong table survived: it was only ever asked about oxygen.

  * **`geometry/distances.py`** is the new node.  `d(M, D) = base(M) + delta(element)`,
    where `base` IS the old `D_ML` table (so every O/N-donor geometry ever built is
    bit-identical) and `delta` is a per-element offset **calibrated against typical
    divalent bond lengths, not derived from radii** — radii were tried first and are
    wrong in a consistent direction, because a dative M-O bond is ~0.12 A longer than
    r_cov(M)+r_cov(O) while M-I is about right, so the residual is not a constant.
    Unknown pairs fall back to covalent radii and **carry `source="covalent-radii"`**,
    which travels into the choice vector, the QC notes and the inspector: a QC verdict
    on an uncalibrated pair is evidence about this table, not about the chemistry.
  * The donor's **element is read off the molecule**, never parsed from the donor-type
    name.  `hydrohalide_X` and `halide_I` name the same iodine, and a donor type added
    later must not silently inherit oxygen's bond length.
  * The placer now assigns vertices on **unit** vectors (assignment scores angles, which
    are scale-free) and scales each one by its own donor's distance.  The chelate branch
    takes its bite angle from the law of cosines instead of a symmetric arcsine, so a
    mixed N,O or S,O pocket gets both bonds right; with d1 == d2 it reduces exactly to
    the old expression.  `d_ml=` still overrides everything, which is what the tests
    that pin an exact bond length use.
  * Result: THQ + iodide on Fe(III) builds clean at CN 4 and 6, as `[I-]` and as neutral
    HI.  Cl, Br the same.  *This is not a claim that Zn/I is good chemistry* — it is that
    the refusal now has to come from the chemistry rather than from a table that had
    never heard of iodine.

- *(rev 18)* **The run inspector: a refusal has to say what it refused.**  The console was
  the only place a run's reasoning appeared, and a console that refreshes a progress line
  is not somewhere you can read a rejection.  Everything below is stored on the task row,
  so it is still readable tomorrow.

  * **Structured outcomes.**  `tasks` gains `error_code`, `detail_json`,
    `structure_created`, `geometry_created`.  The message is for reading; the **code** is
    for grouping (200 rejections with one cause should be one line, and free text cannot
    be grouped — "closest 1.40 A" and "closest 1.41 A" are one finding and two strings);
    the **detail** is for diagnosis without re-running anything.
  * **QC reports carry their evidence.**  A clash now records both atoms, their elements,
    which ligand each came from, the measured distance, the limit it was measured
    against, and the overlap; a bad M-L bond records its target distance and where that
    number came from.  `QCReport.code` and `elements_involved()` are what turn a pile of
    refusals into "×212 qc_clash, elements Fe/I".
  * **Regenerated vs skipped, at a glance.**  `structure_created=False` is the D2
    idempotency case — the task ran and wrote nothing because the registry already had
    that identity.  That is a *success*, and it was indistinguishable from building
    something new, so a run of pure duplicates looked exactly like a run of discoveries.
    `attempts > 1` marks a task that had to be re-run.  Both are summary cards and
    per-task badges.
  * **Embed provenance.**  `embed_with_report` records the ETKDG random-coordinates
    retry, whether MMFF or **UFF** parameterised the relaxation, and the exception that
    used to be swallowed by a bare `except: pass`.  A UFF-relaxed conformer is not
    comparable to an MMFF one and a QC clash after a retried embed is as likely to be
    about the embed as about the chemistry — neither was visible.
  * **`co_ligand` with no perceivable donor is an ANSWER.**  This was `perceive(co)[0]`
    and an `IndexError` reported as a machinery failure, for what is a plain statement
    about the co-ligand.  It is now `co_ligand_no_donor` with a hint naming the fix.
  * **`/runs`** renders it: summary cards, causes grouped with counts (click one to filter
    the task list), the planner's never-queued skips beside the workers' rejections, and
    per-task rows saying what was attempted, what came back, and which flags fired.  It
    reads through a read-only connection — the builder writes specs and queues, the
    inspector does not write at all.
  * Gates: `tests/test_distances.py`, `tests/test_run_inspector.py`.

- *(rev 18)* **A latent migration bug, surfaced by adding four columns at once.**
  `_reconcile_columns` stamped every added column with the same migration version
  (`SCHEMA_VERSION * 1000 + len(added)`), which is a UNIQUE violation on
  `migrations.version` as soon as one revision adds more than one column — and the whole
  migration then rolls back.  No revision ever had, so it sat there until rev 18 added
  four to `tasks`; the failure mode would have been every existing registry refusing to
  open after a `git pull`.  One row and one version per column now, verified against a
  real 125-task `runs.db` and gated by `test_several_columns_can_be_added_in_one_revision`.

- *(rev 17)* **Labile-pattern audit: four of twelve new SMARTS could never match.**
  Twelve donor patterns were added to `sites/perception.py` (hydroxamate, oxime, imide,
  acyl sulfonamide, thiophosphate, selenol, peroxy, halide, and three that are gone
  again).  Probing each one against a molecule it must recognise found that **four never
  fired and one fired on the wrong atom**, and that nothing anywhere would have said so.

  The cause is a convention the module documents and the patterns did not follow: the
  donor is the **last atom of the match**.  Spelled the obvious way, an imide
  (`[CX3](=[OX1])[NX3H1][CX3]=[OX1]`) ends on a carbonyl oxygen that carries no hydrogen,
  so the match is made and then silently discarded by the `h_idx` lookup.  Four patterns
  failed exactly this way; `peroxy_O` as `[OX2H1][OX2][#6]` ended on **carbon** and
  duly offered to deprotonate the methyl of CH3OOH to a carbanion.
  * Fixed by putting the donor last, with `$(...)` recursive SMARTS where the environment
    has to be described around it.  `nhydroxy_N` renamed `nhydroxy_O` (its donor is the
    hydroxyl O, not the N).  `acyl_sulfonamide_N` moved **above** `sulfonamide_N`, which
    was claiming the same nitrogen first.  `phosphinate_O` dropped as redundant —
    `phosphonate_O` already matches it — and `dicarbonyl_CH` / `nitro_CH` dropped because
    carbon is not a donor type this package can place; their real donor is the
    delocalised enolate/nitronate oxygen, which `enol_O` already covers.
  * **Halides.**  A halide almost never belongs in a *labile* list: as a ligand it is
    already X- (a whole ligand, i.e. a `co_ligand` such as `[Cl-]`), and on carbon it is
    a leaving group.  The one protic case is free HX, so the pattern is named
    `hydrohalide_X` for what it does.  It exposed a real gap underneath:
    `_classify_anionic` handled O, S and N and returned `None` for everything else, so
    the chloride produced by activating HCl **perceived no donor at all on recall** —
    activation destroying the site it creates, which is the opposite of perceive-once
    (D5).  A `halide_X` branch closes it.
  * `tests/test_donor_patterns.py` is the gate: every pattern must have a probe molecule,
    must fire on it, and must land on the element its own name promises.  A SMARTS added
    without a probe fails the suite, which is the check that was missing.

- *(rev 16)* **M7 first half: energy backends, and a reference scheme that refuses bad
  subtractions.**  `energy/backends.py` puts xTB (tblite), MACE and a test double behind
  one protocol; `energy/relax.py` gets the body its signature has been waiting for since
  M3; `energy/reference.py` is the new part.

  **The thing that was actually wrong is not what §11 said it was.**  The risk register
  described the archived formation energies as "gas-phase GFN2-xTB on isolated highly
  charged anions", and the fix was assumed to be a *balanced* equation.  Writing the
  balance checker disproved that: the legacy scheme

      E_form = E(EBU, q) - E(M^q+) - SUM E(free ligand anion)

  **balances**.  Atoms and charge both.  `docs/reports/ni_btc_report.md` even calls it
  "charge-conserving", correctly.  A balance check alone would have passed the equation
  whose Ni/EDTA verdict later reversed from "-8.5 eV, a deep trap" to "-0.3 eV, near
  parity" the moment a medium was added.  Balance is necessary and **not sufficient**,
  and a scheme that only checked it would have been a fix in name.

  So there are two gates, not one.  The second asks whether the two sides are comparable
  chemistry, and fails the legacy equation on three counts it can measure off the typed
  graph: a **bare ion** (a lone charged metal atom, the `E(M^q+)` term), a **naked
  polyanion** (a metal-free species with |charge| > 1, the `E(BTC^3-)` term), and
  **coordination change** — the count of metal-donor dative bonds must be equal on both
  sides.  That last one is the sharp instrument: it is what separates a ligand-exchange
  equation, where the metal is coordinated throughout and the errors cancel, from a
  formation-from-free-ions equation, where six bonds appear out of nothing.  A test
  asserts the archived equation balances *and* is still refused, so the distinction
  cannot quietly collapse back into a single check.
  * `strict=False` reproduces a legacy number on purpose; the result carries
    `isodesmic=False` and `store_reaction_energy` refuses it without `force=True`,
    stamping `NOT-ISODESMIC:bare_ion,...` into the note when forced.  The caveat travels
    with the number instead of living in a report nobody re-reads.
  * Mixed levels of theory in one equation are refused; a **charge-blind** backend (MACE
    sees elements and positions, not formal charge) is refused on any charged equation,
    which is a property of the `methods` row, not a docstring.  A species with no energy
    is **named**, never skipped — dropping a term from a balanced equation is how you get
    a confidently wrong number.
  * The equation's fidelity is the **weakest** rung in it.  An xTB product minus a
    raw-construct reagent is a raw-construct number.
  * **The spin convention is code now.**  The archived runs used minimal spin for Ni(II)
    and a fixed sextet for Fe(III), passed in through `METAL_MULT`.  Both come out of
    `high_spin_multiplicity` / `minimal_multiplicity`, and which one was used lands in
    the stored `MethodSpec`.  `check_spin` refuses a multiplicity the electron count
    cannot reach — xTB will happily return a number for an even-electron doublet.
  * `MethodSpec` moved to `_types` beside `Fidelity` (one definition, no energy->registry
    dependency; re-exported so every existing import still works).  `xtb_go` joins
    `RUN_MODES` and `SPEC_VERSION` is 3 — no field changed, but a v3 spec must be refused
    by a v2 build with "newer than I understand" rather than a bare run_mode error.
  * `available_modes()` now asks the backends whether they import, so the UI's disabled
    options track the machine.  DFT stays a `NotBuiltYet` — there is no external code
    wired up, and an option that silently ran xTB instead is the failure ground rule 8
    exists for.  `EnergyBackendUnavailable` ("install it") is deliberately a different
    exception from `NotBuiltYet` ("it is not written"): the laptop and the workstation
    give different answers to the first and the same answer to the second.

  **Exit gate: half done, and the other half is blocked on missing data.**
  `scripts/regress_m7.py --refs` recomputes the four archived Fe(III) reference energies
  from `legacy/fe_btc_refs.json` — same species, same method — which checks the backend
  is the calculator that produced the archived numbers.  The ranking regression the plan
  asks for cannot run: `ni_btc_report.md` cites `ni_btc_outputs/` and
  `ni_btc_xtb_results.csv`, neither is in the repo, and `legacy/ni_btc_results.csv` has an
  **empty `formation_eV` column**.  The archived ranking currently exists only as a table
  in a markdown report.  Stage 2 is a declared stub that says exactly this rather than
  regressing against numbers that are not there.

- *(rev 15)* **Migrations that actually migrate.**  Submitting a run against a registry
  created one revision earlier failed with `no such column: runs.diagnostics_json`:
  `migrate()` only ran `CREATE TABLE IF NOT EXISTS`, which does nothing to an existing
  table, so no column added after a database was first written had ever reached it.  It
  now builds the target schema in an in-memory database (letting SQLite parse its own
  schema rather than parsing SQL by hand), adds missing columns via `ALTER TABLE`,
  recreates views — `CREATE VIEW IF NOT EXISTS` silently keeps a stale definition, which
  had left `v_structures` out of date on older files — and records each change in the
  `migrations` table.  Destructive differences are reported and left alone.  Verified
  against the user's own `data/runs.db`, which upgraded in place and then ran their
  anthrarufin/Cu spec to 11 Cu complexes.  **Ground rule 7 added** to say so, since this
  was the second schema change in two revisions to silently break existing data (spec
  files were the first).
- *(rev 14)* **A run that builds nothing now says why**, plus a spec browser in the UI.
  Reported case: anthrarufin + Cu perceived its sites and produced three ligand records
  and no complexes at all.  The spec was correct.  The planner skipped every coordination
  candidate because `co_ligand` was null while the coordination list asked for CN 4 and 6 —
  a bidentate chelate uses two sites and there was nothing to fill the rest — and it did so
  **silently**, so "your filter matched nothing" and "it worked" were indistinguishable
  from outside.  `plan()` now records a diagnostic for every skipped candidate with a
  reason and a concrete hint, stored on `runs.diagnostics_json`, returned by the submit
  endpoint and shown in the page.
  * `allow_unsaturated` added (default off): builds the coordinatively unsaturated product
    instead of skipping.  Off by default because quietly building something smaller than
    asked for is the failure this project keeps rediscovering.
  * Coordination geometry now follows the CN **actually** built, not the requested one —
    after collapsing to an unsaturated product the old code asked for tetrahedral with two
    sites, which the placer rightly refused.
  * A placer refusal is a **rejection, not a failure**.  "A bidentate ligand cannot span a
    linear two-coordinate centre" is chemistry answering; recording it as `failed` makes a
    run full of sound chemistry look like a run full of bugs.
  * **Specs are versioned.**  Adding `run_mode`/`allow_unsaturated` broke every spec
    already on disk, including the tracked reference one; `SPEC_VERSION` is 2 with a
    migration from 1, unknown fields are refused rather than ignored, and a FUTURE version
    is refused rather than guessed at.
  * UI: `/api/specs` lists seed specs from `data/reference/` and `data/specs/` with their
    molecules, metals and degree; one click loads a spec back into the form; non-spec JSON
    (the molecule library) is filtered out rather than listed as broken.
- *(rev 13)* **Spec builder UI, and the first real concurrency bug.**
  `mofsbu/ui/builder.py` + `static/builder.html` at `/builder`: SMILES entry with a live
  RDKit 2D depiction and a read-out of what the pipeline sees (donors, chelate pockets with
  their structural descriptors, how many distinct activation states), a named molecule
  library persisted to `data/reference/molecule_library.json` (tracked, so it travels), the
  run's molecule and metal lists, construction options, and run submission.  Runs execute on
  a background thread and progress is polled from the task table — the page can be closed
  and reopened without losing anything, because the state was never in the page.

  * **Metal centres are optional and first-class.**  A spec with no metals is a molecular
    construction (COF, cage), not a degenerate case; `BuildSpec.__post_init__` requires
    molecules and never metals.  Metal-free runs currently produce activation states and
    sites; joining molecule to molecule is still M5, the same stub that gates degree >= 2.
  * **Unimplemented options are labelled by the BACKEND.**  `/api/capabilities` reports
    which run modes can execute (`energy.relax.available_modes`) and the page renders from
    it, so a disabled option cannot drift out of date against a hard-coded string in the
    markup.  ML and DFT optimisation are shown, disabled, tagged "not implemented — M7",
    and refused server-side at plan time with the reason.
  * The viewer stays read-only.  The builder writes spec files, library entries, and
    run/task rows — never a structure, geometry or label.  Those still only arrive through
    `registry.api` from a worker.

  **The bug the UI found:** submitting two runs at once failed about two times in five with
  `IntegrityError` on the identity constraint.  `put_structure` reads before it writes, so
  two workers both miss and both insert.  *Idempotency is not concurrency safety* — the
  earlier claim that idempotent writes made parallel workers safe was untested, and wrong.
  The UNIQUE constraint is what makes the race safe rather than silently duplicating, so
  the fix is to let it arbitrate and re-read, not to add locking; `put_geometry` and
  `method_id` had the same shape and got the same treatment.  Two deterministic tests now
  interleave two connections explicitly rather than relying on the flake.
- *(rev 12)* **Work is queued, not called** — the decoupling that makes MPI/CUDA a
  configuration change rather than a rewrite.  `mofsbu/spec.py` (a build is serialisable
  data with a content digest), `runs` + `tasks` tables with an atomic IMMEDIATE claim,
  `registry/jobs.py`, `mofsbu/runner.py` (plan / work split), `scripts/run_spec.py`, and a
  worked spec at `data/reference/spec_zn_thq.json` — **that file is the input interface:
  edit it and re-run**.  Nothing assumes planning and execution share a process, a machine
  or a moment, so an MPI launcher or GPU worker becomes another consumer of the task table.
  `rejected` is distinguished from `failed`: a candidate that fails QC is an answer, not a
  crash, and conflating them makes "produced nothing" indistinguishable from "broke".
  Honest limit recorded: SQLite's locking is fine for a workstation with local disk and not
  dependable on a shared cluster filesystem, so no signature in `jobs.py` mentions SQLite.

  **Ground rule 7 added** (agreed with the user, plus one clause): a missing function gets
  a stub with its real signature rather than reshaping the pipeline around what is currently
  implementable — *and the stub must fail loudly*, never return a plausible default, because
  a stub that quietly answers is indistinguishable from a working function until the answers
  are wrong.  Applied immediately: `assembly.join` (`compatible`, `join`, `grow`,
  `BuildingBlock.open_sites`) and `geometry.placer.place_multicentre` are settled interfaces
  that all raise `NotBuiltYet` naming their milestone.  `degree > 1` therefore refuses with
  an explanation instead of silently building a mononuclear complex.

  **Ground rule 8 added:** parallelism is opt-in and the laptop is the default.  Nothing
  spawns a process unless `MOFSBU_PROFILE=workstation` or `MOFSBU_WORKERS=N` is set; an
  unconfigured machine runs one in-process worker.  A test asserts that a run on an
  undeclared machine never touches `multiprocessing.Process`.
- *(rev 11)* **Reference cases became an artifact instead of a claim.**  The 13-molecule
  validation in rev 10 ran in a throwaway shell one-liner: nothing was written, nothing was
  re-runnable, and the only evidence it happened was a paragraph saying so.  Now
  `data/reference/ligand_cases.tsv` holds the cases (tracked — it is the one part of
  `data/` git keeps), `scripts/check_cases.py` runs them and exits non-zero on
  disagreement, and `tests/test_reference_cases.py` parametrises them into the suite so a
  case added there is a permanent regression guard.  **The TSV is the interface for new
  test cases**: add a row, no code changes.

  It earned its keep on the first run by failing two cases — both of which were MY
  expectations, not the code: oxalic acid has 3 protomers, not 4 (its two mono
  deprotonations are symmetry-equivalent), and BTC does chelate, through a single
  carboxylate's own two oxygens as a 4-ring.  Both rows are annotated EXPECTATION
  CORRECTED with the reasoning, so a later reader can see the expectation moved and why —
  rather than finding a number that silently matches whatever the code emits.

  Also noted: generated `.db` files and xyz folders written under `data/` had disappeared
  from disk between runs while `data/store/` survived, which is the OneDrive hazard flagged
  in rev 3.  Set `MOFSBU_DATA` outside the synced tree; the reference cases are unaffected
  because they are tracked rather than generated.
- *(rev 10)* **General driver (`scripts/build.py`) + five defects it exposed.**  Any SMILES
  can now be run through the whole pipeline; `--dry-run` reports protomers, donors and
  pockets without writing.  Enumeration is liberal and identity deduplicates; rejections
  are printed with reasons rather than dropped.  Pointing it at ordinary chelators found:

  * **`pocket_type` was another taxonomy** ("dioxolene", "enolate-carbonyl", ...).  Pockets
    are now described structurally — chelate ring size, how many donors are anionic, which
    donor types — and selected by predicate (`find_pocket(ring_size=5, n_anionic=2)`).
  * **Only 4- and 5-membered chelates were findable.**  Pockets were matched by "attachment
    atoms bonded or identical", which silently misses every 6-ring — including salicylate,
    whose 6-membered chelate is the reason anyone uses it.  Now found by shortest path.
  * **`binding_modes` gated which donors could chelate** — a per-donor-type table that
    denied bipyridine and ethylenediamine any pocket at all.  Whether two donors can
    chelate is a property of the PAIR, so the gate is gone.
  * **Rotatable-bond counting included the bonds at the donors.**  A donor lies on its own
    bond axis, so rotating that bond cannot move it; counting them made every hydroxyl look
    flexible and let a 7-ring "chelate" across a benzene pass as reachable.  Interior bonds
    only (this is what `legacy._interior_rotatable_bonds` did; the port had lost it).
  * **The site frame picked an arbitrary SIDE.**  In tilted mode the axis was rotated about
    a plane normal whose sign came from whichever substituent RDKit listed first, so
    adjacent ring donors could be handed axes 155 degrees apart and a good pocket looked
    impossible.  An sp2 donor has two real lone pairs: the side IS the torsion well (D13),
    so `site_frame(..., well=)` selects it and the pocket test searches both.

  **Feasibility is calibrated, not guessed.**  A rigid pocket must put a metal at a sane
  bite angle and M-D distance without standing inside the ligand.  The M-D cut (1.35) was
  measured: rigid pockets that must be found sit at 1.08-1.14, rigid pairs that must not at
  1.55-1.63.  The predecessor was a raw "predictions within 1.6 A" test that rejected
  doubly-deprotonated THQ's own chelate at 1.61 — a hundredth of an angstrom deciding
  whether the project's central pocket exists.  Verified against 13 reference molecules,
  positive and negative: catechol, salicylate, acac, bipy, en, oxalate, phenanthroline,
  biphenyl-2,2'-diol found; resorcinol, benzene-1,4-diol, phenol, benzene, ethanol not.
- *(rev 9)* **Configuration taxonomy deleted; descriptors derived from Aut(parent).**
  `classify_configuration` was a three-entry lookup on ring distance
  (`{1: ortho, 2: meta, 3: para}`) — a promise to enumerate chemistry in advance that had
  nothing to say about a five-membered ring, two sites on different rings, fac/mer, or any
  scaffold nobody wrote down.  It was also never load-bearing: the 16 -> 7 collapse came
  entirely from L1 hashes, and the names were painted on afterwards.

  * **Automorphisms promoted out of the canonical-labelling search** (`atom_orbits`,
    `automorphism_generators`, `automorphism_group`).  They were already being discovered
    and used for pruning, then discarded — the same move as promoting the site frame out
    of the placer (D13).  Symmetry can now be *asked about*, not just exploited internally.
  * Two activation selections are the same configuration exactly when a symmetry of the
    **parent** maps one onto the other.  The descriptor is the canonical orbit
    representative (`cfg(0,2)`), generated and unique within the family, plus a
    **separation profile** — sorted pairwise topological distances — as the invariant a
    human can read (`sep(3)` is the adjacent pair).  No chemistry is encoded anywhere.
  * Verified general: BTC's three carboxylates are mutually equivalent so every
    doubly-deprotonated selection is ONE configuration (a ring-distance table would have
    called them all "meta"); a four-membered ring, which has no ortho/meta/para at all,
    still gets correct distinct descriptors.
  * **Selection is structural too.**  The build no longer asks for "the ortho protomer";
    it asks for the q-2 protomer that *presents a dioxolene pocket*, and which descriptor
    that turns out to be is a result.  A test asserts exactly one of the three qualifies,
    and independently that it is the one with the smallest separation.
  * A test asserts the orbit partition and the L1 partition agree — two independent routes
    to the same answer, so a disagreement means one of them is wrong.

  **Where configurations belong.**  For THQ protomers, ortho/meta/para is an **L1**
  distinction and always was: the protons sit on different oxygens, so the graphs differ
  and L1 separates them exactly.  The descriptor explains an existing distinction; it does
  not create one.  `?L2` remains for the genuinely different case — same graph, different
  spatial arrangement (cis/trans at a metal centre, fac/mer, Delta/Lambda) — which is M5.
  When it lands it uses this same machinery: L2 is orbits of the coordination sphere under
  Aut(complex), resolved by geometry.  Familiar words stay available as display-only
  aliases keyed on a descriptor, in exactly the way fragments get names, and nothing in the
  pipeline depends on anyone having supplied one.
- *(rev 8)* **Labels became a projection; search became structural.**  `display_label` was
  composed from `TypedGraph.name`, a hand-typed string — so authored examples read
  perfectly and anything the enumerator generated collapsed to a bare formula that
  collides across isomers.  Worse, the viewer's free-text search matched on it, which
  makes a query depend on whether someone had named the thing.

  * `mofsbu/naming.py` composes a label from the graph alone: metal centres plus ligand
    fragments with counts and bridging, e.g. `Zn[H2O]2[THQ-2H(ortho)] q0`.  No author name
    is consulted anywhere.
  * `structure_fragments` (fragment L1 hash, formula, count, bridging) makes "which
    structures contain this ligand" an **indexed join on fragment identity**.  `find()`
    gains `fragment_l1` / `fragment_formula`; the viewer's `q` searches whole formula, L0,
    and fragment formula — never a label.
  * `fragment_aliases` names a **fragment once**, not each structure containing it; aliases
    are display-only, so adding or changing one moves no query.  `relabel_all` refreshes
    the cached labels.
  * `set_display_label` is **removed** — there is no longer any way to write a label by
    hand — and `display_label` is out of the sortable set.  `verify()` fails if a stored
    label or fragment decomposition disagrees with the graph.
  * Tests assert the invariant directly: a string that exists only in a label
    (`"paddlewheel"`) must return zero rows.

  **Placeholders are explicit and ugly on purpose.**  Every label carries `?L2` because
  configurational isomers are not distinguished yet (M5) — so ortho / meta / para siblings
  visibly share a label rather than appearing to be properly named.  `n_open_sites` stays
  NULL and a separate `n_perceived_donors` holds what was actually counted, because open
  vs occupied needs `site_state`, which does not exist yet.  Nothing incomplete is dressed
  up to look finished.
- *(rev 7)* **Viewer thread bug fixed.**  Every request 500'd under a real server while the
  whole test suite stayed green: FastAPI runs a synchronous `yield` dependency in one
  threadpool thread and the endpoint that consumes it in another, so the per-request SQLite
  connection crossed threads and sqlite3 refused it.  Connections are now opened with
  `check_same_thread=False` — safe here because one request owns a connection from open to
  close and it is read-only.  The regression test asserts the property directly (open on one
  thread, execute on another) rather than firing concurrent requests at TestClient, which
  reuses a single worker and passes either way; the useless version was written first and
  discarded after checking that reverting the fix did not make it fail.  **Any test claimed
  as a regression guard should be verified to fail without the fix.**
- *(rev 6)* **Activation states enumerated instead of chosen** (`sites/protomers.py`).  The
  first Zn/THQ build hard-coded "deprotonate one, deprotonate two" and picked hydroxyls by
  list index, which is not activation-at-viable-sites, it is a script author guessing.  THQ
  has four acidic hydroxyls: 16 subsets, which the L1 hash collapses to **seven distinct
  protomers** — nothing is told the molecule's symmetry, the identity layer discovers it.
  The dianion splits into **ortho / meta / para**, and that distinction is chemistry, not
  labelling: `sites/model.pocket_type` shows only the ortho dianion presents a **dioxolene**
  pocket (two anionic oxygens on adjacent carbons), so meta and para THQ(2-) carry the same
  formula and charge and cannot chelate at all.  Coordination is now built on the pocket
  asked for by name (`find_pocket(kind="dioxolene")`) rather than whichever came first;
  the singly-activated enolate-carbonyl complex is kept as a transient.  Deprotonation is
  recorded as a reaction edge with ammonia as the base, so the counterion's charge
  accounting is provenance rather than assertion.

  Known limitation: placement is rigid, so a chelate keeps the FREE ligand's bite angle —
  the dioxolene pocket comes out near 97 degrees where a bound catecholate sits nearer 85,
  because the free dianion's two O(-) repel.  Correct at `raw_construct` fidelity; it wants
  the relaxation of M7, and the geometry ladder is what records the difference.
- *(rev 5)* **Molecular input, sites, and the N=1 placer landed** on the Zn / THQ /
  ammonium test system.  `graph/from_mol.py` (SMILES -> typed graph, with D15 charge
  placement: localised charge stays on an ammonium N, delocalised charge goes to the graph),
  `sites/perception.py` (ported from legacy), `sites/frames.py` (D13 promoted out of the
  placer), `sites/model.py`, `geometry/embed.py` (ETKDG + FF), `geometry/qc.py` (ported),
  `geometry/placer.py` (frame-directed single-centre placement), `registry.put_sites`.
  Four real Zn complexes built, QC-clean, registered with sites — `scripts/build_zn_thq.py`.

  Four defects found and fixed on the way, three of them inherited:
  * legacy perception saw **nothing** on THQ — its ring is not aromatic so `phenolate_O`
    missed the hydroxyls, and C=O was never a donor at all.  Added `enol_O` and
    `carbonyl_O`; the pocket a metal actually binds is the enolate plus its neighbouring
    carbonyl.
  * `GetTotalNumHs()` counts only IMPLICIT hydrogens, so on hydrogens-explicit molecules
    every H-count test in the classifier read zero — water was not a donor and an N-H
    heteroaromatic would have passed as pyridyl.
  * perception rejected any charged atom, so **deprotonating a donor destroyed the site it
    created**; a structure recalled from the registry is normally already activated.
  * the placer handed chelates polyhedron vertices in table order, which on an octahedron
    is a **trans** pair no five-membered ring can span; and forcing a rigid chelate onto
    fixed vertices strained Zn-O from 1.85 to 2.39 A.  Vertices now choose the sector and
    the ligand's own donor-donor distance sets the bite angle — every Zn-O is 2.000 A.
- *(rev 4)* **Integrity checking added** (`registry/verify.py`, `scripts/verify_registry.py`,
  wired into `scripts/check.sh`).  Ground rule 1 said the registry is written only through
  `registry.api`; the demo seed wrote raw INSERTs and produced 40 rows describing impossible
  molecules — a paddlewheel with no M-M bond, `muN` on a two-metal node, 8 dative bonds recorded
  as 3, one geometry shared by four structures, and every `typed_graph_hash` pointing at a blob
  that was never stored.  Nothing caught it because nothing compared a column to its graph.
  `verify()` re-derives every derivable column from the stored typed graph and reports
  disagreements, separating **error** (contradicts its graph) from **stale** (older algorithm
  version — recompute, not corrupt).  Ten tests break something on purpose and assert it is
  reported, because a verifier that only says "fine" is worse than none.  Example structures moved
  into the package (`mofsbu/examples.py`) so tests, scripts and the demo registry share one source
  of real chemistry; the seed now builds through the API.  Placeholder geometry
  (`geometry/layout.py`, spring embedding) is labelled `spring-embed-3d` and synthetic energies
  `synthetic-demo`, so neither can be mistaken for a calculation in the viewer.  `depth` on
  reaction edges reclassified live-from-M3, not reserved-for-M8.
- *(rev 3)* **M3 core landed** — schema (14 live + reserved tables, `v_structures` read view),
  content-addressed `BlobStore`, `Registry` with migrations and algo-version bookkeeping, and the
  single write surface `put_structure`/`put_geometry`/`put_reaction`.  Identity, the canonical
  order and the derived query columns are persisted; duplicate inserts collapse to one row with
  two provenance edges.  **M3.5 viewer landed** in parallel (FastAPI + 3Dmol.js, read-only).
  Schema amendments from viewer feedback: `structure_metals` child table, `min_depth` on the
  view, a stated delimiter convention for set columns.  Added `mofsbu/config.py` (`MOFSBU_DATA`)
  so the registry can live outside a synced folder, and a journal-mode probe in `Registry` that
  degrades WAL -> TRUNCATE -> MEMORY instead of failing at commit with a bare I/O error.
  Legacy-corpus ingest is blocked on a charge/multiplicity decision — see M3.
- *(rev 2)* Added §2.5 **interface architecture**: inputs and outputs split, inputs move to a
  declarative YAML spec (making every UI a thin editor of it, per D13), outputs get a read-only
  FastAPI + 3Dmol.js viewer.  Inserted **M3.5** before M4 so the site-frame work can be debugged
  visually.  M3 gains denormalised query columns for all four filter axes (composition,
  connectivity/binding mode, energy/fidelity, provenance) and a derived non-identifying
  `display_label`.  Recorded that a desktop GUI and any run-launching UI are deferred, and named
  the seams they would attach to.
- *(rev 1)* Initial implementation plan: ground rules; package map with API signatures; milestone
  spine M0–M9 with exit gates; functionality ledger; test strategy; salvage ledger for `legacy/`;
  open checkpoints C2/C5/C6/C7/C8 placed as decision gates; risk register with declared escape
  hatches (notably the M6 template fallback).
