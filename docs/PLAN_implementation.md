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
| **Done** | M0 rails · M1 descriptors · M2 graph + identity · M3 registry · M3.5 viewer · M4 sites |
| **Partly done** | **M7** — backends, `MethodSpec`, reference scheme and the relax runner all shipped; the regression exit gate has a harness and no recorded result |
| **Next** | **M5** recursive assembly — the critical path, and the first milestone that manufactures conformers |
| **Then** | M6 placer (headline cost) → M8 pathways (the actual contribution) · M9 folded in opportunistically |

```
[M0 · M1 · M2 · M3 · M3.5 · M4 — done] ─► M5 assembly ─► M6 placer ─┐
                                                                    ├─► M8 pathways
                       M7 energy — built, exit gate never run ──────┘
```

**Open decision gates:** C2 (forced by M5), C6 and C7 (forced by M8). See §3. Nothing is
pre-resolved here; the design doc's leanings stand.

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

# templates.py  — plan-B / accelerator for M6 (see §4 risk)
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

**Plan B, declared in advance (see §4):** if constrained multi-center placement stalls, fall back to
`geometry/templates.py` — place from a stored reference node geometry and graft ligands onto it
(`ebu_tools_v2.PaddlewheelBuilder` is the seed). Honest, much cheaper, and it keeps M8 reachable.
Templates are a *geometry source*, not a second identity path — they enter through the same API.

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
| **C2** — L3 thresholds (θ_geom, energy window) | **M5**, at the first conformer generation | clustering can't run without numbers | the pairwise core-RMSD distribution over the M5 fixture set — calibrate, don't guess |
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
| **M6 placer overruns** (the known headline cost) | two weeks in and the paddlewheel still won't converge to sane M–M distances | switch to `geometry/templates.py` (grafting onto stored reference nodes). Decide by a pre-set date, not by mood — write the date down when M6 starts. |
| **Two machines disagree about identity** | a fixture hash differs laptop vs. workstation | already mitigated by §0.7 (WL primary, nauty as verifier). The golden-hash test in M2 is what catches it — do not skip it. |
| **L3 conformer explosion** | thousands of near-identical rows per (L1, L2) | choice-vector dedup runs *before* geometric clustering; cap conformers per (L1,L2); `TORSION_FREE` sites never branch (that tag is the guard) |
| **xTB numbers can't carry route claims** | M7 regression reproduces rankings but absolute ΔG look implausible | keep M8 claims *relative and within-metal*; the reaction-balanced reference scheme is the gate on any quantitative statement |
| **Identity retrofit pressure** | a temptation in M5/M6 to "just add a flag" to L1 | the M2 discrimination table is the contract; changing it means a version bump and a re-hash of the corpus, which is exactly the cost that should make you think twice |
| **Scope creep into periodic frameworks** | topology/net questions start appearing in tickets | out of scope for v1 (design §11) — record them as future work; MOFid/MOFkey is the reference oracle when you get there |

---

## 5. What the tool can do after each remaining milestone

| After | You can… |
|---|---|
| *(today)* | store, dedupe, query and **see** every structure built; ask a stored structure what sites it has, which are open, and which are worth spending QM on |
| M5 | build a structure from stored blocks, and rebuild it exactly from its provenance |
| M6 | build real SBUs — paddlewheels, µ₃-oxo trimers — not just mononuclear nodes |
| M7 *(gate)* | trust the relative energies attached to any of it, because they reproduce the archived results |
| M8 | compare two synthesis routes to the same product and say which is more viable, and why |

---

| M8 | compare two synthesis routes to the same product and say which is more viable, and why |

---

## 6. Test strategy (what "done" means, per layer)

**Fixture set — build it once in M2, reuse everywhere.** `tests/fixtures/`:
hand-written typed graphs for Cu paddlewheel, Fe₃-µ₃-oxo (both valence patterns), Zn₄O,
mononuclear Zn/BDC; molecules BTC, BDC, bipy, EDTA, anthrarufin; the cis/trans anthrarufin–Cu pair.

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
