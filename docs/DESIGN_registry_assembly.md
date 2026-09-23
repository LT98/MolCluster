# Design: Recursive Assembly + Persistent Structure Registry

**Status:** living design doc — decisions and open checkpoints tracked at the bottom.
**Scope:** the realignment of the EBU screening tool toward the north-star goal.
**How to use:** edit freely. The *Decision Ledger* (§7) is the source of truth for what's
locked; the *Open Checkpoints* (§8) are what still needs a call. When you change a decision,
add a revision to [`archive/DESIGN_history.md`](archive/DESIGN_history.md) saying why — a
decision that changed without a recorded argument is one nobody can re-examine later.

---

## 1. North-star goal

A method to **predict and guide synthesis design** by building structures from molecules,
where:

- molecules (and structures) can be selectively **activated** at viable sites;
- structures are **infinitely expandable** — a built structure re-exposes its own remaining
  sites and can accept the next building block, including other structures;
- the system can **detail the construction** of a structure from existing molecules/structures,
  **evaluate the viability of the synthesis route** (reagents → product, with conditions and
  intermediates), and **store all of this** for recall and further construction.

The intellectual contribution is the **route-design layer**: not just ranking which product is
most stable, but evaluating and comparing *pathways* (e.g. bimetal-nucleus-first vs. sequential
metal addition for a Cu paddlewheel) and the effect of intermediate chelating/modulating agents.

---

## 2. Current state (honest assessment)

The existing tool (`ebu_core.py`, `ebu_toolkit.py`, `donor_perception.py`, `geometry_qc.py`) is a
**mononuclear node screening pipeline**:

- `LigandBuilder._activate` — generic SMARTS donor perception + deprotonation. **This is the
  "activate at viable sites" primitive and the most goal-aligned piece that already exists.**
- `EBUBuilder` + `GeometryPlacer` — place N ligands around **one** metal to a valid CN geometry
  (planar / tetrahedral / octahedral), with real QC (clash / bond / ring-planarity).
- `autonomous_enumerate_sbus` — search copy-counts × dentations for that one metal.
- `formation_energies` / `solvation_sweep` — GFN2-xTB ranking + the nucleation-drive /
  sequestration-margin framing. Ni/Fe/BTC + EDTA/EDDA study reports are solid.

**Gaps vs. the goal:** no recursive/expandable structure object; no polynuclear (no M–M, no
µ-bridges — can't build a single canonical SBU); no persistent registry (storage is `.xyz` +
CSV + an in-memory dict); no route/reaction layer; endpoint thermodynamics only; gas-phase
charged-ion energetics are screening-grade at best. Code hygiene: multiple parallel script
versions, 2 git commits — consolidate to one lineage before building the new layer.

---

## 3. Core architectural principle

### 3.1 Recursive `BuildingBlock`

Reframe the central object from "an EBU" to a **recursive BuildingBlock** with one interface,
whether it is a molecule, a metal, or a large assembled fragment:

```
BuildingBlock = { geometry, open_sites[], net_charge, provenance }
open_site     = { type, frame, live_dof_tag, binding_modes[], occupied? }   # frame, not a lone vector — see §6.3
```

**Assembly** = match an open site on A with a compatible open site on B → join → return a new
BuildingBlock that recomputes its remaining open sites. Molecule+metal, metal+metal-via-µ-oxo,
and structure+structure all become the same operation. This is what makes "infinite
expandability" and "structures built from structures" real.

### 3.2 Separate three things that the old filename conflated

The filename length problem was a symptom, not the disease. Split:

1. **Identifier** — a semantic, canonical key answering "are these the same structure?"
   (build-order-invariant, isomer-aware). Stored *inside* the record.
2. **Address** — where the bytes live on disk. A content hash (fixed length). No meaning encoded.
3. **Provenance** — how the structure was made. A DAG **edge**, not part of identity or address.

Consequence: identity lives on the node, sequence lives on the edges. A structure reachable by
two build sequences is **one node with two incoming edges**.

---

## 4. Identity model — a layered key

Do not look for one magic string. Build a **hierarchical key** and match at whatever resolution a
query needs:

```
BLOCK_ID = L0_composition / L1_graph_hash / L2_isomer_tag / L3_conformer_id
```

| Level | Meaning | Invariances / discriminations |
|---|---|---|
| **L0** | composition + charge + spin (`Fe1_C9O6_q-3_s6`) | coarse bucket; separates protomers/charge states |
| **L1** | canonical hash of the **typed molecular graph** | build-order invariant; separates coordination isomers with different connectivity |
| **L2** | isomer / configuration tag (cis/trans, fac/mer, Δ/Λ) | separates same-graph spatial isomers L1 can't see |
| **L3** | conformer / distinct 3D minimum | separates "same route, substantially different structure" that is not a discrete isomer |

**L2 is confirmed in scope, and the discriminator is downstream-assembly-relevance, NOT the
energy gap.** Two configurational isomers can be near-degenerate and still be distinct, persistent,
separately-addressable species because a **barrier** separates them (you can't interconvert without
breaking/reforming bonds or a large hindered rotation). Operational test for L2:

> Does the arrangement change the geometry/vector of the **open sites the next block will see**,
> and is it barrier-separated into a persistent species? If yes → L2.

Motivating case: **anthrarufin–Cu**. Anthrarufin (1,5-dihydroxyanthraquinone) presents a peri-chelate
pocket at each end (1-OH/9-C=O and 5-OH/10-C=O); the cis/trans disposition of what each Cu presents
onward decides whether the extended assembly runs linear or kinks — configurational, barrier-separated,
and decisive for later building-block assembly. An energy-based cut would wrongly merge near-degenerate
cis/trans pairs, which is exactly why the criterion is downstream-relevance, not ΔE.

**Why not SMILES / InChIKey:** SMILES isn't canonical across toolkit versions and mishandles
dative bonds. InChIKey *is* fixed-length (solves filename length) but standard InChI **disconnects
metal–ligand bonds by default**, so it silently merges the very coordination isomers we need to
tell apart. Both fail on the identity axis that matters here.

### 4.1 Canonicalizing the metal graph (the L1 hash)

Build an **explicit typed graph** — nodes labeled `(element, formal_charge, oxidation_state,
spin_class)`, edges typed `covalent | dative(→M) | metal–metal`. **Do not route through
SMILES/InChI.**

Bridging is **not** an edge type: a donor carrying dative edges to two distinct metals *is* µ2,
and labelling it again would give one chemistry two encodings that hash differently. µ2/µ3 are
derived on demand (`TypedGraph.bridge_class`) — that is D14.

The key and its two supporting structures (D16):

1. **Canonical certificate → sha256** (`graph.canon.certificate_digest`) — individualisation-
   refinement with branch-and-bound and automorphism pruning. **This is L1**, and it also fixes
   the atom ordering used as the coordinate system for site annotations (see §6).
2. **Weisfeiler–Lehman hash** (`wl_hash`, surfaced as `identity.wl_index`) — a fast **bucket
   index**, not a key. 1-WL cannot separate µ2-bridging from chelating, which is precisely the
   distinction this project exists to draw.
3. **VF2 isomorphism** (`networkx`) — collision resolver, run only when two records share an L1.

Pure Python throughout, so both machines agree. **nauty/Traces is not used** — its canonical form
differs from ours and would produce different keys; the original plan to use it as a verifier was
dropped with D16.

Polynuclear (M–M, µ-bridges) is handled natively by this typed graph — a key reason to prefer it
over SMILES, which is exactly where those break.

### 4.2 "Substantial difference" is a policy (L3)

You cannot hash floating-point coordinates. Assign L3 at insert time by **clustering** against
existing entries with the same L1/L2. **Whole-molecule RMSD is a bad sole trigger** — a distant
floppy tail gives huge RMSD but is chemically irrelevant, while a small twist that parks a site
open/blocked barely moves RMSD. Use a **compound rule** tied to the same principle as L2:

> New L3 if ( RMSD on the **rigid core + coordination sphere only** > θ_geom ) **OR** the pose
> flips any **open-site accessibility flag** — gated by ( energy within a sane window of the
> representative, so a bad geometry isn't mistaken for a real minimum ).

Optionally use a rotation/translation/permutation-invariant **shape fingerprint** (USR/USRCAT, or a
sorted internal-distance histogram) as a fast pre-filter before RMSD. The open-site-accessibility
trigger is the important half: it makes the L2 and L3 boundaries key off the *same* question
("does this change the open-site picture for the next step?"), at different barrier regimes.

**What counts as an L3 conformer — worked examples (calibrates C2):**

- *Trivial — do NOT spawn L3:* rotation of a non-coordinating methyl/terminal group; the O–H torsion
  of a free (non-binding) –COOH or –OH; libration of a coordinated linker about the M–O axis that
  keeps the same donor bound.
- *Substantial — SHOULD spawn L3:* a functionalized linker tail that in one minimum folds over an
  open axial site and **blocks** it vs. points away and leaves it **open**; a pendant donor that
  self-satisfies via an intramolecular H-bond (deactivated) vs. is exposed; a chelate arm swung
  in/out so a metal can/can't reach the next pocket.
- *Borderline — C2 must legislate:* flexible-chelate ring puckers (EDTA λ/δ backbone; envelope vs
  half-chair); syn/anti rotamer of a bridging carboxylate shifting a bridge angle a few degrees.

Pin the thresholds *and the algorithm version* in the record — a canonical hash is only
reproducible if its recipe is pinned. (Open checkpoint C2: set θ_geom and the energy window.)

---

## 5. Data model (schema)

```
structures   : id, L0(formula,charge,spin), L1_graph_hash, L2_isomer_tag,
               typed_graph_blob, best_geometry_id, algo_versions

geometries   : id, structure_id, L3_conformer_id, fidelity(enum), method,
               coords_hash, energy, converged, relaxed_from, created

site_catalog : id, structure_id, canonical_atom_idx, donor_type            # perceive once
site_state   : id, site_id, geometry_id, status(open/occupied/blocked),
               pka, fukui, buried_vol, marginal_dE, ease_scalar, fidelity   # per-geometry

reactions    : id, product_structure_id, reagent_ids[], intermediate?,
               conditions, dG, barrier_proxy, atom_map, method, fidelity

pathways     : id, name, ordered reaction_ids[], path_score                 # A-vs-B compare
```

- `structures`↔`geometries` = the fidelity/optimization-status ladder.
- `site_catalog`↔`site_state` = perceive-once + refresh-accessibility.
- `reactions` + `pathways` = the polynuclear route-design layer; `atom_map` keeps site identity
  coherent across joins.

**Storage mechanics:** geometry coordinate sets are **content-addressed** (SHA-256, sharded
`store/ab/cdef….xyz`); the DB row references the hash. Meaning lives in DB rows, never in
filenames — the length limit disappears. **Backend:** SQLite single-file registry + the
content-addressed blob store (lean default). ASE-db is the least-new-code alternative given the
existing ASE stack; a graph DB is held in reserve for when the *route network itself* becomes the
object of study.

---

## 6. Subsystems

### 6.1 Structures + typed graph + polynuclear

The typed graph (§4.1) is the identity substrate. Adding metal–metal edges and bridging-ligand
node/edge types (µ2-O, µ3-OH, aqua, syn-syn bridging carboxylates) is all that polynuclear needs
— no new identity machinery. A Cu paddlewheel `Cu2(µ-O2CR)4` stores as one node: two
symmetry-related Cu, four occupied bridging sites, two open axial sites.

**C4 DECIDED: polynuclear-native from v1 — a single metal is just N=1, never a second standard.**
"Done right the first time" locks these commitments *before any hash is computed* (they touch
identity, so they cannot be retrofitted):

- Metals are ordinary graph nodes; **M–M edges and µ-bridge nodes are first-class from the first
  schema**, even while early test cases are mononuclear. No `if one_metal` branch anywhere.
- **Per-center oxidation state and spin are baked into node labels**, not one global charge/spin —
  required for Fe₃-µ₃-oxo (mixed FeII/FeIII), heterometallic nodes, mixed-valence Cu. This changes
  L0, so it is decided now.
- **Bridging-ligand representation is frozen now:** a µ₂-carboxylate = **one ligand node with two
  typed dative edges** (each carrying which metal + which donor atom), never split/duplicated.
- ~~**The geometry placer generalizes to "a set of coordination centers with inter-center
  constraints" — the headline engineering cost of the whole plan.**~~ **Superseded by D20**,
  measured against the M6 battery: a polynuclear node *emerges* from joins, the multi-centre
  placer reconciles only the edges where two determinants disagree, and the cost turned out
  not to be on the critical path. The argument it replaces is in
  [`archive/DESIGN_history.md`](archive/DESIGN_history.md). The rest of this bullet stands:
  identity (the canonical certificate) and site inheritance (atom map) get polynuclear for
  free.

### 6.2 Geometries + fidelity ladder

**Fidelity is a property of a geometry, not of a structure.** One structure identity, many
geometry realizations:

- `fidelity` ordered enum: `raw_construct → FF → ML → xTB → DFT`.  A rung is **not** a
  theory: ML is served by two MACE foundation models (MP-0, charge- and spin-blind;
  OMOL-0, charge- and spin-aware), whose energies share a rung and nothing else.
  `method_id` says which, and no query may order by energy across method rows.
- `method` records functional/basis/dispersion/solvent/spin + code+version (reproducibility).
- `relaxed_from` chain *is* the optimization-status history: `raw → ML-relaxed → DFT-relaxed`.
- Maintain a materialized **best-available geometry** pointer per structure for fast recall.
- A basin-crossing relaxation (DFT lands in a different minimum than the raw build) becomes a
  **new L3** with a `derived_from` edge to the coarse one — nothing is lost, and the cheap
  geometry is recorded as having correctly pointed at it.

Prior-art pattern: ASE-db KV / Materials-Project `TaskDoc` provenance — but SQLite is sufficient
at this scale.

### 6.3 Sites: perceive once, refresh accessibility

Store every site by its **canonical atom index** (from §4.1), not the build-order index — so
sites are stable across recall and across fidelity. Two tiers:

- **`site_catalog`** (geometry-independent): which donor sites exist, by canonical index and type.
  Perceived once per structure.
- **`site_state`** (per-geometry): open / occupied / sterically-blocked, plus the activation
  numbers — because accessibility genuinely changes when a buried site opens after a higher-fidelity
  relax. Recomputed only when fidelity improves; tagged with that fidelity.

**A site is a FRAME, not a vector.** A single outward vector is provably insufficient to specify a
join — it cannot represent the **torsion** (roll) about the new bond, and only partially the bend.
Fixing the relative pose of two rigid fragments across a bond needs a **local frame** (position +
two orthogonal reference directions), plus a symmetry tag saying which angular DOF are live. The
placer already computes this (`_donor_placement_frame`) and already searches the torsion
(`_best_azimuthal_rotation`), but throws both away — the move is to **promote them into the stored
site model** so geometry is reproducible, not re-derived by stochastic search. Site geometry
therefore carries:

- **`frame`** — donor position + outward (dative) axis + a reference direction (e.g. the carboxylate
  O–C–O plane, the pyridyl ring plane) that defines the torsion zero.
- **`live_dof_tag`** — which angular coordinates actually matter, a property of the functional group:
  carboxylate O and pyridyl N have a plane → **torsion is live** (syn/anti, cis/trans); aqua O and
  primary amine N are near axially symmetric → **torsion is a "don't care."** This tag is what tells
  the enumerator whether to branch torsion (spawn rotamer pairs) or ignore it — it is the guard
  against combinatorial explosion.
- **binding modes** — a functional group maps to a *set* of frames (carboxylate: monodentate /
  bidentate-chelate / bridging), not one frame.

The **bend angle** is mostly determined by donor electronics (sp² carboxylate O ⇒ ~120° lone-pair
direction) — derived, not stored, unless strain scoring later wants a nominal + stiffness. The
**realized join** (bond length + chosen torsion) lives on the assembly/reaction edge, and the
torsion is stored as a **discrete well index**, not a raw float — because that index *is* the
conformer-generating coordinate (see D11). Prior art for the connection-site-as-frame abstraction:
reticular builders (pormake, ToBaCCo) align building blocks by connection-point frames.

**Site compatibility** (can A join B?) then becomes a **frame-alignment test** at low strain, not
just "both open" — the assembly step (roadmap) needs this.

**Incremental perception:** when a new structure is built by joining stored blocks, inherit the
parents' sites through the assembly's **atom map** — minus sites consumed by the new bond, plus
any newly exposed ones. Perception is incremental, never from scratch. The atom map (carried on
each `reactions` edge) also gives site provenance for free ("this axial site descends from the Cu
that entered at step 1").

### 6.4 Activation-ease model

"Activation" = deprotonating a donor to expose it **and** that donor's readiness to form the next
bond. The honest quantity is `activation_ease(site, partner, conditions)`, not a bare scalar.
Store **named components** at whatever fidelity is available, plus a derived scalar for ranking;
never discard the parts:

| Component | Meaning | Cheap → rigorous |
|---|---|---|
| **deprotonation** (primary) | pKa or ΔE_deprot (implicit solvent) — the "activate = deprotonate" primitive, numeric | group-contribution/Hammett pKa → xTB/DFT ΔG_deprot |
| **electronic donor strength** | condensed Fukui f⁻ / partial charge / frontier-orbital coefficient on the donor | one xTB single point → DFT |
| **steric accessibility** (per-geometry) | cone/solid angle or buried volume %V_bur around the site | geometric, cheap |
| **marginal bond energy** | incremental ΔE of adding the next partner at that site — couples into pathway scoring | reuse formation-energy machinery |

Keep an **intrinsic scalar** (pKa/nucleophilicity, partner-free) as the default sort key, with an
optional **per-partner correction** layered on (a site easy for Cu²⁺ can be hard for Zr⁴⁺). Tag
every activation number with the fidelity/method that produced it.

### 6.5 Reactions + pathways (the route-design layer)

Model each candidate mechanism as a **path of stored intermediates** over structure nodes. Cu
paddlewheel example:

- **Path A (nucleus-first):** 2 Cu²⁺ → bridged Cu₂ pre-nucleus → carboxylates add → paddlewheel.
- **Path B (sequential):** Cu²⁺ + carboxylate → mono-Cu carboxylate → second Cu adds → paddlewheel.

Each **edge** is a `Reaction` (reagents, product, ΔG, optional barrier proxy, conditions,
atom_map); each **node** is a registry structure, so intermediates are first-class, storable,
reusable, and their stored open sites say what can add next.

**What writes those edges today.** A run with `spec.pathways` set plans, for every product it
was asked for, the rung one ligand below it — the coordinatively unsaturated intermediate, not
the co-ligand-filled complex, because only the former is reachable *by addition* — and performs
the step with `assembly.join`. So a ligand-count sweep (`ligands_per_metal: "1~3"`) comes back
as a chain of stored intermediates with real edges between them, which is Path B's shape for a
mononuclear centre. Under `co_ligand_counts: range` (D26) the co-ligand is stepped the same way
— `M(H2O)2 + H2O → M(H2O)3`, the free co-ligand a reagent — inside a window of empty vertices,
so the seeds of one metal/CN form one ladder rooted at the lowest co-ligand state in the window.
Under D2 the step usually lands on the identity the direct construction already built: one node, two routes, and the second route is the one that says where it came
from. Nothing here scores anything — the scoring below is still M8.

**Pathway viability** = a function over the edge set:

- rate-limiting step (max barrier / least-favorable ΔG),
- cumulative ΔG,
- **sink detection**: any intermediate so stable it stalls the route (the sequestration idea,
  generalized from endpoints to intermediates).

Comparing Path A vs. B = comparing these path-level scores. **Intermediate chelators** enter as
alternate branches (the chelator is a reagent later displaced); score "does chelator X open a
lower-max-barrier path to the same product." Mental-model anchors: chemical-reaction-network /
microkinetic graphs; cluster-nucleation pathway modeling.

### 6.6 Shared descriptor layer — the substrate for C5/C6/C7 (leaning, not locked)

Key observation: **hardness/softness (HSAB), formal charge, ionic radius, and ligand-exchange
lability recur across all three of the activation/kinetics/partner questions.** Build ONE small
descriptor layer — mostly tabulated, near-zero cost — and C5/C6/C7 all draw from it instead of each
inventing its own numbers:

- **per-donor descriptors:** donor_type, predicted pKa, HSAB class (hard/borderline/soft), Fukui
  f⁻/partial charge (when a wavefunction exists), denticity, in-pocket/H-bonded flag.
- **per-metal descriptors:** charge, ionic radius, HSAB class, preferred CN/geometry, high-spin
  d-count (already in `_D_ELECTRONS`), **water/ligand-exchange lability class** (Cu²⁺ fast,
  Cr³⁺/Al³⁺ slow — tabulated).

**C5 — activation-ease floor (leaning).** Floor = a zero-QM **heuristic-tier** record so no
`raw_construct` structure has an empty ease field: predicted pKa from a donor-type table (+ optional
additive/Hammett substituent nudge) + HSAB tag, stored as `fidelity='heuristic'` with an explicit
confidence. Higher tiers add Fukui/charge/marginal-ΔE. Donors sitting in a chelate/H-bond pocket
(e.g. the anthrarufin peri-OH, whose pKa is shifted by the quinone H-bond) **self-flag as
"provisional — promote to xTB"** — the heuristic tier tells you which sites are worth spending QM on.

**C6 — barrier proxy (leaning).** ΔG + sink-detection alone is kinetically blind, and the
paddlewheel question *is* kinetic. v1 = ΔG + sink-detection **plus two tabulated, zero-extra-calc
`barrier_proxy` components** that answer it directly: (a) **concurrent-bond-change count**
(many-body nucleus-first penalty) and (b) **metal exchange-lability class** (desolvation penalty).
Optional cheap add-ons: Coulomb approach penalty (two like-charged reagents), BEP relation
(barrier ≈ α·ΔG + β). Hook reserved for a mid-tier **1D relaxed scan along the forming-bond distance
at xTB** when a specific edge deserves a real number.

**C7 — partner dependence (leaning).** **Factorize, don't store a matrix.** Keep intrinsic site
descriptors + a per-metal partner-descriptor row; compute `ease(site, partner)` at query time via
an **HSAB-match term + charge/size/geometry-fit**; promote a specific (site, partner) pair to a
stored measured marginal-ΔE only when actually computed. Cost of adding a metal = one descriptor
row → multi-metal screening is O(sites + metals), not O(sites × metals). Reproduces the
anthrarufin O,O′ (borderline) + Cu²⁺ (borderline) = good-match → high-ease case for free.

Through-line: **build the descriptor layer first; it is mostly a lookup table and it is the shared
substrate all three sit on.** None of C5/C6/C7 should be a bespoke calculation.

### 6.7 Construction as an explicit decision tree (how conformer pairs are *generated*)

Conformers are not just detected after the fact — they are **manufactured at branch points** in the
construct call. Three kinds:

- **Kind A — stochastic/continuous** (ETKDG embed, distance-geometry, azimuth search): generate
  conformers *by accident*; these are **noise to dedup**, and are the real job of RMSD clustering.
- **Kind B — discrete/enumerable** (donor selection, ligand→CN-vertex assignment, dentation
  partition, geometry pick, **torsion well**): generate conformers *deliberately*; these **are** the
  L2/L3 distinctions and should be enumerated and labeled by which branch produced them.
- **Kind C — inference** (CN, coordination geometry, spin/charge, protonation): when ambiguous, a
  branch in its own right — **enumerate the alternatives as separate construct calls, don't silently
  default.**

Treat construction as a tree: Kind-C/Kind-B choices are internal nodes, leaves are candidate
structures, and branch depth maps to identity level (charge/spin ⇒ L0; connectivity/dentation ⇒ L1;
vertex-assignment/stereo ⇒ L2; torsion/rotamer ⇒ L3). "Conformer pairs" are sibling leaves — or
leaves that converge/diverge under relaxation. The path to each leaf (the choice-vector) is stored
as generative provenance, reusing the same edge/atom-map machinery pointed at construction. This is
what makes stored conformers **reproducible and regenerable at higher fidelity**, not just frozen
coordinates.

---

## 7. Decision Ledger (locked)

- **D1** Central object = recursive `BuildingBlock`; assembly re-exposes open sites.
- **D2** Split identifier / address / provenance. Address = content hash; identity = composite
  L0–L3 key in a DB record; provenance = DAG edges. Identity on node, sequence on edges.
- **D3 (producer superseded by D16)** L1 = canonical hash of an **explicit typed graph**. Never
  via SMILES/InChI — that half stands and is the part this decision exists for. The *producer*
  was originally specified as WL-primary with nauty canonical and VF2 tiebreak; **D16 replaces
  it** with an individualisation-refinement certificate, after 1-WL turned out to be unable to
  separate a µ2-bridging carboxylate from a chelating one. Of the original three, WL survives as
  a bucket index, VF2 survives as the collision resolver, and nauty was dropped entirely.
- **D4** Fidelity/status/coords in a child `geometries` table; structure identity is
  fidelity-invariant. `relaxed_from` chain = optimization history; keep a best-geometry pointer.
- **D5** Sites stored against **canonical atom indices** → perceive once; split into
  geometry-free `site_catalog` + per-geometry `site_state`; inherit sites across joins via atom map.
- **D6** Activation ease = **named components** (deprotonation, electronic, steric, marginal-ΔE)
  + derived scalar, tagged by fidelity; intrinsic scalar default + optional per-partner correction.
- **D7** Polynuclear via typed graph (M–M + µ-bridge edges) — no new identity machinery.
- **D8** Pathways = reaction DAG over structure nodes; viability = path-level score (max barrier +
  cumulative ΔG + sink detection); chelators = alternate branches; `atom_map` on every edge.
- **D9** Backend = SQLite registry + content-addressed blob store (ASE-db fallback; graph DB in
  reserve).
- **D10** L2 is in scope; discriminator = **downstream-assembly-relevance + barrier-separated
  persistence**, NOT the energy gap (near-degenerate cis/trans are still distinct nodes). Justifying
  case: anthrarufin–Cu. (Resolves C1 in favour of L2-aware identity.)
- **D11** L3 trigger = compound rule: RMSD on **rigid core + coordination sphere** OR a flip of any
  **open-site accessibility flag**, gated by an energy window. Open-site flag is the primary trigger.
  **L3 identity is provenance-primary, geometry-verifier:** the conformer's label is the discrete
  construction **choice-vector** (Kind-B/C branches: donor selection, vertex assignment, torsion-well
  index), and geometric clustering only (a) collapses stochastic embedding duplicates and (b)
  reconciles divergence (one choice → many wells) / convergence (many choices → one well). Near-
  degenerate cis/trans survive threshold-merging because they carry different choice-vectors.
- **D13** A site's geometry is a **frame + live-DOF tag + binding-mode set**, not a lone outward
  vector (§6.3). Torsion about a join is stored on the edge as a **discrete well index** = the
  conformer coordinate. `construct` is a deterministic function of (choice-vector, seed) and **emits
  its choice-vector**; the inference layer (CN/geometry/protonation) **branches explicitly** on
  ambiguity instead of silently defaulting.
- **D14** **Bridging is DERIVED, not an edge type.** §4.1 originally typed `bridging(µ2,µ3)`
  as an edge alongside `dative`, but a donor atom carrying dative edges to two distinct metals *is* µ2 —
  labelling it again gives two encodings of one chemistry, and the two hash differently.
  `EdgeType` is therefore `{COVALENT, DATIVE, METAL_METAL}` and `TypedGraph.bridge_class()`
  derives µ2/µ3 on demand. Dative *direction* is derived the same way (donor is the non-metal),
  and validated on insert.
- **D15** **Net charge is a graph-level field, not a sum over atoms.** Writing a
  carboxylate's −1 onto one of its two oxygens makes those oxygens inequivalent, so a
  paddlewheel would hash differently depending on which way round four chemically identical
  bridges happened to be written. Delocalised charge is not given a home it does not have.
  Per-atom `formal_charge` survives for genuinely localised charge (an ammonium N) as a label
  that enters identity but not the total. This is safe *because* hydrogens are explicit nodes:
  protomers stay distinct through the H count, not through where the charge was written.
  Bond order is excluded from the hash for the same reason (Kekulé forms, C=O/C–O resonance).
- **D16** **L1 = sha256 of the canonical certificate, not the WL hash.** D3 specified WL as
  the primary stored key. It cannot be: **1-WL does not separate a µ2-bridging carboxylate from
  a chelating one** — an 8-membered M–O–C–O–M–O–C–O ring versus two 4-membered chelate rings
  give every atom the same local environment, so the colours are stable from the first iteration
  and more iterations do not help. That is a binding-mode distinction this project exists to
  draw. WL is demoted to a fast bucket index; L1 comes from an individualisation-refinement
  canonical certificate (pure Python, so both machines agree — pynauty is dropped rather than
  kept as a verifier, since its canonical form differs and would produce different keys). Cost measured on the fixture
  set: 120 ms total, worst case 37 ms for [Fe(H₂O)₆]²⁺ (|Aut| = 46080) with branch-and-bound
  plus automorphism pruning. Demonstration lives in `tests/test_canon.py`.
- **D17** **An energy difference requires an isodesmic equation, not merely a
  balanced one.** Raised while implementing M7. Balance in atoms and charge is necessary
  and demonstrably not sufficient: the archived `E(EBU) − E(M^q+) − Σ E(anion)` scheme
  satisfies it and still reversed its own qualitative verdict once a medium was added.
  `energy.reference` therefore refuses, by default, any equation containing a bare metal
  ion, a metal-free species with |charge| > 1, or a net change in metal–donor bond count.
  Reproducing a legacy number needs `strict=False`, and the resulting value carries
  `isodesmic=False` into anything that stores it. The cost of this decision is that some
  equations a user considers reasonable will be refused; the alternative is a table of
  numbers that all look equally good.
- **D12 (C4 resolved)** **Polynuclear-native from v1**; mononuclear = N=1. Per-center oxidation
  state/spin in node labels; µ-carboxylate = one node + two typed dative edges; multi-center geometry
  placer is the headline build cost.
- **D18 (C5 resolved)** **The activation-ease floor is the zero-QM heuristic tier, and an
  absent component is absent — never defaulted.** Every perceived site gets a record from
  the M1 tables alone, so no structure carries an empty ease field; a site the model cannot
  score is reported *unscored*, never *hard*, because a missing number that renders as a low
  one is the failure this layer exists to prevent. The scalar is a weighted mean over the
  components actually present (`deprotonation` 0.55, `steric` 0.25, `electronic` 0.10,
  `marginal_de` 0.10), renormalised by those weights: nothing is filled in with 0, or 0.5, or
  "neutral", since a defaulted component is indistinguishable from a computed one.
  `confidence = coverage x sharpness` answers *how much of the model ran*, not *how right is
  it* — coverage from the weights present, sharpness from the table's own `pka_sigma`.
  `provisional` marks sites where **the table value is the wrong question** — not sites that
  merely scored badly — making the floor's most useful output a QM work list rather than a
  verdict. Two things sharpen it from the §6.6 leaning, both learned by measuring what it
  flagged: the pocket must be **inter-group** (a chelate ring ≥ 5), because a carboxylate
  closes a 4-ring through its own two oxygens and "is in a pocket" therefore flags every
  carboxylate in existence — and a work list containing everything is no work list; and a
  donor with **no pKa never self-flags**, since promoting a pyridyl N to xTB to refine a
  deprotonation it does not undergo spends QM to learn nothing. The rule now separates
  salicylate and catechol (pKa₁ 2.97 vs benzoic 4.20; 9.25 vs phenol 9.99 — real neighbour
  effects) from benzoate, BTC and *para*-hydroxybenzoate (where the table value is simply
  correct). A neutral donor's empty `pka` scores 1.0 with a note: "no activation step" is an
  answer, not a gap.
  **The rung above the floor is MACE-OMOL-0, and that is what made this callable now.**
  Deprotonation is charge-changing (`A-H -> A(-) + H(+)`), so a charge-blind potential is
  blind to it in principle, and the ladder for the primary component used to run
  `table -> (nothing) -> xTB` — the ML rung, the one that makes screening affordable, could
  not serve the component the model is mostly made of. Cost of this decision: the scalar's
  pKa anchor (centre 10.0, width 3.0) is a stated convention, not a calibration, which is
  precisely why the components are the record and the scalar is only a sort key (D6).
  **C7 is untouched**: `hsab_match` still raises and `activation_ease(partner=...)` refuses
  rather than returning the partner-free number under a partner-shaped call.
- **D19** **A stored identity is never re-derived under a new recipe version; it keeps the
  answer its own version gave.** Forced by M5, which fills in `l2_isomer_tag` and so would
  otherwise **split** every identity written while L2 was a stub. Backfilling is refused:
  re-deriving L2 for the existing rows would rewrite stored identities and every `reactions`
  edge pointing at them, which is the one thing in this registry that is not regenerable
  (D2). So `ALGO_VERSIONS["l2_isomer_tag"]` bumps from `0-stub`, rows written under the stub
  keep `l2_isomer_tag = ''`, and anything built from M5 onward carries a real tag.

  **What makes this honest rather than a silent fork is that the generation is already on
  the row** — `structures.algo_l2` records which recipe produced each value, so a `''` is
  readable as *"this predates L2"* rather than as *"this has no isomerism"*. The cost is
  real and is accepted: the same species built before and after M5 can occupy two rows.
  That is a visible duplicate with its cause recorded, which is strictly better than an
  invisible one — and it follows ground rule 6, which already says a version bump marks
  rows stale rather than re-labelling them.

  **The narrower case goes the other way, and the difference is the point.** `site_catalog`
  IS re-derived on a version bump (`perception/1` → `/2`), because a catalog is a derived
  annotation that nothing points at, so rewriting it costs only the state rows underneath
  it. An identity is pointed at by the provenance DAG. The rule is therefore not "always
  re-derive" or "never" but: **re-derive what is only an annotation; version what is an
  address.**

- **D20 (supersedes §6.1's last bullet)** **A polynuclear node EMERGES from joins; the
  multi-centre placer is for RECONCILIATION only.** §6.1 called the multi-centre placer "the
  headline engineering cost of the whole plan". Measured against the M6 battery, it is not:
  a Cu paddlewheel builds QC-clean and hashes to its M2 fixture at **both L0 and L1** from
  machinery that already existed — the sp2 donor's second lone-pair well, `kabsch`, the
  curated donor-distance table, the existing `qc` — with no solver anywhere. So **M···M is an
  output to validate, not an input to impose**, and `data/reference/node_cases.tsv` carries
  the literature window a built node is measured *against* rather than a target a placer is
  driven *towards*.

  `place_multicentre` is kept, for a measured reason rather than a blanket one: where two
  determinants fix the same edge they disagree. A µ₃-oxo asks 3.291 Å of Fe₃ while syn-syn
  formate offers 2.696 (0.595 Å apart); µ₄-oxo asks 3.168 Å of Zn₄O against 2.677 (0.491 Å).
  That is **one scalar on one edge**, not a general constrained optimisation, and the
  difference matters: a determined skeleton is *built*, an under-determined one **raises**
  naming what is missing, and neither falls into a minimiser.

  The residual is reported as strain and closed by relaxation. Averaging the two determinants
  is specifically refused — it stores a number neither determinant asked for, and the
  disagreement is a real property of a rigid-ligand model, not noise to be smoothed. Cost of
  this decision, stated up front: the model is quantifiably wrong by ~0.5 Å on oxo-centred
  clusters and says so, rather than being plausibly wrong and silent.

  The one place a distance is a legitimate *input* is a **declared nucleus**, where the caller
  states it on purpose — which is a declaration, not an inference, and so is the same rule.

- **D21 (C2, geometric half)** **θ_geom = 0.15 Å, calibrated and not chosen.** The threshold
  sits in the valley between two measured populations of core-RMSD over the rigid core plus
  coordination sphere, 12 xTB-relaxed structures over 4 choice vectors × 3 embedding seeds:

  | population (xTB-relaxed) | n | min | median | max |
  |---|---|---|---|---|
  | Kind A — one choice vector, different embedding seed | 12 | 0.0000 | 0.0000 | **0.0316** |
  | Kind B — different choice vectors | 54 | **0.6435** | 0.8320 | 1.1871 |

  **RAW constructs could not set this number, and that is the first result rather than an
  obstacle.** A join is a deterministic function of its choice vector and absorbs the ligand's
  embedding noise completely, so at RAW the Kind-A spread is *identically zero* (105 pairs,
  max 0.0000) — there is no stochastic peak for a threshold to sit above, and a number
  calibrated there would be calibrated against no noise. Relaxation re-introduces it.

  0.15 is the **geometric** mean of the two bounds, not the arithmetic one, because these are
  ratios of distances: it puts 4.7× margin above the widest duplicate and 4.3× below the
  closest real branch, where an arithmetic midpoint would sit 20× above one and 1.5× below the
  other. Both bounds ship beside the constant as `CALIBRATION_KIND_A_MAX` /
  `CALIBRATION_KIND_B_MIN`, and a test asserts the threshold stays between them with margin,
  so the number cannot drift away from the data that set it.

  **C2's energy half stays open, deliberately.** Over the same set the Kind-A *energy* spread
  reached 16.6 kcal/mol between samples whose cores agreed to 0.03 Å — all of it motion
  outside the core, because the rigid-core rule cuts a delocalised carboxylate C–O as if it
  were rotatable ([B9](BUGS.md#b9)). A window set from that data would bake a known defect
  into a stored threshold. `DEFAULT_ENERGY_WINDOW is None`, the mechanism is built and tested,
  and the number waits for the core fix.

- **D22 (C2's storage corollary; narrows D19)** **A tag is derived where the coordinates are,
  and a caller that shares a node with an untagged path says so explicitly.** `put_structure`
  computes L2 from a graph alone, and a `structures` row exists before its geometry does — so
  only the builder holds coordinates at insert time. `store_block` therefore derives the tag
  by default and takes `l2=` as an override for the one case that needs it: the run pipeline
  does not tag isomers, so a route that tagged its own product would file it apart from the
  node every other route reached ([B2](BUGS.md#b2)). It passes `l2=""` **on purpose**, and
  closing B2 moves both paths together or neither.

  No backfill, in either direction. Of the 39 structures stored under the stub, 18 would gain
  a real tag and 21 would still tag `""` — and re-deriving them would rewrite stored
  identities and every `reactions` edge pointing at them, which is the one thing here that is
  not regenerable (D2). `structures.algo_l2` records which recipe produced each value, so a
  `""` reads as *"this predates L2"* rather than as *"this has no isomerism"*. The
  choice-vector digest goes the other way and for the same rule: it is an **annotation**
  nothing points at, so `ix_geometries_choice` is repopulated by rewriting, not versioned.

- **D23** **`MAX_BITE_MISMATCH_DEG = 40°`, and what it separates is pinned rather than the
  threshold itself.** A chelate's verdict is its bite angle against the angular separation of
  the two vertices it is asked to span. Against an octahedral *cis* pair (90°) acetate
  measures 59.7°, a 30.3° mismatch, and it is the strained end of what really forms — the
  common chelators (acac ~92°, en ~85°, bipy ~78°) sit within 12°. Against a *trans* pair
  (180°) every one of them misses by 88° or more, acetate by 120°. Anything from ~35 to ~60
  draws the same line; 40 is taken from the low half so the strained-but-real case passes with
  margin while nothing comes near trans.

  The number is one named constant in one module, and the *populations* are what the test
  asserts (`test_the_bite_angle_populations_stay_far_apart`) — so a donor type whose geometry
  moves cannot quietly cross the line, which is the failure mode "pick a tolerance from a
  figure" always has.

- **D24 (C9 resolved)** **A polynuclear node's multiplicity is STATED, never combined.** Two
  centres' unpaired electrons add only if they are independent, and whether two d⁹ Cu(II) give
  a singlet or a triplet is exchange coupling — not recoverable from the centres, and so not
  derivable by the rule that derives it. The two ground-truth fixtures make the point without
  argument: `examples.cu_paddlewheel` declares multiplicity **1** (AF-coupled d⁹–d⁹) where
  `energy.backends.combined_multiplicity(2, 2)` gives **3**, while `examples.fe3_mu3_oxo`
  declares 16 and the additive rule agrees.

  Multiplicity is in L0, so this decides identity: a derived 3 would file the paddlewheel at
  `Cu2_C4H4O8_q0_s3` and miss its own fixture. The multi-centre path therefore **requires** a
  multiplicity and raises otherwise — which is already what `from_rdkit` does, and which makes
  both fixtures correct instead of one of them wrong. `combined_multiplicity` keeps its job on
  the mononuclear path, where "one metal plus a closed-shell ligand set" genuinely is additive.

- **D25 (C10 resolved)** **A METAL–METAL edge is DECLARED, never inferred from distance.** At
  2.673 Å two Cu are bonded and at 5.516 Å they are not, so a threshold looks available — and
  taking it would be exactly the silent inference ground rule 5 forbids, applied to the one
  field that is in the certificate. `EdgeType.METAL_METAL` is hashed into L1, so a guessed
  edge is a guessed *identity*.

  The battery settles it as data rather than as argument: the Fe₃ trimer at **3.29 Å has no
  M–M edge** and the Cu₂ paddlewheel at **2.62 Å has one**, the two fixtures differ at L1 by
  precisely that, and `node_cases.tsv` carries `mm_bond` as its own column saying it is a
  chemical decision not derivable from `d_mm`. So: branch on it or declare it; never default.
  `InterCentreConstraint.metal_metal_bond` is the input side of the same rule.

- **D26** **A new spec builds co-ligand counts in a window below saturation and steps between
  them; an old spec replays the run it planned.** `BuildSpec.co_ligand_counts` is `fill` (a
  co-ligand on every vertex the ligands leave — the only behaviour before v8) or `range`, the
  v8 default: `n_co ∈ [max(0, full − w), full]` with `w = co_ligand_window`, default **2**,
  the uncovered vertices left **empty in the requested polyhedron**. It applies to whatever
  co-ligand or solvent the spec names; water is only the default `co_ligand`. CN is never
  collapsed and never invented: a count *above* a CN's full needs a higher CN, and exists only
  where the spec's coordination list has one (CN 4 and 6 give a one-ligand Ni 1–3 and 3–5
  waters). `allow_unsaturated` does not gate these products — it answers "what if nothing can
  fill the vertices", and `range` with a co-ligand is itself the request for them.

  Why the default moved: with `fill`, the MVP registry (`data/mvp_ni_thq_cl.db`, built from
  `spec_ni_thq_cl_slice.json`) bottoms out at two unrelated seeds, Ni(H2O)2 and Ni(H2O)3, with
  no incoming edge and no common ancestor, and "the same complex with one water fewer" is
  neither a node nor an edge. Solvent gain and loss is what a solvated centre actually does, so
  it belongs in what a run produces unasked.

  With `pathways`, the window also bounds the ladder: **no rung may leave more than `w`
  vertices empty**, for ligand steps and co-ligand steps alike, and a co-ligand step joins every
  pair of in-window rungs one co-ligand apart. The seeds of one metal/CN/polyhedron are then one
  connected ladder whose root is the lowest co-ligand state inside the window — Ni(H2O)3 ←
  Ni(H2O)2 + H2O at CN 4 — and **no bare-metal row is created**. A composition's lowest
  in-window state is a root of its own (nothing below it is in the window) and is reached from
  the ligand-free chain through the rungs above it. `w = 0` allows no empty vertex, so it plans
  no step and says so. The free co-ligand is its own task (`co_ligand`), stored and relaxed at
  the run's fidelity like a ligand, so a co-ligand edge cites `[rung below, free co-ligand]` —
  the ligand step's shape — and prices. `place` edges still cite nothing (C15 unchanged).

  Size, measured on the reference slice re-read at v8: `fill` 722 tasks; `range` with `w` = 1
  1141; **`w` = 2 (the default) 2041**; `w` = 3 2521; CN `4,6` at `w` = 2, 3059. `MAX_PATHWAY_TASKS`
  was raised from 2000 to 4000 so the default window fits CN `4,6` uncut.
  The cap is reported with the uncapped size rather than truncating silently. Migration: every
  `spec_version ≤ 7` reads `co_ligand_counts: fill` (window 2, inert under `fill`) and plans
  byte-identical tasks — pinned against the pre-v8 planner by digest
  (`tests/test_co_ligand_counts.py`), not by comparing the migration with itself. A charged
  co-ligand is not stepped until B20 is fixed; the planner says so.

Every entry above is locked and has a test that fails if it is reversed. Revision
history — how each one was argued and what it cost — is in
[`archive/DESIGN_history.md`](archive/DESIGN_history.md).

## 8. Open Checkpoints (need a call)

*Resolved: C1 → D10 (L2-aware). C4 → D12 (polynuclear-native). C5 → D18 (heuristic floor).
C8 → curated tables (`data/reference/*.tsv`, each row carrying `source` + `source_version`).*

Each open checkpoint is scheduled at the milestone where code first forces the call — see
[`PLAN_implementation.md`](PLAN_implementation.md) §3 for that placement.

- **C2 — L3 thresholds (still open):** the rule is set (D11); still need the **numbers** — θ_geom
  RMSD cutoff on the rigid core + coordination sphere, and the energy window that gates a real
  minimum vs. a bad geometry. Also governs whether a basin-crossing DFT relax spawns a new L3.
- **C6 — Barrier proxy (leaning, §6.6):** confirm v1 = ΔG + sink-detection + concurrent-bond-change
  + exchange-lability proxies (with the 1D-scan hook), vs. thermodynamics-only. Input constraint
  already fixed: a **pivot** node on a composed route (`pathways.route`, `pivot: true`) is
  bookkeeping, and its `y` must never be read as a barrier or an intermediate by any proxy.
- **C7 — Partner-dependence (leaning, §6.6):** confirm the factorized HSAB-match model (descriptor
  vectors combined at query time) over a stored ease matrix.
- **C15 — what tells a placed structure from a grown one (open, found by measurement):** the
  ladder currently distinguishes them by **whether the edge has reagents**: a `place` edge cites
  nothing, a `grow` edge cites the rung below and the ligand added. `test_the_step_reaches_the_
  structure_the_place_task_built` asserts both kinds arrive at one node, and
  `test_the_ladder_reaches_down_to_the_bare_centre` walks down by following only reagent-bearing
  edges.
  That discriminator is load-bearing and it is also the reason a `place` edge cannot be priced:
  with no reagents, `reaction_terms` injects the product alone and the balance report is its
  entire composition. Giving `place` its real reagents (the ion, the ligands, the co-ligands —
  the runner knows all three) was **built and reverted** on this branch, because it makes every
  edge reagent-bearing and the ladder walk then descends into the bare ion. The patch is small;
  the decision is not.
  *Options:* discriminate on `kind` instead of on reagent presence, and let the ladder reach the
  literal bare centre · keep `place` citing nothing and answer "what is this made of" by
  derivation only (`energy.routes.decompositions`, which needs no rebuild) · record the pieces
  under a role that the ladder walk ignores.
  *Note:* `put_reaction` can now express stoichiometry, which was the blocking prerequisite
  either way — a diaqua complex has to be able to say *two* waters.
*(C8 — descriptor-layer sourcing — was here; resolved in M1 in favour of curated tables with
per-row `source` + `source_version`. See `archive/PLAN_completed.md` rev 19.)*

## 9. Suggested phased roadmap

Per D12 the schema is **polynuclear-native from step 1** — no mononuclear-only interim standard.
The multi-center geometry placer is sequenced where its cost is unavoidable, but the *identity/schema*
carries M–M + µ-bridges + per-center labels from the start.

1. **Consolidate** the code to one lineage; freeze the current pipeline as the "raw_construct +
   xTB-screen" fidelity levels.
2. **Descriptor layer** (§6.6): the tabulated per-donor / per-metal descriptor tables — cheap, and
   the shared substrate for C5/C6/C7. Small, do it early.
3. **Registry v1:** typed graph with **M–M + µ-bridge edge/node types and per-center oxidation-state
   /spin labels baked in**, L0/L1 keys (L2 tag from D10), SQLite + content-addressed store,
   `structures` + `geometries` tables. Migrate existing `.xyz`/`MODELS`. (L3 slot reserved, stubbed.)
4. **Sites v1:** `site_catalog` from existing donor perception, keyed to canonical indices, each
   site carrying a **frame + live-DOF tag + binding-mode set** (D13, promote `_donor_placement_frame`
   out of the placer); `site_state` with the heuristic-tier activation-ease floor (C5).
5. **Recursive assembly (mononuclear = N=1 path of the general operation):** join by **frame
   alignment** on compatible sites, torsion chosen from the live-DOF well set and recorded as a
   choice-vector coordinate; atom-map carry; `construct` emits its choice-vector. Validate a
   structure built two ways lands on one L1 node, and that cis/trans survive as distinct L3s on
   provenance.
6. **Multi-center geometry placer** — the headline cost (D12): coordination centers with inter-center
   constraints (M–M distance, bridge bite angle) + per-center local geometry. Build + store a Cu
   paddlewheel and an Fe₃-µ₃-oxo trimer as ground-truth targets.
7. **Pathway layer:** `reactions` + `pathways` with the barrier proxies (C6); run the paddlewheel
   A-vs-B comparison as the first real route-design result.

## 10. Prior art / tools to evaluate (don't reinvent)

- **MOFid / MOFkey** (Snurr group, *Cryst. Growth Des.* 2019) — an InChI-like canonical identifier
  *for MOFs* (node + linker + topology). Reference design / validation oracle; targets periodic
  frameworks, not arbitrary partial assemblies. (Check current maintenance state.)
- **pymatgen** `StructureMatcher` / `MoleculeMatcher` — battle-tested structural-equivalence/dedup
  logic to borrow for the L3 comparison.
- ~~**pynauty** (nauty/Traces) — canonical graph labeling.~~ **Evaluated and rejected (D16):**
  its canonical form differs from ours, and a key producer that is installed on only one machine
  breaks two-machine parity. **networkx** — WL bucket index + VF2 collision resolver; both used.
- **RDKit** — USR/USRCAT shape fingerprints; `CanonicalRankAtoms` for organic sub-parts.
- **ASE-db** — atomistic KV store; the least-new-code registry backend given the existing stack.
- Buried-volume / %V_bur descriptors (SambVca-style) for steric accessibility.

## 11. Known risks / weaknesses to keep in view

- ~~Gas-phase GFN2-xTB on isolated highly-charged anions: only within-metal *relative* rankings are
  trustworthy; a better reference scheme (reaction-balanced, consistent-charge, protonation-aware)
  is needed before energies drive route viability quantitatively.~~ **Addressed in M7, and the
  diagnosis above was half wrong.** "Consistent-charge" was already true: the archived equation
  `E(EBU) − E(M^q+) − Σ E(ligand anion)` balances in atoms *and* in charge, which is why nothing
  caught it. What makes it untrustworthy is that its reference species are not comparable to its
  product — a bare cation with no ligand field, a polyanion with nowhere to put its charge, and
  **six metal–donor bonds appearing out of nothing**. `energy.reference` therefore checks two
  things: balance (necessary) and an isodesmic condition (sufficient in practice) whose sharpest
  rule is that the dative-bond count must be equal on both sides. Absolute energies at this level
  of theory are still not claimable; what is now enforced is that a *difference* is only computed
  between species the method describes comparably.
- No topology/tiling check yet — a correct node+linker doesn't prove it forms a periodic net.
- Canonical-hash reproducibility depends on pinning the algorithm + method versions (store them).
- Activation ease as a single scalar is lossy — the named-component design mitigates but the
  partner/condition dependence (C7) is real.

## 12. Changelog

Moved to [`archive/DESIGN_history.md`](archive/DESIGN_history.md) — every revision, newest
first. This file carries the *current* reasoning; that one carries how it got here.
