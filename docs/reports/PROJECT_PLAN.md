# MolCluster — project plan

## Why this project

### Where computational MOF research is strong
Computational work on metal–organic frameworks is heavily weighted toward **applications and
structure screening**, and it is mathematically elegant:
- **Reticular chemistry** treats a framework as a net: nodes (SBUs) joined by linkers [62, 63].
- **Topological generators** enumerate hypothetical frameworks by placing known SBUs and linkers onto
  RCSR nets — hMOF [64], ToBaCCo [10], AuToGraFS [11], PORMAKE [9]. The result is databases of 10⁵–10⁶
  structures.
- Those structures are **screened for properties**: gas uptake by GCMC, electronic structure by DFT
  (QMOF [65]), stability and more, increasingly by ML. The experimental CoRE MOF [66] and ARC-MOF [67]
  collections serve the same pipeline.

### Where it is weak: the synthesis end
Every one of these pipelines **starts from the building unit, and assumes it forms**. Whether a given
metal salt and ligand, in a given solvent, actually assemble into that SBU — rather than into a
different cluster, a coordination polymer of the wrong connectivity, an amorphous solid or nothing —
is outside their scope. That is the step where syntheses succeed or fail, and it is governed by
solution chemistry that the topology does not see:
- **Metal speciation.** Aqua and hydroxo complexes, and hydrolysis, depend on pH.
- **Ligand protonation.** A carboxylic acid or polyphenol must lose protons to bind, and something
  must take them.
- **Competition from counterions, modulators and solvent** for the same coordination sites.
  Modulators (acetate, formate, benzoate) are routinely used to steer phase, crystallinity and defects
  [68, 69]. The salt anion — chloride, nitrate or acetate — changes the product.
- **Several near-degenerate clusters** at similar energy, where conditions tip the balance.

So MOF synthesis is still optimised by **empirical grids** of salt, solvent, temperature, modulator
and pH, and a successful synthesis is often a product of chance. Current approaches treat the
problem from the outside:
- **Statistical models of synthesis conditions** mined from the literature [70, 71] learn which
  conditions tend to work, without the molecular chemistry of why.
- **Framework free-energy calculations** [72] judge whether a *finished* framework is plausible,
  not whether its node can be reached.
- **In situ nucleation and growth studies** [73] reveal the mechanism case by case,
  experimentally.

### Working hypotheses about the synthesis problem
Four hypotheses shape what the project computes. Each is stated with what the project computes for
it, what would test it, and where it stops. H0 is the premise the other three rest on.

**H0 — What forms in solution carries over into the framework, in some systems and not in others.**
The whole approach assumes that the species favoured before the solid appears tell you about the node
that ends up in it. That holds where the node, or a recognisable piece of it, already exists among the
precursors. It need not hold where the node only assembles during nucleation or at the growing
surface, and the solution species is just feedstock.
- *What the project computes:* whether the target node, or a fragment of it with the same bridging
  pattern, appears among the low-lying species for a given recipe.
- *Test:* systems where the node is known to pre-exist should show it among the favoured species.
  Copper(II) acetate, already a paddlewheel dimer, is the reference case. Systems where the framework
  forms from quite different precursors should not, and the plan should then say its prediction does
  not apply.
- *Stated limit:* this is a condition of validity, not a result. Where it fails, H1–H3 say nothing
  about the framework.

**H1 — Monomer metastability decides whether a regular, crystalline framework can form.**
The first metal–ligand species to assemble (the monomer, and the small clusters built from it) has
to sit in a *metastable* window:
- **Too stable**, and it is a sink. The metal stays locked in a discrete complex and never extends
  into a framework.
- **Too unstable**, and the ligand never holds the metal long enough to build on. Alternatively,
  bonds form irreversibly as soon as they can, and the product is disordered.

To be testable the window needs numbers fixed in advance. Two quantities define it:
- **formation**: the free energy of forming the monomer from its ingredients, per new metal–ligand
  bond;
- **extension**: the free energy of the next step up, monomer to dimer or to the first bridged
  species.

The hypothesis is that both are downhill, and neither by more than a bound. The bound is set on one
half of a set of systems with known outcomes, and tested on the other half.
- *What the project computes:* both quantities for every monomer, in every hydration and
  protonation state.
- *Test:* on the held-out half, crystalline products fall inside the window; discrete complexes sit
  beyond the formation bound; amorphous or failed syntheses miss on extension. Mononuclear copper
  intermediates have been observed during milling of a copper framework [77], so the species this
  hypothesis is about are real and observable.
- *Stated limit:* metastability and reversibility are kinetic properties, and these are
  thermodynamic depths. Depth is used as a proxy for how readily bonds exchange. The proxy is
  reasonable for labile metal ions and is not expected to hold for inert ones.

**H2 — The synthesis recipe is the object to design. It sets the driving force, the height of the
worst intermediate, and the competing sinks.**
This is the project's main focus. A recipe — which salt, which base or modulator, which solvent, in
what order — fixes which species are present and the energy level of each. The aim is to **regulate
those levels so that intermediate species form as smoothly as possible**. "Smooth" is made of three
quantities, each of which the project can defend without computing a reaction barrier:

| Quantity | What it says | Why it is legitimate |
|---|---|---|
| **driving force** | how far downhill the whole recipe is, ingredients to target | a difference between two states; independent of path |
| **highest intermediate** | the highest point any route must pass through | a barrier cannot be lower than the intermediate it passes, so this is a *lower bound* on the kinetic cost |
| **competing sinks** | species off the route that lie below the target | where the metal ends up at equilibrium if they exist |

A smooth recipe is strongly downhill overall, has no high intermediate, and has no sink below the
target.
- *What the project computes:* all three for a route under a stated recipe, so that two recipes for
  the same target can be compared. The first study is an instance of the first quantity: giving the
  released proton to acetate instead of water changes the net reaction, and with it the driving
  force. It does not smooth an otherwise unchanged path.
- *Test:* recipes that score better on the three should be the ones that work at the bench (§5.4).
  The salt is known to change the product for a copper framework [78], which gives cases to check.
- *Stated limits:*
  - Where the chemistry is fast and reversible, only the driving force and the sinks decide the
    outcome; the height of an intermediate matters only where the reaction is slow.
  - A lower bound on a barrier is not the barrier.
  - The steps on a computed route are the order in which species were *built*. They form a valid
    thermodynamic cycle, but they are not necessarily the sequence the chemistry follows.

**H3 — Synthesis conditions can be approximated by the degree of solvation shell in the model.**
How much solvent surrounds the reacting species is the main thing that separates one synthesis
method from another, and the model has a matching setting for each:

| Model treatment | Synthesis method it approximates | What the species see |
|---|---|---|
| gas phase, no solvent | **dry high-energy ball milling** of anhydrous reagents | each other only; nothing screens a charge |
| gas phase + one explicit solvent shell | **wet milling** (liquid-assisted grinding), and dry milling of *hydrated* salts, which bring their own first-shell water | a few solvent molecules on the metal's open sites, and no bulk liquid |
| explicit shell + continuum | **wet-phase (solution) synthesis** | a first shell plus the screening of the bulk solvent |

So one enumeration of species can be priced three ways, and the same recipe compared across methods.
A recipe (H2) that scores poorly in solution may score well in the mill, or the reverse.
Mechanochemical routes to MOFs are established experimentally [74–76], and one copper framework has
been made by dry milling, by wet milling and from solution, with the added liquid changing the
pathway [77]. That makes the comparison testable.
- *What the project computes:* every route under each treatment. The first study already holds all
  three for the Ni–tetrahydroxyquinone system.
- *Test:* for one target made by all three methods, the treatment matched to each method should rank
  the recipes the way the bench does. The calibration against measured pKa and stability constants
  (§5.2) checks the solution end.
- *Stated limits:*
  - **A mill is not a gas.** Reagents and products in a mill are solids, and ions in a solid are
    stabilised by the lattice. The gas-phase treatment has no such term, so it over-penalises any step
    that creates or separates ions. Until a solid-state reference for the reagents is added, the
    no-solvent treatment is read only for steps whose by-product leaves as a neutral molecule (water,
    a carboxylic acid).
  - Local heating and pressure in the mill are also absent.
  - The no-solvent row is the roughest of the three.

pH, base, counterion and modulator enter in all three as species in the same enumeration.
Temperature and concentration enter later, through free-energy and activity terms.

**One system to test all four.** Copper(II) carboxylate frameworks are the natural validation case.
The paddlewheel node pre-exists in copper acetate (H0). Mononuclear intermediates are observed (H1).
The copper salt changes the product (H2). The same framework is reported from solution, wet milling
and dry milling (H3). Copper(II) is also free of the oxidation-state ambiguity of the nickel
system studied so far. It needs polynuclear nodes, which are in progress (§8).

### Hence: MolCluster
MolCluster models the step **before** the topology — **the molecular chemistry of node formation in
solution**. Given the ingredients of a synthesis, it asks:
> Which metal–ligand clusters (candidate SBUs and their precursors) can form, by which routes, and
> which recipe makes the route to a target one smoothest?

This turns the unpredictable part of synthesis into an explicit, auditable graph of species and
reactions:
- **Every species is identified exactly.** Bridging vs chelating, protonation state and hydration
  number are part of identity.
- **Every step is a balanced reaction**, including the proton and counterion bookkeeping that
  decides real syntheses. Energies are priced with a stated level of theory and refused when the
  equation would mislead.
- **Conditions are inputs:** salt anion, solvent, pH and modulator. The first study already
  recovered a textbook synthesis rule from first principles: the acetate salt succeeds largely
  because acetate is a built-in base. It also showed that with a common base the chloride complex
  forms more readily.

**The link to the existing field is direct.** Topological generators *consume* SBUs; MolCluster
*produces* them, together with the conditions under which they are accessible. Coupled, the two give
**synthesis-aware screening**: hypothetical frameworks whose nodes are reachable from real
ingredients, ranked beside their predicted properties.

### Why now
- **MLIPs trained on large DFT datasets** (MACE-OMOL-0 / OMol25, ωB97M-V) relax a coordination
  species in about 5 s at near-DFT geometry quality. That makes it affordable to enumerate thousands
  of candidate species and hydration states, a combinatorial space that DFT alone could not cover.
- **The registry and identity layer** keep that enumeration from becoming noise. Each species is
  stored once, every route is recorded, and every number carries its provenance, so screening
  results can be escalated to DFT and checked against experiment (§5).

### What it does not claim
- It does not model nucleation, crystal growth or framework crystallisation kinetics.
- It predicts the **solution-phase precursor chemistry**: which species dominate and how they
  interconvert. That is a prior on which conditions favour a target node, not a guarantee of a
  crystal.
- It covers labile metals, where thermodynamics governs (§2.6). Inert centres are a stated boundary.

---

## 0. The pitch

MolCluster takes the **ingredients of a coordination synthesis** — a metal salt, ligands in the
form they take in the medium, a solvent — and returns **the species that can form and the routes
between them**, each priced with a stated level of theory and error bar.

It runs as two lines of work with opposite cost profiles:

| | Line 1 — structure generation, classification, identification | Line 2 — energy evaluation and pathway determination |
|---|---|---|
| question | what can exist, and what *is* each thing? | which of those matter, and how do they interconvert? |
| cost per structure | milliseconds to seconds | seconds (screen) to hours (refine) |
| scales with | the size of the enumeration | the number of questions asked |
| bound by | **storage and bookkeeping** | **compute** |
| role | the product's memory | the product's judgement |

The lines meet at the registry: Line 2 reads only species Line 1 has identified, and writes numbers
back against them.

Two further layers close the loop:
- **Learning (§4)** turns what Line 1 stores into compute Line 2 doesn't have to spend. The registry
  is training data for models that say *what to compute next*, not models that replace computing.
- **Experiment (§5)** is the only external judge. Every claim either line makes is scored against a
  measured observable.

```
ingredients ─► 1.1 ligands ─► 1.2 sites ─► 1.3 spheres ─► 1.4 identity ─► 1.5 registry ◄──┐
                                                                              │            │
                                                                              ▼            │
routes ◄─ 2.6 routes ◄─ 2.5 calibrate ◄─ 2.4 refine ◄─ 2.3 select ◄─ 2.2 screen ◄─ 2.1 relax
   │                        ▲                           ▲                                   │
   │                        │                           └──── 4 learning (triage, Δ) ◄──────┘
   └──► 5 experiment ───────┘  (calibration, validation)
```

## 1. Measured demand today

Two runs (NiCl₂ and Ni(OAc)₂ with tHQ, ≤1 deprotonation) into one registry.

| | Line 1 (storage) | Line 2 (compute) |
|---|---|---|
| volume | 1,123 structures · 15,722 geometries · 19,762 reaction edges · 126,235 site-state rows | 14,766 relax tasks, 3,600 MACE geometries; 3,600 continuum corrections |
| size / time | 97 MB database + 234 MB blob store (56k blobs, ~4 KB each) → **~90 KB database per identified structure** | construction 2.2 CPU-h (place median 0.8 s, grow 0.16 s); **MACE 5.0 h (~5 s each)**; xTB-ALPB screen ~7 s per geometry |
| reuse | builds landing on an existing identity write nothing new | 1,364 relaxations skipped as already computed (acetate run) |
| concurrency | 32 workers, 0 lock failures | 8 workers, 1 GPU; wall time 2 h 24 min + 1 h 19 min |

**Scaling.**
- **Line 1** grows with ligands × protomers × co-ligand window × CN × geometry. The primary runs (tHQ
  ≤2 deprotonations) are 9–11k tasks each, ~5× the fallback. At 10⁵ structures the database is
  ~10 GB: SQLite copes, but the single writer and the site-state row count are the pressure points.
- **Line 2** grows with what is **selected**. The screen is linear in geometries at ~12 s each; DFT
  is spent on ~15 species per route pair (~45 jobs). Cost follows the questions, not the space.

---

## 2. Line 1 — structure generation, classification, identification

### 1.1 Ligands and protonation states
**MolCluster.** SMILES entered as the species exists in the medium. Protomers are enumerated over
labile sites, deduplicated by identity, and named by the parent's symmetry orbit — tHQ's 16 subsets
become 7 protomers with no hand-written symmetry.

| Existing | Approach | Relative to MolCluster |
|---|---|---|
| Dimorphite-DL [1] | rule-based protonation states in a pH window | pH-aware; MolCluster is not yet, and enumerates all states up to *k* instead |
| Epik [2] (commercial) | predicted pKa → populated states | predicts populations; MolCluster computes them later, in Line 2 |
| RDKit / OpenEye enumerators | tautomers | different target; tautomers are a **gap** here |
| molSimplify ligand library [3] | curated ligands with charges | fixed library; MolCluster derives charge states from the SMILES |

**Benchmark.** Protomer sets against Dimorphite-DL on 50 tmQM ligands [15]; orbit counts against
hand-derived values for symmetric polyacids.

### 1.2 Binding-site perception
**MolCluster.** Resonance-invariant donor groups (a carboxylate is one group); a site frame
(origin, axis, reference); lone-pair lobes; torsion as a discrete well. A geometry-free catalog plus
a per-geometry state.

| Existing | Approach | Relative to MolCluster |
|---|---|---|
| molSimplify [3] | connecting atoms from a library or the user | manual; MolCluster perceives |
| Architector [4] | coordinating atoms given with the input | manual; far wider element coverage (s–f block) |
| DENOPTIM [5] | attachment-point classes on fragments | the closest in spirit; no lone-pair or torsion model |
| CSD / Mogul | empirical coordination statistics | data, not perception — the validation source |

**Known gaps:** pyrazolate N typing, oxide vs hydroxide, steric blocking of sp³ amines.
**Benchmark.** Donor-set recall and precision against CSD Ni(II)/Cu(II) complexes: are the donors
seen coordinated in crystals among the donors perceived?

### 1.3 Coordination-sphere construction — the model generator
**MolCluster.**
- Polyhedral vertex assignment, with bites limited to 55–115°.
- Reserved vacancies under the co-ligand window. The window's ±2 waters **are the explicit
  first-shell solvent**.
- Frame-onto-frame joins, choosing lobe, well and roll on a deterministic grid; every choice is
  kept in a replayable choice vector.
- Whole-sphere placement and one-ligand growth. The growth ladder *is* the route graph.

| Existing | Scope | Relative to MolCluster |
|---|---|---|
| molSimplify [3] | mononuclear TMCs, high-throughput | more mature geometry; no route graph |
| Architector [4] | mononuclear, s–f block, xTB conformer sampling | the strongest geometry comparator; no provenance or routes |
| Molassembler [6] | any polyhedron, stereopermutations | best at stereo; MolCluster's isomer level is behind |
| AaronTools [7] | ligand swaps on TS templates | catalysis and TS focus |
| stk [8] | cage and supramolecular assembly | polynuclear assembly without chemical donor perception |
| PORMAKE [9], ToBaCCo [10], AuToGraFS [11] | MOFs from SBUs + topology | *consumers* of SBUs; MolCluster aims to *produce* them |

**Positioning.** MolCluster does not win on breadth (Architector) or stereo (Molassembler). What
is different is that every construction is an **edge** — recorded, replayable and priceable — and
hydration states come from the window. Identity is polynuclear-native; construction has so far been
demonstrated on mononuclear Ni.

**Benchmark.**
1. **Coverage:** 100 CSD/tmQM Ni(II) complexes with supported ligands. Is the crystal isomer
   produced within 0.5 Å heavy-atom RMSD after relaxation? Run side by side with Architector and
   molSimplify.
2. **Yield per rejection code.**
3. **Byte-identical replay.**

### 1.4 Classification and identity
**MolCluster.**
- **The typed graph:** explicit H; covalent, dative and M–M edges; charge on the graph; bond order
  not hashed.
- **L0:** composition plus per-centre oxidation state.
- **L1:** sha256 of a canonical certificate. The WL hash is only a bucket, because 1-WL cannot tell
  µ2-bridging from chelating.
- **L2:** isomer tag. **L3:** choice vector plus 0.15 Å RMSD.
- Labels are derived for display and never searched.

| Existing | Approach | Relative to MolCluster |
|---|---|---|
| InChI 1.07 [12] | canonical string; inorganic and organometallic support new in 1.07, molecular inorganics as a prototype | interoperable standard; historically lossy for dative and bridging |
| RDKit canonical SMILES | canonical ranking, dative bonds allowed | not stable for organometallics across versions |
| nauty / Traces [13] | canonical labelling — the gold standard | exact and faster; MolCluster's pure-Python version is to be **checked against it** |
| WL graph hash | 1-WL colour refinement | cannot separate µ2 from chelate |
| tmQM / tmQMg [15, 16] | 86k TMCs; graph representations | a ready test corpus |
| MOFid / MOFkey [17] | MOF → SBUs + topology | framework-level identity; MolCluster's is SBU-level — the two can be joined |
| Molassembler [6] | graph + stereopermutation identity | the only comparator with stereo in identity; a candidate L2 engine |

**Benchmark.**
1. **Invariance:** 10⁴ random relabellings of tmQMg graphs give one L1 each.
2. **Discrimination:** a hand set (µ2/chelate, κ¹/κ², protomers, linkage isomers), with a
   collision table per method (L1, InChI, SMILES, WL).
3. **Correctness:** L1 equivalence agrees with nauty across tmQMg.
4. **Speed:** ms per structure against nauty.
5. **L2:** cis/trans and fac/mer against Molassembler.

### 1.5 Registry, provenance and search
**MolCluster.**
- SQLite plus a content-addressed blob store, with one write surface.
- One row per species: composition, connectivity and isomer together are unique.
- A reaction DAG in which one species is one node with many incoming edges.
- Versioned algorithms, soft delete, structural filters.

| Existing | Approach | Relative to MolCluster |
|---|---|---|
| SCINE Database / Chemoton [18] | MongoDB of structures, calculations, elementary steps | the **closest analogue**; server-based and built for exploration with transition states |
| QCArchive [19] | result store with provenance, distributed compute | scales further; molecule-level, no identity hierarchy |
| AiiDA [20] | full provenance graph of calculations | stronger provenance of *calculations*; MolCluster records provenance of *chemistry* |
| ASE db | structures + energies | no identity or reactions |
| CSD | curated experimental structures | the validation source, not a peer |

**Benchmark.**
1. **Scaling:** synthetic registries at 10⁴, 10⁵ and 10⁶ structures — bytes per structure, filter
   and route query latency, write throughput against worker count. Baseline: 32 workers, 0 lock
   failures, 141 s.
2. **Deduplication yield.**
3. A decision point on moving from SQLite to Postgres, taken on the 10⁶ numbers.

---

## 3. Line 2 — energy evaluation and pathway determination

The governing idea: **the machine-learned potential decides what to compute; it does not compute the
answer.** It relaxes, ranks within an identity, triages across routes, and accelerates Hessians and
sampling. Reported numbers come only from the refinement stage (2.4).

### 2.1 Relaxation and first evaluation — the MLIP accelerator
**MolCluster.** MACE-OMOL-0 [21] (trained on OMol25 [22], ωB97M-V/def2-TZVPD) relaxes every
construct today, ~5 s on the GPU. After each relaxation the bonding is re-derived from the
coordinates and compared with the intended structure, so a relaxation that breaks or forms a bond
is flagged and kept out of the energies.

A model-comparison study against DFT has been run on six Ni species:
- **MACE-MH-1** was rejected as the relaxer. It is blind to charge and spin, and DFT prefers
  OMOL-0's Ni minima.
- **MACE-POLAR-1-L** is DFT-equivalent to OMOL-0 on five of the six and better on the sixth. Run as
  a small → large funnel it costs ~0.8× all-OMOL-0. Whether it replaces OMOL-0 waits on the
  geometry validation in 2.4.

**Planned:** a partial-relaxation funnel that stops unpromising relaxations early.

| Existing | Coverage | Relative to MACE-OMOL-0 |
|---|---|---|
| GFN2-xTB [23] | whole periodic table, continuum solvation | the semi-empirical baseline (Architector uses it); relaxes a MACE Ni minimum by −1.9 eV |
| g-xTB [24] | H–Lr, targets ωB97M-V at TB cost; ~half GFN2's errors | the direct challenger to test |
| UMA [25] | OMol25-trained, multi-domain | the sibling model; benchmark head to head first |
| AIMNet2 [26] | 14 main-group elements; no transition metals (a Pd variant exists) | out of scope for Ni |
| ANI-2x, MACE-OFF | organic | out of scope |

**Benchmark.** MACE-OMOL-0, MACE-POLAR-1, MACE-MH, UMA, GFN2-xTB and g-xTB against ωB97M-V on the escalated species,
MOR41 [27], WCCR10 [28] and a tmQM Ni subset. Metrics:
- geometry RMSD;
- **Spearman ρ of conformer ranking** — the property MACE is actually used for;
- reaction-energy MAE;
- seconds per structure.

### 2.2 Screening in the medium
**MolCluster.**
- MACE gas + xTB-ALPB [29] ΔG_solv on the same geometry.
- The co-ligand window supplies explicit first-shell water.
- Isodesmic and charge-separation checks refuse or flag.
- A proton sink turns deprotonation into a proton exchange with a reference acid.

| Existing | Relative to MolCluster |
|---|---|
| SMD [30], CPCM, COSMO-RS [31] | standard continua; SMD is the refine-tier target |
| cluster-continuum [32, 33] | the same idea; MolCluster generates the explicit shell rather than placing it by hand |
| pKa methodology [34], SAMPL challenges [35] | the error levels to expect: 1–2 pKa units relative, far worse direct |

**Where it stands:**
- Direct deprotonation to water: free-tHQ pKa ≈ 25.
- Proton exchange with acetate: +0.25 eV. Its experimental reference is pending a sourced tHQ pKa.
- Status: **screening only**.

**Benchmark.** A pKa calibration set, as RMSE in pKa units for three schemes: direct (H₃O⁺),
relative (proton sink) and cluster-continuum (H₃O⁺(H₂O)₃); plus the ALPB–GBSA–SMD spread per
species.

### 2.3 Selection — what earns compute
**MolCluster.** A candidate route escalates as a closed equation: every term, the reference
species, and the top-k conformers. One recipe is pinned per equation.

| Existing | Relative to MolCluster |
|---|---|
| Chemoton exploration heuristics [18] | autonomous choice of what to explore next |
| active learning for TMCs (Kulik group) | ML-guided DFT selection |
| Δ-learning [36] | correct a cheap model with a few expensive points |

MolCluster's selection is **question-driven** — routes a user walks, or the top-N by screen. That is
auditable, but not exhaustive.
**Benchmark.** *Selection recall:* on a small system refined exhaustively at DFT, do the top-N
screened routes contain the top-M DFT routes?

### 2.4 Refinement — DFT and thermochemistry
**MolCluster (planned).**
- DFT single point at the MACE geometry. At OMol25's level, E_DFT − E_MACE is the per-species MACE
  error.
- SMD single point on the same geometry.
- A MACE Hessian → quasi-RRHO G [37].
- Spin-state and oxidation-state checks.
- Backend: **GPU4PySCF** (decided 2026-09-29), in-process on the workstation GPU.

| Existing | Relative to MolCluster |
|---|---|
| ORCA, Gaussian, Q-Chem, Psi4, PySCF / GPU4PySCF | consumed, not replaced |
| quasi-RRHO [37], GoodVibes [38], xTB single-point Hessians | thermochemistry; MolCluster's contribution is the MACE Hessian as accelerator |
| r²SCAN-3c [39], ωB97X-3c [40], DLPNO-CCSD(T) | cheaper recipes, and reference points |
| SSE17 [41] | experimentally derived spin-state energetics, Ni²⁺ included |

**Benchmark.**
1. G_RRHO from the MACE Hessian against a DFT Hessian on ~10 species.
2. DFT at the MACE geometry against the DFT minimum — the accelerator claim.
3. Ni(II) spin ordering on the route species.
4. DLPNO-CCSD(T) on 3–5 small equations to bound the DFT error.

### 2.5 Calibration and conditions
**MolCluster (planned).** Standard state, water activity and pH as named terms; a linear fit of
computed ΔG against experiment gives the reported uncertainty.

| Existing | Relative to MolCluster |
|---|---|
| NIST SRD 46 [42], IUPAC SC-Database | the reference constants |
| HySS [43], PHREEQC [44], Visual MINTEQ | speciation codes that **consume** measured constants; MolCluster **computes** constants for species that have none |

**Benchmark.**
- Ni(II) stepwise log β for Cl⁻, acetate and one catecholate against NIST 46 (RMSE, log units).
- The first hydrolysis of Ni(H₂O)₆²⁺.
- The **Irving–Williams order** [45] (Mn < Fe < Co < Ni < Cu > Zn) for one ligand across metals — the
  qualitative sanity check a reviewer expects.

### 2.6 Pathway determination and reporting
**MolCluster.**
- A route graph of construction edges, deprotonation edges and derived splits, walked hop by hop.
- Reverse legs compose exchanges.
- A basis token decides when two routes may be compared.
- The net equation is judged whole.
- Every number carries its recipe and caveats.
- Next: speciation from priced equilibria. No transition states (a barrier proxy is deferred).

| Existing | Approach | Relative to MolCluster |
|---|---|---|
| Chemoton + KiNetX [18] | automated elementary steps with TSs, then kinetics | **far more complete kinetics**; MolCluster is thermodynamic, curated, composition-level |
| AutoMeKin [46], YARP [47], GSM [48] | automated reaction discovery with TSs | kinetics-first; mostly organic or gas phase |
| RMG [49] | rate-based mechanism generation | gas and combustion |
| energetic span model [50] | turnover-determining states | a candidate for the barrier proxy |
| HySS / PHREEQC | speciation from constants | MolCluster's intended output form |

**Domain boundary — state it.** Ni(II) exchanges water at ~3×10⁴ s⁻¹ [51], so ligand substitution is
under thermodynamic control, and composition-level routes plus speciation answer the synthesis
question without transition states. For inert centres (Cr(III), Co(III)) that fails, and
Chemoton-style kinetics would be needed.

**Benchmark.**
1. **End to end:** aqueous Ni(II)–acetate speciation against pH, from ingredients alone.
2. **The salt question** under a common proton acceptor (screen: Cl route −0.66 eV, OAc route −0.08
   eV), re-priced at the refine tier.
3. **Sign stability** across recipes (ALPB/GBSA/SMD; MACE/DFT).

---

## 4. Learning layer — the registry as training data

**The idea.** Every structure Line 1 identifies, every construction it rejects and every equation
Line 2 prices is a labelled example. A model trained on them can **nudge computation toward viable
routes** — which species to relax, which routes to refine — and, with the right labels, predict a
little beyond what was computed. Line 1's storage cost pays for Line 2's compute savings.

**The limit, stated first.** MACE-OMOL-0 and UMA are already the "large model", trained on ~10⁸
ωB97M-V calculations (OMol25). A model trained on MACE-labelled registry rows is a distillation of
MACE onto a narrower domain: faster, never more accurate. Prediction *beyond* explicit calculation
needs labels MACE lacks — DFT from refinement (2.4) and experiment (§5).

### What the registry already holds

| Asset | Salt study | Learnable signal | Limitation |
|---|---|---|---|
| nodes | 1,123 typed graphs, L0–L2, charge, spin | species energy, stability, hydration preference | one metal, three ligands |
| geometries | 15,722 (3,600 MACE) | conformer priors: which choice vector reaches the minimum | near-duplicates |
| energies | MACE + ALPB | — | labels at screening fidelity |
| edges | 19,762 with roles, stoichiometry and atom maps | reaction energies; a condensed graph of reaction is directly buildable | ΔE depends on proton sink, medium and recipe |
| task outcomes | 33,829 with rejection codes | construction feasibility | partly the constructor's limits, not chemistry |

### Models, ranked by value against risk

| # | Model | Labels | Data needed | Use | Existing comparator |
|---|---|---|---|---|---|
| L1 | **Δ-learning** E_DFT − E_MACE and ΔG_SMD − ΔG_ALPB, with uncertainty (GP or ensemble on RAC descriptors or graph kernels) | refine stage | 10²–10³ | corrects screening numbers; a "MACE trust" map (e.g. redox non-innocent THQ) | Δ-ML [36]; RACs for TMCs [52] |
| L2 | **Active-learning escalation** — acquisition (expected improvement or UCB) over the route graph from screen mean plus L1 uncertainty; replaces the fixed selection rule | refine stage, iteratively | grows as it runs | **the nudge**: what to refine next | uncertainty-driven TMC discovery [53]; Chemoton heuristics [18] |
| L3 | **Construction-feasibility classifier** — predicts clashes, refused placements and chelates that cannot span, before building | task outcomes | 10⁴ already | pruning at polynuclear scale, where enumeration explodes | — (no comparator records rejections) |
| L4 | **Geometry-free species-energy surrogate** — a GNN on the typed graph | MACE + corrections over a diverse generated corpus | 10⁵ | pruning 10⁶-scale enumerations before relaxing | tmQMg GNNs [16]; Chemprop [54] |
| L5 | **Foundation fine-tune** — MACE-OMOL or UMA fine-tuned on refine-stage DFT | refine stage | 10²–10³ | the realistic "larger model": better screening *for this chemistry* | MACE foundation fine-tuning [55] |
| L6 | **Experimental priors** — "is this species observed / stable in solution?" | CSD existence [56]; NIST log β [42] | 10³–10⁵ | the only prediction genuinely beyond calculation | — |

**Not planned:** a from-scratch large model trained on the registry to replace calculation. The data
is too narrow, the labels are the wrong fidelity, and OMol25-scale models already fill that niche.

### Rules it inherits from the design
- **A learned number is its own tier.** The fidelity ladder already has a rung for numbers with
  no structure-specific compute. A learned number carries a record of the model version and
  training-set hash. The pinned-theory rule refuses any equation mixing learned and computed terms.
  Learned numbers are never reported as results.
- **Species-level, not edge-level, energies** for thermodynamics. Edge ΔE then follows by Hess's
  law, so every cycle closes and balance holds by construction. An edge-level (CGR) model can
  predict cycles that don't close; keep those for barriers later.
- **Split by chemistry** (ligand, metal, composition), never at random. Enumerated neighbours are
  near-duplicates, and a random split overstates accuracy.
- **Uncertainty is mandatory** for anything that steers computation.

### Data campaign
The generator is the data engine. At ~5 s per MACE relaxation, 10⁵ structures across the Mn–Zn
series and a few ligand families is ~140 GPU-h — about 6 days on one workstation GPU — and ~9 GB
of registry. The same corpus serves the Irving–Williams benchmark (§5).

### Benchmarks
| Model | Metric | Baseline |
|---|---|---|
| L1 | held-out MAE of corrected screening vs DFT; calibration of the uncertainty (coverage of 1σ/2σ) | uncorrected screening |
| L2 | **selection recall** — top-M DFT routes found per DFT job spent | the fixed selection rule |
| L3 | precision and recall of rejection; construction time saved | build everything |
| L4 | MAE vs MACE on held-out metals and ligands; enrichment of top-k | random pruning |
| L5 | MAE vs DFT on held-out species vs MACE-OMOL-0 | MACE-OMOL-0 |
| L6 | AUROC for CSD existence; log β RMSE | a composition-only model |

**First step (cheap, useful now):** a dataset exporter. Nodes, edges with atom maps and roles,
outcomes with rejection codes and energies per recipe, split by chemistry, with a datasheet. It also
makes the registry citable as a data contribution.

---

## 5. Experimental benchmarking

Computation benchmarked only against more computation can be consistently wrong. Each claim below is
paired with a measured observable, the model output it is scored against, and where the data
comes from. **Status:** *exists* = published data to collect; *source* = the value still needs a
citation; *measure* = needs a wet-lab experiment.

### 5.1 Structure and identity (Line 1)
| Observable | Data | Compared with | Metric | Status |
|---|---|---|---|---|
| coordination geometry and M–L distances of Ni(II) complexes | CSD [56] | relaxed geometries of the same species | heavy-atom RMSD; M–O distance error | exists |
| Ni(H₂O)₆²⁺ structure (Ni–O ≈ 2.05–2.07 Å, CN 6) | X-ray/neutron diffraction, EXAFS [57] | the aqua ion as built and relaxed | Ni–O error; preferred hydration number | exists |
| which isomer crystallises (cis/trans, fac/mer) | CSD | lowest-energy L2 isomer | hit rate | exists (needs working L2) |
| donor atoms actually used by each ligand | CSD | perceived donor sets | recall and precision | exists |
| metal–THQ coordination motifs (bridging vs chelating) | CSD survey of metal–THQ and related polymers | predicted dominant motif | agreement | **source** (survey to do) |

### 5.2 Solution thermodynamics (Line 2)
| Observable | Data | Compared with | Metric | Status |
|---|---|---|---|---|
| pKa of acetic acid, phenol, catechol, H₃O⁺ | IUPAC / standard compilations | proton-exchange ΔG (three schemes) | RMSE in pKa units | exists |
| pKa₁, pKa₂ of tHQ | literature (the analogue DHBQ measures 2.95 / 5.25) | same | error | **source** |
| stepwise log β of Ni²⁺ with Cl⁻, acetate, catecholate/salicylate | NIST SRD 46 [42] | ligand-exchange ΔG | RMSE in log units | exists |
| first hydrolysis constant of Ni(H₂O)₆²⁺ | Brown & Ekberg [58] | deprotonation of the aqua ion | error in log units | exists |
| absolute hydration free energy of Ni²⁺, Cl⁻ | Marcus [59] | solvation scheme on bare ions | error in eV | exists |
| Irving–Williams order (Mn < Fe < Co < Ni < Cu > Zn) | Irving & Williams [45] | the same ligand across the Mn–Zn series | rank order | exists (needs the multi-metal corpus) |

**Comparing like with like.**
- **Temperature:** constants are reported at 25 °C.
- **Ionic strength:** often I = 0.1–1 M, so correct to I = 0 with the Davies equation [60], or model
  the supporting electrolyte.
- **Standard state:** 1 M solutes, pure liquid water.
- **Hydrolysis:** it competes at high pH, so compare over the stated pH range.

### 5.3 Electronic structure (Line 2)
| Observable | Data | Compared with | Metric | Status |
|---|---|---|---|---|
| Ni(II) spin ground state (octahedral triplet; μ_eff ≈ 2.9–3.3 μB) | magnetic susceptibility literature; SSE17 [41] | spin ordering at the refine stage | correct ground state; splitting error | exists |
| d–d band positions of Ni(II) species (UV-vis) | literature; own measurements | TD-DFT, or ligand-field ordering of predicted species | band shift direction on substitution | exists / measure |
| metal oxidation state in THQ complexes | EPR, XANES where reported | spin density on Ni | agreement | **source** |

### 5.4 The project's own question — proposed experiments
These test the salt result (acetate vs chloride as the Ni source for THQ complexation) directly:
1. **Potentiometric (pH-metric) titration** of Ni²⁺ + THQ in NiCl₂ and in Ni(OAc)₂ media at fixed
   ionic strength. Fit with Hyperquad [61] to get log β and the species distribution against pH, and
   compare with the predicted speciation. This is the decisive experiment.
2. **UV-vis titration** of the same systems. The appearance of the complex band against pH is an
   independent speciation readout.
3. **Isothermal titration calorimetry** gives ΔH and TΔS separately, which tests the entropy terms
   that screening lacks.
4. **Synthesis outcome.** Does the solid obtained (PXRD, single-crystal structure) contain the SBU
   predicted as the dominant solution species? It's a weaker link, since crystallisation selects, but
   it is the question MOF synthesis actually asks.

### 5.5 Kinetic domain check
Water-exchange rates [51] set where thermodynamic control holds: fast for Ni(II), slow for Cr(III)
and Co(III). A benchmark system with an inert centre (for example Co(III) ammines) is a deliberate
**negative control**. There, the thermodynamic prediction should *disagree* with the observed
kinetic product, and the plan should say so rather than claim coverage.

---

## 6. Where MolCluster sits

| Step | Best existing | MolCluster |
|---|---|---|
| protomers | Epik, Dimorphite-DL | exhaustive, symmetry-named; pH to come |
| site perception | mostly manual elsewhere | **automatic, frame-based — differentiator** |
| sphere construction | Architector, molSimplify | adequate; the difference is recorded choices as edges |
| stereo (L2) | Molassembler | **behind** — adopt or benchmark |
| identity | nauty, InChI | **exact for bridging and dative — differentiator**, to verify against nauty |
| registry + provenance | SCINE, AiiDA, QCArchive | provenance of chemistry (routes as edges) — differentiator |
| MLIP | MACE-OMOL, UMA, g-xTB | consumer — benchmark, don't build |
| solvation | SMD, COSMO-RS, cluster-continuum | standard, with the explicit first shell generated |
| DFT + thermo | ORCA / PySCF + qRRHO | consumer; MACE Hessians as accelerator |
| pathways | Chemoton (kinetics) | thermodynamic, curated; **refusal of non-isodesmic and charge-separating equations — differentiator** |
| speciation | HySS, PHREEQC | computes the constants they take as input |
| learning | Δ-ML, Chemprop, MACE fine-tuning, TMC active learning | consumer of methods; **the data — balanced edges with recorded rejections — is the differentiator** |
| validation | CSD, NIST SRD 46, titration | uses standard data; adds its own titration of the salt question |

**In one line:** the contribution is integration, identity and refusal discipline — not a new
method at any single step. Every step that is not a differentiator uses, or benchmarks against, the
field's standard tool.

## 7. Benchmark roadmap

| # | Benchmark | Line | Cost | Establishes |
|---|---|---|---|---|
| 1 | identity: invariance, discrimination, nauty agreement on tmQMg | 1 | CPU-hours | the identity claim |
| 2 | pKa calibration, three schemes | 2 | minutes (screen) | honest error bars |
| 3 | MACE vs UMA / GFN2 / g-xTB vs DFT on route species | 2 | ~45 DFT single points | the accelerator claim |
| 4 | construction coverage vs Architector / molSimplify on CSD Ni(II) | 1 | a day | the generator claim |
| 5 | MACE Hessian vs DFT Hessian G_RRHO | 2 | ~10 DFT frequency jobs | free energies |
| 6 | Ni(II) log β vs NIST 46; Irving–Williams | 2 | ~a week of DFT | the chemistry claim |
| 7 | registry scaling to 10⁶ | 1 | a day, synthetic | the storage claim |
| 8 | end-to-end Ni–acetate speciation | 1 + 2 | after 1–6 | the project claim |
| 9 | dataset exporter + chemistry-split datasheet | 4 | a day | learning is possible at all |
| 10 | Δ-learning with calibrated uncertainty (L1) | 4 | after 3 | the "MACE trust" map |
| 11 | active-learning selection recall vs the fixed rule (L2) | 4 | after 10 | the nudging claim |
| 12 | CSD geometry and Ni aqua-ion structure (5.1) | 5 | days | structures match experiment |
| 13 | pKa, log β, hydrolysis, hydration energies (5.2) | 5 | with 2 and 6 | numbers match experiment |
| 14 | Ni–THQ titration in chloride vs acetate media (5.4) | 5 | wet lab, weeks | **the project's question, answered by experiment** |

Benchmark 3 has started (six Ni species, three MACE variants against DFT). Next in reach are 1, 2
and 9. Benchmark 14 needs a collaborator with a titrator. Status of every item is in §8.

---

## 8. Progress ledger — done against planned

Status on 13 October 2026. **Done** = built and in use; **Partly** = built in part, with the
remainder named; **In progress** = being built now; **Planned** = not started.

### Line 1 — structure generation, classification, identification
| Step | What was planned | Status | What remains |
|---|---|---|---|
| 1.1 | ligand input from SMILES; protomers enumerated and named by symmetry | **Done** | pH-aware enumeration; tautomers |
| 1.2 | automatic binding-site perception with frames, lone pairs and torsion wells | **Done** | three known typing gaps (pyrazolate, oxide, sp³ amines) |
| 1.3 | mononuclear coordination spheres: whole-sphere placement, one-ligand growth, reserved vacancies | **Done** | shown on Ni(II), tetrahedral and octahedral |
| 1.3 | explicit first-shell solvent through the co-ligand window | **Done** | — |
| 1.3 | polynuclear nodes (bridged and metal–metal bonded) | **In progress** | bridging and metal–metal joins are built; validated full nodes are not |
| 1.4 | identity by composition and exact connectivity | **Done** | verification against nauty (benchmark 1) |
| 1.4 | isomer level (cis/trans, fac/mer) | **Partly** | the classifier exists; production runs do not apply it yet |
| 1.4 | conformer level | **Done** | the energy window is still open |
| 1.5 | registry with provenance: one species once, every route an edge | **Done** | scaling test (benchmark 7) |
| 1.5 | each geometry linked atom-by-atom to its structure | **Done** | — |
| — | interface: recipe builder, run inspector, registry browser, energy-level graph | **Done** | — |

### Line 2 — energy evaluation and pathway determination
| Step | What was planned | Status | What remains |
|---|---|---|---|
| 2.1 | relaxation of every construct with a machine-learned potential | **Done** | — |
| 2.1 | bonding re-checked after each relaxation | **Done** | — |
| 2.1 | comparison of candidate potentials against DFT | **Partly** | done on six Ni species; the choice of relaxer is open; UMA, g-xTB and external sets remain |
| 2.1 | partial-relaxation funnel | **Planned** | — |
| 2.2 | screening in a continuum (water) on the same geometry | **Done** | error bar from a second continuum; pH and water-reference options; a flag for hydration changes |
| 2.2 | proton exchange with a reference acid | **Done** | — |
| 2.2 | refusal or flagging of unbalanced, non-isodesmic and charge-separating equations | **Done** | — |
| 2.3 | selection of whole equations for refinement | **Planned** | — |
| 2.4 | DFT refinement | **Partly** | single points with spin and oxidation-state checks are built; DFT-level solvation and the spin ladder are not |
| 2.4 | thermochemistry from machine-learned Hessians | **Planned** | — |
| 2.5 | calibration against pKa and stability constants; standard state, pH | **Planned** | the pKa of tetrahydroxyquinone still needs a source |
| 2.6 | routes walked hop by hop, forward and reverse, priced as whole equations | **Done** | — |
| 2.6 | automatic route enumeration and ranking | **Planned** | needs a barrier proxy |
| 2.6 | speciation from priced equilibria | **Planned** | — |
| — | hydrogen-bonded solvent shells for free ions and the proton | **Planned** | — |

### Studies and hypotheses
| Item | Status | Note |
|---|---|---|
| Ni–tetrahydroxyquinone, chloride vs acetate salt | **Done at screening level** | re-pricing at the refinement level is planned |
| H0 solution species carry over into the framework | **Planned** | needs systems where the node is known to pre-exist, and ones where it is not |
| H1 monomer metastability | **Planned** | the window's bounds are not yet set; needs systems with known crystalline and failed outcomes |
| H2 recipe sets driving force, highest intermediate and sinks | **Partly** | one worked instance of the driving force (the proton acceptor in the salt study); the other two quantities are not yet reported |
| H3 solvation shell as synthesis method | **Partly** | all three treatments exist for the salt study; no solid-state reference, and no comparison with milling results yet |
| copper carboxylate validation system (tests H0–H3) | **Planned** | needs polynuclear nodes |

### Learning layer and experiments
| Item | Status |
|---|---|
| dataset exporter | **Planned** |
| Δ-learning, active-learning selection, feasibility classifier, surrogate, fine-tuning | **Planned** |
| benchmark 3 (potentials against DFT) | **In progress** |
| benchmarks 1, 2, 4–13 | **Planned** |
| benchmark 14 (titration in chloride and acetate media) | **Planned** — needs a wet-lab partner |

**Reading the ledger.** Line 1 is largely built for mononuclear species; its open front is
polynuclear nodes and the isomer level. Line 2 is built up to screening. Every number produced so
far is a screening number, and the first reportable ones wait on the rest of the refinement stage
and on calibration.

---

## References

Check volume and page numbers before any external use; DOIs are given where the source was
checked for this document.

1. Ropp et al., *J. Cheminform.* **11**, 14 (2019) — Dimorphite-DL.
2. Shelley et al., *J. Comput.-Aided Mol. Des.* **21**, 681 (2007) — Epik.
3. Ioannidis, Gani & Kulik, *J. Comput. Chem.* **37**, 2106 (2016) — molSimplify.
4. Taylor et al., *Nat. Commun.* **14**, 2786 (2023) — Architector.
5. Foscato, Venkatraman & Jensen, *J. Chem. Inf. Model.* **59**, 4077 (2019) — DENOPTIM.
6. Sobez & Reiher, *J. Chem. Inf. Model.* **60**, 3884 (2020) — Molassembler.
7. Ingman et al., *WIREs Comput. Mol. Sci.* (2021) — AaronTools.
8. Turcani, Berardo & Jelfs, *J. Comput. Chem.* **39**, 1931 (2018) — stk.
9. Lee et al., *ACS Appl. Mater. Interfaces* **13**, 23647 (2021) — PORMAKE.
10. Colón, Gómez-Gualdrón & Snurr, *Cryst. Growth Des.* **17**, 5801 (2017) — ToBaCCo.
11. Addicoat et al., *J. Chem. Theory Comput.* **10**, 880 (2014) — AuToGraFS.
12. InChI 1.07 (IUPAC / InChI Trust, 2024), iupac.org/inchi-1-07-available-on-github; Blanke et al.,
    *Faraday Discuss.* (2025), doi:10.1039/d4fd00145a.
13. McKay & Piperno, *J. Symb. Comput.* **60**, 94 (2014) — nauty and Traces.
14. Heller et al., *J. Cheminform.* **7**, 23 (2015) — InChI, the worldwide chemical structure identifier standard.
15. Balcells & Skjelstad, *J. Chem. Inf. Model.* **60**, 6135 (2020) — tmQM.
16. Kneiding et al., *Digital Discovery* **2**, 618 (2023) — tmQMg.
17. Bucior et al., *Cryst. Growth Des.* **19**, 6682 (2019) — MOFid.
18. Unsleber et al., *J. Chem. Theory Comput.* (2022) — SCINE Chemoton 2.0 and its database.
19. Smith et al., *WIREs Comput. Mol. Sci.* **11**, e1491 (2021) — QCArchive.
20. Huber et al., *Sci. Data* **7**, 300 (2020) — AiiDA 1.0.
21. Batatia et al., *NeurIPS* (2022) — MACE; MACE-OMOL-0 model release (2025).
22. Levine et al., arXiv (2025) — OMol25.
23. Bannwarth, Ehlert & Grimme, *J. Chem. Theory Comput.* **15**, 1652 (2019) — GFN2-xTB.
24. Grimme group, ChemRxiv (23 June 2025), doi:10.26434/chemrxiv-2025-bjxvt — g-xTB.
25. Meta FAIR, arXiv (2025) — UMA.
26. Anstine, Zubatyuk & Isayev, *Chem. Sci.* **16**, 10228 (2025) — AIMNet2.
27. Dohm et al., *J. Chem. Theory Comput.* **14**, 2596 (2018) — MOR41.
28. Weymuth et al., *J. Chem. Theory Comput.* **10**, 3092 (2014) — WCCR10.
29. Ehlert et al., *J. Chem. Theory Comput.* **17**, 4250 (2021) — ALPB.
30. Marenich, Cramer & Truhlar, *J. Phys. Chem. B* **113**, 6378 (2009) — SMD.
31. Klamt, *J. Phys. Chem.* **99**, 2224 (1995) — COSMO-RS.
32. Pliego & Riveros, *J. Phys. Chem. A* **105**, 7241 (2001) — cluster-continuum.
33. Bryantsev, Diallo & Goddard, *J. Phys. Chem. B* **112**, 9709 (2008) — water-cluster reference.
34. Ho & Coote, *Theor. Chem. Acc.* **125**, 3 (2010) — pKa methods.
35. Işık et al., *J. Comput.-Aided Mol. Des.* (2021) — SAMPL6 pKa.
36. Ramakrishnan et al., *J. Chem. Theory Comput.* **11**, 2087 (2015) — Δ-learning.
37. Grimme, *Chem. Eur. J.* **18**, 9955 (2012) — quasi-RRHO.
38. Luchini et al., *F1000Research* **9**, 291 (2020) — GoodVibes.
39. Grimme et al., *J. Chem. Phys.* **154**, 064103 (2021) — r²SCAN-3c.
40. Müller et al., *J. Chem. Phys.* **158**, 014103 (2023) — ωB97X-3c.
41. Radoń et al., *Chem. Sci.* (2024), PMC11577268 — SSE17.
42. NIST SRD 46, *Critically Selected Stability Constants of Metal Complexes*.
43. Alderighi et al., *Coord. Chem. Rev.* **184**, 311 (1999) — HySS.
44. Parkhurst & Appelo, USGS Techniques and Methods 6-A43 (2013) — PHREEQC 3.
45. Irving & Williams, *J. Chem. Soc.* 3192 (1953).
46. Martínez-Núñez, *J. Comput. Chem.* **36**, 222 (2015) — AutoMeKin.
47. Zhao & Savoie, *Nat. Comput. Sci.* **1**, 479 (2021) — YARP.
48. Zimmerman, *J. Chem. Theory Comput.* **9**, 3043 (2013) — growing string method.
49. Gao et al., *Comput. Phys. Commun.* **203**, 212 (2016) — RMG.
50. Kozuch & Shaik, *Acc. Chem. Res.* **44**, 101 (2011) — energetic span.
51. Helm & Merbach, *Chem. Rev.* **105**, 1923 (2005) — water exchange rates.
52. Janet & Kulik, *J. Phys. Chem. A* **121**, 8939 (2017) — revised autocorrelation (RAC) descriptors.
53. Janet, Duan, Yang, Nandy & Kulik, *Chem. Sci.* **10**, 7913 (2019) — uncertainty-controlled
    discovery for TMCs.
54. Heid & Green, *J. Chem. Inf. Model.* **62**, 2101 (2022) — Chemprop condensed graph of reaction.
55. Batatia et al., arXiv:2401.00096 (2023) — a foundation model for atomistic materials chemistry
    (MACE fine-tuning).
56. Groom, Bruno, Lightfoot & Ward, *Acta Cryst. B* **72**, 171 (2016) — the Cambridge Structural
    Database.
57. Ohtaki & Radnai, *Chem. Rev.* **93**, 1157 (1993) — structure of hydrated ions.
58. Brown & Ekberg, *Hydrolysis of Metal Ions* (Wiley-VCH, 2016).
59. Marcus, *J. Chem. Soc., Faraday Trans.* **87**, 2995 (1991) — ion hydration free energies.
60. Davies, *Ion Association* (Butterworths, 1962) — the Davies activity equation.
61. Gans, Sabatini & Vacca, *Talanta* **43**, 1739 (1996) — Hyperquad.
62. Yaghi et al., *Nature* **423**, 705 (2003) — reticular synthesis.
63. O'Keeffe, Peskov, Ramsden & Yaghi, *Acc. Chem. Res.* **41**, 1782 (2008) — the RCSR.
64. Wilmer et al., *Nat. Chem.* **4**, 83 (2012) — large-scale hypothetical MOF screening (hMOF).
65. Rosen et al., *Matter* **4**, 1578 (2021) — QMOF database.
66. Chung et al., *Chem. Mater.* **26**, 6185 (2014); *J. Chem. Eng. Data* **64**, 5985 (2019) — CoRE MOF.
67. Burner et al., *Chem. Mater.* **35**, 900 (2023) — ARC-MOF.
68. Tsuruoka et al., *Angew. Chem. Int. Ed.* **48**, 4739 (2009) — coordination modulation.
69. Schaate et al., *Chem. Eur. J.* **17**, 6643 (2011) — modulated synthesis of Zr MOFs.
70. Moosavi et al., *Nat. Commun.* **10**, 539 (2019) — capturing chemical intuition in MOF synthesis.
71. Luo et al., *Angew. Chem. Int. Ed.* **61**, e202200242 (2022) — MOF synthesis prediction from mined
    data.
72. Anderson & Gómez-Gualdrón, *Chem. Sci.* **11**, 4164 (2020) — free energies toward synthetic
    likelihood.
73. Van Vleet, Weng, Li & Schmidt, *Chem. Rev.* **118**, 3681 (2018) — in situ studies of MOF
    nucleation and growth.
74. Pichon, Lazuen-Garay & James, *CrystEngComm* **8**, 211 (2006) — solvent-free mechanochemical
    synthesis of a microporous MOF.
75. Friščić et al., *CrystEngComm* **11**, 418 (2009) — liquid-assisted grinding and the role of the
    added liquid.
76. James et al., *Chem. Soc. Rev.* **41**, 413 (2012) — mechanochemistry: opportunities for new and
    cleaner synthesis.
77. Stolar et al., *Inorg. Chem.* **56**, 6599 (2017) — in situ monitoring of the mechanosynthesis of
    HKUST-1; the effect of liquid additives; mononuclear copper intermediates.
78. "Insights into HKUST-1 metal–organic framework's morphology and physicochemical properties
    induced by changing the copper(II) salt precursors" (2025), PMC11819983.
