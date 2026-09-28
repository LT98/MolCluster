# MolCluster — project plan

*Presentation-oriented: what the project does, how each step compares with existing tools, and how
each claim will be benchmarked. The code-level plans are `docs/PLAN_implementation.md` and
`docs/WORKPLAN_energy.md`. Written 2026-09-28; demand figures are measured on `salt_study_k1.db`.*

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

```
ingredients ─► 1.1 ligands ─► 1.2 sites ─► 1.3 spheres ─► 1.4 identity ─► 1.5 registry
                                                                              │
routes ◄─ 2.6 routes ◄─ 2.5 calibrate ◄─ 2.4 refine ◄─ 2.3 select ◄─ 2.2 screen ◄─ 2.1 relax
```

## 1. Measured demand today

Two runs (NiCl₂ and Ni(OAc)₂ with tHQ, ≤1 deprotonation) into one registry.

| | Line 1 (storage) | Line 2 (compute) |
|---|---|---|
| volume | 1,123 structures · 15,722 geometries · 19,762 reaction edges · 126,235 site-state rows | 14,766 relax tasks, 3,600 MACE geometries; 3,600 continuum corrections |
| size / time | 97 MB database + 234 MB blob store (56k blobs, ~4 KB each) → **~90 KB database per identified structure** | construction 2.2 CPU-h (place median 0.8 s, grow 0.16 s); **MACE 5.0 h (~5 s each)**; xTB-ALPB screen ~7 s per geometry |
| reuse | builds landing on an existing identity write nothing new | 1,364 relaxations skipped as already computed (acetate run) |
| concurrency | 32 workers, 0 lock failures, after the write-lock fix | 8 workers, 1 GPU; wall time 2 h 24 min + 1 h 19 min |

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

**Known gaps:** pyrazolate N typing (B14), oxide vs hydroxide (B15), steric blocking of sp³ amines
(B8).
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
| Molassembler [6] | any polyhedron, stereopermutations | best at stereo; MolCluster's L2 lags (B2) |
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
- `UNIQUE(l0, l1, l2)`.
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

The governing idea (`WORKPLAN_energy.md`): **MACE decides what to compute; it does not compute the
answer.** It relaxes, ranks within an identity, triages across routes, and accelerates Hessians and
sampling. Reported numbers are stage-5 recipe results.

### 2.1 Relaxation and first evaluation — the MLIP accelerator
**MolCluster.** MACE-OMOL-0 [21] (trained on OMol25 [22], ωB97M-V/def2-TZVPD) relaxes every
construct, ~5 s on the GPU. **Missing:** the connectivity check after relaxation (gate E1).

| Existing | Coverage | Relative to MACE-OMOL-0 |
|---|---|---|
| GFN2-xTB [23] | whole periodic table, continuum solvation | the semi-empirical baseline (Architector uses it); relaxes a MACE Ni minimum by −1.9 eV |
| g-xTB [24] | H–Lr, targets ωB97M-V at TB cost; ~half GFN2's errors | the direct challenger to test |
| UMA [25] | OMol25-trained, multi-domain | the sibling model; benchmark head to head first |
| AIMNet2 [26] | 14 main-group elements; no transition metals (a Pd variant exists) | out of scope for Ni |
| ANI-2x, MACE-OFF | organic | out of scope |

**Benchmark (E-1).** MACE-OMOL-0, UMA, GFN2-xTB and g-xTB against ωB97M-V on the escalated species,
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

**Benchmark (E-2).** A pKa calibration set, as RMSE in pKa units for three schemes: direct (H₃O⁺),
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
- Backend: ORCA or GPU4PySCF — **open decision**.

| Existing | Relative to MolCluster |
|---|---|
| ORCA, Gaussian, Q-Chem, Psi4, PySCF / GPU4PySCF | consumed, not replaced |
| quasi-RRHO [37], GoodVibes [38], xTB `bhess` | thermochemistry; MolCluster's contribution is the MACE Hessian as accelerator |
| r²SCAN-3c [39], ωB97X-3c [40], DLPNO-CCSD(T) | cheaper recipes, and reference points |
| SSE17 [41] | experimentally derived spin-state energetics, Ni²⁺ included |

**Benchmark (E-3/E-4).**
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

**Benchmark (E-6).**
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
- Next: speciation from priced equilibria. No transition states (barrier proxy C6 is deferred).

| Existing | Approach | Relative to MolCluster |
|---|---|---|
| Chemoton + KiNetX [18] | automated elementary steps with TSs, then kinetics | **far more complete kinetics**; MolCluster is thermodynamic, curated, composition-level |
| AutoMeKin [46], YARP [47], GSM [48] | automated reaction discovery with TSs | kinetics-first; mostly organic or gas phase |
| RMG [49] | rate-based mechanism generation | gas and combustion |
| energetic span model [50] | turnover-determining states | a candidate for the C6 barrier proxy |
| HySS / PHREEQC | speciation from constants | MolCluster's intended output form |

**Domain boundary — state it.** Ni(II) exchanges water at ~3×10⁴ s⁻¹ [51], so ligand substitution is
under thermodynamic control, and composition-level routes plus speciation answer the synthesis
question without transition states. For inert centres (Cr(III), Co(III)) that fails, and
Chemoton-style kinetics would be needed.

**Benchmark (E-7).**
1. **End to end:** aqueous Ni(II)–acetate speciation against pH, from ingredients alone.
2. **The salt question** under a common proton acceptor (screen: Cl route −0.66 eV, OAc route −0.08
   eV), re-priced at the refine tier.
3. **Sign stability** across recipes (ALPB/GBSA/SMD; MACE/DFT).

---

## 4. Where MolCluster sits

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

**In one line:** the contribution is integration, identity and refusal discipline — not a new
method at any single step. Every step that is not a differentiator uses, or benchmarks against, the
field's standard tool.

## 5. Benchmark roadmap

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

Feasible before an external review: 1 and 2. Benchmark 3 needs the DFT backend (E4).

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
