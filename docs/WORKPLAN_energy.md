# Energy work plan — from ingredients to priced routes, MACE as accelerator not endpoint

**Every energy slice cites this file.** It makes explicit what the other plans assumed: the
builder is a *model generator*, MACE is an *intermediate evaluator and accelerator*, and a
number reported as a result is a **recipe result** from stage 5, never a MACE number. The
explicit-solvent half of cluster-continuum is the co-ligand window the builder already has (§2).
Solvation detail stays in `WORKPLAN_solvation.md` (C16–C21); this file places it in the whole.

Delete this file when the slices in §5 land; what each settles goes to the ledger under one
D-number (D-TBD: "MACE is a screening and accelerator tier; reported numbers are recipe results").

## Why now

The ion/proton check of 2026-09-28 on `salt_study_k1.db` showed today's number — MACE gas energy
plus an xTB-ALPB correction, at pH 0, without thermal terms — is a screening number:

| Equation | gas | ALPB (eV) | caveat |
|---|---|---|---|
| free tHQ → tHQ⁻, H⁺ to water | +7.36 | +1.60 | charge_separation |
| free tHQ → tHQ⁻, **H⁺ to acetate** | −0.64 | **+0.25** | — |
| NiCl₂ route (#132←#94←#28→#68), H⁺ to water | +7.45 | +0.69 | charge_separation |
| NiCl₂ route, **H⁺ to acetate** | −0.55 | **−0.66** | — |
| Ni(OAc)₂ route (#328…#301), H⁺ to water | +7.95 | +1.28 | charge_separation |
| Ni(OAc)₂ route, **H⁺ to acetate** | −0.04 | **−0.08** | — |

- Direct deprotonation to water gives free-tHQ pKa ≈ 25 (1.60 / 0.059 − log 55.5). The ligand is a
  moderately strong acid. **Its experimental pKa₁ is not yet sourced**: the analogue
  2,5-dihydroxy-1,4-benzoquinone measures 2.95 / 5.25. Stage 6 needs the value before the
  free-ligand row can serve as a calibration point.
- A caveat-free equation is not a correct one: the free-ligand proton exchange moves 0.9 eV between
  gas and continuum, because a small localised anion and a large delocalised one are solvated with
  different errors.
- With the **same** proton acceptor on both routes, chloride is the easier route (−0.66 vs
  −0.08 eV). "Acetate helps" in the salt study is mostly *acetate is a base*, not a Ni–OAc vs Ni–Cl
  bonding effect.

---

## 0. The principle — MACE decides what to compute; it does not compute the answer

The generator produces thousands of candidates (salt study: 1,123 structures, 15,722 geometries, 3,600
MACE-relaxed). DFT on all of them is unaffordable and pointless — most never sit on a route
anyone asks about. So the pipeline is a **funnel**; each stage is cheaper per structure than the
next and hands it a smaller set.

| Stage | Does | Evaluator | Structures (salt-study scale) | Cost each |
|---|---|---|---|---|
| 1 Generate | enumerate, build, identity, clash QC | none (RAW) | 10³–10⁴ | ms |
| 2 Relax + triage | geometry, **connectivity check**, conformer ranking and pruning | **MACE-OMOL-0** (E-MH, §3a: MACE-MH-1 is charge-blind and DFT prefers OMOL-0's minima) | 10³ | ~5 s GPU |
| 3 Screen in medium | continuum correction, route pricing, route pruning | MACE + xTB ALPB | 10³ | ~7 s |
| 4 Select | close candidate routes' equations into an escalation set | rules | 10–50 species × top-k | — |
| 5 Refine | DFT single points, DFT solvation, thermochemistry | DFT, with **MACE Hessians and sampling** as accelerators | 10² | min–h |
| 6 Calibrate + condition | pKa / log β fit, standard state, pH, concentration | experiment | reference set | — |
| 7 Report | routes carrying recipe, conditions and error bar | — | — | — |

MACE's four jobs — none is "final energy":

1. **Geometry accelerator.** OMol25 labels are ωB97M-V/def2-TZVPD, so a MACE minimum is used as the
   DFT geometry; stage 5a measures whether that holds. Relaxation moves to MACE-MH (§3a), whose
   head and training level differ, so the argument is re-earned by study E-MH and then by 5a.
2. **Ranker within one identity.** Conformers and L2 isomers of one species: same atoms and
   charge, where systematic error cancels best.
3. **Triage evaluator.** Comparisons across species and routes at stage 3 decide *what escalates*,
   never *what is reported*.
4. **Thermochemistry and sampling accelerator.** Analytic MACE Hessians for ZPE, thermal terms and
   entropy; short MACE MD for aqua-H orientations. These are the expensive parts of a free energy
   that tolerate a cheaper model.

**Rule:** a number is at the recipe of its equation's weakest term, and only stage-5 recipes are
reported as results. Mixing tiers inside one equation stays refused (`energy/reference`
pins the first term's theory already).

---

## 1. Stage 0 — the medium is an ingredient

Ingredients are unchanged: metal salt, ligands entered as they exist in the medium (D27), the
co-ligand, CN and geometries, maximum deprotonations, and the co-ligand window (D26). **New:** the
spec names the medium once, and the continuum, the default co-ligand, the proton carrier and the
water reference derive from it. Today these are three separate choices (spec `co_ligand`, the
`/graph` medium select, the proton-sink select), which lets an equation fill vertices with water
and price in another continuum.

---

## 2. Explicit solvent — what the co-ligand already provides

**The co-ligand window is explicit first-shell solvent at identity level.** A water on a vertex
is a DATIVE edge, so `Ni(tHQ⁻)Cl(H₂O)₂` and `…(H₂O)₃` are different L1 species. The ±2 window
generates the hydration states and the grow ladder connects them. That is the explicit half of
cluster-continuum (`WORKPLAN_solvation.md`, option (b)) for **everything bound to a metal**,
including aqua→donor hydrogen bonds inside a complex, which MACE relaxation finds.

It does **not** cover:

- **Free ions and the proton carrier** (Cl⁻, acetate, tHQ⁻, H₃O⁺). They have no vertices and stay
  bare in the continuum — the worst terms in the check above. Remedy: H-bonded shells (C19–C20,
  solvation S4–S5: Cl⁻(H₂O)ₙ, H₃O⁺(H₂O)₃). Until then, price proton transfers relative to an acid
  (the proton sink).
- **The second shell** around a complex — same remedy, lower priority.

The explicit waters force bookkeeping that is missing today:

- **Hydration equilibria** (Ni(L)(H₂O)ₙ ⇌ Ni(L)(H₂O)ₙ₋₁ + H₂O) change the number of free species,
  so they need entropy (5d) and a **liquid-water reference**:
  - the *monomer* reference, H₂O(aq) + RT ln 55.5 (0.10 eV per water);
  - the *cluster* reference, (H₂O)ₙ (Bryantsev, Diallo & Goddard 2008).

  Offer both beside the proton sink, pinned and named in the net equation.
- **`hydration_change`**, a new caveat code: the count of free solvent molecules changes across the
  arrow. Screening has no TΔS (~0.3–0.5 eV per released water at 298 K), so these comparisons are
  flagged until stage 5.
- **Aqua-H orientations are L3 conformers**, sampled at stage 2 and pruned by MACE within one
  identity.

---

## 3. Stage by stage

### Stage 1 — Generate *(built)*
Construction, identity, clash QC, choice vectors (`runner`, `assembly/`, `geometry/placer.py`).

### Stage 2 — Relax and triage *(mostly built)*
- **Built:** `runner._execute_relax` → `energy/relax.relax_geometry`, reuse on (source geometry,
  method), fidelity per geometry.
- **Missing — gate E1: check connectivity after relaxation.** `_execute_relax` states the gap: a
  relaxation that forms or breaks a bond (a water drifting off, a proton moving) is filed under
  the old identity. Derive bonding from the relaxed coordinates with the pair distances in
  `geometry/distances.py`, compare with the stored graph, and on mismatch mark the geometry
  `connectivity_changed` and exclude it from energy selection (a D11 `derived_from` record comes
  later). **No number is reported before this lands.**
- **Missing — pruning within identity.** Keep conformers within 0.25 eV (~10 kT; proposal) of the
  identity minimum. The rest are kept, below the escalation line. This settles C2 for ranking.

### Stage 3 — Screen in the medium *(built)*
- MACE gas + xTB-ALPB ΔG_solv on the same geometry (`energy/solvation.correct`,
  `solvation_corrections`, C17), pinned per equation; `pathways/route.price_path` with quality
  codes and proton sink.
- **Relabelled as the screening tier** in the UI, with the ALPB–GBSA spread stored beside it as the
  error bar (one more xTB single point per geometry).
- **Route pruning:** a route more than 2 × spread above the best route to the same target falls
  below the escalation line — shown greyed, never deleted.

### Stage 4 — Select; the escalation unit is the equation
One recipe per equation means escalation cannot be per species. A candidate route escalates as a
**closed set**:
- every node, and every reagent, solvent and leaving term of each leg;
- the proton-sink couple and the water-reference species;
- the top-k conformers of each species by MACE (k = 3, proposal), so the DFT minimum does not rest on
  MACE's ordering alone.

For the two demo routes this is ~15 species and ~45 jobs. New:
`energy/escalate.escalation_set(reg, routes, k)` and a `refine` task kind queued at the lowest
priority.

### Stage 5 — Refine
A composite per species × conformer, each term its own `MethodSpec` (ground rule: no bare floats):

```
G(aq) = E_DFT(gas, g)                                 5b
      + ΔG_solv,DFT(g)          SMD / CPCM            5c
      + G_RRHO(g)               MACE Hessian, qRRHO   5d
      + standard-state, activity, pH terms — named beside ΔG, never stored   (stage 6)
```

- **5a Geometry.** MACE geometry by default, validated on ~10 species re-optimised at DFT
  (r²SCAN-3c): RMSD, and ΔE(single point at the MACE geometry − DFT minimum). Under 0.05 eV mean
  (proposal) and the accelerator claim stands. Not xTB geometries: xTB relaxes a MACE Ni minimum
  by −1.9 eV (`WORKPLAN_solvation.md` §1).
- **5b Electronic energy.** ωB97M-V/def2-TZVPD, MACE's training level, so E_DFT − E_MACE is a
  per-species MACE error ("MACE trust"). ωB97X-3c or r²SCAN-3c as a cheaper recipe of its own.
  - **Spin check:** the triplet and the next multiplicity for every Ni species.
  - **Oxidation-state check:** spin density on Ni. If it disagrees with L0, refuse with
    `oxidation_state_mismatch` — THQ is redox non-innocent.
- **5c Solvation.** SMD (water) single point on the same geometry. Its difference from ALPB is the
  continuum error bar at this tier.
- **5d Thermochemistry.** An autograd MACE Hessian → quasi-RRHO G at 298 K, validated against DFT
  or xTB Hessians on the 5a sample. Imaginary modes send the geometry back as `not_a_minimum`.
  DFT Hessians are the costliest term and the most tolerant of a cheaper model — MACE's biggest
  accelerator win.
- **Backend: GPU4PySCF** (called 2026-09-29; in-process, on the RTX 4000 Ada). `energy/dft.py`'s
  `DFTBackend` serves `Fidelity.DFT` with `single_point` and `relax`; no `hessian` (5d uses MACE).
  What was verified on install is in §3a.

### Stage 6 — Calibrate, then condition
- **Terms named beside ΔG (C21):**
  - 1 atm → 1 M: +0.08 eV per solute; cancels when the solute count does not change;
  - water activity (monomer or cluster reference, §2);
  - **pH as an input**: −0.059 × pH eV per released H⁺ relative to 1 M H₃O⁺. That 1 M state is
    what "pH 0" means today.
- **Calibration set, same recipe.**
  - pKa: acetic acid 4.76, phenol, catechol, H₂O/H₃O⁺, and tHQ once sourced.
  - Ni(II) stepwise log β: acetate, Cl⁻, and the first hydrolysis of Ni(H₂O)₆²⁺ (NIST SRD 46).

  Fit ΔG_calc against 2.303 RT × experiment. The slope, intercept and scatter are the reported
  uncertainty.
- **Conditioning, the goal.** At stated pH, concentration and water activity, the route legs
  become equilibrium constants. "Which salt forms the complex" becomes a **speciation** — the
  fraction of Ni in each form — solved by mass balance over the escalated species (Newton on
  log-concentrations).

### Stage 7 — Report
Every number shown carries:
- its recipe — `screen` = MACE+ALPB; `refine` = DFT//MACE + SMD + qRRHO(MACE) + calibration;
- the stage-6 conditions;
- an error bar — the continuum spread at `screen`, the calibration scatter at `refine`.

Refine totals appear where the whole equation reached stage 5; screening totals appear greyed
elsewhere.

---

## 3a. Models and backends — called 2026-09-29; E-MH run and E4's backend built 2026-10-02

Three calls by the user. None has a body yet; each takes a D-number when the slice that
implements it lands with a test.

### Relaxation moves to MACE-MH — study E-MH first

MACE-OMOL-0 relaxes today; the call is to relax with **MACE-MH** and measure how well it
reconciles with OMOL-0 before anything depends on the switch. What is known without running
anything: `mace-torch` 0.3.16 in `ebu` lists `mh-0` and `mh-1` among `mace_mp`'s models, and
unlike `mace_omol` (which pins `head="omol"`) that call does not choose a head. The code has no
MH backend: `config.ML_BACKENDS` is MP-0 and OMOL-0 only, and the default is MP-0.

**E-MH — run 2026-10-02 (simple form), results below the plan.** The plan as written:

0. **Zero-compute questions first**, answered from the model card and the loaded model's
   metadata: which heads `mh-1` carries and at what level of theory each was trained; whether the
   chosen head takes **total charge and spin** (OMOL-0 does; a charge-blind head cannot relax
   Ni(tHQ⁻)Cl(H₂O)ₙ⁺ honestly, and would set `charge_aware = False` on its backend exactly as
   MP-0 does). A charge-blind head ends the study here: it is a relaxer for neutral species only.
1. **Code, committed rather than throwaway:** a `MACEMHBackend` (the head is a constructor
   argument and part of the `MethodSpec`, so two heads never share a `methods` row), a config
   alias `mace-mh-1`, and a comparison script over a registry **copy** — both models' relaxations
   stored as ordinary geometries under their own method rows, so the comparison is a registry
   query, not a side file. E1 (connectivity after relaxation) lands first, because metric (a)
   needs it.
2. **Sample, ~40 species × up to 3 conformers**, stratified from `salt_study_k1.db` (read
   through a copy): neutral and charged Ni(II) complexes across the co-ligand window, both salts'
   routes (#132←#94←#28→#68 and #328…#301), free ligands and anions (tHQ, tHQ⁻, Cl⁻, OAc⁻),
   H₂O / H₃O⁺. Both models start from the **same** raw construct; a second pass starts MH from the
   OMOL-0 minimum, which separates "different basin" from "different start".
3. **Metrics, each with a proposed threshold** (proposals, to be argued before the run):
   - (a) connectivity changes after relaxation, per model — MH adds none that OMOL-0 does not;
   - (b) heavy-atom RMSD and Ni–donor distances between the two minima — median < 0.1 Å;
   - (c) OMOL-0 single point on the MH minimum minus the OMOL-0 minimum — median < 0.05 eV,
     the same test 5a applies against DFT;
   - (d) within-identity conformer ranking (Kendall τ, top-1 agreement) — ≥ 80 % top-1;
   - (e) the six equations of the ion/proton check re-priced on each model's geometries — no
     sign or ordering change;
   - (f) failures (non-convergence, charge refusals), wall time and GPU memory per relaxation.
4. **Cost:** ~40 × 3 × 2 relaxations at ~5 s ≈ 20 min of GPU, plus the cross-start pass.
5. **Outcomes to choose between:** MH relaxes and OMOL-0 stays the stage-2/3 energy evaluator
   (a single point on the MH geometry — recorded as two methods, like DFT//MACE); MH does both;
   or OMOL-0 stays for charged species. Stage 5a then re-tests whichever relaxer won against DFT.

**E-MH results.** `MACEMHBackend` (head in the `MethodSpec`, alias `mace-mh-1`) and
`scripts/compare_relaxers.py` are committed. The study read `salt_study_k1.db` read-only rather than
storing into a copy, and skipped metric (e) (route re-pricing). Sample: 39 structures, 93 starts
(9 free species, 43 neutral and 41 charged Ni(II) starts, the demo-route nodes included), both
models from the same start, LBFGS at fmax 0.05 eV/Å (0.20 run also kept).

- **Step 0 ends the plan's main branch: every MH-1 head is charge- and spin-blind.** The `omol`
  head returns the same energy for water at charge 0, ±1 and as a triplet; OMOL-0 moves by tens
  of meV. `MACEMHBackend` carries `charge_aware = spin_aware = False`, so `energy.reference`
  refuses its energies on a charged equation, exactly as MP-0's.

| Metric (fmax 0.05, `omol` head) | free | neutral Ni | charged Ni | proposed |
|---|---|---|---|---|
| (b) median heavy-atom RMSD, MH vs OMOL-0 min | 0.004 Å | 0.39 Å | 0.47 Å | < 0.1 Å |
| (b) median largest M–donor change | — | 0.10 Å | 0.14 Å | |
| cross-start: MH from the OMOL-0 min, median RMSD | 0.002 Å | 0.08 Å | 0.17 Å | |
| (c) OMOL-0 at MH min − OMOL-0 min, median | 0.002 eV | 0.19 eV | 0.54 eV | < 0.05 eV |
| (c) starts over 0.05 eV | 3/9 | 29/43 | 36/41 | |
| (d) top-1 conformer agreement | | 13 / 29 structures (45 %) | | ≥ 80 % |
| (a) relaxations changing connectivity, OMOL-0 / MH | 0 / 0 | 30 / 29 | 32 / 33 | MH adds none |
| (f) median wall time, OMOL-0 / MH (s, shared GPU) | 0.8 / 0.4 | 30 / 14 | 27 / 12 | |

The `spice_wB97M` head gives the same picture (median RMSD 0.40 Å, median (c) 0.39 eV, top-1
11/29). Two findings outside the comparison:

- **Metric (a) indicts today's relaxer as much as MH.** 62 of 84 Ni relaxations (74 %) change
  connectivity under OMOL-0 at fmax 0.05 (52, 62 %, at 0.20), and MH is no different: a carboxylate O moving onto Ni,
  a water leaving to 3–6 Å, a proton crossing an H-bond. E1 is the gate this needs.
- Disagreement with OMOL-0 is not evidence against MH. Only DFT can say which minimum is
  better — `scripts/dft_arbitrate.py`.

**DFT arbitration** (`scripts/dft_arbitrate.py`, ωB97M-V/def2-SVP single point + gradient on both
minima, same starts, fmax 0.05). ΔE = E_DFT(MH min) − E_DFT(OMOL-0 min):

| Structure | class | RMSD | ΔE_DFT | DFT max force OMOL-0 / MH (eV/Å) |
|---|---|---|---|---|
| #132 | neutral | 0.61 Å | **−0.87 eV** | 3.78 / 1.87 |
| #68 | neutral | 0.05 Å | −0.004 eV | 0.25 / 0.34 |
| #569 | neutral | 0.39 Å | +0.21 eV | 0.73 / 1.12 |
| #28 | charged | 0.25 Å | +0.26 eV | 0.38 / 1.35 |
| #94 | charged | 0.33 Å | +0.30 eV | 0.68 / 2.14 |
| #718 | charged | 0.19 Å | +0.69 eV | 1.07 / 3.21 |

DFT prefers the OMOL-0 minimum in 4 of 6 (all 3 charged, by 0.26–0.69 eV), ties on one, and
prefers MH only on #132, where *both* minima sit far from a DFT stationary point (3.8 / 1.9 eV/Å)
— the structure whose relaxation moves a proton; unexplained. Median DFT max force: OMOL-0 0.70,
MH 1.61 eV/Å. Caveats: six species, def2-SVP rather than the training basis, single points
rather than DFT re-optimisation.

**Outcome (2026-10-02): OMOL-0 stays the stage-2 relaxer; MH-1 is selectable, not default.**
The rule agreed beforehand was to flip the default only if DFT backed MH, and it does not. Stage
5a (§3, at def2-TZVPD with DFT re-optimisation) is where #132 and the relaxer get re-tested.

### MACE-POLAR-1 — the same study, three sizes (2026-10-05)

`MACEPolarBackend` (`mace-polar-1-s|m|l`; size in the `MethodSpec`) needs `graph-longrange==0.4.0`
beside mace-torch 0.3.16 — 0.4.3/0.4.4 dropped the `force_pbc_evaluator` argument mace-torch
passes. **Step 0 passes:** every size is given charge and spin and its energy responds — water
loses an electron at +12.65 eV (vertical IP ≈ 12.6 eV), acetate at +3.24 eV.

**OMOL-0's spin response is weak** (`MACEOmolBackend`, same probes): H₂O→H₂O⁺ +0.04 eV, H₂O
triplet −0.04 eV, acetate triplet +0.61 eV, where POLAR-1-L gives +12.65, +6.57, +4.40 eV.
Acetate's electron detachment is right in both (+3.13 / +3.24 eV). Every Ni(II) species here is a
triplet, so this belongs beside any OMOL-0 number that crosses a spin or charge state.

Same 93 starts and settings as E-MH (`compare_relaxers.py --candidate mace_polar_s mace_polar_m
mace_polar_l`):

| vs OMOL-0 | MH-1 | POLAR-S | POLAR-M | POLAR-L |
|---|---|---|---|---|
| neutral Ni: median RMSD / OMOL-0 penalty at the minimum | 0.39 Å / 0.19 eV | 0.27 / 0.12 | 0.15 / 0.035 | **0.11 / 0.004** |
| charged Ni: median RMSD / penalty | 0.47 / 0.54 | 0.32 / 0.30 | 0.24 / 0.083 | **0.21 / 0.028** |
| charged Ni: cross-start RMSD | 0.17 Å | 0.12 | 0.030 | **0.021** |
| top-1 conformer agreement | 13/29 | 11/29 | 13/29 | 17/29 |
| median s per Ni relaxation (OMOL-0 10.5–13.5) | 12–14 | 6–7 | 10–17 | 26–33 |
| relaxations changing connectivity (OMOL-0: 62) | 62 | 58 | 61 | 60 |

Agreement improves monotonically with size; started from OMOL-0's minimum POLAR-L barely moves
(0.015–0.02 Å), so the two share minima and their raw-start differences are mostly basin choice.
**POLAR-S is not a proxy for POLAR-L** — it disagrees about as much as MH-1.

**DFT check, POLAR-L vs OMOL-0** (same six species and protocol as MH-1):

| | #132 | #68 | #569 | #28 | #94 | #718 |
|---|---|---|---|---|---|---|
| ΔE_DFT (POLAR-L − OMOL-0) | −2.69 eV ⚠ | +0.004 | −0.006 | +0.002 | +0.042 | +0.028 |
| DFT max force OMOL-0 / POLAR-L (eV/Å) | 3.78 / 0.59 | 0.25 / 0.23 | 0.64 / 0.66 | 0.38 / 0.32 | 0.68 / 0.68 | 1.07 / 1.37 |

Five of six are DFT-equivalent (|ΔE| ≤ 0.04 eV). ⚠ **#132 — rechecked with the SCF state
recorded:** at OMOL-0's minimum the SCF is unstable — ⟨S²⟩ 2.58 against 2.00 for a triplet, 68
cycles, and three runs at that geometry gave DFT energies 0.56 eV apart (−76242.51 / −76242.93 /
−76243.07 eV). At POLAR-L's minimum it is a clean triplet (⟨S²⟩ 2.003, Ni spin 1.82, 26 cycles)
and reproducible to 1 meV. POLAR-L's minimum is **2.1–2.7 eV lower** than OMOL-0's by DFT, and
the MH minimum (−0.87 eV against its run's OMOL-0 energy) is ~1.4 eV above POLAR-L's. #569 used
different starts in the MH and POLAR checks, so those two rows are each internally fair but not
comparable with each other.

**Where this leaves the relaxer (2026-10-05, not yet called):** POLAR-L is DFT-equivalent to
OMOL-0 on five of six species and much better on the sixth, is given charge and spin with a
physical response, and costs ~2.5× OMOL-0 per relaxation. Run as an S → L funnel (k = 3) its
all-conformer cost is ~0.3 × all-L ≈ 0.8 × all-OMOL-0. Whether POLAR-L (direct or funnelled)
replaces OMOL-0 is the user's call; stage 5a at def2-TZVPD with DFT re-optimisation is the test
that should precede it.

**Small → large on one structure (prototype, 12 Ni starts).** S → L saves 10 % of L-alone time,
S → M → L 25 %; a warm start cuts the L stage's steps only 17 %, because S's minimum is not near
L's. The start also changes the answer: S → L ends in L-alone's minimum in 6/12 starts, S → M →
L in 3/12, in either direction (−0.88 to +1.06 eV). Per call S is 4–6× cheaper than L and M 2×,
at 18–41 atoms. Storage is not a constraint: a geometry is a few kB.

**Small → large across conformers — the funnel** (`scripts/funnel_eval.py`; 12 Ni structures, 6
neutral and 6 charged, 10 raw starts each, fmax 0.05). Every start relaxed by POLAR-L (baseline);
every start by S or M, then that model's k lowest minima refined by L. Regret = best funnel L
energy − best baseline L energy, one theory throughout:

| funnel | k | median regret | worst | missed > 0.05 eV | found lower | time vs all-L |
|---|---|---|---|---|---|---|
| S → L | 1 | −0.08 eV | +0.57 | 4/12 | 6/12 | 0.22 |
| S → L | 2 | −0.13 | +0.34 | 3/12 | 7/12 | 0.26 |
| **S → L** | **3** | **−0.15** | **+0.11** | **2/12** | **7/12** | **0.31** |
| M → L | 1–3 | +0.01 to −0.04 | +0.23 | 3–4/12 | 5–6/12 | 0.46–0.48 |

- **S → L, k = 3 costs 31 % of all-L and is not worse in the median.** It misses the baseline's
  best by at most 0.11 eV (2/12) and finds a lower L minimum in 7/12 — by 1.0–1.1 eV on two charged
  structures (#678, #804).
- **It works as sampling, not as ranking.** S puts L's best conformer first in 2/12; the gain is that
  L refined from S minima reaches basins L from raw constructs does not. The same fact says ten raw
  starts do not converge the conformer search for these complexes under any model.
- **M is not worth a stage:** half of all-L's time, no better ranking than S.
- Timings for #596, #614 and #804 were partly measured while two jobs shared the GPU.

### A smarter funnel — partial relaxation (planned, not built)

The funnel above relaxes every conformer **to convergence** with S, then refines k of them with
L. What it measured decides what "smarter" can mean:

- S's *ranking* is weak (L's best first in 2/12), so pruning harder on S energies alone loses
  answers; the gain was **sampling** — L reaching basins from S minima it does not reach from raw
  constructs. Ten raw starts do not converge the search for any model.
- A warm start cuts L's steps only 17 %: the L stage is paid in L steps, whatever S did.
- Per call S is 4–6× cheaper than L, and its cost barely grows with system size (~50–70 ms at
  18–41 atoms), so S is overhead-bound — conformers of one identity share atoms and could share a
  call.
- 62 of 84 Ni relaxations change connectivity under every model; finishing those relaxations
  buys a geometry filed under the wrong identity.

So the levers are: stop a trajectory as soon as its outcome is known (pruned, duplicate, or
torn apart), spend the saving on *more starts*, batch S, and converge L only where the result is
used. Partial relaxation is the mechanism for all of them. Proposed slices, each measured
before the next is built:

| # | Slice | What it settles | Gate to proceed |
|---|---|---|---|
| F0 | **Trajectory capture.** `funnel_eval.py` records per-step energy, fmax and coordinates (every n steps) for S, M and L, on the 12 structures and on 20 starts each. One run, kept as data | every policy below is then evaluated **offline** by replaying trajectories — no GPU per policy | — |
| F1 | **Policy simulator.** Replays F0 under a policy (S step budgets, prune fraction or energy window, dedup θ, k, L stop criterion) and reports regret and cost against all-L, and against the union of everything found | which schedule; whether pruning on a *partial* S energy keeps L's best | a policy at ≤ 0.3× all-L with regret no worse than k = 3 full-S |
| F2 | **Early dedup.** At each checkpoint, cluster live trajectories with `identity.conformers.cluster` (θ_geom 0.15 Å, calibrated) and keep the lower-energy member. The energy window stays off until [B9](BUGS.md#b9) is fixed — it is uncalibrated | how many of N starts are duplicates after a few dozen S steps | fraction merged, and that no merge removes L's eventual best |
| F3 | **Early stop on connectivity change.** Needs E1's `connectivity_check`; a trajectory whose graph changes is stopped and recorded `connectivity_changed` rather than finished | steps saved, and that the stop never fires on a relaxation that would have returned to the stored graph | E1 landed; false-stop rate on F0 data |
| F4 | **Batched S.** All live conformers of one identity in one MACE call (same atoms, so one graph batch) | S wall time per identity versus N separate calls | measured speedup on the 18–41-atom Ni set |
| F5 | **Staged L.** Refine the survivors to a loose fmax (0.15 eV/Å), re-rank on L, converge only the top 1–2 to 0.05; optionally stop refining when the next candidate's partial L energy is beyond a window of the best | L steps saved without changing the reported minimum | regret ≤ 0.02 eV against full-L refinement on F0 data |
| F6 | **Runner mode `funnel_go`.** Spec fields (small and large model, schedule, k), a spec-version bump so an old spec replays its run (invariant 10), one `relax` task per *identity* instead of per geometry | the pipeline uses it | F1–F5 settled; the user's call on POLAR-L as the relaxer |

**Storage — a decision for F6, proposed here:** store every **converged** S minimum (cheap, a
`methods` row of its own, and discarding paid compute is the expensive mistake) and every L
result; store **no** truncated trajectory point as a geometry — it is at no rung's convergence and
would read as a minimum. The pruning decisions and their reasons go in the task detail, so a
pruned start can be re-run.

**Not in scope:** changing S's ranking (fine-tuning) and the energy window (B9). Stage 5a — DFT
re-optimisation at def2-TZVPD on funnel winners — still comes before POLAR-L becomes the default.

### DFT backend is GPU4PySCF — built (E4, first half)

`energy/dft.DFTBackend`: ωB97M-V/def2-TZVPD by default, RKS for a singlet and UKS otherwise,
density fitting, continua `smd:<solvent>` / `pcm:<solvent>`, geomeTRIC for `relax`, Mulliken
charges and spin populations and ⟨S²⟩ in `EnergyResult.extras`, and `check_oxidation_state`
(`oxidation_state_mismatch`). Functional, basis, grid levels and guess are all in the
`MethodSpec`. Measured in `ebu` on 2026-10-02:

| Check | Result |
|---|---|
| Install | torch 2.13 pins the CUDA **13.0** toolkit (`nvidia/cu13`). `cupy-cuda12x` loads that NVRTC and those headers against its own 12.9 runtime, and `cupy-cuda13x` 14.2 bundles a 13.2 runtime; both fail to compile CuPy's CUB reductions (`cuda_fp8.hpp`), which PCM/SMD hit — and both "work" only when torch is imported first, because torch's 13.0 `libcudart` then wins. **`cupy-cuda13x==14.0.1` + `gpu4pyscf-cuda13x==1.8.1`** has a 13.0 runtime and passes everything without torch. `DFTBackend.available()` compiles a reduction, so the broken pairing reports as unavailable |
| ωB97M-V incl. VV10 on GPU | yes |
| UKS triplet | yes (O₂: ⟨S²⟩ 2.006) |
| SMD and IEF-PCM on GPU | yes, gas and solvated single points |
| Gradients / geomeTRIC | yes |
| **VV10 cost** | the NLC grid defaulted to the full XC grid (306k points on a 23-atom Ni complex) and VV10 is quadratic in it: 20.5 s of a 22 s Fock build. NLC level 1 (104k points) costs 4.5 s and moves the energy by 1.1 µEh; level 0 by 0.17 mEh. Default `nlc_grids_level=1` |
| **Ni(II) SCF** | from the minao guess the triplet oscillates (ΔE of hundreds of Eh in early cycles) and does not converge in 150 cycles; ωB97X also failed in 60. Started from a converged PBE density it converges in 31 cycles. Default `guess_xc="pbe"` |
| Memory, ~40-atom Ni at def2-TZVPD | not yet measured — the first attempt ran before the two fixes above |
| GPU sharing | DFT and a MACE study on one card slowed both several-fold; DFT should run exclusively (an E5 `refine` queue concern) |

### xTB — what it does today, and whether it stays (to be decided)

| Role today | Where | Replaceable by |
|---|---|---|
| **Screening continuum** — ΔG_solv(ALPB) on each MACE geometry, C17's composite, and the `alpb:water` corrections `finalise_run` writes after every run | `energy/solvation.correct`, `runner.finalise_run` | SMD/PCM in GPU4PySCF, at DFT cost per geometry — orders of magnitude more at 10³ geometries |
| The screening **error bar** — ALPB–GBSA spread (E2) | planned | a DFT continuum spread, at the same cost |
| `xtb_go` relaxation run mode | `spec.RUN_MODES`, `relax.MODE_FIDELITY` | MACE; §3 5a already rejects xTB geometries for Ni (−1.9 eV relaxation of a MACE minimum) |
| A charge-aware rung for deprotonation ease | `descriptors/ease.py` | MACE-OMOL-0, already |
| M7's calculator-identity check against the archive | `scripts/regress_m7.py` | nothing — the archive *is* GFN2-xTB |
| Hessian cross-check for 5d | §3 5d | GPU4PySCF on the 5a sample |

**What the decision turns on:** the solvent-reference rule in `WORKPLAN_solvation.md` §3a
(item 4) adds the continuum "on the same rung" — and at the screening rung the only continuum
the stack has is xTB's. Dropping xTB leaves screening with no continuum at all, so every
screening equation falls back to cluster energies and C17's composite goes. Retiring the
relaxation mode costs nothing, but the run mode cannot be deleted: an old spec replays the run it
planned (architecture invariant 10). **Recommendation:** keep xTB narrowly, as the screening
continuum provider and the archive calculator; stop offering `xtb_go` for new specs; revisit when
g-xTB or a GPU continuum at screening cost is benchmarked (`reports/PROJECT_PLAN.md` E-1).

---

## 4. Storage and invariants

**Already enforced:** fidelity per geometry; one theory per equation (`_energy_row` pins the
first term); corrections never borrowed across geometries (C17); absent ≠ zero; named caveats.

**Added:**
- **Recipe.** A named, versioned composite (electronic, solvation, thermal, reference choices).
  An equation is "at recipe R" when every term has all of R's components on one geometry.
  `energy/reference` selects by recipe instead of fidelity rank alone, and `/graph`'s medium
  select becomes a recipe select.
- **Tables.** `thermo_corrections` (geometry, method, zpe, h_corr, g_rrho, n_imag); DFT energies as
  ordinary rows at `Fidelity.DFT`; SMD reuses `solvation_corrections` through C16's
  `model:solvent` token (`smd:water`).
- **Codes.** `connectivity_changed`, `hydration_change`, `oxidation_state_mismatch`,
  `not_a_minimum`.

---

## 5. Slices (each its own PR)

| # | Slice | Size | Unblocks | Milestone |
|---|---|---|---|---|
| E0 | this plan; the Q2 reframing recorded with the salt-study results | S | — | — |
| E1 | connectivity check after relaxation + `connectivity_changed` | S–M | every reported number; E-MH metric (a) | M7 |
| E-MH | MACE-MH and MACE-POLAR-1 backends + the §3a comparisons with MACE-OMOL-0 — **done** (2026-10-05); the relaxer call is open | S–M | the stage-2 relaxer | M7 |
| F0–F6 | the partial-relaxation funnel (§3a, "A smarter funnel") | M | the stage-2 relaxer at screening scale | M7 |
| E2 | screening tier labelled; ALPB+GBSA spread; `hydration_change`; pH and water reference beside the proton sink (C21); the solvation §3a pricing rules | S | honest screening | M7 |
| E3 | MACE Hessian → `thermo_corrections`; recipe selection in `energy/reference` | M | 5d | M7R |
| E4 | GPU4PySCF `DFTBackend` + spin and oxidation-state checks — **backend and checks built**; the `smd:water` producer and the spin ladder remain | M | 5b, 5c | M7R |
| E5 | `escalation_set` + `refine` task kind | M | stage 4 | M7R |
| E6 | calibration specs + fit script | S | stage 6 | M7R |
| E7 | speciation solver; `/graph` shows refine totals | M | the goal | M8 |
| — | H-bonded shells for ions and H₃O⁺ (solvation S4–S6) | M+ | after E6 shows whether the relative scheme suffices | after M7R |

```
E0 ─► E1 ─► E-MH ─► E2 ─► E3 ─┐
       │            E4 ───────┼─► E5 ─► E6 ─► E7
       └─► F3   F0 ─► F1 ─► F2, F4, F5 ─► F6
```

E-MH sits before E2 because E2's error bar and pricing are computed on the relaxer's
geometries; E4 needs only the GPU4PySCF install and can start in parallel. F0–F2, F4 and F5
need no other slice; F3 needs E1's connectivity check, and F6 the relaxer call.

Where the slices touch: `runner._execute_relax` (E1); `energy/backends.backend_for`,
`energy/relax.MODE_FIDELITY` (E4); `energy/reference` `_energy_row`, `_corrected_row`,
`check_reference_quality` (E2, E3); `energy/solvation.correct` (E2); `pathways/route._apply_proton_sink`
(E2).

## 6. Risks

| Risk | Signal | Response |
|---|---|---|
| A screening number is quoted as a result | a MACE+ALPB ΔE beside a pKa or log β | recipe label on every number; only `refine` is a result |
| Triage drops the right answer | a DFT-refined route beats every screened-in route | measure selection recall on a small system refined exhaustively |
| MACE geometry is not a DFT minimum | 5a ΔE > threshold, imaginary modes | re-optimise at DFT for that chemistry; record it on the recipe |
| The L0 oxidation state is fiction | spin density disagrees | `oxidation_state_mismatch`, refused, never priced |
| Hydration comparisons ride on missing entropy | a `hydration_change` step decides a route | caveat until 5d exists |
| The MACE-MH head is charge-blind | E-MH step 0 | MH relaxes neutral species only; OMOL-0 keeps the charged ones |
| One GPU, three tenants | a DFT job out of memory while MACE or the local LLM server holds the card | DFT runs exclusively; the queue declares it (architecture invariant 3: declared, never detected) |
| Dropping xTB removes the screening continuum | xTB retired before a replacement is benchmarked | §3a: keep it as the continuum provider until one is |
