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
| 2 Relax + triage | geometry, **connectivity check**, conformer ranking and pruning | **MACE-MH** relaxes (decided, not built — §3a); MACE-OMOL-0 until study E-MH | 10³ | ~5 s GPU |
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
- **Backend: GPU4PySCF** (called 2026-09-29; in-process, on the RTX 4000 Ada). A `DFTBackend`
  with `single_point` / `relax` / `hessian` behind the existing contract, filling the
  `Fidelity.DFT` slot that `energy/backends.backend_for` answers with `NotBuiltYet` today. Not
  installed in `ebu` yet (nor PySCF); what must be verified before E4 relies on it is in §3a.

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

## 3a. Models and backends — called 2026-09-29, not built

Three calls by the user. None has a body yet; each takes a D-number when the slice that
implements it lands with a test.

### Relaxation moves to MACE-MH — study E-MH first

MACE-OMOL-0 relaxes today; the call is to relax with **MACE-MH** and measure how well it
reconciles with OMOL-0 before anything depends on the switch. What is known without running
anything: `mace-torch` 0.3.16 in `ebu` lists `mh-0` and `mh-1` among `mace_mp`'s models, and
unlike `mace_omol` (which pins `head="omol"`) that call does not choose a head. The code has no
MH backend: `config.ML_BACKENDS` is MP-0 and OMOL-0 only, and the default is MP-0.

**E-MH — planned, not run.** Nothing below runs until this plan is agreed.

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

### DFT backend is GPU4PySCF

Replaces the open "ORCA or GPU4PySCF". Neither `gpu4pyscf` nor `pyscf` is installed in `ebu`.
To verify on install, before E4 is sized: SMD (5c needs it; otherwise PCM with SMD as a gap) on
GPU; ωB97M-V's non-local (VV10) term on GPU; unrestricted Kohn–Sham for the Ni(II) triplet;
memory for a ~40-atom Ni complex at def2-TZVPD on 20 GB; and **GPU sharing** — MACE relaxations
and the local LLM server use the same card, so a DFT job is scheduled exclusively. DFT Hessians
stay out of scope: 5d's validation uses MACE vs DFT on the 5a sample only.

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
| E-MH | MACE-MH backend + the §3a comparison with MACE-OMOL-0 | S–M | the stage-2 relaxer | M7 |
| E2 | screening tier labelled; ALPB+GBSA spread; `hydration_change`; pH and water reference beside the proton sink (C21); the solvation §3a pricing rules | S | honest screening | M7 |
| E3 | MACE Hessian → `thermo_corrections`; recipe selection in `energy/reference` | M | 5d | M7R |
| E4 | GPU4PySCF `DFTBackend` + spin and oxidation-state checks | M | 5b, 5c | M7R |
| E5 | `escalation_set` + `refine` task kind | M | stage 4 | M7R |
| E6 | calibration specs + fit script | S | stage 6 | M7R |
| E7 | speciation solver; `/graph` shows refine totals | M | the goal | M8 |
| — | H-bonded shells for ions and H₃O⁺ (solvation S4–S6) | M+ | after E6 shows whether the relative scheme suffices | after M7R |

```
E0 ─► E1 ─► E-MH ─► E2 ─► E3 ─┐
                    E4 ───────┼─► E5 ─► E6 ─► E7
```

E-MH sits before E2 because E2's error bar and pricing are computed on the relaxer's
geometries; E4 needs only the GPU4PySCF install and can start in parallel.

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
