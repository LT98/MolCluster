# Solvation work plan — cluster-continuum, one medium per equation

**Decided by the user, not reopened here.** Solvent effects on reaction energies are modelled as
**cluster-continuum**: (c) an implicit continuum on every energy, plus (b) a small number of
**explicit** first-shell solvent molecules where they matter — held as a **count attached to an
existing species**, with their arrangements treated as conformers (L3), never as new
identity-level species. The explicit count is limited to **±2 around a reference count**, for
water or any solvent a spec names. The goal: coordination of a free ligand, of a leaving Cl⁻ and
of the proton carrier H₃O⁺ by solvent becomes modellable as part of a pathway.

Delete this file when the slices below land; what each decision settles goes to the ledger.
**Ledger number: D-TBD.** D26 is taken (feat/variable-coordination, spec v8) and other branches
are adding entries, so this plan claims checkpoint numbers only (C16–C21, next free after C15)
and leaves the D-number to whoever merges it.

---

## 1. Measured, before designing

All on `data/mvp_ni_thq_cl.db`, read through a backup copy in `/tmp` (the live file was never
opened for writing), on its stored MACE-OMOL-0 best geometries. `ebu`, tblite 0.7.0,
mace-torch 0.3.16, 4 CPU threads for xTB, RTX 4000 Ada for MACE.

### What runs

| Model | Through | Runs in `ebu` today |
|---|---|---|
| ALPB (named solvent, parametrised CDS + shift) | tblite `alpb-solvation` | **yes** — already what `XTBBackend(solvent=…)` uses |
| GBSA (named solvent) | tblite `gbsa-solvation` | **yes** |
| GBε / GB (dielectric only, non-empirical) | tblite `gbe-` / `gb-solvation` | **yes** |
| solution-state option `gsolv` / `bar1mol` / `reference` | ALPB, GBSA 3rd argument | **yes** — H₂O differs by 0.082 eV (`bar1mol`) and 0.186 eV (`reference`) from `gsolv` |
| ddCOSMO / ddPCM / CPCM | tblite `ddX-solvation` | **no** — "ddX solvation model support is not available in this build of tblite" |
| any continuum inside MACE-OMOL-0 | — | **no** — the model has none; `MACEBackend` already refuses `solvent=` |

### Cost

GFN2-xTB single point on a stored Ni complex, median of 5:

| structure | atoms | gas | ALPB water | GBSA water |
|---|---|---|---|---|
| 17 `Ni(tHQ⁻)₂(H₂O)₂` q0 s3 | 37 | 0.08 s | 0.07 s | 0.08 s |
| 15 `Ni(tHQ)₂(H₂O)₂` q+2 s3 | 39 | 0.08 s | 0.06 s | 0.08 s |

Free species are ≤ 0.01 s. A continuum on every stored geometry of this registry (637 with
energies) is two single points each — well under two minutes. Cost is not the constraint.

### What a continuum does to the two equations that prompted this

ΔG_solv = E(xTB, continuum) − E(xTB, gas), **same coordinates** (the stored MACE geometry).

| equation | MACE gas | xTB gas | xTB ALPB | ΔΔG_solv ALPB | ΔΔG_solv GBSA | **MACE + ALPB** | MACE + GBSA |
|---|---|---|---|---|---|---|---|
| R1 `tHQ + H₂O → tHQ⁻ + H₃O⁺` (reaction 633) | +7.356 | +6.501 | +0.747 | −5.754 | −5.400 | **+1.602** | +1.956 |
| R2 `NiCl₂(H₂O)₂ + 2 tHQ + 2 H₂O → Ni(tHQ⁻)₂(H₂O)₂ + 2 Cl⁻ + 2 H₃O⁺` | +17.221 | +14.155 | +0.280 | −13.874 | −12.917 | **+3.347** | +4.304 |

Per species, ΔG_solv ALPB / GBSA (eV): tHQ −0.68/−0.51 · H₂O −0.39/−0.26 · tHQ⁻ −2.42/−2.37 ·
H₃O⁺ −4.41/−3.81 · Cl⁻ −3.10/−2.98 · NiCl₂(H₂O)₂ −1.09/−0.89 · Ni(tHQ⁻)₂(H₂O)₂ −2.11/−1.79.

Geometry sensitivity — ΔG_solv at the MACE geometry vs at each phase's own xTB minimum:
tHQ −0.682/−0.687 · tHQ⁻ −2.422/−2.443 · H₃O⁺ and H₂O identical · **Ni complex 17: −2.110 vs
−1.874** (and xTB relaxes that MACE geometry by −1.9 eV in gas: the two theories do not share a
minimum for the Ni complex).

### Reading

1. The gas-phase numbers are dominated by charge separation in vacuum: the continuum removes
   5.8 eV of R1's 7.4 and 13.9 eV of R2's 17.2. `check_reference_quality` passed both as
   isodesmic with no caveat. That is the defect §3 C18 closes.
2. The continuum does **not** make them right. R1 at +1.6 eV corresponds to a pKa near 27; the
   ligand is a moderately strong acid (a pKa₁ near 5 would mean ~+0.4 eV, *to be sourced before
   anyone quotes it*). The residual is the continuum's known weakness on small ions — H₃O⁺
   above all — plus missing thermal terms. That is what (b), explicit solvation of H₃O⁺ and the
   anions, exists to fix (C20, C21).
3. **Model spread is the error bar**: ALPB vs GBSA on the same geometries is 0.35 eV on R1 and
   0.96 eV on R2. A medium therefore has to name its model, not just its solvent (C16).
4. **ML energy + xTB ΔG_solv on the same geometry (the composite) is sound enough for
   screening, and only on the same geometry.** It is a Δ-correction: the xTB gas-phase error
   largely cancels inside ΔG_solv, which is why MACE+ALPB and xTB-ALPB disagree by 0.9 eV on R1
   while their gas-phase numbers disagree by the same 0.9 eV. It is not a free energy of
   solvation of the ML Hamiltonian, and on the Ni complex the correction moved 0.24 eV between
   two geometries — so a correction is stored against one geometry and never borrowed from a
   sibling geometry of the same structure (C17).

---

## 2. What is built, what is missing

| Piece | State |
|---|---|
| `MethodSpec.solvent` / `methods.solvent`, in the `methods` UNIQUE key | built; **every stored row is NULL** (157 + 8781 + 637 energy rows across the three local registries) |
| `XTBBackend(solvent=…)` → ALPB | built; writes a **bare** solvent name (`'water'`), so the model is implied, not recorded |
| `MACEBackend` refuses `solvent=` | built |
| `energy.reference` selects `m.solvent IS ?` for every term, `same_theory` compares solvent | built — so a direct-in-solvent equation already cannot mix media |
| `routes.price_reaction(solvent=…)` | threads it through; nothing stores a solvated energy for it to find |
| a continuum on an **ML** energy | **missing** — this branch (S1) |
| a reference-quality rule for charge separation | **missing** — this branch (S1) |
| a producer that computes corrections over a registry | missing (S2) |
| explicit solvent count on a species | missing; the identity layer is graph-derived (covalent + dative), and an H-bonded adduct has no edge to hash (S4) |
| a proton reference in a continuum | missing (S6) |

---

## 3. Decisions (C16–C21)

**C16 — where a medium is recorded.** *Options:* (a) `MethodSpec.solvent` only — a solvated
energy is a new `geometries` row at the same coordinates under a solvated method row;
(b) a `geometries.solvent` column; (c) a separate table of corrections keyed on a geometry.
*Recommendation:* **(a) for energies computed in the medium, (c) for a continuum added to an
energy that was not**, and never (b): a geometry is not solvated, an energy is, and a column on
`geometries` would put two media under `_refresh_best_geometry`'s ordering. The medium token is
`'model:solvent'` (`'alpb:water'`); a solution state other than tblite's default goes in the
correction method's `extras`. **`methods.solvent IS NULL` means gas phase, stated, not
unknown** — true of every row today (every producer either writes a solvent or refuses one).
*Sub-decision:* `XTBBackend` still writes bare `'water'`. Move it to `'alpb:water'` under an
`energy_backends` bump; old rows keep `'water'` and simply stop matching new queries (D19 —
never relabel). Not done on this branch.

**C17 — is ML + xTB-continuum correction an acceptable medium?** *Recommendation:* yes, as a
screening-grade medium with these rules, all enforced in S1: the correction is on the **same
geometry** as the energy it corrects; one correction theory for every term (pinned like the
base theory); an equation is met entirely by corrected energies or entirely by energies
computed in the medium, never a mix; the stored `methods` row for the result says so
(`extras.solvation_correction`). Error bar to quote: the ALPB–GBSA spread (§1).

**C18 — charge separation as a reference-quality issue.** New code `charge_separation`: the
stoich-weighted count of free charged species changes across the arrow (either direction —
ion pairing is the same error with the other sign). *Recommendation:* it is **not** a bond-type
rule, so `isodesmic` keeps its meaning (the proton-transfer edges stay isodesmic, as
`protons.py` says); instead `ReferenceQuality.acceptable` is false while it is **blocking**.
Blocking in gas phase — `strict=True` refuses, `store_reaction_energy` needs `force` — and a
**non-blocking caveat** when the equation is weighed in a continuum, because the continuum
screens most of the Coulomb term and not all of it (R1 still +1.6 eV). It never disappears.
The rule reads only the terms it is given, so it is the same rule on a stored reaction, a
derived decomposition, or a route's composed net equation (feat/reversible-hops).
*Consequence to accept:* every gas-phase deprotonation edge now carries a blocking caveat,
and the legacy formation equation gains a third code.

**C19 — how an explicit shell attaches to a species.** *Options:* (a) new `structures` rows
for S·(H₂O)ₙ; (b) extra atoms on `geometries` rows of S; (c) a `solvent_shells` table —
(structure_id, solvent_structure_id, n, coords, method, energy, conformer label) — with a
`Term` gaining `shell=(solvent_id, n)`. *Recommendation:* **(c)**. (a) is ruled out by the
user's decision and by the graph (no H-bond edge to hash). (b) breaks the invariant every
geometry reader relies on — atom count equals the structure's — in a dozen places.
In (c) balance adds `n × solvent` to the term's composition; the dative count is untouched
(H-bonds are not DATIVE edges); arrangements of one (S, solvent, n) are L3 conformers.

**C20 — the reference count, the ±2 window, and which arrangement counts.** *Options for the
reference:* a spec field · a curated `data/reference/solvent_shells.tsv` (species pattern,
solvent, n_ref, source). *Recommendation:* the curated table — it is chemistry, not run data,
and must carry a source per row; the api refuses n outside `[n_ref−2, n_ref+2] ∩ ℕ`. Within one
(S, solvent, n, theory), the **lowest** arrangement, with the number sampled recorded so a
one-arrangement shell reads as unsampled rather than as converged (absent ≠ zero). Boltzmann
weighting waits for thermochemistry.

**C21 — the proton reference and standard state.** *Recommendation:* the proton carrier in a
continuum is H₃O⁺(H₂O)ₙ as a C19 shell on H₃O⁺ (n_ref = 3, the Eigen cation, to be confirmed
against the ALPB residual in §1); water as a reagent is at its liquid activity, a
RT·ln(55.3) ≈ 0.10 eV term per water consumed; other solutes at 1 M. These go into
`ReactionEnergy` as **named terms beside dE**, never folded into stored energies, and wait for
the thermochemistry slice.

---

## 4. Slices

### S0 — call C16–C21 · S
**Done when:** each has a line in the ledger under one D-number (D-TBD) and this file says so.

### S1 — a medium on a stored energy, and the charge-separation rule · S · *this branch*
1. `solvation_corrections` table (geometry, correction method, e_gas, e_solv, dG_solv) and a
   migration log line for a table new to an existing database.
2. `registry.put_solvation_correction` / `solvation_corrections`: refuses a bare solvent name,
   a correction on an energy already in a medium, an unstated or mismatched charge/spin, and a
   repeat that disagrees with the stored value.
3. `energy.reference`: `solvent='model:solvent'` met by direct or by corrected energies —
   one route per equation, loud refusal on a mix; the correction theory is pinned.
4. `charge_separation` per C18; `reference_scheme` → `balanced2`.
**Done when:** tests cover each refusal, the composite arithmetic, and the gas/continuum
behaviour of the new code; `scripts/check.sh` passes.

### S2 — the producer · S
`energy/solvation.py`: `correct(reg, geometry_id, medium)` runs the xTB pair on the stored
coordinates and writes through S1; `scripts/solvate.py --medium alpb:water` over a registry,
resumable. C16's `XTBBackend` token change rides here.
**Done when:** reaction 633 and R2 price in `alpb:water` from the registry and reproduce §1.

### S3 — readers · S
`routes.price_reaction` passes the medium into `check_reference_quality`; the medium becomes a
selectable input of the `/graph` page (ui branch). **Done when:** a route renders with the
`charge_separation (caveat)` badge in a continuum and blocking in gas.

### S4 — explicit shells: storage and balance · M
C19's table, `Term.shell`, balance and quality over shell terms, the ±2 window from C20.
**Done when:** `tHQ⁻·(H₂O)₂` and `H₃O⁺·(H₂O)₃` can be stored and an equation over them balances.

### S5 — explicit shells: producer · M
H-bond-directed placement of n solvent molecules on donor/acceptor sites (the site catalog
already knows the acceptors), arrangement sampling, L3 labels, relaxation in the continuum.

### S6 — proton reference and standard-state terms · S
C21. **Done when:** R1 is re-priced with H₃O⁺(H₂O)₃ and the water-activity term, and the
result is compared with a sourced pKa — the first external check this plan has.

### S7 — the archive · S
Re-derive `legacy/solvation_raw.csv`'s per-species ALPB single points for the species this
registry can build; the archive's scheme stays non-isodesmic (WORKPLAN_M7 §6).

```
S0 ─► S1 ─► S2 ─► S3
       └──► S4 ─► S5 ─► S6        S7 after S2
```

---

## 5. Where this sits

- **B19** (two heights for one structure): the 7.356 eV gap between the routes *is* R1 — one
  gas-phase charge separation. S1 puts `charge_separation` on the step that carries it, and
  in `alpb:water` the gap shrinks to ~1.6 eV but does not vanish: the routes still differ by a
  proton in their basis, so B19's labelling fix is still needed. A common-basis view becomes
  less dangerous in a continuum, not safe.
- **M7**: "the same for the solvation sweep" in the exit gate has had nothing to run against;
  S2 + S7 are that. The medium is one more dimension `same_theory` already guards.
- **M8**: the secondary exit gate (sink detection reproducing the EDTA sequestration finding)
  is a solvation result — it reversed once a medium was added — so S1–S2 are prerequisites
  of it, and S4–S6 of any pathway in which a ligand, Cl⁻ or H₃O⁺ is solvated along the way.

---

## 6. Risks

| Risk | Signal | What to do |
|---|---|---|
| A continuum number is read as a free energy | a dE in `alpb:water` quoted against a pKa | it is an electronic energy + ΔG_solv; C21's terms are separate and named |
| The composite hides a geometry mismatch | xTB relaxes an ML geometry by eV (Ni complex: −1.9 eV) | the correction is per geometry, never borrowed (C17); S2 records e_gas so the mismatch is visible |
| Media mixed inside an equation | a corrected term beside a direct one | refused by name (S1) |
| Model spread mistaken for chemistry | two media compared as if one were right | quote ALPB–GBSA spread as the error bar |
| The caveat becomes wallpaper | every continuum step badged | badge only when the ion count changes; it is one code, not a flag on everything |
| Shell counts explode | n_ref ± 2 over every species in a route | the window is per species and a curated reference; S5 samples, never enumerates all |
| A bare `'water'` row matches a `'water'` query and means ALPB by accident | old xTB rows | C16's sub-decision; the corrected route requires `model:solvent` already |
