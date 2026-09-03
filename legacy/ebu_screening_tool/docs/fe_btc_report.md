# Fe(III) + BTC validation run — and Ni(II) vs Fe(III) comparison

**System:** Fe³⁺ + benzene-1,3,5-tricarboxylic acid (BTC), ± EDTA / EDDA additives.
Fe(III) trimesate is the canonical BTC framework (MIL-100(Fe)), so Fe³⁺ is the relevant
oxidation state. Same `ebu_core` enumeration + GFN2-xTB formation energies as the Ni run.

**Spin treatment:** Fe(III) with carboxylate/O-donors is high-spin d⁵, so every Fe EBU and the
Fe³⁺ reference were computed as a **sextet (5 unpaired e⁻)** — not the minimal-spin convention
used for Ni(II). A fixed, consistent choice; it shifts all Fe formation energies by the same
reference offset and does not affect their relative ranking.

**One enumeration note:** Fe required a slightly looser clash cutoff (`clash_scale` 0.60 vs the
0.65 default used for Ni) to keep the Fe–BTC framework node. At 0.65 that node was rejected
purely because Fe's van der Waals radius is ~0.1 Å larger than Ni's, pushing an innocuous
Fe···H contact (2.11 Å, on the far side of the ring) just past the cutoff — it relaxes away
immediately. This is a QC-threshold artifact, not chemistry.

---

## Headline

- **Most stable species overall:** Fe–EDTA hexadentate octahedral chelate, **−62.8 eV**.
- **Fe–BTC framework node** (target): **−59.0 eV** — much deeper than Ni's −36.4 eV, as expected
  for a +3 vs +2 cation binding anionic carboxylates. Additive-insensitive (−57.7 to −59.0 eV
  across all scenarios).
- Per bond, the BTC carboxylate contacts are far stronger than EDTA's (−29 vs −10 eV/bond);
  EDTA only wins on *total* energy by making six bonds. Same pattern as Ni, more pronounced.

---

## Landscape per scenario (target Fe–BTC node = −59.0 eV)

| Scenario | # EBUs | span (eV) | states **below** node (out-compete) | gap of deepest state below node |
|---|---|---|---|---|
| BTC only | 3 | 8.4 | none | — |
| BTC + EDTA | 8 | 12.2 | 4 (down to −62.8) | **−3.8 eV** |
| BTC + EDDA | 14 | 19.2 | none | — |
| BTC + EDTA + EDDA | 24 | 35.0 | 4 (down to −62.3) | −3.4 eV |

---

## The key Ni vs Fe difference — EDTA behaves very differently

The additive *shape* is the same for both metals (EDDA = smoothing ladder above the node; EDTA =
states below the node). **What changes is how deep the EDTA sink is relative to the framework
node:**

| | Ni(II) | Fe(III) |
|---|---|---|
| Framework node | −36.4 eV | −59.0 eV |
| Deepest EDTA-chelate state | −44.9 eV | −62.8 eV |
| **Gap below the framework node** | **−8.5 eV (deep trap)** | **−3.8 eV (near-parity)** |

- **Ni(II) + EDTA = sequestration.** The Ni–EDTA chelate sits 8.5 eV below the framework node —
  a deep well that locks Ni up before it can build framework. EDTA poisons Ni–BTC crystallization.
- **Fe(III) + EDTA = genuine modulation.** The Fe–EDTA chelate is only ~3.8 eV below the (very
  deep) Fe–BTC node — near-thermodynamic-parity. EDTA competes reversibly for Fe rather than
  irreversibly trapping it. This is the well-behaved "strong modulator" regime, and it matches
  why Fe(III) carboxylate MOFs (MIL-100/101) are routinely and successfully grown *with* strong
  competing acids/modulators. Fe(III) binds BTC so strongly that even EDTA can't create a runaway
  trap.
- **EDDA** smooths for both metals: a graded ladder of competing states at/just above the node,
  none below it — the mild-modulator signature, on both Ni and Fe.

Also note Fe's **bare-BTC** landscape is already a little rougher than Ni's (span 8.4 vs 0.6 eV):
Fe(III) supports a mono-BTC bidentate node (−50.6 eV) alongside the bridged node, so Fe has some
intrinsic competing-motif structure even without additives.

---

## Direct answers

- **Most thermodynamically favourable EBU:** Fe–EDTA hexadentate octahedral (−62.8 eV);
  the framework-relevant target is the bridged Fe–BTC node (−59.0 eV), very deep and
  additive-insensitive.
- **Do additives smooth competing EBUs / aid crystallization?**
  - **Fe(III): EDTA is now viable as a modulator** (near-parity competition, not a trap) — this is
    the main change from Ni. EDDA remains a mild smoothing modulator.
  - Best synthesis translation for a Fe–BTC (MIL-100-type) attempt: try **EDTA as a strong
    modulator** and/or **EDDA as a mild one**; both are usable, EDTA is the more aggressive lever.
    Contrast with Ni, where EDTA should be avoided (or kept far sub-stoichiometric) and EDDA is
    the safe choice.

---

## Caveats

Same as the Ni run: screening-grade GFN2-xTB, gas-phase, no solvent/counterions (real syntheses
compress these gaps), fixed high-spin Fe(III), relative ranking only. Many large multi-ligand
EBUs did not fully relax within the step cap (`converged=False`) — their energies are upper
bounds and they are high-lying states, so the landscape conclusions are unaffected. The Ni run
used a 90→50-step relaxation cap and Fe a 50→35-step cap for throughput; this only affects the
crowded high-lying states, not the low-lying framework/chelate nodes that drive the conclusions.

## Files

- `fe_btc_xtb_results.csv` — all 49 Fe EBUs (energies, charge, convergence).
- `ni_vs_fe_landscape.png` — side-by-side landscape, each normalized to its own framework node.
- `fe_btc_outputs/` — Fe EBU `.xyz` structures.
- `run_metal.py`, `xtb_energy_metal.py` — metal-parametrized drivers (set `METAL_SYMBOL`,
  `METAL_CHARGE`, `METAL_MULT`) — reuse for the next metal.
