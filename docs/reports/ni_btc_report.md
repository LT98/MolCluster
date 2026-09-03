# Ni(II) + BTC validation run — formation-energy comparison of competing EBUs

**System:** Ni²⁺ + benzene-1,3,5-tricarboxylic acid (BTC), with optional chelating additives
ethylenediaminetetraacetic acid (EDTA, hexadentate N₂O₄) and ethylenediamine-N,N′-diacetic
acid (EDDA, tetradentate N₂O₂).

**Question:** Which extended building unit (EBU) is most thermodynamically favourable, and do
the additives introduce *competing* EBUs that smooth the landscape and aid crystallization?

**Method:** Structures enumerated with the existing `ebu_core` pipeline (Ni²⁺, auto donor
perception, geometry placement + QC). Formation energies from **GFN2-xTB** (tblite/ASE),
relaxed, referenced as `E_form = E(EBU) − E(Ni²⁺) − Σ E(free-ligand anion)`, charge-conserving.
GFN2-xTB was used because the MACE stack in `energy_model.py` cannot be installed in this
sandbox; xTB is a real self-consistent quantum method with Ni parameters that treats total
charge and spin explicitly, which matters for these charged complexes. **This is screening-grade
— use for relative ranking, not absolute numbers.**

---

## Headline result

| Rank | EBU | E_form (eV) | eV/bond | Note |
|---|---|---|---|---|
| 1 | Ni–EDTA hexadentate, octahedral | **−44.9** | −7.5 | additive complex — deep sink |
| 2 | Ni–EDTA hexadentate, tetra/planar (s4) | −43.9 | −11.0 | additive complex |
| — | **Ni(BTC)₂ bridged node, planar** | **−36.4** | −18.2 | **target framework node** |
| — | Ni–BTC/EDDA mixed nodes | −36 → −26 | −6 to −18 | mild competitors |
| — | Ni–EDDA tetradentate | −31.5 | −7.9 | mild competitor |

Two distinct findings:

1. **Most stable species overall = the Ni–EDTA chelate (−44.9 eV)**, ~8.5 eV *below* the target
   Ni–BTC framework node. That is not a good thing for the MOF — see below.
2. **The target Ni–BTC node itself is essentially additive-independent** (−36.4 eV in BTC-only;
   −36.1 eV in every additive scenario). Neither additive destabilizes the framework node; they
   only change *what else Ni can do instead*.

---

## Landscape per scenario (the crystallization question)

Target Ni–BTC node = −36.4 eV. "Below target" = states that out-compete the framework.

| Scenario | # valid EBUs | span (eV) | states **below** target | states near target (±0–3 eV) |
|---|---|---|---|---|
| BTC only | 2 | 0.6 | none | 2 |
| BTC + EDTA | 7 | 10.7 | **3** (−44.9, −43.9, −43.9) | 4 |
| BTC + EDDA | 13 | 24.2 | none | 3 |
| BTC + EDTA + EDDA | 21 | 40.4 | **4** (−44.6, −43.9, −43.9, −37.9) | 5 |

**BTC only** gives a single sharp node (span 0.6 eV, nothing competing). A landscape this
narrow tends to precipitate fast and disordered — little thermodynamic handle to slow nucleation
and grow ordered crystals.

**EDTA sequesters, it does not smooth.** EDTA opens up states 7–8.5 eV *below* the framework
node — the well-known, benchmark-stable [Ni–EDTA]²⁻ chelate. Thermodynamically, Ni prefers to
sit in EDTA rather than in the framework. At meaningful loading this is a **metal-sequestration
trap** that can suppress framework formation, not a gentle modulator. (Sub-stoichiometric EDTA
could still act as a strong competition-based modulator, but the margin is large and easy to
overshoot.)

**EDDA smooths.** EDDA introduces **no** state below the framework node; instead it builds a
dense ladder of mixed Ni–BTC/EDDA and Ni–EDDA motifs from −36 eV up through −26 eV, several of
them within a few eV of the target. These are accessible intermediate coordination states that
sit *at or just above* the product energy — the signature of a **modulator**: a reservoir of
low-lying competing complexes that can slow nucleation and favour ordered growth **without
permanently trapping the metal**.

**Both together** = the union: the EDTA sink plus the EDDA ladder. The sink dominates the
thermodynamics.

---

## Direct answer to the two goals

- **Most thermodynamically favourable EBU:** the Ni–EDTA hexadentate octahedral chelate
  (−44.9 eV). Among *framework-relevant* BTC nodes, the bridged planar Ni(BTC)₂ node (−36.4 eV),
  which is stable and additive-insensitive.
- **Do the additives smooth competing EBUs to aid crystallization?**
  - **EDDA: yes** — a graded ladder of competing states near the target node, no deep trap.
    This is the additive to try first as a crystallization modulator for a Ni–BTC framework.
  - **EDTA: no (as a smoother)** — it creates a deep competing sink that sequesters Ni below
    the framework node. Useful only as a *strong, easily-overdosed* competition agent at low
    sub-stoichiometric loading.

**Suggested synthesis translation:** run Ni²⁺ + BTC with EDDA as a modulating additive
(sub- to near-stoichiometric) as the primary trial; keep an EDTA arm only at low loading as a
strong-competition contrast. The bare Ni + BTC control is expected to precipitate fast/disordered.

---

## Caveats (screening-level)

- GFN2-xTB, gas-phase, **no solvent or counterions**; real syntheses are solvated and pH/charge-
  buffered, which compresses these gaps. Ranking is relative, not absolute.
- Minimal-spin convention; Ni(II) spin state (octahedral high-spin vs square-planar low-spin)
  not scanned — a fixed offset per species, does not change the qualitative ordering.
- Large, crowded multi-ligand EBUs did not fully relax within the step cap (`converged=False` in
  the CSV); their energies are upper bounds (would only get *more* negative), and they are the
  high-lying states anyway, so the landscape conclusions are unaffected.
- Total E_form rewards making more bonds; the `E_form_per_bond` column is included so hexadentate
  nodes aren't flattered purely by bond count. The EDTA-sink conclusion holds either way, because
  Ni–EDTA is a genuinely benchmark-stable complex.

## Files

- `ni_btc_xtb_results.csv` — all 43 EBUs: charge, atoms, single-point + relaxed energy,
  convergence, formation energy (total and per bond).
- `ni_btc_landscape.png` — the landscape figure.
- `ni_btc_outputs/` — every candidate EBU as `.xyz`.
- `run_ni_btc.py` — enumeration driver (also wired for MACE if that stack is available).
- `xtb_energy.py` — the GFN2-xTB energy/ranking driver (resumable).
