# Solvation shell — cluster-continuum model (water / ethanol / DMF)

## How it fits in (the neat scheme)

Solvation enters as a **cluster-continuum** correction that reuses the whole existing pipeline:

1. **Explicit first shell** — the metal's octahedral coordination sphere is completed with real
   solvent molecules (H₂O, EtOH, DMF as neutral O-donors) on whatever sites BTC/chelator leave
   open. Solvent thus competes for *the same sites the framework would grow into*.
2. **Implicit bulk** — the whole cluster is then evaluated in a GFN2-xTB/ALPB dielectric
   continuum for the chosen solvent.
3. **Mixtures (1:1, 1:1:1)** — explicit molecules split as evenly as possible across the chosen
   solvents; the continuum energy is the equal-weight average of the component solvents' ALPB
   single points. This makes every 1:1 / 1:1:1 system fall out of just three single points per
   structure (water, ethanol, DMF), so the 7-system sweep is nearly free.

Charge-conserving formation energy in solvent system *S*:

```
E_form,S(cluster) = E_S(cluster) − E_S(M^q+) − Σ E_S(ligand anion) − Σ E_S(free solvent)
```

all terms in the same continuum, each = mean over *S*'s components. Species per metal (octahedral,
CN 6): **ion** [M(solv)₆] (dissolved baseline), **node** [M(BTC)₂(solv)₄] (solvated framework
node), **edta** [M(EDTA)] (saturated chelate). Geometry: light gas relaxation, then ALPB single
points (standard screening cluster-continuum).

Two decision metrics come out:
- **Nucleation driving force** = E_form(node) − E_form(ion). Negative → BTC beats coordinated
  solvent, metal leaves solution to form framework.
- **Sequestration margin** = E_form(EDTA) − E_form(node). Negative → EDTA out-competes the
  framework (traps the metal).

---

## Headline: solvation flips the EDTA verdict

| | Ni(II) seq. margin | Fe(III) seq. margin |
|---|---|---|
| **Gas phase (earlier run)** | **−8.5 eV — deep trap** | −3.8 eV |
| **With solvation (any of the 7 systems)** | **−0.3 to +0.5 eV — near-parity** | **+1.6 to +2.2 eV — framework wins** |

The gas-phase calculation badly over-stated EDTA's advantage. It compared a compact, highly
charged Ni–EDTA chelate against the node with no medium to screen that charge. Once the solvent
is present — explicitly capping the node's open sites *and* screening the charge through the
continuum — the −4 EDTA chelate loses almost all of its apparent edge:

- **Ni(II):** EDTA and the solvated framework node are within a few tenths of an eV in every
  solvent — EDTA behaves as a **modulator, not a trap**, and the balance is solvent-tunable.
- **Fe(III):** the framework node is favoured by ~2 eV in every solvent — EDTA is a pure
  modulator, never a sink.

This is exactly why solvation had to be included: it changes the qualitative recommendation for
Ni from "avoid EDTA (it sequesters)" to "EDTA is a usable, tunable modulator in solution."

---

## Solvent-specific readout

**Nucleation driving force (top panel):**
- **Water gives the weakest drive for both metals (~−0.5 eV)** — the aqua-ion is so well
  stabilized that BTC barely wins. Slow, reversible, controllable nucleation (good for large
  single crystals, risky for yield).
- **Ethanol and DMF give the strongest drive** (Ni −1.9 to −2.0; Fe −2.1 to −3.4 eV). DMF is the
  classic MOF solvent and shows a strong, clean driving force here — consistent with practice.
- Mixtures interpolate; **W:E:D and E:D** give a strong drive while keeping EDTA at bay.

**Sequestration margin (bottom panel):**
- **Ni:** water, DMF and W:E:D put the node slightly *below* EDTA (framework favoured);
  ethanol and E:D tip marginally toward EDTA. So for Ni, **avoid ethanol-rich mixtures if EDTA is
  present; prefer DMF or aqueous/DMF blends.**
- **Fe:** every solvent keeps the framework ~2 eV below EDTA — solvent choice here is about
  nucleation control, not sequestration risk.

---

## Suggested synthesis translation

- **Ni–BTC + EDTA modulator:** viable in solution (contrary to the gas-phase verdict). Use **DMF
  or a water/DMF blend** for a clean framework-favoured window with a moderate nucleation drive;
  water alone for slow controllable growth; **avoid ethanol-dominant mixtures** with EDTA.
- **Fe–BTC (MIL-100-type):** framework is robustly favoured over EDTA in all solvents. Pick the
  solvent by nucleation control — **water** for slow growth, **DMF/ethanol** for strong drive.

## Caveats

Screening-grade throughout: GFN2-xTB energies, ALPB continuum (a simple dielectric, not a full
solvation free energy; charged-anion solvation from ALPB is approximate and these species carry
large charges), mixtures approximated as equal-weight averages of pure-solvent single points,
one representative explicit shell per composition, light relaxation. The **direction and rough
magnitude** of the effects (solvation collapses the EDTA advantage; water = weak drive; DMF =
clean strong drive) are the trustworthy output — treat individual numbers as ±several tenths of
an eV.

## Files
- `solvation_summary.csv` — per metal × 7 solvent systems: E_form(node/ion/edta), nucleation
  driving force, sequestration margin.
- `solvation_raw.csv` — raw ALPB single-point energies (water/ethanol/DMF) for every cluster and
  reference, for re-analysis at any mixture ratio.
- `solvation_landscape.png` — the two-metric figure.
- `solvation_model.py` — the cluster-continuum driver (resumable; extend species/metals here).
