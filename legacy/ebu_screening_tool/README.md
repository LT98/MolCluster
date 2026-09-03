# EBU Screening Tool

Screening-level pipeline for metal–organic **extended building units (EBUs)**: enumerate the
extended building units a metal centre forms with a linker and optional chelating additives,
rank them by GFN2-xTB formation energy (gas or implicit solvent), model the solvation shell
(explicit + continuum), and view any structure in 3D. Configurable by metal, linker, additive
and solvent.

## Quick start (fresh Linux machine)

```bash
unzip ebu_screening_tool.zip && cd ebu_screening_tool
bash setup.sh                 # makes .venv, installs everything, runs the self-test
source .venv/bin/activate
jupyter lab ebu_tool.ipynb    # <- the tool
```

`setup.sh` ends by running `selftest.py`; if it prints **ALL CORE CHECKS PASSED** you're ready.

Prefer conda? `conda env create -f environment.yml && conda activate ebu` instead of `setup.sh`.

## The notebook (`ebu_tool.ipynb`)

Run top-to-bottom; edit the **CONFIG** cell (section 1):

```python
METALS       = ['Ni', ('Fe', 3)]   # symbol (default charge) or (symbol, charge)
LINKERS      = ['BTC']             # framework linker(s)
ADDITIVES    = ['EDTA', 'EDDA']    # chelating / modulating additives
ADDITIVE_MODE = 'subsets'          # 'subsets' | 'all' | 'none'
MAX_COPIES   = 2                   # max copies of each ligand per metal
CLASH_SCALE  = 0.65
SPIN_MODE    = 'highspin'          # 'highspin' (d-electron based) | 'auto'
```

Sections: **2** enumerate → **3** rank by formation energy (+ plot vs *chelation extent*) →
**4** solvation sweep (water/ethanol/DMF combinations) → **5** interactive 3D viewer.

Add ligands/solvents on the fly: `LIGAND_LIBRARY['NDC'] = 'OC(=O)c1ccc2cc(C(=O)O)ccc2c1'`, or pass
`'Name:SMILES'` directly in LINKERS/ADDITIVES.

## Use as a library

```python
import ebu_toolkit as T
df     = T.enumerate_models(metals=[('Cu',2)], linkers=['BTC'], additives=['EDTA'])
ranked = T.formation_energies(solvent='dmf', spin_mode='highspin')   # or solvent=None (gas)
solv   = T.solvation_sweep([('Fe',3)], linker='BTC', chelator='EDTA')
T.view_model(T.list_models()[0])                                     # 3D
```

## Files

| file | purpose |
|---|---|
| `ebu_tool.ipynb`      | the interactive tool (start here) |
| `ebu_toolkit.py`      | high-level API: `enumerate_models`, `formation_energies`, `solvation_sweep`, `view_model` |
| `ebu_core.py`         | geometry engine: ligand perception, EBU assembly, placement, enumeration |
| `geometry_qc.py`      | post-placement QC (clash / bond / ring-planarity checks) |
| `donor_perception.py` | generic donor-atom detection + deprotonation |
| `energy_model.py`     | optional MACE-MP-0 backend (alternative to xTB; `pip install torch mace-torch`) |
| `ebu_viz.py`          | py3Dmol 3D rendering helpers |
| `run_metal.py`        | batch enumeration driver (env-configurable) |
| `xtb_energy_metal.py` | batch formation-energy ranking driver |
| `solvation_model.py`  | standalone cluster-continuum solvation model |
| `dft_predict.py`      | MACE-MP-0 energies straight from `.xyz` files |
| `docs/`               | method notes + the Ni / Fe / solvation study reports |

## Notes

- **Energies are screening-grade GFN2-xTB** — trust relative rankings and the direction of
  solvent/additive effects, not absolute numbers.
- `tblite` installs from pip on Linux/macOS. (On Windows only, use `conda install -c conda-forge
  tblite-python`.)
- Optional higher-fidelity path: install MACE (`pip install torch mace-torch`) and use
  `energy_model.py` / `dft_predict.py`.
