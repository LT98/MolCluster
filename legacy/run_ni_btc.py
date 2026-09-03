"""
run_ni_btc.py -- validation run: Ni(2+) + 1,3,5-benzenetricarboxylic acid (BTC),
with optional chelating additives EDTA and/or ethylenediamine-N,N'-diacetate
(EDDA).

Goal
----
1. Enumerate candidate EBUs (extended building units) for four ligand sets:
      A) BTC only            -- the baseline target linker
      B) BTC + EDTA          -- add hexadentate N2O4 chelator
      C) BTC + EDDA          -- add tetradentate N2O2 chelator
      D) BTC + EDTA + EDDA    -- both additives present
2. For each surviving (geometry-valid) EBU, estimate an approximate MACE-MP-0
   formation energy so candidates can be RANKED against each other. This tells
   us (i) which EBU is most thermodynamically favourable, and (ii) whether the
   additives introduce lower-lying competing EBUs that could smooth the
   crystallization energy landscape.

The energy step needs the MLIP stack (torch/ase/mace). If it is not installed
in the current environment, the script still does the full geometry
enumeration + QC and writes every candidate to disk; the energy columns are
left blank and a note is printed. Install the stack with:

    pip install torch --index-url https://download.pytorch.org/whl/cpu
    pip install mace-torch ase

and re-run to fill in the formation energies.
"""
import os, sys, csv, warnings, argparse
warnings.filterwarnings('ignore')

from ebu_core import (LigandBuilder, EBUBuilder, autonomous_enumerate_sbus,
                      CN_GEOMETRIES)

# --- Ligand definitions ------------------------------------------------------
# BTC: benzene-1,3,5-tricarboxylic acid -> tricarboxylate (3 x carboxylate_O)
BTC_SMILES  = 'OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1'
# EDTA: ethylenediaminetetraacetic acid -> hexadentate N2O4
EDTA_SMILES = 'C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O'
# EDDA: ethylenediamine-N,N'-diacetic acid -> tetradentate N2O2
EDDA_SMILES = 'OC(=O)CNCCNCC(=O)O'


def make_ligands():
    btc  = LigandBuilder(BTC_SMILES,  name='BTC')
    edta = LigandBuilder(EDTA_SMILES, name='EDTA')
    edda = LigandBuilder(EDDA_SMILES, name='EDDA')
    return btc, edta, edda


def try_energy_model():
    """Return an MLEnergyModel instance if the MLIP stack is importable,
    else None (so the run degrades to geometry-only)."""
    try:
        from energy_model import MLEnergyModel
        model = MLEnergyModel()
        # force the lazy stack check now so we fail fast / clearly
        model._get_calculator()
        print('[energy] MACE-MP-0 stack available -- formation energies WILL be computed.')
        return model
    except Exception as e:
        print('[energy] MACE-MP-0 stack NOT available -- geometry + QC only.')
        print(f'         reason: {type(e).__name__}: {str(e).splitlines()[0]}')
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(os.getcwd(), 'ni_btc_outputs'))
    ap.add_argument('--max-copies', type=int, default=2,
                    help='max copies per ligand species (default 2)')
    ap.add_argument('--csv', default=os.path.join(os.getcwd(), 'ni_btc_results.csv'))
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    btc, edta, edda = make_ligands()
    ni = EBUBuilder('Ni', 2)
    model = try_energy_model()

    scenarios = {
        'A_BTC_only':      [btc],
        'B_BTC_EDTA':      [btc, edta],
        'C_BTC_EDDA':      [btc, edda],
        'D_BTC_EDTA_EDDA': [btc, edta, edda],
    }

    all_rows = []
    for tag, species in scenarios.items():
        print('\n' + '=' * 70)
        print(f'SCENARIO {tag}: {" + ".join(l.name for l in species)}')
        print('=' * 70)
        res = autonomous_enumerate_sbus(
            builder                = ni,
            ligand_species          = species,
            out_dir                 = args.out,
            label_prefix            = f'Ni2+_{tag}',
            max_copies_per_species  = args.max_copies,
            energy_model            = model,
            relax_with_energy_model = True,
            validate                = True,
            discard_invalid         = True,
        )
        for label, tup in res.items():
            mol, mol3d, qc, formation = tup
            all_rows.append({
                'scenario'     : tag,
                'label'        : label,
                'n_atoms'      : mol3d.GetNumAtoms(),
                'qc_valid'     : bool(qc.is_valid) if qc is not None else '',
                'min_clash'    : (round(qc.min_clash_ratio, 3)
                                  if (qc is not None and qc.min_clash_ratio is not None) else ''),
                'formation_eV' : (round(formation.formation_energy_eV, 4)
                                  if formation is not None else ''),
            })

    # write CSV
    with open(args.csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['scenario', 'label', 'n_atoms',
                                           'qc_valid', 'min_clash', 'formation_eV'])
        w.writeheader()
        w.writerows(all_rows)
    print(f'\n[done] {len(all_rows)} candidate EBUs -> {args.csv}')

    # global ranking if energies present
    ranked = [r for r in all_rows if r['formation_eV'] != '']
    if ranked:
        ranked.sort(key=lambda r: r['formation_eV'])
        print('\nGLOBAL RANKING by approximate formation energy (lowest = most stable):')
        for r in ranked:
            print(f"  {r['formation_eV']:11.4f} eV  [{r['scenario']}]  {r['label']}")
    else:
        print('\n[note] No formation energies computed (MLIP stack absent). '
              'Structures + QC are complete; install torch/mace and re-run to rank.')


if __name__ == '__main__':
    main()
