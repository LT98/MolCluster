"""
xtb_energy.py -- screening-level formation-energy ranking of the Ni+BTC(+EDTA
+EDDA) EBUs using GFN2-xTB (semi-empirical tight binding, via tblite/ASE).

MACE (energy_model.py) needs a torch/CUDA stack that will not fit in this
sandbox. GFN2-xTB is a legitimate screening alternative: a real self-consistent
quantum method with Ni parameters that treats total charge and spin EXPLICITLY
(important for these charged complexes). Screening-grade only -> relative ranking.

Formation energy (charge-conserving):
    E_form = E(EBU, q_ebu) - E(Ni2+) - sum_i E(free ligand_i anion, q_i)
each ligand fully deprotonated as ebu_core built it (BTC -3, EDTA -4, EDDA -2),
q_ebu = +2 + sum q_i. Minimal-spin convention (mult 1 if even electrons else 2;
Ni2+ reference = triplet d8). Resumable + time-budgeted for chunked running.
"""
import os, csv, glob, re, json, warnings, time
warnings.filterwarnings('ignore')
from ase.io import read
from ase.optimize import LBFGS
from ase import Atoms
from tblite.ase import TBLite
from rdkit import Chem
from ebu_core import LigandBuilder, EBUExporter

OUT_DIR   = os.path.join(os.getcwd(), 'ni_btc_outputs')
CSV       = os.path.join(os.getcwd(), 'ni_btc_xtb_results.csv')
REF_CACHE = os.path.join(os.getcwd(), 'ni_btc_refs.json')
FMAX, MAXSTEP = 0.20, 50
TIME_BUDGET = float(os.environ.get('XTB_TIME_BUDGET', '38'))

LIG_SMILES = {
    'BTC':  'OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1',
    'EDTA': 'C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O',
    'EDDA': 'OC(=O)CNCCNCC(=O)O',
}

def minimal_multiplicity(atoms, charge):
    ne = int(sum(atoms.get_atomic_numbers())) - charge
    return 1 if ne % 2 == 0 else 2

def calc_for(charge, mult):
    return TBLite(method='GFN2-xTB', charge=charge, multiplicity=mult,
                  verbosity=0, max_iterations=400, accuracy=1.0)

def relax(atoms, charge, mult, fmax=FMAX, steps=MAXSTEP):
    atoms.calc = calc_for(charge, mult)
    e0 = atoms.get_potential_energy()
    opt = LBFGS(atoms, logfile=None)
    conv = opt.run(fmax=fmax, steps=steps)
    return e0, float(atoms.get_potential_energy()), bool(conv), opt.get_number_of_steps()

def ligand_ref_energy(name):
    lig = LigandBuilder(LIG_SMILES[name], name=name)
    q = Chem.GetFormalCharge(lig.mol)
    tmp = os.path.join(OUT_DIR, f'_ref_{name}.xyz')
    open(tmp, 'w').write(EBUExporter.to_xyz(lig.mol))
    a = read(tmp); mult = minimal_multiplicity(a, q)
    _, e, conv, n = relax(a, q, mult)
    print(f'[ref] {name} q{q:+d} mult{mult}: E={e:.3f} eV (conv={conv}, {n} steps)', flush=True)
    return e

def ni_ref_energy():
    a = Atoms('Ni', positions=[[0, 0, 0]]); a.calc = calc_for(2, 3)
    e = a.get_potential_energy()
    print(f'[ref] Ni2+ mult3: E={e:.3f} eV', flush=True)
    return e

def parse_charge(b):
    m = re.search(r'_q([+-]\d+)', b); return int(m.group(1)) if m else 0

def parse_ligand_counts(b):
    c = {'BTC': 0, 'EDTA': 0, 'EDDA': 0}
    for tok in re.findall(r'(BTC|EDTA|EDDA)\[s\d+\]', b): c[tok] += 1
    return c

def n_donor_bonds(b):
    return sum(int(s) for s in re.findall(r'\[s(\d+)\]', b))

def get_references():
    if os.path.exists(REF_CACHE):
        d = json.load(open(REF_CACHE)); print('[ref] loaded cache', flush=True)
        return d['Ni'], d['lig']
    e_ni = ni_ref_energy()
    e_lig = {n: ligand_ref_energy(n) for n in ('BTC', 'EDTA', 'EDDA')}
    json.dump({'Ni': e_ni, 'lig': e_lig}, open(REF_CACHE, 'w'))
    return e_ni, e_lig

def main():
    fields = ['label', 'scenario', 'n_atoms', 'charge', 'n_bonds', 'E_sp_eV',
              'E_relax_eV', 'converged', 'steps', 'E_form_eV', 'E_form_per_bond_eV']
    done = set()
    if os.path.exists(CSV):
        for r in csv.DictReader(open(CSV)):
            if r.get('E_form_eV') not in (None, ''): done.add(r['label'])
    else:
        csv.DictWriter(open(CSV, 'w', newline=''), fieldnames=fields).writeheader()

    files = sorted(glob.glob(os.path.join(OUT_DIR, 'Ni2+_*.xyz')))
    remaining = [f for f in files if os.path.basename(f)[:-4] not in done]
    print(f'{len(files)} total, {len(done)} done, {len(remaining)} remaining (budget {TIME_BUDGET}s).', flush=True)

    e_ni, e_lig = get_references()
    t0 = time.time()
    for f in remaining:
        if time.time() - t0 > TIME_BUDGET:
            print('[chunk] budget reached.', flush=True); break
        base = os.path.basename(f)[:-4]
        sc = re.match(r'Ni2\+_([A-D])_', base); sc = sc.group(1) if sc else '?'
        q = parse_charge(base); counts = parse_ligand_counts(base); nb = n_donor_bonds(base)
        try:
            a = read(f); mult = minimal_multiplicity(a, q)
            e_sp, e_rx, conv, nstep = relax(a, q, mult)
            e_form = e_rx - (e_ni + sum(counts[k] * e_lig[k] for k in counts))
            row = dict(label=base, scenario=sc, n_atoms=len(a), charge=q, n_bonds=nb,
                       E_sp_eV=round(e_sp, 3), E_relax_eV=round(e_rx, 3), converged=conv,
                       steps=nstep, E_form_eV=round(e_form, 3),
                       E_form_per_bond_eV=round(e_form / nb, 3) if nb else '')
            print(f'  {base[:52]:52s} q{q:+d} Ef={e_form:8.3f} ({e_form/nb:6.3f}/bond) {nstep}st conv={conv}', flush=True)
        except Exception as ex:
            row = dict(label=base, scenario=sc, n_atoms='', charge=q, n_bonds=nb,
                       E_sp_eV='', E_relax_eV='', converged='FAIL', steps='',
                       E_form_eV='', E_form_per_bond_eV='')
            print(f'  {base[:52]:52s} FAILED: {str(ex)[:55]}', flush=True)
        csv.DictWriter(open(CSV, 'a', newline=''), fieldnames=fields).writerow(row)

    rows = [r for r in csv.DictReader(open(CSV)) if r.get('E_form_eV') not in (None, '')]
    print(f'\n{len(rows)}/{len(files)} EBUs evaluated.', flush=True)
    if len(rows) >= len(files):
        rows.sort(key=lambda r: float(r['E_form_eV']))
        print('ALL DONE. Ranking (lowest E_form = most stable):', flush=True)
        for r in rows:
            print(f"  {float(r['E_form_eV']):9.3f} eV  [{r['scenario']}] {r['label']}", flush=True)

if __name__ == '__main__':
    main()
