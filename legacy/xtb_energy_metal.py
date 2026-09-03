"""xtb_energy_metal.py -- GFN2-xTB formation-energy ranking for M(q)+BTC(+EDTA/EDDA).
Env: METAL_SYMBOL, METAL_CHARGE, OUT_TAG, METAL_MULT (fixed EBU+ion multiplicity;
if unset -> minimal spin), XTB_TIME_BUDGET. Resumable, chunked.
E_form = E(EBU,q) - E(M ion) - sum E(free ligand anion). Screening-grade."""
import os, csv, glob, re, json, warnings, time
warnings.filterwarnings('ignore')
from ase.io import read
from ase.optimize import LBFGS
from ase import Atoms
from tblite.ase import TBLite
from rdkit import Chem
from ebu_core import LigandBuilder, EBUExporter

SYM = os.environ.get('METAL_SYMBOL','Ni'); CHG=int(os.environ.get('METAL_CHARGE','2'))
TAG = os.environ.get('OUT_TAG', SYM.lower())
FIX_MULT = os.environ.get('METAL_MULT')  # str or None
FIX_MULT = int(FIX_MULT) if FIX_MULT else None
OUT_DIR = os.path.join(os.getcwd(), f'{TAG}_btc_outputs')
CSV     = os.path.join(os.getcwd(), f'{TAG}_btc_xtb_results.csv')
REF_CACHE = os.path.join(os.getcwd(), f'{TAG}_btc_refs.json')
FMAX, MAXSTEP = 0.20, 35
TIME_BUDGET = float(os.environ.get('XTB_TIME_BUDGET','26'))
LIG_SMILES = {'BTC':'OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1',
              'EDTA':'C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O',
              'EDDA':'OC(=O)CNCCNCC(=O)O'}

def minimal_mult(atoms,charge):
    ne=int(sum(atoms.get_atomic_numbers()))-charge; return 1 if ne%2==0 else 2
def ebu_mult(atoms,charge):
    return FIX_MULT if FIX_MULT else minimal_mult(atoms,charge)
def calc_for(charge,mult):
    return TBLite(method='GFN2-xTB',charge=charge,multiplicity=mult,verbosity=0,max_iterations=400,accuracy=1.0)
def relax(atoms,charge,mult):
    atoms.calc=calc_for(charge,mult); e0=atoms.get_potential_energy()
    opt=LBFGS(atoms,logfile=None); conv=opt.run(fmax=FMAX,steps=MAXSTEP)
    return e0,float(atoms.get_potential_energy()),bool(conv),opt.get_number_of_steps()
def lig_ref(name):
    lig=LigandBuilder(LIG_SMILES[name],name=name); q=Chem.GetFormalCharge(lig.mol)
    tmp=os.path.join(OUT_DIR,f'_ref_{name}.xyz'); open(tmp,'w').write(EBUExporter.to_xyz(lig.mol))
    a=read(tmp); _,e,conv,n=relax(a,q,minimal_mult(a,q))
    print(f'[ref] {name} q{q:+d}: E={e:.3f} ({n}st conv={conv})',flush=True); return e
def metal_ref():
    a=Atoms(SYM,positions=[[0,0,0]]); mult=FIX_MULT or (1 if (a.get_atomic_numbers()[0]-CHG)%2==0 else 2)
    a.calc=calc_for(CHG,mult); e=a.get_potential_energy()
    print(f'[ref] {SYM}{CHG}+ mult{mult}: E={e:.3f}',flush=True); return e
def get_refs():
    if os.path.exists(REF_CACHE):
        d=json.load(open(REF_CACHE)); print('[ref] cache',flush=True); return d['M'],d['lig']
    eM=metal_ref(); el={n:lig_ref(n) for n in ('BTC','EDTA','EDDA')}
    json.dump({'M':eM,'lig':el},open(REF_CACHE,'w')); return eM,el
def pq(b):
    m=re.search(r'_q([+-]\d+)',b); return int(m.group(1)) if m else 0
def pc(b):
    c={'BTC':0,'EDTA':0,'EDDA':0}
    for t in re.findall(r'(BTC|EDTA|EDDA)\[s\d+\]',b): c[t]+=1
    return c
def nb(b): return sum(int(s) for s in re.findall(r'\[s(\d+)\]',b))

def main():
    fields=['label','scenario','n_atoms','charge','n_bonds','E_sp_eV','E_relax_eV','converged','steps','E_form_eV','E_form_per_bond_eV']
    done=set()
    if os.path.exists(CSV):
        for r in csv.DictReader(open(CSV)):
            if r.get('E_form_eV') not in (None,''): done.add(r['label'])
    else:
        csv.DictWriter(open(CSV,'w',newline=''),fieldnames=fields).writeheader()
    files=sorted(glob.glob(os.path.join(OUT_DIR,f'{SYM}{CHG}+_*.xyz')))
    rem=[f for f in files if os.path.basename(f)[:-4] not in done]
    print(f'{len(files)} total, {len(done)} done, {len(rem)} left (budget {TIME_BUDGET}s).',flush=True)
    eM,el=get_refs(); t0=time.time()
    for f in rem:
        if time.time()-t0>TIME_BUDGET: print('[chunk] budget.',flush=True); break
        base=os.path.basename(f)[:-4]
        sc=re.match(rf'{re.escape(SYM)}{CHG}\+_([A-D])_',base); sc=sc.group(1) if sc else '?'
        q=pq(base); counts=pc(base); n=nb(base)
        try:
            a=read(f); e_sp,e_rx,conv,ns=relax(a,q,ebu_mult(a,q))
            ef=e_rx-(eM+sum(counts[k]*el[k] for k in counts))
            row=dict(label=base,scenario=sc,n_atoms=len(a),charge=q,n_bonds=n,E_sp_eV=round(e_sp,3),
                     E_relax_eV=round(e_rx,3),converged=conv,steps=ns,E_form_eV=round(ef,3),
                     E_form_per_bond_eV=round(ef/n,3) if n else '')
            print(f'  {base[:50]:50s} q{q:+d} Ef={ef:8.3f} ({ef/n:6.3f}/bond) {ns}st conv={conv}',flush=True)
        except Exception as ex:
            row=dict(label=base,scenario=sc,n_atoms='',charge=q,n_bonds=n,E_sp_eV='',E_relax_eV='',
                     converged='FAIL',steps='',E_form_eV='',E_form_per_bond_eV='')
            print(f'  {base[:50]:50s} FAILED: {str(ex)[:50]}',flush=True)
        csv.DictWriter(open(CSV,'a',newline=''),fieldnames=fields).writerow(row)
    rows=[r for r in csv.DictReader(open(CSV)) if r.get('E_form_eV') not in (None,'')]
    print(f'\n{len(rows)}/{len(files)} evaluated.',flush=True)
    if len(rows)>=len(files):
        rows.sort(key=lambda r:float(r['E_form_eV'])); print('ALL DONE.',flush=True)
        for r in rows: print(f"  {float(r['E_form_eV']):9.3f} [{r['scenario']}] {r['label']}",flush=True)

if __name__=='__main__': main()
