"""
solvation_model.py -- cluster-continuum solvation for the M/BTC/chelator EBUs.

NEAT SCHEME (how solvation fits in):
  * EXPLICIT first shell: the metal's octahedral coordination sphere is completed
    with real solvent molecules (H2O / EtOH / DMF, neutral, O-donor) on whatever
    sites BTC/chelator leave open. Solvent thus competes for the SAME sites the
    framework would grow into.
  * IMPLICIT bulk: everything is then evaluated in a GFN2-xTB/ALPB dielectric
    continuum for the chosen solvent.
  * MIXTURES (1:1, 1:1:1): explicit molecules split as evenly as possible across
    the chosen solvents; the continuum energy is the equal-weight average over the
    component solvents' ALPB single points (screening approximation).

Charge-conserving cluster-continuum formation energy in solvent system S:
  E_form,S(cluster) = E_S(cluster)
                      - E_S(M^q+) - sum E_S(ligand anion) - sum E_S(free solvent)
All terms in the same continuum S; each E_S(.) = mean over S's component solvents.

Species per metal (octahedral, CN target 6):
  ion    : [M(solv)6]              -- fully solvated free ion (solution baseline)
  node   : [M(BTC)2(solv)4]        -- solvated framework node (2 BTC, trans)
  edta   : [M(EDTA)]               -- hexadentate chelate (saturated, no solvent)
"""
import os, csv, json, warnings, itertools, time
warnings.filterwarnings('ignore')
import numpy as np
from copy import deepcopy
from ase import Atoms
from ase.optimize import LBFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem import AllChem
from ebu_core import LigandBuilder, EBUExporter

PT = Chem.GetPeriodicTable()
ALPB = {'W':'water','E':'ethanol','D':'dmf'}
SOLV_SMILES = {'W':'O','E':'CCO','D':'O=CN(C)C'}
SYSTEMS = {'W':['W'],'E':['E'],'D':['D'],'WE':['W','E'],'WD':['W','D'],'ED':['E','D'],'WED':['W','E','D']}
OCTA = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
D_ML = 2.05

def rot(a,b):
    a=a/np.linalg.norm(a); b=b/np.linalg.norm(b); v=np.cross(a,b); c=float(np.dot(a,b)); s=np.linalg.norm(v)
    if s<1e-8: return np.eye(3) if c>0 else -np.eye(3)
    K=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    return np.eye(3)+K+K@K*(1-c)/s**2

def solvent_frag(key):
    m=Chem.AddHs(Chem.MolFromSmiles(SOLV_SMILES[key])); AllChem.EmbedMolecule(m,AllChem.ETKDGv3())
    try: AllChem.MMFFOptimizeMolecule(m)
    except Exception: pass
    o=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==8][0]
    syms=[a.GetSymbol() for a in m.GetAtoms()]; pos=m.GetConformer().GetPositions()
    return syms,pos,o,0  # charge 0

def btc_frag():
    lig=LigandBuilder(SOLV_SMILES.get('__',''),name='x') if False else LigandBuilder('OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1',name='BTC')
    m=lig.mol; o=lig.donor_indices[0]
    syms=[a.GetSymbol() for a in m.GetAtoms()]; pos=m.GetConformer().GetPositions()
    return syms,pos,o,Chem.GetFormalCharge(m)

def place(frag, site_hat):
    syms,pos,o,q=frag; pos=np.array(pos,float)
    heavy=np.array([p for s,p in zip(syms,pos) if s!='H'])
    body=heavy.mean(0)-pos[o]
    if np.linalg.norm(body)<1e-3: body=pos.mean(0)-pos[o]
    R=rot(body, site_hat)
    pos2=(R@(pos-pos[o]).T).T + site_hat*D_ML
    return syms,pos2,q

def assemble(metal, mcharge, occupants):
    """occupants: list of frag for each octahedral site used (in OCTA order)."""
    syms=[metal]; pos=[[0,0,0]]; q=mcharge
    for i,frag in enumerate(occupants):
        s,p,fq=place(frag, OCTA[i]); syms+=s; pos+=[list(x) for x in p]; q+=fq
    return syms,np.array(pos,float),q

def to_atoms(syms,pos): a=Atoms(symbols=syms,positions=pos); a.pbc=False; return a

def mult_for(syms,charge,metal,fixmult):
    ne=sum(PT.GetAtomicNumber(s) for s in syms)-charge
    if fixmult: return fixmult
    return 1 if ne%2==0 else 2

def energy_S(syms,pos,charge,metal,fixmult,comps,relax=True):
    """mean ALPB energy over component solvents; relax under first component."""
    es=[]; relaxed_pos=pos
    for j,c in enumerate(comps):
        a=to_atoms(syms,relaxed_pos); mult=mult_for(syms,charge,metal,fixmult)
        a.calc=TBLite(method='GFN2-xTB',charge=charge,multiplicity=mult,verbosity=0,
                      max_iterations=400,accuracy=1.0,solvation=("alpb",ALPB[c]))
        if relax and j==0:
            LBFGS(a,logfile=None).run(fmax=0.20,steps=int(os.environ.get('SOLV_STEPS','30')))
            relaxed_pos=a.get_positions()
        es.append(float(a.get_potential_energy()))
    return float(np.mean(es)), relaxed_pos

# ---------------- driver ----------------
CSV = os.path.join(os.getcwd(),'solvation_raw.csv')
FIELDS=['kind','metal','name','comp','charge','n_btc','n_edta','n_W','n_E','n_D','E_W','E_E','E_D']

def distribute(n_sites, comps):
    """round-robin assign solvent types (keys) to n_sites -> list of keys."""
    return [comps[i % len(comps)] for i in range(n_sites)]

def solv_counts(keys):
    return {k:keys.count(k) for k in ('W','E','D')}

def three_alpb(syms,pos,charge,metal,fixmult, relax_gas=True):
    """gas-relax (light) once, then ALPB single points under W,E,D. Returns dict."""
    # gas relax
    a=to_atoms(syms,pos); mult=mult_for(syms,charge,metal,fixmult)
    a.calc=TBLite(method='GFN2-xTB',charge=charge,multiplicity=mult,verbosity=0,max_iterations=400,accuracy=1.0)
    if relax_gas:
        LBFGS(a,logfile=None).run(fmax=0.25,steps=int(os.environ.get('SOLV_STEPS','25')))
    gpos=a.get_positions(); out={}
    for c in ('W','E','D'):
        b=to_atoms(syms,gpos)
        b.calc=TBLite(method='GFN2-xTB',charge=charge,multiplicity=mult,verbosity=0,
                      max_iterations=400,accuracy=1.0,solvation=("alpb",ALPB[c]))
        out[c]=float(b.get_potential_energy())
    return out

def load_xyz(path):
    lines=open(path).read().splitlines(); n=int(lines[0]); syms=[]; pos=[]
    for l in lines[2:2+n]:
        p=l.split(); syms.append(p[0]); pos.append([float(p[1]),float(p[2]),float(p[3])])
    return syms,np.array(pos,float)

def main():
    import glob
    metals=[('Ni',2,None),('Fe',3,6)]
    # work list
    done=set()
    if os.path.exists(CSV):
        for r in csv.DictReader(open(CSV)): done.add((r['kind'],r['metal'],r['name'],r['comp']))
    else:
        csv.DictWriter(open(CSV,'w',newline=''),fieldnames=FIELDS).writeheader()
    def emit(row): csv.DictWriter(open(CSV,'a',newline=''),fieldnames=FIELDS).writerow(row)

    t0=time.time(); budget=float(os.environ.get('SOLV_BUDGET','26'))
    def over(): return time.time()-t0>budget

    # ---- solvent-independent free refs (metal-independent): free solvent molecules
    for k in ('W','E','D'):
        if over(): print('[budget]'); return
        if ('ref','-','solv_'+k,'-') in done: continue
        s,p,o,q=solvent_frag(k); e=three_alpb(s,p,0,None,None)
        emit(dict(kind='ref',metal='-',name='solv_'+k,comp='-',charge=0,n_btc=0,n_edta=0,
                  n_W=0,n_E=0,n_D=0,E_W=round(e['W'],4),E_E=round(e['E'],4),E_D=round(e['D'],4)))
        print(f'ref solv {k} done',flush=True)

    for metal,mq,fm in metals:
        # ---- refs: bare ion, BTC anion, EDTA anion
        if not over() and ('ref',metal,'ion','-') not in done:
            e=three_alpb([metal],np.zeros((1,3)),mq,metal,fm,relax_gas=False)
            emit(dict(kind='ref',metal=metal,name='ion',comp='-',charge=mq,n_btc=0,n_edta=0,n_W=0,n_E=0,n_D=0,
                      E_W=round(e['W'],4),E_E=round(e['E'],4),E_D=round(e['D'],4))); print(f'{metal} ion ref',flush=True)
        if not over() and ('ref',metal,'BTC','-') not in done:
            s,p,o,q=btc_frag(); e=three_alpb(s,p,q,metal,None)
            emit(dict(kind='ref',metal=metal,name='BTC',comp='-',charge=q,n_btc=1,n_edta=0,n_W=0,n_E=0,n_D=0,
                      E_W=round(e['W'],4),E_E=round(e['E'],4),E_D=round(e['D'],4))); print(f'{metal} BTC ref',flush=True)
        if not over() and ('ref',metal,'EDTA','-') not in done:
            lig=LigandBuilder('C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O',name='EDTA'); m=lig.mol
            s=[a.GetSymbol() for a in m.GetAtoms()]; p=m.GetConformer().GetPositions(); q=Chem.GetFormalCharge(m)
            e=three_alpb(s,p,q,metal,None)
            emit(dict(kind='ref',metal=metal,name='EDTA',comp='-',charge=q,n_btc=0,n_edta=1,n_W=0,n_E=0,n_D=0,
                      E_W=round(e['W'],4),E_E=round(e['E'],4),E_D=round(e['D'],4))); print(f'{metal} EDTA ref',flush=True)

        # ---- species with explicit shells
        comps_sets={'W':['W'],'E':['E'],'D':['D'],'WE':['W','E'],'WD':['W','D'],'ED':['E','D'],'WED':['W','E','D']}
        for comp_label,comps in comps_sets.items():
            # ion [M(solv)6]
            if not over() and ('cluster',metal,'ion','_'+comp_label) not in done:
                keys=distribute(6,comps); frags=[solvent_frag(k)[:3]+ (0,) for k in keys]
                frags=[ (solvent_frag(k)) for k in keys]
                syms,pos,q=assemble(metal,mq,[ (f[0],f[1],f[2],f[3]) for f in frags])
                e=three_alpb(syms,pos,q,metal,fm); c=solv_counts(keys)
                emit(dict(kind='cluster',metal=metal,name='ion',comp='_'+comp_label,charge=q,n_btc=0,n_edta=0,
                          n_W=c['W'],n_E=c['E'],n_D=c['D'],E_W=round(e['W'],4),E_E=round(e['E'],4),E_D=round(e['D'],4)))
                print(f'{metal} ion {comp_label} q{q}',flush=True)
            # node [M(BTC)2(solv)4]
            if not over() and ('cluster',metal,'node','_'+comp_label) not in done:
                keys=distribute(4,comps); occ=[btc_frag(),btc_frag()]+[solvent_frag(k) for k in keys]
                syms,pos,q=assemble(metal,mq,occ); e=three_alpb(syms,pos,q,metal,fm); c=solv_counts(keys)
                emit(dict(kind='cluster',metal=metal,name='node',comp='_'+comp_label,charge=q,n_btc=2,n_edta=0,
                          n_W=c['W'],n_E=c['E'],n_D=c['D'],E_W=round(e['W'],4),E_E=round(e['E'],4),E_D=round(e['D'],4)))
                print(f'{metal} node {comp_label} q{q}',flush=True)
        # edta cluster (saturated, no solvent) -- one geometry, 3 ALPB SP
        if not over() and ('cluster',metal,'edta','_none') not in done:
            d=os.path.join(os.getcwd(),f'{metal.lower()}_btc_outputs')
            fs=[os.path.join(d,fn) for fn in os.listdir(d) if 'EDTA[s6]_octahedral' in fn and fn.startswith(f'{metal}{mq}+')]
            if fs:
                s,p=load_xyz(fs[0]); q=mq-4
                e=three_alpb(s,p,q,metal,fm)
                emit(dict(kind='cluster',metal=metal,name='edta',comp='_none',charge=q,n_btc=0,n_edta=1,n_W=0,n_E=0,n_D=0,
                          E_W=round(e['W'],4),E_E=round(e['E'],4),E_D=round(e['D'],4)))
                print(f'{metal} edta q{q}',flush=True)
    print('PASS complete (or budget).',flush=True)

if __name__=='__main__': main()
