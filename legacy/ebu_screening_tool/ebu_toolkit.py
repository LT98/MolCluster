"""
ebu_toolkit.py -- one-stop toolkit wrapping the EBU pipeline into a few
customizable functions, for use from ebu_tool.ipynb as a tool with editable
inputs (metal centres, linkers, additives, solvents) and 3D visualization.

Layers (all screening-grade, GFN2-xTB):
  1. enumerate_models(...)      geometry enumeration + QC  (needs rdkit only)
  2. formation_energies(...)    gas or implicit-solvent E_form ranking  (needs tblite/ase)
  3. solvation_sweep(...)       explicit+implicit cluster-continuum over solvents
  4. view_model(...)            interactive 3D of any generated model  (needs py3Dmol)

Nothing here is chunked/time-budgeted -- that machinery was only for the build
sandbox; in a normal Jupyter kernel these just run to completion.
"""
from __future__ import annotations
import os, re, itertools, warnings
from copy import deepcopy
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import AllChem
from ebu_core import (LigandBuilder, EBUBuilder, autonomous_enumerate_sbus,
                      CN_GEOMETRIES, METAL_CHARGE)

PT = Chem.GetPeriodicTable()

# ---------------------------------------------------------------- libraries
LIGAND_LIBRARY = {   # name -> SMILES (neutral / protonated form)
    'BTC':  'C1=C(C=C(C=C1C(=O)O)C(=O)O)C(=O)O',      # benzene-1,3,5-tricarboxylic acid
    'BDC':  'OC(=O)c1ccc(C(=O)O)cc1',              # terephthalic acid
    'EDTA': 'C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O',
    'EDDA': 'OC(=O)CNCCNCC(=O)O',                  # ethylenediamine-N,N'-diacetic acid
    'bipy': 'c1cc(-c2ccncc2)ccn1',                 # 4,4'-bipyridine
}
SOLVENTS   = {'water': 'O', 'ethanol': 'CCO', 'dmf': 'O=CN(C)C'}
ALPB_NAME  = {'water': 'water', 'ethanol': 'ethanol', 'dmf': 'dmf'}

# high-spin d-electron count for common MOF metal ions -> unpaired e- (high spin)
_D_ELECTRONS = {  # (symbol,charge): d-count
    ('Ni',2):8,('Fe',3):5,('Fe',2):6,('Cu',2):9,('Co',2):7,('Mn',2):5,
    ('Zn',2):10,('Cr',3):3,('Cd',2):10,('Zr',4):0,('Hf',4):0,('Al',3):0,
    ('Ti',4):0,('V',3):2,('In',3):0,('Mg',2):0,('Ca',2):0,
}
def _highspin_unpaired(symbol, charge):
    d = _D_ELECTRONS.get((symbol, charge))
    if d is None: return None
    return d if d <= 5 else 10 - d

# module-level registry so the viz cell can find generated structures
MODELS: dict = {}    # label -> dict(mol, mol3d, qc, xyz_path, metal, charge, n_bonds, scenario)


# ---------------------------------------------------------------- helpers
def build_ligand(name, smiles=None, **kw):
    smi = smiles or LIGAND_LIBRARY.get(name)
    if smi is None:
        raise KeyError(f'Unknown ligand {name!r}; pass smiles= or add it to LIGAND_LIBRARY.')
    return LigandBuilder(smi, name=name, **kw)

def _scenarios(additives, mode):
    """mode='subsets' -> linker-only + every non-empty additive subset;
       mode='all'     -> single scenario using all additives together;
       mode='none'    -> linker-only."""
    if mode == 'none' or not additives:
        return {'linker_only': []}
    if mode == 'all':
        return {'+'.join(additives): list(additives)}
    out = {'linker_only': []}
    for r in range(1, len(additives) + 1):
        for combo in itertools.combinations(additives, r):
            out['+'.join(combo)] = list(combo)
    return out

def metal_multiplicity(symbol, charge, atoms_or_syms, spin_mode='auto'):
    if hasattr(atoms_or_syms, 'get_atomic_numbers'):
        ztot = int(sum(atoms_or_syms.get_atomic_numbers()))
    else:
        ztot = int(sum(PT.GetAtomicNumber(s) for s in atoms_or_syms))
    ne = ztot - charge
    if spin_mode == 'highspin':
        u = _highspin_unpaired(symbol, charge)
        if u is not None:
            # match parity of total electron count
            if (ne - u) % 2 != 0: u += 1
            return u + 1
    return 1 if ne % 2 == 0 else 2


# ---------------------------------------------------------------- 1. enumerate
def enumerate_models(metals, linkers=('BTC',), additives=(), additive_mode='subsets',
                     max_copies=2, clash_scale=0.65, geometries=CN_GEOMETRIES,
                     out_root='tool_outputs', verbose=False):
    """
    metals    : list of symbols ['Ni','Fe'] or (symbol,charge) tuples.
    linkers   : framework linker names/SMILES (from LIGAND_LIBRARY or 'Name:SMILES').
    additives : chelating/modulating additive names.
    additive_mode : 'subsets' (default) | 'all' | 'none'.
    Returns a pandas DataFrame; also fills module-level MODELS registry.
    """
    os.makedirs(out_root, exist_ok=True)
    def _resolve(names):
        libs = []
        for n in names:
            if ':' in n: nm, smi = n.split(':', 1); libs.append(build_ligand(nm, smi))
            else: libs.append(build_ligand(n))
        return libs
    linker_ligs = _resolve(linkers)
    add_ligs    = {n: build_ligand(n) for n in additives}
    rows = []
    for m in metals:
        sym, chg = (m if isinstance(m, (tuple, list)) else (m, METAL_CHARGE.get(m, 2)))
        builder = EBUBuilder(sym, chg)
        for scen, adds in _scenarios(list(additives), additive_mode).items():
            species = list(linker_ligs) + [add_ligs[a] for a in adds]
            out_dir = os.path.join(out_root, f'{sym}{chg}')
            res = autonomous_enumerate_sbus(
                builder=builder, ligand_species=species, out_dir=out_dir,
                label_prefix=f'{sym}{chg}+_{scen}', max_copies_per_species=max_copies,
                clash_scale=clash_scale, geometries=geometries,
                validate=True, discard_invalid=True)
            for label, (mol, mol3d, qc, _) in res.items():
                q = _parse_q(label); nb = _n_bonds(label)
                chel = _lig_bonds(label, set(additives))
                xyzp = os.path.join(out_dir, label + '.xyz')
                MODELS[label] = dict(mol=mol, mol3d=mol3d, qc=qc, xyz_path=xyzp,
                                     metal=sym, charge=q, n_bonds=nb, scenario=scen,
                                     chel_bonds=chel, chel_extent=(chel/nb if nb else 0.0))
                rows.append(dict(metal=sym, scenario=scen, label=label,
                                 n_atoms=mol3d.GetNumAtoms(), charge=q, n_bonds=nb,
                                 chelation_extent=round(chel/nb, 3) if nb else 0.0,
                                 qc_valid=bool(qc.is_valid) if qc else None))
    df = pd.DataFrame(rows)
    if verbose: print(f'{len(df)} models generated across {len(metals)} metal(s).')
    return df

def _parse_q(lbl):
    m = re.search(r'_q([+-]\d+)', lbl); return int(m.group(1)) if m else None
def _n_bonds(lbl):
    return sum(int(s) for s in re.findall(r'\[s(\d+)\]', lbl))
def _lig_bonds(lbl, names):
    """Total metal bonds (sum of [sN]) contributed by ligand tokens in `names`."""
    return sum(int(s) for nm, s in re.findall(r'([A-Za-z]+)\[s(\d+)\]', lbl) if nm in names)

def _lig_counts(lbl, names):
    c = {n: 0 for n in names}
    for nm, s in re.findall(r'([A-Za-z]+)\[s(\d+)\]', lbl):
        if nm in c: c[nm] += 1
    return c


# ---------------------------------------------------------------- 2. energies
_TBLITE_MSG = (
    "GFN2-xTB backend 'tblite' is not installed in this kernel.\n"
    "  Linux / macOS : pip install tblite\n"
    "  Windows       : conda install -c conda-forge tblite-python\n"
    "(then restart the Jupyter kernel). Enumeration and 3D viewing work without it;\n"
    "only the energy-ranking and solvation cells need it."
)
def require_tblite():
    try:
        import tblite.ase  # noqa
    except ImportError as e:
        raise RuntimeError(_TBLITE_MSG) from e

def _tblite(charge, mult, solvent=None, method='GFN2-xTB'):
    require_tblite()
    from tblite.ase import TBLite
    kw = dict(method=method, charge=charge, multiplicity=mult, verbosity=0,
              max_iterations=400, accuracy=1.0)
    if solvent: kw['solvation'] = ('alpb', ALPB_NAME[solvent])
    return TBLite(**kw)

def _relax_energy(atoms, charge, mult, solvent=None, relax=True, fmax=0.2, steps=60, method='GFN2-xTB'):
    from ase.optimize import LBFGS
    atoms = atoms.copy(); atoms.calc = _tblite(charge, mult, solvent, method)
    if relax: LBFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
    return float(atoms.get_potential_energy())

def _mol_to_atoms(mol3d):
    from ase import Atoms
    conf = mol3d.GetConformer()
    syms = [a.GetSymbol() for a in mol3d.GetAtoms()]
    pos  = np.array([list(conf.GetAtomPosition(i)) for i in range(mol3d.GetNumAtoms())])
    a = Atoms(symbols=syms, positions=pos); a.pbc = False; return a

_REF_CACHE = {}
def _ligand_ref(name, solvent, spin_mode, relax, method):
    key = ('lig', name, solvent)
    if key in _REF_CACHE: return _REF_CACHE[key]
    lig = build_ligand(name); a = _mol_to_atoms(lig.mol); q = Chem.GetFormalCharge(lig.mol)
    mult = 1 if (int(sum(a.get_atomic_numbers())) - q) % 2 == 0 else 2
    e = _relax_energy(a, q, mult, solvent, relax, method=method)
    _REF_CACHE[key] = (e, q); return e, q
def _metal_ref(sym, charge, solvent, spin_mode, method):
    from ase import Atoms
    key = ('metal', sym, charge, solvent, spin_mode)
    if key in _REF_CACHE: return _REF_CACHE[key]
    a = Atoms(sym, positions=[[0,0,0]]); mult = metal_multiplicity(sym, charge, a, spin_mode)
    e = _relax_energy(a, charge, mult, solvent, relax=False, method=method)
    _REF_CACHE[key] = e; return e

def formation_energies(df=None, solvent=None, spin_mode='auto', relax=True,
                       fmax=0.2, steps=60, method='GFN2-xTB', verbose=True):
    """
    Screening formation energy for each model: E(EBU) - E(metal ion) - sum E(free ligand anion).
    solvent : None (gas) or 'water'/'ethanol'/'dmf' (implicit ALPB).
    spin_mode : 'auto' (minimal) or 'highspin' (use metal high-spin d-count).
    Returns a ranked DataFrame (also usable as the input df to filter which models to run).
    """
    require_tblite()
    labels = list(MODELS) if df is None else list(df['label'])
    lig_names = list(LIGAND_LIBRARY)
    out = []
    for i, lbl in enumerate(labels, 1):
        info = MODELS[lbl]; sym = info['metal']; q = info['charge']; nb = info['n_bonds']
        counts = _lig_counts(lbl, lig_names)
        try:
            a = _mol_to_atoms(info['mol3d']); mult = metal_multiplicity(sym, q, a, spin_mode)
            e_ebu = _relax_energy(a, q, mult, solvent, relax, fmax, steps, method)
            e_ref = _metal_ref(sym, _metal_charge_of(lbl, sym), solvent, spin_mode, method)
            for nm, n in counts.items():
                if n: e_ref += n * _ligand_ref(nm, solvent, spin_mode, relax, method)[0]
            ef = e_ebu - e_ref
            out.append(dict(label=lbl, metal=sym, scenario=info['scenario'], charge=q,
                            n_bonds=nb, chelation_extent=round(info.get('chel_extent',0.0),3),
                            E_form_eV=round(ef,3),
                            E_form_per_bond_eV=round(ef/nb,3) if nb else None))
            if verbose: print(f'[{i}/{len(labels)}] {lbl[:48]:48s} Ef={ef:8.3f} eV')
        except Exception as ex:
            out.append(dict(label=lbl, metal=sym, scenario=info['scenario'], charge=q,
                            n_bonds=nb, chelation_extent=round(info.get('chel_extent',0.0),3),
                            E_form_eV=None, E_form_per_bond_eV=None))
            if verbose: print(f'[{i}/{len(labels)}] {lbl[:48]:48s} FAILED {str(ex)[:40]}')
    res = pd.DataFrame(out).sort_values('E_form_eV', na_position='last').reset_index(drop=True)
    return res

def _metal_charge_of(lbl, sym):
    m = re.match(rf'{re.escape(sym)}(\d+)\+', lbl); return int(m.group(1)) if m else METAL_CHARGE.get(sym, 2)


# ---------------------------------------------------------------- 3. solvation (cluster-continuum)
_OCTA = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], float)
def _rot(a, b):
    a=a/np.linalg.norm(a); b=b/np.linalg.norm(b); v=np.cross(a,b); c=float(np.dot(a,b)); s=np.linalg.norm(v)
    if s<1e-8: return np.eye(3) if c>0 else -np.eye(3)
    K=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]); return np.eye(3)+K+K@K*(1-c)/s**2
def _frag_from_smiles(smi):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m,AllChem.ETKDGv3())
    try: AllChem.MMFFOptimizeMolecule(m)
    except Exception: pass
    o=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==8][0]
    return [a.GetSymbol() for a in m.GetAtoms()], m.GetConformer().GetPositions(), o, Chem.GetFormalCharge(m)
def _frag_from_ligand(name):
    lig=build_ligand(name); m=lig.mol
    return [a.GetSymbol() for a in m.GetAtoms()], m.GetConformer().GetPositions(), lig.donor_indices[0], Chem.GetFormalCharge(m)
def _place(frag, site_hat, d_ml=2.05):
    syms,pos,o,q=frag; pos=np.array(pos,float)
    heavy=np.array([p for s,p in zip(syms,pos) if s!='H']); body=heavy.mean(0)-pos[o]
    if np.linalg.norm(body)<1e-3: body=pos.mean(0)-pos[o]
    R=_rot(body, site_hat); return syms,(R@(pos-pos[o]).T).T+site_hat*d_ml,q
def _assemble(metal, mq, occupants):
    from ase import Atoms
    syms=[metal]; pos=[[0,0,0]]; q=mq
    for i,fr in enumerate(occupants):
        s,p,fq=_place(fr,_OCTA[i]); syms+=s; pos+=[list(x) for x in p]; q+=fq
    return syms, np.array(pos,float), q
def _three_alpb(syms, pos, charge, sym, spin_mode, relax=True, steps=25, method='GFN2-xTB'):
    from ase import Atoms; from ase.optimize import LBFGS
    a=Atoms(symbols=syms,positions=pos); a.pbc=False
    mult=metal_multiplicity(sym,charge,a,spin_mode)
    a.calc=_tblite(charge,mult,None,method)
    if relax: LBFGS(a,logfile=None).run(fmax=0.25,steps=steps)
    g=a.get_positions(); out={}
    for c in ('water','ethanol','dmf'):
        b=Atoms(symbols=syms,positions=g); b.pbc=False; b.calc=_tblite(charge,mult,c,method)
        out[c]=float(b.get_potential_energy())
    return out

def _chelate_cluster(sym, mq, chelator):
    """Build a properly-wrapped [M(chelator)] complex via the pipeline placer."""
    from ebu_core import GeometryPlacer
    lig=build_ligand(chelator); k=min(lig.n_sites,6); lig=lig.with_sites(k)
    geom='octahedral' if k==6 else ('tetrahedral' if k==4 else 'planar')
    builder=EBUBuilder(sym,mq); mol=builder.build([lig],sites_per_ligand=[k])
    mol3d=GeometryPlacer.place(mol,sym,geom,[lig],d_ml=2.05)
    conf=mol3d.GetConformer()
    syms=[a.GetSymbol() for a in mol3d.GetAtoms()]
    pos=np.array([list(conf.GetAtomPosition(i)) for i in range(mol3d.GetNumAtoms())])
    q=mq+Chem.GetFormalCharge(build_ligand(chelator).mol)
    return syms,pos,q

def solvation_sweep(metals, linker='BTC', chelator='EDTA',
                    systems=('water','ethanol','dmf','water+ethanol','water+dmf','ethanol+dmf','water+ethanol+dmf'),
                    spin_mode='highspin', steps=25):
    """
    Cluster-continuum comparison of [M(solv)6] vs [M(linker)2(solv)4] vs [M(chelator)] across
    solvent systems (equal-ratio explicit shell + averaged ALPB continuum). Returns summary df
    with nucleation driving force and sequestration margin per metal x solvent system.
    """
    require_tblite()
    def comps(sysk): return sysk.split('+')
    def distribute(n, cs): return [cs[i % len(cs)] for i in range(n)]
    rows=[]
    for m in metals:
        sym,mq=(m if isinstance(m,(tuple,list)) else (m,METAL_CHARGE.get(m,2)))
        # references (3 pure-solvent single points each)
        eM={c:_three_alpb([sym],np.zeros((1,3)),mq,sym,spin_mode,relax=False)[c] for c in ('water','ethanol','dmf')}
        eLk=_three_alpb(*(_frag_from_ligand(linker)[:2]),_frag_from_ligand(linker)[3],sym,spin_mode,steps=steps)
        eCh=_three_alpb(*(_frag_from_ligand(chelator)[:2]),_frag_from_ligand(chelator)[3],sym,spin_mode,steps=steps)
        eS={s:_three_alpb(*_frag_from_smiles(SOLVENTS[s])[:2],0,sym,spin_mode,steps=steps) for s in SOLVENTS}
        lk_q=_frag_from_ligand(linker)[3]; ch_q=_frag_from_ligand(chelator)[3]
        # chelator cluster (saturated) — reuse its free geometry as proxy shell (hexadentate wraps metal)
        for sysk in systems:
            cs=comps(sysk)
            mean=lambda d: float(np.mean([d[c] for c in cs]))
            # ion [M(solv)6]
            keys=distribute(6,cs); ionc=_three_alpb(*_assemble(sym,mq,[_frag_from_smiles(SOLVENTS[k]) for k in keys])[:2],
                                                    mq,sym,spin_mode,steps=steps)
            # node [M(linker)2(solv)4]
            keys4=distribute(4,cs); occ=[_frag_from_ligand(linker),_frag_from_ligand(linker)]+[_frag_from_smiles(SOLVENTS[k]) for k in keys4]
            nsy,npo,nq=_assemble(sym,mq,occ); nodec=_three_alpb(nsy,npo,nq,sym,spin_mode,steps=steps)
            # chelator [M(chelator)] octahedral wrap: assemble chelator donors? use frag on all sites is wrong;
            # approximate saturated chelate as metal + chelator anion placed (single fragment, 6-dentate) -> use free-lig proxy
            chsy,chpo,chq=_chelate_cluster(sym,mq,chelator)
            chelc=_three_alpb(chsy,chpo,chq,sym,spin_mode,steps=max(steps,40))
            def Ef(clab, nlk, nsolv_keys, nch):
                e=mean(clab)-mean(eM)-nlk*mean(eLk)-nch*mean(eCh)
                for k in set(nsolv_keys): e-=nsolv_keys.count(k)*mean(eS[k])
                return e
            Ef_ion=Ef(ionc,0,keys,0); Ef_node=Ef(nodec,2,keys4,0); Ef_ch=Ef(chelc,0,[],1)
            rows.append(dict(metal=sym,system=sysk,Ef_ion=round(Ef_ion,2),Ef_node=round(Ef_node,2),
                             Ef_chelator=round(Ef_ch,2),nucleation_dG=round(Ef_node-Ef_ion,2),
                             seq_margin=round(Ef_ch-Ef_node,2)))
            print(f'{sym} {sysk:20s} nucleation={Ef_node-Ef_ion:6.2f}  seq_margin={Ef_ch-Ef_node:6.2f}',flush=True)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- 4. visualize
def list_models(metal=None, scenario=None, contains=None):
    labs=list(MODELS)
    if metal: labs=[l for l in labs if MODELS[l]['metal']==metal]
    if scenario: labs=[l for l in labs if MODELS[l]['scenario']==scenario]
    if contains: labs=[l for l in labs if contains in l]
    return labs

def view_model(label, width=560, height=420, spin=False):
    """Interactive 3D of a generated model (needs py3Dmol)."""
    from ebu_viz import view_xyz
    info=MODELS.get(label)
    if info is None: raise KeyError(f'No model {label!r}. Use list_models().')
    return view_xyz(info['xyz_path'], metal_symbols=info['metal'], width=width, height=height, spin=spin)
