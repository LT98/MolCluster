"""run_metal.py -- enumerate M(q) + BTC (+EDTA/+EDDA) EBUs. Metal via env:
   METAL_SYMBOL (default Ni), METAL_CHARGE (default 2), OUT_TAG (default metal).
Geometry-only here; energies handled by xtb_energy_metal.py."""
import os, warnings
warnings.filterwarnings('ignore')
from ebu_core import LigandBuilder, EBUBuilder, autonomous_enumerate_sbus

SYM = os.environ.get('METAL_SYMBOL', 'Ni')
CHG = int(os.environ.get('METAL_CHARGE', '2'))
TAG = os.environ.get('OUT_TAG', SYM.lower())
OUT = os.path.join(os.getcwd(), f'{TAG}_btc_outputs')
os.makedirs(OUT, exist_ok=True)

btc  = LigandBuilder('OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1', name='BTC')
edta = LigandBuilder('C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O', name='EDTA')
edda = LigandBuilder('OC(=O)CNCCNCC(=O)O', name='EDDA')
M = EBUBuilder(SYM, CHG)

scenarios = {'A_BTC_only':[btc], 'B_BTC_EDTA':[btc,edta],
             'C_BTC_EDDA':[btc,edda], 'D_BTC_EDTA_EDDA':[btc,edta,edda]}
total=0
for tag, sp in scenarios.items():
    print('='*60, f'\n{SYM}{CHG}+ {tag}')
    res = autonomous_enumerate_sbus(builder=M, ligand_species=sp, out_dir=OUT,
            label_prefix=f'{SYM}{CHG}+_{tag}', max_copies_per_species=2, clash_scale=float(os.environ.get('CLASH_SCALE','0.65')),
            validate=True, discard_invalid=True)
    total += len(res)
print(f'\nTOTAL kept: {total} -> {OUT}')
