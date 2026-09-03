"""selftest.py -- fast end-to-end check that the toolkit works on this machine.
Runs enumeration (rdkit-only) then a couple of single-point GFN2-xTB energies
(tblite). Prints PASS/FAIL for each stage. Takes ~10-20 s."""
import sys, warnings
warnings.filterwarnings('ignore')

def _ok(msg): print(f'  \033[32mPASS\033[0m  {msg}')
def _bad(msg): print(f'  \033[31mFAIL\033[0m  {msg}')

print('EBU screening tool — self test\n' + '-'*40)

# 1. imports
try:
    import ebu_toolkit as T
    _ok('import ebu_toolkit (+ ebu_core, geometry_qc, donor_perception)')
except Exception as e:
    _bad(f'import failed: {e}'); sys.exit(1)

# 2. enumeration (geometry engine, rdkit only)
try:
    df = T.enumerate_models(metals=['Ni', ('Fe', 3)], linkers=['BTC'],
                            additives=['EDDA'], additive_mode='subsets',
                            max_copies=2, verbose=False)
    assert len(df) > 0 and 'chelation_extent' in df.columns
    nodes = sorted({f"{l.count('BTC[')}BTC" for l in df[df.scenario=='linker_only'].label})
    _ok(f'enumeration: {len(df)} models; linker-only nodes = {nodes or "(need max_copies>=2)"}')
except Exception as e:
    _bad(f'enumeration failed: {e}'); sys.exit(1)

# 3. energies (tblite / GFN2-xTB) — single point, a few models
try:
    ranked = T.formation_energies(df=df.head(4), solvent=None, spin_mode='highspin',
                                  relax=False, verbose=False)
    vals = ranked['E_form_eV'].dropna()
    assert len(vals) > 0
    _ok(f'GFN2-xTB energies: {len(vals)} formation energies, '
        f'range {vals.min():.1f}..{vals.max():.1f} eV')
except Exception as e:
    _bad(f'energy backend (tblite) failed: {e}')
    print('       -> install with: pip install tblite   (Linux) ')
    sys.exit(1)

# 4. viewer availability (optional)
try:
    import py3Dmol  # noqa
    _ok('py3Dmol available (3D viewer will work)')
except ImportError:
    print('  note  py3Dmol not installed — 3D viewer disabled (pip install py3Dmol)')

print('-'*40 + '\nALL CORE CHECKS PASSED — open ebu_tool.ipynb to start.')
