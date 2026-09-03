import nbformat as nbf
nb = nbf.v4.new_notebook()
C=[]
def md(s): C.append(nbf.v4.new_markdown_cell(s))
def code(s): C.append(nbf.v4.new_code_cell(s))

md("""# EBU Screening Tool — metal centres, linkers, additives, solvation & 3D viewer

A configurable front-end to the `ebu_toolkit` pipeline. Edit the **CONFIG** cells, run top-to-bottom, and inspect / visualize the generated extended building units (EBUs).

**Pipeline layers (all GFN2-xTB, screening-grade):**
1. **Enumerate** candidate EBUs (geometry + QC) for your metals × linkers × additives.
2. **Rank** them by formation energy — gas phase or implicit solvent.
3. **Solvation sweep** — explicit-shell + ALPB continuum over water/ethanol/DMF combinations, giving nucleation driving force and additive-sequestration margin.
4. **Visualize** any generated model in 3D.

**Requirements:** `rdkit`, `ase`, `tblite` (GFN2-xTB), `py3Dmol` (3D), `pandas`, `matplotlib`.
```
pip install rdkit ase tblite py3Dmol pandas matplotlib
```
""")

md("## 0 · Setup")
code("""import warnings; warnings.filterwarnings('ignore')
import pandas as pd, matplotlib.pyplot as plt
import ebu_toolkit as T
from ebu_toolkit import LIGAND_LIBRARY, SOLVENTS
pd.set_option('display.max_colwidth', 60); pd.set_option('display.width', 160)
print('Built-in ligand library :', ', '.join(LIGAND_LIBRARY))
print('Built-in solvents        :', ', '.join(SOLVENTS))
print('Add your own any time, e.g.:  LIGAND_LIBRARY[\\'NDC\\'] = \\'OC(=O)c1ccc2cc(C(=O)O)ccc2c1\\'')""")

md("""## 1 · CONFIG — what to build  ✏️

Edit this cell. Metals may be plain symbols (default oxidation state) or `(symbol, charge)` tuples.
Linkers/additives are names from the library above, or `'Name:SMILES'` to define inline.

`additive_mode`: `'subsets'` builds linker-only **plus every additive combination** (best for competition studies) · `'all'` uses all additives together · `'none'` ignores additives.""")
code("""METALS       = ['Ni', ('Fe', 3)]      # e.g. 'Cu', ('Co',2), 'Zr' ...
LINKERS      = ['BTC']                 # framework linker(s)
ADDITIVES    = ['EDTA', 'EDDA']        # chelating / modulating additives
ADDITIVE_MODE = 'subsets'              # 'subsets' | 'all' | 'none'

MAX_COPIES   = 2      # max copies of each ligand species per metal
CLASH_SCALE  = 0.65   # QC clash tolerance (lower = more permissive; Fe often needs 0.60)
SPIN_MODE    = 'highspin'   # 'highspin' (d-electron based) | 'auto' (minimal spin)""")

md("## 2 · Enumerate candidate EBUs")
code("""models_df = T.enumerate_models(
    metals=METALS, linkers=LINKERS, additives=ADDITIVES,
    additive_mode=ADDITIVE_MODE, max_copies=MAX_COPIES,
    clash_scale=CLASH_SCALE, verbose=True)

print(f'\\n{len(models_df)} QC-valid models. Counts by metal × scenario:')
display(models_df.groupby(['metal','scenario']).size().rename('n_models').reset_index())
display(models_df.head(20))""")

md("""## 3 · Rank by formation energy

`solvent=None` → gas phase.  `solvent='water'|'ethanol'|'dmf'` → implicit ALPB continuum.
Lower `E_form` = more stable (thermodynamic sink); `E_form_per_bond` exposes bond *quality* vs bond *count*.""")
code("""SOLVENT_FOR_RANKING = None    # None (gas) | 'water' | 'ethanol' | 'dmf'

ranked = T.formation_energies(solvent=SOLVENT_FOR_RANKING, spin_mode=SPIN_MODE,
                              relax=True, steps=60, verbose=True)
display(ranked.head(25))""")

code("""# Landscape plot: formation energy per scenario, per metal (framework node dashed)
fig, axes = plt.subplots(1, ranked['metal'].nunique(), figsize=(6*ranked['metal'].nunique(),5), squeeze=False)
for ax,(metal,g) in zip(axes[0], ranked.groupby('metal')):
    node = g[g['scenario']=='linker_only']['E_form_eV'].min()
    for i,(scen,gs) in enumerate(g.groupby('scenario')):
        ax.scatter([i]*len(gs), gs['E_form_eV'], s=60, edgecolor='k', alpha=.8)
        ax.annotate(scen, (i,ax.get_ylim()[0]), rotation=90, fontsize=7, va='bottom')
    if node==node: ax.axhline(node, ls='--', color='green', lw=1)
    ax.set_title(metal); ax.set_ylabel('E_form (eV)'); ax.set_xticks([]); ax.invert_yaxis()
plt.tight_layout(); plt.show()""")

md("""## 4 · Solvation sweep (explicit shell + ALPB continuum)

Completes each metal's octahedral sphere with explicit solvent (equal-ratio combinations of water/ethanol/DMF), wraps the bulk in ALPB, and reports two decision metrics per solvent system:

* **nucleation_dG** = E_form(node) − E_form(solvated ion): *negative → framework beats the dissolved ion* (metal leaves solution to build framework).
* **seq_margin** = E_form(chelator complex) − E_form(node): *negative → the additive sequesters the metal* (competes with the framework).""")
code("""SOLV_METALS   = METALS
SOLV_LINKER   = 'BTC'
SOLV_CHELATOR = 'EDTA'
SOLV_SYSTEMS  = ['water','ethanol','dmf','water+ethanol','water+dmf','ethanol+dmf','water+ethanol+dmf']

solv_df = T.solvation_sweep(SOLV_METALS, linker=SOLV_LINKER, chelator=SOLV_CHELATOR,
                            systems=SOLV_SYSTEMS, spin_mode=SPIN_MODE, steps=25)
display(solv_df)""")

code("""# Two-metric bar chart across solvent systems
fig,(a1,a2)=plt.subplots(2,1,figsize=(11,8),sharex=True)
import numpy as np
mets=solv_df['metal'].unique(); x=np.arange(len(SOLV_SYSTEMS)); w=.8/len(mets)
for j,mt in enumerate(mets):
    g=solv_df[solv_df.metal==mt].set_index('system').reindex(SOLV_SYSTEMS)
    a1.bar(x+j*w, g['nucleation_dG'], w, label=mt)
    a2.bar(x+j*w, g['seq_margin'],   w, label=mt)
for a,t in [(a1,'nucleation driving force (more negative = framework favoured)'),
            (a2,'sequestration margin (below 0 = additive traps metal)')]:
    a.axhline(0,color='k',lw=.8); a.set_ylabel(t); a.legend(); a.grid(axis='y',alpha=.25)
a2.set_xticks(x+w*(len(mets)-1)/2); a2.set_xticklabels(SOLV_SYSTEMS,rotation=30,ha='right')
plt.tight_layout(); plt.show()""")

md("""## 5 · 3D viewer  🔬

Pick any generated model to render. Run the cell; use the dropdowns (needs `ipywidgets`) or call `T.view_model('<label>')` directly.
`T.list_models(metal=..., scenario=..., contains=...)` filters the available labels.""")
code("""labels = T.list_models()
print(len(labels), 'models available. Examples:'); [print('  ', l) for l in labels[:8]]

try:
    import ipywidgets as W
    from IPython.display import display, clear_output
    metals = sorted({T.MODELS[l]['metal'] for l in labels})
    dd_metal = W.Dropdown(options=['(all)']+metals, description='Metal:')
    dd_model = W.Dropdown(options=labels, description='Model:', layout=W.Layout(width='70%'))
    out = W.Output()
    def _refresh(*_):
        m=dd_metal.value
        dd_model.options=[l for l in labels if m=='(all)' or T.MODELS[l]['metal']==m]
    def _show(*_):
        with out:
            clear_output(wait=True)
            if dd_model.value: T.view_model(dd_model.value)
    dd_metal.observe(_refresh,'value'); dd_model.observe(_show,'value')
    display(W.VBox([dd_metal, dd_model, out])); _show()
except ImportError:
    print('\\nipywidgets not installed — call directly, e.g.:')
    print(\"   T.view_model('%s')\" % (labels[0] if labels else '<label>'))""")

code("""# Direct call (works without ipywidgets):
T.view_model(T.list_models()[0], spin=True)""")

md("""---
### Notes
* All energies are **screening-grade GFN2-xTB** — trust *relative* rankings and the *direction/magnitude* of solvent effects, not absolute values. See the report `.md` files for full caveats.
* Extend the tool by adding entries to `LIGAND_LIBRARY` / `SOLVENTS`, or new metals to `METALS`. High-spin multiplicities for uncommon ions can be added to `ebu_toolkit._D_ELECTRONS`.
* Generated `.xyz` files are written under `tool_outputs/<metal>/` and can be opened in any molecular viewer.""")

nb['cells']=C
nbf.write(nb,'ebu_tool.ipynb')
print('wrote ebu_tool.ipynb with',len(C),'cells')
