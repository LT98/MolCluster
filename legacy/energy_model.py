"""
energy_model.py -- ML interatomic potential (MLIP) geometry optimization and
energy estimation, replacing the removed UFF relaxation step in ebu_core.py.

Why this exists
----------------
GeometryPlacer.place() used to finish with a UFF relaxation pass, but UFF is
an organic-only force field: it had no real parameters for a transition
metal's coordination sphere, so ebu_core.py disguised the metal atom as a
neutral carbon just to get UFF to accept it. That trick silently broke (a
caught AtomValenceException, since carbon's standard valence caps at 4) for
any coordination number >= 5 and for any anionic (deprotonated) donor once
its dative bond was converted to a real bond -- i.e. for nearly every
candidate this pipeline actually cares about. Removing UFF and doing the
analytic placement correctly (see ebu_core._donor_placement_frame) fixed the
geometry's local bond ANGLES, but there is still no step that actually
relaxes/polishes the full 3D structure or estimates its energy.

A universal ML interatomic potential is the right replacement: unlike UFF it
has real parameters for transition metals (trained on relaxed DFT structures
across the periodic table, including coordination compounds and open-shell
metals), so the metal can be optimized as itself -- no disguise hack needed.
This module wraps MACE-MP-0 (github.com/ACEsuit/mace-foundations), a
foundation MLIP trained on the Materials Project's relaxed DFT trajectories,
via ASE's Calculator interface (this is also what GeometryPlacer.to_ase()'s
docstring was already anticipating). It is intentionally NOT wired into
GeometryPlacer.place() or enumerate_sbus() automatically -- per the decision
to keep the default pipeline geometry-only, this is an explicit, opt-in step
you run on candidates that already passed geometry_qc.

What this gives you (mapped onto the original 4-step pipeline spec)
---------------------------------------------------------------------
  - optimize_geometry(): fast local relaxation of a placed EBU (step 1 polish)
  - single_point_energy(): one energy evaluation, no relaxation
  - estimate_formation_energy(): EBU energy minus isolated metal-ion and
    free-ligand reference energies -- an approximate formation energy
    (step 2 of the original spec). "Approximate" is not a hedge to skip:
    MACE-MP-0 is charge/spin-agnostic (it only sees atomic species and
    positions, not the formal charges RDKit tracks for our own placement
    bookkeeping), so this number is best used for RELATIVE ranking across
    candidate ligand/metal/geometry combinations, not as an absolute,
    publication-grade formation energy. Getting the latter would need
    matched protonation states, explicit counterions, and a real DFT
    cross-check -- a natural next step once the ML screening layer here
    narrows candidates down.

Environment note
-----------------
This module lazy-imports torch/ase/mace so the rest of the pipeline (which
does not need them) keeps working if they are not installed. In the sandbox
this was developed in, installing them was not practical: default PyPI
`torch` pulls a full CUDA/cuda-toolkit dependency stack (multiple GB, no GPU
present here to use it anyway), and the legitimate CPU-only wheel index
(download.pytorch.org/whl/cpu) is blocked by the sandbox's network proxy.
Install with, e.g.:

    pip install torch --index-url https://download.pytorch.org/whl/cpu
    pip install mace-torch ase

on a machine where that index is reachable (this is normal on a typical
laptop/workstation/HPC login node -- the block is specific to this sandbox).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem

MACE_MP_MODEL_SIZES = ('small', 'medium', 'large')
DEFAULT_MODEL_SIZE = 'medium'


def _require_mlip_stack():
    missing = []
    for mod_name in ('torch', 'ase', 'mace'):
        try:
            __import__(mod_name)
        except ImportError:
            missing.append(mod_name)
    if missing:
        raise ImportError(
            f"energy_model requires {missing} which "
            f"{'is' if len(missing) == 1 else 'are'} not installed here. "
            f"Install with (on a machine that can reach the real PyTorch CPU "
            f"wheel index):\n"
            f"  pip install torch --index-url https://download.pytorch.org/whl/cpu\n"
            f"  pip install mace-torch ase\n"
            f"See the module docstring for why this specific sandbox couldn't "
            f"install them directly."
        )


@dataclass
class OptResult:
    mol3d: 'Chem.Mol'
    energy_eV: float
    converged: bool
    n_steps: int
    fmax_reached: float


@dataclass
class FormationEnergyResult:
    ebu_energy_eV: float
    metal_ref_energy_eV: float
    ligand_ref_energies_eV: List[float] = field(default_factory=list)
    formation_energy_eV: float = 0.0
    notes: str = ''


class MLEnergyModel:
    """
    Thin wrapper around a MACE-MP-0 ASE calculator. Instantiating this class
    is cheap on its own (no imports happen until you actually call a method
    that needs the model), so it's safe to construct one at module import
    time in a notebook and only pay the model-loading cost when first used.
    """

    def __init__(self, model_size: str = DEFAULT_MODEL_SIZE, device: str = 'cpu'):
        if model_size not in MACE_MP_MODEL_SIZES:
            raise ValueError(f'model_size must be one of {MACE_MP_MODEL_SIZES}, got {model_size!r}')
        self.model_size = model_size
        self.device = device
        self._calc = None   # lazily constructed on first use

    def _get_calculator(self):
        if self._calc is not None:
            return self._calc
        _require_mlip_stack()
        from mace.calculators import mace_mp
        self._calc = mace_mp(model=self.model_size, device=self.device, default_dtype='float64')
        return self._calc

    # -- conversion helpers --------------------------------------------------

    @staticmethod
    def mol3d_to_ase(mol3d: 'Chem.Mol'):
        """RDKit Mol (with a 3D conformer) -> ASE Atoms. Formal charges are
        NOT carried over -- MACE-MP-0 takes species + positions only."""
        _require_mlip_stack()
        import ase
        if mol3d.GetNumConformers() == 0:
            raise ValueError('mol3d has no 3D conformer -- run GeometryPlacer.place() first.')
        conf = mol3d.GetConformer()
        syms = [a.GetSymbol() for a in mol3d.GetAtoms()]
        pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol3d.GetNumAtoms())])
        atoms = ase.Atoms(symbols=syms, positions=pos)
        atoms.pbc = False
        return atoms

    @staticmethod
    def _ase_to_mol3d(atoms, template_mol3d: 'Chem.Mol') -> 'Chem.Mol':
        """Copy optimized ASE positions back onto a COPY of template_mol3d
        (same atom order/connectivity, only the conformer changes)."""
        from rdkit.Chem import RWMol
        out = RWMol(template_mol3d)
        conf = out.GetConformer()
        pos = atoms.get_positions()
        for i in range(out.GetNumAtoms()):
            conf.SetAtomPosition(i, pos[i].tolist())
        return out.GetMol()

    # -- core capabilities ----------------------------------------------------

    def single_point_energy(self, mol3d: 'Chem.Mol') -> float:
        """One energy evaluation (eV) at the given geometry, no relaxation."""
        atoms = self.mol3d_to_ase(mol3d)
        atoms.calc = self._get_calculator()
        return float(atoms.get_potential_energy())

    def optimize_geometry(
        self,
        mol3d       : 'Chem.Mol',
        fmax        : float = 0.05,
        max_steps   : int   = 200,
        fixed_atoms : Optional[List[int]] = None,
        optimizer   : str   = 'FIRE',
        logfile     : Optional[str] = None,
    ) -> OptResult:
        """
        Relax an EBU's 3D geometry with a real ML potential -- unlike the
        removed UFF step, the metal is optimized AS the metal (no carbon
        disguise), so coordination numbers above 4 are not a problem.

        fixed_atoms: optional atom indices to freeze (e.g. the metal index,
        if you want to check ligand-only relaxation around a fixed center).
        Leave None to let everything relax, which is usually what you want.
        """
        _require_mlip_stack()
        from ase.optimize import FIRE, BFGS, LBFGS
        from ase.constraints import FixAtoms

        atoms = self.mol3d_to_ase(mol3d)
        atoms.calc = self._get_calculator()
        if fixed_atoms:
            atoms.set_constraint(FixAtoms(indices=list(fixed_atoms)))

        opt_cls = {'FIRE': FIRE, 'BFGS': BFGS, 'LBFGS': LBFGS}.get(optimizer)
        if opt_cls is None:
            raise ValueError(f"optimizer must be one of 'FIRE', 'BFGS', 'LBFGS', got {optimizer!r}")

        dyn = opt_cls(atoms, logfile=logfile)
        converged = dyn.run(fmax=fmax, steps=max_steps)
        n_steps = dyn.get_number_of_steps()
        fmax_reached = float(np.sqrt((atoms.get_forces() ** 2).sum(axis=1).max()))
        energy = float(atoms.get_potential_energy())

        optimized_mol3d = self._ase_to_mol3d(atoms, mol3d)
        return OptResult(
            mol3d=optimized_mol3d,
            energy_eV=energy,
            converged=bool(converged),
            n_steps=n_steps,
            fmax_reached=fmax_reached,
        )

    # -- formation energy (approximate, for ranking candidates) --------------

    def estimate_formation_energy(
        self,
        mol3d        : 'Chem.Mol',
        metal_symbol : str,
        ligands      : List,   # List[LigandBuilder] -- typed loosely to avoid a hard ebu_core dependency
        relax_ebu    : bool = True,
        relax_refs   : bool = True,
    ) -> FormationEnergyResult:
        """
        Approximate formation energy = E(assembled EBU) - E(isolated metal
        atom) - sum(E(each free ligand, its own relaxed geometry)).

        This is a screening-level number, not a thermodynamically rigorous
        formation energy: it compares a formally-charged EBU (per RDKit's
        own bookkeeping) against neutral isolated fragments evaluated by a
        charge-agnostic MLIP, ignores solvent/counterion contributions
        entirely (solvation is a separate, later step in the original
        pipeline spec), and the ligand reference geometry is whatever
        LigandBuilder's own free-conformer embedding produced (already
        computed once at LigandBuilder construction time, not necessarily
        MLIP-relaxed unless relax_refs=True). Use it to RANK candidates
        against each other, not as an absolute number.
        """
        _require_mlip_stack()
        import ase

        if relax_ebu:
            ebu_energy = self.optimize_geometry(mol3d).energy_eV
        else:
            ebu_energy = self.single_point_energy(mol3d)

        metal_atoms = ase.Atoms(symbols=[metal_symbol], positions=[[0., 0., 0.]])
        metal_atoms.calc = self._get_calculator()
        metal_ref_energy = float(metal_atoms.get_potential_energy())

        ligand_ref_energies = []
        for lig in ligands:
            if relax_refs:
                e = self.optimize_geometry(lig.mol).energy_eV
            else:
                e = self.single_point_energy(lig.mol)
            ligand_ref_energies.append(e)

        formation_energy = ebu_energy - metal_ref_energy - sum(ligand_ref_energies)

        return FormationEnergyResult(
            ebu_energy_eV=ebu_energy,
            metal_ref_energy_eV=metal_ref_energy,
            ligand_ref_energies_eV=ligand_ref_energies,
            formation_energy_eV=formation_energy,
            notes=('Screening-level estimate: charge/spin-agnostic MLIP energies, '
                   'no solvent/counterion terms. Use for relative ranking.'),
        )
