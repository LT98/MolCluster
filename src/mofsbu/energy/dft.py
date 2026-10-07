"""DFT through GPU4PySCF — the `Fidelity.DFT` rung (WORKPLAN_energy E4, stage 5b/5c).

In-process Kohn–Sham on the declared GPU.  Restricted for a singlet, unrestricted for
anything else; density fitting on; a continuum is a `model:solvent` medium token
(`smd:water`, `pcm:water`), so the `methods` row says which continuum it was (C16).

What this backend refuses rather than approximates:

* **No GPU, no DFT.**  `MOFSBU_DEVICE=cpu` makes it unavailable, never a CPU fallback at
  another cost or another code (invariant 3: declared, never detected).
* **A broken CUDA toolchain is a missing install.**  `available()` compiles one kernel, so
  a CuPy whose NVRTC and headers disagree is reported here rather than as an exception
  in the middle of a run.
* **A geometry the electrons disagree with.**  `check_oxidation_state` compares the
  metal's Mulliken spin population with the d-count its graph claims and raises
  `OxidationStateMismatch` (`oxidation_state_mismatch`) — THQ is redox non-innocent.

Hessians are out of scope (§3a): stage 5d uses MACE Hessians.
"""
from __future__ import annotations

import math
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from mofsbu._types import EnergyBackendUnavailable, Fidelity, MethodSpec, MofsbuError, split_medium

HARTREE_EV = 27.211386245988
BOHR_ANGSTROM = 0.529177210903

#: Solvents a medium may name, with the static dielectric PCM needs (SMD reads its own
#: descriptor table, which has all of these).
DIELECTRIC = {"water": 78.3553, "methanol": 32.613, "ethanol": 24.852, "dmso": 46.826,
              "acetonitrile": 35.688, "dmf": 37.219, "thf": 7.4257}
CONTINUA = ("smd", "pcm")

#: Spin population on a metal further than this from the d-count's unpaired electrons is
#: a different oxidation or spin state from the one the graph records.
SPIN_POPULATION_TOL = 0.5


class OxidationStateMismatch(MofsbuError):
    """The SCF put a different number of unpaired electrons on the metal than L0 claims."""

    code = "oxidation_state_mismatch"


_TOOLCHAIN: dict[str, str | None] = {}


def _toolchain_problem() -> str | None:
    """None if GPU4PySCF imports and CuPy can compile a kernel here; else why not.  Cached."""
    if "problem" in _TOOLCHAIN:
        return _TOOLCHAIN["problem"]
    try:
        import cupy
        import gpu4pyscf  # noqa: F401

        if cupy.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError("no CUDA device visible")
        k = cupy.ElementwiseKernel("float64 x", "float64 y", "y = x * 2", "mofsbu_probe")
        k(cupy.arange(2, dtype=cupy.float64))
        # A reduction goes through CUB, which includes the toolkit's fp8 headers: the part
        # that breaks when CuPy's runtime and the installed headers differ in version.
        (cupy.arange(64) > 3).sum().get()
        problem = None
    except ImportError as exc:
        problem = (f"{exc.name or exc} is not importable: pip install gpu4pyscf-cuda13x "
                   "with a cupy-cuda13x built for the CUDA runtime torch pins")
    except Exception as exc:                                            # noqa: BLE001
        first = str(exc).strip().splitlines()[0][:200] if str(exc).strip() else ""
        problem = (f"CuPy cannot compile a kernel ({type(exc).__name__}: {first}); its "
                   "NVRTC and headers must match its CUDA major — check `cupy.show_config()`")
    _TOOLCHAIN["problem"] = problem
    return problem


class DFTBackend:
    """Kohn–Sham DFT on GPU4PySCF behind the `EnergyBackend` protocol.

    `xc` and `basis` are constructor arguments and part of `method`, so two recipes never
    share a `methods` row.  The default is ωB97M-V/def2-TZVPD, the OMol25 level MACE-OMOL-0
    was trained on, which makes E_DFT − E_MACE a per-species MACE error (5b).
    """

    name = "dft"
    fidelity = Fidelity.DFT
    charge_aware = True
    spin_aware = True
    code = "gpu4pyscf"

    def __init__(self, *, xc: str = "wb97m-v", basis: str = "def2-tzvpd",
                 density_fit: bool = True, grids_level: int = 3, nlc_grids_level: int = 1,
                 conv_tol: float = 1e-9, max_cycle: int = 200, guess_xc: str | None = "pbe",
                 device: str | None = None) -> None:
        from mofsbu.config import compute_device

        self.xc = xc.lower()
        self.basis = basis.lower()
        self.density_fit = density_fit
        self.grids_level = grids_level
        self.nlc_grids_level = nlc_grids_level
        self.conv_tol = conv_tol
        self.max_cycle = max_cycle
        self.guess_xc = None if guess_xc is None else guess_xc.lower()
        self.device = compute_device() if device is None else device
        self.method = f"{self.xc}/{self.basis}"

    # -- availability ---------------------------------------------------

    def available(self) -> bool:
        return self._unavailable_reason() is None

    def _unavailable_reason(self) -> str | None:
        if not self.device.startswith("cuda"):
            return (f"GPU4PySCF runs on a CUDA device and MOFSBU_DEVICE={self.device!r}; "
                    "declare MOFSBU_DEVICE=cuda")
        return _toolchain_problem()

    def install_hint(self) -> str:
        return self._unavailable_reason() or "installed"

    def _require(self) -> None:
        why = self._unavailable_reason()
        if why is not None:
            raise EnergyBackendUnavailable(f"{self.name} backend: {why}")

    def code_version(self) -> str:
        from mofsbu.energy.backends import _installed_version

        return (f"{_installed_version('gpu4pyscf-cuda12x', 'gpu4pyscf-cuda13x', 'gpu4pyscf')}"
                f"/pyscf-{_installed_version('pyscf')}")

    # -- the method row -------------------------------------------------

    def _medium(self, solvent: str | None) -> tuple[str, str] | None:
        if solvent is None:
            return None
        model, name = split_medium(solvent)
        if model not in CONTINUA:
            raise ValueError(f"continuum {model!r} is not wired here; have {list(CONTINUA)}")
        if name not in DIELECTRIC:
            raise ValueError(f"unknown solvent {name!r}; have {sorted(DIELECTRIC)}")
        return model, name

    def method_spec(self, *, charge: int, multiplicity: int,
                    solvent: str | None = None) -> MethodSpec:
        from mofsbu.versions import ALGO_VERSIONS

        medium = self._medium(solvent)
        extras: dict[str, Any] = {
            "algo": ALGO_VERSIONS["energy_backends"], "xc": self.xc, "basis": self.basis,
            "density_fit": self.density_fit, "grids_level": self.grids_level,
            "nlc_grids_level": self.nlc_grids_level, "guess": self.guess_xc or "minao",
            "reference": "RKS" if multiplicity == 1 else "UKS"}
        return MethodSpec(code=self.code, code_version=self.code_version(),
                          method=self.method,
                          solvent=None if medium is None else f"{medium[0]}:{medium[1]}",
                          charge=charge, multiplicity=multiplicity, extras=extras)

    # -- running it -----------------------------------------------------

    def _mol(self, symbols: Sequence[str], positions: Any, charge: int, multiplicity: int):
        from pyscf import gto

        atoms = [(s, tuple(float(c) for c in xyz))
                 for s, xyz in zip(symbols, np.asarray(positions, dtype=float), strict=True)]
        return gto.M(atom=atoms, unit="Angstrom", basis=self.basis, charge=int(charge),
                     spin=int(multiplicity) - 1, verbose=0)

    def _mf(self, mol, multiplicity: int, solvent: str | None):
        from gpu4pyscf import dft

        mf = (dft.RKS if multiplicity == 1 else dft.UKS)(mol, xc=self.xc)
        if self.density_fit:
            mf = mf.density_fit()
        mf.grids.level = self.grids_level
        # VV10 is quadratic in its grid; level 1 moves a Ni complex's energy by 1 µEh and
        # costs a fifth of level 3 (WORKPLAN_energy §3a).
        mf.nlcgrids.level = self.nlc_grids_level
        mf.conv_tol = self.conv_tol
        mf.max_cycle = self.max_cycle
        medium = self._medium(solvent)
        if medium is not None:
            model, name = medium
            if model == "smd":
                mf = mf.SMD()
                mf.with_solvent.solvent = name
            else:
                mf = mf.PCM()
                mf.with_solvent.method = "IEF-PCM"
                mf.with_solvent.eps = DIELECTRIC[name]
        return mf

    def _prepare(self, symbols, positions, charge, multiplicity, solvent):
        from mofsbu.energy.backends import check_spin

        check_spin(list(symbols), charge, multiplicity)
        self._medium(solvent)
        self._require()
        mol = self._mol(symbols, positions, charge, multiplicity)
        return mol, self._mf(mol, multiplicity, solvent)

    def _guess(self, mol, multiplicity: int):
        """Starting density from a converged gas-phase SCF at `guess_xc`, or None (minao).

        From the minao guess a Ni(II) triplet at wB97M-V oscillates and does not converge in
        150 cycles; started from a PBE density it converges in ~30 (WORKPLAN_energy §3a).
        """
        if self.guess_xc is None or self.guess_xc == self.xc:
            return None
        from gpu4pyscf import dft

        pre = (dft.RKS if multiplicity == 1 else dft.UKS)(mol, xc=self.guess_xc)
        if self.density_fit:
            pre = pre.density_fit()
        pre.conv_tol = 1e-6
        pre.max_cycle = self.max_cycle
        pre.kernel()
        return pre.make_rdm1()

    def single_point(self, symbols, positions, *, charge, multiplicity,
                     solvent=None):
        from mofsbu.energy.backends import EnergyResult

        mol, mf = self._prepare(symbols, positions, charge, multiplicity, solvent)
        t0 = time.perf_counter()
        e_tot = float(mf.kernel(dm0=self._guess(mol, multiplicity)))
        extras = _analysis(mf, multiplicity)
        extras["seconds"] = round(time.perf_counter() - t0, 2)
        return EnergyResult(energy=e_tot * HARTREE_EV,
                            method=self.method_spec(charge=charge, multiplicity=multiplicity,
                                                    solvent=solvent),
                            fidelity=self.fidelity, converged=bool(mf.converged),
                            extras=extras)

    def relax(self, symbols, positions, *, charge, multiplicity, solvent=None,
              fmax=0.05, steps=100):
        """geomeTRIC on the SCF gradient.  `fmax` (eV/Å) sets the gradient criteria."""
        from pyscf.geomopt.geometric_solver import kernel as geometric_kernel

        from mofsbu.energy.backends import RelaxResult

        mol, mf = self._prepare(symbols, positions, charge, multiplicity, solvent)
        mf.kernel(dm0=self._guess(mol, multiplicity))  # each later step starts from the last
        gmax = fmax / HARTREE_EV * BOHR_ANGSTROM        # eV/Å -> Eh/Bohr
        trail: list[tuple[float, float]] = []

        def callback(env):
            g = np.asarray(env["gradients"])
            trail.append((float(env["energy"]),
                          float(np.sqrt((g ** 2).sum(axis=1)).max()) * HARTREE_EV / BOHR_ANGSTROM))

        converged, final = geometric_kernel(
            mf, maxsteps=steps, callback=callback, assert_convergence=False,
            convergence_gmax=gmax, convergence_grms=gmax * 2 / 3)
        # The energy returned belongs to the returned coordinates, so it is recomputed
        # there rather than read off the optimiser's last trial step.
        sp = self.single_point([final.atom_pure_symbol(i) for i in range(final.natm)],
                               final.atom_coords(unit="Angstrom"), charge=charge,
                               multiplicity=multiplicity, solvent=solvent)
        return RelaxResult(
            symbols=tuple(symbols), positions=final.atom_coords(unit="Angstrom"),
            energy=sp.energy, method=sp.method, fidelity=self.fidelity,
            converged=bool(converged) and sp.converged, n_steps=len(trail),
            fmax=trail[-1][1] if trail else math.nan,
            initial_energy=trail[0][0] * HARTREE_EV if trail else None)


def _analysis(mf, multiplicity: int) -> dict[str, Any]:
    """Mulliken charges and spin populations per atom, and <S²> for an open shell."""
    mol = mf.mol
    dm = mf.make_rdm1()
    dm = np.asarray(dm.get() if hasattr(dm, "get") else dm)
    s = mol.intor_symmetric("int1e_ovlp")
    if dm.ndim == 2:
        dm_tot, dm_spin = dm, None
    else:
        dm_tot, dm_spin = dm[0] + dm[1], dm[0] - dm[1]
    slices = mol.aoslice_by_atom()[:, 2:4]

    def per_atom(d):
        pop = np.einsum("ij,ji->i", d, s)
        return [float(pop[a:b].sum()) for a, b in slices]

    out: dict[str, Any] = {
        "mulliken_charges": [round(float(z) - p, 4)
                             for z, p in zip(mol.atom_charges(), per_atom(dm_tot), strict=True)],
        "scf_cycles": getattr(mf, "cycles", None)}
    if dm_spin is not None:
        out["spin_populations"] = [round(p, 4) for p in per_atom(dm_spin)]
        s2, _ = mf.spin_square()
        out["s2"] = round(float(s2), 4)
        out["s2_ideal"] = (multiplicity - 1) / 2 * ((multiplicity - 1) / 2 + 1)
    return out


def expected_unpaired(symbol: str, oxidation_state: int, spin_class: str) -> int:
    """Unpaired d electrons the graph claims for one metal centre."""
    from mofsbu.energy.backends import spin_class_multiplicity

    return spin_class_multiplicity(symbol, oxidation_state, spin_class) - 1


def check_oxidation_state(result: Any, metal_index: int, expected: int, *,
                          tol: float = SPIN_POPULATION_TOL) -> float:
    """Return the metal's spin population, or raise `OxidationStateMismatch`.

    A closed-shell result has no spin density, so `expected` must then be zero.
    """
    pops = result.extras.get("spin_populations")
    got = 0.0 if pops is None else float(pops[metal_index])
    if abs(abs(got) - expected) > tol:
        raise OxidationStateMismatch(
            f"oxidation_state_mismatch: atom {metal_index} carries spin population "
            f"{got:+.2f} where its recorded oxidation and spin state give {expected} "
            f"unpaired electrons (tolerance {tol}); the electrons sit elsewhere, so this "
            f"geometry is not the species its graph names")
    return got

