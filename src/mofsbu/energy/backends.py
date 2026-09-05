"""Energy backends behind one protocol (M7).

Three things this module is careful about, because all three have already gone wrong
somewhere in `legacy/`:

* **Charge and spin are inputs, not afterthoughts.**  `legacy/xtb_energy_metal.py` used a
  minimal-spin convention for Ni(II) and a fixed sextet for Fe(III), passed through an
  environment variable.  Here the multiplicity is part of the `MethodSpec` that is stored
  with the number, so two energies can be checked for comparability instead of assumed
  comparable.
* **A backend that cannot see charge says so.**  MACE-MP-0 is a charge- and spin-blind
  potential: it sees elements and positions.  That is fine for ranking neutral isomers and
  wrong for anything the reference scheme does, so `charge_aware` is a property of the
  backend and `energy.reference` refuses to use a blind one on a charged equation.
* **A missing stack is not a missing feature.**  `EnergyBackendUnavailable` means "install
  it"; `NotBuiltYet` means "it isn't written".  Collapsing the two is how you end up
  believing the laptop has no xTB because the workstation was busy.
"""
from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from typing import Any, Protocol, runtime_checkable

from mofsbu._types import EnergyBackendUnavailable, Fidelity, MethodSpec
from mofsbu.versions import ALGO_VERSIONS

def _installed_version(*distributions: str) -> str:
    """The version of whichever distribution is actually installed.

    Ground rule 6 says a stored number carries the version of the recipe that made it,
    and `tblite` exposes no `__version__` attribute at all — asking for one yields
    "unknown", which then becomes a `methods` row that two different builds of the code
    share.  The package metadata is the thing that is always there.
    """
    for name in distributions:
        try:
            return importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            continue
    return "unknown"


# ── results ──────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class EnergyResult:
    """One energy, with the thing that produced it attached.  Ground rule 3."""

    energy: float                       # eV
    method: MethodSpec
    fidelity: Fidelity
    converged: bool = True
    extras: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class RelaxResult:
    """A relaxation: the moved coordinates AND the energy they belong to."""

    symbols: tuple[str, ...]
    positions: Any                      # (n, 3) float array
    energy: float
    method: MethodSpec
    fidelity: Fidelity
    converged: bool
    n_steps: int
    fmax: float
    initial_energy: float | None = None

    @property
    def relaxation_energy(self) -> float | None:
        """How far the geometry fell.  A large drop means the construct was poor."""
        if self.initial_energy is None:
            return None
        return self.energy - self.initial_energy

    def to_xyz(self, comment: str = "") -> str:
        lines = [str(len(self.symbols)), comment]
        for sym, (x, y, z) in zip(self.symbols, self.positions, strict=True):
            lines.append(f"{sym:<2s} {float(x):15.8f} {float(y):15.8f} {float(z):15.8f}")
        return "\n".join(lines) + "\n"


# ── the protocol ─────────────────────────────────────────────────────────────


@runtime_checkable
class EnergyBackend(Protocol):
    """One interface for xTB, an MLIP, and the test double.

    `symbols` + `positions` rather than an ASE `Atoms` object on purpose: ASE is an
    implementation detail of two of the three backends, and the registry stores .xyz.
    """

    name: str
    fidelity: Fidelity
    charge_aware: bool
    spin_aware: bool

    def available(self) -> bool: ...

    def method_spec(self, *, charge: int, multiplicity: int,
                    solvent: str | None = None) -> MethodSpec: ...

    def single_point(self, symbols: Sequence[str], positions: Any, *,
                     charge: int, multiplicity: int,
                     solvent: str | None = None) -> EnergyResult: ...

    def relax(self, symbols: Sequence[str], positions: Any, *,
              charge: int, multiplicity: int, solvent: str | None = None,
              fmax: float = 0.20, steps: int = 250) -> RelaxResult: ...


# ── spin ─────────────────────────────────────────────────────────────────────

# d-electron counts are derived, not tabulated: group number minus oxidation state.
# The first-row transition metals plus the ones this project actually meets.
_GROUP = {
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10,
    "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10,
    "Ag": 11, "Cd": 12,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
}


def d_electrons(symbol: str, charge: int) -> int:
    """d-count for a transition-metal ion.  Raises for anything else."""
    if symbol not in _GROUP:
        raise ValueError(f"{symbol} is not a transition metal this table covers")
    n = _GROUP[symbol] - charge
    if not 0 <= n <= 10:
        raise ValueError(f"{symbol}({charge:+d}) implies d{n}, which is not a real ion")
    return n


def high_spin_multiplicity(symbol: str, charge: int) -> int:
    """Maximum-multiplicity ground state for a d-block ion: 2S+1 with S = unpaired/2.

    Hund's rule on five d orbitals — d0-d5 fill singly, d6-d10 pair up.  This is the
    convention the Fe(III) run used (sextet, five unpaired) written down as a rule
    rather than an environment variable, so the Ni(II) and Fe(III) cases come out of
    the same line of code.  Ground rule 5: the caller chooses high or low spin; nothing
    here guesses which one the chemistry wants.
    """
    n = d_electrons(symbol, charge)
    unpaired = n if n <= 5 else 10 - n
    return unpaired + 1


def minimal_multiplicity(symbols: Sequence[str], charge: int) -> int:
    """Singlet if the electron count is even, doublet if odd.

    The `legacy` default.  Kept because it is what the archived Ni(II) numbers used and
    the M7 regression has to reproduce them, not because it is a good convention for
    an open-shell metal.
    """
    electrons = sum(_atomic_number(s) for s in symbols) - charge
    return 1 if electrons % 2 == 0 else 2


_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
    "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Hf": 72,
    "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83,
}


def _atomic_number(symbol: str) -> int:
    try:
        return _Z[symbol]
    except KeyError:
        raise ValueError(f"unknown element symbol {symbol!r}") from None


def electron_count(symbols: Sequence[str], charge: int) -> int:
    return sum(_atomic_number(s) for s in symbols) - charge


def check_spin(symbols: Sequence[str], charge: int, multiplicity: int) -> None:
    """A multiplicity has to be reachable from the electron count.

    An even electron count cannot be a doublet.  xTB will happily accept the request and
    return a number for it, which is precisely the kind of plausible-looking answer this
    project keeps having to unlearn.
    """
    if multiplicity < 1:
        raise ValueError(f"multiplicity {multiplicity} is not physical")
    electrons = electron_count(symbols, charge)
    unpaired = multiplicity - 1
    if (electrons - unpaired) % 2 != 0:
        raise ValueError(
            f"{electrons} electrons at charge {charge:+d} cannot give multiplicity "
            f"{multiplicity} ({unpaired} unpaired leaves an odd number to pair)")


# ── backends ─────────────────────────────────────────────────────────────────


class _Base:
    name = "base"
    fidelity = Fidelity.RAW
    charge_aware = True
    spin_aware = True
    code = "base"
    method = "base"

    def available(self) -> bool:
        return True

    def code_version(self) -> str:
        return "0"

    def method_spec(self, *, charge: int, multiplicity: int,
                    solvent: str | None = None) -> MethodSpec:
        extras: dict[str, Any] = {"algo": ALGO_VERSIONS["energy_backends"]}
        if not self.charge_aware:
            # Recorded in the method row, so a stored number carries the caveat with it
            # rather than relying on whoever reads it later remembering.
            extras["charge_blind"] = True
        return MethodSpec(code=self.code, code_version=self.code_version(),
                          method=self.method, solvent=solvent,
                          charge=charge, multiplicity=multiplicity, extras=extras)

    def _require(self) -> None:
        if not self.available():
            raise EnergyBackendUnavailable(
                f"{self.name} backend: {self.install_hint()}")

    def install_hint(self) -> str:
        return "not installed"


class NullBackend(_Base):
    """A deterministic, deliberately non-physical backend for tests.

    Ground rule 8 says a stub must never return a plausible value.  This returns a
    number, so it is fenced instead: `code='null'` lands in the `methods` table, the
    energies are absurd on their face (whole eV per electron), and `energy.reference`
    refuses to score a reaction with it unless the caller passes `allow_null=True`.
    Nothing that reaches the registry from a real run can quietly be one of these.
    """

    name = "null"
    fidelity = Fidelity.RAW
    charge_aware = True
    spin_aware = True
    code = "null"
    method = "null-not-physical"

    def code_version(self) -> str:
        return "1"

    def single_point(self, symbols, positions, *, charge, multiplicity,
                     solvent=None) -> EnergyResult:
        check_spin(list(symbols), charge, multiplicity)
        energy = -1.0 * electron_count(list(symbols), charge) - 0.5 * (multiplicity - 1)
        return EnergyResult(energy=energy,
                            method=self.method_spec(charge=charge, multiplicity=multiplicity,
                                                    solvent=solvent),
                            fidelity=self.fidelity, converged=True)

    def relax(self, symbols, positions, *, charge, multiplicity, solvent=None,
              fmax=0.20, steps=250) -> RelaxResult:
        sp = self.single_point(symbols, positions, charge=charge,
                               multiplicity=multiplicity, solvent=solvent)
        return RelaxResult(symbols=tuple(symbols), positions=positions,
                           energy=sp.energy, method=sp.method, fidelity=self.fidelity,
                           converged=True, n_steps=0, fmax=0.0,
                           initial_energy=sp.energy)


class XTBBackend(_Base):
    """GFN2-xTB through tblite's ASE calculator.

    The port of `legacy/xtb_energy.py` + `xtb_energy_metal.py` + the ALPB half of
    `solvation_model.py`, with the environment variables replaced by arguments and the
    reference arithmetic removed — that now lives in `energy.reference`, where it can be
    checked for balance before it is believed.
    """

    name = "xtb"
    fidelity = Fidelity.XTB
    charge_aware = True
    spin_aware = True
    code = "tblite"
    method = "GFN2-xTB"

    #  ALPB solvent keys tblite accepts, mapped from the names this project uses.
    SOLVENTS = {"water": "water", "ethanol": "ethanol", "dmf": "dmf",
                "methanol": "methanol", "acetonitrile": "acetonitrile",
                "thf": "thf", "dmso": "dmso", "chloroform": "chcl3",
                "toluene": "toluene", "hexane": "hexane"}

    def __init__(self, *, accuracy: float = 1.0, max_iterations: int = 400) -> None:
        self.accuracy = accuracy
        self.max_iterations = max_iterations

    def available(self) -> bool:
        try:
            import ase  # noqa: F401
            import tblite.ase  # noqa: F401
        except Exception:
            return False
        return True

    def install_hint(self) -> str:
        return "conda install -c conda-forge tblite-python ase"

    def code_version(self) -> str:
        return _installed_version("tblite")

    def method_spec(self, *, charge, multiplicity, solvent=None) -> MethodSpec:
        spec = super().method_spec(charge=charge, multiplicity=multiplicity,
                                   solvent=self._solvent_key(solvent))
        return replace(spec, extras={**spec.extras, "accuracy": self.accuracy,
                                     "max_iterations": self.max_iterations})

    def _solvent_key(self, solvent: str | None) -> str | None:
        if solvent is None:
            return None
        key = solvent.strip().lower()
        if key not in self.SOLVENTS:
            raise ValueError(
                f"unknown solvent {solvent!r}; ALPB knows {sorted(self.SOLVENTS)}")
        return key

    def _calc(self, charge: int, multiplicity: int, solvent: str | None):
        from tblite.ase import TBLite

        kwargs: dict[str, Any] = dict(method="GFN2-xTB", charge=charge,
                                      multiplicity=multiplicity, verbosity=0,
                                      accuracy=self.accuracy,
                                      max_iterations=self.max_iterations)
        key = self._solvent_key(solvent)
        if key is not None:
            kwargs["solvation"] = ("alpb", self.SOLVENTS[key])
        return TBLite(**kwargs)

    def _atoms(self, symbols, positions):
        from ase import Atoms

        atoms = Atoms(symbols=list(symbols), positions=positions)
        atoms.pbc = False
        return atoms

    def single_point(self, symbols, positions, *, charge, multiplicity,
                     solvent=None) -> EnergyResult:
        self._require()
        check_spin(list(symbols), charge, multiplicity)
        atoms = self._atoms(symbols, positions)
        atoms.calc = self._calc(charge, multiplicity, solvent)
        energy = float(atoms.get_potential_energy())
        return EnergyResult(energy=energy,
                            method=self.method_spec(charge=charge, multiplicity=multiplicity,
                                                    solvent=solvent),
                            fidelity=self.fidelity, converged=True)

    def relax(self, symbols, positions, *, charge, multiplicity, solvent=None,
              fmax=0.20, steps=250) -> RelaxResult:
        self._require()
        check_spin(list(symbols), charge, multiplicity)
        from ase.optimize import LBFGS

        atoms = self._atoms(symbols, positions)
        atoms.calc = self._calc(charge, multiplicity, solvent)
        e0 = float(atoms.get_potential_energy())
        opt = LBFGS(atoms, logfile=None)
        converged = bool(opt.run(fmax=fmax, steps=steps))
        forces = atoms.get_forces()
        reached = float(max(math.sqrt(float(f @ f)) for f in forces))
        return RelaxResult(
            symbols=tuple(symbols), positions=atoms.get_positions(),
            energy=float(atoms.get_potential_energy()),
            method=self.method_spec(charge=charge, multiplicity=multiplicity, solvent=solvent),
            fidelity=self.fidelity, converged=converged,
            n_steps=int(opt.get_number_of_steps()), fmax=reached, initial_energy=e0)


class MACEBackend(_Base):
    """MACE-MP-0, a universal MLIP — fast, and blind to charge and spin.

    `charge_aware = False` is the whole point of the flag.  MACE sees elements and
    positions; the formal charges this project tracks are invisible to it.  That makes
    it useful for polishing a construct and for ranking species of the SAME charge, and
    unusable for the reference scheme, which is exactly what `energy.reference` enforces
    rather than leaving to a docstring nobody reads at 2 a.m.
    """

    name = "mace"
    fidelity = Fidelity.ML
    charge_aware = False
    spin_aware = False
    code = "mace"
    method = "MACE-MP-0"

    def __init__(self, *, model: str = "medium", device: str = "cpu",
                 default_dtype: str = "float64") -> None:
        self.model = model
        self.device = device
        self.default_dtype = default_dtype
        self._calc_cache: Any = None

    def available(self) -> bool:
        try:
            import ase  # noqa: F401
            import mace  # noqa: F401
            import torch  # noqa: F401
        except Exception:
            return False
        return True

    def install_hint(self) -> str:
        return ("pip install torch --index-url https://download.pytorch.org/whl/cpu "
                "&& pip install mace-torch ase")

    def code_version(self) -> str:
        return f"{_installed_version('mace-torch', 'mace')}/{self.model}"

    def method_spec(self, *, charge, multiplicity, solvent=None) -> MethodSpec:
        if solvent is not None:
            raise ValueError("MACE-MP-0 has no solvation model; use xtb for a continuum")
        spec = super().method_spec(charge=charge, multiplicity=multiplicity, solvent=None)
        return replace(spec, extras={**spec.extras, "model": self.model})

    def _calculator(self):
        if self._calc_cache is None:
            from mace.calculators import mace_mp

            self._calc_cache = mace_mp(model=self.model, device=self.device,
                                       default_dtype=self.default_dtype)
        return self._calc_cache

    def _atoms(self, symbols, positions):
        from ase import Atoms

        atoms = Atoms(symbols=list(symbols), positions=positions)
        atoms.pbc = False
        return atoms

    def single_point(self, symbols, positions, *, charge, multiplicity,
                     solvent=None) -> EnergyResult:
        self._require()
        atoms = self._atoms(symbols, positions)
        atoms.calc = self._calculator()
        return EnergyResult(
            energy=float(atoms.get_potential_energy()),
            method=self.method_spec(charge=charge, multiplicity=multiplicity, solvent=solvent),
            fidelity=self.fidelity, converged=True)

    def relax(self, symbols, positions, *, charge, multiplicity, solvent=None,
              fmax=0.05, steps=250) -> RelaxResult:
        self._require()
        from ase.optimize import LBFGS

        atoms = self._atoms(symbols, positions)
        atoms.calc = self._calculator()
        e0 = float(atoms.get_potential_energy())
        opt = LBFGS(atoms, logfile=None)
        converged = bool(opt.run(fmax=fmax, steps=steps))
        forces = atoms.get_forces()
        reached = float(max(math.sqrt(float(f @ f)) for f in forces))
        return RelaxResult(
            symbols=tuple(symbols), positions=atoms.get_positions(),
            energy=float(atoms.get_potential_energy()),
            method=self.method_spec(charge=charge, multiplicity=multiplicity, solvent=solvent),
            fidelity=self.fidelity, converged=converged,
            n_steps=int(opt.get_number_of_steps()), fmax=reached, initial_energy=e0)


# ── selection ────────────────────────────────────────────────────────────────

_BACKENDS: dict[str, Any] = {"xtb": XTBBackend, "mace": MACEBackend, "null": NullBackend}

# Which backend serves each rung of the ladder.  FF is RDKit's MMFF, which lives in
# `geometry.embed` and is not an EnergyBackend — it produces geometries, not comparable
# energies.  DFT has no backend: there is no external code wired up, and inventing one
# that silently ran xTB instead would be the exact failure ground rule 8 exists for.
_BY_FIDELITY = {Fidelity.ML: "mace", Fidelity.XTB: "xtb"}


def get_backend(name: str, **kwargs: Any) -> EnergyBackend:
    try:
        cls = _BACKENDS[name]
    except KeyError:
        raise ValueError(f"unknown backend {name!r}; have {sorted(_BACKENDS)}") from None
    return cls(**kwargs)


def backend_for(fidelity: Fidelity, **kwargs: Any) -> EnergyBackend:
    """The backend that produces geometries at this rung, or a reason why not."""
    if fidelity not in _BY_FIDELITY:
        from mofsbu.assembly.join import NotBuiltYet

        raise NotBuiltYet(
            f"no energy backend serves {fidelity.name}. "
            + ("DFT needs an external code wired up (M7 leaves the slot open on purpose)."
               if fidelity is Fidelity.DFT else
               f"{fidelity.name} geometries come from the constructor, not from a backend."))
    return get_backend(_BY_FIDELITY[fidelity], **kwargs)


def available_backends() -> dict[str, bool]:
    """What this machine can actually run.  The UI reads this rather than guessing."""
    return {name: get_backend(name).available() for name in _BACKENDS}
