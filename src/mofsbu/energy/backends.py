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

    def __init__(self, *, fidelity: Fidelity = Fidelity.RAW) -> None:
        # A test that exercises the fidelity LADDER needs a double that can sit on a rung
        # above RAW.  Safe to allow because `code` stays "null": the reference scheme
        # refuses these energies by code, not by rung, so a masquerading fidelity cannot
        # smuggle a non-physical number into a reaction energy.
        self.fidelity = fidelity

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

    `charge_aware = False` is the whole point of the flag.  MACE-MP-0 sees elements and
    positions; the formal charges this project tracks are invisible to it.  That makes
    it useful for polishing a construct and for ranking species of the SAME charge, and
    unusable for the reference scheme, which is exactly what `energy.reference` enforces
    rather than leaving to a docstring nobody reads at 2 a.m.

    `MACEOmolBackend` below is the same code path with a different foundation model and
    the flags flipped.  Everything the two share lives here; everything that differs is
    a class attribute, so "does this model see charge?" is answered in one place and
    travels into the `methods` row with every number.
    """

    name = "mace"
    fidelity = Fidelity.ML
    charge_aware = False
    spin_aware = False
    code = "mace"
    method = "MACE-MP-0"

    #: the factory in `mace.calculators` that loads this foundation model
    loader = "mace_mp"
    #: default checkpoint size for that factory
    default_model = "medium"
    #: what the model was trained on — recorded in the method row, because two MLIPs on
    #: the same rung of the ladder are still two different theories
    training_set = "Materials Project (MPtrj)"
    min_mace_version = "0.3.6"

    def __init__(self, *, model: str | None = None, device: str | None = None,
                 default_dtype: str = "float64") -> None:
        from mofsbu.config import compute_device

        # `device="cpu"` used to be the hard default, so a workstation with a GPU ran the
        # MLIP on its CPU and the only symptom was a card that never warmed up.  The
        # device is DECLARED (`MOFSBU_DEVICE`), never detected, for the same reason
        # parallelism is opt-in.
        device = compute_device() if device is None else device
        self.model = self.default_model if model is None else model
        self.device = device
        self.default_dtype = default_dtype
        self._calc_cache: Any = None

    # -- availability ---------------------------------------------------

    def available(self) -> bool:
        try:
            import ase  # noqa: F401
            import torch  # noqa: F401
            from mace import calculators
        except Exception:
            return False
        # An installed mace-torch that predates this foundation model is NOT the same
        # failure as no mace-torch at all, and collapsing the two is how "MACE is
        # installed" becomes "MACE-OMOL-0 will run".  The loader either exists or it
        # does not, and `install_hint` says which upgrade fixes it.
        return hasattr(calculators, self.loader)

    def install_hint(self) -> str:
        return ("pip install torch --index-url https://download.pytorch.org/whl/cpu "
                f"&& pip install 'mace-torch>={self.min_mace_version}' ase")

    def code_version(self) -> str:
        return f"{_installed_version('mace-torch', 'mace')}/{self.model}"

    # -- the method row -------------------------------------------------

    def method_spec(self, *, charge, multiplicity, solvent=None) -> MethodSpec:
        if solvent is not None:
            raise ValueError(
                f"{self.method} has no solvation model; use xtb for a continuum")
        spec = super().method_spec(charge=charge, multiplicity=multiplicity, solvent=None)
        extras = {**spec.extras, "model": self.model, "device": self.device,
                  "training_set": self.training_set}
        if not self.spin_aware:
            # Same reasoning as `charge_blind`: the caveat is stored WITH the number.
            # Without it, an MP-0 energy for a sextet and one for a singlet are the same
            # float and nothing in the database says why.
            extras["spin_blind"] = True
        return replace(spec, extras=extras)

    # -- running it -----------------------------------------------------

    def _calculator(self):
        if self._calc_cache is None:
            from mace import calculators

            factory = getattr(calculators, self.loader, None)
            if factory is None:                     # pragma: no cover - guarded by available()
                raise EnergyBackendUnavailable(
                    f"{self.name} backend: {self.install_hint()}")
            self._calc_cache = factory(model=self.model, device=self.device,
                                       default_dtype=self.default_dtype)
        return self._calc_cache

    def _atoms(self, symbols, positions, *, charge: int, multiplicity: int):
        """Build the ASE object, handing the model the electronic state IF it takes one.

        This is the one place the two MACE backends genuinely differ at run time.  A
        charge-aware model reads `atoms.info["total_charge"]` and
        `atoms.info["total_spin"]` (MACE maps those onto its internal `charge`/`spin`
        config keys; OMol25's `spin` is the MULTIPLICITY, 2S+1, not the unpaired count).
        Setting them on a charge-blind model would be worse than useless: the keys are
        ignored, and the run would look configured.
        """
        from ase import Atoms

        atoms = Atoms(symbols=list(symbols), positions=positions)
        atoms.pbc = False
        if self.charge_aware:
            atoms.info["total_charge"] = int(charge)
        if self.spin_aware:
            atoms.info["total_spin"] = int(multiplicity)
        return atoms

    def _prepare(self, symbols, positions, charge: int, multiplicity: int, solvent):
        # The REQUEST is validated before the machine is: "you asked for a multiplicity
        # this electron count cannot reach" is true on the laptop and on the workstation,
        # and hiding it behind "mace-torch is not installed here" turns a permanent bug
        # into an environment problem.
        if solvent is not None:
            raise ValueError(
                f"{self.method} has no solvation model; use xtb for a continuum")
        if self.spin_aware:
            # A model that is GIVEN the multiplicity has to be given a reachable one.
            # MP-0 is not asked, so there is nothing to check and pretending otherwise
            # would reject perfectly good MP-0 work on a number the model never sees.
            check_spin(list(symbols), charge, multiplicity)
        self._require()
        atoms = self._atoms(symbols, positions, charge=charge, multiplicity=multiplicity)
        atoms.calc = self._calculator()
        return atoms

    def single_point(self, symbols, positions, *, charge, multiplicity,
                     solvent=None) -> EnergyResult:
        atoms = self._prepare(symbols, positions, charge, multiplicity, solvent)
        return EnergyResult(
            energy=float(atoms.get_potential_energy()),
            method=self.method_spec(charge=charge, multiplicity=multiplicity, solvent=solvent),
            fidelity=self.fidelity, converged=True)

    def relax(self, symbols, positions, *, charge, multiplicity, solvent=None,
              fmax=0.05, steps=250) -> RelaxResult:
        from ase.optimize import LBFGS

        atoms = self._prepare(symbols, positions, charge, multiplicity, solvent)
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


class MACEOmolBackend(MACEBackend):
    """MACE-OMOL-0 — the same architecture, trained on OMol25, and NOT charge-blind.

    Why this class exists at all, given `MACEBackend` already runs MACE: OMol25 labels
    carry total charge and spin multiplicity, and the model takes both as inputs.  That
    single difference moves the model across the line `energy.reference` draws.  MP-0
    cannot be used on the charged species this project is made of — a bare Ni(2+), a
    carboxylate anion — and OMOL-0 can, at ML cost instead of xTB cost.

    Three things it is NOT:

    * **not a drop-in replacement for MP-0 in stored data.**  Different training set,
      different reference (wB97M-V/def2-TZVPD total energies, eV), different absolute
      scale.  An MP-0 energy and an OMOL-0 energy share a rung on the fidelity ladder
      and nothing else; `MethodSpec.same_theory` already refuses to mix them, and the
      registry's best-geometry rule no longer compares them by magnitude.
    * **not solvated.**  Gas phase, like MP-0.  A solvent request still raises.
    * **not exempt from the spin convention.**  It is handed a multiplicity, so it is
      handed a REACHABLE one: `check_spin` runs here and does not for MP-0.
    """

    name = "mace_omol"
    charge_aware = True
    spin_aware = True
    method = "MACE-OMOL-0"

    loader = "mace_omol"
    default_model = "extra_large"
    training_set = "OMol25 (wB97M-V/def2-TZVPD)"
    #: `mace_omol` landed in mace-torch 0.3.14; an older install imports fine and has no
    #: OMOL loader, which `available()` reports as "not installed" rather than crashing
    #: halfway through a 500-structure run.
    min_mace_version = "0.3.14"


# ── selection ────────────────────────────────────────────────────────────────

_BACKENDS: dict[str, Any] = {"xtb": XTBBackend, "mace": MACEBackend,
                             "mace_omol": MACEOmolBackend, "null": NullBackend}

#: Backend keys that serve the ML rung.  More than one, which is the whole point: the
#: ladder says how good a number is, not which theory produced it.
ML_BACKENDS = ("mace", "mace_omol")

# Which backend serves each rung of the ladder.  FF is RDKit's MMFF, which lives in
# `geometry.embed` and is not an EnergyBackend — it produces geometries, not comparable
# energies.  DFT has no backend: there is no external code wired up, and inventing one
# that silently ran xTB instead would be the exact failure ground rule 8 exists for.
#
# ML is deliberately absent from this table.  It has two backends and picking between
# them is a DECLARED choice (`MOFSBU_ML_MODEL`, or `BuildSpec.ml_model`), not a lookup —
# see `ml_backend_key` below.
_BY_FIDELITY = {Fidelity.XTB: "xtb"}


def ml_backend_key(ml_model: str | None = None) -> str:
    """Which ML backend a request means: the one named, else the one declared.

    `None` is not "whichever" — it is "whatever this machine declares", which
    `config.ml_backend()` answers and which lands in the `methods` row either way, so a
    number never loses the name of the model that made it.
    """
    from mofsbu.config import resolve_ml_backend

    return resolve_ml_backend(ml_model)


def get_backend(name: str, **kwargs: Any) -> EnergyBackend:
    try:
        cls = _BACKENDS[name]
    except KeyError:
        raise ValueError(f"unknown backend {name!r}; have {sorted(_BACKENDS)}") from None
    return cls(**kwargs)


def backend_for(fidelity: Fidelity, *, ml_model: str | None = None,
                **kwargs: Any) -> EnergyBackend:
    """The backend that produces geometries at this rung, or a reason why not."""
    if fidelity is Fidelity.ML:
        return get_backend(ml_backend_key(ml_model), **kwargs)
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
