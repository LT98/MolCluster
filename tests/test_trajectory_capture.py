"""F0: a relaxation reports every optimiser point, and funnel_eval keeps it as a trajectory file."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mofsbu.energy.backends import NullBackend, _Base, _lbfgs

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from funnel_eval import Trajectories  # noqa: E402

# Three argon atoms off their Lennard-Jones minimum: a few dozen LBFGS steps, no model needed.
AR3 = (["Ar", "Ar", "Ar"], [[0.0, 0.0, 0.0], [3.3, 0.0, 0.0], [0.0, 4.2, 0.0]])


class LJBackend(_Base):
    """Real LBFGS on ASE's Lennard-Jones calculator, through the shared `_lbfgs`."""

    name = code = method = "lj-test"

    def relax(self, symbols, positions, *, charge, multiplicity, solvent=None,
              fmax=0.05, steps=250, on_step=None):
        from ase import Atoms
        from ase.calculators.lj import LennardJones

        atoms = Atoms(symbols=list(symbols), positions=positions)
        atoms.calc = LennardJones(epsilon=0.0104, sigma=3.4, rc=10.0)
        return _lbfgs(self, atoms, symbols, charge, multiplicity, solvent, fmax, steps, on_step)


def _graph():
    return SimpleNamespace(charge=0, multiplicity=1)


def test_on_step_sees_the_start_every_step_and_the_returned_geometry():
    seen: list[tuple[float, float, np.ndarray]] = []
    r = LJBackend().relax(*AR3, charge=0, multiplicity=1, fmax=0.01, steps=200,
                          on_step=lambda e, f, x: seen.append((e, f, x)))
    assert r.converged and r.n_steps > 3
    assert len(seen) == r.n_steps + 1
    np.testing.assert_allclose(seen[0][2], AR3[1])
    assert seen[0][0] == pytest.approx(r.initial_energy)
    assert seen[-1][0] == pytest.approx(r.energy)
    assert seen[-1][1] == pytest.approx(r.fmax) and seen[-1][1] < 0.01
    np.testing.assert_allclose(seen[-1][2], r.positions)
    # each point is a copy, not a view of the moving atoms
    assert not np.allclose(seen[0][2], seen[-1][2])


def test_on_step_left_out_relaxes_the_same():
    a = LJBackend().relax(*AR3, charge=0, multiplicity=1, fmax=0.01, steps=200)
    b = LJBackend().relax(*AR3, charge=0, multiplicity=1, fmax=0.01, steps=200,
                          on_step=lambda *_: None)
    assert a.energy == b.energy and a.n_steps == b.n_steps
    np.testing.assert_array_equal(a.positions, b.positions)


def test_null_backend_reports_its_single_point():
    seen = []
    NullBackend().relax(*AR3, charge=0, multiplicity=1, on_step=lambda *p: seen.append(p))
    assert len(seen) == 1 and seen[0][1] == 0.0


def test_trajectory_file_keeps_every_energy_and_strided_frames_with_the_last(tmp_path):
    traj = Trajectories(tmp_path / "t", every=4)
    r = traj.relax(LJBackend(), "lj", "base", 7, 42, AR3[0], np.array(AR3[1]), _graph(),
                   0.01, 200)
    z = np.load(tmp_path / "t" / "7_42_base_lj.npz")
    n = r["steps"] + 1
    assert z["energy"].shape == (n,) and z["fmax"].shape == (n,)
    assert z["energy"][-1] == pytest.approx(r["energy"])
    idx = list(z["frame_index"])
    assert idx[:-1] == list(range(0, n, 4))[:len(idx) - 1] and idx[-1] == n - 1
    assert z["coords"].shape == (len(idx), 3, 3)
    np.testing.assert_allclose(z["coords"][-1], r["xyz"])
    np.testing.assert_allclose(z["start"], AR3[1])
    meta = json.loads(str(z["meta"]))
    assert meta | {"symbols": None} == meta | {
        "sid": 7, "start_geometry": 42, "stage": "base", "model": "lj", "method": "lj-test",
        "converged": True, "every": 4, "symbols": None}


def test_without_a_directory_nothing_is_written(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r = Trajectories(None, every=4).relax(LJBackend(), "lj", "base", 7, 42, AR3[0],
                                          np.array(AR3[1]), _graph(), 0.01, 200)
    assert "energy" in r and not list(tmp_path.iterdir())
