"""Two cards, one per process: `MOFSBU_DEVICE=cuda:0,cuda:1` and the worker pinning it drives."""
from __future__ import annotations

import sys

import pytest

from mofsbu import runner
from mofsbu.config import (
    compute_device,
    compute_devices,
    device_note,
    physical_device,
    relax_workers,
)


@pytest.fixture(autouse=True)
def _clean(monkeypatch):
    for var in ("MOFSBU_DEVICE", "MOFSBU_DEVICE_PHYSICAL", "MOFSBU_RELAX_WORKERS"):
        monkeypatch.delenv(var, raising=False)


def test_a_pool_is_declared_cards_and_a_process_uses_the_first(monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", " CUDA:0, cuda:1 ")
    assert compute_devices() == ["cuda:0", "cuda:1"]
    assert compute_device() == "cuda:0"
    assert "cuda:0,cuda:1" in device_note() and "uses cuda:0" in device_note()


@pytest.mark.parametrize("value, n", [("cpu", 0), ("cuda", 1), ("cuda:1", 1),
                                      ("cuda:0,cuda:1", 2)])
def test_one_relax_worker_per_declared_card(monkeypatch, value, n):
    monkeypatch.setenv("MOFSBU_DEVICE", value)
    assert relax_workers() == n
    monkeypatch.setenv("MOFSBU_RELAX_WORKERS", "3")
    assert relax_workers() == 3


def test_relax_workers_take_the_cards_in_turn_and_see_only_their_own(monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:0,cuda:1")
    envs = [runner.worker_env("relax", i, split=True) for i in range(3)]
    assert [e["CUDA_VISIBLE_DEVICES"] for e in envs] == ["0", "1", "0"]
    assert {e["MOFSBU_DEVICE"] for e in envs} == {"cuda:0"}
    assert [e["MOFSBU_DEVICE_PHYSICAL"] for e in envs] == ["cuda:0", "cuda:1", "cuda:0"]


def test_a_single_indexed_card_is_pinned_too(monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:1")
    assert runner.worker_env("relax", 0, split=True) == {
        "CUDA_VISIBLE_DEVICES": "1", "MOFSBU_DEVICE": "cuda:0",
        "MOFSBU_DEVICE_PHYSICAL": "cuda:1"}


@pytest.mark.parametrize("value", ["cuda", "cpu", "mps"])
def test_nothing_to_pin_is_left_alone(monkeypatch, value):
    monkeypatch.setenv("MOFSBU_DEVICE", value)
    assert runner.worker_env("relax", 0, split=True) == {}


def test_build_workers_see_no_card_only_when_relax_has_its_own(monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:0,cuda:1")
    assert runner.worker_env("build", 0, split=True) == {"CUDA_VISIBLE_DEVICES": ""}
    assert runner.worker_env("build", 0, split=False) == {}
    monkeypatch.setenv("MOFSBU_DEVICE", "cpu")
    assert runner.worker_env("build", 0, split=True) == {}


def test_a_pinned_worker_records_its_real_card(monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:0")
    assert physical_device() == "cuda:0"
    monkeypatch.setenv("MOFSBU_DEVICE_PHYSICAL", "cuda:1")
    assert physical_device() == "cuda:1"


def test_two_cards_plan_two_relax_workers_and_keep_a_builder(monkeypatch):
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:0,cuda:1")
    plan = runner.plan_workers(6, relaxes=True)
    assert (plan.build, plan.relax) == (4, 2)
    plan = runner.plan_workers(2, relaxes=True)
    assert (plan.build, plan.relax) == (1, 1)


def test_pin_process_confines_this_process_to_its_card(monkeypatch):
    import os

    from mofsbu.config import pin_process

    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.delitem(sys.modules, "torch", raising=False)   # as in a fresh process
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:0,cuda:1")
    assert pin_process(1) == "cuda:1"
    assert os.environ["CUDA_VISIBLE_DEVICES"] == "1"
    assert compute_device() == "cuda:0"
    # already confined: a second call changes nothing
    assert pin_process(0) == "cuda:1" and os.environ["CUDA_VISIBLE_DEVICES"] == "1"


def test_pin_process_refuses_once_cuda_has_started(monkeypatch):
    from types import SimpleNamespace

    from mofsbu.config import pin_process

    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    fake = SimpleNamespace(cuda=SimpleNamespace(is_initialized=lambda: True))
    monkeypatch.setitem(sys.modules, "torch", fake)
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:1")
    with pytest.raises(RuntimeError, match="call it first"):
        pin_process()


def test_pin_process_leaves_cpu_alone(monkeypatch):
    import os

    from mofsbu.config import pin_process

    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.setenv("MOFSBU_DEVICE", "cpu")
    assert pin_process() == "cpu" and "CUDA_VISIBLE_DEVICES" not in os.environ


@pytest.mark.parametrize("visible, n", [("0,1", 2), ("1", 1), ("", 0)])
def test_visible_cards_reads_the_declaration(monkeypatch, visible, n):
    from mofsbu.config import visible_cards

    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", visible)
    assert visible_cards() == n


def test_dft_refuses_a_process_that_sees_two_cards(monkeypatch):
    from mofsbu.energy.dft import DFTBackend

    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0,1")
    b = DFTBackend(device="cuda:0")
    assert not b.available()
    assert "#70" in b.install_hint() and "CUDA_VISIBLE_DEVICES" in b.install_hint()
