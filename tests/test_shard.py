"""Splitting a measurement script over cards (`scripts/_shard.py`) and joining its shards."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import _shard  # noqa: E402
import dft_reopt  # noqa: E402
import funnel_eval  # noqa: E402


@pytest.mark.parametrize("text, want", [(None, None), ("0/2", (0, 2)), (" 1/2 ", (1, 2))])
def test_parse(text, want):
    assert _shard.parse(text) == want


@pytest.mark.parametrize("text", ["2/2", "1", "-1/2", "a/b"])
def test_parse_refuses(text):
    with pytest.raises(SystemExit):
        _shard.parse(text)


def test_shards_cover_the_pick_once():
    items = list(range(7))
    parts = [_shard.take(items, (i, 3)) for i in range(3)]
    assert sorted(x for p in parts for x in p) == items
    assert _shard.take(items, None) == items


def _funnel_row(sid, regret, seconds):
    return {"sid": sid, "charge": 0, "n_atoms": 20,
            "baseline": {"seconds": 10.0, "stored": 2, "failed": 0},
            "s": {"per_k": {"1": {"regret": regret, "seconds": seconds, "stored": 3}},
                  "seconds_small": seconds, "rank_of_large_best": 1}}


def _write(path, rows, **head):
    path.write_text(json.dumps({"large": "L", "small": ["s"], "fmax": 0.05, "starts": 2,
                                "k_max": 1, "seed": 0, **head, "results": rows}))
    return path


def test_funnel_merge_summarises_every_shard(tmp_path):
    a = _write(tmp_path / "f.shard0of2.json", [_funnel_row(1, 0.0, 2.0)])
    b = _write(tmp_path / "f.shard1of2.json", [_funnel_row(2, 0.2, 4.0)])
    out = tmp_path / "f.json"
    assert funnel_eval.main(["--out", str(out), "--merge", str(a), str(b)]) == 0
    doc = json.loads(out.read_text())
    assert [r["sid"] for r in doc["results"]] == [1, 2]
    k1 = doc["summary"]["s"]["k=1"]
    assert k1["median_regret"] == pytest.approx(0.1) and k1["worst_regret"] == 0.2
    assert k1["missed_over_0.05eV"] == 1 and k1["time_vs_all_large"] == pytest.approx(0.3)


def test_merge_refuses_shards_of_different_runs(tmp_path):
    a = _write(tmp_path / "a.json", [_funnel_row(1, 0.0, 2.0)])
    b = _write(tmp_path / "b.json", [_funnel_row(2, 0.0, 2.0)], fmax=0.2)
    with pytest.raises(SystemExit, match="fmax"):
        funnel_eval.main(["--out", str(tmp_path / "f.json"), "--merge", str(a), str(b)])


def test_merge_refuses_a_structure_in_two_shards(tmp_path):
    a = _write(tmp_path / "a.json", [_funnel_row(1, 0.0, 2.0)])
    b = _write(tmp_path / "b.json", [_funnel_row(1, 0.0, 2.0)])
    with pytest.raises(SystemExit, match="more than one shard"):
        funnel_eval.main(["--out", str(tmp_path / "f.json"), "--merge", str(a), str(b)])


def _reopt(path, sid, de):
    rec = {"ml": {"energy": -1.0, "converged": True},
           "dft_at_ml_min": {"energy": -2.0 + de, "converged": True},
           "dft_min": {"energy": -2.0, "converged": True}, "dft_relax": {"converged": True},
           "dE_ml_min_minus_dft_min": de, "heavy_rmsd_ml_dft": 0.1}
    path.write_text(json.dumps({
        "dft": "pbe/def2-svp", "contraction_engine": "x", "relaxers": {"r": "R"},
        "fmax": 0.05, "steps": 500, "dft_fmax": 0.05, "dft_steps": 150,
        "acceptance_eV": 0.05, "records": [{"sid": sid, "relaxers": {"r": rec}}]}))
    return path


def test_dft_reopt_merge_summarises_every_shard(tmp_path):
    a = _reopt(tmp_path / "r.shard0of2.json", 1, 0.02)
    b = _reopt(tmp_path / "r.shard1of2.json", 2, 0.04)
    out = tmp_path / "r.json"
    assert dft_reopt.main(["--out", str(out), "--merge", str(a), str(b)]) == 0
    doc = json.loads(out.read_text())
    assert [r["sid"] for r in doc["records"]] == [1, 2]
    r = doc["summary"]["r"]
    assert r["n"] == r["n_comparable"] == 2 and r["mean_dE"] == pytest.approx(0.03)
    assert r["within_acceptance"] is True
    assert doc["dft"] == "pbe/def2-svp" and "contraction_engine" in doc


def test_a_preset_card_is_reported_not_silently_shared(tmp_path, monkeypatch, capsys):
    import argparse

    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "1")
    monkeypatch.setenv("MOFSBU_DEVICE", "cuda:0,cuda:1")
    a = argparse.Namespace(shard="0/2", out=tmp_path / "f.json")
    assert _shard.start(a) == (0, 2)
    assert a.out.name == "f.shard0of2.json"
    assert "already set and wins" in capsys.readouterr().out
