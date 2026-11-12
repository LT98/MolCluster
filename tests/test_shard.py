"""Splitting a measurement script over cards (`scripts/_shard.py`) and joining its shards."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import _shard  # noqa: E402
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
