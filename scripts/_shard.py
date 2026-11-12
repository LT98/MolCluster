r"""Split one long measurement over several processes, one card each, and join the outputs.

    MOFSBU_DEVICE=cuda:0,cuda:1 python scripts/funnel_eval.py ... --shard 0/2 &
    MOFSBU_DEVICE=cuda:0,cuda:1 python scripts/funnel_eval.py ... --shard 1/2 &
    python scripts/funnel_eval.py --out funnel.json \
        --merge funnel.shard0of2.json funnel.shard1of2.json

Shard i takes items i, i+n, i+2n ... of the same seeded pick and runs on card i mod n of the
declared pool (`config.pin_process`), writing `<out stem>.shard<i>of<n><suffix>`.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def add_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--shard", default=None, metavar="I/N",
                   help="run every N-th item starting at I, on card I mod (pool size)")
    p.add_argument("--merge", nargs="+", type=Path, default=None, metavar="SHARD_JSON",
                   help="join shard outputs into --out and recompute the summary; runs nothing")


def parse(text: str | None) -> tuple[int, int] | None:
    if text is None:
        return None
    m = re.fullmatch(r"(\d+)/(\d+)", text.strip())
    if m is None or not 0 <= int(m.group(1)) < int(m.group(2)):
        raise SystemExit(f"--shard {text!r}: expected I/N with 0 <= I < N")
    return int(m.group(1)), int(m.group(2))


def start(a: argparse.Namespace) -> tuple[int, int] | None:
    """Pin this process to its shard's card and point `a.out` at the shard file.

    Call before any backend is built: the card has to be chosen before CUDA starts.
    """
    from mofsbu.config import pin_process

    shard = parse(a.shard)
    i = 0 if shard is None else shard[0]
    device = pin_process(i)
    if shard is not None:
        a.out = a.out.with_name(f"{a.out.stem}.shard{shard[0]}of{shard[1]}{a.out.suffix}")
        print(f"shard {shard[0]}/{shard[1]} on {device} -> {a.out}", flush=True)
    return shard


def take(items: list, shard: tuple[int, int] | None) -> list:
    return items if shard is None else items[shard[0]::shard[1]]


def load(paths: list[Path], rows_key: str, same: tuple[str, ...]) -> tuple[dict, list]:
    """The first file's header and every file's rows; refuses shards of different runs."""
    docs = [json.loads(p.read_text()) for p in paths]
    for key in same:
        values = {json.dumps(d.get(key), sort_keys=True) for d in docs}
        if len(values) > 1:
            raise SystemExit(f"--merge: the shards disagree on {key!r}; they are not one run")
    rows = [r for d in docs for r in d[rows_key]]
    ids = [r["sid"] for r in rows]
    if len(ids) != len(set(ids)):
        raise SystemExit("--merge: a structure appears in more than one shard")
    return docs[0], rows
