"""Regenerate tests/golden_hashes.json.

Run deliberately, never automatically: a changed golden file means every stored key
in the registry is stale, so it must land in the same commit as the ALGO_VERSIONS
bump that explains it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "tests")]

import fixtures as fx                     # noqa: E402
from mofsbu.identity import identity      # noqa: E402
from mofsbu.versions import ALGO_VERSIONS  # noqa: E402

out = {"algo_versions": ALGO_VERSIONS, "l0": {}, "l1": {}}
for name, build in fx.ALL.items():
    ident = identity(build())
    out["l0"][name] = ident["l0"]
    out["l1"][name] = ident["l1"]

path = ROOT / "tests" / "golden_hashes.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(f"wrote {path} ({len(out['l1'])} fixtures)")
