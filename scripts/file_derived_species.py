"""File the species that connectivity-changed relaxations became (D31) in an existing registry.

    python scripts/file_derived_species.py --db data/salt_study_k1.db              # dry run
    python scripts/file_derived_species.py --db data/salt_study_k1.db --apply

Every geometry E1 flagged (`connectivity_changed = 1`) that has not been filed yet: the dry
run says what would be filed and what would be refused, by reason; `--apply` files it
(`assembly.derive.file_derived`).  A separation is counted as refused: it waits for the
second-shell probe (D31).  New relax tasks file derived species on their own.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.assembly.derive import DeriveRefused, derive_species, file_derived  # noqa: E402
from mofsbu.registry import BlobStore, Registry, geometry_graph  # noqa: E402

NO_ENERGY = ("SELECT COUNT(*) FROM structures s WHERE EXISTS (SELECT 1 FROM geometries g "
             "WHERE g.structure_id=s.id AND g.energy IS NOT NULL) AND NOT EXISTS (SELECT 1 "
             "FROM geometries g WHERE g.structure_id=s.id AND g.energy IS NOT NULL AND "
             "g.connectivity_changed IS NOT 1)")


def _reason(text: str) -> str:
    """A refusal reason without the atom indices, so like refusals count together."""
    words = text.split(": ", 1)[-1] if ": the skeleton" in text else text
    return words.split(";")[0][:90]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--db", type=Path, required=True)
    p.add_argument("--store", type=Path, default=None,
                   help="blob store (default: 'store' beside the database)")
    p.add_argument("--apply", action="store_true", help="write; default is a dry run")
    a = p.parse_args(argv)
    with Registry(a.db, store=BlobStore(a.store) if a.store else None) as reg:
        reg.migrate("file_derived_species")
        todo = [int(r[0]) for r in reg.conn.execute(
            "SELECT g.id FROM geometries g WHERE g.connectivity_changed = 1 AND NOT EXISTS "
            "(SELECT 1 FROM geometries d WHERE d.derived_from = g.id) ORDER BY g.id")]
        out: dict = {"dry_run": not a.apply, "candidates": len(todo),
                     "no_selectable_energy_before": reg.conn.execute(NO_ENERGY).fetchone()[0]}
        filed = collections.Counter()
        refused = collections.Counter()
        for k, gid in enumerate(todo, 1):
            if a.apply:
                res = file_derived(reg, gid)
                if "refused" in res:
                    refused[_reason(res["refused"])] += 1
                    continue
                filed["single"] += 1
                filed["new structures"] += sum(res["created"])
                filed["existing structures"] += len(res["created"]) - sum(res["created"])
                if k % 200 == 0:
                    reg.conn.commit()
            else:
                verdict = json.loads(reg.conn.execute(
                    "SELECT connectivity_json FROM geometries WHERE id=?",
                    (gid,)).fetchone()[0])
                try:
                    derived = derive_species(geometry_graph(reg, gid), verdict)
                except DeriveRefused as exc:
                    refused[_reason(str(exc))] += 1
                    continue
                if len(derived.pieces) > 1:
                    refused["separation — second shell or leaving, pending"] += 1
                else:
                    filed["single"] += 1
        if a.apply:
            reg.conn.commit()
            out["no_selectable_energy_after"] = reg.conn.execute(NO_ENERGY).fetchone()[0]
        out["filed" if a.apply else "would_file"] = dict(filed)
        out["refused"] = dict(refused.most_common())
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
