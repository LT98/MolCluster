"""Link geometries stored before D30 to their structure's graph (B28).

    python scripts/recover_atom_links.py --db data/registry.db                  # dry run
    python scripts/recover_atom_links.py --db data/registry.db --apply
    python scripts/recover_atom_links.py --db data/registry.db --apply --connectivity

Dry run by default: reports how many geometries are in the stored numbering, how many in
another build's, and which cannot be matched, and writes nothing.  `--apply` records each
link (coordinates are never touched) and drops site states that were keyed through the wrong
map.  `--connectivity` then runs the E1 backfill on rows that still have no verdict.  See
`mofsbu.registry.recover`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.registry import BlobStore, Registry, backfill_connectivity   # noqa: E402
from mofsbu.registry.recover import recover_atom_links                    # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--db", type=Path, required=True)
    p.add_argument("--store", type=Path, default=None,
                   help="blob store (default: 'store' beside the database)")
    p.add_argument("--apply", action="store_true", help="write; default is a dry run")
    p.add_argument("--connectivity", action="store_true",
                   help="after --apply, run the E1 backfill")
    a = p.parse_args(argv)
    store = BlobStore(a.store) if a.store else None
    with Registry(a.db, store=store) as reg:
        reg.migrate("recover_atom_links")
        report = recover_atom_links(reg, dry_run=not a.apply)
        report["unmatched"] = {"n": len(report["unmatched"]),
                               "geometry_ids": report["unmatched"][:50]}
        out = {"dry_run": not a.apply, "recover": report}
        if a.apply and a.connectivity:
            out["connectivity"] = backfill_connectivity(reg)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
