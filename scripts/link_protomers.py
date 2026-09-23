"""Link every protomer in a registry to its parent, with a balanced deprotonation edge.

    python scripts/link_protomers.py --db data/mvp_ni_thq_cl.db
    python scripts/link_protomers.py --db data/mvp_ni_thq_cl.db --dry-run

The registry stores tHQ and tHQ(-1H) as unrelated rows. They are not unrelated, and the
reference scheme can price the difference — but only once the proton has somewhere to go,
so this first ensures the H2O / H3O+ couple exists with an energy at the same level of
theory, then writes `AH + H2O -> A- + H3O+` for every pair exactly one proton apart.
A run with `pathways` does this itself when it finishes (`runner.finalise_run`); this
script is for a registry built before that.

Those edges are **isodesmic**: no metal-donor bond changes across the arrow, so unlike
every assembly edge in this project they survive `strict=True`.

Relaxing the couple needs the ML stack; everything else is bookkeeping.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu._types import Fidelity                                       # noqa: E402
from mofsbu.config import data_root, store_root                         # noqa: E402
from mofsbu.energy.protons import (                                     # noqa: E402
    deprotonation_pairs, ensure_couple, link_protomers,
)
from mofsbu.registry import BlobStore, Registry                         # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--db", type=Path, default=data_root() / "registry.db")
    p.add_argument("--store", type=Path, default=None)
    p.add_argument("--ml-model", default="mace_omol")
    p.add_argument("--dry-run", action="store_true",
                   help="report the pairs and write nothing")
    a = p.parse_args(argv)

    store = BlobStore(a.store or store_root())
    with Registry(a.db, store) as reg:
        pairs = deprotonation_pairs(reg)
        print(f"{a.db}: {len(pairs)} protomer pairs")
        if a.dry_run:
            for protonated, deprotonated, n in pairs[:20]:
                print(f"  x{n}  {protonated} -> {deprotonated}")
            return 0

        water, hydronium = ensure_couple(reg, target=Fidelity.ML, ml_model=a.ml_model)
        print(f"couple: water={water} hydronium={hydronium}")

        written = link_protomers(reg, water_id=water, hydronium_id=hydronium)
        reg.conn.commit()

    ok = [w for w in written if w["reaction_id"] is not None]
    print(f"edges written: {len(ok)} of {len(written)} attempted")
    refused: dict[str, int] = {}
    for w in written:
        if w["why_not"]:
            refused[w["why_not"][:70]] = refused.get(w["why_not"][:70], 0) + 1
    for reason, n in sorted(refused.items(), key=lambda kv: -kv[1])[:5]:
        print(f"  [{n}] {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
