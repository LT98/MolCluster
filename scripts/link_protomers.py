"""Link every protomer in a registry to its parent, with a balanced deprotonation edge.

    python scripts/link_protomers.py --db data/mvp_ni_thq_cl.db
    python scripts/link_protomers.py --db data/mvp_ni_thq_cl.db --dry-run

The registry stores tHQ and tHQ(-1H) as unrelated rows. They are not unrelated, and the
reference scheme can price the difference — but only once the proton has somewhere to go,
so this first ensures the H2O / H3O+ couple exists with an energy at the same level of
theory, then writes `AH + n H2O -> A(n-) + n H3O+` for every pair it finds.

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

from mofsbu._types import Fidelity, MethodSpec                          # noqa: E402
from mofsbu.config import data_root, store_root                         # noqa: E402
from mofsbu.energy.protons import deprotonation_pairs, link_protomers   # noqa: E402
from mofsbu.energy.relax import relax_geometry                          # noqa: E402
from mofsbu.geometry.embed import embed_molecule, to_xyz                # noqa: E402
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles           # noqa: E402
from mofsbu.registry import (                                           # noqa: E402
    BlobStore, Registry, find, put_geometry, put_structure,
)

FF = MethodSpec(code="rdkit", code_version="2026.03", method="ETKDGv3+MMFF")

#: The couple. Charge is declared, not inferred: `[OH3+]` is only a proton carrier
#: because we say it carries one.
COUPLE = (("water", "O", 0), ("hydronium", "[OH3+]", 1))


def ensure_species(reg, name: str, smiles: str, charge: int, *, ml_model: str | None,
                   relax: bool) -> int:
    """Register the species if absent, and give it an energy if it has none."""
    mol = embed_molecule(mol_from_smiles(smiles), seed=7)
    graph = from_rdkit(mol, charge=charge, multiplicity=1, name=name)
    put = put_structure(reg, graph, tags=[name, "reference"])
    xyz_text = to_xyz(mol, name)
    put_geometry(reg, put.id, xyz_text, fidelity=Fidelity.FF, method=FF)
    if not relax:
        return put.id
    has_ml = reg.conn.execute(
        "SELECT 1 FROM geometries WHERE structure_id=? AND fidelity>=? AND energy IS NOT NULL",
        (put.id, int(Fidelity.ML))).fetchone()
    if has_ml:
        return put.id
    lines = xyz_text.splitlines()[2:]
    symbols = [line.split()[0] for line in lines if line.strip()]
    coords = [[float(x) for x in line.split()[1:4]] for line in lines if line.strip()]
    result = relax_geometry(coords, symbols, charge=charge, multiplicity=1,
                            target=Fidelity.ML, ml_model=ml_model)
    put_geometry(reg, put.id, result.to_xyz(name), fidelity=result.fidelity,
                 method=result.method, energy=result.energy,
                 converged=result.converged)
    return put.id


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

        ids = {name: ensure_species(reg, name, smiles, charge,
                                    ml_model=a.ml_model, relax=True)
               for name, smiles, charge in COUPLE}
        print(f"couple: water={ids['water']} hydronium={ids['hydronium']}")

        written = link_protomers(reg, water_id=ids["water"],
                                 hydronium_id=ids["hydronium"])
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
