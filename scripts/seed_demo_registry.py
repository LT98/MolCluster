"""Build `data/demo_registry.db` from REAL structures, through the real write surface.

Why this exists: the viewer needs content before the multi-centre placer (M6) can
produce any.  Why it looks like this: the first version of this script wrote rows with
raw INSERTs, and produced a database full of molecules that cannot exist — a
paddlewheel with no metal-metal bond, `muN` bridging on a two-metal node, one geometry
reused for four different structures, and every `typed_graph_hash` pointing at a blob
that was never stored.  Nothing caught it because nothing compared a column to its
graph.

So: every structure here is a real typed graph from `mofsbu.examples`, and every row is
written by `registry.api`.  The columns cannot disagree with the chemistry because
nothing computes them by hand.  `scripts/verify_registry.py` proves it.

Two things ARE synthetic, and are labelled so in the data rather than hidden:

* **coordinates** — a spring embedding of the graph (`spring-embed-3d`), not a relaxed
  structure.  Right atoms, right count, right connectivity, meaningless bond angles.
* **energies** — invented, recorded under the method `synthetic-demo` so that no one
  reading the viewer can mistake them for a calculation.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu import examples as ex                                    # noqa: E402
from mofsbu.config import data_root                                  # noqa: E402
from mofsbu.geometry.layout import to_xyz                            # noqa: E402
from mofsbu.graph import TypedGraph                                  # noqa: E402
from mofsbu.registry import (                                        # noqa: E402
    BlobStore, MethodSpec, Provenance, Registry, put_geometry, put_structure,
)
from mofsbu.registry.verify import summarise, verify                 # noqa: E402
from mofsbu._types import Fidelity                                   # noqa: E402

LAYOUT = MethodSpec(code="demo", code_version="1", method="spring-embed-3d")
SYNTH = MethodSpec(code="demo", code_version="1", method="synthetic-demo")


def catalogue() -> list[tuple[TypedGraph, list[str]]]:
    """40 distinct real structures, with tags."""
    out: list[tuple[TypedGraph, list[str]]] = []

    # paddlewheels — mu2, metal-metal bonded, one per metal plus charge variants
    for metal in ("Cu", "Zn", "Ni", "Co", "Fe", "Mn", "Cr", "Ru"):
        out.append((ex.paddlewheel(metal), ["sbu", "paddlewheel"]))
    out.append((ex.paddlewheel("Cu", charge=-2, multiplicity=3), ["sbu", "paddlewheel", "anionic"]))
    out.append((ex.paddlewheel("Cu", oxidation_state=1, spin_class="ls", charge=-2),
                ["sbu", "paddlewheel", "mixed"]))

    # mu3-oxo trimers — the mixed-valence case
    for valences in ((3, 3, 3), (2, 3, 3), (2, 2, 3)):
        out.append((ex.fe3_mu3_oxo(valences), ["sbu", "trimer"]))

    # mononuclear aqua complexes
    for metal, ox, spin, mult in (("Fe", 2, "hs", 5), ("Fe", 2, "ls", 1), ("Fe", 3, "hs", 6),
                                  ("Ni", 2, "hs", 3), ("Co", 2, "hs", 4), ("Zn", 2, "ls", 1),
                                  ("Mn", 2, "hs", 6), ("Cr", 3, "hs", 4)):
        out.append((ex.hexaaqua(metal, ox, spin, mult), ["mononuclear", "aqua"]))

    # mononuclear carboxylate + aqua — terminal/chelate binding
    for metal, n_aqua in (("Zn", 4), ("Cu", 4), ("Ni", 2), ("Co", 3), ("Mg", 4)):
        out.append((ex.aqua_carboxylate(metal, n_aqua=n_aqua), ["mononuclear", "carboxylate"]))

    # chelates — the fac/mer and Delta/Lambda cases L2 will have to separate (M5)
    for metal, ox, spin, n in (("Fe", 2, "ls", 3), ("Ru", 2, "ls", 3), ("Ni", 2, "hs", 2)):
        out.append((ex.metal_bipy(metal, ox, spin, n, multiplicity=1 if spin == "ls" else 3),
                    ["chelate", "bipy"]))

    # the discrimination pairs the identity work exists to draw
    out.append((ex.zn2_bridged_formates(), ["pair", "bridging"]))
    out.append((ex.zn2_chelated_formates(), ["pair", "chelating"]))
    out.append((ex.pt_ammine_dichloride(), ["pair", "cis-trans"]))

    # free ligands and organics
    for build, tags in ((ex.water, ["ligand"]), (ex.formate, ["ligand"]),
                        (ex.formic_acid, ["ligand", "protomer"]),
                        (lambda: ex.btc(True), ["linker"]),
                        (lambda: ex.btc(False), ["linker", "protomer"]),
                        (ex.bipyridine, ["linker"]),
                        (ex.cyclohexane, ["organic"]),
                        (ex.two_cyclopropanes, ["organic"])):
        out.append((build(), tags))
    return out



def reset_files(paths) -> None:
    """Clear stale database files, tolerating mounts that forbid deletion.

    Some mounts (sandboxes, network shares) allow writes but refuse unlink, and SQLite
    reports the consequences much later as a bare I/O error.  Truncating is equivalent
    for our purposes and always permitted.
    """
    for path in paths:
        try:
            path.unlink(missing_ok=True)
        except OSError:
            try:
                path.open("wb").close()
            except OSError:
                pass


def seed(db_path: Path, store_root: Path, n: int | None = None) -> None:
    """`n` truncates the catalogue; the default seeds all of it."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    reset_files(db_path.parent.glob(db_path.name + "*"))

    entries = catalogue()[:n] if n else catalogue()
    with Registry(db_path, BlobStore(store_root)) as reg:
        reg.migrate("demo registry")
        ids: list[int] = []
        for k, (g, tags) in enumerate(entries):
            put = put_structure(reg, g, tags=[*tags, "demo"])
            if not put.created:
                continue
            ids.append(put.id)

            # every third structure keeps no geometry, so the viewer's empty state is real
            if k % 7 == 3:
                continue

            put_geometry(reg, put.id, to_xyz(g, seed=k), fidelity=Fidelity.RAW, method=LAYOUT)
            if k % 3 == 0:
                put_geometry(reg, put.id, to_xyz(g, seed=k + 100), fidelity=Fidelity.XTB,
                             method=SYNTH, energy=-40.0 * len(g) + k, converged=True)
            if k % 5 == 0:
                put_geometry(reg, put.id, to_xyz(g, seed=k + 200), fidelity=Fidelity.DFT,
                             method=SYNTH, energy=-41.0 * len(g) + k, converged=k % 2 == 0)

        # provenance: pretend the later structures were assembled from earlier ones
        for k, sid in enumerate(ids):
            if k % 4 == 0 and k >= 8:
                put_reaction_edge(reg, sid, ids[k - 8], ids[k - 4], depth=1 + k % 3)

        # one geometry whose blob has been lost — a real failure mode, not a fake row
        orphan = reg.conn.execute(
            "SELECT id, coords_hash FROM geometries ORDER BY id DESC LIMIT 1").fetchone()
        orphaned = False
        if orphan is not None:
            try:
                reg.store.path(orphan["coords_hash"]).unlink(missing_ok=True)
                orphaned = True
            except OSError:
                pass          # some mounts forbid delete; the 404 test skips in that case

        n_struct, n_geom = reg.count("structures"), reg.count("geometries")
        problems = [p for p in verify(reg) if p.check != "coords_blob_missing"]

    print(f"seeded {db_path}: {n_struct} structures, {n_geom} geometries")
    print(f"  blobs in {store_root}")
    print(f"  integrity: {summarise(problems)}")
    if not orphaned:
        print("  NOTE: could not remove a blob to exercise the missing-coordinates path")
    print("  NOTE: coordinates are spring-embedding placeholders and energies are "
          "synthetic;\n        both are labelled as such in the `methods` table.")


def put_reaction_edge(reg, product_id: int, *reagent_ids: int, depth: int) -> None:
    from mofsbu.registry import put_reaction

    put_reaction(reg, product_id, Provenance(
        kind="assembly", reagent_ids=tuple(reagent_ids), depth=depth,
        note="synthetic demo route"))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--db", type=Path, default=data_root() / "demo_registry.db")
    p.add_argument("--store", type=Path, default=data_root() / "store")
    a = p.parse_args(argv)
    seed(a.db, a.store)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
