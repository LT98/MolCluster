"""Build the Zn / THQ / ammonium test system and register it.

THQ = tetrahydroxy-1,4-benzoquinone.  Four acidic hydroxyls, so its activation state is
an ENUMERATED BRANCH, not a choice someone makes in a script: 16 subsets, which the
identity layer collapses to seven distinct protomers, with the doubly-deprotonated state
splitting into ortho, meta and para.  Only the ortho protomer presents a dioxolene
pocket — two anionic oxygens on adjacent carbons — and that is the pocket a metal binds.
Meta and para have none, which is why the configuration matters and is not cosmetic.

Ammonium is the counterion: it never binds zinc, it carries the charge the linker sheds.
Deprotonation is recorded as a reaction edge with ammonia as the base, so the charge
accounting is visible as provenance rather than asserted.

Everything goes through the real pipeline: SMILES -> typed graph -> protomer enumeration
-> donor perception -> site frames -> frame-directed placement -> QC -> registry.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.config import data_root                                        # noqa: E402
from mofsbu.geometry.embed import embed_molecule, to_xyz                   # noqa: E402
from mofsbu.geometry.placer import (                                       # noqa: E402
    LigandPlacement, place_mononuclear, to_rdkit,
)
from mofsbu.graph._types import TypedGraph                                 # noqa: E402
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles              # noqa: E402
from mofsbu.naming import decompose                                        # noqa: E402
from mofsbu.registry import (                                              # noqa: E402
    BlobStore, MethodSpec, Provenance, Registry, alias_fragment, put_geometry, put_sites,
    put_structure, relabel_all,
)
from mofsbu.sites.frames import BindingMode                                # noqa: E402
from mofsbu.sites.model import chelate_pockets, find_pocket, perceive  # noqa: E402
from mofsbu.sites.protomers import enumerate_protomers                     # noqa: E402
from mofsbu._types import Fidelity                                         # noqa: E402

THQ_SMILES = "OC1=C(O)C(=O)C(O)=C(O)C1=O"
FF = MethodSpec(code="rdkit", code_version="2026.03", method="ETKDGv3+MMFF")
BUILD = MethodSpec(code="mofsbu", code_version="0.0.1", method="frame-directed-placement")
SEED = 7


def register_molecule(reg, mol, *, name, charge, multiplicity, tags, provenance=None):
    g = from_rdkit(mol, charge=charge, multiplicity=multiplicity, name=name)
    put = put_structure(reg, g, tags=tags, provenance=provenance)
    put_geometry(reg, put.id, to_xyz(mol, name), fidelity=Fidelity.FF, method=FF)
    put_sites(reg, put.id, perceive(mol))
    return put.id


def build(db_path: Path, store_root: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    with Registry(db_path, BlobStore(store_root)) as reg:
        reg.migrate("zn/thq test system")

        # ── the base and its conjugate acid: the modifier ───────────────────
        nh3 = embed_molecule(mol_from_smiles("N"), seed=3)
        nh4 = embed_molecule(mol_from_smiles("[NH4+]"), seed=3)
        nh3_id = register_molecule(reg, nh3, name="ammonia", charge=0, multiplicity=1,
                                   tags=["modifier", "base", "demo"])
        nh4_id = register_molecule(reg, nh4, name="ammonium", charge=1, multiplicity=1,
                                   tags=["modifier", "counterion", "demo"])
        w = embed_molecule(mol_from_smiles("O"), seed=3)
        w_id = register_molecule(reg, w, name="water", charge=0, multiplicity=1,
                                 tags=["ligand", "demo"])

        zn = TypedGraph(charge=2, multiplicity=1, name="Zn2+")
        zn.add_atom("Zn", oxidation_state=2, spin_class="ls")
        zn_id = put_structure(reg, zn, tags=["metal", "demo"]).id

        # ── every activation state of THQ, enumerated and deduplicated ──────
        protomers = enumerate_protomers(mol_from_smiles(THQ_SMILES))
        lines.append(f"THQ protonation states: 16 subsets -> {len(protomers)} distinct protomers")
        neutral_id = None
        by_config: dict[str, tuple[int, object]] = {}
        for proto in protomers:
            mol = embed_molecule(proto.mol, seed=SEED)
            sites = perceive(mol)
            kinds = sorted({pk.descriptor for pk in chelate_pockets(mol, sites)})
            name = f"THQ{proto.label if proto.n_deprotonated else ''}".strip()
            provenance = None
            if proto.n_deprotonated and neutral_id is not None:
                # THQ + n NH3 -> THQ(n-) + n NH4+ : the charge the linker sheds is carried
                provenance = Provenance(
                    kind="reaction", reagent_ids=(neutral_id, nh3_id), depth=1,
                    note=f"deprotonation x{proto.n_deprotonated} ({proto.configuration or 'n/a'})")
            sid = register_molecule(reg, mol, name=name, charge=proto.charge, multiplicity=1,
                                    tags=["linker", "thq", f"q{proto.charge}", "demo"],
                                    provenance=provenance)
            if proto.n_deprotonated == 0:
                neutral_id = sid
            # Name the FRAGMENT once, not each structure that contains it.  "THQ" is the
            # one author-supplied token in the whole run — it is the molecule the user
            # typed — and the suffix comes from the enumeration, not from a person.
            for frag in decompose(proto.graph):
                alias_fragment(reg, frag.l1, f"THQ{proto.label if proto.n_deprotonated else ''}"
                               .replace(" ", ""), source="build_zn_thq")
            by_config[proto.l1] = (sid, mol, proto)
            lines.append(
                f"  q{proto.charge:+d} {proto.label:16s} id={sid:<3d} "
                f"{len(proto.selections)} equivalent selection(s)  pockets: {', '.join(kinds)}")

        # ── coordination, on the pocket that actually binds ─────────────────
        def ligand_from(mol, predicate: dict, label: str) -> LigandPlacement:
            sites = perceive(mol)
            pocket = find_pocket(mol, sites, **predicate)
            if pocket is None:
                raise SystemExit(f"{label}: no pocket matching {predicate}")
            by_idx = {s.atom_idx: s for s in sites}
            return LigandPlacement(mol=mol, donor_idxs=tuple(pocket.donors),
                                   donor_types=tuple(by_idx[i].donor_type
                                                     for i in pocket.donors),
                                   mode=BindingMode.CHELATE, name=label)

        def aqua() -> LigandPlacement:
            return LigandPlacement(mol=embed_molecule(mol_from_smiles("O"), seed=3),
                                   donor_idxs=(0,), donor_types=("aqua_O",), name="H2O")

        # Select the ligand by what it CAN DO, never by what it is called.  The protomer
        # that chelates is the one carrying a dioxolene pocket; which descriptor that
        # turns out to be is a result, not an assumption.
        def pick(charge: int, **predicate):
            for sid, mol, proto in by_config.values():
                if proto.charge != charge:
                    continue
                if find_pocket(mol, perceive(mol), **predicate) is not None:
                    return sid, mol, proto
            raise SystemExit(f"no q{charge} protomer has a pocket matching {predicate}")

        CHELATING = {"ring_size": 5, "n_anionic": 2}       # both donors activated
        TRANSIENT = {"ring_size": 5, "n_anionic": 1}       # only one activated
        ortho_id, ortho_mol, ortho_proto = pick(-2, **CHELATING)
        mono_id, mono_mol, _mono_proto = pick(-1, **TRANSIENT)
        lines.append("")
        lines.append(f"chelating q-2 protomer selected structurally: "
                     f"{ortho_proto.label}  (of "
                     f"{len([p for p in protomers if p.charge == -2])} q-2 configurations, "
                     f"this is the one with a fully-activated 5-ring pocket)")

        complexes = [
            ("[Zn(THQ-2H_chelating)(H2O)2]", [ligand_from(ortho_mol, CHELATING, "THQ-2H")]
             + [aqua(), aqua()], "tetrahedral", 0, [zn_id, ortho_id, w_id]),
            ("[Zn(THQ-2H_chelating)(H2O)4]", [ligand_from(ortho_mol, CHELATING, "THQ-2H")]
             + [aqua() for _ in range(4)], "octahedral", 0, [zn_id, ortho_id, w_id]),
            ("[Zn(THQ-2H_chelating)2]2-", [ligand_from(ortho_mol, CHELATING, "THQ-2H"),
                                       ligand_from(ortho_mol, CHELATING, "THQ-2H")],
             "tetrahedral", -2, [zn_id, ortho_id]),
            ("[Zn(THQ-1H)(H2O)2]+", [ligand_from(mono_mol, TRANSIENT, "THQ-1H")]
             + [aqua(), aqua()], "tetrahedral", 1, [zn_id, mono_id, w_id]),
            ("[Zn(H2O)4]2+", [aqua() for _ in range(4)], "tetrahedral", 2, [zn_id, w_id]),
        ]

        for label, ligands, geometry, charge, reagents in complexes:
            result = place_mononuclear("Zn", ligands, geometry=geometry)
            mol = to_rdkit("Zn", ligands, result)
            g = from_rdkit(mol, charge=charge, multiplicity=1, name=label,
                           oxidation_states={0: 2}, spin_classes={0: "ls"})
            note = label + (" — transient, singly activated" if "1H" in label else "")
            put = put_structure(reg, g, tags=["complex", "zn", "demo"],
                                provenance=Provenance(kind="assembly", depth=2,
                                                      reagent_ids=tuple(reagents), note=note))
            put_geometry(reg, put.id, result.to_xyz(label), fidelity=Fidelity.RAW,
                         method=BUILD, choice_vector=result.choice_vector, seed=0,
                         qc=result.report.to_dict())
            put_sites(reg, put.id, perceive(mol))
            if charge < 0:                       # the counterion carries the balance
                from mofsbu.registry import put_reaction
                put_reaction(reg, put.id, Provenance(
                    kind="reaction", reagent_ids=(nh4_id,), depth=3,
                    note=f"{abs(charge)} x NH4+ counterion balances q{charge}"))
            (out_dir / f"{label}.xyz").write_text(result.to_xyz(label))
            pockets = ", ".join(sorted({lig.mode.value for lig in ligands}))
            lines.append(f"  {label:26s} id={put.id:<3d} q{charge:+d} CN="
                         f"{sum(l.denticity for l in ligands):<2d} {geometry:12s} {result.report}")

        relabel_all(reg)          # labels are a projection; refresh after aliasing
        lines.append("")
        lines.append("labels, composed from the graph + fragment aliases (never from a name):")
        for row in reg.conn.execute(
                "SELECT display_label FROM structures ORDER BY id"):
            lines.append(f"    {row['display_label']}")

        counts = (reg.count("structures"), reg.count("geometries"),
                  reg.count("site_catalog"), reg.count("reactions"))
    print("\n".join(lines))
    print(f"\n  registry: {counts[0]} structures, {counts[1]} geometries, "
          f"{counts[2]} catalogued sites, {counts[3]} provenance edges")
    print(f"  xyz written to {out_dir}")


def reset_files(paths) -> None:
    for path in paths:
        try:
            path.unlink(missing_ok=True)
        except OSError:
            try:
                path.open("wb").close()
            except OSError:
                pass


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--db", type=Path, default=data_root() / "zn_thq.db")
    p.add_argument("--store", type=Path, default=data_root() / "store")
    p.add_argument("--out", type=Path, default=ROOT / "data" / "zn_thq_xyz")
    a = p.parse_args(argv)
    reset_files(a.db.parent.glob(a.db.name + "*"))
    build(a.db, a.store, a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
