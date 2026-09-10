"""Build structures from your own SMILES.  The general driver.

    # look before you build: protomers, donors, pockets — nothing is written
    python scripts/build.py --smiles "OC1=C(O)C(=O)C(O)=C(O)C1=O" --name THQ --dry-run

    # ligand only: enumerate activation states and register them
    python scripts/build.py --smiles "OC(=O)c1ccccc1O" --name salicylic --multiplicity 1

    # with a metal: enumerate coordination and place it
    python scripts/build.py --smiles "OC1=C(O)C(=O)C(O)=C(O)C1=O" --name THQ \\
        --metal Zn --oxidation-state 2 --spin ls --multiplicity 1 \\
        --coordination 4 6 --ligands-per-metal 1 2 --db data/my_run.db

Enumeration is deliberately liberal: every protomer x every chelate pocket or monodentate
donor x every requested coordination number and geometry is attempted, and IDENTITY does
the deduplication.  The run reports how many candidates collapsed onto how many distinct
structures, and every rejection is printed with its reason — a candidate that fails QC is
reported, not silently dropped, because "it built nothing" and "it built garbage" need to
look different.

Nothing is named by hand.  Configurations are orbit descriptors, labels are composed from
the graph, and the only author-supplied token is `--name`, which aliases the ligand
fragment for display.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.config import data_root                                        # noqa: E402
from mofsbu.energy.backends import (                                       # noqa: E402
    combined_multiplicity, spin_class_multiplicity,
)
from mofsbu.geometry.embed import embed_molecule, to_xyz                   # noqa: E402
from mofsbu.geometry.placer import (                                       # noqa: E402
    GEOMETRIES, LigandPlacement, place_mononuclear, to_rdkit,
)
from mofsbu.graph._types import TypedGraph                                 # noqa: E402
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles              # noqa: E402
from mofsbu.naming import decompose                                        # noqa: E402
from mofsbu.registry import (                                              # noqa: E402
    BlobStore, MethodSpec, Provenance, Registry, alias_fragment, put_geometry, put_sites,
    put_structure, relabel_all,
)
from mofsbu.sites.frames import BindingMode                                # noqa: E402
from mofsbu.sites.model import chelate_pockets, find_pockets, perceive      # noqa: E402
from mofsbu.sites.protomers import enumerate_protomers                     # noqa: E402
from mofsbu._types import Fidelity, MofsbuError                            # noqa: E402

FF = MethodSpec(code="rdkit", code_version="2026.03", method="ETKDGv3+MMFF")
BUILD = MethodSpec(code="mofsbu", code_version="0.0.1", method="frame-directed-placement")

# Coordination geometries to try for each coordination number.
GEOMETRY_BY_CN = {
    2: ["linear"], 3: ["trigonal"], 4: ["tetrahedral", "square_planar"],
    5: ["trigonal_bipyramidal", "square_pyramidal"], 6: ["octahedral"],
}


def describe(mol, label: str) -> list[str]:
    """What the perception layer sees.  Printed by --dry-run."""
    sites = perceive(mol)
    by_idx = {s.atom_idx: s for s in sites}
    lines = [f"    donors: " + (", ".join(f"{s.atom_idx}:{s.donor_type}" for s in sites) or "none")]
    pockets = chelate_pockets(mol, sites)
    if pockets:
        lines.append("    chelate pockets:")
        for pk in pockets:
            lines.append(f"      {pk.donors} {pk.descriptor}  [{pk.feasible_by}]")
    else:
        lines.append("    chelate pockets: none")
    return lines


def candidate_bindings(mol, *, modes: list[str]):
    """Every way this ligand could bind: chelate pockets and monodentate donors."""
    sites = perceive(mol)
    by_idx = {s.atom_idx: s for s in sites}
    out = []
    if "chelate" in modes:
        for pocket in chelate_pockets(mol, sites):
            out.append((tuple(pocket.donors),
                        tuple(by_idx[i].donor_type for i in pocket.donors),
                        BindingMode.CHELATE, pocket.descriptor))
    if "mono" in modes:
        for s in sites:
            out.append(((s.atom_idx,), (s.donor_type,), BindingMode.MONODENTATE, s.donor_type))
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("--smiles", required=True, help="the ligand, as SMILES")
    p.add_argument("--name", default="L", help="display alias for the ligand fragment")
    p.add_argument("--multiplicity", type=int,
                   help="spin multiplicity — required, never guessed (ground rule 5)")
    p.add_argument("--max-deprotonations", type=int, default=None,
                   help="cap the activation enumeration (default: all acidic protons)")
    p.add_argument("--metal", help="coordination centre; omit to build the ligand only")
    p.add_argument("--oxidation-state", type=int, default=2)
    p.add_argument("--spin", default="ls", help="per-centre spin class: hs | ls | is")
    p.add_argument("--coordination", type=int, nargs="+", default=[4, 6],
                   help="coordination numbers to attempt")
    p.add_argument("--geometry", nargs="+", default=None,
                   help=f"restrict geometries (default: all for each CN). {sorted(GEOMETRIES)}")
    p.add_argument("--ligands-per-metal", type=int, nargs="+", default=[1],
                   help="how many copies of the ligand to place")
    p.add_argument("--co-ligand", default="O", help="SMILES filling the remaining sites")
    p.add_argument("--binding", nargs="+", default=["chelate", "mono"],
                   choices=["chelate", "mono"])
    p.add_argument("--pocket-ring", type=int, default=None,
                   help="only chelate through pockets closing this ring size")
    p.add_argument("--pocket-anionic", type=int, default=None,
                   help="only chelate through pockets with this many anionic donors")
    p.add_argument("--db", type=Path, default=data_root() / "build.db")
    p.add_argument("--store", type=Path, default=data_root() / "store")
    p.add_argument("--out", type=Path, default=None, help="also write .xyz here")
    p.add_argument("--dry-run", action="store_true",
                   help="enumerate and report; write nothing")
    p.add_argument("--save-spec", type=Path, default=None,
                   help="write the resolved run spec as JSON, so the run is reproducible")
    a = p.parse_args(argv)

    if a.multiplicity is None and not a.dry_run:
        p.error("--multiplicity is required (it is never guessed); pass 1 for a closed shell")
    multiplicity = a.multiplicity or 1

    try:
        base = mol_from_smiles(a.smiles)
    except MofsbuError as exc:
        p.error(str(exc))

    protomers = enumerate_protomers(base, max_deprotonations=a.max_deprotonations,
                                    multiplicity=multiplicity)
    n_subsets = sum(len(x.selections) for x in protomers)
    print(f"ligand {a.name!r}  {a.smiles}")
    print(f"  activation states: {n_subsets} subsets -> {len(protomers)} distinct protomers")
    for proto in protomers:
        print(f"    q{proto.charge:+d}  {proto.label:34s} "
              f"{len(proto.selections)} equivalent selection(s)")
        if a.dry_run:
            for line in describe(embed_molecule(proto.mol, seed=7), proto.label):
                print(line)

    if a.save_spec:
        a.save_spec.write_text(json.dumps({
            k: (str(v) if isinstance(v, Path) else v) for k, v in vars(a).items()
        }, indent=2, sort_keys=True) + "\n")
        print(f"  spec written to {a.save_spec}")

    if a.dry_run:
        print("\n  --dry-run: nothing written")
        return 0

    attempted = built = duplicates = rejected = 0
    rejections: list[str] = []
    a.db.parent.mkdir(parents=True, exist_ok=True)

    with Registry(a.db, BlobStore(a.store)) as reg:
        reg.migrate(f"build {a.name}")
        co_mol = embed_molecule(mol_from_smiles(a.co_ligand), seed=3)
        co_sites = perceive(co_mol)
        co_donor = co_sites[0] if co_sites else None

        ligand_ids: dict[str, int] = {}
        for proto in protomers:
            mol = embed_molecule(proto.mol, seed=7)
            g = from_rdkit(mol, charge=proto.charge, multiplicity=multiplicity,
                           name=f"{a.name}{proto.label if proto.n_deprotonated else ''}")
            put = put_structure(reg, g, tags=[a.name, "ligand"])
            put_geometry(reg, put.id, to_xyz(mol, a.name), fidelity=Fidelity.FF, method=FF)
            put_sites(reg, put.id, perceive(mol))
            ligand_ids[proto.l1] = put.id
            for frag in decompose(g):
                alias_fragment(reg, frag.l1,
                               f"{a.name}{proto.label if proto.n_deprotonated else ''}"
                               .replace(" ", ""), source="build.py")

        if a.metal:
            # The bare ion's own multiplicity comes from its spin_class, not a literal
            # 1 — a Cu(II) "ls" ion is a doublet (d9, one unpaired electron); no
            # singlet is reachable, and MACE-OMOL-0's check_spin will say so.
            metal_multiplicity = spin_class_multiplicity(
                a.metal, a.oxidation_state, a.spin)
            metal_g = TypedGraph(charge=a.oxidation_state, multiplicity=metal_multiplicity,
                                 name=f"{a.metal}{a.oxidation_state:+d}")
            metal_g.add_atom(a.metal, oxidation_state=a.oxidation_state, spin_class=a.spin)
            metal_id = put_structure(reg, metal_g, tags=["metal"]).id

            print(f"\ncoordination candidates ({a.metal} {a.oxidation_state:+d}):")
            for proto in protomers:
                mol = embed_molecule(proto.mol, seed=7)
                bindings = candidate_bindings(mol, modes=a.binding)
                if a.pocket_ring or a.pocket_anionic is not None:
                    keep = find_pockets(mol, perceive(mol),
                                        ring_size=a.pocket_ring,
                                        n_anionic=a.pocket_anionic)
                    allowed = {p.donors for p in keep}
                    bindings = [b for b in bindings
                                if len(b[0]) == 1 or tuple(b[0]) in allowed]
                for donors, types, mode, kind in bindings:
                    for n_lig, cn in itertools.product(a.ligands_per_metal, a.coordination):
                        used = len(donors) * n_lig
                        n_co = cn - used
                        if n_co < 0 or (n_co and co_donor is None):
                            continue
                        geometries = a.geometry or GEOMETRY_BY_CN.get(cn, [])
                        for geometry in geometries:
                            attempted += 1
                            ligands = [
                                LigandPlacement(mol=mol, donor_idxs=donors, donor_types=types,
                                                mode=mode, name=a.name)
                                for _ in range(n_lig)
                            ] + [
                                LigandPlacement(mol=co_mol, donor_idxs=(co_donor.atom_idx,),
                                                donor_types=(co_donor.donor_type,),
                                                name=a.co_ligand)
                                for _ in range(n_co)
                            ]
                            tag = (f"q{proto.charge}/{kind}/{n_lig}x/CN{cn}/{geometry}")
                            try:
                                result = place_mononuclear(a.metal, ligands, geometry=geometry)
                            except (ValueError, MofsbuError) as exc:
                                rejected += 1
                                rejections.append(f"    {tag}: {exc}")
                                continue
                            if not result.ok:
                                rejected += 1
                                rejections.append(f"    {tag}: {result.report}")
                                continue
                            charge = a.oxidation_state + proto.charge * n_lig
                            complex_mol = to_rdkit(a.metal, ligands, result)
                            # Unpaired electrons add, multiplicities don't: the metal
                            # centre's own multiplicity combines with the ligand's
                            # rather than the ligand's overwriting it (that was the bug
                            # — a Cu(II) complex requested as a singlet).
                            complex_multiplicity = combined_multiplicity(
                                metal_multiplicity, multiplicity)
                            cg = from_rdkit(complex_mol, charge=charge,
                                            multiplicity=complex_multiplicity,
                                            oxidation_states={0: a.oxidation_state},
                                            spin_classes={0: a.spin})
                            cput = put_structure(
                                reg, cg, tags=[a.name, "complex", a.metal],
                                provenance=Provenance(
                                    kind="assembly", depth=1, note=tag,
                                    reagent_ids=(metal_id, ligand_ids[proto.l1])))
                            put_geometry(reg, cput.id, result.to_xyz(tag),
                                         fidelity=Fidelity.RAW, method=BUILD,
                                         choice_vector=result.choice_vector, seed=0,
                                         qc=result.report.to_dict())
                            put_sites(reg, cput.id, perceive(complex_mol))
                            if cput.created:
                                built += 1
                            else:
                                duplicates += 1
                            if a.out:
                                a.out.mkdir(parents=True, exist_ok=True)
                                (a.out / f"{cput.id:04d}_{tag.replace('/', '_')}.xyz"
                                 ).write_text(result.to_xyz(tag))

        relabel_all(reg)
        print()
        for row in reg.conn.execute(
                "SELECT display_label FROM v_structures ORDER BY n_metals, id"):
            print("   ", row["display_label"])
        counts = (reg.count("structures"), reg.count("geometries"))

    print(f"\n  attempted {attempted} coordination candidates: "
          f"{built} new, {duplicates} collapsed onto an existing structure, {rejected} rejected")
    if rejections:
        print("  rejections (reported, not hidden):")
        for line in rejections[:20]:
            print(line)
        if len(rejections) > 20:
            print(f"    ... and {len(rejections) - 20} more")
    print(f"  registry {a.db}: {counts[0]} structures, {counts[1]} geometries")
    if a.out:
        print(f"  xyz written to {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
