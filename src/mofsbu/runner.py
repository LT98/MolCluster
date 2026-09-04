"""Turn a spec into tasks, and execute tasks.  The two halves are deliberately separate.

`plan` writes tasks and returns; `work` claims and executes them.  Nothing here assumes
they happen in the same process, or on the same machine, or at the same time — which is
what lets one in-process worker on a laptop and N processes on the workstation (and, later,
an MPI launcher or a GPU worker) all be the same code.

Ground rule 8: `work` runs IN-PROCESS unless the machine is declared a workstation.  It
will not spawn anything on an unconfigured laptop.
"""
from __future__ import annotations

import itertools
import traceback
from pathlib import Path
from typing import Any

from mofsbu.config import max_workers, parallel_enabled
from mofsbu.geometry.embed import embed_molecule, to_xyz
from mofsbu.geometry.placer import GEOMETRIES, LigandPlacement, place_mononuclear, to_rdkit
from mofsbu.graph._types import TypedGraph
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles
from mofsbu.naming import decompose
from mofsbu.registry import (
    MethodSpec, Provenance, Registry, alias_fragment, put_geometry, put_sites, put_structure,
)
from mofsbu.registry.jobs import (
    add_task, claim_task, complete_task, create_run, fail_task, finish_run,
    set_diagnostics, task_counts,
)
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.model import chelate_pockets, find_pockets, perceive
from mofsbu.sites.protomers import enumerate_protomers
from mofsbu.spec import BuildSpec
from mofsbu._types import Fidelity, MofsbuError

FF = MethodSpec(code="rdkit", code_version="2026.03", method="ETKDGv3+MMFF")
BUILD = MethodSpec(code="mofsbu", code_version="0.0.1", method="frame-directed-placement")

GEOMETRY_BY_CN = {2: ["linear"], 3: ["trigonal"], 4: ["tetrahedral", "square_planar"],
                  5: ["trigonal_bipyramidal", "square_pyramidal"], 6: ["octahedral"]}


# ── planning ─────────────────────────────────────────────────────────────────

def plan(reg: Registry, spec: BuildSpec) -> tuple[int, int]:
    """Create the run and its tasks.  Returns (run_id, n_tasks)."""
    if spec.degree > 1:
        from mofsbu.assembly.join import grow

        grow(None, (), degree=spec.degree)      # raises NotBuiltYet, loudly and with why

    if spec.run_mode != "construct":
        from mofsbu._types import Fidelity
        from mofsbu.energy.relax import relax_geometry

        # Refuse before any work is queued rather than building structures and failing at
        # the optimisation step, which would leave a half-done run to interpret.
        relax_geometry(None, [], charge=0, multiplicity=1,
                       target=Fidelity.ML if spec.run_mode == "ml_go" else Fidelity.DFT)

    run_id = create_run(reg, spec)
    n = 0
    skipped: dict[str, dict[str, Any]] = {}

    def skip(reason: str, hint: str) -> None:
        entry = skipped.setdefault(reason, {"reason": reason, "hint": hint, "count": 0})
        entry["count"] += 1
    for mol_ix, molecule in enumerate(spec.molecules):
        base = mol_from_smiles(molecule.smiles)
        protomers = enumerate_protomers(base, max_deprotonations=molecule.max_deprotonations,
                                        multiplicity=molecule.multiplicity)
        for proto in protomers:
            n += 1
            add_task(reg, run_id, "ligand", {
                "molecule": mol_ix, "selection": list(proto.representative),
                "charge": proto.charge, "label": proto.label,
            }, priority=10)                      # ligands first: complexes reference them

        if not spec.metals:
            continue
        for proto, metal_ix in itertools.product(protomers, range(len(spec.metals))):
            mol = embed_molecule(proto.mol, seed=spec.seed or 7)
            sites = perceive(mol)
            bindings: list[tuple[tuple[int, ...], str]] = []
            if "chelate" in spec.binding:
                allowed = find_pockets(mol, sites, **spec.pocket.as_kwargs())
                bindings += [(tuple(p.donors), BindingMode.CHELATE.value) for p in allowed]
            if "mono" in spec.binding:
                bindings += [((s.atom_idx,), BindingMode.MONODENTATE.value) for s in sites]

            for donors, mode in bindings:
                for n_lig, cn in itertools.product(spec.ligands_per_metal, spec.coordination):
                    used = len(donors) * n_lig
                    n_co = cn - used
                    if n_co < 0:
                        skip(f"{used} donor sites exceed CN {cn}",
                             f"{n_lig} x {len(donors)}-dentate needs CN >= {used}")
                        continue
                    if n_co and not spec.co_ligand:
                        if not spec.allow_unsaturated:
                            skip(f"CN {cn} leaves {n_co} site(s) unfilled and no co-ligand is set",
                                 f"set a co-ligand (e.g. O for water), add {used} to the "
                                 f"coordination list, or enable allow_unsaturated")
                            continue
                        n_co, cn = 0, used      # build the unsaturated product instead
                    # Geometries must match the coordination number ACTUALLY being built.
                    # Taking them from the requested CN after collapsing to an unsaturated
                    # product asks for e.g. tetrahedral with two sites, which the placer
                    # rightly refuses — and a refusal there is a crash, not a chemistry
                    # answer, so it must not be reachable from a legal spec.
                    candidates = [g for g in (spec.geometries or GEOMETRY_BY_CN.get(cn, []))
                                  if g in GEOMETRY_BY_CN.get(cn, [])]
                    if not candidates:
                        skip(f"no coordination geometry for CN {cn}",
                             f"known CNs: {sorted(GEOMETRY_BY_CN)}"
                             + (f"; requested geometries {list(spec.geometries)} do not "
                                f"apply to CN {cn}" if spec.geometries else ""))
                        continue
                    for geometry in candidates:
                        n += 1
                        add_task(reg, run_id, "place", {
                            "molecule": mol_ix, "metal": metal_ix,
                            "selection": list(proto.representative),
                            "charge": proto.charge, "label": proto.label,
                            "donors": list(donors), "mode": mode, "n_ligands": n_lig,
                            "n_co": n_co, "cn": cn, "geometry": geometry,
                        })
    if not spec.metals:
        skip("no metal centres in the spec",
             "molecular-only construction: activation states and sites are produced; "
             "joining molecule to molecule is M5")
    set_diagnostics(reg, run_id, list(skipped.values()))
    return run_id, n


# ── execution ────────────────────────────────────────────────────────────────

def _protomer_mol(spec: BuildSpec, payload: dict[str, Any]):
    from mofsbu.sites.perception import deprotonate
    from mofsbu.sites.protomers import labile_sites

    molecule = spec.molecules[payload["molecule"]]
    base = mol_from_smiles(molecule.smiles)
    selection = payload["selection"]
    if selection:
        sites = labile_sites(base)
        base, _ = deprotonate(base, [sites[i] for i in selection])
    return embed_molecule(base, seed=spec.seed or 7), molecule


def execute(reg: Registry, task, spec: BuildSpec) -> tuple[int | None, int | None]:
    """Run one task.  Returns (structure_id, geometry_id)."""
    payload = task.payload
    mol, molecule = _protomer_mol(spec, payload)
    name = f"{molecule.name}{payload['label'] if payload['selection'] else ''}"

    if task.kind == "ligand":
        g = from_rdkit(mol, charge=payload["charge"], multiplicity=molecule.multiplicity,
                       name=name)
        put = put_structure(reg, g, tags=[molecule.name, "ligand"])
        geom = put_geometry(reg, put.id, to_xyz(mol, name), fidelity=Fidelity.FF, method=FF)
        put_sites(reg, put.id, perceive(mol))
        for frag in decompose(g):
            alias_fragment(reg, frag.l1, name.replace(" ", ""), source="runner")
        return put.id, geom.id

    if task.kind == "place":
        metal = spec.metals[payload["metal"]]
        sites = perceive(mol)
        by_idx = {s.atom_idx: s for s in sites}
        donors = tuple(payload["donors"])
        ligands = [LigandPlacement(
            mol=mol, donor_idxs=donors,
            donor_types=tuple(by_idx[i].donor_type for i in donors),
            mode=BindingMode(payload["mode"]), name=molecule.name)
            for _ in range(payload["n_ligands"])]
        if payload["n_co"]:
            co = embed_molecule(mol_from_smiles(spec.co_ligand), seed=3)
            co_site = perceive(co)[0]
            ligands += [LigandPlacement(mol=co, donor_idxs=(co_site.atom_idx,),
                                        donor_types=(co_site.donor_type,),
                                        name=spec.co_ligand)
                        for _ in range(payload["n_co"])]

        try:
            result = place_mononuclear(metal.symbol, ligands, geometry=payload["geometry"])
        except ValueError as exc:
            # The placer refusing a request is an ANSWER, not a breakage: a bidentate
            # ligand cannot span a linear two-coordinate centre, and saying so is the
            # correct behaviour.  Reporting it as `failed` would make a run full of sound
            # chemistry look like a run full of bugs.
            raise _Rejected(str(exc)) from exc
        if not result.ok:
            raise _Rejected(str(result.report))
        complex_mol = to_rdkit(metal.symbol, ligands, result)
        charge = metal.oxidation_state + payload["charge"] * payload["n_ligands"]
        g = from_rdkit(complex_mol, charge=charge, multiplicity=metal.multiplicity,
                       oxidation_states={0: metal.oxidation_state},
                       spin_classes={0: metal.spin_class}, name="")
        put = put_structure(reg, g, tags=[molecule.name, "complex", metal.symbol],
                            provenance=Provenance(kind="assembly", depth=1,
                                                  note=payload["geometry"]))
        geom = put_geometry(reg, put.id, result.to_xyz(name), fidelity=Fidelity.RAW,
                            method=BUILD, choice_vector=result.choice_vector,
                            seed=spec.seed, qc=result.report.to_dict())
        put_sites(reg, put.id, perceive(complex_mol))
        return put.id, geom.id

    raise MofsbuError(f"no executor for task kind {task.kind!r}")


class _Rejected(MofsbuError):
    """The task was well formed and the answer was no.  Not a failure of the machinery."""


def work(reg: Registry, spec: BuildSpec, run_id: int, *, limit: int | None = None) -> int:
    """Claim and execute tasks until none remain.  In-process."""
    done = 0
    while limit is None or done < limit:
        task = claim_task(reg, run_id)
        if task is None:
            break
        try:
            structure_id, geometry_id = execute(reg, task, spec)
            complete_task(reg, task.id, structure_id=structure_id, geometry_id=geometry_id)
        except _Rejected as exc:
            fail_task(reg, task.id, str(exc), rejected=True)
        except Exception as exc:                                   # noqa: BLE001
            fail_task(reg, task.id, f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}")
        done += 1
    return done


def _worker_process(db_path: str, store_path: str, spec_json: str, run_id: int) -> None:
    from mofsbu.registry import BlobStore

    with Registry(Path(db_path), BlobStore(Path(store_path))) as reg:
        work(reg, BuildSpec.from_json(spec_json), run_id)


def run(reg: Registry, spec: BuildSpec, *, workers: int | None = None) -> dict[str, Any]:
    """Plan and execute a spec.  Serial on a laptop; parallel only where declared."""
    run_id, n_tasks = plan(reg, spec)
    reg.conn.commit()
    n_workers = workers if workers is not None else max_workers()

    if n_workers <= 1 or not parallel_enabled():
        work(reg, spec, run_id)
    else:
        import multiprocessing as mp

        db, store = str(reg.db_path), str(reg.store.root)
        procs = [mp.Process(target=_worker_process, args=(db, store, spec.to_json(None), run_id))
                 for _ in range(n_workers)]
        for proc in procs:
            proc.start()
        for proc in procs:
            proc.join()

    from mofsbu.registry import relabel_all

    relabel_all(reg)          # labels are a projection; refresh once the fragments exist
    status = finish_run(reg, run_id)
    return {"run_id": run_id, "tasks": n_tasks, "status": status,
            "counts": task_counts(reg, run_id), "workers": n_workers}
