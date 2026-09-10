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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mofsbu.config import compute_device, device_note, max_workers, parallel_enabled
from mofsbu.energy.backends import combined_multiplicity, spin_class_multiplicity
from mofsbu.geometry.embed import embed_molecule, embed_with_report, to_xyz
from mofsbu.geometry.placer import GEOMETRIES, LigandPlacement, place_mononuclear, to_rdkit
from mofsbu.graph._types import TypedGraph
from mofsbu.graph.from_mol import from_rdkit, mol_from_smiles
from mofsbu.naming import decompose
from mofsbu.registry import (
    MethodSpec, Provenance, Registry, alias_fragment, catalog_drift, put_geometry,
    put_site_state, put_sites, put_structure,
)
from mofsbu.registry.jobs import (
    add_task, cancel_requested, claim_task, complete_task, create_run, fail_task,
    finish_run,
    outcome_summary, set_diagnostics, task_counts,
)
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.model import find_pockets, perceive, shifting_pocket_donors
from mofsbu.sites.protomers import enumerate_protomers
from mofsbu.sites.state import refresh_state
from mofsbu.spec import BuildSpec
from mofsbu._types import EnergyBackendUnavailable, Fidelity, MofsbuError

@dataclass(frozen=True)
class Outcome:
    """What a task produced, and whether any of it was NEW.

    `structure_created=False` means the registry already had this identity: the task
    succeeded and wrote nothing.  Under D2 that is the expected outcome for a large part
    of any enumeration, and it has to be visible or a run of pure duplicates looks
    exactly like a run of discoveries.
    """

    structure_id: int | None
    geometry_id: int | None
    structure_created: bool | None = None
    geometry_created: bool | None = None
    detail: dict[str, Any] = field(default_factory=dict)


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
        from mofsbu.energy.relax import mode_status

        # Refuse before any work is queued rather than building structures and failing at
        # the optimisation step, which would leave a half-done run to interpret.  Asking
        # the backend whether it is installed is not the same question as whether the
        # code exists, and the two get different answers on the laptop.
        status = mode_status(spec.ml_model).get(spec.run_mode)
        if status is None or not status["available"]:
            note = (status or {}).get("note") or f"run mode {spec.run_mode!r} cannot execute"
            # An unwired mode is a missing BODY, not a missing install, and the two get
            # different exceptions (ground rule 8).  Getting this backwards is what let a
            # 500-structure ml_go run report success while never leaving RAW.
            if status is not None and status.get("backend_available") and not status.get("wired"):
                from mofsbu.assembly.join import NotBuiltYet

                raise NotBuiltYet(f"{spec.run_mode}: {note}")
            raise EnergyBackendUnavailable(f"{spec.run_mode}: {note}")

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
    mol, embed_report = embed_with_report(base, seed=spec.seed or 7)
    return mol, molecule, embed_report


RELAX_PRIORITY = -10        # constructs drain first; relaxes are the expensive tail


def resolve_ml_model(spec: BuildSpec) -> str:
    """The ML model this run means, as a name, resolved once.

    `BuildSpec.ml_model` may be None ("use what this machine declares").  Resolving it
    at queue time and storing the ANSWER is the difference between a run that is
    reproducible and one that depends on an environment variable at the moment a worker
    happened to claim a task.
    """
    from mofsbu.config import resolve_ml_backend

    return resolve_ml_backend(spec.ml_model)


def relax_method_spec(reg: Registry, structure_id: int, target: Fidelity,
                      solvent: str | None = None, ml_model: str | None = None):
    """The MethodSpec a relaxation WOULD produce, without running it.

    This is what makes "has this already been computed?" answerable before spending the
    compute rather than after.  The backend can describe itself — code, pinned version,
    solvent, charge, multiplicity, device — from the structure's own charge and spin,
    which is everything the `methods` row is keyed on.
    """
    from mofsbu.energy.relax import backend_for

    row = reg.conn.execute(
        "SELECT net_charge, multiplicity FROM structures WHERE id=?",
        (structure_id,)).fetchone()
    if row is None:
        raise MofsbuError(f"no structure {structure_id}")
    return backend_for(target, ml_model=ml_model).method_spec(
        charge=int(row["net_charge"]), multiplicity=int(row["multiplicity"]),
        solvent=solvent)


def existing_relaxation(reg: Registry, source_geometry_id: int, target: Fidelity,
                        structure_id: int, solvent: str | None = None,
                        ml_model: str | None = None) -> int | None:
    """The geometry this exact relaxation already produced, or None.

    Keyed on (source geometry, method row).  Deliberately NOT on the structure alone:
    two raw constructs of one identity are two different starting points and may relax
    into two different minima, and collapsing them here would silently discard the
    second one.  What this catches is the genuine duplicate — the same starting
    geometry, at the same level of theory, relaxed again, which is what a re-run of an
    unchanged spec produces for every structure it already has.
    """
    from mofsbu.registry import find_method_id

    try:
        spec = relax_method_spec(reg, structure_id, target, solvent, ml_model)
    except Exception:                                                   # noqa: BLE001
        # Cannot describe the method (backend gone?).  Say "unknown", not "no": the
        # executor re-checks, and queueing a task that turns out redundant is cheap
        # while skipping one that was needed is not.
        return None
    mid = find_method_id(reg, spec)
    if mid is None:
        return None
    # `fidelity` is in the key as well as `method_id`, which looks redundant and is not:
    # two rungs served by backends that describe themselves identically would otherwise
    # collide, and a request for an xTB relaxation would be answered with an ML one that
    # happened to share a method row.  A test asserts exactly this.
    #
    # The converse now also matters: ONE rung served by two backends.  An MP-0 and an
    # OMOL-0 relaxation of the same construct are both fidelity=ML and are different
    # method rows, so this correctly reports "not done yet" for the second — a re-run
    # that switches model recomputes instead of handing back the other model's answer.
    row = reg.conn.execute(
        "SELECT id FROM geometries WHERE relaxed_from=? AND method_id=? AND "
        "structure_id=? AND fidelity=? AND energy IS NOT NULL LIMIT 1",
        (source_geometry_id, mid, structure_id, int(target))).fetchone()
    return None if row is None else int(row["id"])


def read_xyz(text: str) -> tuple[list[str], list[list[float]]]:
    """Parse a stored .xyz back into symbols and coordinates."""
    lines = text.strip().splitlines()
    n = int(lines[0])
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2:2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if len(symbols) != n:
        raise MofsbuError(f"xyz claims {n} atoms and carries {len(symbols)}")
    return symbols, coords


def _execute_relax(reg: Registry, task, spec: BuildSpec) -> Outcome:
    """Optimise a stored geometry and store the RESULT AS A SECOND ROW.

    The raw construct is never replaced.  `geometries.relaxed_from` chains the two, which
    is the fidelity ladder of §6.2: one identity, several realisations, and the cheap one
    recorded as having correctly pointed at the better one.  It also means a relax that
    fails costs nothing already built — the construct is still there and only this task
    is marked failed.

    KNOWN GAP, stated rather than hidden: the relaxed coordinates are stored against the
    SAME structure identity without re-deriving the typed graph from them.  A relaxation
    that breaks or forms a bond should become a new L3 with a `derived_from` edge (D11,
    §6.2), and detecting that needs a coordinates -> graph path this package does not
    have yet.  Until it does, a relax that tears a node apart will be filed under the
    identity of the node it destroyed.  `qc` on the stored row is the partial guard.
    """
    from mofsbu.energy.relax import relax_geometry
    from mofsbu.registry import geometry_xyz, get_graph

    payload = task.payload
    structure_id = int(payload["structure_id"])
    source_gid = int(payload["geometry_id"])
    target = Fidelity(int(payload["target"]))

    graph = get_graph(reg, structure_id)
    if graph.charge is None or graph.multiplicity is None:
        raise MofsbuError(
            f"structure {structure_id} has no recorded charge/multiplicity; a relaxation "
            f"cannot be set up without both, and guessing either is how the archived runs "
            f"ended up with a spin convention nobody could reconstruct")

    # Checked again here, not only at queue time: between queueing and claiming, another
    # worker may have done this exact relaxation.  The check is one indexed SELECT and
    # the thing it guards is an hour of xTB, so it is worth doing twice.
    # The MODEL comes from the task payload, not from the environment, and not from the
    # spec re-read at execution time.  A queue filled on the workstation and drained on
    # the laptop must produce the same theory on both, and `MOFSBU_ML_MODEL` is per
    # machine.  `spec.ml_model` is the fallback for tasks queued before this field
    # existed; None then resolves to the declared default, as it always did.
    ml_model = payload.get("ml_model", spec.ml_model)
    already = existing_relaxation(reg, source_gid, target, structure_id,
                                  payload.get("solvent"), ml_model)
    if already is not None:
        return Outcome(structure_id, already, None, False,
                       {"relax": {"skipped": "already relaxed at this level of theory",
                                  "from_geometry": source_gid, "existing_geometry": already,
                                  "target": target.name}})

    symbols, coords = read_xyz(geometry_xyz(reg, source_gid))
    result = relax_geometry(coords, symbols, charge=int(graph.charge),
                            multiplicity=int(graph.multiplicity), target=target,
                            solvent=payload.get("solvent"), ml_model=ml_model)

    geom = put_geometry(
        reg, structure_id, result.to_xyz(graph.name or ""), fidelity=result.fidelity,
        method=result.method, energy=result.energy, converged=result.converged,
        relaxed_from=source_gid)
    detail = {"relax": {"from_geometry": source_gid, "target": target.name,
                        "converged": result.converged, "steps": result.n_steps,
                        "fmax": round(result.fmax, 4), "energy": result.energy,
                        "relaxation_energy": result.relaxation_energy,
                        "device": compute_device(),
                        "method": result.method.describe()}}
    return Outcome(structure_id, geom.id, None, geom.created, detail)


def _record_sites(reg: Registry, structure_id: int, geometry_id: int, mol,
                  graph: TypedGraph, fidelity: Fidelity) -> dict[str, Any]:
    """Perceive once, then derive state from that one perception.

    Both halves of D5 are written here, in this order, deliberately: `refresh_state` is
    handed the very list that went into `site_catalog`, so the two tiers cannot disagree
    about which atoms are donors.  Perception happening anywhere else in a build is the
    failure mode `tests/test_sites_state.py::test_perception_runs_once_per_structure`
    exists to catch.
    """
    sites = perceive(mol)
    drift = catalog_drift(reg, structure_id, sites)
    put_sites(reg, structure_id, sites)
    conf = mol.GetConformer()
    coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
               conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    pocket_donors = shifting_pocket_donors(mol, sites)
    states = refresh_state(sites, coords, graph=graph, symbols=symbols,
                           pocket_donors=pocket_donors, fidelity=fidelity)
    n_stored = put_site_state(reg, structure_id, geometry_id, states, fidelity=fidelity)
    report: dict[str, Any] = {
        "n_sites": len(sites),
        "n_open": sum(1 for s in states if s.is_open),
        "n_provisional": sum(1 for s in states
                             if s.ease is not None and s.ease.provisional)}
    if drift:
        # Reported, never swallowed.  See `registry.api.catalog_drift`.  The case this
        # used to fire on — two routes to one identity rendering the same chemistry as
        # different resonance forms — is closed in perception; anything that reaches here
        # now is a fresh disagreement about which atoms are donors, and it is still a
        # finding rather than a task failure, because the build is fine and the first
        # catalog stands.
        report["catalog_drift"] = drift
        report["n_state_dropped"] = len(states) - n_stored
    return report


def execute(reg: Registry, task, spec: BuildSpec) -> Outcome:
    """Run one task.  Returns what it produced AND what it took to produce it."""
    payload = task.payload
    if task.kind == "relax":
        # Handled first: a relax payload carries a structure and a geometry, not a
        # molecule index, so `_protomer_mol` has nothing to work with.
        return _execute_relax(reg, task, spec)

    mol, molecule, embed_report = _protomer_mol(spec, payload)
    detail: dict[str, Any] = {"embed": embed_report}
    name = f"{molecule.name}{payload['label'] if payload['selection'] else ''}"

    if task.kind == "ligand":
        g = from_rdkit(mol, charge=payload["charge"], multiplicity=molecule.multiplicity,
                       name=name)
        put = put_structure(reg, g, tags=[molecule.name, "ligand"])
        geom = put_geometry(reg, put.id, to_xyz(mol, name), fidelity=Fidelity.FF, method=FF)
        detail["sites"] = _record_sites(reg, put.id, geom.id, mol, g, Fidelity.FF)
        for frag in decompose(g):
            alias_fragment(reg, frag.l1, name.replace(" ", ""), source="runner")
        return Outcome(put.id, geom.id, put.created, geom.created, detail)

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
            co, co_embed = embed_with_report(mol_from_smiles(spec.co_ligand), seed=3)
            detail["co_ligand_embed"] = co_embed
            co_sites = perceive(co)
            if not co_sites:
                # This used to be `perceive(co)[0]` and an IndexError, reported as a
                # crash with a traceback — for what is a plain statement about the
                # co-ligand: nothing in it can bind.  `[SiH4]` and `CC` reach here, and
                # so does any donor type perception does not yet cover.
                raise _Rejected(
                    f"co-ligand {spec.co_ligand!r} has no perceivable donor atom, so it "
                    f"cannot occupy a coordination site",
                    code="co_ligand_no_donor",
                    detail={"co_ligand": spec.co_ligand,
                            "hint": "give the co-ligand as the species that actually "
                                    "binds — '[I-]' rather than 'I2', 'O' for water, "
                                    "'[OH-]' for hydroxide — or add its donor type to "
                                    "sites.perception"})
            co_site = co_sites[0]
            detail["co_ligand"] = {"smiles": spec.co_ligand,
                                   "donor_type": co_site.donor_type,
                                   "donor_element":
                                       co.GetAtomWithIdx(co_site.atom_idx).GetSymbol(),
                                   "n": payload["n_co"]}
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
            raise _Rejected(str(exc), code="placer_refused",
                            detail={**detail, "geometry": payload["geometry"],
                                    "cn": payload["cn"]}) from exc
        detail["distances"] = [x.to_dict() for x in result.donor_distances]
        if not result.ok:
            # The QC report is stored STRUCTURED, not just stringified: which atoms,
            # which elements, how far inside which limit, and where each target M-L
            # distance came from.  That is the difference between "2 clash(es), closest
            # 1.40 A" — which is what a whole run used to collapse into — and a finding
            # you can group, sort and act on.
            raise _Rejected(str(result.report), code=result.report.code,
                            detail={**detail, "qc": result.report.to_dict(),
                                    "geometry": payload["geometry"],
                                    "cn": payload["cn"]})
        complex_mol = to_rdkit(metal.symbol, ligands, result)
        charge = metal.oxidation_state + payload["charge"] * payload["n_ligands"]
        # The complex's multiplicity is the metal centre's own (from its spin_class,
        # not a stale literal) combined with the ligand's — unpaired electrons add,
        # multiplicities don't.  See `spin_class_multiplicity`/`combined_multiplicity`.
        metal_multiplicity = spin_class_multiplicity(
            metal.symbol, metal.oxidation_state, metal.spin_class)
        multiplicity = combined_multiplicity(metal_multiplicity, molecule.multiplicity)
        g = from_rdkit(complex_mol, charge=charge, multiplicity=multiplicity,
                       oxidation_states={0: metal.oxidation_state},
                       spin_classes={0: metal.spin_class}, name="")
        put = put_structure(reg, g, tags=[molecule.name, "complex", metal.symbol],
                            provenance=Provenance(kind="assembly", depth=1,
                                                  note=payload["geometry"]))
        geom = put_geometry(reg, put.id, result.to_xyz(name), fidelity=Fidelity.RAW,
                            method=BUILD, choice_vector=result.choice_vector,
                            seed=spec.seed, qc=result.report.to_dict())
        detail["sites"] = _record_sites(reg, put.id, geom.id, complex_mol, g, Fidelity.RAW)
        detail["qc"] = result.report.to_dict()
        return Outcome(put.id, geom.id, put.created, geom.created, detail)

    raise MofsbuError(f"no executor for task kind {task.kind!r}")


class _Rejected(MofsbuError):
    """The task was well formed and the answer was no.  Not a failure of the machinery.

    Carries a `code` for grouping and a `detail` dict for diagnosis.  A rejection whose
    only content is a sentence is a rejection nobody can count.
    """

    def __init__(self, message: str, *, code: str = "rejected",
                 detail: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.code = code
        self.detail = detail or {}


def queue_relax(reg: Registry, spec: BuildSpec, task, out: "Outcome") -> int | None:
    """Emit the follow-up relaxation for a task that just built something.

    A separate task rather than a step inside the build, because relaxation is the
    expensive, resumable half: it retries without rebuilding, it survives a closed
    laptop, and it lets one machine drain a queue another machine filled — which is the
    entire point of the two-machine setup.  Lower priority than the constructs, so the
    cheap enumeration finishes first and the queue that remains is all GPU work.

    Ligands are relaxed too, not just complexes.  `energy.reference` refuses to subtract
    energies computed at different levels of theory, so a relaxed complex and an
    unrelaxed free ligand could never appear in the same equation — which would leave the
    reference scheme with nothing to say.
    """
    from mofsbu.energy.relax import MODE_FIDELITY

    if spec.run_mode == "construct" or task.kind == "relax":
        return None
    if out.structure_id is None or out.geometry_id is None:
        return None
    target = MODE_FIDELITY.get(spec.run_mode)
    if target is None:
        return None
    if existing_relaxation(reg, out.geometry_id, target, out.structure_id,
                           None, spec.ml_model) is not None:
        # Re-running an unchanged spec rebuilds the same constructs, recognises them by
        # identity (D2), and hands back the geometry ids it already had.  Without this,
        # every one of them was queued for relaxation again -- the same starting
        # geometry, the same level of theory, the same answer, at full price.
        return None
    return add_task(reg, task.run_id, "relax", {
        "structure_id": out.structure_id, "geometry_id": out.geometry_id,
        "target": int(target), "solvent": None,
        # Resolved at QUEUE time and carried in the payload, so the theory is fixed by
        # the run rather than by whichever machine happens to claim the task.
        "ml_model": resolve_ml_model(spec) if target is Fidelity.ML else None,
    }, priority=RELAX_PRIORITY)


def work(reg: Registry, spec: BuildSpec, run_id: int, *, limit: int | None = None) -> int:
    """Claim and execute tasks until none remain.  In-process."""
    done = 0
    while limit is None or done < limit:
        if cancel_requested(reg, run_id):
            # Cooperative: the task in hand has already finished, and nothing new is
            # claimed.  Everything computed so far is in the registry.
            break
        task = claim_task(reg, run_id)
        if task is None:
            break
        try:
            out = execute(reg, task, spec)
            complete_task(reg, task.id, structure_id=out.structure_id,
                          geometry_id=out.geometry_id,
                          structure_created=out.structure_created,
                          geometry_created=out.geometry_created,
                          detail=out.detail)
            queue_relax(reg, spec, task, out)
        except _Rejected as exc:
            fail_task(reg, task.id, str(exc), rejected=True,
                      code=exc.code, detail=exc.detail)
        except Exception as exc:                                   # noqa: BLE001
            # An unexpected exception is a bug, and the type is the most useful thing to
            # group by — twenty tasks dying of one IndexError is one problem, not twenty.
            fail_task(reg, task.id, f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
                      code=type(exc).__name__,
                      detail={"traceback": traceback.format_exc()[-4000:],
                              "task_kind": task.kind, "payload": task.payload})
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
    if spec.run_mode != "construct":
        # Printed before any work starts: an accelerator sitting idle for a whole run is
        # otherwise only visible in nvidia-smi, an hour later, by accident.
        print(f"mofsbu run {run_id}: mode={spec.run_mode}  {device_note()}  "
              f"workers={n_workers}")

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
            "counts": task_counts(reg, run_id), "workers": n_workers,
            "summary": outcome_summary(reg, run_id)}
