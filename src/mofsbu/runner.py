"""Turn a spec into tasks, and execute tasks.  The two halves are deliberately separate.

`plan` writes tasks and returns; `work` claims and executes them.  Nothing here assumes
they happen in the same process, or on the same machine, or at the same time — which is
what lets one in-process worker on a laptop and N processes on the workstation (and, later,
an MPI launcher or a GPU worker) all be the same code.

Ground rule 8: `work` runs IN-PROCESS unless the machine is declared a workstation.  It
will not spawn anything on an unconfigured laptop.

A run holds two kinds of work and they do not want the same hardware.  Construction
(`BUILD_KINDS`) is CPU work that scales with cores; `relax` is the optimiser, which on a
declared accelerator is one device's worth of work however many processes ask for it.
`execute_run` therefore divides the declared workers between the two queues rather than
letting one pool claim either — see `plan_workers`.
"""
from __future__ import annotations

import itertools
import json
import time
import traceback
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Sequence

from rdkit import Chem

from mofsbu.config import (
    compute_device, device_note, max_workers, parallel_enabled, relax_workers,
)
from mofsbu.energy.backends import combined_multiplicity, spin_class_multiplicity
from mofsbu.geometry.embed import embed_molecule, embed_with_report, to_xyz
from mofsbu.geometry import qc as qc_mod
from mofsbu.geometry.placer import (
    GEOMETRIES, LigandPlacement, cis_vertices, place_mononuclear, to_rdkit)
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
    outcome_summary, set_diagnostics, sweep_interrupted, task_counts, touch_run,
)
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.model import (
    find_pockets, perceive, shifting_pocket_donors, vacancy_sites,
)
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

#: Hard ceiling on compositions per (metal, n_ligands, CN).  A guard rail, not a policy:
#: `max_distinct_ligands` is the knob a user turns, and this only stops a spec that turns
#: it too far from filling a queue with tens of thousands of tasks before anyone notices.
#: Hitting it is reported through the run diagnostics rather than silently truncating.
MAX_COMPOSITIONS = 4000


def _compositions(kinds: list[dict[str, Any]], total: int,
                  max_distinct: int) -> list[list[dict[str, Any]]]:
    """Every multiset of `total` ligand pieces drawn from `kinds`, at most `max_distinct`
    of them different.

    A "kind" is one (molecule, protomer, donor set, binding mode).  Two protomers of one
    molecule are therefore two kinds, which is deliberate — a partially deprotonated set
    on one centre is real chemistry for a polyprotic linker, and refusing to enumerate it
    would be a chemistry decision smuggled in as a data-structure limit.

    `max_distinct=1` returns exactly the homoleptic compositions, which is what the
    planner did before this function existed, so an unchanged spec plans an unchanged run.
    """
    out: list[list[dict[str, Any]]] = []
    for d in range(1, min(max_distinct, total) + 1):
        for chosen in itertools.combinations(range(len(kinds)), d):
            # Ordered compositions of `total` into exactly `d` positive parts.
            for cuts in itertools.combinations(range(1, total), d - 1):
                bounds = (0,) + cuts + (total,)
                counts = [bounds[i + 1] - bounds[i] for i in range(d)]
                out.append([{**kinds[k], "count": c}
                            for k, c in zip(chosen, counts) if c])
                if len(out) > MAX_COMPOSITIONS:
                    return out[:MAX_COMPOSITIONS]
    return out


#: Ceiling on the intermediates `pathways` adds.  Same kind of guard rail as
#: `MAX_COMPOSITIONS`: the ladder of a homoleptic sweep is short, the sub-multiset lattice
#: of a mixed one is not, and hitting the ceiling is reported through the diagnostics
#: rather than silently truncating the chain.
MAX_PATHWAY_TASKS = 2000

#: A step runs after every rung is built (`place` is 0) and before the relaxations: it
#: builds nothing the place tasks do not build, and what it adds is the EDGE.
PATHWAY_PRIORITY = -5


@dataclass(frozen=True)
class PlannedTask:
    """One task, before it has an id.

    `task_refs` names payload keys whose value is an INDEX into the plan rather than a
    task id — the only thing enumeration cannot know, because ids exist once the rows do.
    `plan` rewrites them on the way in, and every reference points BACKWARD, so a task's
    parent is always already written when its turn comes.
    """

    kind: str
    payload: dict[str, Any]
    priority: int = 0
    task_refs: tuple[str, ...] = ()


@dataclass(frozen=True)
class PlannedRun:
    """What a spec would queue, and why it would not queue more.

    Enumeration is separated from writing so that the page can ASK — "how big is this
    run?" is a question you want answered while still editing the spec, and answering it
    by planning a throwaway run would put rows in the registry for a question.  `plan` and
    `estimate` are then the same enumeration, which is the only way the number the page
    shows can be trusted to be the number you get.
    """

    tasks: tuple[PlannedTask, ...] = ()
    diagnostics: tuple[dict[str, Any], ...] = ()
    n_kinds: int = 0

    def by_kind(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for task in self.tasks:
            out[task.kind] = out.get(task.kind, 0) + 1
        return out


def _payload_key(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True)


def _used_vertices(components: Sequence[dict[str, Any]]) -> int:
    return sum(int(c["denticity"]) * int(c["count"]) for c in components)


def _parent_payloads(payload: dict[str, Any], co_ligand: str | None = None,
                     ) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    """`(parent payload, the component one copy of which the join adds)`, one per kind.

    The parent is this product with one ligand copy removed and the vertices that copy
    occupied left EMPTY — not filled with a co-ligand.  That is what makes the step an
    addition rather than a substitution, and it is the only version of the parent a
    `join` can actually turn into this product.

    `co_ligand` (set under `co_ligand_counts="range"`, D26) adds one more parent: this
    product with one co-ligand fewer, its vertex left empty.  Listed last, so the ligand
    parents keep the order a `fill` spec plans them in.
    """
    out = []
    for i, comp in enumerate(payload["components"]):
        rest = [dict(c) for j, c in enumerate(payload["components"]) if j != i]
        if int(comp["count"]) > 1:
            rest.insert(i, dict(comp, count=int(comp["count"]) - 1))
        parent = dict(payload, components=rest)
        parent["n_vacant"] = int(payload["cn"]) - int(payload["n_co"]) - _used_vertices(rest)
        out.append((parent, dict(comp, count=1)))
    if co_ligand and int(payload["n_co"]) > 0:
        parent = dict(payload, n_co=int(payload["n_co"]) - 1,
                      n_vacant=int(payload.get("n_vacant", 0)) + 1)
        out.append((parent, _co_component(co_ligand)))
    return out


#: Why a charged co-ligand is not stepped (see `_charged`).
CO_LIGAND_CHARGED = ("B20: `place` leaves a co-ligand's charge out of the complex's net "
                     "charge, so a join-built co-ligand step would land on a different "
                     "identity from the rung it should reach")


def _charged(smiles: str) -> bool:
    """A co-ligand carrying a formal charge — which `place` does not count (B20)."""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None and Chem.GetFormalCharge(mol) != 0


def _co_component(smiles: str) -> dict[str, Any]:
    """The co-ligand as the component a `grow` step adds: one monodentate copy."""
    return {"co_ligand": smiles, "denticity": 1, "count": 1,
            "mode": BindingMode.MONODENTATE.value}


def enumerate_plan(spec: BuildSpec) -> PlannedRun:
    """Every task this spec implies, with no registry in sight.

    Pure enumeration: the same function answers "plan this run" and "how big would this
    run be", so the estimate on the page cannot drift from what submitting actually does.
    """
    tasks: list[PlannedTask] = []
    skipped: dict[str, dict[str, Any]] = {}

    def skip(reason: str, hint: str, **extra: Any) -> None:
        entry = skipped.setdefault(reason, {"reason": reason, "hint": hint, "count": 0,
                                            **extra})
        entry["count"] += 1

    kinds: list[dict[str, Any]] = []
    ligand_task: dict[tuple[int, tuple[int, ...]], int] = {}
    for mol_ix, molecule in enumerate(spec.molecules):
        base = mol_from_smiles(molecule.smiles)
        protomers = enumerate_protomers(base, max_deprotonations=molecule.max_deprotonations,
                                        multiplicity=molecule.multiplicity)
        for proto in protomers:
            ligand_task[(mol_ix, tuple(proto.representative))] = len(tasks)
            tasks.append(PlannedTask("ligand", {
                "molecule": mol_ix, "selection": list(proto.representative),
                "charge": proto.charge, "label": proto.label,
            }, priority=10))                     # ligands first: complexes reference them

        if not spec.metals:
            continue
        # Every way this molecule could occupy vertices, flattened into ONE list across
        # all molecules.  That flattening is the whole change: the enumeration used to
        # run per molecule, so a coordination sphere could only ever hold copies of one
        # of them and [Mg(dtBK)(Cl)] was not expressible by any spec.
        for proto in protomers:
            mol = embed_molecule(proto.mol, seed=spec.seed or 7)
            sites = perceive(mol)
            bindings: list[tuple[tuple[int, ...], str]] = []
            if "chelate" in spec.binding:
                allowed = find_pockets(mol, sites, **spec.pocket.as_kwargs())
                bindings += [(tuple(p.donors), BindingMode.CHELATE.value) for p in allowed]
            if "mono" in spec.binding:
                bindings += [((s.atom_idx,), BindingMode.MONODENTATE.value) for s in sites]
            for donors, mode in bindings:
                kinds.append({
                    "molecule": mol_ix, "selection": list(proto.representative),
                    "charge": proto.charge, "label": proto.label,
                    "donors": list(donors), "mode": mode,
                    "denticity": len(donors),
                })

    # Under `range` the co-ligand is a reagent a step consumes, so it is built and
    # relaxed like a free ligand — without it a co-ligand step balances one water short.
    vary_co = bool(spec.metals and spec.co_ligand and spec.co_ligand_counts == "range")
    if vary_co:
        ligand_task[("co_ligand", spec.co_ligand)] = len(tasks)
        tasks.append(PlannedTask("co_ligand", {"smiles": spec.co_ligand}, priority=10))

    place_at: dict[str, int] = {}
    place = _placer(tasks, place_at)

    requested: list[int] = []
    for metal_ix in range(len(spec.metals)):
        for total_n, cn in itertools.product(spec.ligands_per_metal, spec.coordination):
            compositions = _compositions(kinds, total_n, spec.max_distinct_ligands)
            if not compositions and kinds:
                skip(f"no ligand composition of {total_n} piece(s) from "
                     f"{len(kinds)} kind(s)",
                     "every combination exceeded max_distinct_ligands "
                     f"({spec.max_distinct_ligands})")
                continue
            for combo in compositions:
                used = _used_vertices(combo)
                n_co = cn - used
                if n_co < 0:
                    skip(f"{used} donor sites exceed CN {cn}",
                         f"this composition needs CN >= {used}")
                    continue
                n_vacant = 0
                if n_co and not spec.co_ligand:
                    if not spec.allow_unsaturated:
                        skip(f"CN {cn} leaves {n_co} site(s) unfilled and no co-ligand is set",
                             f"set a co-ligand (e.g. O for water), add {used} to the "
                             f"coordination list, or enable allow_unsaturated")
                        continue
                    # The requested CN is KEPT and the surplus vertices are left
                    # empty.  This used to be `n_co, cn = 0, used` — rebuild it at
                    # whatever CN the ligands could fill — and that is a different
                    # molecule: two ketones on a tetrahedral Mg is not a linear
                    # 2-coordinate Mg, and only the first is what the next assembly
                    # step attaches to.  `allow_unsaturated`'s own docstring warned
                    # against "quietly building something smaller than you asked
                    # for", which is exactly what the collapse did.
                    n_co, n_vacant = 0, n_co
                # Geometries come from the requested CN, which is now also the CN
                # actually built.
                candidates = [g for g in (spec.geometries or GEOMETRY_BY_CN.get(cn, []))
                              if g in GEOMETRY_BY_CN.get(cn, [])]
                if not candidates:
                    skip(f"no coordination geometry for CN {cn}",
                         f"known CNs: {sorted(GEOMETRY_BY_CN)}"
                         + (f"; requested geometries {list(spec.geometries)} do not "
                            f"apply to CN {cn}" if spec.geometries else ""))
                    continue
                # `fill` builds the one count that covers every free vertex; `range`
                # (D26) builds each count down to `co_ligand_window` fewer, the vertices
                # left uncovered empty in the same polyhedron — CN is never collapsed.
                lowest = max(0, n_co - spec.co_ligand_window)
                co_counts = range(n_co, lowest - 1, -1) if vary_co else (n_co,)
                for geometry in candidates:
                    for k in co_counts:
                        requested.append(place({
                            "metal": metal_ix, "components": combo,
                            "n_co": k, "cn": cn, "geometry": geometry,
                            "n_vacant": n_vacant + n_co - k,
                        }))

    if spec.pathways and spec.metals:
        step_co = spec.co_ligand if vary_co else None
        if step_co and _charged(step_co):
            skip(f"co-ligand steps not planned: {step_co!r} is charged",
                 CO_LIGAND_CHARGED + ". The lower counts are still built; they are not "
                 "linked to each other by co-ligand steps")
            step_co = None
        if vary_co and spec.co_ligand_window == 0:
            skip("co_ligand_window 0 plans no pathway step",
                 "every step adds onto an empty vertex, and a window of 0 allows no rung "
                 "one; widen the window, or use co_ligand_counts fill for the ligand ladder")
        walk = {"co_ligand": step_co,
                "max_vacant": spec.co_ligand_window if vary_co else None}
        shadow, shadow_at = list(tasks), dict(place_at)
        if _plan_pathways(tasks, place, requested, ligand_task, cap=MAX_PATHWAY_TASKS,
                          **walk):
            # The cap bit.  The same walk run uncapped on a copy says how much was cut, so
            # the refusal carries a size and not only the fact of one; nothing is queued
            # from the copy.
            _plan_pathways(shadow, _placer(shadow, shadow_at), requested, ligand_task,
                           cap=None, **walk)
            skip(f"pathway intermediates stopped at {MAX_PATHWAY_TASKS} tasks",
                 "the ladder below this product is not planned; narrow the "
                 "composition (max_distinct_ligands) or the ligand-count range"
                 + ("; lower co_ligand_window, or set co_ligand_counts to fill"
                    if vary_co else ""),
                 cap=MAX_PATHWAY_TASKS, uncapped_tasks=len(shadow))
    if not spec.metals:
        skip("no metal centres in the spec",
             "molecular-only construction: activation states and sites are produced; "
             "joining molecule to molecule is M5")
    return PlannedRun(tuple(tasks), tuple(skipped.values()), len(kinds))


def _placer(tasks: list[PlannedTask], place_at: dict[str, int]):
    def place(payload: dict[str, Any]) -> int:
        """Queue a coordination sphere once, however many ladders reach it."""
        key = _payload_key(payload)
        index = place_at.get(key)
        if index is None:
            index = place_at[key] = len(tasks)
            tasks.append(PlannedTask("place", payload))
        return index
    return place


def _plan_pathways(tasks: list[PlannedTask], place, requested: list[int],
                   ligand_task: dict[tuple[Any, Any], int], *,
                   co_ligand: str | None = None, max_vacant: int | None = None,
                   cap: int | None) -> bool:
    """Add the intermediates each product is reached from, and the steps between them.

    Walked breadth first from the products the spec asked for, one ligand copy at a time,
    down to the centre carrying none.  Every rung is a real, coordinatively unsaturated
    species — the thing `allow_unsaturated` names — and the `grow` task between two rungs
    is the addition itself: it re-places the parent, joins one ligand onto the vertices
    that copy will occupy, and stores the product with the parent as a reagent.

    The chain is walked even when the spec asked for only one count, because "what is this
    complex assembled from" is the same question at every rung.

    Under `co_ligand_counts="range"` (D26) `max_vacant` is the window: no rung with more
    empty vertices is planned, so the walk stops at the lowest co-ligand state inside it
    rather than at the bare centre.  `co_ligand` set also steps the co-ligand, which is
    what joins the ligand-free rungs of one metal and CN into one chain.

    Returns True when `cap` stopped the walk before it reached the bottom.
    """
    seen: set[tuple[int, int]] = set()
    frontier = list(requested)
    while frontier:
        nxt: list[int] = []
        for child_ix in frontier:
            child = tasks[child_ix].payload
            for parent_payload, component in _parent_payloads(child, co_ligand):
                if max_vacant is not None and parent_payload["n_vacant"] > max_vacant:
                    continue
                if cap is not None and len(tasks) >= cap:
                    return True
                parent_ix = place(parent_payload)
                if (parent_ix, child_ix) in seen:
                    continue
                seen.add((parent_ix, child_ix))
                nxt.append(parent_ix)
                tasks.append(PlannedTask("grow", {
                    "parent_task": parent_ix, "child_task": child_ix,
                    "ligand_task": ligand_task.get(
                        ("co_ligand", component["co_ligand"]) if "co_ligand" in component
                        else (component["molecule"], tuple(component["selection"]))),
                    "component": component, "metal": child["metal"],
                    "cn": child["cn"], "geometry": child["geometry"],
                }, priority=PATHWAY_PRIORITY,
                    task_refs=("parent_task", "child_task", "ligand_task")))
        frontier = [ix for ix in nxt if tasks[ix].payload["components"]
                    or (co_ligand and tasks[ix].payload["n_co"])]
    return False


def estimate(spec: BuildSpec) -> dict[str, Any]:
    """How big this run would be, without queueing any of it.

    The same enumeration `plan` uses, so the number is the number — an estimate computed
    by a second, simpler formula is a number that goes wrong exactly when a spec gets
    interesting, which is when it is being read.

    A spec this machine could not execute is reported rather than raised: the page asks
    this question WHILE the spec is being edited, and half-finished is the normal state of
    the thing being measured.
    """
    try:
        _refuse_unrunnable(spec)
    except MofsbuError as exc:
        return {"ok": False, "reason": str(exc), "refusal": type(exc).__name__}
    planned = enumerate_plan(spec)
    by_kind = planned.by_kind()
    builds = sum(by_kind.get(k, 0) for k in BUILD_KINDS)
    # One relaxation per thing built, at most: `queue_relax` skips a structure this
    # registry has already relaxed at this level of theory, and this function has no
    # registry to ask.  Stated as a ceiling rather than quietly counted as certain.
    relaxations = builds if spec.run_mode != "construct" else 0
    return {
        "ok": True, "tasks": len(planned.tasks), "by_kind": by_kind,
        "relaxations": relaxations, "total": len(planned.tasks) + relaxations,
        "ligand_kinds": planned.n_kinds,
        "diagnostics": [dict(d) for d in planned.diagnostics],
        "run_mode": spec.run_mode, "pathways": spec.pathways,
        "co_ligand_counts": spec.co_ligand_counts,
        "co_ligand_window": spec.co_ligand_window,
        # A cap that bit, stated on its own rather than left as one diagnostic among
        # several: the run it describes is incomplete, not merely smaller.
        "capped": [dict(d) for d in planned.diagnostics if "cap" in d],
    }


def _refuse_unrunnable(spec: BuildSpec) -> None:
    """Refuse a spec this build or this machine cannot execute, BEFORE anything is queued.

    Both refusals are about the run as a whole, so they belong ahead of enumeration:
    building half a run and failing at the optimisation step leaves a half-done run to
    interpret, which is worse than a refusal with a reason.
    """
    if spec.degree > 1:
        from mofsbu._types import NotBuiltYet

        # This used to CALL `grow(None, (), degree=...)` purely to borrow the exception it
        # raised.  That was fine while grow was a stub and became a lie the moment M5/S4
        # built it: the tripwire would have thrown an AttributeError on the None seed
        # instead of explaining anything.  The refusal is stated directly now, and it
        # refuses something different from what it used to — not "growth is not written"
        # but "the RUN PIPELINE does not drive it yet".
        raise NotBuiltYet(
            f"degree {spec.degree}: the branch tree and `assembly.grow` are built (M5/S4), "
            "but planning a multi-step construction as TASKS is not — that is S4.1, where "
            "the enumerator becomes what the builder drives. Call "
            "`assembly.construct.enumerate_constructions` directly in the meantime. A run "
            "that quietly queued degree-1 work instead would be the wrong answer wearing "
            "the right count.")

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


def plan(reg: Registry, spec: BuildSpec) -> tuple[int, int]:
    """Create the run and its tasks.  Returns (run_id, n_tasks)."""
    _refuse_unrunnable(spec)
    planned = enumerate_plan(spec)
    run_id = create_run(reg, spec)
    ids: list[int] = []
    for task in planned.tasks:
        payload = task.payload
        if task.task_refs:
            # Plan-time indices become task ids here.  Every reference points backward,
            # so the row it names is already written; a forward one would silently store
            # an id that does not exist yet.
            payload = dict(payload)
            for key in task.task_refs:
                index = payload.get(key)
                payload[key] = ids[index] if isinstance(index, int) else None
        ids.append(add_task(reg, run_id, task.kind, payload, priority=task.priority))
    set_diagnostics(reg, run_id, list(planned.diagnostics))
    return run_id, len(ids)


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


def _components(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """The ligand pieces of a `place` payload, old shape or new.

    A payload used to name ONE molecule and a count (`molecule`/`n_ligands`); it now
    carries a `components` list so a coordination sphere can hold more than one ligand
    kind.  Both are read, because a queue filled by an earlier build may still be
    draining — tasks outlive the process that wrote them, which is the entire premise of
    `plan` and `work` being separate.
    """
    if "components" in payload:
        return payload["components"]
    return [{"molecule": payload["molecule"], "selection": payload["selection"],
             "charge": payload["charge"], "label": payload["label"],
             "donors": payload["donors"], "mode": payload["mode"],
             "denticity": len(payload["donors"]),
             "count": payload["n_ligands"]}]


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

    # A relaxation MOVES ATOMS, so whatever QC said about the construct is a statement
    # about coordinates that no longer exist.  Storing the result unchecked is how a
    # geometry built despite a marginal clash would end up in the registry looking exactly
    # like one that never had a problem — and after `MARGINAL_OVERLAP` this is no longer
    # hypothetical, it is the intended path for hundreds of structures.
    after = _clash_qc(graph, result.symbols, result.positions)
    geom = put_geometry(
        reg, structure_id, result.to_xyz(graph.name or ""), fidelity=result.fidelity,
        method=result.method, energy=result.energy, converged=result.converged,
        relaxed_from=source_gid, qc=after.to_dict())
    detail = {"relax": {"from_geometry": source_gid, "target": target.name,
                        "converged": result.converged, "steps": result.n_steps,
                        "fmax": round(result.fmax, 4), "energy": result.energy,
                        "relaxation_energy": result.relaxation_energy,
                        "device": compute_device(),
                        "method": result.method.describe()}}
    if not after.ok:
        # Stored, flagged, and NOT rejected.  The optimisation was paid for and its
        # result is a real answer — "these ligands still do not fit once relaxed" is a
        # finding about the chemistry, and rev 21's lesson was that discarding compute
        # after it has been spent is the expensive mistake.  What it must never do is
        # look clean: the QC report is on the geometry row.
        detail["qc_failed_after_relax"] = after.to_dict()
    return Outcome(structure_id, geom.id, None, geom.created, detail)


def _clash_qc(graph: TypedGraph, symbols, positions) -> qc_mod.QCReport:
    """Clash check on coordinates nothing measured a bond length against.

    Clashes only, deliberately, and for the same reason at both call sites.  After a
    relaxation the bond lengths are the optimiser's answer rather than a target that was
    missed, and flagging them would report a converged minimum as a construction defect.
    After a join they were SET from `geometry.distances` — the donor is placed at its
    metal-donor distance by construction — so a bond-length check there asks a question
    whose answer the placement already guaranteed, while the overlaps between the incoming
    ligand and the ones already on the centre are the real question and are what this
    reports.
    """
    import numpy as np

    bonded = {(min(i, j), max(i, j)) for i, j, _ in graph.edges()}
    report = qc_mod.QCReport()
    report.clashes = qc_mod.check_clashes(list(symbols), np.asarray(positions), bonded)
    report.ok = not report.clashes
    return report


def _record_sites(reg: Registry, structure_id: int, geometry_id: int, mol,
                  graph: TypedGraph, fidelity: Fidelity,
                  vacancies: tuple = (), metal_idx: int = 0) -> dict[str, Any]:
    """Perceive once, then derive state from that one perception.

    Both halves of D5 are written here, in this order, deliberately: `refresh_state` is
    handed the very list that went into `site_catalog`, so the two tiers cannot disagree
    about which atoms are donors.  Perception happening anywhere else in a build is the
    failure mode `tests/test_sites_state.py::test_perception_runs_once_per_structure`
    exists to catch.

    `vacancies` are the coordination vertices the placer left empty
    (`PlacementResult.vacancies`).  They are NOT perceived — perception looks at a
    molecule and finds donors, and an empty vertex is not in the molecule.  They come from
    the construction, which is the only thing that knows the polyhedron was bigger than
    the ligand set, and they join the same site list so `open_sites()` returns both kinds.
    """
    sites = perceive(mol)
    drift = catalog_drift(reg, structure_id, sites)
    conf = mol.GetConformer()
    coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
               conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    if vacancies:
        sites = sites + vacancy_sites(metal_idx, coords[metal_idx], vacancies)
    put_sites(reg, structure_id, sites)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    pocket_donors = shifting_pocket_donors(mol, [s for s in sites if not s.is_vacancy])
    states = refresh_state(sites, coords, graph=graph, symbols=symbols,
                           pocket_donors=pocket_donors, fidelity=fidelity)
    n_stored = put_site_state(reg, structure_id, geometry_id, states, fidelity=fidelity)
    report: dict[str, Any] = {
        "n_sites": len(sites),
        "n_vacancies": sum(1 for s in sites if s.is_vacancy),
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


@dataclass(frozen=True)
class _Built:
    """One placed coordination sphere, before anything is stored.

    Shared by the two tasks that need one.  `place` stores what comes out of here;
    `grow` needs the PARENT rung as something to join onto, and re-places it rather
    than reading its coordinates back — construction is a deterministic function of the
    payload and the seed (D13), so the rebuilt parent IS the stored parent, and the
    product of joining onto it is therefore the product `place` builds for the next rung:
    one identity, two routes (D2), which is the whole point of recording the step.
    """

    metal: Any
    mol: Any
    graph: TypedGraph
    result: Any
    name: str
    molecule_names: list[str]
    components: list[dict[str, Any]]
    charge: int
    multiplicity: int
    detail: dict[str, Any]


def _build_sphere(spec: BuildSpec, payload: dict[str, Any]) -> _Built:
    """Place one coordination sphere from a `place` payload.  Writes nothing."""
    metal = spec.metals[payload["metal"]]
    detail = {}
    components = _components(payload)
    # One embed per COMPONENT, and every copy of a component shares that one
    # conformer.  Two copies of a ligand are two placements of the same molecule, not
    # two independent embeddings of it — re-embedding each copy would make a
    # homoleptic complex's ligands silently non-identical and put stochastic
    # conformer noise inside a single structure.
    ligands: list[LigandPlacement] = []
    parts: list[str] = []
    ligand_charge = 0
    multiplicities: list[int] = []
    for comp in components:
        cmol, cmolecule, cembed = _protomer_mol(spec, comp)
        detail.setdefault("embed", cembed)
        csites = perceive(cmol)
        by_idx = {s.atom_idx: s for s in csites}
        cdonors = tuple(comp["donors"])
        missing = [i for i in cdonors if i not in by_idx]
        if missing:
            # The donor indices were chosen at PLAN time against a conformer embedded
            # then; if perception no longer sees them the two halves have diverged and
            # placing atom 9 because the payload says 9 would bind whatever now
            # happens to sit there.
            raise _Rejected(
                f"donor index/indices {missing} are no longer perceived on "
                f"{cmolecule.name}; the plan and this executor disagree about the "
                f"molecule", code="donor_index_stale",
                detail={"component": comp, "perceived": sorted(by_idx)})
        cname = f"{cmolecule.name}{comp['label'] if comp['selection'] else ''}"
        parts.append(f"{cname}x{comp['count']}" if comp["count"] > 1 else cname)
        ligand_charge += int(comp["charge"]) * int(comp["count"])
        multiplicities += [cmolecule.multiplicity] * int(comp["count"])
        ligands += [LigandPlacement(
            mol=cmol, donor_idxs=cdonors,
            donor_types=tuple(by_idx[i].donor_type for i in cdonors),
            mode=BindingMode(comp["mode"]), name=cname)
            for _ in range(int(comp["count"]))]
    name = "+".join(parts)
    molecule_names = sorted({spec.molecules[c["molecule"]].name for c in components})
    if payload["n_co"]:
        co, co_embed, co_site = _co_ligand_mol(spec.co_ligand)
        detail["co_ligand_embed"] = co_embed
        detail["co_ligand"] = {"smiles": spec.co_ligand,
                               "donor_type": co_site.donor_type,
                               "donor_element":
                                   co.GetAtomWithIdx(co_site.atom_idx).GetSymbol(),
                               "n": payload["n_co"]}
        ligands += [LigandPlacement(mol=co, donor_idxs=(co_site.atom_idx,),
                                    donor_types=(co_site.donor_type,),
                                    name=spec.co_ligand)
                    for _ in range(payload["n_co"])]

    # A rung's vacancies are what the NEXT rung binds to, so their arrangement is this
    # step's business: left to fill order they come back trans, and a ~90 deg chelate is
    # then refused with `chelate_cannot_span` (measured: strain 2.156 trans against 0.094
    # cis, on one centre and one ligand).  Below two there is no arrangement to choose,
    # and those builds are byte-identical to a `placement 2` one.
    cn_asked = payload.get("cn")
    n_vacant = 0 if cn_asked is None else int(cn_asked) - sum(l.denticity for l in ligands)
    reserve = (cis_vertices(payload["geometry"], int(cn_asked), n_vacant)
               if n_vacant >= 2 else None)
    try:
        # `cn` is passed explicitly so an unsaturated centre keeps the polyhedron it
        # was asked for and reports its empty vertices, instead of being silently
        # rebuilt as a smaller, differently-shaped complex.
        result = place_mononuclear(metal.symbol, ligands,
                                   geometry=payload["geometry"],
                                   cn=cn_asked, reserve=reserve)
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
        # A NEAR MISS is not a refusal.  A rigid placement that lands a hydrogen
        # 0.06 A inside its limit has not made a chemical mistake — it has made a
        # geometric one that an optimiser undoes in a few steps, and throwing the
        # construct away means never finding that out.  So a marginal report is
        # built and stored, and the relaxation this run is going to do anyway
        # becomes the thing that decides.  `geometry.qc.MARGINAL_OVERLAP` documents
        # where the threshold came from.
        #
        # Two guards on that leniency.  It only applies when there IS a relaxation
        # coming: under `construct` there is no optimiser to appeal to, so the
        # construct stays rejected and the code says a GO mode would have retried
        # it.  And a marginal geometry is never stored looking clean — its QC
        # report travels with it and `_execute_relax` re-checks it afterwards.
        if not (result.report.marginal and spec.run_mode != "construct"):
            # The QC report is stored STRUCTURED, not just stringified: which atoms,
            # which elements, how far inside which limit, and where each target M-L
            # distance came from.  That is the difference between "2 clash(es),
            # closest 1.40 A" — which is what a whole run used to collapse into —
            # and a finding you can group, sort and act on.
            code = result.report.code
            extra = {}
            if result.report.marginal:
                code = "qc_clash_marginal_no_go"
                extra = {"hint": "this is a near miss, within "
                                 f"{qc_mod.MARGINAL_OVERLAP} A. A run mode that "
                                 "relaxes (ml_go, xtb_go) would have built it and "
                                 "let the optimiser try to resolve the overlap."}
            raise _Rejected(str(result.report), code=code,
                            detail={**detail, "qc": result.report.to_dict(),
                                    **extra,
                                    "geometry": payload["geometry"],
                                    "cn": payload["cn"]})
        detail["built_despite_qc"] = {
            "reason": "marginal clash, retried under a relaxation",
            "worst_overlap": result.report.to_dict()["worst_overlap"],
            "threshold": qc_mod.MARGINAL_OVERLAP}
    complex_mol = to_rdkit(metal.symbol, ligands, result)
    # Summed over components, so a mixed sphere gets the charge it actually carries:
    # [Mg(dtBK)(Cl)] is +1, not the +2 a neutral-ligand assumption would give or the
    # 0 that two chlorides would.  Each component contributes `charge x count`.
    charge = metal.oxidation_state + ligand_charge
    # The complex's multiplicity is the metal centre's own (from its spin_class,
    # not a stale literal) combined with EVERY ligand's — unpaired electrons add,
    # multiplicities don't.  See `spin_class_multiplicity`/`combined_multiplicity`.
    metal_multiplicity = spin_class_multiplicity(
        metal.symbol, metal.oxidation_state, metal.spin_class)
    multiplicity = combined_multiplicity(metal_multiplicity, *multiplicities)
    g = from_rdkit(complex_mol, charge=charge, multiplicity=multiplicity,
                   oxidation_states={0: metal.oxidation_state},
                   spin_classes={0: metal.spin_class}, name="")
    detail["composition"] = {
        "ligands": [{"molecule": spec.molecules[c["molecule"]].name,
                     "label": c["label"], "charge": c["charge"],
                     "mode": c["mode"], "count": c["count"]} for c in components],
        "n_distinct": len(components), "ligand_charge": ligand_charge,
        "n_co": payload["n_co"], "n_vacant": payload.get("n_vacant", 0)}

    return _Built(metal=metal, mol=complex_mol, graph=g, result=result, name=name,
                  molecule_names=molecule_names, components=components, charge=charge,
                  multiplicity=multiplicity, detail=detail)


def _co_ligand_mol(smiles: str):
    """The co-ligand embedded once, and the donor site every copy of it binds through.

    `place`, a co-ligand `grow` step and the free co-ligand all read this, so the vertex
    filler, the reagent a step consumes and the species priced against it are one
    molecule bound through one atom.
    """
    co, co_embed = embed_with_report(mol_from_smiles(smiles), seed=3)
    co_sites = perceive(co)
    if not co_sites:
        # A plain statement about the co-ligand rather than an IndexError: nothing in it
        # can bind.  `[SiH4]` and `CC` reach here, as does an uncovered donor type.
        raise _Rejected(
            f"co-ligand {smiles!r} has no perceivable donor atom, so it "
            f"cannot occupy a coordination site",
            code="co_ligand_no_donor",
            detail={"co_ligand": smiles,
                    "hint": "give the co-ligand as the species that actually "
                            "binds — '[I-]' rather than 'I2', 'O' for water, "
                            "'[OH-]' for hydroxide — or add its donor type to "
                            "sites.perception"})
    return co, co_embed, co_sites[0]


def _sphere_block(built: _Built, structure_id: int | None = None):
    """A placed sphere AS a building block: donors, empty vertices, and their state.

    Everything `assembly.join` needs and nothing the registry needs.  The vacancies come
    from the placer, which is the only thing that knows the polyhedron was bigger than the
    ligand set; the states come from `refresh_state`, because a block whose state is
    unknown refuses to answer "what is open" rather than guessing (and it is right to).
    """
    import numpy as np

    from mofsbu.assembly.join import BuildingBlock

    conf = built.mol.GetConformer()
    coords = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                        conf.GetAtomPosition(i).z]
                       for i in range(built.mol.GetNumAtoms())])
    sites = perceive(built.mol)
    if built.result.vacancies:
        sites = sites + vacancy_sites(built.result.metal_idx,
                                      coords[built.result.metal_idx],
                                      built.result.vacancies)
    symbols = [atom.GetSymbol() for atom in built.mol.GetAtoms()]
    states = refresh_state(sites, coords, graph=built.graph, symbols=symbols,
                           fidelity=Fidelity.RAW)
    return BuildingBlock(graph=built.graph, sites=tuple(sites), geometry=coords,
                         state={(s.atom_idx, s.slot): s for s in states},
                         structure_id=structure_id)


def _ligand_for(spec: BuildSpec, comp: dict[str, Any], structure_id: int | None = None):
    """The component as a joinable block, plus the donor sites the plan chose."""
    from mofsbu.assembly.construct import ligand_block

    cmol, cmolecule, embed_report = _protomer_mol(spec, comp)
    name = f"{cmolecule.name}{comp['label'] if comp['selection'] else ''}"
    block = ligand_block(cmol, charge=int(comp["charge"]), name=name,
                         multiplicity=cmolecule.multiplicity)
    # The blocks name themselves, so the choice vector records WHICH stored structures
    # were joined rather than two nulls.
    block = replace(block, structure_id=structure_id)
    by_idx = {s.atom_idx: s for s in block.sites}
    missing = [i for i in comp["donors"] if i not in by_idx]
    if missing:
        # Same disagreement `place` refuses: the donor indices were chosen at plan time
        # against a conformer embedded then, and binding atom 9 because the payload says
        # 9 would bind whatever now happens to sit there.
        raise _Rejected(
            f"donor index/indices {missing} are no longer perceived on {name}; the plan "
            f"and this executor disagree about the molecule",
            code="donor_index_stale",
            detail={"component": comp, "perceived": sorted(by_idx)})
    return block, [by_idx[i] for i in comp["donors"]], name


def _co_ligand_for(smiles: str, structure_id: int | None = None):
    """The co-ligand as a joinable block, bound through the site `place` binds it by."""
    from mofsbu.assembly.construct import ligand_block

    co, _, co_site = _co_ligand_mol(smiles)
    block = replace(ligand_block(co, charge=Chem.GetFormalCharge(co), name=smiles),
                    structure_id=structure_id)
    by_idx = {s.atom_idx: s for s in block.sites}
    if co_site.atom_idx not in by_idx:
        raise _Rejected(
            f"co-ligand {smiles!r}: the donor `place` binds (atom {co_site.atom_idx}) is "
            f"not a site of the joinable block, so the step would bind a different atom",
            code="donor_index_stale",
            detail={"co_ligand": smiles, "perceived": sorted(by_idx)})
    return block, [by_idx[co_site.atom_idx]], smiles


def _rung(reg: Registry, task_id: Any, which: str) -> tuple[int, dict[str, Any]]:
    """The structure a previous task built, and the payload that built it."""
    row = None if task_id is None else reg.conn.execute(
        "SELECT status, structure_id, payload_json FROM tasks WHERE id=?",
        (task_id,)).fetchone()
    if row is None:
        raise _Rejected(f"this step records no {which} rung to start from",
                        code="pathway_parent_missing", detail={"task": task_id})
    if row["status"] != "done" or row["structure_id"] is None:
        raise _Rejected(
            f"the {which} rung was not built ({row['status']}), so there is nothing to "
            f"join onto. The step is not wrong — the rung below it is missing, and "
            f"whatever rejected that task says why",
            code="pathway_parent_missing",
            detail={"task": task_id, "status": row["status"]})
    return int(row["structure_id"]), json.loads(row["payload_json"])


def _vertex_pair(block, donors: Sequence[Any], vacancies: Sequence[Any], metal: str,
                 elements: Sequence[str]):
    """The two vertices a chelating ligand should span, and the verdict that chose them.

    Every pair is judged and the least strained feasible one wins — which is how cis
    beats trans without anything naming either: a ligand whose donors sit 3 A apart
    subtends about 90 degrees at its bonds' own length and cannot stretch across 180.
    Ranked by the SAME verdict that will place it (`chelate_reach`), because ranking on
    one criterion and placing under another is how a step picks the pair it then cannot
    build.  Ties break on slot order so a replay picks the same pair.
    """
    from mofsbu.assembly.join import chelate_reach

    best, best_verdict, worst = None, None, None
    for first, second in itertools.combinations(vacancies, 2):
        verdict = chelate_reach(donors, [first, second], partner=metal,
                                donor_elements=elements)
        if verdict.feasible and (best_verdict is None or verdict.strain < best_verdict.strain):
            best, best_verdict = (first, second), verdict
        if worst is None:
            worst = verdict
    return best, (best_verdict or worst)


def _execute_grow(reg: Registry, task, spec: BuildSpec) -> Outcome:
    """Add one ligand to the rung below and record the step as a route to the product.

    This is the pathway half of a run.  `place` builds M(L) and M(L)2 as two independent
    constructions and nothing in the registry says the second is the first plus a ligand;
    this performs that addition — `assembly.join` onto the parent's own empty vertices —
    so the relationship is a `reactions` edge that was actually carried out rather than an
    inference from two formulas.

    Under D2 the product is usually an identity the run already has, and that is the
    RESULT, not a wasted task: one node, two routes, and the second route is the one that
    explains where the node came from.  The structure row is shared, a second geometry
    hangs off it, and the edge carries the choice vector that regenerates the step.
    """
    import numpy as np

    from mofsbu.assembly.join import IncompatibleJoin, join, join_chelate
    from mofsbu.assembly.persist import store_block

    payload = task.payload
    parent_sid, parent_payload = _rung(reg, payload.get("parent_task"), "parent")
    # The free ligand is a REAGENT of this step, not a decoration on it: an edge written
    # without it balances short by exactly one ligand and is unpriceable for a reason
    # that has nothing to do with chemistry.  Guarded like the rung below it.
    ligand_sid = None
    if payload.get("ligand_task") is not None:
        row = reg.conn.execute("SELECT status, structure_id FROM tasks WHERE id=?",
                               (payload["ligand_task"],)).fetchone()
        if row is None or row["structure_id"] is None:
            raise _Rejected(
                "the ligand this step adds was not built, so the step cannot record "
                "what it consumed. The step is not wrong — the ligand task is missing, "
                "and whatever rejected that task says why",
                code="pathway_ligand_missing",
                detail={"task": payload["ligand_task"],
                        "status": None if row is None else row["status"]})
        ligand_sid = int(row["structure_id"])

    parent = _build_sphere(spec, parent_payload)
    block = _sphere_block(parent, structure_id=parent_sid)
    comp = payload["component"]
    if "co_ligand" in comp:
        if _charged(comp["co_ligand"]):
            # The planner does not queue this; a hand-built or older queued task is refused.
            raise _Rejected(f"co-ligand {comp['co_ligand']!r} is charged: {CO_LIGAND_CHARGED}",
                            code="co_ligand_charged", detail={"co_ligand": comp["co_ligand"]})
        ligand, donors, name = _co_ligand_for(comp["co_ligand"], structure_id=ligand_sid)
    else:
        ligand, donors, name = _ligand_for(spec, comp, structure_id=ligand_sid)
    metal = parent.metal
    elements = [ligand.graph.label(s.atom_idx).element for s in donors]
    detail: dict[str, Any] = {"step": {
        "adds": name, "mode": comp["mode"], "from_structure": parent_sid,
        "from_task": payload.get("parent_task"), "ligand_structure": ligand_sid}}

    vacancies = sorted(block.open_vacancies(), key=lambda s: s.slot)
    if len(vacancies) < len(donors):
        raise _Rejected(
            f"the rung below offers {len(vacancies)} open vertex/vertices and this step "
            f"needs {len(donors)}; the parent is coordinatively saturated, so reaching "
            f"this product from it is a SUBSTITUTION, not an addition, and a join cannot "
            f"express one",
            code="pathway_no_open_vertex", detail=detail)
    try:
        if len(donors) == 2:
            pair, verdict = _vertex_pair(block, donors, vacancies, metal.symbol, elements)
            if pair is None:
                raise _Rejected(verdict.reason, code="chelate_cannot_span",
                                detail={**detail, "strain": round(verdict.strain, 4)})
            result = join_chelate(block, ligand, list(pair), donors, seed=spec.seed)
        elif len(donors) == 1:
            result = join(block, ligand, vacancies[0], donors[0], mode=comp["mode"],
                          seed=spec.seed)
        else:
            raise _Rejected(
                f"{len(donors)} donors at once: one donor is `join`, two are "
                f"`join_chelate`, and three (fac/mer tridentate) is not built",
                code="pathway_denticity_unbuilt", detail=detail)
    except IncompatibleJoin as exc:
        # A refusal from the site pair is an ANSWER about the chemistry, the same way the
        # placer's is, and it is reported as a rejection rather than a failure.
        raise _Rejected(exc.verdict.reason, code="join_refused",
                        detail={**detail, "strain": round(exc.verdict.strain, 4)}) from exc

    graph = result.block.graph
    coords = np.asarray(result.block.geometry)
    symbols = [graph.label(i).element for i in graph.nodes()]
    report = _clash_qc(graph, symbols, coords)
    if not report.ok and not (report.marginal and spec.run_mode != "construct"):
        raise _Rejected(str(report), code=report.code,
                        detail={**detail, "qc": report.to_dict()})

    added = [] if "co_ligand" in comp else [spec.molecules[comp["molecule"]].name]
    # Pieces on the product; co-ligands count only where they are stepped (D26).
    stepped_co = int(parent_payload.get("n_co", 0)) if spec.co_ligand_counts == "range" else 0
    stored = store_block(
        reg, result.block, choice_vector=result.choice_vector,
        # The reagents ARE the pathway: the rung below and the free ligand that was added.
        reagent_ids=[i for i in (parent_sid, ligand_sid) if i is not None],
        depth=sum(int(c["count"]) for c in parent.components) + stepped_co + 1,
        note=f"+{name} onto structure {parent_sid}",
        tags=[*sorted({*parent.molecule_names, *added}), "complex", metal.symbol],
        seed=spec.seed, qc=report.to_dict(),
        # The L2 tag is left where the `place` path leaves it — unset.  Computing one here
        # would put this product in a DIFFERENT structures row from the rung the run's own
        # place task built, and the edge would then point at a node nothing else reached
        # (B2 is the seam; both paths move when it is closed, together).
        l2="")
    detail["join"] = {
        "strain": round(result.strain, 4), "mode": comp["mode"],
        "reaction": stored.reaction_id, "choice_vector": stored.choice_digest,
        "product": stored.structure_id,
        "reached_existing": not stored.structure_created,
        "note": ("the step reached a structure this registry already had — one node, two "
                 "routes (D2); the edge is the new information"
                 if not stored.structure_created else
                 "this product is new: no place task built it")}
    detail["qc"] = report.to_dict()
    return Outcome(stored.structure_id, stored.geometry_id, stored.structure_created,
                   stored.geometry_created, detail)


def execute(reg: Registry, task, spec: BuildSpec) -> Outcome:
    """Run one task.  Returns what it produced AND what it took to produce it."""
    payload = task.payload
    if task.kind == "relax":
        # Handled first: a relax payload carries a structure and a geometry, not a
        # molecule index, so `_protomer_mol` has nothing to work with.
        return _execute_relax(reg, task, spec)

    if task.kind == "ligand":
        mol, molecule, embed_report = _protomer_mol(spec, payload)
        detail: dict[str, Any] = {"embed": embed_report}
        name = f"{molecule.name}{payload['label'] if payload['selection'] else ''}"
        g = from_rdkit(mol, charge=payload["charge"], multiplicity=molecule.multiplicity,
                       name=name)
        put = put_structure(reg, g, tags=[molecule.name, "ligand"])
        geom = put_geometry(reg, put.id, to_xyz(mol, name), fidelity=Fidelity.FF, method=FF)
        detail["sites"] = _record_sites(reg, put.id, geom.id, mol, g, Fidelity.FF)
        for frag in decompose(g):
            alias_fragment(reg, frag.l1, name.replace(" ", ""), source="runner")
        return Outcome(put.id, geom.id, put.created, geom.created, detail)

    if task.kind == "co_ligand":
        # The free co-ligand, stored and relaxed like a ligand: it is the reagent a
        # co-ligand step consumes, and an edge citing it prices only if it has an energy.
        mol, co_embed, _ = _co_ligand_mol(payload["smiles"])
        detail = {"embed": co_embed}
        g = from_rdkit(mol, multiplicity=1, name=payload["smiles"])
        put = put_structure(reg, g, tags=[payload["smiles"], "ligand", "co-ligand"])
        geom = put_geometry(reg, put.id, to_xyz(mol, payload["smiles"]),
                            fidelity=Fidelity.FF, method=FF)
        detail["sites"] = _record_sites(reg, put.id, geom.id, mol, g, Fidelity.FF)
        return Outcome(put.id, geom.id, put.created, geom.created, detail)

    if task.kind == "place":
        built = _build_sphere(spec, payload)
        detail = built.detail
        put = put_structure(reg, built.graph,
                            tags=[*built.molecule_names, "complex", built.metal.symbol],
                            provenance=Provenance(kind="assembly", depth=1,
                                                  note=payload["geometry"]))
        geom = put_geometry(reg, put.id, built.result.to_xyz(built.name),
                            fidelity=Fidelity.RAW, method=BUILD,
                            choice_vector=built.result.choice_vector,
                            seed=spec.seed, qc=built.result.report.to_dict())
        # The placer's empty vertices travel with the complex: metal at index 0,
        # which is how `to_rdkit` lays the centre out.
        detail["sites"] = _record_sites(reg, put.id, geom.id, built.mol, built.graph,
                                        Fidelity.RAW, vacancies=built.result.vacancies,
                                        metal_idx=built.result.metal_idx)
        detail["qc"] = built.result.report.to_dict()
        return Outcome(put.id, geom.id, put.created, geom.created, detail)

    if task.kind == "grow":
        return _execute_grow(reg, task, spec)

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


#: The two halves of a run, and the reason a worker pool is not homogeneous.  Everything
#: in `BUILD_KINDS` is CPU work — embed, perceive, place, hash — and scales with cores.
#: `relax` is the optimiser, which on a declared accelerator is one device's worth of work
#: however many processes ask for it.
BUILD_KINDS = ("ligand", "co_ligand", "place", "grow")
RELAX_KINDS = ("relax",)


def work(reg: Registry, spec: BuildSpec, run_id: int, *, limit: int | None = None,
         kinds: Sequence[str] | None = None,
         exclude_kinds: Sequence[str] | None = None,
         wait_while: Any = None, poll: float = 0.25) -> int:
    """Claim and execute tasks until none remain.  In-process.

    `kinds` / `exclude_kinds` restrict this worker to part of the queue (see
    `BUILD_KINDS`).  `wait_while` is a zero-argument predicate: when the queue hands back
    nothing and it returns True, the worker sleeps `poll` seconds and asks again instead
    of exiting.  That is what a relax worker needs and a build worker does not — relax
    tasks are queued BY the build tasks as they finish, so an empty relax queue early in
    a run means "not yet", while an empty build queue means "never again".
    """
    done = 0
    while limit is None or done < limit:
        if cancel_requested(reg, run_id):
            # Cooperative: the task in hand has already finished, and nothing new is
            # claimed.  Everything computed so far is in the registry.
            break
        task = claim_task(reg, run_id, kinds=kinds, exclude_kinds=exclude_kinds)
        if task is None:
            if wait_while is not None and wait_while():
                time.sleep(poll)
                continue
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
        # After the task, not only at the claim: a worker part-way through a long queue
        # is the case a heartbeat exists to distinguish from one that died at the first.
        touch_run(reg, run_id)
        done += 1
    return done


def _worker_process(db_path: str, store_path: str, spec_json: str, run_id: int,
                    kinds: Sequence[str] | None = None,
                    exclude_kinds: Sequence[str] | None = None,
                    more_coming: Any = None) -> None:
    """One worker, in its own process.  Arguments only — nothing is inherited.

    `more_coming` is an `mp.Event` the parent sets once the build workers have all
    exited: while it is clear, an empty relax queue means the builders have not caught up
    yet, and the relax worker waits rather than finishing a run that is not over.
    """
    from mofsbu.registry import BlobStore

    with Registry(Path(db_path), BlobStore(Path(store_path))) as reg:
        work(reg, BuildSpec.from_json(spec_json), run_id, kinds=kinds,
             exclude_kinds=exclude_kinds,
             wait_while=(None if more_coming is None else
                         lambda: not more_coming.is_set()))


@dataclass(frozen=True)
class WorkerPlan:
    """How a run's declared worker count is divided between the two queues.

    Separated from the spawning so the page and the CLI can state the division before a
    run starts, and so it is testable without starting a process.
    """

    total: int = 1
    build: int = 1
    relax: int = 0

    @property
    def split(self) -> bool:
        return self.relax > 0

    def describe(self) -> str:
        if self.total <= 1:
            return "1 worker, in-process"
        if not self.split:
            return f"{self.total} workers, any task"
        return f"{self.build} building + {self.relax} relaxing"


def plan_workers(workers: int | None = None, *, relaxes: bool = True) -> WorkerPlan:
    """Divide the declared workers between construction and relaxation.

    An explicit `workers` is a declaration in its own right (a `--workers` flag, the
    number typed on the builder page), so it is honoured as given; `None` means "whatever
    this machine declares", which on an unconfigured laptop is one in-process worker and
    no subprocesses at all (ground rule 8).

    The relax share comes from `config.relax_workers` and is capped so at least one
    worker is left building: a pool that is ALL relax workers would leave the queue that
    feeds them empty and the machine with one busy core, which is the shape this split
    exists to undo.  `relaxes=False` — a `construct` run — has no second queue at all, and
    reserving a worker for it would idle a core for the length of the run.
    """
    total = workers if workers is not None else max_workers()
    total = max(1, int(total))
    if total <= 1 or (workers is None and not parallel_enabled()):
        return WorkerPlan(1, 1, 0)
    n_relax = min(max(0, relax_workers()) if relaxes else 0, total - 1)
    return WorkerPlan(total, total - n_relax, n_relax)


def execute_run(reg: Registry, spec: BuildSpec, run_id: int, *,
                workers: int | None = None) -> WorkerPlan:
    """Drain a planned run with whatever parallelism this machine declares.

    The one place a worker pool is built.  Every caller that drains a queue — the CLI,
    the builder page's background thread — comes through here, so no entry point can be
    parallel while another is quietly serial.

    Processes are started with the **spawn** context, not forked.  A forked child cannot
    re-initialise CUDA, so on the one configuration this split is for — a GPU relaxing
    while the cores build — fork produces workers that die at their first relaxation.
    Spawn also keeps a fork out of the web server, which is a threaded process.

    A spawned worker starts from nothing, so exactly three things reach it: the DATABASE
    it is pointed at, the SPEC as JSON, and the ENVIRONMENT (which is where the device,
    the model and the worker count live — `config` reads them at execution time, so a
    declaration made in this process does reach the children).  In-process state does
    not, which is worth knowing before substituting a backend in memory and expecting a
    pool to use it.
    """
    # A `construct` run has no relax queue, so it gets no relax worker: the division
    # follows the work that exists, not the hardware alone.
    pool = plan_workers(workers, relaxes=spec.run_mode != "construct")
    if pool.total <= 1:
        work(reg, spec, run_id)
        return pool

    import multiprocessing as mp

    ctx = mp.get_context("spawn")
    args = (str(reg.db_path), str(reg.store.root), spec.to_json(indent=None), run_id)
    # Build workers never claim a relax task when a relax worker exists: one shared queue
    # puts construction behind the optimiser, and on a GPU it also means several
    # processes on one card.
    builders = [ctx.Process(target=_worker_process, args=args,
                            kwargs={"exclude_kinds": RELAX_KINDS if pool.split else None},
                            name=f"mofsbu-build-{run_id}-{i}")
                for i in range(pool.build)]
    builds_done = ctx.Event() if pool.split else None
    relaxers = [ctx.Process(target=_worker_process, args=args,
                            kwargs={"kinds": RELAX_KINDS, "more_coming": builds_done},
                            name=f"mofsbu-relax-{run_id}-{i}")
                for i in range(pool.relax)]
    for proc in [*builders, *relaxers]:
        proc.start()
    try:
        for proc in builders:
            proc.join()
    finally:
        if builds_done is not None:
            # Set even if the join raised: a relax worker waiting on this event is the
            # one thing that would not stop on its own.
            builds_done.set()
    for proc in relaxers:
        proc.join()

    # A worker that DIED is not a task that failed: it recorded nothing, so its share of
    # the queue is still pending.  `finish_run` refuses to call that done; this says which
    # process it was and what killed it.
    dead = [f"{p.name} (exit {p.exitcode})"
            for p in [*builders, *relaxers] if p.exitcode not in (0, None)]
    if dead:
        raise MofsbuError(
            f"{len(dead)} of {pool.total} workers died before the queue was drained: "
            + ", ".join(dead) + ". Their tasks are still pending; the run can be "
            "resumed once the cause is fixed")
    return pool


def run(reg: Registry, spec: BuildSpec, *, workers: int | None = None) -> dict[str, Any]:
    """Plan and execute a spec.  Serial on a laptop; parallel only where declared."""
    # Starting work on a database is the moment to notice that the last process to work
    # on it never came back — a run killed with its shell writes no ending, and until
    # something says so it reads as still going.  Reported, because a run that silently
    # changed status between two commands is worse than one that says it did.
    for swept in sweep_interrupted(reg):
        print(f"mofsbu: run {swept['run_id']} is {swept['now']} — {swept['reason']}"
              + (f"; {swept['returned_claims']} task(s) returned to the queue"
                 if swept["returned_claims"] else ""))
    run_id, n_tasks = plan(reg, spec)
    reg.conn.commit()
    pool = plan_workers(workers, relaxes=spec.run_mode != "construct")
    if spec.run_mode != "construct":
        # Printed before any work starts: an accelerator sitting idle for a whole run is
        # otherwise only visible in nvidia-smi, an hour later, by accident.  The division
        # is printed for the same reason — a GPU fed by a single core is the other way to
        # own the hardware and not use it.
        print(f"mofsbu run {run_id}: mode={spec.run_mode}  {device_note()}  "
              f"workers={pool.total} ({pool.describe()})")

    try:
        execute_run(reg, spec, run_id, workers=workers)
    finally:
        # Closed out even when a worker died, so the run row says what happened instead
        # of staying `pending` for ever while the exception travels.
        from mofsbu.registry import relabel_all

        relabel_all(reg)      # labels are a projection; refresh once the fragments exist
        status = finish_run(reg, run_id)
    return {"run_id": run_id, "tasks": n_tasks, "status": status,
            "counts": task_counts(reg, run_id), "workers": pool.total,
            "build_workers": pool.build, "relax_workers": pool.relax,
            "summary": outcome_summary(reg, run_id)}
