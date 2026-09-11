"""Specs, the task queue, and the two ground rules they exist to serve."""
from __future__ import annotations

import os

import pytest

from mofsbu import config
from mofsbu.registry import BlobStore, Registry, find
from mofsbu.registry.jobs import (
    add_task, claim_task, complete_task, create_run, fail_task, finish_run,
    reset_stale_claims, task_counts,
)
from mofsbu.runner import plan, run
from mofsbu.spec import BuildSpec, MetalSpec, MoleculeSpec, PocketPredicate


@pytest.fixture()
def reg(tmp_path):
    with Registry(tmp_path / "r.db", BlobStore(tmp_path / "s")) as r:
        r.migrate()
        yield r


def catechol_spec(**over) -> BuildSpec:
    base = dict(
        molecules=(MoleculeSpec("catechol", "Oc1ccccc1O", 1),),
        metals=(MetalSpec("Zn", 2, "ls"),), degree=1, coordination=(4,),
        binding=("chelate",), pocket=PocketPredicate(ring_size=5, n_anionic=2),
        ligands_per_metal=(1,), seed=7)
    base.update(over)
    return BuildSpec(**base)


# ── a build is data ──────────────────────────────────────────────────────────

def test_spec_round_trips_and_hashes(tmp_path):
    spec = catechol_spec(note="round trip")
    path = spec.save(tmp_path / "spec.json")
    assert BuildSpec.load(path) == spec
    assert BuildSpec.load(path).digest == spec.digest


def test_identical_specs_share_a_digest_and_different_ones_do_not():
    assert catechol_spec().digest == catechol_spec().digest
    assert catechol_spec().digest != catechol_spec(coordination=(6,)).digest


# ── the queue ────────────────────────────────────────────────────────────────

def test_a_task_is_claimed_exactly_once(reg):
    """Two workers pulling at once must not get the same row."""
    run_id = create_run(reg, catechol_spec())
    for i in range(6):
        add_task(reg, run_id, "place", {"i": i})
    claimed = []
    while (task := claim_task(reg, run_id)) is not None:
        claimed.append(task.id)
        complete_task(reg, task.id)
    assert len(claimed) == len(set(claimed)) == 6
    assert claim_task(reg, run_id) is None


def test_rejected_is_not_failed(reg):
    """A candidate that fails QC is an answer, not a crash.

    Conflating them makes 'the run produced nothing' and 'the run broke' look identical.
    """
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    add_task(reg, run_id, "place", {})
    a, b = claim_task(reg, run_id), claim_task(reg, run_id)
    fail_task(reg, a.id, "QC said no", rejected=True)
    complete_task(reg, b.id)
    assert task_counts(reg, run_id) == {"rejected": 1, "done": 1}
    assert finish_run(reg, run_id) == "done"          # rejections do not fail the run


def test_a_dead_worker_returns_its_task_to_the_queue(reg):
    run_id = create_run(reg, catechol_spec())
    add_task(reg, run_id, "place", {})
    claim_task(reg, run_id)                            # claimed, then the worker dies
    assert claim_task(reg, run_id) is None
    assert reset_stale_claims(reg, run_id) == 1
    assert claim_task(reg, run_id) is not None


# ── ground rule 8: the laptop is the default ─────────────────────────────────

def test_parallelism_is_opt_in(monkeypatch):
    monkeypatch.delenv("MOFSBU_PROFILE", raising=False)
    monkeypatch.delenv("MOFSBU_WORKERS", raising=False)
    assert config.machine_profile() == config.LAPTOP
    assert config.max_workers() == 1
    assert not config.parallel_enabled()

    monkeypatch.setenv("MOFSBU_WORKERS", "4")
    assert config.max_workers() == 4 and config.parallel_enabled()

    monkeypatch.delenv("MOFSBU_WORKERS")
    monkeypatch.setenv("MOFSBU_PROFILE", "workstation")
    assert config.machine_profile() == config.WORKSTATION
    assert config.max_workers() == max(1, (os.cpu_count() or 2) - 1)


def test_a_run_on_an_undeclared_machine_stays_in_process(reg, monkeypatch):
    monkeypatch.delenv("MOFSBU_PROFILE", raising=False)
    monkeypatch.delenv("MOFSBU_WORKERS", raising=False)
    import multiprocessing as mp

    def explode(*_a, **_k):                      # pragma: no cover - must never run
        raise AssertionError("spawned a process on an undeclared machine")

    monkeypatch.setattr(mp, "Process", explode)
    summary = run(reg, catechol_spec())
    assert summary["workers"] == 1
    assert summary["counts"].get("done", 0) > 0


# ── ground rule 7: missing functions are stubs that shout ────────────────────

def test_degree_above_one_refuses_rather_than_doing_something_smaller(reg):
    from mofsbu.assembly.join import NotBuiltYet

    with pytest.raises(NotBuiltYet, match="M5"):
        plan(reg, catechol_spec(degree=2))


def test_the_assembly_interfaces_exist_and_all_raise():
    """Signatures are settled so callers can be written; bodies are scheduled work."""
    from mofsbu.assembly.join import NotBuiltYet, compatible, grow, join
    from mofsbu.geometry.placer import place_multicentre

    for call in (lambda: compatible(None, None),
                 lambda: join(None, None, None, None),
                 lambda: grow(None, (), degree=2),
                 lambda: place_multicentre([1, 2], [], [])):
        with pytest.raises(NotBuiltYet):
            call()


# ── end to end ───────────────────────────────────────────────────────────────

def test_a_spec_produces_structures(reg):
    summary = run(reg, catechol_spec(ligands_per_metal=(1, 2)))
    assert summary["status"] == "done"
    labels = [r["display_label"] for r in find(reg, limit=100)]
    assert any(lbl.startswith("Zn[") for lbl in labels)
    assert any("catechol" in lbl for lbl in labels)     # fragment alias applied on relabel


# ── a run that builds nothing must say why ───────────────────────────────────

def anthrarufin_spec(**over) -> BuildSpec:
    """The user's actual failing case: a bidentate ligand, no co-ligand, CN 4 and 6."""
    base = dict(
        molecules=(MoleculeSpec("ATF", "C1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C(=CC=C3)O", 1),),
        metals=(MetalSpec("Cu", 2, "ls"),), degree=1, coordination=(4, 6),
        binding=("chelate", "mono"), co_ligand=None, ligands_per_metal=(1,))
    base.update(over)
    return BuildSpec(**base)


def test_a_run_that_queues_nothing_explains_itself(reg):
    """Silence was the bug: a chelate needs 2 more sites at CN 4 and there is no
    co-ligand, so every candidate was skipped and the run looked like a success that
    simply produced no complexes."""
    from mofsbu.registry.jobs import get_diagnostics
    from mofsbu.runner import plan

    run_id, n_tasks = plan(reg, anthrarufin_spec())
    diagnostics = get_diagnostics(reg, run_id)
    assert diagnostics, "a run that builds no complexes must record why"
    assert any("co-ligand" in d["reason"] for d in diagnostics)
    assert all(d["hint"] for d in diagnostics)
    # ligands are still built; it is the coordination step that had nothing to do
    kinds = {r["kind"] for r in reg.conn.execute(
        "SELECT DISTINCT kind FROM tasks WHERE run_id=?", (run_id,))}
    assert kinds == {"ligand"}


def test_setting_a_co_ligand_builds_the_complexes(reg):
    from mofsbu.registry import find

    summary = run(reg, anthrarufin_spec(co_ligand="O"))
    assert summary["counts"].get("done", 0) > 0
    assert find(reg, metal="Cu", limit=50), "Cu complexes should now exist"


def test_the_placer_refusing_is_a_rejection_not_a_failure(reg):
    """A bidentate ligand cannot span a linear two-coordinate centre.

    Saying so is correct behaviour; recording it as `failed` would make a run full of
    sound chemistry look like a run full of bugs.
    """
    summary = run(reg, anthrarufin_spec(allow_unsaturated=True))
    assert summary["counts"].get("failed", 0) == 0
    assert summary["counts"].get("rejected", 0) > 0


# ── specs are files, so the schema needs migrations ──────────────────────────

def test_a_version_1_spec_still_loads(tmp_path):
    """Changing the dataclass without a migration silently breaks every spec on disk."""
    import json

    v1 = {"spec_version": 1, "molecules": [{"name": "x", "smiles": "O", "multiplicity": 1,
                                            "max_deprotonations": None}],
          "metals": [], "degree": 1, "coordination": [4], "geometries": None,
          "ligands_per_metal": [1], "co_ligand": "O", "binding": ["chelate"],
          "pocket": {"ring_size": None, "n_anionic": None, "min_anionic": None,
                     "rigid": None},
          "seed": 0, "relax_to": None, "note": "old"}
    path = tmp_path / "old.json"
    path.write_text(json.dumps(v1))
    spec = BuildSpec.load(path)
    assert spec.run_mode == "construct" and spec.allow_unsaturated is False
    # v4 added `ml_model`.  A v1 spec never expressed a choice of ML potential, so it
    # comes forward as None ("use what this machine declares") rather than being
    # back-filled with "mace-mp-0", which would invent a decision nobody made.
    assert spec.ml_model is None


def test_a_version_3_spec_comes_forward_without_inventing_a_model(tmp_path):
    spec = BuildSpec.from_dict({
        "spec_version": 3, "run_mode": "ml_go",
        "molecules": [{"name": "x", "smiles": "O", "multiplicity": 1,
                       "max_deprotonations": None}]})
    assert spec.ml_model is None
    assert spec.to_dict()["spec_version"] == 6
    # v6 adds mixed-ligand enumeration.  A v3 spec meant one molecule per coordination
    # sphere, so it comes forward homoleptic — defaulting it to anything else would
    # silently multiply the size of every run already on disk.
    assert spec.max_distinct_ligands == 1


def test_a_version_4_spec_drops_the_stale_metal_multiplicity():
    """v5 removes `MetalSpec.multiplicity` — it was never kept in sync with
    `oxidation_state`/`spin_class`, which is exactly how a Cu(II) centre ended up
    stored as a singlet.  A v4 spec's stored value is dropped, not trusted."""
    spec = BuildSpec.from_dict({
        "spec_version": 4,
        "molecules": [{"name": "x", "smiles": "O", "multiplicity": 1,
                       "max_deprotonations": None}],
        "metals": [{"symbol": "Cu", "oxidation_state": 2, "spin_class": "ls",
                    "multiplicity": 1}]})
    assert not hasattr(spec.metals[0], "multiplicity")
    assert spec.metals[0].spin_class == "ls"


def test_a_spec_may_name_the_ml_potential_and_a_typo_is_refused():
    """`ml_go` is a rung, not a theory; the spec is where the theory is recorded."""
    mol = [{"name": "x", "smiles": "O", "multiplicity": 1, "max_deprotonations": None}]
    spec = BuildSpec.from_dict({"spec_version": 4, "molecules": mol,
                                "run_mode": "ml_go", "ml_model": "MACE-OMOL-0"})
    assert spec.ml_model == "MACE-OMOL-0"
    # Two specs differing only in the model are two different runs, not one.
    other = BuildSpec.from_dict({"spec_version": 4, "molecules": mol,
                                 "run_mode": "ml_go", "ml_model": "mace-mp-0"})
    assert spec.digest != other.digest
    with pytest.raises(ValueError, match="unknown ML model"):
        BuildSpec.from_dict({"spec_version": 4, "molecules": mol,
                             "run_mode": "ml_go", "ml_model": "mace-omol-1"})


def test_a_future_spec_version_is_refused_not_guessed():
    with pytest.raises(ValueError, match="newer than this build"):
        BuildSpec.from_dict({"spec_version": 999, "molecules": [
            {"name": "x", "smiles": "O", "multiplicity": 1, "max_deprotonations": None}]})


def test_an_unknown_field_is_refused():
    with pytest.raises(ValueError, match="unknown spec field"):
        BuildSpec.from_dict({"spec_version": 2, "molecules": [
            {"name": "x", "smiles": "O", "multiplicity": 1, "max_deprotonations": None}],
            "wishful_thinking": True})


def test_non_spec_json_is_not_mistaken_for_a_spec():
    assert not BuildSpec.looks_like_spec({"molecules": [{"name": "a", "smiles": "O"}]})
    assert BuildSpec.looks_like_spec({"spec_version": 2, "molecules": []})


# ── heteroleptic enumeration: more than one ligand kind on a centre ──────────

def _place_compositions(reg, spec):
    from mofsbu.runner import plan
    import json as _json

    run_id, _ = plan(reg, spec)
    out = []
    for r in reg.conn.execute(
            "SELECT payload_json FROM tasks WHERE run_id=? AND kind='place'", (run_id,)):
        p = _json.loads(r["payload_json"])
        out.append(tuple(sorted(
            (spec.molecules[c["molecule"]].name, bool(c["selection"]), c["count"])
            for c in p["components"])))
    return set(out)


def _two_ligand_spec(**over):
    base = dict(
        molecules=(MoleculeSpec("A", "CC(=O)C", 1), MoleculeSpec("B", "Cl", 1)),
        metals=(MetalSpec("Mg", 2, "ls"),), coordination=(4,),
        ligands_per_metal=(2,), binding=("mono",), co_ligand=None,
        allow_unsaturated=True, geometries=("tetrahedral",))
    base.update(over)
    return BuildSpec(**base)


def test_one_distinct_ligand_is_homoleptic_and_is_the_default(reg):
    """The default must reproduce what every spec written before this field meant.

    Two molecules in a spec used to be enumerated independently and could never share a
    coordination sphere; that is still what happens unless you ask for otherwise, so no
    stored run changes size or content underneath you.
    """
    spec = _two_ligand_spec()
    assert spec.max_distinct_ligands == 1
    combos = _place_compositions(reg, spec)
    assert all(len(c) == 1 for c in combos), f"a default spec mixed ligands: {combos}"


def test_raising_the_cap_mixes_two_molecules_on_one_centre(reg):
    combos = _place_compositions(reg, _two_ligand_spec(max_distinct_ligands=2))
    mixed = {c for c in combos if len({name for name, _, _ in c}) > 1}
    assert mixed, "asked for 2 distinct ligand kinds and got only homoleptic spheres"
    assert any(len(c) == 1 for c in combos), "the cap is a maximum, not a target"


def test_two_protomers_of_one_molecule_are_two_kinds(reg):
    """HCl and Cl- may share a metal: a partially deprotonated set is real chemistry
    for a polyprotic linker, and this is the decision that allows it."""
    combos = _place_compositions(reg, _two_ligand_spec(max_distinct_ligands=2))
    same_mol_mixed = {
        c for c in combos
        if len(c) == 2 and len({name for name, _, _ in c}) == 1
        and len({deprot for _, deprot, _ in c}) == 2}
    assert same_mol_mixed, "no complex mixed a protonated and deprotonated form"


def test_a_composition_never_exceeds_the_requested_ligand_count(reg):
    for combo in _place_compositions(reg, _two_ligand_spec(max_distinct_ligands=2)):
        assert sum(count for _, _, count in combo) == 2


def test_max_distinct_ligands_must_be_at_least_one():
    with pytest.raises(ValueError, match="at least 1"):
        _two_ligand_spec(max_distinct_ligands=0)
