"""Pinned algorithm versions, and which checkout is answering right now.

Two different questions, in one file because both are "what version", and worth keeping
apart because they behave oppositely:

* `ALGO_VERSIONS` is what produced a STORED ROW.  Ground rule 6
  (docs/PLAN_implementation.md §0): every stored value carries the version of the recipe
  that produced it.  Bump one here when its algorithm changes; never silently re-label
  existing rows.  These travel with the data and outlive the process.
* `build_info()` is which CODE is running.  It is never stored and never enters an
  identity — it answers "what am I looking at", which is a question you ask of a server,
  not of a row.  The pages show it so that a branch under test is distinguishable from
  `main` at a glance.
"""
from __future__ import annotations

import functools
import subprocess

ALGO_VERSIONS: dict[str, str] = {
    "graph_schema":    "1",       # NodeLabel/EdgeType layout + what enters a hash
    "l0_composition":  "1",       # L0 string format
    "l1_certificate":  "cert1",   # sha256 over the IR canonical certificate — the L1 key
    "wl_index":        "wl3-nx",  # networkx WL, 3 iterations — fast bucket index only
    "canonical_order": "ir1",     # individualisation-refinement canonical labelling
    # What counts as an arrangement, how donors are grouped into classes, and where
    # the cis/trans split sits.  Structures tagged under a different version are not
    # comparable at L2 and are never re-labelled (ground rule 6).
    "l2_isomer_tag":   "iso1",    # M5/S5 — cis/trans, fac/mer, Delta/Lambda (D10)
    # The L3 label is provenance-primary (D11): the choice vector names the conformer and
    # geometry only verifies.  So this version covers both halves — what `core_atoms` calls
    # the rigid core, and theta_geom, the clustering threshold that collapses stochastic
    # duplicates (D21).  Rows written under `0-stub` carry no label at all rather than a
    # different one, so they are not comparable and are never re-derived (ground rule 6).
    "l3_conformer_id": "conf1",   # M5/S6 — choice-vector label + geometric verifier
    # How a choice vector becomes a key: what is dropped (seed, a self-referential
    # digest), how floats are quantised, how it serialises.  It prefixes the digest it
    # produces rather than living in a column, because a choice digest is not part of an
    # address — see `assembly.choice`.  Bump it and every later digest stops comparing
    # equal to the ones already stored, which is the intended and only honest outcome.
    # 2: a join records WHICH LONE PAIR the metal bound.  An sp2 donor has two in-plane
    #    lobes and they are 2.8 A apart in the M...M they imply, so the same recorded path
    #    under `cv1` is a path that did not say which one it took — it took lobe 0 because
    #    that was the only one selectable, which replays correctly and is not the same
    #    statement as choosing it.  Every digest moves, which is the honest outcome: a key
    #    over a larger set of coordinates is a different key.
    "choice_vector":   "cv2",     # M6/S1 — the L3 label's producer (D11, D13)
    # 2: the MACE backends became a family (MP-0 / OMOL-0).  A method row now records
    #    `training_set`, and `spin_blind` alongside `charge_blind`, so a stored ML number
    #    says which foundation model made it instead of only "mace".  Old rows keep
    #    algo=1 and are never re-labelled (ground rule 6); they are simply not the same
    #    method as anything computed from here on, which is correct — they are not.
    "energy_backends":  "2",      # backend protocol + how a MethodSpec is filled in
    "reference_scheme": "balanced2",   # + charge_separation, corrected media (C16-C18)
    "spin_convention":  "hs1",    # high_spin_multiplicity: d-count table + charge
    "descriptor_tables": "1",     # curated donor table + generated metal table (M1, C8)
    # The C5 floor as ratified in D18: which components exist, how they are weighted,
    # what confidence and provisional mean.  Bumping this is what marks stored ease rows
    # as products of a different policy — the numbers are only comparable within a
    # version, because the weights ARE the model.
    "ease_model":       "floor1",  # pKa-table floor + cone occlusion (M4, C5/D18)
    # Which atoms are donors and what type each one is.  Was an unpinned literal inside
    # `put_sites`; pinned here because a catalog is only comparable to one written by the
    # same recipe, and because the staleness check needs something to compare against.
    # 2: coordination is not constitution.  `1` counted a metal as an ordinary heavy
    #    neighbour, so `aqua_O`, `ether_O` and `carbonyl_O` were perceived while free and
    #    NOT perceived once bound — a donor left the catalog at the moment it became
    #    occupied, and the bond that formed had no row to be recorded against.  A `1`
    #    catalog on a metal-bearing structure is therefore incomplete rather than merely
    #    old, which is why `put_sites` replaces one instead of keeping it.
    # 3: a site records EVERY lone pair, not just the first.  An sp2 donor has two
    #    in-plane lobes and `site_frame` always knew it; `perceive` kept one, so a stored
    #    carboxylate pointed its metal at the syn lobe and the anti lobe existed nowhere.
    #    Same donors and same types as `2` — what changed is the frame record, and it
    #    changed in the one way that is not safely readable as "absent": a `2` frame has
    #    no `lone_pairs` key and therefore reads as a donor with exactly ONE lobe, which
    #    is the correct answer for an aqua and the wrong one for a carboxylate.  A bridge
    #    judged against a `2` catalog would refuse, for a reason that is an artifact of
    #    when the row was written rather than a fact about the chemistry.  Old rows keep
    #    their own answer and are rewritten the next time anything touches the structure
    #    (D19), which is how the `2` bump behaved too.
    "perception":       "3",       # which atoms are donors, of what type, and their lobes
    "site_state":       "1",       # status rule + buried-volume convention (M4)
    # 2: pinning a donor and its outward axis leaves two rotations undetermined — the
    #    ligand's spin about the M-L axis, and the metal's swing out of the donor's
    #    plane — and `1` left both wherever the alignment arithmetic dropped them.  They
    #    are now searched over fixed grids and the winning indices are recorded, so the
    #    coordinates of every monodentate placement differ from a `1` one.  Old rows keep
    #    algo=1 and are never re-labelled (ground rule 6): they are not a worse version
    #    of this recipe, they are a different one, and only comparable among themselves.
    #    The grid sizes and the out-of-plane cap are part of the recipe, because an index
    #    only means anything against a grid of a known size.
    # 3: WHICH vertices a centre leaves open is now chosen, not left to fill order.  A
    #    centre with two or more vacancies reserves a mutually-cis set for them, so the
    #    ligands take the remaining vertices and the empty ones are where the next step
    #    can reach: a CN-6 rung carrying co-ligands used to hand back a TRANS pair, and a
    #    ~90 deg chelate was refused across it (measured, one centre and one ligand:
    #    strain 2.156 trans against 0.094 cis).  Ligands therefore sit on different
    #    vertices than a `2` build put them on, so the coordinates differ wherever a
    #    centre has two or more vacancies.  A saturated centre, and one with a single
    #    vacancy, have no arrangement to choose and are byte-identical to `2`.  Old rows
    #    keep algo=2 and are never re-labelled (ground rule 6 / D19).
    "placement":        "3",       # how a ligand is oriented on a coordination vertex
}


def version_block() -> str:
    return "\n".join(f"{k:18s} {v}" for k, v in sorted(ALGO_VERSIONS.items()))


def _git(*args: str) -> str:
    """One git question, answered from the checkout this package lives in, or "".

    Never raises.  A checkout is a convenience here, not a requirement — an installed
    wheel has no `.git` and is a perfectly good way to run this — so a missing git, a
    missing repository and a timeout all mean the same thing: no branch to report.
    """
    from mofsbu.config import REPO_ROOT

    try:
        out = subprocess.run(("git", *args), cwd=REPO_ROOT, capture_output=True,
                             text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return ""
    return out.stdout.strip() if out.returncode == 0 else ""


@functools.lru_cache(maxsize=1)
def build_info() -> dict[str, str]:
    """Which code this process is running: package version, branch, commit, checkout.

    **Read once, at first ask, and cached for the life of the process** — deliberately.
    The question being answered is "what am I running", and what a running server is
    running is what it imported at start-up; re-reading git would make the stamp track
    the working tree instead of the process, so a page could report a commit whose code
    is not the code answering the request.  Restarting the server is what changes it,
    which is also when it actually changes.

    `dirty` is the one soft edge: it is the working tree's state at start-up, and it is
    reported because "0.0.1 on main" means something quite different with uncommitted
    changes under it.
    """
    from mofsbu import __version__

    commit = _git("rev-parse", "HEAD")
    branch = _git("rev-parse", "--abbrev-ref", "HEAD")
    from mofsbu.config import REPO_ROOT

    return {
        "version": __version__,
        "branch": branch,
        "commit": commit,
        "short": commit[:7],
        "dirty": "1" if _git("status", "--porcelain") else "",
        "committed_at": _git("log", "-1", "--format=%cs"),
        "checkout": str(REPO_ROOT),
        "source": "git" if commit else "installed package (no checkout)",
    }


def build_line() -> str:
    """The stamp as one line, for a console banner and the page's tooltip."""
    info = build_info()
    where = (f"{info['branch']} @ {info['short']}{'*' if info['dirty'] else ''}"
             if info["commit"] else info["source"])
    when = f", {info['committed_at']}" if info["committed_at"] else ""
    return f"mofsbu {info['version']} — {where}{when}"
