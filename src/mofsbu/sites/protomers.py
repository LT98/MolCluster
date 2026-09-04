"""Protonation states as an enumerated branch, not a hard-coded choice.

"Activate at viable sites" is the project's founding primitive, and a molecule with
several acidic protons has a *set* of activation states, not one.  Four hydroxyls give
16 subsets, and they are not 16 structures: by the molecule's own symmetry they collapse
to seven.

**Configurations are DERIVED, not named.**  An earlier version classified the
doubly-deprotonated cases as ortho / meta / para from a three-entry lookup on ring
distance.  That is a taxonomy, and a taxonomy is a promise to enumerate chemistry in
advance — it says nothing about a five-membered ring, about two sites on different
rings, about fac/mer, or about anything nobody wrote down.  Instead:

* which selections are equivalent is decided by **Aut(parent)**, the automorphism group
  of the un-deprotonated molecule.  Two selections are the same configuration exactly
  when a symmetry of the parent maps one onto the other.  This needs no chemistry
  knowledge and is correct for any scaffold.
* the configuration's *name* is the **canonical orbit representative** — the
  lexicographically smallest set of canonical site ordinals in that orbit.  It is
  generated, unique within the family, and meaningless-looking on purpose.
* a **separation profile** (sorted pairwise topological distances between the chosen
  sites) is carried alongside as an invariant humans can read.

Familiar words like "ortho" can be attached later as display-only aliases keyed on the
descriptor, in the same way fragments get names.  Nothing in the pipeline depends on
anyone having supplied one, and a scaffold no one anticipated still gets a correct,
distinct descriptor.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import networkx as nx
from rdkit import Chem

from mofsbu.graph._types import TypedGraph
from mofsbu.graph.canon import automorphism_group, canonical_index_map
from mofsbu.graph.from_mol import from_rdkit
from mofsbu.identity import l1_graph_hash
from mofsbu.sites.perception import DonorSite, deprotonate, find_donor_sites


@dataclass(frozen=True)
class Protomer:
    """One deprotonation state, with the selections that all lead to it."""

    n_deprotonated: int
    charge: int
    l1: str
    configuration: str                  # derived orbit representative, e.g. "0,1"
    separation: tuple[int, ...]         # sorted pairwise topological distances
    selections: tuple[tuple[int, ...], ...]
    representative: tuple[int, ...]     # the selection used to build `mol`
    mol: Chem.Mol
    graph: TypedGraph

    @property
    def label(self) -> str:
        """Generated, not named.  `cfg` is the orbit representative; `sep` the profile."""
        if not self.n_deprotonated:
            return "neutral"
        tag = f"-{self.n_deprotonated}H"
        if not self.configuration:
            return tag
        sep = ",".join(str(d) for d in self.separation)
        return f"{tag}cfg({self.configuration})sep({sep})"


def site_ordinals(parent: TypedGraph, site_atoms: list[int]) -> dict[int, int]:
    """Number the activation sites canonically, so ordinals do not depend on input order."""
    cmap = canonical_index_map(parent)
    ranked = sorted(site_atoms, key=lambda a: cmap[a])
    return {atom: k for k, atom in enumerate(ranked)}


def configuration_key(parent: TypedGraph, site_atoms: list[int],
                      selection: tuple[int, ...]) -> str:
    """The canonical representative of this selection's orbit under Aut(parent).

    Two selections related by a symmetry of the parent molecule are the same
    configuration and get the same key.  No chemistry is encoded: the group decides.
    """
    if not selection:
        return ""
    # A SINGLE deprotonation needs a descriptor too.  On a molecule whose acidic sites
    # are symmetry-equivalent (four identical hydroxyls) every choice gives one class and
    # the descriptor is constant; on one whose sites are inequivalent (a carboxylic acid
    # and a phenol) the choices are different species, and returning "" for both labelled
    # two distinct structures identically.
    ordinals = site_ordinals(parent, site_atoms)
    chosen = [site_atoms[i] for i in selection]
    best: tuple[int, ...] | None = None
    for perm in automorphism_group(parent):
        mapped = tuple(sorted(ordinals[perm[a]] for a in chosen if perm[a] in ordinals))
        if len(mapped) != len(chosen):
            continue                       # symmetry moved a site off the site set
        if best is None or mapped < best:
            best = mapped
    return ",".join(str(i) for i in (best or ()))


def separation_profile(parent: TypedGraph, atoms: list[int]) -> tuple[int, ...]:
    """Sorted pairwise topological distances — an invariant a human can read.

    Not used for equivalence (the orbit key decides that); carried because "these two
    sites are one bond apart" is the thing a chemist actually wants to know, and it
    generalises to any scaffold and any number of sites.
    """
    if len(atoms) < 2:
        return ()
    nxg = parent.to_networkx()
    lengths = dict(nx.all_pairs_shortest_path_length(nxg))
    return tuple(sorted(lengths[a].get(b, -1)
                        for a, b in itertools.combinations(sorted(atoms), 2)))


def labile_sites(mol: Chem.Mol, donor_type: str | None = None) -> list[DonorSite]:
    sites = [s for s in find_donor_sites(mol) if s.labile]
    if donor_type is not None:
        sites = [s for s in sites if s.donor_type == donor_type]
    return sites


def enumerate_protomers(
    mol: Chem.Mol,
    *,
    donor_type: str | None = None,
    max_deprotonations: int | None = None,
    multiplicity: int = 1,
) -> list[Protomer]:
    """Every distinct protonation state, deduplicated by identity.

    Selections that give the same L1 hash are the same structure — that is how the
    symmetry-equivalent choices collapse without anyone writing down the symmetry.
    """
    sites = labile_sites(mol, donor_type)
    limit = len(sites) if max_deprotonations is None else min(max_deprotonations, len(sites))
    # Equivalence of selections is decided on the PARENT: the un-deprotonated molecule,
    # whose symmetry is what makes two choices the same choice.
    parent = from_rdkit(mol, charge=0, multiplicity=multiplicity)
    site_atoms = [s.idx for s in sites]

    found: dict[str, dict] = {}
    for k in range(limit + 1):
        for selection in itertools.combinations(range(len(sites)), k):
            variant = mol
            if selection:
                variant, _ = deprotonate(mol, [sites[i] for i in selection])
            graph = from_rdkit(variant, charge=-k, multiplicity=multiplicity)
            key = l1_graph_hash(graph)
            entry = found.get(key)
            if entry is None:
                chosen_atoms = [sites[i].idx for i in selection]
                found[key] = {
                    "n": k, "charge": -k, "l1": key,
                    "configuration": configuration_key(parent, site_atoms, selection),
                    "separation": separation_profile(parent, chosen_atoms),
                    "selections": [selection], "representative": selection,
                    "mol": variant, "graph": graph,
                }
            else:
                entry["selections"].append(selection)

    return [
        Protomer(n_deprotonated=e["n"], charge=e["charge"], l1=e["l1"],
                 configuration=e["configuration"], separation=e["separation"],
                 selections=tuple(e["selections"]), representative=e["representative"],
                 mol=e["mol"], graph=e["graph"])
        for e in sorted(found.values(), key=lambda e: (e["n"], e["configuration"]))
    ]
