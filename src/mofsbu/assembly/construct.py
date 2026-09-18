"""Construction as an explicit decision tree (§6.7), and its deterministic replay.

A structure is not "built"; it is *chosen*.  §6.7 names the three kinds of choice and this
module is where they stop being prose:

* **Kind B — discrete, enumerable** (which donor, which vertex, which torsion well).  These
  are branches to walk, and `enumerate_constructions` walks them.  Each leaf carries the
  path that reached it, which is D11's conformer label: L3 identity is provenance-primary.
* **Kind C — inference** (CN, coordination geometry, protonation, spin).  When ambiguous
  these branch explicitly or raise — never default.  `resolve_geometry` is the one with
  teeth here: **CN 4 does not determine a polyhedron.**  Tetrahedral and square planar are
  both CN 4 and they are different structures with different chemistry, so a builder that
  quietly picks one has invented a result.  CN 6 *is* determined, and says so.
* **Kind A — stochastic** (the embedding, the azimuth search) is deliberately absent from
  this module.  It is noise to be deduplicated, it lives behind `seed`, and collapsing it
  is geometric clustering's job in S6, not the enumerator's.

**What this module refuses to do, and why that is not a gap.**  Enumerating the six
vertices of an octahedron produces six leaves for the first ligand, and they are the same
structure — one L1 with six routes.  It would be easy to collapse them here by hashing the
product, and it would be wrong: near-degenerate cis and trans *also* share an L1 and are
exactly what D10 insists must stay distinct.  Telling a symmetry duplicate from a real
branch needs L2 and the θ_geom clustering that S5/S6 deliver; until then the honest output
is every branch, a cap, and a `capped` flag that says the cap was reached.  A silently
truncated tree is worse than a big one.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterator, Sequence

import numpy as np

from mofsbu._types import AmbiguousSpecError
from mofsbu.assembly.choice import ChoiceVector, digest_of
from mofsbu.assembly.join import BuildingBlock, IncompatibleJoin, JoinResult, join
from mofsbu.geometry.placer import GEOMETRIES, site_vectors
from mofsbu.graph._types import TypedGraph
from mofsbu.graph.from_mol import from_rdkit
from mofsbu.sites.frames import BindingMode
from mofsbu.sites.model import Site, perceive, vacancy_sites
from mofsbu.sites.state import refresh_state

#: {coordination number: the polyhedra that have it}.  Derived from the placer's own
#: table rather than retyped, so a geometry added there cannot go missing here.
GEOMETRIES_FOR_CN: dict[int, tuple[str, ...]] = {
    cn: tuple(sorted(name for name, cns in GEOMETRIES.items() if cn in cns))
    for cn in sorted({cn for cns in GEOMETRIES.values() for cn in cns})
}


# ── Kind C: inference that branches or refuses ───────────────────────────────

def geometry_branches(cn: int, geometry: str | None = None) -> tuple[str, ...]:
    """The polyhedra a centre of this CN could adopt — the Kind-C branch list.

    Naming a geometry collapses the branch to one, which is a caller making the choice
    explicitly rather than the model guessing on their behalf.
    """
    if geometry is not None:
        if geometry not in GEOMETRIES:
            raise ValueError(
                f"unknown coordination geometry {geometry!r}; have {sorted(GEOMETRIES)}")
        if cn not in GEOMETRIES[geometry]:
            raise ValueError(
                f"{geometry} is not a CN-{cn} geometry; it has "
                f"{sorted(GEOMETRIES[geometry])}")
        return (geometry,)
    options = GEOMETRIES_FOR_CN.get(int(cn), ())
    if not options:
        raise ValueError(
            f"no coordination geometry for CN {cn}; have {sorted(GEOMETRIES_FOR_CN)}")
    return options


def resolve_geometry(cn: int, geometry: str | None = None) -> str:
    """The ONE geometry for this centre, or a refusal naming the alternatives.

    Ground rule 5 at its sharpest.  A CN-4 centre is tetrahedral or square planar, and the
    difference is not cosmetic: it changes the site frames, the accessible vertices, every
    downstream join, and — once S5 lands — the L2 tag.  Defaulting to either one would
    make the model's most consequential guess its most invisible.
    """
    options = geometry_branches(cn, geometry)
    if len(options) > 1:
        raise AmbiguousSpecError(
            f"CN {cn} does not determine a coordination geometry: it could be "
            f"{' or '.join(options)}, and those are different structures. Name one, or "
            f"call geometry_branches({cn}) and build both — this is a Kind-C branch "
            f"(design doc §6.7), not a missing default.")
    return options[0]


# ── blocks ───────────────────────────────────────────────────────────────────

def _with_state(graph: TypedGraph, sites: Sequence[Site], coords: np.ndarray,
                ) -> BuildingBlock:
    states = refresh_state(sites, coords, graph=graph,
                           symbols=[graph.label(i).element for i in graph.nodes()])
    return BuildingBlock(graph=graph, sites=tuple(sites), geometry=coords,
                         state={(s.atom_idx, s.slot): s for s in states})


def metal_block(symbol: str, oxidation_state: int, *, cn: int,
                geometry: str | None = None, spin_class: str = "hs",
                charge: int | None = None, multiplicity: int | None = None,
                ) -> BuildingBlock:
    """One metal centre with every coordination vertex empty, at the origin.

    Raises rather than guessing when `cn` admits more than one polyhedron — see
    `resolve_geometry`.  `metal_block_branches` is the same call in branching form.

    Multiplicity is DERIVED from the d-count and the spin class, never defaulted to 1: a
    high-spin Fe(III) is a sextet, and a module whose whole argument is "do not guess"
    handing back a singlet for it would be the loudest possible contradiction.

    The metal's per-atom `formal_charge` is **0** — charge is graph-level (D15), and the
    oxidation state travels in its own label field.  Writing it twice made this producer
    disagree with `from_rdkit` and `examples.py` about the same species (B12).
    """
    from mofsbu.energy.backends import spin_class_multiplicity

    chosen = resolve_geometry(cn, geometry)
    if multiplicity is None:
        multiplicity = spin_class_multiplicity(symbol, oxidation_state, spin_class)
    g = TypedGraph(charge=oxidation_state if charge is None else charge,
                   multiplicity=multiplicity, name=f"{symbol}({chosen})")
    g.add_atom(symbol, formal_charge=0,
               oxidation_state=oxidation_state, spin_class=spin_class)
    directions = [tuple(float(x) for x in v) for v in site_vectors(chosen, cn, 1.0)]
    return _with_state(g, vacancy_sites(0, (0.0, 0.0, 0.0), directions), np.zeros((1, 3)))


def metal_block_branches(symbol: str, oxidation_state: int, *, cn: int,
                         geometry: str | None = None,
                         **kw: Any) -> tuple[BuildingBlock, ...]:
    """Every polyhedron this CN admits, as separate blocks — Kind C enumerated."""
    return tuple(metal_block(symbol, oxidation_state, cn=cn, geometry=g, **kw)
                 for g in geometry_branches(cn, geometry))


def ligand_block(mol: Any, *, charge: int, name: str = "",
                 multiplicity: int = 1) -> BuildingBlock:
    """A perceived, framed ligand from an embedded RDKit molecule.

    Perception happens HERE, once, and never again: everything downstream inherits
    through atom maps (D5).  The molecule must already carry a conformer — a site without
    a frame cannot be joined, and embedding is a Kind-A step that belongs to whoever owns
    the seed.
    """
    from mofsbu.geometry.embed import coordinates

    if not mol.GetNumConformers():
        raise AmbiguousSpecError(
            f"ligand {name or '<unnamed>'} has no conformer, so its donors have no "
            "frames and nothing can be aligned to them. Embed it first "
            "(geometry.embed.embed_molecule) — which conformer you embedded is a Kind-A "
            "choice and travels with the seed, not with the identity.")
    graph = from_rdkit(mol, charge=charge, multiplicity=multiplicity, name=name)
    return _with_state(graph, perceive(mol), coordinates(mol))


# ── the tree ─────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Construction:
    """One leaf: the product, and the ordered path of choices that reached it."""

    block: BuildingBlock
    steps: tuple[dict[str, Any], ...]
    strain: float
    #: Atom count of the block this path started from.  The seed is always `join`'s first
    #: argument, so its atoms keep indices 0..seed_atoms-1 in every product along the
    #: path — which is what lets a recorded step address a vertex by index at any depth.
    seed_atoms: int = 0

    @property
    def choice_vector(self) -> dict[str, Any]:
        """The path AS a choice vector — what D11 makes the conformer's label."""
        return {"op": "construct", "n_steps": len(self.steps), "steps": list(self.steps)}

    @property
    def digest(self) -> str:
        return digest_of(self.choice_vector)

    @property
    def key(self) -> frozenset[str]:
        """Order-independent identity of the path: the SET of steps it took.

        Two orderings of the same joins reach the same product — the anchor never moves,
        so each ligand lands on its own vertex regardless of when — and this is what lets
        the enumerator emit that product once instead of once per permutation. It is the
        "one node, two routes" claim made at the level where the routes are generated.
        """
        return frozenset(digest_of(step) for step in self.steps)


@dataclass(frozen=True)
class BranchTree:
    """The leaves, plus what happened to everything that is not a leaf."""

    leaves: tuple[Construction, ...] = ()
    capped: bool = False
    #: Why a branch was refused, as `{reason: count}`.  A tree with no leaves has to be
    #: able to say why — silence was a real bug once (`test_a_run_that_queues_nothing_
    #: explains_itself`), and an enumerator is exactly where it would recur.
    refusals: dict[str, int] = field(default_factory=dict)

    def __iter__(self) -> Iterator[Construction]:
        return iter(self.leaves)

    def __len__(self) -> int:
        return len(self.leaves)

    def explain(self) -> str:
        if self.leaves:
            note = f"{len(self.leaves)} construction(s)"
            return note + (" (capped — there are more)" if self.capped else "")
        if not self.refusals:
            return ("no constructions: nothing was offered to join (no open donors, or "
                    "no open vertices to put them on)")
        worst = sorted(self.refusals.items(), key=lambda kv: -kv[1])
        return "no constructions; every branch was refused: " + "; ".join(
            f"{reason} (x{n})" for reason, n in worst)


def _candidate_steps(block: BuildingBlock, partners: Sequence[BuildingBlock],
                     modes: Sequence[str],
                     ) -> tuple[list[tuple[int, Site, Site, str, int]], dict[str, int]]:
    """Every (partner, donor, vertex, mode, well) this block could try, in a fixed order.

    **Both directions.**  A partner's donor can go on this block's vertex — adding a
    ligand — and this block's donor can go on a PARTNER's vertex, which is adding a metal.
    Only offering the first would make polynuclear growth unreachable through the
    enumerator while looking like it had simply found nothing, and "the paddlewheel branch
    silently does not exist" is a worse answer than the M6 refusal that direction actually
    produces.

    Deterministic ordering is not tidiness: the output order is what `max_products`
    truncates, so an unstable order would make a capped tree depend on dict iteration.

    Skipped candidates are COUNTED, not dropped.  A mode the donor does not offer is a
    real reason a tree came back empty, and an enumerator that filtered it out silently
    would hand back "nothing was offered to join" — true, unhelpful, and the exact shape
    of the silence `test_a_run_that_queues_nothing_explains_itself` exists to prevent.
    """
    out: list[tuple[int, Site, Site, str, int]] = []
    skipped: dict[str, int] = {}
    own_vacancies = sorted(block.open_vacancies(), key=lambda s: (s.atom_idx, s.slot))
    own_donors = sorted(block.open_donors(), key=lambda s: s.atom_idx)
    for p, partner in enumerate(partners):
        # (site on THIS block, site on the partner) — the order `join` expects, so the
        # growing block stays its first argument and keeps its atom numbering.
        pairs = [(v, d) for d in sorted(partner.open_donors(), key=lambda s: s.atom_idx)
                 for v in own_vacancies]
        pairs += [(d, v) for d in own_donors
                  for v in sorted(partner.open_vacancies(), key=lambda s: (s.atom_idx, s.slot))]
        for site_a, site_b in pairs:
            donor = site_b if site_a.is_vacancy else site_a
            for mode in modes:
                if mode not in donor.binding_modes:
                    reason = (f"{donor.donor_type} does not bind {mode!r}; it offers "
                              f"{', '.join(donor.binding_modes)}")
                    skipped[reason] = skipped.get(reason, 0) + 1
                    continue
                for well in range(len(torsion_wells_for(donor, mode))):
                    out.append((p, site_a, site_b, mode, well))
    return out, skipped


def torsion_wells_for(donor: Site, mode: str) -> tuple[float, ...]:
    """The wells this donor branches over in this mode.  `TORSION_FREE` gives exactly one.

    That single-well answer is the declared guard against conformer explosion, and it is
    why `aqua` never doubles the tree.
    """
    from mofsbu.sites.frames import torsion_wells

    return torsion_wells(donor.donor_type, BindingMode(mode))


def enumerate_constructions(
    seed: BuildingBlock,
    partners: Sequence[BuildingBlock],
    *,
    degree: int = 1,
    modes: Sequence[str] = (BindingMode.MONODENTATE.value,),
    max_products: int = 1000,
    with_geometry: bool = True,
) -> BranchTree:
    """Walk the Kind-B branch tree `degree` additions deep.

    Breadth first, so a cap truncates the deepest level rather than an arbitrary subtree,
    and so `degree=2` always contains every `degree=1` product's children.

    The growing block is always the FIRST argument to `join`, which makes its atom
    numbering stable across the whole path: a vertex recorded at step 1 still means that
    vertex at step 3, and a replayed path can therefore address sites by index without
    carrying a chain of atom maps.
    """
    if degree < 1:
        raise ValueError(f"degree {degree} builds nothing; the smallest construction is 1")
    if degree > 1 and not with_geometry:
        raise ValueError(
            "degree > 1 needs geometry: the next step's candidates are the product's OPEN "
            "sites, and openness is a per-geometry fact (`sites.state`). A graph-only "
            "product cannot say which of its vertices are still reachable, and treating "
            "unknown as open is how a second ligand gets joined onto a buried donor.")
    refusals: dict[str, int] = {}
    n_seed = len(seed.graph)
    frontier = [Construction(seed, (), 0.0, n_seed)]
    leaves: dict[frozenset[str], Construction] = {}
    capped = False

    for _level in range(degree):
        nxt: dict[frozenset[str], Construction] = {}
        for node in frontier:
            candidates, skipped = _candidate_steps(node.block, partners, modes)
            for reason, n in skipped.items():
                refusals[reason] = refusals.get(reason, 0) + n
            for p, site_a, site_b, mode, well in candidates:
                try:
                    result = join(node.block, partners[p], site_a, site_b,
                                  mode=mode, torsion_well=well,
                                  with_geometry=with_geometry)
                except IncompatibleJoin as exc:
                    refusals[exc.verdict.reason] = refusals.get(exc.verdict.reason, 0) + 1
                    continue
                step = dict(result.choice_vector, partner=p)
                child = Construction(result.block, node.steps + (step,),
                                     max(node.strain, result.strain), n_seed)
                if child.key in nxt:
                    continue                     # same joins, different order
                if len(leaves) + len(nxt) >= max_products:
                    capped = True
                    break
                nxt[child.key] = child
            if capped:
                break
        leaves.update(nxt)
        frontier = list(nxt.values())
        if capped or not frontier:
            break

    ordered = tuple(sorted(leaves.values(), key=lambda c: (len(c.steps), c.digest)))
    return BranchTree(ordered, capped=capped, refusals=refusals)


def construct(seed: BuildingBlock, partners: Sequence[BuildingBlock],
              path: Sequence[dict[str, Any]] | Construction | dict[str, Any], *,
              with_geometry: bool = True) -> Construction:
    """Rebuild a construction from its recorded path.  A pure function of its inputs.

    This is the half of D13 that makes a stored conformer a regenerable object rather than
    frozen coordinates: hand back the same seed, the same partners and the same path, and
    the same structure comes out — at whatever fidelity you now want to spend on it.

    Every step is addressed exactly as it was recorded, and a step that under-determines
    its join raises rather than filling the gap in.  A replay that quietly substituted a
    default would return a *different structure wearing the original's provenance*, which
    is the single worst failure available to this layer.
    """
    steps = _steps_of(path)
    block, strain, taken = seed, 0.0, []
    n_seed = len(seed.graph)
    for n, step in enumerate(steps, start=1):
        donor_spec, vacancy_spec = step.get("donor"), step.get("vacancy")
        if not isinstance(donor_spec, dict) or not isinstance(vacancy_spec, dict):
            raise AmbiguousSpecError(
                f"step {n} records no donor/vacancy pair, so there is nothing to replay. "
                f"A path is the exact sequence of joins, not a description of one.")
        try:
            p = int(step["partner"])
            partner = partners[p]
        except (KeyError, TypeError, ValueError):
            raise AmbiguousSpecError(
                f"step {n} does not say which partner it joined; with {len(partners)} "
                f"available there is no way to choose one that is not a guess.") from None
        except IndexError:
            raise AmbiguousSpecError(
                f"step {n} names partner {step['partner']}, but only {len(partners)} were "
                f"supplied. Replay needs the SAME partners, in the same order.") from None

        # Which side held which role is recorded, not re-derived: `order` says whether the
        # growing block brought the donor (adding a metal) or the vertex (adding a ligand),
        # and guessing it from the blocks would silently repair a corrupted path into a
        # different structure wearing the original's provenance.
        if step.get("order") == "donor-block-first":
            site_a = _find_site(block, donor_spec.get("atom"), 0, f"donor (step {n})")
            site_b = _find_site(partner, vacancy_spec.get("atom"),
                                vacancy_spec.get("slot"), f"vacancy (step {n})")
        else:
            site_a = _find_site(block, vacancy_spec.get("atom"),
                                vacancy_spec.get("slot"), f"vacancy (step {n})")
            site_b = _find_site(partner, donor_spec.get("atom"), 0, f"donor (step {n})")
        # `lone_pair` defaults to 0 rather than raising on an older path, and that is safe
        # for the one reason a replay default ever is: 0 IS what those paths did — it is
        # the stored `frame`, and the lobes were not selectable when they were written.
        result = join(block, partner, site_a, site_b,
                      mode=step.get("mode", BindingMode.MONODENTATE.value),
                      torsion_well=int(step.get("torsion_well", 0)),
                      lone_pair=int(donor_spec.get("lone_pair", 0)),
                      seed=int(step.get("seed", 0)), with_geometry=with_geometry)
        block, strain = result.block, max(strain, result.strain)
        taken.append(dict(result.choice_vector, partner=p))
    return Construction(block, tuple(taken), strain, n_seed)


def _steps_of(path: Sequence[dict] | Construction | dict) -> list[dict]:
    if isinstance(path, Construction):
        return list(path.steps)
    if isinstance(path, dict):
        steps = path.get("steps")
        if steps is None:
            raise AmbiguousSpecError(
                "a construction choice vector must carry its `steps`; without them the "
                "digest identifies a path nobody can walk.")
        return list(steps)
    return list(path)


def _find_site(block: BuildingBlock, atom: Any, slot: Any, what: str) -> Site:
    """The site a recorded step refers to, or a refusal that says what it looked for."""
    if atom is None:
        raise AmbiguousSpecError(f"{what} records no atom index")
    for site in block.sites:
        if site.atom_idx == atom and (slot is None or site.slot == slot):
            return site
    raise AmbiguousSpecError(
        f"{what}: no site at atom {atom}"
        + (f" slot {slot}" if slot is not None else "")
        + ". Replay addresses sites by the index they had when the path was recorded, so "
        "this means the seed or the partners are not the ones it was recorded against.")


def constructions_as_joins(tree: BranchTree) -> list[JoinResult]:
    """`BranchTree` in the shape `grow`'s settled signature promised.

    The path, not the last join, is what a multi-step growth means — so the choice vector
    carried here is the whole construction's and `strain` is the worst step on it, the
    same "the path is as viable as its worst point" reading M8 will use for barriers.
    """
    return [JoinResult(block=c.block, atom_map={i: i for i in range(c.seed_atoms)},
                       choice_vector=c.choice_vector, strain=c.strain)
            for c in tree.leaves]


def replayable(construction: Construction) -> ChoiceVector:
    """The construction's choice vector, keyed — what a `geometries` row would store."""
    return ChoiceVector(construction.choice_vector)
