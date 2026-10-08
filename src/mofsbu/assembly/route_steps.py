"""Execute a `spec.RouteSpec`: an explicit, label-addressed sequence of assembly steps.

A route is how a polynuclear node is written down as data (D20: a node is what a sequence
of joins produces).  Every step goes through the same assembly calls a hand-written build
uses — `join`, `join_bridge`, `bridge_compatible` — so a route replays bit-identically
from its spec, and nothing here places an atom the assembly layer did not.

Labels, never indices.  A step names atoms by label and the executor carries every label
through each join's atom maps (`atom_map` for the first argument, `partner_atom_map` for
the second), because a join renumbers both blocks:

  ``cu1``        a metal's atom (``metal`` step)
  ``f1``         every atom of a ligand (``ligand`` step); ``f1.d0``, ``f1.d1`` its open
                 donors in atom order
  ``cu2.cap``    every atom of the co-ligands a ``cap`` step put on ``cu2``
  ``<reserve>``  vertices held open by a ``reserve`` step, for a later ``bridge``

Steps (required keys in `spec.ROUTE_OPS`).  The first three are the ground-up path: every
reagent is built and STORED by the pipeline's own `place` / `ligand` code, and a join only
ever acts on blocks read back from the registry by `load` (`persist.load_block`):

  sphere      {as, metal: {symbol, oxidation_state?, spin_class?}, cn, geometry,
               ligands?: [{molecule, count?, donor?}], co_ligands?, co_ligand?}
              one coordination sphere through `runner._build_sphere`, stored as `place`
              stores it; empty vertices are the ones the placer leaves (a cis set at >= 2)
  free_ligand {as, molecule}                     the molecule as written, stored as the
                                                 `ligand` task stores it
  load        {as, species | block_id, geometry?}  a stored structure as a block; labels
              `as`, `as.m<k>` metals, `as.d<k>` open donors, `as.L<k>` ligand fragments,
              `as.<formula>` (e.g. `as.H2O`), all in canonical order

and the in-memory construction steps:

  metal    {as, symbol, cn, geometry?, oxidation_state?, spin_class?}
  ligand   {as, molecule, seed?}
  join     {donor, onto, slot?, lone_pair?, torsion?}
                                                donor's block onto a metal's vacancy;
                                                `lone_pair` an index or {toward|away: label},
                                                the lobe whose axis points most at/from it;
                                                `torsion` a well index or {bring, toward}: the
                                                well that brings one label nearest the other
  reserve  {as, for, metals: [m1, m2]}          the vertex pair `bridge_compatible` ranks
                                                best for ligand `for`
           {as, metal, nearest}                 the vacancy on `metal` nearest an atom
  cap      {block, molecule, keep?, seed?}      one co-ligand on every open vertex of the
                                                block holding `block`, except `keep`
  bridge   {ligand, vacancies, donors?, lone_pairs?}
  release  {ligand, as, multiplicity?}          remove the ligand fragment holding the
                                                `ligand` atoms (`assembly.release`); its
                                                vertex opens, the fragment becomes block
                                                `as` — store it and name it in `leaving`
  store    {as, block, tags?, from?, leaving?, watch?}
                                                persist the block holding `block`; `from`
                                                are reagents, `leaving` released species
  perturb  {as, of, move, from, to, by}         a start geometry of stored `of` with the
                                                `move` atoms shifted `by` Å along the
                                                centroid(from) -> centroid(to) direction

A `perturb` geometry is a starting point for a relaxation, not a construct: it is stored on
its structure with the displacement in its choice vector and its clash QC, and it is never
offered as a built geometry.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from mofsbu._types import MofsbuError

__all__ = ["RouteError", "RouteResult", "StoredSpecies", "execute_route"]


class RouteError(MofsbuError):
    """A route step that cannot be carried out as written.  Says which step and why."""


@dataclass
class StoredSpecies:
    """One `store` or `perturb` result: where it went and how to read its rows back."""

    name: str
    structure_id: int
    geometry_id: int
    kind: str                                   # "built" | "perturbed"
    rows: dict[str, list[int]]                  # label -> xyz row indices of this geometry
    watch: dict[str, list[str]] = field(default_factory=dict)
    reserved: dict[str, list[list[Any]]] = field(default_factory=dict)
    perturbation: dict[str, Any] | None = None
    created: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "structure_id": self.structure_id,
                "geometry_id": self.geometry_id, "kind": self.kind, "rows": self.rows,
                "watch": self.watch, "reserved": self.reserved,
                "perturbation": self.perturbation, "created": self.created}


@dataclass
class RouteResult:
    route: str
    stored: list[StoredSpecies]


class _Workspace:
    """Blocks under construction and where every label currently points."""

    def __init__(self, route: str) -> None:
        self.route = route
        self.blocks: dict[int, Any] = {}
        self.next_id = 0
        self.where: dict[str, tuple[int, list[int]]] = {}      # label -> (block, atoms)
        self.reserved: dict[str, list[tuple[str, int]]] = {}   # name -> [(metal label, slot)]
        self.cv: list[dict[str, Any]] = []
        self.stored: dict[str, StoredSpecies] = {}
        self.blocks_of_stored: dict[str, Any] = {}
        self.loaded: dict[str, int] = {}                        # load label -> structure id

    def add(self, block: Any) -> int:
        bid = self.next_id
        self.next_id += 1
        self.blocks[bid] = block
        return bid

    def atoms(self, ref: str, step: int) -> tuple[int, list[int]]:
        if ref not in self.where:
            raise RouteError(f"route {self.route!r} step {step}: no label {ref!r}; have "
                             f"{sorted(self.where)}")
        return self.where[ref]

    def atom(self, ref: str, step: int) -> tuple[int, int]:
        bid, atoms = self.atoms(ref, step)
        if len(atoms) != 1:
            raise RouteError(f"route {self.route!r} step {step}: {ref!r} names "
                             f"{len(atoms)} atoms where one is needed")
        return bid, atoms[0]

    def merge(self, a: int, b: int, result: Any) -> int:
        """Replace blocks `a` and `b` by a join's product, carrying every label across."""
        new = self.add(result.block)
        for label, (bid, atoms) in list(self.where.items()):
            if bid == a:
                self.where[label] = (new, [result.atom_map[i] for i in atoms])
            elif bid == b:
                self.where[label] = (new, [result.partner_atom_map[i] for i in atoms])
        del self.blocks[a], self.blocks[b]
        cv = getattr(result.choice_vector, "data", result.choice_vector)
        self.cv.append(dict(cv) if cv else {})
        return new

    def vacancy(self, metal: str, slot: int | None, step: int):
        bid, m = self.atom(metal, step)
        held = {(self.where[lab][1][0], s) for v in self.reserved.values() for lab, s in v}
        opts = [v for v in self.blocks[bid].open_vacancies() if v.atom_idx == m]
        if slot is None:
            opts = [v for v in opts if (v.atom_idx, v.slot) not in held]
            if not opts:
                raise RouteError(f"route {self.route!r} step {step}: {metal!r} has no open "
                                 f"vertex that is not reserved")
            return bid, opts[0]
        hit = [v for v in opts if v.slot == int(slot)]
        if not hit:
            raise RouteError(f"route {self.route!r} step {step}: {metal!r} vertex {slot} is "
                             f"not open; open: {sorted(v.slot for v in opts)}")
        return bid, hit[0]

    def rows(self, bid: int) -> dict[str, list[int]]:
        """Labels on block `bid` as xyz rows (graph node order, as `to_xyz` writes)."""
        order = {n: k for k, n in enumerate(self.blocks[bid].graph.nodes())}
        return {lab: [order[i] for i in atoms]
                for lab, (b, atoms) in self.where.items() if b == bid}


def _molecule(spec: Any, name: str):
    return next(m for m in spec.molecules if m.name == name)


def _ligand(spec: Any, name: str, seed: int):
    from rdkit import Chem

    from mofsbu.assembly.construct import ligand_block
    from mofsbu.geometry.embed import embed_molecule
    from mofsbu.graph.from_mol import mol_from_smiles

    m = _molecule(spec, name)
    mol = mol_from_smiles(m.smiles)
    charge = sum(a.GetFormalCharge() for a in Chem.AddHs(mol).GetAtoms())
    return ligand_block(embed_molecule(mol, seed=seed), charge=charge, name=name,
                        multiplicity=m.multiplicity)


def labels_for(reg: Any, structure_id: int, block: Any, prefix: str) -> dict[str, list[int]]:
    """Labels for a block read from the registry, in CANONICAL order so they are portable.

    ``p`` every atom; ``p.m0``… the metals; ``p.d0``… the OPEN donors; ``p.L0``… each ligand
    fragment; ``p.<formula>`` every ligand fragment of that formula (``p.H2O``: the waters).
    """
    from mofsbu.identity.keys import hill_formula
    from mofsbu.registry.api import canonical_map

    cmap = canonical_map(reg, structure_id)
    g = block.graph
    out: dict[str, list[int]] = {prefix: list(g.nodes())}
    for k, m in enumerate(sorted(g.metals(), key=cmap.get)):
        out[f"{prefix}.m{k}"] = [m]
    for k, d in enumerate(sorted(block.open_donors(), key=lambda d: cmap[d.atom_idx])):
        out[f"{prefix}.d{k}"] = [d.atom_idx]
    frags = sorted((tuple(f) for f in g.ligand_fragments()), key=lambda f: min(cmap[a] for a in f))
    for k, frag in enumerate(frags):
        out[f"{prefix}.L{k}"] = list(frag)
        counts: dict[str, int] = {}
        for a in frag:
            el = g.label(a).element
            counts[el] = counts.get(el, 0) + 1
        out.setdefault(f"{prefix}.{hill_formula(counts)}", []).extend(frag)
    return out


def _lobe(ws: _Workspace, site: Any, spec_lp: Any, step: int) -> int | None:
    """`lone_pair` as written: an index, or {toward|away: label} resolved to one."""
    from mofsbu.sites.model import frame_lobes

    if spec_lp is None or isinstance(spec_lp, int):
        return spec_lp
    if not isinstance(spec_lp, dict) or len(spec_lp) != 1 or \
            next(iter(spec_lp)) not in ("toward", "away"):
        raise RouteError(f"route {ws.route!r} step {step}: lone_pair must be an index or "
                         f"{{'toward'|'away': label}}, got {spec_lp!r}")
    (how, ref), = spec_lp.items()
    bid, atoms = ws.atoms(ref, step)
    blk = ws.blocks[bid]
    order = {n: k for k, n in enumerate(blk.graph.nodes())}
    xyz = np.asarray(blk.geometry, dtype=float)
    lobes = frame_lobes(site.frame)
    if not lobes:
        raise RouteError(f"route {ws.route!r} step {step}: the donor has no frame to choose "
                         f"a lobe from")
    target = xyz[[order[a] for a in atoms]].mean(axis=0) - np.asarray(lobes[0]["origin"])
    target /= np.linalg.norm(target)
    dots = [float(np.dot(np.asarray(lb["axis"], dtype=float), target)) for lb in lobes]
    return int(np.argmax(dots) if how == "toward" else np.argmin(dots))


def _torsion_well(ws: _Workspace, a: int, b: int, donor: Any, vac: Any, kw: dict,
                  spec_t: Any, step: int) -> int:
    """`torsion` as {bring: label, toward: label}: the donor's well that brings one near the other.

    A well index is not a side of anything — a vacancy's torsion reference is not consistent
    across vertices (B27) — so a route that means "the same face as that ligand" says so, and
    this picks the well that puts the `bring` atoms closest to the `toward` atoms.  Ties go to
    the lower well.
    """
    from mofsbu.assembly.join import join
    from mofsbu.sites.frames import BindingMode, torsion_wells

    if not isinstance(spec_t, dict) or set(spec_t) != {"bring", "toward"}:
        raise RouteError(f"route {ws.route!r} step {step}: torsion must be a well index or "
                         f"{{'bring': label, 'toward': label}}, got {spec_t!r}")
    (ba, bring), (ta, toward) = ws.atoms(spec_t["bring"], step), ws.atoms(spec_t["toward"], step)
    if {ba, ta} - {a, b}:
        raise RouteError(f"route {ws.route!r} step {step}: torsion labels must be on the two "
                         f"blocks being joined")
    best = None
    for w in range(len(torsion_wells(donor.donor_type, BindingMode("mono")))):
        res = join(ws.blocks[a], ws.blocks[b], donor, vac, **kw, torsion_well=w)
        moved = {a: res.atom_map, b: res.partner_atom_map}
        order = {n: k for k, n in enumerate(res.block.graph.nodes())}
        xyz = np.asarray(res.block.geometry, dtype=float)
        p = xyz[[order[moved[ba][x]] for x in bring]].mean(axis=0)
        q = xyz[[order[moved[ta][x]] for x in toward]].mean(axis=0)
        d = round(float(np.linalg.norm(p - q)), 6)
        if best is None or d < best[0]:
            best = (d, w)
    return best[1]


def _release_target(ws: _Workspace, spec_lig: Any, step: int) -> tuple[int, list[int]]:
    """A `release` step's ligand: a label, or {on: metal, nearest: label}.

    The second form picks, among the ligand fragments bound to `on`, the one whose donor
    is nearest the centroid of `nearest` (which may be a released block's atoms — its
    coordinates stay in the same frame).  Ties go to the lowest atom index.
    """
    from mofsbu.graph._types import EdgeType

    if isinstance(spec_lig, str):
        return ws.atoms(spec_lig, step)
    if not isinstance(spec_lig, dict) or not {"on", "nearest"} <= set(spec_lig) <= {
            "on", "nearest", "formula"}:
        raise RouteError(f"route {ws.route!r} step {step}: release ligand must be a label "
                         f"or {{'on': metal, 'nearest': label, 'formula'?: str}}, got "
                         f"{spec_lig!r}")
    from mofsbu.identity.keys import hill_formula
    bid, metal = ws.atom(spec_lig["on"], step)
    blk = ws.blocks[bid]
    g = blk.graph
    xyz = np.asarray(blk.geometry, dtype=float)
    order = {n: k for k, n in enumerate(g.nodes())}
    tb, tatoms = ws.atoms(spec_lig["nearest"], step)
    txyz = np.asarray(ws.blocks[tb].geometry, dtype=float)
    torder = {n: k for k, n in enumerate(ws.blocks[tb].graph.nodes())}
    target = txyz[[torder[a] for a in tatoms]].mean(axis=0)
    frag_of = {a: f for f in g.ligand_fragments() for a in f}

    def formula(frag) -> str:
        counts: dict[str, int] = {}
        for a in frag:
            counts[g.label(a).element] = counts.get(g.label(a).element, 0) + 1
        return hill_formula(counts)

    want = spec_lig.get("formula")
    donors = sorted(d for d in g.neighbors(metal, EdgeType.DATIVE)
                    if not (tb == bid and d in tatoms)
                    and (want is None or formula(frag_of[d]) == want))
    if not donors:
        raise RouteError(f"route {ws.route!r} step {step}: {spec_lig['on']!r} has no bound "
                         f"ligand{' of formula ' + want if want else ''} to release")
    best = min(donors, key=lambda d: (round(float(np.linalg.norm(xyz[order[d]] - target)), 6), d))
    return bid, [best]


def _donor_site(ws: _Workspace, ref: str, step: int):
    bid, atom = ws.atom(ref, step)
    site = next((s for s in ws.blocks[bid].open_donors() if s.atom_idx == atom), None)
    if site is None:
        raise RouteError(f"route {ws.route!r} step {step}: {ref!r} is not an open donor")
    return bid, site


def _do(ws: _Workspace, spec: Any, i: int, s: dict[str, Any], reg: Any,
        stored_cb: Any) -> None:
    from mofsbu.assembly.construct import metal_block
    from mofsbu.assembly.join import bridge_compatible, join, join_bridge

    op = s["op"]
    if op in ("sphere", "free_ligand"):
        # Built and STORED by the runner's own `place` / `ligand` paths; usable in a later
        # step only after a `load`, so every join operand is read back from the registry.
        ws.stored[s["as"]] = stored_cb(ws, i, s, None)
    elif op == "load":
        from mofsbu.assembly.persist import load_block

        if reg is None:
            raise RouteError(f"route {ws.route!r} step {i}: load needs a registry")
        if "species" in s:
            sp = ws.stored.get(s["species"])
            if sp is None or sp.kind != "built":
                raise RouteError(f"route {ws.route!r} step {i}: no stored species "
                                 f"{s['species']!r} to load; have "
                                 f"{sorted(n for n, x in ws.stored.items() if x.kind == 'built')}")
            sid, gid = sp.structure_id, sp.geometry_id
        else:
            row = reg.conn.execute("SELECT id FROM structures WHERE block_id=?",
                                   (s["block_id"],)).fetchone()
            if row is None:
                raise RouteError(f"route {ws.route!r} step {i}: no structure "
                                 f"{s['block_id']!r} in this registry")
            sid, gid = int(row["id"]), s.get("geometry")
        blk = load_block(reg, sid, gid)
        bid = ws.add(blk)
        for lab, atoms in labels_for(reg, sid, blk, s["as"]).items():
            ws.where[lab] = (bid, atoms)
        ws.loaded[s["as"]] = sid
    elif op == "release":
        from mofsbu.assembly.release import release

        bid, atoms = _release_target(ws, s["ligand"], i)
        res = release(ws.blocks[bid], atoms, multiplicity=int(s.get("multiplicity", 1)))
        new, lig = ws.add(res.block), ws.add(res.ligand)
        for label, (b, ats) in list(ws.where.items()):
            if b != bid:
                continue
            rest = [res.atom_map[a] for a in ats if a in res.atom_map]
            gone = [res.ligand_map[a] for a in ats if a in res.ligand_map]
            if rest:                        # a label spanning both keeps what remains
                ws.where[label] = (new, rest)
            else:                           # wholly on the ligand: it leaves with it
                ws.where[label] = (lig, gone)
        ws.where[s["as"]] = (lig, list(res.ligand.graph.nodes()))
        del ws.blocks[bid]
        ws.cv.append(res.choice_vector)
    elif op == "metal":
        blk = metal_block(s["symbol"], int(s.get("oxidation_state", 2)), cn=int(s["cn"]),
                          geometry=s.get("geometry"), spin_class=s.get("spin_class", "hs"))
        ws.where[s["as"]] = (ws.add(blk), [blk.graph.metals()[0]])
    elif op == "ligand":
        blk = _ligand(spec, s["molecule"], int(s.get("seed", 7)))
        bid = ws.add(blk)
        ws.where[s["as"]] = (bid, list(blk.graph.nodes()))
        for k, d in enumerate(sorted(blk.open_donors(), key=lambda d: d.atom_idx)):
            ws.where[f"{s['as']}.d{k}"] = (bid, [d.atom_idx])
    elif op == "join":
        a, donor = _donor_site(ws, s["donor"], i)
        b, vac = ws.vacancy(s["onto"], s.get("slot"), i)
        if a == b:
            raise RouteError(f"route {ws.route!r} step {i}: donor and vertex are on one "
                             f"block; that is a chelate or a bridge, not a join")
        lobe = _lobe(ws, donor, s.get("lone_pair"), i)
        kw = {} if lobe is None else {"lone_pair": int(lobe)}
        if isinstance(s.get("torsion"), int):
            kw["torsion_well"] = int(s["torsion"])
        elif s.get("torsion") is not None:
            kw["torsion_well"] = _torsion_well(ws, a, b, donor, vac, kw, s["torsion"], i)
        ws.merge(a, b, join(ws.blocks[a], ws.blocks[b], donor, vac, **kw))
    elif op == "reserve":
        if "for" in s:
            lig_bid, _ = ws.atoms(s["for"], i)
            donors = sorted(ws.blocks[lig_bid].open_donors(), key=lambda d: d.atom_idx)
            m1, m2 = s["metals"]
            b1, a1 = ws.atom(m1, i)
            b2, a2 = ws.atom(m2, i)
            if b1 != b2:
                raise RouteError(f"route {ws.route!r} step {i}: {m1!r} and {m2!r} are not "
                                 f"in one block yet")
            vac = ws.blocks[b1].open_vacancies()
            pairs = [(x, y) for x in vac if x.atom_idx == a1 for y in vac if y.atom_idx == a2]
            if not pairs:
                raise RouteError(f"route {ws.route!r} step {i}: no open vertex pair across "
                                 f"{m1!r} and {m2!r}")
            elems = [ws.blocks[lig_bid].graph.label(d.atom_idx).element for d in donors[:2]]
            partner = ws.blocks[b1].graph.label(a1).element
            best = min(pairs, key=lambda p: bridge_compatible(
                donors[:2], list(p), partner=partner, donor_elements=elems).strain)
            ws.reserved[s["as"]] = [(m1, best[0].slot), (m2, best[1].slot)]
        elif "metal" in s and "nearest" in s:
            bid, m = ws.atom(s["metal"], i)
            tb, target = ws.atoms(s["nearest"], i)
            blk = ws.blocks[bid]
            xyz = np.asarray(blk.geometry, dtype=float)
            order = {n: k for k, n in enumerate(blk.graph.nodes())}
            p = xyz[[order[t] for t in target]].mean(axis=0)
            opts = [v for v in blk.open_vacancies() if v.atom_idx == m]
            axis = lambda v: np.asarray((v.frame or {})["axis"], dtype=float)  # noqa: E731
            best = min(opts, key=lambda v: float(np.linalg.norm(xyz[order[m]] + axis(v) - p)))
            ws.reserved[s["as"]] = [(s["metal"], best.slot)]
        else:
            raise RouteError(f"route {ws.route!r} step {i}: reserve needs either "
                             f"'for' + 'metals' or 'metal' + 'nearest'")
    elif op == "cap":
        bid, _ = ws.atoms(s["block"], i)
        keep = {(lab, slot) for name in s.get("keep", ()) for lab, slot in ws.reserved[name]}
        seed = int(s.get("seed", 11))
        while True:
            blk = ws.blocks[bid]
            metals = set(blk.graph.metals())
            metal_label = {atoms[0]: lab for lab, (b, atoms) in ws.where.items()
                           if b == bid and "." not in lab and len(atoms) == 1
                           and atoms[0] in metals}
            held = {(ws.where[lab][1][0], slot) for lab, slot in keep}
            open_v = [v for v in blk.open_vacancies() if (v.atom_idx, v.slot) not in held]
            if not open_v:
                break
            v = open_v[0]
            lig = _ligand(spec, s["molecule"], seed)
            seed += 1
            lb = ws.add(lig)
            donor = sorted(lig.open_donors(), key=lambda d: d.atom_idx)[0]
            tag = f"{metal_label.get(v.atom_idx, 'metal')}.cap"
            tmp = f"__cap{lb}"
            ws.where[tmp] = (lb, list(lig.graph.nodes()))
            bid = ws.merge(lb, bid, join(lig, blk, donor, v))
            prev = ws.where.get(tag, (bid, []))[1]
            ws.where[tag] = (bid, prev + ws.where.pop(tmp)[1])
    elif op == "bridge":
        lig_bid, _ = ws.atoms(s["ligand"], i)
        names = s.get("donors") or [f"{s['ligand']}.d0", f"{s['ligand']}.d1"]
        donors = [_donor_site(ws, n, i)[1] for n in names]
        held = ws.reserved.get(s["vacancies"])
        if held is None:
            raise RouteError(f"route {ws.route!r} step {i}: no reservation "
                             f"{s['vacancies']!r}")
        target = [ws.vacancy(lab, slot, i) for lab, slot in held]
        bid = target[0][0]
        if any(t[0] != bid for t in target):
            raise RouteError(f"route {ws.route!r} step {i}: reserved vertices span blocks")
        del ws.reserved[s["vacancies"]]
        lps = tuple(int(x) for x in s.get("lone_pairs", (0, 0)))
        ws.merge(lig_bid, bid, join_bridge(ws.blocks[lig_bid], ws.blocks[bid], donors,
                                           [t[1] for t in target], lone_pairs=lps))
    elif op == "store":
        bid, _ = ws.atoms(s["block"], i)
        sp = stored_cb(ws, i, s, bid)
        ws.stored[s["as"]] = sp
        ws.blocks_of_stored[s["as"]] = (ws.blocks[bid], ws.rows(bid))
    elif op == "perturb":
        if s["of"] not in ws.stored:
            raise RouteError(f"route {ws.route!r} step {i}: {s['of']!r} is not stored yet")
        ws.stored[s["as"]] = stored_cb(ws, i, s, None)


def perturbed_coords(block: Any, rows: dict[str, list[int]], s: dict[str, Any],
                     route: str = "", step: int = -1) -> np.ndarray:
    """`block`'s coordinates with the `move` atoms shifted along centroid(from)->centroid(to)."""
    xyz = np.array(block.geometry, dtype=float)

    def pick(refs):
        out = []
        for r in ([refs] if isinstance(refs, str) else refs):
            if r not in rows:
                raise RouteError(f"route {route!r} step {step}: no label {r!r} on "
                                 f"{s['of']!r}")
            out += rows[r]
        return sorted(set(out))

    move, a, b = pick(s["move"]), pick(s["from"]), pick(s["to"])
    d = xyz[b].mean(axis=0) - xyz[a].mean(axis=0)
    n = float(np.linalg.norm(d))
    if n < 1e-6:
        raise RouteError(f"route {route!r} step {step}: 'from' and 'to' coincide")
    xyz[move] += float(s["by"]) * d / n
    return xyz


def execute_route(spec: Any, route: Any, stored_cb: Any, reg: Any = None) -> RouteResult:
    """Run every step of `route`; `stored_cb(ws, step, step_dict, block_id)` persists."""
    ws = _Workspace(route.name)
    for i, step in enumerate(route.steps):
        _do(ws, spec, i, step, reg, stored_cb)
    return RouteResult(route.name, list(ws.stored.values()))
