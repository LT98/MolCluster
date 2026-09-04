"""The typed molecular graph — the identity substrate (design doc D3, D7, D12).

Nodes are ATOMS (hydrogens explicit).  Edges are typed bonds.  Everything that
distinguishes one structure from another must be visible here, and nothing that
does NOT distinguish structures may be visible here.

Two deliberate exclusions from identity, both to stop chemically-identical inputs
hashing differently (see PLAN §0 rule 2):

* **Bond order is stored but NOT hashed.**  Otherwise the two Kekule forms of
  benzene, and the C=O / C-O resonance forms of a carboxylate, would be different
  structures.  For a carboxylate that would additionally make its two oxygens
  inequivalent, which is chemically wrong and would poison every syn/anti question
  downstream.
* **Bridging (mu2/mu3) is DERIVED, not an edge type.**  A donor atom carrying
  dative edges to two distinct metals *is* mu2; storing that fact again as an edge
  label would give two encodings of one chemistry, and the two would hash
  differently.  See `TypedGraph.bridge_class`.  (Deviation from design doc §4.1's
  literal wording — proposed as D14.)

* **Net charge is a graph-level field, not a sum over atoms.**  Writing the -1 of a
  carboxylate onto one of its two oxygens makes those oxygens inequivalent, so a
  paddlewheel would hash differently depending on which way round four chemically
  identical bridges happened to be written.  Delocalised charge has no single home,
  so it is not given one.  `NodeLabel.formal_charge` remains, for charge that really
  is localised (an ammonium N), as a LABEL that enters identity but not the total.
  Protomers stay distinct because hydrogens are explicit nodes.
  (Proposed as D15.)

Dative direction is likewise derived, not stored: a DATIVE edge always runs
non-metal donor -> metal acceptor, and that is validated on insert.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from enum import Enum
from typing import Iterable, Iterator

import networkx as nx

from mofsbu._types import AmbiguousSpecError, GraphValidationError

# Elements treated as coordination centres.  Deliberately explicit rather than a
# periodic-table rule: metalloid edge cases should be a decision, not an accident.
METALS: frozenset[str] = frozenset("""
Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh Pd
Ag Cd In Sn Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir
Pt Au Hg Tl Pb Bi
""".split())


class EdgeType(str, Enum):
    """Bond types that enter identity."""

    COVALENT = "cov"
    DATIVE = "dat"          # non-metal donor -> metal acceptor
    METAL_METAL = "mm"


class BridgeClass(str, Enum):
    """Derived, never stored: how many distinct metals one donor atom serves."""

    NONE = "none"           # not coordinated
    TERMINAL = "terminal"   # one metal
    MU2 = "mu2"
    MU3 = "mu3"
    MU_N = "muN"            # four or more


_BRIDGE_RANK = {BridgeClass.NONE: 0, BridgeClass.TERMINAL: 1, BridgeClass.MU2: 2,
                BridgeClass.MU3: 3, BridgeClass.MU_N: 4}


def _bridge_class(n: int) -> BridgeClass:
    return {0: BridgeClass.NONE, 1: BridgeClass.TERMINAL,
            2: BridgeClass.MU2, 3: BridgeClass.MU3}.get(n, BridgeClass.MU_N)


@dataclass(frozen=True, slots=True)
class NodeLabel:
    """An atom's identity-bearing label.

    `oxidation_state` and `spin_class` are PER-CENTRE (D12) and belong on metals.
    They are part of the hashed label, so a mixed-valence node is a different
    structure from its homovalent analogue.
    """

    element: str
    formal_charge: int = 0
    oxidation_state: int | None = None
    spin_class: str | None = None       # "hs" | "ls" | "is" | None

    @property
    def is_metal(self) -> bool:
        return self.element in METALS

    def key(self) -> str:
        """The string fed to the hasher.  Stable across runs and platforms."""
        if self.is_metal:
            ox = "?" if self.oxidation_state is None else f"{self.oxidation_state:+d}"
            sp = self.spin_class or "?"
            return f"{self.element}/{self.formal_charge:+d}/{ox}/{sp}"
        return f"{self.element}/{self.formal_charge:+d}"

    def __str__(self) -> str:      # pragma: no cover - convenience
        return self.key()


class TypedGraph:
    """An explicit typed molecular graph.

    `charge` and `multiplicity` are both required before an L0 key can be computed,
    and both default to None so that construction never guesses (ground rule 5).
    """

    __slots__ = ("_g", "_next", "charge", "multiplicity", "name")

    def __init__(
        self,
        *,
        charge: int | None = None,
        multiplicity: int | None = None,
        name: str = "",
    ) -> None:
        self._g = nx.Graph()
        self._next = 0
        self.charge = charge
        self.multiplicity = multiplicity
        self.name = name

    # -- construction ----------------------------------------------------

    def add_atom(
        self,
        element: str,
        *,
        formal_charge: int = 0,
        oxidation_state: int | None = None,
        spin_class: str | None = None,
    ) -> int:
        label = NodeLabel(element, formal_charge, oxidation_state, spin_class)
        if not label.is_metal and (oxidation_state is not None or spin_class is not None):
            raise GraphValidationError(
                f"per-centre oxidation state / spin class is for metals only, got {element!r}"
            )
        idx = self._next
        self._next += 1
        self._g.add_node(idx, label=label)
        return idx

    def add_bond(self, i: int, j: int, etype: EdgeType, *, order: float = 1.0) -> None:
        if i == j:
            raise GraphValidationError(f"self-bond on atom {i}")
        for n in (i, j):
            if n not in self._g:
                raise GraphValidationError(f"unknown atom index {n}")
        mi, mj = self.is_metal(i), self.is_metal(j)
        if etype is EdgeType.DATIVE:
            if mi == mj:
                raise GraphValidationError(
                    f"DATIVE must run non-metal donor -> metal acceptor; got "
                    f"{self.label(i).element}-{self.label(j).element}"
                )
        elif etype is EdgeType.METAL_METAL:
            if not (mi and mj):
                raise GraphValidationError("METAL_METAL must join two metals")
        elif etype is EdgeType.COVALENT:
            if mi and mj:
                raise GraphValidationError(
                    "two metals joined COVALENT; use METAL_METAL so identity sees it"
                )
        self._g.add_edge(i, j, etype=etype, order=float(order))

    # -- accessors -------------------------------------------------------

    def __len__(self) -> int:
        return self._g.number_of_nodes()

    def nodes(self) -> list[int]:
        return sorted(self._g.nodes)

    def edges(self) -> list[tuple[int, int, EdgeType]]:
        return sorted((min(u, v), max(u, v), d["etype"]) for u, v, d in self._g.edges(data=True))

    def label(self, i: int) -> NodeLabel:
        return self._g.nodes[i]["label"]

    def is_metal(self, i: int) -> bool:
        return self.label(i).is_metal

    def metals(self) -> list[int]:
        return [i for i in self.nodes() if self.is_metal(i)]

    def neighbors(self, i: int, etype: EdgeType | None = None) -> list[int]:
        out = []
        for j in self._g.neighbors(i):
            if etype is None or self._g.edges[i, j]["etype"] is etype:
                out.append(j)
        return sorted(out)

    def edge_type(self, i: int, j: int) -> EdgeType:
        return self._g.edges[i, j]["etype"]

    def adjacency(self) -> dict[int, list[tuple[str, int]]]:
        """{atom: [(edge_type_value, neighbour), ...]} — the hashing view."""
        return {
            i: sorted((self._g.edges[i, j]["etype"].value, j) for j in self._g.neighbors(i))
            for i in self.nodes()
        }

    # -- derived chemistry -----------------------------------------------

    def bridging_metals(self, i: int) -> tuple[int, ...]:
        """Distinct metals this atom donates to.  Derived — never stored."""
        if self.is_metal(i):
            return ()
        return tuple(sorted(j for j in self.neighbors(i, EdgeType.DATIVE)))

    def bridge_class(self, i: int) -> BridgeClass:
        """How many metals this ATOM donates to.  A single-atom mu3-oxo is MU3 here."""
        return _bridge_class(len(self.bridging_metals(i)))

    def ligand_fragments(self) -> list[tuple[int, ...]]:
        """Connected components of the non-metal atoms, over covalent bonds only.

        Metals are the things fragments are joined *through*, so they are removed
        before the components are found — otherwise every ligand on a node would come
        back as one fragment.
        """
        seen: set[int] = set()
        out: list[tuple[int, ...]] = []
        for start in self.nodes():
            if start in seen or self.is_metal(start):
                continue
            stack, comp = [start], []
            seen.add(start)
            while stack:
                v = stack.pop()
                comp.append(v)
                for u in self.neighbors(v, EdgeType.COVALENT):
                    if u not in seen and not self.is_metal(u):
                        seen.add(u)
                        stack.append(u)
            out.append(tuple(sorted(comp)))
        return sorted(out)

    def fragment_bridge_class(self, fragment: Iterable[int]) -> BridgeClass:
        """How many distinct metals a LIGAND spans — the chemically meaningful one.

        A paddlewheel's carboxylate is mu2 even though each of its oxygens is terminal:
        the bridge is a property of the ligand, not of one donor atom.  This is what
        binding-mode filtering must key on.
        """
        metals: set[int] = set()
        for i in fragment:
            metals.update(self.bridging_metals(i))
        return _bridge_class(len(metals))

    def max_bridge_class(self) -> BridgeClass:
        classes = [self.fragment_bridge_class(f) for f in self.ligand_fragments()]
        return max(classes, key=_BRIDGE_RANK.__getitem__, default=BridgeClass.NONE)

    def net_charge(self) -> int:
        """The graph-level total.  Not a sum over atoms — see the module docstring."""
        if self.charge is None:
            raise AmbiguousSpecError(
                f"{self!r} has no net charge; set it explicitly "
                "(TypedGraph(charge=...) or .with_charge(...))"
            )
        return self.charge

    def localised_charge(self) -> int:
        """Sum of the per-atom formal-charge labels.  Diagnostic, not the total."""
        return sum(self.label(i).formal_charge for i in self.nodes())

    def element_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for i in self.nodes():
            el = self.label(i).element
            counts[el] = counts.get(el, 0) + 1
        return counts

    # -- transformation ---------------------------------------------------

    def relabel(self, mapping: dict[int, int]) -> TypedGraph:
        """Return a copy with atom indices remapped.  Identity must be invariant."""
        if sorted(mapping) != self.nodes() or sorted(mapping.values()) != self.nodes():
            raise GraphValidationError("relabel mapping must be a permutation of the atom indices")
        out = TypedGraph(charge=self.charge, multiplicity=self.multiplicity, name=self.name)
        out._next = self._next
        for i in self.nodes():
            out._g.add_node(mapping[i], label=self.label(i))
        for u, v, et in self.edges():
            out._g.add_edge(mapping[u], mapping[v], etype=et, order=self._g.edges[u, v]["order"])
        return out

    def permuted(self, seed: int) -> TypedGraph:
        """A deterministically shuffled relabelling — the invariance-test workhorse."""
        import random

        nodes = self.nodes()
        shuffled = list(nodes)
        random.Random(seed).shuffle(shuffled)
        return self.relabel(dict(zip(nodes, shuffled)))

    def subgraph(self, idxs: Iterable[int]) -> TypedGraph:
        keep = sorted(set(idxs))
        out = TypedGraph(charge=None, multiplicity=None, name=f"{self.name}:sub")
        remap = {}
        for i in keep:
            remap[i] = out.add_atom(
                self.label(i).element,
                formal_charge=self.label(i).formal_charge,
                oxidation_state=self.label(i).oxidation_state,
                spin_class=self.label(i).spin_class,
            )
        for u, v, et in self.edges():
            if u in remap and v in remap:
                out.add_bond(remap[u], remap[v], et, order=self._g.edges[u, v]["order"])
        return out

    def with_multiplicity(self, multiplicity: int) -> TypedGraph:
        out = self.relabel({i: i for i in self.nodes()})
        out.multiplicity = multiplicity
        return out

    def with_charge(self, charge: int) -> TypedGraph:
        out = self.relabel({i: i for i in self.nodes()})
        out.charge = charge
        return out

    # -- interchange -------------------------------------------------------

    def to_networkx(self) -> nx.Graph:
        """A plain nx.Graph carrying only the hashed labels."""
        g = nx.Graph()
        for i in self.nodes():
            g.add_node(i, lbl=self.label(i).key())
        for u, v, et in self.edges():
            g.add_edge(u, v, elbl=et.value)      # NB: bond order deliberately absent
        return g

    def to_dict(self) -> dict:
        return {
            "schema": "typed_graph/1",
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "name": self.name,
            "atoms": [
                {
                    "i": i,
                    "el": self.label(i).element,
                    "q": self.label(i).formal_charge,
                    "ox": self.label(i).oxidation_state,
                    "spin": self.label(i).spin_class,
                }
                for i in self.nodes()
            ],
            "bonds": [
                {"i": u, "j": v, "t": et.value, "order": self._g.edges[u, v]["order"]}
                for u, v, et in self.edges()
            ],
        }

    def to_json(self) -> bytes:
        return json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":")).encode()

    @classmethod
    def from_dict(cls, d: dict) -> TypedGraph:
        g = cls(charge=d.get("charge"), multiplicity=d.get("multiplicity"),
                name=d.get("name", ""))
        remap: dict[int, int] = {}
        for a in d["atoms"]:
            remap[a["i"]] = g.add_atom(
                a["el"], formal_charge=a.get("q", 0),
                oxidation_state=a.get("ox"), spin_class=a.get("spin"),
            )
        for b in d["bonds"]:
            g.add_bond(remap[b["i"]], remap[b["j"]], EdgeType(b["t"]), order=b.get("order", 1.0))
        return g

    @classmethod
    def from_json(cls, blob: bytes) -> TypedGraph:
        return cls.from_dict(json.loads(blob))

    def __repr__(self) -> str:      # pragma: no cover - convenience
        return (f"TypedGraph({self.name!r}, {len(self)} atoms, "
                f"{self._g.number_of_edges()} bonds, q={self.charge}, "
                f"mult={self.multiplicity})")
