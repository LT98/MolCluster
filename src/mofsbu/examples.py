"""Hand-written typed graphs — real chemistry, used as both test ground truth and demo data.

Identity needs no 3D.  These polynuclear nodes are written out by hand so that the
D12 commitments (per-centre oxidation state and spin, M-M edges, mu-bridges as one
ligand with two dative edges) are tested before any geometry placer exists.  When
the placer lands in M6, its output must hash to exactly these values.

These live in the package rather than under `tests/` so that scripts, notebooks and the
demo registry all draw their example structures from the same place.  A demo database
seeded from real graphs cannot drift into chemistry that does not exist, which is what
happened when the seed invented column values directly.
"""
from __future__ import annotations

from mofsbu.graph import EdgeType as E
from mofsbu.graph import TypedGraph

COV, DAT, MM = E.COVALENT, E.DATIVE, E.METAL_METAL


def _formate(g: TypedGraph) -> tuple[int, int, int, int]:
    """HCOO: one carbon, two equivalent oxygens, one hydrogen.

    The oxygens carry NO formal charge — the -1 is delocalised and lives on the
    graph, which is what keeps the two O's interchangeable (see D15).
    """
    c = g.add_atom("C")
    o1 = g.add_atom("O")
    o2 = g.add_atom("O")
    h = g.add_atom("H")
    g.add_bond(c, o1, COV)
    g.add_bond(c, o2, COV)
    g.add_bond(c, h, COV)
    return c, o1, o2, h


def _water(g: TypedGraph) -> tuple[int, int, int]:
    o = g.add_atom("O")
    h1 = g.add_atom("H")
    h2 = g.add_atom("H")
    g.add_bond(o, h1, COV)
    g.add_bond(o, h2, COV)
    return o, h1, h2


# -- free ligands --------------------------------------------------------------

def water() -> TypedGraph:
    g = TypedGraph(charge=0, multiplicity=1, name="water")
    _water(g)
    return g


def formate() -> TypedGraph:
    g = TypedGraph(charge=-1, multiplicity=1, name="formate")
    _formate(g)
    return g


def formic_acid() -> TypedGraph:
    """The protomer of `formate` — distinct at L1 because H is an explicit node."""
    g = TypedGraph(charge=0, multiplicity=1, name="formic_acid")
    _c, o1, _o2, _h = _formate(g)
    g.add_bond(o1, g.add_atom("H"), COV)
    return g


def btc(deprotonated: bool = True) -> TypedGraph:
    """Benzene-1,3,5-tricarboxylate (BTC3-) or its triacid."""
    g = TypedGraph(charge=-3 if deprotonated else 0, multiplicity=1,
                   name="btc" if deprotonated else "btc_h3")
    ring = [g.add_atom("C") for _ in range(6)]
    for k in range(6):
        g.add_bond(ring[k], ring[(k + 1) % 6], COV)
    for k in range(6):
        if k % 2 == 0:                       # 1,3,5 positions bear the carboxylates
            c = g.add_atom("C")
            o1, o2 = g.add_atom("O"), g.add_atom("O")
            g.add_bond(ring[k], c, COV)
            g.add_bond(c, o1, COV)
            g.add_bond(c, o2, COV)
            if not deprotonated:
                g.add_bond(o1, g.add_atom("H"), COV)
        else:
            g.add_bond(ring[k], g.add_atom("H"), COV)
    return g


# -- polynuclear nodes ---------------------------------------------------------

def cu_paddlewheel() -> TypedGraph:
    """Cu2(mu-O2CH)4 — the canonical paddlewheel.

    One node: two symmetry-related Cu joined by an M-M edge, four bridging
    carboxylates each contributing one ligand fragment with two dative edges to
    different metals, two open axial sites (absent here, hence unoccupied).
    """
    g = TypedGraph(charge=0, multiplicity=1, name="cu_paddlewheel")
    cu_a = g.add_atom("Cu", oxidation_state=2, spin_class="hs")
    cu_b = g.add_atom("Cu", oxidation_state=2, spin_class="hs")
    g.add_bond(cu_a, cu_b, MM)
    for _ in range(4):
        _c, o1, o2, _h = _formate(g)
        g.add_bond(o1, cu_a, DAT)
        g.add_bond(o2, cu_b, DAT)
    return g


def fe3_mu3_oxo(valences: tuple[int, int, int] = (3, 3, 3)) -> TypedGraph:
    """Fe3(mu3-O)(mu-O2CH)6 — the basic carboxylate trimer.

    `valences` sets the per-centre oxidation states, so the mixed-valence
    Fe(II)/Fe(III)/Fe(III) node and the all-Fe(III) node are different structures
    at L0 while sharing a connectivity (D12).
    """
    charge = sum(valences) - 2 - 6
    spins = {2: 5, 3: 6}                       # high-spin d6 / d5, per centre
    total_s = sum((spins[v] - 1) / 2 for v in valences)
    tag = "".join(str(v) for v in valences)
    g = TypedGraph(charge=charge, multiplicity=int(2 * total_s + 1),
                   name=f"fe3_mu3_oxo_{tag}")
    fes = [g.add_atom("Fe", oxidation_state=v, spin_class="hs") for v in valences]
    mu3 = g.add_atom("O")
    for fe in fes:
        g.add_bond(mu3, fe, DAT)
    for a, b in ((0, 1), (1, 2), (2, 0)):
        for _ in range(2):
            _c, o1, o2, _h = _formate(g)
            g.add_bond(o1, fes[a], DAT)
            g.add_bond(o2, fes[b], DAT)
    return g


# -- connectivity discrimination pair (same composition, different binding) -----

def zn2_bridged_formates() -> TypedGraph:
    """Zn2(mu-O2CH)2 — both carboxylates bridge the two zincs."""
    g = TypedGraph(charge=2, multiplicity=1, name="zn2_bridged")
    za = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    zb = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    for _ in range(2):
        _c, o1, o2, _h = _formate(g)
        g.add_bond(o1, za, DAT)
        g.add_bond(o2, zb, DAT)
    return g


def zn2_chelated_formates() -> TypedGraph:
    """Same atoms as `zn2_bridged_formates`, each formate chelating one zinc."""
    g = TypedGraph(charge=2, multiplicity=1, name="zn2_chelated")
    za = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    zb = g.add_atom("Zn", oxidation_state=2, spin_class="ls")
    for zn in (za, zb):
        _c, o1, o2, _h = _formate(g)
        g.add_bond(o1, zn, DAT)
        g.add_bond(o2, zn, DAT)
    return g


# -- L1 must NOT separate these (that is L2's job, M5) -------------------------

def pt_ammine_dichloride(order: str = "a") -> TypedGraph:
    """Pt(NH3)2Cl2 — cis and trans share this connectivity exactly.

    `order` only changes the sequence atoms are added in, so the two builds are the
    same graph written two ways.  L1 must return one hash for both.
    """
    g = TypedGraph(charge=0, multiplicity=1, name=f"pt_ammine_dichloride_{order}")
    pt = g.add_atom("Pt", oxidation_state=2, spin_class="ls")

    def ammine() -> None:
        n = g.add_atom("N")
        for _ in range(3):
            g.add_bond(n, g.add_atom("H"), COV)
        g.add_bond(n, pt, DAT)

    def chloride() -> None:
        g.add_bond(g.add_atom("Cl"), pt, DAT)

    steps = [ammine, ammine, chloride, chloride] if order == "a" else \
            [chloride, ammine, chloride, ammine]
    for step in steps:
        step()
    return g


# -- spin as identity ----------------------------------------------------------

def fe_hexaaqua(spin_class: str = "hs") -> TypedGraph:
    """[Fe(H2O)6]2+ high-spin vs low-spin — same connectivity, different species."""
    mult = 5 if spin_class == "hs" else 1
    g = TypedGraph(charge=2, multiplicity=mult, name=f"fe_hexaaqua_{spin_class}")
    fe = g.add_atom("Fe", oxidation_state=2, spin_class=spin_class)
    for _ in range(6):
        o, _h1, _h2 = _water(g)
        g.add_bond(o, fe, DAT)
    return g


# -- a genuine 1-WL collision --------------------------------------------------

def cyclohexane() -> TypedGraph:
    g = TypedGraph(charge=0, multiplicity=1, name="cyclohexane")
    cs = [g.add_atom("C") for _ in range(6)]
    for k in range(6):
        g.add_bond(cs[k], cs[(k + 1) % 6], COV)
        for _ in range(2):
            g.add_bond(cs[k], g.add_atom("H"), COV)
    return g


def two_cyclopropanes() -> TypedGraph:
    """Same formula and same local environments as cyclohexane, not isomorphic.

    1-WL cannot tell a 6-ring from two 3-rings, so this pair exercises the
    collision-resolution path that VF2 / the canonical certificate exist for.
    """
    g = TypedGraph(charge=0, multiplicity=1, name="two_cyclopropanes")
    for _ring in range(2):
        cs = [g.add_atom("C") for _ in range(3)]
        for k in range(3):
            g.add_bond(cs[k], cs[(k + 1) % 3], COV)
            for _ in range(2):
                g.add_bond(cs[k], g.add_atom("H"), COV)
    return g


# ── parameterised families, for populating a registry with real variety ──────

def paddlewheel(metal: str = "Cu", oxidation_state: int = 2, spin_class: str = "hs",
                charge: int | None = None, multiplicity: int = 1) -> TypedGraph:
    """M2(mu-O2CH)4 for any metal — the paddlewheel family."""
    q = 2 * oxidation_state - 4 if charge is None else charge
    g = TypedGraph(charge=q, multiplicity=multiplicity,
                   name=f"{metal}2(mu-O2CH)4 paddlewheel")
    a = g.add_atom(metal, oxidation_state=oxidation_state, spin_class=spin_class)
    b = g.add_atom(metal, oxidation_state=oxidation_state, spin_class=spin_class)
    g.add_bond(a, b, MM)
    for _ in range(4):
        _c, o1, o2, _h = _formate(g)
        g.add_bond(o1, a, DAT)
        g.add_bond(o2, b, DAT)
    return g


def hexaaqua(metal: str = "Fe", oxidation_state: int = 2, spin_class: str = "hs",
             multiplicity: int = 5) -> TypedGraph:
    g = TypedGraph(charge=oxidation_state, multiplicity=multiplicity,
                   name=f"[{metal}(H2O)6]{oxidation_state:+d}")
    m = g.add_atom(metal, oxidation_state=oxidation_state, spin_class=spin_class)
    for _ in range(6):
        o, _h1, _h2 = _water(g)
        g.add_bond(o, m, DAT)
    return g


def aqua_carboxylate(metal: str = "Zn", oxidation_state: int = 2, spin_class: str = "ls",
                     n_aqua: int = 4, multiplicity: int = 1) -> TypedGraph:
    """A mononuclear metal with one chelating formate and some aqua ligands."""
    g = TypedGraph(charge=oxidation_state - 1, multiplicity=multiplicity,
                   name=f"[{metal}(O2CH)(H2O){n_aqua}]")
    m = g.add_atom(metal, oxidation_state=oxidation_state, spin_class=spin_class)
    _c, o1, o2, _h = _formate(g)
    g.add_bond(o1, m, DAT)
    g.add_bond(o2, m, DAT)
    for _ in range(n_aqua):
        o, _h1, _h2 = _water(g)
        g.add_bond(o, m, DAT)
    return g


def bipyridine() -> TypedGraph:
    """2,2'-bipyridine — two pyridyl N donors, torsion-live."""
    g = TypedGraph(charge=0, multiplicity=1, name="bipy")
    joins = []
    for _ring in range(2):
        n = g.add_atom("N")
        cs = [g.add_atom("C") for _ in range(5)]
        ring = [n, *cs]
        for k in range(6):
            g.add_bond(ring[k], ring[(k + 1) % 6], COV)
        for c in cs[1:]:
            g.add_bond(c, g.add_atom("H"), COV)
        joins.append(cs[0])
    g.add_bond(joins[0], joins[1], COV)
    return g


def metal_bipy(metal: str = "Fe", oxidation_state: int = 2, spin_class: str = "ls",
               n_bipy: int = 3, multiplicity: int = 1) -> TypedGraph:
    """[M(bipy)n] — a tris-chelate, the fac/mer and Delta/Lambda case for L2 (M5)."""
    g = TypedGraph(charge=oxidation_state, multiplicity=multiplicity,
                   name=f"[{metal}(bipy){n_bipy}]{oxidation_state:+d}")
    m = g.add_atom(metal, oxidation_state=oxidation_state, spin_class=spin_class)
    for _ in range(n_bipy):
        joins, ns = [], []
        for _ring in range(2):
            n = g.add_atom("N")
            cs = [g.add_atom("C") for _ in range(5)]
            ring = [n, *cs]
            for k in range(6):
                g.add_bond(ring[k], ring[(k + 1) % 6], COV)
            for c in cs[1:]:
                g.add_bond(c, g.add_atom("H"), COV)
            joins.append(cs[0])
            ns.append(n)
        g.add_bond(joins[0], joins[1], COV)
        for n in ns:
            g.add_bond(n, m, DAT)
    return g


ALL: dict[str, "callable[[], TypedGraph]"] = {
    "water": water,
    "formate": formate,
    "formic_acid": formic_acid,
    "btc": lambda: btc(True),
    "btc_h3": lambda: btc(False),
    "cu_paddlewheel": cu_paddlewheel,
    "fe3_mu3_oxo_333": lambda: fe3_mu3_oxo((3, 3, 3)),
    "fe3_mu3_oxo_233": lambda: fe3_mu3_oxo((2, 3, 3)),
    "zn2_bridged": zn2_bridged_formates,
    "zn2_chelated": zn2_chelated_formates,
    "pt_ammine_dichloride_a": lambda: pt_ammine_dichloride("a"),
    "pt_ammine_dichloride_b": lambda: pt_ammine_dichloride("b"),
    "fe_hexaaqua_hs": lambda: fe_hexaaqua("hs"),
    "fe_hexaaqua_ls": lambda: fe_hexaaqua("ls"),
    "cyclohexane": cyclohexane,
    "two_cyclopropanes": two_cyclopropanes,
}
