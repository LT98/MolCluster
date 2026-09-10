"""Reaction-balanced energies — the fix for the charged-ion reference problem (M7).

**What was wrong, precisely.**  `legacy/xtb_energy_metal.py` scored every candidate as

    E_form = E(EBU, q) - E(M^q+) - SUM E(free ligand anion)

and the report calls that equation "charge-conserving", which it is.  **That is the trap.**
It balances in atoms and in charge and it is still not a number you can subtract, because
balance is a NECESSARY condition and not a sufficient one.  What breaks it is that the
reference species are nothing like the product: a bare Ni(2+) with no ligand field at all,
and a BTC(3-) whose three negative charges have nothing to sit on, are exactly the two
species gas-phase GFN2-xTB describes worst, and the product is a saturated complex where
neither error appears.  Nothing cancels.  The equation destroys six dative bonds and
creates none.

`docs/reports/solvation_report.md` is the measurement of the damage: adding a medium moved
the Ni/EDTA sequestration margin from -8.5 eV to -0.3 eV, and the qualitative conclusion
("EDTA poisons Ni-BTC crystallisation") reversed.  A balance check alone would have passed
that equation without a murmur, which is why this module checks two things and not one.

**What replaces it.**  Nothing here invents a better absolute energy — that is not
available at this level of theory.  What it does is refuse to subtract two numbers that
are not comparable, and make the comparison the caller is actually entitled to:

* the equation must **balance in atoms and in charge** — necessary, and checked first;
* it must be **isodesmic enough to be worth subtracting**: no bare metal ion standing in
  for a coordinated one, no naked polyanion, and the same number of metal-donor bonds on
  both sides.  A ligand-EXCHANGE equation satisfies this; a formation-from-free-ions
  equation does not, and that difference is the whole of M7;
* every species must come from the **same level of theory** — same code, version, method
  and medium — since the difference of two energies from two theories is a difference of
  theories;
* a **charge-blind backend** is refused on any equation that involves charge at all.
  "MACE" is not the test — MACE-MP-0 sees only elements and positions and is refused,
  while MACE-OMOL-0 is given total charge and spin multiplicity and is not.  The flag
  lives on the backend and is written into the `methods` row (`charge_blind`), so the
  refusal reads the stored number's own provenance rather than its code name;
* the non-physical test backend is refused unless the caller opts in by name.

Each refusal raises with the specific imbalance named.  None of them returns the number
anyway with a warning attached, because a warning is not visible in a CSV six months
later and the number looks exactly like a good one.

**The shape of a balanced equation.**  A ligand-exchange form gets there without ever
isolating a highly charged species:

    [M(H2O)6]^q+  +  n H(k)L  ->  [M L(n) (H2O)(6-m)]^(q-nk)  +  m H2O  +  nk H3O^+

Both sides carry the same atoms and the same total charge; the proton released by the
ligand is carried by a solvent molecule rather than emitted as a bare H+.  `check_balance`
is what verifies that, and it does not care how the equation was arrived at — a caller
that has a better one gets the same guarantee.
"""
from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

from mofsbu._types import Fidelity, MethodSpec, ReferenceSchemeError
from mofsbu.graph import EdgeType
from mofsbu.registry.api import RegistryError, get_graph
from mofsbu.registry.db import Registry, utcnow
from mofsbu.versions import ALGO_VERSIONS

# `leaving` is the one role that sits on the RIGHT of the arrow.  Displaced water, an
# expelled counter-ion, the proton acceptor that carries away an H+ — without a
# product-side role the equation can never balance and the checker would reject every
# real exchange reaction.
PRODUCT_SIDE_ROLES = frozenset({"leaving"})
REAGENT_SIDE_ROLES = frozenset({"reagent", "chelator", "solvent"})


# ── the equation ─────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Term:
    """One species in the equation, with the sign the arrow gives it."""

    structure_id: int
    stoich: int                 # always positive
    role: str
    side: str                   # 'reagent' | 'product'
    label: str = ""

    @property
    def sign(self) -> int:
        return +1 if self.side == "product" else -1


@dataclass(frozen=True)
class BalanceReport:
    """Why an equation is or is not usable.  Returned, not raised, so a UI can show it."""

    balanced: bool
    element_delta: dict[str, int]           # products - reagents, non-zero entries only
    charge_delta: int
    terms: tuple[Term, ...] = ()

    def describe(self) -> str:
        if self.balanced:
            return "balanced in atoms and charge"
        bits = []
        if self.element_delta:
            over = ", ".join(f"{el}{d:+d}" for el, d in sorted(self.element_delta.items()))
            bits.append(f"atoms {over} (products minus reagents)")
        if self.charge_delta:
            bits.append(f"charge {self.charge_delta:+d}")
        return "unbalanced: " + "; ".join(bits)


@dataclass(frozen=True)
class ReactionEnergy:
    """A reaction energy that knows what it is allowed to claim."""

    dE: float                               # eV, products minus reagents
    method: MethodSpec
    fidelity: Fidelity
    balance: BalanceReport
    quality: ReferenceQuality
    all_converged: bool
    contributions: tuple[tuple[int, int, float], ...] = ()   # (structure_id, signed n, E)
    scheme: str = field(default_factory=lambda: ALGO_VERSIONS["reference_scheme"])

    def describe(self) -> str:
        flags = []
        if not self.all_converged:
            flags.append("unconverged geometries")
        if not self.quality.isodesmic:
            flags.append("NOT isodesmic — relative comparisons only")
        note = f"  [{'; '.join(flags)}]" if flags else ""
        return f"dE = {self.dE:+.3f} eV  ({self.method.describe()}, {self.fidelity.name}){note}"


# ── reading the equation out of the registry ─────────────────────────────────


def reaction_terms(reg: Registry, reaction_id: int) -> tuple[Term, ...]:
    """Every species in the reaction, each on the side its role puts it."""
    row = reg.conn.execute(
        "SELECT product_structure_id FROM reactions WHERE id=?", (reaction_id,)).fetchone()
    if row is None:
        raise RegistryError(f"no reaction {reaction_id}")

    terms = [Term(structure_id=int(row["product_structure_id"]), stoich=1,
                  role="product", side="product")]
    for r in reg.conn.execute(
            "SELECT structure_id, stoich, role FROM reaction_reagents WHERE reaction_id=? "
            "ORDER BY role, structure_id", (reaction_id,)):
        role = str(r["role"])
        if role in PRODUCT_SIDE_ROLES:
            side = "product"
        elif role in REAGENT_SIDE_ROLES:
            side = "reagent"
        else:
            raise ReferenceSchemeError(
                f"reaction {reaction_id}: role {role!r} is on neither side of the arrow. "
                f"Reagent-side roles are {sorted(REAGENT_SIDE_ROLES)}, product-side "
                f"{sorted(PRODUCT_SIDE_ROLES)}.")
        terms.append(Term(structure_id=int(r["structure_id"]), stoich=int(r["stoich"]),
                          role=role, side=side))
    return tuple(terms)


def check_balance(reg: Registry, reaction_id: int) -> BalanceReport:
    """Does this equation conserve atoms and charge?

    Uses the typed graphs, not the formulas or the labels: the graph is the only thing
    identity is allowed to come from (ground rule 2), and a formula string parsed back
    into counts is a second encoding of the same fact waiting to disagree.
    """
    terms = reaction_terms(reg, reaction_id)
    elements: Counter[str] = Counter()
    charge = 0
    for term in terms:
        g = get_graph(reg, term.structure_id)
        if g.charge is None:
            raise ReferenceSchemeError(
                f"structure {term.structure_id} has no recorded charge; a charge balance "
                f"cannot be checked and so no reaction energy can be computed")
        weight = term.sign * term.stoich
        for element, n in g.element_counts().items():
            elements[element] += weight * n
        charge += weight * int(g.charge)

    delta = {el: n for el, n in elements.items() if n}
    return BalanceReport(balanced=not delta and charge == 0,
                         element_delta=delta, charge_delta=charge, terms=terms)


# ── is this equation worth subtracting? ──────────────────────────────────────
#
# Balance is necessary and not sufficient.  The legacy formation equation balances and is
# still worthless, so the second check asks whether the two sides describe comparable
# chemistry.  Three rules, each one a species the archived runs actually used as a
# reference and each one measurable off the typed graph.

BARE_ION = "bare_ion"
NAKED_POLYANION = "naked_polyanion"
COORDINATION_CHANGE = "coordination_change"


@dataclass(frozen=True)
class QualityIssue:
    code: str
    detail: str
    hint: str
    structure_id: int | None = None


@dataclass(frozen=True)
class ReferenceQuality:
    """Whether the equation's two sides are comparable enough to subtract."""

    isodesmic: bool
    dative_delta: int                       # metal-donor bonds, products minus reagents
    issues: tuple[QualityIssue, ...] = ()

    def describe(self) -> str:
        if self.isodesmic:
            return "isodesmic: same metal-donor bond count, no isolated reference species"
        return "; ".join(f"{i.code}: {i.detail}" for i in self.issues)


def _dative_count(g: Any) -> int:
    return sum(1 for _, _, etype in g.edges() if etype is EdgeType.DATIVE)


def check_reference_quality(reg: Registry, reaction_id: int, *,
                            max_naked_charge: int = 1) -> ReferenceQuality:
    """Are the reference species comparable to the product, or merely balanced?

    * **bare ion** — a lone metal atom carrying a charge, the `E(M^q+)` term.  Gas-phase
      xTB on a naked cation has no ligand field to get right or wrong, so its error has
      nothing in common with the complex it is subtracted from.
    * **naked polyanion** — a metal-free species with |charge| > `max_naked_charge`, the
      `E(BTC^3-)` term.  Three negative charges with nowhere to go is the single worst
      case for a semi-empirical gas-phase method.
    * **coordination change** — the count of metal-donor (dative) bonds must be the same
      on both sides.  This is the rule that separates a ligand-exchange equation, where
      the metal stays coordinated throughout and the errors cancel, from a
      formation-from-free-ions equation, where six bonds appear out of nothing.
    """
    terms = reaction_terms(reg, reaction_id)
    issues: list[QualityIssue] = []
    dative = 0
    for term in terms:
        g = get_graph(reg, term.structure_id)
        dative += term.sign * term.stoich * _dative_count(g)
        charge = int(g.charge or 0)
        metals = g.metals()
        if len(g) == 1 and metals and charge != 0:
            issues.append(QualityIssue(
                BARE_ION, structure_id=term.structure_id,
                detail=f"{g.label(metals[0]).element}{charge:+d} appears as a bare ion",
                hint="reference it as a solvated complex — [M(H2O)n]^q+ — so the metal "
                     "has a ligand field on both sides of the arrow"))
        if not metals and abs(charge) > max_naked_charge:
            issues.append(QualityIssue(
                NAKED_POLYANION, structure_id=term.structure_id,
                detail=f"{g.name or 'species'} carries charge {charge:+d} with no metal "
                       f"to sit on",
                hint="use the protonated acid plus an explicit proton acceptor "
                     "(H3O+, protonated solvent) instead of the free polyanion"))
    if dative:
        issues.append(QualityIssue(
            COORDINATION_CHANGE,
            detail=f"metal-donor bonds change by {dative:+d} across the arrow",
            hint="write it as a LIGAND EXCHANGE — the incoming donor replaces a solvent "
                 "molecule rather than binding a metal that started with nothing"))
    return ReferenceQuality(isodesmic=not issues, dative_delta=dative,
                            issues=tuple(issues))


# ── the energies behind the equation ─────────────────────────────────────────


def _method_of(reg: Registry, method_id: int | None) -> MethodSpec:
    if method_id is None:
        raise ReferenceSchemeError("a stored energy has no method row (ground rule 3)")
    row = reg.conn.execute("SELECT * FROM methods WHERE id=?", (method_id,)).fetchone()
    if row is None:
        raise ReferenceSchemeError(f"no method {method_id}")
    return MethodSpec(code=row["code"], code_version=row["code_version"],
                      method=row["method"], solvent=row["solvent"],
                      charge=row["charge"], multiplicity=row["multiplicity"],
                      extras=json.loads(row["extras_json"] or "{}"))


def _energy_row(reg: Registry, structure_id: int, *, fidelity: Fidelity | None,
                solvent: str | None, method: str | None = None) -> Any:
    """The energy this species brings to the equation.

    Highest fidelity with an energy, or exactly the requested rung.  A species with no
    energy at all is named rather than skipped: dropping a term from a balanced equation
    is how you get a number that is confidently wrong.

    `method` pins the THEORY, and exists because a rung is not one.  Since MACE-OMOL-0
    joined MACE-MP-0 on the ML rung, "the best ML energy for this species" can mean two
    numbers on two scales, and the old ordering picked whichever was numerically lower —
    i.e. whichever model had the deeper reference, for every species, consistently
    wrong.  `reaction_balanced_energy` pins it from the first term it resolves, so the
    whole equation is answered in one theory or refuses; the tie-break below only
    decides which theory that first term offers, and it prefers a charge-aware method,
    never a smaller float.
    """
    sql = ("SELECT g.id, g.energy, g.converged, g.fidelity, g.method_id "
           "FROM geometries g JOIN methods m ON m.id = g.method_id "
           "WHERE g.structure_id=? AND g.energy IS NOT NULL AND m.solvent IS ?")
    args: list[Any] = [structure_id, solvent]
    if fidelity is not None:
        sql += " AND g.fidelity=?"
        args.append(int(fidelity))
    if method is not None:
        sql += " AND m.method=?"
        args.append(method)
    sql += (" ORDER BY g.fidelity DESC, (g.converged IS 1) DESC, "
            "         COALESCE(json_extract(m.extras_json, '$.charge_blind'), 0) ASC, "
            "         m.id ASC, g.energy ASC LIMIT 1")
    return reg.conn.execute(sql, args).fetchone()


def reaction_balanced_energy(
    reg: Registry,
    reaction: int,
    *,
    fidelity: Fidelity | None = None,
    solvent: str | None = None,
    allow_null: bool = False,
    allow_unconverged: bool = True,
    strict: bool = True,
) -> ReactionEnergy:
    """dE for a balanced reaction, or a refusal that says which rule was broken.

    `solvent` selects the medium the whole equation is evaluated in (None = gas phase);
    every species must have an energy in that same medium, which is what stops a solvated
    product being compared against a gas-phase reagent.

    `strict` (default on) additionally requires the equation to be isodesmic — see
    `check_reference_quality`.  Turning it off is how you reproduce a legacy number on
    purpose: the value comes back with the quality report attached and `isodesmic=False`,
    so the caveat travels with the number instead of living in a report nobody re-reads.
    """
    balance = check_balance(reg, reaction)
    if not balance.balanced:
        raise ReferenceSchemeError(
            f"reaction {reaction} is {balance.describe()}. A formation energy from an "
            f"unbalanced equation is the charged-ion reference problem by another name — "
            f"add the missing species (displaced solvent, the proton acceptor, a "
            f"counter-ion) as `leaving` reagents so both sides carry the same atoms and "
            f"the same charge.")

    quality = check_reference_quality(reg, reaction)
    if strict and not quality.isodesmic:
        lines = "\n".join(f"  - {i.code}: {i.detail}\n    {i.hint}" for i in quality.issues)
        raise ReferenceSchemeError(
            f"reaction {reaction} balances but is not isodesmic, so subtracting its two "
            f"sides compares different chemistry, not different arrangements of the same "
            f"chemistry:\n{lines}\n"
            f"This is the archived E(EBU) - E(M^q+) - SUM E(ligand anion) scheme, whose "
            f"Ni/EDTA verdict reversed once a medium was added "
            f"(docs/reports/solvation_report.md). Pass strict=False to reproduce it "
            f"deliberately; the result will carry isodesmic=False.")

    method: MethodSpec | None = None
    total = 0.0
    all_converged = True
    rungs: list[int] = []
    contributions: list[tuple[int, int, float]] = []

    for term in balance.terms:
        # `method` is None for the first term and pinned for every one after it, so the
        # equation cannot be assembled out of two theories that happen to share a rung.
        # The `same_theory` check below still stands: it catches a mismatch in code
        # version or medium that the method NAME alone would let through.
        row = _energy_row(reg, term.structure_id, fidelity=fidelity, solvent=solvent,
                          method=None if method is None else method.method)
        if row is None:
            medium = solvent or "gas phase"
            rung = f" at {fidelity.name}" if fidelity is not None else ""
            theory = f" from {method.method}" if method is not None else ""
            raise ReferenceSchemeError(
                f"structure {term.structure_id} ({term.role}) has no energy{rung}{theory} "
                f"in {medium}; every species in the equation needs one, because dropping "
                f"a term silently changes what the number means"
                + (f". It may have one from another method — a rung is not a theory, and "
                   f"the equation was pinned to {method.method} by its first term"
                   if method is not None else ""))
        spec = _method_of(reg, row["method_id"])

        if spec.code == "null" and not allow_null:
            raise ReferenceSchemeError(
                f"structure {term.structure_id} was scored by the null backend, which is "
                f"deliberately non-physical; pass allow_null=True if this is a test")
        if spec.extras.get("charge_blind"):
            charges = {int(get_graph(reg, t.structure_id).charge or 0) for t in balance.terms}
            if charges != {0}:
                raise ReferenceSchemeError(
                    f"structure {term.structure_id} was scored by a charge-blind backend "
                    f"({spec.method}) but the equation involves charged species "
                    f"({sorted(charges)}). That backend cannot see the difference between "
                    f"them; use MACE-OMOL-0 (MOFSBU_ML_MODEL=mace-omol-0, or "
                    f"ml_model='mace-omol-0' on the spec) at ML cost, or xtb, for "
                    f"anything the reference scheme has to weigh.")
        if method is None:
            method = spec
        elif not method.same_theory(spec):
            raise ReferenceSchemeError(
                f"mixed levels of theory in one equation: {method.describe()} vs "
                f"{spec.describe()}. The difference of two energies from two theories is "
                f"a difference of theories.")

        rungs.append(int(row["fidelity"]))
        converged = row["converged"]
        if converged is not None and not converged:
            all_converged = False
        signed = term.sign * term.stoich
        total += signed * float(row["energy"])
        contributions.append((term.structure_id, signed, float(row["energy"])))

    if method is None:                      # unreachable: a reaction always has a product
        raise ReferenceSchemeError(f"reaction {reaction} has no species")
    if not all_converged and not allow_unconverged:
        raise ReferenceSchemeError(
            f"reaction {reaction} draws on an unconverged geometry; its energy is an "
            f"upper bound, not a minimum")

    # The equation is only as good as its weakest term: a dE mixing an xTB product with
    # a raw-construct reagent is a raw-construct number, and labelling it xTB is how a
    # screening result gets quoted as if it were a relaxation.
    rung = Fidelity(min(rungs))
    return ReactionEnergy(dE=total, method=method, fidelity=rung, balance=balance,
                          quality=quality, all_converged=all_converged,
                          contributions=tuple(contributions))


# ── writing it back ──────────────────────────────────────────────────────────


def store_reaction_energy(reg: Registry, reaction_id: int, energy: ReactionEnergy, *,
                          force: bool = False) -> None:
    """Record dE on the reaction, with the method and fidelity that produced it.

    The one writer of `reactions.dG` / `method_id` / `fidelity`.  The column is named
    dG for the pathway layer that will eventually put a free energy in it; what M7 puts
    there is an electronic dE, and the method row says so rather than the name implying
    a thermochemistry that has not been done.
    """
    from mofsbu.registry.api import method_id as ensure_method

    if not energy.quality.isodesmic and not force:
        raise ReferenceSchemeError(
            f"refusing to store a non-isodesmic dE on reaction {reaction_id} "
            f"({energy.quality.describe()}); pass force=True and the caveat is recorded "
            f"in the note alongside it")
    mid = ensure_method(reg, energy.method)
    caveats = "" if energy.quality.isodesmic else (
        " NOT-ISODESMIC:" + ",".join(i.code for i in energy.quality.issues))
    if not energy.all_converged:
        caveats += " UNCONVERGED"
    reg.conn.execute(
        "UPDATE reactions SET dG=?, method_id=?, fidelity=?, note=? WHERE id=?",
        (energy.dE, mid, int(energy.fidelity),
         f"balanced-dE {ALGO_VERSIONS['reference_scheme']} @ {utcnow()}{caveats}",
         reaction_id))


def put_balanced_reaction(
    reg: Registry,
    product_structure_id: int,
    *,
    reagents: Iterable[tuple[int, int]] = (),
    leaving: Iterable[tuple[int, int]] = (),
    chelators: Iterable[tuple[int, int]] = (),
    solvents: Iterable[tuple[int, int]] = (),
    kind: str = "reaction",
    note: str = "",
    depth: int | None = None,
) -> int:
    """Store a reaction and REFUSE it if it does not balance.

    `put_reaction` stays as it is — an assembly provenance edge is a record of what was
    built from what, and it neither claims nor needs a balanced equation.  This is the
    stricter door, for the reactions that will carry energies: an equation that cannot
    hold a number is rejected at write time rather than at read time, when whoever hits
    it no longer remembers what it was meant to say.
    """
    cur = reg.conn.execute(
        "INSERT INTO reactions (product_structure_id, kind, note, depth, created_at) "
        "VALUES (?,?,?,?,?)", (product_structure_id, kind, note, depth, utcnow()))
    rid = int(cur.lastrowid)
    rows = ([(sid, n, "reagent") for sid, n in reagents]
            + [(sid, n, "chelator") for sid, n in chelators]
            + [(sid, n, "solvent") for sid, n in solvents]
            + [(sid, n, "leaving") for sid, n in leaving])
    for sid, n, role in rows:
        if n <= 0:
            raise ReferenceSchemeError(f"stoichiometry must be positive, got {n} for {sid}")
        reg.conn.execute(
            "INSERT INTO reaction_reagents (reaction_id, structure_id, stoich, role) "
            "VALUES (?,?,?,?)", (rid, sid, n, role))

    report = check_balance(reg, rid)
    if not report.balanced:
        reg.conn.execute("DELETE FROM reactions WHERE id=?", (rid,))
        raise ReferenceSchemeError(
            f"refusing to store this reaction: {report.describe()}. "
            f"An equation that cannot carry an energy should not be in the table that "
            f"energies are read from.")
    return rid
