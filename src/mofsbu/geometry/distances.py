"""How far a donor sits from the metal — a property of the PAIR, not of the metal.

`placer.D_ML` gave one distance per metal and applied it to every donor, so a Zn(II)
centre placed an iodide at the same 2.00 A as a water oxygen.  Zn-I is about 2.60 A.
Half an angstrom of error is not a cosmetic problem: the QC clash check compares against
van der Waals radii, iodine's is 1.98 A against oxygen's 1.52, and the two errors compound
in the same direction.  The result was that **every halide co-ligand was rejected by
construction** with "QC FAILED: 2 clash(es)" — a message that names the symptom and hides
the cause, and which reads identically to a genuine steric refusal.

The model here is deliberately small:

    d(M, D) = base(M) + delta(donor element)

`base` is the metal's M-O distance — the old `D_ML` table, unchanged, so every geometry
this package has ever built with an O or N donor comes out bit-identical.  `delta` is a
per-element offset relative to oxygen, calibrated against typical divalent first-row
transition-metal bond lengths rather than derived from covalent radii.  Radii were tried
first and are wrong in a consistent direction: a dative M-O bond is ~0.12 A longer than
r_cov(M) + r_cov(O), while M-I is about right, so a radii sum underestimates the light
donors and the correction is not a constant.

Everything the model does not know falls back to covalent radii and **says so**.  A
distance carries its `source`, the placer writes it into the choice vector, and the run
inspector shows it, because "this bond length is an estimate" is exactly the kind of fact
that has to travel with the number instead of being rediscovered from a clash report.
"""
from __future__ import annotations

from dataclasses import dataclass

from mofsbu._types import AmbiguousSpecError

#: M-O distance per metal, angstrom.  This IS the old `placer.D_ML` table, kept as the
#: base of the model so existing geometries and their tests do not move.
BASE_MO: dict[str, float] = {
    "Zn": 2.00, "Cu": 1.98, "Ni": 2.06, "Co": 2.08, "Fe": 2.10,
    "Mn": 2.18, "Cr": 2.00, "Mg": 2.10, "Ca": 2.40, "Cd": 2.28,
}
DEFAULT_BASE_MO = 2.05

#: Offset from the metal's M-O distance, angstrom, by donor element.
#:
#: Calibrated against typical divalent first-row transition-metal bond lengths, not
#: against radii.  Check the arithmetic on Zn (base 2.00): Zn-N 2.05, Zn-S 2.32,
#: Zn-Cl 2.25, Zn-Br 2.40, Zn-I 2.60 — each within a few hundredths of the CSD-typical
#: value.  On Cu (base 1.98): Cu-N 2.03, Cu-Cl 2.23, Cu-S 2.30.  Same on Fe(II).
#:
#: Known weak spot, recorded rather than hidden: phosphine on a low-spin square-planar
#: d8 centre is genuinely shorter than this (Ni-P ~2.2, not 2.5).  Spin state is not an
#: input to this model.  A P donor therefore reports `source="element-offset"` like the
#: rest, and if that case starts mattering the fix is a (metal, spin, donor) table, not
#: a fudge here.
DELTA_BY_ELEMENT: dict[str, float] = {
    "O": 0.00,      # the reference donor
    "N": 0.05,
    "F": -0.05,
    "C": 0.05,      # carbene / cyanide carbon
    "P": 0.45,
    "S": 0.32,
    "Cl": 0.25,
    "Se": 0.45,
    "Br": 0.40,
    "As": 0.50,
    "I": 0.60,
    "Te": 0.62,
}

#: Cordero (2008) covalent radii, angstrom.  Used ONLY by the fallback.
COVALENT_RADII: dict[str, float] = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Te": 1.38, "I": 1.39,
    "Li": 1.28, "Na": 1.66, "K": 2.03, "Mg": 1.41, "Ca": 1.76, "Al": 1.21,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.61, "Fe": 1.52,
    "Co": 1.50, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Zr": 1.75, "Ru": 1.46,
    "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
}
DEFAULT_COVALENT = 1.40

#: A dative M-L bond is longer than the covalent-radii sum.  Measured off the pairs the
#: calibrated model does cover, this is the mean residual — it is a correction to a
#: fallback, and it is why the fallback is labelled an estimate.
DATIVE_LENGTHENING = 0.10


@dataclass(frozen=True)
class Distance:
    """A metal-donor distance and where it came from.

    `source` is not decoration.  `"element-offset"` is a calibrated number;
    `"covalent-radii"` is an estimate for a pair nobody has checked; `"override"` means
    a caller supplied it.  The difference decides whether a QC failure is evidence about
    the chemistry or evidence about this table.
    """

    value: float
    source: str
    metal: str
    donor_element: str

    @property
    def estimated(self) -> bool:
        return self.source == "covalent-radii"

    def to_dict(self) -> dict:
        return {"d": round(self.value, 3), "source": self.source,
                "metal": self.metal, "donor": self.donor_element,
                "estimated": self.estimated}


def base_distance(metal: str) -> float:
    """The metal's M-O distance: the base of the model."""
    return BASE_MO.get(metal, DEFAULT_BASE_MO)


def metal_donor_distance(metal: str, donor_element: str,
                         *, override: float | None = None) -> Distance:
    """Ideal M-D distance for this pair.

    `donor_element` is the ELEMENT SYMBOL of the donor atom, read off the molecule —
    never parsed out of a donor-type name.  `hydrohalide_X` and `halide_I` name the same
    iodine; `carboxylate_O` and `enolate_O` name the same oxygen at different distances
    only because of what they are attached to, which this model does not claim to know.
    Taking the element from the atom keeps a new donor type from silently inheriting
    oxygen's bond length.
    """
    if override is not None:
        return Distance(float(override), "override", metal, donor_element)
    delta = DELTA_BY_ELEMENT.get(donor_element)
    if delta is not None:
        return Distance(base_distance(metal) + delta, "element-offset", metal, donor_element)
    d = (COVALENT_RADII.get(metal, DEFAULT_COVALENT)
         + COVALENT_RADII.get(donor_element, DEFAULT_COVALENT)
         + DATIVE_LENGTHENING)
    return Distance(d, "covalent-radii", metal, donor_element)


def donor_elements(mol, donor_idxs) -> tuple[str, ...]:
    """Element symbols of a ligand's donor atoms, from the molecule itself."""
    return tuple(mol.GetAtomWithIdx(int(i)).GetSymbol() for i in donor_idxs)


def metal_metal_distance(m1: str, m2: str, *, motif: str | None = None,
                         override: float | None = None) -> Distance:
    """M-M distance for a DECLARED nucleus.  The caller states it or gets an estimate.

    `Distance.donor_element` carries the second metal here; the dataclass is not widened
    for one caller, and no M-M distance reaches a stored choice vector as a dict.

    There is no curated M-M table and this function will not stand in for one.  A
    quadruply-bonded Mo2 at 2.09 A and a Cu2 paddlewheel at 2.62 A are the same question
    to every input this function has — the bond order is the answer and the element
    symbols do not carry it.  So `motif` exists to be REFUSED: naming one is the moment a
    caller believes a table is being consulted, and raising there is what keeps a made-up
    number out of a geometry that is afterwards measured against literature.

    `node_cases.tsv` is not that table either.  Its `d_mm` values are what a built node is
    checked AGAINST; driving the placer from them would make the M6 exit gate compare a
    number with itself.

    A metal-metal bond is not dative, so `DATIVE_LENGTHENING` is not added: the fallback
    is the plain covalent-radii sum, and it reports `estimated` because for a bond whose
    order it does not know that is the honest label (D20).
    """
    if override is not None:
        return Distance(float(override), "override", m1, m2)
    if motif is not None:
        raise AmbiguousSpecError(
            f"there is no M-M distance table to look {motif!r} up in, and {m1}-{m2} does "
            f"not determine one: a quadruply-bonded Mo2 sits at 2.09 A and a Cu2 "
            f"paddlewheel at 2.62 A, and nothing in the two element symbols tells them "
            f"apart. A nucleus is DECLARED, so pass override=<distance in A>; call "
            f"without motif= to accept the covalent-radii estimate instead.")
    d = (COVALENT_RADII.get(m1, DEFAULT_COVALENT)
         + COVALENT_RADII.get(m2, DEFAULT_COVALENT))
    return Distance(d, "covalent-radii", m1, m2)
