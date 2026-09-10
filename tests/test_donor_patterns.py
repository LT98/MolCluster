"""Every donor pattern must actually fire, on the atoms it is named for.

This file exists because four patterns were added that could never match anything.  The
donor is the LAST atom of a SMARTS match (`_find_labile`), so the natural spelling of a
multi-atom acid — `[CX3](=[OX1])[NX3H1][CX3]=[OX1]` for an imide — ends on a carbonyl
oxygen with no hydrogen.  The pattern matches the molecule, the `h_idx` lookup then
throws the match away, and the list gains an entry that reads like coverage and finds
nothing.  Nothing failed; the donor simply was not there.

A pattern with no probe molecule is a test failure too.  Adding a SMARTS without adding
the molecule it is supposed to recognise is how the four dead ones got in.

`DELOCALISED_GROUPS` is held to the same discipline, plus one it does not share: a group
has to perceive the SAME donor set in every resonance form of itself, because that is the
entire reason it is matched as a group instead of an atom at a time.
"""
from __future__ import annotations

from collections import Counter

import pytest
from rdkit import Chem

from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.sites.perception import (
    DELOCALISED_GROUPS, LABILE_DONOR_PATTERNS, deprotonate, find_donor_sites,
)
from mofsbu.sites.protomers import labile_sites

# One molecule per pattern that the pattern MUST recognise.
PROBES: dict[str, str] = {
    "carboxylate_O":       "CC(=O)O",              # acetic acid
    "phenolate_O":         "Oc1ccccc1",            # phenol
    "enol_O":              "OC=C",                 # vinyl alcohol / hydroxyquinone OH
    "sulfonate_O":         "CS(=O)(=O)O",          # methanesulfonic acid
    "sulfinate_O":         "CS(=O)O",              # methanesulfinic acid
    "phosphonate_O":       "CP(=O)(O)O",           # methylphosphonic acid
    "boronate_O":          "CB(O)O",               # methylboronic acid
    "acyl_sulfonamide_N":  "CC(=O)NS(=O)(=O)C",    # N-acylmethanesulfonamide
    "sulfonamide_N":       "CS(=O)(=O)NC",         # N-methylmethanesulfonamide
    "azolate_N":           "c1cc[nH]c1",           # pyrrole
    "hydroxamate_O":       "CC(=O)NO",             # acetohydroxamic acid
    "oxime_O":             "CC(C)=NO",             # acetone oxime
    "nhydroxy_O":          "CN(C)O",               # N,N-dimethylhydroxylamine
    "imide_N":             "O=C1CCC(=O)N1",        # succinimide
    "thiophosphate_S":     "COP(=O)(OC)S",         # O,O-dimethyl thiophosphate
    "selenol_Se":          "C[SeH]",               # methaneselenol
    "peroxy_O":            "COO",                  # methyl hydroperoxide
    "thiolate_S":          "CS",                   # methanethiol
    "alkoxide_O":          "CCO",                  # ethanol
    "hydrohalide_X":       "Cl",                   # hydrogen chloride
}

# Element the donor atom must be, taken from the pattern's own name.
HALOGENS = {"F", "Cl", "Br", "I"}


def _expected_element(name: str) -> set[str]:
    suffix = name.rsplit("_", 1)[-1]
    return HALOGENS if suffix == "X" else {suffix}


def test_every_pattern_has_a_probe_molecule():
    """A SMARTS with no example is untested by construction."""
    named = [name for name, _, _ in LABILE_DONOR_PATTERNS]
    assert set(named) == set(PROBES), (
        "patterns without a probe: " + str(sorted(set(named) - set(PROBES)))
        + "; probes without a pattern: " + str(sorted(set(PROBES) - set(named))))


def test_no_pattern_parses_to_none():
    for name, smarts, _ in LABILE_DONOR_PATTERNS:
        assert Chem.MolFromSmarts(smarts) is not None, f"{name}: unparseable SMARTS"


@pytest.mark.parametrize("name", list(PROBES))
def test_pattern_fires_on_its_probe(name):
    """The pattern finds a site, and the site is on the element the name promises."""
    smiles = PROBES[name]
    sites = [s for s in find_donor_sites(mol_from_smiles(smiles)) if s.donor_type == name]
    assert sites, (
        f"{name} matched nothing in {smiles}. Either the SMARTS is wrong, an earlier "
        f"pattern claimed the atom first, or — the usual cause — the LAST atom of the "
        f"match carries no hydrogen, so `_find_labile` discarded it. Put the donor last, "
        f"with `$(...)` for its environment if need be.")
    mol = mol_from_smiles(smiles)
    for site in sites:
        element = mol.GetAtomWithIdx(site.idx).GetSymbol()
        assert element in _expected_element(name), (
            f"{name} put its site on {element}{site.idx} in {smiles}; the name says "
            f"{sorted(_expected_element(name))}. A mistyped donor silently picks the "
            f"wrong entry out of IDEAL_MDA_ANGLE and the wrong binding modes.")


def test_no_carbon_donors():
    """A carbanion has no lone pair to point at a metal, and no frame model here.

    1,3-diketones and nitroalkanes are real C-H acids; their DONOR is the delocalised
    oxygen (enolate, nitronate), which `enol_O` already covers on the enol tautomer.
    """
    for name, _, _ in LABILE_DONOR_PATTERNS:
        assert not name.endswith("_C") and not name.endswith("_CH"), (
            f"{name}: carbon is not a donor type this package can place")


def test_peroxy_takes_the_hydroxyl_oxygen_not_the_carbon():
    """The regression that motivated the element check.

    Spelled `[OX2H1][OX2][#6]` the match ends on carbon, and methyl hydroperoxide came
    back with its METHYL carbon as the activatable site — deprotonating CH3 to a
    carbanion.
    """
    sites = find_donor_sites(mol_from_smiles("COO"))
    peroxy = [s for s in sites if s.donor_type == "peroxy_O"]
    assert peroxy and all(
        mol_from_smiles("COO").GetAtomWithIdx(s.idx).GetSymbol() == "O" for s in peroxy)


def test_a_halide_survives_its_own_activation():
    """HX -> X- must still perceive a donor, or activation destroys the site (D5).

    `_classify_anionic` handled O, S and N and returned None for everything else, so the
    chloride produced by deprotonating HCl was invisible the moment it was read back.
    """
    hcl = mol_from_smiles("Cl")
    sites = labile_sites(hcl)
    assert sites, "HCl has no labile site"
    anion, _ = deprotonate(hcl, sites)
    assert Chem.MolToSmiles(anion) == "[Cl-]"
    recalled = find_donor_sites(anion)
    assert recalled, "the chloride left by activation perceives no donor at all"
    assert recalled[0].donor_type == "halide_Cl"
    assert recalled[0].labile is False


def test_a_free_halide_co_ligand_is_a_donor():
    """`[Cl-]` handed straight in as a co-ligand is the normal way a halide arrives."""
    sites = find_donor_sites(mol_from_smiles("[Cl-]"))
    assert [s.donor_type for s in sites] == ["halide_Cl"]


# ── the delocalised groups ───────────────────────────────────────────────────
#
# One anionic probe per group, and how many oxygens the whole group must present.  The
# count is the assertion that matters: a group that reports two of its three oxygens has
# not been matched as a group at all, it has been classified an atom at a time and the
# third one fell through whichever valence rule its bond order happened to miss.
GROUP_PROBES: dict[str, tuple[str, int]] = {
    "carboxylate_O": ("CC(=O)[O-]",      2),   # acetate
    "sulfonate_O":   ("CS(=O)(=O)[O-]",  3),   # methanesulfonate
    "sulfinate_O":   ("CS(=O)[O-]",      2),   # methanesulfinate
    "phosphonate_O": ("CP(=O)([O-])[O-]", 3),  # methylphosphonate
    "nitro_O":       ("C[N+](=O)[O-]",   2),   # nitromethane
    "boronate_O":    ("CB([O-])[O-]",    2),   # methylboronate
}

#: Spellings of one anion that differ ONLY in where the bond orders and the minus sign
#: went.  D15 hashes these identically, so perception has to type them identically.
RESONANCE_PAIRS: list[tuple[str, str, str]] = [
    ("acetate",      "CC(=O)[O-]",         "CC([O-])=O"),
    ("benzoate",     "[O-]C(=O)c1ccccc1",  "O=C([O-])c1ccccc1"),
    ("sulfonate",    "CS(=O)(=O)[O-]",     "CS([O-])(=O)=O"),
    ("phosphonate",  "CP(=O)([O-])[O-]",   "CP([O-])([O-])=O"),
    ("nitro",        "C[N+](=O)[O-]",      "C[N+]([O-])=O"),
    ("oxalate",      "[O-]C(=O)C(=O)[O-]", "O=C([O-])C([O-])=O"),
]


def test_every_group_has_a_probe_molecule():
    named = [name for name, _, _ in DELOCALISED_GROUPS]
    assert set(named) == set(GROUP_PROBES), (
        "groups without a probe: " + str(sorted(set(named) - set(GROUP_PROBES)))
        + "; probes without a group: " + str(sorted(set(GROUP_PROBES) - set(named))))


def test_no_group_smarts_parses_to_none():
    for name, smarts, _ in DELOCALISED_GROUPS:
        assert Chem.MolFromSmarts(smarts) is not None, f"{name}: unparseable SMARTS"


@pytest.mark.parametrize("name", list(GROUP_PROBES))
def test_a_group_types_every_one_of_its_oxygens_the_same(name):
    smiles, n_oxygens = GROUP_PROBES[name]
    found = Counter(s.donor_type for s in find_donor_sites(mol_from_smiles(smiles)))
    assert found == {name: n_oxygens}, (
        f"{smiles} perceived {dict(found)}, want {n_oxygens} x {name}. The oxygens of a "
        f"delocalised group are equivalent; a stray `carbonyl_O` or `alkoxide_O` here "
        f"means the group was classified one atom at a time, off its bond orders.")


@pytest.mark.parametrize("name", list(GROUP_PROBES))
def test_a_group_carries_one_charge_across_all_its_oxygens(name):
    """Charge belongs to the group, not to whichever oxygen the minus sign landed on."""
    smiles, _ = GROUP_PROBES[name]
    sites = find_donor_sites(mol_from_smiles(smiles))
    assert len({s.charge_after for s in sites}) == 1, (
        f"{smiles}: {[(s.idx, s.charge_after) for s in sites]} — the oxygens of one group "
        f"disagree about their charge, so the stored catalog depends on the spelling.")


@pytest.mark.parametrize("name,left,right", RESONANCE_PAIRS,
                         ids=[c[0] for c in RESONANCE_PAIRS])
def test_two_resonance_forms_perceive_the_same_donors(name, left, right):
    """The invariance itself, at the free-ligand level.

    Indices are not comparable between the two spellings (the atoms are written in a
    different order), so the claim is on the multiset of types and charges: same donors,
    same many, same charge, whichever form RDKit is holding.
    """
    def described(smiles: str) -> Counter:
        return Counter((s.donor_type, s.charge_after, s.labile)
                       for s in find_donor_sites(mol_from_smiles(smiles)))

    assert described(left) == described(right), (
        f"{name}: {left} and {right} are one identity under D15 and perceive differently")


@pytest.mark.parametrize("smiles,forbidden", [
    ("CC(=O)OC",     "carboxylate_O"),   # methyl acetate: -OR is not a group oxygen
    ("COP(=O)(OC)S", "phosphonate_O"),   # phosphate diester: only one terminal O left
    ("CC(C)=O",      "carboxylate_O"),   # acetone: one oxygen is not a carboxylate
    ("CS(=O)C",      "sulfinate_O"),     # DMSO is a sulfoxide, not a sulfinate
])
def test_an_esterified_oxygen_does_not_pull_its_neighbour_into_a_group(smiles, forbidden):
    """The terminal-oxygen rule, which is what keeps the group table from over-reaching.

    An -OR oxygen has a second heavy neighbour, so it is not part of the delocalised set
    and does not count towards the group's oxygen tally.  Without that, every ester would
    perceive as a carboxylate and its alkyl oxygen would be handed a charge it has not got.
    """
    found = {s.donor_type for s in find_donor_sites(mol_from_smiles(smiles))}
    assert forbidden not in found, f"{smiles} was swept into a group: {sorted(found)}"
