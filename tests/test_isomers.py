"""M5/S5: L2 separates what L1 must not, and only when it has the geometry to do it.

The structures under test are BUILT, not hand-written — S4's enumerator puts ligands on
named vertices, so "cis" here means two ligands that really are 90 degrees apart in real
coordinates rather than a fixture asserting its own label.

Amine ligands are conspicuously absent and that is not an oversight: `occlusion` currently
reports every sp3 amine as sterically blocked, so `open_donors()` excludes them and the
classic Pt(NH3)2Cl2 pair cannot be assembled at all. Aqua and chloride stand in. See the
note in the M5 work plan.
"""
from __future__ import annotations

import numpy as np
import pytest

from mofsbu.assembly.construct import ligand_block, metal_block
from mofsbu.assembly.join import join
from mofsbu.geometry.embed import embed_molecule
from mofsbu.graph._types import EdgeType, TypedGraph
from mofsbu.graph.from_mol import mol_from_smiles
from mofsbu.identity.isomers import AMBIGUOUS, CIS_TRANS_SPLIT_DEG
from mofsbu.identity.keys import l1_graph_hash, l2_isomer_tag
from mofsbu.versions import ALGO_VERSIONS

#: Octahedral vertex order from `site_vectors`: +x +y -x -y +z -z.
#: Square planar: +x +y -x -y.  Slots index into those, so 0/1 are adjacent and 0/2 face off.
CIS_SLOTS, TRANS_SLOTS = (0, 1), (0, 2)


@pytest.fixture(scope="module")
def aqua():
    return ligand_block(embed_molecule(mol_from_smiles("O")), charge=0, name="aqua")


@pytest.fixture(scope="module")
def chloride():
    return ligand_block(embed_molecule(mol_from_smiles("[Cl-]")), charge=-1, name="chloride")


def occupy(metal, pairs):
    """Put each ligand on a named vertex, in order."""
    block = metal
    for slot, ligand in pairs:
        vertex = next(s for s in block.open_vacancies() if s.slot == slot)
        block = join(block, ligand, vertex, ligand.open_donors()[0]).block
    return block


def tag(block) -> str:
    return l2_isomer_tag(block.graph, block.geometry)


# ── cis / trans ──────────────────────────────────────────────────────────────

def test_cis_and_trans_share_an_l1_and_split_at_l2(aqua, chloride):
    """The D10 claim, tested: one connectivity, two structures.

    L1 must NOT separate these — a finer L1 would also split things that are genuinely one
    structure — so the whole distinction rests on this tag.
    """
    pt = metal_block("Pt", 2, cn=4, geometry="square_planar", spin_class="ls")
    cis = occupy(pt, [(0, aqua), (1, aqua), (2, chloride), (3, chloride)])
    trans = occupy(pt, [(0, aqua), (2, aqua), (1, chloride), (3, chloride)])

    assert l1_graph_hash(cis.graph) == l1_graph_hash(trans.graph)
    assert tag(cis) != tag(trans)
    assert "cis" in tag(cis) and "trans" in tag(trans)


def test_both_ligand_classes_are_reported(aqua, chloride):
    """MA2B2 says something about the A's and about the B's, and they can differ."""
    pt = metal_block("Pt", 2, cn=4, geometry="square_planar", spin_class="ls")
    cis = occupy(pt, [(0, aqua), (1, aqua), (2, chloride), (3, chloride)])
    assert tag(cis).count("=") == 2
    assert "Cl1=" in tag(cis) and "O3=" in tag(cis)


# ── fac / mer ────────────────────────────────────────────────────────────────

def test_fac_and_mer_share_an_l1_and_split_at_l2(aqua, chloride):
    ru = metal_block("Ru", 3, cn=6, spin_class="ls")
    fac = occupy(ru, [(0, aqua), (1, aqua), (4, aqua),
                      (2, chloride), (3, chloride), (5, chloride)])
    mer = occupy(ru, [(0, aqua), (2, aqua), (1, aqua),
                      (3, chloride), (4, chloride), (5, chloride)])

    assert l1_graph_hash(fac.graph) == l1_graph_hash(mer.graph)
    assert "fac" in tag(fac) and "mer" in tag(mer)


# ── Delta / Lambda ───────────────────────────────────────────────────────────

VERTEX = {"+x": (1, 0, 0), "-x": (-1, 0, 0), "+y": (0, 1, 0),
          "-y": (0, -1, 0), "+z": (0, 0, 1), "-z": (0, 0, -1)}
#: A true C3 propeller: every arm spans one vertex of the {+x,+y,+z} face and one of the
#: opposite face, so the three arms sit at 120 degrees in the plane perpendicular to the
#: (1,1,1) axis.  Built by hand rather than joined because a CHELATE join is not built —
#: S3's `join` places one donor, and `chelate_compatible` only judges whether two could be
#: placed together.
PROPELLER = [("+x", "-y"), ("+y", "-z"), ("+z", "-x")]


def tris_chelate(pairs=PROPELLER, *, mirror: bool = False):
    """A minimal tris(bidentate) octahedron: three O-C-C-O arms on one metal."""
    g = TypedGraph(charge=-3, multiplicity=2, name="tris-chelate")
    metal = g.add_atom("Fe", formal_charge=3, oxidation_state=3, spin_class="ls")
    coords = [np.zeros(3)]
    for first, second in pairs:
        a, b = np.array(VERTEX[first], float) * 2.0, np.array(VERTEX[second], float) * 2.0
        o1, c1 = g.add_atom("O", formal_charge=-1), g.add_atom("C")
        c2, o2 = g.add_atom("C"), g.add_atom("O", formal_charge=-1)
        for u, v in ((o1, c1), (c1, c2), (c2, o2)):
            g.add_bond(u, v, EdgeType.COVALENT)
        g.add_bond(o1, metal, EdgeType.DATIVE)
        g.add_bond(o2, metal, EdgeType.DATIVE)
        outward = (a + b) / float(np.linalg.norm(a + b))
        coords += [a, a + outward * 1.3, b + outward * 1.3, b]
    arr = np.array(coords)
    return g, (arr * np.array([-1.0, 1.0, 1.0]) if mirror else arr)


def test_a_tris_chelate_and_its_mirror_image_are_two_structures():
    """Enantiomers share every achiral descriptor there is, including L1."""
    g, right = tris_chelate()
    mirrored_g, left = tris_chelate(mirror=True)
    assert l1_graph_hash(g) == l1_graph_hash(mirrored_g)
    assert {l2_isomer_tag(g, right), l2_isomer_tag(mirrored_g, left)} == \
           {"6:O4=Delta", "6:O4=Lambda"}


def test_the_handedness_does_not_depend_on_which_end_you_look_from():
    """A chirality that moved when the molecule was rotated would be a point of view.

    The C3 axis is built from a cross product whose sign depends on an arbitrary ordering,
    and that is safe precisely because flipping it swaps which donor of each arm counts as
    "upper" at the same time — the two sign changes cancel.
    """
    g, coords = tris_chelate()
    reference = l2_isomer_tag(g, coords)
    for rotation in (np.diag([-1.0, 1.0, -1.0]),                      # 180 about y
                     np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]]),  # 90 about z
                     np.array([[0, 0, 1.0], [1, 0, 0], [0, 1, 0]])):  # cyclic x->y->z
        assert l2_isomer_tag(g, coords @ rotation.T) == reference


def test_a_chelates_own_bite_is_not_an_isomer():
    """One instance of a class is a ligand binding the way it binds, not a choice made."""
    g = TypedGraph(charge=-1, multiplicity=1, name="one-chelate")
    metal = g.add_atom("Zn", formal_charge=2, oxidation_state=2, spin_class="ls")
    o1, c1 = g.add_atom("O", formal_charge=-1), g.add_atom("C")
    c2, o2 = g.add_atom("C"), g.add_atom("O", formal_charge=-1)
    for u, v in ((o1, c1), (c1, c2), (c2, o2)):
        g.add_bond(u, v, EdgeType.COVALENT)
    g.add_bond(o1, metal, EdgeType.DATIVE)
    g.add_bond(o2, metal, EdgeType.DATIVE)
    coords = np.array([[0.0, 0.0, 0.0],      # Zn
                       [2.0, 0.0, 0.0],      # O
                       [2.6, 1.0, 0.0],      # C
                       [1.4, 1.6, 0.0],      # C
                       [0.0, 2.0, 0.0]])     # O — a 90 degree bite, and still not an isomer
    assert l2_isomer_tag(g, coords) == ""


# ── the boundaries of what it will claim ─────────────────────────────────────

def test_no_geometry_means_no_tag(aqua, chloride):
    """Not a degraded answer: nothing in the graph distinguishes the isomers."""
    pt = metal_block("Pt", 2, cn=4, geometry="square_planar", spin_class="ls")
    cis = occupy(pt, [(0, aqua), (1, aqua), (2, chloride), (3, chloride)])
    assert l2_isomer_tag(cis.graph) == ""
    assert l2_isomer_tag(cis.graph, None) == ""


def test_a_geometry_in_the_wrong_order_is_refused(aqua, chloride):
    pt = metal_block("Pt", 2, cn=4, geometry="square_planar", spin_class="ls")
    cis = occupy(pt, [(0, aqua), (1, aqua), (2, chloride), (3, chloride)])
    with pytest.raises(ValueError, match="label it with confidence"):
        l2_isomer_tag(cis.graph, cis.geometry[:-1])


def test_a_distorted_geometry_declines_to_call_it(aqua, chloride):
    """Near the split the honest output is no tag, not a coin flip."""
    pt = metal_block("Pt", 2, cn=4, geometry="square_planar", spin_class="ls")
    cis = occupy(pt, [(0, aqua), (1, aqua), (2, chloride), (3, chloride)])
    donors = [i for i in cis.graph.nodes()
              if cis.graph.is_metal(0) and i in cis.graph.neighbors(0, EdgeType.DATIVE)]
    metal_pos = cis.geometry[0]
    bent = cis.geometry.copy()
    # Swing one aqua until its angle at the metal sits on the cis/trans boundary.
    axis_ref = bent[donors[0]] - metal_pos
    radius = float(np.linalg.norm(axis_ref))
    theta = np.radians(CIS_TRANS_SPLIT_DEG)
    bent[donors[1]] = metal_pos + radius * np.array(
        [np.cos(theta), np.sin(theta), 0.0])
    assert AMBIGUOUS not in l2_isomer_tag(cis.graph, bent)      # never leaks into a key
    assert "O3=" not in l2_isomer_tag(cis.graph, bent)          # the class simply drops out


def test_the_tag_survives_an_atom_order_shuffle(aqua, chloride):
    """The invariance L1 is already held to, and the reason a tag may be part of a key."""
    pt = metal_block("Pt", 2, cn=4, geometry="square_planar", spin_class="ls")
    cis = occupy(pt, [(0, aqua), (1, aqua), (2, chloride), (3, chloride)])
    nodes = cis.graph.nodes()
    mapping = {old: new for old, new in zip(nodes, reversed(nodes))}
    shuffled = cis.graph.relabel(mapping)
    coords = np.empty_like(cis.geometry)
    for old, new in mapping.items():
        coords[new] = cis.geometry[old]
    assert l2_isomer_tag(shuffled, coords) == tag(cis)


def test_the_recipe_version_moved_off_the_stub():
    assert ALGO_VERSIONS["l2_isomer_tag"] not in ("0-stub", "")
