# The formation-energy graph — `/graph`

Pick a structure, walk its provenance backwards a step at a time, and read what each route
cost on one chart. The provenance panel on the registry page answers *what was this made
from, and what did that step cost*; this page answers the question that needs several steps
composed — **which of the routes to this structure is the better one**.

Read-only, like the rest of the viewer. Nothing on this page writes.

## Getting here

- **"graph this →"** under the provenance heading on the registry page, with a row selected.
- The **energy graph** tab, then a structure id in the box at the top left.
- A link: `/graph#target=17`.

## The walk

The **standing on** box is the node the picker is relative to: its label, id, formula, charge,
multiplicity, how many edges reach it, and a turnable 3D view of it.

**Precursor steps** lists every recorded edge that reaches that node, collapsed so that one
chemistry is one row — same reagents, same thing added — with a count for the edges behind it.
Columns are `dE · added · source · rung · caveat · count`, ordered by dE with unpriceable rows
next and `place` edges last. Clicking a row steps back to the **source**, the piece being built
on; clicking a **species** steps into that species instead, so you can follow the ligand rather
than the metal it was added to. `↩16` after a species says sixteen edges reach it, `·end` says
none do and the row will not move.

Which reagent is the source is a fact about the species rather than a reading of its label: the
registry knows which one carries the metal and how many atoms each has. Where two precursors
would print the same short name — two protomers differ only in which oxygen lost its proton, and
the short form drops exactly that — the id is shown, and only there.

The **trail** above is the route so far, target first. Click any crumb, or any level on the
chart, to move the picker there.

## Several routes, and forking

Standing at the end of a route, a step extends it. Standing in the middle of one — having
clicked back to a crumb or a level — a step that differs from the one already there **forks** a
new route that keeps everything to its right. That shared right-hand side is exactly the
condition under which two routes can be read against each other.

Colour follows the branching rather than a fixed palette: the routes form a tree rooted at the
target, each divergence splits its parent's arc of the hue circle among its children, and a
route takes the middle of the arc it lands in. So how close two colours are *is* how late the
two routes parted. With **colour** off they keep dash patterns instead, which is also what makes
the chart survive being printed.

Each route has a line in **routes**: its swatch, where it bottoms out, its total, its rung, its
worst single step, and a × to take it off the chart. **drop node** on the *standing on* heading
removes a node and with it every route that goes through it — a route is a claim about how the
target is reached, and leaving the remainder of one behind would draw a route nobody built.

## Reading the chart

The **target sits at zero**, on the right, and routes reach leftwards. Every route ends at the
one node they all share, which is what pins them to a common point; their left ends are wherever
each bottoms out and are not required to agree. `y` falls leftwards by each step's dE, so the
precursors of an exothermic assembly sit *above* the thing they fall to.

Because every step balances, at each node

```
node + pieces not yet consumed  ==  target + what has left
```

so a `y` is the energy of a system with exactly the target's atoms, relative to the target in
the same accounting. That is what lets it be read down a column.

**Two routes may be compared only where that accounting matches** — same target, and the same
things shed on the way to it. The identity above then forces the unconsumed pieces to agree too,
so the shed multiset alone decides it. It is reported per node as `basis`, and the legend badges
a route whose basis differs from the others. An assembly step sheds nothing and a deprotonation
sheds `n H3O+`, so what this catches in practice is a protonated route set beside a deprotonated
one — where the vertical gap between the curves is not a comparison of anything. See
[B19](../BUGS.md#b19) for a worked case and for what the chart still does not say well.

A node with no number gets a **dashed column** and the words *no number*, not a position:
everything to the left of an unpriceable step has no `y`, and putting it at the last known
height would be inventing data.

## Derived steps

A `place` edge records a construction rather than a reaction — the runner built the whole
coordination sphere at once — so it cites no reagents and cannot be priced (**C15**). Those rows
read `— (built whole)` and are not steps.

**derived too**, beside the *precursor steps* heading, additionally offers splits found by
composition arithmetic over the registry: pairs of existing structures whose atoms and charge
add up to this one. They are dashed, badged `inferred`, and never merged into the recorded list.
A node whose recorded edges all cite nothing asks for them automatically, because derivation is
then the only answer to what it is made of. A derived row that happens to match a recorded edge
says so.

## Controls

| | |
|---|---|
| **dE** | the step's energy change, on every connector |
| **level value** | each level's own coordinate on its bar — reading a height off the axis by eye is an estimate |
| **added / shed** | what each step builds in and what it loses, beside the connector |
| **colour** | off = grey plus dash patterns |
| **zero at root** | each route re-zeroed on its own left end. Offered and badged: routes bottom out at different species, so those zeros are different systems |
| **text** | S…XXL; every label on the chart and the margins with them |
| **derived too** | see above |

The three panes resize by dragging the line between them — double-click to reset, arrow keys
when focused. Pane sizes and control settings are per browser and deliberately **not** in the
link: a shared link should give the other reader their own chart of the same routes.

## What is in the link

`#target=17&p=851-33.339-74&p=956-38.392-77` — the target, then one `p` per route, each a chain
of legs. A leg is `<reaction>-<node>`, or `d<other part>-<node>` along a derived split. Both
halves are needed: the edge, because an id pair does not say *which* edge joins them; the node,
because one edge names two reagents and stepping into either is a different walk.

Replaying a link rebuilds each route against the registry in front of it rather than trusting
it. A leg that no longer holds stops that route and says which one, instead of quietly
producing a shorter walk that would look deliberate.

## Where the numbers come from

`energy/routes.py` prices one edge. `pathways/route.py` composes them — the running sum, the
check that consecutive steps actually join, the shed and unconsumed multisets, and the basis.
That is in Python rather than in the page because `scripts/check.sh` cannot execute a line of
JavaScript, and a cumulative sum no test can reach is the kind of number this project refuses
to show. `tests/test_path_energy.py` is where it is pinned.
