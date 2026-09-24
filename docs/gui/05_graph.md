# The formation-energy graph — `/graph`

Pick a structure, walk its provenance a step at a time — back to what it was made from, or
forward to what consumes a node — and read what each route
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

Two tables offer the next hop, one per direction.

**made from · step back** lists every recorded edge that reaches that node, collapsed so that
one chemistry is one row — same reagents, same products — with a count for the edges behind it.
Columns are `dE · added · source · rung · caveat · count`, ordered by dE with unpriceable rows
next and `place` edges last. Clicking a row steps back to the **source**, the piece being built
on; clicking a **species** steps into that species instead, so you can follow the ligand rather
than the metal it was added to.

**consumed by · step forward** lists every recorded edge that takes the node as a **reagent**,
columns `dE (edge) · with · makes · rung · caveat · count`. Clicking a row steps forward to
what the edge **makes**; a co-reagent in the *with* column is not a step, because the edge
leads only to its product. The dE shown is the edge's own, forward as recorded — walked this
way it enters the route with the opposite sign (see below). Only `reagent` terms count: the
water and H3O⁺ that carry a proton through a deprotonation edge are not built on, so they are
never offered as a way forward and would otherwise be a hub leading to every deprotonated
species in the registry.

Neither table lists the edge you arrived by — stepping back along it is the crumb to your
right. After a species, `↩16` says sixteen other edges make it, `↪3` that three other edges
consume it, and `·end` that no other edge does either and the row leads nowhere further.

Which reagent is the source is a fact about the species rather than a reading of its label: the
registry knows which one carries the metal and how many atoms each has. Where two precursors
would print the same short name — two protomers differ only in which oxygen lost its proton, and
the short form drops exactly that — the id is shown, and only there.

The **trail** above is the route so far, target first. `←` before a crumb means the node before
it was made from it, `→` that the node before it was consumed by the edge that made it. Click
any crumb, or any level on the chart, to move the picker there.

## Walking forward: composing an exchange

A route may change direction. Walking back from a target to a shared parent and then forward
again composes a **ligand exchange** out of recorded assemblies. On `data/mvp_ni_thq_cl.db`:

```
17 ← 73 ← 81 → 80 → 24     via 306, 616, c631, c316
```

reads, in the direction the chemistry runs, NiCl₂(H₂O)₂ [24] loses Cl⁻ twice to Ni(H₂O)₂²⁺
[81], which takes up tHQ⁻ twice to Ni(tHQ⁻)₂(H₂O)₂ [17]. A consumed-by leg is the edge **run in
reverse** along the route, so it contributes **−dE**, and what the edge added is released on the
way while anything it shed is taken up. Here the legs are −4.707, −12.573, +12.884, +6.906 eV,
and the route comes to **+2.509 eV**.

**A route is judged on its net equation, not on its legs.** Summing every leg's terms in the
route's direction and cancelling whatever is on both sides gives
`NiCl₂(H₂O)₂ + 2 tHQ⁻ → Ni(tHQ⁻)₂(H₂O)₂ + 2 Cl⁻` — balanced, **isodesmic**, no caveats, and
priced whole it is the same +2.509 eV. Every leg on its own carries `coordination_change`; the
exchange as a whole does not, and that is the point of composing it. The legs keep their own
diagnostics; the route's `isodesmic` and caveats come from the net equation. The spectator
tally is netted the same way, so a ligand taken up on one leg and given back on another cancels
rather than appearing on both sides of `basis`.

**A pivot is bookkeeping, not an intermediate.** Where the walk turns — #81 above — the node is
drawn as a **hollow bar** tagged *pivot · not a barrier*. Its `y` is the energy of the bare
parent plus every ligand on both sides, which is how the composition is booked, not a height the
exchange climbs: nothing here says the exchange is dissociative. A turn the other way (forward
to a shared product, then back) is a pivot too. **A pivot's `y` must never feed a barrier proxy
(C6)**, and the legend reports no *worst step* for a route with one, because the legs either
side of it are halves of one exchange.

The legend badges such a route **composed**: it is not one recorded edge and not a derived split,
so it is neither `recorded` nor `inferred`. Its tooltip names the reaction ids it was composed
from (the witness), the net equation, and the net equation's verdict.

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
worst single step (none for a route with a pivot), a **composed** badge where it turns, and a ×
to take it off the chart. **drop node** on the *standing on* heading
cuts every route that passes through the node: each keeps its part from the target up to the node
before it, and the node and everything beyond it go. What is left is still a route that was walked
hop by hop, only shorter, and it is priced again. Routes the cut makes identical merge into one; a
route left as the target alone is removed. The target itself cannot be dropped.

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

**derived too**, beside the *made from* heading, additionally offers splits found by
composition arithmetic over the registry: pairs of existing structures whose atoms and charge
add up to this one. They are dashed, badged `inferred`, and never merged into the recorded list.
A node whose recorded edges all cite nothing asks for them automatically, because derivation is
then the only answer to what it is made of. A derived row that happens to match a recorded edge
says so.

## Medium, and where a released proton goes

Two selects in the toolbar are inputs to **every** route on the chart at once, because routes
priced under different ones compare nothing:

- **medium** — `gas`, or a continuum the registry holds corrections in (`alpb:water` once
  `runner.finalise_run` or `scripts/solvate.py` has run). A medium is always `model:solvent`;
  each energy is `E + dG_solv` on its own geometry, and a route with any term lacking a
  correction is refused rather than half-solvated. The legend says `in <medium>`.
- **H⁺ to** — `water (H3O+)`, or a free base whose conjugate acid is also in the registry
  (acetate → acetic acid in the Ni(OAc)₂ run; the chloride run offers none, since HCl is not a
  species). Each step that releases n H₃O⁺ also runs n × (H₃O⁺ + A⁻ → H₂O + HA), priced in the
  same medium — exact by Hess's law, and visible: the step's terms, the basis (`HOAc` instead
  of `H3O+`) and the net equation all carry it, and the legend badges the route `H⁺ → <base>`.

The route's caveats — `charge_separation` among them — are listed on its legend line whenever
there are any, isodesmic or not.

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

`#target=17&p=851-33.339-74&p=306-73.616-81.c631-80.c316-24` — the target, then one `p` per
route, each a chain of legs. A leg is `<reaction>-<node>` made from a recorded edge,
`c<reaction>-<node>` consumed by one, or `d<other part>-<node>` along a derived split. A bare
reaction id still means *made from*, so every link written before the `c` form reads as it did.
Both halves are needed: the edge, because an id pair does not say *which* edge joins them; the node,
because one edge names two reagents and stepping into either is a different walk.

Replaying a link rebuilds each route against the registry in front of it rather than trusting
it. A leg that no longer holds stops that route and says which one, instead of quietly
producing a shorter walk that would look deliberate.

**The walk survives leaving the page.** The tabs link to a bare `/graph`, so the page also keeps
the link (and the chosen medium and H⁺ destination) in this browser's `localStorage`, **per
database**: a structure id means a different species in another registry. Opening `/graph` with no
hash puts back the last walk for the active database; a link that carries its own hash always
wins. It is per browser and per viewer address: a different machine, a private window or a
different forwarded port (`localhost:8000` vs `localhost:8013`) starts clean.

## Where the numbers come from

`energy/routes.py` prices one edge. `pathways/route.py` composes them — the running sum with
each leg's sign, the check that consecutive steps actually join, the netted spectator tally, the
basis, the pivots, and the net equation with its quality.
That is in Python rather than in the page because `scripts/check.sh` cannot execute a line of
JavaScript, and a cumulative sum no test can reach is the kind of number this project refuses
to show. `tests/test_path_energy.py` is where it is pinned.
