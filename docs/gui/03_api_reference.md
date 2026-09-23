# HTTP API reference

Every page is built on these. They are usable directly — `curl`, a notebook, a script —
which is the point of keeping run state in the task table rather than in the page.

Two routers: `ui/app.py` (the read-only viewer) and `ui/builder.py` (the only writer).

| Method | Path | Area | Notes |
|---|---|---|---|
| GET | `/` | Registry | the viewer page |
| GET | `/builder` | Builder | the spec builder page |
| GET | `/runs` | Runs | the run inspector page |
| GET | `/graph` | Graph | the formation-energy graph page |
| GET | `/api/capabilities` | Builder | what the pipeline can do; the page renders from this |
| GET | `/api/compute` | Setup | devices present, and the one declared |
| POST | `/api/compute` | Setup | declare the device for this server |
| GET | `/api/databases` | Setup | list registries |
| POST | `/api/databases` | Setup | create one |
| POST | `/api/database/select` | Setup | switch the whole app to one |
| POST | `/api/molecule/preview` | Builder | embed + perceive a SMILES (≤300 atoms) |
| GET | `/api/library` | Builder | saved molecules |
| POST | `/api/library` | Builder | save one |
| GET | `/api/specs` | Builder | list named specs |
| GET | `/api/specs/{source}/{name}` | Builder | load one |
| POST | `/api/spec` | Builder | write one |
| POST | `/api/estimate` | Builder | **the same enumeration submitting performs** |
| POST | `/api/runs` | Runs | queue a run |
| GET | `/api/runs` | Runs | list runs, with computed liveness |
| GET | `/api/runs/{run_id}` | Runs | one run |
| POST | `/api/runs/{run_id}/cancel` | Runs | cancel |
| POST | `/api/runs/{run_id}/resume` | Runs | resume |
| GET | `/api/runs/{run_id}/tasks` | Runs | tasks; `status=`, `code=`, `limit=`, `offset=` |
| GET | `/api/structures` | Registry | listing; filters + `include_hidden=1` |
| GET | `/api/structures/{id}` | Registry | one record |
| GET | `/api/structures/{id}/spec` | Registry | the spec that built it |
| GET | `/api/structures/{id}/origin` | Registry | provenance |
| GET | `/api/structures/{id}/xyz` | Registry | best geometry's coordinates, for the hover preview |
| GET | `/api/structures/{id}/routes` | Graph | incoming edges priced (`routes`), outgoing ones (`consumed_by`), `species`, and `derived` with `inferred=1` |
| GET | `/api/paths/price` | Graph | a whole walk costed; `nodes=` target first, `via=` one edge per gap (`<id>`, `c<id>`, `d<part>`) |
| POST | `/api/structures/{id}/rerun` | Registry | re-run it (**builder router**) |
| POST | `/api/structures/{id}/hide` | Registry | soft-delete (**builder router**) |
| GET | `/api/geometries/{id}/xyz` | Registry | plain-text coordinates, for 3Dmol |
| GET | `/api/filters` | Registry | live filter domains + `reserved_inactive` |
| GET | `/api/build` | Registry | version, branch, commit, dirty-at-start-up |
| GET | `/api/meta` | Registry | db/store paths, read-only flag, `algo_versions` |

Note the two `POST /api/structures/...` rows: they are served by the *builder* router even
though they act on a registry entry, because the viewer has no write path at all.

`GET /api/meta` is the one to reach for when a result looks wrong for its age — it returns
the `algo_versions` rows as stored, which is how a `placement 2` geometry is told apart
from a `placement 3` one.

## The two graph endpoints, in more detail

`GET /api/structures/{id}/routes` answers the same question as the `routes` key of
`GET /api/structures/{id}` and is the one to use in a loop: the full record also lists every
geometry and asks the blob store about each one, which a walk does not want once per hop.
Every entry in `routes` carries `direction: "made_from"`.

`consumed_by` is the other direction: every edge that takes the structure as a **reagent**
(`registry.api.outgoing_routes`), priced the same way, each with `product_structure_id` and
`direction: "consumed_by"`. Only `role = 'reagent'` counts — a structure that appears in an
edge only as its `solvent` or `leaving` term (the H₂O/H₃O⁺ couple of a deprotonation) does not
list it, because that couple carries a proton rather than being built on.

`species` has one row per structure the terms mention — label, formula, charge, `n_metals`,
`n_incoming_routes`, `n_outgoing_routes` — so a candidate step can be named, and its own
reachability in both directions shown, without a request per candidate. `structure` carries the
same columns for the node itself.

`inferred=1` additionally returns `derived`: two-part splits `energy.routes.decompositions`
found by composition arithmetic. They arrive in a key of their own and every entry carries
`origin: "inferred"`, plus `recorded_as` when an edge happens to record the same split. The
separation is not cosmetic — conflating derived with recorded is the thing the reference
scheme exists to prevent.

`GET /api/paths/price?nodes=17,33,74&via=851,339` costs a whole walk. `nodes` is target-first
and `via` has one entry per gap:

| `via` | leg | the edge |
|---|---|---|
| `851` | made from | produces `nodes[i]` from `nodes[i+1]` |
| `c631` | consumed by | produces `nodes[i+1]` from `nodes[i]`, which it takes as a `reagent`; enters the route at **−dE** |
| `d44` | derived | a split of `nodes[i]` into `nodes[i+1]` and structure 44, found by composition |

Both parts of a derived leg are needed because composition fixes the complement only up to its
formula and charge, and the registry holds protomers that share those. A bare id keeps meaning
*made from*, so links from before the `c` form still resolve.

It returns a `y` per node with the target at zero, the `shed` and `free_pieces` multisets that
say what that `y` is an energy *of* (one netted tally: a piece taken up and later given back
cancels), and a `basis` token — two nodes with the same basis carry energies of systems with
the same atoms, and only those may be compared. Each step carries `direction`, `dE` in the
route's direction and `reaction_dE` as the edge records it, plus its own `isodesmic` and
`caveats`.

On the route: `net_equation` (every leg's terms summed in the route's direction with anything on
both sides cancelled — `terms`, `balanced`, `isodesmic`, `caveats`, `dative_delta`, `dE` priced
whole, `agrees_with_steps`, `why_not`); the route's `isodesmic` and `caveats` are the net
equation's, not the union of the steps'. `pivots` lists the node indices where the walk changes
direction, and those nodes carry `pivot: true` and `pivot_kind` (`parent` or `product`); a
route with a pivot has no `worst_step`, and `worst_step_why` says so. `origin` is `recorded`
only when the route is one made-from edge, `inferred` when any leg is derived, and otherwise
`composed`, with `witness` listing each reaction id and the direction it was walked.

A malformed walk is a 400 with a sentence — including an edge that does not join the nodes it
is given for, a consumed-by leg through a proton carrier, and a leg that walks the previous
edge straight back. A step the reference scheme refuses comes back as data, with `y` absent
from that step leftwards rather than zeroed.
