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
| POST | `/api/estimate` | Builder | **the same enumeration submitting performs**; `capped` lists a cap that cut the ladder, with its uncapped size |
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
| GET | `/api/structures/{id}/routes` | Graph | incoming edges priced, `species`, and `derived` with `inferred=1` |
| GET | `/api/paths/price` | Graph | a whole walk costed; `nodes=` target first, `via=` one edge per gap |
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
geometry and asks the blob store about each one, which a backwards walk does not want once
per hop. It adds `species`, one row per structure the terms mention — label, formula, charge,
`n_metals`, `n_incoming_routes` — so a candidate step can be named, and its own reachability
shown, without a request per candidate.

`inferred=1` additionally returns `derived`: two-part splits `energy.routes.decompositions`
found by composition arithmetic. They arrive in a key of their own and every entry carries
`origin: "inferred"`, plus `recorded_as` when an edge happens to record the same split. The
separation is not cosmetic — conflating derived with recorded is the thing the reference
scheme exists to prevent.

`GET /api/paths/price?nodes=17,33,74&via=851,339` costs a whole walk. `nodes` is target-first
and `via` has one entry per gap: a reaction id, or `d<other part>` for a step along a derived
split. Both parts of a derived leg are needed because composition fixes the complement only up
to its formula and charge, and the registry holds protomers that share those.

It returns a `y` per node with the target at zero, the `shed` and `free_pieces` multisets that
say what that `y` is an energy *of*, and a `basis` token — two nodes with the same basis carry
energies of systems with the same atoms, and only those may be compared. A malformed walk is a
400 with a sentence; a step the reference scheme refuses comes back as data, with `y` absent
from that step leftwards rather than zeroed.
