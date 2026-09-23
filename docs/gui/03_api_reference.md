# HTTP API reference

Every page is built on these. They are usable directly — `curl`, a notebook, a script —
which is the point of keeping run state in the task table rather than in the page.

Two routers: `ui/app.py` (the read-only viewer) and `ui/builder.py` (the only writer).

| Method | Path | Area | Notes |
|---|---|---|---|
| GET | `/` | Registry | the viewer page |
| GET | `/builder` | Builder | the spec builder page |
| GET | `/runs` | Runs | the run inspector page |
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
