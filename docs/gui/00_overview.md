# MolCluster's HTML GUI — overview

Four pages served by one FastAPI app. Everything below is what the code actually exposes
as of the formation-energy graph; anything not yet reachable is in `04_limits.md` rather
than described as if it worked.

## Launching it

```bash
conda run -n ebu python -m mofsbu.ui                 # http://127.0.0.1:8000
conda run -n ebu python -m mofsbu.ui --db other.db --port 8001
```

Flags: `--host` (defaults to `127.0.0.1`, localhost only), `--port` (8000), `--reload`
for development, `--log-level` (defaults to `warning` so an open run does not bury the
console). `launch/mofsbu.sh`, `launch/mofsbu.bat` and
`launch/install-desktop-entry.sh` wrap this for a double-click start.

On an interactive terminal, start-up asks two questions before serving: whether to use
CUDA for this server's `ml_go`/`xtb_go` runs, and whether to allow more than one
in-process worker. Neither is inferred from the hardware — the device is *declared*, and
the page shows which declaration is in force.

## The four pages

| Page | Path | What it is for |
|---|---|---|
| **Registry** | `/` | Browse, filter and inspect everything the pipeline has built |
| **Builder** | `/builder` | Compose a spec and queue a run. Never edit JSON by hand |
| **Run inspector** | `/runs` | Watch a run, and read why anything was refused |
| **Energy graph** | `/graph` | Walk a structure's provenance backwards and compare the routes ([05](05_graph.md)) |

They share one tab strip (`chrome.css` / `chrome.js`), and each page names the database it
is pointed at, so two servers on two registries are not confusable. The registry page and the
graph also share `routes.js` — how to read a provenance edge, what to call a species, when two
edges are one chemistry — because two answers to any of those would be two accounts of one
registry.

## The one boundary worth knowing

**The viewer cannot write, by construction.** `ui/app.py` opens every connection
`file:...?mode=ro` with `PRAGMA query_only`, and reads the `v_structures` view rather than
base tables. There is no write path in that module at all.

**The builder writes exactly three things:** spec files and library entries on disk, and
rows in `runs`/`tasks` when you submit. It never writes a structure, a geometry or a
label — those only ever arrive through `registry.api`, from a worker.

That split is what lets a page queue an hour-long job without becoming stateful: run state
lives in the task table, so you can close the page, reopen it, or replace it with a
notebook and lose nothing.

A footer stamp (`/api/build`) reports version, branch, commit and whether the working tree
was dirty at start-up — served from the app rather than baked into the static HTML, so it
describes the process answering the request rather than whenever the page was written.
