# Production branch

This branch is for **production runs**: campaigns whose results are meant to be used. Tool
building and optimisation happen on `main`. This file exists only here.

## Rules

1. **One direction.** `production` takes merges from `main`. Nothing is merged from `production`
   into `main`.
2. **Fixes go to `main` first.** A defect found during a campaign is fixed on `main` through the
   usual PR, then merged down here. No code is edited on this branch.
3. **Update between campaigns, never during one.** To take new code: merge `main`, run the gate,
   tag, then start.

   ```bash
   git fetch origin
   git merge origin/main
   conda run -n ebu bash scripts/check.sh
   git tag prod-YYYY-MM-<campaign>-<n>
   ```

4. **Own checkout, own data.** This branch lives in its own worktree with its own `data/` folder.
   The scripts put this checkout's `src` on the path, so a run here uses this code and writes only
   here. Give the viewer its own port so it does not collide with a development viewer:

   ```bash
   ./launch/mofsbu.sh --gpu --port 8100 --db data/campaigns/<campaign>/<campaign>.db
   ```

5. **What is committed here:** this file and the campaign log below. Databases, blob stores and
   results are never committed (`/data/*` is ignored).
6. **Specs for unpublished targets are not committed while the repository is public.** Keep them
   with the campaign's data under `data/campaigns/<campaign>/`, which git ignores. A spec is
   committed under `data/reference/` only once its target may be named publicly.

## Starting a campaign

1. Merge `main`, run the gate, tag (rule 3).
2. Create `data/campaigns/<campaign>/` and put the spec there.
3. Run into a database of its own in that folder.
4. Add a row to the log.

## Campaign log

| Tag | Date | Campaign | Commit | Database | Notes |
|---|---|---|---|---|---|
| — | — | — | — | — | none yet |
