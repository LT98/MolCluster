# Experiment log (`/experiments`)

Wet-lab records: what was made, how, what it turned out to be, and how it performed in a
cell. Kept beside the registry, never in it. The page writes only its own folder,
`<data root>/experiments/`, which is git-ignored — the records stay on the machine that
holds them while the code is public.

## What a value can be

| In a field | Means |
|---|---|
| blank | **still to record** — highlighted amber and listed on the right |
| `n/a` | does not apply to this experiment |
| `unknown` | cannot be recovered (never written down at the time) |
| anything else | the value; a wrong kind of value (text where a number goes) is a **problem**, red |

A whole group (a characterization technique, the cells, the washes) can be marked `n/a` or
`unknown` from the selector beside its heading. Nothing is ever filled in for you, and blank
is never read as zero. While an entry is `planned` or `abandoned`, product,
characterization and cells are not asked for.

## How a recipe is laid out

- **Preparations** — each solution, suspension or solid mixture made before the reaction:
  its components (solutes *and* solvents, with mass, volume or mmol, and the planned molar
  ratio), how it was made up (stirring, sonication, …, rpm, time, temperature).
- **Steps** — in order: what was added to what, how (dropwise over 30 min), how long, at
  what temperature, with what mixing and speed, under what atmosphere.
- **Apparatus** — what the method needs: the mill, jar and balls and the grinding liquid for
  ball milling; vessel and fill for solvothermal; thimble and flask for Soxhlet.
- **Workup, product, characterization, cells** — separation and washes (solvent, volume,
  times), drying; yield; PXRD phase, FTIR, CHN, TGA, BET, microscopy, oxidation state; one
  item per coin cell with electrode, cell, protocol and results.

## What is computed (never stored)

With a formula in the reagent table, the right-hand panel shows for every component its
mmol and where it came from: mass, volume × density (× assay for a solution sold by wt%),
or volume × molarity. From those: equivalents against the first solute, flagged red when
they differ from the planned ratio by more than 5 %; the concentration of each solute in its
preparation and in the whole reaction mixture; and η (µL/mg) for liquid-assisted grinding.

## The other views

- **reagents / phases** — the two tables entries refer to. A code typed into an entry that
  is not in the table is a problem. The reagent table shows the molar mass computed from the
  formula next to it.
- **overview** — one row per entry; the same table is rewritten to `INDEX.md` in the log
  folder on every save.
- **compare** — the fields that differ between two entries. Preparations pair by position
  and components by role, so a swapped additive is one row.
- **copy as markdown** (on an entry) — the saved entry as a readable card.

## Saving

`save` (or Ctrl+S) writes the entry. Each save names the version it started from and is
refused if the file changed in between — another tab, another person — rather than
overwriting it; reload and redo the change. **remove** moves the file to
`entries/.deleted/`; nothing is destroyed.
