-- mofsbu registry, schema version 1.
--
-- Tables for M3 are LIVE.  Tables and columns for M4-M8 are RESERVED: created now,
-- empty until their milestone, so that later work is a fill-in rather than a migration.
-- A reserved column is nullable or defaulted and is documented with the milestone that
-- populates it.  Nothing reads a reserved column before then.
--
-- Ground rules honoured here (docs/PLAN_implementation.md §0):
--   * identity is the (L0, L1, L2) key — never a filename, never a build order
--   * every computed number references a `methods` row (fidelity + method + version)
--   * the canonical atom order is STORED, not recomputed, so an algorithm-version bump
--     cannot silently re-key existing site annotations

PRAGMA foreign_keys = ON;

-- ── meta ─────────────────────────────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS migrations (
    version     INTEGER PRIMARY KEY,
    applied_at  TEXT    NOT NULL,
    note        TEXT    NOT NULL DEFAULT ''
);

-- Which algorithm versions this database's rows were produced with.  A mismatch
-- against mofsbu.versions.ALGO_VERSIONS means stored keys are stale.
CREATE TABLE IF NOT EXISTS algo_versions (
    name        TEXT PRIMARY KEY,
    version     TEXT NOT NULL,
    recorded_at TEXT NOT NULL
);

-- Every energy, ease value and barrier points at one of these.  No bare floats.
--
-- `code` is the PACKAGE and `method` is the THEORY, and the difference carries weight:
-- MACE-MP-0 and MACE-OMOL-0 are both code='mace' and are not the same theory.  MP-0 is
-- trained on Materials Project relaxations and is blind to charge and spin; OMOL-0 is
-- trained on OMol25 (wB97M-V/def2-TZVPD) and takes total charge and spin multiplicity
-- as inputs.  Their absolute energies are on different scales and must never be
-- subtracted from one another.  Nothing here needed migrating for the second model to
-- arrive: `method`, `code_version` (which pins the checkpoint) and `extras_json` (which
-- carries `training_set`, and `charge_blind`/`spin_blind` when they apply) are all in
-- the UNIQUE key already, so the two models land in two rows on their own.  What DID
-- need changing is every query that ordered geometries by energy across method rows —
-- see `_refresh_best_geometry` and `energy.reference._energy_row`.
CREATE TABLE IF NOT EXISTS methods (
    id           INTEGER PRIMARY KEY,
    code         TEXT NOT NULL,                 -- 'tblite' | 'mace' | 'heuristic' | 'legacy'
    code_version TEXT NOT NULL,                 -- pins the checkpoint too: '0.3.14/extra_large'
    method       TEXT NOT NULL,                 -- 'GFN2-xTB' | 'MACE-MP-0' | 'MACE-OMOL-0'
    solvent      TEXT,
    charge       INTEGER,
    multiplicity INTEGER,
    extras_json  TEXT NOT NULL DEFAULT '{}',
    UNIQUE (code, code_version, method, solvent, charge, multiplicity, extras_json)
);

-- ── structures: identity (LIVE, M3) ──────────────────────────────────────────

CREATE TABLE IF NOT EXISTS structures (
    id                   INTEGER PRIMARY KEY,

    -- identity: the composite L0-L3 key.  L2/L3 are stubs until M5 but the columns
    -- and the UNIQUE constraint are final, so filling them in is not a migration.
    l0_composition       TEXT NOT NULL,
    l1_graph_hash        TEXT NOT NULL,         -- sha256 of the canonical certificate (D16)
    l2_isomer_tag        TEXT NOT NULL DEFAULT '',   -- RESERVED M5
    block_id             TEXT NOT NULL,
    wl_index             TEXT NOT NULL,         -- fast bucket index, NOT identity

    -- the graph and the coordinate system its annotations key against
    typed_graph_hash     TEXT NOT NULL,         -- blob digest of TypedGraph.to_json()
    canonical_order_json TEXT NOT NULL,         -- frozen at insert; see ground rule 6

    -- the recipe versions that produced the keys above
    algo_l0              TEXT NOT NULL,
    algo_l1              TEXT NOT NULL,
    algo_l2              TEXT NOT NULL,
    algo_canon           TEXT NOT NULL,

    -- denormalised query columns ─ axis 1: composition / metal / charge / spin
    formula              TEXT    NOT NULL,
    metals               TEXT    NOT NULL DEFAULT '',   -- ',Cu,Cu,' — LIKE-searchable
    n_metals             INTEGER NOT NULL DEFAULT 0,
    n_atoms              INTEGER NOT NULL,
    net_charge           INTEGER NOT NULL,
    multiplicity         INTEGER NOT NULL,
    ox_states            TEXT    NOT NULL DEFAULT '',

    -- denormalised query columns ─ axis 2: connectivity / binding mode
    max_bridge_class     TEXT    NOT NULL DEFAULT 'none',   -- none|terminal|mu2|mu3|muN
    n_dative_bonds       INTEGER NOT NULL DEFAULT 0,
    has_metal_metal      INTEGER NOT NULL DEFAULT 0,
    -- Set columns use the bracketed-comma convention, as `metals` does: ',a,b,'
    -- so that `LIKE '%,carboxylate_O,%'` cannot match a prefix of another value.
    donor_types          TEXT             DEFAULT NULL,     -- populated by put_sites (M4)
    binding_modes        TEXT             DEFAULT NULL,     -- RESERVED M5
    -- Number of donors perceived on the free structure.  NOT the number of OPEN sites:
    -- that needs site_state, which is not populated yet, so n_open_sites stays NULL
    -- rather than being filled with a number that means something else.
    n_perceived_donors   INTEGER          DEFAULT NULL,
    n_open_sites         INTEGER          DEFAULT NULL,     -- RESERVED: needs site_state

    -- denormalised query columns ─ axis 3: energy / fidelity
    -- Maintained together by registry.api._refresh_best_geometry and nowhere else;
    -- tests assert they never drift from the geometries table they cache.
    best_geometry_id     INTEGER REFERENCES geometries(id) ON DELETE SET NULL,
    best_fidelity        INTEGER DEFAULT NULL,

    -- Display only: DERIVED from the graph (mofsbu.naming.compose_label), non-unique,
    -- regenerable, and NEVER an input to a query.  No query in this codebase may filter,
    -- join or sort on it: searching is structural (see structure_fragments below), and a
    -- label is composed only after rows have been retrieved.  A cached copy lives here so
    -- listings do not have to re-parse every graph; `registry.relabel_all` recomputes it
    -- and `registry.verify` fails if it has drifted from what the graph derives.
    display_label        TEXT NOT NULL DEFAULT '',

    -- Soft delete, and deliberately NOT a real one.  Dropping a structure row cascades
    -- to its geometries, site_catalog, site_state and its `reactions` edges — and
    -- `reactions` is a DAG in which most products have more than one incoming edge, so
    -- deleting "one entry" routinely severs some other route's history.  The registry is
    -- regenerable by construction; cascaded-away provenance is not.  So a mistake is
    -- HIDDEN (filtered out of every listing) rather than removed, and can be brought
    -- back.  Nothing in the identity layer reads this column: a hidden structure still
    -- occupies its L0/L1/L2 key, so re-deriving it is still recognised under D2 rather
    -- than silently inserted a second time.
    hidden               INTEGER NOT NULL DEFAULT 0,
    hidden_at            TEXT             DEFAULT NULL,
    hidden_reason        TEXT    NOT NULL DEFAULT '',

    created_at           TEXT NOT NULL,
    UNIQUE (l0_composition, l1_graph_hash, l2_isomer_tag)
);

CREATE INDEX IF NOT EXISTS ix_structures_l1        ON structures (l1_graph_hash);
CREATE INDEX IF NOT EXISTS ix_structures_wl        ON structures (wl_index);
CREATE INDEX IF NOT EXISTS ix_structures_l0        ON structures (l0_composition);
CREATE INDEX IF NOT EXISTS ix_structures_formula   ON structures (formula);
CREATE INDEX IF NOT EXISTS ix_structures_nmetals   ON structures (n_metals, net_charge);
CREATE INDEX IF NOT EXISTS ix_structures_bridge    ON structures (max_bridge_class);
CREATE INDEX IF NOT EXISTS ix_structures_fidelity  ON structures (best_fidelity);

-- ── geometries: the fidelity ladder (LIVE, M3) ───────────────────────────────

CREATE TABLE IF NOT EXISTS geometries (
    id                   INTEGER PRIMARY KEY,
    structure_id         INTEGER NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    l3_conformer_id      TEXT    NOT NULL DEFAULT '',   -- RESERVED M5
    coords_hash          TEXT    NOT NULL,              -- blob digest of the .xyz
    n_atoms              INTEGER NOT NULL,
    -- A RUNG, not a theory: one rung can be served by more than one method (ML is
    -- served by both MACE models).  `method_id` says which, and is part of the UNIQUE
    -- key below, so the same construct relaxed by both models is two rows, not a clash.
    fidelity             INTEGER NOT NULL,              -- 0 raw 1 FF 2 ML 3 xTB 4 DFT
    method_id            INTEGER REFERENCES methods(id),
    energy               REAL,
    converged            INTEGER,
    relaxed_from         INTEGER REFERENCES geometries(id) ON DELETE SET NULL,

    -- RESERVED M5: construction is a deterministic function of (choice_vector, seed),
    -- so a geometry can be regenerated rather than only replayed from coordinates.
    choice_vector_digest TEXT    DEFAULT NULL,
    choice_vector_json   TEXT    DEFAULT NULL,
    seed                 INTEGER DEFAULT NULL,

    qc_json              TEXT    DEFAULT NULL,          -- RESERVED M6 (clash/bond/inter-centre)
    created_at           TEXT    NOT NULL,
    UNIQUE (structure_id, coords_hash, method_id)
);

CREATE INDEX IF NOT EXISTS ix_geometries_structure ON geometries (structure_id, fidelity);
CREATE INDEX IF NOT EXISTS ix_geometries_energy    ON geometries (energy);
CREATE INDEX IF NOT EXISTS ix_geometries_choice    ON geometries (choice_vector_digest);

-- ── sites: perceive once, refresh accessibility (RESERVED, M4) ───────────────

-- A site is a place a BOND CAN FORM, and that is two things, not one: a donor atom on a
-- ligand, and a vacant coordination vertex on a metal.  `role` tells them apart.  Both are
-- frames (D13) -- origin, outward axis, reference direction -- which is what lets
-- `assembly.compatible(Site, Site)` ask one question of a donor/vacancy pair in M5,
-- instead of needing a second type and a second code path for the metal side.
--
-- `slot` exists because a metal carries SEVERAL vacancies on ONE atom, which the old
-- UNIQUE (structure_id, canonical_idx) could not hold.  A donor is always slot 0 and
-- therefore keeps exactly its previous uniqueness: one donor row per atom.
CREATE TABLE IF NOT EXISTS site_catalog (
    id              INTEGER PRIMARY KEY,
    structure_id    INTEGER NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    canonical_idx   INTEGER NOT NULL,          -- keyed to structures.canonical_order_json
    role            TEXT    NOT NULL DEFAULT 'donor',   -- donor | vacancy
    slot            INTEGER NOT NULL DEFAULT 0, -- 0 for a donor; 0..n-1 per metal vacancy
    donor_type      TEXT    NOT NULL,          -- '' for a vacancy: it is not a donor
    labile          INTEGER NOT NULL DEFAULT 0,
    charge_after    INTEGER NOT NULL DEFAULT 0,
    live_dof        TEXT,                      -- 'live' | 'free'      (D13)
    binding_modes   TEXT,                      -- comma-separated set  (D13)
    frame_json      TEXT,                      -- origin + axis + ref  (D13)
    algo_perception TEXT,
    UNIQUE (structure_id, canonical_idx, slot)
);

CREATE TABLE IF NOT EXISTS site_state (
    id                   INTEGER PRIMARY KEY,
    site_id              INTEGER NOT NULL REFERENCES site_catalog(id) ON DELETE CASCADE,
    geometry_id          INTEGER NOT NULL REFERENCES geometries(id) ON DELETE CASCADE,
    status               TEXT    NOT NULL,     -- open | occupied | blocked
    pka                  REAL,
    fukui                REAL,
    buried_vol           REAL,
    marginal_de          REAL,
    ease_scalar          REAL,
    ease_components_json TEXT,
    confidence           REAL,
    provisional          INTEGER DEFAULT 0,    -- in-pocket donor: promote to xTB (C5)
    fidelity             INTEGER,
    method_id            INTEGER REFERENCES methods(id),
    UNIQUE (site_id, geometry_id)
);

CREATE INDEX IF NOT EXISTS ix_site_catalog_type ON site_catalog (donor_type);
CREATE INDEX IF NOT EXISTS ix_site_state_status ON site_state (status, ease_scalar);

-- ── provenance and routes (edges LIVE from M3, scores RESERVED M8) ───────────
-- Identity lives on the node, sequence on the edges (D2).  One structure reachable
-- by two routes is one row here with two incoming edges.

CREATE TABLE IF NOT EXISTS reactions (
    id                   INTEGER PRIMARY KEY,
    product_structure_id INTEGER NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    kind                 TEXT    NOT NULL DEFAULT 'assembly',  -- ingest|assembly|reaction
    intermediate         INTEGER NOT NULL DEFAULT 0,
    atom_map_json        TEXT,                                 -- RESERVED M5
    conditions_json      TEXT,                                 -- RESERVED M8
    dG                   REAL,                                 -- RESERVED M8
    barrier_proxy        REAL,                                 -- RESERVED M8 (C6)
    barrier_proxy_json   TEXT,                                 -- RESERVED M8
    method_id            INTEGER REFERENCES methods(id),
    fidelity             INTEGER,
    choice_vector_digest TEXT,                                 -- RESERVED M5
    depth                INTEGER,                              -- steps from a reagent
    note                 TEXT NOT NULL DEFAULT '',
    created_at           TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS reaction_reagents (
    reaction_id  INTEGER NOT NULL REFERENCES reactions(id) ON DELETE CASCADE,
    structure_id INTEGER NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    stoich       INTEGER NOT NULL DEFAULT 1,
    role         TEXT    NOT NULL DEFAULT 'reagent',   -- reagent|chelator|solvent|leaving
    PRIMARY KEY (reaction_id, structure_id, role)
);

CREATE TABLE IF NOT EXISTS pathways (               -- RESERVED M8
    id         INTEGER PRIMARY KEY,
    name       TEXT NOT NULL,
    score_json TEXT,
    created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS pathway_steps (          -- RESERVED M8
    pathway_id  INTEGER NOT NULL REFERENCES pathways(id) ON DELETE CASCADE,
    step_ix     INTEGER NOT NULL,
    reaction_id INTEGER NOT NULL REFERENCES reactions(id) ON DELETE CASCADE,
    PRIMARY KEY (pathway_id, step_ix)
);

CREATE INDEX IF NOT EXISTS ix_reactions_product ON reactions (product_structure_id);
CREATE INDEX IF NOT EXISTS ix_reactions_depth   ON reactions (depth);
CREATE INDEX IF NOT EXISTS ix_reagents_struct   ON reaction_reagents (structure_id);

-- ── runs and tasks: work is queued, not called ───────────────────────────────
-- Ground rule 8's substrate.  A build is a stored spec that fans out into independently
-- claimable tasks, so the thing that EXECUTES is separable from the thing that ASKS.
-- One in-process worker on a laptop and N processes on the workstation use the same
-- table and the same claim; an MPI launcher or a GPU worker becomes another consumer
-- rather than a rewrite.  State lives here and not in any interface, which is why the
-- viewer can stay read-only and disposable while jobs run for hours.

CREATE TABLE IF NOT EXISTS runs (
    id          INTEGER PRIMARY KEY,
    spec_digest TEXT NOT NULL,
    spec_json   TEXT NOT NULL,
    status      TEXT NOT NULL DEFAULT 'pending',   -- pending|running|done|failed
    note        TEXT NOT NULL DEFAULT '',
    host        TEXT NOT NULL DEFAULT '',
    -- Which accelerator this run was submitted under.  Ground rule 6 applied to
    -- hardware: the device is DECLARED, and a stored run must still say which one
    -- produced it — otherwise "this took four hours" loses the only fact that explains
    -- it.  Empty on runs recorded before the column existed, which is honest: they
    -- never expressed a choice, and back-filling 'cpu' would invent provenance.
    device      TEXT NOT NULL DEFAULT '',
    workers     INTEGER          DEFAULT NULL,
    -- Why candidates were NOT queued.  A run that plans zero tasks is a legitimate
    -- outcome and must never be silent: "your filter matched nothing" and "it worked"
    -- have to look different from the outside.
    diagnostics_json TEXT NOT NULL DEFAULT '[]',
    created_at  TEXT NOT NULL,
    finished_at TEXT
);

CREATE TABLE IF NOT EXISTS tasks (
    id           INTEGER PRIMARY KEY,
    run_id       INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    kind         TEXT    NOT NULL,                 -- place | join (M5) | relax (M7)
    payload_json TEXT    NOT NULL,
    status       TEXT    NOT NULL DEFAULT 'pending',
                                        -- pending|claimed|done|failed|rejected
    priority     INTEGER NOT NULL DEFAULT 0,
    claimed_by   TEXT,
    claimed_at   TEXT,
    attempts     INTEGER NOT NULL DEFAULT 0,
    structure_id INTEGER REFERENCES structures(id) ON DELETE SET NULL,
    geometry_id  INTEGER REFERENCES geometries(id) ON DELETE SET NULL,
    error        TEXT,
    -- A stable machine-readable reason, so a run's outcomes can be GROUPED.  Free text
    -- cannot be: "2 clash(es), closest 1.40 A" and "...1.41 A" are one finding and two
    -- strings.  Set for rejected and failed tasks alike (`qc_clash`, `placer_refused`,
    -- `IndexError`, ...).
    error_code   TEXT,
    -- Everything the outcome was actually made of: the QC report with atoms, elements,
    -- measured distances and the limits they were measured against; the M-L distance
    -- each donor was placed at and whether that number was calibrated or estimated; the
    -- embed report (retry / force field); and whether the structure and geometry were
    -- newly written or recognised as already present.  A run's diagnosis has to be
    -- readable AFTER the console has scrolled away.
    detail_json  TEXT,
    -- Did this task write a new registry row, or recognise one that already existed?
    -- Idempotency on identity is the design (D2), which means a task doing nothing new
    -- is the SUCCESS case — and until now it was indistinguishable from one that built
    -- something.  NULL = not applicable / not recorded.
    structure_created INTEGER,
    geometry_created  INTEGER,
    created_at   TEXT NOT NULL,
    finished_at  TEXT
);

CREATE INDEX IF NOT EXISTS ix_tasks_claimable ON tasks (run_id, status, priority, id);
CREATE INDEX IF NOT EXISTS ix_tasks_status    ON tasks (status);
CREATE INDEX IF NOT EXISTS ix_runs_digest     ON runs (spec_digest);

-- ── descriptor tables (RESERVED, M1) ─────────────────────────────────────────
-- The shared substrate for C5/C6/C7.  Mostly lookup; every row carries its source.

CREATE TABLE IF NOT EXISTS donor_descriptors (
    donor_type        TEXT PRIMARY KEY,
    pka               REAL,
    pka_sigma         REAL,
    hsab              TEXT,          -- hard | borderline | soft
    default_denticity INTEGER,
    live_dof          TEXT,
    source            TEXT NOT NULL DEFAULT '',
    source_version    TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS metal_descriptors (
    symbol            TEXT    NOT NULL,
    charge            INTEGER NOT NULL,
    ionic_radius      REAL,
    hsab              TEXT,
    preferred_cn      TEXT,
    d_electrons       INTEGER,
    exchange_lability TEXT,          -- fast | moderate | slow
    source            TEXT NOT NULL DEFAULT '',
    source_version    TEXT NOT NULL DEFAULT '',
    PRIMARY KEY (symbol, charge)
);

-- Ligand fragments, so that "which structures contain this ligand" is an indexed join
-- rather than a string match on a label.  This is the searchable half of naming: queries
-- go through here, labels are composed afterwards from what comes back.
CREATE TABLE IF NOT EXISTS structure_fragments (
    structure_id   INTEGER NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    fragment_l1    TEXT    NOT NULL,        -- identity of the fragment on its own
    formula        TEXT    NOT NULL,
    count          INTEGER NOT NULL DEFAULT 1,
    bridge_class   TEXT    NOT NULL DEFAULT 'none',
    n_metals_bound INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (structure_id, fragment_l1)
);
CREATE INDEX IF NOT EXISTS ix_fragments_l1      ON structure_fragments (fragment_l1);
CREATE INDEX IF NOT EXISTS ix_fragments_formula ON structure_fragments (formula);

-- Human names for fragments.  DISPLAY ONLY: an alias may be absent, wrong or renamed and
-- no query changes, because nothing searches on it.  A fragment with no alias falls back
-- to its formula, which is honest rather than pretty.
CREATE TABLE IF NOT EXISTS fragment_aliases (
    fragment_l1 TEXT PRIMARY KEY,
    name        TEXT NOT NULL,
    source      TEXT NOT NULL DEFAULT '',
    created_at  TEXT NOT NULL
);

-- Child table for the metal vocabulary: makes "which metals exist" and "structures
-- containing Cu" index lookups instead of string splitting.  Maintained alongside
-- structures.metals, which stays for cheap single-row display.
CREATE TABLE IF NOT EXISTS structure_metals (
    structure_id    INTEGER NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    symbol          TEXT    NOT NULL,
    oxidation_state INTEGER,
    spin_class      TEXT,
    count           INTEGER NOT NULL DEFAULT 1,
    PRIMARY KEY (structure_id, symbol, oxidation_state, spin_class)
);
CREATE INDEX IF NOT EXISTS ix_structure_metals_symbol ON structure_metals (symbol);

-- ── free-form tagging (LIVE) ─────────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS structure_tags (
    structure_id INTEGER NOT NULL REFERENCES structures(id) ON DELETE CASCADE,
    tag          TEXT    NOT NULL,
    PRIMARY KEY (structure_id, tag)
);

-- ── the read surface consumers should use ────────────────────────────────────
-- The viewer reads this view, not the base tables, so the underlying layout can
-- change without breaking it.

CREATE VIEW IF NOT EXISTS v_structures AS
SELECT
    s.*,
    g.energy      AS best_energy,
    g.fidelity    AS best_geom_fidelity,
    g.converged   AS best_converged,
    g.coords_hash AS best_coords_hash,
    g.n_atoms     AS best_n_atoms,
    (SELECT COUNT(*) FROM geometries gg WHERE gg.structure_id = s.id)          AS n_geometries,
    (SELECT COUNT(*) FROM reactions  r  WHERE r.product_structure_id = s.id)   AS n_incoming_routes,
    (SELECT MIN(r.depth) FROM reactions r WHERE r.product_structure_id = s.id) AS min_depth,
    (SELECT GROUP_CONCAT(t.tag, ',') FROM structure_tags t WHERE t.structure_id = s.id) AS tags
FROM structures s
LEFT JOIN geometries g ON g.id = s.best_geometry_id;
