# Design history — the changelog of `DESIGN_registry_assembly.md`

**Archive.** The active design doc carries the *current* reasoning: the locked decision
ledger, the subsystems as they now stand, and the checkpoints still open. This file carries
how it got there — every revision, newest first.

Go here when a decision looks arbitrary and you want the argument that produced it, or when
you are about to reverse one and need to know what it cost last time.

---


- *(M6/S0)* **D20, D24 and D25 — the polynuclear calls, made before the placer exists.**
  Three decisions M6 could not start without, and the first of them reverses the plan's own
  estimate of where the milestone's cost lives.

  * **D20 — the multi-centre placer is not the headline cost.** §6.1 had said it was, in the
    same breath as "done right the first time", so it had been carried unexamined since the
    schema was drawn. Measured against the battery, a Cu paddlewheel builds QC-clean and
    matches its M2 fixture at both L0 and L1 with **no solver anywhere** — the sp2 donor's
    second lone-pair well, `kabsch`, the curated distance table and the existing `qc` are
    the whole of it. What the battery *did* turn up is narrower and was not in the plan at
    all: there are **two bridge mechanisms** and only the multi-atom one was supported, and
    for oxo-centred clusters two independent determinants fix the same M···M and disagree by
    about half an ångström. So `place_multicentre` survives for a reason that can be
    measured (Fe₃ 0.595 Å, Zn₄O 0.491 Å) rather than for a reason that was assumed, and it
    is one scalar per edge, not a constrained optimisation. The superseded text is the last
    bullet of §6.1, struck there rather than deleted.

    The reframing has a cost and it is stated in the ledger rather than discovered later:
    the rigid-ligand model is **wrong by ~0.5 Å** on oxo-centred clusters. Real structures
    close that gap in the ligand — the O–C–O angle opens — and a rigid one cannot, so the
    residual is reported as strain and left for relaxation. Averaging the two determinants
    was specifically refused: it would store a number neither of them asked for, and it
    would make a known limitation invisible.

  * **D24 (C9) — multiplicity is stated.** Found by the two ground-truth fixtures
    disagreeing: `cu_paddlewheel` declares 1 where `combined_multiplicity` gives 3, and
    `fe3_mu3_oxo` declares 16 where it agrees. The tempting read is that one fixture is
    wrong. The right read is that **the rule has a precondition** — unpaired electrons add
    when the centres are independent — and two d⁹ Cu(II) 2.6 Å apart are not. Coupling is
    not recoverable from the centres, so it is a Kind-C declaration, and the multicentre
    path requires it the way `from_rdkit` already does. Both fixtures then come out correct.

  * **D25 (C10) — an M–M edge is declared.** Reached twice, independently, which is the
    strongest evidence a gate is real: from the code side, because `EdgeType.METAL_METAL` is
    in the certificate and a threshold applied silently is a guessed identity; and from the
    data side, by the parallel session that built `node_cases.tsv` and gave `mm_bond` its own
    column with the note that it is "a chemical decision, NOT derivable from d_mm". The Fe₃
    trimer has no edge at 3.29 Å and the paddlewheel has one at 2.62 Å, and the fixtures
    differ at L1 by exactly that.

  Two things about the battery worth keeping. Its targets are **typical of the named
  compound class, not a refinement of a deposited structure** — no CIF was consulted — which
  is why every row carries a window and the tests treat them as a shape check; narrowing one
  against the CSD is a separate, per-row job. And the reference table and the `@m6` gates
  arrived from a parallel session *before* any placer existed, which is what made D20
  callable on measurements instead of on the estimate it replaced.

- *(M6/S0)* **B12 closed: `metal_block` labelled a metal differently from every other
  producer.** Details in [`BUGS_resolved.md`](BUGS_resolved.md). The part worth recording
  here is what it says about coverage rather than about the bug: the fix moved **no test and
  no golden hash**, because every stored hash descends from `examples.py` and nothing had
  ever asserted on an identity that came out of the enumerator. D2's claim — one node, two
  routes — was being tested along one route. The fix therefore ships with a
  producer-agreement test pinning the label string itself, which is the assertion that would
  have caught it.

- *(M5 closeout)* **D21, D22, D23 — the three calls M5 made in code and never wrote down.**
  Recorded when M6 started, because M6 adds ledger entries and a ledger that is behind the
  code is worse than no ledger.

  * **D21 — θ_geom = 0.15 Å.** The calibration's first result is that it **could not be done
    where the plan said to do it**: at RAW the Kind-A spread is identically zero over 105
    pairs, because a join is a deterministic function of its choice vector and absorbs the
    ligand's embedding noise completely. There is no stochastic peak to sit above, and a
    threshold placed there would be calibrated against no noise. Relaxation re-introduces it
    and the valley is then 0.61 Å wide with nothing in it. The number is the geometric mean
    of the two bounds, and both bounds ship as named constants beside it so the threshold
    cannot drift away from the data that set it.

    **The calibration paid for itself by finding a defect instead of a number.** The first
    run put Kind A at a median of 1.20 Å — wider than most of Kind B, which would have made
    the exercise meaningless. One product pair had its pyridine ring carbons 2.3 Å apart
    while the coordinating N moved 0.15 Å: the ligand was bound in the right place and
    *rotated about the M–N axis*. `_place_donor_block` was matching axes and discarding the
    frames' `ref` vectors, so the roll was settled by `rotation_between`'s minimal rotation,
    which depends on how the ligand happened to be oriented in its own coordinate file.
    Re-embedding rotates a ligand rigidly, so one choice vector gave different rolls. That is
    precisely what `sites/frames.py` opens by saying a lone outward vector cannot do — the
    join was using the vector half of the frame it was handed. Frame onto frame through an
    orthonormal triad: Kind-A spread 1.56 Å → 0.0000 Å.

    C2's energy half stays open on purpose. The Kind-A *energy* spread reached 16.6 kcal/mol
    between samples whose cores agreed to 0.03 Å, all of it outside the core, because the
    rigid-core rule cuts a delocalised carboxylate C–O as rotatable (B9). Setting a window
    from that would bake a known defect into a stored threshold.

  * **D22 — where a tag is derived, and what is never backfilled.** The blocking version of
    this question dissolved on inspection: `put_structure` calls `l2_isomer_tag(g)` with no
    coordinates and nothing in `src/` passed `l2=`, so building the classifier split
    nothing and the corpus was untouched. What survives is smaller and real — the builder
    holds the coordinates, so it derives the tag, and `store_block`'s `l2=` override exists
    for the run pipeline, which passes `""` deliberately so its products land on the same
    node every other route reaches. Also folded in: the choice-vector digest is **not**
    backfilled either, and the contrast with L2 is the point — a digest is an annotation
    nothing points at, an identity is an address the provenance DAG points at. Re-derive
    annotations; version addresses.

  * **D23 — `MAX_BITE_MISMATCH_DEG = 40°`.** What makes this not folklore is that the two
    populations it separates are far apart and are *themselves* pinned by a test: the common
    chelators sit within 12° of an octahedral cis pair and every one of them misses a trans
    pair by 88° or more. Anything from ~35 to ~60 draws the same line, so the constant is
    reported as taken from the low half of a wide interval rather than as a measurement.

- *(M4 follow-up)* **Perception is resonance-invariant; the `site_catalog` seam below is
  closed.** The filed job was to stop perception reading bond order, and the shape of the
  fix is that **a delocalised oxo-acid is one donor, not several oxygens**:
  `sites.perception.DELOCALISED_GROUPS` matches the group's central atom, takes every
  terminal oxygen on it — ignoring metal neighbours, since coordination is not
  constitution — and gives all of them one donor type and one charge. No bond order is
  read anywhere in that path, which is what makes it invariant rather than merely
  patched. Three things worth recording:

  * The bug was **wider than the acetate case that surfaced it.** A sulfonate's two S=O
    oxygens were being typed `carbonyl_O` and only its anionic one `sulfonate_O`; nitro
    came out as one `carbonyl_O` and one `alkoxide_O`. Those are not near-misses, and they
    were route-dependent for the same reason acetate was.
  * **Charge is a group property here too**, exactly as in D15. All of a group's oxygens
    now report the same `charge_after` — anionic (−1) or not — instead of whichever one
    the resonance form parked the minus sign on. Per-site −1 is already what
    `LABILE_DONOR_PATTERNS` says for a diprotic acid, so this is the existing convention,
    not a new one.
  * The taxonomy is **deliberately coarse**: carbonate is `carboxylate_O`, sulfate is
    `sulfonate_O`, a phosphate diester is `phosphonate_O`. A name per oxo-acid is a
    promise to have anticipated every one of them — the same promise `Pocket` refuses to
    make. `nitro_O` is the one genuinely new type.

  `catalog_drift` stays as a guard rather than a known finding, and
  `tests/test_sites_state.py::test_no_build_route_drifts_from_the_stored_catalog` is the
  gate: the build routes in `tests/build_routes.py` reach one identity through every
  resonance form of it, and the assertion is that the registry cannot tell them apart.

  Still open, and NOT this job: the per-atom valence rules count a metal as an ordinary
  heavy neighbour, so a coordinated aqua oxygen is perceived as no donor at all. That is
  a `SiteStatus.OCCUPIED` row that never gets written, and it is a different fix with a
  different blast radius (`n_perceived_donors` moves for every assembled structure).

- *(M4 second half)* **C5 called, as D18.** The floor was ratified essentially as §6.6
  leaned, with two things the leaning did not say. First, the rule that turned out to
  matter most is not *what* the components are but that **an absent one stays absent**:
  the scalar renormalises over the components present, so a geometry-free record and a
  fully-populated one are distinguishable instead of both landing somewhere plausible.
  Second, **MACE-OMOL-0 is what made the call safe to make now.** Deprotonation is a
  charge change, so the ML rung was structurally unable to serve the model's primary
  component until a charge-aware potential existed — the floor was carrying the whole
  model with no affordable rung above it, which is a bad position from which to ratify a
  floor. §6.4's "cheap → rigorous" column for deprotonation should now read
  `pKa table → MACE-OMOL-0 ΔE → xTB → DFT`.

  Two findings from building it, both about seams rather than bugs:

  * **`site_catalog` is not a pure function of identity.** D15 excludes bond order from
    the L1 hash so C=O/C–O⁻ resonance forms hash identically; perception reads bond order.
    A monodentate acetate bound through either oxygen is one identity whose two build
    routes perceive different donor sets. The first catalog stands, `catalog_drift`
    reports the disagreement, and making perception resonance-invariant is filed as its
    own job — it needs a fixture set, not a patch inside a registry write. *(Done — see
    the M4 follow-up entry at the top of this changelog.)*
  * **`put_sites` deleting before inserting cascaded into `site_state`.** Under D2,
    re-deriving an identity the registry already has is the *expected* outcome for most of
    an enumeration, so a structure built twice kept state only on its second geometry and
    `n_open_sites` counted against a best geometry that no longer had any. The catalog is
    geometry-independent (§6.3); it is now written once per structure and kept.

- *(M7 implementation session)* Built the energy backends and the reference scheme. The
  session's finding is a correction to §11's first bullet: **the archived formation-energy
  equation is charge-balanced, and that is exactly why it went undetected.** A
  reaction-balance check — the fix §11 proposed — would have passed it. The rule that
  actually catches it is conservation of metal–donor bonds across the arrow, which
  distinguishes a ligand-exchange equation (errors cancel) from a formation-from-free-ions
  equation (they do not). See `docs/PLAN_implementation.md` rev 16. **D17 (proposed)**
  below follows from it.


- *(M2 implementation session)* Built the typed graph, canonicalisation and L0/L1 keys against
  hand-written polynuclear fixtures (paddlewheel, Fe₃-µ₃-oxo in two valence patterns, a
  bridging/chelating pair, cis/trans-Pt(NH₃)₂Cl₂, HS/LS hexaaqua). Raised **D14/D15/D16** as
  proposals in §7.1 — all three were forced by the code, and D16 in particular contradicts D3's
  choice of WL as the primary key. L2/L3 ship as stubs with final signatures. See
  `docs/PLAN_implementation.md` for the milestone context.

- *(this session, rev 1)* Initial spec: recursive BuildingBlock; identifier/address/provenance
  split; layered L0–L3 identity with typed-graph canonicalization; fidelity-laddered geometries;
  perceive-once site catalog + per-geometry state; named-component activation-ease model;
  polynuclear + reaction/pathway route-design layer; SQLite + content-addressed store. Ledger
  D1–D9 locked; checkpoints C1/C2/C4/C5/C6/C7 open.
- *(this session, rev 3)* Added §6.7 **construction-as-decision-tree** (Kind A/B/C branch points;
  conformer pairs are generated, not just detected). Revised **D11**: L3 identity is
  provenance-primary (choice-vector), geometry-verifier. Added **D13**: site = frame + live-DOF tag
  + binding-mode set (not a lone vector); join torsion stored as a discrete well index = conformer
  coordinate; `construct` emits its choice-vector; inference layer branches explicitly. Updated
  §3.1 `open_site`, §6.3 site model, roadmap steps 4–5.
- *(this session, rev 2)* Resolved **C1 → D10**: L2 in scope, discriminator is
  downstream-assembly-relevance + barrier (NOT energy gap); anthrarufin–Cu justifying case. Added
  **D11**: compound L3 trigger (rigid-core RMSD OR open-site-flag flip). Resolved **C4 → D12**:
  polynuclear-native from v1 (per-center labels; µ-carboxylate = one node + two dative edges;
  multi-center placer = headline cost). Added §6.6 **shared descriptor layer** with leaning
  proposals for C5 (heuristic pKa+HSAB floor), C6 (ΔG + sink + concurrent-bond-change +
  exchange-lability proxies), C7 (factorized HSAB-match, not a matrix). Added C2 worked examples,
  new **C8** (descriptor sourcing/provenance). Roadmap re-sequenced polynuclear-native.

---

- *(pre-M5)* **D19 — re-derive what is only an annotation; version what is an address.**

  Two questions arrived together and looked like one. Perception had a bug (a coordinated
  donor was perceived as no donor — see `BUGS_resolved.md`), so `ALGO_VERSIONS["perception"]`
  had to bump and 34 stored catalogs were known-wrong. Separately M5 is about to fill in
  `l2_isomer_tag`, which would split every identity written while L2 was a stub.

  Both are "a recipe changed; what happens to rows written under the old one?" — and the
  answers are opposite, which is the useful part:

  * **`site_catalog` is re-derived.** It is a derived annotation that nothing points at, so
    rewriting it costs only the `site_state` rows beneath it, and those were computed against
    a catalog now known to be incomplete. Absent beats stale (ground rule 9). Done lazily —
    a `perception/1` catalog is replaced the next time anything touches its structure — so
    there is no bulk migration and no window where the registry silently claims the old rows
    are current.
  * **An identity is versioned, never re-derived.** It is pointed at by the provenance DAG,
    and provenance is the one thing in this registry that is not regenerable (D2). Backfilling
    L2 would rewrite stored identities *and* every `reactions` edge pointing at them.

  What makes the second one honest rather than a silent fork is that `structures.algo_l2`
  already stamps the generation on every row, so `''` reads as "this predates L2" rather than
  "this has no isomerism". The accepted cost, stated rather than discovered later: the same
  species built before and after M5 can occupy two rows. A visible duplicate with its cause
  recorded is strictly better than an invisible one.

  Also fixed while pinning this down: `put_sites` took its version as a hardcoded
  `"perception/1"` literal, so the recipe was unpinned in violation of ground rule 6 and there
  was nothing for a staleness check to compare against. It now reads `ALGO_VERSIONS`.
