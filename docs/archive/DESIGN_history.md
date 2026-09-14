# Design history — the changelog of `DESIGN_registry_assembly.md`

**Archive.** The active design doc carries the *current* reasoning: the locked decision
ledger, the subsystems as they now stand, and the checkpoints still open. This file carries
how it got there — every revision, newest first.

Go here when a decision looks arbitrary and you want the argument that produced it, or when
you are about to reverse one and need to know what it cost last time.

---


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
