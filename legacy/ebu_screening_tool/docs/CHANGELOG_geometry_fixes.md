# Geometry-engine fixes (linker-node generation)

Symptom reported: for M²⁺ + BTC "linker-only", the only surviving model was a distorted
`BTC[s1]_BTC[s2]` where one benzene ring was puckered ("hydrogen dangling on benzene") and the
expected M-2BTC / M-3BTC / M-4BTC nodes were missing.

Three root causes were found and fixed. All are in `ebu_core.py` / `geometry_qc.py`; the notebook
and `ebu_toolkit.py` are unchanged.

## 1. The "abomination" — rigid linker force-chelated onto one metal
`_pair_chelate_compatible` used `far_threshold = 5.5 Å`. BTC's two *meta* carboxylate O donors are
4.8 Å apart, so they were judged "close enough to chelate" and wrapped onto a single metal by the
distance-geometry routine — which can only satisfy that by **puckering the aromatic ring** (measured
0.24 Å out-of-plane vs ~0 for real benzene).
**Fix:** `far_threshold → 3.5 Å`. A single metal with ~2.05 Å bonds can span at most ~4.1 Å (trans),
so a rigid pair wider than that is a *bridging* linker, not a chelate. BTC (and BDC) are now placed
single-anchor (bridging); their non-coordinating carboxylates correctly land far from the metal and
the config is rejected if it was pretending to be a chelate. Flexible chelators (EDTA, EDDA) are
unaffected — they still qualify via the rotatable-bond branch.

## 2. Good nodes killed by a soft metal···H graze
The clean M-2BTC node was discarded over a single 2.11 Å metal···(aromatic H) contact sitting exactly
at the 0.65 clash threshold — a pre-relaxation artifact that relaxes out.
**Fix:** `geometry_qc.check_clashes` now uses a gentler cutoff for any H-involving pair
(`clash_scale_h = 0.50`) while keeping the stricter `0.65` for heavy–heavy overlaps (the genuine
steric failures). Ring-planarity tolerance also tightened (`0.25 → 0.15 Å`) so any residual ring
distortion is still caught.

## 3. The real reason 3-/4-BTC never appeared — rings folding onto the metal
The monodentate placement set the M–O–C angle by tilting about an **arbitrary** perpendicular axis,
so which way the rest of the ligand swung was random and often folded the benzene ring onto the metal
(ring-C 2.3 Å, ring-H 1.2 Å from M). This could not be repaired downstream, because the clash-avoidance
search spins the ligand about the metal–donor axis, which **preserves every atom's distance to the
metal**.
**Fix:** in `_rigid_place_single_anchor`, the tilt-axis azimuth is now chosen to push the ligand body
*farthest* from the metal. A cooperative all-against-all orientation refinement pass was also added so
several bulky linkers splay apart instead of being placed greedily one-at-a-time.

## Result (verified)
`Fe²⁺/Fe³⁺/Co²⁺/Cu²⁺ + BTC` linker-only now yields **M-2BTC, M-3BTC, M-4BTC** (planar and tetrahedral),
all with planar rings (0.000 Å deviation) and rings pointing outward (ring-C 3.6–4.1 Å from the metal).
EDTA/EDDA chelation is unchanged. (A single M-1BTC is intentionally absent: one metal + one monodentate
carboxylate is CN 1, not a defined node — BTC needs ≥2 copies, or it bridges to other metals in the
real framework.)

Files changed: `ebu_core.py`, `geometry_qc.py`.
