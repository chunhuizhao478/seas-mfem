# P4a Stage 2 — cross-cluster (diff-1) rank seams: impl + review — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed).  Design:
`PLAN_p4a_stage2_diff1_2026-07-19.md`.  Builds on Stage 1 + 2a.

## What Stage 2 delivers
The fault-free BULK multi-rate stepper now handles CROSS-CLUSTER (diff-1) rank
seams at np>1 — completing Phase 4a for the bulk path.  A diff-1 seam joins a
FINE local element (cluster c) to a COARSE element (cluster c+1) on another rank;
the coarse in-tray is remote, so each rank fills its OWN element's contribution
from exchanged ghost data + the closed-form `[a,b]` (rank-identical), with no
cross-rank buffer messages.

## Implementation (commits aa48f41, 62778bf, ffbe2d1)
- **2a — D(k) ghost exchange** (`ExchangeClusterProviderDkGhost`): packs a
  cluster's provider Taylor stacks D(k) (raw, dt-independent) and exchanges them
  (scalar `ghost_gf_`, NUM_STATE·order per-component), caching `ghost_dk_[c]`.
- **2c — corrector** (`ComputeADERClusterSeamFaceFluxRHS`, 3 modes): diff-0
  (Stage 1); fine-side (local finer) integrates the GHOST coarse forecast over
  `[a,b]` → its own rhs; coarse-side (local coarser) integrates its LOCAL forecast
  over `[a,b]` against the ghost fine `I` → a LOCAL seam in-tray
  (`seam_coarse_buf_`), consumed at its own correct.  Both feed
  `SharedInteriorFaceFlux_` the SAME (fine I, coarse forecast) ⇒ bit-identical
  flux ⇒ conservation exact.
- **seam-coarse D(k) retention** (`EnsureSeamCoarseElems_` /
  `PrepareSeamCoarseForecast`): a coarse element at a rank seam is not a layout
  provider, so an extra element-local CK retains its raw D(k) operator-side
  (exact — D(k) is dt-independent).  Avoids a driver/layout ripple.
- **matched collectives**: the D(k) exchange fires for every predicting cluster
  c≥1 (rank-uniform gate) behind `meta.exchange_bulk_provider_dk`;
  `n_collectives += NUM_STATE·ader_order` per such predict.  Per-cluster caches
  sized to the clusters a rank TOUCHES (local + ghost), and an untouched cluster
  still exchanges (all-(-1) map) so no peer hangs.
- **`BuildLtsLayout`**: contiguity assert relaxed — a rank may own a subset of the
  global clusters at np>1.

## Local gates (all GREEN)
| Gate | Result |
|---|---|
| **diff-1 np=2 == np=1** (Nc=2; np=1 uses the VALIDATED Phase-2 LOCAL diff-1 buffer path) | 74/74; Δenergy=1.3e-16, Δcentroid=**0.0**, Δmax=**0.0** |
| **diff-1 np=3 == np=1** (Nc=3 CHAIN — middle cluster provider+consumer across two seams; B.11) | 111/111; Δenergy=3e-16, Δcentroid=**0.0**, Δmax=**0.0** |
| gate sensitivity (drop diff-1) | 99.8% energy divergence — the gate bites |
| per-sync ghost-exchange counter (P-007) | matched (incl. the D(k) term), all ranks |
| SeamCoarseBuffersZero + local BuffersZero at sync | green |
| Stage-1 diff-0 seam np=2 | 60/60 (unchanged) |
| LTS unit regression | predictor 31 · reorder 14 · scheduler 31 · layout 102 · tpv102-ADER 4 |

The Nc=3 np=3 run also exercises the **untouched-cluster matched-exchange** path
(rank 0 exchanges for cluster 2 with an empty map) — green.

## Scope / notes
- 4a is "simple exchange, correctness first."  Known 4a inefficiencies (4b/Phase-5
  optimize): the D(k) exchange is per-component (NUM_STATE·order) not batched; it
  fires for every c≥1 predict (not only true cross-rank providers); and
  seam-coarse D(k) is a recompute (extra CK) rather than folded into the main
  predict.  All are perf, not correctness.
- Fault half at np>1 is still separate (the fault stepper keeps its `nprocs==1`
  guard); the bulk seam path does not touch it.
- Production-mesh physics fidelity (TPV104 arrivals / SAFS breakout) + the
  performance measurement remain Expanse-staged (Phase 5).

## Adversarial review — PASS-WITH-FIXES (no CRITICAL/MODERATE live bug)
A reviewer audited the Stage-2c diff-1 numerics (the 3-mode corrector, the D(k)
exchange/retention, the step-2d consume, the stepper wiring, `BuildTickTable`,
and `lts_layout.cpp`), compared the seam `[a,b]` against the validated local
consumer, and traced the tick schedule + every `AdvanceADERClusterBulk` caller.
Verified correct: matched collectives at scale (incl. untouched-cluster ranks),
`[a,b]` identical to the local consumer, the fault path unaffected
(`lts_fault_stepper` omits `cluster_id_for_seam` ⇒ np=1 byte-exact), and the
index/reassembly/slot maps.  Four LOW hardening items — ALL FIXED:
- **[LOW] stale-D(k) guard (GAP-A3) absent** — `ghost_dk_epoch_` was written, never
  read.  FIXED: the mode-1 reader asserts `ghost_dk_epoch_[c_coarse] ==
  floor(tick/2^c_coarse)*2^c_coarse` (the coarse's current-step predict tick).
- **[LOW] maxdiff>1 at a seam silently skipped** — FIXED: the `else` branch now
  `MFEM_VERIFY(|c1-c_nbr|<=1)` before `continue` (fail loud on a clustering bug).
- **[LOW] ragged final interval (T_actual<T_s) untested at np>1** — FIXED: the
  diff-1 gate now makes its LAST sync ragged (`T_s - dt_base`); Nc=2 + Nc=3 still
  np>1==np=1 (d<7e-16).
- **[LOW] seam in-tray fill-count invariant not asserted** — FIXED: the step-2d
  consume asserts `fill ∈ [1,2]` (maxdiff-1 ⇒ ≤2 fine sub-steps; ragged ⇒ 1).

Post-fix gates GREEN: diff-1 Nc=2 74/74 + Nc=3 chain 111/111 (ragged final,
d<7e-16); diff-0 seam 60/60; predictor 31 · reorder 14 · scheduler 31 · layout
102 · tpv102-ADER 4.
