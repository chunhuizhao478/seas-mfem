# P4a — fault half at np>1: impl + review — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed).  Commit 3543026.
Flow: scope → implement → review → fix → verify.

## Insight that made it small
D-2 (lts ⇒ fault-locality partition, P-017) keeps EVERY fault face
**rank-interior** at np>1.  So the fault half is purely local — the interior-only
per-cluster fault-QP eval (`EvaluateBulkAtFaultQPsCanonicalRange`, guarded on
`GetNumSharedFaultQPs()==0`, which HOLDS at np>1 under fault-locality) adds ZERO
collectives.  The only np>1 concern is the BULK half's rank seams, which reuse the
already-validated Stage-1/2 seam machinery.

## Implementation
- `LtsFaultSyncStepper`: `Predict` runs `PrepareSeamCoarseForecast` (bulk diff-1
  coarse-D(k) exchange, c≥1) under `SetExchangeProviderDk`; `Correct` forwards
  `cluster_id_for_seam` + the schedule to `AdvanceADERClusterBulk` (bulk rank
  seams); the ctor assert relaxed so a fault-FREE rank at np>1 is legal (fault
  half no-op, bulk-only).
- driver `lts_fault_stepping`: gate made PURELY rank-uniform
  (`LtsEnabled() && num_fault_global>0 && optin`); dropped `nprocs==1`;
  `exchange_bulk_provider_dk` + `SetExchangeProviderDk` = `(nprocs>1 && Nc>1)`;
  per-sync ghost-exchange counter + `SeamCoarseBuffersZero` in the fault loop.
  V_max is already an `MPI_Allreduce(MAX)`; the paraview/station output reuses the
  np>1-ready GTS machinery.

## Local gate (all GREEN)
| Gate | Result |
|---|---|
| **np=2 == np=1** fault interleave (TPV104-1000m, Nc=6, 2 syncs) | max\|Δfield\|=**1.67e-27**; 9/9 SCEC stations BYTE-IDENTICAL |
| **np=2 == np=1** fault interleave (8 syncs, ragged growth) | max\|Δfield\|=2.26e-26; stations byte-identical |
| clean completion (first-ever np>1 fault run) | no hang/abort; LTS-aware partition injected; matched-collective counter green |
| np=1 fault interleave | unchanged (new paths no-op) |
| regression | bulk seam 60/60 + 111/111 · predictor 31 · reorder 14 · scheduler 31 · layout 102 · tpv102-ADER 4 |

The differential is self-relative (np=1 vs np=2 of the SAME binary/deck), so it is
locally valid (the reorder-canary precedent), NOT a physics-vs-reference reproducer.

## Adversarial review — PASS-WITH-FIXES (no wrong-answer bug on the np=2 path)
Verified: `num_fault_global`/`num_shared_global` are real Allreduces; the fault
half is exchange-free (interior-only eval, confirmed by grep); the counter
accounts for bulk seams only; `exchange_bulk_provider_dk` is set identically on
`lts_meta` and the stepper; the iterator's `SetSubSteps` is unconditional (valid on
a fault-free rank); V_max is a matched Allreduce.
- **[MODERATE] FIXED** — the gate's per-rank clause `(FaultFacesReordered() ||
  num_fault_total==0)` could DIVERGE (a rank owning only SHARED fault faces →
  false while clean ranks → true → deadlock) if D-2 were ever violated.  Fixed:
  the gate is now purely rank-uniform; D-2 is enforced fail-loud + uniformly by
  `num_shared_global==0` (before the block) + the stepper ctor asserts (on all
  ranks).
- **[LOW] FIXED** — stale "np=1" wording in the ctor D-2 assert message.
- **[LOW] noted** — the fault-FREE-rank path (np≥4, a rank owns zero fault faces)
  is defensively handled (relaxed ctor + skip-fault-half) but UNTESTED locally
  (the np=2 TPV104 gate has fault on both ranks); verify at np≥4 on Frontera.

## np-scaling validation (local, np≤10 per the laptop-MPI rule) — GREEN
Fault-interleave differential np ∈ {2,4,10} vs np=1, TPV104-1000m, Nc=6, over 15
sync intervals (≥12-sync gate), `SEAS_LTS_FAULT_INTERLEAVE=1`:

| np | fault field vs np=1 | SCEC stations | shared fault QPs | run |
|---|---|---|---|---|
| 2 | max\|Δ\| ~3.9e-26 | 9/9 byte-identical | 0 (D-2) | clean, 15 syncs |
| 4 | max\|Δ\| ~3.6e-26 | 9/9 byte-identical | 0 | clean, 15 syncs |
| 10 | max\|Δ\| ~3.6e-25 | 9/9 byte-identical | 0 | clean, 15 syncs |

`num_fault_global` printed `shared = 0` at every np ⇒ the fault-locality partition
kept all fault faces rank-interior (D-2 held).  No hang / abort / matched-collective
violation at any np.  **np=10 exercised the fault-FREE-rank path** (far-field ranks
own zero fault faces) and completed cleanly — closing the review's one untested-path
LOW item.  (A pre-existing blanket "rate_state at np>1 … shared psi 1st-order"
WARNING fires on nprocs>1 regardless of the shared count; moot here since shared=0 —
a candidate one-line refinement is to gate it on `num_shared_global>0`.)

## Scope / remaining
- Active-rupture fidelity (T_nuc=1.0s ⇒ ~72 syncs, too slow locally) + np≥4
  fault-free-rank coverage + the SAFS/TPV104 physics acceptance stay Expanse/
  Frontera-staged (Phase 5), as does the performance measurement.
