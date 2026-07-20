# P3-4 — Fault-half interleave (np=1): implementation + review record — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed). Flow: implement →
review → fix → verify locally.

## What P3-4 is
Each fault-bearing cluster advances its fault QPs at its own rate (dt_base·2^c)
instead of every fault QP every substep. The primitives (reorder P-006, range
friction `Advance` A.7, absolute range nucleation D-3) landed earlier; P3-4 is the
driver assembly + the two wave-operator pieces that were missing.

## Stages (all implemented)

**Stage A — extraction (behavior-preserving, user-chosen approach).**
The ~918-line per-face body of `WaveOperator::ComputeADERFaceFluxRHS` (fault /
boundary / interior flux, per QP) was extracted as a PURE MOVE into a private
helper `ProcessADERFaceToRHS_(int f, const real_t* I_data, real_t dt, Vector& rhs,
FaultEvalStageAvgMode eval_avg_mode)`. `ComputeADERFaceFluxRHS` now loops calling
it (after the interior cached fast-path); the 3 face-level `continue`s became
`return`s; `enum FaultEvalStageAvgMode` was hoisted to a member. Scripted move
(`extract_face_helper.py`), asserted anchors, brace-balanced.
- **Byte-exact gate PASSED:** a GTS TPV104 fault run (fault VTU c0,c68) is
  BYTE-IDENTICAL pre vs post the extraction.

**Stage B — corrector fault-face flux.**
`AdvanceADERClusterBulk` gained optional `const int* fault_face_ids, int
n_fault_faces`. When non-null, after the non-fault face sweep and before the mass
inverse it applies the cluster's fault-face flux via `ProcessADERFaceToRHS_`
(consuming the imposed states the friction range Advance wrote). nullptr = a
perfect no-op (Phase-2 bulk unchanged). This is the "reuse the EXACT GTS fault
flux" payoff of Stage A — no divergence between the LTS and GTS fault flux.

**Stage C — fault-aware stepper `dynamic/lts_fault_stepper.hpp`.**
`LtsFaultSyncStepper` (ILtsClusterStepper) mirrors the fault-free bulk stepper +
the fault half. Per due cluster c: Predict(c) rescales the iterator sub-steps to
dt_step, computes tau_nodes = scaled sub-step midpoints, runs the cluster
predictor, evaluates the bulk traction at cluster c's fault QPs, and runs the
friction RANGE `Advance` over the cluster's cluster-contiguous QP range with an
ABSOLUTE range nuc callback; the imposed states land in a persistent global
buffer. Correct(c) runs the bulk+fault corrector. Per-cluster QP range is derived
from the reordered `GetFaultInteriorFaces()` + `clusters[c].fault_faces`, with a
contiguity assertion.

**Stage D — driver wiring (`drivers/spatial_dyn_driver.cpp`).**
`lts_fault_stepping = LtsEnabled() && num_fault_total>0 && nprocs==1 &&
FaultFacesReordered() && SEAS_LTS_FAULT_INTERLEAVE-set`. OPT-IN (default OFF ⇒ the
still-GTS reorder canary + unchanged default science, per D-5). New sync-interval
block mirrors the bulk block + per-sync fault/station output + slip_rate reset +
V_max reduction + NaN tripwire + V2 checkpoint (with canonical perm). GTS loop and
final-checkpoint guards updated to exclude `lts_fault_stepping`.

## Local validation (all GREEN)
| Gate | Result |
|---|---|
| Stage-A extraction byte-exact (GTS fault VTU pre/post) | c0, c68 **byte-identical** |
| byte-exact-off regression (GTS run, current binary vs pre-Stage-BCD) | c0, c68 **byte-identical** |
| **Nc=1 interleave == GTS (correctness gate)** | max\|Δfield\| = **8.08e-27**; 9/9 stations byte-identical; VERDICT PASS |
| Nc=6 multi-cluster interleave smoke | runs clean, `done: 5 syncs, t=0.06s, V_max=1.02e-16`, no NaN / no contiguity or buffer errors |
| test-wave-operator | 25/0 |
| test-ader-tpv102-smoke (fault path byte gate) | 4/0 |
| test-lts-predictor | 31/31 |
| test-lts-fault-reorder / friction-range / nucleation-absolute | 14/0, 10/0, 62/0 |
| test-spatial-nucleation / lts-partition / lts-clustering / lts-layout | 73/73, 19/0, 4043/4043, 102/102 |

**Robustness hardening (proactive):** `LtsFaultSyncStepper::Predict` captures the
ORIGINAL configured friction sub-steps once in the ctor (`cfg_dT0_`/`cfg_w0_`) and
rescales from those each cluster, rather than reading the just-mutated shared
iterator — order-independent + obviously correct (was subtly-correct via
proportional-rescale invariance).

**Bug caught+fixed by the local smoke (implement→verify→fix worked):** the friction
range `Advance` asserts `Σ deltaT == dt_step`; the stepper must rescale the shared
iterator's sub-steps to each cluster's dt_step (and use the scaled sub-step
midpoints as the predictor tau_nodes), exactly as the GTS macro-step does. Fixed in
`LtsFaultSyncStepper::Predict`.

## Scope / what stays Frontera-staged
Multi-cluster PHYSICS fidelity (TPV104 station arrivals within 1%, SAFS breakout
within 2%) is production-mesh — Frontera, per the plan. The local gates prove the
plumbing + the Nc=1 reduction to GTS; the interleave is OPT-IN until that Frontera
gate + the D-5 default flip.

## Adversarial review findings
(Three parallel reviewers: stepper / corrector+extraction / driver wiring.)

### Reviewer B — corrector + extraction (wave_operator.inl/.hpp)
Extraction confirmed a FAITHFUL pure move; corrector logic confirmed correct
(fault flux before mass inverse; role-4 skip ⇒ no double-processing; I_cluster
valid by D-1; nullptr = byte-exact no-op; no np=1 MPI hazard).
- **R1 [MED] FIXED** — corrector now `MFEM_VERIFY`s the imposed-state buffer is set
  before the fault loop (else the inline EvaluateADER re-solve would double-advance
  the fault state).
- **R3 [MED] FIXED (corrected after re-verify)** — the first attempt asserted each
  cluster fault face is tagged role `Fault` *in the sweep* (face_ids); the re-verify
  ABORTED on it because `BuildLtsLayout` tracks fault faces ONLY in `cl.fault_faces`
  and never adds them to `cl.face_ids`.  The assertion was corrected to the actual
  invariant: no cluster fault face may appear in the bulk sweep `face_ids` (so the
  fault loop is its sole processor).  [Caught by the review→fix→verify loop.]
- **R4 [LOW] FIXED** — removed the dead caller-scope `bulk_bg_scaled` in
  ComputeADERFaceFluxRHS (the helper recomputes it).
- **R6 [LOW/Phase-4] FIXED** — corrector aborts loud if a listed fault face is a
  partition seam (shared) instead of the helper silently dropping it.
- **R7 [NIT] FIXED** — docstring `AdvanceADERCluster` → `AdvanceADERClusterBulk`.
- **R2 [MED] — test coverage:** the corrector fault loop has no dedicated UNIT
  test (the byte gates are fault-free). COVERED by the Nc=1 interleave==GTS
  INTEGRATION gate (interleave-with-fault reproduces GTS to 8e-27 on the real
  mesh). A small 2-element fault-fixture unit test is a noted follow-up.
- **R5 [LOW] — has_bulk_bg_ consensus:** the R-1505 fail-loud/collective guard in
  GTS ComputeADERFaceFluxRHS is not run on the LTS corrector path. PRE-EXISTING
  (shared with the fault-free bulk stepper); the LTS driver sets the background on
  all ranks before the first Correct, and LTS is np=1. Noted; not P3-4-specific.

### Reviewer C — driver wiring (spatial_dyn_driver.cpp)
Wiring confirmed a faithful, correctly-wired mirror of the bulk LTS block: all
`lts_stepping` uses updated for the fault mode; stepper ctor args match; checkpoint
`qp_canon_ptr` correct + symmetric on restart; station/paraview pairing right;
V_max valid; opt-in safe against half-activation; LSW+LTS ok; dt_base correct.
- **C1 [LOW] FIXED** — the opt-in now parses the VALUE (`0`/`false`/empty ⇒ OFF),
  not mere presence, since it changes the science.
- **C2 [MED] FIXED** — seed `sync = step0` (both LTS blocks) so output cycles + logs
  continue across `--restart` instead of colliding pre-restart ParaView frames.
- **C3 [MED, staged] FIXED** — rank-0 WARNING at Nc>1 that fault/station output is
  coarse (per-sync) cadence (sub-sync sampling is Phase 4; Nc=1 local gate unaffected).
- **C4 [LOW]** — NaN handling aborts (MFEM_VERIFY) vs the GTS loop's graceful exit;
  CONSISTENT with the bulk LTS block, left as-is (noted).
- **C5 [LOW/verify]** — `interior_flux="matrix"` (bimaterial) + lts + fault is
  reachable and NOT locally tested (TPV104 is scalar).  The cluster methods DO route
  through the bimaterial virtuals (Phase-2 validated bimaterial single-cluster==GTS;
  the extracted fault flux uses the same InteriorFaceFlux_/FluxForElem_ virtuals as
  GTS), so no restrictive guard added — flagged as a Frontera verification item for
  the fault+matrix combination.

### Reviewer A — stepper (lts_fault_stepper.hpp)
Deepest review — VERIFIED CORRECT all the load-bearing areas: tick schedule /
t_step_start, imposed-state non-staleness across ticks (disjoint per-cluster ranges),
QP↔offset coupling, buffer lifecycle, the sub-step rescale (now from captured
originals), contiguity + role-tag consistency, D-3 nucleation, and the
slip_rate_substep_max contract.  No np=1 wrong-answer defect.
- **A2 [MED] FIXED** — ctor now asserts `GetNumSharedFaultQPs()==0` (D-2, np=1).
- **A3 [LOW] FIXED** — added `<cmath>` + `<algorithm>` (were transitive-only).
- **A4 [LOW] FIXED** — ctor asserts `nbf>0 && n_total_fault_qps_>0` (fault layout built).
- **A5 [LOW] FIXED** — `Qpw_*` sized to exactly O (RunSubSteps_ requires `size()==O`).
- **A1 [MED] — unrestricted `EvaluateBulkAtFaultQPsCanonical` per cluster — DEFERRED
  (documented follow-up), NOT an np=1 correctness bug.** Predict evaluates ALL fault
  QPs from the cluster's predictor state (only the cluster's range is consumed).
  Consequences are all Frontera-scope: (a) ~2× GTS fault-trace cost ⇒ the fault-half
  speedup is not yet realized (perf is a Frontera measurement); (b) at Nc>1 it reads
  non-cluster (uninitialized/stale) predictor blocks — harmless in practice (unused;
  the Nc=6 smoke ran clean) but a UBSan/FP-trap nicety; (c) at np>1 the per-cluster
  ExchangeFaceNbrData violates the matched-collective contract (P-007) — a Phase-4
  blocker. **At Nc=1 (the local gate) cluster 0 owns every element, so there is NO
  uninitialized read and NO perf gap** — the interleave is correct + validated.  The
  proper fix is a face-subset / no-exchange `EvaluateBulkAtFaultQPs` overload; it is
  the TOP P3-4 follow-up, aligned with the plan's Frontera staging of multi-cluster
  perf + the Phase-4 np>1 work.

## Follow-ups
1. **Restricted per-cluster fault-QP eval** (Reviewer A1) — ✅ **IMPLEMENTED** (2026-07-19):
   extracted the per-interior-fault-face eval body into `EvalBulkAtFaultQPsForInteriorFace_`
   (PURE MOVE) and added `EvaluateBulkAtFaultQPsCanonicalRange(Q, fi_begin, fi_end, ...)`
   (interior-only, no-exchange); the stepper's Predict now evaluates only the cluster's
   contiguous face block `[qp_begin/nbf, qp_end/nbf)`.  Removes the ~2× whole-mesh cost,
   the Nc>1 uninitialized read (write range == read range), and the per-cluster collective
   (unblocks Phase 4).  Verified: extraction byte-exact (GTS fault VTU pre/post) + Nc=1==GTS
   preserved + Nc=6 clean.
2. **Dedicated corrector-fault-flux unit test** (Reviewer B2) — small 2-element fault
   fixture, `AdvanceADERClusterBulk(...,fault_face_ids)` vs GTS `AdvanceADER`.
3. **has_bulk_bg_ consensus on the LTS corrector path** (Reviewer B5) — pre-existing,
   shared with the bulk stepper.
4. **fault + interior_flux="matrix" (bimaterial)** (Reviewer C5) — reachable + should
   work via the bimaterial virtuals (Phase-2 validated bulk); Frontera verification item.
