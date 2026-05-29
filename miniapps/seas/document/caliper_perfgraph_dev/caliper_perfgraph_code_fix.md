# Fix Report: Caliper perfgraph IMPLEMENTATION — round 2 (code review fixes)

Source review: `caliper_perfgraph_code_check.md` (R-001, R-002, R-003).

## Summary
- Findings addressed: **3 of 3**
- Files modified:
  - `dynamic/wave_operator.inl` (4 new `MFEM_PERF_SCOPE` — R-001 ×3, R-003 ×1)
  - `document/caliper_perfgraph_dev/README_caliper_build_and_run.md` (R-002)
- Tests added: n/a in the unit sense (annotations are no-ops on the default
  build); the runtime verification is the README's ADER-region grep (Caliper
  build only). Local regression gate re-run below.
- Test suite: **PASS** — force-clean rebuild from the worktree against the
  main-repo MFEM; 6/6 wave-op/ADER tests pass, `seas_spatial_dyn_driver` links,
  0 compiler errors (`/tmp/caliper_verify3.log`, `VERIFY_RC=0`).

## Changes Made
1. **R-001 (MODERATE) — ADER face-flux cost unprofiled for TPV31.** Added
   `MFEM_PERF_SCOPE` to the three ADER-path compute methods that `AdvanceADER`
   actually calls on the ADER path (`wave_operator.inl:5327-5332`):
   - `ComputeADERFaceFluxRHS` (`:3466`) → `"seas::WaveOperator::ComputeADERFaceFluxRHS"`
   - `ComputeADERSharedFaceFluxRHS` (`:4486`) → `"seas::WaveOperator::ComputeADERSharedFaceFluxRHS"`
   - `ComputeADERVolumeUpdate` (`:2260`) → `"seas::WaveOperator::ComputeADERVolumeUpdate"`
     (the review's "optional but recommended" item — added for a symmetric ADER
     volume/face breakdown; it wraps the already-annotated `ComputeVolumeRHS`).
   The pre-existing `Mult`/`ComputeFaceFluxRHS`/`ComputeSharedFaceFluxRHS`
   annotations were KEPT (correct for the RK path); they are simply inactive on
   the ADER/TPV31 path, which is now covered by the new regions.
2. **R-002 (LOW) — README nesting wrong.** Rewrote the region table (added the 4
   new regions; tagged each region `ADER` / `RK only` / `both`) and replaced the
   nesting paragraph with an accurate ADER call tree
   (`step → AdvanceADERWithSubStep → {ComputeADERSubStepStates(→ApplySpatialDerivative),
   friction_substep, AdvanceADER(→ComputeADERVolumeUpdate→ComputeVolumeRHS,
   ComputeADERFaceFluxRHS, ComputeADERSharedFaceFluxRHS)}`). Also fixed the
   stale **Sanity-check grep**, which looked for `seas::WaveOperator::Mult` — a
   region that never fires for TPV31 (RK-only), so it would have falsely reported
   "MISSING"; it now greps `seas::WaveOperator::ComputeADERFaceFluxRHS`. Updated
   the per-face-hooks note to reference `Compute{,ADER}FaceFluxRHS`.
3. **R-003 (LOW/POSSIBLE) — predictor depth.** Verified `ApplySpatialDerivative`
   (`:892`) is called inside a `(node × CK-iter × direction)` loop
   (`wave_operator.inl:1151, 1273`) — O(nodes·(order−1)·3) ≈ a handful of calls
   per sub-step, NOT per-element (the element loop is INSIDE the method). Safe
   per the granularity rule, so annotated it:
   `"seas::WaveOperator::ApplySpatialDerivative"`.

## Deviations from suggested fixes
- R-001: also annotated `ComputeADERVolumeUpdate` (the review marked it
  "optional but recommended"). Done for ADER symmetry; not a contradiction.
- R-003 was "POSSIBLE — annotate only if not per-element." Confirmed not
  per-element, so applied. No reviewer suggestion was wrong.

## Unresolved Findings
- None.

## New Tests
- The README runtime sanity check now asserts an ADER region
  (`seas::WaveOperator::ComputeADERFaceFluxRHS`) appears in `runtime-report` —
  the executable check that R-001 is fixed (Caliper build, Frontera).
- Local regression: the 6-test force-clean gate (re-run) is the guard that the
  4 new annotations are still no-ops on the default build.

## Annotations-only diff (post-fix)
- 16 insertions / 0 removals across the 3 source files; every added line is an
  `MFEM_PERF_*` macro (verified via `git diff | grep`).
- `wave_operator.inl`: 11, `bimaterial_wave_operator.inl`: 1,
  `spatial_dyn_driver.cpp`: 4.

## Ready for Re-Review: YES
TPV31 (ADER) perfgraph now resolves volume, interior-face, shared-face, friction,
and predictor cost as distinct regions; the default build is byte-unchanged
(6/6 tests pass, 0 compiler errors).
