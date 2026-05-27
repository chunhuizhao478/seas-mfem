# Fix Report — [DIAG-SIGN] speckle instrumentation review (2026-05-26)

Companion to `REVIEW.md` (speckle-instrument review). Addresses R-001…R-006.

## Summary
- Findings addressed: 5 of 6 (R-005 verified as a false finding — see below)
- Files modified:
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` (R-001, R-003, R-004, R-006)
  - `miniapps/seas/dynamic/fault_face_flux.hpp` (R-004 + R-003 doc)
  - `miniapps/seas/jobs/safs/spatial_dyn_slipweakening_nocap_triq_8N_400r_dev_2hr_safs.sbatch` (R-002)
- Tests: 164/164 tpv104-checkpoint (bit-exact round-trip), 48/48 phaseh
  ComputeMaxDt parity, 18/18 tpv102 nuc-callback parity — byte-exactness +
  serialization preserved (new field is transient/non-serialized).
- Driver rebuild: clean (relinked, no errors/warnings).

## Changes Made
1. **R-001 (multi-spot localization)** → replaced the single global-`MINLOC`
   `[DIAG-SIGN-DOF]` dump with a per-rank loop that prints EVERY local DOF whose
   interval sub-step σ_n is below the floor (or tensile when cap off), capped at
   32/rank. All concurrent speckle spots are now localized, not just the worst.
2. **R-002 (no opening/bulk cause)** → nocap sbatch now sets
   `SEAS_DIAG_SLIP=1`, `SEAS_DIAG_SLIP_VTHR=15.0` so the existing per-sub-step
   `sn_vjump`/`sn_sterm` decomposition fires on a speckle DOF (above the ~13 m/s
   front) — captures WHY (predictor/ghost opening vs tensile bulk), not just
   WHERE.
3. **R-003 (pre-onset sampling gap)** → `sigma_n_substep_min` is no longer reset
   every macro-step; it now ACCUMULATES the most-tensile sub-step over the whole
   diag interval and is reset inside the diag block only AFTER it is reported
   (globally consistent gate → no divergent collective). A tensile transient on
   a step the diag does not sample is captured at the next print. The
   `[DIAG-SIGN]` SUB-STEP line is relabelled `SUB-STEP(interval)`.
4. **R-004 (magic sentinel)** → replaced `1.0e300`/`1.0e299` with
   `std::numeric_limits<real_t>::max()` (local `SN_UNSET` in the driver; field
   default in the header). `<limits>` already included in both TUs.
5. **R-006 (log bloat)** → the `[DIAG-SIGN]` aggregate line and per-DOF dump now
   print only when `step%100==0 || n_ss_tensile_g>0`, so the block does not emit
   an all-compressive line every step once V>10.

## Unresolved Findings
- **R-005 (shared-QP pre-average σ_n)** — NOT a bug. Verified that
  `fault_face_flux.cpp:182` (`sigma_n_total = σ_n0+σ_n_nuc+sigma_n_trial`) lives
  in `CompleteFromTrial`, which runs AFTER the shared-QP trial averaging in
  `wave_operator.inl`. So `states[qq].sigma_n_total` is already recomputed from
  the averaged trial and is consistent with the end-of-step `data.sigma_n_corr`
  (`:322`). The suggested "fix" would be redundant (and wrong if a non-Trial
  avg mode is ever used). Left unchanged per fix-agent rule 5.

## New Tests
- None added in-tree: the instrumentation is env-gated diagnostic output in an
  MPI driver with no unit-test harness for the diag block. Validation is the
  byte-exact/serialization regression suite above (confirms no physics/format
  change) + the runtime checks listed in REVIEW.md "Unreviewed Areas", to be
  done on the first nocap+diag Frontera run.

## Ready for Re-Review: YES
