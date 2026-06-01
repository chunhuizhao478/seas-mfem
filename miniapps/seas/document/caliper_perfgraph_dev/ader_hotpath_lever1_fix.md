# Fix Report: Lever-1 review R-001 + R-002

## Summary
- Findings addressed: **2 of 2** requested (R-001, R-002). R-003/R-004/R-005
  deferred per instruction.
- Files modified: `tests/unit/test_wave_operator_spatial_derivative.cpp`,
  `drivers/spatial_dyn_driver.cpp`.
- Tests added: 4 assertions (scalar + bimaterial recursion equivalence × orders {2,3}).
- Test suite: **PASS** — `seas_test_wave_operator_spatial_derivative` 38/38
  (was 34); regression set all green incl. `bimaterial_wave_operator_parity` 9/9
  (default path still byte-exact). Driver builds + links; `--deriv-cache` present
  in the binary.

## Changes Made
1. **R-001 (recursion-level equivalence)** →
   `tests/unit/test_wave_operator_spatial_derivative.cpp`:
   - New `RecursionCachedVsOnTheFlyMaxRelErr(op, N, order, seed)` — sets the
     absorbing background (AdvanceADER's `has_bulk_bg_` precondition), then runs
     **`AdvanceADER`** and **`ComputeADERSubStepStates`** in both DerivModes and
     returns the max relative error.
   - Called in `TestCachedEquivalence` for the scalar `WaveOperator` AND the
     `BimaterialWaveOperator` (the production matrix path), orders {2,3}.
   - Result: max rel **8.9e-16 / 7.9e-16** (order 2), **1.4e-15 / 1.5e-15**
     (order 3) — the recursion is even tighter than the single call (dominated by
     the identical `D(0)=Q` term). Gate tol 1e-12.
2. **R-002 (production activation path)** → `drivers/spatial_dyn_driver.cpp`
   (right after the operator is bound, `WaveOperator<ParMesh> &wave = *wave_ptr;`):
   - `const bool use_deriv_cache = HasFlag(argc, argv, "--deriv-cache");`
   - `if (use_deriv_cache) { wave.SetDerivMode(DerivMode::Cached); /*root print*/ }`
   - Works for the scalar and bimaterial operators (cache is geometry-only).
     **Default (flag absent) stays `OnTheFly` → byte-identical**; the flag opts
     into the machine-eps cached path (R-002). R-004 budget guard still fires
     inside `SetDerivMode` if the cache is too large.

## Verification
- [x] **R-001:** recursion (AdvanceADER + ComputeADERSubStepStates) Cached vs
      OnTheFly ≤ 1.5e-15 (scalar + bimaterial, orders 2/3). The production CK path
      is now correctness-gated, not just timed.
- [x] **R-002:** `seas_spatial_dyn_driver` compiles + links; `--deriv-cache` is in
      the binary; default-off ⇒ `SetDerivMode` not called ⇒ byte-identical
      (verified by inspection + the unchanged operator code + green parity test).
      End-to-end driver run not exercised locally (no-local-full-mesh rule); the
      flag activates the cache on Frontera.
- [x] Regression unchanged: `wave_operator` 23/23, `bimaterial_parity` 9/9,
      `ader_ck_predictor` 11/11, `ader_linear_wave_equivalence` 4/4.

## Unresolved Findings (deferred per instruction)
- **R-003** (ParMesh parallel coverage) — not added; the cached kernel is
  element-local and the np=1 recursion test now covers the algorithm, but a
  parallel test remains the right belt-and-suspenders.
- **R-004** (homogeneous-order guard in the Cached branch) — not added.
- **R-005** (`SetDerivMode` idempotent early-out) — not added. Note: the driver
  calls `SetDerivMode(Cached)` exactly ONCE (outside the time loop), so the
  rebuild-per-call footgun is not triggered by the R-002 wiring.

## New Tests
- `RecursionCachedVsOnTheFlyMaxRelErr` + 4 `TEST_LE` assertions in
  `test_wave_operator_spatial_derivative.cpp` — covers R-001.

## Production note (not a code change)
To activate on Frontera, add `--deriv-cache` to the `ibrun ./seas_spatial_dyn_driver`
line in the tpv*_spatial sbatch jobs. This is a deliberate opt-in (changes results
at round-off); left to the user. Compare a `--deriv-cache` run's on-fault traces
against a baseline run within a physical tolerance (R-002 end-to-end acceptance).

## Ready for Re-Review: YES (R-001, R-002)
