# Fix Report — REVIEW.md round 3 (BP5 Phase 6 — post safety-net cleanup)

**Date:** 2026-05-16 (round 3)
**Branch:** `feature/paraview-compaction`
**Review document:** `REVIEW.md` (round 3, 4 findings: R-201, R-202, R-203, R-204)

## Summary

- Findings addressed: **4 of 4**
- Files modified: **4** (Makefile, sbatch, paraview_output.hpp, test_bp5_petsc_ts_zero_fault_rank.cpp)
- Tests added: **1 new sub-test** (R-204 driver-source grep), appended to existing R-003 test binary
- Test suite: **PASS** locally — schedule-cap 23/23, R-003+R-204 8/8 (was 7/7), kinematics 31/31, R-104 PASS under `mpirun -np 2`
- All four production drivers rebuild clean.

## Changes Made

### R-201 [MODERATE] — Wire new tests into `make test` umbrella

**File:** `miniapps/seas/Makefile:1973-1979`

Added the three serial regression tests as prerequisites of the seas `test:` target. The MPI-only test stays out (per the reviewer's explicit guidance — `mpirun` is not guaranteed available in every dev environment) with an inline comment directing users to invoke it manually.

```diff
       test-ader-tpv102-smoke \
       test-R001-driver-defaults-to-ader \
-      test-R002-driver-init-total-q
+      test-R002-driver-init-total-q \
+      test-paraview-schedule-cap \
+      test-kinematics-field-set \
+      test-bp5-petsc-ts-zero-fault-rank
+# Note: `test-paraview-rank0-warning-gate` is INTENTIONALLY excluded
+# from the `test:` umbrella above — it requires `mpirun -np 2` which
+# is not guaranteed available in every dev environment.  Run it
+# manually as a pre-submit check:
+#     make test-paraview-rank0-warning-gate
```

Verification:
```
$ cd miniapps/seas && make --dry-run test | grep -c "./seas_test_paraview_schedule_cap"
1   # was 0 pre-fix
$ make --dry-run test | grep -c "./seas_test_bp5_petsc_ts_zero_fault_rank"
1   # was 0 pre-fix
$ make --dry-run test | grep -c "rank0_warning_gate"
0   # MPI-only test correctly excluded
```

### R-202 [MODERATE] — Fix sbatch bench_out caveat (wrong file/line/cadence)

**File:** `miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch:82-93`

Replaced the misleading IMPORTANT CAVEAT block with the corrected file path (`bp5_benchmark_output.hpp`, no inline line number — line numbers go stale), and the actual cadences (0.1 s coseismic / 0.1 s nucleation / 0.1 yr interseismic from `BP5BenchmarkOutput::OutputInterval`).

```diff
 # IMPORTANT CAVEAT.  --paraview-dt-co/-nu/-inter-yr only affect ParaView
 # fault/volume snapshots.  The station probe files (*_fltst_*.txt the
 # SCEC benchmark consumes) go through BP5BenchmarkOutput::OutputInterval
-# which constructs a FRESH AdaptiveSchedule from hardcoded defaults
-# and ignores these CLI flags (paraview_output.hpp:1589-1593).  The
-# bench_out probes therefore still write at 1 yr / 1 s / 0.01 s for
-# SCEC compliance.  Only the ParaView .vtkhdf size is reduced here.
+# which has its OWN hard-coded thresholds and ignores these CLI flags
+# (miniapps/seas/io/bp5_benchmark_output.hpp, the
+# `BP5BenchmarkOutput::OutputInterval` static method).  The bench_out
+# probes therefore still write at:
+#     coseismic    V > 1e-3 m/s         every 0.1 s    (SCEC spec)
+#     nucleation   1e-6 < V <= 1e-3     every 0.1 s
+#     interseismic V <= 1e-6 m/s        every 0.1 yr   (SCEC spec)
+# regardless of the --paraview-dt-* flags.  Only the ParaView .vtkhdf
+# size is reduced by those flags.
```

The new comment is verifiable against the actual code:
```
$ grep -A 12 "static real_t OutputInterval" miniapps/seas/io/bp5_benchmark_output.hpp
   static real_t OutputInterval(real_t V_max)
   {
      if (V_max > 1e-3)         return 0.1;                          // coseismic
      else if (V_max > 1e-6)    return 0.1;                          // nucleation
      else                      return 0.1 * BP5Params::seconds_per_year;  // interseismic
   }
```

### R-203 [LOW] — Drop dead `regime` parameter from `RecomputeIntervalForCap`

**File:** `miniapps/seas/io/paraview_output.hpp:275-276` (function signature) + `316-336` (warning text) + `1893-1894` (call site)

The parameter was structurally dead — the only caller is `SnapshotCapAwareInterval` which early-exits for non-interseismic regimes, so `regime` always arrived as `0`. Removed the parameter, removed the dead `regime_name` lookup, and hardcoded "interseismic" in the warning text with an explanatory clause for the reader.

```diff
       real_t RecomputeIntervalForCap(real_t time_to_end,
-                                     int snapshots_so_far,
-                                     int regime = 0) const
+                                     int snapshots_so_far) const
       {
          ...
                {
-                  const char *regime_name = (regime == 2) ? "coseismic"
-                                          : (regime == 1) ? "nucleation"
-                                                          : "interseismic";
                   mfem::out << "ParaViewOutput: max_total_snapshots="
                             << max_total_snapshots
-                            << " exhausted (detected in "
-                            << regime_name
-                            << " regime; snapshots_so_far="
+                            << " exhausted in the interseismic regime "
+                               "(coseismic / nucleation regimes always "
+                               "use their natural cadence and never reach "
+                               "this branch; snapshots_so_far="
                             << snapshots_so_far
                             ...

       const real_t time_to_end = total_run_time_ - time;
       return adaptive_.RecomputeIntervalForCap(time_to_end,
-                                               total_snapshots_written_,
-                                               regime);
+                                               total_snapshots_written_);
```

Schedule-cap test still 23/23 PASS after the change. The new warning text appears in test output and accurately describes the only branch that can reach the warning.

### R-204 [LOW] — Add driver-source grep sub-test to R-003 test

**File:** `miniapps/seas/tests/unit/test_bp5_petsc_ts_zero_fault_rank.cpp:34-37` (new `#include`s) + `107-141` (new sub-test)

The existing test verifies the underlying MFEM `Vector::SetSize` invariant the workaround depends on, but didn't verify that the BP5 driver itself still uses the workaround. Added a static-grep sub-test that opens `tests/verification/bp5_verification_full.cpp` and asserts the literal `padded = std::max(actual, 1)` and `state.SetSize(actual)` substrings are present.

The sub-test gracefully skips (with an INFO message, not a failure) if the driver source isn't reachable from the current working directory — so the test remains runnable in arbitrary contexts. The seas `make test-bp5-petsc-ts-zero-fault-rank` invokes from `miniapps/seas/` where the relative path resolves correctly.

```
$ ./seas_test_bp5_petsc_ts_zero_fault_rank
  ...
  PASSED: R-204: BP5 driver must still contain the R-003 padded-shrink workaround
    at tests/verification/bp5_verification_full.cpp;
    a refactor that removed it would crash on Frontera on the next zero-fault-DOF rank.
    (has_padded=1 has_shrink=1)
=== Summary: 8 / 8 passed; 0 failed ===
```

If a future commit removes the workaround from the driver, the test fails with `has_padded=0 has_shrink=0` and a pointer back to the original R-003 fix discussion.

## New Tests

| Test                                                              | Purpose                                          | How to run                                       |
|-------------------------------------------------------------------|--------------------------------------------------|--------------------------------------------------|
| R-204 sub-test in `seas_test_bp5_petsc_ts_zero_fault_rank`         | Verify BP5 driver still applies the R-003 padded-shrink workaround | `make test-bp5-petsc-ts-zero-fault-rank` (also runs as part of `make test` post-R-201) |

No standalone new binary; the sub-test is appended to the existing R-003 test.

## Build / Test Verification

```
$ cd miniapps/seas
$ make seas_bp5_full seas_tpv102_driver seas_tpv104_driver seas_tpv205_driver -j4
  ... all four build clean.

$ make seas_test_paraview_schedule_cap seas_test_kinematics_field_set \
       seas_test_bp5_petsc_ts_zero_fault_rank seas_test_paraview_rank0_warning_gate -j4
  ... all four test binaries link.

$ ./seas_test_paraview_schedule_cap
  === Summary: 23 / 23 passed; 0 failed ===

$ ./seas_test_bp5_petsc_ts_zero_fault_rank
  === Summary: 8 / 8 passed; 0 failed ===   (was 7/7 pre-R-204)

$ ./seas_test_kinematics_field_set
  === Summary: 31 / 31 passed; 0 failed ===

$ mpirun --oversubscribe -np 2 ./seas_test_paraview_rank0_warning_gate
  PASS [R-104]: cap-exhausted warning fires exactly once across 2 ranks, on rank 0.

$ make --dry-run test | grep -c "./seas_test_paraview_schedule_cap"
  1   (R-201 verified — was 0 pre-fix)

$ make --dry-run test | grep -c "./seas_test_bp5_petsc_ts_zero_fault_rank"
  1   (R-201 verified — was 0 pre-fix)
```

## Unresolved Findings

None.

## Deviations from the Review's Suggested Fixes

None substantive. The R-201 fix added the comment block recommended by the reviewer to flag the deliberately-excluded MPI-only test. The R-204 fix uses a slightly different grep pattern than the reviewer sketched (the reviewer's `"padded = std::max(actual, 1)"` matches the literal in the driver — kept as-is), but adds a clearer error-message diagnostic line that prints which of the two substrings (`has_padded` / `has_shrink`) is missing.

## Notes for Reviewer Re-Review

- The sbatch comment in R-202 deliberately omits a specific line number for `BP5BenchmarkOutput::OutputInterval` to avoid the same "line numbers go stale" problem that the original comment had (line 1589 cited, actual was 794, then this round 794). The file path + function name + signature are enough to grep for.

- The R-203 docstring update notes that the deleted `regime` parameter was REPLACED by an explanatory clause in the warning text ("coseismic / nucleation regimes always use their natural cadence and never reach this branch"). A future reader who wonders why the warning hardcodes "interseismic" has the answer in the warning itself.

- The R-204 sub-test is a static-grep, not a structural / AST-based check. A pathological refactor (e.g., renaming `padded` to `_padded` or splitting the `std::max(actual, 1)` across lines) would defeat it. This is intentional — the goal is to catch innocent removals of the workaround, not to bullet-proof against every conceivable refactor. The R-105 inline comment in the driver still documents the workaround clearly enough that a maintainer who edits that block will see they're touching load-bearing code.

- The `make test` umbrella now picks up 3 new test binaries when invoked from `miniapps/seas/`. From the project root, `make test` runs MFEM's top-level test target which does not recurse into the seas miniapp's `test:` — so a CI that builds from the root will still miss the new tests. Wiring up the project-root umbrella to delegate to the seas one is OUT OF SCOPE for this fix (the seas miniapp has never been part of MFEM's top-level test set); flagged here so the reviewer can decide whether to track it as a separate task.

## Ready for Re-Review: YES
