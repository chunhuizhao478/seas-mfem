# Code Review: post-revert audit of feature/paraview-compaction (2026-05-16, round 3)

**Date:** 2026-05-16 (round 3)

**Context:** Fresh adversarial audit of the changes that landed in this branch since `main`, with particular focus on what was added or modified during the morning R-001..R-006 fix pass, the afternoon hard-cap revert, and the R-104..R-106 round-2 fixes. The fix agent reported PASS for all of R-101..R-106 (with R-101..R-103 marked OBSOLETE-by-revert). This round walks the resulting code with no benefit of the doubt to look for *new* bugs introduced by all the back-and-forth, and to verify that the regression tests actually catch the bugs they claim to guard.

**Verdict preview:** the code itself is now correct — the BP5 sbatch can be resubmitted to Frontera. But two issues materially weaken the post-fix safety net: (1) every one of the four new regression tests is unwired from `make test`, so a future PR can re-introduce any of the R-001..R-006 / R-104 bugs and CI will pass; (2) the sbatch's "bench_out probes therefore still write at 1 yr / 1 s / 0.01 s for SCEC compliance" caveat is *factually wrong* — bench_out actually writes at 0.1 yr / 0.1 s / 0.1 s, which means an operator following that comment to estimate disk size will be 10× off.

## Review Scope
- Branch: `feature/paraview-compaction` (now ~14 commits ahead of `main`).
- Files re-reviewed (changed since the round-2 audit):
  - `miniapps/seas/io/paraview_output.hpp` (hard cap reverted; `RecomputeIntervalForCapWithFloor` deleted; `SnapshotCapAwareInterval` restored to interseismic-only soft cap; stale doc-comment fixed at line 296-304).
  - `miniapps/seas/tests/verification/bp5_verification_full.cpp` (R-105 comment expansion at 1571-1591; R-106 env-var decoupling at 1997-2009).
  - `miniapps/seas/Makefile` (R-104 test build rules at lines 231, 445, 1070-1072, 1656-1658, 2250-2251).
  - `miniapps/seas/tests/unit/test_paraview_schedule_cap.cpp` (R-002 test re-targeted to verify SOFT cap at line 476-510).
  - `miniapps/seas/tests/unit/test_paraview_rank0_warning_gate.cpp` (NEW MPI-only R-104 test).
  - `miniapps/seas/tests/unit/test_bp5_petsc_ts_zero_fault_rank.cpp` (NEW R-003 test).
  - `miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch` (cap flag dropped at line 207-210; comment block rewritten).
- Plan / context: previous round's `FIX_REPORT.md` and `REVIEW.md`, plus `regime_budget_2026-05-16.md`.
- Local build verification: all four production drivers (BP5 + TPV102/104/205) build clean; schedule-cap tests 23/23, R-003 7/7, R-104 PASS under `mpirun -np 2`, kinematics 31/31.

## Findings

### [R-201] [MODERATE] [Makefile:1931] — Four new regression tests are unwired from `make test`; future commits can silently regress any of R-001..R-006 / R-104 / R-003 / kinematics

**Category:** EDGE_CASE (test coverage gap)

**Description:**
The `make test` umbrella target (line 1931 onward) lists the test binaries to build and run as routine CI. Four tests added during this branch's lifetime are NOT included:

| Test                                       | Guards finding | Standalone target          |
|--------------------------------------------|----------------|----------------------------|
| `seas_test_paraview_schedule_cap`          | R-001/R-005    | `test-paraview-schedule-cap` |
| `seas_test_kinematics_field_set`           | split-bulk R-001 (from May 12) | `test-kinematics-field-set` |
| `seas_test_bp5_petsc_ts_zero_fault_rank`   | R-003          | `test-bp5-petsc-ts-zero-fault-rank` |
| `seas_test_paraview_rank0_warning_gate`    | R-001 / R-104  | `test-paraview-rank0-warning-gate` |

Verified by:
```
$ grep -E "test-paraview-schedule-cap|test-kinematics-field-set|test-bp5-petsc-ts-zero-fault-rank|test-paraview-rank0-warning-gate" miniapps/seas/Makefile
2238:test-paraview-schedule-cap: seas_test_paraview_schedule_cap
2241:test-kinematics-field-set: seas_test_kinematics_field_set
2245:test-bp5-petsc-ts-zero-fault-rank: seas_test_bp5_petsc_ts_zero_fault_rank
2250:test-paraview-rank0-warning-gate: seas_test_paraview_rank0_warning_gate
```

All four appear only at their own definition lines. None is a prerequisite of `test:` at line 1931.

This was already noted as a deliberate deviation in the May-12 `STATUS_2026-05-12.md` ("`make test` doesn't run `seas_test_kinematics_field_set` automatically — existing PV unit tests aren't in `SEQ_MINIAPPS` either"), but that note is **not** a justification for adding three MORE unwired tests on top of it. The whole point of writing regression tests is to catch regressions *automatically*. A guard that only fires when someone remembers to invoke it manually is no guard at all.

**Trigger:**
Any future PR that breaks one of the gated invariants — the rank-0 warning gate, the padded-shrink Memory invariant, the kinematics-only volume PV registration, the soft-cap interseismic-only contract. None of these would be caught by `make test`.

**Actual behavior:**
`make test` passes; the regression ships; production breaks on Frontera.

**Expected behavior:**
`make test` invokes the four new tests (or at minimum the three that catch a CRITICAL or MODERATE issue: schedule-cap, bp5-petsc-ts-zero-fault-rank, paraview-rank0-warning-gate; kinematics-field-set is a MODERATE deviation guard).

**Suggested fix:**

```diff
 # miniapps/seas/Makefile, line 1931 onwards — extend the `test` umbrella.
 test: test-friction test-state-evolution test-antiplane test-fault-operator \
       ...
       test-tpv102-setup test-tpv102-local \
       test-fault-face-flux-frame \
       test-fault-face-flux-frame-and-flux \
+      test-paraview-schedule-cap \
+      test-kinematics-field-set \
+      test-bp5-petsc-ts-zero-fault-rank \
       ...
```

The MPI-only `test-paraview-rank0-warning-gate` should stay out of `test:` (because `make test` runs serially on the dev box and the MPI launch may not be available everywhere), but should be added to a parallel umbrella if one exists, OR documented in the project README as a required pre-submit step. The current state — "run if you remember" — is not adequate.

**Test case:**
This is a Makefile-wiring fix, not a code fix. Verification is mechanical:
```bash
$ cd miniapps/seas
$ make --dry-run test 2>&1 | grep -c seas_test_paraview_schedule_cap
# Pre-fix: 0.  Post-fix: 1 (the umbrella's prereq build step).
$ make --dry-run test 2>&1 | grep -c "./seas_test_bp5_petsc_ts_zero_fault_rank"
# Pre-fix: 0.  Post-fix: 1 (the umbrella's run step).
```

If the Makefile change is applied correctly, `make test` should depend on (and run) the three serial tests.

---

### [R-202] [MODERATE] [bp5_phase6_paraview_zfp_normal_48hr.sbatch:82-88] — Sbatch "IMPORTANT CAVEAT" block cites the wrong file, wrong line numbers, AND wrong cadences for `BP5BenchmarkOutput::OutputInterval`

**Category:** DEVIATION (documentation says one thing, code does another)

**Description:**
The sbatch comment block at lines 82-88 says:

```
# IMPORTANT CAVEAT.  --paraview-dt-co/-nu/-inter-yr only affect ParaView
# fault/volume snapshots.  The station probe files (*_fltst_*.txt the
# SCEC benchmark consumes) go through BP5BenchmarkOutput::OutputInterval
# which constructs a FRESH AdaptiveSchedule from hardcoded defaults
# and ignores these CLI flags (paraview_output.hpp:1589-1593).  The
# bench_out probes therefore still write at 1 yr / 1 s / 0.01 s for
# SCEC compliance.  Only the ParaView .vtkhdf size is reduced here.
```

Three independent factual errors:

1. **Wrong file path.** `BP5BenchmarkOutput::OutputInterval` lives in
   `miniapps/seas/io/bp5_benchmark_output.hpp`, NOT
   `miniapps/seas/io/paraview_output.hpp`. The function defined at
   `paraview_output.hpp:1635` is `ParaViewOutput::OutputInterval`
   (a different static helper, used by `test_io.cpp`).
2. **Wrong line numbers.** Even within `paraview_output.hpp`, line
   1589-1593 does NOT have `OutputInterval`; the file's actual line
   1635 has the unrelated `ParaViewOutput::OutputInterval` helper.
3. **Wrong cadences.** The actual `BP5BenchmarkOutput::OutputInterval`
   at `bp5_benchmark_output.hpp:794-808` is:
   ```cpp
   if      (V_max > 1e-3) return 0.1;                              // coseismic
   else if (V_max > 1e-6) return 0.1;                              // nucleation
   else                   return 0.1 * BP5Params::seconds_per_year; // interseismic
   ```
   So bench_out probes write at **0.1 s / 0.1 s / 0.1 yr**, not
   "1 yr / 1 s / 0.01 s" as the sbatch claims.

The user-visible consequence: an operator reading the sbatch comment to estimate the size of `*_fltst_*.txt` station files will be **10× off** in interseismic (0.1 yr cadence produces 10× more rows than 1 yr) and **10×–100×** off in coseismic (0.1 s cadence vs the claimed 0.01 s). Across a 250 yr / 1 event run, the actual probe-file size estimate is dominated by interseismic at 0.1 yr × 250 yr × 10 stations × ~8 columns × ~25 bytes/value ≈ 5 MB per station × 10 stations = 50 MB — not catastrophic in absolute terms, but the claim of "1 yr" cadence implies 5 MB total which the user might use to under-provision disk for the broader job.

**Trigger:**
Reading the sbatch comment to estimate probe-file disk usage.

**Actual behavior:**
The cited file, line numbers, and cadences are all wrong. The reader is misled about which code path produces probe files and what cadences it uses.

**Expected behavior:**
The cited file path matches reality; line numbers are accurate or omitted; cadences match the actual `OutputInterval` implementation.

**Suggested fix:**

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

**Test case:**
This is a doc fix; verification is by code-grep:
```bash
$ grep -n "OutputInterval" miniapps/seas/io/bp5_benchmark_output.hpp
794:   static real_t OutputInterval(real_t V_max)
$ grep -A 10 "static real_t OutputInterval" \
        miniapps/seas/io/bp5_benchmark_output.hpp | grep -E "return"
# Should show: return 0.1; (×2), return 0.1 * BP5Params::seconds_per_year;
```

The sbatch comment should be edited so its claimed cadences match the actual `return` values.

---

### [R-203] [LOW] [paraview_output.hpp:275-277, 322-324] — R-006 `regime` parameter on `RecomputeIntervalForCap` is structurally dead code; warning will always say "interseismic"

**Category:** QUALITY (dead parameter, misleading API)

**Description:**
`RecomputeIntervalForCap` (paraview_output.hpp:275-277) takes a defaulted
`int regime = 0` parameter, plumbed in by the R-006 fix so the warning
text can name the regime where the cap was exhausted. The warning at
line 322-324 looks it up:

```cpp
const char *regime_name = (regime == 2) ? "coseismic"
                        : (regime == 1) ? "nucleation"
                                        : "interseismic";
```

But there is now only ONE caller of `RecomputeIntervalForCap`, the
internal `SnapshotCapAwareInterval` at line 1875-1896. And that caller
has an early-exit at line 1879:

```cpp
if (regime != 0 || adaptive_.max_total_snapshots <= 0)
{
   return base;
}
// ... only here does RecomputeIntervalForCap get called ...
return adaptive_.RecomputeIntervalForCap(time_to_end,
                                         total_snapshots_written_,
                                         regime);
```

So `regime` is provably `== 0` at the call site (the `regime != 0`
branch returned before getting here). Therefore the `regime_name`
lookup will always pick the `else` branch and produce "interseismic".
The "coseismic" / "nucleation" cases are unreachable.

The previous-round `RecomputeIntervalForCapWithFloor` had a separate
purpose for the parameter (the floor variant called the original
method on cap-exhaustion with the actual regime). With that function
deleted, the parameter has no real consumer.

This is dead-code-style, not a bug per se. But a future reader will
spend time wondering why the lookup exists if the only path that gets
here is interseismic, and may incorrectly conclude there's a bug in
either the caller or the regime classification.

**Trigger:**
A reader auditing `RecomputeIntervalForCap` to understand the warning
text. They will trace the `regime` parameter back to its only call site
and conclude the parameter is vestigial.

**Actual behavior:**
The `regime_name` lookup compiles into a chain of comparisons that
always picks "interseismic".

**Expected behavior:**
Either the parameter and the lookup are removed (since only
interseismic can reach this code path), OR the lookup is moved into
`SnapshotCapAwareInterval` so it's adjacent to the regime decision,
OR the comment explicitly notes that the `coseismic` / `nucleation`
branches are reachable only by external callers that bypass the
soft-cap early-exit.

**Suggested fix:**
Simplest — drop the parameter and hardcode "interseismic":

```diff
       real_t RecomputeIntervalForCap(real_t time_to_end,
-                                     int snapshots_so_far,
-                                     int regime = 0) const
+                                     int snapshots_so_far) const
       {
          if (max_total_snapshots <= 0) { return dt_interseismic; }
          const int remaining_budget = max_total_snapshots - snapshots_so_far;
          if (remaining_budget <= 0)
          {
             ...
             if (_rank == 0)
 #endif
             {
-               const char *regime_name = (regime == 2) ? "coseismic"
-                                       : (regime == 1) ? "nucleation"
-                                                       : "interseismic";
                mfem::out << "ParaViewOutput: max_total_snapshots="
                          << max_total_snapshots
-                         << " exhausted (detected in "
-                         << regime_name
-                         << " regime; snapshots_so_far="
+                         << " exhausted in the interseismic regime "
+                         "(coseismic / nucleation regimes always use "
+                         "their natural cadence and never reach this "
+                         "branch; snapshots_so_far="
                          << snapshots_so_far
                          << ...
```

And update the single caller:

```diff
       const real_t time_to_end = total_run_time_ - time;
       return adaptive_.RecomputeIntervalForCap(time_to_end,
-                                               total_snapshots_written_,
-                                               regime);
+                                               total_snapshots_written_);
```

No test case — this is purely a dead-code cleanup; no observable behaviour change.

---

### [R-204] [LOW] [test_bp5_petsc_ts_zero_fault_rank.cpp:68-87] — R-003 regression test re-implements the padded-shrink pattern inline; a future driver refactor that drops the trick won't fail the test

**Category:** EDGE_CASE (incomplete test coverage)

**Description:**
The R-003 regression test at `test_bp5_petsc_ts_zero_fault_rank.cpp:68-87`
re-implements the padded-shrink pattern locally:

```cpp
const int actual = 0;
Vector state;
const int padded = std::max(actual, 1);
state.SetSize(padded);
state.SetSize(actual);
TEST_ASSERT(!state.GetMemory().Empty(), ...);
```

This verifies that the MFEM `Vector::SetSize` semantics the fix depends
on are still in effect — useful for catching an MFEM upgrade that
breaks the contract. But it does NOT verify that the BP5 driver itself
still applies the workaround. If a future commit reverts the driver to
the original `Vector state(fault_op.StateSize())` line, this test would
still pass (because the test does its own padded-shrink inline).

The crash would then manifest only on Frontera at the end of TSSolve,
on the next production submission — exactly the scenario the test was
written to catch.

**Trigger:**
A future commit that refactors `bp5_verification_full.cpp:1572-1591`
back to the naive `Vector state(fault_op.StateSize())` line, without
realising the workaround is load-bearing.

**Actual behavior:**
The local-MFEM test passes (MFEM `SetSize` semantics unchanged); the
driver-side regression is not caught.

**Expected behavior:**
The test directly invokes the driver code path, OR static-greps the
driver file for the expected workaround pattern.

**Suggested fix:**
Add a second sub-test that statically asserts the driver file contains
the workaround (a grep-style test, NOT compiling the driver):

```cpp
// In test_bp5_petsc_ts_zero_fault_rank.cpp, after the existing sub-tests:

// -- R-204 (REVIEW.md round 3): static check that the BP5 driver still
//    applies the padded-shrink workaround.  Without this check, a
//    refactor that reverts to `Vector state(fault_op.StateSize())`
//    would crash on Frontera but pass this unit test.
{
   const std::string driver_path =
      "miniapps/seas/tests/verification/bp5_verification_full.cpp";
   std::ifstream driver(driver_path);
   if (driver.is_open())
   {
      std::stringstream buf;
      buf << driver.rdbuf();
      const std::string src = buf.str();
      // Look for the canonical padded-shrink pattern.  Any sequence
      // that allocates >=1 element then shrinks back is acceptable;
      // the test is conservative and just requires the literal
      // padded/actual idiom from the R-003 comment.
      const bool has_padded = src.find("padded = std::max(actual, 1)")
                              != std::string::npos;
      const bool has_shrink = src.find("state.SetSize(actual)")
                              != std::string::npos;
      TEST_ASSERT(has_padded && has_shrink,
                  "R-204: BP5 driver must still contain the R-003 "
                  "padded-shrink workaround at "
                  << driver_path << "; a future refactor that removed "
                  "it would crash on Frontera on the next zero-fault-"
                  "DOF rank.");
   }
   else
   {
      std::cout << "  INFO: R-204 driver-source check skipped — "
                   "driver file not found at " << driver_path
                << " (test run from unexpected CWD?)\n";
   }
}
```

(Path resolution against CWD is fragile but acceptable for a local
dev check; the test gracefully skips if the file isn't reachable.)

**Test case:**
The sub-test above IS the test. If applied, the existing test binary
gains an 8th assertion. The assertion FAILS if the driver workaround
is removed.

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-201, R-202)
- Low issues: 2 (R-203, R-204)
- Plan compliance: **FULL** — R-001..R-006 are all resolved (R-002 by revert per user direction; R-001, R-004, R-005, R-006 by code change; R-003 by driver workaround). R-104, R-105, R-106 from round 2 are all addressed. The branch is functionally ready.
- Verdict: **PASS WITH FIXES** — R-201 (wire the new tests into `make test`) and R-202 (fix the misleading sbatch comment) should be addressed before the next Frontera resubmit, because they affect the post-fix safety net and operator documentation respectively. Neither blocks the run itself. R-203 and R-204 are quality cleanups for a follow-up commit.

## Unreviewed Areas
- The MFEM-side `linalg/petsc.cpp:899` patch (the alternative R-003 fix
  that would let `pdata.Empty()` pass for zero-local-size aliases) is
  STILL pending. The driver-side workaround at
  `bp5_verification_full.cpp:1572-1591` is in place, the R-105 comment
  flags the fragile dependency, and `test_bp5_petsc_ts_zero_fault_rank`
  pins the underlying invariant — but R-204 (above) notes the driver
  refactor isn't directly tested. The right long-term fix is still the
  upstream MFEM patch; tracking that is out of scope for this review.
- The R-104 MPI test was verified by running under `mpirun -np 2`
  locally. The behaviour under np >= 3 was not exercised. The test
  uses `MPI_Allreduce(MPI_SUM)` over a binary local, so it scales
  trivially to any rank count, but a paranoid follow-up could
  parameterise the recommended `make test-paraview-rank0-warning-gate`
  target to also try np=4 / np=8.
- The hard-cap revert removed code without removing the historical
  context from `paraview_output.hpp` (R-002 fix description is gone
  but the round-2 revert comment at lines 1857-1866 of
  `SnapshotCapAwareInterval` references it). Not flagged as a finding
  because the existing comments document the round-2 revert clearly
  enough.
- `bp5_verification_full.cpp` was scanned for OTHER places that
  construct a `Vector` sized by `fault_op.StateSize()` or
  `fault_op.SlipSize()`; only `state` at line 1572 is so constructed.
  All other paths (`pv_local_slip`, `pv_local_state`, etc.) use
  explicit `SetSize(N)` with N >= 2 * num_local_dofs which is
  positive on zero-fault-DOF ranks (zero, but SetSize(0) on a
  freshly-default-constructed Vector is safe-but-empty just like
  Vector(0) — same crash if those vectors were ever
  PlaceMemory'd into a PETSc Vec). They are NOT placed into PETSc
  via PetscODESolver, so they're fine. Confirmed by grep: only
  `state` is passed to `petsc_ode->Run(...)`.
