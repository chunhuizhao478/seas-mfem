# Code Review: PETSc TS restart — round 5 (post round-4 fixes)

**Date:** 2026-05-17
**Reviewer:** fresh adversarial pass on the round-4 fix work
**Objective:** Find new bugs the round-4 fixes introduced or that round-4 missed.  Treat this as a first-look review — do NOT just check off round-4 items.

## Review Scope

- Plan: `miniapps/seas/debug_document/paraview_output_debug_document/petsc_ts_restart_plan_2026-05-16.md`
- Files reviewed (current working tree):
  - `miniapps/seas/io/petsc_ts_checkpoint.hpp` (R-006 probe added)
  - `miniapps/seas/tests/unit/test_bp5_petsc_ts_restart.cpp` (Sub-test 3 + MPI init + R-004/R-005 changes)
  - `miniapps/seas/tests/verification/bp5_verification_full.cpp` (R-002 doc, R-003 scale fix, R-008 `v2_authoritative_dt`, R-009 comment)
  - `miniapps/seas/io/paraview_output.hpp` (R-004 setters; unchanged from round 4)
  - `setup_mfem.sh` (R-007 hint)
- Cross-checked: `linalg/petsc.cpp:4357-4394` (PetscODESolver::Run), `linalg/petsc.hpp:949-980` (PetscODESolver API), `linalg/operator.hpp:343-500` (TimeDependentOperator + ExplicitMult signature).
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`, prior REVIEW.md round 4 (now overwritten).

---

## Findings

### [R-001] [MODERATE] [test_bp5_petsc_ts_restart.cpp:main + all sub-tests] — Test races / corrupts files under `mpirun -np N>=2`

**Category:** BUG (latent — does not fire under the current Makefile target, but the plan explicitly mandated MPI launch)

**Description:**
The plan §"Testing Strategy" / build wiring originally specified `mpirun -np 2 ./seas_test_bp5_petsc_ts_restart`.  The implementer documented a deviation: the Makefile target now invokes `./seas_test_bp5_petsc_ts_restart` directly (serial).  Round-4 added `MPI_Init` to `main()` for Sub-test 3's PetscODESolver init, so the test is now half-MPI-aware: it INITIALISES MPI but never CHECKS rank/size.

Every sub-test passes `mpi=nullptr` to `WriteCheckpoint` / `WritePetscTSCheckpoint`, which hard-codes `rank = mpi ? mpi->Rank() : 0`.  Result: every rank writes to `*_checkpoint_r0.txt`.  Under `mpirun -np 2`:
- Both ranks call `MakeTmpDir("subtest1")` → `mkdir` race (probably harmless).
- Both ranks call `WriteCheckpoint("/tmp/.../subtest1/v1_only", ...)` → write to `/tmp/.../subtest1/v1_only_checkpoint_r0.txt` SIMULTANEOUSLY.  File contents corrupted / interleaved.
- Sub-test 2 / 3 / 4 / 7 / 8 / 9 all have the same race.
- Sub-test 3 additionally creates a PetscODESolver on `MPI_COMM_SELF` on each rank; each rank tries to write to the same V2 file via `WritePetscTSCheckpoint(..., nullptr)`.

The fact that the Makefile target doesn't use mpirun masks the bug, but:
1. Anyone running the test under `mpirun -np 2 ./seas_test_bp5_petsc_ts_restart` (matching the plan) hits silent file corruption.
2. The plan's `test-bp5-petsc-ts-restart` target was originally `$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./...` — the implementer changed it.  A future regression to the planned target reintroduces the race.

**Trigger:**
`mpirun -np 2 ./seas_test_bp5_petsc_ts_restart`

**Actual behavior:**
File-format sub-tests (1, 2, 4, 7, 8, 9, 10, 11, 12) race on `/tmp/seas_test_bp5_petsc_ts_restart/subtestN/*_checkpoint_r0.txt`.  Sub-test 3 races on its V2 file.  Asserts may pass or fail nondeterministically; in the worst case, the binary segfaults inside `std::ifstream` parsing of a half-written file.

**Expected behavior:**
Either:
(a) Detect `MPI_Comm_size > 1` at `main()` start and SKIP all sub-tests on non-rank-0 (with a clear "test is serial-only; running on rank 0 only" message); OR
(b) Use rank-suffixed `/tmp` paths so each rank has its own sandbox.

Option (a) is cleaner since the sub-tests are designed as serial unit tests.

**Suggested fix:**
```diff
 int main(int argc, char *argv[])
 {
 #ifdef MFEM_USE_PETSC
    int already_inited = 0;
    MPI_Initialized(&already_inited);
    if (!already_inited) { MPI_Init(&argc, &argv); }
    mfem::MFEMInitializePetsc(&argc, &argv, NULL, NULL);
+
+   // R-001 (REVIEW.md round 5): the sub-tests use file I/O with
+   // mpi=nullptr (hard-coding rank=0 in the filename), so under
+   // mpirun -np N>=2 every rank would race on the same /tmp paths.
+   // Detect parallel launch and run sub-tests only on rank 0; the
+   // rest spin idle until MPI_Finalize.
+   int mpi_size = 1, mpi_rank = 0;
+   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
+   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
+   if (mpi_size > 1 && mpi_rank != 0)
+   {
+      mfem::MFEMFinalizePetsc();
+      if (!already_inited) { MPI_Finalize(); }
+      return 0;
+   }
+   if (mpi_size > 1 && mpi_rank == 0)
+   {
+      std::cout << "INFO: launched under mpirun -np " << mpi_size
+                << "; sub-tests run on rank 0 only (file I/O is "
+                << "serial-only and would race on shared /tmp paths "
+                << "otherwise).  Use np=1 for the full suite.\n";
+   }
 #else
    (void)argc; (void)argv;
 #endif
```

**Test case:**
```bash
# Without fix:
mpirun -np 2 ./seas_test_bp5_petsc_ts_restart
# Expect: file corruption, asserts may pass/fail nondeterministically,
# OR segfault parsing half-written files.

# With fix:
mpirun -np 2 ./seas_test_bp5_petsc_ts_restart
# Expect: rank 0 runs all sub-tests as if serial (51/51 pass), rank 1
# prints "running on rank 0 only", both ranks exit 0.
```

---

### [R-002] [MODERATE] [test_bp5_petsc_ts_restart.cpp:Subtest3] — Test mirrors the driver's restart sequence but never invokes the actual driver source — deletions in the driver's V2 block go undetected

**Category:** DEVIATION (round-4 R-001 partial fix — Sub-test 3 proves PetscODESolver+file format work together, but does NOT prove the DRIVER's V2 restart block calls them correctly)

**Description:**
Sub-test 3 hand-rolls the V2 restart sequence:
```cpp
// Lines 894-897:
TSSetStepNumber(ts, static_cast<PetscInt>(r_step));
real_t t_post  = r_t;
real_t dt_post = r_dt;
ode.Run(y_B, t_post, dt_post, T_full);
```

This proves the V2 restart MACHINERY (file format + PetscODESolver semantics) works.  But it does NOT prove the DRIVER at `bp5_verification_full.cpp:2441-2563` contains the corresponding calls.

Concretely: if I delete these load-bearing lines from the driver
- `t = ts_t;` (line 2530)
- `current_dt = ts_dt_next;` (line 2533)
- `TSSetStepNumber(ts, static_cast<PetscInt>(ts_step));` (line 2522)

then Sub-test 3 still passes (because it has its own copy of those lines), but BP5 production restart is broken.  The other "grep" sub-tests (7, 8, 9) catch text strings but NOT these specific lines.

Round-4 R-001's intent was to verify end-to-end restart works.  The current Sub-test 3 verifies that the BUILDING BLOCKS work; it does NOT verify that the driver assembles them correctly.

**Trigger:**
A refactor that accidentally deletes line 2530 / 2533 / 2522 from the driver.  Sub-test 3 still passes.  BP5 production restart silently uses wrong t / dt / step.

**Actual behavior:**
Sub-test 3 passes regardless of what the driver does.

**Expected behavior:**
The test should fail if the driver's V2 restart block loses its load-bearing assignments.  Either add grep-style source checks (matching R-204 / R-007 / R-008 / R-009 patterns), OR refactor the driver's V2 block into a callable helper that the test can invoke directly.

**Suggested fix (grep-style, minimal):**
Add to Sub-test 3 (after the trajectory assertion, before the `}`):

```diff
    TEST_ASSERT(diff < bound,
                "Sub-test 3e (R-001 end-to-end): fresh-vs-restart "
                ...
                "does NOT actually work.");
+
+   // R-002 (REVIEW.md round 5): grep the driver source for the
+   // three load-bearing assignments in the V2 restart block.  Sub-test
+   // 3 above proves the building blocks work; this check proves the
+   // driver assembles them correctly.  Without this, a refactor that
+   // accidentally deletes any of `t = ts_t`, `current_dt = ts_dt_next`,
+   // or `TSSetStepNumber(ts, ...)` from the driver's V2 restart block
+   // would silently break BP5 production restart while sub-test 3
+   // continues to pass.
+   const std::string driver_path =
+      "tests/verification/bp5_verification_full.cpp";
+   std::ifstream driver(driver_path);
+   if (driver.is_open())
+   {
+      std::stringstream buf; buf << driver.rdbuf();
+      const std::string src = buf.str();
+      const bool has_t_assign  =
+         src.find("t          = ts_t;") != std::string::npos
+         || src.find("t = ts_t;") != std::string::npos;
+      const bool has_dt_assign =
+         src.find("current_dt = ts_dt_next;") != std::string::npos;
+      const bool has_setstep   =
+         src.find("TSSetStepNumber(ts, static_cast<PetscInt>(ts_step))")
+         != std::string::npos;
+      TEST_ASSERT(has_t_assign && has_dt_assign && has_setstep,
+                  "Sub-test 3f (R-002 round 5): driver V2 restart block "
+                  "must contain `t = ts_t`, `current_dt = ts_dt_next`, "
+                  "and `TSSetStepNumber(ts, static_cast<PetscInt>"
+                  "(ts_step))`.  has_t_assign=" << has_t_assign
+                  << " has_dt_assign=" << has_dt_assign
+                  << " has_setstep=" << has_setstep);
+   }
+   else
+   {
+      std::cout << "  INFO: Sub-test 3f driver-grep skipped — could "
+                   "not open " << driver_path << " (CWD?)\n";
+   }
```

**Test case:**
```cpp
// Demonstration: delete `current_dt = ts_dt_next;` from
// bp5_verification_full.cpp:2533 and re-run.
//
// Without the fix: ./seas_test_bp5_petsc_ts_restart → 51/51 PASS
//                  (Sub-test 3 doesn't catch it).
// With the fix:    ./seas_test_bp5_petsc_ts_restart → Sub-test 3f
//                  FAILS with "has_dt_assign=0".
```

---

### [R-003] [LOW] [petsc_ts_checkpoint.hpp:WritePetscTSCheckpoint probe] — V1 probe catches missing files but not empty / truncated ones

**Category:** EDGE_CASE

**Description:**
The R-006 round-4 fix added a probe before opening in `ios::app` mode:
```cpp
{
   std::ifstream probe(filename);
   MFEM_VERIFY(probe.good(),
               "WritePetscTSCheckpoint: prerequisite V1 checkpoint "
               "file " << filename << " does not exist; ...");
}
```

`probe.good()` returns true if the file exists and was opened successfully — regardless of contents.  If the file is EMPTY (0 bytes) or TRUNCATED (writer crashed mid-V1), the probe passes and `WritePetscTSCheckpoint` appends the V2 trailing block to a header-less file.  Downstream `ReadCheckpoint` then aborts with "expected 'SEAS_CHECKPOINT_V1', got '...'", pointing at the reader rather than at the writer that produced the half-formed file.

Not a critical bug — the system fails loudly, just with a confusing message — but the round-4 fix didn't fully address the failure mode it documented.

**Trigger:**
Any scenario where the V1 file exists but is empty or truncated — e.g., disk full mid-write, MPI_Abort during WriteCheckpoint, manual `> file` to zero-out for testing.

**Actual behavior:**
`WritePetscTSCheckpoint` appends V2 to an empty file.  `ReadCheckpoint` later reports confusing parse error.

**Expected behavior:**
Probe should also verify the file ends with a recognizable V1 trailer.  Cheapest check: stat for `size > min_v1_size` where `min_v1_size` is the size of an empty-vector V1 file (header + 0-element vectors ≈ 250 bytes).

**Suggested fix:**
Either accept this as a known limitation (cheapest), or add a size check:

```diff
    {
       std::ifstream probe(filename);
       MFEM_VERIFY(probe.good(),
                   "WritePetscTSCheckpoint: prerequisite V1 checkpoint "
                   "file " << filename << " does not exist; call "
                   "WriteCheckpoint on the same prefix first.");
+      // R-003 (REVIEW.md round 5): also verify the file is non-empty.
+      // probe.good() is true for empty files too, so without this
+      // check a truncated V1 file (e.g., from a crashed write) would
+      // have V2 appended to it and later confuse ReadCheckpoint.
+      probe.seekg(0, std::ios::end);
+      const auto file_size = probe.tellg();
+      MFEM_VERIFY(file_size > 0,
+                  "WritePetscTSCheckpoint: prerequisite V1 checkpoint "
+                  "file " << filename << " exists but is empty "
+                  "(prior WriteCheckpoint crashed mid-stream?).  "
+                  "Cannot safely append the V2 trailing block.");
    }
```

**Test case:**
N/A (would require simulating a crashed write; downgrade or drop if not worth the complexity).

---

### [R-004] [LOW] [paraview_output.hpp:SetLastCommittedCycle] — No input validation; arbitrary negative or sentinel-collision values accepted

**Category:** ASSUMPTION

**Description:**
The R-004 round-4 fix added:
```cpp
void SetLastCommittedCycle(int cycle) { last_committed_cycle_ = cycle; }
```

This setter accepts ANY int.  The CommitSchedule dedup logic at `paraview_output.hpp:1607-1609` uses:
```cpp
const bool same_step_as_last_commit =
   (last_committed_cycle_ != std::numeric_limits<int>::min())
   && (std::abs(time - last_write_time_) <= tol);
```

The first conjunct (`!= INT_MIN`) is the only thing distinguishing "fresh" from "post-commit" state.  If a V2 checkpoint accidentally stored `INT_MIN` (e.g., from a never-committed pre-checkpoint run that somehow wrote V2), the post-restart dedup would falsely treat the state as "no prior commit ever happened", over-bumping the counter on the first CommitSchedule.

A V2-aware writer always sets `paraview_last_committed_cycle` to either a real cycle number or `INT_MIN` (from the `mon->pv_out ? ... : INT_MIN` default), so a real `INT_MIN` value can legitimately appear in a V2 file (paraview disabled scenario).  Restoring `INT_MIN` is correct in that case.

But the setter does NOT distinguish "you intended INT_MIN" from "the V2 read returned garbage that happens to be INT_MIN".  For other corruption modes (e.g., a negative-but-not-INT_MIN value), the setter accepts it silently.

Compare R-007 round-4 (regime clamp) — that fix DID add validation.  The dedup-cycle setter omits it.

**Trigger:**
A corrupted V2 file that contains a negative `paraview_last_committed_cycle` other than `INT_MIN`.  Or any value that violates the implicit invariant "monotonically non-decreasing cycle".

**Actual behavior:**
Setter accepts silently.  Subsequent dedup behavior depends on whether the corrupted value happens to be `INT_MIN` or not.

**Expected behavior:**
Either clamp negative values to `INT_MIN` (treating "any negative" as "no prior commit"), or document explicitly that the setter is unvalidated and the V2 producer is trusted.

**Suggested fix:**
```diff
-   void SetLastCommittedCycle(int cycle) { last_committed_cycle_ = cycle; }
+   void SetLastCommittedCycle(int cycle)
+   {
+      // R-004 (REVIEW.md round 5): clamp negative values other than
+      // the documented INT_MIN sentinel to INT_MIN.  This treats any
+      // corrupted-but-negative V2 value as "no prior commit", which
+      // is the safest default — falsely identifying "fresh state"
+      // leads to one extra bump on the first CommitSchedule (minor
+      // over-count), whereas a corrupted positive value can cause
+      // legitimate dedups to misfire (silent under-count for the
+      // remainder of the run).
+      last_committed_cycle_ = (cycle < 0)
+                              ? std::numeric_limits<int>::min()
+                              : cycle;
+   }
```

**Test case:**
```cpp
// Add to sub-test 12 after the existing assertions:
pv_post.SetLastCommittedCycle(-42);
TEST_EQ(pv_post.GetLastCommittedCycle(),
        std::numeric_limits<int>::min(),
        "Sub-test 12 (R-004 round 5): negative cycle key clamped to "
        "INT_MIN (any negative is treated as 'no prior commit')");
```

---

### [R-005] [LOW] [test_bp5_petsc_ts_restart.cpp:Subtest3] — Mid-run `mid_dt_next` may be exactly 0 in unusual scenarios; test doesn't exercise R-003 driver fallback

**Category:** EDGE_CASE / QUALITY

**Description:**
Sub-test 3 saves `mid_dt_next` from `TSGetTimeStep` after a Run that ended at `T_mid` via `TS_EXACTFINALTIME_MATCHSTEP`.  Per PETSc source, this can produce a near-zero `mid_dt_next` because the controller may shrink dt at the final-step match.  In a typical run, `mid_dt_next` is some small positive value (the test logged it implicitly via the trajectory check).

But the test does NOT explicitly verify what `mid_dt_next` is, nor exercise the R-003 driver-side fallback (`if (current_dt <= 0.0) current_dt = dt_init;` at `bp5_verification_full.cpp:2539-2548`).  Sub-test 8's grep covers the source-text existence; nothing covers actual runtime behavior.

If a future PETSc upgrade changes the post-MATCHSTEP `TSGetTimeStep` behavior to return exactly 0 (instead of a small positive), the driver's fallback would fire and the test would silently continue, but a separate code path is now exercised that we have no regression test for.

**Trigger:**
A PETSc version change OR a checkpoint taken from a TSSetConvergedReason(TS_DIVERGED_*) state, where TSGetTimeStep returns 0.

**Actual behavior:**
The driver falls back to `dt_init` silently; no test catches misbehavior in the fallback path.

**Expected behavior:**
Either explicitly test the dt_next=0 case end-to-end, or accept the grep-test coverage as sufficient.

**Suggested fix:**
Add a small variant of Sub-test 3 that hand-writes V2 with `dt_next = 0`, drives the restart sequence (using `dt_init` as the fallback the driver would compute), and asserts trajectory continuity:

```cpp
// In Sub-test 3, after the main test:
{
   // Verify the R-003 driver-side fallback path produces a stable
   // restart even when V2 ts_dt_next == 0.  We simulate the driver's
   // fallback explicitly here (test doesn't invoke driver).
   const std::string prefix_zero = MakeTmpDir("subtest3_dtzero") + "/v2";
   WriteMinimalV1(prefix_zero);
   WritePetscTSCheckpoint(prefix_zero, T_mid, /*dt_next=*/0.0,
                          static_cast<int>(mid_step),
                          static_cast<int>(mid_rej),
                          0, -1e30, 0.0, 0,
                          std::numeric_limits<int>::min(), -1e30,
                          nullptr);
   real_t r_t2 = 0, r_dt2 = 0, r_lw2 = 0, r_vmax2 = 0, r_voltime2 = 0;
   int r_step2 = 0, r_rej2 = 0, r_snap2 = 0, r_regime2 = 0, r_commit2 = 0;
   const bool ok2 = ReadPetscTSCheckpoint(prefix_zero, r_t2, r_dt2,
                                          r_step2, r_rej2, r_snap2,
                                          r_lw2, r_vmax2, r_regime2,
                                          r_commit2, r_voltime2, nullptr);
   TEST_ASSERT(ok2, "Sub-test 3g: V2 read with dt_next=0 succeeds");
   TEST_DOUBLE_EQ(r_dt2, 0.0, "Sub-test 3g: dt_next=0 round-trips");

   // Apply the driver's R-003 fallback explicitly.
   real_t dt_post2 = r_dt2;
   if (dt_post2 <= 0.0) { dt_post2 = dt_init; }
   TEST_ASSERT(dt_post2 > 0.0,
               "Sub-test 3g (R-005 round 5): R-003 driver-side "
               "fallback produces a positive dt when V2 returns 0");

   // Run with the fallback dt; trajectory should still land at T_full.
   mfem::Vector y_B2(1); y_B2(0) = final_B;  // continue from y_B
   DecayOp op_B3;
   mfem::PetscODESolver ode2(MPI_COMM_SELF, "");
   ode2.Init(op_B3, mfem::PetscODESolver::ODE_SOLVER_GENERAL);
   mfem::petsc::TS ts2 = ode2;
   TSSetType(ts2, TSRK); TSRKSetType(ts2, TSRK5DP);
   TSAdapt tsad2; TSGetAdapt(ts2, &tsad2);
   TSAdaptSetType(tsad2, TSADAPTBASIC);
   ode2.SetAbsTol(1e-7); ode2.SetRelTol(1e-10);
   real_t t_post2 = T_full;  // already at T_full from main Sub-test 3
   // Push to T_full + epsilon to exercise one more step
   ode2.Run(y_B2, t_post2, dt_post2, T_full + 1e-3);
   TEST_ASSERT(std::isfinite(y_B2(0)),
               "Sub-test 3g: R-003 fallback produces finite state");
}
```

OR, accept this as a documented gap (the grep test in Sub-test 8 verifies the source text exists).

**Test case:** As above.

---

### [R-006] [LOW] [bp5_verification_full.cpp:R-008 guard] — Assertion at Run() site uses bit-exact `==` on real_t; may false-positive if any intervening code legitimately recomputes dt

**Category:** POSSIBLE / ASSUMPTION

**Description:**
The R-008 round-4 fix added:
```cpp
if (v2_authoritative_dt > 0.0)
{
   MFEM_VERIFY(current_dt == v2_authoritative_dt,
               "R-008: current_dt (" << current_dt
               << ") was modified between the V2 restart block "
               "and the Run() call ...");
}
```

This uses bit-exact `==` on `real_t`.  The assertion's intent is "no code between V2 block and Run() modified current_dt".  Today, that's true (verified by grep).

But a LEGITIMATE future change could introduce a no-op-looking modification that nonetheless changes the bit pattern.  E.g., `current_dt = std::min(current_dt, dt_max);` where dt_max > current_dt always — semantically a no-op, but the assignment may or may not produce the same bit pattern depending on compiler.  The assertion would fire on what is actually a correct change.

A `std::abs(current_dt - v2_authoritative_dt) < 1e-15 * v2_authoritative_dt` check would be more robust, but ALSO masks real bugs (a CFL clamp at the second ULP could be a real concern in some configurations).

Trade-off: bit-exact is stricter (catches more, but false-positives on harmless code motion).  Reviewer judgment: keep the bit-exact check for now since it's exactly what we want today.  Flag as POSSIBLE so the next reviewer sees the concern.

**Trigger:**
A future change inserts a value-preserving but bit-changing assignment to `current_dt` between V2 block and Run().

**Actual behavior:**
Assertion fires; restart aborts.

**Expected behavior:**
For value-preserving changes, no abort.  But this is a value judgment.

**Suggested fix:**
Document the trade-off in the assertion comment:

```diff
       if (v2_authoritative_dt > 0.0)
       {
+         // R-006 (REVIEW.md round 5): bit-exact `==` is intentional —
+         // any modification of current_dt between the V2 block and
+         // here, even a value-preserving one (e.g., `current_dt =
+         // std::min(current_dt, dt_max)`), would change the bit
+         // pattern under some compilers and fire this assertion.
+         // That's by design: any insertion HERE deserves a deliberate
+         // re-examination of whether V2 is still the authoritative
+         // source of post-restart dt.  If you legitimately need to
+         // clamp post-V2 dt, update `v2_authoritative_dt` in the
+         // same statement, OR widen this check to a relative
+         // tolerance with a documented bound.
          MFEM_VERIFY(current_dt == v2_authoritative_dt,
                      "R-008: current_dt (" << current_dt
```

**Test case:** N/A.

---

### [R-007] [LOW] [test_bp5_petsc_ts_restart.cpp:Subtest3] — `final_A` and `final_B` not asserted to be within atol of the analytic solution `exp(-T_full)`

**Category:** QUALITY (missing sanity check)

**Description:**
Sub-test 3 asserts `|final_A - final_B| < atol + rtol·|final_A|`.  This catches fresh-vs-restart divergence.  But if BOTH are equally wrong (e.g., a PETSc bug that consistently produces 2.0 instead of 0.368), the test passes silently.

Adding an absolute check against the analytic solution `exp(-1)` would catch this class of "consistently wrong" bugs.

**Suggested fix:**
```diff
    TEST_ASSERT(diff < bound,
                "Sub-test 3e (R-001 end-to-end): fresh-vs-restart "
                ...
                "does NOT actually work.");
+
+   // R-007 (REVIEW.md round 5): sanity check against the analytic
+   // solution.  Catches the failure mode "fresh and restart agree
+   // with each other but are both wrong" (e.g., PETSc adapter bug,
+   // RK type mismatch).  exp(-T_full) for T_full=1 is 0.367879...
+   const real_t exact = std::exp(-T_full);
+   TEST_ASSERT(std::abs(final_A - exact) < 1e-5,
+               "Sub-test 3e (R-007 round 5): fresh-run final state "
+               "agrees with the analytic solution exp(-T_full) — "
+               "y_A(T_full)=" << final_A << ", exact=" << exact
+               << ", diff=" << std::abs(final_A - exact)
+               << ".  Catches 'both runs equally wrong' bugs that "
+               "the relative comparison alone would miss.");
```

**Test case:** N/A.

---

## Summary

- Critical issues: **0**
- Moderate issues: **2** (R-001 mpirun race, R-002 driver-grep gap)
- Low issues: **5** (R-003, R-004, R-005, R-006, R-007)
- Plan compliance: **PARTIAL** — Sub-test 3 satisfies the round-4 spirit of R-001 but with the qualifications in R-002 (driver-internal calls unverified) and the unfixed parallel-launch race from R-001.
- Verdict: **PASS WITH FIXES** — the implementation correctly proves the PetscODESolver round-trip works; restart works **for the components tested**.  The two MODERATE findings (R-001, R-002) close the remaining "is the driver actually wired up correctly + does the test handle parallel launch" gaps.  None of the findings indicate the existing implementation is wrong; they identify coverage gaps and a parallel-launch race.

## Unreviewed Areas

- **The `seas_bp5_full --restart --petsc-ts` workflow under actual BP5 conditions** — only the building blocks (file format, PetscODESolver round-trip) are unit-tested; the BP5 driver's V2 block has not been exercised end-to-end.  The plan AC #14 (sbatch re-submission with `--restart`) is the gate for this.  Out of scope for the unit-test review round.
- **The `setup_mfem.sh` external-build path** — only in-tree builds were exercised this round.  No regressions reported.
- **R-008 round 4** — the bit-exact guard at the Run() site is documented but not exercised by any test (it only fires on broken futures).  Acceptable as a tripwire; flagged in R-006 above as a possible source of false-positive aborts.
- **Phase 2 / Phase 3 of the plan** — still deferred per the plan §"Dependencies".  No new findings.
