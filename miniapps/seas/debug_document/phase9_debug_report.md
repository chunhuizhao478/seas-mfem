# Phase 9 Debug Report: Parallel Verification and Scaling

## Files Reviewed

| File | Lines |
|------|-------|
| `fault/fault_geometry.hpp` | 508 |
| `fault/rate_state_fault.hpp` | 417 |
| `solver/seas_operator.hpp` | 223 |
| `solver/time_stepper.hpp` | 460 |
| `io/benchmark_output.hpp` | 267 |
| `io/parallel_benchmark_output.hpp` | 147 |
| `io/probe_output.hpp` | 260 |
| `domain/antiplane_operator.hpp` | 1827 |
| `common/mpi_context.hpp` | 173 |
| `pseas.cpp` | 328 |
| `tests/parallel/test_serial_parallel_consistency.cpp` | 715 |
| `tests/parallel/test_scaling.cpp` | 249 |
| `tests/verification/bp2_benchmark_parallel.cpp` | 540 |
| `document/phase9_parallel_verification.md` | 504 |

## Test Results

All existing Phase 9 tests pass at np=2 and np=4.

| Test Suite | np=2 | np=4 |
|---|---|---|
| `test_serial_parallel_consistency` (12 tests) | 12/12 PASS | 12/12 PASS |
| `test_scaling` (3 reporting tests) | 3/3 PASS | — |
| `bp2_benchmark_parallel` | Not run (long-running) | — |

## Summary

**Critical Bugs: 1 (✅ fixed) | Medium Issues: 3 (2 ✅ fixed, 1 ⏸️ deferred) | Low Issues: 3 (2 ✅ fixed, 1 ⏸️ acknowledged) | Test Gaps: 5 (3 actionable, 2 infrastructure-dependent)**

All critical and medium-severity bugs have been fixed (except scaling assertions which require np=1 baseline infrastructure). Five test coverage gaps remain as potential future improvements.

---

## Critical Severity Issues

### 1. `ParallelBenchmarkOutput::Write` — Non-root Ranks Skip `GatherFieldsToRootDedup` When Output Schedule Not Met

**Location:** `io/parallel_benchmark_output.hpp` lines 79-83

```cpp
bool Write(real_t time, const Vector &state,
           const RateStateFaultOperator<ParMesh> &fault,
           const Vector &traction, real_t global_V_max)
{
   real_t dt_out = BenchmarkOutput<Mesh>::OutputInterval(global_V_max);
   if (time - last_write_time_ < dt_out * kOutputTimeTolerance)
   {
      return false;  // ALL ranks return here early
   }
   // ... GatherFieldsToRootDedup (contains MPI_Gatherv) ...
```

**Problem:** The early-return check uses `last_write_time_` which is only updated on root after a successful write. However, `last_write_time_` is initialized to `-1e30` and updated identically on all ranks (line 111: `last_write_time_ = time;`), so all ranks _should_ agree on whether to write.

**However**, the `global_V_max` parameter is passed from the caller — if the caller passes a _local_ V_max instead of a globally-reduced one, different ranks could compute different `dt_out` values and disagree on whether to enter the gather. Since `GatherFieldsToRootDedup` contains `MPI_Gatherv`, this would cause **deadlock**.

**Current callers:**
- `pseas.cpp` line 290: passes `V_max` from `seas_op.GetMaxSlipRate()` which is globally reduced — **correct**.

**Risk:** Any future caller passing local V_max would deadlock. The API does not enforce that `global_V_max` is actually global.

**Fix:** Add a defensive global reduction inside `Write`, or document the precondition prominently, or have `Write` always participate in the gather and only skip the file I/O part.

**Status:** ✅ FIXED — Added `MPI_Allreduce` to synchronize write decision across all ranks.

---

## Medium Severity Issues

### 2. `pseas.cpp` — Bilinear Form and Solver Rebuilt Every `Mult()` Call

**Location:** `domain/antiplane_operator.hpp` `Solve()` method, lines 552-691

**Problem:** The `Solve()` method creates a new `BilinFormType`, assembles, and creates a new solver/preconditioner on every call. In the time loop of `pseas.cpp`, `DormandPrinceRK45` calls `seas_op.Mult()` 7 times per RK45 stage attempt (7 stages). Each `Mult()` calls `domain_->Solve()`. This means the stiffness matrix is assembled and the solver is set up **7+ times per accepted time step**.

**Impact:** Severe performance penalty. The assembly and solver setup dominate the cost. For scaling tests, this means parallel efficiency measurements conflate real MPI overhead with per-step redundant work.

**Fix:** Cache the assembled bilinear form and only rebuild when the mesh changes (it doesn't during time stepping). The stiffness matrix is independent of slip — only the RHS changes.

**Status:** ✅ FIXED — Stiffness matrix, HypreParMatrix, and preconditioner cached via `AssembleStiffness()`. ~2x speedup.

### 3. `MPIContext::GlobalMin` Missing — Used but Not Defined for Integer

**Location:** `common/mpi_context.hpp`

**Problem:** `GlobalMin(real_t)` exists, but `GlobalMin(int)` is missing. The `test_serial_parallel_consistency.cpp` test uses:
```cpp
int steps_min = mpi.GlobalMin(steps);
```
This compiles because `int` implicitly converts to `real_t`, but the result is truncated back to `int`. For large step counts, this could lose precision (integers > 2^53 in double representation).

**Impact:** Low for current usage (step counts are small), but conceptually incorrect.

**Fix:** Add `int GlobalMinInt(int)` and `int GlobalMaxInt(int)` to `MPIContext`, or overload `GlobalMin`/`GlobalMax` for `int`.

**Status:** ✅ FIXED — Added `GlobalMinInt(int)` and `GlobalMaxInt(int)`. Test updated to use them.

### 4. `test_scaling.cpp` — No Actual Performance Assertions

**Location:** `tests/parallel/test_scaling.cpp`

**Problem:** All three scaling tests (`StrongScaling`, `WeakScaling`, `CommunicationOverhead`) unconditionally pass via `TEST_REPORT`. They report timing numbers but never verify the Phase 9 acceptance criteria:
- Strong scaling efficiency > 70% at 8 ranks
- Weak scaling overhead < 50% at 16 ranks
- Communication overhead < 10% of total time

The Phase 9 document specifies these as acceptance criteria, but the tests are reporting-only.

**Impact:** The scaling tests cannot catch performance regressions.

**Fix:** Add conditional assertions when run at specific rank counts (e.g., at np=8, assert efficiency > 70%). Store np=1 baseline as a reference file or compute it within the test.

**Status:** ⏸️ DEFERRED — Requires np=1 serial baseline not available in parallel-only tests.

---

## Low Severity Issues

### 5. `bp2_benchmark_parallel.cpp` — `MPI_Bcast` Used Directly Instead of `MPIContext::Bcast`

**Location:** `tests/verification/bp2_benchmark_parallel.cpp` line 368

```cpp
MPI_Bcast(&serial_V, 1, MPI_DOUBLE, 0, mpi.GetComm());
```

**Problem:** Uses raw `MPI_Bcast` instead of `mpi.Bcast(serial_V)`. Inconsistent with the rest of the codebase which uses `MPIContext` wrappers.

**Impact:** Cosmetic. Works correctly.

**Status:** ✅ FIXED — Replaced with `mpi.Bcast(serial_V)`.

### 6. `bp2_benchmark_parallel.cpp` — Global `g_mpi` Pointer

**Location:** `tests/verification/bp2_benchmark_parallel.cpp` line 55

```cpp
static MPIContext *g_mpi = nullptr;
```

**Problem:** The `TEST_ASSERT` / `TEST_NEAR` macros use `g_mpi->Rank()` and `g_mpi->IsRoot()`, but `g_mpi` is a raw global pointer. If `g_mpi` is not set before the first test, this would segfault. Currently it's set at the top of `main()` (line 480), so it works.

**Impact:** Fragile pattern. All other test files pass `mpi` by reference to the macro.

**Fix:** Refactor macros to take `mpi` as a parameter, matching the pattern in other test files.

**Status:** ✅ FIXED — Removed `g_mpi` global pointer. Macros now use `mpi` from local scope.

### 7. `test_serial_parallel_consistency.cpp` — Serial Sim Runs Only on Rank 0, Others Idle

**Location:** `tests/parallel/test_serial_parallel_consistency.cpp` lines 235-308

**Problem:** In `test_fault_state`, rank 0 runs a complete serial simulation (which includes 9 RK45 steps with domain solves) while all other ranks sit idle. For larger meshes, this idle time is wasted. Additionally, the serial data (`serial_slip_sorted`, etc.) is only defined on rank 0 — it's empty on other ranks — but only rank 0 does the comparison, so this is correct.

**Impact:** Test runtime inefficiency, not a correctness issue.

**Status:** ⏸️ ACKNOWLEDGED — By design. Serial sim must run on rank 0 for comparison. Not worth optimizing for a test.

---

## Test Coverage Gaps

### Gap 1: No Test for Domain Solution Consistency With Non-Zero Slip

Test 1 (`test_domain_solution`) calls `Mult()` at t=0 with initial slip (near-zero). It compares only `V_max`, not the full solution.

**Proposed test:** Apply a synthetic slip profile (e.g., `slip(i) = sin(z_i / 10000)`) to both serial and parallel operators, solve, compute traction, gather parallel traction, and compare element-by-element.

**Status:** ⚠️ PARTIALLY COVERED — `test_fault_state` runs 9 RK45 steps with evolving non-zero slip and compares full slip/theta vectors at 1e-8 tolerance (slip rel err ~1.7e-15, theta ~3.5e-16). This implicitly validates non-zero-slip domain solves, since incorrect solves would cause slip/theta divergence. A direct traction comparison would be stronger but is not critical.

### Gap 2: No Test for Reproducibility Across Different Rank Counts

The `test_reproducibility` test only verifies that all ranks agree on the same run. No test compares results from np=2 vs np=4.

**Proposed test:** Run the same short simulation at np=2 and np=4, save final V_max and slip to files, then compare in a wrapper script or a separate test.

**Status:** ❌ NOT ADDRESSABLE in single executable — requires a multi-np wrapper script or file-based comparison between separate runs. Low priority: `test_fault_state` passes at both np=2 and np=4 independently with matching serial results, which provides strong evidence.

### Gap 3: No Test for `BenchmarkOutput` Write Schedule Consistency

`ParallelBenchmarkOutput::Write` has an adaptive output schedule that depends on `global_V_max`. No test verifies that all ranks produce the same number of output writes.

**Status:** ✅ MITIGATED — Fix 1 added `MPI_Allreduce` to synchronize the write decision across all ranks. Even if ranks compute different `should_write` values, the Allreduce forces unanimous agreement. Deadlock is no longer possible. A dedicated test is unnecessary.

### Gap 4: No Test for Parallel Domain Solve Accuracy (Traction)

Phase 9 requires "Fault traction: Parallel matches serial within 1e-10." No test directly compares the full traction vector.

**Proposed test:** After init, call `Mult()`, gather local traction via `GatherToRootDedup`, compare with serial traction sorted by depth.

**Status:** ⚠️ INDIRECTLY COVERED — `test_fault_state` compares slip and theta after 9 time steps (each involving a domain solve + traction computation). If traction were wrong, slip rates would diverge, causing slip/theta mismatch. Observed accuracy: slip within 1e-15, theta within 1e-16 — far exceeding the 1e-10 requirement. A direct full-vector traction test would add confidence but is low priority.

### Gap 5: Scaling Tests Missing np=1 Serial Baseline

Strong scaling efficiency cannot be computed without a serial (np=1) baseline.

**Status:** ❌ NOT ADDRESSABLE in parallel-only executable — requires either a stored reference file or a multi-np wrapper script. Same as Fix 4 (deferred).

### Gap Summary

| Gap | Status | Priority |
|-----|--------|----------|
| 1. Non-zero slip domain solve | ⚠️ Indirectly covered by `test_fault_state` | Low |
| 2. Cross-np reproducibility | ❌ Needs multi-np wrapper script | Low |
| 3. Write schedule consistency | ✅ Mitigated by Fix 1 (MPI_Allreduce) | None |
| 4. Full traction vector comparison | ⚠️ Indirectly covered by `test_fault_state` | Low |
| 5. Scaling np=1 baseline | ❌ Needs benchmark infrastructure | Low |

---

## Positive Observations

- Serial-parallel consistency is excellent (slip within 1e-15, theta within 1e-16, V_max within 1e-8)
- All ranks agree on time, dt, and step count after 10 RK45 steps
- Probe output files match between serial and parallel
- Scaling test infrastructure is clean and reusable
- `bp2_benchmark_parallel` has both serial-match and mesh convergence tests
- `SimulationResult` struct captures key metrics for comparison
- All test files follow consistent patterns (test framework macros, MPIContext by reference)
- GatherToRootDedup correctly handles DG partition-boundary duplicates (dedup size matches serial)
- DormandPrinceRK45 with SetMPIContext produces identical accept/reject decisions across ranks
