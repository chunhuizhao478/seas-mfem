# Phase 8 Debug Report: Parallel Fault and I/O

## Files Reviewed

| File | Lines |
|------|-------|
| `fault/fault_geometry.hpp` | ~470 |
| `fault/rate_state_fault.hpp` | 416 |
| `io/parallel_benchmark_output.hpp` | 147 |
| `io/benchmark_output.hpp` | 267 |
| `io/probe_output.hpp` | 260 |
| `solver/seas_operator.hpp` | 223 |
| `solver/time_stepper.hpp` | ~450 |
| `pseas.cpp` | 314 |
| `tests/parallel/test_parallel_fault.cpp` | 1011 |
| `document/phase8_parallel_fault_io.md` | 497 |

## Test Results

All Phase 8 tests pass at np=2 and np=4.

| Test Suite | np=2 | np=4 |
|---|---|---|
| `test_parallel_fault` (27 tests) | 27/27 PASS | 27/27 PASS |

## Summary

**Critical Bugs: 0 (1 FIXED) | Medium Issues: 0 (2 FIXED) | Low Issues: 2 | Investigated & OK: 1 | Test Gaps: 0 (5 FILLED)**

All critical and medium bugs have been fixed and verified by new tests. All previously identified test gaps have been filled. The remaining low issues are cosmetic/documentation only.

---

## Critical Severity Issues

### 1. DormandPrinceRK45 Error Norm Is Local — Ranks Can Disagree on Accept/Reject — FIXED

**Location:** `solver/time_stepper.hpp` lines 294-300

**Problem:** `DormandPrinceRK45::Step()` computed a local L-infinity error norm. Each rank made independent accept/reject decisions, causing potential divergence and deadlock.

**Resolution:** Added `MPIContext *mpi_ctx_` member with `SetMPIContext()` setter. After computing local `err_norm`, `mpi_ctx_->GlobalMax(err_norm)` ensures all ranks agree. `pseas.cpp` calls `ode_solver.SetMPIContext(&mpi)`.

**Verified by:** Test 6 (`test_rk45_sync`) and Test 10 (`test_dt_agreement`) — both confirm all ranks agree on `t` and `dt` after 5 RK45 steps at np=2 and np=4.

**Status:** FIXED.

---

## Medium Severity Issues

### 2. Gathered Fault Data Contains Duplicate DOFs at Partition Boundaries — FIXED

**Location:** `fault/fault_geometry.hpp` `GatherToRootDedup()` and `GatherFieldsToRootDedup()`

**Problem:** In parallel DG, shared fault faces at partition boundaries produced duplicate DOFs in gathered data. At np=4: 52 gathered DOFs vs 50 serial. `ProbeInterpolator` could pick inconsistent field values at duplicates.

**Resolution:** Added `GatherToRootDedup()` and `GatherFieldsToRootDedup()` methods with `DeduplicateByDepth()` / `DeduplicateMultipleByDepth()` static helpers. These sort by depth and average duplicate entries within a 1.0 m tolerance. `ParallelBenchmarkOutput::Write()` now uses `GatherFieldsToRootDedup()`. `pseas.cpp` uses `GatherToRootDedup()` for initial fault depth setup.

**Verified by:** Test 7 (`test_compute_rhs_consistency`) confirms dedup size matches serial (50==50). Test 8 (`test_dedup_probe_interpolation`) confirms interpolated probe values match serial within 0.01.

**Status:** FIXED.

### 3. Diagnostic Output from All Ranks in DormandPrinceRK45 — FIXED

**Location:** `solver/time_stepper.hpp` lines 333 and 364-365

**Problem:** `[RK45 dt_min]` and `[RK45 STUCK]` diagnostics were printed by all ranks.

**Resolution:** Guarded with `(!mpi_ctx_ || mpi_ctx_->IsRoot())` so only root prints.

**Status:** FIXED.

---

## Low Severity Issues

### 4. `pseas.cpp` Uses `std::system()` to Delete Old Output Files

**Location:** `pseas.cpp` lines 196-198

```cpp
std::string rm_cmd = "rm -f " + output_prefix + "_z*.txt";
std::system(rm_cmd.c_str());
```

**Problem:** Uses shell command for file deletion. Not portable; `std::system()` return value ignored.

**Impact:** Low — prefix is a simple string in practice.

### 5. Doc vs Implementation: No Separate `parallel_fault_data.hpp` or `parallel_time_stepper.hpp`

**Document (tasks 8.1, 8.3):** Specifies these as separate deliverables.

**Implementation:** Parallel fault data is integrated into `FaultGeometry<ParMesh>`. Parallel time stepping is handled by `DormandPrinceRK45` with `SetMPIContext()`.

**Status:** Integration into existing classes is cleaner than separate files. No action needed.

### 6. `FaultGeometry::GatherToRoot` Uses `SEAS_USE_MPI` Guard, Not `MFEM_USE_MPI` — INVESTIGATED, CORRECT

**Location:** `fault/fault_geometry.hpp` line 108

**Investigation result:** `SEAS_USE_MPI` is the correct guard. It controls `MPIContext::GetComm()` availability, which is SEAS-specific. `MFEM_USE_MPI` is always defined when MFEM was built with MPI — even for serial SEAS targets. Using `MFEM_USE_MPI` would cause compile errors in serial targets that include `FaultGeometry<Mesh>`.

The two macros serve different purposes:
- `MFEM_USE_MPI` → guards MFEM parallel types (`ParMesh`, `ParFiniteElementSpace`)
- `SEAS_USE_MPI` → guards SEAS parallel utilities (`MPIContext::GetComm()`, MPI calls)

**Status:** No action needed. The existing guard is correct.

---

## Test Coverage — ALL GAPS FILLED

### Test 6: RK45 Time Step Synchronization — IMPLEMENTED, PASSING

Takes 5 accepted RK45 steps with `SetMPIContext`, verifies `t` and `dt` agree across all ranks after each step. Confirms Fix #1.

### Test 7: ComputeRHS Parallel Consistency — IMPLEMENTED, PASSING

Calls `Mult()` on serial and parallel operators, gathers parallel rates via `GatherToRootDedup`, compares max(V) and dedup size to serial. Confirms rates match and dedup produces correct count.

### Test 8: Deduplicated Probe Interpolation vs Serial — IMPLEMENTED, PASSING

Creates synthetic field `sin(z/10000)` on local depths, gathers with dedup, interpolates on root at 4 probe depths. Compares to serial interpolation. Confirms Fix #2.

### Test 9: Parallel Benchmark Output File Content — IMPLEMENTED, PASSING

Runs serial and parallel `ForceWrite` at t=0, reads output files on root, compares time and `log10(V)` values. Confirms `ParallelBenchmarkOutput` produces correct file content.

### Test 10: Time Stepper dt Agreement Across Ranks — IMPLEMENTED, PASSING

Takes 5 accepted RK45 steps, checks `GlobalMin(dt) == GlobalMax(dt)` after each step. Confirms Fix #1 from a dt-specific angle.

---

## Positive Observations

- Clean template-based parallel dispatch (`FaultGeometry<ParMesh>`, `RateStateFaultOperator<ParMesh>`)
- Correct `MPI_Gatherv` pattern with precomputed `recv_counts_`/`recv_displs_`
- `GetGlobalMaxSlipRate()` correctly reduces across ranks
- `SEASQuasiDynamicOperator::SetInitialCondition()` correctly reduces stress equilibrium error across ranks
- Adaptive output schedule uses globally-reduced `V_max`
- `ParallelBenchmarkOutput` correctly uses `GatherFieldsToRootDedup` to handle partition-boundary duplicates
- FSAL optimization in RK45 is correctly implemented
- `DormandPrinceRK45` correctly synchronizes error norm across ranks via `GlobalMax`
- `DeduplicateByDepth` / `DeduplicateMultipleByDepth` sort by depth and merge duplicates with averaging
- Comprehensive test suite (27 tests) covering DOF distribution, gather, V_max reduction, serial-parallel consistency, probe interpolation, RK45 sync, ComputeRHS consistency, dedup interpolation, benchmark output, and dt agreement
