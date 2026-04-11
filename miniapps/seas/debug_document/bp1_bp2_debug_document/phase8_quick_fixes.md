# Phase 8 Quick Fixes

## Fix 1: DormandPrinceRK45 — Global Error Norm for Parallel Consistency [CRITICAL] — FIXED

**File:** `solver/time_stepper.hpp`

**Problem:** Each rank computed a local error norm and made independent accept/reject decisions. Ranks could diverge, causing deadlock in subsequent `op.Mult()` calls (which contain MPI collectives).

**Resolution:** Added `MPIContext *mpi_ctx_` member with `SetMPIContext()` setter. After computing local `err_norm`, `mpi_ctx_->GlobalMax(err_norm)` ensures all ranks agree. Diagnostic output (`[RK45 dt_min]`, `[RK45 STUCK]`) is now guarded with `(!mpi_ctx_ || mpi_ctx_->IsRoot())`. `pseas.cpp` calls `ode_solver.SetMPIContext(&mpi)`.

**Status:** FIXED.

---

## Fix 2: Deduplicate Gathered Fault Data [MEDIUM] — FIXED

**File:** `fault/fault_geometry.hpp`, `io/parallel_benchmark_output.hpp`

**Problem:** At partition boundaries, DG fault DOFs are counted by both ranks. `GatherToRoot` produced a vector with duplicates. The `ProbeInterpolator` could bracket a probe depth between a duplicate pair with inconsistent field values.

**Resolution:** Added `GatherToRootDedup()` and `GatherFieldsToRootDedup()` methods with `DeduplicateByDepth()` / `DeduplicateMultipleByDepth()` static helpers. These sort by depth (descending, surface first) and average duplicate entries within a 1.0 m tolerance. `ParallelBenchmarkOutput::Write()` uses `GatherFieldsToRootDedup()`. `pseas.cpp` uses `GatherToRootDedup()` for initial fault depth setup.

**Status:** FIXED.

---

## ~~Fix 3: Consistent MPI Guard Macros~~ [INVESTIGATED — NO FIX NEEDED]

**File:** `fault/fault_geometry.hpp` line 108

**Initial concern:** Uses `SEAS_USE_MPI` while `antiplane_operator.hpp` uses `MFEM_USE_MPI`.

**Resolution:** `SEAS_USE_MPI` is correct. It guards `MPIContext::GetComm()` which is only available in parallel SEAS targets. `MFEM_USE_MPI` is always defined when MFEM was built with MPI (even for serial SEAS executables), so using it here would cause compile errors in serial targets like `seas.cpp`. The two macros serve different purposes:
- `MFEM_USE_MPI` → guards MFEM parallel types (`ParMesh`, `ParFiniteElementSpace`)
- `SEAS_USE_MPI` → guards SEAS parallel utilities (`MPIContext::GetComm()`, MPI calls)

**Status:** No action needed.

---

## Fix Priority

| Fix | Severity | Status |
|-----|----------|--------|
| 1. RK45 global error norm | Critical | FIXED |
| 2. Deduplicate gathered data | Medium | FIXED |
| ~~3. MPI guard consistency~~ | ~~Low~~ | Investigated, correct as-is |
