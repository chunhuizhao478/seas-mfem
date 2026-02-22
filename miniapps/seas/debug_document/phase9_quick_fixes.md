# Phase 9 Quick Fixes

## Fix 1: ParallelBenchmarkOutput::Write — Potential MPI Deadlock on Output Schedule Disagreement [CRITICAL]

**File:** `io/parallel_benchmark_output.hpp`

**Problem:** `Write()` returns early (line 83) before reaching the `MPI_Gatherv` inside `GatherFieldsToRootDedup()` if the output schedule check fails. The check depends on `global_V_max` which is passed by the caller. If a future caller passes a local (not globally-reduced) V_max, different ranks compute different `dt_out`, causing some ranks to skip the gather while others enter it — deadlock.

Currently safe because `pseas.cpp` passes a globally-reduced V_max, but the API is fragile.

**Proposed fix (option A — defensive reduction):**
```cpp
bool Write(real_t time, const Vector &state,
           const RateStateFaultOperator<ParMesh> &fault,
           const Vector &traction, real_t global_V_max)
{
   // Ensure all ranks agree on the output decision
   real_t dt_out = BenchmarkOutput<Mesh>::OutputInterval(global_V_max);
   int should_write = (time - last_write_time_ >= dt_out * kOutputTimeTolerance) ? 1 : 0;
   // Synchronize decision
   int global_should_write = should_write;
   MPI_Allreduce(MPI_IN_PLACE, &global_should_write, 1, MPI_INT, MPI_MAX,
                 mpi_ctx_.GetComm());
   if (!global_should_write) { return false; }
   // ... proceed with gather ...
```

**Proposed fix (option B — documentation):** Add a prominent comment/assert:
```cpp
// PRECONDITION: global_V_max must be globally reduced (same on all ranks)
// Otherwise ranks may disagree on output schedule, causing MPI_Gatherv deadlock.
```

**Status:** ✅ FIXED. Added `MPI_Allreduce(MPI_IN_PLACE, &should_write, 1, MPI_INT, MPI_MAX, ...)` to synchronize the write decision across all ranks before entering the gather.

---

## Fix 2: Stiffness Matrix Rebuilt Every Mult() — Major Performance Issue [MEDIUM]

**File:** `domain/antiplane_operator.hpp`, `Solve()` method

**Problem:** `Solve()` assembles a new `BilinFormType`, creates a new preconditioner, and sets up the solver on every call. With RK45 calling `Mult()` 7 times per step attempt, this dominates the cost. The stiffness matrix doesn't change between calls — only the RHS (slip BC) changes.

**Proposed fix:** Cache the assembled stiffness matrix and solver/preconditioner. Only rebuild when first called or when mesh changes. Move assembly from `Solve()` to a new `AssembleStiffness()` called once during construction or lazily on first `Solve()`.

```cpp
// In AntiplaneDomainOperator:
mutable bool stiffness_assembled_ = false;
mutable std::unique_ptr<BilinFormType> a_;  // Cached bilinear form
// ...
void Solve(...) {
   if (!stiffness_assembled_) {
      AssembleStiffness();
      stiffness_assembled_ = true;
   }
   // Only rebuild RHS and solve
}
```

**Status:** ✅ FIXED. Extracted `AssembleStiffness()` method; stiffness matrix, HypreParMatrix, and preconditioner are cached. Only RHS is rebuilt per `Solve()`. ~2x speedup observed (np=2: 16s→7s for 50 steps).

---

## Fix 3: MPIContext Missing Integer GlobalMin/GlobalMax [MEDIUM]

**File:** `common/mpi_context.hpp`

**Problem:** `GlobalMin(real_t)` and `GlobalMax(real_t)` exist, but integer variants are missing. `test_serial_parallel_consistency.cpp` calls `mpi.GlobalMin(steps)` which implicitly converts `int → real_t → int`, potentially losing precision for large values.

**Proposed fix:**
```cpp
int GlobalMinInt(int local_val) const
{
#ifdef SEAS_USE_MPI
   int global_val;
   MPI_Allreduce(&local_val, &global_val, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   return global_val;
#else
   return local_val;
#endif
}

int GlobalMaxInt(int local_val) const
{
#ifdef SEAS_USE_MPI
   int global_val;
   MPI_Allreduce(&local_val, &global_val, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   return global_val;
#else
   return local_val;
#endif
}
```

**Status:** ✅ FIXED. Added `GlobalMinInt(int)` and `GlobalMaxInt(int)` to `MPIContext`. Updated `test_serial_parallel_consistency.cpp` to use them.

---

## Fix 4: Scaling Tests Have No Performance Assertions [MEDIUM]

**File:** `tests/parallel/test_scaling.cpp`

**Problem:** All tests unconditionally pass (`TEST_REPORT`). Phase 9 acceptance criteria require:
- Strong scaling efficiency > 70% at 8 ranks
- Weak scaling overhead < 50%
- Communication overhead < 10%

None of these are checked.

**Proposed fix:** Store a serial baseline time (either from np=1 run or from a reference file). Add conditional assertions:
```cpp
// In test_strong_scaling:
if (mpi.Size() == 8 && serial_baseline > 0)
{
   double efficiency = (serial_baseline / wall_time) / mpi.Size();
   TEST_ASSERT(efficiency > 0.70,
               "Strong scaling efficiency > 70% at np=8");
}
```

Note: Real performance assertions require a serial baseline. For now, the tests are useful as reporting tools. Full assertions can be added when benchmark infrastructure is available.

**Status:** ⏸️ DEFERRED. Requires a serial (np=1) baseline to compute efficiency. The tests remain reporting-only, which is appropriate — real performance assertions need dedicated benchmark infrastructure or a multi-np wrapper script.

---

## Fix Priority

| Fix | Severity | Status |
|-----|----------|--------|
| 1. ParallelBenchmarkOutput deadlock risk | Critical (latent) | ✅ Fixed |
| 2. Stiffness matrix caching | Medium (performance) | ✅ Fixed |
| 3. MPIContext integer min/max | Medium | ✅ Fixed |
| 4. Scaling test assertions | Medium | ⏸️ Deferred (needs np=1 baseline) |
