# Phase 6 Quick Fixes

**Version 2.0** - 2026-02-15 (all fixes completed)

This document provides actionable fixes for discrepancies identified in the Phase 6 debug report.

---

## Status Summary

| Fix | Priority | Status |
|-----|----------|--------|
| Fix 1: `AllRanksAgree` generic template uses wrong MPI type in doc | LOW | COMPLETED (doc already updated) |
| Fix 2: `BroadcastFromRoot` uses `MPI_BYTE` in doc | LOW | COMPLETED (doc already updated) |
| Fix 3: Gather/Scatter API uses `mfem::Vector` in doc, `std::vector` in impl | LOW | COMPLETED (doc already updated) |
| Fix 4: Add `GetComm()` to `MPIContext` in doc | LOW | COMPLETED (doc already updated) |
| Fix 5: Add `Bcast` methods to `MPIContext` in doc | LOW | COMPLETED (doc already updated) |
| Fix 6: `DomainOperator::GetComm()` uses `SEAS_USE_MPI` instead of `MFEM_USE_MPI` | LOW | COMPLETED (code already uses `MFEM_USE_MPI`) |
| Fix 7: Remove or update unused `SEASSolver`/`SEASPreconditioner` aliases | LOW | COMPLETED (removed from doc; never existed in impl) |

---

## Fix 1: `AllRanksAgree` Generic Template Uses Wrong MPI Type (LOW)

**File:** `document/phase6_parallel_infrastructure.md`

**Current (lines 259-264):**
```cpp
template <typename T>
bool AllRanksAgree(T local_value, MPI_Comm comm = MPI_COMM_WORLD) {
    T min_val, max_val;
    MPI_Allreduce(&local_value, &min_val, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_value, &max_val, 1, MPI_DOUBLE, MPI_MAX, comm);
    return (min_val == max_val);
}
```

Uses `MPI_DOUBLE` for all `T`, which is wrong for `int`.

**Fix:** Update doc to show explicit specializations matching the implementation:
```cpp
template <typename T>
bool AllRanksAgree(T local_val, MPI_Comm comm = MPI_COMM_WORLD);

template <>
bool AllRanksAgree<int>(int local_val, MPI_Comm comm) {
    int global_min, global_max;
    MPI_Allreduce(&local_val, &global_min, 1, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(&local_val, &global_max, 1, MPI_INT, MPI_MAX, comm);
    return global_min == global_max;
}

template <>
bool AllRanksAgree<real_t>(real_t local_val, MPI_Comm comm) {
    real_t global_min, global_max;
    MPI_Allreduce(&local_val, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_val, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    return global_min == global_max;
}
```

---

## Fix 2: `BroadcastFromRoot` Uses `MPI_BYTE` in Doc (LOW)

**File:** `document/phase6_parallel_infrastructure.md`

**Current (lines 253-255):**
```cpp
template <typename T>
void BroadcastFromRoot(T &value, MPI_Comm comm = MPI_COMM_WORLD) {
    MPI_Bcast(&value, sizeof(T), MPI_BYTE, 0, comm);
}
```

`MPI_BYTE` with `sizeof(T)` bypasses MPI type checking.

**Fix:** Update doc to show explicit specializations:
```cpp
template <typename T>
void BroadcastFromRoot(T &value, MPI_Comm comm = MPI_COMM_WORLD);

template <>
void BroadcastFromRoot<int>(int &value, MPI_Comm comm) {
    MPI_Bcast(&value, 1, MPI_INT, 0, comm);
}

template <>
void BroadcastFromRoot<real_t>(real_t &value, MPI_Comm comm) {
    MPI_Bcast(&value, 1, MPI_DOUBLE, 0, comm);
}
```

---

## Fix 3: Gather/Scatter API Container Type Mismatch (LOW)

**File:** `document/phase6_parallel_infrastructure.md`

**Current (lines 238-250):** Doc shows API using `mfem::Vector`, `mfem::Array<int>`, and `MPI_Comm` as last parameter:
```cpp
void GatherVectorToRoot(const Vector &local_data,
                        Vector &global_data,
                        const Array<int> &recv_counts,
                        const Array<int> &displacements,
                        MPI_Comm comm = MPI_COMM_WORLD);
```

**Implementation (`parallel_utils.hpp`):** Uses `std::vector<real_t>` and computes counts/displacements internally:
```cpp
void GatherVectorToRoot(const std::vector<real_t> &local,
                        std::vector<real_t> &global,
                        MPI_Comm comm = MPI_COMM_WORLD);
```

**Fix:** Update doc to match the simpler implementation API that handles counts/displacements internally. The implementation's API is superior — callers don't need to precompute displacements.

---

## Fix 4: Add `GetComm()` to `MPIContext` in Doc (LOW)

**File:** `document/phase6_parallel_infrastructure.md`

The doc's `MPIContext` class (lines 137-213) is missing `GetComm()`.

**Fix:** Add to the parallel branch of the `MPIContext` class definition in the doc:
```cpp
#ifdef SEAS_USE_MPI
    // ... existing constructor/destructor ...

    /// Get MPI communicator
    MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#else
    // ... no GetComm in serial ...
#endif
```

---

## Fix 5: Add `Bcast` Methods to `MPIContext` in Doc (LOW)

**File:** `document/phase6_parallel_infrastructure.md`

The doc's `MPIContext` only shows `GlobalMax`, `GlobalMin`, `GlobalSum`, `GlobalSumInt`, `Barrier`. The implementation adds two `Bcast` methods.

**Fix:** Add to the `MPIContext` class definition in the doc:
```cpp
    /// Broadcast a scalar from root to all processes
    void Bcast(real_t &value) const {
#ifdef SEAS_USE_MPI
        MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }

    /// Broadcast a vector from root to all processes
    void Bcast(Vector &vec) const {
#ifdef SEAS_USE_MPI
        int n = vec.Size();
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!IsRoot()) { vec.SetSize(n); }
        MPI_Bcast(vec.GetData(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }
```

---

## Fix 6: `DomainOperator::GetComm()` Guard Inconsistency (LOW)

**File:** `domain/domain_operator.hpp`

**Current (lines 110-113):**
```cpp
#ifdef SEAS_USE_MPI
   virtual MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#endif
```

But the convenience alias (lines 118-120) uses `MFEM_USE_MPI`:
```cpp
#ifdef MFEM_USE_MPI
using ParallelDomainOperator = DomainOperator<ParMesh>;
#endif
```

`MPI_Comm` comes from MFEM's MPI headers, so availability depends on `MFEM_USE_MPI`, not `SEAS_USE_MPI`.

**Fix:** Change `SEAS_USE_MPI` to `MFEM_USE_MPI` on line 110:
```cpp
#ifdef MFEM_USE_MPI
   virtual MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#endif
```

**Note:** In practice this never causes a build failure because `SEAS_USE_MPI` is only defined when `MFEM_USE_MPI` is also set.

---

## Fix 7: Remove or Update Unused Solver/Preconditioner Aliases (LOW)

**File:** `common/seas_types.hpp`

**Current (lines 38-39, 74-75):**
```cpp
// Parallel
using SEASSolver = HyprePCG;
using SEASPreconditioner = HypreBoomerAMG;

// Serial
using SEASSolver = CGSolver;
using SEASPreconditioner = GSSmoother;
```

These are never used anywhere in the codebase. The actual domain operator creates `CGSolver` + `GSSmoother` (serial) or `CGSolver` + `HypreSmoother` (parallel) directly. The parallel aliases are also misleading — `HyprePCG`/`HypreBoomerAMG` were abandoned because AMG fails on singular (all-Neumann) DG systems.

**Fix (Option A — remove):** Delete the aliases entirely since they're unused and misleading.

**Fix (Option B — update):** Change to match actual usage:
```cpp
// Parallel
using SEASSolver = CGSolver;
using SEASPreconditioner = HypreSmoother;

// Serial
using SEASSolver = CGSolver;
using SEASPreconditioner = GSSmoother;
```

**Recommendation:** Option A (remove). Dead code that doesn't reflect reality should be deleted.
