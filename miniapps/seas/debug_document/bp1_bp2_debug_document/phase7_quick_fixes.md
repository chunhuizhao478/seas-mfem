# Phase 7 Quick Fixes

## Status Summary

| Fix | Priority | Status |
|-----|----------|--------|
| Fix 1: Add BR2 parallel domain tests | HIGH | COMPLETE |
| Fix 2: Add BR2 parallel MMS test | HIGH | COMPLETE |
| Fix 3: Strengthen serial-parallel traction test | HIGH | COMPLETE |
| Fix 4: Add fault depth consistency test | MEDIUM | COMPLETE |

---

## Fix 1: Add BR2 Parallel Domain Tests

**File:** `tests/parallel/test_parallel_domain.cpp`

Add two new tests mirroring `test_parallel_solve` and `test_parallel_uniform_slip` but with `DGMethod::BR2`:

```cpp
/// Test: Parallel BR2 solve with zero slip converges
bool test_parallel_solve_br2(MPIContext &ctx)
{
   // Same as test_parallel_solve but with DGMethod::BR2
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf, DGMethod::BR2);
   // ... zero slip, verify solution near zero
}

/// Test: Parallel BR2 uniform slip gives max_u ~ 0.5
bool test_parallel_uniform_slip_br2(MPIContext &ctx)
{
   // Same as test_parallel_uniform_slip but with DGMethod::BR2
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, bp2.mu(), bp2.Vp, Wf, DGMethod::BR2);
   // ... uniform slip=1, verify max_u in [0.35, 0.65]
}
```

---

## Fix 2: Add BR2 Parallel MMS Test

**File:** `tests/parallel/mms_antiplane_parallel.cpp`

Run the MMS convergence test twice: once with IP, once with BR2. Both should achieve ~2nd order convergence.

```cpp
// After the existing IP test loop:
// Run again with BR2
for (int level = 0; level < 4; level++)
{
   AntiplaneDomainOperator<ParMesh> op(pmesh, 1, mu, Vp, Wf, DGMethod::BR2);
   // ... same MMS verification
}
```

---

## Fix 3: Strengthen Serial-Parallel Traction Test

**File:** `tests/parallel/test_parallel_domain.cpp`

Replace or supplement `test_serial_parallel_consistency` with a non-uniform slip test:

```cpp
/// Test: Serial and parallel traction match with non-uniform slip
bool test_traction_nonuniform_slip(MPIContext &ctx)
{
   // Use slip(i) = sin(pi * depth(i) / Lz) for depth-varying slip
   // Compare gathered parallel traction against serial traction pointwise
   // Tolerance: relative error < 1e-6
}
```

Key steps:
1. Compute serial traction with depth-dependent slip
2. Broadcast serial traction to all ranks
3. Compute parallel traction, gather to root
4. Sort both by depth and compare pointwise

---

## Fix 4: Add Fault Depth Consistency Test

**File:** `tests/parallel/test_parallel_domain.cpp`

```cpp
/// Test: Gathered parallel fault depths match serial fault depths
bool test_fault_depth_consistency(MPIContext &ctx)
{
   // 1. Get serial fault depths (sorted)
   // 2. Get parallel fault depths on each rank, gather to root
   // 3. Sort gathered depths
   // 4. Compare sorted sets within tolerance
}
```
