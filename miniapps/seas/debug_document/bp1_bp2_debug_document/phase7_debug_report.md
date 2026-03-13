# Phase 7 Debug Report: Parallel Domain Operator

## Files Reviewed

| File | Lines |
|------|-------|
| `domain/antiplane_operator.hpp` | 1776 |
| `tests/parallel/test_parallel_domain.cpp` | 368 |
| `tests/parallel/mms_antiplane_parallel.cpp` | 119 |
| `document/phase7_parallel_domain.md` | 408 |

## Test Results

All existing tests pass at np=2 and np=4.

| Test Suite | np=2 | np=4 |
|---|---|---|
| `test_parallel_domain` (10 tests) | 10/10 PASS | 10/10 PASS |
| `mms_antiplane_parallel` (IP + BR2) | PASS (IP ~1.95, BR2 ~1.95) | PASS (IP ~1.95, BR2 ~1.95) |

## Summary

**Critical Bugs: 0 | Medium Issues: 2 (ALL FIXED) | Low Issues: 3 | Test Gaps: 4 (ALL FIXED)**

The parallel domain operator is well-implemented. All 10 domain tests and both MMS convergence tests (IP + BR2) pass at np=2 and np=4.

---

## Medium Severity Issues

### 1. No Parallel Tests for BR2 Method — FIXED

**Location:** `tests/parallel/test_parallel_domain.cpp`, `tests/parallel/mms_antiplane_parallel.cpp`

**Resolution:** Added `test_parallel_solve_br2` and `test_parallel_uniform_slip_br2` to domain tests. Added BR2 convergence loop to MMS test. All pass at np=2 and np=4. BR2 MMS convergence rate ~1.95 matches IP.

### 2. Serial-Parallel Traction Comparison Is Weak — FIXED

**Location:** `test_traction_nonuniform_slip` (line 429-559)

**Resolution:** Added `test_traction_nonuniform_slip` which uses depth-varying slip `sin(π(z+Lz)/Lz)`, gathers parallel traction to root, matches by depth, and compares pointwise. Results: max_rel_err=1.76e-11 (np=2), 7.84e-12 (np=4) — essentially machine precision.

---

## Low Severity Issues

### 3. Doc vs Implementation: Solver Type Mismatch (INFO)

**Document (line 56-57):** Shows `HyprePCG` / `HypreBoomerAMG` for parallel.

**Implementation (line 466-478):** Uses `CGSolver(comm)` / `HypreSmoother`. This is correct — AMG fails on singular all-Neumann DG systems. The document table was already updated externally.

**Status:** No action needed.

### 4. Doc vs Implementation: No Separate `antiplane_operator_impl.hpp` (INFO)

**Document (line 104):** Shows separate `antiplane_operator_impl.hpp` for serial/parallel specializations.

**Implementation:** Everything is in `antiplane_operator.hpp` using `if constexpr` dispatch. This is a cleaner approach.

**Status:** No action needed.

### 5. BR2 One-Sided Lifting Approximation May Lose Accuracy

**Location:** `antiplane_operator.hpp` lines 1696-1713

```cpp
// Double the contribution (one-sided approximation)
sum += n_j * (2.0 * contrib1);
```

For shared faces, only Elem1's lifting is computed and doubled. This approximation assumes the two elements are geometrically similar. On non-uniform meshes or at partition boundaries where elements differ in size, this may introduce O(h) errors in the BR2 penalty.

**Severity:** Low for current use (uniform meshes), but worth verifying with a non-uniform mesh test.

---

## Test Coverage Gaps

### Gap 1: BR2 Parallel MMS Convergence — FIXED

Added BR2 convergence loop to `mms_antiplane_parallel.cpp`. Rate ~1.95 at np=2 and np=4.

### Gap 2: Non-Uniform Slip in Parallel — FIXED

Added `test_traction_nonuniform_slip` with depth-varying slip. Pointwise traction matches serial within 1e-11.

### Gap 3: Fault Depth Consistency Across Ranks — FIXED

Added `test_fault_depth_consistency`. Gathered parallel depths (deduplicated) match serial depths exactly.

### Gap 4: Traction Computation on Shared Faces — FIXED

Covered by `test_traction_nonuniform_slip` which exercises shared-face traction with non-trivial slip. At np=4 all fault faces are shared (0 interior, 16 shared), so the shared-face traction path is fully exercised.

---

## Recommended Additional Tests

| Test | Purpose | Priority | Status |
|------|---------|----------|--------|
| BR2 parallel MMS convergence | Verify BR2 convergence in parallel | High | DONE |
| BR2 parallel solve (zero/uniform slip) | Verify BR2 shared face assembly | High | DONE |
| Non-uniform slip serial-parallel comparison | Verify slip indexing + traction | High | DONE |
| Fault depth gathering across ranks | Verify depth consistency | Medium | DONE |

---

## Positive Observations

- Clean `if constexpr` dispatch for serial/parallel code paths
- Correct DG property exploited: TrueVSize == VSize (no shared DOFs)
- Proper face-neighbor data exchange before shared face operations
- Fault DOF accounting test correctly handles double-counted shared faces
- MMS convergence rate ~1.95 matches expected 2nd order
- Solver handles singular all-Neumann system correctly
