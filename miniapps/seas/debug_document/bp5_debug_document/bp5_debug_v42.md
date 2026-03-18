# BP5 Debug v42: IP Penalty `c_N_1` Hardcoded to 1.0 — p≥2 Fault Locking Fix

**Date**: 2026-03-18
**Status**: Code fix applied, all 528 tests pass, ready for TACC runs
**Previous**: v41 (parametric study — p≥2 fault locking identified as blocking bug)
**Branch**: `feature/elasticity`

---

## 1. Problem Statement

v41 parametric study revealed that **IP p=2 and p=4 produce zero earthquakes** over 1800 yr.
The fault locks up completely: V decays from 1e-9 to 1e-13 m/s, tau_strike drops from
13.27 to 9.2 MPa (p=2) or 8.5 MPa (p=4). This prevents using p-refinement to resolve
the p=1 discretization artifacts (dip mismatch, 20% strike deficit).

## 2. Root Cause: `c_N_1` Hardcoded to 1.0 in All IP RHS Paths

### 2.1 The Bug

The IP penalty constant `c_N_1` (inverse trace inequality constant) is used in the penalty
formula: `η_F = (D+1) × c_N_1 × (A/V) × (c₁²/c₀)`.

The **bilinear form** (`DGElasticityIPPenaltyIntegrator`, line 113) and **ComputeTraction()**
(lines 2985, 3344) correctly use the order-dependent formula:

```cpp
real_t c_N_1 = order_ * (order_ + dim - 1.0) / dim;
// p=1: 1.0,  p=2: 2.667,  p=4: 8.0,  p=6: 16.0
```

But **all 5 RHS assembly paths** had `c_N_1` hardcoded to `1.0`:

```cpp
real_t p0 = (dim + 1) * 1.0 * (nl_q / detJ1) * (c1_mat * c1_mat / c0_mat);
//                      ^^^
//                      Should be c_N_1 = order_ * (order_ + dim - 1.0) / dim
```

### 2.2 Affected Locations (All in `elasticity_operator.hpp`)

| # | Function | Lines | Purpose |
|---|----------|-------|---------|
| 1 | `AssembleSlipContributionIP()` | 1007-1008 | Interior fault face slip RHS |
| 2 | `AssembleSlipContributionIPShared()` | 1437-1438 | Shared fault face slip RHS |
| 3 | `AssembleDirichletLoading()` | 1790 | Boundary Dirichlet RHS |
| 4 | `AssembleDirichletLoading()` | 2062-2064 | Interior Dirichlet RHS (Y=0 faces) |
| 5 | `AssembleDirichletLoading()` | 2441-2444 | Shared Dirichlet RHS (parallel Y=0) |

### 2.3 Why p=1 Was Unaffected

At p=1: `c_N_1 = 1 × (1+3-1)/3 = 1.0`. The hardcoded `1.0` equals the correct value.
The bug is **invisible at p=1** — which is why v40 worked perfectly.

### 2.4 Impact at p≥2

| Order | c_N_1 (correct) | RHS used | LHS/RHS ratio | Effect |
|-------|-----------------|----------|---------------|--------|
| p=1 | 1.0 | 1.0 | **1.0×** | No effect |
| p=2 | 2.667 | 1.0 | **2.67×** | RHS under-penalized by 2.67× |
| p=4 | 8.0 | 1.0 | **8.0×** | RHS under-penalized by 8× |
| p=6 | 16.0 | 1.0 | **16×** | RHS under-penalized by 16× |

The bilinear form (LHS) has `c_N_1`× more penalty than the RHS at p≥2. This means:

1. **Slip assembly**: The prescribed fault slip `delta_u` is enforced with penalty 1.0,
   but the bilinear form penalizes jumps with penalty `c_N_1`. The net effect: the DG
   solution has `[[u]] ≈ delta_u / c_N_1` instead of `[[u]] ≈ delta_u`. The effective slip
   is **divided by c_N_1**, reducing tectonic loading by 2.67× (p=2) or 8× (p=4).

2. **Dirichlet loading**: Similarly, the boundary plate loading `u_D = sgn(Y) × Vp×t/2`
   is enforced 2.67-8× too weakly relative to the bilinear form constraint.

3. **Combined effect**: The fault receives dramatically insufficient tectonic loading.
   Stress decays instead of building up. The fault locks — no nucleation possible.

### 2.5 Tandem Reference

Tandem's formula (from `InverseInequality.h`):
```cpp
constexpr static double trace_constant(unsigned N) {
    return (N + 1) * (N + D) / static_cast<double>(D);
}
```
Called with `N = PolynomialDegree - 1`:
- p=1: trace_constant(0) = 1×3/3 = 1.0
- p=2: trace_constant(1) = 2×4/3 = 2.667
- p=4: trace_constant(3) = 4×6/3 = 8.0

This matches our formula `c_N_1 = p*(p+D-1)/D` exactly. Tandem uses this consistently
in **all** penalty computations (bilinear form, RHS, and traction) through the single
`penalty_[fctNo]` precomputed value. There is no opportunity for the LHS and RHS to
use different penalties.

### 2.6 Why ComputeTraction Was Already Correct

The traction extraction (`ComputeTraction()`) was rewritten in v35-v36 with per-quad-point
evaluation. During that rewrite, the IP penalty was correctly computed with `c_N_1`:

```cpp
// ComputeTraction, line 2985:
real_t c_N_1 = order_ * (order_ + dim - 1.0) / dim;  // ← Correct!
```

The RHS paths (`AssembleSlipContribution*` and `AssembleDirichletLoading`) were written
earlier (v28-v34) when only p=1 was tested, and the `1.0` was never updated.

## 3. The Fix

**One-line change at each of 5 locations**: Replace `(dim + 1) * 1.0 *` with
`(dim + 1) * c_N_1 *`, adding the computation `c_N_1 = order_ * (order_ + dim - 1.0) / dim`
before each penalty formula.

### Example (Location 1, lines 1007-1010):

```cpp
// BEFORE (v41):
real_t p0 = (dim + 1) * 1.0 * (nl_q / detJ1) * (c1_mat * c1_mat / c0_mat);
real_t p1 = (dim + 1) * 1.0 * (nl_q / detJ2) * (c1_mat * c1_mat / c0_mat);

// AFTER (v42):
real_t c_N_1 = order_ * (order_ + dim - 1.0) / dim;
real_t p0 = (dim + 1) * c_N_1 * (nl_q / detJ1) * (c1_mat * c1_mat / c0_mat);
real_t p1 = (dim + 1) * c_N_1 * (nl_q / detJ2) * (c1_mat * c1_mat / c0_mat);
```

### Consistency Verification

After fix, all 7 IP penalty locations use the same `c_N_1`:

| Location | Lines | `c_N_1` | Status |
|----------|-------|---------|--------|
| DGElasticityIPPenaltyIntegrator (bilinear form) | integrator line 113 | `p*(p+dim_-1)/dim_` | ✓ Already correct |
| AssembleSlipContributionIP (interior) | 1008-1009 | `order_*(order_+dim-1)/dim` | ✓ **Fixed** |
| AssembleSlipContributionIPShared (shared) | 1441-1442 | `order_*(order_+dim-1)/dim` | ✓ **Fixed** |
| AssembleDirichletLoading (boundary) | 1796 | `order_*(order_+dim-1)/dim` | ✓ **Fixed** |
| AssembleDirichletLoading (interior Dirichlet) | 2070-2072 | `order_*(order_+dim-1)/dim` | ✓ **Fixed** |
| AssembleDirichletLoading (shared Dirichlet) | 2451-2453 | `order_*(order_+dim-1)/dim` | ✓ **Fixed** |
| ComputeTraction (interior) | 2999-3001 | `order_*(order_+dim-1)/dim` | ✓ Already correct |
| ComputeTraction (shared) | 3358-3360 | `order_*(order_+dim-1)/dim` | ✓ Already correct |

## 4. Why This Is the Complete Fix

1. **No other structural issues found**: The IP path correctly uses per-quad-point
   evaluation in all three functions (ComputeTraction, AssembleSlipContribution,
   AssembleDirichletLoading). The centroid-only bug that plagued BR2 in v33-v34
   does not exist in the current IP path.

2. **BR2 is unaffected**: BR2 uses `penalty = dim+1` (number of faces per element),
   which is order-independent. The `c_N_1` formula is IP-specific.

3. **p=1 regression**: `c_N_1 = 1.0` at p=1, so results are **identical** to v40/v41a.

## 5. Verification

### Build
```
conda activate mfem-dev && make -j8
```
Build succeeds (only benign duplicate -rpath linker warnings).

### Serial Unit Tests
All test suites pass:

| Suite | Tests | Status |
|-------|-------|--------|
| seas_test_elasticity_operator | 84 | PASS |
| seas_test_elasticity_br2 | 46 | PASS |
| seas_test_fault_basis | 187 | PASS |
| seas_test_domain_interface | 32 | PASS |
| seas_test_bp5_params | 96 | PASS |

### Parallel Unit Tests (8 MPI ranks)
All parallel suites pass:

| Suite | Tests | Status |
|-------|-------|--------|
| seas_test_parallel_domain | 11 | PASS |
| seas_test_parallel_elasticity | 6 | PASS |
| seas_test_parallel_fault | 27 | PASS |
| seas_test_parallel_utils | 9 | PASS |
| seas_test_serial_parallel_consistency | 12 | PASS |
| seas_test_br2_consistency | 7 | PASS |
| seas_test_bp5_parallel_smoke | 13 | PASS |

Total: 528 tests (serial + parallel), all pass.

## 6. Expected Impact

| Config | v41 (bug) | v42 (expected) | Tandem |
|--------|-----------|----------------|--------|
| p=1, 1000m | 8 events, ~250 yr | **Same** (c_N_1=1.0) | 8 events, ~240 yr |
| p=2, 1000m | 0 events (locked) | **7-8 events, ~240 yr** | 8 events, ~240 yr |
| p=4, 2500m | 0 events (locked) | **7-8 events, ~240 yr** | 8 events, ~240 yr |

At p=2/p=4:
- Dip quantities should be near-zero (as v41c/v41d showed before locking)
- Strike slip deficit should be ~4% (p=2) or ~0.1% (p=4) instead of ~20% (p=1)
- Recurrence should converge toward Tandem's 240 yr

## 7. Next Steps

1. **Submit v42 TACC runs**:
   - v42a: IP p=1, 1000m (regression check, should match v41a exactly)
   - v42b: IP p=2, 1000m (critical test — should nucleate)
   - v42c: IP p=4, 2500m (match Tandem uphoff.2 — definitive convergence test)

2. **If p=2 nucleates and matches Tandem**: The formulation is verified correct. The
   remaining dip mismatch and strike deficit at p=1 are confirmed as discretization
   artifacts. Use p=2 or higher for benchmark submissions.

3. **If p=2 still locks**: Other p≥2-specific issues may exist (unlikely given the
   clear root cause analysis). Would need detailed traction diagnostics at p=2.

## 8. Files Changed

| File | Change |
|------|--------|
| `domain/elasticity_operator.hpp` | Replace `1.0` with `c_N_1` at 5 IP RHS penalty locations |

## 9. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| v30-v40 | All previous fixes | Done |
| v41 | Parametric study: identified p≥2 fault locking | Analysis done |
| **v42** | **IP `c_N_1` fix: 1.0 → p(p+D-1)/D at 5 RHS locations** | **Applied, tests pass** |
