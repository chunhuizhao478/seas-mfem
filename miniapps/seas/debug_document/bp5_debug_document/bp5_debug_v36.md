# BP5 Debug v36: Interior Dirichlet BR2 Skeleton Fix + Per-Quadrature-Point Traction

**Date**: 2026-03-16
**Status**: Code applied, all tests pass, TACC runs pending
**Previous**: v34 (current baseline), v35 (reverted — same code changes, wrong analysis)
**Branch**: `feature/elasticity`

---

## 1. Problem Statement

BP5 v34 results show three related issues:

| Issue | Observation | Tandem Reference |
|-------|-------------|-----------------|
| Growing antisymmetric tau_dip | tau_dip grows linearly with time, antisymmetric along strike | tau_dip small and stable (0.01-0.2 MPa) |
| Recurrence too long (p=1) | ~271 yr | ~197 yr |
| Fault locking at p>=2 | p=2: ~396 yr; p=4, p=6: NO earthquakes | All p converge to ~197 yr |

## 2. What Happened with v35

v35 identified 3 bugs and applied fixes (commit `a9f8ebb`). It was reverted to v34 (`da6a6c7`) because the analysis concluded that Fixes 1+2 (per-quadrature-point traction) are "exact for p=1" and Fix 3 (interior Dirichlet BR2 skeleton) is "moderate" — so v35 was expected to change nothing at p=1.

**This analysis was wrong.** Re-analysis shows Fix 3 affects p=1 because it changes the boundary loading magnitude on non-fault Y=0 faces by a factor of 2x. These are the very faces whose introduction in v34 caused the growing tau_dip.

v36 re-applies the identical code from v35 with corrected understanding of why it matters at p=1.

## 3. Root Cause Analysis

### 3.1 Bug A: Interior Dirichlet BR2 Skeleton — Missing Cross-Element Lifts + 2x Penalty (affects ALL p)

**File**: `domain/elasticity_operator.hpp`, `AssembleDirichletLoading()`, interior Dirichlet BR2 path

**The bug**: v34 computed the BR2 lifted function separately for each element:

```cpp
// v34 code — elem1 sees ONLY its own lift
real_t ev = 0.0;
for (int m = 0; m < ndof1; m++)
   ev += shapes1(m, q) * f_lifted1(u * dim + s, m);
sum += tn * ev;
...
fl_q1(i, q) = sum;  // NO 0.5 factor
```

But the bilinear form (`AssembleSlipContributionBR2`) uses **both** elements' lifts with a **0.5** factor:

```cpp
// Bilinear form — cross-element evaluation
sum += tn * (eval1 + eval2);  // BOTH elements combined
f_lifted_q(i, q) = 0.5 * sum;  // 0.5 factor
```

Two sub-bugs:
1. **Missing cross-element terms**: Each element's loading only saw its own lift, not the neighbor's. The bilinear form's stabilization matrix uses both lifts.
2. **Missing 0.5 factor**: The loading had `fl_q1(i,q) = sum` vs the bilinear form's `f_lifted_q(i,q) = 0.5 * sum`. The penalty contribution to the RHS was **2x too large**.

**Why this matters at p=1**: The interior Dirichlet faces are the source of the tau_dip growth (confirmed by v33->v34 comparison: v33 had no interior Dirichlet BC and had stable tau_dip; v34 introduced them and tau_dip started growing). The 2x penalty overshoot creates excessive stress at fault edges where the prescribed tectonic loading u_D = +/-Vp*t meets the fault slip. The overshoot manifests as the growing antisymmetric dip stress pattern:

| Station | v33 tau_dip @50yr | v34 @50yr | v34 @250yr |
|---------|-------------------|-----------|------------|
| strk+00dp+00 | 0.031 | 0.064 | **0.553** |
| strk+36dp+00 | -0.289 | -0.264 | **-0.795** |
| strk-36dp+00 | 0.347 | 0.315 | **0.848** |

The doubling of the v33 values at 50yr (0.031 -> 0.064) and continued linear growth are consistent with a 2x penalty overshoot on the loading side.

### 3.2 Bug B: ComputeTraction Evaluates at Centroid Only (affects p>=2)

**File**: `domain/elasticity_operator.hpp`, `ComputeTraction()`

The average stress traction {sigma.n} was evaluated at the face centroid only:
- For p=1: gradient is constant per element, so centroid evaluation is exact
- For p>=2: gradient varies across the face as O(h^{p-1}), and a single centroid sample is NOT the face average

This explains the progressive fault locking with increasing p: the traction fed to the rate-and-state friction law was increasingly wrong.

### 3.3 Bug C: BR2 Correction Uses avg_shapes (affects p>=2)

**File**: `domain/elasticity_operator.hpp`, `ComputeTraction()` BR2 path

The lifted function was evaluated using face-averaged shapes (`avg_shape1`, `avg_shape2`) instead of per-quadrature-point shapes. For p>=2, the lifted function is a polynomial of degree p whose spatial variation is lost when projected through avg_shapes.

The Minv eigenvalues grow as O(p^3/h^3), so high-degree modes in f_lifted are massively amplified — and avg_shapes completely misrepresents these modes.

### 3.4 Why Tandem Doesn't Have These Issues

Tandem evaluates traction at each face quadrature point natively:
```python
traction_q = 0.5 * (traction(0, n_unit_q) + traction(1, n_unit_q))
           + c0 * (E_q[0] * u[0] - E_q[1] * u[1] - f_q)
```

Both {sigma.n} and the BR2 penalty correction are per-quadrature-point. The result is L2-projected to fault DOFs, producing accurate face-averaged values at any polynomial order.

For the interior Dirichlet faces, Tandem processes them through the same skeleton flux path as fault faces, with `flux_u_add_bc` and `flux_sigma_add_bc` kernels that correctly use `0.5*sign` decomposition. No separate loading assembly needed.

## 4. Fixes Applied

All three fixes are re-applied from v35 commit `a9f8ebb` via `git diff da6a6c7..a9f8ebb`.

### Fix 1: Interior Dirichlet BR2 Skeleton — Cross-Element Lifts + 0.5 Factor

Replaced separate `fl_q1` (elem1-only) and `fl_q2` (elem2-only) with combined `f_lifted_q`:

```cpp
// v36: Cross-element evaluation at each quadrature point
real_t eval1 = 0.0, eval2 = 0.0;
for (int m = 0; m < ndof1; m++)
   eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
for (int m = 0; m < ndof2; m++)
   eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
sum += tn * (eval1 + eval2);  // BOTH elements
...
f_lifted_q(i, q) = 0.5 * sum;  // Matches bilinear form
```

The same `f_lifted_q` is used for both elem1 (+sign) and elem2 (-sign) RHS contributions.

### Fix 2: Per-Quadrature-Point Stress Traction in ComputeTraction

Rewrote ComputeTraction (both local interior and shared parallel faces) to:
1. Loop over face quadrature points (not centroid)
2. At each q: compute grad1_q, grad2_q -> avg_grad_q -> strain_q -> stress_q -> T_stress_q
3. Face-average: `T = sum_q(T_q * w_q) / sum_q(w_q)`

### Fix 3: Per-Quadrature-Point BR2 Correction in ComputeTraction

Evaluate the lifted function at each face quadrature point using per-point shapes with cross-element evaluation:

```cpp
// v36: Per-point evaluation with cross-element lifting
real_t eval1 = 0.0, eval2 = 0.0;
for (int m = 0; m < ndof1; m++)
   eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
for (int m = 0; m < ndof2; m++)
   eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
sum += tn * (eval1 + eval2);
...
correction_q[i] = br2_penalty * 0.5 * sum;
```

This replaces the previous `avg_shape1(m)` / `avg_shape2(m)` evaluation.

## 5. Impact Summary

| Fix | What it corrects | p=1 effect | p>=2 effect |
|-----|-----------------|------------|-------------|
| Fix 1 | Interior Dirichlet BR2 loading 2x too large + missing cross-element terms | **PRIMARY** — eliminates tau_dip growth, should shorten recurrence | Also improves accuracy |
| Fix 2 | Centroid-only stress traction | Exact for p=1 (no change) | **CRITICAL** — enables correct traction at higher p |
| Fix 3 | avg_shapes BR2 correction | Negligible for p=1 | **CRITICAL** — prevents fault locking at p>=2 |

## 6. Files Modified

| File | Lines | Change |
|------|-------|--------|
| `domain/elasticity_operator.hpp` | ~2185-2330 | Fix 1: Interior Dirichlet BR2 skeleton |
| `domain/elasticity_operator.hpp` | ~2625-2989 | Fix 2+3: ComputeTraction per-quad-point (interior faces) |
| `domain/elasticity_operator.hpp` | ~3003-3350 | Fix 2+3: ComputeTraction per-quad-point (shared faces) |

Total: 414 insertions, 430 deletions (mostly structural rewrite of ComputeTraction).

## 7. Verification

### 7.1 Build
```
conda activate mfem-dev && make -j8
```
Build succeeds with no errors (only benign duplicate -rpath linker warnings).

### 7.2 Serial Unit Tests
All test suites pass:

| Suite | Tests | Status |
|-------|-------|--------|
| seas_test_antiplane | all | PASS |
| seas_test_bp2_short | all | PASS |
| seas_test_bp5_mesh | 25 | PASS |
| seas_test_bp5_output | 45 | PASS |
| seas_test_bp5_params | 96 | PASS |
| seas_test_br2_consistency | all | PASS |
| seas_test_checkpoint | 21 | PASS |
| seas_test_domain_interface | 32 | PASS |
| seas_test_elasticity_br2 | 46 | PASS |
| seas_test_elasticity_operator | 84 | PASS |
| seas_test_fault_basis | 187 | PASS |
| seas_test_fault_detection | 7 | PASS |
| seas_test_fault_operator | all | PASS |
| seas_test_friction_law | all | PASS |
| seas_test_io | all | PASS |
| seas_test_mesh_regions | 7 | PASS |
| seas_test_mpi_context | 15 | PASS |
| seas_test_parallel_domain | 11 | PASS |
| seas_test_parallel_elasticity | 6 | PASS |
| seas_test_parallel_fault | all | PASS |
| seas_test_parallel_utils | 9 | PASS |
| seas_test_psi_state | all | PASS |
| seas_test_quasi_dynamic | all | PASS |
| seas_test_scaling | 3 | PASS |
| seas_test_serial_parallel_consistency | all | PASS |
| seas_test_state_evolution | all | PASS |
| seas_test_vector_friction | 45 | PASS |
| seas_test_bp5_parallel_smoke | 13 | PASS |

### 7.3 Parallel Unit Tests (8 MPI ranks)
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

## 8. TACC Test Plan

### 8.1 Full Runs (1800 yr each)

| Job | Config | Nodes | Ranks | Key Question |
|-----|--------|-------|-------|-------------|
| `bp5_v36_1000m_br2_p1_full` | p=1, 1000m, BR2 | 8 | 400 | Does tau_dip stop growing? Recurrence shorter than 271 yr? |
| `bp5_v36_1000m_br2_p2_full` | p=2, 1000m, BR2 | 16 | 800 | Does p=2 converge toward p=1 (not away from it)? |
| `bp5_v36_2500m_br2_p4_full` | p=4, 2500m, BR2 | 16 | 800 | Do earthquakes nucleate? (v34: locked for 1800 yr) |
| `bp5_v36_4000m_br2_p6_full` | p=6, 4000m, BR2 | 16 | 800 | Do earthquakes nucleate? Match Tandem's uphoff submission? |

### 8.2 Expected Outcomes

**p=1 (primary test for Fix 1)**:
- tau_dip should stop growing linearly (the 2x penalty overshoot is eliminated)
- Recurrence should shorten from ~271 yr toward Tandem's ~197 yr
- The magnitude of improvement depends on how much the interior Dirichlet overshoot was delaying nucleation

**p=2 (tests Fix 1 + Fix 2 + Fix 3 combined)**:
- Should converge TOWARD p=1 recurrence, not away from it (v34 had p=2 at 396 yr vs p=1 at 271 yr)
- Per-quadrature-point traction (Fix 2) now captures gradient variation
- Per-point BR2 correction (Fix 3) now captures lifted function variation

**p=4, p=6 (critical tests for Fixes 2+3)**:
- Earthquakes should nucleate (v34: complete fault locking)
- These configs match Tandem's SCEC submissions (uphoff.2 at p=4/2.5km, uphoff at p=6/4km)
- If earthquakes nucleate with recurrence ~197 yr, confirms formulation is correct
- Any remaining mismatch at p=4/p=6 would indicate a different bug

### 8.3 Diagnostic Checks

When results arrive, check:

1. **tau_dip time series** at strk+00dp+00, strk+36dp+00, strk-36dp+00:
   - v34: grows linearly (0.06 MPa at 50yr -> 0.55 MPa at 250yr)
   - v36 expected: stable or oscillating around small value

2. **V_max time series**: earthquake nucleation timing and peak slip rate

3. **tau_strike along-strike profile**: should be symmetric (antisymmetric tau_dip was the bug signature)

4. **Convergence with p**: recurrence intervals should converge as p increases:
   p=1 >= p=2 >= p=4 ~= p=6 ~= 197 yr (Tandem reference)

## 9. Relationship to v35

v36 applies the **identical code changes** as v35 (commit `a9f8ebb`). The difference is:

| | v35 | v36 |
|--|-----|-----|
| Analysis of p=1 impact | "Fixes 1+2 exact for p=1, Fix 3 moderate" | "Fix 3 is the primary p=1 fix (2x penalty overshoot)" |
| Why reverted | "didn't fix anything at p=1" (expected no change) | N/A — now understand it SHOULD change p=1 |
| TACC results | Submitted but reverted before results analyzed | Pending |

The revert was premature. Fix 3 directly affects the interior Dirichlet loading magnitude, which is the source of the tau_dip growth that appeared in v34. The "exact for p=1" reasoning only applies to Fixes 1+2 (ComputeTraction), not to Fix 3 (AssembleDirichletLoading).

## 10. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H32 | Previous fixes (see v21-v31 docs) | Done |
| v32 | h-refinement study (250m, 500m) | Done |
| v33 | p-refinement (p=2 general order support) | Done |
| v34 | Non-fault Y=0 interior Dirichlet BC (jump from element centroids) | Done |
| v35 | Per-quad-point traction + BR2 skeleton fix (reverted prematurely) | Reverted |
| **v36** | **Re-apply v35 fixes with corrected analysis** | **Applied, TACC pending** |
