# BP5 Debug v35: Per-Quadrature-Point Traction for General Polynomial Order

## Problem

v34 results show progressive fault locking with increasing polynomial order:
- **p=1**: Earthquakes at ~271 yr recurrence (works, but longer than Tandem's ~197 yr)
- **p=2**: Earthquakes at ~396 yr (significantly delayed)
- **p=4**: NO earthquakes in 1800 yr (fault locked, 17-20% of expected interseismic slip)
- **p=6**: NO earthquakes in 1291 yr (fault locked, similar to p=4)

This indicates a formulation bug that grows with polynomial order, not just discretization error.

## Root Cause Analysis

### Bug 1 (CRITICAL): ComputeTraction evaluates {sigma.n} at face centroid only

**Location**: `domain/elasticity_operator.hpp`, ComputeTraction (interior + shared faces)

The average stress traction {sigma.n} was computed by:
1. Evaluating gradients (dshape) at the face centroid from each side
2. Averaging the two gradients
3. Computing strain -> stress -> traction at that single point

For p=1 (linear elements): gradient is constant per element, centroid is exact.
For p>=2: gradient varies across the face. Centroid evaluation gives ONE sample of
a polynomial that varies as O(h^{p-1}) across the face. This is NOT the face average.

**Impact**: Wrong mean stress traction at the fault. Error grows with p because higher-order
polynomials have more variation across the face.

### Bug 2 (CRITICAL): BR2 correction uses avg_shapes evaluation

**Location**: `domain/elasticity_operator.hpp`, ComputeTraction BR2 path

The BR2 lifting correction evaluated the lifted function using "face-averaged shapes":
```cpp
avg_shape1(m) = sum_q w_q * shapes1(m,q) / sum_w
correction[i] = br2_penalty * 0.5 * sum(TN * avg_shape1(m) * f_lifted1(u*dim+s, m))
```

The lifted function f_lifted is a polynomial of degree p. Evaluating it at avg_shapes
(a single fixed set of shape function values) instead of at each quadrature point
loses spatial variation information.

For p=1: f_lifted is linear, avg_shape evaluation is reasonable.
For p>=2: f_lifted varies across the face, evaluation at avg_shapes introduces errors
proportional to the polynomial's variation. At p=4-6, f_lifted has high-degree modes
whose Minv-amplified values are incorrectly captured by avg_shapes.

**Impact**: Wrong BR2 stabilization correction. Combined with Bug 1, the total fault
traction is wrong, causing the rate-and-state friction law to produce incorrect slip rates.

### Bug 3 (MODERATE): Interior Dirichlet BR2 skeleton missing cross-element lifts + wrong scaling

**Location**: `domain/elasticity_operator.hpp`, AssembleDirichletLoading BR2 skeleton path

Two issues in the interior Dirichlet skeleton (non-fault Y=0 faces) BR2 path:

**3a: Missing cross-element lifting terms**
The v34 code computed the lifted function separately for each element:
```cpp
fl_q1(i,q) = TN * shapes1(m,q) * f_lifted1(m)  // elem1 self-lift ONLY
fl_q2(i,q) = TN * shapes2(m,q) * f_lifted2(m)  // elem2 self-lift ONLY
```

The bilinear form (BR2 integrator) uses full cross-element lifting:
```cpp
sum += tn * (eval1 + eval2);  // BOTH elements combined
```
where eval1 evaluates f_lifted1 at shapes1, eval2 evaluates f_lifted2 at shapes2.

**3b: Missing 0.5 factor on combined evaluation**
The bilinear form computes L_q = 0.5 * tn * (eval1 + eval2), with Lift *= 0.5 (for
the {psi} average). The interior Dirichlet code had the 0.5 in face_int (matching
Lift *= 0.5) but was missing the second 0.5 factor on f_lifted_q, making the penalty
contribution 2x too large.

**Impact**: Incorrect BC loading at non-fault Y=0 faces. These faces are far from the
main fault region, so the effect on traction is indirect (through the displacement
solution), but contributes to inaccuracy at higher p.

### Why Tandem doesn't have these issues

Tandem's traction computation (`Elasticity.cpp` line 972, `elasticity.py` line 242):
```python
traction_q = 0.5 * (traction(0, n_unit_q) + traction(1, n_unit_q))
           + c0 * (E_q[0] * u[0] - E_q[1] * u[1] - f_q)
```

1. **Per-quadrature-point evaluation**: Both {sigma.n} and the penalty correction
   are evaluated at each face quadrature point q, not at the centroid.
2. **Simple scalar penalty**: Uses c0 = -NumFacets = -4 for BR2 (from `penalty()`
   in Elasticity.h), which is dimensionless. At p=6 the jump residual is O(h^6),
   making the penalty term negligible regardless of dimensional consistency.
3. **L2 projection**: Per-point traction is projected to fault DOFs via
   `ElasticityAdapter::traction()`, producing accurate face-averaged values.

For our low-order (p=1, p=2) simulations, the penalty term is NOT negligible,
so we need the BR2 lifting (which IS dimensionally correct in our code).

## Fixes Applied

### Fix 1+2: Per-Quadrature-Point Traction (ComputeTraction, interior + shared faces)

Rewrote `ComputeTraction` for both local interior faces and shared (parallel) faces:

1. **Loop over face quadrature points** (not centroid):
   - At each q: compute grad1_q, grad2_q from both elements
   - At each q: avg_grad_q -> strain_q -> stress_q -> T_stress_q = stress * n_hat

2. **BR2 correction per quadrature point**:
   - Face integrals use per-point shapes (correct, from v34)
   - Evaluate lifted function at each q using per-point shapes (NOT avg_shapes):
     ```
     correction_q = br2_penalty * 0.5 * sum_us(TN * (shapes1(m,q)*f_lifted1 + shapes2(m,q)*f_lifted2))
     ```
   - This matches the bilinear form's per-point evaluation

3. **Face-average** the per-point traction:
   ```
   T = sum_q T_q * w_q * |J_q| / sum_q w_q * |J_q|
   ```
   For flat tet faces: |J_q| = const, simplifies to weighted average.

### Fix 3: Interior Dirichlet BR2 Skeleton Cross-Element Lifts

Replaced separate self-lift `fl_q1` (elem1 only) and `fl_q2` (elem2 only)
with a combined `f_lifted_q` that evaluates both elements' lifted functions:

```cpp
// Cross-element evaluation: sum over BOTH elements
real_t eval1 = 0.0, eval2 = 0.0;
for (int m = 0; m < ndof1; m++)
   eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
for (int m = 0; m < ndof2; m++)
   eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
sum += tn * (eval1 + eval2);
...
f_lifted_q(i, q) = 0.5 * sum;  // Match bilinear form's L_q = 0.5 * sum
```

The same combined `f_lifted_q` is used for both elem1 (+sign) and elem2 (-sign)
RHS contributions, matching the bilinear form's structure.

## Files Modified

| File | Change |
|------|--------|
| `domain/elasticity_operator.hpp` | ComputeTraction interior faces: per-quad-point {sigma.n} + BR2 correction |
| `domain/elasticity_operator.hpp` | ComputeTraction shared faces: same per-quad-point rewrite |
| `domain/elasticity_operator.hpp` | Interior Dirichlet BR2 skeleton: cross-element lifts + 0.5 factor |

## Verification

### Completed
1. **Build**: `conda activate mfem-dev && make -j8` — ✅ No errors/warnings
2. **Serial unit tests**: 22/22 test suites pass (regression check) — ✅
3. **Parallel unit tests**: 8/8 test suites pass with 8 MPI ranks — ✅

### Pending
4. Run BP5 p=1 on 8 ranks: verify traction values similar to v34
5. Run BP5 p=2 on 8 ranks: verify improved (shorter) recurrence interval
6. TACC runs at p=1, p=2, p=4: compare with v34 results
