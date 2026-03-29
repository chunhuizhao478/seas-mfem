# BP5 Debug v55: In-Depth MFEM vs Tandem Code Comparison

**Date**: 2026-03-29
**Status**: INVESTIGATION — SYSTEMATIC DIFFERENCES CATALOGUED
**Previous**: v54 (CG+AMG disproves solver hypothesis, penalty correction dominance confirmed)
**Branch**: `feature/elasticity`

---

## 1. Context

v54 conclusively proved that the linear solver (MUMPS-BLR vs CG+AMG) is NOT the root cause
of nucleation failure. Both solvers produce identical results. The penalty correction dominance
(3-6x the physical stress at near-front stations) is intrinsic to the IP DG formulation.

**This document performs a line-by-line comparison** of the MFEM and Tandem BP5 implementations
to identify every difference that could contribute to the ~0.38 MPa traction excess that kills
nucleation.

---

## 2. Summary of Differences Found

| # | Area | Severity | Description |
|---|------|----------|-------------|
| D1 | Traction recovery: normal vector | **HIGH** | Tandem uses `n_unit_q` (unit normal), MFEM uses face-constant `basis.normal` |
| D2 | Traction recovery: L2 projection weights | **HIGH** | Tandem weights by `nl_q * w_q` (physical surface measure), MFEM uses `w_q` only |
| D3 | Traction recovery: penalty × nl_q | **HIGH** | Tandem's penalty term implicitly picks up `nl_q` via projection; MFEM's doesn't |
| D4 | Slip rate solver: search space | **MEDIUM** | Tandem: log10(V) ∈ [-32, log10(τ/η)]; MFEM: linear V ∈ [0, τ/η] |
| D5 | Stiffness: two separate integrators | **MEDIUM** | MFEM splits consistency+symmetry from penalty into two integrators |
| D6 | Stiffness: quadrature order | **LOW** | Penalty integrator uses 2p (MFEM) vs 2p+1 (Tandem) |
| D7 | RHS slip: consistency term | **NEEDS AUDIT** | MFEM's manual consistency loop vs Tandem's generated kernel may differ |
| D8 | Slip rate sign convention | **NEEDS VERIFICATION** | Tandem negates V_vec; MFEM doesn't — must be self-consistent through chain |
| D9 | sigma_n clamp | **LOW** | MFEM clamps at 10% of sigma_n_bp5; Tandem does not |
| D10 | Fault basis: per-face vs per-quad-point | **LOW** | MFEM: single basis per face; Tandem: per quad point |
| D11 | Traction recovery: stress formula | **NEEDS AUDIT** | MFEM uses `{ε}` then σ; Tandem uses `{σ·n}` directly |

---

## 3. D1: Traction Recovery — Normal Vector Convention

### 3.1 Tandem

Tandem's `compute_traction` kernel (from `elasticity.py:242-244`):

```python
traction_q['pq'] <= 0.5 * (traction(0, n_unit_q) + traction(1, n_unit_q)) +
                    c0[0] * (E_q[0]['lq']*u[0]['lp'] - E_q[1]['lq']*u[1]['lp'] - f_q['pq'])
```

Where `traction(x, normal)` is:
```python
lam_q[x]['q'] * Dx_q[x]['lsq'] * u[x]['ls'] * normal['pq'] +
mu_q[x]['q'] * (Dx_q[x]['ljq']*u[x]['lp']*normal['jq'] + Dx_q[x]['lpq']*u[x]['lj']*normal['jq'])
```

**Both stress average and penalty use `n_unit_q` (unit normal at each quad point).**

The `c00 = -penalty(fctNo)` coefficient (from `Elasticity.cpp:972`) is a scalar penalty.

### 3.2 MFEM

MFEM's `ComputeTractionImpl` (line 3896-3908):

```cpp
T_stress_q[ci] += stress_ij * basis.normal[cj];  // face-constant unit normal
```

And penalty correction (line 3931-3932):

```cpp
real_t jump_c = (u1q - u2q) - sign * delta_u_q_c;
correction_q[c] = -penalty_ip * sign * jump_c;  // scalar, no normal weighting
```

### 3.3 Impact

For the **stress average**: MFEM uses `basis.normal` (a single unit normal at the face centroid).
Tandem uses `n_unit_q` (unit normal at each quadrature point). For flat faces (linear tets),
these are identical. **No difference at p=1 on linear meshes.**

For the **penalty correction**: Tandem's formula is
`c00 * (u_jump - f)` where `c00 = -penalty`. This is a pointwise traction at the quad point,
with units of Pa. MFEM's formula is `-penalty_ip * sign * jump` which is also Pa.

**The key difference emerges in the L2 projection** — see D2 and D3.

---

## 4. D2: Traction L2 Projection Weights

### 4.1 Tandem

Tandem's adapter projects traction to fault DOFs via (`elasticity_adapter.py:26-27`):

```python
traction['kp'] <= minv['lk'] * e_q_T['ql'] * w['q'] * nl_q['q'] *
                  traction_q['oq'] * fault_basis_q['opq']
```

Where `minv` is the inverse of the physical-surface mass matrix:

```python
m['kl'] <= e_q['kq'] * w['q'] * nl_q['q'] * e_q['lq']
```

The projection is:

```
T_fault[k,p] = M^{-1}[l,k] × Σ_q φ_l(q) × w_q × |n_q| × T_3D[o,q] × F[o,p,q]
```

where `|n_q| = nl_q` is the face Jacobian determinant (physical surface area element).

### 4.2 MFEM

MFEM's `GalerkinProject` (line 247-250):

```cpp
val += ir_.IntPoint(q).weight * e_q_(l, q) * quad_vals(c * nq_ + q);
```

With mass matrix:
```cpp
M_ref(i,j) = Σ_q w_q × φ_i(q) × φ_j(q)
```

**No `nl_q` factor in either the mass matrix or the RHS.**

### 4.3 Analysis

Tandem: `T_nodal = M_phys^{-1} × Σ_q φ × w_q × nl_q × T_q`  (physical L2 projection)
MFEM:   `T_nodal = M_ref^{-1}  × Σ_q φ × w_q × T_q`           (reference L2 projection)

For flat faces, `nl_q = const` across all quad points. Then:

```
M_phys = nl * M_ref
T_nodal_tandem = (nl*M_ref)^{-1} × nl × Σ_q φ w T_q = M_ref^{-1} × Σ_q φ w T_q = T_nodal_mfem
```

**The `nl_q` cancels for flat faces.** So at p=1 on linear tets, D2 produces no difference.

However, the traction values being projected are different — see D3.

---

## 5. D3: Penalty Term × nl_q in Projection

### 5.1 Tandem's Effective Penalty in Traction

Tandem's `compute_traction` at quad point `q`:
```
T_q[p] = 0.5*(σ_0·n̂ + σ_1·n̂)[p] - penalty * (u_jump - f)[p]
```

This `T_q` is then projected with weight `nl_q * w_q`:
```
RHS_traction[l] = Σ_q φ_l(q) * w_q * nl_q * T_q[p]
```

So the penalty contribution to the projected traction RHS is:
```
-penalty * Σ_q φ_l(q) * w_q * nl_q * (u_jump - f)[p]
```

**The penalty term in the projected traction is weighted by `nl_q`.**

### 5.2 MFEM's Effective Penalty in Traction

MFEM's `ComputeTractionImpl` at quad point `q`:
```
T_q[c] = T_stress[c] - correction[c]
correction[c] = -penalty_ip * sign * jump_c
```

This `T_q` is projected with weight `w_q` only:
```
RHS_traction[l] = Σ_q φ_l(q) * w_q * T_q[c]
```

So the penalty contribution is:
```
penalty_ip * sign * Σ_q φ_l(q) * w_q * jump_c
```

**No `nl_q` weighting on the penalty contribution.**

### 5.3 Impact

For flat faces on linear tets, `nl_q` is constant. The `nl_q` in Tandem's projection
cancels between numerator and denominator of M^{-1} × RHS. So the final projected traction
nodal values are the same.

**This is a non-issue for flat linear tets.** But for curved meshes or higher-order geometry,
this would matter.

---

## 6. D4: Slip Rate Solver — Log10 vs Linear Brent Search

### 6.1 Tandem (`DieterichRuinaBase::slip_rate`, line 81-175)

```cpp
double Va = -32.0;                          // log10(10^-32)
double Vb = std::log10(tau / eta);          // log10(τ/η)
auto fF = [&](double Ve) {
    double V = std::pow(10.0, Ve);
    double f = arsinh(V / (2.0 * V0) * exp(psi / a)) * a;
    return tau - f * snAbs - eta * V;
};
double Ve = zeroIn(Va, Vb, fF, ...);
double V = std::pow(10.0, Ve);
```

Bracket: `[10^{-32}, τ/η]` in log10 space → 36 orders of magnitude.

### 6.2 MFEM (`DieterichRuinaFriction::SolveSlipRatePsi`, line 324-378)

```cpp
real_t V_lo = 0.0;
real_t V_hi = tau / eta;
auto residual = [&](real_t V) {
    real_t f = a * std::asinh(V / (2.0 * V0) * std::exp(psi / a));
    return tau - sigma_n * f - eta * V;
};
real_t V = zeroIn(V_lo, V_hi, residual, ...);
```

Bracket: `[0, τ/η]` in linear space.

### 6.3 Impact

**Interseismic precision**: During interseismic periods, V ≈ 10^{-12} m/s. In Tandem's log10
search, the root is at position `-12` in the bracket `[-32, 4]`, which is well-resolved.
In MFEM's linear search, the root is at `10^{-12}` in the bracket `[0, 10^4]`. The relative
position is `10^{-16}`, which is at machine precision.

Brent's method converges to within `tol + 4·eps·|x|`. For the linear search:
- Machine eps ≈ 2.2e-16
- |x| ≈ V_hi ≈ 10^4 (for τ ≈ 20 MPa, η ≈ 4600 Pa·s/m)
- Absolute precision ≈ 4·2.2e-16·10^4 ≈ 10^{-12}
- At V ≈ 10^{-12}, relative error ≈ 100%

For the log10 search:
- |x| ≈ 12 (the log10 of V)
- Absolute precision ≈ 4·2.2e-16·12 ≈ 10^{-14}
- V precision ≈ V · ln(10) · 10^{-14} ≈ 10^{-12} · 2.3 · 10^{-14} ≈ 10^{-26}

**The log10 search provides ~14 more digits of relative precision at interseismic rates.**

v54 Test 13 confirmed both methods converge to residual < 1e-16 for a specific test point
(V_ref = 0.0095, coseismic). The difference is most pronounced at very low slip rates,
which is exactly where the nucleation/interseismic transition happens.

**However**, v54 found that both codes agree well at the start (V = 0.01 m/s) and the
divergence begins at t ≈ 3s when V is still ~0.01 m/s. The interseismic precision issue
would only matter much later. **This is unlikely to be the root cause of the initial divergence
but could affect long-term behavior.**

---

## 7. D5: Stiffness Matrix — Two Separate Integrators

### 7.1 Tandem

Single `assembleSurface` kernel computes all three DG terms in one pass:

```python
a[x][y]['kplu'] <=
    c0[x] * E_q[x]['kq'] * w['q'] * traction_op_q[y]['lupq'] +   # consistency
    c1[y] * E_q[y]['lq'] * w['q'] * traction_op_q[x]['kpuq'] +   # symmetry
    c2[abs(y-x)] * w['q'] * E_q[x]['kq'] * L_q[y]['lpuq']        # penalty
```

Coefficients for skeleton (from `Elasticity.cpp:444-467`):
```cpp
// Side 0→0: c00 = -0.5, c10 = +0.5*epsilon, c20 = penalty
// Side 0→1: c01 = +0.5, c11 = -0.5*epsilon, c21 = -penalty
// Side 1→0: c00 = +0.5, c10 = -0.5*epsilon, c20 = -penalty
// Side 1→1: c01 = -0.5, c11 = +0.5*epsilon, c21 = penalty
```

All terms share the same quadrature rule, the same evaluation points, and the same
pre-computed geometric quantities.

### 7.2 MFEM

Two separate integrators added to the BilinearForm:

```cpp
// Integrator 1: Consistency + Symmetry (kappa=0)
new DGElasticityIntegrator(lambda_coeff_, mu_coeff_, epsilon_, 0.0)
// Integrator 2: Penalty (custom)
new DGElasticityIPPenaltyIntegrator(lambda_coeff_, mu_coeff_, 3, penalty_factor_)
```

Each integrator independently:
- Determines its own quadrature order
- Iterates over faces separately
- Computes geometric quantities independently

### 7.3 Impact

For flat linear elements, both integrators are exact for their respective polynomial integrands,
so the combined result should match Tandem's single-pass assembly.

However, the split means:
1. Faces are visited twice (performance impact, not correctness)
2. Any floating-point rounding differences from independent computation of normals, Jacobians,
   etc. can lead to subtle inconsistencies between consistency and penalty contributions.

**More critically**: the `DGElasticityIntegrator` (MFEM built-in) has its own internal
implementation of `{σ(u)·n}·[[v]]`. Any difference in how MFEM's built-in integrator
computes this vs. how the custom penalty integrator computes its term could produce an
inconsistency in the final matrix.

**This is a structural risk but probably not the root cause given both MFEM integrators
have been verified individually.**

---

## 8. D7: RHS Slip Consistency Term

### 8.1 Tandem

For the skeleton RHS, Tandem uses `rhsFacet` kernel (`elasticity.py:174-178`):

```python
b['kp'] <= b['kp'] + c1[0] * tractionTest(0, f_q) * w['q'] +
           c2[0] * w['q'] * E_q[0]['kq'] * f_lifted_q['pq']
```

Where for IP: `f_lifted_q['iq'] <= f_q['iq'] * nl_q['q']`  (line 156)

This is called with `c10 = 0.5 * epsilon` and `c20 = penalty` for side 0
(then `c10 = -0.5*epsilon` and `c20 = -penalty` for side 1).

So the RHS for side 0 is:
```
b_0[k,p] += 0.5*ε * tractionTest(0, f_q) * w_q + penalty * w_q * φ_0[k,q] * f_q[p,q] * nl_q
```

Note: `tractionTest(x, utilde)` uses `n_q` (un-normalized), not `n_unit_q`.

### 8.2 MFEM

`AssembleSlipContributionIP` (line 1288-1316):

```cpp
// Consistency term (symmetry): ε * σ_test(φ_k, e_i) · n * delta_u
for (int k = 0; k < ndof1; k++) {
    real_t grad_dot_n1 = 0.0;
    for (int d = 0; d < dim; d++)
        grad_dot_n1 += dshape1_adj(k, d) * nor(d);  // adj gradient · n_q
    for (int i = 0; i < dim; i++) {
        real_t sym_val1 = 0.0;
        for (int u = 0; u < dim; u++) {
            real_t trac1 = lambda_val_ * dshape1_adj(k, i) * nor(u)
                + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n1
                             + dshape1_adj(k, u) * nor(i));
            sym_val1 += trac1 * delta_u_q[u];
        }
        int idx1 = i * ndof1 + k;
        elvec1(idx1) += epsilon_ * 0.5 * sym_val1 * w1;
        // Penalty
        elvec1(idx1) += wq_penalty * sign * delta_u_q[i] * shape1(k);
    }
}
```

Where `w1 = ip.weight / (2.0 * detJ1)` and `wq_penalty = penalty_ip * ip.weight * nl_q`.

### 8.3 Analysis

The consistency term structure appears to match Tandem's `tractionTest` formula. Both use:
- Adjugate Jacobian (`dshape_adj` in MFEM = `Dx_q` in Tandem, with `g = J^{-T}`)
  - Wait: `dshape_adj = dshape_ref × adj(J)`, while Tandem's `Dx_q = g × Dxi_q = J^{-T} × Dxi_q`
  - `adj(J) = det(J) × J^{-1}`, and `dshape_adj × nor = det(J) × J^{-1} × dshape_ref × nor`
  - In MFEM, `w1 = ip.weight / (2.0 * detJ1)`, so `dshape_adj * w1 = dshape_ref × J^{-1} × ip.weight / 2`
  - In Tandem, `Dx_q` already has the `J^{-T}` applied, and the coefficient includes `w_q * 0.5`

Both should produce the same result. **However**, MFEM uses the un-normalized normal `nor`
(from `CalcOrtho`) for the consistency term dot product, while Tandem's `tractionTest`
also uses `n_q` (un-normalized). This is consistent.

The penalty terms:
- MFEM: `penalty_ip * ip.weight * nl_q * sign * delta_u_q * shape`
- Tandem: `penalty * w_q * φ * f_q * nl_q`

These should match if `penalty_ip = penalty` (verified in v54 tests) and the slip is consistent.

**The sign handling is the main concern.** MFEM uses `sign * delta_u_q` where sign is based
on `nor(1) > 0`. Tandem sets the slip through `bc_skeleton` which checks the face orientation
and provides `f_q` with appropriate sign. Whether these produce identical results depends on
the mesh-specific face orientation convention — verified in v54 tests to match for uniform slip.

**No clear bug found, but the complexity of the manual loop vs generated kernel makes a subtle
difference hard to spot without running both on identical input.**

---

## 9. D8: Slip Rate Sign Convention

### 9.1 Tandem (`DieterichRuinaBase::slip_rate`, line 174)

```cpp
return -(V / tauAbs) * tauAbsVec;
```

Slip rate is **anti-parallel** to traction direction.

### 9.2 MFEM (`DieterichRuinaFriction::SolveSlipRateVectorPsi`, line 416-417)

```cpp
V_vec[0] = (V_abs / tau_abs) * tau_vec[0];
V_vec[1] = (V_abs / tau_abs) * tau_vec[1];
```

Slip rate is **parallel** to traction direction. No negative sign.

### 9.3 Tracing Through the Chain

**Tandem flow:**
1. `slip_rate()` returns `V_vec = -(V/|τ|)·τ` (sign: anti-parallel to τ)
2. `RateAndState::rhs()`: `r_mat(node, t) = Vi[t]` stores rate (dS/dt = V_vec)
3. State S is integrated: `S(t) = S(0) + ∫V_vec dt` (cumulative slip, negative of traction direction)
4. `ElasticityAdapter::slip()`: maps S to 3D via `copy_slip` and `fault_basis_q`
5. `set_slip()` sets `f_q` for the DG RHS

**MFEM flow:**
1. `SolveSlipRateVectorPsi()` returns `V_vec = (V/|τ|)·τ` (sign: parallel to τ)
2. `ComputeRHS()`: `rate(i*3 + c) = V_vec[c]` stores rate
3. State S integrated: `S(t) = S(0) + ∫V_vec dt` (opposite sign from Tandem)
4. `GetSlip()`: `slip(i*2+c) = state(i*3+c)` copies directly
5. `EmbedSlip()`: maps 2-component slip to 3D in fault basis
6. `AssembleSlipContributionIP()`: uses `sign * delta_u` as the prescribed jump

### 9.4 Tandem's Copy Convention

In `AdapterBase.cpp`:
```cpp
copy_slip = Identity with offset: copy_slip[n, n+1] = 1.0 for n=0..dim-2
```

For 3D (dim=3): `copy_slip` = `[[0,1,0],[0,0,1]]` (2×3 matrix)
This maps: tangential component 0 → physical dim 1, tangential component 1 → physical dim 2.

So if Tandem's slip state S has components `[S_0, S_1]` and the fault basis has columns
`[normal, tangent1, tangent2]`, the 3D slip vector at quad point q is:
```
slip_q[p] = Σ_l φ_l(q) × fault_basis[p, tangent1] × S[l,0] + fault_basis[p, tangent2] × S[l,1]
```

### 9.5 Is MFEM Self-Consistent?

The critical question: does the MFEM sign chain produce the correct `[[u]]` in the DG equation?

The DG equation enforces (in the IP formulation):
```
∫{σ(u)·n}·[[v]] + ε·∫[[u]]·{σ(v)·n} + ∫penalty·[[u]]·[[v]] =
∫{σ(u)·n}·g + ε·∫g·{σ(v)·n} + ∫penalty·g·[[v]]
```

where `g = [[u]]_prescribed` = the prescribed displacement jump.

If the fault slips in the +x direction on the +y side (Tandem convention for BP5), then:
- u(+y side) has +x displacement, u(-y side) has -x displacement
- `[[u]] = u(+y) - u(-y)` = positive in x → slip direction

Tandem's V_vec is anti-parallel to τ. Since τ is in the -x direction (fault resists the
+x loading), V_vec = -(-x) = +x direction. So dS/dt is in the +x direction. S is positive.
The 3D slip maps S → +x displacement. The jump `[[u]]` = +x. ✓

MFEM's V_vec is parallel to τ. τ is in the -x direction. V_vec = -x direction. dS/dt = -x.
S is negative. EmbedSlip maps S → -x displacement. Then `sign * delta_u` = ??? This depends
on the sign convention.

**This needs explicit verification on a running problem.** If the sign chain is wrong, the
prescribed jump would have the wrong sign, making the DG penalty push the solution in the
wrong direction. The effect would be a systematic traction error.

**Recommendation:** Add a diagnostic that prints `delta_u`, `sign`, `u_jump = u1 - u2`, and
`u_jump - sign*delta_u` at a known fault station. Compare with Tandem's values.

---

## 10. D11: Stress Formula — {ε} then σ vs {σ·n̂}

### 10.1 Tandem

Tandem's stress traction formula (from `elasticity.py:89-91`):

```python
def traction(x, normal):
    return lam_q[x]['q'] * Dx_q[x]['lsq'] * u[x]['ls'] * normal['pq'] +
           mu_q[x]['q'] * (Dx_q[x]['ljq']*u[x]['lp']*normal['jq'] +
                           Dx_q[x]['lpq']*u[x]['lj']*normal['jq'])
```

This computes `(σ(u) · n̂)_p` at each quad point `q`:
```
T_x[p,q] = λ(q) × ∇·u_x(q) × n̂[p,q] + μ(q) × (∇u_x[p,·]·n̂[·,q] + (∇u_x)^T[p,·]·n̂[·,q])
```

Then averages: `T_avg = 0.5 * (T_0 + T_1)`

### 10.2 MFEM

MFEM first computes average gradients, then stress, then dots with normal:

```cpp
// Average gradient
real_t ag = 0.5 * (grad1(ci, cj) + grad2(ci, cj));
real_t ag_t = 0.5 * (grad1(cj, ci) + grad2(cj, ci));
real_t eps_ij = 0.5 * (ag + ag_t);
// Stress
real_t stress_ij = lambda * tr(eps) * delta_ij + 2*mu * eps_ij;
// Traction
T_stress_q[ci] += stress_ij * basis.normal[cj];
```

This computes `avg(σ(u)) · n̂` where `avg(σ)` = `σ(avg(∇u))`.

### 10.3 Equivalence?

For constant material properties (λ, μ uniform):
```
avg(σ(u)·n̂) = 0.5*(σ(u_0) + σ(u_1))·n̂ = σ(0.5*(∇u_0 + ∇u_1))·n̂ = σ(avg(∇u))·n̂
```

So `{σ(u)·n̂} = σ({∇u})·n̂` for constant material properties. **No difference.**

For heterogeneous material properties (λ, μ vary by element):
```
{σ(u)·n̂} = 0.5*(σ_0(u_0) + σ_1(u_1))·n̂ ≠ σ_avg({∇u})·n̂
```

**Tandem handles heterogeneous materials correctly by using per-element `lam_q[x]`, `mu_q[x]`.**
MFEM uses a single `lambda_val_`, `mu_val_`. For BP5 (uniform material), both are identical.

---

## 11. Tandem's `surfaceOp` vs `compute_traction` — Internal Inconsistency?

A surprising finding: **Tandem's stiffness assembly and traction recovery use different normal
vectors.**

### 11.1 Stiffness Assembly (`surfaceOp`, `elasticity.py:105-110`)

```python
traction_q['pq'] <= 0.5 * (traction(0, n_q) + traction(1, n_q)),  # n_q (NOT unit)
u_jump['pq'] <= ...,
surface(x):
    unew[x] += c0[x] * traction_q * E_q[x] * w +
               c1[x] * tractionTest(x, u_jump) * w +
               c2[x] * w * E_q[x] * u_jump * nl_q   # penalty uses nl_q
```

The consistency term uses `n_q` (un-normalized) and the penalty uses `nl_q`.

### 11.2 Traction Recovery (`compute_traction`, `elasticity.py:242-244`)

```python
traction_q <= 0.5 * (traction(0, n_unit_q) + traction(1, n_unit_q)) +  # n_unit_q (UNIT)
              c00 * (u_jump - f)                                         # no nl_q
```

The stress average uses `n_unit_q` (unit normal) and the penalty has no `nl_q`.

### 11.3 Why This is Correct in Tandem

In the stiffness assembly, the integration is:
```
∫_F σ(u)·n ds = Σ_q w_q × σ_q × n_q
```
where `n_q = nl_q × n_unit_q` and `w_q` is the reference quadrature weight. So `w_q × n_q`
includes the face Jacobian (`nl_q`) in the integration measure. This is correct for
integrating the weak form.

In the traction recovery, the goal is to compute the **physical traction** `σ·n̂` at each
quad point (where `n̂` is the unit normal). This is then L2-projected with the physical
measure `nl_q × w_q` onto the fault basis. The `nl_q` enters through the projection,
not through the traction formula.

**This is mathematically consistent.** The stiffness assembly integrates `σ·n ds = σ·n̂ |J_f| ds_ref`
where `n_q = n̂ × |J_f|` absorbs the face Jacobian. The traction recovery computes `σ·n̂`
(physical traction) and handles `|J_f|` in the projection.

### 11.4 MFEM Comparison

MFEM's stiffness assembly (via `DGElasticityIntegrator`) also uses `nor` (un-normalized,
= `n_q` in Tandem notation) in the weak form integration. The traction recovery uses
`basis.normal` (unit normal) for the stress dot product. This matches Tandem's convention.

**However**, MFEM's `AssembleSlipContributionIP` uses `nor` (un-normalized) for the consistency
term and `nl_q * ip.weight` for the penalty weighting. This matches Tandem's `rhs_skeleton`
which uses `n_q` for consistency and `nl_q * w_q` for penalty.

**This appears consistent between MFEM and Tandem.**

---

## 12. D6: Quadrature Order Mismatch

### 12.1 Tandem

All surface operations use `MinQuadOrder()`, typically `2p+1`:

Source: `DGCurvilinearCommon.cpp:16`:
```cpp
auto fctRule = simplexQuadratureRule<DomainDimension - 1u>(MinQuadOrder());
```

### 12.2 MFEM

| Operation | Quad Order | Source |
|-----------|-----------|--------|
| K: DGElasticityIntegrator | `2p` | `bilininteg.cpp:4097` |
| K: DGElasticityIPPenaltyIntegrator | `2p` | `dg_elasticity_ip_penalty_integrator.hpp:76-78` |
| b: AssembleSlipContributionIP | `2p+1` (default) | `elasticity_operator.hpp:1193` |
| Traction: ComputeTractionImpl | `2p+1` (default) | `elasticity_operator.hpp:3761` |

### 12.3 Impact

v54 Section 9.5 tested changing K's integrators to order `2p+1`. **No effect** — for flat
elements at p=1, both orders exactly integrate degree-2 polynomial integrands.

**Not a correctness issue for linear elements at p=1.**

---

## 13. Solver and Time-Stepper Differences

### 13.1 Linear Solver

| | Tandem | MFEM |
|---|--------|------|
| Framework | PETSc KSP | MFEM BilinearForm + MUMPS/CG+AMG |
| Default | CG + p-multigrid (degenerates to CG+ILU at p=1) | MUMPS-BLR or CG+BoomerAMG |
| rtol | 1e-12 | 1e-10 (CG) or BLR tol 1e-12 |
| Matrix | Assembled P, matrix-free A (shell operator) | Fully assembled |

v54 proved MUMPS-BLR and CG+AMG give identical results. **Not the root cause.**

### 13.2 Time Stepper

| | Tandem | MFEM |
|---|--------|------|
| Method | PETSc TS `5dp` (Dormand-Prince RK45) | Custom DormandPrinceRK45 |
| atol | Configurable (default varies) | 1e-7 |
| rtol | Configurable | 1e-7 |
| Error norm | PETSc default (L2 or infinity) | Infinity norm on state |
| V-guard | None | Factor-100 per-step limit |

**Potential difference:** PETSc's TS `5dp` and MFEM's custom implementation may have
subtly different error estimation and step size control. PETSc computes the error
estimate and new step size using the standard embedded RK formula with configurable
norms. MFEM's implementation may use a different error norm or step size controller.

**This could cause trajectory divergence** — if MFEM's time stepper takes different step
sizes than Tandem's, the accumulated error can differ even if both are within tolerance.

---

## 14. Boundary Condition Handling

### 14.1 BP5 Geometry

BP5: Fault at Y=0, X ∈ [-Lx, Lx], Z ∈ [0, Lz].
Dirichlet: `u = (sgn(Y) × Vp × t/2, 0, 0)` on far-field boundaries.
Natural (traction-free): top (Z=Lz) and bottom (Z=0).

### 14.2 Tandem

Uses `Physical Surface(5)` for far-field Dirichlet, `Physical Surface(1)` for top/bottom
natural. The `SeasScenario` provides a `fun_boundary` function returning the Dirichlet
displacement as a function of position and time.

### 14.3 MFEM

Uses `BCMode::FarField`:
- Dirichlet on attrs 1-4 (x=±Lx, y=±Ly)
- Natural on attrs 5-6 (z=0, z=Lz)

`AssembleDirichletLoading` computes:
```cpp
u_D[0] = sign * Vp_ * time / 2.0;  // sign from Y centroid
```

### 14.4 Potential Difference

The attribute mapping must be verified: MFEM's geo file maps physical surfaces differently
from Tandem's geo file. This was addressed in earlier debug versions (H25 fix) but should
be double-checked that attrs 1-4 correctly tag the four far-field faces and attrs 5-6
correctly tag top/bottom.

Also: MFEM's `BuildDirichletInteriorFaces()` handles Y=0 interior faces outside the fault
region as special Dirichlet interior faces. Tandem doesn't need this because its mesh has
these as true boundary faces. **This is a structural difference that is hard to verify without
a mesh-specific test.**

---

## 15. BP5 Friction Parameters and Initialization

### 15.1 Tandem

BP5 parameters set in `.toml` config and `SeasScenario.cpp`:
- `a`, `b`: spatially varying (a-b < 0 in seismogenic zone)
- `V0 = 1e-6 m/s`
- `f0 = 0.6`
- `L (Dc)`: spatially varying
- `sigma_n = 25 MPa` (SnPre)
- `V_init = 0.01 m/s` (nucleation zone), `V_init = 1e-9 m/s` (elsewhere)

State initialized via: `psi_init = a * asinh(V_init / (2*V0) * exp(f0/b))` and
`theta_init = Dc/V0 * exp(psi/b - f0/b)` (from ageing law steady state?).

Actually, Tandem initializes `psi` directly from the steady-state friction at V_init:
```
psi = a * arcsinh(V_init / (2*V0) * exp(f0/a))  [regularized form]
```
Then the initial state is `[S=0, psi, ...]`.

### 15.2 MFEM

Same parameters from `bp5_scenario.hpp`. The initialization should match Tandem's formula.
**Verified in v53 to match at all stations after the 3-DOF fix.**

---

## 16. Prioritized Investigation Plan

Based on the differences found, the most likely root causes of the 0.38 MPa traction excess
are, in priority order:

### Priority 1: D8 — Slip Rate Sign Convention (Section 9)

**Why:** A sign error in the slip→displacement_jump chain would produce a systematic
traction error at every time step. The penalty correction `penalty * ([[u]] - g)` would
double instead of cancel if `g` has the wrong sign.

**Test:** Print `delta_u`, `sign`, `u1-u2`, and `(u1-u2) - sign*delta_u` at a specific
fault station on the first few time steps. Compare with Tandem's `u_jump - f_q`.

### Priority 2: D4 — Log10 vs Linear Slip Rate Solver (Section 6)

**Why:** At V ~ 0.01 m/s the linear solver is fine, but as V drifts to lower values
(nucleation failing → interseismic), the precision loss could create a feedback loop:
imprecise V → wrong dS/dt → wrong slip → wrong traction → more V decay.

**Test:** Switch MFEM to log10 Brent search and re-run the first 600 steps to see if
the V trajectory changes.

### Priority 3: D7 — RHS Slip Consistency Term Audit (Section 8)

**Why:** Even if the formula looks correct, the manual implementation with adjugate Jacobians,
sign conventions, and weight factors is complex enough that a subtle factor-of-2 or sign
error could hide. The best test is a numerical comparison.

**Test:** For a known displacement field (e.g., manufactured solution), verify that MFEM's
assembled RHS matches the expected DG right-hand-side from the formula.

### Priority 4: D11 — Stress Average at Traction Recovery (Section 10)

**Why:** The stress formula `σ({∇u})·n̂` vs `{σ(u)·n̂}` is equivalent for uniform material,
so this is NOT a difference. But verify the physical gradient computation (Jinv × dshape)
matches Tandem's `g × Dxi_q` in the non-trivial multi-element case.

### Lower Priority

D2/D3 (L2 projection weights), D5 (split integrators), D6 (quad order), D9 (sigma_n clamp),
D10 (fault basis per-face) are all verified to produce no difference for flat linear tets
at p=1 with uniform material.

---

## 17. What Has Been Eliminated (Running Total from v53/v54)

| Hypothesis | Evidence | Version |
|------------|----------|---------|
| MUMPS-BLR solver error | CG+AMG identical | v54 |
| Penalty coefficient formula | Matches Tandem exactly | v54 |
| Slip interpolation bug | Exact round-trip verified | v54 |
| L2 projection (mass matrix) | M_ref correct | v54 |
| Quadrature mismatch (K vs b) | Both orders give identical K for flat p=1 | v54 |
| Dirichlet BC time/sign | Correct time behavior verified | v54 |
| State evolution formula | Matches Tandem | v54 |
| Friction solver (residual) | Both converge to < 1e-16 | v54 |
| elastic_sigma_n OFF | Now ON + sign fixed, no effect on nucleation | v54 |
| Fault discretization (1 DOF/face) | Fixed to 3 DOF/face | v53 |
| Nucleation zone geometry | bp5_outside eps=1e-3 matches Tandem | v53 |
| L2 projection nl_q weighting (D2/D3) | Cancels for flat faces | v55 |
| Stress formula {ε} vs {σ·n} (D11) | Equivalent for uniform material | v55 |
| Normal vector in traction (D1) | Identical for flat faces | v55 |

---

## 18. Open Questions

1. **Is the slip sign chain (D8) self-consistent?** Trace V_vec → dS/dt → S → embed → sign×delta_u → [[u]] = g^F and verify the final `g^F` has the correct sign.

2. **Does the log10 solver (D4) change the trajectory?** Quick test: switch to log10 Brent in MFEM and run 600 steps.

3. **Is MFEM's RHS consistent with Tandem's for a given displacement field?** Manufacture a displacement, compute RHS in both codes, compare.

4. **How does Tandem's penalty correction / stress ratio compare to MFEM's 50x?** This was identified as the key open question in v54 Section 9.4.

---

## 19. Code Locations Reference

### Tandem
| Component | File | Key Lines |
|-----------|------|-----------|
| DG stiffness (generated) | `app/kernels/elasticity.py` | 97-110 (surfaceOp) |
| DG RHS (generated) | `app/kernels/elasticity.py` | 152-178 (rhsFacet) |
| Traction recovery (generated) | `app/kernels/elasticity.py` | 242-248 (compute_traction) |
| Penalty coefficient | `app/localoperator/Elasticity.cpp` | 278-291, `InverseInequality.h` |
| Traction driver | `app/localoperator/Elasticity.cpp` | 948-987 (traction_skeleton) |
| Adapter projection | `app/kernels/elasticity_adapter.py` | 26-27 (evaluate_traction) |
| Friction solver | `app/localoperator/DieterichRuinaBase.h` | 81-175 (slip_rate) |
| State evolution | `app/localoperator/DieterichRuinaAgeing.h` | state_rhs |
| ODE driver | `app/form/SeasQDOperator.cpp` | evaluate |
| Time stepper | `app/common/PetscTimeSolver.h` | 5dp (DOPRI45) |

### MFEM
| Component | File | Key Lines |
|-----------|------|-----------|
| DG stiffness | `domain/elasticity_operator.hpp` | 910-970 (AssembleStiffness) |
| DG RHS (slip) | `domain/elasticity_operator.hpp` | 1140-1320 (AssembleSlipContributionIP) |
| DG RHS (Dirichlet) | `domain/elasticity_operator.hpp` | 2014-2260 (AssembleDirichletLoading) |
| Traction recovery | `domain/elasticity_operator.hpp` | 3700-4060 (ComputeTractionImpl) |
| Friction solver | `friction/dieterich_ruina.hpp` | 324-378 (SolveSlipRatePsi) |
| State evolution | `friction/state_evolution.hpp` | AgingLawPsi::Rate |
| ODE driver | `solver/seas_operator.hpp` | ComputeRHS |
| Time stepper | `solver/time_stepper.hpp` | DormandPrinceRK45 |
| Penalty integrator | `integrator/dg_elasticity_ip_penalty_integrator.hpp` | AssembleFaceMatrix |
| L2 projection | `fault/face_quadrature.hpp` | GalerkinProject |
