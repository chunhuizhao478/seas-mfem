# BP5 Debug v60: Normal Stress Oscillation at VW/VS Boundary

## 1. Observation

ParaView fault surface output (v60 ParaView run) reveals **+-5 MPa oscillatory normal stress** at the VW/VS parameter boundary at t = 2.6e7 s (~ 0.82 yr). Key facts:

- **Tangential fields** (slip, slip rate, traction dip/strike) are **smooth** — no oscillation
- **Normal stress only** shows the +-5 MPa checkerboard pattern at the VW/VS boundary
- **Tandem reference** shows < 0.3 MPa normal stress variation even at 25 yr
- The VW/VS boundary is where friction parameter `a` transitions from 0.004 (VW) to 0.04 (VS) over a 2 km zone
- This normal stress oscillation likely feeds into the friction law via sigma_n_eff = 25 MPa + normal_traction, eventually causing the ~25.27 yr blowup

## 2. Complete Difference Catalog: MFEM vs Tandem

### 2.1 Penalty Parameter Computation

**Tandem** (`Elasticity.cpp:487-499`):
```cpp
auto const p = [&](int side) {
    const auto [c0, c1] = stiffness_tensor_bounds(info.up[side]);
    constexpr double c_N_1 = InverseInequality<Dim>::trace_constant(PolynomialDegree - 1);
    return (Dim + 1) * c_N_1 * (area_[fctNo] / volume_[info.up[side]]) * (c1 * c1 / c0);
};
penalty_[fctNo] = (p(0) + p(1)) / 4.0;  // single scalar per face
```

**MFEM** (`dg_elasticity_ip_combined_integrator.hpp:827-846`):
```cpp
int p = std::max(fe1.GetOrder(), fe2.GetOrder());
real_t c_N_1 = p * (p + dim_ - 1.0) / dim_;
real_t c0 = 2.0 * mu_val;   // <-- from lam1t, mu1t (Element 1 only!)
real_t c1 = dim_ * lam + 2.0 * mu_val;
real_t ratio = c1 * c1 / c0;
real_t p0 = (dim_ + 1) * c_N_1 * (real_t(dim_) * nl_q / detJ1) * ratio;
real_t p1 = (dim_ + 1) * c_N_1 * (real_t(dim_) * nl_q / detJ2) * ratio;
return penalty_factor_ * (p0 + p1) / 4.0;  // recomputed per QP
```

#### Differences:

| # | Aspect | Tandem | MFEM | Impact for BP5 order 1 uniform mat |
|---|--------|--------|------|------------------------------------|
| D1 | Penalty evaluation | Single precomputed scalar per face | Recomputed at each face QP | Same for affine tets (constant Jacobians) |
| D2 | Material for p(side) | `stiffness_tensor_bounds(info.up[side])` — each side uses **its own** element's material | `lam1t, mu1t` — **both sides** use Element 1's material | Same for uniform lambda, mu |
| D3 | Material values | min(2*mu) / max(D*lam+2*mu) over all element QPs | Single point value at face QP | Same for uniform lambda, mu |
| D4 | Area/volume ratio | `area_[fctNo] / volume_[info.up[side]]` — integrated physical area / volume | `dim * nl_q / detJ` at single QP | Same for affine tets: dim * nl_q/detJ = area_phys/vol_phys |
| D5 | penalty_factor_ | Not present | Multiplier, default 1.0 (CLI `--penalty-factor`) | Same at default 1.0 |
| D6 | Trace constant arg | `trace_constant(PolynomialDegree - 1)` | `p*(p+dim-1)/dim` where p=fe.GetOrder() | Same: both = (N+1)(N+D)/D with N=p-1, giving 1.0 for p=1 |

### 2.2 Traction at Quadrature Points (ComputeTractionAtQuadPoints)

**Tandem** (`elasticity.py:242-244`, `Elasticity.cpp:1246-1284`):
```python
traction_q <= 0.5 * (traction(0, n_unit_q) + traction(1, n_unit_q))
            + c00 * (E_q[0]*u[0] - E_q[1]*u[1] - f_q)
# c00 = -penalty(fctNo)  (single scalar per face)
```

**MFEM** (`dg_elasticity_ip_combined_integrator.hpp:500-538`):
```cpp
T_p = 0.5 * (sig1 + sig2) * n_hat(j);          // stress average with unit normal
traction_q(c*nq+q) += (-penalty) * (u1q - u2q - f_q);  // penalty per QP
```

#### Differences:

| # | Aspect | Tandem | MFEM | Impact |
|---|--------|--------|------|----|
| D7 | Stress formula | Generated kernel: `lam_q*Dx_q*u * n_unit_q + mu_q*(Dx_q*u*n_unit_q + ...)` | Explicit loop: sig_ij = lam*tr(eps)*delta_ij + 2*mu*eps_ij, then T_p = 0.5*(sig1+sig2)*n_hat(j) | Mathematically identical |
| D8 | Penalty value | `c00 = -penalty(fctNo)` — precomputed, single per face | `(-penalty)` recomputed at each QP via ComputePenalty | Same for affine tets, uniform mat (see D1) |
| D9 | Side 0/1 convention | `E_q[0]*u[0] - E_q[1]*u[1]` where up[0]/up[1] from FacetInfo | `u1q - u2q` where Elem1/Elem2 from FaceElementTransformations | Verified matching at first step in v59 (traction matches 8-11 digits) |
| D10 | Physical gradients | `Dx_q = g * Dxi_q` where g = J_inv precomputed per face | `dshape_phys = dshape_ref * Jinv` computed per QP | Mathematically identical for affine tets |

### 2.3 L2 Projection (ProjectTractionToFaultDOFs)

**Tandem** (`elasticity_adapter.py:24-27`, `AdapterBase.cpp:87-96`):
```python
# Mass matrix: M[i,j] = sum_q w[q] * nl_q[q] * e_q[i,q] * e_q[j,q]
# Projection: traction[k,p] = Minv[l,k] * e_q_T[q,l] * w[q] * nl_q[q] * traction_q[o,q] * fault_basis_q[o,p,q]
```

**MFEM** (`dg_elasticity_ip_combined_integrator.hpp:714-804`):
```cpp
M(i,j) += w * nl_q(q) * e_q(i,q) * e_q(j,q);       // same mass matrix
rhs(l,t) += wn * e_q(l,q) * T_local;                 // T_local = traction_q . basis_vec[t]
traction_local(t*nbf+k) = sum_l Minv(k,l) * rhs(l,t);
```

#### Differences:

| # | Aspect | Tandem | MFEM | Impact |
|---|--------|--------|------|----|
| D11 | Projection formula | Single fused kernel: Minv * E_q_T * w * nl * T_q * fault_basis | Explicit loops: build M, invert, build RHS, multiply | Mathematically identical |
| D12 | Mass matrix symmetry | Uses `EigenMap(m).inverse()` | Uses `DenseMatrixInverse(M).GetInverseMatrix(Minv)` | Both compute M^{-1} of same matrix |
| D13 | Fault basis vectors | `fault_basis_q[o,p,q]` — precomputed per-QP from Curvilinear::facetBasis() | `FaultBasisQPData::normal/tangent1/tangent2` — per-QP from BuildFaultGeometry | Verified matching in v59 (bbae289 sign fix) |

### 2.4 Stiffness Matrix Assembly

**Tandem** (`elasticity.py:105-110`):
```python
# surfaceOp uses n_q (un-normalized normal) for traction
traction_q <= 0.5 * (traction(0, n_q) + traction(1, n_q))
# surface(x) uses c0[x], c1[x], c2[x] per-side coefficients
```

**MFEM** (`dg_elasticity_ip_combined_integrator.hpp:115-166`):
```cpp
// AssembleFaceBlock uses nor (un-normalized) with dshape_adj / detJ
// Explicit symmetrization of face matrix (lines 173-181)
```

#### Differences:

| # | Aspect | Tandem | MFEM | Impact |
|---|--------|--------|------|----|
| D14 | Consistency term normal | `traction(x, n_q)` — un-normalized | `dshape_adj * nor / detJ` — equivalent to `dshape_phys * n_hat * nl` | Mathematically identical |
| D15 | Explicit symmetrization | None — SIPG is symmetric by construction | `elmat(i,j) = 0.5*(elmat(i,j)+elmat(j,i))` at lines 173-181 | Ensures exact machine-precision symmetry. Tandem doesn't need this (uses MUMPS direct solver). MFEM also uses MUMPS but symmetrizes anyway |
| D16 | Penalty in matrix | `penalty_[fctNo]` — single per face | `ComputePenalty(...)` per QP | Same for affine tets, uniform mat (see D1) |

### 2.5 Quadrature Rules

**Tandem**: `fctRule = simplexQuadratureRule<2>(MinQuadOrder())` where `MinQuadOrder() = 2*PolynomialDegree+1 = 3` for p=1.

**MFEM**: `IntRules.Get(TRIANGLE, 2*p+1)` where `p = max(fe1.GetOrder(), fe2.GetOrder()) = 1`, so order 3.

| # | Aspect | Tandem | MFEM | Impact |
|---|--------|--------|------|----|
| D17 | Quadrature order | `2*PolynomialDegree + 1 = 3` | `2*p + 1 = 3` | Same order |
| D18 | QP locations/weights | Tandem's simplexQuadratureRule (custom library) | MFEM's IntRules (Dunavant/Stroud) | **Potentially different point locations and weights for the same order.** Both integrate polynomials of degree 3 exactly, but the specific points differ between implementations. For affine tets with polynomial basis, this doesn't affect the integral. But if there are non-polynomial features (e.g., oscillations), different sampling points yield different intermediate values. |

### 2.6 Traction Projection Quadrature (Adapter)

**Tandem**: The projection `evaluate_traction` kernel uses the SAME `fctRule` QPs as the traction computation.

**MFEM**: `ProjectTractionToFaultDOFs` receives the same `IntegrationRule &ir` used by `ComputeTractionAtQuadPoints`.

| # | Aspect | Tandem | MFEM | Impact |
|---|--------|--------|------|----|
| D19 | Projection QP rule | Same as traction computation | Same as traction computation | Consistent in both codes |

### 2.7 PETSc Time Stepping Configuration

**Tandem** (`examples/tandem/3d/bp5-qd/rk45.cfg`):
No `-ts_adapt_dt_max` option.

**MFEM** (`tests/verification/petsc_ts_rk45_tandem.cfg`):
`-ts_adapt_dt_max 3.15576e6` (= 0.1 yr)

| # | Aspect | Tandem | MFEM | Impact |
|---|--------|--------|------|----|
| D20 | dt_max constraint | None — PETSc default (no maximum) | 0.1 yr max step size | **Could affect how the time integrator resolves the normal stress evolution. Tandem may take larger steps during interseismic period.** |

## 3. Assessment

For BP5 with order-1 affine tets and uniform elastic material (lambda, mu constant everywhere):

- **Differences D1-D6, D8, D10, D14, D16**: All involve per-QP vs per-face evaluation, or per-side vs single-side material. For affine tets with uniform material, these are **mathematically equivalent** — but the code paths differ.

- **Difference D15 (explicit symmetrization)**: MFEM averages `elmat(i,j)` and `elmat(j,i)`. For exact SIPG with uniform material, `elmat` should already be symmetric, so this averaging changes nothing beyond O(machine epsilon). Could accumulate over very many steps.

- **Difference D18 (QP locations)**: Different quadrature implementations may use different triangle QP sets for the same polynomial order. Both integrate the bilinear form exactly, so K and b are the same. But intermediate per-QP values (used in traction computation) are sampled at different spatial locations.

- **Difference D20 (dt_max)**: This is a **configuration difference** that could affect the time-stepping behavior. Tandem may take larger interseismic steps, potentially "stepping over" short-lived oscillations that MFEM resolves due to the 0.1 yr cap.

## 4. Key Question: Where Does the Oscillation Enter?

Since all formulas are mathematically equivalent, the oscillation must arise from one of:

1. **A subtle numerical difference** that gets amplified through the feedback loop:
   normal jump error --> penalty amplification --> sigma_n oscillation --> friction response --> slip correction --> displacement update --> larger normal jump error

2. **The dt_max difference (D20)**: Tandem may take ~1 yr steps during interseismic, skipping over oscillation growth. MFEM with dt_max=0.1 yr takes smaller steps, allowing the oscillation to develop.

3. **An actual code bug** in a path not covered by the first-step comparison (e.g., a path only active after many steps, or a path specific to the normal component).

## 5. Proposed Diagnostics

### 5.1 Remove dt_max (test D20)

Remove `-ts_adapt_dt_max 3.15576e6` from `petsc_ts_rk45_tandem.cfg` and rerun. If the oscillation disappears or is reduced, the dt_max constraint is forcing MFEM to resolve dynamics that Tandem skips.

### 5.2 Dump normal traction evolution at VW/VS boundary

Add diagnostic output of the **normal component** of `traction_q` (at QPs, BEFORE L2 projection) for a specific fault face at the VW/VS boundary, at steps 0, 1, 10, 100, 1000. Compare:
- Stress average part: `0.5*(sig1+sig2) . n_hat` projected onto fault normal
- Penalty correction: `-penalty * (u1_n - u2_n)` (the normal displacement jump)

This identifies whether the oscillation enters through the stress average or the penalty term.

### 5.3 Compare penalty values

At the first time step, dump the penalty value for each fault face and compare with Tandem's precomputed `penalty_[fctNo]`. Even though the formula should give the same result, verify numerically.

### 5.4 Reduce penalty (test feedback amplification)

Run with `--penalty-factor 0.5` to halve the penalty. If the normal stress oscillation is reduced proportionally, the penalty term is the amplifier.
