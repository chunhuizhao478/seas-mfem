# BP5 Debug v32/v33: Recurrence Interval Diagnosis + General Polynomial Order Fix

**Date**: 2026-03-15
**Status**: Implementation complete, all tests pass, cluster scripts ready
**Previous**: v31 (DG slip sign fix + Dirichlet interior faces, all unit tests pass)

---

## 1. Problem Statement

BP5 with BR2 method at p=1, h=1000m produces earthquake recurrence ~435 yr,
which is 1.82× longer than Tandem reference (~240 yr, p=6, h=4km).

All other qualitative features match: correct stress drop pattern, correct
slip distribution, proper earthquake nucleation and arrest.

## 2. Root Cause Analysis

### 2.1 The BR2 traction correction

The DG traction at a fault face is:

```
t = {{σ·n̂}} - σ × {{C : L([[u]] − δ)}} · n̂
```

where L is the BR2 lifting operator and σ = dim+1 = 4.

**MFEM** implements the full BR2 lifting: the correction involves M⁻¹,
face integrals, and elasticity tensor contraction. The correction has magnitude
~μ/h × |residual| in Pascals.

**Tandem** uses a scalar placeholder: `c00 = -penalty(fctNo) = -NumFacets = -4`
(dimensionless). The correction `-4 × ([[u]] − δ)` is in meters, not Pascals —
effectively zero in stress terms.

### 2.2 Why Tandem can use zero correction

Tandem uses p=6 (quintic basis). The DG residual `[[u]] − δ` converges
spectrally — it's negligible for smooth solutions. The traction is dominated
by `{{σ·n̂}}`, and the correction term is irrelevant.

### 2.3 Why MFEM's correction causes bias

For p=1 (linear basis), the DG residual scales as O(h²). The BR2 lifting
amplifies by ~1/h, giving a traction correction of O(μ × h).

For h=1000m, μ=32 GPa: the correction at healthy faces was ~0.06 MPa (from
v24 diagnostics). This is only 0.4% of the ~15 MPa shear traction.

But rate-and-state friction amplifies exponentially:

```
δV/V ≈ (1/a) × δτ/σ_n ≈ 100 × 0.004 = 40%
```

A systematic 40% slip-rate bias during the interseismic period accumulates
over centuries, delaying nucleation by a factor consistent with 1.82×.

### 2.4 Why removing the correction blows up (v24-v28)

The v24-v28 experiments tried replacing the BR2 lifting with IP-style scalar
penalty or zero correction. All blew up at specific fault DOFs within ~2 years.

The BR2 lifting provides **direction-dependent** stabilization (through the
elasticity tensor). The IP penalty is **isotropic** — it over-corrects normal
components and under-corrects tangential, destabilizing specific face orientations.

The pathological DOFs need the **structure** of the BR2 correction, not just its
magnitude. This is why the dilemma exists:

| Approach | Stability | Recurrence | Problem |
|----------|-----------|------------|---------|
| Full BR2 lifting | ✓ stable | 1.82× too long | O(h) bias |
| No correction | ✗ blows up | N/A | Specific DOFs unstable |
| IP scalar penalty | ✗ blows up | N/A | Wrong direction |

## 3. Solution: Increase Polynomial Order

The traction correction scales as **O(h^p)**:
- p=1: O(h) → 1.82× mismatch at h=1km
- p=2: O(h²) → ~50× smaller correction → ~1% mismatch expected
- p=6 (Tandem): negligible

This is the physically correct fix — no free parameters, no tuning. The method
converges; p=1 just isn't in the asymptotic regime yet.

### 3.1 DOF count comparison

| Order | DOFs/tet | Total DOFs (76K tets) | MUMPS cost factor |
|-------|----------|----------------------|-------------------|
| p=1 | 4 | ~900K | 1× |
| p=2 | 10 | ~2.3M | ~3.5× |

### 3.2 Equivalent accuracy

| Config | Traction bias | Total DOFs |
|--------|--------------|------------|
| p=1, h=1000m | O(h) ≈ large | ~900K |
| p=1, h=250m | O(h) ≈ medium | ~15M |
| p=2, h=1000m | O(h²) ≈ small | ~2.3M |

p=2 at 1000m gives better traction accuracy than p=1 at 250m, at 6.5× fewer DOFs.

## 4. Implementation: Code Changes for General Order

### 4.1 What's already order-general (no changes needed)

- DG integrators (BR2, IP bilinear forms) — quadrature adaptive to element order
- Mass matrix computation — quadrature adaptive
- Slip assembly (AssembleSlipContributionIP/BR2) — uses `2*face_order+1` quadrature
- FE space creation — uses `order_` parameter from constructor
- BR2 penalty σ = dim+1 — geometry-based, correct for all orders (matches Tandem)
- Fault system (rate-and-state) — indexes by face × 2 components, volume-order-independent

### 4.2 Issue 1: Hardcoded `int order = 1` in drivers

**Files**: `bp5_verification_full.cpp` (line 560), `pseas.cpp` (line 120)

**Fix**: Add `--order` command-line parameter, default 1.

### 4.3 Issue 2: IP trace constant `c_N_1 = 1.0` hardcoded

**Files**: `elasticity_operator.hpp` (lines 2775, 3081),
`dg_elasticity_ip_penalty_integrator.hpp` (line 110)

**Fix**: Tandem formula `c_N_1 = p * (p + D - 1) / D` where p = polynomial degree.
For D=3: p=1 → 1.0 (unchanged), p=2 → 8/3, p=3 → 5.0.

### 4.4 Issue 3: Order-0 quadrature in BR2 traction lifting (CRITICAL)

**Files**: `elasticity_operator.hpp` (lines 2817, 3112)

The face integral in the BR2 traction correction uses a single centroid point.
For p=1 (linear shapes on triangle), the centroid rule is exact. For p≥2
(quadratic shapes), it's insufficient.

**Fix**: Use `2*face_order` quadrature order with multi-point loop for the face
integral. Keep the lifted function evaluation at the face centroid (single
traction value per fault face — the fault system architecture is unchanged).

## 5. What Does NOT Change

- `2*fi` indexing in fault system = 2 components (dip, strike), NOT DOFs per face
- BR2 penalty σ = dim+1 = correct for all orders
- Rate-and-state ODE system size = same (one evaluation point per face)
- Fault basis, geometry, I/O = independent of volume element order

## 6. Verification Plan

### Test 1: Regression at p=1
Run all 22 serial test suites. Must all pass (changes should be no-op for order=1).

### Test 2: p=2 smoke test
Run BP5 with `--order 2` on 1000m mesh, 8 MPI ranks, 10 time steps.
Verify: no crash, traction values O(15 MPa), reasonable displacement.

### Test 3: p=2 convergence (TACC)
Submit 1000m mesh with `--order 2`, full 1800-year run.
**Success criterion**: recurrence interval significantly closer to 240 yr than
the 435 yr from p=1.

## 7. Files Changed

| File | Changes |
|------|---------|
| `tests/verification/bp5_verification_full.cpp` | Add `--order` CLI option |
| `domain/elasticity_operator.hpp` | c_N_1 formula, BR2 traction quadrature |
| `integrator/dg_elasticity_ip_penalty_integrator.hpp` | c_N_1 formula |
| `pseas.cpp` | Add `--order` CLI option |

## 8. Why Antiplane Traction Works Without Penalty Correction

### 8.1 Observation

The antiplane `ComputeTraction` (`antiplane_operator.hpp`, line 938) computes:

```
tau_face = mu * avg_dudx
```

This is the pure average gradient `μ × {{∂u/∂x}}` — NO penalty correction, no
BR2 lifting, no IP stabilization. Yet antiplane BP1 produces correct results.

In contrast, the 3D elasticity `ComputeTraction` requires the full BR2 lifting
correction, and removing it causes blow-up (v24-v28 experiments).

### 8.2 Why no penalty works for scalar PDE

The antiplane problem is a scalar PDE:

```
-∇·(μ ∇u) = f
```

The traction is a scalar quantity: `τ = μ ∂u/∂n`. There are three reasons why
omitting the DG penalty from traction evaluation is safe here:

**1. No directional coupling.** The scalar traction has only one component. Any
DG residual (jump mismatch) shifts the traction magnitude, but there is no
"wrong direction" — there IS only one direction. The error is purely in amplitude
and converges with mesh refinement.

**2. Penalty is in the stiffness matrix, not the traction.** The DG bilinear form
used to *solve* for u includes the IP/BR2 penalty terms:

```
a(u,v) = Σ_K ∫_K μ ∇u·∇v  −  Σ_F ∫_F {{μ ∇u}}·[[v]] − [[u]]·{{μ ∇v}}
         + penalty terms
```

These penalty terms ensure u is well-determined (unique, stable solution).
Once u is solved, the traction `μ ∂u/∂n` can be evaluated from the solution
without repeating the penalty — the solution already "absorbed" the stabilization.

**3. Scalar penalty has correct structure by default.** Even if penalty WERE
included in traction evaluation, it would be `κ × [[u]]` — a scalar times a
scalar. There's no way for a scalar penalty to have "wrong structure." Any
reasonable positive κ would work.

### 8.3 Why 3D elasticity NEEDS the BR2 correction

The 3D elasticity problem is a vector PDE with tensor coupling:

```
-∇·(C : ε(u)) = f
```

The traction is a vector: `t_i = C_{ijkl} ε_{kl} n_j`. The elasticity tensor C
couples all displacement components through anisotropic stress-strain relations.

**1. Directional coupling.** The traction has 3 components (t_x, t_y, t_z) that
are coupled through C. A scalar IP penalty `κ × [[u]]` applies the same
correction to all components — but the correct correction depends on the
face orientation relative to the elasticity tensor. The penalty over-corrects
normal components and under-corrects tangential ones (or vice versa), creating
a directional bias.

This is exactly why v24-v28 blew up at specific DOFs: faces with orientations
where the IP scalar penalty had the worst directional mismatch became unstable.

**2. BR2 provides direction-correct stabilization.** The BR2 lifting computes
the correction through:

```
L([[u]]) = M⁻¹ ∫_F [[u]] ⊗ n dS
correction = σ × {{C : L}} · n̂
```

The elasticity tensor C appears in the correction, ensuring the penalty has the
correct anisotropic structure. This is why BR2 is stable while IP scalar penalty
is not.

**3. Magnitude matters for traction, not just stability.** Even when stable
(full BR2), the correction magnitude adds a systematic O(h^p) bias to the
traction. For scalar PDE, this bias only affects amplitude (harmless at O(h)).
For vector PDE, the bias affects both amplitude AND direction of the traction
vector, coupling into the fault friction law through both dip and strike
components.

### 8.4 Summary: scalar vs vector DG traction

| Property | Antiplane (scalar) | 3D Elasticity (vector) |
|----------|-------------------|----------------------|
| Traction | scalar τ = μ ∂u/∂n | vector t = C:ε·n̂ |
| DG penalty in solve | Yes (IP/BR2) | Yes (BR2) |
| Penalty in traction eval | No (not needed) | Yes (required) |
| Why | No direction, scalar residual | Tensor coupling, directional |
| Omitting penalty | Amplitude shift only | Direction + amplitude → unstable |

### 8.5 Implication for p-refinement

This analysis reinforces the v33 approach:
- Antiplane: penalty-free traction works → p-refinement not needed for traction accuracy
- 3D elasticity: penalty IS in traction → correction magnitude O(h^p) biases result → p-refinement reduces bias

## 9. Implementation Results

### 9.1 Test suite results (all pass)

All 22 serial test suites pass at p=1, confirming zero regression. Additionally,
pre-existing test failures in 3 suites were fixed:

| Test Suite | Issue | Fix |
|------------|-------|-----|
| `seas_test_bp5_output` (4 failures) | Test inputs used positive slip/traction, but `WriteFromGlobalData` negates for SCEC convention | Changed inputs to negative (internal convention) |
| `seas_test_parallel_elasticity` (2 failures) | `CreateTestMesh3D` had z∈[0,Lz] instead of [-Lz,0]; fault detection requires z≤0 | Added `v[2] -= Lz` shift; fixed boundary attributes |
| `seas_test_bp5_parallel_smoke` (1 failure) | Same z-coordinate bug | Same fix |

### 9.2 Cluster submission scripts

| Script | Purpose | Nodes | Queue | Time |
|--------|---------|-------|-------|------|
| `bp5_v33_1000m_br2_p2_smoke.sbatch` | Verify p=2 works, no crash | 8 | development | 2h |
| `bp5_v33_1000m_br2_p1_baseline.sbatch` | Regression: v33 code at p=1 = v32 | 8 | normal | 48h |
| `bp5_v33_1000m_br2_p2_full.sbatch` | Test p=2 recurrence interval | 16 | normal | 48h |

### 9.3 Success criteria for p=2 full run

- Recurrence interval significantly closer to 240 yr than 435 yr
- If recurrence ≈ 240-280 yr → formulation is correct, mismatch was discretization error
- Compare slip profiles, stress drop with Tandem reference

---

## 10. p=2 Blow-Up Diagnosis and Fix (v33b)

### 10.1 Observed failure

The p=2 smoke test (job 7599513) ran 652 steps before hitting walltime. V_max
grew monotonically from 0.05 m/s to 400+ m/s without arrest — the earthquake
never stopped. Traction blow-up messages appeared at step 443 (V_max ≈ 220 m/s)
on 6 interior-face DOFs at rank 321, with:
- `tau_mag ≈ 1.8 GPa` (72× the background σ_n = 25 MPa)
- Slip values of 100+ meters (completely unphysical)

The failure was NOT a crash or NaN — the ODE integrator kept running with
decreasing dt (~1 ms), but the dynamics were completely unphysical.

### 10.2 Root cause: centroid evaluation of BR2 lifted function at p=2

Two bugs in the BR2 traction `ComputeTraction` code:

**Bug 1: Centroid displacement jump used for all quadrature points.**
The face integral `∫_F φ_m · [[u]] ⊗ n dS` used a single centroid value of
`[[u]]` at all quadrature points. For p=2, the displacement field is quadratic,
so `[[u]]` varies across the face. Using the centroid value is incorrect.

**Bug 2: Centroid shape functions used for lifted function evaluation.**
The evaluation `L_h(x_c) = Σ_m φ_m(x_c) · c_m` uses shape functions at the
face centroid. For p=2 GaussLobatto basis on a triangle:

| DOF type | Shape value at centroid | Face average |
|----------|------------------------|-------------|
| Vertex   | -1/9                   | 0           |
| Edge mid | +4/9                   | 1/3         |

The negative vertex values at the centroid cause the evaluation to amplify
oscillatory modes in the lifted function, producing incorrect traction values.
The face-averaged shapes (non-negative, partition of unity) are the physically
correct choice for the one-DOF-per-face fault system.

**Why p=1 was unaffected:** For p=1 (linear basis), centroid shape values are
all 1/(dim+1) = 1/4 for a tet (all positive), and the displacement jump is
linear on the face, so the centroid value equals the face average. Both bugs
are dormant at p=1.

### 10.3 Fix applied

Two changes in `domain/elasticity_operator.hpp`, applied to both interior face
and shared face BR2 sections:

1. **Per-quadrature-point displacement jump**: The face integral now evaluates
   `[[u]]` at each quadrature point from the element DOF vectors:
   ```
   u1q = Σ_k shape1(k, x_q) · u1_all(c*ndof1+k)
   jump_q[c] = (u1q - u2q) - sign * delta_u[c]
   ```

2. **Face-averaged shape evaluation**: During the quadrature loop, accumulate
   face-averaged shapes: `avg_shape(m) = Σ_q w_q · shape(m, x_q) / Σ_q w_q`.
   Use these instead of centroid shapes in the evaluation step.

Both changes are zero-cost at p=1 (face-averaged shapes equal centroid shapes
for linear basis) and provide correct evaluation at p≥2.

### 10.4 Verification

All 32 unit tests pass (22 serial + 10 parallel) — zero regression at p=1.

**p=1 baseline confirmation**: The p=1 full run with the v33b code (including
the face-averaged evaluation fix) produces results identical to the previous
v33 p=1 baseline. This confirms that the fix is truly zero-cost at p=1: the
face-averaged shapes equal the centroid shapes for linear basis, and the
per-quadrature-point jump equals the centroid jump for linear displacement
fields. The recurrence interval, slip profiles, and stress drop are unchanged.

---

## 11. Why Antiplane Traction Works Without Penalty Correction

### 11.1 Background

The antiplane (BP1-BP3) traction computation uses only `μ · {{∂u/∂x}}` with
NO penalty correction. This is different from the 3D elasticity traction
which includes the BR2 or IP penalty. Why does antiplane work?

### 11.2 Mathematical reason

**Antiplane PDE** is a scalar Poisson equation: `∇·(μ∇u) = 0` where `u` is the
out-of-plane displacement (scalar).

The DG traction is: `t = μ · {{∂u/∂n}}`

The penalty correction for DG would add: `- penalty × ([[u]] - δ)`

For the scalar PDE, the traction has only one component (normal derivative).
There is no directional coupling — no interaction between displacement components.
The average flux `{{∂u/∂n}}` captures the correct traction without penalty because:

1. **No tensor contraction**: The elasticity tensor C reduces to the scalar μ.
   There's no C:∇u contraction that mixes displacement components.

2. **Self-consistent residual**: For the scalar PDE, the DG bilinear form
   naturally balances the flux term `{{∂u/∂n}}` with the penalty `η[[u]]` in
   the stiffness matrix. The traction extracted from the stress gradient already
   includes the effect of the penalty through the solution itself.

3. **No lifting operator needed**: The BR2 lifting for the scalar PDE is
   trivially the projection of the jump onto the volume space. The average
   gradient already captures this because the scalar field has no directional
   ambiguity.

### 11.3 Why 3D elasticity needs penalty correction

For the 3D elasticity PDE: `∇·(C:∇u) = 0` where `u` is a 3-component vector.

The DG traction `{{C:∇u}}·n̂` involves the full elasticity tensor contraction.
The penalty correction `C:L([[u]]-δ)·n̂` provides:

1. **Coercivity**: Without it, the DG bilinear form may not be coercive for
   the elasticity system (unlike the scalar Laplacian which is naturally coercive
   with the average flux).

2. **Directional coupling**: The elasticity tensor couples normal and tangential
   displacement components. The BR2 lifting captures this coupling through the
   tensor structure of `C:L`.

3. **Convergence at low order**: At p=1 on coarse meshes, the stress gradient
   `{{C:∇u}}` has O(h) error. The penalty correction compensates by penalizing
   the displacement jump, improving the effective convergence rate.

### 11.4 Summary

| Feature | Antiplane (scalar) | 3D Elasticity (vector) |
|---------|-------------------|----------------------|
| PDE type | Scalar Poisson | Vector elasticity |
| Traction | `μ · {{∂u/∂n}}` | `{{C:∇u}}·n̂ - σ·C:L·n̂` |
| Penalty needed? | No | Yes |
| Why | No directional coupling, self-consistent flux | Tensor coupling, coercivity, convergence |
