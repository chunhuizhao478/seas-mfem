# BP5 Debug v32/v33: Recurrence Interval Diagnosis + General Polynomial Order Fix

**Date**: 2026-03-15
**Status**: Implementation in progress
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
