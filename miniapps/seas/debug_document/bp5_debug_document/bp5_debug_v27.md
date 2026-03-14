# BP5 Debug v27: Matching Tandem's Traction Exactly

**Date**: 2026-03-14
**Status**: Detailed comparison, precise differences identified
**Previous**: v26 (IP blows up worse, need deeper investigation)

---

## 1. Precise Traction Formula Comparison

### Tandem's traction on a fault face (skeleton, two-sided)

From `elasticity.py` line 242-244, `Elasticity.cpp` line 972:

```
T_p(q) = 0.5 * (σ₁·n̂ + σ₂·n̂)_p  +  c0 * (u₁_p - u₂_p - slip_p)
```

where:
- `c0 = -penalty(fctNo)`
- For BR2: `penalty = NumFacets = 4` (dimensionless)
- For IP: `penalty = (p(0)+p(1))/4` where `p(side) = (D+1)*c_N_1*(A/V)*(c₁²/c₀)`

**The penalty correction is SCALAR × VECTOR**: a single scalar `c0` multiplied by
the displacement jump vector. Each component `correction_p = c0 * jump_p`.
No elasticity tensor in the correction.

### MFEM's traction on a fault face (IP path)

From `elasticity_operator.hpp` lines 2141-2169:

```
correction_i = ip_coeff * Σ_{u,s} [λ δ_{us} n̂_i + μ(δ_{iu} n̂_s + δ_{is} n̂_u)] * n̂_s * jump_u
```

where `ip_coeff = κ * |n|² * (1/(2 detJ₁) + 1/(2 detJ₂))`

**The penalty correction is TENSOR-COUPLED**: the elasticity tensor C contracted
with the normal vector, applied to the jump. This is NOT a scalar × vector —
it's `ip_coeff * (C : n̂⊗n̂) · jump`.

### MFEM's traction on a fault face (BR2 path, current H26 fix)

```
correction_c = penalty_val * jump_c    (penalty_val = 4 for tet)
```

This matches Tandem's BR2 formula: scalar × vector.

---

## 2. The Structural Difference

| Aspect | Tandem | MFEM (IP) | MFEM (BR2 + H26) |
|--------|--------|-----------|-------------------|
| Average stress | 0.5*(σ₁+σ₂)·n̂ | 0.5*(σ₁+σ₂)·n̂ | 0.5*(σ₁+σ₂)·n̂ |
| Correction type | **scalar × jump** | **tensor × jump** | **scalar × jump** |
| Correction value (BR2) | -4 * jump | N/A | -4 * jump |
| Correction value (IP) | -penalty * jump | -ip_coeff * (C:n̂⊗n̂)·jump | N/A |
| Penalty units (BR2) | dimensionless | N/A | dimensionless |
| Penalty units (IP) | Pa/m (scalar) | 1/m² × Pa (tensor) | N/A |

The key difference for the IP method: **Tandem uses scalar penalty, MFEM uses
tensor-coupled penalty**. The tensor coupling means different jump directions
get different effective penalties, which can amplify certain components at
fault edge faces.

---

## 3. What MFEM's IP Traction Actually Computes

Expanding the MFEM IP correction for a fault with normal n̂ = (1,0,0):

```
correction_i = ip_coeff * Σ_{u,s} C_{ius,j} n̂_j n̂_s * jump_u
```

With n̂ = (1,0,0):
```
correction_1 = ip_coeff * [(λ+2μ) * jump_1]                    (normal component)
correction_2 = ip_coeff * [μ * jump_2]                          (strike component)
correction_3 = ip_coeff * [μ * jump_3]                          (dip component)
```

The normal-direction penalty is (λ+2μ)/μ = 3× stronger than the tangential
penalty. This anisotropy doesn't exist in Tandem's scalar penalty.

For BP5 with λ=μ=32 GPa: the normal penalty is 96 GPa * ip_coeff, the
tangential is 32 GPa * ip_coeff. The 3:1 ratio means normal jumps are penalized
3× more than tangential jumps.

### Why this matters at fault edges

At fault edge faces where element geometry transitions, the displacement jump
may have components in all three directions. The tensor penalty amplifies the
normal component 3× more than tangential. If the normal component is large
(from element mismatch), the correction overshoots, creating an even larger
jump in the next time step → positive feedback → blowup.

Tandem's scalar penalty treats all directions equally, avoiding this
directional amplification.

---

## 4. Why Tandem's Approach Works

### For BR2:
- Penalty = 4 (dimensionless) → correction ≈ 0 for any jump
- Traction = average stress only
- Stable because the BR2 bilinear form already ensures [[u]] ≈ δ
- No directional amplification

### For IP:
- Penalty = material-dependent scalar (Pa/m)
- Correction = penalty * jump (isotropic, all directions treated equally)
- Stable because the penalty magnitude is properly scaled:
  - Strong enough to prevent large jumps (Pa/m units)
  - No directional bias (scalar × vector)
  - Consistent with the IP bilinear form penalty

---

## 5. The Fix: Match Tandem's Traction Formula Exactly

### For BR2 method:

Keep the H26 fix: `correction = 4 * jump` (scalar × vector).
This matches Tandem exactly. The instability we saw was NOT from this formula —
it was from the combination of BR2 bilinear form + this traction on certain
mesh faces. This needs further investigation (mesh quality or time stepping).

### For IP method:

**Replace the tensor-coupled penalty with Tandem's scalar penalty**:

```cpp
// CURRENT (MFEM, tensor-coupled):
correction[i] = ip_coeff * Σ_{u,s} C_{ius,j} n̂_j n̂_s * jump[u]

// FIX (match Tandem, scalar):
real_t penalty = compute_tandem_ip_penalty(FTr);  // material-dependent, Pa/m
for (int c = 0; c < dim; c++)
    correction[c] = penalty * jump[c];
```

where the penalty matches Tandem's formula:
```cpp
real_t compute_tandem_ip_penalty(FaceElementTransformations *FTr) const
{
    // Tandem: penalty = (p(0) + p(1)) / 4
    // p(side) = (D+1) * c_N_1 * (area/volume) * (c1²/c0)
    int dim = 3;
    real_t c0 = 2.0 * mu_val_;                                // min eigenvalue of C
    real_t c1 = dim * lambda_val_ + 2.0 * mu_val_;            // max eigenvalue of C
    real_t c_N_1 = 1.0;  // trace constant for p=0 (order 1 DG, PolynomialDegree-1=0)

    real_t area = FTr->Face->Weight();  // face area (det of face Jacobian)
    // Note: need to integrate properly for curved elements

    real_t vol1 = FTr->Elem1->Weight();
    real_t vol2 = FTr->Elem2->Weight();

    real_t p0 = (dim + 1) * c_N_1 * (area / vol1) * (c1 * c1 / c0);
    real_t p1 = (dim + 1) * c_N_1 * (area / vol2) * (c1 * c1 / c0);

    return (p0 + p1) / 4.0;
}
```

### The same formula should be used in the IP bilinear form

The MFEM `DGElasticityIntegrator` uses `kappa * |n|²/detJ` as the penalty
in the bilinear form. For full consistency with Tandem, this should also use
the material-dependent scalar penalty. However, this is a larger change —
it affects the stiffness matrix assembly, not just the traction.

For now, we can test the traction fix alone (scalar penalty in traction)
while keeping the existing bilinear form. If the traction is stable, we can
investigate the bilinear form penalty later.

---

## 6. Open Questions

### Why does BR2 + scalar traction (H26) blow up at ~2 years?

With BR2 bilinear form and `correction = 4 * jump` in traction, specific DOFs
still blow up at t≈2yr on both old and new mesh. This can't be explained by
the traction formula alone (it matches Tandem exactly for BR2).

Possible causes:
1. **The bilinear form penalty (BR2, η_e=4)** may not be strong enough on
   certain faces, leading to large [[u]] that the negligible traction penalty
   can't control
2. **The time stepper** may handle the feedback loop differently from Tandem's
   PETSc TS (Tandem uses PETSc's TSRK 5dp with line search)
3. **The slip contribution to the RHS** in the bilinear form may differ
   between Tandem and MFEM, affecting how well [[u]] tracks δ

### Does Tandem ACTUALLY blow up with the same mesh?

We assumed Tandem is stable, but we haven't verified this. Tandem might also
have specific DOFs with large traction, but PETSc's time stepper might handle
them (step rejection, line search) without diverging.

---

## 7. Test Plan

### Test A: IP method + Tandem scalar penalty in traction + new mesh

Replace MFEM's tensor-coupled IP traction penalty with Tandem's scalar penalty.
Keep the existing IP bilinear form (MFEM's kappa-based). Run uniform test.

### Test B: BR2 method + original BR2 traction + new mesh

Revert ALL traction changes (restore original BR2 lifting in traction).
Use the new uniform mesh. This tests whether the better mesh alone reduces
the VS lockup enough to be acceptable.

### Test C: Add diagnostic output for ALL fault face tractions

Before fixing anything, add output that prints the traction magnitude at every
fault face at step 1 and when blowup occurs. This reveals which faces have
anomalous traction from the very start vs which develop it over time.

---

## 8. Implementation: H28 — Tandem Scalar IP Traction Penalty

### Code changes

File: `domain/elasticity_operator.hpp`, `ComputeTraction()`, both interior and
shared face IP paths.

**Before (tensor-coupled)**:
```cpp
real_t ip_coeff = kappa * nor_sq * (1/(2*detJ1) + 1/(2*detJ2));
// correction[i] = ip_coeff * Σ_{u,s} C_{ius,j} n̂_j n̂_s * jump[u]
// → 3:1 directional ratio (normal vs tangential)
```

**After (Tandem scalar)**:
```cpp
real_t c0_mat = 2.0 * mu_val_;                       // min eigenvalue of C
real_t c1_mat = dim * lambda_val_ + 2.0 * mu_val_;   // max eigenvalue of C
real_t c_N_1 = 1.0;                                  // trace constant for p=0
real_t face_area = nor.Norml2();
real_t vol1 = FTr->Elem1->Weight();
real_t vol2 = FTr->Elem2->Weight();
real_t p0 = (dim+1) * c_N_1 * (face_area/vol1) * (c1_mat*c1_mat/c0_mat);
real_t p1 = (dim+1) * c_N_1 * (face_area/vol2) * (c1_mat*c1_mat/c0_mat);
real_t penalty_ip = (p0 + p1) / 4.0;
// correction[c] = penalty_ip * jump[c]
// → isotropic, all directions penalized equally
```

### Why this should work

1. **Isotropic**: No directional amplification at fault edge faces
2. **Properly scaled**: Has units Pa/m (face_area/volume ~ 1/h, c1²/c0 ~ Pa)
3. **Material-dependent**: Scales with elastic moduli
4. **Matches Tandem exactly**: Same formula as `prepare_penalty` in Elasticity.cpp:280-290

### Note on bilinear form vs traction penalty

The IP **bilinear form** in MFEM uses `jmatcoef = kappa * |n|² * wLM` where
`wLM` includes `(λ+2μ)/detJ` (line 4161 of bilininteg.cpp). This is
tensor-coupled in the assembly but effectively computes a similar scalar
magnitude.

Tandem's IP **bilinear form** also uses tensor-coupled L_q (lines 103-109 of
elasticity.py). But the **traction** uses scalar penalty in both codes (now).

This means the bilinear form and traction use different penalty structures —
tensor vs scalar. This is intentional in Tandem. The bilinear form needs the
full tensor for coercivity. The traction is post-processing where isotropic
scalar penalty is more stable.

---

## 9. Test Plan

### Test 1: IP + scalar traction + uniform fault (no earthquake)

```
sbatch: bp5_v27_test1_ip_scalar.sbatch
Flags:  --dg-method IP --delta-tau-factor 0 --V-nuc 1e-9
Mesh:   bp5_1000m.msh (uniform res_f=1)
```

PASS: V/Vp in [0.5, 2.0] at z=22km for 300 years, no traction blowup.

### Test 2: IP + scalar traction + Tandem init (if Test 1 passes)

```
Flags:  --dg-method IP --psi-init-mode tandem --V-nuc 0.01
```

PASS: First event ~150yr, VS maintains V≈Vp, cycling.

### Test 3: IP + scalar traction + SCEC init (if Test 2 passes)

```
Flags:  --dg-method IP (defaults)
```

PASS: Immediate earthquake, VS recovers, cycling.

---

## 10. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H25 | Previous fixes | Done |
| H26 | Traction: scalar penalty=4 for BR2 (matches Tandem BR2) | Done, blows up at ~2yr |
| H27 | Mesh: uniform fault resolution (Tandem-style) | Done |
| **H28** | **Traction: Tandem scalar IP penalty (isotropic, Pa/m)** | **Done, testing** |
