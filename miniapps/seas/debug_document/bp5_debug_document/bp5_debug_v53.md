# BP5 Debug v53: p=1 Fault Discretization Mismatch with Tandem

**Date**: 2026-03-28
**Status**: ROOT CAUSE IDENTIFIED, NO CODE CHANGES YET
**Previous**: v52 (traction coherence, BLR residual amplification, dip contamination)
**Branch**: `feature/elasticity`

---

## 1. Motivation

After the recent reverts, the previously suspected `dt_init` / v53 startup-path issue is no
longer the leading explanation for the BP5 `p=1`, near-fault 1000 m mismatch.

Two earlier hypotheses were also checked and are no longer the primary suspects:

1. **The IP `3x` penalty factor is wrong at `p=1`**  
   This is **false**. A term-by-term code comparison shows the `3x` factor is the geometric
   normalization required for MFEM's tetrahedral `CalcOrtho/Weight()` convention to match
   Tandem's physical `A/V` penalty coefficient.

2. **MFEM boundary loading misses Tandem's center-line `Vp*t` branch**  
   This is also **false**. MFEM has separate non-fault `Y=0` interior/shared Dirichlet logic
   that reproduces the Tandem center-line jump construction.

The remaining question was: if penalty and boundary loading are not the source, what still
makes MFEM `p=1` disagree with Tandem `p=1`?

The answer is now clear:

> **Current MFEM `p=1` does not solve the same discrete fault problem as Tandem `p=1`.**

MFEM collapses each triangular fault face to **one constant fault DOF** at `p=1`, while
Tandem uses **three nodal fault DOFs per triangular face** at `p=1`.

That changes the sampled BP5 problem data (`a`, `Dc/L`, `V_init`, `tau_pre`) near the sharp
nucleation-zone boundary and makes the comparison non-like-for-like.

---

## 2. Main Finding

### 2.1 Tandem `p=1` Uses Three Fault Nodes per Triangle

Tandem's rate-and-state fault space is always built from
`NodalRefElement<DomainDimension - 1>(PolynomialDegree)`.

Source:

```cpp
auto RateAndStateBase::Space() -> NodalRefElement<DomainDimension - 1u> {
    return NodalRefElement<DomainDimension - 1u>(
        PolynomialDegree, WarpAndBlendFactory<DomainDimension - 1u>(), ALIGNMENT);
}
```

For 3D elasticity, the fault is 2D, so at `PolynomialDegree = 1` a triangular fault face has
**three nodal basis functions / three nodal samples**.

Tandem then assigns friction/problem data at **every fault node**:

```cpp
for (std::size_t index = 0; index < num_nodes; ++index) {
    auto params = pfun(fault_.storage()[index].template get<Coords>());
    law_.set_params(index, params);
}
```

So Tandem `p=1` is already a **multi-DOF fault discretization**.

### 2.2 MFEM `p=1` Uses One Face-Averaged Fault DOF

MFEM explicitly hard-codes:

- `p=1` -> `face_order = 0`
- `face_order = 0` -> `nbf = 1`
- `nbf = 1` means one constant fault DOF per face

Source:

```cpp
// p=1: face_order=0, nbf=1 -> backward compatible (face average)
int face_fe_order = (method_ == DGMethod::IP && order_ >= 2) ? order_ : 0;
face_quad_ = std::make_unique<FaceQuadrature>(face_fe_order, ...);
nbf_per_face_ = face_quad_->NumBasisFunctions();
```

and

```cpp
///   - p=1 (vol_order=1): face_order=0, nbf=1 -> proven stable, backward compatible
///   - p≥2 (vol_order≥2): face_order=vol_order, nbf=(p+1)(p+2)/2 -> matches Tandem
```

So current MFEM only matches Tandem's fault discretization at `p>=2`, not at `p=1`.

### 2.3 Consequence for BP5 Initialization

MFEM evaluates BP5 parameter fields at its fault DOFs. At `p=1`, that means one sample per
face, effectively at the face centroid:

```cpp
// At nbf=1 (p=1): single centroid point -> same as old code
const IntegrationRule &nir = face_quad_->GetNodalRule();
...
fault_x2_(i * nbf + kk) = coords(0);
fault_x3_(i * nbf + kk) = -coords(2);
```

Then:

```cpp
for (int i = 0; i < num_fault_dofs_; i++)
{
   real_t x2 = coords_x2_(i);
   real_t x3 = coords_x3_(i);

   a_values_(i) = bp5_params_.a_of_x2_x3(x2, x3);
   dc_values_(i) = bp5_params_.Dc_of_x2_x3(x2, x3);
   bp5_params_.tau0_vec(x2, x3, tau);
   bp5_params_.V_init_vec(x2, x3, Vi);
}
```

Therefore, on each `p=1` fault triangle:

- Tandem samples BP5 data at **3 face nodes**
- MFEM samples BP5 data at **1 face centroid**

This is not a small implementation detail. It changes the discrete nucleation patch itself.

---

## 3. Why This Matters for the 1000 m BP5 `p=1` Case

BP5 has a **sharp rectangular nucleation zone**:

- depth: `hs + ht <= x3 <= hs + ht + H`
- strike: `-l/2 <= x2 <= -l/2 + w`

In Tandem, those tests are applied to each fault node. In current MFEM `p=1`, they are
applied only at the face centroid.

This means a face that straddles the nucleation boundary can be treated as:

- **partially inside** in Tandem (`3` nodal samples, potentially mixed values), but
- **fully inside or fully outside** in MFEM (`1` centroid sample only)

That directly changes:

1. `Dc = L_nuc` vs `L0`
2. `V_init = V_nuc` vs `Vp`
3. `tau_pre`
4. local friction state through `psi` equilibrium initialization

The result is a different discrete seed geometry and different initial traction/slip-rate
balance near the nucleation front.

For the 1000 m mesh, this is a strong candidate for why Tandem `p=1` nucleates while current
MFEM `p=1` does not.

---

## 4. Secondary Configuration Difference: `bp5_outside`

There is also a smaller, but real, BP5 configuration mismatch:

- Tandem's stock BP5 case uses `scenario = "bp5_outside"`
- `bp5_outside = BP5.new({eps=1e-3})`
- MFEM currently uses exact inclusion logic with no `eps` offset

This affects points lying exactly on the nucleation-zone boundary.

However, this is **not** the primary issue found here. The larger inconsistency is the fault
discretization mismatch itself:

- Tandem `p=1`: 3 nodal samples / face
- MFEM `p=1`: 1 centroid sample / face

---

## 5. Conclusion

The current BP5 `p=1` mismatch is best explained by a **discrete model mismatch**, not by a
penalty-coefficient mismatch.

### 5.1 What Is Now Disproved

1. **`3x` penalty is wrong for `p=1`**  
   Disproved. The `3x` factor is required to match Tandem's IP penalty coefficient.

2. **MFEM is missing Tandem's `Y=0` center-line loading branch**  
   Disproved. That logic exists in the interior/shared Dirichlet assembly.

3. **The reverted `bp5_v53_p1_exact_7618935.out` startup path is still the active blocker**  
   No longer supported after the reverts.

### 5.2 What Is Now the Leading Root Cause

> **MFEM `p=1` fault discretization is not Tandem `p=1` fault discretization.**

Current MFEM uses one face-averaged fault DOF per triangle at `p=1`.  
Tandem uses three nodal fault DOFs per triangle at `p=1`.

As a result, MFEM and Tandem are not sampling the same BP5 fields on the same discrete fault
space, so the observed nucleation mismatch is unsurprising.

---

## 6. Proposed Fix

### 6.1 Goal

Achieve **exact Tandem-style fault sampling at `p=1`**:

- triangular fault face
- polynomial degree `p=1`
- **3 nodal fault DOFs per face**
- BP5 parameters evaluated at those 3 face nodes
- slip interpolated from those 3 DOFs to quadrature points
- traction projected back onto those 3 DOFs

In short:

> **MFEM `p=1` should use the same nodal fault space as Tandem `p=1`, not a face average.**

### 6.2 Concrete Implementation Target

For the IP elasticity path:

1. Remove the special-case reduction
   - current: `p=1 -> face_order=0`
   - target: `p=1 -> face_order=1`

2. Construct `FaceQuadrature` with `face_order = order_` for `p=1` as well
   - for triangles, this gives `nbf = (1+1)(1+2)/2 = 3`

3. Keep per-DOF BP5 parameter evaluation exactly as already written
   - once `nbf_per_face_ = 3`, the existing `GetFaultCoords2D()` and
     `ComputeBP5Params()` machinery will naturally sample all three nodes

4. Keep the fault traction/slip machinery in nodal form
   - interpolate nodal slip to quad points
   - compute traction at quad points
   - Galerkin project traction back to nodal DOFs

This is already the design of the multi-DOF fault path. The required fix is mainly to make
that path active for `p=1`, not only for `p>=2`.

### 6.3 Exact Matching Target

To match Tandem as closely as possible, the intended final state is:

- `p=1` fault faces use **3 nodal samples per face**
- node locations follow the same reference-face nodal set as Tandem's
  `NodalRefElement<2>(1)`
- BP5 parameter evaluation is nodal, not centroidal
- `bp5_outside`-style inclusion tolerance should also be considered after the nodal fix

The nodal fault-space change is the first-priority correction.

---

## 7. Recommended Validation Plan After the Fix

After implementing the `p=1 -> 3 DOF / face` change:

1. Print `nbf_per_face_` at startup for `p=1`
   - expected: `3`

2. Dump a few fault-face coordinates at the nucleation boundary
   - confirm three distinct nodal points per face, not one centroid

3. Compare nucleation-zone classification counts
   - number of DOFs with `Dc = 0.13`
   - number of DOFs with `V_init = 0.01`
   - compare with Tandem

4. Re-run BP5 near-fault 1000 m, `p=1`
   - compare station nucleation timing against Tandem

5. Only after this, revisit secondary differences
   - `bp5_outside` vs exact inclusion
   - any remaining strike/dip traction differences

---

## 8. Bottom Line

The current BP5 `p=1` comparison is not apples-to-apples.

The critical mismatch is:

- **Tandem `p=1`: 3 fault nodes per triangular face**
- **MFEM `p=1`: 1 face-average DOF per triangular face**

The proposed fix is to make MFEM `p=1` use the same nodal fault discretization as Tandem,
with **three node samples per element face exactly**.
