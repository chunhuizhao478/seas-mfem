# BP5 Debug v11: BR2 Traction Correction — Normal Inconsistency Audit

## Problem

The BR2 penalty correction in `ComputeTraction` uses two different normals:
- **face_int** (lifting construction): `nor(s)` — unnormalized normal from `CalcOrtho` (= |J_F| * n̂)
- **TestNormal** (traction evaluation): `basis.normal[s]` — unit normal n̂

This inconsistency exists in both interior faces (lines ~2000,2004) and shared faces (lines ~2232,2234).

## Document Conventions

| Symbol | Code Variable | Description |
|--------|--------------|-------------|
| n̂ | `basis.normal` | Unit outward normal at fault face |
| nor | `nor` from `CalcOrtho` | Unnormalized face normal = \|J_F\| · n̂ |
| [[u]] | `u_jump` | Displacement jump u₁ - u₂ |
| δ | `delta_u` | Prescribed slip embedded in 3D |
| M⁻¹ | `Minv1`, `Minv2` | Element mass matrix inverse |
| φ_m | `shape1(m)`, `shape2(m)` | Shape function values |
| η_F | `br2_penalty` | BR2 penalty parameter (2·dim for hex, dim+1 for tet) |
| C_{iu,sq} | `TestNormal` | Elasticity tensor contracted with normal |

## Term-by-Term Comparison

### Bilinear Form Integrator (`dg_elasticity_br2_integrator.hpp`)

In the stiffness matrix assembly, the BR2 integrator uses **unnormalized** normals consistently:

```
// Lifting (line 261):
face_int[...] = E_src * nor_all[q*dim+s] * wq    // unnormalized

// TestNormal (line 291):
n_q(d) = nor_all[q*dim+d]                         // unnormalized
```

Comment at lines 225-229 explicitly states:
> "the lifting uses the unnormalized normal (from CalcOrtho), and L_q contracts
> with unnormalized normal again. This gives the correct face-area scaling.
> The test_normal also uses the unnormalized normal."

This is correct for the **bilinear form** because the quadrature weight `w_q` (reference-face weight) combined with `|J_F|²` from the two appearances of `nor` gives the correct area-weighted integration.

### ComputeTraction (point evaluation)

In `ComputeTraction`, we perform a **point evaluation** of the BR2 correction at the face centroid — no integration. The correction should be:

  T_correction_i = η_F · 0.5 · Σ_{u,s} C_{iu,s}(n̂) · [Σ_m φ_m(x_c) · (M⁻¹ · face_int)_{us,m}]

where `face_int_{us,m} = φ_m(x_c) · jump_u · n̂_s` (using unit normal).

**Current code (BUGGY)**:
```cpp
// face_int: uses nor(s) = |J_F| · n̂_s  [UNNORMALIZED]
face_int1(u * dim + s, m) = shape1(m) * jump[u] * nor(s);

// TestNormal: uses basis.normal[s] = n̂_s  [UNIT]
real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * basis.normal[i]
   + mu_val_ * ((i == u ? 1.0 : 0.0) * basis.normal[s]
                + (i == s ? 1.0 : 0.0) * basis.normal[u]);
```

This mixes `|J_F| · n̂` in the lifting with `n̂` in the traction, producing:
  T ∝ |J_F| · (C : (jump ⊗ n̂) · n̂)

The spurious `|J_F|` factor scales the correction by the face area, which is mesh-dependent and wrong for a point evaluation.

### Expected Impact

- On uniform meshes with unit-sized faces, `|J_F| ≈ 1`, so the bug has minimal effect
- On graded meshes (like BP5 with 1km → 40km elements), `|J_F|` varies widely, causing:
  - Oversized corrections on large faces → traction noise
  - Undersized corrections on small faces → insufficient stabilization
- This explains some of the DG traction noise observed in v10 debug

## Fix

**Change `nor(s)` → `basis.normal[s]`** in the `face_int` computation in both interior and shared face loops of `ComputeTraction`.

### Interior faces (~lines 2000, 2004):
```cpp
// BEFORE:
face_int1(u * dim + s, m) = shape1(m) * jump[u] * nor(s);
face_int2(u * dim + s, m) = shape2(m) * jump[u] * nor(s);

// AFTER:
face_int1(u * dim + s, m) = shape1(m) * jump[u] * basis.normal[s];
face_int2(u * dim + s, m) = shape2(m) * jump[u] * basis.normal[s];
```

### Shared faces (~lines 2232, 2234):
```cpp
// BEFORE:
face_int1(u * dim + s, m) = shape1(m) * jump[u] * nor(s);
face_int2(u * dim + s, m) = shape2(m) * jump[u] * nor(s);

// AFTER:
face_int1(u * dim + s, m) = shape1(m) * jump[u] * basis.normal[s];
face_int2(u * dim + s, m) = shape2(m) * jump[u] * basis.normal[s];
```

### Rationale

- `ComputeTraction` is a **point evaluation**, not a face integral
- The bilinear form uses unnormalized normals because `|J_F|²` provides the correct area element for integration
- For point evaluation: T = C : R_h(jump) · n̂, all normals should be **unit** n̂
- The `TestNormal` (line 2025-2027) already uses `basis.normal` (unit) — the fix makes `face_int` consistent

## Tests Added

1. **BR2 traction correction consistency**: Verify that the traction correction is independent of face area scaling (same correction on a 2x-scaled mesh)
2. **Patch test**: Constant strain field → exact stress recovery with BR2 traction
3. **Slip sign convention**: Positive strike slip → negative traction (stress drop)
