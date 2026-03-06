# Discontinuous Galerkin Formulation for 3D Linear Elasticity

This document provides the complete mathematical formulation for the 3D DG elasticity solver used in the BP5 SEAS benchmark, with clear mapping to implementation files. It parallels the antiplane theory document (`document/archives/dg_antiplane_theory.md`).

---

## Table of Contents

**Part I: Mathematical Formulation**
1. [Strong Form (Governing Equation)](#1-strong-form-governing-equation)
2. [Computational Domain and Boundary Conditions](#2-computational-domain-and-boundary-conditions)
3. [DG Notation](#3-dg-notation)
4. [BR2 Method: Complete Formulation](#4-br2-method-complete-formulation)
5. [IP Method: Complete Formulation](#5-ip-method-complete-formulation)
6. [Traction Computation](#6-traction-computation)
7. [Fault Coordinate Transformation](#7-fault-coordinate-transformation)

**Part II: Implementation Mapping**
8. [File-to-Formula Mapping](#8-file-to-formula-mapping)

**Appendix: Tandem Reference**
- [A. Tandem Implementation Details](#appendix-a-tandem-implementation-details)
- [B. Tandem Source File Reference](#appendix-b-tandem-source-file-reference)

---

# Part I: Mathematical Formulation

---

## 1. Strong Form (Governing Equation)

### 1.1 Physical Problem: 3D Linear Elasticity

**Displacement field** (full 3D vector):

```
u = (u₁(x), u₂(x), u₃(x)),    x = (x₁, x₂, x₃) ∈ Ω ⊂ ℝ³
```

where (x₁, x₂, x₃) correspond to (fault-normal, along-strike, depth) in the SCEC convention.

**Strain tensor** (symmetric gradient):

```
ε_ij(u) = (1/2)(∂u_i/∂x_j + ∂u_j/∂x_i)
```

**Stress tensor** (isotropic linear elasticity, Hooke's law):

```
σ_ij(u) = λ tr(ε) δ_ij + 2μ ε_ij
         = λ (∂u_k/∂x_k) δ_ij + μ (∂u_i/∂x_j + ∂u_j/∂x_i)
```

where:
- `λ` = first Lame parameter [Pa]
- `μ` = shear modulus [Pa]
- For BP5: `ρ = 2670 kg/m³`, `cs = 3464 m/s`, `ν = 0.25`
  - `μ = ρ cs² = 32.04 GPa`
  - `λ = 2νμ/(1-2ν) = μ` (since ν = 0.25)

### 1.2 Strong Form PDE

**Equilibrium** (quasi-static, no body forces):

```
┌───────────────────────────────────────────────────┐
│  -∇·σ(u) = 0      in Ω                            │
│                                                    │
│  i.e., -∂σ_ij/∂x_j = 0    for i = 1,2,3          │
└───────────────────────────────────────────────────┘
```

### 1.3 Boundary Conditions

**Dirichlet BC** (prescribed displacement):
```
u = g_D      on Γ_D
```

**Neumann BC** (zero traction, natural):
```
σ(u) · n = 0      on Γ_N
```

**Interior fault condition** (prescribed slip jump, no-opening):
```
[[u]] = u⁺ - u⁻ = δ      on Γ_F
```

where the slip `δ` is constrained to lie in the fault tangent plane:
```
δ = δ_dip · t₁ + δ_strike · t₂      (no opening: δ · n = 0)
```

### 1.4 Comparison with Antiplane (Scalar) Case

| Aspect | Antiplane (BP1/BP2) | Full Elasticity (BP5) |
|--------|--------------------|-----------------------|
| PDE | -∇·(K∇u) = 0 (Laplace) | -∇·σ(u) = 0 (Navier) |
| Unknown | Scalar u(x,z) | Vector u = (u₁, u₂, u₃) |
| Flux | K∇u | σ(u) = λ tr(ε)I + 2με |
| Traction | τ = μ ∂u/∂x | T = σ·n (vector) |
| Slip | Scalar δ | Vector δ = (δ_dip, δ_strike) |
| DG space | Scalar DG | Vector DG (vdim=3) |
| Spatial dimension | 2D | 3D |

The scalar Laplace case is recovered when the displacement is purely out-of-plane, `u = (0, 0, u₃(x₁, x₃))`, reducing the stress to `σ₁₃ = μ ∂u₃/∂x₁`, `σ₃₃ = μ ∂u₃/∂x₃`.

---

## 2. Computational Domain and Boundary Conditions

### 2.1 Domain Geometry (BP5)

```
         x₃ = 0 (free surface)
         Neumann: σ·n = 0
    ┌─────────────────────────────────────────────────────────────┐
    │                              │                              │
    │                              │                              │
    │      Left half               │       Right half             │
    │      (x₁ < 0)               │       (x₁ > 0)              │
    │                              │                              │
x₁=-Lx                             │x₁=0                      x₁=+Lx
Neumann                             │                          Neumann
σ·n = 0                          FAULT                        σ·n = 0
    │                    [[u]] = δ  │                              │
    │                              │                              │
    │                              │                              │
    └─────────────────────────────────────────────────────────────┘
         x₃ = Lz (bottom)
         Neumann: σ·n = 0

    x₂ = -Ly  ────────────────────────────────────────  x₂ = +Ly
    Dirichlet: u₂ = -Vp·t/2                            Dirichlet: u₂ = +Vp·t/2
               u₁ = u₃ = 0                                        u₁ = u₃ = 0
```

The domain is `Ω = [-Lx, Lx] × [-Ly, Ly] × [0, Lz]` with the fault at `x₁ = 0`.

### 2.2 Boundary Classification (BP5)

| Boundary | Location | Attr | BC Type | Value | Implementation |
|----------|----------|------|---------|-------|----------------|
| Fault-normal left | x₁ = -Lx | 1 | Neumann | σ·n = 0 | Natural (no integrator) |
| Fault-normal right | x₁ = +Lx | 2 | Neumann | σ·n = 0 | Natural (no integrator) |
| Strike+ wall | x₂ = +Ly | 3 | Dirichlet | u₂ = +Vp·t/2 | `AssembleDirichletLoading()` |
| Strike- wall | x₂ = -Ly | 4 | Dirichlet | u₂ = -Vp·t/2 | `AssembleDirichletLoading()` |
| Free surface | x₃ = 0 | 5 | Neumann | σ·n = 0 | Natural (no integrator) |
| Bottom | x₃ = Lz | 6 | Neumann | σ·n = 0 | Natural (no integrator) |
| Fault | x₁ = 0 | — | Jump | [[u]] = δ | `AssembleSlipContribution*()` |

**Key difference from antiplane**: Plate loading is on the ±x₂ walls (not the bottom), and only the x₂ component of displacement is prescribed.

### 2.3 Coordinate Convention (SCEC BP5)

```
x₁ = fault-normal direction     (mesh x-axis)
x₂ = along-strike direction     (mesh y-axis)
x₃ = depth (positive downward)  (mesh z-axis)
```

Fault at x₁ = 0, extending from x₃ = 0 to x₃ = Wf in depth and |x₂| ≤ lf/2 in strike.

---

## 3. DG Notation

### 3.1 Mesh Decomposition

- **Elements**: Ω = ∪_K K (hexahedral or tetrahedral)
- **Interior faces**: F_I = interior faces shared by two elements, excluding fault faces
- **Fault faces**: F_F = interior faces with prescribed slip (F_F ⊂ interior faces, F_F ∩ F_I = ∅)
- **Dirichlet faces**: F_D = faces on Γ_D (±x₂ walls)
- **Neumann faces**: F_N = faces on Γ_N (all other boundaries)

**Note**: F_I and F_F are disjoint subsets of all interior faces. Following Arnold et al. (2002) eq. 3.24, the bilinear form's interior face terms sum over **F_I ∪ F_F** (all interior faces). The fault RHS shifts the penalty enforcement from [[u]] = 0 to [[u]] = δ.

### 3.2 Jump and Average Operators (Vector-Valued)

For an interior face `e` shared by elements K⁺ and K⁻, with outward normal `n` pointing from K⁻ to K⁺:

**Jump** (vector):
```
[[v]] = v⁺ - v⁻
```

**Average** (vector):
```
{{v}} = (v⁺ + v⁻) / 2
```

**Average of stress flux** (tensor contracted with normal):
```
{{σ(v)·n}} = (σ(v⁺)·n + σ(v⁻)·n) / 2
```

For boundary faces: `[[v]] = v`, `{{v}} = v`.

### 3.3 DG FE Space

**Vector DG space** (3 components):
```
V_h = {v ∈ [L²(Ω)]³ : v|_K ∈ [P^p(K)]³  ∀K ∈ T_h}
```

where `P^p(K)` is the polynomial space of degree ≤ p on element K.

**Implementation**: `DG_FECollection(order, 3, BasisType::GaussLobatto)` with `FiniteElementSpace(&mesh, fec, 3, Ordering::byNODES)`.

### 3.4 Traction Operator

The traction operator `σ(v)·n` applied to a single basis function `φ_k` with component `a` gives:

```
[σ(φ_k e_a) · n]_j = λ (∂φ_k/∂x_a) n_j + μ (δ_ja ∇φ_k · n + (∂φ_k/∂x_j) n_a)
```

where `e_a` is the unit vector in direction `a`. This is the vector elasticity generalization of the scalar flux `K ∇φ_k · n`.

---

## 4. BR2 Method: Complete Formulation

The **Bassi-Rebay 2 (BR2)** method is the default DG method, matching Tandem. It uses lifting operators with elasticity tensor coupling.

### 4.1 Lifting Operator Definition (Vector Elasticity)

For an interior face `e` shared by elements K⁻ and K⁺, the lifting operator `r_e` maps the jump `[[v]]` to a tensor-valued function:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  ∫_{K⁻ ∪ K⁺} r_e([[v]]) : w dx = ∫_e {{w · n}} · [[v]] ds    ∀ tensor w   │
└──────────────────────────────────────────────────────────────────────────────┘
```

**Local computation** (per element, per component):

For the scalar case, the lifting is `r_e|_{K±} = M_K⁻¹ ∫_e φ ⊗ n [[v]] ds`. For vector elasticity, the lifting structure is identical but applied component-by-component:

```
Lift[elem][l, s, m] = (1/2) · M_K⁻¹[m, o] · Σ_q E[elem][o, q] · E[source][l, q] · n[s, q] · w[q]
```

where:
- `l` = source DOF index
- `s` = normal-direction index
- `m` = destination DOF index
- `M_K⁻¹` = scalar element mass matrix inverse
- `E[elem][m, q]` = shape function m evaluated at quadrature point q
- `n[s, q]` = face normal component s at quadrature point q

### 4.2 The `test_normal` Operator (Elasticity Tensor Coupling)

The key difference from scalar BR2 is the **elasticity tensor coupling** in the penalty term. Following Tandem's `elasticity.py` (lines 118-120), the `test_normal` operator is:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  test_normal_{iu,s} = λ · δ_us · n_i + μ · (δ_iu · n_s + δ_is · n_u)      │
└──────────────────────────────────────────────────────────────────────────────┘
```

This is the contraction of the elasticity tensor `C_{ijrs}` with the face normal: `test_normal_{iu,s} = C_{i_s_u_s'} n_{s'} δ_{ss'}` ... more precisely, it encodes the traction operator `[σ(·)·n]` in index form.

**Physical meaning**: `test_normal_{iu,s}` gives the `i`-th component of the traction on a face with normal `n`, when the displacement gradient has a unit entry in the `(u, s)` position.

**Reduction to scalar**: For the scalar Laplace case with diffusion coefficient K, `test_normal` reduces to `K · n_i`. This recovers the scalar `BR2InteriorFaceIntegrator` in `dg_br2_integrator.hpp`.

### 4.3 BR2 Lifted Flux (with Elasticity Coupling)

The combined lifted flux, incorporating the elasticity tensor, is:

```
L_q[source][l, i, u, q] = (1/2) · Σ_{elem∈{K⁻,K⁺}} Σ_{s,m} test_normal(elem)_{iu,s}
                           · E[elem][m, q] · Lift[elem][l, s, m]
```

This produces a `(ndof_source × dim × dim × nqp)` tensor. Indices:
- `l` = source DOF
- `i` = test function component (row of face matrix)
- `u` = trial function component (column of face matrix)
- `q` = quadrature point

### 4.4 BR2 Penalty Parameter

```
┌──────────────────────────────────────────────────────────────┐
│  σ^BR2 = number of faces per element                         │
│        = 2·dim = 6    (hexahedra)                            │
│        = dim+1 = 4    (tetrahedra)                           │
└──────────────────────────────────────────────────────────────┘
```

This follows Tandem's convention. Compared to the antiplane case (σ = dim+1 = 3 for triangles), the 3D hex case uses σ = 6.

**Implementation**: `sigma_ = 2 * dim_` in `DGElasticityBR2Integrator` constructor, with runtime check for tet elements.

### 4.5 BR2 Bilinear Form

The complete BR2 bilinear form `a(u_h, v_h)`:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  a(u_h, v_h) = a^vol(u_h, v_h) + a^int(u_h, v_h) + a^Dir(u_h, v_h)        │
└──────────────────────────────────────────────────────────────────────────────┘
```

where `a^int` sums over **all** interior faces (F_I ∪ F_F), including fault faces.

---

#### 4.5.1 Volume Term

```
┌──────────────────────────────────────────────────────────────┐
│  a^vol(u_h, v_h) = ∫_Ω σ(u_h) : ε(v_h) dx                  │
│                   = ∫_Ω [λ tr(ε(u_h)) tr(ε(v_h))            │
│                        + 2μ ε(u_h) : ε(v_h)] dx              │
└──────────────────────────────────────────────────────────────┘
```

| Term | File | Class |
|------|------|-------|
| Volume integral | `domain/elasticity_operator.hpp` | `ElasticityIntegrator` (MFEM built-in) |

---

#### 4.5.2 Interior Face Terms (Including Fault Faces)

For each interior face `e ∈ F_I ∪ F_F`:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  a^int(u_h, v_h) = Σ_{e ∈ F_I ∪ F_F} [a^cons_e + a^sym_e + a^lift_e]      │
└──────────────────────────────────────────────────────────────────────────────┘
```

**Consistency term** (ensures convergence):
```
a^cons_e = -∫_e {{σ(u_h) · n}} · [[v_h]] ds
```

**Symmetry term** (SIPG, ε = -1, makes bilinear form symmetric):
```
a^sym_e = ε · (-∫_e {{σ(v_h) · n}} · [[u_h]] ds)
        = +∫_e {{σ(v_h) · n}} · [[u_h]] ds       (for ε = -1)
```

**BR2 lifting term** (stabilization via lifted flux with elasticity coupling):
```
a^lift_e = σ^BR2 · ∫_e Σ_{i,u} [Σ_q w_q E_x[k,q] L_q[y][l,i,u,q]] ds
```

More explicitly, for the face matrix block `(x_elem, y_elem)`, component `(i, u)`:
```
a[(x,k,i), (y,l,u)] = c₀ · ∫ [σ_i(φ_x^k) · n]_u · φ_y^l ds     (consistency)
                     + c₁ · ∫ [σ_u(φ_y^l) · n]_i · φ_x^k ds     (symmetry)
                     + c₂ · Σ_q w_q E_x[k,q] L_q[y][l,i,u,q]   (BR2 penalty)
```

**Sign conventions** (following Tandem):

| Block | c₀ (consistency) | c₁ (symmetry) | c₂ (BR2 penalty) |
|-------|-------------------|----------------|-------------------|
| (0,0): elem1←elem1 | -0.5 | ε·0.5 = +0.5 | +σ^BR2 |
| (0,1): elem1←elem2 | +0.5 | ε·0.5 = +0.5 | -σ^BR2 |
| (1,0): elem2←elem1 | -0.5 | ε·(-0.5) = +0.5 | -σ^BR2 |
| (1,1): elem2←elem2 | +0.5 | ε·(-0.5) = +0.5 | +σ^BR2 |

Note: Block (1,0) consistency uses `-0.5` (not `+0.5`) because elem2 has normal sign `+1`, giving the same contribution form as block (0,0). Block (1,1) symmetry uses `ε·(-0.5) = +0.5` because ε = -1 and elem2's test function has the opposite sign convention.

| Term | File | Class |
|------|------|-------|
| Interior face integrator | `integrator/dg_elasticity_br2_integrator.hpp` | `DGElasticityBR2Integrator` |
| test_normal operator | `integrator/dg_elasticity_br2_integrator.hpp` | `TestNormal()` method |

---

#### 4.5.3 Dirichlet Boundary Terms

For each Dirichlet face `e ∈ F_D`:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  a^Dir(u_h, v_h) = Σ_{e ∈ F_D} [a^Dir,cons_e + a^Dir,sym_e + a^Dir,lift_e] │
└──────────────────────────────────────────────────────────────────────────────┘
```

**Dirichlet consistency**:
```
a^Dir,cons_e = -∫_e (σ(u_h) · n) · v_h ds
```

**Dirichlet symmetry**:
```
a^Dir,sym_e = ε · (-∫_e (σ(v_h) · n) · u_h ds)
```

**Dirichlet BR2 lifting** (full factor, no 0.5 averaging for boundary):
```
a^Dir,lift_e = σ^BR2 · ∫_K C : r_e(u_h) : r_e(v_h) dx
```

| Term | File | Class |
|------|------|-------|
| Boundary face integrator | `integrator/dg_elasticity_br2_integrator.hpp` | `DGElasticityBR2BoundaryIntegrator` |

---

### 4.6 BR2 Linear Form (Right-Hand Side)

The complete BR2 linear form:

```
┌──────────────────────────────────────────────────┐
│  L(v_h) = L^Dir(v_h) + L^Fault(v_h)              │
└──────────────────────────────────────────────────┘
```

(Neumann BC with σ·n = 0 contributes nothing.)

---

#### 4.6.1 Dirichlet Contribution

For prescribed displacement `g_D = (0, ±Vp·t/2, 0)` on Γ_D:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  L^Dir(v_h) = Σ_{e ∈ F_D} [ ε·(-∫_e (σ(v_h) · n) · g_D ds)               │
│                            + σ^BR2 · ∫_K C : r_e(g_D) : r_e(v_h) dx ]     │
└──────────────────────────────────────────────────────────────────────────────┘
```

Expanding component-wise for each test function `φ_k` with component `i`:

```
L^Dir_k,i = ε · Σ_u [σ_i(φ_k) · n]_u · g_D_u · w
          + σ^BR2 · φ_k · f_lifted_q[i, q] · w
```

where `f_lifted_q` is the BR2 lifted Dirichlet data with elasticity tensor coupling:

```
f_lifted_q[i, q] = Σ_{u,s,m} test_normal_{iu,s} · E[m, q] · Minv[m, o] · face_int[u·dim+s, o]
face_int[u·dim+s, o] = Σ_q E[o, q] · g_D[u] · n[s, q] · w[q]
```

Note: For boundary faces, the lifting uses full factor (no 0.5), unlike interior faces.

| Term | File | Function |
|------|------|----------|
| Dirichlet loading RHS | `domain/elasticity_operator.hpp` | `AssembleDirichletLoading()` |

---

#### 4.6.2 Fault Slip Contribution

For prescribed slip `δ` on fault faces F_F, the slip is first embedded from local to global frame:

```
δ_global = δ_dip · t₁ + δ_strike · t₂
```

Then the RHS contribution:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  L^Fault(v_h) = Σ_{e ∈ F_F} [ ε·(-∫_e {{σ(v_h) · n}} · δ ds)             │
│                              + σ^BR2 · ∫_{K⁻∪K⁺} C:r_e(δ)·r_e([[v_h]]) ] │
└──────────────────────────────────────────────────────────────────────────────┘
```

Expanding component-wise for each test function `φ_k` with component `i`:

**Element 1 (K⁻)**:
```
b^fault_{k,i} = c₁ · Σ_u [σ_i(φ_k) · n]_u · sign · δ_u · w / det(J₁)
              + σ^BR2 · φ_k · f_lifted_q[i, q] · w
```

**Element 2 (K⁺)**:
```
b^fault_{k,i} = c₁ · Σ_u [σ_i(φ_k) · n]_u · sign · δ_u · w / det(J₂)
              + (-σ^BR2) · φ_k · f_lifted_q[i, q] · w
```

where `c₁ = ε · 0.5 = +0.5` (SIPG symmetry coefficient), `sign` accounts for normal orientation relative to the fault convention, and `f_lifted_q` is the BR2 lifted slip with elasticity tensor coupling.

The BR2 lifted slip is computed as:

```
Step 1: face_int[elem][u·dim+s, m] = Σ_q E[m,q] · sign · δ_u · n[s,q] · w[q]

Step 2: f_lifted[elem][u·dim+s, m] = 0.5 · Minv[m, o] · face_int[elem][u·dim+s, o]

Step 3: f_lifted_q[i, q] = 0.5 · Σ_{elem,u,s,m} test_normal_{iu,s} · E[elem][m,q]
                            · f_lifted[elem][u·dim+s, m]
```

| Term | File | Function |
|------|------|----------|
| Fault slip RHS (BR2) | `domain/elasticity_operator.hpp` | `AssembleSlipContributionBR2()` |
| Fault slip RHS (IP) | `domain/elasticity_operator.hpp` | `AssembleSlipContributionIP()` |
| Slip embedding | `fault/fault_basis.hpp` | `FaultBasis::EmbedSlip()` |

---

### 4.7 BR2 Summary: Complete Discrete System

**Find** u_h ∈ V_h **such that**:

```
┌──────────────────────────────────────────────────┐
│  a(u_h, v_h) = L(v_h)      ∀ v_h ∈ V_h          │
└──────────────────────────────────────────────────┘
```

**Expanding all terms** (interior face terms sum over F_I ∪ F_F, following Arnold et al. eq. 3.24):

```
BILINEAR FORM (Left-Hand Side):
===============================

[Volume]
  ∫_Ω σ(u_h) : ε(v_h) dx

[Interior consistency — all interior faces including fault]
- Σ_{e ∈ F_I ∪ F_F} ∫_e {{σ(u_h) · n}} · [[v_h]] ds

[Interior symmetry — all interior faces including fault]    (ε = -1)
+ Σ_{e ∈ F_I ∪ F_F} ∫_e {{σ(v_h) · n}} · [[u_h]] ds

[Interior BR2 lifting — all interior faces including fault]
+ Σ_{e ∈ F_I ∪ F_F} σ^BR2 ∫ C : r_e([[u_h]]) : r_e([[v_h]]) dx

[Dirichlet consistency]
- Σ_{e ∈ F_D} ∫_e (σ(u_h) · n) · v_h ds

[Dirichlet symmetry]
+ Σ_{e ∈ F_D} ∫_e (σ(v_h) · n) · u_h ds

[Dirichlet BR2 lifting]
+ Σ_{e ∈ F_D} σ^BR2 ∫ C : r_e(u_h) : r_e(v_h) dx


LINEAR FORM (Right-Hand Side):
==============================

[Dirichlet RHS (symmetry)]
+ Σ_{e ∈ F_D} ∫_e (σ(v_h) · n) · g_D ds

[Dirichlet RHS (BR2 lifting)]
+ Σ_{e ∈ F_D} σ^BR2 ∫ C : r_e(g_D) : r_e(v_h) dx

[Fault RHS (symmetry)]
+ Σ_{e ∈ F_F} ε · (-∫_e {{σ(v_h) · n}} · δ ds)

[Fault RHS (BR2 lifting)]
+ Σ_{e ∈ F_F} σ^BR2 ∫ C : r_e(δ) : r_e([[v_h]]) dx
```

**How fault enforcement works**: Same mechanism as antiplane. On fault faces `e ∈ F_F`, the bilinear form contributes LHS terms with `[[u_h]]` and the RHS contributes terms with `δ`. The net effect is:

- Symmetry: `+{{σ(v)·n}}·[[u_h]]` (LHS) vs `−{{σ(v)·n}}·δ` (RHS) → enforces `[[u_h]] → δ`
- Lifting: `σ r_e([[u_h]])·r_e([[v_h]])` (LHS) vs `σ r_e(δ)·r_e([[v_h]])` (RHS) → penalizes `[[u_h]] - δ`
- Consistency: `-{{σ(u)·n}}·[[v_h]]` (LHS only) → provides flux coupling across fault

---

## 5. IP Method: Complete Formulation

The **Interior Penalty (IP)** method is the alternative DG method, using MFEM's built-in `DGElasticityIntegrator`.

### 5.1 IP Penalty Parameter

```
┌─────────────────────────────────┐
│  κ = (p + 1)²                    │
└─────────────────────────────────┘
```

where p = polynomial degree. MFEM's `DGElasticityIntegrator` handles the h-dependent scaling internally.

### 5.2 IP Bilinear Form

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  a(u_h, v_h) = a^vol + a^int + a^Dir                                        │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

#### 5.2.1 Volume Term

Same as BR2:
```
a^vol(u_h, v_h) = ∫_Ω σ(u_h) : ε(v_h) dx
```

---

#### 5.2.2 Interior Face Terms

For each interior face `e ∈ F_I ∪ F_F`:

**Consistency** (same as BR2):
```
a^cons_e = -∫_e {{σ(u_h) · n}} · [[v_h]] ds
```

**Symmetry** (same as BR2):
```
a^sym_e = ε · (-∫_e {{σ(v_h) · n}} · [[u_h]] ds)
```

**IP penalty** (replaces BR2 lifting):
```
a^pen_e = κ · |n|² · (w₁ + w₂) · ∫_e [[u_h]] · [[v_h]] ds
```

where `w_i = 1/(2·det(J_i))` and `|n|² = n · n` (unnormalized face normal).

| Term | File | Class |
|------|------|-------|
| Interior face integrator | MFEM built-in | `DGElasticityIntegrator` |

---

#### 5.2.3 Dirichlet Boundary Terms

Same structure as interior, with full penalty (single element):

```
a^Dir,pen_e = κ · |n|² · w · ∫_e u_h · v_h ds
```

| Term | File | Class |
|------|------|-------|
| Boundary face integrator | MFEM built-in | `DGElasticityIntegrator` |

---

### 5.3 IP Linear Form (Right-Hand Side)

#### 5.3.1 Dirichlet Contribution

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  L^Dir(v_h) = Σ_{e ∈ F_D} [ ε · (-∫_e (σ(v_h) · n) · g_D ds)             │
│                            + κ · |n|² · w · ∫_e g_D · v_h ds ]            │
└──────────────────────────────────────────────────────────────────────────────┘
```

#### 5.3.2 Fault Slip Contribution

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  L^Fault(v_h) = Σ_{e ∈ F_F} [ ε · (-∫_e {{σ(v_h) · n}} · δ ds)           │
│                              + κ · |n|² · (w₁+w₂) · ∫_e δ · [[v_h]] ds ] │
└──────────────────────────────────────────────────────────────────────────────┘
```

| Term | File | Function |
|------|------|----------|
| Fault slip RHS (IP) | `domain/elasticity_operator.hpp` | `AssembleSlipContributionIP()` |

---

### 5.4 IP Summary: Complete Discrete System

```
BILINEAR FORM (Left-Hand Side):
===============================

[Volume]
  ∫_Ω σ(u_h) : ε(v_h) dx

[Interior consistency — all interior faces including fault]
- Σ_{e ∈ F_I ∪ F_F} ∫_e {{σ(u_h) · n}} · [[v_h]] ds

[Interior symmetry — all interior faces including fault]
+ Σ_{e ∈ F_I ∪ F_F} ∫_e {{σ(v_h) · n}} · [[u_h]] ds

[Interior IP penalty — all interior faces including fault]
+ Σ_{e ∈ F_I ∪ F_F} κ |n|² (w₁+w₂) ∫_e [[u_h]] · [[v_h]] ds

[Dirichlet consistency]
- Σ_{e ∈ F_D} ∫_e (σ(u_h) · n) · v_h ds

[Dirichlet symmetry]
+ Σ_{e ∈ F_D} ∫_e (σ(v_h) · n) · u_h ds

[Dirichlet IP penalty]
+ Σ_{e ∈ F_D} κ |n|² w ∫_e u_h · v_h ds


LINEAR FORM (Right-Hand Side):
==============================

[Dirichlet RHS (symmetry)]
+ Σ_{e ∈ F_D} ∫_e (σ(v_h) · n) · g_D ds

[Dirichlet RHS (penalty)]
+ Σ_{e ∈ F_D} κ |n|² w ∫_e g_D · v_h ds

[Fault RHS (symmetry)]
+ Σ_{e ∈ F_F} ε · (-∫_e {{σ(v_h) · n}} · δ ds)

[Fault RHS (penalty)]
+ Σ_{e ∈ F_F} κ |n|² (w₁+w₂) ∫_e δ · [[v_h]] ds
```

---

## 6. Traction Computation

### 6.1 Physical Traction on Fault

The traction vector on the fault is:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  T_i = σ_ij · n_j                                                           │
│      = [λ tr(ε) δ_ij + 2μ ε_ij] · n_j                                      │
│      = λ (∂u_k/∂x_k) n_i + μ (∂u_i/∂x_j + ∂u_j/∂x_i) n_j                │
└──────────────────────────────────────────────────────────────────────────────┘
```

where:
- `σ` is computed from the **averaged** gradient across the fault:
  `{{∇u}} = (∇u⁺ + ∇u⁻) / 2`
- `n` is the fault normal (from `FaultBasis`)

### 6.2 Traction Computation Steps

1. **Evaluate ∇u on both sides**: For each fault face, compute the displacement gradient from each adjacent element using physical derivatives:
   ```
   (∇u)_cd = Σ_k (dφ_k/dx_d) · u_c,k
   ```

2. **Average**: `{{∇u}}_{cd} = (∇u⁺_{cd} + ∇u⁻_{cd}) / 2`

3. **Strain**: `ε_{ij} = ({{∇u}}_{ij} + {{∇u}}_{ji}) / 2`

4. **Stress**: `σ_{ij} = λ tr(ε) δ_{ij} + 2μ ε_{ij}`

5. **Global traction**: `T_i^{global} = Σ_j σ_{ij} · n_j`

6. **Project to local frame**: `τ_local = FaultBasis::ProjectTraction(T_global)`
   - `τ_dip = T_global · t₁`
   - `τ_strike = T_global · t₂`

### 6.3 Traction Output Layout

The traction vector has size `2N` (two components per fault DOF):

```
traction = [τ_dip_0, τ_strike_0, τ_dip_1, τ_strike_1, ..., τ_dip_{N-1}, τ_strike_{N-1}]
```

### 6.4 Comparison with Antiplane Traction

| Aspect | Antiplane | Full Elasticity |
|--------|-----------|-----------------|
| Formula | `τ = μ {{∂u/∂x}}` | `T = {{σ(u)}} · n` |
| Components | Scalar | Vector (2: dip, strike) |
| Normal needed | Implicit (fault at x=0) | From `FaultBasis` |
| Coordinate transform | None needed | `ProjectTraction()` |
| Output size | `N` | `2N` |

| Component | File | Function |
|-----------|------|----------|
| Traction computation | `domain/elasticity_operator.hpp` | `ComputeTraction()` |
| Gradient evaluation | `domain/elasticity_operator.hpp` | Inline (CalcDShape + CalcInverse) |
| Traction projection | `fault/fault_basis.hpp` | `FaultBasis::ProjectTraction()` |

---

## 7. Fault Coordinate Transformation

### 7.1 The Two Coordinate Frames

**Global frame (x₁, x₂, x₃)**: Defined by the 3D mesh. All PDE quantities live here — displacement `u`, stress `σ`, displacement jumps `Δu`.

**Fault-local frame (n, t₁, t₂)**: Defined per face on the fault surface. All fault physics quantities live here — tangential traction `τ[2]`, slip rate `V[2]`, accumulated slip `s[2]`, friction parameters.

### 7.2 Basis Construction (FaultBasis)

For each fault face, following Tandem's `Curvilinear::facetBasis()`:

```
1. Raw normal:  n_raw = CalcOrtho(J_face)
2. Orient:      if (n_raw · ref_normal < 0) → n_raw = -n_raw
3. Normalize:   n = n_raw / |n_raw|
4. Strike:      t₂ = normalize(up × n)
5. Dip:         t₁ = t₂ × n
```

For BP5 setup: `ref_normal = (1, 0, 0)`, `up = (0, 0, -1)`.
Result: `n ≈ (1,0,0)`, `t₁ ≈ (0,0,1)` (dip), `t₂ ≈ (0,1,0)` (strike).

### 7.3 The Two Transforms

**Traction projection** (global → local):
```
τ_dip    = T_global · t₁
τ_strike = T_global · t₂
```

**Slip embedding** (local → global):
```
Δu = s_dip · t₁ + s_strike · t₂      (no-opening: Δu · n = 0)
```

### 7.4 Data Flow

```
                    GLOBAL FRAME                      LOCAL FRAME
                   (mesh x₁,x₂,x₃)                (n, t₁=dip, t₂=strike)

  ┌───────────────────────────────────────────────────────────────────────┐
  │ ElasticityDomainOperator::Solve()                                     │
  │                                                                       │
  │   slip_local[2] ──EmbedSlip──→ Δu[3] ──→ DG face RHS assembly        │
  └───────────────────────────────────────────────────────────────────────┘

  ┌───────────────────────────────────────────────────────────────────────┐
  │ ElasticityDomainOperator::ComputeTraction()                           │
  │                                                                       │
  │   ∇u → ε → σ → σ·n = T_global[3] ──ProjectTraction──→ τ_local[2]    │
  └───────────────────────────────────────────────────────────────────────┘

  ┌───────────────────────────────────────────────────────────────────────┐
  │ RateStateFaultOperator::ComputeRHS()  (entirely in local frame)       │
  │                                                                       │
  │   τ_local[2] + τ₀[2] → SolveSlipRateVectorPsi → V_local[2], dψ/dt   │
  └───────────────────────────────────────────────────────────────────────┘
```

The domain operator is the **only place** that touches both frames. Everything above the fault interface (friction, state evolution, parameters) works entirely in the local frame.

| Component | File | Function |
|-----------|------|----------|
| Basis computation | `fault/fault_basis.hpp` | `FaultBasis::Compute()` |
| Traction projection | `fault/fault_basis.hpp` | `FaultBasis::ProjectTraction()` |
| Slip embedding | `fault/fault_basis.hpp` | `FaultBasis::EmbedSlip()` |
| Normal stress | `fault/fault_basis.hpp` | `FaultBasis::NormalStress()` |

---

### 7.5 Rate-and-State Friction (Vector Formulation)

In the fault-local frame, the rate-and-state friction law extends to 2D:

**Scalar equation** (solved by Newton's method):
```
|τ| = σ_n · f(|V|, ψ) + η · |V|

where f(|V|, ψ) = a · asinh[(|V| / 2V₀) · exp(ψ/a)]
```

**Vector decomposition**:
```
1. |τ| = sqrt(τ_dip² + τ_strike²)
2. Solve scalar: |V| = SolveSlipRatePsi(|τ|, ψ, σ_n, η, a)
3. Project: V = -(|V| / |τ|) · τ      (anti-parallel to traction)
```

**State evolution** (uses scalar |V|):
```
dψ/dt = (b V₀ / Dc) [exp((f₀ - ψ)/b) - |V|/V₀]
```

| Component | File | Function |
|-----------|------|----------|
| Vector slip rate | `friction/dieterich_ruina.hpp` | `SolveSlipRateVectorPsi()` |
| Scalar Newton solve | `friction/dieterich_ruina.hpp` | `SolveSlipRatePsi()` |
| State evolution | `friction/state_evolution.hpp` | `AgingLawPsi::Rate()` |

---

### 7.6 Mass Matrix Inverse (for BR2)

The vector DG space `V_h` has `vdim = 3` copies of a scalar DG space. The scalar mass matrix `M_K` is the same for all components. Therefore, the BR2 lifting uses the **scalar** mass matrix inverse, precomputed once per element.

```
M_K[m, o] = ∫_K φ_m · φ_o dx

M_K⁻¹ is computed and stored for each element.
```

| Component | File | Function |
|-----------|------|----------|
| Mass inverse precomputation | `domain/elasticity_operator.hpp` | `PrecomputeMassInverse()` |
| Scalar FE space | `domain/elasticity_operator.hpp` | `scalar_fes_` |

---

# Part II: Implementation Mapping

---

## 8. File-to-Formula Mapping

### 8.1 Source File Overview

```
miniapps/seas/
├── config/
│   ├── bp2_params.hpp                # BP2 benchmark parameters
│   └── bp5_params.hpp                # BP5 benchmark parameters (2D a, L, tau0, V_init)
├── common/
│   ├── seas_types.hpp                # Shared type definitions, DG method enum
│   └── mpi_context.hpp               # MPI context utilities
├── domain/
│   ├── domain_operator.hpp           # Abstract DomainOperator interface
│   ├── antiplane_operator.hpp        # AntiplaneDomainOperator (scalar, BP1/BP2)
│   └── elasticity_operator.hpp       # ElasticityDomainOperator (vector, BP5)
├── integrator/
│   ├── dg_br2_integrator.hpp         # Scalar BR2 face integrators (Laplace)
│   └── dg_elasticity_br2_integrator.hpp  # Vector BR2 face integrators (elasticity)
├── fault/
│   ├── fault_basis.hpp               # FaultBasis: per-face (n, t1, t2) coordinate frames
│   ├── fault_geometry.hpp            # FaultGeometry: spatial parameters (BP2 1D, BP5 2D)
│   └── rate_state_fault.hpp          # RateStateFaultOperator (scalar + vector)
├── friction/
│   ├── friction_law.hpp              # Abstract FrictionLaw interface
│   ├── dieterich_ruina.hpp           # DieterichRuinaFriction (scalar + vector solver)
│   └── state_evolution.hpp           # AgingLaw, SlipLaw, AgingLawPsi
├── solver/
│   └── seas_operator.hpp             # SEASQuasiDynamicOperator (time integration coupling)
└── tests/unit/
    ├── test_elasticity_br2.cpp       # BR2 integrator tests (14 tests)
    ├── test_elasticity_operator.cpp   # Elasticity operator tests (32 tests)
    ├── test_fault_basis.cpp           # FaultBasis tests (113 tests)
    ├── test_bp5_params.cpp            # BP5 parameter tests (77 tests)
    └── test_vector_friction.cpp       # Vector friction tests (45 tests)
```

### 8.2 BR2 Method: Term-to-Code Mapping

| Formula Term | File | Class/Function | Status |
|--------------|------|----------------|--------|
| **Bilinear Form** | | | |
| Volume: ∫ σ(u):ε(v) | `elasticity_operator.hpp` | `AssembleStiffness()` → `ElasticityIntegrator` | Done |
| Interior consistency: -∫ {{σ(u)·n}}·[[v]] | `dg_elasticity_br2_integrator.hpp` | `DGElasticityBR2Integrator` | Done |
| Interior symmetry: ε·(-∫ {{σ(v)·n}}·[[u]]) | `dg_elasticity_br2_integrator.hpp` | `DGElasticityBR2Integrator` | Done |
| Interior BR2 lifting: σ ∫ C:r_e([[u]]):r_e([[v]]) | `dg_elasticity_br2_integrator.hpp` | `DGElasticityBR2Integrator` | Done |
| test_normal operator | `dg_elasticity_br2_integrator.hpp` | `TestNormal()` | Done |
| Dirichlet terms | `dg_elasticity_br2_integrator.hpp` | `DGElasticityBR2BoundaryIntegrator` | Done |
| **Linear Form** | | | |
| Dirichlet loading RHS | `elasticity_operator.hpp` | `AssembleDirichletLoading()` | Done |
| Fault slip RHS (BR2) | `elasticity_operator.hpp` | `AssembleSlipContributionBR2()` | Done |
| **Supporting** | | | |
| Mass matrix inverse | `elasticity_operator.hpp` | `PrecomputeMassInverse()` | Done |
| Scalar FE space (for Minv) | `elasticity_operator.hpp` | `scalar_fes_` | Done |
| Penalty (σ = 2·dim) | `dg_elasticity_br2_integrator.hpp` | Constructor | Done |
| Fault face detection | `elasticity_operator.hpp` | `SetupFaultInfo()` | Done |
| FaultBasis computation | `fault_basis.hpp` | `FaultBasis::Compute()` | Done |

### 8.3 IP Method: Term-to-Code Mapping

| Formula Term | File | Class/Function | Status |
|--------------|------|----------------|--------|
| **Bilinear Form** | | | |
| Volume: ∫ σ(u):ε(v) | `elasticity_operator.hpp` | `ElasticityIntegrator` | Done |
| Interior consistency + symmetry + penalty | MFEM built-in | `DGElasticityIntegrator` | Built-in |
| Dirichlet terms | MFEM built-in | `DGElasticityIntegrator` | Built-in |
| **Linear Form** | | | |
| Dirichlet loading RHS | `elasticity_operator.hpp` | `AssembleDirichletLoading()` | Done |
| Fault slip RHS (IP) | `elasticity_operator.hpp` | `AssembleSlipContributionIP()` | Done |
| **Supporting** | | | |
| Penalty κ = (p+1)² | `elasticity_operator.hpp` | `AssembleStiffness()` | Done |

### 8.4 Cross-Cutting Components

| Component | File | Description |
|-----------|------|-------------|
| DG method selection | `antiplane_operator.hpp` | `enum class DGMethod { BR2, IP }` |
| Vector DG FE space | `elasticity_operator.hpp` | `DG_FECollection`, vdim=3, byNODES |
| Coordinate transforms | `fault_basis.hpp` | `FaultBasis` class |
| Traction computation | `elasticity_operator.hpp` | `ComputeTraction()` |
| Solver setup | `elasticity_operator.hpp` | `SetupSolver()` (CG + GS preconditioner) |
| Parallel support | `elasticity_operator.hpp` | Template `<MeshType>` |

### 8.5 Scalar-to-Vector Generalization

| Scalar (Antiplane) | Vector (Elasticity) | Key Difference |
|--------------------|---------------------|----------------|
| `DiffusionIntegrator` | `ElasticityIntegrator` | Stress tensor coupling |
| `DGDiffusionIntegrator` | `DGElasticityIntegrator` | Vector traction operator |
| `BR2InteriorFaceIntegrator` | `DGElasticityBR2Integrator` | `test_normal` coupling |
| `K · n_i` | `test_normal_{iu,s}` | Elasticity tensor vs scalar coeff |
| Scalar lifting | Component-wise lifting + tensor coupling | dim² more terms per DOF |
| `AssembleSlipContribution()` | `AssembleSlipContributionBR2/IP()` | 3-component δ, traction operator |
| `ComputeTraction()` → scalar τ | `ComputeTraction()` → (τ_dip, τ_strike) | Stress tensor + FaultBasis |

---

# Appendix: Tandem Reference

---

## Appendix A: Tandem Implementation Details

### A.1 BR2 Penalty (from Tandem Elasticity)

**From `Elasticity.h`**:
```cpp
double penalty(std::size_t fctNo) const {
    if (method_ == DGMethod::BR2) {
        return NumFacets;  // = 2*Dim for hex, Dim+1 for tet
    }
    return penalty_[fctNo];
}
```

For 3D hexahedra: `NumFacets = 6 = 2*3`.

### A.2 test_normal Operator (from Tandem)

**From `elasticity.py` lines 118-120**:
```python
test_normal[0][x]['iu_sq'] <= lam[0][x]['q'] * unit['us'] * n_unit[x]['iq'] + \
    mu[0][x]['q'] * (unit['iu'] * n_unit[x]['sq'] + unit['is'] * n_unit[x]['uq'])
test_normal[1][x]['iu_sq'] <= lam[1][x]['q'] * unit['us'] * n_unit[x]['iq'] + \
    mu[1][x]['q'] * (unit['iu'] * n_unit[x]['sq'] + unit['is'] * n_unit[x]['uq'])
```

This gives:
```
test_normal_{iu,sq} = λ · δ_us · n_i + μ · (δ_iu · n_s + δ_is · n_u)
```

### A.3 Lifting Operator (from Tandem Elasticity)

**From `elasticity.py` lines 87-93**:
```python
# Scalar part of lifting (same structure as Poisson)
Lift[0]['liu'] <= 0.5 * Minv[0]['us'] * E_q[0]['sq'] * E_q[x]['lq'] * n_q['iq'] * w['q']
Lift[1]['liv'] <= 0.5 * Minv[1]['vs'] * E_q[1]['sq'] * E_q[x]['lq'] * n_q['iq'] * w['q']

# Combined lifted flux with elasticity tensor coupling
L_q[x]['lq'] <= 0.5 * (test_normal[0][0]['iu_sq'] * E_q[0]['uq'] * Lift[0]['lis'] +
                        test_normal[1][0]['iu_sq'] * E_q[1]['vq'] * Lift[1]['lis'] +
                        test_normal[0][1]['iu_sq'] * E_q[0]['uq'] * Lift[0]['lis'] +
                        test_normal[1][1]['iu_sq'] * E_q[1]['vq'] * Lift[1]['lis'])
```

### A.4 Assembly Coefficients (from Tandem Elasticity)

**From `Elasticity.cpp`**:
```cpp
// Same sign conventions as Poisson
assemble.c00 = -0.5;              // flux from side 0 on test from side 0
assemble.c01 = +0.5;              // flux from side 1 on test from side 0
assemble.c10 = epsilon * 0.5;     // = -0.5 for SIPG
assemble.c11 = +0.5;
assemble.c20 = penalty(fctNo);    // = 6 for hex
assemble.c21 = -penalty(fctNo);   // opposite sign for cross-element
```

### A.5 Fault Slip Embedding (from Tandem)

**From `AdapterBase.cpp`**:
```cpp
// Tandem stores basis as D×D matrix (columns = n, t1, t2)
// Slip embedding: delta_u = basis * (0, s_dip, s_strike)^T
// Traction projection: tau_local = basis^T * tau_global
```

### A.6 BP5 Dirichlet Loading (from Tandem)

**From `bp5.lua`**:
```lua
boundary = {
    {5, "fault", {}},
    {1, "dirichlet", {0, Vh, 0}},    -- +x2 wall
    {2, "dirichlet", {0, -Vh, 0}},   -- -x2 wall
    {3, "natural", {}},               -- all others
}
```

where `Vh = Vp * t / 2`.

---

## Appendix B: Tandem Source File Reference

### B.1 Key Source Files

```
tandem/
├── app/
│   ├── localoperator/
│   │   ├── Elasticity.h           # Elasticity operator header
│   │   ├── Elasticity.cpp         # Full implementation
│   │   └── DieterichRuinaAgeing.h # Friction law (scalar + vector)
│   ├── kernels/
│   │   ├── elasticity.py          # Tensor kernel definitions (test_normal, lifting)
│   │   └── poisson.py             # Scalar kernel definitions (for comparison)
│   └── form/
│       ├── BC.h                   # Boundary condition enum
│       └── FacetInfo.h            # Face information structure
└── src/
    ├── form/
    │   ├── DGCurvilinearCommon.h  # DG method enum
    │   └── InverseInequality.h    # Penalty constants
    └── geometry/
        └── Curvilinear.cpp        # facetBasis() — per-face coordinate frames
```

### B.2 Critical Implementation Locations

| Component | File | Lines |
|-----------|------|-------|
| BR2 penalty | `Elasticity.h` | penalty() method |
| IP penalty | `Elasticity.cpp` | prepare_penalty() |
| test_normal | `elasticity.py` | 118-120 |
| Lifting operator | `elasticity.py` | 87-93 |
| Volume kernel | `elasticity.py` | 30-35 |
| Skeleton assembly | `Elasticity.cpp` | assemble_skeleton() |
| Boundary assembly | `Elasticity.cpp` | assemble_boundary() |
| Fault RHS | `elasticity.py` | 140+ |
| Coordinate frame | `Curvilinear.cpp` | facetBasis() |

---

## References

1. **Arnold, D.N., Brezzi, F., Cockburn, B., and Marini, L.D.** (2002). "Unified Analysis of Discontinuous Galerkin Methods for Elliptic Problems." *SIAM J. Numer. Anal.* 39(5), 1749-1779.

2. **Bassi, F. and Rebay, S.** (1997). "A high-order accurate discontinuous finite element method for the numerical solution of the compressible Navier-Stokes equations." *J. Comput. Phys.* 131, 267-279.

3. **Uphoff, C., May, D.A., Gabriel, A.-A.** (2023). "A discontinuous Galerkin method for sequences of earthquakes and aseismic slip on multiple faults using unstructured curvilinear grids." *Geophys. J. Int.* 233(1), 586-626.

4. **Tandem Code**: https://github.com/TEAR-ERC/tandem

5. **SCEC SEAS Benchmarks**: https://strike.scec.org/cvws/seas/

6. **MFEM Examples**: `ex17.cpp` (DG elasticity), `ex14.cpp` (DG diffusion)

7. **Dieterich, J.H.** (1979). "Modeling of rock friction: 1. Experimental results and constitutive equations." *J. Geophys. Res.* 84, 2161-2168.

8. **Ruina, A.** (1983). "Slip instability and state variable friction laws." *J. Geophys. Res.* 88, 10359-10370.
