# Discontinuous Galerkin Formulation for Antiplane Shear

This document provides the complete mathematical formulation for the SEAS miniapp DG solver, with clear mapping to implementation files.

---

## Table of Contents

**Part I: Mathematical Formulation**
1. [Strong Form (Governing Equation)](#1-strong-form-governing-equation)
2. [Computational Domain and Boundary Conditions](#2-computational-domain-and-boundary-conditions)
3. [DG Notation](#3-dg-notation)
4. [BR2 Method: Complete Formulation](#4-br2-method-complete-formulation)
5. [IP Method: Complete Formulation](#5-ip-method-complete-formulation)
6. [Traction Computation](#6-traction-computation)

**Part II: Implementation Mapping**
7. [File-to-Formula Mapping](#7-file-to-formula-mapping)

**Appendix: Tandem Reference**
- [A. Tandem Implementation Details](#appendix-a-tandem-implementation-details)
- [B. Tandem Source File Reference](#appendix-b-tandem-source-file-reference)

---

# Part I: Mathematical Formulation

---

## 1. Strong Form (Governing Equation)

### 1.1 Physical Problem: Antiplane Shear

**Displacement field** (Mode III / out-of-plane shear):

```
u = (0, 0, u(x, z))
```

where `u(x,z)` is the out-of-plane (y-direction) displacement.

**Stress-strain relations** (isotropic elastic, shear modulus μ):

```
σ_xy = μ ∂u/∂x
σ_zy = μ ∂u/∂z
```

**Equilibrium** (quasi-static, no body forces):

```
∂σ_xy/∂x + ∂σ_zy/∂z = 0
```

### 1.2 Strong Form PDE

Substituting Hooke's law into equilibrium gives the **Laplace equation**:

```
┌─────────────────────────────────────┐
│  -∇·(K ∇u) = 0      in Ω           │
└─────────────────────────────────────┘
```

where:
- `K = 1` (unit diffusion coefficient; physical μ applied to traction computation)
- `Ω ⊂ ℝ²` is the computational domain

### 1.3 Boundary Conditions

**Dirichlet BC** (prescribed displacement):
```
u = g_D      on Γ_D
```

**Neumann BC** (prescribed traction):
```
K ∇u · n = g_N      on Γ_N
```

**Interior fault condition** (prescribed slip jump):
```
[[u]] = u⁺ - u⁻ = δ(z, t)      on Γ_F
```

---

## 2. Computational Domain and Boundary Conditions

### 2.1 Domain Geometry

```
        z = 0 (free surface)
        Neumann: K∇u · n = 0
    ┌─────────────────────────────────────────────────┐
    │                        │                        │
    │                        │                        │
    │     Left half          │      Right half        │
    │     (x < 0)            │      (x > 0)           │
    │                        │                        │
x=-Lx                        │x=0                   x=+Lx
Neumann                      │                     Neumann
K∇u·n = 0                 FAULT                   K∇u·n = 0
    │                  [[u]] = δ                     │
    │                        │                        │
    │                        │                        │
    └─────────────────────────────────────────────────┘
        z = -Lz (bottom)
        Dirichlet: u = sign(x)·Vp·t/2
```

### 2.2 Boundary Classification (BP2)

| Boundary | Location | BC Type | Value | Implementation |
|----------|----------|---------|-------|----------------|
| Free surface | z = 0 | Neumann | g_N = 0 | `antiplane_operator.hpp` |
| Far-field left | x = -L_x | Neumann | g_N = 0 | `antiplane_operator.hpp` |
| Far-field right | x = +L_x | Neumann | g_N = 0 | `antiplane_operator.hpp` |
| Bottom | z = -L_z | Dirichlet | g_D = sign(x)·V_p·t/2 | `antiplane_operator.hpp` |
| Fault | x = 0, z ∈ [-W_f, 0] | Jump | [[u]] = δ | `rate_state_fault.hpp` |

---

## 3. DG Notation

### 3.1 Mesh Decomposition

- **Elements**: Ω = ∪_K K (triangulation)
- **Interior faces**: F_I = interior faces shared by two elements, excluding fault faces
- **Fault faces**: F_F = interior faces with prescribed slip (F_F ⊂ interior faces, F_F ∩ F_I = ∅)
- **Dirichlet faces**: F_D = faces on Γ_D
- **Neumann faces**: F_N = faces on Γ_N

**Note**: F_I and F_F are disjoint subsets of all interior faces. Following Arnold et al. (2002)
eq. 3.24, the bilinear form's interior face terms sum over **F_I ∪ F_F** (all interior faces).
The fault contribution to the RHS provides the inhomogeneous shift to enforce [[u]] = δ.

### 3.2 Jump and Average Operators

For an interior face `e` shared by elements K⁺ and K⁻:

**Jump**:
```
[[v]] = v⁺ - v⁻
```

**Average**:
```
{{v}} = (v⁺ + v⁻) / 2
```

**Average of flux**:
```
{{K ∇v · n}} = (K⁺ ∇v⁺ · n + K⁻ ∇v⁻ · n) / 2
```

For boundary faces: `[[v]] = v`, `{{v}} = v`.

---

## 4. BR2 Method: Complete Formulation

The **Bassi-Rebay 2 (BR2)** method uses lifting operators instead of direct penalty terms.

### 4.1 Lifting Operator Definition

For an interior face `e` shared by elements K⁻ and K⁺, the lifting operator `r_e : L²(e) → [Pᵖ(K⁻)]ᵈ × [Pᵖ(K⁺)]ᵈ` satisfies:

```
┌──────────────────────────────────────────────────────────────────────────┐
│  ∫_{K⁻ ∪ K⁺} r_e([[v]]) · w dx = ∫_e {{w · n}} [[v]] ds    ∀w          │
└──────────────────────────────────────────────────────────────────────────┘
```

**Local computation**:
```
r_e|_{K±} = M_{K±}⁻¹ ∫_e φ ⊗ n± [[v]] ds
```

where `M_K` is the element mass matrix.

**Implementation**: `domain/antiplane_operator.hpp` → `ComputeLifting()`

### 4.2 BR2 Penalty Parameter

```
┌─────────────────────────┐
│  σ^BR2 = D + 1          │
└─────────────────────────┘
```

where D = spatial dimension:
- **2D**: σ = 3
- **3D**: σ = 4

**Implementation**: `domain/antiplane_operator.hpp` → `SetupPenalty()`

### 4.3 BR2 Bilinear Form

The complete BR2 bilinear form `a(u_h, v_h)` consists of:

```
┌─────────────────────────────────────────────────────────────────────────┐
│  a(u_h, v_h) = a^vol(u_h, v_h) + a^int(u_h, v_h) + a^Dir(u_h, v_h)     │
└─────────────────────────────────────────────────────────────────────────┘
```

where `a^int` sums over **all** interior faces (F_I ∪ F_F), including fault faces.

---

#### 4.3.1 Volume Term

```
┌─────────────────────────────────────────────────────────────┐
│  a^vol(u_h, v_h) = ∫_Ω K ∇u_h · ∇v_h dx                    │
└─────────────────────────────────────────────────────────────┘
```

| Term | File | Function/Class |
|------|------|----------------|
| Volume integral | `domain/antiplane_operator.hpp` | `SetupBilinearForm()` |
| MFEM integrator | MFEM built-in | `DiffusionIntegrator` |

---

#### 4.3.2 Interior Face Terms (Including Fault Faces)

For each interior face `e ∈ F_I ∪ F_F` (all interior faces, including fault):

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  a^int(u_h, v_h) = Σ_{e ∈ F_I ∪ F_F} [ a^cons_e + a^sym_e + a^lift_e ]             │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

**Note**: Fault faces F_F receive the same bilinear form treatment as regular interior
faces F_I. This is consistent with Arnold et al. (2002) eq. 3.24, where the face terms
sum over Γ (all element boundaries), and with Tandem's implementation where
`assemble_skeleton` is called unconditionally for all interior faces regardless of BC type.
The fault RHS terms (Section 4.4.2) then shift the penalty enforcement from [[u]] = 0
to [[u]] = δ.

**Consistency term** (ensures convergence):
```
a^cons_e = -∫_e {{K ∇u_h · n}} [[v_h]] ds
```

**Symmetry term** (makes bilinear form symmetric):
```
a^sym_e = -∫_e {{K ∇v_h · n}} [[u_h]] ds
```

**BR2 lifting term** (stabilization via lifting operator):
```
a^lift_e = σ ∫_{K⁻ ∪ K⁺} K r_e([[u_h]]) · r_e([[v_h]]) dx
```

| Term | File | Function/Class |
|------|------|----------------|
| Interior face loop | `domain/antiplane_operator.hpp` | `AssembleInteriorFaces()` |
| Consistency term | `domain/dg_integrators.hpp` | `BR2ConsistencyIntegrator` |
| Symmetry term | `domain/dg_integrators.hpp` | `BR2SymmetryIntegrator` |
| Lifting term | `domain/dg_integrators.hpp` | `BR2LiftingIntegrator` |
| Mass matrix inverse | `domain/antiplane_operator.hpp` | `SetupMassInverse()` |

---

#### 4.3.3 Dirichlet Boundary Terms

For each Dirichlet face `e ∈ F_D`:

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│  a^Dir(u_h, v_h) = Σ_{e ∈ F_D} [ a^Dir,cons_e + a^Dir,sym_e + a^Dir,lift_e ]        │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

**Dirichlet consistency**:
```
a^Dir,cons_e = -∫_e (K ∇u_h · n) v_h ds
```

**Dirichlet symmetry**:
```
a^Dir,sym_e = -∫_e (K ∇v_h · n) u_h ds
```

**Dirichlet BR2 lifting**:
```
a^Dir,lift_e = σ ∫_K K r_e(u_h) · r_e(v_h) dx
```

| Term | File | Function/Class |
|------|------|----------------|
| Dirichlet BC loop | `domain/antiplane_operator.hpp` | `AssembleDirichletFaces()` |
| All three terms | `domain/dg_integrators.hpp` | `BR2DirichletIntegrator` |

---

### 4.4 BR2 Linear Form (Right-Hand Side)

The complete BR2 linear form `L(v_h)` consists of:

```
┌─────────────────────────────────────────────────────────┐
│  L(v_h) = L^Dir(v_h) + L^Fault(v_h)                    │
└─────────────────────────────────────────────────────────┘
```

(Neumann BC with g_N = 0 contributes nothing.)

---

#### 4.4.1 Dirichlet Contribution

For prescribed displacement `g_D` on Γ_D:

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  L^Dir(v_h) = Σ_{e ∈ F_D} [ -∫_e (K ∇v_h · n) g_D ds                               │
│                            + σ ∫_K K r_e(g_D) · r_e(v_h) dx ]                       │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

| Term | File | Function/Class |
|------|------|----------------|
| Dirichlet RHS | `domain/antiplane_operator.hpp` | `AssembleDirichletRHS()` |
| Bottom BC value | `domain/antiplane_operator.hpp` | `BottomBCCoefficient` |

---

#### 4.4.2 Fault Slip Contribution

For prescribed slip `δ` on fault faces F_F:

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  L^Fault(v_h) = Σ_{e ∈ F_F} [ -∫_e {{K ∇v_h · n}} δ ds                             │
│                              + σ ∫_{K⁻ ∪ K⁺} K r_e(δ) · r_e([[v_h]]) dx ]          │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

**Note**: The slip `δ` enters the RHS (not the matrix) because it is time-dependent.

| Term | File | Function/Class |
|------|------|----------------|
| Fault RHS assembly | `domain/antiplane_operator.hpp` | `AddSlipContribution()` |
| Slip data | `fault/rate_state_fault.hpp` | `GetSlip()` |
| Fault face identification | `domain/antiplane_operator.hpp` | `IdentifyFaultFaces()` |

---

### 4.5 BR2 Summary: Complete Discrete System

**Find** u_h ∈ V_h **such that**:

```
┌─────────────────────────────────────────────────────────┐
│  a(u_h, v_h) = L(v_h)      ∀ v_h ∈ V_h                 │
└─────────────────────────────────────────────────────────┘
```

**Expanding all terms** (following Arnold et al. eq. 3.24, interior face terms sum
over F_I ∪ F_F):

```
BILINEAR FORM (Left-Hand Side):
===============================

[Volume]
  ∫_Ω K ∇u_h · ∇v_h dx

[Interior consistency — all interior faces including fault]
- Σ_{e ∈ F_I ∪ F_F} ∫_e {{K ∇u_h · n}} [[v_h]] ds

[Interior symmetry — all interior faces including fault]
- Σ_{e ∈ F_I ∪ F_F} ∫_e {{K ∇v_h · n}} [[u_h]] ds

[Interior BR2 lifting — all interior faces including fault]
+ Σ_{e ∈ F_I ∪ F_F} σ ∫ K r_e([[u_h]]) · r_e([[v_h]]) dx

[Dirichlet consistency]
- Σ_{e ∈ F_D} ∫_e (K ∇u_h · n) v_h ds

[Dirichlet symmetry]
- Σ_{e ∈ F_D} ∫_e (K ∇v_h · n) u_h ds

[Dirichlet BR2 lifting]
+ Σ_{e ∈ F_D} σ ∫ K r_e(u_h) · r_e(v_h) dx


LINEAR FORM (Right-Hand Side):
==============================

[Dirichlet RHS (symmetry)]
- Σ_{e ∈ F_D} ∫_e (K ∇v_h · n) g_D ds

[Dirichlet RHS (lifting)]
+ Σ_{e ∈ F_D} σ ∫ K r_e(g_D) · r_e(v_h) dx

[Fault RHS (symmetry)]
- Σ_{e ∈ F_F} ∫_e {{K ∇v_h · n}} δ ds

[Fault RHS (lifting)]
+ Σ_{e ∈ F_F} σ ∫ K r_e(δ) · r_e([[v_h]]) dx
```

**How fault enforcement works**: On fault faces e ∈ F_F, the bilinear form contributes
LHS terms with [[u_h]] and the RHS contributes terms with δ. The net effect is:

- Symmetry: -{{K∇v·n}}[[u_h]] (LHS) vs -{{K∇v·n}}δ (RHS) → enforces [[u_h]] → δ
- Lifting: σ r([[u_h]])·r([[v_h]]) (LHS) vs σ r(δ)·r([[v_h]]) (RHS) → penalizes [[u_h]] - δ
- Consistency: -{{K∇u·n}}[[v_h]] (LHS only) → provides flux coupling across fault

This is analogous to how Dirichlet BCs are handled: the bilinear form penalizes u toward 0,
and the RHS shifts the target to g_D. Here, the penalty enforces [[u_h]] toward 0, and the
fault RHS shifts it to δ.

---

## 5. IP Method: Complete Formulation

The **Interior Penalty (IP)** method uses direct penalty terms instead of lifting operators.

### 5.1 IP Penalty Parameter

For each face `e`, compute:

```
┌─────────────────────────────────────────┐
│  σ^IP_e = (1/4) (p_0 + p_1)            │
└─────────────────────────────────────────┘
```

where for each adjacent element K_i:

```
p_i = (D + 1) · c_N · (|e| / |K_i|) · (K_max² / K_min)
```

**Inverse inequality constant**:
```
c_N = (N + 1)(N + D) / D
```

where N = polynomial degree, D = spatial dimension.

**Example values** (2D, uniform K = 1):

| Order p | c_N | Notes |
|---------|-----|-------|
| 1 | 3 | Linear elements |
| 2 | 6 | Quadratic elements |
| 3 | 10 | Cubic elements |

For boundary faces: σ^IP_e = p_0 (single element).

**Implementation**: `domain/antiplane_operator.hpp` → `SetupIPPenalty()`

### 5.2 IP Bilinear Form

The complete IP bilinear form `a(u_h, v_h)` consists of:

```
┌─────────────────────────────────────────────────────────────────────────┐
│  a(u_h, v_h) = a^vol(u_h, v_h) + a^int(u_h, v_h) + a^Dir(u_h, v_h)     │
└─────────────────────────────────────────────────────────────────────────┘
```

where `a^int` sums over **all** interior faces (F_I ∪ F_F), including fault faces.

---

#### 5.2.1 Volume Term

```
┌─────────────────────────────────────────────────────────────┐
│  a^vol(u_h, v_h) = ∫_Ω K ∇u_h · ∇v_h dx                    │
└─────────────────────────────────────────────────────────────┘
```

**(Same as BR2)**

| Term | File | Function/Class |
|------|------|----------------|
| Volume integral | `domain/antiplane_operator.hpp` | `SetupBilinearForm()` |
| MFEM integrator | MFEM built-in | `DiffusionIntegrator` |

---

#### 5.2.2 Interior Face Terms (Including Fault Faces)

For each interior face `e ∈ F_I ∪ F_F` (all interior faces, including fault):

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  a^int(u_h, v_h) = Σ_{e ∈ F_I ∪ F_F} [ a^cons_e + a^sym_e + a^pen_e ]              │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

**Note**: Same as BR2 — fault faces receive the same bilinear form treatment as regular
interior faces (see Section 4.3.2 note).

**Consistency term** (same as BR2):
```
a^cons_e = -∫_e {{K ∇u_h · n}} [[v_h]] ds
```

**Symmetry term** (same as BR2):
```
a^sym_e = -∫_e {{K ∇v_h · n}} [[u_h]] ds
```

**IP penalty term** (replaces BR2 lifting):
```
a^pen_e = (σ_e / h_e) ∫_e [[u_h]] [[v_h]] ds
```

where h_e = face mesh size.

| Term | File | Function/Class |
|------|------|----------------|
| Interior face loop | `domain/antiplane_operator.hpp` | `AssembleInteriorFaces()` |
| All three terms | MFEM built-in | `DGDiffusionIntegrator` |
| Penalty parameter | `domain/antiplane_operator.hpp` | `SetupIPPenalty()` |

---

#### 5.2.3 Dirichlet Boundary Terms

For each Dirichlet face `e ∈ F_D`:

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│  a^Dir(u_h, v_h) = Σ_{e ∈ F_D} [ a^Dir,cons_e + a^Dir,sym_e + a^Dir,pen_e ]         │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

**Dirichlet consistency** (same as BR2):
```
a^Dir,cons_e = -∫_e (K ∇u_h · n) v_h ds
```

**Dirichlet symmetry** (same as BR2):
```
a^Dir,sym_e = -∫_e (K ∇v_h · n) u_h ds
```

**Dirichlet IP penalty**:
```
a^Dir,pen_e = (σ_e / h_e) ∫_e u_h v_h ds
```

| Term | File | Function/Class |
|------|------|----------------|
| Dirichlet BC | MFEM built-in | `DGDiffusionIntegrator` (boundary mode) |

---

### 5.3 IP Linear Form (Right-Hand Side)

The complete IP linear form `L(v_h)` consists of:

```
┌─────────────────────────────────────────────────────────┐
│  L(v_h) = L^Dir(v_h) + L^Fault(v_h)                    │
└─────────────────────────────────────────────────────────┘
```

---

#### 5.3.1 Dirichlet Contribution

For prescribed displacement `g_D` on Γ_D:

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  L^Dir(v_h) = Σ_{e ∈ F_D} [ -∫_e (K ∇v_h · n) g_D ds                               │
│                            + (σ_e / h_e) ∫_e g_D v_h ds ]                           │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

| Term | File | Function/Class |
|------|------|----------------|
| Dirichlet RHS | `domain/antiplane_operator.hpp` | `AssembleDirichletRHS()` |

---

#### 5.3.2 Fault Slip Contribution

For prescribed slip `δ` on fault faces F_F:

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│  L^Fault(v_h) = Σ_{e ∈ F_F} [ -∫_e {{K ∇v_h · n}} δ ds                             │
│                              + (σ_e / h_e) ∫_e δ [[v_h]] ds ]                       │
└─────────────────────────────────────────────────────────────────────────────────────┘
```

| Term | File | Function/Class |
|------|------|----------------|
| Fault RHS assembly | `domain/antiplane_operator.hpp` | `AddSlipContribution()` |

---

### 5.4 IP Summary: Complete Discrete System

**Find** u_h ∈ V_h **such that**:

```
┌─────────────────────────────────────────────────────────┐
│  a(u_h, v_h) = L(v_h)      ∀ v_h ∈ V_h                 │
└─────────────────────────────────────────────────────────┘
```

**Expanding all terms** (following Arnold et al. eq. 3.24, interior face terms sum
over F_I ∪ F_F):

```
BILINEAR FORM (Left-Hand Side):
===============================

[Volume]
  ∫_Ω K ∇u_h · ∇v_h dx

[Interior consistency — all interior faces including fault]
- Σ_{e ∈ F_I ∪ F_F} ∫_e {{K ∇u_h · n}} [[v_h]] ds

[Interior symmetry — all interior faces including fault]
- Σ_{e ∈ F_I ∪ F_F} ∫_e {{K ∇v_h · n}} [[u_h]] ds

[Interior IP penalty — all interior faces including fault]
+ Σ_{e ∈ F_I ∪ F_F} (σ_e / h_e) ∫_e [[u_h]] [[v_h]] ds

[Dirichlet consistency]
- Σ_{e ∈ F_D} ∫_e (K ∇u_h · n) v_h ds

[Dirichlet symmetry]
- Σ_{e ∈ F_D} ∫_e (K ∇v_h · n) u_h ds

[Dirichlet IP penalty]
+ Σ_{e ∈ F_D} (σ_e / h_e) ∫_e u_h v_h ds


LINEAR FORM (Right-Hand Side):
==============================

[Dirichlet RHS (symmetry)]
- Σ_{e ∈ F_D} ∫_e (K ∇v_h · n) g_D ds

[Dirichlet RHS (penalty)]
+ Σ_{e ∈ F_D} (σ_e / h_e) ∫_e g_D v_h ds

[Fault RHS (symmetry)]
- Σ_{e ∈ F_F} ∫_e {{K ∇v_h · n}} δ ds

[Fault RHS (penalty)]
+ Σ_{e ∈ F_F} (σ_e / h_e) ∫_e δ [[v_h]] ds
```

**How fault enforcement works**: Same mechanism as BR2 (see Section 4.5 explanation).
On fault faces, the IP penalty (σ/h)[[u_h]][[v_h]] on the LHS combined with
(σ/h)δ[[v_h]] on the RHS enforces [[u_h]] → δ.

---

## 6. Traction Computation

### 6.1 Physical Traction on Fault

The shear traction on the fault is:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  τ = μ · {{∂u/∂x}} = μ · (1/2) (∂u⁺/∂x + ∂u⁻/∂x)                           │
└─────────────────────────────────────────────────────────────────────────────┘
```

where μ = shear modulus (physical, not the normalized K=1).

### 6.2 Implementation

| Component | File | Function/Class |
|-----------|------|----------------|
| Traction computation | `domain/antiplane_operator.hpp` | `ComputeTraction()` |
| Gradient evaluation | MFEM | `GradientGridFunctionCoefficient` |
| Fault face loop | `domain/antiplane_operator.hpp` | `fault_interior_faces_` |
| Traction output | `fault/rate_state_fault.hpp` | `GetShearStress()` |

### 6.3 State Variable Notation

The rate-and-state friction literature uses two equivalent representations for the state variable:

**ψ (psi) form - Natural logarithm representation:**
```
f(V, ψ) = a · asinh[(V / 2V₀) · exp(ψ/a)]
dψ/dt = (b·V₀/Dc) · [exp((f₀ - ψ)/b) - V/V₀]
```

**θ (theta) form - Physical age representation:**
```
f(V, θ) = a · asinh[(V / 2V₀) · exp((f₀ + b·ln(V₀θ/Dc)) / a)]
dθ/dt = 1 - V·θ/Dc
```

**Equivalence:** The two forms are related by:
```
ψ = f₀ + b·ln(V₀·θ/Dc)
θ = (Dc/V₀) · exp((ψ - f₀)/b)
```

The SEAS miniapp implementation uses the **θ form** because:
1. θ has physical meaning (contact age) with units of [seconds]
2. The evolution law dθ/dt = 1 - V·θ/Dc is simpler
3. Steady state θ_ss = Dc/V is intuitive

This matches Tandem's internal implementation despite literature often using ψ notation.

---

# Part II: Implementation Mapping

---

## 7. File-to-Formula Mapping

### 7.1 Source File Overview

```
miniapps/seas/
├── config/
│   └── bp2_params.hpp              # BP2 benchmark parameters
├── common/
│   ├── seas_types.hpp              # Shared type definitions
│   └── mpi_context.hpp             # MPI context utilities
├── domain/
│   ├── domain_operator.hpp         # Abstract DomainOperator interface
│   ├── antiplane_operator.hpp      # AntiplaneDomainOperator class + DGMethod enum
│   └── bp2_mesh.hpp                # BP2 mesh generation + boundary attributes
├── integrator/
│   └── dg_br2_integrator.hpp       # BR2 face integrators (interior + boundary)
├── fault/
│   ├── fault_geometry.hpp          # FaultGeometry class
│   └── rate_state_fault.hpp        # RateStateFaultOperator class
├── friction/
│   ├── friction_law.hpp            # Abstract FrictionLaw interface
│   ├── dieterich_ruina.hpp         # DieterichRuinaFriction class
│   └── state_evolution.hpp         # AgingLaw/SlipLaw classes
└── tests/unit/
    └── test_antiplane.cpp          # 24 unit tests for domain operator
```

### 7.2 BR2 Method: Term-to-Code Mapping

| Formula Term | File | Class/Function | Status |
|--------------|------|----------------|--------|
| **Bilinear Form** | | | |
| Volume: ∫ K ∇u · ∇v | `antiplane_operator.hpp` | `Solve()` | ✓ Done |
| Interior consistency: -∫ {{K ∇u · n}} [[v]] | `dg_br2_integrator.hpp` | `BR2InteriorFaceIntegrator` | ✓ Done |
| Interior symmetry: -∫ {{K ∇v · n}} [[u]] | `dg_br2_integrator.hpp` | `BR2InteriorFaceIntegrator` | ✓ Done |
| Interior lifting: σ ∫ r([[u]]) · r([[v]]) | `dg_br2_integrator.hpp` | `BR2InteriorFaceIntegrator` | ✓ Done |
| Dirichlet terms | `dg_br2_integrator.hpp` | `BR2BoundaryFaceIntegrator` | ✓ Done |
| **Linear Form** | | | |
| Dirichlet RHS | `antiplane_operator.hpp` | `DGDirichletLFIntegrator` | ✓ Done |
| Fault slip RHS | `antiplane_operator.hpp` | `AssembleSlipContributionBR2()` | ✓ Done |
| **Supporting** | | | |
| Mass matrix inverse | `antiplane_operator.hpp` | `PrecomputeMassMatrixInverses()` | ✓ Done |
| Lifting operator | `dg_br2_integrator.hpp` | Inline in `AssembleFaceMatrix()` | ✓ Done |
| Penalty (σ = D+1) | `antiplane_operator.hpp` | `br2_penalty_ = 3.0` | ✓ Done |
| Fault face ID | `antiplane_operator.hpp` | `SetupFaultInfo()` | ✓ Done |

### 7.3 IP Method: Term-to-Code Mapping

| Formula Term | File | Class/Function | Status |
|--------------|------|----------------|--------|
| **Bilinear Form** | | | |
| Volume: ∫ K ∇u · ∇v | `antiplane_operator.hpp` | `Solve()` | ✓ Done |
| Interior consistency | MFEM | `DGDiffusionIntegrator` | ✓ Built-in |
| Interior symmetry | MFEM | `DGDiffusionIntegrator` | ✓ Built-in |
| Interior penalty: (σ/h) ∫ [[u]] [[v]] | MFEM | `DGDiffusionIntegrator` | ✓ Built-in |
| Dirichlet terms | MFEM | `DGDiffusionIntegrator` | ✓ Built-in |
| **Linear Form** | | | |
| Dirichlet RHS | `antiplane_operator.hpp` | `DGDirichletLFIntegrator` | ✓ Done |
| Fault slip RHS | `antiplane_operator.hpp` | `AssembleSlipContributionIP()` | ✓ Done |
| **Supporting** | | | |
| IP penalty κ = (p+1)² | MFEM | `DGDiffusionIntegrator` (h-scaling internal) | ✓ Built-in |

### 7.4 Cross-Cutting Components

| Component | File | Description |
|-----------|------|-------------|
| DG method selection | `antiplane_operator.hpp` | `enum class DGMethod { BR2, IP }` |
| DG FE space | `antiplane_operator.hpp` | `DG_FECollection`, `L2_FECollection` |
| Solver setup | `antiplane_operator.hpp` | `SetupSolver()` |
| Parallel support | `antiplane_operator.hpp` | Template `<MeshType>` |

---

# Appendix: Tandem Reference

---

## Appendix A: Tandem Implementation Details

### A.1 BR2 Penalty (from Tandem)

**From `Poisson.h` lines 112-117**:
```cpp
double penalty(std::size_t fctNo) const {
    if (method_ == DGMethod::BR2) {
        return NumFacets;  // = Dim + 1
    }
    return penalty_[fctNo];
}
```

### A.2 IP Penalty (from Tandem)

**From `Poisson.cpp` lines 209-223**:
```cpp
auto const p = [&](int side) {
    auto Kfield = material[info.up[side]].get<K>().data();
    auto k0 = *std::min_element(Kfield, Kfield + numBasisFunctions);
    auto k1 = *std::max_element(Kfield, Kfield + numBasisFunctions);
    constexpr double c_N = InverseInequality<Dim>::trace_constant(PolynomialDegree - 1);
    return (Dim + 1) * c_N * (area_[fctNo] / volume_[info.up[side]]) * (k1 * k1 / k0);
};

if (info.up[0] != info.up[1]) {
    penalty_[fctNo] = (p(0) + p(1)) / 4.0;  // Interior face
} else {
    penalty_[fctNo] = p(0);  // Boundary face
}
```

### A.3 Inverse Inequality Constant (from Tandem)

**From `InverseInequality.h` lines 27-29**:
```cpp
constexpr static double trace_constant(unsigned N) {
    return (N + 1) * (N + D) / static_cast<double>(D);
}
```

### A.4 Lifting Operator (from Tandem)

**From `poisson.py` lines 59-65**:
```python
# Lifting operator on each side
Lift[0]['liu'] <= 0.5 * Minv[0]['us'] * E_q[0]['sq'] * E_q[x]['lq'] * n_q['iq'] * w['q']
Lift[1]['liv'] <= 0.5 * Minv[1]['vs'] * E_q[1]['sq'] * E_q[x]['lq'] * n_q['iq'] * w['q']

# Combined lifted flux
L_q[x]['lq'] <= 0.5 * n_q['iq'] * (K_q[0]['q'] * E_q[0]['uq'] * Lift[0]['liu'] +
                                    K_q[1]['q'] * E_q[1]['vq'] * Lift[1]['liv'])
```

### A.5 Assembly Coefficients (from Tandem)

**From `Poisson.cpp` lines 305-310**:
```cpp
assemble.c00 = -0.5;              // flux from side 0 on test from side 0
assemble.c01 = +0.5;              // flux from side 1 on test from side 0
assemble.c10 = epsilon * 0.5;     // = -0.5 (symmetric)
assemble.c11 = +0.5;              // opposite sign
assemble.c20 = penalty(fctNo);    // = Dim + 1 for BR2
assemble.c21 = -penalty(fctNo);   // opposite sign for cross-element
```

### A.6 Boundary Condition Enumeration (from Tandem)

**From `BC.h` line 6**:
```cpp
enum class BC : int { None = 0, Natural = 1, Fault = 3, Dirichlet = 5 };
```

### A.7 Fault RHS Assembly (from Tandem)

**From `poisson.py` lines 102-104**:
```python
b['k'] <= b['k'] + c1[0] * w['q'] * K_Dx_q[0]['kiq'] * n_q['iq'] * f_q['q'] +
                   c2[0] * w['q'] * E_q[0]['kq'] * f_lifted_q['q']
```

---

## Appendix B: Tandem Source File Reference

### B.1 Key Source Files

```
tandem/
├── app/
│   ├── localoperator/
│   │   ├── Poisson.h           # Poisson operator header
│   │   ├── Poisson.cpp         # Full implementation
│   │   └── DieterichRuinaAgeing.h  # Friction law
│   ├── kernels/
│   │   └── poisson.py          # Tensor kernel definitions
│   └── form/
│       ├── BC.h                # Boundary condition enum
│       └── FacetInfo.h         # Face information structure
└── src/
    └── form/
        ├── DGCurvilinearCommon.h   # DG method enum
        ├── DGOperator.h            # Base operator class
        └── InverseInequality.h     # Penalty constants
```

### B.2 Critical Implementation Locations

| Component | File | Lines |
|-----------|------|-------|
| BR2 penalty | `Poisson.h` | 112-117 |
| IP penalty | `Poisson.cpp` | 209-223 |
| Volume kernel | `poisson.py` | 30-32 |
| Surface assembly | `poisson.py` | 71-75 |
| Lifting operator | `poisson.py` | 59-65 |
| Skeleton BC | `Poisson.cpp` | 406-417 |
| Boundary BC | `Poisson.cpp` | 418-432 |
| Fault RHS | `poisson.py` | 102-104 |

---

## References

1. **Arnold, D.N., et al.** (2002). "Unified Analysis of Discontinuous Galerkin Methods for Elliptic Problems." *SIAM J. Numer. Anal.* 39(5), 1749-1779.

2. **Bassi, F. and Rebay, S.** (1997). "A high-order accurate discontinuous finite element method..." *J. Comput. Phys.* 131, 267-279.

3. **Uphoff, C., May, D.A., Gabriel, A.-A.** (2023). "A discontinuous Galerkin method for sequences of earthquakes and aseismic slip..."

4. **Tandem Code**: https://github.com/TEAR-ERC/tandem

5. **SCEC SEAS Benchmarks**: https://strike.scec.org/cvws/seas/

6. **MFEM Examples**: `ex14.cpp` (DG diffusion), `ex17.cpp` (DG elasticity)
