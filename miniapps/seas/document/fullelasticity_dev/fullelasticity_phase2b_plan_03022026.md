# Phase 2b: DG Elasticity BR2 Integrator + ElasticityDomainOperator

**Date**: 2026-03-02
**Status**: COMPLETE (14 BR2 tests + 32 elasticity operator tests passing)
**Prerequisite**: Phase 1 complete (BP5Params, vector friction, FaultBasis, DomainOperator interface)

---

## Overview

Phase 2b implements the two core components for 3D DG elasticity:

1. **`DGElasticityBR2Integrator`** — Vector elasticity generalization of the scalar BR2 face integrator
2. **`ElasticityDomainOperator`** — Full 3D DG elasticity solver with fault slip BC

Both support BR2 (default, matching Tandem) and IP (alternative, using MFEM built-in) DG methods.

---

## Part 1: DGElasticityBR2Integrator

**File**: `integrator/dg_elasticity_br2_integrator.hpp` (NEW, ~850 lines)

### 1.1 Mathematical Formulation

For each interior face, the DG bilinear form has three terms:

```
a(u,v) = -∫_F {{σ(u)·n}} · [[v]] ds              (consistency)
        - ε∫_F {{σ(v)·n}} · [[u]] ds              (symmetry, ε=-1 for SIPG)
        + σ_BR2 ∫_F C:R_h(u) : R_h(v) ds          (BR2 lifting penalty)
```

where:
- `σ(u) = λ·tr(ε(u))·I + 2μ·ε(u)` is the stress tensor
- `R_h([[u]])` is the lifted jump (tensor-valued)
- `σ_BR2 = num_faces_per_element` (6 for hex, 4 for tet)
- `C` is the elasticity tensor

### 1.2 The `test_normal` Operator

Following Tandem's `elasticity.py` (lines 118-120), the elasticity tensor coupling in the BR2 lifting is:

```
test_normal(x)_{iu,sq} = λ(x)·δ_us·n_i + μ(x)·(δ_iu·n_s + δ_is·n_u)
```

This is `C_ijrs·n_s` contracted with the face normal — the traction operator applied to the lifting. For scalar Laplace, this reduces to `K·n_i` (just diffusion coefficient times normal), recovering the scalar BR2 integrator.

### 1.3 Lifting Computation

Generalizes the scalar BR2 pattern from `integrator/dg_br2_integrator.hpp`:

**Step 1: Scalar lifting** (same structure as existing scalar BR2):
```
Lift[elem][l, s, m] = 0.5 · Minv[elem][m,o] · Σ_q E[elem][o,q] · E[source][l,q] · n[s,q] · w[q]
```

**Step 2: Apply elasticity tensor coupling** (new for vector case):
```
L_q[source][l, i, u, q] = 0.5 · Σ_elem test_normal(elem)_{iu,sq} · E[elem][m,q] · Lift[elem][l, s, m]
```

This produces a `(ndof_source × dim × dim × nqp)` tensor. Indices: `l`=source DOF, `i`=test component, `u`=trial component, `q`=quadrature point.

### 1.4 Face Matrix Assembly

The face matrix has block structure `(dim·ndof1 + dim·ndof2) × (dim·ndof1 + dim·ndof2)`. For each block pair `(x_elem, y_elem)`:

```
a[(x,i), (y,u)] = c0[y] · ∫ σ_i(φ_x)·n · φ_y^u ds     (consistency)
                 + c1[x] · ∫ σ_u(φ_y)·n · φ_x^i ds     (symmetry)
                 + c2[|x-y|] · ∫ φ_x^i · L_q[y][l,i,u,q] ds  (BR2 penalty)
```

Sign conventions: `c0=-0.5`, `c1=ε·0.5`, `c2[same]=+σ_BR2`, `c2[cross]=-σ_BR2`.

### 1.5 Class Interface

```cpp
class DGElasticityBR2Integrator : public BilinearFormIntegrator
{
public:
   DGElasticityBR2Integrator(Coefficient &lambda, Coefficient &mu,
                              real_t epsilon,
                              const std::vector<DenseMatrix> &elem_mass_inv,
                              int dim = 3);

   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override;
protected:
   Coefficient &lambda_, &mu_;
   real_t epsilon_, sigma_;
   const std::vector<DenseMatrix> &elem_mass_inv_;
   int dim_;
};

class DGElasticityBR2BoundaryIntegrator : public BilinearFormIntegrator
{
   // Same interface, for Dirichlet boundary faces
};
```

### 1.6 Key Implementation Details

- **Mass matrix inverse**: Per-component (scalar), not the full vector mass matrix. The vector DG space is `vdim` copies of the same scalar DG space, so the scalar mass matrix `M_K` is identical for all components. Stored per-element as `DenseMatrix` in `elem_mass_inv_`.
- **BR2 penalty**: `σ = 2·dim` for hexahedra (6 faces), `dim+1` for tetrahedra (4 faces).
- **SIPG sign**: `ε = -1` produces symmetric face matrices.
- **Symmetry check**: Uses Frobenius-norm-relative check (robust for large material values like BP5's μ ≈ 32 GPa).

---

## Part 2: ElasticityDomainOperator

**File**: `domain/elasticity_operator.hpp` (NEW, ~900 lines)

### 2.1 Class Structure

```cpp
template <typename MeshType = Mesh>
class ElasticityDomainOperator : public DomainOperator<MeshType>
{
public:
   ElasticityDomainOperator(MeshType &mesh, int order,
                             real_t lambda, real_t mu,
                             real_t Vp, real_t Wf, real_t lf,
                             DGMethod method = DGMethod::BR2);

   int NumComponents() const override { return 3; }
   int Dimension() const override { return 3; }
   int NumSlipComponents() const override { return 2; }

   void Solve(real_t time, const Vector &slip_bc, GridFuncType &u) override;
   void ComputeTraction(const GridFuncType &u, const Vector &slip_bc,
                        Vector &traction) override;

   const FaultBasis *GetFaultBasis() const override;
   void GetFaultCoords2D(Vector &x2, Vector &x3) const override;
   int GetNumFaultDOFs() const override;
   void GetFaultDepths(Vector &depths) const override;
   real_t GetShearModulus() const override;
};
```

### 2.2 FE Space Setup

Two DG FE spaces:
- **Vector space** (`vdim=3`, `Ordering::byNODES`): For the actual elasticity problem
- **Scalar space** (`vdim=1`): For computing BR2 mass matrix inverses

```cpp
fec_ = std::make_unique<DG_FECollection>(order_, 3, BasisType::GaussLobatto);
fes_ = std::make_unique<FESpaceType>(&mesh_, fec_.get(), 3, Ordering::byNODES);
scalar_fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
```

### 2.3 Stiffness Matrix Assembly

Both BR2 (default) and IP are supported:

```cpp
void AssembleStiffness() {
   // Volume term (same for both methods):
   cached_a_->AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff_, mu_coeff_));

   if (method_ == DGMethod::BR2) {
      PrecomputeMassInverse();
      cached_a_->AddInteriorFaceIntegrator(
         new DGElasticityBR2Integrator(...));
      cached_a_->AddBdrFaceIntegrator(
         new DGElasticityBR2BoundaryIntegrator(...), dirichlet_bdr_marker_);
   } else {
      real_t kappa = (order_ + 1) * (order_ + 1);
      cached_a_->AddInteriorFaceIntegrator(
         new DGElasticityIntegrator(lambda_coeff_, mu_coeff_, -1.0, kappa));
      cached_a_->AddBdrFaceIntegrator(
         new DGElasticityIntegrator(...), dirichlet_bdr_marker_);
   }
   cached_a_->Assemble();
   cached_a_->Finalize();
}
```

### 2.4 Fault Face Detection (3D)

Interior faces at `x1=0` within the fault zone (`|x2| <= lf/2`, `0 <= x3 <= Wf`):

```cpp
bool IsFaultFace3D(const Vector &center) const {
   return std::abs(center(0)) < tol
       && std::abs(center(1)) <= lf_/2.0 + tol
       && center(2) >= -tol
       && center(2) <= Wf_ + tol;
}
```

### 2.5 Slip BC Assembly

Slip is received in fault-local frame `(s_dip, s_strike)` and embedded into global frame via `FaultBasis::EmbedSlip`. Both BR2 and IP slip RHS assembly are implemented, selected by `method_`.

### 2.6 Traction Computation

After solving for displacement `u`, extract stress on fault faces and project to local frame:

1. Evaluate `∇u` on both sides of each fault face, average
2. Compute stress tensor: `σ_ij = λ·tr(ε)·δ_ij + 2μ·ε_ij`
3. Compute global traction: `T = σ·n`
4. Project to local frame via `FaultBasis::ProjectTraction` → `(τ_dip, τ_strike)`

### 2.7 Boundary Conditions

| Face | Attribute | BC Type | Treatment |
|------|-----------|---------|-----------|
| x₁=-Lx | 1 | Natural (zero traction) | No integrator |
| x₁=+Lx | 2 | Natural (zero traction) | No integrator |
| x₂=+Ly (strike+) | 3 | Dirichlet: u₂=+Vp·t/2 | Boundary face integrator |
| x₂=-Ly (strike-) | 4 | Dirichlet: u₂=-Vp·t/2 | Boundary face integrator |
| x₃=0 (free surface) | 5 | Natural (zero traction) | No integrator |
| x₃=Lz (bottom) | 6 | Natural (zero traction) | No integrator |

### 2.8 MFEM `mutable` Pattern

MFEM's `Coefficient::Eval()` and `AddBdrFaceIntegrator()` take non-const references, but assembly is logically const. Solution: `lambda_coeff_`, `mu_coeff_`, and `dirichlet_bdr_marker_` are declared `mutable`.

---

## Tests

### BR2 Integrator Tests (`tests/unit/test_elasticity_br2.cpp`) — 14 tests

1. **Face matrix symmetry** (SIPG ε=-1): Verified symmetric for interior faces
2. **Face matrix non-zero**: Verified norm > 0
3. **Boundary face symmetry**: Verified for boundary faces
4. **Full bilinear form assembly**: Complete stiffness matrix assembled, symmetric, correct size
5. **Order 0**: Simplest case (6×6 face matrix), verified symmetric
6. **BP5 material params**: λ=μ≈32 GPa, verified symmetry (Frobenius-norm-relative check)

### Elasticity Operator Tests (`tests/unit/test_elasticity_operator.cpp`) — 32 tests

1. **Construction**: NumComponents=3, Dimension=3, NumSlipComponents=2, FaultBasis non-null
2. **Fault detection at x1=0**: Interior fault faces found, depths and 2D coords correct
3. **FaultBasis properties**: Normal ≈ (1,0,0), orthonormal tangent vectors, dot products ≈ 0
4. **Zero slip equilibrium**: Zero slip + zero loading → zero displacement
5. **Traction extraction**: Zero slip → near-zero traction
6. **Stiffness assembly**: IP and BR2 both assemble successfully
7. **FaultGeometry 3D**: BP5Params constructor, DOF count matches, a/eta/Dc/tau_pre sizes correct

---

## Key Design Decisions

1. **BR2 as default method**: Matches Tandem's default. MFEM has no built-in BR2 for elasticity, so we implement `DGElasticityBR2Integrator` following the existing scalar BR2 pattern with elasticity tensor coupling.

2. **Scalar mass matrix inverse for BR2**: The vector DG space is `vdim` copies of a scalar DG space, so the component-wise mass matrix is identical. Stored per-element.

3. **Both slip assembly methods**: `AssembleSlipContributionBR2` and `AssembleSlipContributionIP` are both implemented, selected by `method_`. Mirrors the antiplane operator exactly.

4. **Traction in local frame**: `ComputeTraction` returns 2-component local-frame traction `(τ_dip, τ_strike)`. The domain operator is the only place touching both coordinate frames.

---

## Files Created/Modified

| File | Action | Lines |
|------|--------|-------|
| `integrator/dg_elasticity_br2_integrator.hpp` | NEW | ~850 |
| `domain/elasticity_operator.hpp` | NEW | ~900 |
| `tests/unit/test_elasticity_br2.cpp` | NEW | ~406 |
| `tests/unit/test_elasticity_operator.cpp` | NEW | ~375 |
| `Makefile` | EXTENDED | New targets + header groups |

## Phase 2a Carry-Forward Items

The following items were identified during Phase 2a review
(`fullelasticity_phase2a_check_03032026.md`) as deferred to Phase 2b:

### CF-1: Rotated (non-axis-aligned) fault FaultBasis test

All existing FaultBasis tests use axis-aligned faults (normal along x-axis).
A test with a non-axis-aligned 3D fault would exercise the full generality
of the strike/dip decomposition. Now that `ElasticityDomainOperator` provides
realistic 3D mesh geometry, this test can be written using the operator's
fault face detection on a rotated mesh or by constructing a simple mesh with
a tilted interior face.

### CF-2: DG two-sided sign handling for fault faces

`FaultBasis` computes one consistent basis per face (oriented by `ref_normal`).
Tandem additionally flips the entire basis matrix for the "minus" side via a
`SignFlipped` flag in `AdapterBase::prepare()`. SEAS-MFEM handles the plus/minus
distinction in `ElasticityDomainOperator` itself (sign applied when computing
jumps/averages in the DG numerical flux), not in `FaultBasis`. Verify that the
sign handling in `ComputeTraction` and `AssembleSlipContribution*` is correct
for both sides of each fault face.

### CF-3: Mass matrix and higher-order DG on fault

For higher-order DG elements (order > 0), traction projection from quadrature
points to fault DOFs requires a fault mass matrix and its inverse, as Tandem
does. The current direct dot-product approach in `FaultBasis::ProjectTraction`
works for DG0 (one DOF per face). If higher-order elements are needed on the
fault for BP5, the traction computation would need to integrate over quadrature
points with a fault mass matrix. Assess whether DG0 is sufficient or if
higher-order is needed.

### CF-4: Per-face vs per-quadrature-point basis granularity

`FaultBasis` stores one basis per face (constant normal per planar face).
For higher-order elements on curved faults, per-DOF or per-quadrature-point
basis would be needed. For BP5 with a planar fault, per-face is sufficient.
Only relevant if curved faults are considered in the future.

---

## Tandem Reference

| SEAS-MFEM | Tandem Reference |
|-----------|-----------------|
| `DGElasticityBR2Integrator` | `app/localoperator/Elasticity.cpp` + `app/kernels/elasticity.py` |
| `ElasticityDomainOperator` | `app/localoperator/Elasticity.h` |
| `test_normal` operator | `elasticity.py` lines 118-120 |
| BR2 lifting | `lift_skeleton` kernel in `Elasticity.cpp` |
