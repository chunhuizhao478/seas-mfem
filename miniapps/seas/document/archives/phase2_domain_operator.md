# Phase 2: Domain Operator (Serial)

## Overview

This phase implements the domain operator for solving the antiplane shear problem (Laplace equation) with appropriate boundary conditions. **Following Tandem's approach**, we implement a DG solver with BR2 (default) or IP method support.

## Reference Implementation

**Tandem source files** (in `/Users/chunhuizhao/projects/tandem/`):
- `app/localoperator/Poisson.h` - Poisson operator header
- `app/localoperator/Poisson.cpp` - Full DG implementation with BR2/IP
- `app/kernels/poisson.py` - Tensor kernel definitions
- `src/form/BC.h` - Boundary condition enumeration

**Theory document**: See `dg_antiplane_theory.md` for complete mathematical derivation:
- **Section 1**: Strong form (Laplace equation)
- **Section 2**: Domain and boundary conditions
- **Section 3**: DG notation (jump/average operators)
- **Section 4**: BR2 method complete formulation
- **Section 5**: IP method complete formulation
- **Section 6**: Traction computation

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 2.1 Abstract domain interface | `domain_operator.hpp` | - |
| 2.2 DG Antiplane implementation | `antiplane_operator.hpp` | `test_antiplane.cpp` |
| 2.3 BR2/IP integrators | `dg_br2_integrator.hpp` | In `test_antiplane.cpp` |
| 2.4 Traction computation | In `antiplane_operator.hpp` | In `test_antiplane.cpp` |
| 2.5 MMS verification | In `antiplane_operator.hpp` (`SolveMMS`) | In `test_antiplane.cpp` |
| 2.6 BP2 mesh construction | `bp2_mesh.hpp` | In `test_antiplane.cpp` |
| 2.7 Fault face identification | In `antiplane_operator.hpp` | In `test_antiplane.cpp` |

## Verification Checkpoint

MMS test shows expected convergence order (p+1 in L² norm for polynomial order p).

---

## Physical Background

### BP2 Governing Equations

**Equilibrium equation** (quasi-static):
```
0 = ∂σ_xy/∂x + ∂σ_yz/∂z
```

**Hooke's law** (antiplane shear):
```
σ_xy = μ ∂u/∂x
σ_yz = μ ∂u/∂z
```

This reduces to the **Laplace equation**:
```
∇²u = ∂²u/∂x² + ∂²u/∂z² = 0
```

### Computational Domain (Full Domain)

```
Computational Domain: x ∈ [-Lx, +Lx], z ∈ [-Lz, 0]

                    z = 0 (free surface)
    <--------------------------------------------------->
    |                    |                              |
    |   Left half        |    Right half                |
    |   (x < 0)          |    (x > 0)                   |
    |                    |                              |
x=-Lx      u → -Vp·t/2   |0   u → +Vp·t/2            x=+Lx
    |                    |                              |
    |                  FAULT                            |
    |                (slip = δ)                         |
    |                    |                              |
    <--------------------------------------------------->
                    z = -Lz (bottom)
```

### Boundary Conditions (Full Domain)

| Boundary | Location | Condition | Mathematical Form |
|----------|----------|-----------|-------------------|
| Fault | x = 0, z > -Wf | Interior Jump | [[u]] = δ |
| Below fault | x = 0, z ≤ -Wf | Natural | K∇u · n = 0 (free slip) |
| Far-field left | x = -Lx | Natural | K∇u · n = 0 |
| Far-field right | x = +Lx | Natural | K∇u · n = 0 |
| Free surface | z = 0 | Natural | K∇u · n = 0 |
| Bottom | z = -Lz | Dirichlet | u = sign(x)·Vp·t/2 |

**Note**: Following Tandem's approach, plate loading is applied at the **bottom boundary** (z = -Lz), NOT at far-field boundaries.

---

## DG Method Selection (Following Tandem)

### Tandem's DG Methods

From `src/form/DGCurvilinearCommon.h`:
```cpp
enum class DGMethod { BR2, IP };
```

**Tandem default: BR2**

### BR2 Method (Bassi-Rebay 2)

BR2 uses lifting operators instead of direct penalty terms:

```
a(u,v) = ∫_Ω K∇u·∇v dx                           (volume)
       - ∫_{F_I} {{K∇u·n}} [[v]] ds               (consistency)
       - ∫_{F_I} {{K∇v·n}} [[u]] ds               (symmetry, ε=-1)
       + ∫_{F_I} σ r_e([[u]])·r_e([[v]]) dx       (BR2 lifting)
```

**BR2 Penalty**: σ = Dim + 1 (3 for 2D triangles, 4 for 3D tetrahedra)

**Lifting operator** (from Tandem `poisson.py` lines 59-65):
```python
Lift[i] = 0.5 * M_inv[i] * E_q[i]^T * n_q * w
L_q = 0.5 * n · (K_0 * E_0 * Lift[0] + K_1 * E_1 * Lift[1])
```

### IP Method (Interior Penalty)

The alternative SIPG method with direct penalty:

```
a(u,v) = ∫_Ω K∇u·∇v dx
       - ∫_{F_I} {{K∇u·n}} [[v]] ds
       - ∫_{F_I} {{K∇v·n}} [[u]] ds
       + σ/h ∫_{F_I} [[u]][[v]] ds
```

**IP Penalty** (from Tandem `Poisson.cpp` lines 209-223):
```
σ = (1/4)(p_0 + p_1)

p_i = (D+1) · c_N · (|e|/|E_i|) · (K_max²/K_min)

c_N = (N+1)(N+D)/D   (inverse inequality constant)
```

---

## Detailed Component Design

### 2.1 Abstract Domain Operator Interface

```cpp
// domain/domain_operator.hpp

#include "common/seas_types.hpp"

namespace mfem {
namespace seas {

/// Abstract base class for domain operators
/// Template parameter controls serial vs parallel MFEM types
/// Type aliases use FESpaceForMesh<>, GridFunctionForMesh<>, etc.
/// from common/seas_types.hpp for compile-time serial/parallel switching.
template <typename MeshType = Mesh>
class DomainOperator {
public:
    // Type aliases derived from mesh type via seas_types.hpp
    using FESpaceType = FESpaceForMesh<MeshType>;
    using GridFuncType = GridFunctionForMesh<MeshType>;
    using BilinFormType = BilinearFormForMesh<MeshType>;
    using LinFormType = LinearFormForMesh<MeshType>;

    virtual ~DomainOperator() = default;

    /// Get number of displacement components (1 for antiplane)
    virtual int NumComponents() const = 0;

    /// Get spatial dimension (2 for 2D)
    virtual int Dimension() const = 0;

    /// Solve domain problem with given fault slip BC
    virtual void Solve(real_t time, const Vector &slip_bc,
                       GridFuncType &displacement) = 0;

    /// Compute traction at fault from displacement
    /// Takes slip_bc as input for penalty/consistency terms
    virtual void ComputeTraction(const GridFuncType &displacement,
                                  const Vector &slip_bc,
                                  Vector &traction) = 0;

    /// Get reference to finite element space
    virtual FESpaceType &GetFESpace() = 0;
    virtual const FESpaceType &GetFESpace() const = 0;

    /// Get the underlying mesh
    virtual MeshType &GetMesh() = 0;
    virtual const MeshType &GetMesh() const = 0;

    /// Get shear modulus
    virtual real_t GetShearModulus() const = 0;

    /// Get number of fault DOFs
    virtual int GetNumFaultDOFs() const = 0;

    /// Get depth at each fault DOF
    virtual void GetFaultDepths(Vector &depths) const = 0;

    /// Get fault DOF indices
    virtual const Array<int> &GetFaultDOFs() const = 0;

#ifdef SEAS_USE_MPI
    /// Get MPI communicator (parallel only)
    virtual MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#endif
};

// Convenience type aliases
using SerialDomainOperator = DomainOperator<Mesh>;
#ifdef MFEM_USE_MPI
using ParallelDomainOperator = DomainOperator<ParMesh>;
#endif

} // namespace seas
} // namespace mfem
```

### 2.2 DG Method Enumeration (Following Tandem)

```cpp
// domain/dg_method.hpp

namespace mfem {
namespace seas {

/// DG method selection (following Tandem)
enum class DGMethod {
    BR2,    ///< Bassi-Rebay 2 (default, Tandem default)
    IP      ///< Interior Penalty (SIPG with ε = -1)
};

} // namespace seas
} // namespace mfem
```

### 2.3 Antiplane Operator Implementation

```cpp
// domain/antiplane_operator.hpp

namespace mfem {
namespace seas {

/// Solves ∇²u = 0 with slip BC on fault
/// Following Tandem's approach with BR2/IP method selection
template <typename MeshType = Mesh>
class AntiplaneDomainOperator : public DomainOperator<MeshType> {
public:
    using Base = DomainOperator<MeshType>;
    using FESpaceType = typename Base::FESpaceType;
    using GridFuncType = typename Base::GridFuncType;
    using BilinFormType = typename std::conditional<
        std::is_same<MeshType, ParMesh>::value,
        ParBilinearForm, BilinearForm>::type;
    using LinFormType = typename std::conditional<
        std::is_same<MeshType, ParMesh>::value,
        ParLinearForm, LinearForm>::type;

    /// Constructor with DG method selection
    /// @param mesh The computational mesh
    /// @param order Polynomial order
    /// @param mu Shear modulus [Pa]
    /// @param Vp Plate rate [m/s] for far-field boundary condition
    /// @param Wf Fault depth [m] - slip applied from z=0 to z=-Wf
    /// @param method DG method (IP default)
    AntiplaneDomainOperator(MeshType &mesh, int order, real_t mu,
                            real_t Vp, real_t Wf = 40.0e3,
                            DGMethod method = DGMethod::IP);

    int NumComponents() const override { return 1; }
    int Dimension() const override { return mesh_.Dimension(); }

    void Solve(real_t time, const Vector &slip_bc,
               GridFuncType &displacement) override;

    /// Solve with custom Dirichlet BC for MMS verification
    void SolveMMS(Coefficient &dirichlet_bc, GridFuncType &displacement);

    /// Compute traction including slip_bc for penalty/consistency terms
    void ComputeTraction(const GridFuncType &displacement,
                          const Vector &slip_bc,
                          Vector &traction) override;

    FESpaceType &GetFESpace() override { return *fes_; }
    const FESpaceType &GetFESpace() const override { return *fes_; }
    MeshType &GetMesh() override { return mesh_; }
    const MeshType &GetMesh() const override { return mesh_; }

    real_t GetShearModulus() const override { return mu_; }
    int GetNumFaultDOFs() const override { return num_fault_dofs_; }
    void GetFaultDepths(Vector &depths) const override;
    const Array<int> &GetFaultDOFs() const override { return fault_dofs_; }

    /// Get fault interior face indices
    const Array<int> &GetFaultInteriorFaces() const;

    /// Get plate rate
    real_t GetPlateRate() const { return Vp_; }

    /// Get fault depth
    real_t GetFaultDepth() const { return Wf_; }

    /// Get polynomial order
    int GetOrder() const { return order_; }

    /// Get DG method
    DGMethod GetMethod() const { return method_; }

private:
    MeshType &mesh_;
    int order_;
    real_t mu_;
    real_t Vp_;         // Plate rate [m/s]
    real_t Wf_;         // Fault depth [m]
    DGMethod method_;

    // DG parameters (following Tandem)
    real_t sigma_;      // Symmetry sign: -1 for SIPG/BR2
    real_t br2_penalty_; // BR2 penalty = Dim + 1 = 3 (2D)

    std::unique_ptr<DG_FECollection> fec_;
    std::unique_ptr<FESpaceType> fes_;

    // Solver and preconditioner
    std::unique_ptr<Solver> solver_;
    std::unique_ptr<Solver> prec_;
    std::unique_ptr<GSSmoother> serial_prec_;  // Lifetime management for serial

    // Fault face information
    Array<int> fault_interior_faces_;
    Array<int> fault_dofs_;
    Vector fault_depths_;
    int num_fault_dofs_;

    // For BR2: precomputed element mass matrix inverses
    std::vector<DenseMatrix> elem_mass_inv_;

    // Work vectors
    Vector X_, B_;

    // Setup methods
    void SetupFESpace();
    void SetupBoundaryMarkers();
    void SetupFaultInfo();        // Identify interior faces at x≈0
    void SetupSolver();

    // Method-specific slip assembly
    void AssembleSlipContributionIP(...);
    void AssembleSlipContributionBR2(...);

    // BR2-specific methods
    void PrecomputeMassMatrixInverses();
    const DenseMatrix &GetMassMatrixInverse(int elem);
    void ComputeLiftingOperator(...);

    // Slip interpolation to fault DOF depths
    real_t InterpolateSlipBC(real_t z, const Vector &slip_bc);
};

extern template class AntiplaneDomainOperator<Mesh>;
#ifdef SEAS_USE_MPI
extern template class AntiplaneDomainOperator<ParMesh>;
#endif

using SerialAntiplaneOperator = AntiplaneDomainOperator<Mesh>;
using ParallelAntiplaneOperator = AntiplaneDomainOperator<ParMesh>;

} // namespace seas
} // namespace mfem
```

### 2.4 Penalty Parameters

Penalty parameters are set in the constructor and used during `Solve()`:

```cpp
// In constructor:
// SIPG sign parameter (always -1 for symmetric interior penalty)
sigma_ = -1.0;

// BR2 penalty: σ = D + 1 = 3 for 2D
br2_penalty_ = 3.0;

// In Solve():
if (method_ == DGMethod::IP)
{
   // IP: kappa = (order + 1)^2 (dimensionless, matches MFEM examples)
   // MFEM's DGDiffusionIntegrator handles h-scaling internally via |n|²/det(J)
   real_t rep_kappa = (order_ + 1) * (order_ + 1);
   a.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(one, sigma_, rep_kappa));
}
else  // BR2
{
   // BR2: custom integrator with lifting operators
   // Requires precomputed element mass matrix inverses
   a.AddInteriorFaceIntegrator(
      new BR2InteriorFaceIntegrator(one, sigma_, elem_mass_inv_, mesh_.Dimension()));
}
```

### 2.5 BR2 Lifting Operator and Mass Matrix Inverses

The BR2 method requires precomputed element mass matrix inverses for the lifting operator:

```cpp
/// Precompute M_K^{-1} for all elements (called once in constructor)
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::PrecomputeMassMatrixInverses() const
{
   int ne = mesh_.GetNE();
   elem_mass_inv_.resize(ne);

   for (int i = 0; i < ne; i++)
   {
      const FiniteElement *fe = fes_->GetFE(i);
      ElementTransformation *T = mesh_.GetElementTransformation(i);

      int ndof = fe->GetDof();
      DenseMatrix M(ndof);
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder());
      Vector shape(ndof);
      M = 0.0;

      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir.IntPoint(j);
         T->SetIntPoint(&ip);
         fe->CalcShape(ip, shape);
         real_t w = ip.weight * T->Weight();
         for (int k = 0; k < ndof; k++)
            for (int l = 0; l < ndof; l++)
               M(k, l) += w * shape(k) * shape(l);
      }

      // Invert mass matrix
      elem_mass_inv_[i].SetSize(ndof);
      DenseMatrixInverse M_inv(M);
      M_inv.GetInverseMatrix(elem_mass_inv_[i]);
   }
   mass_inv_computed_ = true;
}

/// Compute lifting operator for one element adjacent to a face
/// r_e|_K = M_K^{-1} * ∫_e φ_K ⊗ n [[v]] ds
/// For interior faces, use factor 0.5 on each side
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::ComputeLiftingOperator(
   int elem, FaceElementTransformations &FTr, real_t jump, Vector &lift) const
{
   const FiniteElement *fe = fes_->GetFE(elem);
   int ndof = fe->GetDof();
   lift.SetSize(ndof);
   lift = 0.0;

   bool is_elem1 = (elem == FTr.Elem1No);
   int face_order = fe->GetOrder();
   const IntegrationRule &ir = IntRules.Get(FTr.FaceGeom, 2 * face_order + 1);

   Vector shape(ndof);
   Vector int_result(ndof);
   int_result = 0.0;

   for (int p = 0; p < ir.GetNPoints(); p++)
   {
      const IntegrationPoint &ip = ir.IntPoint(p);
      FTr.SetAllIntPoints(&ip);
      const IntegrationPoint &eip = is_elem1 ?
         FTr.GetElement1IntPoint() : FTr.GetElement2IntPoint();
      fe->CalcShape(eip, shape);
      real_t face_weight = ip.weight * FTr.Face->Weight();
      real_t factor = (FTr.Elem2No >= 0) ? 0.5 : 1.0;
      for (int k = 0; k < ndof; k++)
         int_result(k) += factor * shape(k) * jump * face_weight;
   }

   // Apply mass matrix inverse: lift = M^{-1} * int_result
   const DenseMatrix &Minv = GetMassMatrixInverse(elem);
   Minv.Mult(int_result, lift);
}
```

### 2.6 Slip Contribution (Method-Specific)

Slip contribution is assembled separately for IP and BR2 methods. Each fault face has exactly one DOF (the midpoint), so `slip_bc(i)` gives the slip directly for face `i` without depth interpolation.

**IP method** (`AssembleSlipContributionIP`):

RHS for imposed jump `[[u]] = δ`: `σ⟨δ, {∇v·n}⟩_F + κ⟨δ, [[v]]⟩_F`

Key implementation detail: the penalty (jump) term uses the **combined** weight from both elements (`wq_penalty = kappa * |nor|² * (w1 + w2)`), matching MFEM's `DGDiffusionIntegrator` which applies the same combined penalty to all four blocks. Using separate weights gives only half the penalty and under-constrains the slip BC.

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::AssembleSlipContributionIP(
   Vector &rhs, const Vector &slip_bc) const
{
   real_t kappa = (order_ + 1) * (order_ + 1);

   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      // Get DOFs, shape functions, gradients for both elements
      Array<int> dofs1, dofs2;
      fes_->GetElementDofs(FTr->Elem1No, dofs1);
      fes_->GetElementDofs(FTr->Elem2No, dofs2);
      // ... (shape/gradient evaluation)

      // Each fault face has one DOF - use slip directly
      real_t slip_phys = slip_bc(i);
      if (std::abs(slip_phys) < 1e-15) { continue; }

      for (int p = 0; p < ir.GetNPoints(); p++)
      {
         // Normal, shape functions, gradient transforms
         // ...
         real_t slip_imposed = (nor(0) > 0) ? -slip_phys : slip_phys;

         // Element-specific weights for symmetry term
         real_t w1 = ip.weight / (2.0 * detJ1);
         real_t w2 = ip.weight / (2.0 * detJ2);
         // Combined weight for penalty term (matches DGDiffusionIntegrator)
         real_t wq_penalty = kappa * nor_sq * (w1 + w2);

         for (int k = 0; k < dofs1.Size(); k++)
         {
            elvec1(k) += sigma_ * dn1(k) * slip_imposed * w1;
            elvec1(k) += wq_penalty * slip_imposed * shape1(k);
         }
         for (int k = 0; k < dofs2.Size(); k++)
         {
            elvec2(k) += sigma_ * dn2(k) * slip_imposed * w2;
            elvec2(k) -= wq_penalty * slip_imposed * shape2(k);
         }
      }
      // Assemble to global RHS
   }
}
```

**BR2 method** (`AssembleSlipContributionBR2`):

Follows Tandem's `rhs_lift_skeleton` and `rhsFacet`:
1. Consistency: `ε * 0.5 * ∫_F K∇φ·n * δ ds` (same as IP)
2. BR2 lifting: `penalty * ∫_F φ * f_lifted_q ds`

The lifted flux `f_lifted_q` is computed using the BR2 formula:
```
f_lifted[i] = 0.5 * Minv[i] * ∫_F E[i] * δ * n * w ds
f_lifted_q = 0.5 * n · (K * E[0] * f_lifted[0] + K * E[1] * f_lifted[1])
```

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::AssembleSlipContributionBR2(
   Vector &rhs, const Vector &slip_bc) const
{
   int dim = mesh_.Dimension();

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      // ... (DOF/shape setup, slip value from slip_bc(fi))

      // Compute face integrals: ∫_F E * δ * n * w ds for each element
      DenseMatrix face_int1(ndof1, dim), face_int2(ndof2, dim);
      // ... (quadrature loop accumulating face integrals)

      // Apply mass inverse: f_lifted = 0.5 * Minv * face_int
      // ... (for each dimension j)

      // Compute f_lifted_q at quadrature points
      // f_lifted_q[q] = 0.5 * n · (E[0] * f_lifted[0] + E[1] * f_lifted[1])

      // Assemble RHS:
      //   consistency: c1 * ∫_F K∇φ·n * δ ds (c1 = sigma_ * 0.5 = -0.5)
      //   BR2 lifting: penalty * ∫_F φ * f_lifted_q ds
      //   Element 1: +penalty, Element 2: -penalty
      real_t penalty = br2_penalty_;  // D+1 = 3
      // ... (quadrature loop, assemble to global RHS)
   }
}
```

---

## Traction Computation

Traction is computed using the DG numerical flux at interior fault faces. The method takes `slip_bc` as a parameter (needed for penalty corrections in future, currently uses averaged gradient only):

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::ComputeTraction(
   const GridFuncType &displacement, const Vector &slip_bc,
   Vector &traction)
{
   // τ = μ * {{∂u/∂x}} (average of x-gradient from both sides)
   //
   // For IP method: gradient only (penalty residual is too large on coarse mesh)
   // For BR2 method: gradient only (BR2 lifting penalty correction is safe but
   //   currently omitted for simplicity)

   traction.SetSize(num_fault_dofs_);
   traction = 0.0;

   GradientGridFunctionCoefficient grad_u(&displacement);

   int idx = 0;
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      // Single midpoint evaluation, consistent with GetFaultDepths()
      IntegrationPoint ip;
      ip.x = 0.5;

      FTr->SetAllIntPoints(&ip);
      const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
      const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

      // Evaluate gradient in both adjacent elements
      Vector grad1(mesh_.Dimension()), grad2(mesh_.Dimension());
      grad_u.Eval(grad1, *FTr->Elem1, eip1);
      grad_u.Eval(grad2, *FTr->Elem2, eip2);

      // Average gradient (for symmetric flux)
      real_t avg_dudx = 0.5 * (grad1(0) + grad2(0));
      traction(idx++) = mu_ * avg_dudx;
   }
}
```

---

## Unit Test Specifications

All tests are in `tests/unit/test_antiplane.cpp` using a custom test framework (macros `TEST_ASSERT`, `TEST_ASSERT_NEAR`, `RUN_TEST`).

### Mesh Generation Tests

```cpp
test_mesh_creation()          // Mesh dimensions: 2*nx*nz elements, dim=2
test_boundary_attributes()    // Correct attribute counts per boundary
test_test_mesh()              // CreateTestMesh produces 32 elements
test_mms_mesh()               // CreateMMSMesh at different refinement levels
```

### Operator Construction Tests

```cpp
test_operator_construction()  // NumComponents=1, Dimension=2, mu, Wf, fault DOFs > 0
test_fault_dofs()             // 1 DOF per face (midpoint), correct count
test_fault_depths()           // Depths in valid range [-Lz, 0]
```

### Laplace Solution Tests

```cpp
test_zero_slip_solution()        // Zero slip + t=0 → solution ≈ 0
test_uniform_slip_solution()     // Slip=1m → non-trivial, bounded solution
test_full_depth_fault_solution() // Full depth: max_u ≈ 0.5 (slip/2)
test_solution_antisymmetry()     // Sum of u ≈ 0 (antisymmetric about fault)
test_slip_jump_verification()    // Estimated jump ([[u]] ≈ 2*max_u) within 25% of imposed
```

### Traction Computation Tests

```cpp
test_traction_computation()      // Uniform slip → τ ≈ 0; localized slip → τ > 1 MPa
test_stress_kernel_localized()   // Unit slip at one DOF, verify adjacent DOF traction,
                                 // nx/nz resolution studies with graded mesh
test_stress_kernel()             // Source + adjacent DOF traction with analytical comparison
```

### MMS Verification Tests

```cpp
test_mms_convergence()       // Default method: p+1 convergence rate (> 1.5)
test_mms_convergence_ip()    // IP method: 4-level convergence study
test_mms_convergence_br2()   // BR2 method: 4-level convergence study
```

MMS exact solution: `u(x,z) = sin(π(x+Lx)/(2Lx)) · exp(πz/Lz)` (satisfies `∇²u = 0` when `Lx = Lz`)

### DG Method Comparison Tests

```cpp
test_dg_method_selection()    // Enum correctly set (IP/BR2)
test_ip_penalty_formula()     // IP produces non-zero bounded solution
test_br2_method()             // BR2 produces non-zero bounded solution
test_ip_br2_consistency()     // IP and BR2 within 20% relative difference
test_br2_integrator_lifting() // BR2 with lifting produces valid solution
test_br2_ip_convergence()     // Both converge, difference decreases with refinement
test_tandem_ip_penalty_value() // Tandem penalty formula verification
```

---

## File Organization

```
domain/
├── domain_operator.hpp           # Abstract interface (template <typename MeshType>)
├── antiplane_operator.hpp        # Antiplane DG operator (includes implementation)
├── bp2_mesh.hpp                  # BP2MeshGenerator class with boundary attributes

integrator/
└── dg_br2_integrator.hpp         # BR2InteriorFaceIntegrator, BR2BoundaryFaceIntegrator

common/
└── seas_types.hpp                # Type aliases (FESpaceForMesh, IsParallelMesh, etc.)
```

Notes:
- The `DGMethod` enum is defined in `antiplane_operator.hpp` (not a separate file)
- Implementation is inline in the header (template class)
- No separate `antiplane_operator_impl.hpp` or `penalty_calculator.hpp`

---

---

## Mesh Creation and Fault Face Identification

### BP2 Boundary Attributes

All outer boundaries use Natural BC (zero traction). Tectonic loading is applied via prescribed `Vp` on fault faces below `Wf` in the fault operator.

```cpp
// domain/bp2_mesh.hpp

struct BP2BoundaryAttributes
{
   static constexpr int FARFIELD_LEFT = 1;   ///< Left far-field (x = -Lx)
   static constexpr int FARFIELD_RIGHT = 2;  ///< Right far-field (x = +Lx)
   static constexpr int FREE_SURFACE = 3;    ///< Free surface (z = 0)
   static constexpr int BOTTOM = 4;          ///< Bottom (z = -Lz)
   // Note: The fault is NOT a boundary - it's interior faces at x = 0
};
```

### BP2 Mesh Generator

The mesh is created by the `BP2MeshGenerator` class, which provides static factory methods:

```cpp
// domain/bp2_mesh.hpp

class BP2MeshGenerator
{
public:
   /// Parameters for mesh generation
   struct Parameters
   {
      real_t Lx = 100.0e3;   ///< Domain half-width [m]
      real_t Lz = 100.0e3;   ///< Domain depth [m]
      real_t Wf = 40.0e3;    ///< Rate-state fault depth [m]
      int nx = 100;          ///< Elements per half (2*nx total in x)
      int nz = 100;          ///< Elements in z
      real_t grading_x = 1.0; ///< Grading in x (1.0 = uniform)
      real_t grading_z = 1.0; ///< Grading in z (1.0 = uniform)
   };

   /// Create uniform mesh with default parameters
   static std::unique_ptr<Mesh> Create();

   /// Create uniform mesh with specified parameters
   static std::unique_ptr<Mesh> Create(const Parameters &params);

   /// Create graded mesh (sinh/tanh grading near fault and surface)
   static std::unique_ptr<Mesh> CreateGraded(const Parameters &params);

   /// Create small test mesh (8x4, [-10km, +10km] x [-10km, 0])
   static std::unique_ptr<Mesh> CreateTestMesh();

   /// Create MMS mesh for convergence studies ([-1, +1] x [-1, 0])
   static std::unique_ptr<Mesh> CreateMMSMesh(int refinement_level);

   /// Check if point is on fault (x ≈ 0, z > -Wf)
   static bool IsOnFault(real_t x, real_t z, real_t Wf, real_t tol = 1e-10);

private:
   static void ShiftMeshCoordinates(Mesh &mesh, real_t Lx, real_t Lz);
   static void SetBoundaryAttributes(Mesh &mesh, const Parameters &params);
   static void ApplyGrading(Mesh &mesh, const Parameters &params);
};
```

The `Create()` method builds a Cartesian mesh with `2*nx` elements in x and `nx` in z, then shifts coordinates so `x ∈ [-Lx, +Lx]` and `z ∈ [-Lz, 0]`. This ensures element boundaries align at `x = 0` for a clean fault interface.

`CreateGraded()` applies sinh-based grading in x (concentrating elements near the fault at `x = 0`) and tanh-based grading in z (concentrating elements near the surface at `z = 0`).

### Fault Face Identification

The fault is NOT a boundary - it's an **interior face** where we prescribe a jump `[[u]] = δ`. The `SetupFaultInfo()` method identifies all interior faces at `x ≈ 0` (the full fault interface from `z = 0` to `z = -Lz`):

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::SetupFaultInfo()
{
   fault_interior_faces_.SetSize(0);

   int num_faces = mesh_.GetNumFaces();
   for (int f = 0; f < num_faces; f++)
   {
      FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }  // Boundary face

      if (IsFaultFace(f))
      {
         fault_interior_faces_.Append(f);
      }
   }

   // One fault DOF per face, evaluated at the face midpoint.
   // Using (order+1) Gauss-Lobatto points would place DOFs at shared
   // vertices between adjacent faces. With DG order 1 every DOF sits
   // at a shared vertex, so two DOFs at the same depth evolve
   // independently — causing spurious single-DOF nucleation.
   // A single midpoint DOF per face eliminates the duplication.
   num_fault_dofs_ = fault_interior_faces_.Size();

   fault_dofs_.SetSize(num_fault_dofs_);
   for (int i = 0; i < num_fault_dofs_; i++)
   {
      fault_dofs_[i] = i;
   }
}

template <typename MeshType>
bool AntiplaneDomainOperator<MeshType>::IsFaultFace(int face) const
{
   real_t x, z;
   GetFaceCenter(face, x, z);
   // Full fault interface at x ≈ 0 (rate-state for z > -Wf; Vp for z < -Wf)
   const real_t tol = 1e-10 * std::max(Wf_, 1.0);
   return std::abs(x) < tol;
}
```

### Fault Depth Computation

Fault depths are lazily computed and cached. They are eagerly precomputed in the constructor to avoid invalidating `FaceElementTransformations` pointers during slip assembly (since `GetInteriorFaceTransformations` returns a pointer to a shared internal object):

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::GetFaultDepths(Vector &depths) const
{
   if (!fault_depths_computed_)
   {
      fault_depths_.SetSize(num_fault_dofs_);
      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         IntegrationPoint ip;
         ip.x = 0.5;  // Face midpoint
         FTr->SetAllIntPoints(&ip);
         Vector coords(mesh_.Dimension());
         FTr->Face->Transform(ip, coords);
         fault_depths_(i) = coords(1);  // z-coordinate
      }
      fault_depths_computed_ = true;
   }
   depths = fault_depths_;
}
```

### Slip Interpolation

`InterpolateSlipBC` provides linear interpolation of slip values at arbitrary depths along the fault, using fault DOF depths as knot points:

```cpp
template <typename MeshType>
real_t AntiplaneDomainOperator<MeshType>::InterpolateSlipBC(
   real_t z, const Vector &slip_bc) const
{
   // Linear interpolation provides O(h²) accuracy for smooth slip profiles.
   // Finds the bracketing interval and interpolates, with nearest-neighbor
   // fallback for extrapolation.
   // ...
}
```

---

## Boundary Condition Implementation

### All Natural (Zero-Traction) BCs

All **four outer boundaries** use Natural BC (zero traction). There is no Dirichlet BC at the bottom boundary. Tectonic loading is applied by prescribing slip rate `Vp` on the fault interface below `Wf`, handled by the fault operator (Phase 3).

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::SetupBoundaryMarkers()
{
   // All boundaries have Natural BC (zero traction):
   // - FARFIELD_LEFT (x=-Lx): zero traction at far-field
   // - FARFIELD_RIGHT (x=+Lx): zero traction at far-field
   // - FREE_SURFACE (z=0): zero traction at free surface
   // - BOTTOM (z=-Lz): zero traction
   //
   // No Dirichlet boundary markers are needed.
   // In DG with all-Natural BCs, no boundary face integrators are added
   // during Solve() — only interior face integrators contribute.
}
```

### MMS Boundary Conditions

For Method of Manufactured Solutions testing, `SolveMMS()` applies Dirichlet BCs on **all** boundaries using `DGDiffusionIntegrator` (boundary faces) and `DGDirichletLFIntegrator` (RHS):

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::SolveMMS(
   Coefficient &dirichlet_bc, GridFuncType &displacement)
{
   // Mark all boundaries for Dirichlet
   Array<int> all_bdr_marker(mesh_.bdr_attributes.Max());
   all_bdr_marker = 1;

   BilinFormType a(fes_.get());
   ConstantCoefficient one(1.0);
   a.AddDomainIntegrator(new DiffusionIntegrator(one));

   // Boundary face: standard DG penalty (both IP and BR2 use same bdr treatment)
   real_t bdr_kappa = (order_ + 1) * (order_ + 1);

   // Interior face: method-specific integrator
   if (method_ == DGMethod::IP)
      a.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(one, sigma_, bdr_kappa));
   else
      a.AddInteriorFaceIntegrator(
         new BR2InteriorFaceIntegrator(one, sigma_, elem_mass_inv_, mesh_.Dimension()));

   // Boundary face integrator for Dirichlet
   a.AddBdrFaceIntegrator(
      new DGDiffusionIntegrator(one, sigma_, bdr_kappa), all_bdr_marker);

   // RHS: exact solution as Dirichlet BC
   LinFormType b(fes_.get());
   b.AddBdrFaceIntegrator(
      new DGDirichletLFIntegrator(dirichlet_bc, one, sigma_, bdr_kappa),
      all_bdr_marker);
   // ... assemble and solve
}
```

---

## Implementation Notes

### DG Method Implementation

Both IP and BR2 methods are fully implemented:

1. **IP method** (default):
   - Uses MFEM's built-in `DGDiffusionIntegrator` with penalty `κ = (p+1)²`
   - MFEM handles h-scaling internally via `|n|²/det(J)`
   - Interior face + boundary face integrators

2. **BR2 method**:
   - Custom `BR2InteriorFaceIntegrator` class (in `integrator/dg_br2_integrator.hpp`)
   - Implements full lifting operator following Tandem's `poisson.py` kernels
   - Fixed penalty `σ = D + 1 = 3` for 2D
   - Precomputed element mass matrix inverses stored as `std::vector<DenseMatrix>`
   - Also includes `BR2BoundaryFaceIntegrator` for boundary face assembly

### BR2 Integrator Classes

```cpp
// integrator/dg_br2_integrator.hpp

/// Interior face integrator with BR2 lifting operators
class BR2InteriorFaceIntegrator : public BilinearFormIntegrator
{
public:
   BR2InteriorFaceIntegrator(Coefficient &K, real_t epsilon,
                             const std::vector<DenseMatrix> &elem_mass_inv,
                             int dim = 2);

   void AssembleFaceMatrix(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Trans,
                           DenseMatrix &elmat) override;
   // Assembly follows Tandem's assembleSurface:
   // a[x][y] = c0[y]*∫K∇φ_x·n φ_y + c1[x]*∫K∇φ_y·n φ_x + c2[|x-y|]*∫φ_x L_q[y]
   // c0 = -0.5, c1 = ε*0.5, c2[0] = +σ (same elem), c2[1] = -σ (cross elem)
};

/// Boundary face integrator with BR2 lifting
class BR2BoundaryFaceIntegrator : public BilinearFormIntegrator
{
   // Similar structure for boundary faces (single element, full factor)
};
```

### Solver Setup

The solver is configured based on serial/parallel context using `if constexpr`:

```cpp
if constexpr (IsParallelMesh<MeshType>::value)
{
   // HypreBoomerAMG + HyprePCG
}
else
{
   // CGSolver + GSSmoother (created per-solve, stored as member for lifetime)
}
```

The serial preconditioner (`GSSmoother`) must be stored as a member because `CGSolver::SetPreconditioner` stores a pointer, not a copy.

---

## Acceptance Criteria

| Component | Tests | Pass Criteria |
|-----------|-------|---------------|
| Mesh generation | 4 tests | Correct dimensions, boundary attributes, element counts |
| Operator construction | 3 tests | Correct properties, fault DOF count and depths |
| Solution accuracy | 5 tests | Zero slip → zero, full depth → slip/2, antisymmetry, jump |
| Traction | 3 tests | Uniform → near zero, localized → significant, stress kernel |
| MMS convergence | 3 tests (default, IP, BR2) | Rate > 1.5 for linear elements |
| DG method comparison | 7 tests | BR2 ≈ IP within 20%, both converge |
| **Total** | **25 tests** | All pass |
