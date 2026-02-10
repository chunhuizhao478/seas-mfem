# Phase Field Fracture (PFF) Implementation Plan for MFEM

## Overview

This document outlines the plan to implement a Phase Field Fracture method in MFEM, transferring functionality from the MOOSE/FARMS implementation. The implementation targets a three-point bending test for pure solid mechanics.

### Source Reference
- **MOOSE Source**: `/Users/chunhuizhao/projects/farms_cdms/threepointbending/code/case_wohole/`
- **Key Files**: `elasticity.i`, `fracture.i`

### Key Differences from MOOSE Implementation
- **No AD (Automatic Differentiation)**: Manual Jacobian computation
- **MFEM patterns**: Use NonlinearForm, custom Operators, and Integrators
- **Staggered scheme**: Fixed-point iteration between elasticity and damage subproblems

---

## Mathematical Formulation

### Governing Equations

**1. Elasticity Sub-problem (Momentum Balance)**

Given damage field `d`, find displacement `u` such that:
```
∇ · σ = 0   in Ω
```

where the degraded stress is:
```
σ = g(d) · σ₀(ε)
σ₀ = λ tr(ε) I + 2μ ε    (isotropic linear elasticity)
ε = ½(∇u + ∇uᵀ)          (small strain)
```

With spectral decomposition for tension-compression split:
```
ψ = g(d) · ψ⁺(ε) + ψ⁻(ε)
ψ⁺ = ½λ⟨tr(ε)⟩₊² + μ tr(ε⁺²)   (tensile)
ψ⁻ = ½λ⟨tr(ε)⟩₋² + μ tr(ε⁻²)   (compressive)
```

**2. Damage Sub-problem (Phase Field Evolution)**

Given strain energy history `H`, find damage `d` such that:
```
∂ψ/∂d + ∂φ/∂d = 0
```

Which yields the variational form:
```
∫_Ω (Gc·l/c₀)∇d·∇δd dΩ + ∫_Ω (2·Gc/(c₀·l)·d + g'(d)·H)δd dΩ = 0
```

**3. Material Parameters**
```
E = 40 MPa           (Young's modulus)
ν = 0.25             (Poisson's ratio)
K = E/(3(1-2ν))      (Bulk modulus)
G = E/(2(1+ν))       (Shear modulus)
l = 5×10⁻⁵ m         (Regularization length)
σₜ = 6.43 MPa        (Tensile strength)
Gc = 8lσₜ²/(3E)      (Fracture toughness, derived)
```

**4. Degradation Function (AT2 model)**
```
g(d) = (1-d)²(1-η) + η
α(d) = d²                (crack geometric function)
c₀ = 2                   (normalization constant for AT2)
η = 10⁻⁶                 (residual stiffness)
```

**5. Boundary Conditions**
- Top loading: prescribed vertical displacement `u_y = -10⁻⁴ t`
- Bottom supports: fixed `u_y = 0`
- Damage: homogeneous Neumann (∂d/∂n = 0)

**6. Irreversibility Constraint**
```
d(t) ≥ d(t - Δt)   ∀ t
```

---

## Implementation Architecture

### Directory Structure
```
miniapps/pff/
├── pff.cpp                           # Main driver
├── pff_solver.hpp                    # Staggered solver class
├── pff_solver.cpp
├── operators/
│   ├── elasticity_operator.hpp       # Degraded elasticity operator
│   ├── elasticity_operator.cpp
│   ├── damage_operator.hpp           # Phase field damage operator
│   └── damage_operator.cpp
├── materials/
│   ├── pff_material.hpp              # Material parameters container
│   ├── degradation_function.hpp      # g(d), g'(d)
│   └── spectral_decomposition.hpp    # Strain decomposition
├── integrators/
│   ├── degraded_elasticity_integrator.hpp
│   ├── degraded_elasticity_integrator.cpp
│   ├── damage_diffusion_integrator.hpp
│   ├── damage_diffusion_integrator.cpp
│   ├── damage_source_integrator.hpp
│   └── damage_source_integrator.cpp
├── tests/
│   ├── test_degradation_function.cpp
│   ├── test_spectral_decomposition.cpp
│   ├── test_elasticity_operator.cpp
│   ├── test_damage_operator.cpp
│   └── test_staggered_solver.cpp
├── examples/
│   ├── three_point_bending.cpp       # Main demo
│   └── mesh/
│       └── beam.msh
├── CMakeLists.txt
└── makefile
```

---

## Phase 1: Material Models and Utilities

### Task 1.1: Degradation Function
**File**: `materials/degradation_function.hpp`

Implement the AT2 degradation function and its derivatives:
```cpp
class DegradationFunction {
public:
    DegradationFunction(real_t eta = 1e-6, real_t p = 2.0);

    // g(d) = (1-d)^p * (1-eta) + eta
    real_t Eval(real_t d) const;

    // g'(d) = -p * (1-d)^(p-1) * (1-eta)
    real_t EvalDerivative(real_t d) const;

    // g''(d) = p * (p-1) * (1-d)^(p-2) * (1-eta)
    real_t EvalSecondDerivative(real_t d) const;

private:
    real_t eta_;  // Residual stiffness
    real_t p_;    // Exponent
};
```

**Unit Test**: `tests/test_degradation_function.cpp`
```cpp
TEST(DegradationFunction, EvalAtZero) {
    DegradationFunction g;
    EXPECT_NEAR(g.Eval(0.0), 1.0, 1e-10);
}

TEST(DegradationFunction, EvalAtOne) {
    DegradationFunction g(1e-6);
    EXPECT_NEAR(g.Eval(1.0), 1e-6, 1e-10);
}

TEST(DegradationFunction, DerivativeAtZero) {
    DegradationFunction g;
    EXPECT_NEAR(g.EvalDerivative(0.0), -2.0 * (1.0 - 1e-6), 1e-10);
}

TEST(DegradationFunction, MonotonicallyDecreasing) {
    DegradationFunction g;
    for (real_t d = 0.0; d < 1.0; d += 0.1) {
        EXPECT_LT(g.EvalDerivative(d), 0.0);
    }
}

TEST(DegradationFunction, SecondDerivativePositive) {
    DegradationFunction g;
    for (real_t d = 0.01; d < 0.99; d += 0.1) {
        EXPECT_GT(g.EvalSecondDerivative(d), 0.0);
    }
}
```

---

### Task 1.2: Spectral Decomposition
**File**: `materials/spectral_decomposition.hpp`

Implement strain tensor spectral decomposition for tension-compression split:
```cpp
template <int dim>
class SpectralDecomposition {
public:
    // Decompose strain into positive and negative parts
    // ε = ε⁺ + ε⁻
    static void Decompose(const DenseMatrix &strain,
                          DenseMatrix &strain_pos,
                          DenseMatrix &strain_neg);

    // Compute positive strain energy density
    // ψ⁺ = λ/2 ⟨tr(ε)⟩₊² + μ tr(ε⁺²)
    static real_t PositiveEnergy(const DenseMatrix &strain,
                                  real_t lambda, real_t mu);

    // Compute negative strain energy density
    static real_t NegativeEnergy(const DenseMatrix &strain,
                                  real_t lambda, real_t mu);

    // Compute stress from positive strain only (for degradation)
    static void PositiveStress(const DenseMatrix &strain,
                               real_t lambda, real_t mu,
                               DenseMatrix &stress_pos);

    // Compute stress from negative strain only
    static void NegativeStress(const DenseMatrix &strain,
                               real_t lambda, real_t mu,
                               DenseMatrix &stress_neg);

private:
    // Helper: compute eigenvalues and eigenvectors
    static void ComputeEigen(const DenseMatrix &A,
                             Vector &eigenvalues,
                             DenseMatrix &eigenvectors);
};
```

**Unit Test**: `tests/test_spectral_decomposition.cpp`
```cpp
TEST(SpectralDecomposition, PureTension) {
    // Uniaxial tension: ε = diag(ε₁, -ν*ε₁, -ν*ε₁)
    DenseMatrix strain(3, 3);
    strain = 0.0;
    strain(0, 0) = 0.01;  // 1% tensile strain
    strain(1, 1) = -0.0025;  // Poisson contraction
    strain(2, 2) = -0.0025;

    DenseMatrix strain_pos(3, 3), strain_neg(3, 3);
    SpectralDecomposition<3>::Decompose(strain, strain_pos, strain_neg);

    // Positive strain should capture tensile component
    EXPECT_GT(strain_pos(0, 0), 0.0);
}

TEST(SpectralDecomposition, PureCompression) {
    DenseMatrix strain(3, 3);
    strain = 0.0;
    strain(0, 0) = -0.01;  // Compressive
    strain(1, 1) = 0.0025;
    strain(2, 2) = 0.0025;

    DenseMatrix strain_pos(3, 3), strain_neg(3, 3);
    SpectralDecomposition<3>::Decompose(strain, strain_pos, strain_neg);

    // Negative energy should dominate
    real_t psi_neg = SpectralDecomposition<3>::NegativeEnergy(strain, 1.0, 1.0);
    real_t psi_pos = SpectralDecomposition<3>::PositiveEnergy(strain, 1.0, 1.0);
    EXPECT_GT(psi_neg, psi_pos);
}

TEST(SpectralDecomposition, Symmetry) {
    DenseMatrix strain(3, 3);
    strain.Randomize();
    // Make symmetric
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            strain(j, i) = strain(i, j);
        }
    }

    DenseMatrix strain_pos(3, 3), strain_neg(3, 3);
    SpectralDecomposition<3>::Decompose(strain, strain_pos, strain_neg);

    // Result should be symmetric
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_NEAR(strain_pos(i, j), strain_pos(j, i), 1e-12);
            EXPECT_NEAR(strain_neg(i, j), strain_neg(j, i), 1e-12);
        }
    }
}

TEST(SpectralDecomposition, Additivity) {
    // ε⁺ + ε⁻ = ε
    DenseMatrix strain(3, 3);
    strain.Randomize();
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            strain(j, i) = strain(i, j);
        }
    }

    DenseMatrix strain_pos(3, 3), strain_neg(3, 3);
    SpectralDecomposition<3>::Decompose(strain, strain_pos, strain_neg);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_NEAR(strain_pos(i, j) + strain_neg(i, j),
                        strain(i, j), 1e-12);
        }
    }
}

TEST(SpectralDecomposition, EnergyNonNegative) {
    DenseMatrix strain(3, 3);
    strain.Randomize();
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            strain(j, i) = strain(i, j);
        }
    }

    real_t psi_pos = SpectralDecomposition<3>::PositiveEnergy(strain, 1.0, 1.0);
    real_t psi_neg = SpectralDecomposition<3>::NegativeEnergy(strain, 1.0, 1.0);

    EXPECT_GE(psi_pos, 0.0);
    EXPECT_GE(psi_neg, 0.0);
}
```

---

### Task 1.3: Material Parameters Container
**File**: `materials/pff_material.hpp`

```cpp
struct PFFMaterialParameters {
    // Elastic properties
    real_t E;       // Young's modulus
    real_t nu;      // Poisson's ratio
    real_t lambda;  // First Lamé parameter
    real_t mu;      // Second Lamé parameter (shear modulus)
    real_t K;       // Bulk modulus

    // Fracture properties
    real_t Gc;      // Fracture toughness
    real_t l;       // Regularization length
    real_t c0;      // Normalization constant (= 2 for AT2)

    // Degradation parameters
    real_t eta;     // Residual stiffness
    real_t p;       // Degradation exponent

    // Constructor with derived quantities
    PFFMaterialParameters(real_t E_, real_t nu_, real_t l_,
                          real_t sigma_t, real_t eta_ = 1e-6);

    // Compute Lamé parameters from E, nu
    void ComputeLameParameters();

    // Compute Gc from strength-based formula
    void ComputeFractureToughness(real_t sigma_t);
};
```

**Unit Test**: `tests/test_pff_material.cpp`
```cpp
TEST(PFFMaterial, LameParametersFromENu) {
    PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

    real_t expected_lambda = mat.E * mat.nu / ((1 + mat.nu) * (1 - 2 * mat.nu));
    real_t expected_mu = mat.E / (2 * (1 + mat.nu));

    EXPECT_NEAR(mat.lambda, expected_lambda, 1e-6);
    EXPECT_NEAR(mat.mu, expected_mu, 1e-6);
}

TEST(PFFMaterial, BulkModulus) {
    PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);
    real_t expected_K = mat.E / (3 * (1 - 2 * mat.nu));
    EXPECT_NEAR(mat.K, expected_K, 1e-6);
}

TEST(PFFMaterial, FractureToughness) {
    // Gc = 8*l*sigma_t^2 / (3*E)
    real_t E = 40e6, l = 5e-5, sigma_t = 6.43e6;
    PFFMaterialParameters mat(E, 0.25, l, sigma_t);

    real_t expected_Gc = 8 * l * sigma_t * sigma_t / (3 * E);
    EXPECT_NEAR(mat.Gc, expected_Gc, 1e-10);
}
```

---

## Phase 2: Integrators

### Task 2.1: Degraded Elasticity Integrator
**File**: `integrators/degraded_elasticity_integrator.hpp`

Nonlinear integrator for the weak form:
```
∫_Ω σ(u, d) : ε(v) dΩ
```

where `σ = g(d)σ⁺ + σ⁻` (spectral decomposition).

```cpp
class DegradedElasticityIntegrator : public NonlinearFormIntegrator {
public:
    DegradedElasticityIntegrator(const PFFMaterialParameters &mat,
                                  const GridFunction &damage);

    // Compute element residual
    virtual void AssembleElementVector(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       const Vector &elfun,
                                       Vector &elvect) override;

    // Compute element Jacobian
    virtual void AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &Tr,
                                     const Vector &elfun,
                                     DenseMatrix &elmat) override;

    // Compute strain energy at quadrature points (for damage driving force)
    void ComputeStrainEnergy(const FiniteElement &el,
                             ElementTransformation &Tr,
                             const Vector &elfun,
                             Vector &psi_active);

private:
    const PFFMaterialParameters &mat_;
    const GridFunction &damage_;
    DegradationFunction g_;

    // Helper methods
    void ComputeStrain(const DenseMatrix &grad_u, DenseMatrix &strain);
    void ComputeStress(const DenseMatrix &strain, real_t d,
                       DenseMatrix &stress);
};
```

**Unit Test**: `tests/test_degraded_elasticity_integrator.cpp`
```cpp
TEST(DegradedElasticityIntegrator, ZeroDamageMatchesLinearElasticity) {
    // Setup simple 2D mesh
    Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL);
    H1_FECollection fec(1, 2);
    FiniteElementSpace fes(&mesh, &fec, 2);  // vector field

    // Zero damage field
    L2_FECollection l2_fec(0, 2);
    FiniteElementSpace l2_fes(&mesh, &l2_fec);
    GridFunction damage(&l2_fes);
    damage = 0.0;

    PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);

    // Compare with standard elasticity
    DegradedElasticityIntegrator pff_integ(mat, damage);
    ElasticityIntegrator std_integ(mat.lambda, mat.mu);

    // ... compare element matrices
}

TEST(DegradedElasticityIntegrator, FullDamageGivesResidualStiffness) {
    // With d = 1, stiffness should be eta * original
    Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL);
    H1_FECollection fec(1, 2);
    FiniteElementSpace fes(&mesh, &fec, 2);

    L2_FECollection l2_fec(0, 2);
    FiniteElementSpace l2_fes(&mesh, &l2_fec);
    GridFunction damage(&l2_fes);
    damage = 1.0;

    PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);
    DegradedElasticityIntegrator integ(mat, damage);

    // Element stiffness should be scaled by eta for tensile part
    // ... verify
}

TEST(DegradedElasticityIntegrator, JacobianSymmetric) {
    // Jacobian matrix should be symmetric
    // ... setup and verify symmetry
}

TEST(DegradedElasticityIntegrator, StrainEnergyPositive) {
    // Strain energy should always be non-negative
    // ... verify for various deformation states
}
```

---

### Task 2.2: Damage Diffusion Integrator
**File**: `integrators/damage_diffusion_integrator.hpp`

Bilinear integrator for:
```
∫_Ω (Gc·l/c₀) ∇d · ∇δd dΩ
```

```cpp
class DamageDiffusionIntegrator : public BilinearFormIntegrator {
public:
    DamageDiffusionIntegrator(real_t Gc, real_t l, real_t c0 = 2.0);

    virtual void AssembleElementMatrix(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       DenseMatrix &elmat) override;

private:
    real_t coeff_;  // Gc * l / c0
};
```

**Unit Test**: `tests/test_damage_diffusion_integrator.cpp`
```cpp
TEST(DamageDiffusionIntegrator, MatchesDiffusionIntegrator) {
    // Should be equivalent to DiffusionIntegrator with constant coefficient
    real_t Gc = 1.0, l = 0.1, c0 = 2.0;
    real_t expected_coeff = Gc * l / c0;

    DamageDiffusionIntegrator pff_diff(Gc, l, c0);
    ConstantCoefficient coeff(expected_coeff);
    DiffusionIntegrator std_diff(coeff);

    // Compare element matrices
    // ...
}

TEST(DamageDiffusionIntegrator, MatrixSPD) {
    // Stiffness matrix should be symmetric positive semi-definite
    // ...
}
```

---

### Task 2.3: Damage Source Integrator
**File**: `integrators/damage_source_integrator.hpp`

Nonlinear integrator for:
```
∫_Ω (2·Gc/(c₀·l)·d + g'(d)·H) δd dΩ
```

where `H = max(ψ⁺(τ), τ ≤ t)` is the strain energy history.

```cpp
class DamageSourceIntegrator : public NonlinearFormIntegrator {
public:
    DamageSourceIntegrator(const PFFMaterialParameters &mat,
                           const GridFunction &strain_energy_history);

    virtual void AssembleElementVector(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       const Vector &elfun,
                                       Vector &elvect) override;

    virtual void AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &Tr,
                                     const Vector &elfun,
                                     DenseMatrix &elmat) override;

private:
    const PFFMaterialParameters &mat_;
    const GridFunction &H_;  // Strain energy history
    DegradationFunction g_;
};
```

**Unit Test**: `tests/test_damage_source_integrator.cpp`
```cpp
TEST(DamageSourceIntegrator, ZeroHistoryZeroSource) {
    // With H = 0, only the 2*Gc/(c0*l)*d term contributes
    // ...
}

TEST(DamageSourceIntegrator, JacobianCorrect) {
    // Verify Jacobian with finite differences
    // ...
}

TEST(DamageSourceIntegrator, ConsistentWithEnergy) {
    // Verify that -∂E/∂d matches the residual
    // ...
}
```

---

## Phase 3: Operators

### Task 3.1: Elasticity Operator
**File**: `operators/elasticity_operator.hpp`

```cpp
class ElasticityOperator : public Operator {
public:
    ElasticityOperator(ParFiniteElementSpace &fes,
                       const PFFMaterialParameters &mat,
                       Array<int> &ess_bdr);

    // Set the current damage field
    void SetDamageField(const ParGridFunction &d);

    // R(u) = ∫ σ(u,d) : ε(v) dΩ - f_ext
    virtual void Mult(const Vector &u, Vector &r) const override;

    // K = dR/du
    virtual Operator &GetGradient(const Vector &u) const override;

    // Compute strain energy for damage driving force
    void ComputeStrainEnergyDensity(const Vector &u,
                                    ParGridFunction &psi_active) const;

    // Essential boundary DOFs
    void SetEssentialTrueDofs(const Array<int> &ess_tdof_list);

    // Prescribed displacement
    void SetPrescribedDisplacement(const Vector &u_prescribed);

private:
    ParFiniteElementSpace &fes_;
    const PFFMaterialParameters &mat_;
    mutable ParGridFunction d_gf_;  // Current damage
    mutable ParNonlinearForm *nlf_;
    mutable OperatorPtr K_;
    Array<int> ess_tdof_list_;
};
```

**Unit Test**: `tests/test_elasticity_operator.cpp`
```cpp
TEST(ElasticityOperator, PatchTest) {
    // Linear displacement field should give zero residual (minus BCs)
    ParMesh pmesh(MPI_COMM_WORLD, Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON));
    H1_FECollection fec(1, 3);
    ParFiniteElementSpace fes(&pmesh, &fec, 3);

    PFFMaterialParameters mat(40e6, 0.25, 5e-5, 6.43e6);
    Array<int> ess_bdr(pmesh.bdr_attributes.Max());
    ess_bdr = 0;

    ElasticityOperator op(fes, mat, ess_bdr);

    // Set zero damage
    L2_FECollection l2_fec(0, 3);
    ParFiniteElementSpace l2_fes(&pmesh, &l2_fec);
    ParGridFunction d(&l2_fes);
    d = 0.0;
    op.SetDamageField(d);

    // Linear displacement: u = Ax + b
    ParGridFunction u(&fes);
    // ... set linear displacement

    Vector r(fes.GetTrueVSize());
    Vector u_true;
    u.GetTrueDofs(u_true);
    op.Mult(u_true, r);

    EXPECT_NEAR(r.Norml2(), 0.0, 1e-10);
}

TEST(ElasticityOperator, NewtonConvergence) {
    // Newton solver should converge for simple loading
    // ...
}

TEST(ElasticityOperator, SymmetricJacobian) {
    // Verify K is symmetric
    // ...
}
```

---

### Task 3.2: Damage Operator
**File**: `operators/damage_operator.hpp`

```cpp
class DamageOperator : public Operator {
public:
    DamageOperator(ParFiniteElementSpace &fes,
                   const PFFMaterialParameters &mat);

    // Update strain energy history field
    void UpdateStrainEnergyHistory(const ParGridFunction &psi_active);

    // R(d) = diffusion + source terms
    virtual void Mult(const Vector &d, Vector &r) const override;

    // K = dR/dd
    virtual Operator &GetGradient(const Vector &d) const override;

    // Apply irreversibility constraint: d >= d_old
    void ApplyIrreversibility(Vector &d, const Vector &d_old) const;

private:
    ParFiniteElementSpace &fes_;
    const PFFMaterialParameters &mat_;
    mutable ParGridFunction H_gf_;  // Strain energy history max
    mutable ParNonlinearForm *nlf_;
    mutable OperatorPtr K_;
};
```

**Unit Test**: `tests/test_damage_operator.cpp`
```cpp
TEST(DamageOperator, ZeroEnergyZeroDamage) {
    // With H = 0, equilibrium damage should be d = 0
    // ...
}

TEST(DamageOperator, IrreversibilityConstraint) {
    // d should never decrease
    // ...
}

TEST(DamageOperator, NewtonConvergence) {
    // Newton solver should converge
    // ...
}

TEST(DamageOperator, DamageInRange) {
    // 0 <= d <= 1 always
    // ...
}
```

---

## Phase 4: Staggered Solver

### Task 4.1: PFF Staggered Solver
**File**: `pff_solver.hpp`

```cpp
class PFFSolver {
public:
    PFFSolver(ParMesh &mesh, const PFFMaterialParameters &mat, int order = 1);

    // Setup boundary conditions
    void SetElasticityEssentialBdr(const Array<int> &bdr_attr);
    void SetPrescribedDisplacement(VectorCoefficient &u_bc);

    // Solve one time step with fixed-point iteration
    // Returns number of iterations
    int SolveTimeStep(real_t dt);

    // Advance to next time step
    void AdvanceTime(real_t dt);

    // Get solution fields
    const ParGridFunction &GetDisplacement() const { return u_; }
    const ParGridFunction &GetDamage() const { return d_; }
    const ParGridFunction &GetStrainEnergy() const { return psi_; }

    // Fixed-point iteration parameters
    void SetMaxIterations(int max_it) { max_fp_it_ = max_it; }
    void SetRelativeTolerance(real_t tol) { fp_rel_tol_ = tol; }
    void SetAbsoluteTolerance(real_t tol) { fp_abs_tol_ = tol; }

    // Newton solver parameters
    void SetNewtonRelTol(real_t tol);
    void SetNewtonAbsTol(real_t tol);
    void SetNewtonMaxIter(int max_it);

private:
    // Mesh and spaces
    ParMesh &mesh_;
    H1_FECollection u_fec_;      // Displacement FE collection
    H1_FECollection d_fec_;      // Damage FE collection
    L2_FECollection psi_fec_;    // Strain energy (DG)
    ParFiniteElementSpace u_fes_;
    ParFiniteElementSpace d_fes_;
    ParFiniteElementSpace psi_fes_;

    // Solution fields
    ParGridFunction u_;     // Displacement
    ParGridFunction d_;     // Damage
    ParGridFunction d_old_; // Previous damage (for irreversibility)
    ParGridFunction psi_;   // Strain energy density
    ParGridFunction H_;     // Strain energy history (max)

    // Operators
    std::unique_ptr<ElasticityOperator> elasticity_op_;
    std::unique_ptr<DamageOperator> damage_op_;

    // Newton solvers
    NewtonSolver newton_u_;
    NewtonSolver newton_d_;

    // Linear solvers
    CGSolver cg_u_, cg_d_;

    // Material
    const PFFMaterialParameters &mat_;

    // Solver parameters
    int max_fp_it_ = 20;
    real_t fp_rel_tol_ = 1e-8;
    real_t fp_abs_tol_ = 1e-10;

    // Time
    real_t t_ = 0.0;

    // Internal methods
    void SolveElasticity();
    void SolveDamage();
    void UpdateStrainEnergyHistory();
    bool CheckConvergence(const Vector &d_new, const Vector &d_old) const;
};
```

**Unit Test**: `tests/test_staggered_solver.cpp`
```cpp
TEST(PFFSolver, FixedPointConvergence) {
    // Solver should converge for moderate loading
    // ...
}

TEST(PFFSolver, DamageIrreversible) {
    // Damage should never decrease across time steps
    // ...
}

TEST(PFFSolver, EnergyDecreasing) {
    // Total energy should decrease during damage evolution
    // (for quasi-static loading at fixed displacement)
    // ...
}

TEST(PFFSolver, NoLoadNoDamage) {
    // Without loading, no damage should evolve
    // ...
}
```

---

## Phase 5: Main Driver and Examples

### Task 5.1: Main Driver
**File**: `pff.cpp`

```cpp
int main(int argc, char *argv[]) {
    // Initialize MPI
    Mpi::Init(argc, argv);

    // Parse command line options
    OptionsParser args(argc, argv);
    // ... mesh file, order, material params, time stepping

    // Load/create mesh
    Mesh mesh(mesh_file);
    ParMesh pmesh(MPI_COMM_WORLD, mesh);

    // Material parameters
    PFFMaterialParameters mat(E, nu, l, sigma_t);

    // Create solver
    PFFSolver solver(pmesh, mat, order);

    // Setup boundary conditions
    solver.SetElasticityEssentialBdr(ess_bdr);
    solver.SetPrescribedDisplacement(u_bc);

    // Time stepping loop
    real_t t = 0.0;
    while (t < t_final) {
        // Solve
        int fp_its = solver.SolveTimeStep(dt);

        // Output
        if (step % vis_steps == 0) {
            // ParaView output
        }

        // Advance time
        solver.AdvanceTime(dt);
        t += dt;
    }

    return 0;
}
```

---

### Task 5.2: Three-Point Bending Example
**File**: `examples/three_point_bending.cpp`

Complete example matching the MOOSE three-point bending test:
- 3D beam geometry
- Top loading (prescribed displacement)
- Bottom supports (fixed)
- Track force-displacement curve
- Visualize damage evolution

**Unit Test**: `tests/test_three_point_bending.cpp`
```cpp
TEST(ThreePointBending, InitialStiffness) {
    // Linear response before damage
    // ...
}

TEST(ThreePointBending, DamageInitiation) {
    // Damage should initiate at critical load
    // ...
}

TEST(ThreePointBending, SymmetricDamage) {
    // Damage pattern should be symmetric for symmetric geometry
    // ...
}
```

---

## Phase 6: Testing Infrastructure

### Task 6.1: Test Framework Setup
**File**: `tests/CMakeLists.txt`

Setup Google Test integration:
```cmake
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Unit tests
add_executable(pff_tests
    test_degradation_function.cpp
    test_spectral_decomposition.cpp
    test_pff_material.cpp
    test_degraded_elasticity_integrator.cpp
    test_damage_diffusion_integrator.cpp
    test_damage_source_integrator.cpp
    test_elasticity_operator.cpp
    test_damage_operator.cpp
    test_staggered_solver.cpp
    test_three_point_bending.cpp
)

target_link_libraries(pff_tests
    pff_lib
    mfem
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(pff_tests)
```

---

## Implementation Order and Dependencies

```
Phase 1: Material Models (no MFEM dependencies)
├── 1.1 Degradation Function ──┐
├── 1.2 Spectral Decomposition ├── Independent, can parallelize
└── 1.3 Material Parameters ───┘

Phase 2: Integrators (depends on Phase 1)
├── 2.1 Degraded Elasticity Integrator ─── depends on 1.1, 1.2, 1.3
├── 2.2 Damage Diffusion Integrator ────── depends on 1.3
└── 2.3 Damage Source Integrator ───────── depends on 1.1, 1.3

Phase 3: Operators (depends on Phase 2)
├── 3.1 Elasticity Operator ────── depends on 2.1
└── 3.2 Damage Operator ────────── depends on 2.2, 2.3

Phase 4: Staggered Solver (depends on Phase 3)
└── 4.1 PFF Solver ─────────────── depends on 3.1, 3.2

Phase 5: Driver and Examples (depends on Phase 4)
├── 5.1 Main Driver ────────────── depends on 4.1
└── 5.2 Three-Point Bending ────── depends on 5.1

Phase 6: Testing (parallel with all phases)
└── 6.1 Test Framework ─────────── independent
```

---

## Verification and Validation

### Verification Tests
1. **Patch test**: Linear displacement gives zero internal forces
2. **Convergence study**: Mesh refinement convergence rates
3. **Symmetry**: Symmetric loading gives symmetric response
4. **Energy conservation**: Total energy balance

### Validation Against MOOSE
1. **Force-displacement curve**: Match MOOSE results
2. **Damage pattern**: Compare spatial distribution
3. **Critical load**: Match damage initiation threshold
4. **Post-peak response**: Match softening behavior

---

## Performance Considerations

1. **Parallel scalability**: All operators should work with ParMesh
2. **Linear solver**: Use algebraic multigrid (BoomerAMG via HYPRE)
3. **Memory**: Reuse vectors/matrices where possible
4. **Quadrature**: Use appropriate quadrature order for nonlinear terms

---

## References

1. Miehe, C., Welschinger, F., & Hofacker, M. (2010). Thermodynamically consistent phase-field models of fracture: Variational principles and multi-field FE implementations. IJNME, 83(10), 1273-1311.
2. Bourdin, B., Francfort, G. A., & Marigo, J. J. (2000). Numerical experiments in revisited brittle fracture. JMPS, 48(4), 797-826.
3. MFEM Documentation: https://mfem.org/
4. MOOSE Phase Field Module: https://mooseframework.inl.gov/modules/phase_field/
