# Phase 7: Parallel Domain Operator

## Overview

This phase extends the domain operator to support parallel (MPI) execution using ParMesh, ParFiniteElementSpace, and HYPRE solvers.

## Dependencies

- **Phase 2** (`phase2_domain_operator.md`): Serial domain operator implementation
- **Phase 6** (`phase6_parallel_infrastructure.md`): Type abstractions and MPI context

## Reference

- **dg_antiplane_theory.md**: DG formulation (same math, different linear algebra)
- **MFEM Parallel Examples**: `ex1p.cpp`, `ex14p.cpp`

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 7.1 Template domain operator | `domain_operator.hpp` (templated) | - |
| 7.2 Parallel antiplane impl. | `antiplane_operator_impl.hpp` | `test_parallel_domain.cpp` |
| 7.3 HYPRE solver setup | In antiplane impl. | - |
| 7.4 Parallel MMS test | - | `mms_antiplane_parallel.cpp` |

## Verification Checkpoint

Parallel MMS shows same convergence as serial.

---

## Parallel Domain Design

### Template-Based Approach

The domain operator uses C++ templates to support both serial and parallel MFEM types:

```cpp
template <typename MeshType = SEASMesh>
class DomainOperator;

// Specializations:
// - DomainOperator<Mesh>     -> Serial
// - DomainOperator<ParMesh>  -> Parallel
```

### Key Differences: Serial vs Parallel

| Aspect | Serial | Parallel |
|--------|--------|----------|
| Mesh | `Mesh` | `ParMesh` |
| FE Space | `FiniteElementSpace` | `ParFiniteElementSpace` |
| DOFs | All DOFs local | True DOFs distributed |
| Bilinear Form | `BilinearForm` | `ParBilinearForm` |
| Linear Algebra | Dense/Sparse | HypreParMatrix, HypreParVector |
| Solver | CGSolver | HyprePCG |
| Preconditioner | GSSmoother | HypreBoomerAMG |

---

## Detailed Component Design

### 7.1 Templated Domain Operator

```cpp
// domain/domain_operator.hpp

template <typename MeshType = SEASMesh>
class DomainOperator {
public:
    // Type aliases derived from mesh type
    using FESpaceType = typename std::conditional<
        std::is_same<MeshType, ParMesh>::value,
        ParFiniteElementSpace, FiniteElementSpace>::type;
    using GridFuncType = typename std::conditional<
        std::is_same<MeshType, ParMesh>::value,
        ParGridFunction, GridFunction>::type;
    using BilinFormType = typename std::conditional<
        std::is_same<MeshType, ParMesh>::value,
        ParBilinearForm, BilinearForm>::type;

    virtual ~DomainOperator() = default;

    virtual int NumComponents() const = 0;
    virtual void Solve(real_t time, const Vector &slip_bc,
                       GridFuncType &displacement) = 0;
    virtual void ComputeTraction(const GridFuncType &displacement,
                                  Vector &traction) = 0;
    virtual FESpaceType &GetFESpace() = 0;
    virtual MeshType &GetMesh() = 0;

#ifdef SEAS_USE_MPI
    virtual MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#endif
};
```

### 7.2 Parallel Antiplane Operator

```cpp
// domain/antiplane_operator_impl.hpp

// Parallel-specific solver setup
template <>
void AntiplaneDomainOperator<ParMesh>::SetupSolver() {
    // Parallel solver setup using HYPRE
    auto *amg = new HypreBoomerAMG();
    amg->SetPrintLevel(0);

    auto *pcg = new HyprePCG(GetComm());
    pcg->SetTol(1e-12);
    pcg->SetMaxIter(500);
    pcg->SetPrintLevel(0);
    pcg->SetPreconditioner(*amg);

    prec_.reset(amg);
    solver_.reset(pcg);
}

// Serial-specific solver setup
template <>
void AntiplaneDomainOperator<Mesh>::SetupSolver() {
    auto *gs = new GSSmoother();
    auto *cg = new CGSolver();
    cg->SetRelTol(1e-12);
    cg->SetMaxIter(500);
    cg->SetPrintLevel(0);
    cg->SetPreconditioner(*gs);

    prec_.reset(gs);
    solver_.reset(cg);
}
```

### 7.3 Parallel Solve Implementation

```cpp
template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::Solve(
    real_t time, const Vector &slip_bc, GridFuncType &displacement)
{
    // Apply boundary conditions (works for both serial and parallel)
    ApplyFarFieldBC(time, displacement);
    ApplyFaultSlipBC(slip_bc, displacement);

    // Get true DOFs
    // For ParFiniteElementSpace, this extracts the owned true DOFs
    // For FiniteElementSpace, this is the full DOF vector
    fes_->GetEssentialTrueDofs(ess_bdr_, ess_tdof_list_);

    // Form linear system
    // ParBilinearForm creates HypreParMatrix automatically
    // BilinearForm creates SparseMatrix
    a_form_->FormLinearSystem(ess_tdof_list_, displacement, *b_, A_, X_, B_);

    // Solve
    solver_->SetOperator(*A_);
    solver_->Mult(B_, X_);

    // Recover solution
    a_form_->RecoverFEMSolution(X_, *b_, displacement);
}
```

### Parallel Mesh Distribution

```cpp
// Mesh creation and distribution
void CreateParallelMesh(const BP2Params &params, int refinement_level) {
    // Create serial mesh on rank 0
    Mesh serial_mesh = Mesh::MakeCartesian2D(
        params.nx, params.nz,
        Element::QUADRILATERAL,
        true,  // generate_edges
        2.0 * params.Lx, params.Lz  // domain size
    );

    // Shift to center domain at x=0
    // ...

    // Distribute to all ranks
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

    // Parallel mesh now distributes elements across ranks
    // Ghost elements are automatically created at partition boundaries
}
```

---

## Parallel Unit Tests

### 7.2 Parallel Domain Tests

```cpp
// tests/parallel/test_parallel_domain.cpp

class ParallelDomainTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
    }
    int rank_, size_;
};

/// Test 1: Parallel mesh distribution
TEST_F(ParallelDomainTest, MeshDistribution) {
    // Create serial mesh
    Mesh serial_mesh = Mesh::MakeCartesian2D(64, 64, Element::QUADRILATERAL);

    // Distribute to parallel
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

    // Each rank should have portion of elements
    int local_ne = pmesh.GetNE();
    int global_ne;
    MPI_Allreduce(&local_ne, &global_ne, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(global_ne, 64 * 64);
    EXPECT_GT(local_ne, 0);  // Each rank has some elements
}

/// Test 2: Parallel solve produces same result as serial
TEST_F(ParallelDomainTest, ParallelSerialConsistency) {
    BP2Params params;
    real_t mu = 32.04e9;
    int N = 32;

    // Serial solve (on rank 0 only)
    Vector serial_solution;
    if (rank_ == 0) {
        Mesh serial_mesh = CreateBP2Mesh(params, N);
        AntiplaneDomainOperator<Mesh> serial_domain(serial_mesh, 1, mu);

        Vector slip(serial_mesh.GetNBE());
        slip = 0.0;  // Zero slip for comparison

        GridFunction u(&serial_domain.GetFESpace());
        serial_domain.Solve(0.0, slip, u);
        serial_solution = u;
    }

    // Parallel solve
    Mesh serial_mesh_for_par = CreateBP2Mesh(params, N);
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh_for_par);
    AntiplaneDomainOperator<ParMesh> par_domain(pmesh, 1, mu);

    // Setup and solve
    // ...

    // Gather parallel solution to rank 0 and compare
    if (rank_ == 0) {
        for (int i = 0; i < serial_solution.Size(); i++) {
            EXPECT_NEAR(serial_solution(i), gathered_solution(i),
                        1e-10 * std::abs(serial_solution(i)));
        }
    }
}

/// Test 3: DOF distribution is correct
TEST_F(ParallelDomainTest, DOFDistribution) {
    Mesh serial_mesh = Mesh::MakeCartesian2D(32, 32, Element::QUADRILATERAL);
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

    FiniteElementCollection *fec = new H1_FECollection(1, 2);
    ParFiniteElementSpace pfes(&pmesh, fec);

    // Global true DOFs should match serial case
    int local_tdof = pfes.GetTrueVSize();
    int global_tdof;
    MPI_Allreduce(&local_tdof, &global_tdof, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // Account for shared DOFs (approximately)
    int expected_serial_dof = 33 * 33;  // (N+1)^2 for linear elements
    EXPECT_LE(global_tdof, expected_serial_dof);
    EXPECT_GE(global_tdof, expected_serial_dof * 0.9);

    delete fec;
}

/// Test 4: Boundary conditions in parallel
TEST_F(ParallelDomainTest, BoundaryConditions) {
    // Verify BCs are correctly applied on partitioned boundaries
    // Each rank applies BCs only on its local boundary faces
}
```

### 7.4 Parallel MMS Verification

```cpp
// tests/verification/mms_antiplane_parallel.cpp

/// Parallel MMS test: verify convergence in parallel execution
TEST(ParallelMMS, AntiplaneConvergence) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    real_t mu = 1.0;  // Simplified for MMS

    std::vector<int> N_values = {16, 32, 64, 128};
    std::vector<real_t> errors;

    for (int N : N_values) {
        // Create serial mesh, then distribute
        Mesh serial_mesh = Mesh::MakeCartesian2D(N, N, Element::QUADRILATERAL);
        ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

        // Setup MMS: u_exact = sin(pi*x) * sin(pi*z)
        // Source: f = 2*pi^2 * sin(pi*x) * sin(pi*z)
        auto u_exact = [](const Vector &x) {
            return sin(M_PI * x(0)) * sin(M_PI * x(1));
        };
        auto f_source = [](const Vector &x) {
            return 2.0 * M_PI * M_PI * sin(M_PI * x(0)) * sin(M_PI * x(1));
        };

        // Create operator and solve
        AntiplaneDomainOperator<ParMesh> domain(pmesh, 1, mu);
        // ... add source term, solve ...

        // Compute local L2 error
        real_t local_L2_sq = 0.0;
        // ... integrate (u - u_exact)^2 over local elements ...

        // Global L2 error
        real_t global_L2_sq;
        MPI_Allreduce(&local_L2_sq, &global_L2_sq, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        real_t global_L2 = std::sqrt(global_L2_sq);

        if (rank == 0) {
            errors.push_back(global_L2);
        }
    }

    // Verify convergence order (on root)
    if (rank == 0) {
        for (size_t i = 1; i < errors.size(); i++) {
            real_t ratio = errors[i-1] / errors[i];
            EXPECT_GT(ratio, 3.5);  // Second order: ratio ~ 4
        }
    }
}
```

---

## Key Implementation Patterns

### Parallel Assembly

```cpp
// ParBilinearForm handles parallel assembly internally
template <>
void AntiplaneDomainOperator<ParMesh>::SetupBilinearForm() {
    a_form_.reset(new ParBilinearForm(fes_.get()));
    a_form_->AddDomainIntegrator(new DiffusionIntegrator());
    a_form_->Assemble();
    // ParBilinearForm::FormSystemMatrix creates HypreParMatrix
}
```

### True DOFs vs Local DOFs

```cpp
// In parallel, we work with true DOFs (owned by this rank)
void ComputeTraction(const ParGridFunction &u, Vector &traction) {
    // Get true DOF values
    const Vector &u_true = u.GetTrueVector();

    // For fault DOFs on this rank, compute traction
    // ...
}
```

### Ghost Exchange

```cpp
// MFEM handles ghost exchange automatically for GridFunctions
ParGridFunction u(&pfes);
u.Distribute();  // Updates ghost values from neighboring ranks
```

---

## File Organization

```
domain/
├── domain_operator.hpp           # Abstract interface (templated)
├── antiplane_operator.hpp        # Antiplane declaration (templated)
├── antiplane_operator_impl.hpp   # Implementation (serial/parallel specific)
└── plane_strain_operator.hpp     # [Future] 2D plane strain
```

---

## Acceptance Criteria

| Test | Pass Criteria |
|------|---------------|
| Mesh distribution | All ranks have elements, global count preserved |
| DOF distribution | Unique DOF ownership, no gaps or duplicates |
| Parallel solver | Converges to same tolerance as serial |
| Serial-parallel match | Parallel solution matches serial within 1e-10 |
| Parallel MMS | Shows same 2nd order convergence as serial |
