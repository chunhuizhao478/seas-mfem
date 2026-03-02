# Phase 10: Extensions and Optimization

## Overview

This phase extends the implementation to support additional problem types, implements optimizations, and prepares for large-scale simulations.

## Dependencies

- **Phase 1-9**: Complete serial and parallel BP2 implementation

## Reference

- **SCEC BP1 Specification**: https://strike.scec.org/cvws/seas/download/SEAS_BP1_QD.pdf
- **SCEC BP3 Specification**: For 3D extension reference

## Tasks and Deliverables

| Task | Deliverable | Tests |
|------|-------------|-------|
| 10.1 Plane strain operator | `plane_strain_operator.hpp` | - |
| 10.2 BP1 configuration | `bp1_params.hpp` | - |
| 10.3 Documentation | README, examples | - |
| 10.4 Performance optimization | Profiling, tuning | Timing tests |
| 10.5 Large-scale testing | - | 3D prototype tests |

---

## Extensions

### 10.1 Plane Strain Operator (2D Vector Problem)

Extends the antiplane scalar problem to 2D plane strain elasticity (vector problem):

```cpp
// domain/plane_strain_operator.hpp

namespace mfem {
namespace seas {

/// Solves 2D plane strain elasticity
/// Displacement u = (u_x, u_y) is a vector field
template <typename MeshType = Mesh>
class PlaneStrainDomainOperator : public DomainOperator<MeshType> {
public:
    using Base = DomainOperator<MeshType>;
    using FESpaceType = typename Base::FESpaceType;
    using GridFuncType = typename Base::GridFuncType;

    PlaneStrainDomainOperator(MeshType &mesh, int order,
                               real_t lambda, real_t mu);

    int NumComponents() const override { return 2; }  // Vector field

    void Solve(real_t time, const Vector &slip_bc,
               GridFuncType &displacement) override;

    void ComputeTraction(const GridFuncType &displacement,
                          Vector &traction) override;

private:
    real_t lambda_;  // First Lame parameter
    real_t mu_;      // Shear modulus

    // Vector FE space for 2D displacement
    std::unique_ptr<FiniteElementCollection> fec_;
    std::unique_ptr<FESpaceType> fes_;  // dim=2 components

    // Elasticity bilinear form
    std::unique_ptr<BilinearForm> a_form_;
};

} // namespace seas
} // namespace mfem
```

### Plane Strain Equations

**Equilibrium** (quasi-static):
```
div(sigma) = 0

sigma_xx = lambda * (eps_xx + eps_yy) + 2*mu*eps_xx
sigma_yy = lambda * (eps_xx + eps_yy) + 2*mu*eps_yy
sigma_xy = 2*mu*eps_xy
```

**Strain-displacement**:
```
eps_xx = du_x/dx
eps_yy = du_y/dy
eps_xy = 0.5 * (du_x/dy + du_y/dx)
```

### 10.2 BP1 Configuration

BP1 is a 1D slip on a 2D domain (simpler than BP2):

```cpp
// config/bp1_params.hpp

namespace mfem {
namespace seas {

struct BP1Params {
    // Material properties (same as BP2)
    real_t rho = 2670.0;
    real_t cs = 3464.0;
    real_t mu = rho * cs * cs;

    // Friction parameters
    real_t sigma_n = 50.0e6;
    real_t a = 0.015;    // Uniform a (unlike BP2)
    real_t b = 0.019;    // Velocity weakening (a < b)
    real_t Dc = 0.008;   // Different from BP2
    real_t V0 = 1.0e-6;
    real_t f0 = 0.6;

    // Loading
    real_t Vp = 1.0e-9;
    real_t V_init = 1.0e-9;

    // Geometry (1D fault in 2D domain)
    real_t L = 100.0e3;   // Fault length
    real_t W = 100.0e3;   // Domain width

    // Simulation
    real_t t_final = 200.0 * 3.15e7;  // 200 years
};

} // namespace seas
} // namespace mfem
```

### Differences from BP2

| Parameter | BP1 | BP2 |
|-----------|-----|-----|
| Dimension | 2D domain, 1D fault | 2D domain, 1D fault |
| Rate-state a | Uniform (0.015) | Depth-dependent (0.010-0.025) |
| Rate-state b | 0.019 | 0.015 |
| Dc | 0.008 m | 0.004 m |
| Fault model | Spring-slider | Continuum |

---

## Performance Optimization

### 10.4 Optimization Strategies

#### 1. Solver Optimization

```cpp
// Optimize HYPRE solver parameters for SEAS problems
void OptimizeSolver(HypreBoomerAMG &amg, HyprePCG &pcg) {
    // AMG settings for elliptic problems
    amg.SetCoarsening(6);  // Falgout coarsening
    amg.SetInterpolation(0);  // Classical interpolation
    amg.SetRelaxType(6);  // Symmetric SOR
    amg.SetNumSweeps(1);
    amg.SetStrongThreshold(0.25);

    // PCG settings
    pcg.SetTol(1e-12);
    pcg.SetMaxIter(200);
    pcg.SetPrintLevel(0);
}
```

#### 2. Matrix Reuse

The stiffness matrix doesn't change during time stepping:

```cpp
// Reuse assembled matrix across time steps
class OptimizedAntiplaneDomain {
    void SetupOperator() {
        a_form_->Assemble();
        a_form_->Finalize();

        // Store LU factorization or AMG setup
        A_assembled_ = true;
    }

    void Solve(...) {
        if (!A_assembled_) {
            SetupOperator();
        }

        // Only update RHS and BCs, not matrix
        UpdateRHS(time, slip_bc);
        solver_->Mult(B_, X_);
    }
};
```

#### 3. State Variable Optimization

```cpp
// Cache depth-dependent parameters
class OptimizedRateStateFault {
    void Initialize() {
        // Precompute a(z), sigma_n, eta for each DOF
        a_values_.SetSize(num_fault_dofs_);
        for (int i = 0; i < num_fault_dofs_; i++) {
            real_t z = depths_(i);
            a_values_(i) = ComputeA(z);
        }

        // Precompute tau0 once
        tau0_ = ComputePreStress();
    }
};
```

#### 4. Vectorized Friction Evaluation

```cpp
// Vectorize friction coefficient computation
void ComputeFrictionVectorized(const Vector &V, const Vector &theta,
                                 const Vector &a, Vector &f) {
    int n = V.Size();
    for (int i = 0; i < n; i++) {
        real_t arg = (V(i) / (2.0 * V0_)) *
                     std::exp((f0_ + b_ * std::log(V0_ * theta(i) / Dc_)) / a(i));
        f(i) = a(i) * std::asinh(arg);
    }
}
```

### Profiling Results Template

| Component | Time (%) | Optimization Opportunity |
|-----------|----------|-------------------------|
| Domain solve | 60% | AMG tuning, matrix reuse |
| Slip rate solve | 25% | Vectorization, Newton tuning |
| State evolution | 5% | Vectorization |
| I/O | 5% | Async output, buffering |
| MPI | 5% | Overlap computation/communication |

---

## Large-Scale Testing

### 10.5 3D Prototype

Preparing for future 3D extension:

```cpp
// domain/elastic3d_operator.hpp (prototype)

template <typename MeshType = Mesh>
class Elastic3DDomainOperator : public DomainOperator<MeshType> {
public:
    int NumComponents() const override { return 3; }  // 3D vector field

    // 3D elasticity solve
    void Solve(real_t time, const Vector &slip_bc,
               GridFuncType &displacement) override;

    // Traction on fault surface (2D manifold in 3D)
    void ComputeTraction(const GridFuncType &displacement,
                          Vector &traction) override;
};
```

### Large-Scale Test Cases

| Test Case | Mesh Size | DOFs | Target Ranks |
|-----------|-----------|------|--------------|
| BP2 fine | 1024x1024 | ~1M | 64-256 |
| 3D prototype | 64x64x64 | ~1M | 64-256 |
| Production 3D | 256x256x256 | ~50M | 1000+ |

### Memory Scaling

```cpp
// Estimate memory requirements
size_t EstimateMemory(int N, int dim, int order) {
    size_t elements = std::pow(N, dim);
    size_t dofs = std::pow(N + 1, dim) * (dim == 3 ? 3 : (dim == 2 ? 2 : 1));

    // Matrix storage (sparse)
    size_t nnz_per_row = std::pow(2 * order + 1, dim);
    size_t matrix_mem = dofs * nnz_per_row * sizeof(double);

    // Vectors (solution, RHS, work vectors)
    size_t vector_mem = dofs * sizeof(double) * 10;

    return matrix_mem + vector_mem;
}
```

---

## Documentation

### 10.3 Documentation Deliverables

#### README.md

```markdown
# SEAS Miniapp

SEAS (Sequences of Earthquakes and Aseismic Slip) simulation using MFEM.

## Features

- BP1 and BP2 benchmark implementations
- Serial and parallel (MPI) execution
- Rate-and-state friction with aging law
- Adaptive time stepping

## Building

### Serial
```bash
cmake .. -DMFEM_USE_MPI=OFF
make seas
```

### Parallel
```bash
cmake .. -DMFEM_USE_MPI=ON
make pseas
```

## Running

### Serial BP2
```bash
./seas -bp bp2 -tf 1e10 -o output/
```

### Parallel BP2 (4 ranks)
```bash
mpirun -np 4 ./pseas -bp bp2 -tf 1e10 -o output/
```

## Examples

See `examples/` directory for complete examples:
- `bp2_qd.cpp`: BP2 quasi-dynamic example
- `bp2_qd_parallel.cpp`: Parallel BP2 example
- `bp1_qd.cpp`: BP1 for comparison
```

#### Example Files

```cpp
// examples/bp2_qd.cpp

/// Example: Run BP2-QD benchmark (serial)
int main(int argc, char *argv[]) {
    // Parse options
    int order = 1;
    real_t t_final = 1e10;
    std::string output_dir = "bp2_output";

    OptionsParser args(argc, argv);
    args.AddOption(&order, "-o", "--order", "Polynomial order");
    args.AddOption(&t_final, "-tf", "--t-final", "Final time");
    args.AddOption(&output_dir, "-out", "--output", "Output directory");
    args.Parse();

    // Setup BP2
    BP2Params params;
    params.t_final = t_final;

    // Create mesh
    Mesh mesh = CreateBP2Mesh(params);

    // Create operators
    AntiplaneDomainOperator<Mesh> domain(mesh, order, params.mu);
    FaultGeometry geom(mesh, 1, params);
    DieterichRuinaFriction friction({params.V0, params.f0, params.b0, params.Dc});
    AgingLaw aging;
    RateStateFaultOperator fault(&geom, &friction, &aging, params);

    // Create SEAS operator
    SEASQuasiDynamicOperator<Mesh> seas_op(&domain, &fault);

    // Time stepping
    AdaptiveTimeStepper stepper;
    RK4Solver ode;
    ode.Init(seas_op);

    Vector state(fault.StateSize());
    seas_op.SetInitialCondition(state);

    // Run
    BenchmarkOutput output(output_dir, params);
    real_t t = 0.0;
    while (t < params.t_final) {
        real_t dt = stepper.GetDt();
        ode.Step(state, t, dt);
        t += dt;

        stepper.SetDt(stepper.ComputeNewDt(seas_op.GetMaxSlipRate(), dt));
        output.Write(t, state, seas_op);
    }

    return 0;
}
```

---

## Future Extensions

### Planned Enhancements

| Extension | Description | Priority |
|-----------|-------------|----------|
| Slip law | Alternative state evolution | Medium |
| 3D elasticity | Full 3D earthquake cycles | High |
| GPU acceleration | CUDA/HIP solvers | Low |
| Dynamic rupture | Inertial effects | Medium |
| Heterogeneous friction | Spatially varying a, b | Low |
| Thermal pressurization | Coupled thermo-poro-mechanics | Low |

### Architecture for Extensions

```cpp
// Future-proof class hierarchy

class DomainOperator { /* abstract */ };
  class AntiplaneDomainOperator : public DomainOperator { /* 2D scalar */ };
  class PlaneStrainDomainOperator : public DomainOperator { /* 2D vector */ };
  class Elastic3DDomainOperator : public DomainOperator { /* 3D vector */ };

class FrictionLaw { /* abstract */ };
  class DieterichRuinaFriction : public FrictionLaw { /* regularized RS */ };
  class ClassicalRSFriction : public FrictionLaw { /* classical */ };

class StateEvolution { /* abstract */ };
  class AgingLaw : public StateEvolution { /* d(theta)/dt = 1 - V*theta/Dc */ };
  class SlipLaw : public StateEvolution { /* d(theta)/dt = -V*theta/Dc * ln(...) */ };
  class CompositeEvolution : public StateEvolution { /* combined laws */ };
```

---

## File Organization

```
miniapps/seas/
├── CMakeLists.txt
├── seas.cpp                    # Serial driver
├── pseas.cpp                   # Parallel driver
│
├── common/
│   ├── seas_types.hpp
│   ├── mpi_context.hpp
│   └── parallel_utils.hpp
│
├── config/
│   ├── seas_config.hpp
│   ├── bp1_params.hpp
│   └── bp2_params.hpp
│
├── domain/
│   ├── domain_operator.hpp
│   ├── antiplane_operator.hpp
│   ├── plane_strain_operator.hpp
│   └── elastic3d_operator.hpp    # [future]
│
├── fault/
│   ├── fault_operator.hpp
│   ├── fault_geometry.hpp
│   ├── rate_state_fault.hpp
│   └── parallel_fault_data.hpp
│
├── friction/
│   ├── friction_law.hpp
│   ├── dieterich_ruina.hpp
│   └── state_evolution.hpp
│
├── solver/
│   ├── seas_operator.hpp
│   ├── quasi_dynamic.hpp
│   ├── time_stepper.hpp
│   └── parallel_time_stepper.hpp
│
├── io/
│   ├── benchmark_output.hpp
│   ├── probe_output.hpp
│   ├── parallel_probe_output.hpp
│   └── checkpoint.hpp
│
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── parallel/
│   └── verification/
│
├── examples/
│   ├── bp1_qd.cpp
│   ├── bp2_qd.cpp
│   └── bp2_qd_parallel.cpp
│
└── document/
    ├── bp2_implementation_plan.md
    ├── phase1_core_infrastructure.md
    ├── phase2_domain_operator.md
    └── ... (other phase documents)
```

---

## Acceptance Criteria

| Component | Criteria |
|-----------|----------|
| Plane strain | Passes MMS verification |
| BP1 | Runs and produces reasonable results |
| Documentation | README complete with examples |
| Performance | Meets scaling targets from Phase 9 |
| Code quality | Clean, well-documented, follows MFEM style |
