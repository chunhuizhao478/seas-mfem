# BP2 Benchmark Implementation Plan for MFEM

## Executive Summary

This document provides a detailed implementation plan for the SCEC SEAS BP2-QD (Benchmark Problem 2, Quasi-Dynamic) benchmark in MFEM. The design emphasizes:
1. **Accurate reproduction** of BP2 benchmark results
2. **Generality** for future extension to 2D plane strain and 3D problems
3. **Comprehensive unit testing** to verify each component
4. **Serial/Parallel duality** - code ready for MPI parallelism from the start, with serial version for concept verification

---

# Part I: BP2 Benchmark Specification

## 1. Problem Definition

### 1.1 Physical Setup

BP2 is a **2D antiplane shear** problem with the following characteristics:

```
Domain: (x, z) ∈ (-∞, ∞) × (0, ∞)  [half-space]

Coordinate system:
  - x: horizontal (perpendicular to fault)
  - y: along-strike (displacement direction)
  - z: depth (positive downward)

Displacement: u = u(x, z, t) in y-direction only (scalar field)
```

### 1.2 Governing Equations

**Equilibrium equation** (quasi-static):
```
0 = ∂σxy/∂x + ∂σyz/∂z
```

**Hooke's law** (antiplane shear):
```
σxy = μ ∂u/∂x
σyz = μ ∂u/∂z
```

This reduces to the **Laplace equation**:
```
∇²u = ∂²u/∂x² + ∂²u/∂z² = 0
```

### 1.3 Computational Domain (Full Domain)

The physical half-space domain is truncated to a **full domain** for computation:

```
Computational Domain: x ∈ [-Lx, +Lx], z ∈ [0, Lz]

                        z = 0 (free surface)
        ←───────────────────────────────────────────→
        │                    │                      │
        │   Left half        │    Right half        │
        │   (x < 0)          │    (x > 0)           │
        │                    │                      │
x = -Lx │      u → -Vp·t/2   │0   u → +Vp·t/2      │ x = +Lx
        │                    │                      │
        │                  FAULT                    │
        │                (slip = δ)                 │
        │                    │                      │
        ←───────────────────────────────────────────→
                        z = Lz (bottom)
```

**Full Domain Formulation:**
- Domain: x ∈ [-Lx, +Lx], z ∈ [0, Lz] with Lx = Lz = 100 km
- Fault at x = 0 (internal interface with prescribed slip jump)
- Opposite displacements on each side of the fault:
  - Left side (x = 0⁻): u = -δ/2
  - Right side (x = 0⁺): u = +δ/2
- Total slip across fault: u(0⁺) - u(0⁻) = δ

### 1.4 Boundary Conditions (Full Domain)

| Boundary | Location | Condition | Mathematical Form |
|----------|----------|-----------|-------------------|
| Fault (left) | x = 0⁻, z < Wf | Dirichlet | u = -δ/2 |
| Fault (right) | x = 0⁺, z < Wf | Dirichlet | u = +δ/2 |
| Creep (left) | x = 0⁻, z ≥ Wf | Natural | σxy = 0 (free slip) |
| Creep (right) | x = 0⁺, z ≥ Wf | Natural | σxy = 0 (free slip) |
| Far-field left | x = -Lx | Natural | σxy = 0 |
| Far-field right | x = +Lx | Natural | σxy = 0 |
| Free surface | z = 0 | Natural | σyz = 0 |
| Bottom | z = Lz | Dirichlet | u = sign(x)·Vp·t/2 |

**Note on Boundary Conditions:**
- Following Tandem's approach, plate loading is applied at the **bottom boundary** (z = Lz)
- The displacement at bottom varies with x-position: positive on right half, negative on left half
- Far-field boundaries (x = ±Lx) have zero traction (Natural BC), not Dirichlet
- See `tandem_boundary_condition_analysis.md` for detailed analysis

### 1.5 Friction Law

**Regularized rate-and-state friction**:
```
f(V, θ) = a · sinh⁻¹[(V / 2V₀) · exp((f₀ + b·ln(V₀θ/Dc)) / a)]

Fault strength:
F(V, θ) = σn · f(V, θ)

Stress balance with radiation damping:
τ = τ⁰ + τqs - η·V = F(V, θ)

where:
  τ⁰ = initial/pre-stress
  τqs = shear stress from quasi-static deformation
  η = μ/(2cs) = radiation damping coefficient
```

**Aging law state evolution**:
```
dθ/dt = 1 - V·θ/Dc
```

### 1.5 BP2 Parameters

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Density | ρ | 2670 | kg/m³ |
| Shear wave speed | cs | 3464 | m/s |
| Shear modulus | μ | ρ·cs² = 32.04 | GPa |
| Normal stress | σn | 50 | MPa |
| Rate-state a (VW zone) | a₀ | 0.010 | - |
| Rate-state a (VS zone) | amax | 0.025 | - |
| Rate-state b | b₀ | 0.015 | - |
| **Critical slip distance** | **Dc** | **0.004** | **m** |
| Plate rate | Vp | 10⁻⁹ | m/s |
| Initial slip rate | Vinit | 10⁻⁹ | m/s |
| Reference slip rate | V₀ | 10⁻⁶ | m/s |
| Reference friction | f₀ | 0.6 | - |
| VW zone depth | H | 15 | km |
| VW-VS transition width | h | 3 | km |
| Rate-state fault width | Wf | 40 | km |
| Final simulation time | tf | 1200 | years |

### 1.6 Depth-Dependent Parameter a(z)

```
a(z) = {
    a₀,                                    0 ≤ z < H
    a₀ + (amax - a₀)(z - H)/h,            H ≤ z < H + h
    amax,                                  H + h ≤ z < Wf
}

where H = 15 km, h = 3 km, a₀ = 0.010, amax = 0.025
```

### 1.7 Initial Conditions

**Initial slip**:
```
δ(z, 0) = 0
```

**Pre-stress** (uniform, at steady-state with Vinit at depth Wf):
```
τ⁰ = σn · amax · sinh⁻¹[(Vinit / 2V₀) · exp((f₀ + b₀·ln(V₀/Vinit)) / amax)] + η·Vinit
```

**Initial state** (depth-dependent, consistent with τ⁰ and Vinit):
```
θ(z, 0) = (Dc/V₀) · exp{(a/b)·ln[(2V₀/Vinit)·sinh((τ⁰ - η·Vinit)/(a·σn))] - f₀/b}
```

### 1.8 Output Requirements

**Time series at 12 stations** (depths in km):
- z = 0.0, 2.4, 4.8, 7.2, 9.6, 12.0, 14.4, 16.8, 19.2, 24.0, 28.8, 36.0

**Output quantities**:
1. Time (s)
2. Slip δ (m)
3. Slip rate V (log₁₀ m/s)
4. Shear stress τ (MPa)
5. State θ (log₁₀ s)

**Available benchmark data** (depths in km):
- z = 0, 4.8, 12, 16.8, 24

---

# Part II: Architecture Design

## 2. Design Principles

### 2.1 Serial/Parallel Design Philosophy

**Key Principle**: Design for MPI parallelism from the start, but implement serial version first for concept verification.

**Implementation Strategy**:
1. **Template-based abstraction**: Use template parameters or abstract interfaces to switch between serial and parallel MFEM classes
2. **Conditional compilation**: Use `#ifdef MFEM_USE_MPI` for parallel-specific code paths
3. **Serial-first verification**: Complete all unit tests with serial version before parallel implementation
4. **Incremental parallelization**: Parallelize components one at a time with verification at each step

**MFEM Serial/Parallel Class Mapping**:
| Serial Class | Parallel Class | Notes |
|--------------|----------------|-------|
| `Mesh` | `ParMesh` | Distributed mesh with ghost elements |
| `FiniteElementSpace` | `ParFiniteElementSpace` | True DOFs vs local DOFs |
| `BilinearForm` | `ParBilinearForm` | Shared face assembly |
| `LinearForm` | `ParLinearForm` | Parallel vector assembly |
| `GridFunction` | `ParGridFunction` | Distributed solution vector |
| `Vector` | `HypreParVector` | Hypre-compatible parallel vectors |
| `CGSolver` | `HyprePCG` | Parallel Krylov solver |
| `GSSmoother` | `HypreBoomerAMG` | Algebraic multigrid preconditioner |

**Fault Parallelization Considerations**:
- Fault DOFs are distributed across MPI ranks owning fault boundary elements
- State variables (θ) are local to each rank's fault portion
- Slip rates (V) require no inter-rank communication (computed from local traction)
- Output stations may span multiple ranks - need gathering for I/O

### 2.2 Generality Requirements

The implementation must support future extensions:

| Current (BP2) | Future Extension | Design Implication |
|---------------|------------------|-------------------|
| 2D antiplane (scalar) | 2D plane strain (vector) | Abstract PDE operator |
| 2D antiplane (scalar) | 3D elasticity (vector) | Dimension-templated classes |
| Boundary fault | Interior fault | Face integrator abstraction |
| Quasi-dynamic | Fully-dynamic | Time operator modes |
| Aging law | Slip law, other laws | Friction law interface |
| Serial | MPI Parallel | Abstract mesh/FE interfaces |

### 2.3 Class Hierarchy Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                      SEASBenchmark (Driver)                                 │
│                      - Handles both serial and MPI execution                │
│                      - Uses MPI_COMM_WORLD when available                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   ┌─────────────────────────┐    ┌─────────────────────────┐               │
│   │ DomainOperator<MESH>    │◄──►│  FaultOperator<MESH>    │               │
│   │ (Template Abstract)     │    │  (Template Abstract)    │               │
│   └───────────┬─────────────┘    └───────────┬─────────────┘               │
│               │                              │                              │
│   ┌───────────┴───────────┐      ┌───────────┴───────────┐                 │
│   │ AntiplaneDomain       │      │ RateStateFault        │                 │
│   │  - Serial: Mesh       │      │  - Serial: local data │                 │
│   │  - Parallel: ParMesh  │      │  - Parallel: dist.    │                 │
│   │ PlaneStrainDomain     │      └───────────────────────┘                 │
│   │ Elastic3DDomain       │                                                │
│   └───────────────────────┘                                                │
│                                                                             │
│   ┌───────────────────────────────────────────────────────────────────────┐│
│   │              SEASOperator : TimeDependentOperator                     ││
│   │  ┌─────────────────────┐  ┌─────────────────────┐                     ││
│   │  │QuasiDynamicMode     │  │FullyDynamicMode     │                     ││
│   │  │ - Serial/Parallel   │  │ - Serial/Parallel   │                     ││
│   │  └─────────────────────┘  └─────────────────────┘                     ││
│   └───────────────────────────────────────────────────────────────────────┘│
│                                                                             │
│   ┌───────────────────┐    ┌───────────────────┐    ┌───────────────────┐  │
│   │  FrictionLaw      │    │  StateEvolution   │    │  MPIContext       │  │
│   │  (Abstract)       │    │  (Abstract)       │    │  (Comm wrapper)   │  │
│   └─────────┬─────────┘    └─────────┬─────────┘    └───────────────────┘  │
│             │                        │                                      │
│   ┌─────────┴─────────┐    ┌─────────┴─────────┐                           │
│   │ DieterichRuina    │    │ AgingLaw          │                           │
│   │ SlipWeakening     │    │ SlipLaw           │                           │
│   └───────────────────┘    └───────────────────┘                           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 2.4 Parallel Execution Model

```
┌───────────────────────────────────────────────────────────────────────────┐
│                         MPI Parallel Execution                            │
├───────────────────────────────────────────────────────────────────────────┤
│                                                                           │
│  Rank 0              Rank 1              Rank 2              Rank N       │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐│
│  │ ParMesh     │    │ ParMesh     │    │ ParMesh     │    │ ParMesh     ││
│  │ (local part)│    │ (local part)│    │ (local part)│    │ (local part)││
│  │  + ghosts   │    │  + ghosts   │    │  + ghosts   │    │  + ghosts   ││
│  └──────┬──────┘    └──────┬──────┘    └──────┬──────┘    └──────┬──────┘│
│         │                  │                  │                  │        │
│  ┌──────┴──────┐    ┌──────┴──────┐    ┌──────┴──────┐    ┌──────┴──────┐│
│  │ Local Fault │    │ Local Fault │    │ Local Fault │    │ Local Fault ││
│  │ DOFs + θ, V │    │ DOFs + θ, V │    │ DOFs + θ, V │    │ DOFs + θ, V ││
│  └──────┬──────┘    └──────┴──────┘    └──────┴──────┘    └──────┬──────┘│
│         │                  │                  │                  │        │
│         └──────────────────┴──────────────────┴──────────────────┘        │
│                                 │                                          │
│                    ┌────────────┴────────────┐                            │
│                    │   HyprePCG + BoomerAMG  │                            │
│                    │   (Parallel solve)      │                            │
│                    └─────────────────────────┘                            │
│                                                                           │
│  Communication Points:                                                    │
│  1. Domain solve: HYPRE handles parallel linear algebra                   │
│  2. Fault state: No communication needed (local computation)              │
│  3. Output: MPI_Gather for probe stations spanning ranks                  │
│  4. Time stepping: MPI_Allreduce for global max slip rate (CFL)           │
│                                                                           │
└───────────────────────────────────────────────────────────────────────────┘
```

## 3. File Organization

```
miniapps/seas/
├── CMakeLists.txt                    # Builds both serial and parallel versions
├── seas.cpp                          # Main driver for all benchmarks (serial)
├── pseas.cpp                         # Parallel driver (MPI version)
│
├── common/
│   ├── seas_types.hpp                # Type aliases for serial/parallel switching
│   ├── mpi_context.hpp               # MPI wrapper (no-op in serial)
│   └── parallel_utils.hpp            # MPI helper functions
│
├── config/
│   ├── seas_config.hpp               # Configuration parsing
│   ├── bp2_params.hpp                # BP2-specific parameters
│   └── bp1_params.hpp                # BP1 parameters (for comparison)
│
├── domain/
│   ├── domain_operator.hpp           # Abstract domain operator (templated)
│   ├── antiplane_operator.hpp        # 2D antiplane (scalar Laplacian)
│   ├── antiplane_operator_impl.hpp   # Implementation (shared serial/parallel)
│   ├── plane_strain_operator.hpp     # 2D plane strain [future]
│   └── elastic3d_operator.hpp        # 3D elasticity [future]
│
├── fault/
│   ├── fault_operator.hpp            # Abstract fault operator
│   ├── rate_state_fault.hpp          # Rate-and-state fault implementation
│   ├── fault_geometry.hpp            # Fault surface geometry
│   ├── fault_output.hpp              # Fault output handling
│   └── parallel_fault_data.hpp       # Distributed fault data management
│
├── friction/
│   ├── friction_law.hpp              # Abstract friction law interface
│   ├── dieterich_ruina.hpp           # Dieterich-Ruina friction
│   └── state_evolution.hpp           # State evolution laws (aging, slip)
│
├── solver/
│   ├── seas_operator.hpp             # TimeDependentOperator for SEAS
│   ├── quasi_dynamic.hpp             # Quasi-dynamic coupling
│   ├── fully_dynamic.hpp             # Fully-dynamic coupling [future]
│   ├── time_stepper.hpp              # Adaptive time stepping
│   └── parallel_time_stepper.hpp     # Parallel-aware time control
│
├── io/
│   ├── benchmark_output.hpp          # SCEC-format output
│   ├── probe_output.hpp              # Time series at probe locations
│   ├── parallel_probe_output.hpp     # MPI-aware probe gathering
│   └── checkpoint.hpp                # Checkpointing support (serial/parallel)
│
├── tests/
│   ├── unit/
│   │   ├── test_friction_law.cpp     # Friction law unit tests
│   │   ├── test_state_evolution.cpp  # State evolution unit tests
│   │   ├── test_slip_rate_solver.cpp # Newton solver unit tests
│   │   ├── test_initial_state.cpp    # Initial condition unit tests
│   │   ├── test_antiplane.cpp        # Antiplane operator unit tests
│   │   └── test_traction.cpp         # Traction computation unit tests
│   ├── integration/
│   │   ├── test_quasi_dynamic.cpp    # QD coupling integration test
│   │   └── test_bp2_short.cpp        # Short BP2 simulation test
│   ├── parallel/                     # MPI-specific tests
│   │   ├── test_parallel_domain.cpp  # Parallel domain operator tests
│   │   ├── test_parallel_fault.cpp   # Distributed fault tests
│   │   ├── test_parallel_comm.cpp    # Communication pattern tests
│   │   └── test_scaling.cpp          # Strong/weak scaling tests
│   └── verification/
│       ├── mms_antiplane.cpp         # Method of manufactured solutions
│       ├── mms_antiplane_parallel.cpp# Parallel MMS verification
│       └── bp2_benchmark.cpp         # Full BP2 verification
│
└── examples/
    ├── bp2_qd.cpp                    # BP2 quasi-dynamic example (serial)
    ├── bp2_qd_parallel.cpp           # BP2 parallel example
    └── bp1_qd.cpp                    # BP1 for comparison
```

### 3.1 CMake Configuration for Serial/Parallel Builds

```cmake
# CMakeLists.txt

# Always build serial version
add_executable(seas seas.cpp)
target_link_libraries(seas mfem)

# Build parallel version if MPI is available
if(MFEM_USE_MPI)
    add_executable(pseas pseas.cpp)
    target_link_libraries(pseas mfem)
    target_compile_definitions(pseas PRIVATE SEAS_USE_MPI)

    # Parallel tests
    add_executable(test_parallel_domain tests/parallel/test_parallel_domain.cpp)
    target_link_libraries(test_parallel_domain mfem gtest)

    # MPI test runner
    add_test(NAME parallel_domain_np2
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                     $<TARGET_FILE:test_parallel_domain>)
    add_test(NAME parallel_domain_np4
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                     $<TARGET_FILE:test_parallel_domain>)
endif()
```

---

# Part III: Detailed Component Design

## 3.5 Serial/Parallel Type Abstraction

The key to supporting both serial and parallel execution is abstracting MFEM types:

```cpp
// common/seas_types.hpp

#ifndef SEAS_TYPES_HPP
#define SEAS_TYPES_HPP

#include "mfem.hpp"

namespace mfem {
namespace seas {

#ifdef SEAS_USE_MPI

// Parallel type aliases
using SEASMesh = ParMesh;
using SEASFiniteElementSpace = ParFiniteElementSpace;
using SEASBilinearForm = ParBilinearForm;
using SEASLinearForm = ParLinearForm;
using SEASGridFunction = ParGridFunction;

// Parallel solver types
using SEASSolver = HyprePCG;
using SEASPreconditioner = HypreBoomerAMG;

// MPI communicator wrapper
inline MPI_Comm GetSEASComm() { return MPI_COMM_WORLD; }
inline int GetSEASRank() {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank;
}
inline int GetSEASSize() {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size); return size;
}
inline bool IsSEASRoot() { return GetSEASRank() == 0; }

#else

// Serial type aliases
using SEASMesh = Mesh;
using SEASFiniteElementSpace = FiniteElementSpace;
using SEASBilinearForm = BilinearForm;
using SEASLinearForm = LinearForm;
using SEASGridFunction = GridFunction;

// Serial solver types
using SEASSolver = CGSolver;
using SEASPreconditioner = GSSmoother;

// Serial stubs for MPI functions
inline int GetSEASRank() { return 0; }
inline int GetSEASSize() { return 1; }
inline bool IsSEASRoot() { return true; }

#endif

} // namespace seas
} // namespace mfem

#endif // SEAS_TYPES_HPP
```

### 3.6 MPI Context Wrapper

```cpp
// common/mpi_context.hpp

#ifndef MPI_CONTEXT_HPP
#define MPI_CONTEXT_HPP

#include "mfem.hpp"

namespace mfem {
namespace seas {

/// RAII wrapper for MPI initialization (no-op in serial)
class MPIContext {
public:
#ifdef SEAS_USE_MPI
    MPIContext(int *argc, char ***argv) {
        MPI_Init(argc, argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
    }
    ~MPIContext() { MPI_Finalize(); }
#else
    MPIContext(int *argc, char ***argv) : rank_(0), size_(1) {}
    ~MPIContext() = default;
#endif

    int Rank() const { return rank_; }
    int Size() const { return size_; }
    bool IsRoot() const { return rank_ == 0; }

    /// Global reduction for maximum value (for CFL condition)
    real_t GlobalMax(real_t local_val) const {
#ifdef SEAS_USE_MPI
        real_t global_val;
        MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);
        return global_val;
#else
        return local_val;
#endif
    }

    /// Global reduction for minimum value
    real_t GlobalMin(real_t local_val) const {
#ifdef SEAS_USE_MPI
        real_t global_val;
        MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);
        return global_val;
#else
        return local_val;
#endif
    }

    /// Barrier synchronization
    void Barrier() const {
#ifdef SEAS_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

private:
    int rank_;
    int size_;
};

} // namespace seas
} // namespace mfem

#endif // MPI_CONTEXT_HPP
```

## 4. Friction Law Implementation

### 4.1 Abstract Interface

```cpp
// friction/friction_law.hpp

class FrictionLaw {
public:
    virtual ~FrictionLaw() = default;

    /// Compute friction coefficient f(V, θ)
    virtual real_t FrictionCoefficient(real_t V, real_t theta,
                                        const Vector &params) const = 0;

    /// Compute fault strength F = σn * f(V, θ)
    virtual real_t FaultStrength(real_t V, real_t theta, real_t sigma_n,
                                  const Vector &params) const = 0;

    /// Solve for slip rate V given stress τ and state θ
    /// Solves: τ = σn * f(V, θ) + η * V
    virtual real_t SolveSlipRate(real_t tau, real_t theta, real_t sigma_n,
                                  real_t eta, const Vector &params) const = 0;

    /// Compute initial state θ(0) from equilibrium
    virtual real_t InitialState(real_t tau0, real_t V_init, real_t sigma_n,
                                 real_t eta, const Vector &params) const = 0;
};
```

### 4.2 Dieterich-Ruina Implementation

```cpp
// friction/dieterich_ruina.hpp

class DieterichRuinaFriction : public FrictionLaw {
public:
    // Global constants
    struct Constants {
        real_t V0 = 1.0e-6;    // Reference slip rate [m/s]
        real_t f0 = 0.6;       // Reference friction coefficient
        real_t b  = 0.015;     // State evolution parameter
        real_t Dc = 0.004;     // Critical slip distance [m] (BP2 value)
    };

    DieterichRuinaFriction(const Constants &c) : cp_(c) {}

    /// f(V, θ) = a * asinh[(V / 2V₀) * exp((f₀ + b*ln(V₀θ/Dc)) / a)]
    real_t FrictionCoefficient(real_t V, real_t theta,
                                const Vector &params) const override {
        real_t a = params(0);  // Rate-and-state parameter a
        real_t arg = (V / (2.0 * cp_.V0)) *
                     std::exp((cp_.f0 + cp_.b * std::log(cp_.V0 * theta / cp_.Dc)) / a);
        return a * std::asinh(arg);
    }

    real_t FaultStrength(real_t V, real_t theta, real_t sigma_n,
                          const Vector &params) const override {
        return sigma_n * FrictionCoefficient(V, theta, params);
    }

    /// Solve τ = σn * f(V, θ) + η * V for V using Newton-Raphson/Brent
    real_t SolveSlipRate(real_t tau, real_t theta, real_t sigma_n,
                          real_t eta, const Vector &params) const override;

    /// Compute θ(0) such that τ⁰ = σn * f(V_init, θ(0)) + η * V_init
    real_t InitialState(real_t tau0, real_t V_init, real_t sigma_n,
                         real_t eta, const Vector &params) const override;

private:
    Constants cp_;

    // Internal helper: compute f and df/dV for Newton solver
    void FrictionAndDerivative(real_t V, real_t theta, real_t a,
                                real_t &f, real_t &df_dV) const;
};
```

### 4.3 State Evolution Laws

```cpp
// friction/state_evolution.hpp

class StateEvolution {
public:
    virtual ~StateEvolution() = default;

    /// Compute dθ/dt given current V and θ
    virtual real_t Rate(real_t V, real_t theta, real_t Dc) const = 0;

    /// Compute steady-state θ for given V
    virtual real_t SteadyState(real_t V, real_t Dc) const = 0;
};

/// Aging law: dθ/dt = 1 - V*θ/Dc
class AgingLaw : public StateEvolution {
public:
    real_t Rate(real_t V, real_t theta, real_t Dc) const override {
        return 1.0 - V * theta / Dc;
    }

    real_t SteadyState(real_t V, real_t Dc) const override {
        return Dc / V;
    }
};

/// Slip law: dθ/dt = -V*θ/Dc * ln(V*θ/Dc)
class SlipLaw : public StateEvolution {
public:
    real_t Rate(real_t V, real_t theta, real_t Dc) const override {
        real_t x = V * theta / Dc;
        return -x * std::log(x);
    }

    real_t SteadyState(real_t V, real_t Dc) const override {
        return Dc / V;
    }
};
```

## 5. Domain Operator Design

### 5.1 Abstract Domain Operator (Templated for Serial/Parallel)

```cpp
// domain/domain_operator.hpp

#include "common/seas_types.hpp"

namespace mfem {
namespace seas {

/// Abstract base class for domain operators
/// Template parameter controls serial vs parallel MFEM types
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

    virtual ~DomainOperator() = default;

    /// Get number of displacement components
    virtual int NumComponents() const = 0;

    /// Solve domain problem with given fault slip BC
    virtual void Solve(real_t time, const Vector &slip_bc,
                       GridFuncType &displacement) = 0;

    /// Compute traction at fault from displacement
    virtual void ComputeTraction(const GridFuncType &displacement,
                                  Vector &traction) = 0;

    /// Get reference to finite element space
    virtual FESpaceType &GetFESpace() = 0;

    /// Get the underlying mesh
    virtual MeshType &GetMesh() = 0;

#ifdef SEAS_USE_MPI
    /// Get MPI communicator (parallel only)
    virtual MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#endif
};

// Convenience type aliases
using SerialDomainOperator = DomainOperator<Mesh>;
using ParallelDomainOperator = DomainOperator<ParMesh>;

} // namespace seas
} // namespace mfem
```

### 5.2 Antiplane Operator (BP2) - Serial Version

```cpp
// domain/antiplane_operator.hpp

namespace mfem {
namespace seas {

/// Solves ∇²u = 0 with slip BC on fault
/// Template parameter: Mesh for serial, ParMesh for parallel
template <typename MeshType = Mesh>
class AntiplaneDomainOperator : public DomainOperator<MeshType> {
public:
    // Inherit type aliases from base
    using Base = DomainOperator<MeshType>;
    using FESpaceType = typename Base::FESpaceType;
    using GridFuncType = typename Base::GridFuncType;

    // Additional type aliases for bilinear form
    using BilinFormType = typename std::conditional<
        std::is_same<MeshType, ParMesh>::value,
        ParBilinearForm, BilinearForm>::type;
    using LinFormType = typename std::conditional<
        std::is_same<MeshType, ParMesh>::value,
        ParLinearForm, LinearForm>::type;

    AntiplaneDomainOperator(MeshType &mesh, int order, real_t mu);

    int NumComponents() const override { return 1; }

    void Solve(real_t time, const Vector &slip_bc,
               GridFuncType &displacement) override;

    void ComputeTraction(const GridFuncType &displacement,
                          Vector &traction) override;

    FESpaceType &GetFESpace() override { return *fes_; }
    MeshType &GetMesh() override { return mesh_; }

private:
    MeshType &mesh_;
    real_t mu_;  // Shear modulus

    std::unique_ptr<FiniteElementCollection> fec_;
    std::unique_ptr<FESpaceType> fes_;
    std::unique_ptr<BilinFormType> a_form_;

    // Solver (type depends on serial/parallel)
    std::unique_ptr<Solver> solver_;
    std::unique_ptr<Solver> prec_;

    // Boundary attributes
    Array<int> fault_bdr_attr_;
    Array<int> farfield_bdr_attr_;
    Array<int> freesurface_bdr_attr_;

    // Work vectors
    mutable Vector rhs_, sol_;

    // Setup methods
    void SetupBilinearForm();
    void SetupSolver();
};

// Explicit instantiation declarations
extern template class AntiplaneDomainOperator<Mesh>;
#ifdef SEAS_USE_MPI
extern template class AntiplaneDomainOperator<ParMesh>;
#endif

// Convenience type aliases
using SerialAntiplaneOperator = AntiplaneDomainOperator<Mesh>;
using ParallelAntiplaneOperator = AntiplaneDomainOperator<ParMesh>;

} // namespace seas
} // namespace mfem
```

### 5.3 Antiplane Operator - Parallel-Specific Implementation Details

```cpp
// domain/antiplane_operator_impl.hpp (partial, showing parallel specifics)

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

template <>
void AntiplaneDomainOperator<Mesh>::SetupSolver() {
    // Serial solver setup
    auto *gs = new GSSmoother();
    auto *cg = new CGSolver();
    cg->SetRelTol(1e-12);
    cg->SetMaxIter(500);
    cg->SetPrintLevel(0);
    cg->SetPreconditioner(*gs);

    prec_.reset(gs);
    solver_.reset(cg);
}

template <typename MeshType>
void AntiplaneDomainOperator<MeshType>::Solve(
    real_t time, const Vector &slip_bc, GridFuncType &displacement)
{
    // Common solve logic for both serial and parallel
    // BilinearForm handles parallel assembly internally via ParBilinearForm

    // Apply far-field BC: u = ±Vp*t/2 at x boundaries
    ApplyFarFieldBC(time, displacement);

    // Apply fault slip BC
    ApplyFaultSlipBC(slip_bc, displacement);

    // Solve the system
    // For ParBilinearForm, this uses HypreParMatrix internally
    Vector &X = displacement.GetTrueDofs();  // True DOFs in parallel
    a_form_->FormLinearSystem(ess_tdof_list_, displacement, *b_, A_, X, B_);
    solver_->SetOperator(*A_);
    solver_->Mult(B_, X);
    a_form_->RecoverFEMSolution(X, *b_, displacement);
}
```

## 6. Fault Operator Design

### 6.1 Rate-State Fault Operator (Serial/Parallel Compatible)

```cpp
// fault/rate_state_fault.hpp

#include "common/seas_types.hpp"
#include "common/mpi_context.hpp"

namespace mfem {
namespace seas {

/// Manages rate-and-state friction computation on the fault
/// Works with both serial and parallel execution
class RateStateFaultOperator {
public:
    RateStateFaultOperator(FaultGeometry *geom,
                           FrictionLaw *friction,
                           StateEvolution *state_law,
                           const BP2Params &params,
                           const MPIContext *mpi_ctx = nullptr);

    /// Number of LOCAL state DOFs (slip + state variable)
    /// In parallel, this is the local portion only
    int LocalStateSize() const { return num_local_fault_dofs_ * 2; }

    /// Number of GLOBAL state DOFs
    int GlobalStateSize() const { return num_global_fault_dofs_ * 2; }

    /// Initialize state: set δ(0) = 0, compute θ(0) from τ⁰
    void Initialize(const Vector &tau0, Vector &state);

    /// Compute RHS: [dδ/dt, dθ/dt] given traction and state
    void ComputeRHS(real_t time, const Vector &traction,
                    const Vector &state, Vector &rate);

    /// Get slip from state vector (for boundary condition)
    void GetSlip(const Vector &state, Vector &slip) const;

    /// Get maximum slip rate (LOCAL, call GlobalMaxSlipRate for global)
    real_t GetLocalMaxSlipRate() const { return V_max_local_; }

    /// Get GLOBAL maximum slip rate (requires MPI reduction in parallel)
    real_t GetGlobalMaxSlipRate() const {
        if (mpi_ctx_) {
            return mpi_ctx_->GlobalMax(V_max_local_);
        }
        return V_max_local_;
    }

    /// Access local fault DOF information
    int NumLocalDOFs() const { return num_local_fault_dofs_; }

    /// Get local-to-global mapping for fault DOFs (parallel only)
    const Array<int> &GetLocalToGlobalMap() const { return local_to_global_; }

private:
    FaultGeometry *geom_;
    FrictionLaw *friction_;
    StateEvolution *state_law_;
    BP2Params params_;
    const MPIContext *mpi_ctx_;

    int num_local_fault_dofs_;   // Local fault DOFs on this rank
    int num_global_fault_dofs_;  // Total fault DOFs across all ranks
    real_t V_max_local_;

    // Precomputed depth-dependent parameters (local portion)
    Vector a_values_;      // a(z) at each local fault DOF
    Vector sigma_n_;       // Normal stress
    Vector eta_;           // Radiation damping coefficient
    real_t tau0_;          // Pre-stress

    // Parallel data distribution
    Array<int> local_to_global_;  // Mapping from local to global fault DOF
    Array<int> fault_dof_owner_;  // Which rank owns each global fault DOF
};

} // namespace seas
} // namespace mfem
```

### 6.2 Parallel Fault Data Distribution

```cpp
// fault/parallel_fault_data.hpp

#include "common/seas_types.hpp"

namespace mfem {
namespace seas {

/// Manages distribution of fault data across MPI ranks
class ParallelFaultData {
public:
    /// Construct from parallel mesh and fault boundary attribute
    ParallelFaultData(ParMesh &pmesh, int fault_attr);

    /// Get number of fault DOFs owned by this rank
    int NumLocalDOFs() const { return num_local_; }

    /// Get total number of fault DOFs globally
    int NumGlobalDOFs() const { return num_global_; }

    /// Get local-to-global DOF mapping
    const Array<int> &LocalToGlobal() const { return l2g_; }

    /// Get global-to-local DOF mapping (-1 if not local)
    int GlobalToLocal(int global_dof) const;

    /// Check if a global DOF is owned by this rank
    bool IsOwned(int global_dof) const;

    /// Get depth (z-coordinate) for each local fault DOF
    const Vector &GetLocalDepths() const { return local_depths_; }

    /// Gather local data to global array on root rank
    /// Used for I/O operations
    void GatherToRoot(const Vector &local_data, Vector &global_data) const;

    /// Scatter global data from root to local arrays
    void ScatterFromRoot(const Vector &global_data, Vector &local_data) const;

    /// Interpolate global data to specific probe depths (on root only)
    void InterpolateToProbes(const Vector &global_data,
                              const Array<real_t> &probe_depths,
                              Vector &probe_values) const;

private:
    ParMesh &pmesh_;
    int fault_attr_;

    int num_local_;
    int num_global_;
    Array<int> l2g_;
    Array<int> g2l_;
    Array<int> owner_;
    Vector local_depths_;

    // For MPI communication
    Array<int> recv_counts_;
    Array<int> recv_displs_;

    void SetupDistribution();
};

} // namespace seas
} // namespace mfem
```

### 6.3 Fault DOF Distribution Strategy

In parallel execution, fault DOFs are distributed based on mesh partitioning:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    Fault DOF Distribution Example                       │
│                    (Mesh partitioned by depth)                          │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Depth (z)    0 km                                            40 km    │
│     │         │                                                 │       │
│     ▼         ▼                                                 ▼       │
│  ┌──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐   │
│  │  Rank 0  │  Rank 1  │  Rank 2  │  Rank 3  │  Rank 4  │  Rank 5  │   │
│  │ (0-7km)  │ (7-14km) │(14-20km) │(20-26km) │(26-33km) │(33-40km) │   │
│  │  ~10 DOFs│  ~10 DOFs│  ~10 DOFs│  ~10 DOFs│  ~10 DOFs│  ~10 DOFs│   │
│  └──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘   │
│                                                                         │
│  Key Points:                                                            │
│  - Each rank owns fault DOFs on its local mesh portion                  │
│  - State variables (θ) are computed locally (no communication)          │
│  - Slip rates (V) are computed locally from local traction              │
│  - Global max V requires MPI_Allreduce for time step control            │
│  - Output stations may span ranks (need gather for I/O)                 │
│                                                                         │
│  Output Stations (may not align with rank boundaries):                  │
│  z = 0 km → Rank 0 (local)                                              │
│  z = 4.8 km → Rank 0 (local)                                            │
│  z = 12 km → Rank 1 (local)                                             │
│  z = 16.8 km → Rank 2 (local)                                           │
│  z = 24 km → Rank 3 (local)                                             │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## 7. SEAS Operator (Time Integration)

### 7.1 Quasi-Dynamic Operator (Serial/Parallel Compatible)

```cpp
// solver/seas_operator.hpp

#include "common/seas_types.hpp"
#include "common/mpi_context.hpp"

namespace mfem {
namespace seas {

/// SEAS quasi-dynamic time integration operator
/// Template parameter controls serial vs parallel execution
template <typename MeshType = SEASMesh>
class SEASQuasiDynamicOperator : public TimeDependentOperator {
public:
    using DomainOpType = DomainOperator<MeshType>;
    using GridFuncType = typename DomainOpType::GridFuncType;

    SEASQuasiDynamicOperator(DomainOpType *domain,
                              RateStateFaultOperator *fault,
                              const MPIContext *mpi_ctx = nullptr);

    /// Initialize state at t=0
    void SetInitialCondition(Vector &state);

    /// Compute d(state)/dt = RHS(t, state)
    /// state = [δ₁, θ₁, δ₂, θ₂, ...]  (interleaved)
    void Mult(const Vector &state, Vector &rate) const override;

    /// Get current displacement solution
    const GridFuncType &GetDisplacement() const { return *u_gf_; }

    /// Get traction at fault (local portion)
    const Vector &GetTraction() const { return traction_; }

    /// Get maximum slip rate (GLOBAL in parallel)
    real_t GetMaxSlipRate() const {
        return fault_->GetGlobalMaxSlipRate();
    }

    /// Get local maximum slip rate (for debugging)
    real_t GetLocalMaxSlipRate() const {
        return fault_->GetLocalMaxSlipRate();
    }

private:
    DomainOpType *domain_;
    RateStateFaultOperator *fault_;
    const MPIContext *mpi_ctx_;

    std::unique_ptr<GridFuncType> u_gf_;

    mutable Vector slip_;
    mutable Vector traction_;

    /// Internal: solve domain and update traction
    void UpdateDomainAndTraction(const Vector &state) const;
};

// Convenience type aliases
using SerialSEASOperator = SEASQuasiDynamicOperator<Mesh>;
using ParallelSEASOperator = SEASQuasiDynamicOperator<ParMesh>;

} // namespace seas
} // namespace mfem
```

### 7.2 Parallel Time Stepper with Adaptive Control

```cpp
// solver/parallel_time_stepper.hpp

#include "common/mpi_context.hpp"

namespace mfem {
namespace seas {

/// Adaptive time stepper aware of parallel execution
/// Synchronizes time step across all MPI ranks
class ParallelAdaptiveTimeStepper {
public:
    ParallelAdaptiveTimeStepper(const MPIContext *mpi_ctx = nullptr);

    /// Configuration
    void SetDtMin(real_t dt_min) { dt_min_ = dt_min; }
    void SetDtMax(real_t dt_max) { dt_max_ = dt_max; }
    void SetTargetSlipRateMax(real_t V_target) { V_target_ = V_target; }

    /// Compute new time step based on slip rate
    /// In parallel, uses global max slip rate
    real_t ComputeNewDt(real_t V_max_local, real_t dt_current) const {
        // Get global max slip rate
        real_t V_max = mpi_ctx_ ? mpi_ctx_->GlobalMax(V_max_local)
                                : V_max_local;

        // Adaptive formula based on slip rate
        real_t dt_new;
        if (V_max < 1e-12) {
            dt_new = dt_max_;
        } else if (V_max > V_target_) {
            // Rapid slip: reduce time step
            dt_new = dt_current * V_target_ / V_max;
        } else {
            // Slow slip: can increase time step
            dt_new = std::min(dt_current * 1.5, dt_max_);
        }

        // Clamp to bounds
        dt_new = std::max(dt_min_, std::min(dt_max_, dt_new));

        return dt_new;
    }

    /// Get current time step (same on all ranks)
    real_t GetDt() const { return dt_; }

    /// Update time step (must be called synchronously by all ranks)
    void SetDt(real_t dt) { dt_ = dt; }

private:
    const MPIContext *mpi_ctx_;
    real_t dt_min_ = 1e-6;          // Minimum time step (seconds)
    real_t dt_max_ = 3.15e7;        // Maximum time step (1 year)
    real_t V_target_ = 1e-6;        // Target max slip rate for dt control
    real_t dt_ = 1e3;               // Current time step
};

} // namespace seas
} // namespace mfem
```

### 7.3 Main Driver Structure (Serial/Parallel)

```cpp
// seas.cpp (serial) / pseas.cpp (parallel)

#include "common/seas_types.hpp"
#include "common/mpi_context.hpp"
#include "domain/antiplane_operator.hpp"
#include "fault/rate_state_fault.hpp"
#include "solver/seas_operator.hpp"
#include "io/benchmark_output.hpp"

int main(int argc, char *argv[]) {
    // Initialize MPI context (no-op in serial)
    mfem::seas::MPIContext mpi_ctx(&argc, &argv);

    // Parse command line
    mfem::OptionsParser args(argc, argv);
    // ... add options ...

    // Load mesh (serial or parallel)
#ifdef SEAS_USE_MPI
    mfem::ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
#else
    mfem::Mesh &mesh = serial_mesh;
#endif

    // Create operators using type aliases
    mfem::seas::AntiplaneDomainOperator<mfem::seas::SEASMesh>
        domain(mesh, order, mu);

    mfem::seas::RateStateFaultOperator fault(
        &geom, &friction, &aging_law, params, &mpi_ctx);

    mfem::seas::SEASQuasiDynamicOperator<mfem::seas::SEASMesh>
        seas_op(&domain, &fault, &mpi_ctx);

    // Time integration
    mfem::seas::ParallelAdaptiveTimeStepper stepper(&mpi_ctx);
    mfem::RK4Solver ode_solver;
    ode_solver.Init(seas_op);

    Vector state(fault.LocalStateSize());
    seas_op.SetInitialCondition(state);

    real_t t = 0.0, t_final = params.t_final;
    while (t < t_final) {
        real_t dt = stepper.GetDt();
        ode_solver.Step(state, t, dt);
        t += dt;

        // Update time step based on global max slip rate
        real_t V_max = seas_op.GetMaxSlipRate();
        stepper.SetDt(stepper.ComputeNewDt(V_max, dt));

        // Output (only on root in parallel)
        if (mpi_ctx.IsRoot()) {
            output.Write(t, state);
        }
    }

    return 0;
}
```

---

# Part IV: Unit Test Design

## 8. Test Strategy

### 8.1 Testing Levels

| Level | Purpose | Coverage |
|-------|---------|----------|
| Unit | Test individual functions | Friction, state evolution, solver |
| Integration | Test component interactions | Domain-fault coupling |
| Verification | Test against known solutions | MMS, analytical |
| Validation | Test against benchmark | BP2 data comparison |

### 8.2 Critical Functions to Test

Based on the fault interface implementation guide, these functions are critical:

1. **Friction coefficient**: f(V, θ)
2. **Slip rate solver**: Solve τ = σn·f(V,θ) + η·V for V
3. **Initial state**: θ(0) from τ⁰, V_init
4. **State evolution**: dθ/dt = 1 - V·θ/Dc
5. **Pre-stress**: τ⁰ computation
6. **Traction**: T = ∂u/∂x · μ at fault

## 9. Unit Test Specifications

### 9.1 Friction Law Tests

```cpp
// tests/unit/test_friction_law.cpp

/// Test 1: Friction coefficient at reference state
/// At V = V₀, θ = Dc/V₀ (steady state), f should equal f₀
TEST(FrictionLaw, ReferenceState) {
    DieterichRuinaFriction::Constants c;
    c.V0 = 1e-6; c.f0 = 0.6; c.b = 0.015; c.Dc = 0.004;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.015;  // a = b (neutral)
    real_t V = c.V0;
    real_t theta = c.Dc / c.V0;  // steady state

    real_t f = law.FrictionCoefficient(V, theta, params);

    // At reference: f = a * asinh[(V₀/2V₀) * exp((f₀ + b*ln(1))/a)]
    //             = a * asinh[0.5 * exp(f₀/a)]
    real_t expected = params(0) * std::asinh(0.5 * std::exp(c.f0 / params(0)));
    EXPECT_NEAR(f, expected, 1e-12);
}

/// Test 2: Friction increases with slip rate (direct effect)
TEST(FrictionLaw, DirectEffect) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.010;
    real_t theta = 1000.0;  // Fixed state

    real_t V1 = 1e-9, V2 = 1e-6;
    real_t f1 = law.FrictionCoefficient(V1, theta, params);
    real_t f2 = law.FrictionCoefficient(V2, theta, params);

    EXPECT_GT(f2, f1);  // Higher V → higher f (at fixed θ)
}

/// Test 3: Friction decreases with state (evolution effect)
TEST(FrictionLaw, EvolutionEffect) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.010;
    real_t V = 1e-9;

    real_t theta1 = 100.0, theta2 = 10000.0;
    real_t f1 = law.FrictionCoefficient(V, theta1, params);
    real_t f2 = law.FrictionCoefficient(V, theta2, params);

    EXPECT_GT(f2, f1);  // Higher θ → higher f (at fixed V)
}

/// Test 4: Velocity-weakening in nucleation zone (a < b)
TEST(FrictionLaw, VelocityWeakening) {
    DieterichRuinaFriction::Constants c;
    c.b = 0.015;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.010;  // a < b

    // At steady state: θ_ss = Dc/V
    // f_ss = a * asinh[(V/2V₀) * exp((f₀ + b*ln(V₀·Dc/(V·Dc)))/a)]
    //      = a * asinh[(V/2V₀) * exp((f₀ + b*ln(V₀/V))/a)]

    real_t V1 = 1e-10, V2 = 1e-8;
    real_t theta1_ss = c.Dc / V1;
    real_t theta2_ss = c.Dc / V2;

    real_t f1 = law.FrictionCoefficient(V1, theta1_ss, params);
    real_t f2 = law.FrictionCoefficient(V2, theta2_ss, params);

    // For a < b, steady-state friction decreases with V
    EXPECT_GT(f1, f2);
}

/// Test 5: Velocity-strengthening in creeping zone (a > b)
TEST(FrictionLaw, VelocityStrengthening) {
    DieterichRuinaFriction::Constants c;
    c.b = 0.015;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.025;  // a > b

    real_t V1 = 1e-10, V2 = 1e-8;
    real_t theta1_ss = c.Dc / V1;
    real_t theta2_ss = c.Dc / V2;

    real_t f1 = law.FrictionCoefficient(V1, theta1_ss, params);
    real_t f2 = law.FrictionCoefficient(V2, theta2_ss, params);

    // For a > b, steady-state friction increases with V
    EXPECT_LT(f1, f2);
}
```

### 9.2 Slip Rate Solver Tests

```cpp
// tests/unit/test_slip_rate_solver.cpp

/// Test 1: Solver recovers known slip rate
TEST(SlipRateSolver, RecoverKnownRate) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.015;
    real_t sigma_n = 50e6;  // 50 MPa
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);  // μ/(2cs)

    // Pick a V and θ, compute τ, then recover V
    real_t V_true = 1e-8;
    real_t theta = 5000.0;

    real_t f = law.FrictionCoefficient(V_true, theta, params);
    real_t tau = sigma_n * f + eta * V_true;

    real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, params);

    EXPECT_NEAR(V_solved, V_true, V_true * 1e-10);
}

/// Test 2: Solver handles zero radiation damping
TEST(SlipRateSolver, ZeroRadiationDamping) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.015;
    real_t sigma_n = 50e6;
    real_t eta = 0.0;  // No radiation damping

    real_t V_true = 1e-7;
    real_t theta = 3000.0;

    real_t f = law.FrictionCoefficient(V_true, theta, params);
    real_t tau = sigma_n * f;

    real_t V_solved = law.SolveSlipRate(tau, theta, sigma_n, eta, params);

    EXPECT_NEAR(V_solved, V_true, V_true * 1e-10);
}

/// Test 3: Solver handles extreme slip rates
TEST(SlipRateSolver, ExtremeRates) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.015;
    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);

    // Very slow (interseismic)
    real_t V_slow = 1e-12;
    real_t theta_slow = 1e8;
    real_t f_slow = law.FrictionCoefficient(V_slow, theta_slow, params);
    real_t tau_slow = sigma_n * f_slow + eta * V_slow;
    real_t V_solved_slow = law.SolveSlipRate(tau_slow, theta_slow, sigma_n, eta, params);
    EXPECT_NEAR(V_solved_slow, V_slow, V_slow * 1e-8);

    // Fast (coseismic)
    real_t V_fast = 1.0;  // 1 m/s
    real_t theta_fast = 0.01;
    real_t f_fast = law.FrictionCoefficient(V_fast, theta_fast, params);
    real_t tau_fast = sigma_n * f_fast + eta * V_fast;
    real_t V_solved_fast = law.SolveSlipRate(tau_fast, theta_fast, sigma_n, eta, params);
    EXPECT_NEAR(V_solved_fast, V_fast, V_fast * 1e-8);
}

/// Test 4: Solver convergence count
TEST(SlipRateSolver, ConvergenceEfficiency) {
    DieterichRuinaFriction::Constants c;
    DieterichRuinaFriction law(c);

    // Should converge in < 50 iterations for any reasonable input
    // (Implementation should track iteration count)
}
```

### 9.3 Initial State Tests

```cpp
// tests/unit/test_initial_state.cpp

/// Test 1: Initial state satisfies stress balance
TEST(InitialState, StressBalance) {
    DieterichRuinaFriction::Constants c;
    c.Dc = 0.004;  // BP2 value
    DieterichRuinaFriction law(c);

    Vector params(1); params(0) = 0.015;  // a = b (neutral)
    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);
    real_t V_init = 1e-9;

    // Compute pre-stress τ⁰ (BP2 formula at amax)
    real_t amax = 0.025;
    Vector params_max(1); params_max(0) = amax;
    real_t theta_ss = c.Dc / V_init;
    real_t f_init = law.FrictionCoefficient(V_init, theta_ss, params_max);
    real_t tau0 = sigma_n * f_init + eta * V_init;

    // Compute initial state
    real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, params);

    // Verify: τ⁰ = σn * f(V_init, θ₀) + η * V_init
    real_t f_check = law.FrictionCoefficient(V_init, theta0, params);
    real_t tau_check = sigma_n * f_check + eta * V_init;

    EXPECT_NEAR(tau_check, tau0, tau0 * 1e-10);
}

/// Test 2: Initial state varies with depth (a varies)
TEST(InitialState, DepthDependence) {
    DieterichRuinaFriction::Constants c;
    c.Dc = 0.004;
    DieterichRuinaFriction law(c);

    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);
    real_t V_init = 1e-9;

    // Compute τ⁰ at amax
    real_t amax = 0.025;
    Vector params_max(1); params_max(0) = amax;
    real_t theta_ss = c.Dc / V_init;
    real_t f_max = law.FrictionCoefficient(V_init, theta_ss, params_max);
    real_t tau0 = sigma_n * f_max + eta * V_init;

    // Compute θ₀ at different depths (different a values)
    std::vector<real_t> a_values = {0.010, 0.015, 0.020, 0.025};
    std::vector<real_t> theta0_values;

    for (real_t a : a_values) {
        Vector params(1); params(0) = a;
        real_t theta0 = law.InitialState(tau0, V_init, sigma_n, eta, params);
        theta0_values.push_back(theta0);

        // Verify stress balance
        real_t f_check = law.FrictionCoefficient(V_init, theta0, params);
        real_t tau_check = sigma_n * f_check + eta * V_init;
        EXPECT_NEAR(tau_check, tau0, tau0 * 1e-10);
    }

    // θ₀ should vary with a
    EXPECT_NE(theta0_values[0], theta0_values[3]);
}

/// Test 3: BP2-specific initial state values
TEST(InitialState, BP2Values) {
    // Compare with benchmark data initial values
    // At z=0: state_log10 = 3.602... → θ ≈ 4000 s
    // At z=12km: state_log10 = 3.602... → θ ≈ 4000 s

    DieterichRuinaFriction::Constants c;
    c.V0 = 1e-6; c.f0 = 0.6; c.b = 0.015; c.Dc = 0.004;
    DieterichRuinaFriction law(c);

    real_t sigma_n = 50e6;
    real_t eta = 0.5 * std::sqrt(32.04e9 * 2670.0);
    real_t V_init = 1e-9;

    // τ⁰ from benchmark: shear_stress ≈ 26.546 MPa
    real_t tau0_expected = 26.546e6;

    // θ₀ from benchmark: state_log10 ≈ 3.602 → θ ≈ 10^3.602 ≈ 4000 s
    real_t theta0_expected = std::pow(10.0, 3.602);

    // Compute using our formulas
    real_t amax = 0.025;
    Vector params_max(1); params_max(0) = amax;
    real_t theta_ss_max = c.Dc / V_init;
    real_t f_max = law.FrictionCoefficient(V_init, theta_ss_max, params_max);
    real_t tau0_computed = sigma_n * f_max + eta * V_init;

    // τ⁰ should match benchmark
    EXPECT_NEAR(tau0_computed / 1e6, tau0_expected / 1e6, 0.01);

    // Compute θ₀ at a = a₀ = 0.010 (surface)
    Vector params_a0(1); params_a0(0) = 0.010;
    real_t theta0_computed = law.InitialState(tau0_computed, V_init, sigma_n, eta, params_a0);

    // θ₀ should be in reasonable range
    EXPECT_GT(theta0_computed, 1000.0);
    EXPECT_LT(theta0_computed, 10000.0);
}
```

### 9.4 State Evolution Tests

```cpp
// tests/unit/test_state_evolution.cpp

/// Test 1: Aging law at steady state
TEST(StateEvolution, AgingLawSteadyState) {
    AgingLaw aging;
    real_t Dc = 0.004;

    real_t V = 1e-9;
    real_t theta_ss = aging.SteadyState(V, Dc);

    // At steady state: dθ/dt = 0
    real_t rate = aging.Rate(V, theta_ss, Dc);
    EXPECT_NEAR(rate, 0.0, 1e-20);
}

/// Test 2: Aging law increases θ when V*θ < Dc
TEST(StateEvolution, AgingLawIncreasing) {
    AgingLaw aging;
    real_t Dc = 0.004;

    real_t V = 1e-9;
    real_t theta = 100.0;  // V*θ = 1e-7 < Dc = 0.004

    real_t rate = aging.Rate(V, theta, Dc);
    EXPECT_GT(rate, 0.0);  // θ should increase
}

/// Test 3: Aging law decreases θ when V*θ > Dc
TEST(StateEvolution, AgingLawDecreasing) {
    AgingLaw aging;
    real_t Dc = 0.004;

    real_t V = 1.0;        // 1 m/s (coseismic)
    real_t theta = 1000.0; // V*θ = 1000 >> Dc = 0.004

    real_t rate = aging.Rate(V, theta, Dc);
    EXPECT_LT(rate, 0.0);  // θ should decrease rapidly
}

/// Test 4: Slip law at steady state
TEST(StateEvolution, SlipLawSteadyState) {
    SlipLaw slip;
    real_t Dc = 0.004;

    real_t V = 1e-9;
    real_t theta_ss = slip.SteadyState(V, Dc);

    real_t rate = slip.Rate(V, theta_ss, Dc);
    EXPECT_NEAR(rate, 0.0, 1e-20);
}

/// Test 5: Evolution timescale
TEST(StateEvolution, Timescale) {
    AgingLaw aging;
    real_t Dc = 0.004;
    real_t V = 1e-9;

    // Characteristic timescale τ = Dc/V
    real_t timescale = Dc / V;
    EXPECT_NEAR(timescale, 4e6, 1.0);  // ~46 days
}
```

### 9.5 Traction Computation Tests

```cpp
// tests/unit/test_traction.cpp

/// Test 1: Traction from uniform gradient
TEST(Traction, UniformGradient) {
    // Create simple mesh with known geometry
    // Set u = α*x (linear in x)
    // τ = μ * ∂u/∂x = μ * α

    // ... mesh setup ...

    real_t mu = 32.04e9;
    real_t alpha = 1e-6;  // strain

    // ... set up GridFunction with u = α*x ...
    // ... compute traction ...

    real_t tau_expected = mu * alpha;
    // EXPECT_NEAR(tau_computed, tau_expected, tau_expected * 1e-10);
}

/// Test 2: Traction at fault boundary
TEST(Traction, FaultBoundary) {
    // Test that traction is correctly computed at x=0
    // for a mesh with fault at left boundary
}

/// Test 3: Zero traction for uniform displacement
TEST(Traction, UniformDisplacement) {
    // u = constant everywhere → τ = 0
}
```

### 9.6 Antiplane Operator Tests

```cpp
// tests/unit/test_antiplane.cpp

/// Test 1: Solve Laplace equation with known solution
TEST(AntiplaneDomain, LaplaceSolution) {
    // u = x * y satisfies ∇²u = 0
    // Verify solver recovers this solution with appropriate BCs
}

/// Test 2: Free surface boundary condition
TEST(AntiplaneDomain, FreeSurface) {
    // At z=0: ∂u/∂z = 0 (zero traction)
    // Verify this is satisfied
}

/// Test 3: Far-field boundary condition
TEST(AntiplaneDomain, FarField) {
    // At x → ±∞: u → ±Vp*t/2
    // Verify displacement matches at computational boundary
}

/// Test 4: Method of Manufactured Solutions
TEST(AntiplaneDomain, MMS) {
    // Choose u_exact with non-zero ∇²u
    // Add source term f = -∇²u_exact
    // Verify solver recovers u_exact
}
```

## 10. Integration Tests

### 10.1 Domain-Fault Coupling Test

```cpp
// tests/integration/test_quasi_dynamic.cpp

/// Test 1: Steady-state slip at plate rate
TEST(QuasiDynamic, SteadyStateSlip) {
    // Initialize at steady state
    // Run for short time
    // Verify V ≈ V_init everywhere

    BP2Params params;
    // ... setup ...

    // After short time, slip rate should be near V_init
    real_t V_max = seas_op.GetMaxSlipRate();
    EXPECT_NEAR(V_max, params.V_init, params.V_init * 0.01);
}

/// Test 2: Stress balance maintained
TEST(QuasiDynamic, StressBalance) {
    // At each time step: τ = σn * f(V, θ) + η * V
    // Verify this is satisfied at all fault points
}

/// Test 3: Conservation of slip
TEST(QuasiDynamic, SlipConservation) {
    // Total slip should equal ∫ V dt
}
```

### 10.2 Short BP2 Simulation Test

```cpp
// tests/integration/test_bp2_short.cpp

/// Test: Run BP2 for 1 year, verify qualitative behavior
TEST(BP2Short, OneYear) {
    BP2Params params;
    params.t_final = 3.15e7;  // 1 year in seconds

    // ... run simulation ...

    // Verify:
    // 1. Slip rate starts near V_init
    // 2. Slip rate remains bounded
    // 3. State variable remains positive
    // 4. Stress is in reasonable range
}
```

## 11. Verification Tests

### 11.1 Method of Manufactured Solutions

```cpp
// tests/verification/mms_antiplane.cpp

/// MMS test for antiplane operator
/// Choose exact solution, compute source term, verify convergence
TEST(MMS, AntiplaneConvergence) {
    // u_exact = sin(π*x) * sin(π*z)
    // f = -∇²u_exact = 2π² * sin(π*x) * sin(π*z)

    std::vector<int> N_values = {8, 16, 32, 64};
    std::vector<real_t> errors;

    for (int N : N_values) {
        // Create N×N mesh, solve, compute error
        // ...
        errors.push_back(L2_error);
    }

    // Verify second-order convergence
    for (size_t i = 1; i < errors.size(); i++) {
        real_t ratio = errors[i-1] / errors[i];
        EXPECT_GT(ratio, 3.5);  // Should be ~4 for 2nd order
    }
}
```

### 11.2 BP2 Benchmark Comparison

```cpp
// tests/verification/bp2_benchmark.cpp

/// Compare simulation output with benchmark data
TEST(BP2Benchmark, MatchReferenceData) {
    // Run full BP2 simulation
    BP2Params params;
    // ... run simulation ...

    // Load benchmark data
    auto bench_z0 = LoadBenchmarkData("bp2-qd-z0km-res.txt");

    // Compare at selected times
    std::vector<real_t> compare_times = {
        1e6, 1e7, 1e8, 1e9, 1e10  // seconds
    };

    for (real_t t : compare_times) {
        // Interpolate both to time t
        real_t slip_sim = InterpolateSlip(simulation, t, 0.0);
        real_t slip_bench = InterpolateSlip(bench_z0, t);

        // Allow 5% relative error
        EXPECT_NEAR(slip_sim, slip_bench, std::abs(slip_bench) * 0.05);
    }
}
```

## 12. Parallel Unit Tests

### 12.1 Parallel Domain Operator Tests

```cpp
// tests/parallel/test_parallel_domain.cpp

#include <gtest/gtest.h>
#include <mpi.h>
#include "domain/antiplane_operator.hpp"

/// Test fixture for parallel domain tests
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
    // This test compares parallel solution with serial reference

    // Create identical problems
    BP2Params params;
    real_t mu = 32.04e9;

    // Serial solve (on rank 0 only)
    Vector serial_solution;
    if (rank_ == 0) {
        Mesh serial_mesh = CreateBP2Mesh(params, 32);
        AntiplaneDomainOperator<Mesh> serial_domain(serial_mesh, 1, mu);
        // ... solve ...
        serial_solution = serial_domain.GetSolution();
    }

    // Parallel solve
    Mesh serial_mesh_for_par = CreateBP2Mesh(params, 32);
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh_for_par);
    AntiplaneDomainOperator<ParMesh> par_domain(pmesh, 1, mu);
    // ... solve ...

    // Gather parallel solution to rank 0 and compare
    Vector gathered_solution;
    par_domain.GatherSolutionToRoot(gathered_solution);

    if (rank_ == 0) {
        for (int i = 0; i < serial_solution.Size(); i++) {
            EXPECT_NEAR(serial_solution(i), gathered_solution(i),
                        1e-10 * std::abs(serial_solution(i)));
        }
    }
}

/// Test 3: Traction consistency across ranks
TEST_F(ParallelDomainTest, TractionConsistency) {
    // Verify that fault traction is computed correctly at shared faces
    // (no duplication or missing contributions)
}

/// Test 4: Boundary conditions in parallel
TEST_F(ParallelDomainTest, BoundaryConditions) {
    // Verify that BCs are correctly applied on partitioned boundaries
}
```

### 12.2 Parallel Fault Operator Tests

```cpp
// tests/parallel/test_parallel_fault.cpp

/// Test 1: Fault DOF distribution
TEST_F(ParallelFaultTest, DOFDistribution) {
    // Create mesh and distribute
    ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

    // Create fault data manager
    ParallelFaultData fault_data(pmesh, fault_attr);

    // Verify global count equals sum of local counts
    int local_dofs = fault_data.NumLocalDOFs();
    int global_dofs;
    MPI_Allreduce(&local_dofs, &global_dofs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(global_dofs, fault_data.NumGlobalDOFs());
}

/// Test 2: Local-to-global mapping is unique
TEST_F(ParallelFaultTest, UniqueMapping) {
    ParallelFaultData fault_data(pmesh, fault_attr);

    // Gather all local-to-global maps
    std::set<int> all_global_dofs;
    // ... gather from all ranks ...

    // Each global DOF should be owned by exactly one rank
    EXPECT_EQ(all_global_dofs.size(), fault_data.NumGlobalDOFs());
}

/// Test 3: State consistency after time step
TEST_F(ParallelFaultTest, StateConsistency) {
    // Verify that fault state evolution produces same result
    // as serial execution
}

/// Test 4: Global max slip rate reduction
TEST_F(ParallelFaultTest, GlobalMaxSlipRate) {
    // Set known slip rates on each rank
    // Verify GlobalMaxSlipRate returns correct global max

    // Rank 0: V_max = 1e-9
    // Rank 1: V_max = 1e-6 (should be global max)
    // Rank 2: V_max = 1e-8

    real_t local_V_max = (rank_ == 1) ? 1e-6 : 1e-9;
    fault_op.SetLocalMaxSlipRate(local_V_max);  // Test helper

    real_t global_V_max = fault_op.GetGlobalMaxSlipRate();

    EXPECT_NEAR(global_V_max, 1e-6, 1e-12);
}
```

### 12.3 Parallel Communication Tests

```cpp
// tests/parallel/test_parallel_comm.cpp

/// Test 1: Gather to root
TEST_F(ParallelCommTest, GatherToRoot) {
    ParallelFaultData fault_data(pmesh, fault_attr);

    // Create local data
    Vector local_data(fault_data.NumLocalDOFs());
    for (int i = 0; i < local_data.Size(); i++) {
        // Set value = global DOF index
        local_data(i) = fault_data.LocalToGlobal()[i];
    }

    // Gather to root
    Vector global_data;
    fault_data.GatherToRoot(local_data, global_data);

    // On root, verify all values
    if (rank_ == 0) {
        EXPECT_EQ(global_data.Size(), fault_data.NumGlobalDOFs());
        for (int i = 0; i < global_data.Size(); i++) {
            // Each value should equal its index
            EXPECT_NEAR(global_data(i), i, 1e-12);
        }
    }
}

/// Test 2: Scatter from root
TEST_F(ParallelCommTest, ScatterFromRoot) {
    // Test reverse operation: root distributes data to ranks
}

/// Test 3: Probe interpolation in parallel
TEST_F(ParallelCommTest, ProbeInterpolation) {
    // Verify that output at probe stations is correctly gathered
    // from potentially different ranks
}
```

### 12.4 Scaling Tests

```cpp
// tests/parallel/test_scaling.cpp

/// Strong scaling test: fixed problem size, increasing processors
TEST(Scaling, StrongScaling) {
    // Problem: 128x128 mesh, 1-year simulation
    // Run with np = 1, 2, 4, 8, 16
    // Record wall time
    // Verify reasonable speedup

    // Target: >0.7 parallel efficiency at 8 processes
}

/// Weak scaling test: problem size scales with processors
TEST(Scaling, WeakScaling) {
    // Each rank handles 32x32 mesh elements
    // np = 1: 32x32
    // np = 4: 64x64
    // np = 16: 128x128
    // Verify wall time remains approximately constant

    // Target: <1.5x increase in time at 16 processes
}

/// Communication overhead test
TEST(Scaling, CommunicationOverhead) {
    // Measure time spent in MPI operations
    // Should be <10% of total time for reasonable problem sizes
}
```

### 12.5 Parallel Verification Tests

```cpp
// tests/verification/mms_antiplane_parallel.cpp

/// Parallel MMS test: verify convergence in parallel execution
TEST(ParallelMMS, AntiplaneConvergence) {
    // Same as serial MMS but using ParMesh and ParFiniteElementSpace

    std::vector<int> N_values = {16, 32, 64, 128};
    std::vector<real_t> errors;

    for (int N : N_values) {
        // Create serial mesh, then distribute
        Mesh serial_mesh = Mesh::MakeCartesian2D(N, N, Element::QUADRILATERAL);
        ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

        // Solve with parallel operator
        AntiplaneDomainOperator<ParMesh> domain(pmesh, 1, mu);
        // ... setup MMS source term ...
        // ... solve ...

        // Compute global L2 error
        real_t local_L2 = ComputeLocalL2Error(domain, u_exact);
        real_t global_L2;
        MPI_Allreduce(&local_L2, &global_L2, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        global_L2 = std::sqrt(global_L2);

        if (rank == 0) {
            errors.push_back(global_L2);
        }
    }

    // Verify convergence order (on root)
    if (rank == 0) {
        for (size_t i = 1; i < errors.size(); i++) {
            real_t ratio = errors[i-1] / errors[i];
            EXPECT_GT(ratio, 3.5);  // Second order
        }
    }
}
```

---

# Part V: Implementation Roadmap

## 13. Implementation Strategy: Serial First, Then Parallel

The implementation follows a **serial-first** approach for concept verification, followed by parallel extension:

```
Phase 1-5: Serial Implementation (Complete BP2 in serial)
    ↓
Phase 6: Serial Verification (Match benchmark data)
    ↓
Phase 7: Parallel Infrastructure (Type abstraction, MPI context)
    ↓
Phase 8: Parallel Domain (ParMesh, parallel solver)
    ↓
Phase 9: Parallel Fault & I/O (Distributed fault data)
    ↓
Phase 10: Parallel Verification (Serial/parallel consistency)
```

---

## 14. Phase 1: Core Infrastructure (Serial)

| Task | Deliverable | Tests |
|------|-------------|-------|
| 1.1 Friction law interface | `friction_law.hpp` | `test_friction_law.cpp` |
| 1.2 Dieterich-Ruina implementation | `dieterich_ruina.hpp` | Included above |
| 1.3 State evolution laws | `state_evolution.hpp` | `test_state_evolution.cpp` |
| 1.4 Slip rate solver | In `dieterich_ruina.hpp` | `test_slip_rate_solver.cpp` |
| 1.5 Initial state computation | In `dieterich_ruina.hpp` | `test_initial_state.cpp` |

**Verification checkpoint**: All unit tests pass for friction components.

## 15. Phase 2: Domain Operator (Serial)

| Task | Deliverable | Tests |
|------|-------------|-------|
| 2.1 Abstract domain interface | `domain_operator.hpp` | - |
| 2.2 Antiplane implementation | `antiplane_operator.hpp` | `test_antiplane.cpp` |
| 2.3 Traction computation | In `antiplane_operator.hpp` | `test_traction.cpp` |
| 2.4 MMS verification | - | `mms_antiplane.cpp` |

**Verification checkpoint**: MMS test shows expected convergence order.

## 16. Phase 3: Fault Operator (Serial)

| Task | Deliverable | Tests |
|------|-------------|-------|
| 3.1 Fault geometry | `fault_geometry.hpp` | - |
| 3.2 Rate-state fault operator | `rate_state_fault.hpp` | - |
| 3.3 BP2 parameters | `bp2_params.hpp` | `test_initial_state.cpp` (BP2 specific) |

**Verification checkpoint**: Initial state matches benchmark data.

## 17. Phase 4: SEAS Operator (Serial)

| Task | Deliverable | Tests |
|------|-------------|-------|
| 4.1 Quasi-dynamic operator | `seas_operator.hpp` | `test_quasi_dynamic.cpp` |
| 4.2 Time stepper integration | `time_stepper.hpp` | - |
| 4.3 Short simulation test | - | `test_bp2_short.cpp` |

**Verification checkpoint**: 1-year simulation runs without crash, output is qualitatively correct.

## 18. Phase 5: I/O and Serial Validation

| Task | Deliverable | Tests |
|------|-------------|-------|
| 5.1 SCEC output format | `benchmark_output.hpp` | - |
| 5.2 Probe output | `probe_output.hpp` | - |
| 5.3 Full BP2 run (serial) | - | `bp2_benchmark.cpp` |
| 5.4 Comparison scripts | `compare_results.py` | - |

**Validation checkpoint**: Serial output matches benchmark data within tolerance.

---

## 19. Phase 6: Parallel Infrastructure

| Task | Deliverable | Tests |
|------|-------------|-------|
| 6.1 Type abstraction | `common/seas_types.hpp` | - |
| 6.2 MPI context wrapper | `common/mpi_context.hpp` | `test_mpi_context.cpp` |
| 6.3 Parallel utilities | `common/parallel_utils.hpp` | `test_parallel_utils.cpp` |
| 6.4 CMake parallel config | `CMakeLists.txt` updates | Build test |

**Verification checkpoint**: Code compiles with and without `MFEM_USE_MPI`.

## 20. Phase 7: Parallel Domain Operator

| Task | Deliverable | Tests |
|------|-------------|-------|
| 7.1 Template domain operator | `domain_operator.hpp` (templated) | - |
| 7.2 Parallel antiplane impl. | `antiplane_operator_impl.hpp` | `test_parallel_domain.cpp` |
| 7.3 HYPRE solver setup | In antiplane impl. | - |
| 7.4 Parallel MMS test | - | `mms_antiplane_parallel.cpp` |

**Verification checkpoint**: Parallel MMS shows same convergence as serial.

## 21. Phase 8: Parallel Fault and I/O

| Task | Deliverable | Tests |
|------|-------------|-------|
| 8.1 Parallel fault data | `fault/parallel_fault_data.hpp` | `test_parallel_fault.cpp` |
| 8.2 Distributed state mgmt | In `rate_state_fault.hpp` | - |
| 8.3 Parallel time stepper | `solver/parallel_time_stepper.hpp` | - |
| 8.4 Parallel probe output | `io/parallel_probe_output.hpp` | `test_parallel_comm.cpp` |
| 8.5 Parallel main driver | `pseas.cpp` | - |

**Verification checkpoint**: Parallel code runs on multiple ranks.

## 22. Phase 9: Parallel Verification and Scaling

| Task | Deliverable | Tests |
|------|-------------|-------|
| 9.1 Serial-parallel consistency | - | `test_serial_parallel_consistency.cpp` |
| 9.2 BP2 parallel validation | - | `bp2_benchmark_parallel.cpp` |
| 9.3 Strong scaling test | - | `test_scaling.cpp` |
| 9.4 Weak scaling test | - | Included above |
| 9.5 Communication profiling | - | Timing reports |

**Verification checkpoint**:
- Parallel results match serial within tolerance
- Reasonable parallel efficiency (>70% at 8 ranks)

## 23. Phase 10: Extensions and Optimization

| Task | Deliverable | Tests |
|------|-------------|-------|
| 10.1 Plane strain operator | `plane_strain_operator.hpp` | - |
| 10.2 BP1 configuration | `bp1_params.hpp` | - |
| 10.3 Documentation | README, examples | - |
| 10.4 Performance optimization | Profiling, tuning | Timing tests |
| 10.5 Large-scale testing | - | 3D prototype tests |

---

# Part VI: Acceptance Criteria

## 24. Serial Unit Test Criteria

| Component | Required Tests | Pass Criteria |
|-----------|----------------|---------------|
| Friction coefficient | 5+ tests | All pass with tolerance 1e-10 |
| Slip rate solver | 4+ tests | Convergence in <50 iterations |
| Initial state | 3+ tests | Stress balance within 1e-10 |
| State evolution | 5+ tests | All pass |
| Antiplane operator | 4+ tests | MMS shows 2nd order convergence |

## 25. Serial Integration Test Criteria

| Test | Pass Criteria |
|------|---------------|
| Steady-state | V stays within 1% of V_init for 1 year |
| Stress balance | τ = F(V,θ) + ηV within 0.1% at all times |
| Short BP2 | Completes without crash, output in expected ranges |

## 26. Serial Validation Criteria

| Metric | BP2 Benchmark | Tolerance |
|--------|---------------|-----------|
| Initial τ | 26.546 MPa | ±0.01 MPa |
| Initial θ (z=0) | 10^3.602 s | ±10% |
| Slip at t=100 years | ~3 m | ±20% |
| First earthquake time | ~90 years | ±10 years |

## 27. Parallel Test Criteria

### 27.1 Parallel Unit Tests

| Test | Pass Criteria |
|------|---------------|
| Mesh distribution | All ranks have elements, global count preserved |
| DOF distribution | Unique DOF ownership, no gaps or duplicates |
| Parallel solver | Converges to same tolerance as serial |
| Communication | Gather/scatter operations preserve data |

### 27.2 Parallel Consistency Tests

| Test | Pass Criteria |
|------|---------------|
| Solution match | Parallel solution matches serial within 1e-10 relative error |
| Traction match | Parallel traction matches serial within 1e-10 |
| State evolution | Parallel state matches serial after N time steps |
| Output match | Probe values identical between serial and parallel |

### 27.3 Parallel Scaling Criteria

| Metric | Target |
|--------|--------|
| Strong scaling efficiency (8 ranks) | > 70% |
| Weak scaling overhead (16 ranks) | < 50% increase in time |
| Communication overhead | < 10% of total time |
| Memory per rank | Decreases linearly with rank count |

### 27.4 Parallel Validation Criteria

| Metric | Pass Criteria |
|--------|---------------|
| BP2 parallel results | Match serial results exactly |
| BP2 parallel vs benchmark | Same tolerance as serial (above) |
| Reproducibility | Same result with different rank counts (1,2,4,8) |

---

# Appendix A: BP2 Quick Reference

## A.1 Key Formulas

**Pre-stress**:
```
τ⁰ = σn·amax·sinh⁻¹[(V_init/2V₀)·exp((f₀+b·ln(V₀/V_init))/amax)] + η·V_init
   ≈ 26.546 MPa
```

**Initial state**:
```
θ(z,0) = (Dc/V₀)·exp{(a(z)/b)·ln[(2V₀/V_init)·sinh((τ⁰-η·V_init)/(a(z)·σn))] - f₀/b}
```

**Radiation damping**:
```
η = μ/(2cs) = 32.04e9 / (2 × 3464) ≈ 4.625e6 Pa·s/m
```

**Depth profile of a**:
```
a(z) = {
    0.010,                              z < 15 km
    0.010 + 0.005·(z-15)/3,            15 km ≤ z < 18 km
    0.025,                              z ≥ 18 km
}
```

## A.2 Expected Behavior

- **Interseismic**: V ≈ V_init = 10⁻⁹ m/s, θ slowly increasing
- **Nucleation**: V accelerates, θ drops rapidly
- **Coseismic**: V ≈ 1 m/s, θ very small
- **Postseismic**: V decays, θ recovers
- **Cycle period**: ~90 years between large events

---

# Appendix B: Serial/Parallel Quick Reference

## B.1 Type Mapping Summary

| Serial | Parallel | Header |
|--------|----------|--------|
| `Mesh` | `ParMesh` | `mfem.hpp` |
| `FiniteElementSpace` | `ParFiniteElementSpace` | `mfem.hpp` |
| `BilinearForm` | `ParBilinearForm` | `mfem.hpp` |
| `LinearForm` | `ParLinearForm` | `mfem.hpp` |
| `GridFunction` | `ParGridFunction` | `mfem.hpp` |
| `CGSolver` | `HyprePCG` | `mfem.hpp` |
| `GSSmoother` | `HypreBoomerAMG` | `mfem.hpp` |

## B.2 Key Parallel Patterns

**Mesh Distribution:**
```cpp
// Serial mesh creation, then parallel distribution
Mesh serial_mesh = Mesh::MakeCartesian2D(nx, nz, ...);
ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);
```

**Parallel Solve:**
```cpp
// Using HYPRE for parallel linear algebra
HypreBoomerAMG amg;
HyprePCG pcg(MPI_COMM_WORLD);
pcg.SetPreconditioner(amg);
pcg.SetOperator(A);  // HypreParMatrix
pcg.Mult(B, X);      // HypreParVector
```

**Global Reduction:**
```cpp
// For CFL condition: need global max slip rate
real_t local_V_max = ComputeLocalMaxSlipRate();
real_t global_V_max;
MPI_Allreduce(&local_V_max, &global_V_max, 1, MPI_DOUBLE, MPI_MAX, comm);
```

**Output Gathering:**
```cpp
// Gather probe data to root for I/O
if (mpi_ctx.IsRoot()) {
    output_file << t << " " << global_slip << " " << global_V << "\n";
}
```

## B.3 Build Commands

```bash
# Serial build
mkdir build-serial && cd build-serial
cmake .. -DMFEM_USE_MPI=OFF
make seas

# Parallel build
mkdir build-parallel && cd build-parallel
cmake .. -DMFEM_USE_MPI=ON
make pseas

# Run serial
./seas -bp bp2 -tf 1e10

# Run parallel (4 ranks)
mpirun -np 4 ./pseas -bp bp2 -tf 1e10
```

## B.4 Testing Commands

```bash
# Serial unit tests
./test_friction_law
./test_antiplane

# Parallel unit tests
mpirun -np 2 ./test_parallel_domain
mpirun -np 4 ./test_parallel_domain

# Scaling test
for np in 1 2 4 8; do
    mpirun -np $np ./test_scaling
done
```

---

*Document Version: 2.0*
*Updated: February 5, 2026*
*Changes: Added MPI parallelism support, serial/parallel type abstraction, parallel test strategy*
*For use with: fault_interface_implementation_guide.md, tandem_to_mfem_porting_report.md*
