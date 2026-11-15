# Phase Field Fracture (PFF) Refactoring Plan

## Overview

This document outlines a plan to improve the maintainability, readability, and extensibility of the PFF miniapp, inspired by MOOSE's architecture.

## Current Issues

1. **Hardcoded parameters in pff.cpp**: Material properties, solver settings, and boundary conditions are all defined in the main driver file
2. **Monolithic solver class**: `pff_solver.hpp` contains ~1200 lines mixing solver logic, residual computation, boundary conditions, and I/O
3. **Problem-specific code mixed with reusable components**: Mode 2 shear loading is hardcoded; adding new problems requires code changes
4. **No input file system**: Users must recompile to change parameters

## Proposed Architecture

### 1. Input File System

Create a simple input file format (JSON or custom) similar to MOOSE's `.i` files:

```
# example: mode2_shear.pff
[Mesh]
  file = ../mesh/mesh_quad.msh
  refine_levels = 0

[Material]
  E = 2.1e5
  nu = 0.3
  Gc = 2.7
  l = 0.02
  eta = 1e-6
  p = 2

[BoundaryConditions]
  [top]
    type = PrescribedDisplacement
    component = x
    function = "t"  # u_x = t
  [bottom]
    type = Fixed
    components = "x y"

[Solver]
  max_staggered_iter = 20
  fp_rel_tol = 1e-8
  fp_abs_tol = 1e-10
  linear_solver = CG
  preconditioner = AMG

[TimeStep]
  dt = 2e-5
  t_final = 0.02

[Output]
  format = ParaView
  interval = 10
  csv = true
```

**Implementation:**
- Create `PFFInputParser` class to read input files
- Use MFEM's existing JSON support or create simple parser
- Store parameters in `PFFProblemConfig` struct

### 2. Modular Code Structure

Reorganize into separate files:

```
miniapps/pff/
├── core/
│   ├── pff_input.hpp          # Input file parsing
│   ├── pff_config.hpp         # Configuration structs
│   ├── pff_material.hpp       # Material parameters (existing)
│   └── pff_utils.hpp          # Common utilities
│
├── physics/
│   ├── elasticity_physics.hpp # Elasticity sub-problem
│   ├── damage_physics.hpp     # Damage sub-problem
│   └── coupling.hpp           # Staggered coupling logic
│
├── boundary/
│   ├── bc_base.hpp            # Base BC class
│   ├── dirichlet_bc.hpp       # Dirichlet BCs
│   ├── neumann_bc.hpp         # Neumann BCs (future)
│   └── bc_factory.hpp         # BC creation from input
│
├── solvers/
│   ├── staggered_solver.hpp   # Staggered iteration driver
│   ├── linear_solver.hpp      # Linear solver wrapper
│   └── convergence.hpp        # Convergence checking
│
├── postprocess/
│   ├── reaction_force.hpp     # Reaction force computation
│   ├── energy_output.hpp      # Energy quantities
│   └── field_output.hpp       # Field visualization
│
├── problems/
│   ├── mode1_tension.hpp      # Mode 1 problem setup
│   ├── mode2_shear.hpp        # Mode 2 problem setup
│   └── custom_problem.hpp     # User-defined problems
│
├── pff_driver.cpp             # Main driver (minimal)
└── CMakeLists.txt
```

### 3. Base Classes and Interfaces

**Physics Interface:**
```cpp
class PhysicsBase {
public:
    virtual void Setup(const PFFConfig& config) = 0;
    virtual real_t ComputeResidual() = 0;
    virtual void Solve() = 0;
    virtual void UpdateFields() = 0;
};

class ElasticityPhysics : public PhysicsBase { ... };
class DamagePhysics : public PhysicsBase { ... };
```

**Boundary Condition Interface:**
```cpp
class BCBase {
public:
    virtual void Apply(ParGridFunction& gf, real_t time) = 0;
    virtual Array<int> GetEssentialDofs() = 0;
};

class DirichletBC : public BCBase { ... };
class FunctionDirichletBC : public BCBase { ... };
```

**Problem Interface:**
```cpp
class ProblemBase {
public:
    virtual void SetupMesh() = 0;
    virtual void SetupBCs() = 0;
    virtual void SetupMaterial() = 0;
    virtual void Run() = 0;
};
```

### 4. Configuration Structs

```cpp
struct MeshConfig {
    std::string file;
    int refine_levels = 0;
    int order = 1;
};

struct MaterialConfig {
    real_t E, nu, Gc, l, eta, p;
    bool spectral_decomposition = true;
};

struct SolverConfig {
    int max_fp_iter = 20;
    real_t fp_rel_tol = 1e-8;
    real_t fp_abs_tol = 1e-10;
    real_t nl_rel_tol = 1e-8;
    real_t nl_abs_tol = 1e-10;
    std::string linear_solver = "CG";
    std::string preconditioner = "AMG";
};

struct TimeConfig {
    real_t dt, t_final;
};

struct OutputConfig {
    std::string format = "ParaView";
    int interval = 10;
    bool csv = true;
    std::string prefix = "pff_output";
};

struct PFFConfig {
    MeshConfig mesh;
    MaterialConfig material;
    SolverConfig solver;
    TimeConfig time;
    OutputConfig output;
    std::vector<BCConfig> bcs;
};
```

### 5. Simplified Main Driver

```cpp
// pff_driver.cpp - entire file
#include "pff.hpp"

int main(int argc, char *argv[]) {
    Mpi::Init(argc, argv);

    // Parse command line for input file
    const char* input_file = "input.pff";
    OptionsParser args(argc, argv);
    args.AddOption(&input_file, "-i", "--input", "Input file");
    args.Parse();

    // Load configuration
    PFFConfig config = PFFInputParser::Parse(input_file);

    // Create and run problem
    auto problem = ProblemFactory::Create(config);
    problem->Run();

    return 0;
}
```

### 6. Reusable Components

Components that should be generic and reusable:

| Component | Current Location | Proposed |
|-----------|-----------------|----------|
| Residual computation | `pff_solver.hpp` | `physics/residual.hpp` |
| Reaction force | `pff_solver.hpp` | `postprocess/reaction_force.hpp` |
| Convergence check | `pff_solver.hpp` | `solvers/convergence.hpp` |
| MOOSE-style output | `pff_solver.hpp` | `core/moose_output.hpp` |
| Staggered iteration | `pff_solver.hpp` | `solvers/staggered_solver.hpp` |
| Degradation functions | `materials/` | Keep, add more models |

### 7. Implementation Phases

**Phase 1: Input File System (1-2 weeks)**
- Create `PFFInputParser` class
- Define `PFFConfig` structs
- Update main driver to use input files
- Backward compatible: keep command-line options

**Phase 2: Modular Solvers (2-3 weeks)**
- Extract `ElasticityPhysics` class
- Extract `DamagePhysics` class
- Create `StaggeredSolver` class
- Add convergence module

**Phase 3: Boundary Conditions (1-2 weeks)**
- Create BC base class and factory
- Implement common BC types
- Support reading BCs from input file

**Phase 4: Problem Templates (1-2 weeks)**
- Create `ProblemBase` interface
- Implement Mode 1 and Mode 2 as derived classes
- Add custom problem support

**Phase 5: Postprocessing (1 week)**
- Modular reaction force computation
- Energy output utilities
- Flexible field output

### 8. Example: Mode 2 Shear with New Architecture

**Input file (mode2_shear.pff):**
```
[Problem]
  type = Mode2Shear

[Mesh]
  file = mesh_quad.msh

[Material]
  E = 2.1e5
  nu = 0.3
  Gc = 2.7
  l = 0.02

[TimeStep]
  dt = 2e-5
  t_final = 0.02

[Output]
  prefix = mode2_results
```

**Run command:**
```bash
mpirun -np 4 pff -i mode2_shear.pff
```

### 9. Benefits

1. **Easier parameter studies**: Change input file, not code
2. **Better maintainability**: Each module has single responsibility
3. **Extensibility**: Add new physics, BCs, or problems without touching core
4. **Testing**: Unit test individual components
5. **Documentation**: Each module can be documented separately
6. **Collaboration**: Multiple developers can work on different modules

### 10. Backward Compatibility

- Keep command-line options working
- If no input file specified, use defaults
- Existing examples continue to work

## Priority

1. **High**: Input file system (most user-facing benefit)
2. **High**: Modular solver structure (maintainability)
3. **Medium**: BC framework (extensibility)
4. **Medium**: Problem templates (code reuse)
5. **Low**: Postprocessing modules (convenience)

## Dependencies

- MFEM's existing config/option parsing
- Consider using nlohmann/json for JSON parsing (header-only)
- Or MFEM's built-in JSON support if available

## Open Questions

1. Use JSON, YAML, or custom format for input files?
2. How to handle mesh-dependent BCs (like "top" boundary)?
3. Should we support MOOSE `.i` format directly for compatibility?
4. How to handle adaptive mesh refinement in the new structure?
