# SEAS-MFEM Robustness & Refactoring Plan (Draft v1)

**Date:** 2026-04-08
**Status:** First draft for discussion
**Context:** BP5 benchmark results match Tandem at p1/1000m mesh. The codebase is now functionally correct and ready for structural improvements to support long-term development, multi-benchmark generality, and production robustness.

---

## Table of Contents

1. [Current Codebase Summary](#1-current-codebase-summary)
2. [MPI Safety & Parallel Robustness](#2-mpi-safety--parallel-robustness)
3. [File Decomposition: Breaking Up Large Files](#3-file-decomposition-breaking-up-large-files)
4. [Generalizing Beyond BP5](#4-generalizing-beyond-bp5)
5. [Parameter File System](#5-parameter-file-system)
6. [Project Organization & Build System](#6-project-organization--build-system)
7. [Additional Improvements](#7-additional-improvements)
8. [Phased Execution Plan](#8-phased-execution-plan)

---

## 1. Current Codebase Summary

### 1.1 Source Code Inventory

| Directory | Files | LOC | Purpose |
|-----------|-------|-----|---------|
| `domain/` | 6 | 10,687 | Domain PDE operators (elasticity, antiplane) |
| `io/` | 7 | 2,977 | Benchmark output, ParaView, checkpoints |
| `fault/` | 5 | 2,837 | Fault mechanics, geometry, basis functions |
| `integrator/` | 4 | 2,608 | DG bilinear form integrators (IP, BR2) |
| `solver/` | 3 | 1,756 | Time stepping, SEAS quasi-dynamic operator |
| `friction/` | 3 | 1,046 | Dieterich-Ruina friction, state evolution |
| `config/` | 3 | 747 | Benchmark parameter structs (BP1, BP2, BP5) |
| `trace/` | 1 | 670 | Face trace logging |
| `common/` | 3 | 512 | MPI context, parallel utils, type aliases |
| **Total** | **35** | **~24,000** | |

Tests: ~37,000 LOC across 46 files (unit, parallel, verification).

### 1.2 Largest Files (Refactoring Targets)

| File | LOC | Concern |
|------|-----|---------|
| `domain/elasticity_operator.hpp` | 6,505 | Single monolithic class — setup, assembly, solve, traction, debug, verification |
| `fault/rate_state_fault.hpp` | 964 | Moderate but growing — BP2 and BP5 constructors in one class |
| `integrator/dg_elasticity_ip_combined_integrator.hpp` | 919 | Complex but self-contained |
| `io/bp5_benchmark_output.hpp` | 897 | BP5-specific output logic |
| `integrator/dg_elasticity_br2_integrator.hpp` | 826 | Complex but self-contained |

### 1.3 What Already Works Well

- **`DomainOperator<MeshType>` base class** (`domain_operator.hpp`, 217 LOC): Clean abstract interface with `Solve()`, `ComputeTraction()`, mesh-type templating. This is a strong foundation.
- **Template mesh duality**: `FESpaceForMesh<>`, `GridFunctionForMesh<>` in `seas_types.hpp` supports serial/parallel transparently.
- **`SEASQuasiDynamicOperator`**: Domain-agnostic coupling — only talks to `DomainOperator` and `RateStateFaultOperator` through `Vector`-based interfaces.
- **Friction layer**: Fully physics-agnostic. `DieterichRuinaFriction` operates on scalar `V`, `theta`, `a` etc. BP5 feeds `||V||` into the same friction law.
- **DG integrators**: Well-encapsulated IP and BR2 implementations.

---

## 2. MPI Safety & Parallel Robustness

### 2.1 Current State

The codebase has **~80+ direct MPI calls** spread across 10+ files. An `MPIContext` RAII wrapper exists (`common/mpi_context.hpp`, 206 LOC) but is limited to basic reductions and broadcasts on `MPI_COMM_WORLD`.

**Key findings from MPI audit:**

| Pattern | Count | Status |
|---------|-------|--------|
| `MPI_Comm_rank/size` | ~30 | Mixed: some use `mesh_.GetComm()`, some use `MPI_COMM_WORLD` |
| `MPI_Allreduce` | ~25 | Most use `mesh_.GetComm()` (good), some use `MPI_COMM_WORLD` (legacy) |
| `MPI_Allgather/Allgatherv` | ~10 | Used in face-key gathering, fault layout — all use `mesh_.GetComm()` |
| `MPI_Irecv/Isend/Waitall` | 3 | In `ExpandOwnedToLocalFault()` — hardcoded tag `27183` |
| `MPI_Reduce` | ~10 | Root-only reductions for logging/norms |
| `MPI_Gatherv` | ~5 | In `parallel_utils.hpp`, `fault_geometry.hpp` |
| `MPI_Bcast` | ~6 | In `MPIContext` — all on `MPI_COMM_WORLD` |
| `MPI_Barrier` | 2 | In `MPIContext::Barrier()` and `paraview_output.hpp` |

### 2.2 Issues to Address

#### 2.2.1 Inconsistent Communicator Usage

`MPIContext` hardcodes `MPI_COMM_WORLD` for all operations, but domain operations correctly use `mesh_.GetComm()`. If we ever need subcommunicator support (e.g., mesh splitting, ensemble runs), the `MPIContext` wrapper becomes a liability.

**Proposed fix:** `MPIContext` should accept a communicator at construction and store it, defaulting to `MPI_COMM_WORLD` only when none is provided. All reduction/broadcast methods should use the stored communicator.

```cpp
// Current:
MPIContext(int *argc, char ***argv);
MPI_Comm GetComm() const { return MPI_COMM_WORLD; }

// Proposed:
MPIContext(int *argc, char ***argv, MPI_Comm comm = MPI_COMM_WORLD);
MPI_Comm GetComm() const { return comm_; }
```

#### 2.2.2 No MPI Error Checking

Zero MPI return values are checked anywhere. While MFEM internally handles many MPI errors, our raw calls (especially `MPI_Irecv`/`MPI_Isend` in fault communication) should at minimum assert success in debug builds.

**Proposed fix:** Add an `MFEM_SEAS_MPI_CHECK(call)` macro:

```cpp
#ifdef MFEM_DEBUG
#define MFEM_SEAS_MPI_CHECK(call) \
   do { int err = (call); MFEM_ASSERT(err == MPI_SUCCESS, \
        "MPI call failed: " #call " error=" << err); } while(0)
#else
#define MFEM_SEAS_MPI_CHECK(call) (call)
#endif
```

#### 2.2.3 Hardcoded MPI Message Tag

The fault face exchange in `elasticity_operator.hpp:5105-5130` uses hardcoded tag `27183`. This is safe in isolation but fragile if additional point-to-point communication is added later.

**Proposed fix:** Define named tag constants:

```cpp
namespace mpi_tags {
   constexpr int kFaultFaceExchange = 27183;
   // Future tags here...
}
```

#### 2.2.4 Conditional Collective Participation Risk

Several output routines (e.g., `parallel_benchmark_output.hpp`) make write decisions based on a value that the caller must ensure is globally reduced. If a future caller passes a local (not global) value, ranks may disagree on whether to participate in a subsequent collective, causing deadlock.

**Proposed fix:** Defensive `MPI_Allreduce` synchronization before any conditional collective. The existing fix in `parallel_benchmark_output.hpp:84-88` is the right pattern — audit all similar sites.

#### 2.2.5 `parallel_utils.hpp` API Assumptions

`GatherVectorToRoot()` and `ScatterVectorFromRoot()` assume the caller provides correctly sized buffers and matching counts. No runtime validation exists.

**Proposed fix:** Add debug-mode assertions for buffer sizes and count consistency.

### 2.3 Recommended MPI Refactoring Strategy

1. **Extend `MPIContext`** to accept and store a communicator
2. **Consolidate raw MPI calls** behind `MPIContext` methods where possible (especially reductions, gathers, broadcasts)
3. **Add `MFEM_SEAS_MPI_CHECK` macro** for debug-mode error checking
4. **Define tag constants** in a central header
5. **Audit all conditional collectives** for rank-consistent participation
6. **Consider a `FaultCommunicator` class** that encapsulates the Irecv/Isend/Waitall pattern for fault face exchange, hiding the raw MPI from `elasticity_operator.hpp`

---

## 3. File Decomposition: Breaking Up Large Files

### 3.1 `elasticity_operator.hpp` (6,505 lines) — The Primary Target

This single file contains the entire `ElasticityDomainOperator` class: setup, face classification, stiffness assembly, slip RHS assembly, Dirichlet loading, traction computation, fault coordinate extraction, ghost DOF communication, 800 lines of debug utilities, and 400 lines of production verification routines.

**Proposed decomposition into 6 files:**

#### File 1: `elasticity_operator.hpp` (Core — ~1,500 lines)
- Class declaration with public API
- Member variable declarations
- `Solve()` method (main entry point)
- Helper structs: `FaceVertexKey`, `SharedFaultFaceBlock`, `SharedFaultCommBlock`
- `ComputeSkeletonDirichletSign()` helper
- Includes the other 5 implementation files at the bottom (template class pattern)

#### File 2: `elasticity_operator_setup.inl` (~1,200 lines)
- `SetupFESpace()`
- `SetupBoundaryMarkers()`
- `BuildFacetBCTables()` (lines 1819-2202 — the largest single method)
- `ValidateFacetBCTables()`
- `AllgatherKeys()`
- `SetupFaultInfo()`
- `RunStartupFaceAudit()`
- `BuildOwnedFaultLayout()`
- `SetupSolver()`
- `PrecomputeMassInverse()`

#### File 3: `elasticity_operator_assembly.inl` (~1,400 lines)
- `AssembleStiffness()`
- `AssembleSlipContributionIP()` / `AssembleSlipContributionBR2()` (interior faces)
- `AssembleSlipContributionIPShared()` / `AssembleSlipContributionBR2Shared()` (shared faces)
- `AssembleDirichletLoading()`

#### File 4: `elasticity_operator_traction.inl` (~1,200 lines)
- `ComputeTractionImpl()` (lines 5432-6440 — the core traction computation)
- `ComputeTraction()`, `ComputeTractionComponents()`, `ComputeTractionDiagnostics()` wrappers
- `GetFaultDepths()`, `GetFaultCoords2D()`
- `RestrictToOwnedFault()`, `ExpandOwnedToLocalFault()`
- `BuildSlipAtQuadPoints()`

#### File 5: `elasticity_operator_debug.inl` (~800 lines)
- All `DebugDump*()` functions
- `DebugRank()`, `DebugEnabledForTime()`, `DebugPhaseName()`, `DebugFilePath()`
- `DebugTargetLocalElements()`, `DebugShouldDumpFace()`

#### File 6: `elasticity_operator_verify.inl` (~400 lines)
- `VerifyGhostDOFCommunication()`
- `VerifyDirichletSkipSets()`
- `VerifySharedDirichletPerFace()`
- `VerifyRHSNorms()`

**Why `.inl` files?** Since `ElasticityDomainOperator` is a template class, all method definitions must be visible at the point of instantiation. Using `.inl` (inline implementation) files that are `#include`d at the bottom of the main header is the standard pattern for splitting template class implementations without breaking compilation. This is the same approach used by Eigen, Boost, and other template-heavy C++ libraries.

### 3.2 Other Files to Consider

| File | LOC | Action |
|------|-----|--------|
| `rate_state_fault.hpp` | 964 | **Split** BP2 constructor/logic and BP5 constructor/logic into separate `.inl` files or extract the vector-slip RHS computation |
| `bp5_benchmark_output.hpp` | 897 | Fine as-is — self-contained BP5 output |
| `antiplane_operator.hpp` | 1,903 | Consider extracting assembly/traction into `.inl` files if it continues to grow |
| `dg_elasticity_ip_combined_integrator.hpp` | 919 | Fine as-is — single-purpose integrator |
| `dg_elasticity_br2_integrator.hpp` | 826 | Fine as-is — single-purpose integrator |
| `bp5_verification_full.cpp` | ~1,980 | **Refactor** into smaller functions; extract mesh creation, command-line parsing, and time-loop into reusable pieces (see Section 5) |
| `face_trace_logger.hpp` | 670 | Fine as-is |

### 3.3 Include Pattern for Template Decomposition

```cpp
// elasticity_operator.hpp — end of file
// Implementation files for template methods (must be included in header)
#include "elasticity_operator_setup.inl"
#include "elasticity_operator_assembly.inl"
#include "elasticity_operator_traction.inl"
#include "elasticity_operator_debug.inl"
#include "elasticity_operator_verify.inl"
```

Each `.inl` file opens the same namespace and contains only method definitions:

```cpp
// elasticity_operator_setup.inl
namespace mfem { namespace seas {

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::SetupFESpace() { ... }

// ...
}} // namespace mfem::seas
```

---

## 4. Generalizing Beyond BP5

### 4.1 BP5-Specific Code Identified

The following components contain BP5-specific assumptions that prevent direct reuse for other benchmark problems:

#### 4.1.1 Geometry Assumptions in `elasticity_operator.hpp`

| Item | Location | Assumption |
|------|----------|------------|
| Fault plane | `BuildFacetBCTables()` | Fault is Physical Surface 3 (Tandem convention for Y=0 plane) |
| Boundary attrs | `SetupBoundaryMarkers()` | Attr 1 = Natural, Attr 5 = Dirichlet (Tandem BP5 `.geo` convention) |
| Face vertex key | `FaceVertexKey` struct | Assumes tetrahedral faces (3 vertices per face) |
| Coordinate system | Throughout | X = along-strike, Y = fault-normal, Z = depth (Tandem's BP5 mapping) |
| Domain half-widths | `AssembleDirichletLoading()` | `Y > 1000` or `Y < -1000` check for far-field half-rate loading |
| Ref normal | `ref_normal_` member | Hardcoded reference for skeleton Dirichlet orientation matching Tandem BP5 |
| 3D displacement | `NumComponents()` returns 3 | Assumes full 3D elasticity |

#### 4.1.2 Parameter Structures in `config/bp5_params.hpp`

All friction parameters (`a0`, `amax`, `b`, `L0`, `L_nuc`, `V0`, `f0`), geometric parameters (`Wf`, `lf`, `hs`, `ht`, `H`, `l_vw`, `w_nuc`), and stress parameters (`sigma_n`) are hardcoded to BP5 SCEC benchmark values.

The spatial distribution functions `a_of_x2_x3()`, `L_of_x2_x3()`, `V_init_vec()`, `tau0_vec()`, `psi_init()` encode the BP5-specific zoning geometry with hardcoded transition depths.

#### 4.1.3 Output in `io/bp5_benchmark_output.hpp` and `io/bp5_parallel_output.hpp`

- Station names and coordinates are BP5-specific (10 on-fault stations)
- SCEC 8-column output format with BP5 field ordering
- `Probe2DInterpolator` is somewhat general but constructed with BP5-specific stations in the driver

#### 4.1.4 Driver in `tests/verification/bp5_verification_full.cpp`

- ~50 command-line flags with BP5-specific defaults
- Inline mesh creation function `CreateBP5InlineMesh()` with BP5 domain sizes
- BP5-specific initialization sequence
- Hardcoded probe station coordinates

### 4.2 Generalization Strategy

The goal is to separate **framework code** (reusable for any SEAS problem) from **benchmark-specific code** (BP1, BP2, BP5, future problems).

#### Layer 1: Framework (Benchmark-Agnostic)

These should have **zero knowledge** of any specific benchmark:

```
common/          - MPI context, type aliases, parallel utilities
domain/          - DomainOperator base, ElasticityDomainOperator, AntiplaneDomainOperator
fault/           - RateStateFaultOperator, FaultGeometry, FaultBasis
friction/        - FrictionLaw, DieterichRuinaFriction, StateEvolution
integrator/      - DG integrators (IP, BR2)
solver/          - SEASQuasiDynamicOperator, TimeStepper
io/              - GenericBenchmarkOutput, ParaViewOutput, Checkpoint, ProbeOutput
```

**Key changes needed:**
- `ElasticityDomainOperator` constructor should accept boundary attribute configuration (which attrs are Dirichlet, which are Natural, which are Fault) rather than assuming Tandem's BP5 convention
- `SetupBoundaryMarkers()` should be driven by a configuration object, not hardcoded attr numbers
- `AssembleDirichletLoading()` should accept a user-provided displacement function `u_D(x, t)` rather than hardcoding `(Vp*t, 0, 0) * sign`
- `FaceVertexKey` should support both triangular (3-vertex) and quadrilateral (4-vertex) faces
- Reference normal and skeleton Dirichlet sign logic should be configurable

#### Layer 2: Benchmark Scenarios

Each benchmark defines its specific physics:

```
scenarios/
├── scenario_base.hpp       - Abstract scenario interface
├── bp1_scenario.hpp        - BP1: 1D fault, antiplane, depth-dependent a(z)
├── bp2_scenario.hpp        - BP2: 1D fault, antiplane, depth-dependent a(z)
├── bp5_scenario.hpp        - BP5: 2D fault, 3D elasticity, a(x2,x3)
└── custom_scenario.hpp     - User-defined from parameter file
```

A scenario provides:
- Material properties: `mu`, `lambda`, `rho`
- Friction parameters as functions of fault coordinates: `a(x)`, `b(x)`, `L(x)`, `sigma_n(x)`
- Boundary condition specification: which attributes are Dirichlet/Natural/Fault
- Dirichlet displacement function: `u_D(x, t)`
- Initial conditions: `V_init(x)`, `psi_init(x)`, `tau0(x)`
- Domain geometry description
- Output station locations and format

#### Layer 3: Drivers

Thin main programs that wire scenario + framework:

```
drivers/
├── seas_driver.cpp         - Generic driver (reads param file, selects scenario)
├── bp5_driver.cpp          - BP5-specific driver (for backward compatibility / testing)
└── ...
```

### 4.3 Specific Generalization Tasks

#### Task G1: Abstract Boundary Configuration

Create a `BoundaryConfig` struct that replaces hardcoded attribute numbers:

```cpp
struct BoundaryConfig {
   std::set<int> dirichlet_attrs;   // Far-field loading boundaries
   std::set<int> natural_attrs;     // Free surface / natural BC boundaries
   int fault_phys_group;            // Physical group ID for the fault surface
   // Dirichlet displacement function: u_D(x, t) -> Vector
   std::function<void(const Vector &x, real_t t, Vector &u_D)> dirichlet_func;
};
```

#### Task G2: Abstract Fault Zone Configuration

Replace the hardcoded `a_of_x2_x3()`, `L_of_x2_x3()` etc. with a scenario interface:

```cpp
struct FaultZoneConfig {
   std::function<real_t(const Vector &fault_coords)> a;        // Direct effect
   std::function<real_t(const Vector &fault_coords)> b;        // State effect
   std::function<real_t(const Vector &fault_coords)> L;        // Critical slip distance
   std::function<real_t(const Vector &fault_coords)> sigma_n;  // Normal stress
   std::function<real_t(const Vector &fault_coords)> V_init;   // Initial slip rate
   std::function<real_t(const Vector &fault_coords)> psi_init; // Initial state
   // ... etc.
};
```

#### Task G3: Generic Output Configuration

Replace BP5-specific station coordinates with a configurable station list:

```cpp
struct OutputConfig {
   std::vector<ProbeStation> on_fault_stations;    // (name, fault_coord_1, fault_coord_2)
   std::vector<ProbeStation> off_fault_stations;   // (name, x, y, z)
   std::string output_prefix;
   std::string output_format;  // "scec", "csv", "hdf5"
   int checkpoint_interval;
   // ... time-series output frequency, ParaView settings, etc.
};
```

---

## 5. Parameter File System

### 5.1 Current State

Parameters are currently specified through:
1. **C++ structs** (`bp5_params.hpp`) — requires recompilation to change values
2. **Command-line arguments** (~50 flags in `bp5_verification_full.cpp`) — no validation, no documentation, error-prone
3. **Hardcoded defaults** scattered across driver files

### 5.2 Proposed: TOML-Based Configuration

Following Tandem's proven approach (TOML for structure + Lua for parameter functions), but adapted for our needs. We propose TOML as the primary format because:
- Human-readable and writable
- Well-defined specification (no ambiguity like YAML)
- Mature C++ parsing libraries exist (e.g., `toml11`, header-only)
- Tandem already uses it, so our users are familiar with it

#### Example Parameter File: `bp5_1000m.toml`

```toml
[problem]
type = "seas-qd"           # "seas-qd" (quasi-dynamic), "seas-fd" (fully-dynamic)
benchmark = "bp5"           # Pre-defined scenario, or "custom"
dimension = 3

[mesh]
file = "bp5/mesh/p1_1000m/bp5_1000m.msh"
scale = 1000.0              # km -> m conversion factor
format = "gmsh"

[domain]
order = 1                   # Polynomial order
dg_method = "BR2"           # "BR2" or "IP"
face_basis = "GaussLobatto" # "GaussLobatto" or "ClosedUniform"

[solver]
type = "mumps-blr"          # "cg-amg", "mumps", "mumps-blr", "gmres-amg", "superlu"
blr_tolerance = 1e-10       # MUMPS-BLR low-rank tolerance
check_residual = false

[material]
density = 2670.0            # kg/m^3
shear_wave_speed = 3464.0   # m/s
poisson_ratio = 0.25

[friction]
type = "dieterich-ruina"
evolution = "aging-psi"     # "aging", "slip", "aging-psi"
V0 = 1e-6                  # Reference slip rate [m/s]
f0 = 0.6                   # Reference friction coefficient
b = 0.03                   # State effect parameter
L0 = 0.14                  # Critical slip distance (default) [m]
L_nucleation = 0.13        # Critical slip distance (nucleation zone) [m]
a_vw = 0.004               # Direct effect (velocity-weakening)
a_vs = 0.04                # Direct effect (velocity-strengthening)

[loading]
plate_rate = 1e-9           # Far-field plate velocity [m/s]
V_init = 1e-9               # Initial slip rate [m/s]
sigma_n = 25e6              # Effective normal stress [Pa]

[fault.geometry]
depth = 40e3                # Fault depth [m]
length = 100e3              # Fault length [m]
# Zone dimensions (BP5-specific)
shallow_width = 2e3         # Shallow VS zone width [m]
transition_width = 2e3      # VW-VS transition width [m]
uniform_vw_width = 12e3     # Uniform VW zone depth [m]
uniform_vw_length = 60e3    # Uniform VW zone along-strike [m]
nucleation_width = 12e3     # Nucleation zone width [m]

[nucleation]
mode = "smooth"             # "sharp" (SCEC default) or "smooth" (Gaussian taper)
V_nuc = 0.03                # Nucleation slip rate (SCEC: 0.03, Tandem: 0.01)
delta_tau_factor = 1.0      # Pre-stress perturbation multiplier (SCEC: 1.0, Tandem: 0.0)
taper_width = 2000.0        # Gaussian taper width [m] (only for smooth mode)

[time_stepping]
method = "petsc-rk45"       # "dormand-prince", "petsc-rk45"
t_final_years = 1800.0
atol = 1e-7
rtol = 1e-50                # Pure absolute tolerance
dt_min = 1e-6               # Minimum time step [s]
dt_max_years = 0.5          # Maximum time step [years]
dt_init = 0.01              # Initial time step [s]
max_steps = 10000000

[boundary]
# Boundary attribute -> BC type mapping (from mesh physical groups)
dirichlet = [5]             # Far-field loading boundaries
natural = [1]               # Free surface / bottom
fault = [3]                 # Fault surface (interior)

[output]
directory = "."
prefix = "bp5_full"
checkpoint_interval = 5000
print_interval_years = 10.0
write_every_step = false

[output.paraview]
enabled = false
step_interval = 0           # 0 = time-based
dt = 0.0                    # >0 = fixed time interval [s]

[output.stations]
# On-fault stations: name, along-strike [m], depth [m]
on_fault = [
   { name = "strk+00dp+00", x2 = 0.0,    x3 = 0.0 },
   { name = "strk+00dp+10", x2 = 0.0,    x3 = 10e3 },
   { name = "strk+00dp+22", x2 = 0.0,    x3 = 22.5e3 },
   { name = "strk+16dp+00", x2 = 16e3,   x3 = 0.0 },
   { name = "strk+16dp+10", x2 = 16e3,   x3 = 10e3 },
   { name = "strk+36dp+00", x2 = 36e3,   x3 = 0.0 },
   { name = "strk-16dp+00", x2 = -16e3,  x3 = 0.0 },
   { name = "strk-16dp+10", x2 = -16e3,  x3 = 10e3 },
   { name = "strk-24dp+10", x2 = -24e3,  x3 = 10e3 },
   { name = "strk-36dp+00", x2 = -36e3,  x3 = 0.0 },
]
```

### 5.3 Implementation Design

#### Config Parsing Architecture

```
                     +--------------------+
                     |   seas_config.hpp  |
                     |                    |
                     |  SEASConfig struct |
                     |   ├── MeshConfig   |
                     |   ├── DomainConfig |
                     |   ├── SolverConfig |
                     |   ├── MaterialConfig|
                     |   ├── FrictionConfig|
                     |   ├── FaultConfig  |
                     |   ├── TimeConfig   |
                     |   ├── BoundaryConfig|
                     |   └── OutputConfig |
                     +--------+-----------+
                              |
                     +--------v-----------+
                     | seas_config_parser |
                     |                    |
                     |  ParseTOML(file)   |
                     |  -> SEASConfig     |
                     |                    |
                     |  Also supports:    |
                     |  - CLI overrides   |
                     |  - Validation      |
                     |  - Default filling |
                     +--------------------+
```

#### Key Design Decisions

1. **TOML as primary format** — following Tandem's choice
2. **CLI overrides** — any TOML field can be overridden from command line for scripting/batch runs
3. **Validation at parse time** — check ranges, required fields, cross-field consistency
4. **Benchmark presets** — `benchmark = "bp5"` auto-fills all BP5 defaults; user overrides individual values
5. **Backward compatibility** — existing CLI-only invocations continue to work during transition period
6. **No Lua dependency (initially)** — parameter functions (`a(x)`, `sigma_n(x)`) are encoded as zone specifications in TOML rather than arbitrary Lua functions. This avoids a Lua dependency. Can add Lua support later if needed.

#### Zone-Based Parameter Functions (No Lua Required)

Instead of Lua closures, define spatial zones in TOML:

```toml
[friction.zones]
# Zone definitions for a(x2, x3)
# Format: zone type, bounds, value
[[friction.zones.a]]
type = "rectangle"
x2_min = -30e3
x2_max = 30e3
x3_min = 2e3        # Below shallow zone
x3_max = 14e3       # Above deep zone
value = 0.004       # Velocity-weakening a

[[friction.zones.a]]
type = "default"
value = 0.04        # Velocity-strengthening a (everywhere else)

[[friction.zones.a]]
type = "transition"  # Smooth transition between zones
width = 2e3
```

### 5.4 TOML Library Selection

Options for header-only C++ TOML parsing:

| Library | Header-Only | C++ Standard | Stars | Notes |
|---------|-------------|--------------|-------|-------|
| `toml11` | Yes | C++11/17 | 1.9k | Tandem uses this |
| `toml++` | Yes | C++17 | 1.4k | More modern, better error messages |
| `cpptoml` | Yes | C++11 | 500 | Older, less maintained |

**Recommendation:** `toml11` (matches Tandem, C++11 compatible, proven in our reference implementation).

---

## 6. Project Organization & Build System

### 6.1 Current Structure

```
miniapps/seas/
├── common/          (3 files)
├── config/          (3 files)
├── domain/          (6 files)
├── fault/           (5 files)
├── friction/        (3 files)
├── integrator/      (4 files)
├── solver/          (3 files)
├── io/              (7 files)
├── trace/           (1 file)
├── bp1/             (mesh, data, scripts)
├── bp2/             (mesh, data, scripts)
├── bp5/             (mesh, data, scripts)
├── tests/           (unit, parallel, verification)
├── debug_document/  (78 markdown files)
├── tools/           (mesh generation)
└── scripts/         (Python utilities)
```

### 6.2 Proposed Structure

The current directory layout is actually quite reasonable and does not need radical reorganization. The main improvements are:

```
miniapps/seas/
├── common/              # (keep) MPI, types, parallel utils
│   ├── mpi_context.hpp
│   ├── mpi_tags.hpp         # NEW: centralized MPI tag constants
│   ├── parallel_utils.hpp
│   └── seas_types.hpp
│
├── config/              # (keep + extend) Parameter structs
│   ├── bp1_params.hpp
│   ├── bp2_params.hpp
│   ├── bp5_params.hpp
│   ├── seas_config.hpp      # NEW: unified config struct
│   └── seas_config_parser.hpp  # NEW: TOML parser
│
├── domain/              # (keep + decompose)
│   ├── domain_operator.hpp
│   ├── antiplane_operator.hpp
│   ├── antiplane_bdrload_operator.hpp
│   ├── elasticity_operator.hpp          # Slimmed: declarations + includes .inl
│   ├── elasticity_operator_setup.inl    # NEW: setup/init methods
│   ├── elasticity_operator_assembly.inl # NEW: stiffness/RHS assembly
│   ├── elasticity_operator_traction.inl # NEW: traction computation
│   ├── elasticity_operator_debug.inl    # NEW: debug dump utilities
│   ├── elasticity_operator_verify.inl   # NEW: verification routines
│   ├── boundary_config.hpp              # NEW: boundary attribute configuration
│   ├── bp2_mesh.hpp
│   └── seas_boundary_tags.hpp
│
├── fault/               # (keep)
│   ├── rate_state_fault.hpp
│   ├── fault_geometry.hpp
│   ├── fault_basis.hpp
│   ├── face_quadrature.hpp
│   └── fault_nodes.hpp
│
├── friction/            # (keep)
│   ├── friction_law.hpp
│   ├── dieterich_ruina.hpp
│   └── state_evolution.hpp
│
├── integrator/          # (keep)
│   ├── dg_br2_integrator.hpp
│   ├── dg_elasticity_br2_integrator.hpp
│   ├── dg_elasticity_ip_combined_integrator.hpp
│   └── dg_elasticity_ip_penalty_integrator.hpp
│
├── solver/              # (keep)
│   ├── seas_operator.hpp
│   ├── seas_bdrload_operator.hpp
│   └── time_stepper.hpp
│
├── io/                  # (keep + generalize)
│   ├── benchmark_output.hpp
│   ├── bp5_benchmark_output.hpp
│   ├── bp5_parallel_output.hpp
│   ├── paraview_output.hpp
│   ├── probe_output.hpp
│   ├── checkpoint.hpp
│   └── parallel_benchmark_output.hpp
│
├── scenarios/           # NEW: benchmark-specific scenario definitions
│   ├── scenario_base.hpp
│   ├── bp1_scenario.hpp
│   ├── bp2_scenario.hpp
│   └── bp5_scenario.hpp
│
├── drivers/             # NEW: main program entry points
│   ├── seas_driver.cpp          # Generic driver (param file based)
│   └── bp5_quick_check.cpp      # BP5 smoke test driver
│
├── trace/               # (keep)
│   └── face_trace_logger.hpp
│
├── bp1/                 # (keep) benchmark data/mesh/scripts
├── bp2/                 # (keep) benchmark data/mesh/scripts
├── bp5/                 # (keep) benchmark data/mesh/scripts
├── tests/               # (keep)
├── debug_document/      # (keep)
├── tools/               # (keep)
├── scripts/             # (keep)
├── params/              # NEW: example parameter files
│   ├── bp5_1000m.toml
│   ├── bp5_500m.toml
│   ├── bp2_200m.toml
│   └── bp1_100m.toml
│
├── CMakeLists.txt       # (update for new files)
└── Makefile             # (update for new files)
```

### 6.3 Discussion: MOOSE-Style vs. Current Style

MOOSE uses an extremely granular structure with 76+ subdirectories under `framework/src/` (e.g., `actions/`, `auxkernels/`, `bcs/`, `constraints/`, `dgkernels/`, `executioners/`, `functions/`, `ics/`, `indicators/`, `kernels/`, `materials/`, `meshgenerators/`, `outputs/`, `postprocessors/`, `preconditioners/`, `problems/`, `timeintegrators/`, `transfers/`, `userobjects/`, `utils/`, `vectorpostprocessors/`...).

**This level of granularity is overkill for our project.** MOOSE is a general-purpose multi-physics framework with hundreds of physics modules. Our codebase is a focused SEAS simulation miniapp with ~24,000 LOC of framework code.

**Tandem's approach is closer to what we need:** ~10-15 directories with clear functional separation (localoperator, form, common, io, etc.). Our current structure already mirrors this.

**Recommendation:** Keep the current directory layout. The main structural improvement is the `.inl` decomposition of large files and the addition of `scenarios/`, `drivers/`, and `params/` directories.

### 6.4 Build System Updates

The `CMakeLists.txt` needs updates for:
1. New `.inl` files (no compilation changes — they're included by headers)
2. New driver targets (`seas_driver`, `bp5_quick_check`)
3. TOML library integration (header-only, add to include path)
4. New scenario files

No fundamental build system changes are needed. The current CMake structure with phased test registration is well-organized.

---

## 7. Additional Improvements

### 7.1 Error Handling & Assertions

**Current state:** Mix of `MFEM_ASSERT`, `MFEM_ABORT`, raw `std::cerr` + `MPI_Abort`, and unchecked conditions.

**Proposed improvements:**
- Standardize on `MFEM_ASSERT` for debug-mode invariant checks
- Use `MFEM_ABORT` for unrecoverable errors
- Add `MFEM_VERIFY` for user-facing validation (always checked, not just debug)
- Remove raw `MPI_Abort` calls; use MFEM's error handling which handles MPI cleanup

### 7.2 Logging & Diagnostics

**Current state:** `std::cout` / `std::cerr` with `if (rank == 0)` guards scattered throughout. ~800 lines of debug dump utilities in `elasticity_operator.hpp`.

**Proposed improvements:**
- Consider a simple logging wrapper that handles rank-0 filtering automatically
- Controlled by verbosity level (0=errors only, 1=progress, 2=diagnostics, 3=debug dumps)
- Can be as simple as:

```cpp
// In MPIContext or a new Logger class:
void Log(int level, const std::string &msg) const {
   if (IsRoot() && level <= verbosity_) { std::cout << msg << std::endl; }
}
```

### 7.3 Test Organization

**Current state:** 46 test files, ~37,000 LOC. Well-organized into `unit/`, `parallel/`, `verification/`.

**Proposed improvements:**
- Some unit tests are very large (e.g., `test_elasticity_operator.cpp` at 4,594 LOC). Consider splitting into focused test files.
- Add regression test infrastructure: store reference outputs, auto-compare on CI.
- Consider a test parameter file system that mirrors the production parameter files.

### 7.4 Memory Management

**Current state:** Mix of raw pointers (passed by caller) and `std::unique_ptr`. The `SEASQuasiDynamicOperator` takes raw `DomainOpType*` and `FaultOpType*` — caller owns the objects.

**Proposed improvement:** This is fine as-is. The ownership model is clear (caller creates and destroys). No change needed unless we add dynamic object creation from parameter files, in which case the driver should own the objects via `unique_ptr`.

### 7.5 Documentation

**Current state:** 
- `CODEBASE_GUIDE.md` (76K) — comprehensive but may be outdated
- 78 debug documents tracking BP1-BP5 development history
- Inline comments in code are generally good

**Proposed improvements:**
- Update `CODEBASE_GUIDE.md` after refactoring
- Add a `README.md` under `miniapps/seas/` with quick-start instructions
- Parameter file documentation (what each field means, valid ranges)
- Example parameter files with inline comments (the TOML files themselves serve as documentation)

### 7.6 Code Style Consistency

**Not a priority.** The code is internally consistent within files. Minor style variations between files (e.g., brace placement, comment style) are cosmetic and not worth the churn of reformatting.

---

## 8. Phased Execution Plan

### Phase 0: Preparation (No Code Changes)
- [ ] Finalize this plan through discussion
- [ ] Create a `test/refactoring` branch
- [ ] Establish baseline: all existing tests pass

### Phase 1: MPI Safety (Low Risk, High Value)
**Goal:** Harden MPI layer without changing any physics or behavior.

- [ ] 1a. Extend `MPIContext` to accept and store a communicator
- [ ] 1b. Add `MFEM_SEAS_MPI_CHECK` debug macro
- [ ] 1c. Create `common/mpi_tags.hpp` with named tag constants
- [ ] 1d. Audit conditional collectives for rank-consistent participation
- [ ] 1e. Add debug-mode buffer size assertions in `parallel_utils.hpp`
- [ ] 1f. All existing tests must pass with zero behavior change

### Phase 2: File Decomposition (Low Risk, High Value)
**Goal:** Split `elasticity_operator.hpp` without changing any logic.

- [ ] 2a. Create 5 `.inl` files by extracting method definitions
- [ ] 2b. Update `elasticity_operator.hpp` to include the `.inl` files
- [ ] 2c. Update `CMakeLists.txt` and `Makefile` header dependencies
- [ ] 2d. All existing tests must pass identically (byte-for-byte output)
- [ ] 2e. Consider splitting other large files if Phase 2a goes smoothly

### Phase 3: Parameter File Infrastructure (Medium Risk)
**Goal:** Add TOML config parsing alongside existing CLI interface.

- [ ] 3a. Integrate `toml11` as header-only dependency
- [ ] 3b. Define `SEASConfig` struct hierarchy in `config/seas_config.hpp`
- [ ] 3c. Implement `SEASConfigParser` that reads TOML and fills `SEASConfig`
- [ ] 3d. Create example parameter files for BP5 (1000m, 500m)
- [ ] 3e. Add CLI override support (any TOML field overridable from command line)
- [ ] 3f. Add config validation (ranges, required fields, consistency)
- [ ] 3g. Write a new `seas_driver.cpp` that reads param file and runs simulation
- [ ] 3h. Verify new driver produces identical output to existing `bp5_verification_full.cpp`

### Phase 4: Generalization (Medium-High Risk)
**Goal:** Separate framework from benchmark-specific code.

- [ ] 4a. Define `ScenarioBase` interface
- [ ] 4b. Extract BP5-specific logic from `ElasticityDomainOperator` into configurable callbacks
- [ ] 4c. Create `BoundaryConfig` to replace hardcoded boundary attributes
- [ ] 4d. Make `AssembleDirichletLoading` accept user-provided `u_D(x, t)` function
- [ ] 4e. Create `BP5Scenario` that encapsulates all BP5 assumptions
- [ ] 4f. Create `BP1Scenario` and `BP2Scenario`
- [ ] 4g. Verify all three benchmarks produce identical results through the generic driver

### Phase 5: Polish (Low Risk)
- [ ] 5a. Standardize error handling (MFEM_ASSERT/VERIFY/ABORT)
- [ ] 5b. Add logging wrapper with verbosity levels
- [ ] 5c. Update `CODEBASE_GUIDE.md`
- [ ] 5d. Write user-facing documentation for parameter files
- [ ] 5e. Clean up test organization if needed

---

## Open Questions for Discussion

1. **TOML vs. JSON vs. YAML?** We proposed TOML (Tandem uses it). JSON lacks comments. YAML has parsing ambiguities. Are there other preferences?

2. **Lua for parameter functions?** Tandem uses Lua for spatially varying parameters. We proposed zone-based TOML instead (simpler, no external dependency). Is the zone-based approach sufficient for foreseeable use cases, or will we need arbitrary function evaluation?

3. **Scenario abstraction level?** Should scenarios be C++ classes (compiled in) or fully runtime-configurable from parameter files? C++ scenarios are simpler to implement but require recompilation for new problems. Runtime config is more flexible but harder to implement for arbitrary parameter functions.

4. **Header-only template pattern vs. explicit instantiation?** The current header-only approach (all code in `.hpp`) is simple but slow to compile. Explicit template instantiation (`.cpp` files with `template class ElasticityDomainOperator<ParMesh>;`) would speed compilation but adds maintenance burden. Worth doing?

5. **How aggressive on Phase 4 generalization?** We could do a minimal version (just parameterize boundary attributes and Dirichlet function) or a full scenario system. The minimal version is lower risk and may be sufficient for the next several benchmarks.

6. **Test refactoring scope?** Some test files are very large (4,594 LOC). Is splitting them worth the effort, or should we focus on the framework code?

7. **Backward compatibility period?** How long should the old CLI-only interface (`bp5_verification_full.cpp`) coexist with the new param-file driver? Suggest keeping both until the new driver is fully validated, then deprecating (not removing) the old one.
