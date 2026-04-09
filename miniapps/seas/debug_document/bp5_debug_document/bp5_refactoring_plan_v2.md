# SEAS-MFEM Robustness & Refactoring Plan (v2)

**Date:** 2026-04-08
**Status:** Revised draft incorporating critical review feedback
**Predecessor:** `bp5_refactoring_plan_v1.md`
**Context:** BP5 benchmark results match Tandem at p1/1000m. Codebase is functionally correct. This plan addresses structural improvements for robustness, maintainability, and incremental generalization.

---

## Changes from v1

Key revisions based on critical review:

1. **Added Phase 0 (Regression Infrastructure)** — the review's biggest finding. No refactoring code until we have frozen golden outputs, a regression script with hard pass/fail, and a tagged baseline.
2. **Corrected MPI count** — 120 calls across 15 source files (not "~80+"), 58 in `elasticity_operator.hpp` alone.
3. **Removed `ScenarioBase` class hierarchy** — the review correctly identified this as premature abstraction from N=1 working 3D benchmark. Benchmarks are now config presets, not C++ class hierarchies.
4. **Minimized Phase 4** — reduced to BoundaryConfig parameterization (Task G1) and Dirichlet function injection only. Full generalization deferred until a second 3D benchmark concretely exposes what the abstraction needs.
5. **Added concrete acceptance criteria** — defined what "identical output" means for each phase (byte-for-byte for Phase 2, L2 < 1e-12 for Phase 3, L2 < 1e-10 for Phase 4).
6. **Added quick-check regression test** — 100-step BP5 run (~minutes) as a CI-gatable check.
7. **Specified per-phase branches** — no single `test/refactoring` branch. Merge to main sequentially.
8. **Value-type configs** — `BoundaryConfig` and `FaultZoneConfig` are value types with no reference captures to avoid lifetime bugs.
9. **Backward-compatible constructors** — Phase 4 adds overloads, does not remove old signatures.
10. **Noted future dynamic consideration** — `DomainOperator` may need `ExplicitStep()` eventually; don't design for it now but don't close the door.

---

## Table of Contents

0. [Phase 0: Regression Infrastructure (Pre-Requisite)](#0-phase-0-regression-infrastructure)
1. [Current Codebase Summary](#1-current-codebase-summary)
2. [Phase 1: MPI Safety & Parallel Robustness](#2-phase-1-mpi-safety--parallel-robustness)
3. [Phase 2: File Decomposition](#3-phase-2-file-decomposition)
4. [Phase 3: Parameter File System](#4-phase-3-parameter-file-system)
5. [Phase 4: Minimal Generalization](#5-phase-4-minimal-generalization)
6. [Phase 5: Polish](#6-phase-5-polish)
7. [Project Organization & Build System](#7-project-organization--build-system)
8. [Execution Timeline & Branch Strategy](#8-execution-timeline--branch-strategy)
9. [Open Questions](#9-open-questions)

---

## 0. Phase 0: Regression Infrastructure

**This phase must be completed before any refactoring code is written.**

The review identified the plan's biggest gap: the verification tests are informational, not gating. `bp5_verification_full.cpp`'s `RunComparison()` (line 236) prints data point counts but always returns success — `main()` returns 0 regardless. There is no automated regression check.

### 0.1 Freeze Golden Reference Outputs

Run BP5 at p1/1000m with current code (the verified-correct state) and save:
- All 10 on-fault station time series (`bp5_full_fltst_*.txt`)
- Global output (`bp5_full_global.txt`)
- Checkpoint files at step 5000, 10000 (for restart verification)

Store these in `bp5/benchmark_data/golden_p1_1000m/` alongside the existing Tandem reference data.

### 0.2 Write Regression Comparison Script

Create `scripts/regression_check.py` that:
1. Loads simulation output and golden reference (both SCEC 8-column format)
2. Interpolates to common time points (the golden and new runs may have slightly different adaptive time steps)
3. Computes per-station relative L2 error for each field (slip_strike, slip_dip, log10_V_strike, log10_V_dip, tau_strike, tau_dip, log10_state)
4. Reports pass/fail against a specified tolerance
5. Returns nonzero exit code on failure

```
Usage:
  python scripts/regression_check.py \
    --sim-dir ./output \
    --ref-dir bp5/benchmark_data/golden_p1_1000m \
    --prefix bp5_full \
    --tolerance 1e-12 \
    --stations fltst_strk+00dp+00,fltst_strk+00dp+10,...
```

Tolerance levels by phase:
| Phase | Tolerance | Rationale |
|-------|-----------|-----------|
| Phase 2 (code motion) | byte-for-byte diff | Zero logic changes — output must be identical |
| Phase 3 (TOML config) | L2 < 1e-12 | TOML-parsed floats may differ at ULP level from C++ literals |
| Phase 4 (generalization) | L2 < 1e-10 | Code path changes (e.g., `set::count()` vs `if (attr==5)`) may cause FP reordering |

### 0.3 Quick-Check Regression Test

Create a lightweight test that can gate every PR:
- Run BP5 at p1/1000m for **100 accepted time steps** (not full 1800 years)
- Compare against a frozen 100-step reference
- Takes minutes, not hours
- Add to CMakeLists.txt as a registered test (`seas_bp5_regression_quick`)

This bridges the gap between the fast unit tests (seconds) and the full verification run (hours on HPC).

### 0.4 Tag Pre-Refactoring Baseline

```
git tag v1.0-bp5-verified -m "BP5 p1/1000m matches Tandem. Pre-refactoring baseline."
```

All Phase 1+ branches start from this tag.

### 0.5 Phase 0 Checklist

- [ ] 0a. Run BP5 p1/1000m full simulation, collect all outputs
- [ ] 0b. Store golden outputs in `bp5/benchmark_data/golden_p1_1000m/`
- [ ] 0c. Write `scripts/regression_check.py` with hard pass/fail
- [ ] 0d. Run regression script against golden outputs (self-check: tolerance 0)
- [ ] 0e. Create 100-step golden reference for quick-check
- [ ] 0f. Write quick-check test driver or script
- [ ] 0g. Tag current commit as `v1.0-bp5-verified`
- [ ] 0h. Verify byte-for-byte reproducibility: run the same test twice, diff outputs

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

### 1.2 Dependency Graph (Current)

The current dependency DAG is clean with no cycles:

```
common/ ← config/ ← friction/ ← integrator/ ← domain/ ← fault/ ← solver/ ← io/ ← drivers/tests
```

**This must be preserved.** Framework code (`domain/`, `fault/`, `solver/`) must never depend on config-parsing or scenario code.

### 1.3 What Already Works Well

- **`DomainOperator<MeshType>` base class** (217 LOC): Clean abstract interface. `Solve()`, `ComputeTraction()`, `GetFaultDOFs()` all use `Vector`-based interfaces. This is the critical stable API.
- **Template mesh duality**: Serial/parallel transparent via `FESpaceForMesh<>`.
- **`SEASQuasiDynamicOperator`**: Domain-agnostic coupling. Only talks through `Vector`-based interfaces. Phases 1-4 do not touch this.
- **Friction layer**: Physics-agnostic. Operates on scalar `V`, `theta`, `a`. BP5 feeds `||V||` — same law.
- **DG integrators**: Well-encapsulated, single-purpose files.
- **Ownership model**: Caller owns mesh and creates operators with raw pointers. Operators own their internal FE spaces, solvers, and work buffers. Clear and correct.

### 1.4 Largest Files (Refactoring Targets)

| File | LOC | Issue |
|------|-----|-------|
| `domain/elasticity_operator.hpp` | 6,505 | Monolithic: setup + assembly + solve + traction + debug + verify |
| `fault/rate_state_fault.hpp` | 964 | Growing — BP2 and BP5 constructors coexist |
| `io/bp5_benchmark_output.hpp` | 897 | Self-contained, not urgent |

---

## 2. Phase 1: MPI Safety & Parallel Robustness

**Goal:** Harden MPI layer. Zero behavior change. Low risk, high value.
**Branch:** `refactor/phase1-mpi`
**Regression criterion:** Byte-for-byte output vs golden reference.

### 2.1 MPI Call Inventory (Corrected)

Actual count: **120 MPI calls across 15 source files** (not ~80+):

| File | Count | Primary Operations |
|------|-------|--------------------|
| `domain/elasticity_operator.hpp` | 58 | Allreduce, Allgather/v, Irecv/Isend/Waitall, Reduce, Comm_rank/size |
| `common/parallel_utils.hpp` | 14 | Gather/v, Scatter/v, Bcast, Comm_rank/size |
| `common/mpi_context.hpp` | 12 | Allreduce, Bcast, Barrier, Comm_rank/size |
| `tests/parallel/test_parallel_elasticity.cpp` | 10 | Test infrastructure |
| `io/bp5_parallel_output.hpp` | 5 | Allreduce, Comm_rank |
| `solver/seas_operator.hpp` | 4 | Reduce (norm tracking) |
| `tests/parallel/test_parallel_domain.cpp` | 3 | Test infrastructure |
| `tests/parallel/test_br2_consistency.cpp` | 3 | Test infrastructure |
| `common/seas_types.hpp` | 2 | Comm_rank/size |
| `fault/fault_geometry.hpp` | 2 | Allgather, Gatherv |
| Others (5 files) | 7 | Various |

### 2.2 Changes

#### 2.2a Extend `MPIContext` to Store a Communicator

```cpp
// Current (mpi_context.hpp:51):
MPI_Comm GetComm() const { return MPI_COMM_WORLD; }

// Proposed:
MPIContext(int *argc, char ***argv, MPI_Comm comm = MPI_COMM_WORLD);
MPI_Comm GetComm() const { return comm_; }
private:
   MPI_Comm comm_;
```

All 12 MPI calls inside `MPIContext` switch from `MPI_COMM_WORLD` to `comm_`. Default is `MPI_COMM_WORLD` so all existing callers are unaffected.

#### 2.2b Add Debug-Mode MPI Error Check Macro

```cpp
// New file: common/mpi_check.hpp
#ifdef MFEM_DEBUG
#define MFEM_SEAS_MPI_CHECK(call) \
   do { int _mpi_err = (call); \
        MFEM_ASSERT(_mpi_err == MPI_SUCCESS, \
                    "MPI call failed: " #call " error=" << _mpi_err); \
   } while(0)
#else
#define MFEM_SEAS_MPI_CHECK(call) (call)
#endif
```

Apply to the 3 point-to-point calls (`MPI_Irecv`, `MPI_Isend`, `MPI_Waitall` in `elasticity_operator.hpp:5105-5138`) and to the non-blocking calls in `parallel_utils.hpp`. These are the highest-risk MPI calls — collectives rarely fail silently, but p2p can.

#### 2.2c Create `common/mpi_tags.hpp`

```cpp
namespace mfem { namespace seas { namespace mpi_tags {
   constexpr int kFaultFaceExchange = 27183;
}}}
```

Replace hardcoded `27183` in `elasticity_operator.hpp:5106,5130`.

#### 2.2d Audit Conditional Collectives

The review flagged that `parallel_benchmark_output.hpp` had a conditional collective deadlock risk (fixed in phase9). Audit all remaining sites where a collective call is inside a conditional:

Locations to audit:
- `elasticity_operator.hpp`: All `MPI_Reduce` calls in verification/debug methods (guarded by `first_step_debug_` flag — all ranks enter, so this is safe)
- `io/bp5_parallel_output.hpp`: `MPI_Allreduce` at line 237 (already fixed)
- `solver/seas_operator.hpp`: `MPI_Reduce` at lines 500-506 (all ranks participate — safe)

Document the result of each audit as code comments.

#### 2.2e Debug-Mode Buffer Assertions in `parallel_utils.hpp`

Add `MFEM_ASSERT` checks for:
- `GatherVectorToRoot()`: recv buffer on root has correct total size
- `ScatterVectorFromRoot()`: send counts sum to total vector size
- `BroadcastVectorFromRoot()`: all ranks agree on vector size after broadcast

These compile away in release mode — zero performance cost.

### 2.3 Phase 1 Checklist

- [ ] 1a. Extend `MPIContext` to accept/store communicator (default `MPI_COMM_WORLD`)
- [ ] 1b. Add `common/mpi_check.hpp` with `MFEM_SEAS_MPI_CHECK` macro
- [ ] 1c. Create `common/mpi_tags.hpp` with named tag constants
- [ ] 1d. Apply `MFEM_SEAS_MPI_CHECK` to p2p calls in `elasticity_operator.hpp`
- [ ] 1e. Audit all conditional collectives, document results as comments
- [ ] 1f. Add debug-mode buffer assertions in `parallel_utils.hpp`
- [ ] 1g. Run quick-check regression: byte-for-byte match
- [ ] 1h. Run full test suite: all existing tests pass
- [ ] 1i. Merge `refactor/phase1-mpi` to main

---

## 3. Phase 2: File Decomposition

**Goal:** Split `elasticity_operator.hpp` (6,505 lines) into manageable files. Zero logic changes.
**Branch:** `refactor/phase2-decompose`
**Regression criterion:** Byte-for-byte output vs golden reference (same binary — `.inl` is just code motion).

### 3.1 Pre-Check: Verify Byte-for-Byte Reproducibility

Before starting, run the test twice and diff outputs. If the code is not already deterministic (e.g., due to MPI reduction ordering), we need to know now. If outputs differ, the regression criterion must be relaxed.

### 3.2 Decomposition: `elasticity_operator.hpp` into 6 Files

The class is a template (`ElasticityDomainOperator<MeshType>`), so all method definitions must be visible at instantiation. We use `.inl` files `#include`d at the bottom of the header — the standard pattern for splitting template implementations (Eigen, Boost, etc.).

#### File 1: `elasticity_operator.hpp` (Core — ~1,500 lines, down from 6,505)

**Keeps:**
- Class declaration with all public/private member declarations
- Public API methods (inline/short ones: getters, setters)
- `Solve()` method (main entry point, ~140 lines)
- Helper structs: `FaceVertexKey`, `SharedFaultFaceBlock`, `SharedFaultCommBlock`
- `ComputeSkeletonDirichletSign()` helper
- Enums: `SolverType`, `BCMode`, `FacetBC`, `DebugAssemblePhase`
- `#include` directives for the 5 `.inl` files at the bottom

#### File 2: `elasticity_operator_setup.inl` (~1,200 lines)

**Methods moved here (lines 1763-2961 of current file):**
- `SetupFESpace()` — FE collection + spaces creation
- `SetupBoundaryMarkers()` — boundary attribute classification
- `BuildFacetBCTables()` — single source of truth for face BC classification (~380 lines, the largest single method)
- `ValidateFacetBCTables()` — consistency checking
- `AllgatherKeys()` — MPI helper for face key gathering
- `SetupFaultInfo()` — face-neighbor data, quadrature, fault basis
- `RunStartupFaceAudit()` — MPI face classification verification
- `BuildOwnedFaultLayout()` — parallel fault DOF partitioning (~270 lines)
- `GetFaceBC()`, `GetSharedFaceBC()` — trivial accessors
- `SetupSolver()` — linear solver creation
- `PrecomputeMassInverse()` — BR2 element mass inverse caching

#### File 3: `elasticity_operator_assembly.inl` (~1,400 lines)

**Methods moved here (lines 3063-4903 of current file):**
- `AssembleStiffness()` — K matrix assembly with DG integrators (~410 lines)
- `AssembleSlipContributionIP()` — interior fault faces, IP method
- `AssembleSlipContributionBR2()` — interior fault faces, BR2 method (~270 lines)
- `AssembleSlipContributionIPShared()` — shared fault faces, IP method
- `AssembleSlipContributionBR2Shared()` — shared fault faces, BR2 method (~220 lines)
- `AssembleDirichletLoading()` — boundary RHS from Dirichlet BCs (~810 lines)

#### File 4: `elasticity_operator_traction.inl` (~1,200 lines)

**Methods moved here (lines 4909-6440 and 1558-1607 of current file):**
- `ComputeTractionImpl()` — core traction computation (~1,010 lines)
- `ComputeTraction()`, `ComputeTractionComponents()`, `ComputeTractionDiagnostics()` — public wrappers
- `GetFaultDepths()` — fault depth extraction
- `GetFaultCoords2D()` — fault 2D coordinate extraction
- `RestrictToOwnedFault()` — restrict local to owned fault DOFs
- `ExpandOwnedToLocalFault()` — MPI expand owned DOFs with ghost values
- `BuildSlipAtQuadPoints()` — nodal-to-quadrature slip evaluation

#### File 5: `elasticity_operator_debug.inl` (~800 lines)

**Methods moved here (lines 767-1556 of current file):**
- `DebugRank()`, `DebugEnabledForTime()`, `DebugPhaseName()`, `DebugFilePath()` — helper predicates
- `DebugTargetLocalElements()`, `DebugShouldDumpFace()` — filtering
- `DebugDumpFaceData()` — per-face integration details
- `DebugDumpElementVector()` — full element vectors
- `DebugDumpKContributions()` — stiffness matrix element contributions
- `DebugDumpFaultJumpsLocal()`, `DebugDumpFaultJumps()` — fault displacement jumps
- `DebugDumpFaultTractionLocal()`, `DebugDumpFaultTraction()` — fault tractions

#### File 6: `elasticity_operator_verify.inl` (~400 lines)

**Methods moved here (lines 284-722 of current file):**
- `VerifyGhostDOFCommunication()` — MPI ghost DOF expansion test
- `VerifyDirichletSkipSets()` — Dirichlet boundary skip-set audit
- `VerifySharedDirichletPerFace()` — per-face shared Dirichlet diagnostics
- `VerifyRHSNorms()` — ||b_slip||, ||b_dir||, ||b_total|| verification

### 3.3 Include Pattern

At the bottom of `elasticity_operator.hpp`:

```cpp
// Template method implementations — must be in header for instantiation.
#include "elasticity_operator_setup.inl"
#include "elasticity_operator_assembly.inl"
#include "elasticity_operator_traction.inl"
#include "elasticity_operator_debug.inl"
#include "elasticity_operator_verify.inl"
```

Each `.inl` file:

```cpp
// elasticity_operator_setup.inl
// Part of ElasticityDomainOperator — see elasticity_operator.hpp
namespace mfem { namespace seas {

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::SetupFESpace()
{
   // ... (moved verbatim from elasticity_operator.hpp)
}

// ... more methods ...

}} // namespace mfem::seas
```

### 3.4 Build System Updates

- `CMakeLists.txt`: Add `.inl` files to header dependency lists (they're not compiled separately, but IDEs and build systems should track them)
- `Makefile`: Update header dependency lists

### 3.5 Other Files (Deferred)

| File | LOC | Decision |
|------|-----|----------|
| `rate_state_fault.hpp` | 964 | **Defer** — not yet painful enough. Revisit after Phase 4 if BP2+BP5 constructors diverge further |
| `antiplane_operator.hpp` | 1,903 | **Defer** — self-contained, rarely edited |
| `bp5_verification_full.cpp` | ~1,980 | **Defer** — will be addressed by new driver in Phase 3 |

### 3.6 Phase 2 Checklist

- [ ] 2a. Verify byte-for-byte reproducibility (run test twice, diff outputs)
- [ ] 2b. Create 5 `.inl` files by extracting method definitions (pure cut-paste, no logic changes)
- [ ] 2c. Update `elasticity_operator.hpp` to include the `.inl` files
- [ ] 2d. Update `CMakeLists.txt` and `Makefile` header dependencies
- [ ] 2e. Compile — verify zero warnings or errors
- [ ] 2f. Run quick-check regression: byte-for-byte match (if 2a confirmed this is achievable)
- [ ] 2g. Run full test suite: all existing tests pass
- [ ] 2h. Merge `refactor/phase2-decompose` to main

---

## 4. Phase 3: Parameter File System

**Goal:** Add TOML config parsing alongside existing CLI interface. Old driver preserved.
**Branch:** `refactor/phase3-toml`
**Regression criterion:** L2 < 1e-12 between TOML-driven and CLI-driven outputs (TOML-parsed floats may differ at ULP level from C++ literals).

### 4.1 Format Choice: TOML

Rationale (unchanged from v1):
- Tandem uses it — our users are familiar
- Human-readable, supports comments
- Well-defined specification (no YAML ambiguity)
- `toml11` library: header-only, C++11, proven in Tandem

### 4.2 Config Struct Hierarchy

```cpp
// config/seas_config.hpp
namespace mfem { namespace seas {

struct MeshConfig {
   std::string file;
   real_t scale = 1000.0;         // Default: km -> m
   bool inline_mesh = false;
   int nx = 2, ny = 2, nz = 1;   // Inline mesh element counts
   real_t Lx = 200e3, Ly = 100e3, Lz = 100e3;  // Domain half-sizes
};

struct DomainConfig {
   int order = 1;
   std::string dg_method = "BR2";     // "BR2" or "IP"
   int face_basis_type = BasisType::GaussLobatto;
   real_t penalty_factor = 1.0;
};

struct SolverConfig {
   std::string type = "mumps-blr";
   real_t blr_tolerance = 1e-10;
   bool check_residual = false;
};

struct MaterialConfig {
   real_t density = 2670.0;
   real_t shear_wave_speed = 3464.0;
   real_t poisson_ratio = 0.25;
   // Derived:
   real_t mu() const { return density * shear_wave_speed * shear_wave_speed; }
   real_t lambda() const { return 2.0 * poisson_ratio * mu() / (1.0 - 2.0 * poisson_ratio); }
   real_t eta() const { return mu() / (2.0 * shear_wave_speed); }
};

struct FrictionConfig {
   std::string type = "dieterich-ruina";
   std::string evolution = "aging-psi";
   real_t V0 = 1e-6;
   real_t f0 = 0.6;
   real_t b = 0.03;
   real_t L0 = 0.14;
   real_t L_nucleation = 0.13;
   real_t a_vw = 0.004;
   real_t a_vs = 0.04;
};

struct LoadingConfig {
   real_t plate_rate = 1e-9;
   real_t V_init = 1e-9;
   real_t sigma_n = 25e6;
};

struct FaultGeometryConfig {
   real_t depth = 40e3;
   real_t length = 100e3;
   real_t shallow_width = 2e3;
   real_t transition_width = 2e3;
   real_t uniform_vw_width = 12e3;
   real_t uniform_vw_length = 60e3;
   real_t nucleation_width = 12e3;
};

struct NucleationConfig {
   std::string mode = "smooth";
   real_t V_nuc = 0.03;            // SCEC default
   real_t delta_tau_factor = 1.0;   // SCEC default
   real_t taper_width = 2000.0;
   real_t epsilon = 1e-3;
};

struct TimeSteppingConfig {
   std::string method = "petsc-rk45";
   real_t t_final_years = 1800.0;
   real_t atol = 1e-7;
   real_t rtol = 1e-50;
   real_t dt_min = 1e-6;
   real_t dt_max_years = 0.5;
   real_t dt_init = 0.01;
   int max_steps = 10000000;
};

struct BoundaryAttributeConfig {
   std::set<int> dirichlet_attrs = {5};  // Far-field loading
   std::set<int> natural_attrs = {1};    // Free surface / bottom
   int fault_phys_group = 3;             // Fault interior surface
};

struct StationDef {
   std::string name;
   real_t x2;  // Along-strike coordinate
   real_t x3;  // Depth coordinate
};

struct OutputConfig {
   std::string directory = ".";
   std::string prefix = "bp5_full";
   int checkpoint_interval = 5000;
   real_t print_interval_years = 10.0;
   bool write_every_step = false;
   bool paraview_enabled = false;
   int paraview_step_interval = 0;
   real_t paraview_dt = 0.0;
   std::vector<StationDef> on_fault_stations;
};

// Top-level config
struct SEASConfig {
   std::string benchmark;   // "bp1", "bp2", "bp5", "custom"
   MeshConfig mesh;
   DomainConfig domain;
   SolverConfig solver;
   MaterialConfig material;
   FrictionConfig friction;
   LoadingConfig loading;
   FaultGeometryConfig fault_geometry;
   NucleationConfig nucleation;
   TimeSteppingConfig time_stepping;
   BoundaryAttributeConfig boundary;
   OutputConfig output;
};

}} // namespace mfem::seas
```

**Key design decision (from review):** All config sub-structs are **value types**. No reference members, no reference captures in any `std::function`. When config values are passed to operator constructors, they are copied or the struct is passed by `const &` — but the operators must not store references to config members. This eliminates lifetime coupling between config and operators.

### 4.3 Config Parser

```cpp
// config/seas_config_parser.hpp

/// Parse a TOML file into SEASConfig.
/// If benchmark is specified (e.g., "bp5"), fills all defaults first,
/// then overrides with explicit TOML values.
SEASConfig ParseSEASConfig(const std::string &toml_file);

/// Apply command-line overrides to an existing config.
/// Format: --key value where key is dot-separated (e.g., --material.density 2700)
void ApplyCLIOverrides(SEASConfig &config, int argc, char *argv[]);

/// Validate config: check ranges, required fields, cross-field consistency.
/// Returns empty string on success, error description on failure.
std::string ValidateConfig(const SEASConfig &config);
```

**Benchmark presets:** `benchmark = "bp5"` auto-fills all BP5 defaults from the current `BP5Params` struct. User overrides individual values in the TOML file. This means existing `bp5_params.hpp` continues to serve as the authoritative default value source — the parser calls it to fill defaults.

### 4.4 Example Parameter File

See v1 Section 5.2 for the full `bp5_1000m.toml` example — unchanged.

### 4.5 New Generic Driver

```cpp
// drivers/seas_driver.cpp
//
// Usage:
//   mpirun -np N ./seas_driver config.toml [--override.key value ...]
//
// Reads TOML config, constructs domain/fault/solver, runs simulation.
// Replaces benchmark-specific drivers with a single entry point.
```

The driver:
1. Parses TOML → `SEASConfig`
2. Applies CLI overrides
3. Validates config
4. Constructs mesh, domain operator, fault operator, SEAS operator
5. Runs time integration
6. Writes output

### 4.6 Backward Compatibility

The existing `bp5_verification_full.cpp` driver with its ~50 CLI flags is **not modified or removed**. It continues to work as-is. The new `seas_driver.cpp` is an addition, not a replacement.

Validation step (3h): Run both drivers with identical parameters. Compare outputs with L2 < 1e-12 tolerance.

### 4.7 No Lua Dependency

As proposed in v1: spatial parameter functions (`a(x)`, `L(x)`, etc.) are encoded as zone specifications in TOML, not arbitrary Lua closures. This avoids an external dependency. If we later need arbitrary functions (e.g., for non-standard geometries), Lua can be added incrementally without changing the config struct.

### 4.8 Phase 3 Checklist

- [ ] 3a. Integrate `toml11` as header-only dependency (add to include path, no linking)
- [ ] 3b. Define `SEASConfig` struct hierarchy in `config/seas_config.hpp`
- [ ] 3c. Implement `ParseSEASConfig()` and `ValidateConfig()` in `config/seas_config_parser.hpp`
- [ ] 3d. Implement CLI override support (`ApplyCLIOverrides()`)
- [ ] 3e. Create example parameter files: `params/bp5_1000m.toml`, `params/bp5_500m.toml`
- [ ] 3f. Write `drivers/seas_driver.cpp` that reads param file and runs BP5
- [ ] 3g. Write unit tests for config parsing and validation
- [ ] 3h. Run both drivers, compare outputs: L2 < 1e-12
- [ ] 3i. Run full test suite: all existing tests pass
- [ ] 3j. Merge `refactor/phase3-toml` to main

---

## 5. Phase 4: Minimal Generalization

**Goal:** Parameterize boundary configuration and Dirichlet function. Minimum change that unblocks running different problems through the same operator.
**Branch:** `refactor/phase4-boundary-config`
**Regression criterion:** L2 < 1e-10 between old and new code paths.

### 5.1 Scope Reduction (From v1)

The v1 plan proposed a full `ScenarioBase` class hierarchy with C++ scenario classes for BP1, BP2, and BP5. The review correctly identified this as premature abstraction from N=1 working 3D benchmark:

> Building scenarios from one working example guarantees you'll design the wrong abstraction.

**What we do in Phase 4:**
- Parameterize `BoundaryConfig` (which attrs are Dirichlet/Natural/Fault)
- Make `AssembleDirichletLoading()` accept a user-provided displacement function
- Backward-compatible constructor overload

**What we defer:**
- `ScenarioBase` class hierarchy
- `FaultZoneConfig` with `std::function` callbacks for `a(x)`, `L(x)`, etc.
- Generic scenario factory
- Full runtime problem switching

These are deferred until we have a second 3D benchmark (e.g., BP6, or a custom user problem) that concretely reveals what axis of variation matters.

### 5.2 `BoundaryConfig` (Value Type)

```cpp
// domain/boundary_config.hpp

namespace mfem { namespace seas {

/// Configuration for boundary condition assignment.
///
/// This is a VALUE TYPE — no references, no captured state.
/// Copy freely; pass by const reference to constructors.
struct BoundaryConfig {
   /// Boundary attributes that receive Dirichlet (far-field loading) BCs.
   std::set<int> dirichlet_attrs;

   /// Boundary attributes that receive natural (traction-free) BCs.
   std::set<int> natural_attrs;

   /// Physical group ID for the fault interior surface.
   int fault_phys_group;

   /// Static factory for BP5 defaults (Tandem convention).
   static BoundaryConfig BP5Default() {
      return { {5}, {1}, 3 };
   }
};

}} // namespace mfem::seas
```

### 5.3 Dirichlet Function Injection

```cpp
/// Type alias for Dirichlet displacement function.
/// Signature: (x, t) -> u_D
/// x is the spatial point [dim], t is time, u_D is the displacement [dim].
using DirichletFunc = std::function<void(const Vector &x, real_t t, Vector &u_D)>;
```

**Lifetime safety (from review):** The `DirichletFunc` must be a **self-contained closure**. It must capture all needed parameters **by value**, not by reference. Example:

```cpp
// CORRECT: captures Vp by value
real_t Vp = config.loading.plate_rate;
DirichletFunc bp5_dirichlet = [Vp](const Vector &x, real_t t, Vector &u_D) {
   u_D.SetSize(3);
   real_t sign = (x(1) > 0) ? 1.0 : -1.0;
   u_D(0) = sign * 0.5 * Vp * t;
   u_D(1) = 0.0;
   u_D(2) = 0.0;
};

// WRONG: captures config by reference — lifetime bug
DirichletFunc bad = [&config](const Vector &x, real_t t, Vector &u_D) { ... };
```

### 5.4 Constructor Changes

**Backward-compatible overload** (from review recommendation):

```cpp
// NEW constructor (preferred):
ElasticityDomainOperator(MeshType &mesh, int order,
                         real_t lambda, real_t mu, real_t Vp, real_t Wf, real_t lf,
                         const BoundaryConfig &bdr_config,
                         DirichletFunc dirichlet_func,
                         const std::string &method = "BR2",
                         SolverType solver_type = SolverType::CG_AMG,
                         BCMode bc_mode = BCMode::FarField);

// OLD constructor (deprecated but preserved — fills BoundaryConfig from old params):
ElasticityDomainOperator(MeshType &mesh, int order,
                         real_t lambda, real_t mu, real_t Vp, real_t Wf, real_t lf,
                         const std::string &method = "BR2",
                         SolverType solver_type = SolverType::CG_AMG,
                         BCMode bc_mode = BCMode::FarField)
   : ElasticityDomainOperator(mesh, order, lambda, mu, Vp, Wf, lf,
                               BoundaryConfig::BP5Default(),
                               MakeBP5DirichletFunc(Vp),
                               method, solver_type, bc_mode) {}
```

This lets all 7 verification tests and the old driver continue to compile and work without modification. Tests are migrated incrementally, not in a big bang.

### 5.5 Internal Changes

- `SetupBoundaryMarkers()`: Read from `bdr_config_.dirichlet_attrs` instead of hardcoding attr 5
- `BuildFacetBCTables()`: Use `bdr_config_.fault_phys_group` instead of hardcoding physical group 3
- `AssembleDirichletLoading()`: Call `dirichlet_func_(x, t, u_D)` instead of hardcoding `(Vp*t, 0, 0) * sign`
- `ref_normal_`: Passed through `BoundaryConfig` or computed from mesh geometry rather than hardcoded

### 5.6 What This Enables

With just `BoundaryConfig` + `DirichletFunc`:
- Run BP5 with different mesh boundary attribute conventions (not locked to Tandem's `.geo` numbering)
- Apply different far-field loading functions without recompiling
- Prepare for future benchmarks that use the same 3D elasticity but different boundary setups

What this does NOT do:
- Does not generalize friction parameter spatial distributions (still uses `BP5Params::a_of_x2_x3`)
- Does not create a scenario system
- Does not support arbitrary fault orientations (still Y=0)

These are intentional non-goals for Phase 4. They wait for Phase 4+ when a concrete second benchmark drives the design.

### 5.7 Phase 4 Regression

Before starting Phase 4:
- Create frozen regression reference (if not already done in Phase 0)
- Run BP5 at p1/1000m through old constructor → save output A
- After Phase 4: run through new constructor with `BP5Default()` → save output B
- Compare: L2(A, B) < 1e-10

### 5.8 Phase 4 Checklist

- [ ] 4a. Create `domain/boundary_config.hpp` with `BoundaryConfig` value type
- [ ] 4b. Define `DirichletFunc` type alias
- [ ] 4c. Add new constructor overload to `ElasticityDomainOperator` (in `elasticity_operator.hpp`)
- [ ] 4d. Implement deprecated old constructor as delegating constructor
- [ ] 4e. Refactor `SetupBoundaryMarkers()` to use `bdr_config_`
- [ ] 4f. Refactor `BuildFacetBCTables()` to use `bdr_config_.fault_phys_group`
- [ ] 4g. Refactor `AssembleDirichletLoading()` to use `dirichlet_func_`
- [ ] 4h. Verify all tests compile and pass with old constructor (zero-change path)
- [ ] 4i. Update `seas_driver.cpp` to construct `BoundaryConfig` from TOML
- [ ] 4j. Run regression: L2 < 1e-10 between old and new code paths
- [ ] 4k. Merge `refactor/phase4-boundary-config` to main

---

## 6. Phase 5: Polish

**Goal:** Standardize error handling, add logging, update documentation.
**Branch:** `refactor/phase5-polish`
**Regression criterion:** Byte-for-byte match (no computational changes).

### 6.1 Error Handling Standardization

Current state: Mix of `MFEM_ASSERT`, `MFEM_ABORT`, raw `std::cerr` + `MPI_Abort`, and unchecked conditions.

Changes:
- `MFEM_ASSERT`: Debug-mode invariant checks (internal consistency)
- `MFEM_VERIFY`: Always-on user-facing validation (bad config, missing files)
- `MFEM_ABORT`: Unrecoverable errors
- Remove raw `MPI_Abort` calls — use MFEM's error handling which does MPI cleanup

### 6.2 Logging Wrapper

Simple rank-0 filtering with verbosity levels:

```cpp
// In MPIContext or a standalone Logger:
enum class LogLevel { Error = 0, Progress = 1, Diagnostic = 2, Debug = 3 };

void Log(LogLevel level, const std::string &msg) const {
   if (IsRoot() && static_cast<int>(level) <= verbosity_) {
      std::cout << msg << std::endl;
   }
}
```

Replace scattered `if (rank == 0) std::cout << ...` patterns with `mpi.Log(LogLevel::Progress, ...)`.

### 6.3 Documentation Updates

- Update `CODEBASE_GUIDE.md` to reflect new file structure (`.inl` files, `boundary_config.hpp`, etc.)
- Add `README.md` under `miniapps/seas/` with quick-start instructions
- Parameter file documentation (inline comments in example TOML files serve as primary docs)
- Update sbatch scripts and Frontera job configurations for new driver path

### 6.4 Future Dynamic Solver Consideration

The review noted that `DomainOperator::Solve()` assumes an implicit linear solve. A future fully-dynamic solver would need an explicit time step (mass matrix multiply) rather than a stiffness solve. The `DomainOperator` base class (217 LOC) is small and adding a `virtual void ExplicitStep(...)` later is cheap.

**Action:** No change now. Note this as a comment in `domain_operator.hpp`:

```cpp
// NOTE: This interface assumes implicit/quasi-static solves. If fully-dynamic
// (explicit) time integration is added, consider adding ExplicitStep() or
// MassMatrixMult() to this interface.
```

### 6.5 Phase 5 Checklist

- [ ] 5a. Standardize error handling (MFEM_ASSERT/VERIFY/ABORT)
- [ ] 5b. Remove raw `MPI_Abort` calls
- [ ] 5c. Add logging wrapper with verbosity levels
- [ ] 5d. Replace scattered `if (rank == 0) cout` with `Log()` calls
- [ ] 5e. Update `CODEBASE_GUIDE.md`
- [ ] 5f. Add `README.md` with quick-start
- [ ] 5g. Add future-dynamic comment in `domain_operator.hpp`
- [ ] 5h. Update sbatch scripts for new driver
- [ ] 5i. Run regression: byte-for-byte match
- [ ] 5j. Merge `refactor/phase5-polish` to main

---

## 7. Project Organization & Build System

### 7.1 Final Directory Structure (After All Phases)

```
miniapps/seas/
├── common/                          # MPI, types, parallel utils
│   ├── mpi_context.hpp              # (Phase 1: add stored communicator)
│   ├── mpi_check.hpp               # NEW Phase 1: debug MPI error macro
│   ├── mpi_tags.hpp                # NEW Phase 1: named tag constants
│   ├── parallel_utils.hpp           # (Phase 1: add debug assertions)
│   └── seas_types.hpp
│
├── config/                          # Parameters and configuration
│   ├── bp1_params.hpp
│   ├── bp2_params.hpp
│   ├── bp5_params.hpp
│   ├── seas_config.hpp             # NEW Phase 3: unified config struct
│   └── seas_config_parser.hpp      # NEW Phase 3: TOML parser
│
├── domain/                          # Domain PDE operators
│   ├── domain_operator.hpp          # (Phase 5: add future-dynamic comment)
│   ├── antiplane_operator.hpp
│   ├── antiplane_bdrload_operator.hpp
│   ├── elasticity_operator.hpp      # Phase 2: slimmed to ~1500 LOC
│   ├── elasticity_operator_setup.inl     # NEW Phase 2
│   ├── elasticity_operator_assembly.inl  # NEW Phase 2
│   ├── elasticity_operator_traction.inl  # NEW Phase 2
│   ├── elasticity_operator_debug.inl     # NEW Phase 2
│   ├── elasticity_operator_verify.inl    # NEW Phase 2
│   ├── boundary_config.hpp          # NEW Phase 4
│   ├── bp2_mesh.hpp
│   └── seas_boundary_tags.hpp
│
├── fault/                           # (unchanged)
├── friction/                        # (unchanged)
├── integrator/                      # (unchanged)
├── solver/                          # (unchanged)
├── io/                              # (unchanged)
├── trace/                           # (unchanged)
│
├── drivers/                         # NEW Phase 3
│   └── seas_driver.cpp             # Generic TOML-based driver
│
├── params/                          # NEW Phase 3: example config files
│   ├── bp5_1000m.toml
│   ├── bp5_500m.toml
│   └── bp2_200m.toml
│
├── scripts/
│   ├── regression_check.py         # NEW Phase 0
│   └── compare_bp2_first_cycle.py   # (existing)
│
├── bp1/, bp2/, bp5/                 # (unchanged: mesh, data, scripts)
│   └── bp5/benchmark_data/
│       └── golden_p1_1000m/        # NEW Phase 0: frozen golden outputs
│
├── tests/                           # (unchanged structure, updated deps)
├── debug_document/                  # (unchanged)
├── tools/                           # (unchanged)
├── CMakeLists.txt                   # (updated per phase)
└── Makefile                         # (updated per phase)
```

### 7.2 Dependency Direction (Strict Rule)

```
drivers/ → config/seas_config_parser → config/seas_config → config/bp*_params
                                                          ↘
common/ ← friction/ ← integrator/ ← domain/ ← fault/ ← solver/ ← io/
```

**Never let framework code depend on config-parsing or driver code.** The `domain/` directory includes `domain/boundary_config.hpp` (a simple value-type struct with no dependencies), but never `config/seas_config_parser.hpp`.

### 7.3 MOOSE-Style: Not Adopted

Per v1 analysis: MOOSE's 76-directory granularity is overkill for 24K LOC. Our structure already mirrors Tandem's ~10-15 directory layout, which is appropriate for this scale. No change.

---

## 8. Execution Timeline & Branch Strategy

### 8.1 Branch Strategy (From Review)

Per-phase branches merged sequentially to main. No single long-lived refactoring branch.

```
main ──────────────────────────────────────────────────────────►
  │                                                             
  ├─ tag: v1.0-bp5-verified                                     
  │                                                             
  ├─ refactor/phase0-regression ──── merge ─┐                   
  │                                         │                   
  ├─────────────────────────────────────────►│                   
  │                                         │                   
  ├─ refactor/phase1-mpi ───────── merge ───┤                   
  │                                         │                   
  ├─ refactor/phase2-decompose ── merge ────┤                   
  │                                         │                   
  ├─ refactor/phase3-toml ────── merge ─────┤                   
  │                                         │                   
  ├─ refactor/phase4-boundary ── merge ─────┤                   
  │                                         │                   
  └─ refactor/phase5-polish ─── merge ──────┘                   
```

Each branch starts from the latest main (after the previous phase merged). Never rebase — merge only, so history is recoverable.

### 8.2 Phase Summary Table

| Phase | Branch | Risk | Regression Criterion | Key Deliverable |
|-------|--------|------|---------------------|-----------------|
| 0 | `refactor/phase0-regression` | None | Self-check: tolerance 0 | Golden outputs + regression script |
| 1 | `refactor/phase1-mpi` | Low | Byte-for-byte | Hardened MPI layer |
| 2 | `refactor/phase2-decompose` | Low | Byte-for-byte | `elasticity_operator.hpp` → 6 files |
| 3 | `refactor/phase3-toml` | Medium | L2 < 1e-12 | TOML config + generic driver |
| 4 | `refactor/phase4-boundary` | Medium | L2 < 1e-10 | Configurable boundary attrs + Dirichlet func |
| 5 | `refactor/phase5-polish` | Low | Byte-for-byte | Error handling, logging, docs |

### 8.3 Rollback Strategy

- Phase 0: Nothing to roll back (additive only)
- Phase 1: Trivially reversible (macro removal, tag inlining, constructor default)
- Phase 2: `cat *.inl >> elasticity_operator.hpp` and remove includes
- Phase 3: Old driver preserved. Remove TOML parser, example files, new driver. Low cost.
- Phase 4: Old constructor preserved as overload. Can remove new constructor and `boundary_config.hpp`. Tests unchanged.
- Phase 5: Remove logging wrapper. Revert comment changes. Low cost.

---

## 9. Open Questions

### Resolved from v1 (By Review)

| v1 Question | Resolution |
|-------------|------------|
| TOML vs JSON vs YAML? | **TOML** — confirmed. Tandem precedent, comment support, well-defined spec. |
| Lua for parameter functions? | **No Lua initially** — zone-based TOML. Add Lua later only if needed. |
| Scenario abstraction level? | **Skip ScenarioBase.** Use config structs as scenarios. Benchmarks are presets. |
| How aggressive on Phase 4? | **Minimal.** BoundaryConfig + DirichletFunc only. Wait for second 3D benchmark. |
| Backward compatibility period? | **Indefinite.** Old constructor stays as deprecated overload. Old driver stays as-is. |

### Remaining Open Questions

1. **Header-only templates vs explicit instantiation?** The `.inl` pattern keeps everything header-only (current approach). Explicit template instantiation (`template class ElasticityDomainOperator<ParMesh>;` in a `.cpp`) would halve compile times but adds maintenance. Worth doing after Phase 2 if compile times are painful?

2. **Quick-check test: 100 steps or fewer?** Need to determine empirically how many steps are needed for meaningful regression detection while keeping runtime under 5 minutes on 8 MPI ranks.

3. **Test file splitting scope?** Some test files are large (`test_elasticity_operator.cpp`: 4,594 LOC; `test_parallel_elasticity.cpp`: 3,092 LOC). The review asked if splitting is worth the effort. Suggest deferring unless a specific test file becomes a merge-conflict hotspot.

4. **`FaceVertexKey` for non-tet meshes?** Currently hardcoded for 3 vertices (tetrahedra). Hex meshes need 4 vertices per face. Phase 4 does not address this — it's a future generalization triggered when someone needs hex support. Flag for awareness only.

5. **Config validation depth?** How much cross-field validation should `ValidateConfig()` do? Minimum: required fields present, numeric ranges. Maximum: check mesh file exists, verify boundary attributes match mesh, check friction parameter consistency (e.g., `a < b` in VW zone). Suggest starting with minimum and adding validation as users hit errors.
