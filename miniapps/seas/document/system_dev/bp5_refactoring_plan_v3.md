# SEAS-MFEM Robustness & Refactoring Plan (v3)

**Date:** 2026-04-09
**Status:** Final draft — ready for execution
**Predecessors:** `bp5_refactoring_plan_v1.md` (initial), `bp5_refactoring_plan_v2.md` (review round 1)
**Context:** BP5 matches Tandem at p1/1000m. This plan improves code robustness, maintainability, and configurability. It introduces a constitutive model abstraction layer (Phase 3) as active refactoring work, and documents the architectural foundation for future extensions: QD-FD hybrid dynamic rupture, extended rate-and-state friction variants, nonlinear constitutive models, and multi-physics solid-fluid coupling.

---

## Changes from v2

| # | Issue | Resolution |
|---|-------|------------|
| 1 | **Phase order**: TOML before BoundaryConfig | **Reordered.** BoundaryConfig before TOML (simpler, no external dep). |
| 2 | **BCMode / BoundaryConfig redundancy** | **New constructor drops BCMode.** Old constructor delegates. BCMode removed in Phase 6. |
| 3 | **Quick-check needs external .msh** | **Use `CreateBP5InlineMesh()`.** No external mesh dependency. |
| 4 | **Golden outputs not per-rank-count** | **Separate goldens** per rank count. |
| 5 | **`RunComparison()` never fails** | **Patch it**: `--regression-tolerance` flag, nonzero exit on exceed. |
| 6 | **toml11 vendoring unspecified** | **Git submodule** under `extern/toml11`, optional `SEAS_USE_TOML`. |
| 7 | **DirichletFunc lifetime** | **Documented**: capture by value. |
| 8 | **MPI count** | **99 non-test + 21 test = 120 total.** |
| 9 | **Phase 2 byte-for-byte** | Because compiled binary is identical, not because algorithm is deterministic. |
| 10 | **Phase reversibility** | Keep BCMode dual-path until Phase 6. |
| 11 | **Mismatched config attrs** | `MFEM_VERIFY` that attrs exist in mesh. |
| 12 | **Debug doc line staleness** | `.inl` headers note original line ranges. |
| 13 | **`SEASConfig` behavior** | Pure data struct. No factory methods. |
| 14 | **Quick-check step count** | Start with 20 steps. |
| 15 | **v2 open questions** | All resolved. |
| 16 | **Phase 1 MPI too shallow** | **Tandem-inspired encapsulation**: `FaultScatter` class hides Irecv/Isend, type-safe MPI traits, communicator-as-parameter everywhere. Goal: new developers never write raw MPI. |
| 17 | **BoundaryConfig uses numeric attrs** | **Name-based BC identification**: read Gmsh physical surface names from mesh, map to BC roles. User configures by name ("fault", "dirichlet"), not by memorizing attr numbers. |
| 18 | **QD-FD hybrid underspecified** | **Detailed workflow** added: state persistence strategy, data transfer optimization, GPU potential, fault state continuity protocol to prevent wave oscillation. |
| 19 | **Slip-weakening in friction plan** | **Removed.** QD-dynamic requires rate-and-state friction. Focus on R&S variants only. |
| 20 | **Constitutive abstraction is future-only** | **Promoted to Phase 3.** Create `constitutive/` directory now with interface + `LinearElastic` wrapper. Keep elasticity as-is but set up the structure so new constitutive laws are written in one place. |
| 21 | **Damage-breakage called multi-physics** | **Corrected.** Damage-breakage is a constitutive model (Section 11.3), not multi-physics. Multi-physics = solid-fluid coupling with two/three-field formulations (Section 11.4 — noted, not urgent). |
| 22 | **DG method preference** | **SIPG as primary method.** BR2 and others kept as options but SIPG is the default going forward. |

---

## Table of Contents

0. [Phase 0: Regression Infrastructure](#0-phase-0-regression-infrastructure)
1. [Codebase Summary](#1-codebase-summary)
2. [Phase 1: MPI Safety & Encapsulation](#2-phase-1-mpi-safety--encapsulation)
3. [Phase 2: File Decomposition](#3-phase-2-file-decomposition)
4. [Phase 3: Constitutive Model Abstraction](#4-phase-3-constitutive-model-abstraction)
5. [Phase 4: BoundaryConfig](#5-phase-4-boundaryconfig)
6. [Phase 5: TOML Parameter Files](#6-phase-5-toml-parameter-files)
7. [Phase 6: Polish](#7-phase-6-polish)
8. [Project Layout After All Phases](#8-project-layout-after-all-phases)
9. [Execution & Branch Strategy](#9-execution--branch-strategy)
10. [Dependency & Coupling Analysis](#10-dependency--coupling-analysis)
11. [Resolved Questions](#11-resolved-questions)
12. [Future Architecture: Dynamic Rupture, Friction, Multi-Physics](#12-future-architecture)

---

## 0. Phase 0: Regression Infrastructure

**Must complete before any refactoring code.** No exceptions.

**Branch:** `refactor/phase0-regression`
**Risk:** None (additive only).

### 0.1 The Problem

`bp5_verification_full.cpp` `RunComparison()` (line 236) prints station data-point counts, but `main()` always returns 0. Sbatch jobs check `$?` from `ibrun` — they always "pass." There is no automated regression gate.

### 0.2 Patch `RunComparison()` (10-Line Change)

Add a `--regression-tolerance` CLI flag. When provided, `RunComparison()` computes per-station L2 error against reference data and returns nonzero if any station exceeds the tolerance. `main()` propagates the return code. Without the flag, behavior is unchanged (informational only).

### 0.3 Freeze Golden Reference Outputs

| Tier | Ranks | Steps | Mesh | Location | Purpose |
|------|-------|-------|------|----------|---------|
| Local quick-check | 8 | 20 | `CreateBP5InlineMesh(2,2,1, 200e3,100e3,100e3)` | `bp5/benchmark_data/golden_inline_8r_20step/` | Gate every commit locally |
| Frontera serial | 1 | 50 | `bp5_tandem_exact.msh` (p1/1000m) | `bp5/benchmark_data/golden_serial_1r_50step/` | Serial correctness |
| Frontera parallel | 400 | 50 | same | `bp5/benchmark_data/golden_parallel_400r_50step/` | Parallel correctness |
| Full production | 400 | ~1800 yr | same | `bp5/benchmark_data/golden_p1_1000m/` | Major milestones |

**Golden outputs must be per-rank-count** — `MPI_Allreduce(..., MPI_SUM)` on `MPI_DOUBLE` produces different FP rounding at different rank counts.

### 0.4 Regression Script

`scripts/regression_check.py`: load sim + golden, interpolate to common time points, per-station L2 error, hard pass/fail, nonzero exit on failure.

### 0.5 Quick-Check Test

Uses `CreateBP5InlineMesh()` — no external `.msh`, runs anywhere. 8 ranks, 20 steps. Validate detection power by flipping Dirichlet loading sign.

### 0.6 Verify Byte-for-Byte Reproducibility

Run twice at 8 ranks, diff. If not identical, relax Phase 2 to L2 < 1e-15.

### 0.7 Tag Baseline

```
git tag v1.0-bp5-verified -m "BP5 p1/1000m matches Tandem. Pre-refactoring baseline."
```

### 0.8 Phase 0 Checklist

- [ ] 0a. Add `--regression-tolerance` to `bp5_verification_full.cpp`
- [ ] 0b. Run Frontera: 1-rank/50-step + 400-rank/50-step goldens
- [ ] 0c. Run local: 8-rank/20-step inline-mesh golden
- [ ] 0d. Store goldens in `bp5/benchmark_data/golden_*/`
- [ ] 0e. Write `scripts/regression_check.py`
- [ ] 0f. Self-check: tolerance 0 against own golden
- [ ] 0g. Verify byte-for-byte reproducibility
- [ ] 0h. Run Frontera full production golden (background)
- [ ] 0i. Tag `v1.0-bp5-verified`
- [ ] 0j. Merge to `system_update`

---

## 1. Codebase Summary

### 1.1 Module Inventory

| Module | Directory | Files | LOC | Purpose |
|--------|-----------|-------|-----|---------|
| Config | `config/` | 3 | 747 | Benchmark parameter structs (BP1, BP2, BP5) |
| Domain | `domain/` | 6 | 10,687 | DG elasticity: antiplane (scalar) + 3D vector |
| Fault | `fault/` | 5 | 2,837 | Geometry, basis transforms, rate-and-state operator |
| Friction | `friction/` | 3 | 1,046 | Dieterich-Ruina, aging/slip state evolution |
| Integrator | `integrator/` | 4 | 2,608 | DG bilinear forms: BR2 + IP for scalar and 3D |
| Solver | `solver/` | 3 | 1,756 | SEAS QD operator, Dormand-Prince RK45 |
| I/O | `io/` | 7 | 2,977 | SCEC output, ParaView, checkpoint, probes |
| Common | `common/` | 3 | 512 | MPI context, parallel utils, type aliases |
| Trace | `trace/` | 1 | 670 | Face-level traction diagnostics |
| **Total** | | **35** | **~24,000** | |

### 1.2 Execution Flow (One Time Step)

Per `ARCHITECTURE.md`. Hot path per RK stage:

```
Mult(state, rate):
  1. fault_->GetSlip(state, slip_)
  2. domain_->ExpandOwnedToLocalFault(slip_, local_slip_)    [MPI p2p]
  3. domain_->Solve(t, local_slip_, displacement_)            [MUMPS/CG+AMG]
  4. domain_->ComputeTraction(displacement_, local_slip_, traction_)
  5. domain_->RestrictToOwnedFault(local_traction_, traction_)
  6. fault_->ComputeRHS(traction_, state, rate)               [MPI Allreduce]
```

Step 3 dominates wall time. Phase 3 (constitutive) touches how stress is computed in steps 3-4 but does not change the flow.

### 1.3 Stable Interfaces (Untouched by All Phases)

| Interface | LOC | Why stable |
|-----------|-----|------------|
| `DomainOperator<MeshType>` | 217 | Pure virtual base. `Solve()`, `ComputeTraction()` via `Vector`. |
| `SEASQuasiDynamicOperator` | 556 | Couples domain + fault through pointers. |
| `FrictionLaw` / `DieterichRuinaFriction` | 132/635 | Scalar physics. |
| `DormandPrinceRK45` | 701 | Time integrator. |

**Note:** DG integrators (4 files, 2,608 LOC) will be touched by Phase 3 (constitutive abstraction) — they gain a `ConstitutiveModel*` constructor overload. The existing `(lambda, mu)` constructors are preserved.

---

## 2. Phase 1: MPI Safety & Encapsulation

**Goal:** Harden MPI layer AND encapsulate it so future developers never write raw MPI calls.
**Branch:** `refactor/phase1-mpi`
**Regression:** Byte-for-byte match vs golden.
**Risk:** Low.

### 2.1 Tandem's MPI Design (Reference)

Tandem's `src/parallel/` (1,092 LOC across 21 files) provides a layered abstraction that we should learn from:

| Tandem Pattern | What It Does | Our Equivalent |
|----------------|-------------|----------------|
| `MPITraits<T>` | Compile-time C++ type → MPI type mapping | Not needed (we use `MPI_DOUBLE` / `MPI_INT` directly — small surface) |
| `Scatter` + `ScatterPlan` | Non-blocking element exchange with begin/test/wait | **`FaultScatter` — our primary need** |
| `CommPattern` (`AllToAllV`, `GatherV`) | Type-safe collective wrappers | **Extend `MPIContext`** |
| `LocalGhostCompositeView` | Unified local+ghost data access | Already have `ExpandOwnedToLocalFault` |
| Communicator-as-parameter | All classes accept `MPI_Comm` at construction | **Adopt everywhere** |

The key insight from Tandem: **application code never calls MPI directly**. It calls `Scatter::begin_scatter()` / `wait_scatter()`, or `MPIContext::GlobalMax()`. The raw MPI calls are isolated to ~5 infrastructure files.

### 2.2 `FaultScatter` Class (New)

The 3 raw p2p calls in `elasticity_operator.hpp:5090-5140` (Irecv/Isend/Waitall for fault face exchange) are the most dangerous MPI code in the codebase. Encapsulate them:

```cpp
// common/fault_scatter.hpp
class FaultScatter {
public:
   /// Construct from shared fault communication blocks.
   FaultScatter(const std::vector<SharedFaultCommBlock> &blocks, MPI_Comm comm);

   /// Post non-blocking receives and sends for owned->ghost exchange.
   void BeginScatter(const Vector &owned_data, int comps_per_dof);

   /// Non-blocking test — returns true when all p2p is complete.
   bool TestScatter();

   /// Block until scatter completes. Ghost data is written to output.
   void WaitScatter(Vector &local_data);

private:
   MPI_Comm comm_;
   std::vector<SharedFaultCommBlock> blocks_;
   std::vector<MPI_Request> requests_;
   std::vector<Vector> send_bufs_, recv_bufs_;
};
```

`ExpandOwnedToLocalFault()` in `elasticity_operator.hpp` becomes a thin wrapper around `FaultScatter`. All raw `MPI_Irecv`/`MPI_Isend`/`MPI_Waitall` calls move into `FaultScatter` — the operator never touches MPI directly.

### 2.3 Extend `MPIContext`

```cpp
// Additions to MPIContext:
MPIContext(int *argc, char ***argv, MPI_Comm comm = MPI_COMM_WORLD);
MPI_Comm GetComm() const { return comm_; }

// Type-safe collective wrappers (all use stored comm_):
real_t GlobalMax(real_t local) const;
real_t GlobalMin(real_t local) const;
real_t GlobalSum(real_t local) const;
int    GlobalMaxInt(int local) const;
int    GlobalSumInt(int local) const;
void   Barrier() const;
void   Bcast(real_t &value) const;
void   Bcast(Vector &vec) const;
// NEW:
void   GatherToRoot(const Vector &local, Vector &global) const;
void   AllgatherVec(const Vector &local, Vector &global) const;
```

**Goal:** Every MPI operation used by the SEAS code is available through `MPIContext` or `FaultScatter`. No `#include <mpi.h>` needed outside `common/`.

### 2.4 Other Phase 1 Items

- **`common/mpi_check.hpp`**: `MFEM_SEAS_MPI_CHECK(call)` macro for debug-mode error checking.
- **`common/mpi_tags.hpp`**: Named tag constants (`kFaultFaceExchange = 27183`).
- **Audit conditional collectives**: Document safety of each conditional MPI site as code comments.
- **Debug-mode buffer assertions** in `parallel_utils.hpp`.

### 2.5 Phase 1 Checklist

- [ ] 1a. Extend `MPIContext` with stored communicator + new collective methods
- [ ] 1b. Create `common/fault_scatter.hpp` encapsulating p2p fault exchange
- [ ] 1c. Refactor `ExpandOwnedToLocalFault()` to use `FaultScatter`
- [ ] 1d. Add `common/mpi_check.hpp` and `common/mpi_tags.hpp`
- [ ] 1e. Apply `MFEM_SEAS_MPI_CHECK` to all remaining raw MPI calls
- [ ] 1f. Audit conditional collectives, add comments
- [ ] 1g. Add buffer assertions in `parallel_utils.hpp`
- [ ] 1h. Verify: no raw `MPI_*` calls outside `common/` (except MFEM internals)
- [ ] 1i. Run quick-check regression
- [ ] 1j. Run full test suite
- [ ] 1k. Merge to `system_update`

---

## 3. Phase 2: File Decomposition

**Goal:** Split `elasticity_operator.hpp` (6,505 lines) into 6 files. Zero logic changes.
**Branch:** `refactor/phase2-decompose`
**Regression:** Byte-for-byte (compiled binary is literally identical).
**Risk:** Low.

### 3.1 Decomposition Plan

| File | LOC | Content | Original Lines |
|------|-----|---------|----------------|
| `elasticity_operator.hpp` | ~1,500 | Class declaration, `Solve()`, enums, structs, getters | Keeps existing public API |
| `elasticity_operator_setup.inl` | ~1,200 | `SetupFESpace`, `BuildFacetBCTables`, `SetupFaultInfo`, `BuildOwnedFaultLayout`, `SetupSolver`, `PrecomputeMassInverse` | Lines 1763-2961 |
| `elasticity_operator_assembly.inl` | ~1,400 | `AssembleStiffness`, `AssembleSlipContribution{IP,BR2}{,Shared}`, `AssembleDirichletLoading` | Lines 3063-4903 |
| `elasticity_operator_traction.inl` | ~1,200 | `ComputeTractionImpl`, wrappers, fault coords, `RestrictToOwnedFault`, `ExpandOwnedToLocalFault`, `BuildSlipAtQuadPoints` | Lines 4909-6440, 1558-1607 |
| `elasticity_operator_debug.inl` | ~800 | All `DebugDump*` and `Debug*` methods | Lines 767-1556 |
| `elasticity_operator_verify.inl` | ~400 | `VerifyGhostDOFCommunication`, `VerifyDirichletSkipSets`, `VerifySharedDirichletPerFace`, `VerifyRHSNorms` | Lines 284-722 |

Each `.inl` header notes original line range for debug doc traceability.

### 3.2 Phase 2 Checklist

- [ ] 2a. Confirm byte-for-byte reproducibility (Phase 0 step 0.6)
- [ ] 2b. Create 5 `.inl` files (pure cut-paste)
- [ ] 2c. Update `elasticity_operator.hpp` to include them
- [ ] 2d. Add original-line-range comments
- [ ] 2e. Update CMakeLists/Makefile header deps
- [ ] 2f. Compile: zero warnings
- [ ] 2g. Run regression: byte-for-byte
- [ ] 2h. Run full test suite
- [ ] 2i. Merge to `system_update`

---

## 4. Phase 3: Constitutive Model Abstraction

**Goal:** Create a `constitutive/` directory with an abstract interface and a `LinearElastic` implementation that wraps the current hardcoded `(lambda, mu)` stress computation. Update integrators and traction to call through the interface. Keep linear elasticity as the only implementation — but set up the structure so writing a new constitutive law means writing code in one place, and all other code (integrators, traction, assembly) calls it automatically.
**Branch:** `refactor/phase3-constitutive`
**Regression:** L2 < 1e-12 (same physics, different code path through virtual call).
**Risk:** Medium (touches integrators and traction).

### 4.1 The Problem

The stress computation sigma = lambda * tr(epsilon) * I + 2 * mu * epsilon is hardcoded in three separate locations:

1. **`DGElasticityBR2Integrator::TestNormal()`** (`dg_elasticity_br2_integrator.hpp:77-84`): Encodes the elasticity tensor action.
2. **`DGElasticityIPCombinedIntegrator`** (`dg_elasticity_ip_combined_integrator.hpp`): Evaluates stress at quadrature points.
3. **`ComputeTractionImpl()`** (`elasticity_operator_traction.inl`): Computes {sigma.n} on fault faces.

If someone writes a new constitutive law (e.g., viscoelastic, elastoplastic), they would need to find and modify all three locations — easy to miss one, leading to inconsistent physics.

### 4.2 `ConstitutiveModel` Interface

```cpp
// constitutive/constitutive_model.hpp
namespace mfem { namespace seas {

/// Abstract base for material constitutive relations.
///
/// Given strain epsilon and (optionally) internal state variables,
/// compute stress sigma and tangent stiffness C.
///
/// For linear elasticity: sigma = lambda*tr(eps)*I + 2*mu*eps, no internal vars.
/// For plasticity: sigma depends on plastic strain history.
/// For damage-breakage: sigma depends on damage tensor and breakage parameter.
class ConstitutiveModel {
public:
   virtual ~ConstitutiveModel() = default;

   /// Number of internal state variables per quadrature point.
   /// Linear elastic: 0. Plasticity: 7 (6 plastic strain + hardening).
   virtual int NumInternalVars() const = 0;

   /// Whether the model requires Newton iteration (nonlinear stress-strain).
   virtual bool IsNonlinear() const = 0;

   /// Compute stress from strain.
   /// @param epsilon Strain tensor (Voigt: [e_xx, e_yy, e_zz, e_xy, e_yz, e_xz])
   /// @param int_vars Internal state (nullptr if NumInternalVars()==0)
   /// @param sigma Output stress (Voigt: [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz])
   virtual void ComputeStress(const real_t *epsilon, const real_t *int_vars,
                              real_t *sigma) const = 0;

   /// Compute the elasticity tensor action on a normal vector:
   ///   T_{iu} = C_{iuks} * n_k   (used by DG flux)
   /// For isotropic: T_{iu} = lambda * delta_{iu} * n_k * delta_{ks} + mu * (n_i delta_{us} + n_s delta_{iu})
   /// This replaces TestNormal() in the BR2 integrator.
   virtual void TensorDotNormal(const Vector &n, int i, int u, int s,
                                real_t &val) const = 0;

   /// Compute tangent stiffness d(sigma)/d(epsilon) at a quadrature point.
   /// For linear elastic this is constant. For nonlinear it depends on state.
   virtual void ComputeTangent(const real_t *epsilon, const real_t *int_vars,
                               DenseMatrix &C_tang) const = 0;

   /// Update internal state after a converged Newton step.
   /// No-op for elastic. For plasticity: return mapping.
   virtual void UpdateState(const real_t *epsilon, real_t *int_vars) const = 0;

   /// Get Lame parameters (for penalty computation and wave speeds).
   virtual real_t GetLambda() const = 0;
   virtual real_t GetMu() const = 0;
};

}} // namespace mfem::seas
```

### 4.3 `LinearElastic` Implementation

```cpp
// constitutive/linear_elastic.hpp
class LinearElastic : public ConstitutiveModel {
public:
   LinearElastic(real_t lambda, real_t mu) : lambda_(lambda), mu_(mu) {}

   int NumInternalVars() const override { return 0; }
   bool IsNonlinear() const override { return false; }

   void ComputeStress(const real_t *eps, const real_t *, real_t *sig) const override {
      real_t tr = eps[0] + eps[1] + eps[2];
      sig[0] = lambda_ * tr + 2*mu_ * eps[0];  // sigma_xx
      sig[1] = lambda_ * tr + 2*mu_ * eps[1];  // sigma_yy
      sig[2] = lambda_ * tr + 2*mu_ * eps[2];  // sigma_zz
      sig[3] = 2*mu_ * eps[3];                  // sigma_xy
      sig[4] = 2*mu_ * eps[4];                  // sigma_yz
      sig[5] = 2*mu_ * eps[5];                  // sigma_xz
   }

   void TensorDotNormal(const Vector &n, int i, int u, int s, real_t &val) const override {
      val = lambda_ * (u == s ? 1.0 : 0.0) * n(i)
          + mu_ * ((i == u ? 1.0 : 0.0) * n(s) + (i == s ? 1.0 : 0.0) * n(u));
   }
   // ... ComputeTangent, UpdateState (no-op), GetLambda, GetMu
};
```

This is exactly the current inline code, extracted into a single authoritative location.

### 4.4 Integrator Changes

DG integrators gain a new constructor overload:

```cpp
// Existing (preserved for backward compat):
DGElasticityBR2Integrator(Coefficient &lambda, Coefficient &mu, ...);

// New (preferred):
DGElasticityBR2Integrator(const ConstitutiveModel &model, ...);
```

The old constructor internally creates a `LinearElastic` from the coefficients. All existing code compiles unchanged. The `TestNormal()` private method delegates to `model_.TensorDotNormal()`.

**DG method preference:** SIPG (Interior Penalty) is the primary method going forward. BR2 integrators are preserved and continue to work, but new development (constitutive coupling, dynamic rupture) targets SIPG first. The constructor default changes from `"BR2"` to `"IP"` in new drivers.

### 4.5 Traction Changes

`ComputeTractionImpl()` currently constructs a `DGElasticityIPCombinedIntegrator` with `lambda_coeff_` and `mu_coeff_`. After Phase 3, it passes `constitutive_model_*` instead. For linear elastic, the computation is identical.

### 4.6 What This Enables

A developer writing a new constitutive law (e.g., viscoelastic):
1. Creates one file: `constitutive/viscoelastic.hpp`
2. Implements `ConstitutiveModel` interface (5 virtual methods)
3. Passes it to `ElasticityDomainOperator` constructor
4. **Done.** Assembly, traction, and DG flux all use the new law automatically.

No hunting through three separate files to update stress computations.

### 4.7 Phase 3 Checklist

- [ ] 3a. Create `constitutive/constitutive_model.hpp` with abstract interface
- [ ] 3b. Create `constitutive/linear_elastic.hpp` wrapping current `(lambda, mu)` logic
- [ ] 3c. Add `ConstitutiveModel*` constructor overload to DG integrators (IP and BR2)
- [ ] 3d. Old `(lambda, mu)` constructors delegate to new via internal `LinearElastic`
- [ ] 3e. Update `ElasticityDomainOperator` to accept `ConstitutiveModel*`
- [ ] 3f. Old constructor creates `LinearElastic` internally (backward compat)
- [ ] 3g. Update `ComputeTractionImpl()` to use `ConstitutiveModel`
- [ ] 3h. Verify all tests compile and pass with old constructors (zero-change path)
- [ ] 3i. Write unit test constructing operator with explicit `LinearElastic`
- [ ] 3j. Run regression: L2 < 1e-12
- [ ] 3k. Merge to `system_update`

---

## 5. Phase 4: BoundaryConfig

**Goal:** Parameterize boundary conditions using physical surface names from the mesh, not numeric attribute IDs. Users should never memorize which number means which BC.
**Branch:** `refactor/phase4-boundary-config`
**Regression:** L2 < 1e-10.
**Risk:** Medium.

### 5.1 Name-Based BC Identification

**Problem:** Currently `BCMode::FarField` means "attr 5 = Dirichlet, attr 1 = Natural." But Gmsh physical surface numbers are arbitrary — a user's mesh might use attr 7 for the fault and attr 2 for Dirichlet. Worse, inline meshes and Gmsh meshes use different attribute conventions.

**Solution:** Read physical surface names from the Gmsh `.msh` file and map names to BC roles. The user configures by name:

```toml
[boundary]
fault = "fault_surface"          # Physical Surface name in .msh
dirichlet = "far_field"          # Physical Surface name in .msh
natural = ["top_surface", "bottom_surface"]  # Can be multiple
```

For inline meshes (no Gmsh names), fall back to numeric attributes.

### 5.2 `BoundaryConfig` (Value Type)

```cpp
// domain/boundary_config.hpp
namespace mfem { namespace seas {

struct BoundaryConfig {
   // Name-based identification (preferred for Gmsh meshes):
   std::string fault_name;                    // e.g., "fault_surface"
   std::vector<std::string> dirichlet_names;  // e.g., {"far_field"}
   std::vector<std::string> natural_names;    // e.g., {"top", "bottom"}

   // Numeric fallback (for inline meshes or when names not available):
   std::set<int> dirichlet_attrs;
   std::set<int> natural_attrs;
   int fault_attr = -1;  // -1 = use name-based

   /// Resolve names to attribute numbers using mesh physical group names.
   /// Called once at construction. After this, dirichlet_attrs/natural_attrs/fault_attr
   /// are populated and the operator uses numeric attrs internally.
   void ResolveFromMesh(const Mesh &mesh);

   static BoundaryConfig BP5Default() {
      BoundaryConfig bc;
      bc.fault_name = "fault";
      bc.dirichlet_names = {"dirichlet"};
      bc.natural_names = {"natural"};
      // Numeric fallback for Tandem .geo convention:
      bc.dirichlet_attrs = {5};
      bc.natural_attrs = {1};
      bc.fault_attr = 3;
      return bc;
   }
};

}} // namespace mfem::seas
```

`ResolveFromMesh()` reads the Gmsh physical group name-to-attribute mapping (MFEM stores this from Gmsh import), looks up each name, and populates the numeric attr sets. If names are not found (inline mesh), falls back to the numeric attrs.

### 5.3 Multiple Dirichlet Functions

Different boundaries can have different Dirichlet functions:

```cpp
using DirichletFunc = std::function<void(const Vector &x, real_t t, Vector &u_D)>;

struct BoundaryConfig {
   // ...
   /// Per-attribute Dirichlet functions. Key = boundary attr number.
   /// If an attr is in dirichlet_attrs but not in this map, it gets the default.
   std::map<int, DirichletFunc> dirichlet_funcs;

   /// Default Dirichlet function (applied to all Dirichlet attrs not in the map).
   DirichletFunc default_dirichlet_func;
};
```

This allows different loading on different boundaries (e.g., plate motion on far-field X faces, zero on far-field Y faces).

**Lifetime rule:** All `DirichletFunc` lambdas capture parameters by value.

### 5.4 Constructor Changes

Same delegating pattern as v2: new constructor takes `BoundaryConfig`, old constructor delegates with `BoundaryConfig::BP5Default()`. BCMode dual-path kept until Phase 6.

### 5.5 `MFEM_VERIFY` for Mismatched Config

```cpp
// In SetupBoundaryMarkers():
for (int attr : bdr_config_.dirichlet_attrs) {
   MFEM_VERIFY(mesh_has_bdr_attr(attr),
      "BoundaryConfig: Dirichlet attr " << attr << " not found in mesh."
      " Available attrs: " << list_mesh_attrs());
}
```

### 5.6 Phase 4 Checklist

- [ ] 4a. Create `domain/boundary_config.hpp` with name-based `BoundaryConfig`
- [ ] 4b. Implement `ResolveFromMesh()` using MFEM's Gmsh physical group data
- [ ] 4c. Add per-attribute `DirichletFunc` support
- [ ] 4d. Add new constructor overload to `ElasticityDomainOperator`
- [ ] 4e. Implement deprecated old constructor as delegating constructor
- [ ] 4f. Dual-path `SetupBoundaryMarkers()` (BoundaryConfig + old BCMode)
- [ ] 4g. Add `MFEM_VERIFY` for attr existence
- [ ] 4h. Refactor `BuildFacetBCTables()` and `AssembleDirichletLoading()`
- [ ] 4i. Verify all tests pass with old constructor
- [ ] 4j. Write test with name-based BoundaryConfig
- [ ] 4k. Run regression: L2 < 1e-10
- [ ] 4l. Merge to `system_update`

---

## 6. Phase 5: TOML Parameter Files

**Goal:** Add TOML config parsing. Old CLI driver preserved.
**Branch:** `refactor/phase5-toml`
**Regression:** L2 < 1e-12.
**Risk:** Medium (external dependency).

### 6.1 toml11 Integration

Git submodule under `extern/toml11/`. Optional `SEAS_USE_TOML` build flag.

### 6.2 `SEASConfig` Struct

Pure data struct (no factory methods). Mirrors TOML file structure:

```cpp
struct SEASConfig {
   std::string benchmark;   // "bp5", "custom", etc.
   MeshConfig mesh;
   DomainConfig domain;     // order, dg_method ("IP" default), penalty_factor
   SolverConfig solver;
   MaterialConfig material;
   FrictionConfig friction;
   LoadingConfig loading;
   FaultGeometryConfig fault_geometry;
   NucleationConfig nucleation;
   TimeSteppingConfig time_stepping;
   BoundaryConfig boundary;  // Uses name-based from Phase 4
   OutputConfig output;
};
```

**DG method default:** `DomainConfig::dg_method = "IP"` (SIPG). BR2 remains available as `"BR2"`.

### 6.3 Generic Driver

`drivers/seas_driver.cpp`: Parses TOML, applies CLI overrides, validates, constructs operators, runs simulation.

### 6.4 Phase 5 Checklist

- [ ] 5a. Add `extern/toml11` submodule, `SEAS_USE_TOML` flag
- [ ] 5b. Define `SEASConfig` in `config/seas_config.hpp`
- [ ] 5c. Implement parser and validation
- [ ] 5d. Create example TOML files
- [ ] 5e. Write `drivers/seas_driver.cpp`
- [ ] 5f. Run both drivers, compare: L2 < 1e-12
- [ ] 5g. Merge to `system_update`

---

## 7. Phase 6: Polish

**Goal:** Clean up, standardize, document.
**Branch:** `refactor/phase6-polish`
**Regression:** Byte-for-byte.
**Risk:** Low.

- Remove `BCMode` enum and dual-path
- Standardize error handling (`MFEM_ASSERT/VERIFY/ABORT`)
- Add logging wrapper with verbosity levels
- Update `CODEBASE_GUIDE.md`, `ARCHITECTURE.md`
- Add `README.md` with quick-start
- Add future-dynamic comment in `domain_operator.hpp`
- Update sbatch scripts

---

## 8. Project Layout After All Phases

```
miniapps/seas/
├── common/
│   ├── mpi_context.hpp              P1: stored communicator + collective methods
│   ├── fault_scatter.hpp            NEW P1: encapsulated p2p fault exchange
│   ├── mpi_check.hpp               NEW P1: debug error macro
│   ├── mpi_tags.hpp                NEW P1: named tag constants
│   ├── parallel_utils.hpp           P1: debug assertions
│   └── seas_types.hpp
├── constitutive/                    NEW P3
│   ├── constitutive_model.hpp       Abstract interface
│   └── linear_elastic.hpp           Wraps current (lambda, mu)
├── config/
│   ├── bp1_params.hpp
│   ├── bp2_params.hpp
│   ├── bp5_params.hpp
│   ├── seas_config.hpp             NEW P5
│   └── seas_config_parser.hpp      NEW P5 (requires SEAS_USE_TOML)
├── domain/
│   ├── domain_operator.hpp          P6: future-dynamic comment
│   ├── antiplane_operator.hpp
│   ├── antiplane_bdrload_operator.hpp
│   ├── elasticity_operator.hpp      P2: ~1500 LOC
│   ├── elasticity_operator_setup.inl     NEW P2
│   ├── elasticity_operator_assembly.inl  NEW P2
│   ├── elasticity_operator_traction.inl  NEW P2
│   ├── elasticity_operator_debug.inl     NEW P2
│   ├── elasticity_operator_verify.inl    NEW P2
│   ├── boundary_config.hpp          NEW P4
│   ├── bp2_mesh.hpp
│   └── seas_boundary_tags.hpp
├── fault/                           (unchanged)
├── friction/                        (unchanged)
├── integrator/                      P3: gain ConstitutiveModel* overloads
├── solver/                          (unchanged)
├── io/                              (unchanged)
├── trace/                           (unchanged)
├── extern/toml11/                   NEW P5
├── drivers/seas_driver.cpp          NEW P5
├── params/*.toml                    NEW P5
├── scripts/regression_check.py      NEW P0
└── bp5/benchmark_data/golden_*/     NEW P0
```

---

## 9. Execution & Branch Strategy

### 9.1 Phase Order

```
Phase 0 (regression) → Phase 1 (MPI) → Phase 2 (decompose) → Phase 3 (constitutive) → Phase 4 (BoundaryConfig) → Phase 5 (TOML) → Phase 6 (polish)
```

### 9.2 Summary Table

| Phase | Branch | Risk | Regression | Key Deliverable |
|-------|--------|------|-----------|-----------------|
| 0 | `refactor/phase0-regression` | None | Self-check | Golden outputs + regression script |
| 1 | `refactor/phase1-mpi` | Low | Byte-for-byte | `FaultScatter`, extended `MPIContext`, no raw MPI outside `common/` |
| 2 | `refactor/phase2-decompose` | Low | Byte-for-byte | `elasticity_operator.hpp` → 6 files |
| 3 | `refactor/phase3-constitutive` | Medium | L2 < 1e-12 | `constitutive/` directory, `LinearElastic`, integrator overloads |
| 4 | `refactor/phase4-boundary` | Medium | L2 < 1e-10 | Name-based `BoundaryConfig`, per-attr `DirichletFunc` |
| 5 | `refactor/phase5-toml` | Medium | L2 < 1e-12 | TOML config + generic driver |
| 6 | `refactor/phase6-polish` | Low | Byte-for-byte | BCMode removed, logging, docs |

### 9.3 Rollback Cost

| Phase | Rollback |
|-------|----------|
| 0 | Zero (additive) |
| 1 | Remove `FaultScatter`, inline MPI back. Moderate but straightforward. |
| 2 | `cat *.inl >> .hpp`. 5 minutes. |
| 3 | Remove `constitutive/`, remove integrator overloads. Old `(lambda,mu)` constructors still work. |
| 4 | Remove `boundary_config.hpp`. BCMode dual-path still works. |
| 5 | Remove TOML parser + driver + submodule. Old driver untouched. |
| 6 | Restore BCMode. Remove logging. |

---

## 10. Dependency & Coupling Analysis

### 10.1 DAG After All Phases

```
common/ ← constitutive/ ← config/ ← friction/ ← integrator/ ← domain/ ← fault/ ← solver/ ← io/ ← drivers/
                                                      ↑                    ↑
                                            constitutive_model.hpp   boundary_config.hpp
                                            (abstract, no deps)      (value-type, no deps)
```

`constitutive/` depends only on `mfem.hpp`. Integrators depend on `constitutive/constitutive_model.hpp` (abstract interface only). Framework code never depends on config-parsing or drivers.

---

## 11. Resolved Questions

| Question | Resolution |
|----------|------------|
| Explicit template instantiation? | Defer. Only if compile > 2 min per TU. |
| Quick-check step count? | 20. Validate with sign-flip test. |
| Test file splitting? | Defer. |
| `FaceVertexKey` for hex? | Defer. Code comment only. |
| Config validation depth? | Minimum + mesh-attr check. |
| BCMode deprecation? | Phase 4: dual-path. Phase 6: removed. |
| toml11 vendoring? | Git submodule, optional `SEAS_USE_TOML`. |
| Phase order? | 0-1-2-3(constitutive)-4(BoundaryConfig)-5(TOML)-6(polish). |
| `SEASConfig` behavior? | Pure data struct only. |
| `DirichletFunc` lifetime? | Capture by value. Documented. |
| DG method preference? | SIPG primary. BR2 preserved as option. |
| Constitutive in which phase? | Phase 3. Interface + LinearElastic now; nonlinear later. |

---

## 12. Future Architecture: Dynamic Rupture, Friction, Multi-Physics

This section documents planned extensions. **No code for these in Phases 0-6.** The goal is to ensure current decisions don't close doors.

### 12.1 QD-FD Hybrid Dynamic Rupture — Detailed Workflow

**Overview:** A fully-dynamic (FD) elastic wave solver in velocity-stress variables (Dumbser & Kaeser 2006, SeiSol). 9 unknowns per node: Q = (sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz, u_dot, v_dot, w_dot). Explicit time stepping with Riemann solver at interfaces. Coupled to the existing QD system through a slip-rate threshold.

**QD vs FD comparison:**

| Aspect | QD (current) | FD (future) |
|--------|-------------|-------------|
| Unknowns | Displacement u (3/node) | Stress + velocity Q (9/node) |
| Time stepping | Implicit (global sparse solve) | Explicit (element-local, CFL-limited) |
| Interface flux | DG SIPG/BR2 penalty | Riemann solver |
| Radiation damping | eta*V in friction | Not needed (waves propagate) |
| Domain BCs | Dirichlet (far-field) | Absorbing (outgoing waves) |

#### 12.1.1 Regime Switching Protocol

```
QD regime (normal):
  - Monitor V_max via MPI_Allreduce(MPI_MAX) every RK stage
  - When V_max > V_activate (e.g., 1e-3 m/s):
      1. FREEZE QD state: save (displacement, slip, psi) to in-memory snapshot
      2. Convert QD displacement -> FD velocity-stress:
           sigma_ij = ConstitutiveModel->ComputeStress(grad(u))   [element-local]
           v_i = 0 (or finite-difference estimate from recent displacements)
      3. Initialize FD absorbing boundaries
      4. Switch regime = FD

FD regime (earthquake):
  - Time step with CFL constraint: dt_FD ~ h / c_p (much smaller than QD dt)
  - Riemann solver at all element interfaces
  - Friction constraint at fault interface (same FrictionLaw, but eta=0)
  - Monitor V_max
  - When V_max < V_deactivate (e.g., 1e-6 m/s) AND seismic energy has exited:
      1. RESTORE QD snapshot from memory
      2. Update slip: slip_new = slip_frozen + integral(V_FD dt) from FD
      3. Update state variable: psi_new = psi from last FD step
      4. Recompute QD traction from updated slip (one Solve())
      5. Switch regime = QD, resume adaptive RK45
```

#### 12.1.2 State Persistence Strategy

**QD state during FD: keep in memory, do NOT write to disk.**

| Data | Size (p1/1000m, 400 ranks) | Strategy |
|------|---------------------------|----------|
| Displacement `u` | ~3M DOFs * 8 bytes = ~24 MB | In-memory snapshot |
| Slip (owned fault DOFs) | ~50k DOFs * 2 comp * 8 bytes = ~0.8 MB | In-memory |
| State variable psi | ~50k DOFs * 8 bytes = ~0.4 MB | In-memory |
| Stiffness matrix (QD) | ~1 GB sparse | **Keep allocated.** Do not free. Reassembly is expensive. |

Total snapshot: ~25 MB + stiffness matrix. Negligible compared to FD memory requirements (9 DOFs/node vs 3, plus flux matrices). The stiffness matrix stays allocated during FD — freeing and reassembling it when returning to QD is far more expensive than the memory cost.

**Checkpointing:** Write a checkpoint at the QD→FD transition point. If the FD run crashes, restart from this checkpoint rather than re-running the entire QD interseismic period.

#### 12.1.3 Data Transfer Optimization

QD→FD conversion (displacement → velocity-stress) is **element-local**: compute grad(u) per element, then sigma = C : epsilon. No MPI communication needed. This is embarrassingly parallel.

**GPU potential:** The element-local stress conversion is a perfect GPU kernel:
- Each element: read nodal displacements, compute gradient, apply constitutive law, write Q
- No inter-element communication
- Memory-bound, highly parallel
- MFEM already has GPU-accelerated element operations

The FD explicit time stepping is also GPU-friendly (element-local mass inverse, face-local Riemann solve). This is where GPU acceleration has the highest impact — not in QD (dominated by sparse direct solve).

#### 12.1.4 Fault State Continuity — Preventing Wave Pollution

**The critical challenge:** When switching QD→FD, any inconsistency between the stress field and the displacement field creates artificial waves that propagate through the domain. These are numerical artifacts, not physics.

**Prevention protocol:**

1. **Stress-displacement consistency:** Compute sigma from grad(u) using the SAME constitutive model. The Phase 3 `ConstitutiveModel` abstraction is essential here — both QD and FD use the same `ConstitutiveModel::ComputeStress()`, guaranteeing consistency.

2. **Velocity initialization:** Set initial velocity to zero (or a smooth estimate). Do NOT use finite differences from QD time history — the QD adaptive dt varies wildly and FD differences would be noisy.

3. **Ramp-up period:** After QD→FD switch, run FD for a few wave-crossing times (domain_size / c_p) before trusting the FD solution. During ramp-up, absorbing boundaries damp any artificial waves.

4. **FD→QD transition:** The critical invariant is that slip and psi are continuous. The displacement field is reconstructed from the updated slip via one QD Solve(). Any residual stress mismatch is absorbed by the next QD time step's implicit solve.

5. **Fault state is the coupling variable:** Both regimes agree on fault slip and psi. The domain fields (displacement for QD, velocity-stress for FD) are derived quantities. As long as the fault state is continuous, the domain fields will be consistent (up to the radiation damping correction, which is absent in FD).

#### 12.1.5 Proposed Files

```
dynamic/
├── wave_operator.hpp              # DG wave equation operator (velocity-stress)
├── riemann_solver.hpp             # Exact Riemann solver for isotropic elastic waves
├── absorbing_bc.hpp               # Absorbing BCs (outgoing-wave flux)
├── free_surface_bc.hpp            # Free-surface BC via inverse Riemann
├── fault_riemann.hpp              # Fault Riemann solver with friction constraint
├── ader_time_integrator.hpp       # ADER time integration (or explicit RK with CFL)
└── velocity_stress_state.hpp      # Q = (sigma, v) container + conversion utilities

solver/
└── seas_hybrid_operator.hpp       # QD <-> FD regime switching
```

---

### 12.2 Extended Friction Laws (Rate-and-State Variants)

**Scope:** QD-dynamic simulations require rate-and-state friction. Slip-weakening is not compatible with the QD regime (no steady-state slip rate). We focus exclusively on R&S variants.

**Planned variants:**

| Variant | How It Differs from Current | State Variables |
|---------|---------------------------|-----------------|
| Dieterich-Ruina aging (current) | Baseline | psi |
| Dieterich-Ruina slip law (current) | Different state evolution | psi |
| Regularized R&S + flash heating | `f(V,T) = f_RS(V) * [1 - (1-f_w/f_RS) * exp(-V_w*Dc_th / (V*Dc))]` | psi, T |
| R&S + thermal pressurization | `sigma_n_eff = sigma_n - p(T)`, heat diffusion on fault | psi, T, p |

**Interface evolution:** The current `FrictionLaw` base (132 LOC) works for pure R&S. For flash heating and thermal pressurization, the interface needs to accept additional state variables per DOF. Rather than breaking the existing interface, we add a `FaultInterfaceLaw` as a generalized base:

```cpp
class FaultInterfaceLaw {
   virtual int NumStateVars() const = 0;       // 1 for basic R&S, 2 for +heating, 3 for +pressurization
   virtual real_t SolveSlipRate(real_t tau, real_t sigma_n, real_t eta,
                                const real_t *state_vars, const real_t *params) const = 0;
   virtual void StateRates(real_t V, const real_t *state_vars,
                           const real_t *params, real_t *rates) const = 0;
};
```

Existing `FrictionLaw` + `StateEvolution` become a convenience adapter for the basic R&S case. `RateStateFaultOperator` is updated to talk through `FaultInterfaceLaw`.

---

### 12.3 Future Nonlinear Constitutive Models

Phase 3 establishes the `ConstitutiveModel` interface with `LinearElastic`. Future implementations:

| Model | File | Internal Vars | Key Change |
|-------|------|--------------|------------|
| Linear elastic (Phase 3) | `linear_elastic.hpp` | 0 | None — baseline |
| Viscoelastic (Maxwell/Kelvin) | `viscoelastic.hpp` | 6 (viscous strain) | Time-dependent stress relaxation |
| Elastoplastic (J2, Drucker-Prager) | `elastoplastic.hpp` | 7 (plastic strain + hardening) | `Solve()` becomes Newton iteration |
| Damage-breakage | `damage_breakage.hpp` | 2+ (damage D, breakage B) | Stiffness degrades with damage; breakage controls granular transition |

**Damage-breakage is a constitutive model, not multi-physics.** It modifies the stress-strain relation C(D,B) : epsilon. The damage and breakage parameters evolve via ODEs coupled to the strain field. This fits naturally into `ConstitutiveModel::UpdateState()` — no additional PDE solver needed.

**Key architectural point:** When `IsNonlinear()` returns true, `ElasticityDomainOperator::Solve()` must switch from a single linear solve to a Newton iteration loop:

```
while not converged:
   Assemble K_tangent using ConstitutiveModel::ComputeTangent(epsilon_k, int_vars)
   Solve K_tangent * delta_u = residual
   u_{k+1} = u_k + delta_u
   Update epsilon, check convergence
ConstitutiveModel::UpdateState(epsilon_converged, int_vars)
```

This Newton loop is a future addition to `Solve()`, activated only when the constitutive model reports `IsNonlinear()`. The Phase 3 linear path is unaffected.

---

### 12.4 Multi-Physics: Solid-Fluid Coupling (Note)

**This is not urgent. We note the structure for future reference.**

Multi-physics here means coupling solid mechanics with porous flow (Biot's equations) or thermal diffusion. This is fundamentally different from constitutive model changes (Section 12.3) because it introduces **additional PDEs with their own FE spaces and solvers**, plus **coupling terms** that appear in both equations.

**Two-field formulation (u, p):** Displacement + pore pressure. The solid equation gains a `-alpha * grad(p)` coupling term. The fluid equation gains a `alpha * div(du/dt)` coupling term. These are not simply "added terms" — they require a coupled system or a staggered iteration scheme.

**Three-field formulation (u, p, T):** Adds temperature. Shear heating Q = tau*V enters the thermal equation. Thermal expansion enters the solid equation. Pore pressure depends on temperature through thermal pressurization.

**Code structure goal:** Each physics module lives in its own file with its own governing equation. A coupling operator assembles the coupled system or manages staggered iteration:

```
physics/
├── physics_operator.hpp           # Abstract base: Solve(t, coupling_fields)
├── elasticity_physics.hpp         # Wraps ElasticityDomainOperator
├── porous_flow_physics.hpp        # Biot: div(q) + alpha*div(du/dt) = 0
├── thermal_physics.hpp            # Heat eq: rho*cp*dT/dt = div(k*grad(T)) + Q
└── coupled_solver.hpp             # Staggered or monolithic coupling
```

For now, we simply ensure that:
1. The `ConstitutiveModel` interface (Phase 3) can accept auxiliary fields (p, T) as future arguments
2. The `DomainOperator` base class doesn't prevent adding coupling terms later
3. Each physics module is a separate file — no monolithic multi-physics class

---

### 12.5 Design Rules for Future Extensions

1. **One file per model.** Each constitutive model, friction variant, and physics module gets its own `.hpp`. No 6000-line files.

2. **Abstract base → concrete implementations.** New models implement the interface, not modify existing code.

3. **No cross-directory includes between peer modules.** `constitutive/` never includes `friction/`. `dynamic/` never includes `domain/`. Coupling is through interfaces passed by the driver.

4. **Factory pattern for runtime selection.** Models selected by string name in config, constructed by factory.

5. **State variable layout owned by the model.** Each law declares `NumStateVars()` / `NumInternalVars()`. Operators allocate accordingly.

6. **SIPG is the primary DG method.** New development (constitutive integration, dynamic Riemann solver, multi-physics) targets SIPG first. BR2 is preserved for comparison/legacy but is not the default path for new features.

7. **The driver owns the coupling loop.** Regime switching (QD/FD), staggered iteration, convergence checks live in the driver or hybrid operator, not inside individual modules.

### 12.6 What Phases 0-6 Must NOT Do

| Rule | Rationale |
|------|-----------|
| Don't hardcode `dim=3` or `num_slip_comp=2` in interfaces | Future 2D dynamic reuses same abstractions |
| Don't embed `(lambda, mu)` deeper — use `ConstitutiveModel` | Phase 3 addresses this directly |
| Don't merge `RateStateFaultOperator` with friction logic | Fault operator orchestrates; friction computes. Separate. |
| Don't make `SEASQuasiDynamicOperator` aware of specific friction types | Talks through `RateStateFaultOperator` only |
| Don't close `DomainOperator` base to extension | FD wave operator needs sibling base or extended interface |
| Don't store `Vp_` as authoritative Dirichlet source | Phase 4's `DirichletFunc` callback is the right abstraction |
| Keep `FaultBasis`, `FaultGeometry` solver-independent | Both QD and FD use same fault coordinate transforms |
| Don't free stiffness matrix during FD regime | Reassembly cost >> memory cost |

---

## Appendix: Known Limitations

Current status after this plan:

1. **Matrix assembly** — stored, not matrix-free. Memory-bound at high order. *Not addressed.*
2. **Direct solver default** — MUMPS for BP5. Iterative less robust for DG. *Not addressed.*
3. **Single fault plane at Y=0** — multi-fault / non-planar not supported. *Phase 4 BoundaryConfig is a partial step (configurable fault group). Full multi-fault requires future work.*
4. **Quasi-dynamic only** — no fully-dynamic wave propagation. *Addressed architecturally by Section 12.1. Implementation is future work.*
5. **DG noise on unstructured meshes** — BR2 traction noise at sharp element-size transitions. *Not addressed. SIPG as new default may help.*
6. **Linear elasticity only** — no nonlinear constitutive models. *Phase 3 establishes the abstraction. Nonlinear implementations are future work (Section 12.3).*
7. **Single-physics** — no porous flow or thermal coupling. *Noted architecturally in Section 12.4. Not urgent.*
8. **Rate-and-state friction only** — no extended R&S variants (flash heating, thermal pressurization). *Addressed architecturally by Section 12.2. Implementation is future work.*
