# SEAS-MFEM Robustness & Refactoring Plan (v4)

**Date:** 2026-04-09
**Status:** Final — ready for execution
**Predecessors:** v1 (initial), v2 (review round 1), v3 (review round 2 + future architecture)
**Context:** BP5 matches Tandem at p1/1000m. This plan improves robustness, maintainability, and configurability across 7 phases. It introduces a constitutive model abstraction (Phase 3) and documents the architectural foundation for QD-FD hybrid dynamic rupture, extended R&S friction, nonlinear constitutive models, and multi-physics coupling.

---

## Changes from v3

| # | Issue | Resolution |
|---|-------|------------|
| 1 | **B1: MFEM has no Gmsh physical group name API** | **Dropped name-based BC from Phase 4.** Use numeric attrs only. Add companion mapping file (TOML: `{"fault": 3, "dirichlet": 5}`) as future option when MFEM adds the API or when we parse Gmsh `.msh` directly. |
| 2 | **B2: IP-vs-BR2 default change untested** | **Gate default change on explicit comparison.** Phase 5 checklist item: run BP5 p1/1000m with both methods, compare station L2. If IP is worse, keep method-specific defaults. |
| 3 | **D1: Constructor overload stacking** | **Design final constructor signature upfront.** Phase 3 introduces the full signature (ConstitutiveModel + BoundaryConfig, both required). Only 2 overloads ever exist. |
| 4 | **D2: ConstitutiveModel ownership** | **`const &` for new constructor, non-owning `const*` internal.** Old constructor owns `unique_ptr<LinearElastic>` and passes `*owned_model_`. Lifetime rule documented. |
| 5 | **D3: TensorDotNormal couples constitutive to DG** | **Removed from interface.** DG integrators derive it from `ComputeTangent()` by contracting C_{ijks} with n_k. Constitutive model authors only implement stress and tangent — no DG knowledge needed. |
| 6 | **D4: GetLambda/GetMu are isotropic-specific** | **Replaced with `GetPenaltyModulus()` and `GetMaxWaveSpeed()`.** Both well-defined for any elastic-type material (max eigenvalue of C, sqrt(max_eigenvalue/rho)). |
| 7 | **D6: Phase 3 tolerance 1e-12 too loose** | **Target byte-for-byte.** Relax to L2 < 1e-15 only if virtual dispatch causes compiler FP reordering. |
| 8 | **S2: Frontera serial-vs-parallel pair every phase** | **Added to every phase checklist** as explicit item. |
| 9 | **S3: Quick-check detection with 2x Vp** | **Added to Phase 0.** Sign flip may not manifest in 20 interseismic steps; 2x Vp changes loading rate from step 1. |
| 10 | **S5: IP-vs-BR2 comparison before default change** | **Added as Phase 5 gating checklist item.** |
| 11 | **S7: Nonlinear assembly cost note** | **Added to Section 12.3.** Newton iteration changes assembly from ~5% to ~30-50% of wall time. |
| 12 | **S12: Record actual L2 per phase** | **Added to every phase checklist.** Set future tolerances relative to measured baseline. |
| 13 | **T1-T15: 15 new unit tests** | **All added to phase checklists.** ~3,200 LOC of new tests across all phases. |
| 14 | **All v3 review checklist additions** | **Incorporated verbatim** into phase checklists. |
| 15 | **Phase 0 golden: compare all 7 SCEC columns** | **Specified in regression_check.py requirements.** |

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
13. [Phase 7: New Driver Full Simulation & Verification](#13-phase-7-new-driver-full-simulation--verification)

---

## 0. Phase 0: Regression Infrastructure

**Must complete before any refactoring code.** No exceptions.
**Branch:** `refactor/phase0-regression`

### 0.1 Patch `RunComparison()`

Add `--regression-tolerance` CLI flag. `RunComparison()` returns nonzero on exceed. `main()` propagates exit code.

### 0.2 Freeze Golden Reference Outputs

| Tier | Ranks | Steps | Mesh | Location | Status |
|------|-------|-------|------|----------|--------|
| Local quick-check | 8 | 20 | `CreateBP5InlineMesh(...)` | `golden_inline_8r_20step/` | **TODO** |
| Frontera serial | 1 | 50 | bp5_tandem_exact.msh | `golden_serial_1r_50step/` | **DONE** (job 7643780) |
| Frontera parallel | 400 | 50 | same | `golden_parallel_400r_50step/` | **DONE** (job 7643831) |
| Full production | 400 | ~1800 yr | same | `golden_tandem_p1_1000m/` | **DONE** (job 7639398, v62 fix, ~1593 yr) |

Golden outputs are per-rank-count (MPI_Allreduce FP rounding differs by rank count).

**Available golden data (as of 2026-04-09):**

All stored in `bp5/benchmark_data/`:

- **`golden_serial_1r_50step/`**: 10 station files (13 data points each) + 1 global file (52 data points) + 4 trace CSV files. Prefix: `bp5_verify_serial_job7643780`. Frontera job 7643780, 1 rank, 50 accepted steps, t_final = 1.24 s.
- **`golden_parallel_400r_50step/`**: 10 station files (12 data points each) + 1 global file + 48 trace CSV files (12 ranks × 4 trace types). Prefix: `bp5_verify_parallel_job7643831`. Frontera job 7643831, 400 ranks, 50 accepted steps.
- **`golden_tandem_p1_1000m/`**: 10 station files (~86,000 data points each). Prefix: `bp5_v62_fix_job7639398`. Frontera job 7639398, 400 ranks, full production run, t_final ≈ 1593 years. This is the verified-correct v62 output that matches Tandem. No global file in this directory.

All station files use the standard SCEC 8-column format: `time(s), slip_strike(m), slip_dip(m), log10(V_strike), log10(V_dip), tau_strike(MPa), tau_dip(MPa), log10(state)(s)`.

**Remaining:** Local 8-rank/20-step inline-mesh golden (item 0c) is not yet produced — requires running on local machine after Phase 0 code changes.

### 0.3 `scripts/regression_check.py`

Requirements:
- Compares **all 7 non-time SCEC columns** per station (slip_dip, slip_strike, log10_V_dip, log10_V_strike, tau_dip, tau_strike, log10_state) **plus** global output file (V_max, event count)
- Interpolates to common time points
- Hard pass/fail with nonzero exit code
- Reports per-station, per-field L2 errors

### 0.4 Quick-Check Test

Uses `CreateBP5InlineMesh()`. 8 ranks, 20 steps. Validate detection power with **2x Vp perturbation** (not just sign flip — sign flip may not manifest in 20 interseismic steps; 2x Vp changes loading rate from step 1).

### 0.5 Comment Cleanup: Deferred to Phase 6

The codebase contains ~432 references to "Tandem" across ~44 source files. These should be rewritten to cite published literature instead. However, this is a ~200-comment diff touching ~30 files — a large merge-conflict surface that conflicts with Phase 0's "additive only, no risk" principle. **Moved to Phase 6 (polish)** where large non-functional changes belong. See Section 7 for the full cleanup specification.

### 0.6 Phase 0 Checklist

- [ ] 0a. Add `--regression-tolerance` to `bp5_verification_full.cpp`
- [x] 0b. Run Frontera: 1-rank/50-step + 400-rank/50-step goldens *(done: jobs 7643780, 7643831)*
- [ ] 0c. Extract `CreateBP5InlineMesh()` from `bp5_verification_full.cpp` (line 96) into a shared header (e.g., `config/bp5_mesh_utils.hpp`) so both the verification driver and the quick-check test can call it
- [ ] 0c-ii. Run local: 8-rank/20-step inline-mesh golden
- [x] 0d. Store goldens in `bp5/benchmark_data/golden_*/` *(done: serial, parallel, full production)*
- [ ] 0e. Write `scripts/regression_check.py` (all 7 SCEC columns + global output)
- [ ] 0f. Self-check: tolerance 0 against own golden → exit 0
- [ ] 0f-ii. Self-test `regression_check.py`: feed golden against deliberately perturbed copy (multiply slip by 1.001) → nonzero exit; feed golden against itself → zero exit
- [ ] 0f-iii. Self-test `--regression-tolerance`: run bp5_verification_full at 8r/20step with `--regression-tolerance 0` + golden as reference → exit 0; with `--regression-tolerance 1e-30` → nonzero exit
- [ ] 0g. Verify byte-for-byte reproducibility (run twice, diff)
- [ ] 0g-ii. `make test` passes with zero code changes (baseline confirmation)
- [ ] 0h. Validate quick-check detection: run with `Vp *= 2`, confirm test catches it
- [x] 0i. Run Frontera full production golden *(done: job 7639398, v62 fix, ~1593 yr, stored in golden_tandem_p1_1000m/)*
- [ ] 0j. Tag `v1.0-bp5-verified`
- [ ] 0k. Merge to `system_update`

---

## 1. Codebase Summary

### 1.1 Module Inventory

| Module | Directory | Files | LOC |
|--------|-----------|-------|-----|
| Config | `config/` | 3 | 747 |
| Domain | `domain/` | 6 | 10,687 |
| Fault | `fault/` | 5 | 2,837 |
| Friction | `friction/` | 3 | 1,046 |
| Integrator | `integrator/` | 4 | 2,608 |
| Solver | `solver/` | 3 | 1,756 |
| I/O | `io/` | 7 | 2,977 |
| Common | `common/` | 3 | 512 |
| Trace | `trace/` | 1 | 670 |
| **Total** | | **35** | **~24,000** |

### 1.2 Stable Interfaces (Untouched by All Phases)

| Interface | Why stable |
|-----------|------------|
| `DomainOperator<MeshType>` | Pure virtual base. `Solve()`, `ComputeTraction()` via `Vector`. |
| `SEASQuasiDynamicOperator` | Couples domain + fault through pointers. |
| `FrictionLaw` / `DieterichRuinaFriction` | Scalar physics. |
| `DormandPrinceRK45` | Time integrator. |

**DG integrators** (4 files, 2,608 LOC) gain `ConstitutiveModel*` overloads in Phase 3. Existing constructors preserved.

---

## 2. Phase 1: MPI Safety & Encapsulation

**Goal:** Encapsulate raw MPI so future developers never write MPI calls directly. Inspired by Tandem's `src/parallel/` (Scatter + ScatterPlan + CommPattern pattern).
**Branch:** `refactor/phase1-mpi`
**Regression:** Byte-for-byte.

### 2.1 `FaultScatter` Class

Encapsulates the 3 raw p2p calls in `elasticity_operator.hpp:5090-5140`:

```cpp
// common/fault_scatter.hpp
class FaultScatter {
public:
   FaultScatter(const std::vector<SharedFaultCommBlock> &blocks, MPI_Comm comm);
   void BeginScatter(const Vector &owned_data, int comps_per_dof);
   bool TestScatter();       // Non-blocking test
   void WaitScatter(Vector &local_data);
private:
   MPI_Comm comm_;
   std::vector<MPI_Request> requests_;
   std::vector<Vector> send_bufs_, recv_bufs_;
};
```

`ExpandOwnedToLocalFault()` becomes a thin wrapper. No raw `MPI_Irecv`/`MPI_Isend` outside `common/`.

### 2.2 Extended `MPIContext`

Stored communicator + new methods: `GatherToRoot(Vector)`, `AllgatherVec(Vector)`. All existing methods switch to stored `comm_`.

### 2.3 Other Items

- `common/mpi_check.hpp`: `MFEM_SEAS_MPI_CHECK(call)` macro
- `common/mpi_tags.hpp`: Named tag constants
- Audit conditional collectives; document safety as comments
- Debug-mode buffer assertions in `parallel_utils.hpp`

### 2.4 Phase 1 Checklist

- [ ] 1a. Extend `MPIContext` with stored communicator + GatherToRoot, AllgatherVec
- [ ] 1a-ii. Extend `test_mpi_context.cpp`: test GatherToRoot, AllgatherVec at 2/4 ranks; test stored communicator (construct with dup'd comm)
- [ ] 1b. Create `common/fault_scatter.hpp`
- [ ] 1b-ii. Write `test_fault_scatter.cpp` (~300 LOC): known owned data → BeginScatter → WaitScatter → verify ghost data; test at 2/4 ranks; empty blocks; asymmetric block sizes
- [ ] 1c. Refactor `ExpandOwnedToLocalFault()` to use `FaultScatter`
- [ ] 1d. Add `common/mpi_check.hpp` and `common/mpi_tags.hpp`
- [ ] 1d-ii. Test `MFEM_SEAS_MPI_CHECK`: debug build triggers assertion on bad call; release build compiles as no-op
- [ ] 1e. Audit conditional collectives, add comments
- [ ] 1f. Buffer assertions in `parallel_utils.hpp`
- [ ] 1g. Verify: no raw `MPI_Irecv`/`MPI_Isend`/`MPI_Waitall` (p2p) calls outside `common/`. Remaining collectives (`MPI_Allreduce`, `MPI_Comm_rank`, etc.) are migrated to `MPIContext` methods opportunistically — not required for Phase 1 merge.
- [ ] 1h. Run quick-check regression: byte-for-byte
- [ ] 1h-ii. Run specifically: `test_parallel_elasticity`, `test_serial_parallel_consistency`, `test_bp5_parallel_smoke` at 4 ranks
- [ ] 1h-iii. Run Frontera serial-vs-parallel pair, diff `[VERIFY]` lines
- [ ] 1i. Record actual L2 in debug document
- [ ] 1j. Run full test suite
- [ ] 1k. Merge to `system_update`

---

## 3. Phase 2: File Decomposition

**Goal:** Split `elasticity_operator.hpp` (6,505 lines) into 6 files. Zero logic changes.
**Branch:** `refactor/phase2-decompose`
**Regression:** Byte-for-byte (identical binary).

### 3.1 Decomposition Plan

| File | LOC | Original Lines |
|------|-----|----------------|
| `elasticity_operator.hpp` | ~1,500 | Core declarations, `Solve()`, enums, structs |
| `elasticity_operator_setup.inl` | ~1,200 | Lines 1763-2961 |
| `elasticity_operator_assembly.inl` | ~1,400 | Lines 3063-4903 |
| `elasticity_operator_traction.inl` | ~1,200 | Lines 4909-6440, 1558-1607 |
| `elasticity_operator_debug.inl` | ~800 | Lines 767-1556 |
| `elasticity_operator_verify.inl` | ~400 | Lines 284-722 |

### 3.2 Phase 2 Checklist

- [ ] 2a. Confirm byte-for-byte reproducibility (Phase 0 step 0g)
- [ ] 2b. Create 5 `.inl` files (pure cut-paste)
- [ ] 2c. Update `elasticity_operator.hpp` to include them
- [ ] 2d. Add original-line-range comments to each `.inl` header
- [ ] 2e. Update CMakeLists/Makefile header deps
- [ ] 2e-ii. Verify incremental rebuild: edit one line in `elasticity_operator_setup.inl`, run make without clean, confirm test binary relinks with the change
- [ ] 2f. Compile: zero warnings
- [ ] 2g. Run quick-check regression: byte-for-byte
- [ ] 2g-ii. Run `test_elasticity_operator`, `test_cross_verify_tandem`, `test_elasticity_br2` explicitly
- [ ] 2h. Run full test suite
- [ ] 2i. Merge to `system_update`

---

## 4. Phase 3: Constitutive Model Abstraction

**Goal:** Create `constitutive/` directory with abstract interface + `LinearElastic` wrapper. Update integrators and traction to call through the interface. Writing a new constitutive law = one file in one directory.
**Branch:** `refactor/phase3-constitutive`
**Regression:** Target byte-for-byte. Relax to L2 < 1e-15 only if virtual dispatch causes FP reordering.

### 4.1 `ConstitutiveModel` Interface

The interface must handle three tiers of constitutive complexity:

| Tier | Example | State Variables | State Evolution | Non-Local |
|------|---------|----------------|-----------------|-----------|
| 1. Algebraic | Linear elastic | None | None | No |
| 2. Evolving | Damage-breakage (Lyakhovsky & Ben-Zion 2014) | α (damage), B (breakage), ε^p (plastic strain) | dα/dt, dB/dt, dε^p/dt are ODEs integrated alongside fault state | Optional (∇α in stress) |
| 3. Non-local | Damage with gradient regularization | Same as Tier 2 | Same as Tier 2 | Yes — ∇α appears in stress (Eq. 14 of L&BZ 2014) |

**Key insight from the damage-breakage model:** The state variables α and B have their own time-evolution ODEs (Eqs. 25-26 in L&BZ 2014):

```
dα/dt = (1-B) [C_d I₂(ξ - ξ₀) + D∇²α]           (damage evolution)
dB/dt = C_B P(α)(1-B) I₂ × [...]                   (breakage evolution)
dε^(p)_ij/dt = C_g · B^m₁ · τ^(m₂-1) · τᵢⱼ       (plastic strain rate)
```

These are NOT post-convergence updates — they must be **time-integrated by the RK45 stepper** alongside the fault slip and psi. The elastic moduli λ(α), μ(α), γ(α) change with damage (Eq. 12), so the stiffness matrix must be reassembled when α evolves. For non-local models, ∇α appears directly in the stress tensor (Eq. 14), requiring α to have its own FE space.

**Phase 3 implements Tier 1 only** (LinearElastic). But the interface is designed so Tiers 2-3 work without breaking changes — all evolution/non-local methods have default no-op implementations.

```cpp
// constitutive/constitutive_model.hpp
class ConstitutiveModel {
public:
   virtual ~ConstitutiveModel() = default;

   // =========================================================================
   // Tier 1: Core stress-strain interface (required for ALL models)
   // =========================================================================

   /// Internal state variables per quadrature point.
   /// Linear elastic: 0. Damage-breakage: 8 (α, B, 6 plastic strain components).
   virtual int NumInternalVars() const = 0;

   /// Whether Newton iteration is required (nonlinear σ-ε relation).
   virtual bool IsNonlinear() const = 0;

   /// Compute stress from strain (Voigt: [εxx, εyy, εzz, εxy, εyz, εxz]).
   /// For damage-breakage: moduli depend on int_vars[0] = α.
   virtual void ComputeStress(const real_t *epsilon, const real_t *int_vars,
                              real_t *sigma) const = 0;

   /// Compute 6×6 tangent stiffness d(σ)/d(ε) in ENGINEERING Voigt notation.
   /// Index map: 0=xx, 1=yy, 2=zz, 3=xy, 4=yz, 5=xz.
   /// Shear strains are ENGINEERING (γ_xy = 2ε_xy).
   /// For isotropic: C(0,0) = λ+2μ, C(3,3) = μ (NOT 2μ).
   /// DG integrators derive TensorDotNormal from this.
   virtual void ComputeTangent(const real_t *epsilon, const real_t *int_vars,
                               DenseMatrix &C_tang) const = 0;

   /// Update internal state after a converged Newton step (Tier 1 path).
   /// For algebraic plasticity: return mapping. For evolving models: no-op
   /// (state is advanced by ComputeStateRates instead).
   virtual void UpdateState(const real_t *epsilon, real_t *int_vars) const = 0;

   /// Max eigenvalue of C. Used for DG penalty and CFL.
   virtual real_t GetPenaltyModulus() const = 0;

   /// Max wave speed sqrt(max_eig(C)/ρ). Used for CFL and radiation damping.
   virtual real_t GetMaxWaveSpeed(real_t rho) const = 0;

   // =========================================================================
   // Tier 2: Time-evolving state variables (default: no evolution)
   // =========================================================================

   /// Whether state variables evolve via ODEs that must be time-integrated.
   /// When true, the time integrator calls ComputeStateRates() at each RK stage
   /// and advances int_vars alongside the fault slip/psi state.
   /// Linear elastic: false. Damage-breakage: true. Viscoelastic: true.
   virtual bool HasStateEvolution() const { return false; }

   /// Compute d(int_vars)/dt for time integration.
   /// @param epsilon   Current strain at this quadrature point
   /// @param int_vars  Current state [NumInternalVars()]
   /// @param rates     Output: time derivatives [NumInternalVars()]
   virtual void ComputeStateRates(const real_t *epsilon,
                                   const real_t *int_vars,
                                   real_t *rates) const {}

   // =========================================================================
   // Tier 3: Non-local models (default: local only)
   // =========================================================================

   /// Number of state variables that require FE-space gradient computation.
   /// 0 for most models. 1 for non-local damage (∇α in stress, Eq. 14 L&BZ).
   /// When > 0, the domain operator allocates an FE space for these variables
   /// and computes their gradients before calling ComputeStressNonLocal.
   virtual int NumNonLocalVars() const { return 0; }

   /// Compute stress including non-local gradient terms.
   /// Default: delegates to ComputeStress (ignoring gradient).
   /// Damage-breakage override adds: -ϑ ∇ᵢα · ∇ⱼα (structural stresses).
   virtual void ComputeStressNonLocal(const real_t *epsilon,
                                       const real_t *int_vars,
                                       const real_t *grad_state,
                                       real_t *sigma) const {
      ComputeStress(epsilon, int_vars, sigma);
   }
};
```

**What the time integrator does when `HasStateEvolution() == true`:**

```
Mult(state, rate):                        // SEASQuasiDynamicOperator
  1. Extract fault slip from state
  2. Extract constitutive state from state    ← NEW (α, B, ε^p per quad point)
  3. If NonLocal: compute ∇α via FE space     ← NEW
  4. Solve domain (Newton if nonlinear)
  5. Compute traction
  6. Compute fault rates (dslip/dt, dpsi/dt)
  7. Compute constitutive rates               ← NEW: model->ComputeStateRates()
  8. Pack all rates into output vector
```

Steps 2, 3, 7 are no-ops for `LinearElastic` — zero overhead for current code.

**Key design decision (D3):** `TensorDotNormal()` is NOT on the interface. DG integrators derive it from `ComputeTangent()` by contracting `C_{ijks} * n_k`. Constitutive model authors only implement stress and tangent — no DG knowledge required.

**Key design decision (D4):** `GetLambda()`/`GetMu()` replaced by `GetPenaltyModulus()` and `GetMaxWaveSpeed(rho)`. Both well-defined for anisotropic, damaged, or degraded tensors.

**Design principle:** A developer writing damage-breakage creates ONE file (`constitutive/damage_breakage.hpp`) implementing all tiers. The framework (domain operator, time integrator) queries `HasStateEvolution()` and `NumNonLocalVars()` to activate the appropriate code paths. No hunting through multiple files.

### 4.2 `LinearElastic` Implementation

```cpp
// constitutive/linear_elastic.hpp
class LinearElastic : public ConstitutiveModel {
public:
   LinearElastic(real_t lambda, real_t mu) : lambda_(lambda), mu_(mu) {}
   int NumInternalVars() const override { return 0; }
   bool IsNonlinear() const override { return false; }
   void ComputeStress(...) const override { /* sigma = lambda*tr(eps)*I + 2*mu*eps */ }
   void ComputeTangent(...) const override { /* 6x6 isotropic Voigt stiffness */ }
   void UpdateState(...) const override { /* no-op */ }
   real_t GetPenaltyModulus() const override { return lambda_ + 2*mu_; }
   real_t GetMaxWaveSpeed(real_t rho) const override { return std::sqrt((lambda_+2*mu_)/rho); }
   // Convenience (isotropic-only, not on base interface):
   real_t GetLambda() const { return lambda_; }
   real_t GetMu() const { return mu_; }
};
```

### 4.3 Integrator Changes

DG integrators gain a `ConstitutiveModel*` constructor. The integrator derives `TensorDotNormal` from `ComputeTangent()`:

```cpp
// Inside DG integrator — private helper:
// Computed ONCE per face (not per DOF), cached as 27-real table.
void CacheTangentNormalContraction(const ConstitutiveModel &model,
                                    const Vector &n) {
   model.ComputeTangent(nullptr, nullptr, C_cached_);  // 6x6 Voigt
   // Pre-contract: T_{ius} = C_{Voigt(iu), Voigt(ks)} * n_k
   for (int i = 0; i < 3; i++)
     for (int u = 0; u < 3; u++)
       for (int s = 0; s < 3; s++)
         T_cached_(i,u,s) = contract(C_cached_, n, i, u, s);
}
// Then per-DOF call is a single table lookup — O(1):
real_t TestNormal(int i, int u, int s) const { return T_cached_(i,u,s); }
```

**Performance:** Current `TestNormal()` is 7 FLOPs inlined. The generic contraction from `ComputeTangent()` is ~18 FLOPs. But with the per-face cache (27 reals, computed once), the per-DOF call becomes a single array lookup — faster than the current 7-FLOP inline. For `LinearElastic` with homogeneous material, the cache can be computed once globally. Assembly is ~5% of QD wall time, so even without caching the impact is negligible (~0.1%), but for future explicit dynamic where assembly dominates, the cache matters.

**DG method preference:** SIPG is the primary method. BR2 preserved as option. Default in new drivers is `"IP"`, but the change is gated on IP-vs-BR2 comparison (Phase 5, item 5e-ii).

### 4.4 Constructor Design (D1: Final Signature Upfront)

**Problem with raw parameter list:** The current constructor has 11+ parameters including `face_basis_type`, `penalty_factor`, `blr_tol`, `check_residual`, `match_quad_order` — all set via post-construction setters because they don't fit in the constructor. The plan's earlier signature also used `const std::string &method` which silently breaks callers using the `DGMethod` enum.

**Solution: `DomainConfig` struct.** Bundle all discretization/solver parameters into a struct. This absorbs `method`, `solver_type`, `face_basis_type`, `penalty_factor`, `blr_tol`, `check_residual`, and `match_quad_order` — eliminating all post-construction setters. The struct reuses the existing `DGMethod` enum (no string typos). It also aligns with `SEASConfig::DomainConfig` in Phase 5's TOML schema.

```cpp
// domain/domain_config.hpp (created in Phase 3)
enum class DGMethod { IP, BR2 };  // existing enum, moved here

struct DomainConfig {
   DGMethod method = DGMethod::BR2;
   SolverType solver_type = SolverType::CG_AMG;
   int face_basis_type = BasisType::GaussLobatto;
   real_t penalty_factor = 1.0;
   real_t blr_tol = 1e-10;
   bool check_residual = false;
   bool match_quad_order = false;
};
```

**Phase 3 introduces the final constructor signature:**

```cpp
// The ONE new constructor:
ElasticityDomainOperator(MeshType &mesh, int order,
                         const ConstitutiveModel &model,
                         real_t Wf, real_t lf,
                         const BoundaryConfig &bdr_config,   // required — user provides
                         DirichletFunc dirichlet_func,        // required — user provides
                         const DomainConfig &domain_config = {});
```

**Key changes from the current constructor:**
- `Vp` removed: captured inside `DirichletFunc` lambda (the operator never uses Vp directly after Phase 4). During Phase 3 transition, the old constructor passes `Vp` into `MakeBP5DirichletFunc(Vp)`.
- `method` and `solver_type` absorbed into `DomainConfig` (uses enum, not string — no typo risk).
- `face_basis_type`, `penalty_factor`, `blr_tol`, `check_residual`, `match_quad_order` absorbed into `DomainConfig` — no more post-construction setters needed.
- `Wf` and `lf` (fault depth/length) stay as positional params — needed by `SetupFaultInfo()` for geometric face classification.

**Old constructor (deprecated, delegates) — the ONLY place legacy values appear:**

```cpp
ElasticityDomainOperator(MeshType &mesh, int order,
                         real_t lambda, real_t mu, real_t Vp, real_t Wf, real_t lf,
                         DGMethod method = DGMethod::BR2,
                         SolverType solver_type = SolverType::CG_AMG,
                         BCMode bc_mode = BCMode::FarField)
{
   owned_model_ = std::make_unique<LinearElastic>(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3; bc.dirichlet_attrs = {5}; bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);
   DomainConfig dc;
   dc.method = method; dc.solver_type = solver_type;
   // delegate to new constructor...
}
```

**Phase 6 cleanup:** After BCMode removal, also remove `Vp` from the old constructor if `DirichletFunc` has fully replaced it. Verify no internal code reads `Vp_` anymore.

**Ownership (D2):** New constructor takes `const ConstitutiveModel &` (non-owning). Stores `const ConstitutiveModel *model_` internally. Old constructor owns `std::unique_ptr<LinearElastic> owned_model_` and passes `*owned_model_`. Lifetime rule: the model must outlive the operator. Documented in header comment.

### 4.5 Phase 3 Checklist

- [ ] 3a. Create `constitutive/constitutive_model.hpp`
- [ ] 3b. Create `constitutive/linear_elastic.hpp`
- [ ] 3b-ii. Write `test_linear_elastic.cpp` (~300 LOC): ComputeStress (pure tension, shear, hydrostatic, zero strain); ComputeTangent (symmetric, positive definite, correct diagonal); GetPenaltyModulus == lambda+2*mu; GetMaxWaveSpeed; NumInternalVars==0; IsNonlinear==false
- [ ] 3c. Add `ConstitutiveModel*` constructor to DG integrators; derive TensorDotNormal from ComputeTangent
- [ ] 3c-ii. Write `test_constitutive_integrator.cpp` (~400 LOC): **KEY equivalence test** — construct integrators with old `(lambda,mu)` and new `LinearElastic(lambda,mu)`; assemble element matrices on small mesh; compare bit-for-bit. Same for full operator: compare displacement and traction vectors.
- [ ] 3d. Old `(lambda,mu)` integrator constructors delegate via internal `LinearElastic`
- [ ] 3e. Create `domain/boundary_config.hpp` (struct definition + `DirichletFunc` typedef only — no MFEM_VERIFY, no per-attr logic; those come in Phase 4)
- [ ] 3e-ii. Create `domain/domain_config.hpp` with `DomainConfig` struct (bundles method, solver_type, face_basis_type, penalty_factor, blr_tol, check_residual, match_quad_order)
- [ ] 3f. Add new constructor to `ElasticityDomainOperator` (final signature from 4.4; `BoundaryConfig` required, `DomainConfig` has defaults)
- [ ] 3f-ii. Old constructor delegates with `owned_model_` + internally-built `BoundaryConfig` + `DomainConfig` (legacy values for existing test meshes only)
- [ ] 3f-iii. Remove post-construction setters (`SetCheckResidual`, `SetBLRTol`, `SetMatchQuadOrder`, `SetPenaltyFactor`, `SetFaceBasisType`) — all absorbed into `DomainConfig`
- [ ] 3g. Update `ComputeTractionImpl()` to use `ConstitutiveModel`
- [ ] 3h. Verify all tests compile and pass with OLD constructor
- [ ] 3h-ii. Run specifically: `test_cross_verify_tandem`, `test_serial_parallel_consistency`
- [ ] 3i. Write unit test constructing operator with explicit `LinearElastic`
- [ ] 3j. Run quick-check regression: target byte-for-byte; relax to L2 < 1e-15 only if needed
- [ ] 3j-ii. Record actual L2 error in debug document
- [ ] 3j-iii. Run Frontera serial-vs-parallel pair
- [ ] 3k. Run full test suite
- [ ] 3l. Merge to `system_update`

---

## 5. Phase 4: BoundaryConfig

**Goal:** Parameterize boundary attributes and Dirichlet functions. Eliminate hardcoded attr numbers.
**Branch:** `refactor/phase4-boundary-config`
**Regression:** L2 < 1e-10 (record actual L2; set future tolerance to 10x measured).

### 5.1 Design

**No hardcoded numeric presets.** The `BoundaryConfig` struct has no `BP5Default()` or any other preset with hardcoded attribute numbers. Attribute numbers are mesh-specific — a user's mesh may use any numbering convention. The user must always specify boundary attributes explicitly, either through the TOML config file (Phase 5) or as constructor arguments.

**The old deprecated constructor** (from Phase 3) is the only place that internally constructs a `BoundaryConfig` with specific numbers. This is purely for backward compatibility with existing tests and will be removed in Phase 6.

**Numeric attrs for now.** MFEM does not expose Gmsh physical group names. We use numeric attrs, which the user reads from their Gmsh `.geo` file. Future: when MFEM adds a name API (or we parse `.msh` directly), extend `BoundaryConfig` with string-based lookup.

### 5.2 `BoundaryConfig` (Value Type)

```cpp
// domain/boundary_config.hpp
struct BoundaryConfig {
   std::set<int> dirichlet_attrs;   ///< User-specified from TOML/constructor
   std::set<int> natural_attrs;     ///< User-specified from TOML/constructor
   int fault_attr;                  ///< User-specified from TOML/constructor

   /// Per-attribute Dirichlet functions. Key = attr number.
   /// Attrs not in map use default_dirichlet_func.
   std::map<int, DirichletFunc> dirichlet_funcs;
   DirichletFunc default_dirichlet_func;

   /// NO static presets. User always provides values explicitly.
   /// The TOML config file is the intended input mechanism:
   ///
   ///   [boundary]
   ///   dirichlet = [5]
   ///   natural = [1]
   ///   fault = 3
};
```

**DirichletFunc lifetime (D8):** All lambdas capture by value. Documented in `boundary_config.hpp`:

```cpp
/// DirichletFunc captures parameters BY VALUE at construction.
/// Changing Vp later requires creating a new lambda.
using DirichletFunc = std::function<void(const Vector &x, real_t t, Vector &u_D)>;
```

### 5.3 Constructor Changes

Phase 3 introduced the final signature. `BoundaryConfig` is **required** from the start — no default argument:

```cpp
// New constructor (Phase 3+4 combined):
ElasticityDomainOperator(MeshType &mesh, int order,
                         const ConstitutiveModel &model,
                         real_t Vp, real_t Wf, real_t lf,
                         const BoundaryConfig &bdr_config,   // required, no default
                         DirichletFunc dirichlet_func,        // required, no default
                         const std::string &method = "BR2",
                         SolverType solver_type = SolverType::CG_AMG);
```

**Old constructor (deprecated, backward compat only):** Constructs `BoundaryConfig` internally with the numbers that match the existing BP5 test meshes. This is the ONLY place these numbers appear — not in `BoundaryConfig` itself:

```cpp
ElasticityDomainOperator(MeshType &mesh, int order,
                         real_t lambda, real_t mu, real_t Vp, real_t Wf, real_t lf,
                         const std::string &method, SolverType solver_type,
                         BCMode bc_mode)
{
   owned_model_ = std::make_unique<LinearElastic>(lambda, mu);
   // These numbers are legacy — they match existing BP5 test meshes only.
   // New code must pass BoundaryConfig explicitly from user config.
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);
   // delegate to new constructor...
}
```

### 5.4 Internal Changes

- `SetupBoundaryMarkers()`: Dual-path (BoundaryConfig + old BCMode) for reversibility. `MFEM_VERIFY` that all attrs exist in mesh.
- `BuildFacetBCTables()`: Use `bdr_config_.fault_attr`
- `AssembleDirichletLoading()`: Look up `dirichlet_funcs[attr]` or `default_dirichlet_func`

### 5.5 TOML Integration (Phase 5)

When Phase 5 adds TOML, the user specifies boundary attrs in the config file:

```toml
[boundary]
dirichlet = [5]       # User checks their .geo file for these numbers
natural = [1]
fault = 3

[boundary.dirichlet_functions]
# Different loading on different boundaries (optional):
# 5 = "plate_motion"   # future: named function registry
```

The driver reads `[boundary]` from TOML → constructs `BoundaryConfig` → passes to operator. No hardcoded numbers anywhere in the framework.

### 5.6 Phase 4 Checklist

- [ ] 4a. Create `domain/boundary_config.hpp` (no static presets)
- [ ] 4b. Implement per-attribute DirichletFunc support
- [ ] 4b-ii. Write `test_boundary_config.cpp` (~500 LOC): construct with explicit attrs; invalid attr → MFEM_VERIFY; per-attr DirichletFunc dispatch; default func fallback
- [ ] 4c. `BoundaryConfig` is required in new constructor (no default value)
- [ ] 4d. Old deprecated constructor builds BoundaryConfig internally from legacy numbers
- [ ] 4e. Dual-path `SetupBoundaryMarkers()` (BoundaryConfig + old BCMode)
- [ ] 4f. `MFEM_VERIFY` for attr existence in mesh
- [ ] 4g. Refactor `BuildFacetBCTables()` and `AssembleDirichletLoading()`
- [ ] 4h. Verify all tests pass with old constructor
- [ ] 4h-ii. Write `test_boundary_config_operator.cpp` (~400 LOC): **KEY equivalence test** — old constructor vs new constructor with explicit `{5}, {1}, 3`; compare displacement and traction
- [ ] 4h-iii. Write `test_boundary_classification.cpp` (~300 LOC): `BuildFacetBCTables` face-by-face verification
- [ ] 4i. Run quick-check regression: L2 < 1e-10
- [ ] 4i-ii. Record actual L2 in debug document. Set Phase 5+ tolerance to `min(stated_phase_tolerance, 10 × measured)`. E.g., if measured L2 = 1e-14, Phase 5 tolerance = min(1e-12, 1e-13) = 1e-13.
- [ ] 4i-iii. Run Frontera serial-vs-parallel pair
- [ ] 4j. Run full test suite
- [ ] 4k. Merge to `system_update`

---

## 6. Phase 5: TOML Parameter Files

**Goal:** Add TOML config parsing. Old CLI driver preserved.
**Branch:** `refactor/phase5-toml`
**Regression:** L2 < 1e-12.

### 6.1 toml11 Integration

Git submodule `extern/toml11/`. Optional `SEAS_USE_TOML`. Framework compiles without it.

### 6.2 `SEASConfig` Struct

Pure data struct (D9). No factory methods. `DomainConfig::dg_method` default: `"IP"` — **gated on IP-vs-BR2 comparison (item 5e-ii)**.

Includes simulation mode selection:

```cpp
struct SimulationConfig {
   std::string mode = "qd";  // "qd", "dynamic", "hybrid"
   // No hybrid/dynamic subsections yet — those are added when the
   // dynamic/ directory is built. The TOML parser accepts unknown
   // sections gracefully (log warning, don't error).
};
```

The driver uses `SimulationConfig::mode` to select which operators to construct (see Section 12.1.1). For Phase 5, only `"qd"` is functional. `"dynamic"` and `"hybrid"` produce a clear error: `"mode 'dynamic' requires the dynamic solver module (not yet implemented)"`. The `[simulation.hybrid]` and `[simulation.dynamic]` TOML subsections (V_activate, cfl_factor, etc.) are **not defined until the dynamic code exists** — designing config for unimplemented features is premature. The parser ignores unknown subsections under `[simulation]` (log warning, don't error). `test_seas_config_parser.cpp` (T13) tests only `mode = "qd"` parsing; hybrid/dynamic subsection parsing is added when that code arrives.

### 6.3 Phase 5 Checklist

- [ ] 5a. Add `extern/toml11` submodule, `SEAS_USE_TOML` flag
- [ ] 5b. Define `SEASConfig` in `config/seas_config.hpp`
- [ ] 5c. Implement parser and validation
- [ ] 5c-ii. Write `test_seas_config_parser.cpp` (~600 LOC): parse valid TOML (verify every field); benchmark presets (`benchmark="bp5"` fills defaults from `BP5Params`); field overrides (preset + explicit `friction.b = 0.05`); CLI overrides (`--material.density 3000`); validation (negative V0 → error, missing mesh → error, empty Dirichlet set → error, zero tolerance → error); roundtrip (parse → write → parse → compare); bit-level comparison of parsed doubles vs `BP5Params` defaults
- [ ] 5d. Create example TOML files
- [ ] 5e. Write `drivers/seas_driver.cpp`
- [ ] 5e-ii. **GATE: IP-vs-BR2 comparison.** Run BP5 p1/1000m with `--dg-method BR2` (current verified) and `--dg-method IP`. Compare station L2, peak V_max, recurrence interval. If IP L2 > 1e-3, keep method-specific defaults (IP for antiplane, BR2 for 3D elasticity). Document result. Only change default if IP matches.
- [ ] 5e-iii. If default changes, add startup warning: `"[WARNING] DG method default changed from BR2 to IP in v2.0. Use --dg-method BR2 for previous behavior."` Prints once if no explicit `--dg-method`. Remove in Phase 6.
- [ ] 5f. Run both drivers (old CLI + new TOML), compare: L2 < 1e-12
- [ ] 5f-ii. Run Frontera serial-vs-parallel pair with new driver
- [ ] 5g. Run full test suite
- [ ] 5h. Merge to `system_update`

---

## 7. Phase 6: Polish

**Goal:** Clean up, standardize, document.
**Branch:** `refactor/phase6-polish`
**Regression:** Byte-for-byte.

### 7.1 Comment Cleanup: Remove Tandem-Referencing Language

Rewrite ~200 comments across ~30 files to cite published literature instead of Tandem source code. Zero code changes — comments only. Done as a single dedicated commit for clean `git blame`.

| Current Pattern | Replace With |
|----------------|-------------|
| `"following Tandem's lift_skeleton"` | `"BR2 lifting operator (Bassi & Rebay 1997)"` |
| `"Tandem's DieterichRuinaAgeing::slip_rate"` | `"Regularized R&S friction solver (Rice et al. 2001)"` |
| `"Tandem convention"` | `"SEAS convention"` or describe the choice mathematically |
| `"Tandem AdapterBase.cpp:72-73"` | `"Gram-Schmidt with ref-normal alignment"` |
| `"port of Tandem's zeroIn"` | `"Brent's method (Forsythe, Malcolm & Moler 1977)"` |
| `"matches Tandem's assembleSurface"` | `"DG interior face assembly with SIPG flux"` |

**Keep Tandem references for:** BP5 benchmark comparisons in `bp5_params.hpp`, `test_cross_verify_tandem.cpp` (rename to `test_cross_code_verification.cpp`), external reference data directories.

### 7.2 Phase 6 Checklist

- [ ] 6a. Grep all tests for `BCMode` usage
- [ ] 6a-ii. Update each test to use `BoundaryConfig` equivalent; verify all pass
- [ ] 6a-iii. Only then remove `BCMode` enum and dual-path in `SetupBoundaryMarkers()`
- [ ] 6a-iv. After removing BCMode, also remove `Vp` from the new constructor if `DirichletFunc` has fully replaced it. Verify no internal code reads `Vp_` anymore.
- [ ] 6b. Standardize error handling (MFEM_ASSERT/VERIFY/ABORT)
- [ ] 6c. Add logging wrapper with verbosity levels
- [ ] 6c-ii. Write `test_logging.cpp` (~100 LOC)
- [ ] 6d. Replace scattered `if (rank==0) cout` with `Log()`
- [ ] 6e. Remove DG method default change warning (if added in Phase 5)
- [ ] 6f. **Comment cleanup** (Section 7.1): rewrite ~200 Tandem-referencing comments. Single dedicated commit.
- [ ] 6g. Update `CODEBASE_GUIDE.md`, `ARCHITECTURE.md`
- [ ] 6h. Add `README.md` with quick-start
- [ ] 6i. Add future-dynamic comment in `domain_operator.hpp`
- [ ] 6j. Update sbatch scripts
- [ ] 6k. Run regression: byte-for-byte
- [ ] 6k-ii. Run Frontera serial-vs-parallel pair — final confirmation
- [ ] 6l. Run full test suite
- [ ] 6m. Merge to `system_update`

---

## 8. Project Layout After All Phases

```
miniapps/seas/
├── common/
│   ├── mpi_context.hpp              P1: stored communicator + collectives
│   ├── fault_scatter.hpp            NEW P1
│   ├── mpi_check.hpp               NEW P1
│   ├── mpi_tags.hpp                NEW P1
│   ├── parallel_utils.hpp           P1: debug assertions
│   └── seas_types.hpp
├── constitutive/                    NEW P3
│   ├── constitutive_model.hpp       Abstract interface
│   └── linear_elastic.hpp           Current (lambda, mu) wrapped
├── config/
│   ├── bp1_params.hpp
│   ├── bp2_params.hpp
│   ├── bp5_params.hpp
│   ├── seas_config.hpp             NEW P5
│   └── seas_config_parser.hpp      NEW P5
├── domain/
│   ├── domain_operator.hpp          P6: future-dynamic comment
│   ├── elasticity_operator.hpp      P2: ~1500 LOC
│   ├── elasticity_operator_*.inl    NEW P2 (5 files)
│   ├── boundary_config.hpp          NEW P3 (struct+typedef), extended P4 (MFEM_VERIFY, per-attr)
│   ├── domain_config.hpp            NEW P3 (DGMethod, SolverType, face_basis, penalty, etc.)
│   ├── antiplane_operator.hpp
│   ├── antiplane_bdrload_operator.hpp
│   ├── bp2_mesh.hpp
│   └── seas_boundary_tags.hpp
├── fault/                           (unchanged)
├── friction/                        (unchanged)
├── integrator/                      P3: ConstitutiveModel overloads
├── solver/                          (unchanged)
├── io/                              (unchanged)
├── trace/                           (unchanged)
├── extern/toml11/                   NEW P5
├── drivers/seas_driver.cpp          NEW P5
├── params/*.toml                    NEW P5
├── scripts/regression_check.py      NEW P0
└── tests/
    ├── unit/
    │   ├── test_linear_elastic.cpp          NEW P3 (~300 LOC)
    │   ├── test_constitutive_integrator.cpp NEW P3 (~400 LOC)
    │   ├── test_boundary_config.cpp         NEW P4 (~500 LOC)
    │   ├── test_boundary_config_operator.cpp NEW P4 (~400 LOC)
    │   ├── test_boundary_classification.cpp NEW P4 (~300 LOC)
    │   ├── test_seas_config_parser.cpp      NEW P5 (~600 LOC)
    │   ├── test_logging.cpp                 NEW P6 (~100 LOC)
    │   └── ... (existing 24 unit tests)
    ├── parallel/
    │   ├── test_fault_scatter.cpp           NEW P1 (~300 LOC)
    │   └── ... (existing 9 parallel tests)
    └── verification/                        (existing 7)
```

**New test LOC: ~2,900** across 8 new test files.

---

## 9. Execution & Branch Strategy

### 9.1 Summary Table

| Phase | Branch | Risk | Regression | Key Deliverable | New Tests |
|-------|--------|------|-----------|-----------------|-----------|
| 0 | `phase0-regression` | None | Self-check | Golden outputs, regression script | T1, T2, T3 |
| 1 | `phase1-mpi` | Low | Byte-for-byte | FaultScatter, extended MPIContext | T4, T5, T6 |
| 2 | `phase2-decompose` | Low | Byte-for-byte | 6 files from 1 | T7 (incremental rebuild) |
| 3 | `phase3-constitutive` | Medium | Byte-for-byte* | `constitutive/`, integrator overloads | T8, T9 |
| 4 | `phase4-boundary` | Medium | L2 < 1e-10 | BoundaryConfig, per-attr DirichletFunc | T10, T11, T12 |
| 5 | `phase5-toml` | Medium | L2 < 1e-12 | TOML config, generic driver | T13, T14 |
| 6 | `phase6-polish` | Low | Byte-for-byte | BCMode removed, logging, docs | T15 |

*Phase 3: byte-for-byte target; relax to L2 < 1e-15 only if virtual dispatch causes FP reordering.

### 9.2 Branch Strategy

All phase branches fork from and merge to `system_update`. After Phase 6, `system_update` merges to `main`.

```
system_update ─────────────────────────────────────────────────────►  → merge to main
  │
  ├─ tag: v1.0-bp5-verified
  ├─ refactor/phase0-regression ─── merge ─┐
  ├─ refactor/phase1-mpi ────────── merge ─┤
  ├─ refactor/phase2-decompose ──── merge ─┤
  ├─ refactor/phase3-constitutive ─ merge ─┤
  ├─ refactor/phase4-boundary ───── merge ─┤
  ├─ refactor/phase5-toml ──────── merge ──┤
  └─ refactor/phase6-polish ────── merge ──┘
```

Each branch starts from latest `system_update`. Merge only (no rebase).

### 9.3 Cross-Phase Testing Matrix

| Phase | Critical Old Tests | Golden Check | Frontera Pair |
|-------|-------------------|-------------|---------------|
| 0 | `make test` (full baseline) | Self-check at tol 0 | Produce golden |
| 1 | `test_parallel_elasticity`, `test_serial_parallel_consistency`, `test_bp5_parallel_smoke` | Byte-for-byte | Serial-vs-parallel [VERIFY] diff |
| 2 | `test_elasticity_operator`, `test_cross_verify_tandem`, `test_elasticity_br2` | Byte-for-byte | Not needed (binary unchanged) |
| 3 | `test_cross_verify_tandem`, `test_serial_parallel_consistency`, ALL integrator tests | Byte-for-byte target | Run both sbatch |
| 4 | `test_elasticity_operator`, `test_cross_verify_tandem`, `test_bp5_integration` | L2 < 1e-10 | Run both sbatch |
| 5 | All existing (old driver unchanged) | L2 < 1e-12 | New driver on Frontera |
| 6 | All (after BCMode migration) | Byte-for-byte | Final confirmation |

### 9.4 Rollback Cost

| Phase | Rollback Steps | Effort |
|-------|---------------|--------|
| 0 | Zero (additive only) | None |
| 1 | Remove `FaultScatter`, inline MPI calls back into operator. Restore old `MPIContext` constructor. | 30 min |
| 2 | `cat *.inl >> .hpp`, remove includes. | 5 min |
| 3 | (1) Delete `constitutive/` and `domain/boundary_config.hpp`. (2) `git checkout` old integrator constructors. (3) `git checkout` old `ElasticityDomainOperator` constructor. (4) `git checkout` old `ComputeTractionImpl`. Old constructors no longer delegate — they're restored to pre-Phase-3 state. | 30 min with git |
| 4 | Remove per-attr DirichletFunc, MFEM_VERIFY, operator refactor. Old BCMode dual-path still works. | 20 min |
| 5 | Remove TOML parser, driver, submodule. Old CLI driver untouched. | 10 min |
| 6 | Restore BCMode, remove logging wrapper, revert comments. | 30 min |

Phase 3 has the highest rollback cost because it touches integrators, operator constructor, and traction computation. All changes are in tracked files — `git checkout` from pre-Phase-3 tag restores everything.

---

## 10. Dependency & Coupling Analysis

### 10.1 DAG After All Phases

```
common/ ← constitutive/ ← integrator/ ← domain/ ← fault/ ← solver/ ← io/ ← drivers/
              ↑                              ↑
     (abstract, no deps)          boundary_config.hpp (value-type)
                    config/ ← friction/
                      ↑
               toml11 (P5, optional)
```

`constitutive/` depends only on `mfem.hpp`. Framework code never depends on config-parsing or drivers.

---

## 11. Resolved Questions

| Question | Resolution |
|----------|------------|
| Explicit template instantiation? | Defer. Only if compile > 2 min per TU. |
| Quick-check step count? | 20. Validate with 2x Vp perturbation. |
| Test file splitting? | Defer. |
| FaceVertexKey for hex? | Defer. Code comment. |
| Config validation depth? | Minimum + mesh-attr check. |
| BCMode deprecation? | Phase 4: dual-path. Phase 6: removed. |
| toml11 vendoring? | Git submodule, optional `SEAS_USE_TOML`. |
| Phase order? | 0-1-2-3(constitutive)-4(BoundaryConfig)-5(TOML)-6(polish). |
| SEASConfig behavior? | Pure data struct (D9). |
| DirichletFunc lifetime? | Capture by value (D8). |
| DG method preference? | SIPG primary. Gated on IP-vs-BR2 test (B2). |
| Constitutive in which phase? | Phase 3. Mandatory. |
| GetLambda/GetMu on interface? | Replaced by GetPenaltyModulus/GetMaxWaveSpeed (D4). |
| TensorDotNormal on interface? | Removed. Derived from ComputeTangent (D3). |
| Constructor overloads? | Final signature in Phase 3, only 2 overloads (D1). |
| ConstitutiveModel ownership? | const& in, const* stored, old ctor owns unique_ptr (D2). |
| Name-based BC? | Deferred (B1). MFEM has no Gmsh name API. Numeric only for now. |
| IP-vs-BR2 default change? | Gated on explicit comparison (B2, S5). |

---

## 12. Future Architecture: Dynamic Rupture, Friction, Multi-Physics

### 12.1 Three Simulation Modes

The code must support three distinct simulation modes, user-selectable via config:

```toml
[simulation]
mode = "hybrid"   # "qd" | "dynamic" | "hybrid"
```

| Mode | Description | Use Case |
|------|-------------|----------|
| **`qd`** (Quasi-Dynamic) | Current system. Implicit solve, radiation damping, adaptive RK45. | Earthquake cycles, interseismic, slow slip. No wave propagation. |
| **`dynamic`** (Dynamic Rupture) | Velocity-stress system, Riemann solver, explicit time stepping. No QD fallback. | Single-event rupture simulation, wave propagation studies. |
| **`hybrid`** (QD-Dynamic) | Starts in QD. Switches to FD when V_max > threshold. Returns to QD when rupture ends. | Full SEAS cycles with resolved coseismic rupture and wave propagation. |

Each mode uses the same mesh, the same `ConstitutiveModel`, the same `FaultInterfaceLaw`, and the same `FaultBasis`/`FaultGeometry`. The difference is the domain solver and time integration strategy.

#### 12.1.1 Mode-Specific Operator Selection

```cpp
// In the driver (or seas_hybrid_operator.hpp):
if (config.mode == "qd") {
   // Current system — no changes needed
   auto seas_op = SEASQuasiDynamicOperator(domain_qd, fault, mpi);
   time_stepper.Run(seas_op, state, t_final);
}
else if (config.mode == "dynamic") {
   // FD only — no QD operator created, no stiffness matrix
   auto wave_op = WaveOperator(mesh, model, bdr_config);
   auto fd_stepper = ADERTimeStepper(wave_op, fault, cfl);
   fd_stepper.Run(state, t_final);
}
else if (config.mode == "hybrid") {
   // Both systems, with regime switching
   auto hybrid_op = SEASHybridOperator(domain_qd, wave_op, fault, mpi, thresholds);
   hybrid_stepper.Run(hybrid_op, state, t_final);
}
```

For `"dynamic"` mode, the QD stiffness matrix is never assembled — saving significant memory and setup time. For `"qd"` mode, the wave operator is never created. Only `"hybrid"` instantiates both.

#### 12.1.2 Hybrid Mode: Regime Switching Protocol

```
QD regime (interseismic):
  - Monitor V_max via MPI_Allreduce(MPI_MAX) every RK stage
  - When V_max > V_activate (e.g., 1e-3 m/s):
      1. FREEZE QD state: save (displacement, slip, psi) to in-memory snapshot
      2. Convert displacement → velocity-stress (element-local, GPU-friendly):
           sigma = ConstitutiveModel->ComputeStress(grad(u))
           v = 0 (or smooth estimate)
      3. Initialize absorbing boundaries
      4. Checkpoint to disk (crash recovery point)
      5. Switch regime = FD

FD regime (coseismic):
  - Explicit time stepping with CFL: dt_FD ~ h / c_p
  - Riemann solver at interfaces, friction constraint at fault
  - Monitor V_max
  - When V_max < V_deactivate (e.g., 1e-6 m/s) AND seismic energy exited:
      1. RESTORE QD snapshot from memory
      2. Update: slip_new = slip_frozen + integral(V_FD dt)
      3. Update: psi_new from last FD step
      4. One QD Solve() to recompute traction from updated slip
      5. Switch regime = QD, resume adaptive RK45

Dynamic-only mode:
  - No switching. FD runs from t=0 to t_final.
  - Initial conditions: prescribed stress field + nucleation perturbation
  - Absorbing BCs from the start
```

#### 12.1.3 State Persistence During Hybrid Switching

**QD state during FD: keep in memory.**

| Data | Size (p1/1000m, 400 ranks) | Strategy |
|------|---------------------------|----------|
| Displacement u | ~24 MB | In-memory snapshot |
| Slip + psi | ~1.2 MB | In-memory |
| Stiffness matrix | ~1 GB sparse | Keep allocated (reassembly too expensive) |

Checkpoint to disk at QD→FD transition for crash recovery. If FD crashes, restart from checkpoint.

#### 12.1.4 Data Transfer & GPU Potential

QD→FD conversion (displacement → velocity-stress) is element-local: compute grad(u), apply `ConstitutiveModel::ComputeStress()`. No MPI. Embarrassingly parallel — ideal GPU kernel.

FD explicit time stepping is also GPU-friendly: element-local mass inverse, face-local Riemann solve. This is where GPU acceleration has highest impact.

#### 12.1.5 Fault State Continuity

The critical invariant: **fault slip and psi are continuous across regime transitions.**

Prevention of wave pollution at QD→FD switch:
1. Stress-displacement consistency via shared `ConstitutiveModel::ComputeStress()`
2. Initial velocity = 0 (not noisy finite differences)
3. Ramp-up period: absorbing BCs damp artificial waves during first few wave-crossing times
4. FD→QD: slip and psi continuous; displacement reconstructed via one QD Solve()
5. Fault state is the coupling variable — domain fields are derived

#### 12.1.6 Config Parameters (Hybrid Mode)

```toml
[simulation]
mode = "hybrid"

[simulation.hybrid]
V_activate = 1e-3      # V_max threshold to switch QD → FD [m/s]
V_deactivate = 1e-6    # V_max threshold to switch FD → QD [m/s]
wave_exit_time = 0.0    # Extra wait after V < V_deactivate before switching back [s]
                         # 0 = auto (domain_size / c_p)
checkpoint_at_switch = true   # Write checkpoint at each QD→FD transition

[simulation.dynamic]
# Used when mode = "dynamic" or during FD phase of hybrid
cfl_factor = 0.5        # CFL number (dt = cfl * h / c_p)
time_integrator = "ader" # "ader" or "rk4"
```

#### 12.1.7 Proposed Files

```
dynamic/
├── wave_operator.hpp              # Velocity-stress DG operator
├── riemann_solver.hpp             # Exact Riemann solver for elastic waves
├── absorbing_bc.hpp               # Absorbing BCs
├── free_surface_bc.hpp            # Free-surface BC via inverse Riemann
├── fault_riemann.hpp              # Fault Riemann solver + friction constraint
├── ader_time_integrator.hpp       # ADER time integration
└── velocity_stress_state.hpp      # Q = (sigma, v) container + conversions

solver/
└── seas_hybrid_operator.hpp       # Three-mode dispatch (qd / dynamic / hybrid)
```

### 12.2 Extended Friction Laws (R&S Variants)

*(Unchanged from v3 — focus on rate-and-state variants only: flash heating, thermal pressurization. No slip-weakening.)*

### 12.3 Future Nonlinear Constitutive Models

Phase 3 establishes the `ConstitutiveModel` interface with three tiers (algebraic → evolving → non-local). Future implementations use the higher tiers:

| Model | File | Tier | Internal Vars | Key Feature |
|-------|------|------|--------------|-------------|
| Linear elastic (Phase 3) | `linear_elastic.hpp` | 1 | 0 | Baseline |
| Viscoelastic (Maxwell) | `viscoelastic.hpp` | 2 | 6 (viscous strain) | `ComputeStateRates`: dε_v/dt |
| Elastoplastic (J2) | `elastoplastic.hpp` | 1 | 7 (plastic strain + hardening) | `UpdateState`: return mapping |
| Damage-breakage (local) | `damage_breakage.hpp` | 2 | 8 (α, B, 6 ε^p) | `ComputeStateRates`: dα/dt, dB/dt, dε^p/dt (Lyakhovsky & Ben-Zion 2014, Eqs. 25-26, 28b) |
| Damage-breakage (non-local) | `damage_breakage.hpp` | 3 | 8 + FE field for α | `ComputeStressNonLocal`: adds -ϑ∇ᵢα·∇ⱼα structural stresses (Eq. 14) |

**Damage-breakage implementation note:** The damage and breakage parameters have their own evolution equations that the time integrator must advance. The moduli λ(α), μ(α), γ(α) change with damage state (Eq. 12), so the stiffness matrix is reassembled at each Newton step when `IsNonlinear()` is true. For the non-local variant, α requires its own DG FE space for gradient computation — the domain operator detects this via `NumNonLocalVars() > 0` and allocates accordingly.

**Interface stability note:** The Tier 2/3 methods (`HasStateEvolution`, `ComputeStateRates`, `NumNonLocalVars`, `ComputeStressNonLocal`) are designed from the damage-breakage equations (Lyakhovsky & Ben-Zion 2014) but have not been exercised by framework code — `LinearElastic` returns false/0 for all of them. When implementing the first Tier 2 model, expect the interface to need minor adjustments (e.g., `ComputeStateRates` may need access to the full displacement field or neighboring quadrature points, not just local strain). This is acceptable — the interface will stabilize after the first real consumer.

**Performance note (S7):** For Tier 1 nonlinear models (J2 plasticity), `Solve()` becomes a Newton iteration loop. Stiffness re-assembly via `ComputeTangent()` occurs per Newton step per RK stage — changes assembly from ~5% to ~30-50% of wall time. For Tier 2 evolving models (damage-breakage), the constitutive state ODEs add ~N_qp × NumInternalVars additional unknowns to the global ODE system. For a 1000m mesh with ~500k quadrature points and 8 internal vars, that's ~4M additional DOFs in the time integrator state vector. This is significant but manageable — the ODE RHS evaluation (step 7 in the Mult loop) is embarrassingly parallel across quadrature points.

### 12.4 Multi-Physics: Solid-Fluid Coupling (Note)

*(Unchanged from v3 — two/three-field formulations (u,p) or (u,p,T), coupling terms, not urgent, structure noted.)*

### 12.5 Design Rules

1. One file per model.
2. Abstract base → concrete implementations.
3. No cross-directory includes between peer modules.
4. Factory pattern for runtime selection.
5. State variable layout owned by the model.
6. SIPG primary DG method; BR2 preserved for comparison.
7. Driver owns the coupling loop.

### 12.6 What Phases 0-6 Must NOT Do

| Rule | Rationale |
|------|-----------|
| Don't hardcode `dim=3` or `num_slip_comp=2` | Future 2D dynamic reuses same abstractions |
| Don't embed `(lambda,mu)` deeper — use `ConstitutiveModel` | Phase 3 addresses this |
| Don't merge `RateStateFaultOperator` with friction logic | Separate orchestration from computation |
| Don't close `DomainOperator` base to extension | FD wave operator needs sibling interface |
| Don't assume QD is the only mode | Three modes (qd/dynamic/hybrid) must coexist; config selects which operators are instantiated |
| Don't store `Vp_` as authoritative Dirichlet source | Phase 4's DirichletFunc callback replaces it |
| Keep `FaultBasis`, `FaultGeometry` solver-independent | Both QD and FD use same transforms |
| Don't free stiffness matrix during FD | Reassembly cost >> memory cost |

---

## 13. Phase 7: New Driver Full Simulation & Verification

**Date Added:** 2026-04-10
**Goal:** Wire `seas_driver.cpp` to run a complete BP5 simulation using all new code paths (ConstitutiveModel, BoundaryConfig, DomainConfig, TOML config) and verify the new setup produces results identical to the existing verification driver (`bp5_verification_full.cpp`).
**Branch:** `refactor/phase7-driver-verification`
**Regression:** L2 < 1e-12 against existing golden data; target byte-for-byte against old driver at identical parameters.

### 13.0 Overview

After 7 phases of refactoring, the new infrastructure is in place but untested end-to-end:
- **Phase 3** added `ConstitutiveModel` + `LinearElastic` + new operator constructor
- **Phase 4** added `BoundaryConfig` with parameterized boundary attrs and `DirichletFunc`
- **Phase 5** added TOML parsing, `SEASConfig`, and `seas_driver.cpp`
- **Phase 6** added unified BoundaryConfig path, logging cleanup

But `seas_driver.cpp` currently only **parses config and prints it** (line 122: `"Full simulation pipeline will be wired in Phase 6."`). It does NOT construct operators, run time integration, or produce output. This phase completes the driver and proves the new code paths are correct.

### 13.1 Gap Analysis

| Component | Old Path (`bp5_verification_full.cpp`) | New Path (`seas_driver.cpp`) | Gap |
|-----------|---------------------------------------|------------------------------|-----|
| **Config source** | CLI flags + hardcoded `BP5Params` | TOML file → `SEASConfig` | Need `SEASConfig → BP5Params` bridge |
| **Material** | Scalar `lambda, mu` from `BP5Params` | `LinearElastic(lambda, mu)` from `MaterialConfig` | Need construction from config |
| **Boundary** | `BCMode::FarField` enum → legacy internal build | `BoundaryConfig{dirichlet_attrs, natural_attrs, fault_attr, DirichletFunc}` from TOML | Need construction from config |
| **Domain operator** | Legacy constructor (11 params) | New constructor (ConstitutiveModel + BoundaryConfig + DomainConfig) | Need to wire new constructor |
| **DG method** | `DGMethod` enum from CLI `--solver` | String `"IP"` / `"BR2"` from TOML | Need string→enum conversion |
| **Solver type** | `SolverType` enum from CLI `--solver` | String `"mumps"` / `"cg"` from TOML | Need string→enum conversion |
| **Fault geometry** | `FaultGeometry(domain, BP5Params, mpi)` | Same, but params from `SEASConfig` | Need `SEASConfig → BP5Params` bridge |
| **Fault operator** | `RateStateFaultOperator<ParMesh, 2>` with `BP5Params` | Same | Same bridge needed |
| **I/O** | `ParallelBP5BenchmarkOutput` | Same | Need station list + fault coords setup |
| **Time integration** | `DormandPrinceRK45` with CFL-aware dt | Same | Need dt init from config |
| **Checkpoint** | `ReadCheckpoint` / `WriteCheckpoint` | Same | Need restart CLI flag |

**Core blocker:** `FaultGeometry`, `RateStateFaultOperator`, and `ParallelBP5BenchmarkOutput` all accept `BP5Params` as constructor arguments. Rather than refactoring those constructors (high-risk, many callers), we create a `SEASConfig → BP5Params` bridge function.

### 13.2 Constraints

- **Interface constraints:** `FaultGeometry`, `RateStateFaultOperator`, `ParallelBP5BenchmarkOutput`, `DormandPrinceRK45` constructors are untouched.
- **Dependency constraints:** Must use existing `BP5Params` struct for components that require it. No parallel refactoring of fault/IO code.
- **Convention constraints:** Follow existing BP5 verification driver patterns for I/O, earthquake detection, console output.
- **Numerical constraints:** Output must match old driver at identical parameters. L2 < 1e-12 against golden data. Target byte-for-byte.

### 13.3 Phase 7a: SEASConfig → BP5Params Bridge & Enum Mapping

#### Goal
Provide utility functions that convert between `SEASConfig` fields and the types expected by existing components, without modifying those components.

#### Files to Create
- `config/seas_config_bridge.hpp` — utility functions for SEASConfig → BP5Params conversion and string→enum mapping

#### Detailed Requirements

1. **`BP5Params BuildBP5Params(const SEASConfig &config)`**

   Maps SEASConfig fields to BP5Params fields. Field mapping:

   ```cpp
   BP5Params BuildBP5Params(const SEASConfig &config)
   {
      BP5Params p;
      p.rho = config.material.density;
      p.cs = config.material.cs;
      p.nu = config.material.nu;
      p.V0 = config.friction.V0;
      p.f0 = config.friction.f0;
      p.b = config.friction.b;
      p.L0 = config.friction.L0;
      p.L_nuc = config.friction.L_nuc;
      p.a0 = config.friction.a0;
      p.amax = config.friction.amax;
      p.sigma_n = config.friction.sigma_n;
      p.Vp = config.loading.Vp;
      p.V_init = config.loading.V_init;
      p.V_nuc = config.loading.V_nuc;
      p.delta_tau_factor = config.loading.delta_tau_factor;
      p.nucleation_eps = config.loading.nucleation_eps;
      p.smooth_nucleation = config.loading.smooth_nucleation;
      p.Wf = config.fault_geom.Wf;
      p.lf = config.fault_geom.lf;
      p.hs = config.fault_geom.hs;
      p.ht = config.fault_geom.ht;
      p.H = config.fault_geom.H;
      p.l_vw = config.fault_geom.l_vw;
      p.w_nuc = config.fault_geom.w_nuc;
      p.t_final = config.time.t_final;
      return p;
   }
   ```

   Every field in `BP5Params` that has a corresponding `SEASConfig` field must be mapped. Fields not in `SEASConfig` (e.g., `V_zero = 1e-20`, `seconds_per_year`) keep their `BP5Params` defaults.

2. **`DGMethod ParseDGMethod(const std::string &s)`**

   Only accept the same strings that `SEASConfigParser::Validate()` accepts (`"IP"`, `"BR2"`), plus lowercase variants. Do NOT add aliases like `"SIPG"` that Validate() rejects — they would be dead code since Validate() runs upstream and aborts on unrecognized values.

   ```cpp
   DGMethod ParseDGMethod(const std::string &s)
   {
      if (s == "IP" || s == "ip") return DGMethod::IP;
      if (s == "BR2" || s == "br2") return DGMethod::BR2;
      MFEM_ABORT("Unknown DG method: " << s);
      return DGMethod::BR2;
   }
   ```

3. **`SolverType ParseSolverType(const std::string &s)`**

   ```cpp
   SolverType ParseSolverType(const std::string &s)
   {
      if (s == "cg" || s == "CG") return SolverType::CG_AMG;
      if (s == "mumps" || s == "MUMPS") return SolverType::MUMPS;
      if (s == "mumps-blr" || s == "MUMPS-BLR") return SolverType::MUMPS_BLR;
      if (s == "superlu") return SolverType::SUPERLU;
      if (s == "strumpack") return SolverType::STRUMPACK;
      if (s == "gmres") return SolverType::GMRES_AMG;
      MFEM_ABORT("Unknown solver type: " << s);
      return SolverType::CG_AMG;
   }
   ```

4. **`BoundaryConfig BuildBoundaryConfig(const BoundaryTomlConfig &bdr, real_t Vp)`**

   ```cpp
   BoundaryConfig BuildBoundaryConfig(const BoundaryTomlConfig &bdr, real_t Vp)
   {
      BoundaryConfig bc;
      bc.dirichlet_attrs = bdr.dirichlet_attrs;
      bc.natural_attrs = bdr.natural_attrs;
      bc.fault_attr = bdr.fault_attr;
      bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);
      return bc;
   }
   ```

5. **`DomainConfig BuildDomainConfig(const SolverConfig &solver)`**

   ```cpp
   DomainConfig BuildDomainConfig(const SolverConfig &solver)
   {
      DomainConfig dc;
      dc.face_basis_type = solver.face_basis_type;
      dc.penalty_factor = solver.penalty_factor;
      dc.blr_tol = solver.blr_tol;
      dc.check_residual = solver.check_residual;
      return dc;
   }
   ```

#### Interfaces
- All functions are `inline` in the header, no .cpp file needed
- Include `seas_config.hpp`, `bp5_params.hpp`, `boundary_config.hpp`, `domain_config.hpp`, `elasticity_operator.hpp` (for DGMethod, SolverType)

#### Edge Cases
- `SEASConfig` with `benchmark=""` (no preset): all fields must be explicitly set — `BuildBP5Params` still works (copies whatever values are in SEASConfig)
- `SEASConfig` with non-BP5 benchmark: `BuildBP5Params` still works but caller should validate the result makes physical sense
- Empty `dirichlet_attrs` or `fault_attr==0`: caught by `SEASConfigParser::Validate()` upstream, but `BuildBoundaryConfig` does not re-validate

#### Unit Test: `tests/unit/test_seas_config_bridge.cpp` (~200 LOC)

Every previous phase created explicit test files for new code (Phase 1: `test_fault_scatter.cpp`, Phase 3: `test_linear_elastic.cpp`, Phase 4: `test_boundary_config.cpp`, Phase 5: `test_seas_config_parser.cpp`). Phase 7a follows the same pattern.

Test cases:

1. **`BuildBP5Params` round-trip:** Construct `SEASConfig`, call `ApplyBP5Defaults`, then `BuildBP5Params` → every field matches `BP5Params()` default constructor. Compare all 25+ fields individually (not just a few spot-checks).

2. **`BuildBP5Params` with non-preset values:** Construct `SEASConfig` with `material.density=3000`, `friction.b=0.05`, etc. → verify those specific fields in the returned `BP5Params`, and verify unmapped fields (e.g., `V_zero`, `seconds_per_year`) retain their `BP5Params` defaults.

3. **`ParseDGMethod`:** `"IP"` → `DGMethod::IP`, `"ip"` → `DGMethod::IP`, `"BR2"` → `DGMethod::BR2`, `"br2"` → `DGMethod::BR2`. (No "SIPG" — not a valid Validate() input.)

4. **`ParseSolverType`:** All 6 valid strings → correct enum value: `"cg"`, `"mumps"`, `"mumps-blr"`, `"superlu"`, `"strumpack"`, `"gmres"`.

5. **`BuildBoundaryConfig`:** Verify `dirichlet_attrs`, `natural_attrs`, `fault_attr` copied from `BoundaryTomlConfig`. Verify `default_dirichlet_func` is non-null and produces correct displacement at a test point (e.g., `y=2000, t=1e9` → `u(0) = Vp*t/2`).

6. **`BuildDomainConfig`:** Verify all `SolverConfig` fields map to corresponding `DomainConfig` fields: `face_basis_type`, `penalty_factor`, `blr_tol`, `check_residual`.

Add to `Makefile`: `test-config-bridge` target, included in `make test`.

#### Acceptance Criteria
- [ ] `BuildBP5Params(config)` round-trips: `ApplyBP5Defaults` → `BuildBP5Params` → every field matches `BP5Params()` default
- [ ] `ParseDGMethod("IP") == DGMethod::IP`, `ParseDGMethod("BR2") == DGMethod::BR2`
- [ ] `ParseSolverType("mumps") == SolverType::MUMPS`, all 6 solver types covered
- [ ] `BuildBoundaryConfig` produces a `BoundaryConfig` with correct attrs and non-null `default_dirichlet_func`
- [ ] `test_seas_config_bridge` compiles and all tests pass under `make test`

#### Dependencies
- Depends on: Phase 5 (SEASConfig, BoundaryTomlConfig exist), Phase 4 (BoundaryConfig, MakeBP5DirichletFunc exist)
- Required by: Phase 7b

---

### 13.4 Phase 7b: Wire Full Simulation Pipeline in `seas_driver.cpp`

#### Goal
Transform `seas_driver.cpp` from a config-printer into a fully functional BP5 simulation driver using all new code paths.

#### Files to Modify
- `drivers/seas_driver.cpp` — replace the config-print stub with the full simulation pipeline

#### Detailed Requirements

The driver must replicate the exact simulation flow of `bp5_verification_full.cpp` but using the new constructor paths. The pipeline stages are:

**Stage 1: Config & Mesh (lines 35-66 of current driver, extend)**

```cpp
// Parse TOML config (already implemented)
SEASConfig config = SEASConfigParser::ParseFile(config_file);
SEASConfigParser::ApplyCLIOverrides(config, overrides);
SEASConfigParser::Validate(config);

// Build bridge objects
BP5Params params = BuildBP5Params(config);
DGMethod dg_method = ParseDGMethod(config.solver.dg_method);
SolverType solver_type = ParseSolverType(config.solver.solver_type);
BoundaryConfig bdr_config = BuildBoundaryConfig(config.boundary, config.loading.Vp);
DomainConfig domain_config = BuildDomainConfig(config.solver);
```

**Stage 2: Mesh loading + mesh characteristics**

```cpp
// Load mesh from file specified in TOML
auto serial_mesh = std::make_unique<Mesh>(config.mesh.file.c_str(), 1, 1);
if (config.mesh.scale != 1.0)
{
   serial_mesh->Transform([scale=config.mesh.scale](const Vector &x, Vector &p)
   {
      p.SetSize(x.Size());
      for (int i = 0; i < x.Size(); i++) { p(i) = x(i) * scale; }
   });
}
ParMesh pmesh(mpi.GetComm(), *serial_mesh);
serial_mesh.reset();

// Compute mesh characteristics (h_min needed for CFL-aware dt in Stage 7)
real_t h_min, h_max, kappa_min, kappa_max;
pmesh.GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
MFEM_VERIFY(h_min > 0, "Degenerate mesh: h_min = " << h_min);
```

**Stage 3: Domain operator (NEW constructor path)**

```cpp
LinearElastic material(config.material.lambda(), config.material.mu());

ElasticityDomainOperator<ParMesh> domain(
   pmesh, config.mesh.order,
   material,                    // ConstitutiveModel& (Phase 3)
   config.loading.Vp,
   config.fault_geom.Wf, config.fault_geom.lf,
   bdr_config,                  // BoundaryConfig (Phase 4)
   dg_method, solver_type,
   domain_config);              // DomainConfig (Phase 3)
```

This is the **critical line** — it exercises the new constructor that was never called by a real simulation driver before.

**Stage 4: Fault components**

```cpp
FaultGeometry<ParMesh> fault_geom(domain, params, &mpi);

DieterichRuinaFriction::Constants fc;
fc.V0 = params.V0; fc.f0 = params.f0;
fc.b = params.b; fc.Dc = params.L0;
DieterichRuinaFriction friction(fc);
AgingLawPsi aging(params.b, params.V0, params.f0);

RateStateFaultOperator<ParMesh, 2> fault_op(
   &fault_geom, &friction, &aging, params, &mpi);
```

**Stage 5: SEAS operator + initial condition**

```cpp
PBP5SEASOp seas_op(&domain, &fault_op, &mpi);
seas_op.SetElasticSigmaN(true);  // Default ON per Tandem convention

Vector state(fault_op.StateSize());
seas_op.SetInitialCondition(state);
real_t V_init = seas_op.GetMaxSlipRate();
```

**Stage 6: I/O setup**

Use `BP5BenchmarkOutput<MeshType>::DefaultStations()` — do NOT hardcode a station list. The 10 SCEC BP5 on-fault stations (Section 4.1 of spec) are defined in `io/bp5_benchmark_output.hpp:811-824`:

```cpp
// Correct: use the canonical station list from the output class
auto stations = BP5BenchmarkOutput<ParMesh>::DefaultStations();
```

These are (for reference, from `DefaultStations()`):
```
fltst_strk-36dp+00  (-36 km,  0 km)
fltst_strk-16dp+00  (-16 km,  0 km)
fltst_strk+00dp+00  (  0 km,  0 km)
fltst_strk+16dp+00  ( 16 km,  0 km)
fltst_strk+36dp+00  ( 36 km,  0 km)
fltst_strk-24dp+10  (-24 km, 10 km)
fltst_strk-16dp+10  (-16 km, 10 km)
fltst_strk+00dp+10  (  0 km, 10 km)
fltst_strk+16dp+10  ( 16 km, 10 km)
fltst_strk+00dp+22  (  0 km, 22 km)
```

This is exactly what `bp5_verification_full.cpp` uses (line 1034: `auto stations = BP5BenchmarkOutput<Mesh>::DefaultStations()`). The station names must match exactly for golden data comparison — the comparison script matches files by station name suffix.

I/O objects:
```cpp
ParallelBP5BenchmarkOutput bench_out(
   full_prefix, params, stations, fault_geom, mpi,
   local_x2, local_x3, local_tp_dip, local_tp_strike,
   domain.GetNbfPerFace(), face_basis_type);

// Global output
std::unique_ptr<ProbeOutput> global_out;
if (mpi.IsRoot()) {
   global_out = std::make_unique<ProbeOutput>(
      full_prefix + "_global.txt",
      std::vector<std::string>{"time(s)", "log10(Vmax)(m/s)"},
      "BP5-QD global output");
}
```

**Stage 7: Time integration**

**CRITICAL: `tandem_time_stepping` policy.** When `config.time.tandem_time_stepping == true` (the TOML default, `seas_config.hpp:108`), the old driver (`bp5_verification_full.cpp` lines 879-884) overrides three settings:
- `dt_init = 0.01` (fixed, ignoring CFL formula)
- V-guard = OFF
- psi-clamp = OFF (not applicable here — psi-clamp is a fault operator setting)

The new driver must mirror this exactly, otherwise the 7d equivalence test will fail due to different initial conditions.

```cpp
DormandPrinceRK45 ode_solver;
ode_solver.SetMPIContext(&mpi);
ode_solver.SetAbsTol(config.time.atol);
ode_solver.SetRelTol(config.time.rtol);
ode_solver.SetDtMin(1e-6);
ode_solver.SetDtMax(0.1 * BP5Params::seconds_per_year);

// CFL-aware initial dt computation (always computed for logging)
int dim = 3;
real_t c_N_1 = config.mesh.order * (config.mesh.order + dim - 1.0) / dim;
real_t c_N_1_ref = 2.0 * (2.0 + dim - 1.0) / dim;
real_t beta = 4.0 * c_N_1 / c_N_1_ref;
real_t V_max_init = std::max(V_init, params.V_nuc);
real_t dt_V = std::min(1e3, 0.01 * params.L_nuc / std::max(V_max_init, 1e-20));
real_t dt_CFL = 2.0 * params.eta() * h_min / (beta * params.mu());
real_t dt_init = std::min(dt_V, dt_CFL);

// Tandem-style time stepping overrides (mirrors bp5_verification_full.cpp:879-884)
bool use_v_guard = true;   // default: V-guard ON
if (config.time.tandem_time_stepping)
{
   dt_init = 0.01;         // Tandem's fixed startup dt (tandem_dt_init)
   use_v_guard = false;     // Tandem doesn't use V-guard
}

ode_solver.SetDt(dt_init);
ode_solver.SetStatePerNode(3);  // BP5: [slip_dip, slip_strike, psi]
if (use_v_guard)
{
   ode_solver.SetVGuard(100.0);
}
ode_solver.Init(seas_op);
```

**Why this matters for equivalence:** The CFL-based dt at p=1/1000m is typically ~O(1) seconds, while Tandem's dt=0.01 is much smaller. V-guard rejects steps where V grows too fast — with it OFF, the RK45 error control alone governs acceptance. Both differences compound from step 1, making outputs diverge immediately if not matched.

**Stage 8: Time loop**

Standard time loop with:
- Step rejection handling
- Earthquake detection (V_threshold_seismic = 1e-3, V_threshold_interseismic = 1e-6)
- bench_out.Write() + bench_out.Flush() per accepted step
- global_out per step (root only)
- Periodic console output
- Checkpoint writing at configurable intervals
- NaN/Inf detection with early termination

Match the exact logic from `bp5_verification_full.cpp` lines 1849-2100 (the main time loop).

**Stage 9: Final output**

```cpp
bench_out.ForceWrite(t, state, fault_op, seas_op.GetTraction(), V_max);
bench_out.Close();
// Summary statistics
```

#### Additional CLI Overrides for the New Driver

Add support for common CLI flags that `bp5_verification_full.cpp` supports, passed after the TOML file:

```
--override time.t_final=1e8     # TOML field override (already supported)
--max-steps N                   # Shorthand for time.max_steps
--tandem-time-stepping          # Enable Tandem-style dt defaults
--restart PREFIX                # Restart from checkpoint
--write-every-step              # Write I/O at every accepted step
```

These are convenience aliases that translate to `SEASConfig` field modifications before the simulation starts.

#### Interfaces
- `seas_driver.cpp` includes: `seas_config_bridge.hpp`, `elasticity_operator.hpp`, `fault_geometry.hpp`, `rate_state_fault.hpp`, `dieterich_ruina.hpp`, `state_evolution.hpp`, `seas_operator.hpp`, `time_stepper.hpp`, `bp5_parallel_output.hpp`, `probe_output.hpp`, `checkpoint.hpp`
- No new classes; the driver is a `main()` function using existing components via the new constructor paths

#### Edge Cases
- Mesh file not found → MFEM's Mesh constructor throws, caught by driver
- Mesh scale = 1.0 (already in meters) → skip Transform
- `h_min = 0` from degenerate mesh �� MFEM_VERIFY guard on dt_CFL computation
- Zero fault DOFs on a rank → `FaultGeometry` handles this (existing code)
- Restart with different rank count than checkpoint → `ReadCheckpoint` handles or fails cleanly

#### Acceptance Criteria
- [ ] `seas_driver.cpp` compiles with `make seas_driver`
- [ ] `mpirun -np 1 ./seas_driver bp5_example.toml --override time.t_final=100` runs and produces station output files
- [ ] `mpirun -np 8 ./seas_driver bp5_example.toml --override time.max_steps=20` produces output
- [ ] Station file format matches SCEC 8-column spec (time, slip_strike, slip_dip, log10_V_strike, log10_V_dip, tau_strike, tau_dip, log10_state)
- [ ] Global output file produced with time + log10(Vmax)
- [ ] Console output includes step, time, dt, V_max, earthquake count
- [ ] Driver exits cleanly on NaN/Inf with nonzero exit code

#### Dependencies
- Depends on: Phase 7a (bridge functions)
- Required by: Phase 7c, 7d

---

### 13.5 Phase 7c: Verification TOML Config & Build Targets

#### Goal
Create TOML configs and Makefile targets that enable testing the new driver against golden data.

#### Files to Create
- `config/bp5_verification_50step.toml` — matches bp5_verification_full.cpp's 50-step Frontera configuration exactly
- `config/bp5_quick_check.toml` — matches the inline-mesh 8-rank/20-step quick-check configuration

#### Files to Modify
- `Makefile` — add targets for new driver tests

#### Detailed Requirements

1. **`config/bp5_verification_50step.toml`**

   ```toml
   # Matches bp5_verification_full.cpp default configuration
   # for golden data comparison (50 accepted steps, Frontera mesh)

   benchmark = "bp5"

   [mesh]
   file = "bp5/mesh/reference/bp5_tandem_exact.msh"
   scale = 1000
   order = 1

   [boundary]
   dirichlet = [5]
   natural = [1]
   fault = 3

   [solver]
   dg_method = "IP"
   solver_type = "mumps"
   blr_tol = 1e-12

   [time]
   t_final = 56844000000    # ~1800 yr (will override for 50-step test)
   tandem_time_stepping = true
   max_steps = 50

   [output]
   output_dir = "."
   output_prefix = "bp5_driver_verify"
   ```

2. **`config/bp5_quick_check.toml`**

   ```toml
   # Quick-check: inline mesh, 20 steps, local verification
   benchmark = "bp5"

   [mesh]
   file = "__inline__"     # Sentinel: driver uses CreateBP5InlineMesh()
   scale = 1.0             # Inline mesh is already in meters
   order = 1

   [boundary]
   dirichlet = [5]
   natural = [1]
   fault = 3

   [solver]
   dg_method = "IP"
   solver_type = "mumps"
   blr_tol = 1e-12

   [time]
   max_steps = 20
   tandem_time_stepping = true

   [output]
   output_dir = "bp5_driver_quick_check"
   output_prefix = "bp5_driver_quick"
   ```

   **Note:** The `file = "__inline__"` sentinel tells the driver to call `CreateBP5InlineMesh()` instead of loading from file. This requires adding inline-mesh support to the driver (Phase 7b should handle this with an `if (config.mesh.file == "__inline__")` check).

3. **Makefile targets**

   ```makefile
   # New driver quick-check (exercises new constructor path)
   driver-quick-check: seas_driver
       @mkdir -p bp5_driver_quick_check
       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 8 ./seas_driver \
           config/bp5_quick_check.toml
       @echo "driver-quick-check: PASS (no crash)"

   # New driver vs old driver equivalence test
   # IMPORTANT: Both drivers must use the same solver type.
   # New driver: TOML default is solver_type = "mumps" (→ SolverType::MUMPS).
   # Old driver: default is SolverType::MUMPS_BLR. Must pass --solver mumps explicitly.
   driver-equivalence: seas_driver seas_bp5_full
       @mkdir -p bp5_driver_equiv_new bp5_driver_equiv_old
       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 8 ./seas_driver \
           config/bp5_quick_check.toml \
           --override output.output_dir=bp5_driver_equiv_new \
           --override output.output_prefix=bp5_equiv_new
       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 8 ./seas_bp5_full \
           --inline-mesh --max-steps 20 --tandem-time-stepping \
           --solver mumps \
           --output-dir bp5_driver_equiv_old --output-prefix bp5_equiv_old
       python3 scripts/compare_station_files.py \
           bp5_driver_equiv_new bp5_driver_equiv_old --tolerance 1e-12
   ```

#### Acceptance Criteria
- [ ] `make seas_driver` builds successfully
- [ ] `make driver-quick-check` runs without crash, produces station files
- [ ] TOML configs parse without error (validated by `SEASConfigParser::Validate`)

#### Dependencies
- Depends on: Phase 7b (driver runs)
- Required by: Phase 7d

---

### 13.6 Phase 7d: Equivalence Testing — Old Path vs New Path

#### Goal
Prove that the new constructor path (ConstitutiveModel + BoundaryConfig + DomainConfig) produces bit-identical or near-identical results to the legacy constructor path (scalar lambda/mu + BCMode) when given the same parameters.

#### Files to Create
- `scripts/compare_station_files.py` — compares two sets of station output files field-by-field
- `tests/verification/bp5_driver_equivalence.cpp` — C++ test that constructs both paths side-by-side and compares

#### Detailed Requirements

1. **`scripts/compare_station_files.py`**

   ```
   Usage: compare_station_files.py DIR_A DIR_B [--tolerance TOL] [--prefix-a PA] [--prefix-b PB]

   For each station file found in DIR_A:
     1. Find matching station in DIR_B (by station name suffix)
     2. Load both as 8-column SCEC format
     3. Interpolate to common time points
     4. Compute per-field relative L2 error for all 7 non-time columns
     5. Report per-station, per-field errors
     6. Exit nonzero if any field exceeds tolerance
   ```

   This is a lightweight version of the Phase 0 `regression_check.py` focused on comparing two simulation outputs rather than against golden data. Key difference: it auto-detects file prefixes by scanning directory contents, handles different time grids via interpolation, and accepts different file naming prefixes.

   Required columns: `slip_strike, slip_dip, log10_V_strike, log10_V_dip, tau_strike, tau_dip, log10_state`

   Output format:
   ```
   Station: fltst_strk+00dp+00
     slip_strike:    L2 = 1.23e-15  PASS
     slip_dip:       L2 = 4.56e-16  PASS
     ...
   Overall: PASS (max L2 = 1.23e-15 < tolerance 1e-12)
   ```

2. **`tests/verification/bp5_driver_equivalence.cpp`**

   A C++ program that constructs BOTH operator paths and compares outputs directly in memory (no file I/O needed for the comparison). This is the strongest equivalence test:

   **IMPORTANT:** Both paths must use the same solver type for a valid comparison. The TOML default is `"mumps"` → `SolverType::MUMPS`. The legacy path's default is `SolverType::MUMPS_BLR`. Explicitly set both to `SolverType::MUMPS` to avoid false differences from solver numerics.

   ```cpp
   // Both paths use MUMPS (not MUMPS_BLR) — matches TOML default
   const SolverType solver = SolverType::MUMPS;
   const DGMethod method = DGMethod::IP;

   // OLD path: legacy constructor
   ElasticityDomainOperator<ParMesh> domain_old(
      pmesh, order, params.lambda(), params.mu(),
      params.Vp, params.Wf, params.lf,
      method, solver, BCMode::FarField);

   // NEW path: ConstitutiveModel + BoundaryConfig
   LinearElastic material(params.lambda(), params.mu());
   BoundaryConfig bdr_config;
   bdr_config.dirichlet_attrs = {5};
   bdr_config.natural_attrs = {1};
   bdr_config.fault_attr = 3;
   bdr_config.default_dirichlet_func = MakeBP5DirichletFunc(params.Vp);
   DomainConfig domain_config;
   domain_config.blr_tol = 1e-12;

   ElasticityDomainOperator<ParMesh> domain_new(
      pmesh, order, material, params.Vp, params.Wf, params.lf,
      bdr_config, method, solver, domain_config);

   // Compare: stiffness matrix, RHS at t=0, solution at t=0
   // Then run N steps with both, compare displacement and traction
   ```

   Steps to compare:
   - After construction: `domain_old.GetStiffnessMatrix()` vs `domain_new.GetStiffnessMatrix()` → must be identical (same assembly code path)
   - At t=0 with zero slip: solve both → compare displacement vectors → max |diff| < 1e-15
   - At t=1yr with zero slip: solve both → compare displacement and traction → max |diff| < 1e-15
   - After 5 full SEAS steps (domain + fault + SEAS operator): compare state vector → max |diff| < 1e-14

   If the stiffness matrices and displacement at t=0 are bit-identical, the paths are numerically equivalent. Any difference in later steps would come from FP accumulation in the time stepper, which should be bounded by 1e-14.

#### Acceptance Criteria
- [ ] `compare_station_files.py` correctly detects matching station files across two directories
- [ ] `compare_station_files.py` exits 0 when comparing a directory against itself
- [ ] `compare_station_files.py` exits nonzero when comparing against perturbed data (slip * 1.001)
- [ ] `bp5_driver_equivalence` passes at 1/4/8 ranks: max displacement diff < 1e-15 at t=0
- [ ] `bp5_driver_equivalence` passes at 8 ranks: max state vector diff < 1e-14 after 5 steps

#### Dependencies
- Depends on: Phase 7b (driver exists), Phase 7c (TOML configs exist)
- Required by: Phase 7e

---

### 13.7 Phase 7e: Golden Data Regression Test

#### Goal
Run the new driver against golden reference data to confirm it produces correct BP5 results.

#### Files to Modify
- `Makefile` — add regression test target

#### Detailed Requirements

1. **Local regression test (8 ranks, 20 steps, inline mesh)**

   ```makefile
   driver-regression-local: seas_driver
       @mkdir -p bp5_driver_regression
       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 8 ./seas_driver \
           config/bp5_quick_check.toml \
           --override output.output_dir=bp5_driver_regression \
           --override output.output_prefix=bp5_driver_regress
       python3 scripts/compare_station_files.py \
           bp5_driver_regression \
           bp5/benchmark_data/golden_inline_8r_20step \
           --tolerance 0
   ```

   This compares new-driver output against the Phase 0 local golden. Tolerance 0 = byte-for-byte.

2. **Frontera regression test (1 rank, 50 steps, reference mesh)**

   ```bash
   mpirun -np 1 ./seas_driver config/bp5_verification_50step.toml \
       --override output.output_dir=bp5_driver_frontera_1r \
       --override output.output_prefix=bp5_driver_serial

   python3 scripts/compare_station_files.py \
       bp5_driver_frontera_1r \
       bp5/benchmark_data/golden_serial_1r_50step \
       --tolerance 1e-12
   ```

3. **Frontera regression test (400 ranks, 50 steps, reference mesh)**

   Same as above with `np 400` and comparing against `golden_parallel_400r_50step`.

4. **End-to-end sbatch script for Frontera**

   Create `scripts/sbatch_driver_regression.sh` that:
   - Runs new driver at 1 rank + 400 ranks
   - Runs old driver at 1 rank + 400 ranks (as control)
   - Compares all 4 outputs: new-serial vs golden-serial, new-parallel vs golden-parallel, new-serial vs old-serial, new-parallel vs old-parallel
   - Reports all L2 errors in a single summary

#### Acceptance Criteria
- [ ] Local 8-rank/20-step: L2 = 0 (byte-for-byte) against golden
- [ ] If local golden doesn't exist yet: generate it first with old driver, then compare
- [ ] Frontera 1-rank/50-step: L2 < 1e-12 against golden_serial_1r_50step
- [ ] Frontera 400-rank/50-step: L2 < 1e-12 against golden_parallel_400r_50step
- [ ] New driver vs old driver (same rank count): L2 < 1e-15 (near byte-for-byte)

#### Dependencies
- Depends on: Phase 7d (comparison tooling), Phase 0 golden data
- Required by: nothing (final verification)

---

### 13.8 Phase 7 Checklist (All Sub-phases Combined)

**Phase 7a: Bridge & Enum Mapping**
- [ ] 7a-i. Create `config/seas_config_bridge.hpp` with `BuildBP5Params`, `ParseDGMethod`, `ParseSolverType`, `BuildBoundaryConfig`, `BuildDomainConfig`
- [ ] 7a-ii. Write `tests/unit/test_seas_config_bridge.cpp` (~200 LOC): round-trip, non-preset values, all enum variants, BoundaryConfig attrs + DirichletFunc verification, DomainConfig field mapping (see Phase 7a unit test spec for full list)
- [ ] 7a-iii. Add `test-config-bridge` target to Makefile, include in `make test`
- [ ] 7a-iv. `test_seas_config_bridge` compiles and all tests pass

**Phase 7b: Wire Simulation Pipeline**
- [ ] 7b-i. Add includes for all simulation components to `seas_driver.cpp`
- [ ] 7b-ii. Implement Stage 2: mesh loading from TOML file path with scale transform
- [ ] 7b-iii. Implement Stage 3: construct `LinearElastic` + `ElasticityDomainOperator` via NEW constructor
- [ ] 7b-iv. Implement Stage 4: construct `FaultGeometry`, `DieterichRuinaFriction`, `AgingLawPsi`, `RateStateFaultOperator`
- [ ] 7b-v. Implement Stage 5: construct `PBP5SEASOp`, set initial condition
- [ ] 7b-vi. Implement Stage 6: station list, `ParallelBP5BenchmarkOutput`, `ProbeOutput`
- [ ] 7b-vii. Implement Stage 7: `DormandPrinceRK45` with CFL-aware initial dt
- [ ] 7b-viii. Implement Stage 8: time loop with earthquake detection, bench_out.Write, checkpoint
- [ ] 7b-ix. Implement Stage 9: final output, summary statistics
- [ ] 7b-x. Add `__inline__` mesh sentinel support (calls `CreateBP5InlineMesh()`)
- [ ] 7b-xi. Compile test: `make seas_driver` succeeds with zero warnings
- [ ] 7b-xii. Smoke test: `mpirun -np 1 ./seas_driver config/bp5_example.toml --override time.max_steps=2` produces output

**Phase 7c: Verification Configs & Build Targets**
- [ ] 7c-i. Create `config/bp5_verification_50step.toml`
- [ ] 7c-ii. Create `config/bp5_quick_check.toml`
- [ ] 7c-iii. Add Makefile targets: `driver-quick-check`, `driver-equivalence`, `driver-regression-local`
- [ ] 7c-iv. `make driver-quick-check` runs without crash (8 ranks, 20 steps, inline mesh)

**Phase 7d: Equivalence Testing**
- [ ] 7d-i. Write `scripts/compare_station_files.py` with auto-prefix detection, per-field L2, tolerance gating
- [ ] 7d-ii. Self-test: `compare_station_files.py DIR DIR --tolerance 0` → exit 0 (self-comparison)
- [ ] 7d-iii. Self-test: perturbed comparison → exit nonzero
- [ ] 7d-iv. Write `tests/verification/bp5_driver_equivalence.cpp`: side-by-side old-path vs new-path
- [ ] 7d-v. `bp5_driver_equivalence` at 4 ranks: displacement diff < 1e-15 at t=0
- [ ] 7d-vi. `bp5_driver_equivalence` at 8 ranks: state vector diff < 1e-14 after 5 steps
- [ ] 7d-vii. Run `make driver-equivalence`: station-file comparison passes at tolerance 1e-12

**Phase 7e: Golden Data Regression**
- [ ] 7e-i. Generate local 8-rank/20-step golden with OLD driver (if not already done from Phase 0)
- [ ] 7e-ii. `make driver-regression-local`: new driver output matches golden at tolerance 0
- [ ] 7e-iii. Create `scripts/sbatch_driver_regression.sh` for Frontera
- [ ] 7e-iv. Frontera serial (1 rank, 50 steps): L2 < 1e-12 vs golden_serial_1r_50step
- [ ] 7e-v. Frontera parallel (400 ranks, 50 steps): L2 < 1e-12 vs golden_parallel_400r_50step
- [ ] 7e-vi. Record actual L2 in this document
- [ ] 7e-vii. Run full `make test` suite — all existing tests still pass
- [ ] 7e-viii. Merge to `system_update`

### 13.9 Summary Table Update

| Phase | Branch | Risk | Regression | Key Deliverable | New Tests |
|-------|--------|------|-----------|-----------------|-----------|
| 7a | `phase7-driver-verification` | Low | N/A (utility only) | `seas_config_bridge.hpp` | `test_seas_config_bridge.cpp` (~200 LOC) |
| 7b | same | **High** | First run | Full simulation in `seas_driver.cpp` | Smoke test |
| 7c | same | Low | N/A | TOML configs, Makefile targets | Quick-check target |
| 7d | same | Medium | L2 < 1e-12 | Equivalence proof: old == new | `bp5_driver_equivalence.cpp`, `compare_station_files.py` |
| 7e | same | Low | L2 < 1e-12 | Golden regression confirmed | Frontera sbatch |

### 13.10 Rollback Cost

If Phase 7 fails: delete `seas_config_bridge.hpp`, revert `seas_driver.cpp` to the config-printer, remove TOML configs and new Makefile targets. Old driver and all existing tests are untouched. **Rollback: 5 minutes with git.**

### 13.11 Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| New constructor produces different results than old | Low (same assembly code) | High | Phase 7d equivalence test catches it immediately |
| `BoundaryConfig` path classifies faces differently than `BCMode::FarField` | Medium | High | Phase 7d side-by-side test at t=0 catches it before time integration |
| `DomainConfig` defaults differ from post-construction setter values | Medium | Medium | Phase 7a bridge sets `blr_tol` explicitly; Phase 7d catches any drift |
| `BP5Params` bridge misses a field | Low | High | Round-trip unit test in 7a-ii catches it |
| Station file format mismatch | Low | Low | `compare_station_files.py` reports specific column/row differences |
| Inline mesh (from `CreateBP5InlineMesh()`) doesn't have the same boundary attributes as the Gmsh mesh | Medium | High | Verify `__inline__` mesh boundary attributes in Phase 7b before constructing BoundaryConfig |

---

## Appendix: Known Limitations

1. **Matrix assembly** — stored, not matrix-free. *Not addressed.*
2. **Direct solver default** — MUMPS for BP5. *Not addressed.*
3. **Single fault plane at Y=0** — *Phase 4 BoundaryConfig is partial step. Multi-fault is future work.*
4. **Quasi-dynamic only** — *Section 12.1 defines three modes (qd/dynamic/hybrid). Dynamic and hybrid are future work; config infrastructure ready in Phase 5.*
5. **DG noise on unstructured meshes** — *SIPG as new primary may help.*
6. **Linear elasticity only** — *Phase 3 establishes abstraction. Nonlinear is Section 12.3.*
7. **Single-physics** — *Section 12.4 (noted, not urgent).*
8. **R&S friction only** — *Section 12.2 (future work).*
