# Implementation Plan: Heterogeneous (CVM-H) Material Support for SEAS

## Overview

Wire the existing CVM-H sidecar pipeline (DataField3D / FieldCoefficient / sidecar HDF5 in `data_projected/velocity_safs.h5`) into the two consumer paths that currently run with constant material — the quasi-dynamic ElasticityDomainOperator (used by SAFS / BP5 drivers) and the dynamic-rupture WaveOperator (used by tpv102 / tpv104 / tpv205 drivers). The new path is an **additional option** selected by a runtime CLI flag; the existing scalar / `ConstantCoefficient` code path remains the default and must continue to produce byte-identical output for every existing fixture.

The plan deliberately adopts the **`FieldCoefficient` direct path** rather than the `MaterialField::Mode::GridFunction` path sketched in `data_projection_feature_plan_v2.md` Phase 5. The reason was established empirically in the debugging session preceding this plan (recorded in `comparison_old_vs_new_baseline_z1km.png` and the residual-table dumps from `plot_old_vs_new_baseline.py`): P1 linear FE interpolation between corner DOFs cannot represent the sharp basin layering in the upper 3 km (Vs varies from ≈0.3 to ≈3 km/s over <1 km of depth), producing per-tet outlier errors of 1 – 7 km/s in the projected GridFunction.

**Two paths, two granularities (R-006).**  The two consumer operators benefit DIFFERENTLY from `Mode::Coefficient`:

- **Quasi-dynamic `ElasticityDomainOperator`.**  The DG integrators (`DGElasticityBR2Integrator`, `DGElasticityIPPenaltyIntegrator`, etc.) evaluate the `Coefficient` directly at EACH quadrature point inside the bilinear-form assembly.  So `Mode::Coefficient` gives full per-qp sidecar-accurate material here: `λ(qp) = lambda_coef->Eval(T, ip_qp)` = sidecar trilinear at `(x_qp, y_qp, z_qp)`.  The FE-interpolation error vanishes by construction.

- **Dynamic-rupture `WaveOperator`.**  `GodunovFlux` precomputes 9×9 Jacobians from a SINGLE `(λ, μ, ρ)` triple, so the wave path stores per-ELEMENT flat material in this plan — see Phase 3 `GodunovFluxPool` (Detailed Req. 1).  The Coefficient is evaluated ONCE per element at the CENTROID; the per-qp sidecar accuracy is NOT achieved on this path.  This is a per-element upgrade over the previous per-DOMAIN constant; it is NOT the per-qp upgrade the quasi-dynamic path enjoys.  A follow-up plan can promote the wave path to per-qp flux Jacobians at the cost of substantially more memory/compute.

A new `MaterialField::Mode::Coefficient` is added alongside (not replacing) the existing `Mode::Constant` and `Mode::GridFunction`. `Mode::Coefficient` is the recommended runtime mode for new SAFS / CVM-H drivers; `Mode::GridFunction` remains supported for backwards compatibility with code paths that already require a GridFunction handle (e.g. ParaView snapshots of the projected fields).

## Constraints

### Interface constraints — what cannot change
- **Existing scalar constructors must remain bit-equivalent.** `WaveOperator<MeshType>::WaveOperator(MeshType&, int, real_t λ, real_t μ, real_t ρ, const BoundaryConfig&)` and the two existing `ElasticityDomainOperator` constructors (the `ConstitutiveModel`-based one and the legacy scalar `(λ, μ)` one) must accept the same arguments and produce the same outputs they do today.
- **All existing unit tests, including TPV102 / TPV104 / TPV205 / BP5 regression tests, must continue to pass without modification.** Drivers that do not opt into the new flag must not link any sidecar I/O.
- **DG integrator signatures cannot change.** `DGElasticityBR2Integrator`, `DGElasticityIPPenaltyIntegrator`, and `DGElasticityIPCombinedIntegrator` already take `Coefficient&` (verified in `integrator/dg_elasticity_*_integrator.hpp`), so they polymorphically accept any `mfem::Coefficient` subclass. No integrator-side change is required.
- **`FieldCoefficient` and `DataField3D` are NOT to be modified.** The sidecar pipeline is finished and its contract (`DataField3D::Evaluate` is the trilinear interpolant in canonical UTM 11 N coordinates, aborts on out-of-bbox) is the source of truth.
- **`MFEM_ABORT` is the only mechanism for sidecar-related runtime failures.** No silent fallback to scalar material if the sidecar load fails.

### Dependency constraints
- The sidecar HDF5 file must already exist on disk (built by `build_velocity_cvmh.py`). No on-the-fly construction.
- The pre-flight mesh-bbox containment check (`FieldProjector::Project` calls `DataField3D::ContainsBBox` at startup) must run BEFORE any quadrature-point evaluation, regardless of whether `Mode::Coefficient` or `Mode::GridFunction` is selected.
- Density (`density` field in the sidecar) is required in addition to Vp and Vs. The existing sidecar already contains it (`fields/density`).

### Convention constraints
- The CLI flag is named `--material-mode` with three accepted values: `constant` (default), `sidecar-coefficient` (= `Mode::Coefficient`), `sidecar-gridfunction` (= `Mode::GridFunction`). A second flag `--sidecar PATH` supplies the HDF5 path when the mode is one of the two sidecar variants; aborts at startup if mode is sidecar but `--sidecar` is empty.
- New struct fields are placed in the existing `MaterialField` struct in `dynamic/heterogeneous_material.hpp`. No parallel struct.
- Test files follow the project convention: `tests/unit/test_<name>.cpp` driven by `make test-<name>` in `miniapps/seas/Makefile`. The test for `MaterialField::Mode::Coefficient` extends the existing `tests/unit/test_heterogeneous_material.cpp`; new tests for the operator integrations live in new files.
- `--no-touch` zones from `CLAUDE.md`: do NOT modify `bp5/`, `bp1/`, `bp2/`, `domain/` files except `elasticity_operator.hpp` and its `*.inl` companions, and only the additive constructor / additive coefficient-owning members described below. Tandem-derived friction code (`friction/dieterich_ruina.hpp`) is untouched.

### Numerical constraints
- `Mode::Coefficient` material accuracy at each quadrature point must equal sidecar trilinear accuracy. The plan's primary numerical acceptance test (T-COEFF-QP-EQUALITY) requires `|GF(qp) − sidecar_trilinear(qp)| < 1.0 m/s` at every quadrature point of a smoke-test mesh. (Strict machine epsilon is not realistic because `mfem::ElementTransformation::Transform` is FP-arithmetic; 1 m/s out of ≈3000 m/s is well within roundoff.)
- For `Mode::Constant`, every assembled matrix entry must remain bit-identical to the pre-change reference. This is verified by an L∞ comparison against a saved reference (T-CONST-REGRESSION).
- Heterogeneous CFL for the dynamic-rupture path must use **per-element maximum c_p** computed at sidecar quadrature points, then reduced with `MPI_Allreduce(MIN)` over `dt_local = h_e / cp_max_e`. The TPV102 constant-material baseline `Δt` must reproduce to within FP roundoff (T-CFL-CONST-PARITY).

## Phase 1: Extend `MaterialField` with `Mode::Coefficient`

### Goal
Add a third mode to `MaterialField` (in `miniapps/seas/dynamic/heterogeneous_material.hpp`) that holds three `mfem::Coefficient` pointers (one per λ, μ, ρ) and evaluates them at a supplied `(ElementTransformation, IntegrationPoint)`. The struct remains usable in both quasi-dynamic and dynamic-rupture code paths.

### Files to Create
- `miniapps/seas/dynamic/heterogeneous_material.cpp` — non-trivial helpers that cannot live in the header (Coefficient evaluation needs `ElementTransformation` mutation, so the previous header-only design is incomplete).
- `miniapps/seas/io/material_coefficients.hpp` — three small `mfem::Coefficient` subclasses `LambdaFromSidecar`, `MuFromSidecar`, `RhoFromSidecar`. Each holds references to `DataField3D` objects (`vp_field`, `vs_field`, `rho_field`) and computes:
  ```
  ρ(x,y,z) = rho_field.Evaluate(x,y,z)
  μ(x,y,z) = ρ · Vs(x,y,z)²
  λ(x,y,z) = ρ · Vp(x,y,z)² − 2·μ
  ```
  Implementing them as derived `Coefficient` classes (not lambdas via `FunctionCoefficient`) keeps the `Eval` signature stable and means the integrator's `dynamic_cast`-free hot path is preserved.
- `miniapps/seas/io/material_coefficients.cpp`

### Files to Modify
- `miniapps/seas/dynamic/heterogeneous_material.hpp`:
  - Add `Mode::Coefficient = 2` to the `Mode` enum.
  - Add three `mfem::Coefficient*` members (non-owning):
    ```cpp
    mfem::Coefficient* lambda_coef = nullptr;
    mfem::Coefficient* mu_coef     = nullptr;
    mfem::Coefficient* rho_coef    = nullptr;
    ```
    Non-owning because Coefficients depend on `DataField3D` objects that the driver allocates and outlives the operator.
  - Add a factory `MaterialField::MakeCoefficient(mfem::Coefficient* λ, mfem::Coefficient* μ, mfem::Coefficient* ρ)`.
  - Extend the inline accessor with an `EvalAt(int elem, ElementTransformation& T, const IntegrationPoint& ip, real_t& λ_out, real_t& μ_out, real_t& ρ_out) const`. This is the version called from the dynamic-rupture hot path (`Mult`, ADER CK recursion) because integrators in `dynamic/` evaluate material at qpoints inside flux dispatch. The existing `At(int elem, int dof, …)` accessor stays for `Mode::Constant` and `Mode::GridFunction`; in `Mode::Coefficient` it MUST abort (Coefficient mode requires the qpoint context, not a nodal index).
  - Extend `MaxCpInElement(int elem)` to accept an optional `ElementTransformation*` so the Coefficient mode can sample c_p at the element's quadrature points; for `Constant` / `GridFunction` the argument is ignored.
  - **The existing header-only inline body of `MaxCpInElement` MUST BE DELETED** (R-004).  The new version calls `IntRules.Get` in `Mode::Coefficient` and is no longer suitable for inlining.  Phase 1 patch leaves only the declaration in the header (with the `(int elem, mfem::ElementTransformation* T = nullptr) const` signature) and adds a single definition to `dynamic/heterogeneous_material.cpp`.  Leaving the inline body in place AND adding a `.cpp` definition causes a linker ODR error.
  - Tighten the existing `MakeGridFunction` factory (R-006 round-3): add an `MFEM_VERIFY` that all three shared_ptr arguments are non-null (R-004 round-2; full body in Detailed Req. 7).  The factory remains header-inline (the addition is trivial).

### Detailed Requirements
1. **Header additions** (verbatim signatures):
   ```cpp
   enum class Mode : int { Constant = 0, GridFunction = 1, Coefficient = 2 };

   // Non-owning. Caller (driver) owns these and outlives the struct.
   mfem::Coefficient* lambda_coef = nullptr;
   mfem::Coefficient* mu_coef     = nullptr;
   mfem::Coefficient* rho_coef    = nullptr;

   static MaterialField MakeCoefficient(mfem::Coefficient* lambda,
                                        mfem::Coefficient* mu,
                                        mfem::Coefficient* rho);

   // Hot-path accessor for dynamic-rupture (qpoint context).
   inline void EvalAt(int elem,
                      mfem::ElementTransformation& T,
                      const mfem::IntegrationPoint& ip,
                      real_t& lambda_out, real_t& mu_out, real_t& rho_out) const;

   // CFL helper — Coefficient mode needs T to sample qpoints.
   real_t MaxCpInElement(int elem,
                         mfem::ElementTransformation* T = nullptr) const;
   ```

2. **`EvalAt` semantics** (header-inline for the hot path):
   - `Mode::Constant` → branch-free return of `(lambda_const, mu_const, rho_const)`. Same as `At`.
   - `Mode::GridFunction` → MFEM_ABORT (use `At(elem, dof, …)` instead; Coefficient is the only mode that meaningfully consumes a qpoint).
   - `Mode::Coefficient` → `lambda_out = lambda_coef->Eval(T, ip);` etc. **All three pointers must be non-null; abort with a precise message if any is null.**

3. **`MaterialField::At(int elem, int dof, …)`** in `Mode::Coefficient` MUST abort with `MFEM_ABORT("MaterialField::At called in Mode::Coefficient — use EvalAt(elem, T, ip, …)")`. Tests pin this contract.

4. **`MaxCpInElement` Coefficient path** (in `.cpp`, since it walks IntRules):
   - Build the element's natural-order quadrature rule `IntRules.Get(geom, 2 * default_qorder)` with `default_qorder = 2` (R-005 round-2 — `mfem::Coefficient` has no `GetFESpace()` accessor and querying it on every subclass would add dispatch complexity for negligible accuracy improvement; `c_p` varies slowly in space, so 4-qpoint tet sampling is more than enough).  If a future caller needs higher accuracy, add an `int qorder = 2` parameter to `MaxCpInElement` and route it through.
   - At each qpoint, call `EvalAt(elem, T, ip, …)` and compute `cp = sqrt((λ + 2μ) / ρ)`.
   - Return `max(cp)` across all qpoints in the element.

5. **`LambdaFromSidecar`, `MuFromSidecar`, `RhoFromSidecar`** classes — each subclasses `mfem::Coefficient`, holds references to `DataField3D vp_field, vs_field, rho_field` (Coefficient classes only need the fields they actually depend on; `RhoFromSidecar` only stores `rho_field`, etc.). Each `Eval(T, ip)`:
   ```cpp
   real_t Eval(ElementTransformation& T,
               const IntegrationPoint& ip) override
   {
       Vector x(3);
       T.Transform(ip, x);
       const real_t rho = rho_field_.Evaluate(x[0], x[1], x[2]);
       const real_t vs  = vs_field_.Evaluate(x[0], x[1], x[2]);   // MuFromSidecar
       const real_t vp  = vp_field_.Evaluate(x[0], x[1], x[2]);   // LambdaFromSidecar
       // ...
       return /* λ, μ, or ρ formula */;
   }
   ```
   These classes are the only place where Vp/Vs/density → λ/μ/ρ is computed; the rest of the operator code touches only λ/μ/ρ.

6. **`MakeCoefficient` factory**:
   ```cpp
   MaterialField MaterialField::MakeCoefficient(
       mfem::Coefficient* lambda, mfem::Coefficient* mu, mfem::Coefficient* rho)
   {
       MFEM_VERIFY(lambda && mu && rho,
                   "MakeCoefficient: all three Coefficient pointers must be non-null");
       MaterialField m;
       m.mode = Mode::Coefficient;
       m.lambda_coef = lambda;
       m.mu_coef     = mu;
       m.rho_coef    = rho;
       return m;
   }
   ```

7. **`MakeGridFunction` factory hardening (R-004 round-2).**  The existing factory in `heterogeneous_material.hpp` accepts three shared_ptrs but does not validate them.  Tighten it to abort on any null:
   ```cpp
   MaterialField MaterialField::MakeGridFunction(
       std::shared_ptr<mfem::ParGridFunction> rho_gf,
       std::shared_ptr<mfem::ParGridFunction> lambda_gf,
       std::shared_ptr<mfem::ParGridFunction> mu_gf)
   {
       MFEM_VERIFY(rho_gf && lambda_gf && mu_gf,
                   "MakeGridFunction: all three ParGridFunction shared_ptrs "
                   "must be non-null.");
       MaterialField m;
       m.mode      = Mode::GridFunction;
       m.rho_gf    = std::move(rho_gf);
       m.lambda_gf = std::move(lambda_gf);
       m.mu_gf     = std::move(mu_gf);
       return m;
   }
   ```
   This is a strict tightening — every existing caller passes the result of `FieldProjector::ProjectVelocity`, which never returns nulls — so existing tests pass unchanged.

### Interfaces
- `material_coefficients.hpp` exposes:
  ```cpp
  namespace mfem::seas {
  class RhoFromSidecar : public mfem::Coefficient { ... };
  class MuFromSidecar  : public mfem::Coefficient { ... };
  class LambdaFromSidecar : public mfem::Coefficient { ... };
  } // namespace mfem::seas
  ```
- `heterogeneous_material.hpp` exposes `MaterialField::Mode::Coefficient`, `MakeCoefficient`, `EvalAt`.
- No other consumer code needs the new headers in this phase.

### Edge Cases to Handle
- **All three Coefficient pointers must be non-null in `Mode::Coefficient`.** `MakeCoefficient` aborts if any null. `EvalAt` re-asserts in debug builds via `MFEM_ASSERT`.
- **`At(elem, dof, …)` called in `Mode::Coefficient`** → MFEM_ABORT (not silent fallback). Documents that drivers picking `Mode::Coefficient` must route all material lookups through `EvalAt`.
- **`MaxCpInElement(elem, nullptr)` in `Mode::Coefficient`** → MFEM_ABORT (requires the element transformation handle).
- **Vp ≤ √2·Vs at any sidecar voxel** → λ < 0. This is a sidecar correctness issue, not an operator one; the sidecar pre-flight check in `build_velocity_cvmh.py` (`--vp-min-mps`/`--vs-max-mps`) is responsible for catching it. The Coefficient `LambdaFromSidecar::Eval` does NOT clamp.

### Acceptance Criteria
- [ ] `make test-heterogeneous-material` passes, including new tests T-5-4, T-5-5, T-5-6 (see Testing Strategy).
- [ ] Existing tests T-5-1, T-5-2, T-5-3 (constant + GridFunction modes) still pass byte-identically.
- [ ] `material_coefficients.cpp` builds with `-Werror` and is exercised by at least one test that constructs each of the three Coefficient subclasses from a synthetic `DataField3D`.

### Dependencies
- Depends on: nothing (extends an existing struct; no new infrastructure).
- Required by: Phase 2, Phase 3.

---

## Phase 2: ElasticityDomainOperator — accept `MaterialField` with `Mode::Coefficient`

### Goal
Add a new constructor `ElasticityDomainOperator(mesh, order, const MaterialField&, …)` that swaps in user-supplied λ, μ Coefficients while preserving every other behaviour of the existing class. When `material.mode == Mode::Constant` the resulting matrix is bit-identical to today's `ConstantCoefficient` path.

### Files to Modify
- `miniapps/seas/domain/elasticity_operator.hpp`:
  - Add `#include "../dynamic/heterogeneous_material.hpp"` at the top of the header so `MaterialField` is in scope (R-010 round-2).
  - Replace the two `ConstantCoefficient` member fields
    ```cpp
    mutable ConstantCoefficient lambda_coeff_, mu_coeff_;
    ```
    with a **non-owning view pointer** plus an **owned ConstantCoefficient
    used only in the constant path** (R-007 — avoid the `PassthroughCoefficient`
    indirection and its extra virtual dispatch per qpoint):
    ```cpp
    // Owned only in Mode::Constant; default-constructed in the other modes.
    mutable std::unique_ptr<mfem::ConstantCoefficient> owned_const_lambda_;
    mutable std::unique_ptr<mfem::ConstantCoefficient> owned_const_mu_;
    // Non-owning view dereferenced by every integrator.  In Mode::Constant
    // points at owned_const_lambda_.get(); in Mode::Coefficient points
    // directly at the caller's MaterialField coefficients; in
    // Mode::GridFunction points at internally-built GridFunctionCoefficient
    // instances (which need separate ownership — see below).
    mutable mfem::Coefficient* lambda_coef_view_ = nullptr;
    mutable mfem::Coefficient* mu_coef_view_     = nullptr;
    // Mode::GridFunction also needs owned GridFunctionCoefficients.
    mutable std::unique_ptr<mfem::GridFunctionCoefficient> owned_gf_lambda_;
    mutable std::unique_ptr<mfem::GridFunctionCoefficient> owned_gf_mu_;
    ```
    The two existing constructors populate `owned_const_lambda_` /
    `owned_const_mu_` and set the view pointers to their `.get()` values,
    preserving the existing assembly path with EXACTLY ONE virtual
    dispatch per evaluation.
  - Add a new constructor:
    ```cpp
    ElasticityDomainOperator(MeshType& mesh, int order,
                              const MaterialField& material,
                              real_t Vp, real_t Wf, real_t lf,
                              const BoundaryConfig& bdr_config,
                              DGMethod method = DGMethod::BR2,
                              SolverType solver_type = SolverType::MUMPS_BLR,
                              const DomainConfig& config = {});
    ```
    Body branches on `material.mode`:
    - `Mode::Constant`: `owned_const_lambda_ = std::make_unique<ConstantCoefficient>(material.lambda_const);` ditto μ.  `lambda_coef_view_ = owned_const_lambda_.get();` ditto μ.  Set `lambda_val_ = material.lambda_const;` `mu_val_ = material.mu_const;` for the legacy scalar accessors.
    - `Mode::Coefficient`: MFEM_VERIFY that all three `MaterialField` coefficient pointers are non-null; set `lambda_coef_view_ = material.lambda_coef;` ditto μ.  **No `Eval` of element-0-centroid** (R-002): set `lambda_val_ = std::numeric_limits<real_t>::quiet_NaN();` and `mu_val_ = std::numeric_limits<real_t>::quiet_NaN();` as sentinels.  Document on `GetLambda()` / `GetShearModulus()` that they return NaN in `Mode::Coefficient`; any numerical use NaN-propagates rather than silently producing a rank-dependent wrong value.
    - `Mode::GridFunction`: MFEM_VERIFY that `material.rho_gf && material.lambda_gf && material.mu_gf` are all non-null shared_ptrs (R-004 round-2 — a caller building a `MaterialField` by hand instead of via the `MakeGridFunction` factory can leave one slot null; without this check the subsequent `.get()` is `nullptr` and `GridFunctionCoefficient(nullptr)` is UB).  Then `owned_gf_lambda_ = std::make_unique<GridFunctionCoefficient>(material.lambda_gf.get());` ditto μ.  `lambda_coef_view_ = owned_gf_lambda_.get();` ditto μ.  `lambda_val_` / `mu_val_` similarly set to NaN.
  - All assembly sites that today do `new ElasticityIntegrator(lambda_coeff_, mu_coeff_)` change to `new ElasticityIntegrator(*lambda_coef_view_, *mu_coef_view_)`. This is a mechanical search-and-replace; every existing call site (10 hits in `elasticity_operator_assembly.inl`, 2 in `elasticity_operator_traction.inl`, 1 in `elasticity_operator_debug.inl`, 1 in `elasticity_operator_verify.inl`) is updated identically.
- `miniapps/seas/domain/elasticity_operator_assembly.inl`, `elasticity_operator_traction.inl`, `elasticity_operator_debug.inl`, `elasticity_operator_verify.inl` — pure mechanical `lambda_coeff_` → `*lambda_coef_view_` replacement; same for μ.

### Detailed Requirements
1. **`InitOperator()` is unchanged.** All differences are absorbed by which concrete `mfem::Coefficient` lives behind the owning pointers; integrators are agnostic.

2. **No `PassthroughCoefficient` wrapper** (R-007).  The non-owning-view-pointer pattern in §Files to Modify gives each integrator a direct reference to the active `mfem::Coefficient` with exactly ONE virtual dispatch per evaluation.

3. **Two existing constructors stay verbatim** except the body's `lambda_coeff_, mu_coeff_` initializer-list lines are replaced by `owned_const_lambda_` / `owned_const_mu_` construction inside the body:
   ```cpp
   // was: lambda_coeff_(lambda), mu_coeff_(mu),
   // now (inside body, after initializer list):
   owned_const_lambda_ = std::make_unique<ConstantCoefficient>(lambda);
   owned_const_mu_     = std::make_unique<ConstantCoefficient>(mu);
   lambda_coef_view_   = owned_const_lambda_.get();
   mu_coef_view_       = owned_const_mu_.get();
   ```

4. **`ComputeTraction` / `ComputeTractionComponents` / `ComputeTractionDiagnostics`** — these use `lambda_coeff_`, `mu_coeff_` inside `DGElasticityBR2BoundaryIntegrator` and similar.  Same `*lambda_coef_view_` / `*mu_coef_view_` rewrite.  No semantic change.

5. **`GetLambda()` accessor (R-002).**  `real_t GetLambda() const { return lambda_val_; }` and `GetShearModulus() override { return mu_val_; }` continue to return cached scalars.  In `Mode::Constant` they return the constructor-supplied `lambda` / `mu`.  In `Mode::Coefficient` and `Mode::GridFunction` they return `std::numeric_limits<real_t>::quiet_NaN()` as a sentinel — the documentation MUST say "scalar accessor; returns NaN when material is heterogeneous (use the `MaterialField` and `EvalAt(elem, T, ip, ...)` instead)."  Returning NaN makes any numerical misuse propagate as NaN rather than silently use an arbitrary rank-dependent sample.

6. **Stiffness matrix caching** (`stiffness_assembled_`, `cached_a_`, `cached_Ah_`) — UNCHANGED. The stiffness for `Mode::Coefficient` is heterogeneous but still constant-in-time; it is assembled once at construction and reused for every quasi-dynamic time step, exactly like the constant case.

7. **Boundary face integrators** — same mechanical rewrite. The TractionIntegrator path (lines 663 and 1151 in `elasticity_operator_traction.inl`) is the most subtle: `traction_lambda_coeff_, traction_mu_coeff_` pattern stays, just plumbed from `*lambda_coef_view_`.

### Interfaces
- New constructor:
  ```cpp
  ElasticityDomainOperator(MeshType&, int order, const MaterialField&,
                            real_t Vp, real_t Wf, real_t lf,
                            const BoundaryConfig&, DGMethod, SolverType,
                            const DomainConfig& = {});
  ```
- No new public methods.

### Edge Cases to Handle
- **`MaterialField::Mode::GridFunction` passed.** This is supported but emits a one-line `mfem::out` warning at construction: "Note: ElasticityDomainOperator with Mode::GridFunction uses GridFunctionCoefficient at qpoints. For sharp basin layering, prefer Mode::Coefficient (FieldCoefficient direct path)."
- **Constructor receives `MaterialField` with null Coefficient pointers in Mode::Coefficient** → MFEM_ABORT.
- **Material outlives operator?** The operator stores Coefficient pointers but does not own them in `Mode::Coefficient`. Document in the ctor docstring: "Caller must keep the `MaterialField` and its referenced Coefficient objects alive at least until `ElasticityDomainOperator` is destroyed."

### Acceptance Criteria
- [ ] All existing BP5 tests (`make test-bp5-integration`, `make test-bp5-smoke`) pass byte-identically.
- [ ] New test T-EOP-COEFF (see Testing Strategy) constructs `ElasticityDomainOperator` with `Mode::Coefficient` over a small fixture and confirms the assembled stiffness matrix entries match a direct symbolic computation to within FP roundoff.
- [ ] `git diff` of `elasticity_operator_assembly.inl` shows ONLY the mechanical `lambda_coeff_` → `*lambda_coeff_owned_` rewrite, NO semantic changes. (Reviewer-enforced.)

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 4 (driver wiring for quasi-dynamic).

---

## Phase 3: WaveOperator — accept `MaterialField` with `Mode::Coefficient`

### Goal
Add a new `WaveOperator<MeshType>` constructor accepting `const MaterialField&`. In `Mode::Constant` the existing `GodunovFlux` precomputation and per-element CFL are preserved bit-identically. In `Mode::Coefficient` the flux is evaluated **per element** (per-element 9×9 Jacobians stored), and the CFL helper computes `c_p` per element.

This phase is structurally larger than Phase 2 because `GodunovFlux` was designed around a SINGLE constant `(λ, μ, ρ)` triple — its `Ax_`, `Ax_plus_`, `Ax_minus_` matrices (and the three cached `ref_star_[3]` matrices used by the ADER recursion) are precomputed once at construction.  Each `GodunovFlux` instance therefore owns **six** 9×9 `DenseMatrix` objects, each carrying an 81-double heap payload, so the realistic per-instance heap cost is roughly **4 KB** (= 6 × 81 × 8 bytes + small headers — R-001 round-3; the earlier "9×9 doubles = 648 bytes; 1 M elements = 648 MB" estimate was off by 6× because it counted only ONE matrix per instance).  Going heterogeneous means either (a) caching one `GodunovFlux` per element, or (b) computing the 9×9 Jacobians per element on the fly. **Option (a)** is chosen because the uniqueness map (Detailed Req. 2) typically dedups 5–10× across smoothly-varying layers, bringing realistic memory down to **400 MB – 2 GB** for the SAFS 3.2 M-tet fixture, and the lookup is then a flat array index, preserving the hot-path access pattern.  If Option (a) overshoots the budget on real fixtures, fall back to Option (b) — JIT Jacobian build — or strip the cached `ref_star_[3]` triple from heterogeneous-mode pools (it can be re-derived from `Ax_`/`Ay_`/`Az_` on demand at modest CPU cost).

### Files to Create
- `miniapps/seas/dynamic/godunov_flux_pool.hpp`, `godunov_flux_pool.cpp` — a `GodunovFluxPool` container that owns one `GodunovFlux` per element (or shared one per unique material triple, see Detailed Req. 6). Exposes:
  ```cpp
  // R-001 (round-2): template parameter parallels WaveOperator<MeshType>.
  // The MPI face-neighbour exchange (Det. Req. 7a) is only valid when
  // MeshType == mfem::ParMesh.  For the serial specialisation
  // `GodunovFluxPool<mfem::Mesh>` the exchange is skipped and
  // `AtNbr(...)` aborts (serial meshes have no face neighbours).
  //
  // Mesh argument is non-const (R-002 round-3) because Detailed Req. 7a
  // builds a temporary `mfem::ParFiniteElementSpace` over the supplied
  // ParMesh — that ctor requires a non-const `ParMesh*`.  Adding a
  // const_cast at the call site would reintroduce the round-1 R-013
  // const-correctness hazard.
  template <typename MeshType = mfem::Mesh>
  class GodunovFluxPool
  {
  public:
      GodunovFluxPool(const MaterialField& material,
                      MeshType& mesh,
                      int dedup_sig_figs = 6);
      const GodunovFlux& At(int elem_local) const;            // hot-path lookup (local element)
      const GodunovFlux& AtNbr(int nbr_local_idx) const;       // hot-path lookup (face-neighbour, MPI only)
      real_t MaxCp() const;                                    // global maximum, MPI-reduced (parallel only; serial returns local max)
      real_t CpForElement(int elem_local) const;               // per-element c_p
      int NumUniqueFluxes() const;                             // for memory diagnostics
  };
  ```
  In `Mode::Constant`, the pool has a single internal `GodunovFlux`; `elems_to_flux_idx_` is sized `ne` and zero-initialised so every `At(elem)` lookup is well-defined.  The `quadrature_order` parameter previously listed here is REMOVED (R-012); the pool does not need a quadrature rule — `Mode::Coefficient` / `Mode::GridFunction` sample at the element CENTROID only (per Detailed Req. 1).  The new `dedup_sig_figs` parameter (default 6) controls the rounding-bucket size for the uniqueness map (R-011).

### Files to Modify
- `miniapps/seas/dynamic/wave_operator.hpp`:
  - Add `#include "heterogeneous_material.hpp"` and `#include "godunov_flux_pool.hpp"` at the top of the header (R-010).
  - Replace member `GodunovFlux flux_;` with `std::unique_ptr<GodunovFluxPool<MeshType>> flux_pool_;` (matches the WaveOperator's own `MeshType` template parameter; see R-001).  The existing `const GodunovFlux& GetFlux()` accessor changes to `const GodunovFlux& GetFlux(int elem = 0) const { return flux_pool_->At(elem); }` (with default `elem=0` preserving the existing accessor's signature for tests that fetch the "first" flux; this is documented as the `Mode::Constant` accessor).
  - Add a new constructor:
    ```cpp
    WaveOperator(MeshType& mesh, int order,
                 const MaterialField& material,
                 const BoundaryConfig& bc);
    ```
  - The existing `WaveOperator(MeshType&, int, real_t λ, real_t μ, real_t ρ, const BoundaryConfig&)` constructor calls the new one internally with `MaterialField::MakeConstant(λ, μ, ρ)`.
- `miniapps/seas/dynamic/wave_operator.inl`:
  - All references to `flux_.X()` become `flux_pool_->At(elem).X()`. In hot paths (`ComputeVolumeRHS`, `ComputeFaceFluxRHS`, `ComputeSharedFaceFluxRHS`, and the ADER counterparts), the element index is already in scope, so this is a mechanical rewrite. For face fluxes that need BOTH sides' flux Jacobians, both `flux_pool_->At(elem_plus)` and `flux_pool_->At(elem_minus)` are fetched.
  - `flux_.BuildJacobian(d, A_d)` for `d ∈ {0,1,2}` (lines 55-57 in the current ctor) becomes "for each element, store `Ax_e`, `Ay_e`, `Az_e`" inside the `GodunovFluxPool`. The cached `Ax_`, `Ay_`, `Az_` `DenseMatrix` members at WaveOperator scope are removed (they become per-element data living in the pool).
  - `flux_.GetCp()` calls (line 5001 in `ComputeMaxDt`) become per-element c_p with `flux_pool_->CpForElement(e)` inside the existing per-element loop, then reduced with `MPI_Allreduce(MIN, dt_local)`.
- `miniapps/seas/dynamic/godunov_flux.hpp`, `godunov_flux.cpp` — UNTOUCHED. The pool calls existing `GodunovFlux::GodunovFlux(λ, μ, ρ)` per element.

### Detailed Requirements
1. **`GodunovFluxPool` ctor** (in `.cpp`):
   - `Mode::Constant`: allocate a single `GodunovFlux(material.lambda_const, material.mu_const, material.rho_const)` into `flux_storage_` (size 1).  **Size `elems_to_flux_idx_` to `ne` and zero-initialise** (R-001) so `elems_to_flux_idx_[e] == 0` for every local element; every `At(elem)` reads `flux_storage_[0]`.  Extra cost: `ne * sizeof(int)` (≈ 20 MB on a 5 M-element mesh, acceptable).  Alternative tail: `At(elem)` may short-circuit to `flux_storage_.front()` when `mode == Mode::Constant`; either implementation produces a well-defined lookup for every local `elem`.
   - `Mode::GridFunction`: per element, evaluate `lambda_gf`, `mu_gf`, `rho_gf` at the element CENTROID (using `mesh.GetElementTransformation(e)` and the natural reference-element centroid `IntegrationPoint`), build a `GodunovFlux(λ_e, μ_e, ρ_e)`, and store it. **Per-element constant**: the flux is built from a single sampled value, NOT integrated. This matches the existing single-(λ,μ,ρ)-per-flux design.
   - `Mode::Coefficient`: same as GridFunction but use `material.EvalAt(e, T_e, ip_centroid, λ_e, μ_e, ρ_e)`.

2. **Memory optimisation — uniqueness map** (R-011): in `Mode::GridFunction` and `Mode::Coefficient`, many elements share `(λ, μ, ρ)` triples (deep layers where Vs(z) plateaus).  The pool uses a hash map keyed on the rounded triple to deduplicate; `flux_storage_` holds one `GodunovFlux` per unique key, `elems_to_flux_idx_` maps each LOCAL element to its key index.  Rounding precision is controlled by the constructor's `dedup_sig_figs` parameter (default 6).  The 6-sig-fig default works well for high-Vs basement (relative bucket ≈ 1e-6) but may dedup poorly in basin sediments; loosen to `dedup_sig_figs = 4` if the memory budget (see Acceptance Criteria) is exceeded.

   **Rounding algorithm** (R-010 round-3 — pinned to make `NumUniqueFluxes()` deterministic across compilers):
   ```cpp
   inline real_t round_to_sig_figs(real_t v, int sig_figs)
   {
       if (v == 0.0 || std::isnan(v)) { return v; }
       const real_t magnitude = std::pow(real_t{10},
                                         std::floor(std::log10(std::abs(v)))
                                         + real_t{1} - real_t(sig_figs));
       return std::round(v / magnitude) * magnitude;
   }
   ```
   Hash key is `std::tuple<real_t, real_t, real_t>` of the rounded triple; hash is `boost::hash_combine`-style XOR-mixing of element-wise `std::hash<real_t>` results (or equivalent — MFEM does not standardise this helper, so define it locally in `godunov_flux_pool.cpp`).  Do NOT use `snprintf("%g", ...)` for rounding: `%g` and direct numeric rounding can disagree at boundary values (e.g., `3.999995` rounds differently), producing different `NumUniqueFluxes()` across compilers.

3. **`At(int elem)` hot path**:
   ```cpp
   const GodunovFlux& GodunovFluxPool::At(int elem) const
   {
       MFEM_ASSERT(elem >= 0 && elem < static_cast<int>(elems_to_flux_idx_.size()),
                   "GodunovFluxPool::At: elem out of range");
       const int idx = elems_to_flux_idx_[elem];
       return flux_storage_[idx];
   }
   ```
   One indirection, branch-free in Release.  The assert pins R-001's lookup-validity contract in Debug builds.

4. **CFL helper**:
   ```cpp
   template <typename MeshType>
   real_t GodunovFluxPool<MeshType>::CpForElement(int elem) const
   { return At(elem).GetCp(); }

   template <typename MeshType>
   real_t GodunovFluxPool<MeshType>::MaxCp() const
   { return max_cp_; }
   ```
   `max_cp_` is a cached member computed at construction.  In the parallel
   specialisation (R-004 round-3) the ctor runs:
   ```cpp
   real_t local_max = -std::numeric_limits<real_t>::infinity();
   for (const auto& f : flux_storage_) { local_max = std::max(local_max, f.GetCp()); }
   if constexpr (IsParallelMesh<MeshType>::value)
   {
   #ifdef MFEM_USE_MPI
       MPI_Allreduce(&local_max, &max_cp_, 1,
                     MPITypeMap<real_t>::mpi_type, MPI_MAX,
                     mesh.GetComm());
   #endif
   }
   else
   {
       max_cp_ = local_max;
   }
   ```
   so callers can use `MaxCp()` without per-call MPI on either specialisation.

5. **`WaveOperator::ComputeMaxDt(real_t cfl)`** (R-005) changes from
   ```cpp
   return cfl_mixed_flux_factor * cfl * h_min_ / flux_.GetCp();
   ```
   to a per-element walk **followed by an MPI_Allreduce(MIN) of `dt_local`** (the existing code did NOT reduce inside `ComputeMaxDt`; it only reduced `h_min_` once at construction):
   ```cpp
   real_t dt_local = std::numeric_limits<real_t>::max();
   for (int e = 0; e < ne_; ++e)
   {
       const real_t h_e  = h_elem_[e];           // per-element char length (cached at ctor)
       const real_t cp_e = flux_pool_->CpForElement(e);
       dt_local = std::min(dt_local, cfl_mixed_flux_factor * cfl * h_e / cp_e);
   }
   real_t dt_global = dt_local;
   if constexpr (IsParallelMesh<MeshType>::value)
   {
   #ifdef MFEM_USE_MPI
       MPI_Allreduce(&dt_local, &dt_global, 1,
                     MPITypeMap<real_t>::mpi_type, MPI_MIN,
                     static_cast<ParMesh&>(mesh_).GetComm());
   #endif
   }
   return dt_global;
   ```
   The existing `h_min_` scalar becomes a per-element `h_elem_[ne_]` array (populated inside the same ctor loop that already walks every element).  The legacy `h_min_` MPI_Allreduce at lines 111-114 is REMOVED (replaced by the per-call reduce above).  The `h_min_` member variable itself MUST also be removed from the `WaveOperator` class (R-008 round-2 — leaving the member in place uninitialised would let future code silently consume a stale value).  Greppable: after Phase 3, `rg "\bh_min_\b" miniapps/seas/dynamic/wave_operator.*` returns zero hits.

6. **Face-flux dispatch.** Every site in `wave_operator.inl` that calls `flux_.<method>(...)` needs to choose WHICH side's flux to use. For symmetric methods (e.g. `Interior`, `Central`), the existing code uses a single `flux_` for both sides — heterogeneous case needs to handle a bi-material face. **First-pass design (this plan)**: at a bi-material face, use the AVERAGE flux.  `GodunovFlux::Interior(...)` is a NON-STATIC member function (R-007 round-2); the per-element flux is fetched from the pool first, then `.Interior(...)` is dispatched on each instance:
   ```cpp
   real_t F_h_self[9];
   real_t F_h_nbr [9];
   flux_pool_->At(elem_plus ).Interior(nor, Q_self, Q_nbr, F_h_self);
   flux_pool_->At(elem_minus).Interior(nor, Q_nbr, Q_self, F_h_nbr);
   for (int c = 0; c < 9; ++c)
   {
       F_h[c] = 0.5 * (F_h_self[c] + F_h_nbr[c]);
   }
   ```
   This is an approximation valid for small material contrasts; it is REASONABLE for SAFS (the largest contrast at a face is between adjacent tets in the size-graded mesh, ≤ a factor of ~2 in c_p, much less than e.g. soil/rock contrast). A more rigorous Riemann solver for the bi-material case is left for a follow-up plan; **this plan documents the approximation explicitly in `wave_operator.inl` next to the call site and in the new test T-WAVEOP-BIMATERIAL-AVG-FLUX**.

   For `Absorbing` and `FreeSurface` (boundary-only fluxes), there is only one element on each face; use that element's `flux_pool_->At(elem)`.

7. **Fault face flux** — leave the bi-material fault-face logic in `fault_face_flux.cpp` UNCHANGED in this phase. Fault impedances continue to be initialised from the operator's "default" `(λ, μ, ρ)` — i.e. from `flux_pool_->At(elem_plus_of_fault)` for one side and `flux_pool_->At(elem_minus_of_fault)` for the other. The existing impedance-initialisation API (`FaultFaceFlux::InitializeImpedancesFromMaterial`, sketched in `data_projection_feature_plan_v2.md` §Phase 5) is extended in a follow-up plan; here we wire up only what is needed for `Mode::Constant` parity.

7a. **Shared-face neighbour fluxes (parallel only) (R-003).**  Interior-side flux dispatch in §Detailed Req. 6 references `flux_pool_->At(elem_plus)` and `flux_pool_->At(elem_minus)`.  At shared faces (partition boundaries) the NEIGHBOUR element lives on a different MPI rank and has no entry in the local pool.  `GodunovFluxPool` MUST therefore expose a second hot-path lookup `AtNbr(int nbr_local_idx)` that returns the GodunovFlux for the foreign face-neighbour element, indexed by `ParMesh::face_nbr_elements` local index.

   At pool construction (after building the local `flux_storage_`), perform the neighbour-material exchange.  Pick option (a) for simplicity:
   - **(a) ParGridFunction exchange.** Build a **new** `mfem::L2_FECollection` of order 0 (DG0) and a **new** `mfem::ParFiniteElementSpace` over the supplied `ParMesh` with `vdim = 3` (one set of DG0 DOFs per scalar field, three fields).  Allocate a temporary `ParGridFunction` over THAT space.  **DO NOT reuse `WaveOperator::fes_`** (which is `L2_FECollection(order = p, ...)` and has `(p+1)(p+2)(p+3)/6` DOFs per element; reusing it would leave `ndof_per_el − 1` DOFs uninitialised and the subsequent `ExchangeFaceNbrData` would carry garbage λ/μ/ρ values; R-002 round-2).  Set per-element values `(λ_e, μ_e, ρ_e)` at each local element's UNIQUE DG0 DOF.  Call `gf.ExchangeFaceNbrData()` once.  Walk the resulting `face_nbr_data` to obtain `(λ, μ, ρ)` for each face-neighbour element, build `GodunovFlux` objects, and store them in a **separate** array `face_nbr_flux_storage_` keyed by `nbr_local_idx`.  Optionally deduplicate against the local `flux_storage_` via the same rounded-triple hash to avoid double-allocating common triples — when a neighbour's rounded `(λ, μ, ρ)` matches a local entry, point `nbr_local_idx_to_flux_idx_[nbr_local_idx]` at the local `flux_storage_` slot and store an empty entry in `face_nbr_flux_storage_`.  Build `nbr_local_idx_to_flux_idx_` parallel to `elems_to_flux_idx_`.
   - **(b) Coefficient direct.** Alternatively walk `pmesh.GetFaceNbrElementTransformation(...)` and call `material.EvalAt(...)` on each neighbour transformation.  This avoids the ParGridFunction allocation but requires `EvalAt` to accept a Coefficient handle independently of element index — which it does in `Mode::Coefficient`.

   In `Mode::Constant` the entire path is a no-op (single shared flux suffices); `AtNbr` returns `flux_storage_.front()`.

   `ComputeSharedFaceFluxRHS` and `ComputeADERSharedFaceFluxRHS` are updated to fetch local-side flux via `flux_pool_->At(elem_local)` and foreign-side via `flux_pool_->AtNbr(elem_nbr_local_idx)`, then apply the averaged-flux pattern from §Detailed Req. 6.

8. **Existing `Ax_`, `Ay_`, `Az_` cached `DenseMatrix` members at WaveOperator scope** (set at lines 55-57 of the ctor) are REMOVED. All references in `wave_operator.inl` (line 1055 fetches `flux_.GetReferenceStarMatrix(d)`) change to `flux_pool_->At(e).GetReferenceStarMatrix(d)`.

### Interfaces
- `class GodunovFluxPool` (new) with methods `At(elem)`, `CpForElement(elem)`, `MaxCp()`.
- `WaveOperator::WaveOperator(mesh, order, const MaterialField&, const BoundaryConfig&)` (new ctor).
- All `GodunovFlux` calls remain; only their owner changes.

### Edge Cases to Handle
- **Element on a fault face has different material on each side.** This is the bi-material flux case discussed in Detailed Req. 6; using the averaged flux is documented as an approximation, and the test T-WAVEOP-BIMATERIAL-AVG-FLUX asserts that for a homogeneous fixture both sides yield the same `flux_pool_->At(elem)` and the average reduces to the standard flux.
- **MPI element ownership.** The pool indexes by LOCAL element ID. Each MPI rank constructs its own pool from its local elements.  In the **parallel specialisation** `GodunovFluxPool<ParMesh>`, `MaxCp()` does an `MPI_Allreduce(MAX)` over the per-rank local maxima at construction time (so callers can use it without per-call MPI).  In the **serial specialisation** `GodunovFluxPool<Mesh>`, `MaxCp()` returns the local maximum directly with no MPI calls.  The `MIN dt` reduction over `dt_local` is **per-call**, performed inside `WaveOperator::ComputeMaxDt` per §Detailed Req. 5 (the legacy ctor-time reduction of `h_min_` is REMOVED).
- **Empty rank.** If a rank has zero local elements (unusual but possible at high MPI counts), the per-rank `local_max_cp_` initialises to `-std::numeric_limits<real_t>::infinity()` (R-008 round-3 — using `-inf` instead of `0.0` so that any legitimate `c_p = 0` value would still survive `MPI_Allreduce(MAX)`); the reduction recovers the correct global value from non-empty ranks.  An all-empty global mesh is impossible because `FieldProjector::ComputeMeshBBoxParallel` aborts upstream (`field_coefficient.cpp` R-009 guard).

### Acceptance Criteria
- [ ] `make test-tpv102-explicit-rk`, `make test-tpv104-rk` (or whatever current dynamic-rupture tests are named) pass byte-identically with `Mode::Constant`.
- [ ] New test T-WAVEOP-CONST-PARITY constructs `WaveOperator` two ways — old `(λ, μ, ρ)` ctor vs new `MaterialField::MakeConstant(λ, μ, ρ)` ctor — and asserts the L∞ norm of `Mult` outputs on a synthetic state matches to within 1e-12 of the norm.
- [ ] New test T-WAVEOP-COEFF-PARITY constructs `WaveOperator` with `Mode::Coefficient` over a synthetic constant `DataField3D` (every voxel = same value) and asserts the result matches the equivalent `MakeConstant` ctor to within 1e-10 (slight slack because the Coefficient path samples at element centroids instead of using the precomputed constant Jacobians directly).
- [ ] CFL `Δt` for the TPV102 fixture is bit-identical to today's value when using `Mode::Constant`.
- [ ] T-WAVEOP-BIMATERIAL-SHARED-FACE-MPI (R-003) passes at `mpirun -np 2`.
- [ ] **Memory budget (R-011 / R-001 round-3):** on the SAFS z-graded 3.2 M-tet fixture with `Mode::Coefficient` and the 39-slice sidecar, `GodunovFluxPool` allocates AT MOST 2 GB of HEAP-INCLUSIVE memory for `flux_storage_`.  **Do NOT use `sizeof(GodunovFlux)`** in the assertion — that returns only the in-struct layout and ignores the heap-allocated `DenseMatrix::data` arrays inside each instance (R-001 round-3).  The concrete assertion is:
  ```cpp
  // Heap-inclusive per-instance bytes:
  //   6 DenseMatrix * 81 doubles * 8 bytes = 3888 bytes of matrix data
  //   plus headers + 7 scalars ≈ 256 bytes
  constexpr size_t kGodunovFluxBytes = 6 * 81 * sizeof(double) + 256;
  // = 4144 bytes.
  assert(pool.NumUniqueFluxes() * kGodunovFluxBytes <= 2'000'000'000);
  ```
  If exceeded: (a) loosen `dedup_sig_figs` (from 6 to 4), or (b) drop the cached `ref_star_[3]` matrices on heterogeneous-mode paths and revisit Option (b) (JIT Jacobian build).

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 4 (driver wiring for dynamic rupture).

---

## Phase 4: Driver wiring — `--material-mode` and `--sidecar`

### Goal
Add the CLI plumbing and material-field construction logic to both consumer driver families (quasi-dynamic / SAFS BP5; dynamic-rupture TPV102/104/205). Each driver gains exactly one new code path that is reached only when `--material-mode != constant`. Default behaviour is byte-identical to today.

### Files to Create
- `miniapps/seas/io/sidecar_material_factory.hpp`, `.cpp` — small factory:
  ```cpp
  namespace mfem::seas {
  struct SidecarMaterialBundle
  {
      // Owns the DataField3D objects (long-lived; outlives the operator).
      std::unique_ptr<DataField3D> vp_field;
      std::unique_ptr<DataField3D> vs_field;
      std::unique_ptr<DataField3D> rho_field;
      // Owns the Coefficient subclasses, which hold references to the above.
      std::unique_ptr<mfem::Coefficient> lambda_coef;
      std::unique_ptr<mfem::Coefficient> mu_coef;
      std::unique_ptr<mfem::Coefficient> rho_coef;
      // Convenience: produces a MaterialField pointing to the Coefficients.
      MaterialField MakeMaterialField() const;
  };
  SidecarMaterialBundle LoadSidecarMaterialBundle(
      const std::string& sidecar_path,
      mfem::ParMesh& pmesh,                // for pre-flight bbox check
      InterpMode interp = InterpMode::Trilinear);
  } // namespace mfem::seas
  ```
  Implementation:
  1. Construct three `DataField3D` objects (one each for Vp, Vs, density).
  2. Set `InterpMode` on all three.
  3. Compute the mesh bbox via `FieldProjector::ComputeMeshBBoxParallel`.  For EACH of the three loaded fields, if `DataField3D::ContainsBBox(mxmin, mxmax, mymin, mymax, mzmin, mzmax)` returns false, call the **existing** helper `FieldProjector::AbortContainmentFailure(field, mxmin, mxmax, mymin, mymax, mzmin, mzmax)` (in `field_coefficient.cpp`) — that helper already prints both bbox tables and the field's name (via `DataField3D::FieldName()`), so the implementer should NOT roll a new abort path (R-009 round-3).  Checking all three defends against future schema or filesystem corruption that leaves only one field's axes wrong (R-006 round-2).
  4. Construct the three Coefficient subclasses.
  5. Return the bundle.

### Files to Modify
- **Quasi-dynamic driver** — pick `bp5_full_driver.cpp` or whichever the SAFS / BP5 driver is. Add:
  ```cpp
  args.AddOption(&material_mode_str, "-mm", "--material-mode",
                 "constant | sidecar-coefficient | sidecar-gridfunction "
                 "(default constant).");
  args.AddOption(&sidecar_path, "-sc", "--sidecar",
                 "Path to schema-v1 HDF5 sidecar (required when "
                 "--material-mode is sidecar-*).");
  ```
  Then:
  ```cpp
  std::unique_ptr<SidecarMaterialBundle> bundle;   // optional; outlives operator
  MaterialField material;
  if (material_mode_str == "constant") {
      material = MaterialField::MakeConstant(lambda_const, mu_const, rho_const);
  } else if (material_mode_str == "sidecar-coefficient") {
      MFEM_VERIFY(!sidecar_path.empty(),
                  "--sidecar PATH is required when --material-mode=sidecar-coefficient");
      bundle = std::make_unique<SidecarMaterialBundle>(
          LoadSidecarMaterialBundle(sidecar_path, pmesh));
      material = bundle->MakeMaterialField();
  } else if (material_mode_str == "sidecar-gridfunction") {
      // Existing FieldProjector::ProjectVelocity path
      auto vf = FieldProjector::ProjectVelocity(sidecar_path, target_fes, ...);
      material = MaterialField::MakeGridFunction(vf.rho, vf.lambda, vf.mu);
  } else {
      MFEM_ABORT("Unknown --material-mode '" << material_mode_str << "'");
  }

  ElasticityDomainOperator<ParMesh> op(pmesh, order, material,
                                        Vp, Wf, lf, bdr_config, ...);
  ```
- **Dynamic-rupture drivers** — `tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp`: same wiring. The line that constructs `WaveOperator` (`tpv104_driver.cpp:1033-1035`) becomes:
  ```cpp
  WaveOperator<MeshT> wave(pmesh, order, material, bc);
  ```
- `miniapps/seas/Makefile` — link the new `material_coefficients.cpp`, `sidecar_material_factory.cpp`, `godunov_flux_pool.cpp` into the seas common library; add new test targets (Testing Strategy below).

### Detailed Requirements
1. **CLI defaults:**
   - `--material-mode` default = `"constant"`. With this value, NEITHER `DataField3D` NOR `FieldCoefficient` NOR any sidecar I/O is touched. The driver behaviour is byte-identical to today.
   - `--sidecar` default = `""`. Reading from an empty path is forbidden — checked at parse time only if `--material-mode` is sidecar.

2. **Order of operations** in the driver (must match exactly):
   1. Parse CLI.
   2. Load mesh / build ParMesh.
   3. **Before constructing operator**: if sidecar mode, load `SidecarMaterialBundle` (this includes the pre-flight bbox containment check against the full ParMesh).
   4. Construct operator with the chosen `MaterialField`.
   5. Run simulation.
   Critically the `SidecarMaterialBundle` lives in the driver scope (typically in `main`), giving it lifetime ≥ that of the operator. The operator stores only non-owning Coefficient pointers via `MaterialField`.

3. **Pre-flight bbox check semantics.** `LoadSidecarMaterialBundle` aborts (via `MFEM_ABORT`) if any of the three sidecar fields' bbox does not contain the mesh's bbox. Same message format as the existing `FieldProjector::AbortContainmentFailure`. This is identical to the existing pre-flight behaviour in the GridFunction path.

4. **Sidecar paths and interp mode.** `LoadSidecarMaterialBundle` takes `InterpMode` as a parameter (default `Trilinear`). For now hardcode the driver-side to pass the default; a follow-up plan can expose this on the CLI if needed.

5. **No silent fallback.** If `--material-mode` is sidecar but the sidecar load or pre-flight check fails, the driver aborts with the exact error from `MFEM_ABORT` inside `DataField3D` / `FieldProjector`. No "default to constant" path.

6. **Audit `GetShearModulus` / `GetLambda` consumers (R-003 round-2).**  In heterogeneous modes `Mode::Coefficient` and `Mode::GridFunction`, the operator's scalar accessors `GetLambda()` and `GetShearModulus()` return `quiet_NaN()` as a sentinel (see Phase 2 Detailed Req. 5).  Existing SEAS code uses these accessors at multiple sites OUTSIDE the bilinear-form assembly path — chiefly friction-state initialisation, where the bulk impedance `η = √(μ·ρ)` seeds Brent's solver.  If those callers see NaN, the simulation NaN-propagates at step 1 with no informative trace.

   Before enabling heterogeneous mode on any driver, run
   ```
   rg "GetShearModulus|GetLambda" miniapps/seas/
   ```
   and review every match.  Any match that uses the returned value in a numerical formula (friction impedance, pre-stress amplitude, CFL helpers, etc.) MUST be rewritten to evaluate the operator's `MaterialField` at the relevant location (fault QP, element centroid, etc.) when the active mode is `Mode::Coefficient` or `Mode::GridFunction`.  Existing `Mode::Constant` callers continue to work via the same scalar accessors — they return the constructor-supplied real value, not NaN.

   The expected audit set today includes:
   - friction state init in `solver/seas_operator.hpp`
   - pre-stress / nucleation in `drivers/tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp`
   - any utility in `friction/` that hard-codes a Constant-coefficient model

   Practical recipe for fault-DOF material lookup (R-007 round-3 — the most common audit-target context).  Given a fault QP at 3D physical coords `(x, y, z)` and the local bulk element index `elem` that owns it (already in `FaultBasis` / `fault_face_flux.cpp` bookkeeping):
   ```cpp
   ElementTransformation* T = mesh.GetElementTransformation(elem);
   IntegrationPoint ip;
   // Inverse-map physical → reference coords of element `elem`.
   {
       mfem::Vector phys(3);
       phys(0) = x; phys(1) = y; phys(2) = z;
       T->TransformBack(phys, ip);
   }
   real_t lambda, mu, rho;
   material.EvalAt(elem, *T, ip, lambda, mu, rho);
   const real_t eta_p = std::sqrt((lambda + 2 * mu) * rho);
   const real_t eta_s = std::sqrt(mu * rho);
   ```
   Reuse `FaultBasis::GetFaultQPCoords()` (or equivalent) to obtain the `(elem, x, y, z)` tuple per fault QP without re-walking mesh face geometry.

   Document the audited list in the driver's commit log.  A new acceptance test (T-CLI-NO-NAN-IN-COEFF-MODE; see Testing Strategy) verifies that no fault DOF receives a NaN impedance in heterogeneous mode.

### Interfaces
- `LoadSidecarMaterialBundle(sidecar_path, pmesh, interp) -> SidecarMaterialBundle`.
- `SidecarMaterialBundle::MakeMaterialField() -> MaterialField`.
- New CLI flags: `--material-mode`, `--sidecar`.

### Edge Cases to Handle
- **Mesh doesn't fit sidecar bbox** → MFEM_ABORT during pre-flight with both bbox tables printed.
- **`--material-mode=sidecar-...` but no `--sidecar PATH`** → CLI parse error message.
- **Sidecar lacks Vp/Vs/density field** → MFEM_ABORT inside `DataField3D` constructor.
- **Driver passed both `--material-mode=constant` AND `--sidecar PATH`** → log a warning that `--sidecar` is ignored. Do not abort. Tests don't pin this; it is a UX-only rule.

### Acceptance Criteria
- [ ] `mpirun -np 1 bp5_full_driver --mesh bp5/mesh/bp5_1000m.msh --tfinal …` (existing invocation, no new flags) produces byte-identical output to today's reference.
- [ ] `mpirun -np 4 bp5_full_driver --mesh safs_fault_box_nwcut_500m_zgraded.msh --tfinal … --material-mode sidecar-coefficient --sidecar data_projected/velocity_safs.h5` runs to completion without abort.
- [ ] Same two invocations for `tpv102_driver` and `tpv104_driver`.
- [ ] Existing test suite (`make test`) passes with no modifications.

### Dependencies
- Depends on: Phase 2 (quasi-dynamic ctor) and Phase 3 (dynamic-rupture ctor).
- Required by: nothing in this plan; production fault-zone simulations beyond it.

---

## Testing Strategy

All tests live in `miniapps/seas/tests/unit/`. New test targets are added to `miniapps/seas/Makefile` (`test-<name>` and `test-<name>-parallel` as relevant). The Test Catalog below numbers tests per phase.

### Phase 1 tests (extend `test_heterogeneous_material.cpp`)

- **T-5-4 — Mode::Coefficient round-trip.** Build three `ConstantCoefficient` instances (`lambda = 1.0`, `mu = 2.0`, `rho = 3.0`); wrap in `MaterialField::MakeCoefficient`; call `EvalAt(0, T, ip, λ, μ, ρ)` on a synthetic linear-tet element transformation; assert `λ==1.0 && μ==2.0 && ρ==3.0`.
- **T-5-5 — `At` aborts in Mode::Coefficient.** Build a Coefficient-mode `MaterialField` and call `At(0, 0, …)`; assert the call aborts (use a `try/catch`-style harness or the existing `EXPECT_ABORT` macro if available; otherwise wrap in a subprocess test).
- **T-5-6 — `MaxCpInElement` in Mode::Coefficient.** Build a synthetic `DataField3D` whose `vp_field`, `vs_field`, `rho_field` produce known constants. Build `LambdaFromSidecar`, `MuFromSidecar`, `RhoFromSidecar`. Build `MaterialField::MakeCoefficient(...)`. Call `MaxCpInElement(0, T)` and assert it matches `sqrt((λ + 2μ) / ρ)` to FP roundoff.
- **T-5-7 — `LambdaFromSidecar` numerical correctness.** With synthetic `DataField3D`s of Vp = 5000, Vs = 3000, ρ = 2700 (matching the constant-medium BP5 baseline), assert `LambdaFromSidecar::Eval` returns `2700 * 5000² − 2 * 2700 * 3000² = 2700*(25e6 − 18e6) = 1.89e10 Pa` to FP roundoff.

### Phase 2 tests (new `test_elasticity_operator_coefficient_mode.cpp`)

- **T-EOP-CONST-REGRESSION.** Construct `ElasticityDomainOperator` two ways on the same mesh — legacy `(λ, μ)` ctor vs new ctor with `MaterialField::MakeConstant(λ, μ, ρ)` — and assert the assembled stiffness matrix entries match bit-identically (`HypreParMatrix::Read` both, subtract, assert `L∞ < machine_eps`). Fixture: BP5 1000 m mesh, np=1 and np=4.
- **T-EOP-COEFF-CONST-PARITY.** Construct `ElasticityDomainOperator` with `Mode::Coefficient` where the three Coefficients are `ConstantCoefficient` instances wrapping the same scalars as the constant ctor; assert the assembled stiffness matrix entries match to within `1e-10 * |K|` (slight slack because `Coefficient::Eval` goes through `T.Transform(ip, x)` and back, accumulating FP roundoff vs. the direct `ConstantCoefficient` path).
- **T-EOP-COEFF-LAYERED.** Construct a synthetic layered sidecar (Vs = 1000 in z > -3000m, Vs = 3000 in z ≤ -3000m; same Vp scaling; constant ρ). Project velocity onto the mesh in both `Mode::GridFunction` (P1) and `Mode::Coefficient` paths. Solve a Laplace-style problem (constant Dirichlet BC, zero RHS) and compare both solutions. The `Mode::Coefficient` solution should resolve the layer at z=-3000m sharply; the `Mode::GridFunction` solution should show the smearing predicted by the debugging session. **Acceptance** (R-008): max pointwise RELATIVE difference between the two displacement solutions, `||u_coef − u_gf||_∞ / ||u_coef||_∞`, is at least `1e-3` — proves the GridFunction path smears the sharp basin layer in a way the Coefficient path does not.  Both solutions are displacement fields (units: m); the comparison threshold is dimensionless.  The exact threshold is fixture-dependent and is pinned by the fixture builder.
- **T-EOP-MATERIAL-LIFETIME.** Construct an operator with `Mode::Coefficient`; immediately destroy the `MaterialField` (out of scope); call `Solve()` on the operator; assert SEGFAULT or `MFEM_ABORT` is reached. This pins the lifetime contract documented in the constructor comment. (If running this as a unit test is awkward, document it as a manual check.)

### Phase 3 tests (new `test_wave_operator_coefficient_mode.cpp`)

- **T-WAVEOP-CONST-PARITY.** Construct `WaveOperator` legacy and new ctors with same `(λ, μ, ρ)`; apply `Mult` to a random Q; assert L∞ of difference < 1e-12.
- **T-WAVEOP-COEFF-PARITY.** Construct `WaveOperator` with `Mode::Coefficient` over synthetic constant sidecar; assert `Mult` matches `Mode::Constant` to within 1e-10.
- **T-WAVEOP-CFL-CONST-PARITY.** TPV102 fixture; `ComputeMaxDt(0.5)` returns the same value when called via legacy and new ctors. Assert exact equality (the c_p is the same constant).
- **T-WAVEOP-CFL-HETERO.** Build a synthetic two-layer sidecar (low-Vs basin top, high-Vs basement); construct `WaveOperator` with `Mode::Coefficient`. Assert `ComputeMaxDt` returns a value bracketed by `h / cp_top` and `h / cp_bottom` and is less than `h / cp_bottom` (i.e. the basement c_p drives the global dt limit if it is the larger one, or vice versa). Tests pin the correctness of per-element c_p reduction.
- **T-WAVEOP-BIMATERIAL-AVG-FLUX.** Two adjacent elements with different (λ, μ, ρ). Compute the bi-material flux at the shared face two ways: (a) using the average-flux approximation in `wave_operator.inl`, (b) by hand-computing one side's `GodunovFlux::Interior` with one material and the other side's with the other material. Assert (a) equals 0.5*((b1) + (b2)). This is a parity test against the design choice, not a numerical-correctness test against a reference solution.
- **T-WAVEOP-BIMATERIAL-SHARED-FACE-MPI** (R-003).  Two MPI ranks; insert a `(λ, μ, ρ)` discontinuity that crosses the partition boundary so each rank sees different material at adjacent elements.  Assert that `ComputeSharedFaceFluxRHS` does not abort or OOB on either rank, and that the flux from each side respects momentum conservation: `||F_rank0 + F_rank1||_∞ < 1e-10 * ||F_rank0||_∞`.  Pins R-003's neighbour-material exchange.
- **T-WAVEOP-COEFF-CENTROID-PARITY** (R-006 / R-009 replacement for T-WAVEOP-COEFF-QP-EQUALITY).  Build a Coefficient-mode WaveOperator over a small mesh.  For each local element `e`:
  1. Compute the centroid sample `(λ_c, μ_c, ρ_c) = material.EvalAt(e, T_e, ip_centroid, …)`.
  2. Compute `cp_c = sqrt((λ_c + 2 μ_c) / ρ_c)` (m/s).
  3. Assert `|flux_pool_->CpForElement(e) − cp_c| < 1e-12 * cp_c` (relative roundoff, unit-consistent m/s).

  This pins the per-element flat-material design (the per-element flux is built from a single centroid sample, NOT integrated over qpoints).  The test honestly characterises what the dynamic-rupture path delivers: material at each qpoint INSIDE the flux Jacobian's contribution to assembly equals the centroid sample, not the per-qp sidecar value.  See R-006 for why the previous "1.0 m/s" claim was unsatisfiable.
- **T-WAVEOP-COEFF-WITHIN-ELEM-VARIATION** (diagnostic, no hard assertion).  On a basin-straddling tet, report `max_qp |c_p(qp) − c_p(centroid)|` in m/s as a physics diagnostic.  This number quantifies the limitation of the centroid-based design — it is NOT a pass/fail criterion, but a number whose growth over time indicates a follow-up plan should add per-qp flux Jacobians.

### Phase 4 tests (new `test_driver_material_mode_cli.cpp` + integration tests)

- **T-CLI-CONSTANT-REGRESSION.** Spawn `bp5_full_driver --tfinal 100s` (short) without `--material-mode`; capture output `.vtu` files; diff bit-for-bit against a stored reference. Same for `tpv102_driver`. (These reuse the existing CI-style golden test fixtures.)
- **T-CLI-SIDECAR-COEFFICIENT-SMOKE.** Spawn the same drivers with `--material-mode sidecar-coefficient --sidecar tests/fixtures/synthetic_sidecar.h5 --mesh tests/fixtures/small_sidecar_fitted.msh --tfinal 1s`; assert exit code 0 and that the produced `.vtu` files are non-empty.
- **T-CLI-SIDECAR-MISSING-PATH.** Spawn drivers with `--material-mode sidecar-coefficient` and no `--sidecar`; assert MFEM_ABORT with the expected message substring.
- **T-CLI-MESH-OUT-OF-SIDECAR.** Spawn drivers with a mesh whose bbox lies outside the synthetic sidecar; assert MFEM_ABORT with the expected bbox-table message.
- **T-CLI-VALGRIND-LIFETIME.** Run `valgrind --tool=memcheck` against the smoke driver in T-CLI-SIDECAR-COEFFICIENT-SMOKE; assert zero "definitely lost" or "indirectly lost" bytes. This pins the `SidecarMaterialBundle` lifetime contract under real driver flow.
- **T-CLI-NO-NAN-IN-COEFF-MODE (R-003 round-2).**  Drive the SAFS / BP5 quasi-dynamic and TPV dynamic-rupture drivers with `--material-mode sidecar-coefficient` on the small fixture.  After friction-state initialisation (and BEFORE the first time step), assert that every fault DOF's impedance `η_p` / `η_s` is finite (`std::isfinite`).  Pins the Phase-4 audit step.

### Test fixtures to create
- `tests/fixtures/synthetic_sidecar.h5` — a tiny (10×10×10) sidecar built from a known Vp/Vs/ρ field in Python (schema-v1). Reused by T-5-4 onwards.
- `tests/fixtures/synthetic_layered_sidecar.h5` — a 20×20×40 sidecar with a sharp Vs step at z=-3000m. Reused by T-EOP-COEFF-LAYERED.
- `tests/fixtures/small_sidecar_fitted.msh` — a tiny tet mesh that fits inside both sidecar fixtures; reused by smoke tests.

A small Python helper `tests/fixtures/build_synthetic_sidecars.py` constructs the three fixtures from analytic functions; the fixtures themselves are committed to the repo for reproducibility.

---

## Risk Assessment

### High-confidence (low risk)
- **Phase 1 and Phase 2 are mostly mechanical.** The DG integrators already accept `Coefficient&`, so the quasi-dynamic path is a search-and-replace plus a new constructor.  The `MaterialField` struct extension is additive at compile time (R-014): no existing constructor or accessor signature changes.  At runtime the legacy `At(elem, dof, ...)` accessor deliberately aborts in `Mode::Coefficient` — callers that opt into the new mode must audit their accessor calls and move to `EvalAt(elem, T, ip, ...)`.
- **Phase 4 (CLI wiring)** is the standard pattern used by every other SEAS driver.

### Medium-risk
- **GodunovFlux per-element storage memory.** A 1 M tet mesh holds ~5 M elements (incl. fault QP doubling); 5 M × 9² × 8 B = 3.24 GB. With the uniqueness map this drops to maybe 500 MB (one unique flux per unique `(λ, μ, ρ)` triple). Mitigation: the uniqueness map; further mitigation if needed is to keep only the `Ax_*` half of `GodunovFlux` (skipping the cached `ref_star_` triple) per element.
- **Bi-material face flux approximation.** Using the AVERAGED Godunov flux at bi-material faces is a known compromise. Mitigation: T-WAVEOP-BIMATERIAL-AVG-FLUX pins the current behavior so a future correctness improvement is detectable. The plan documents that an exact bi-material Riemann solver is follow-up work.

### Low-confidence (need investigation before implementing)
- **Fault impedance computation (`FaultFaceFlux::InitializeImpedancesFromMaterial`).** Phase 5 of `data_projection_feature_plan_v2.md` sketches a version of this against the `Mode::GridFunction` path. For `Mode::Coefficient`, we want the qpoint-direct version. **This plan does NOT yet specify the fault-side impedance update beyond what `Mode::Constant` does**; a follow-up plan should extend `FaultFaceFlux` to read material at each fault QP via `EvalAt(elem, T, ip, …)`. Investigation needed: confirm that `fault_face_flux.cpp` impedances are EFFECTIVELY constant in our current Mode::Constant production code, so the missing extension does not break TPV102 / TPV104 regression.
- **Element-centroid sampling for `GodunovFlux` per-element triple.** Sampling `(λ, μ, ρ)` at the element CENTROID (Detailed Req. Phase 3.1) is a per-element constant approximation. For tall fault-corridor tets, the element-centroid value may not be representative; a single-point sample at centroid hides the variation that the qpoint-direct Coefficient is designed to expose. This is acceptable for `Mode::Coefficient`'s elastic CFL / Godunov flux (which IS a per-element flat-material method by construction; the heterogeneity is exposed only inside the volume / face quadratures). Document this in the `GodunovFluxPool` constructor.
- **Performance.** Per-qpoint Coefficient evaluation is ~5× slower than `ConstantCoefficient::Eval` (which is constant-fold). A 2× total wall-clock slowdown for `Mode::Coefficient` is expected and acceptable; a 5× slowdown would warrant revisiting. **No specific benchmark target is set in this plan**; the production-mesh smoke tests will reveal actual costs.

### Known tricky areas in existing code
- **`elasticity_operator.hpp` lines 127-135** (constant ctor body) — the `dynamic_cast<const LinearElastic*>(model_)` is needed to extract scalar λ/μ for the legacy `lambda_val_` cache. This stays as-is; the new ctor does not go through `LinearElastic`.
- **`wave_operator.inl` lines 55-57** — the `Ax_, Ay_, Az_` cached matrices are touched in many places. Mitigation: do the rewrite mechanically through Phase 3 and audit by running the existing TPV104 regression test BEFORE moving to Phase 4.
- **`fault_face_flux.cpp`** has its own impedance precomputation that this plan does not modify. Mitigation: tests T-WAVEOP-CONST-PARITY and existing TPV104 regression are the safety net — they fail if fault impedances change unexpectedly.

---

## Out of Scope for This Plan

The following are deliberately left for follow-up plans:

1. **Exact bi-material Riemann solver for the wave operator.** Phase 3 uses an averaged-flux approximation. A rigorous Riemann solver across bi-material faces (Pelties et al. 2012 §3 or de la Puente et al. 2009 §2.3) is a separate planning effort.
2. **`FaultFaceFlux::InitializeImpedancesFromMaterial` Coefficient-mode update.** This plan leaves fault-face impedances in their current `Mode::Constant` form. A follow-up will compute per-fault-qp impedances via `material.EvalAt(elem, T, ip, …)`.
3. **Time-dependent material.** Both modes assume a STATIC sidecar (no time variation). Damage-dependent material is out of scope.
4. **GPU offload.** `MaterialField::EvalAt` in `Mode::Coefficient` cannot run on accelerators (CPU-only HDF5 + trilinear). Any future GPU port will need a device-side material representation that this plan does not cover.
