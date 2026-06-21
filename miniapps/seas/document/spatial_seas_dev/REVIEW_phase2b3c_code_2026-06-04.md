# Code Review: Phase 3c (bp5_analytic prestress) + Phase 2b Stage 1a (operator coefficient infra) — 2026-06-04

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` — §Phase 3c (R-007/OBS-1) and §Phase 2b (heterogeneous material), Stage 1a only.
- Files reviewed:
  - `spatial/code/spatial_stress.hpp` — new `Bp5AnalyticStressSource` functor (+ `config/bp5_params.hpp` include).
  - `tests/unit/test_spatial_seas_bp5_analytic_prestress.cpp` — new parity test.
  - `domain/elasticity_operator.hpp` — new `MaterialCoefficient` adapter, member type switch, heterogeneous `Coefficient&` ctor, NaN-sentinel getters.
  - `Makefile` — `seas_test_spatial_seas_bp5_analytic` target (vars/link/obj).
- Domain context: `CLAUDE.md` (extreme-care `domain/`; sign conventions; no silent fallback), the QD plan (R-002 enum supersession, R-007), the archived `safs/project_7.0_alternative/document/heterogeneous_material_plan.md` (R-002 NaN sentinels, R-003 getter-consumer audit), `fault/fault_geometry_safs_templated.inl` (Cauchy projection convention), `config/bp5_params.hpp` (`tau0_vec`).
- Method: three adversarial passes; build + run exercised (elasticity_operator 461/461, elasticity_br2 46/46, faultgeom_parity np1/np4, bp5_analytic np1/np4 all pass; driver objects recompile clean). Verified the functor projection algebra by hand (`n·S·n=σ_n`, `t1·S·n=τ_dip`, `t2·S·n=τ_strike`).

## Findings

### [R-401] MODERATE [BUG] domain/elasticity_operator.hpp:ElasticityDomainOperator(Coefficient&,...) — heterogeneous ctor ships but RHS assembly + traction recovery still read NaN scalar λ,μ

**Category:** BUG (latent silent-NaN) / DEVIATION (incomplete feature behind a public ctor)

**Description:**
The Stage-1a heterogeneous ctor sets `lambda_val_ = mu_val_ = quiet_NaN()` and delegates λ,μ to external coefficients. But four runtime methods still consume the scalar `lambda_val_/mu_val_`, NOT the coefficient: `AssembleSlipContributionIP` (`assembly.inl:650/726/752`), `AssembleSlipContributionIPShared` (`:977/1039`), `AssembleDirichletLoading` (`:1240/1291/1506/1559/1604/1801/1853`), and `ComputeTractionImpl` (`traction.inl:1189/1206/1661/1675`). For an operator built via the new ctor, every `Solve()` RHS term and every recovered traction becomes **NaN, silently**. `AssembleStiffness` (K) uses only the coefficient, so K itself is fine — which makes the failure mode worse: K assembles cleanly, then the slip/Dirichlet RHS is NaN.

Currently this is **latent**, not live: `spatial_seas_driver.cpp` §4.0 aborts on heterogeneous configs, so the new ctor has no caller yet. **It becomes CRITICAL the moment Stage 2 removes that abort if Stage 1b has not landed.**

**Trigger:** `ElasticityDomainOperator op(mesh, order, lambda_coef, mu_coef, …); op.Solve(t, slip, u);` or `op.ComputeTraction(...)`.

**Actual behavior:** NaN RHS / NaN traction, no diagnostic.

**Expected behavior:** either correct heterogeneous tractions (Stage 1b), or a loud abort until Stage 1b lands.

**Suggested fix:** add a temporary guard in the heterogeneous ctor so the half-wired operator cannot silently produce NaN. Remove it in the Stage 1b commit that generalizes the four functions.
```diff
       lambda_coeff_.SetExternal(&lambda_c);
       mu_coeff_.SetExternal(&mu_c);
+
+      // STAGE 1a GUARD (remove in Stage 1b): the slip/Dirichlet RHS assembly
+      // and traction recovery still read scalar lambda_val_/mu_val_ (now NaN).
+      // Refuse to construct a half-wired heterogeneous operator until Stage 1b
+      // generalizes AssembleSlipContributionIP / ...IPShared /
+      // AssembleDirichletLoading / ComputeTractionImpl to coeff.Eval(T,ip).
+      MFEM_ABORT("ElasticityDomainOperator heterogeneous ctor: traction/RHS "
+                 "assembly not yet generalized to per-qp coefficients "
+                 "(Phase 2b Stage 1b). Do not enable the heterogeneous driver "
+                 "path until Stage 1b lands.");
```
(Alternative if a constructible-but-K-only operator is wanted for an interim test: gate the abort on first `Solve`/`ComputeTraction` instead. Either way, no silent NaN.)

**Test case:**
```cpp
// test_R401_hetero_ctor_no_silent_nan:
//   ConstantCoefficient lc(lambda), mc(mu);   // homogeneous values via the
//                                             // heterogeneous code path
//   build op via the Coefficient& ctor on a small BP5 mesh;
//   EITHER expect a clean abort (Stage 1a guard),
//   OR (post-Stage-1b) Solve + ComputeTraction and assert
//      std::isfinite on every traction entry AND bit-for-bit vs the
//      LinearElastic(lambda,mu) constant ctor.
```

---

### [R-402] MODERATE [BUG] Makefile:$(TEST_SPATIAL_SEAS_BP5_ANALYTIC_OBJ) — obj rule omits spatial_stress.hpp / templated .inl prereqs → stale-test risk

**Category:** BUG (build hygiene — silent stale test)

**Description:**
The new object rule lists `$(SEAS_HEADERS)` as its only header prerequisite. `SEAS_HEADERS` (Makefile:839) is `FRICTION/DOMAIN/FAULT/SOLVER/IO/CHECKPOINT/CONSTITUTIVE/CONFIG/DYNAMIC` headers and does **not** include `spatial/code/spatial_stress.hpp` (the file defining `Bp5AnalyticStressSource`, the unit under test) nor `fault/fault_geometry_safs_templated.inl` (the `ComputeParams<StressSource>` body the test instantiates). So editing the functor or the projection body does **not** trigger a rebuild of the test that validates them — `make seas_test_spatial_seas_bp5_analytic` reports "up to date" against a stale object and the test passes on old code. (This is why the implementer had to `touch` files manually mid-session.)

**Trigger:** edit `spatial/code/spatial_stress.hpp` (e.g., introduce a sign bug in `Bp5AnalyticStressSource::Evaluate`), then `make seas_test_spatial_seas_bp5_analytic` → no rebuild, test still "passes".

**Actual behavior:** stale binary, regression undetected.

**Expected behavior:** edits to the functor / projection body rebuild the test.

**Suggested fix:** add the real prerequisites to the obj rule.
```diff
 $(TEST_SPATIAL_SEAS_BP5_ANALYTIC_OBJ): %.o: $(SRC)%.cpp \
                                       $(SEAS_HEADERS) \
+                                      spatial/code/spatial_stress.hpp \
+                                      fault/fault_geometry_safs_templated.inl \
+                                      config/bp5_params.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
```
(The pre-existing faultgeom_parity rule has the same latent gap for `spatial_friction.hpp`; out of scope for this fix but worth a follow-up.)

**Test case:**
```bash
# test_R402_rebuild_on_functor_edit:
touch spatial/code/spatial_stress.hpp
make seas_test_spatial_seas_bp5_analytic   # MUST recompile the .o, not "up to date"
```

---

### [R-403] LOW [BUG] domain/elasticity_operator.hpp:GetModel() — null-deref in heterogeneous mode

**Category:** BUG (latent) / ASSUMPTION

**Description:**
The heterogeneous ctor sets `model_ = nullptr`. `GetModel()` returns `*model_` unconditionally, so calling it on a heterogeneous operator is undefined behavior (null-ref). It is documented on the ctor but not guarded. No current heterogeneous caller invokes `GetModel()`, so LOW.

**Suggested fix:**
```diff
-   const ConstitutiveModel &GetModel() const { return *model_; }
+   const ConstitutiveModel &GetModel() const
+   {
+      MFEM_VERIFY(model_ != nullptr,
+                  "GetModel(): no scalar ConstitutiveModel in the heterogeneous "
+                  "(Coefficient) operator mode — use the coefficients instead.");
+      return *model_;
+   }
```

**Test case:** (LOW — optional) build a heterogeneous operator and assert `GetModel()` aborts rather than segfaults.

---

### [R-404] LOW [ASSUMPTION] Stage-1 report mischaracterized the λ,μ-scalar site scope (affects Stage 1b completeness)

**Category:** ASSUMPTION / DEVIATION (scope accuracy)

**Description:**
The implementer's Stage-1 report described the `lambda_val_/mu_val_` generalization scope as "the ~20 hand-rolled sites in `ComputeTractionImpl` + the BR2 traction-quadrature path." In fact the 16 sites span **four** functions: `AssembleSlipContributionIP`, `AssembleSlipContributionIPShared`, `AssembleDirichletLoading` (all RHS-assembly, 12 sites) and `ComputeTractionImpl` (4 sites). If Stage 1b generalizes only the traction-recovery sites and misses the three RHS-assembly functions, the heterogeneous K·u=b right-hand side uses NaN λ,μ → wrong displacement → wrong everything (and the 2-layer slab test would expose it only if it covers Dirichlet loading + slip RHS, not traction alone).

**Suggested fix:** Stage 1b checklist must enumerate and cover all four functions:
- `assembly.inl::AssembleSlipContributionIP` (650, 651, 726, 727, 752, 753)
- `assembly.inl::AssembleSlipContributionIPShared` (977, 978, 1039, 1040)
- `assembly.inl::AssembleDirichletLoading` (1240, 1241, 1291, 1292, 1506, 1507, 1559, 1560, 1604, 1605, 1801, 1802, 1853, 1854)
- `traction.inl::ComputeTractionImpl` (1189, 1192, 1206, 1208, 1661, 1664, 1675, 1677)

**Test case:** the Stage-3 slab test must exercise a Dirichlet-loaded heterogeneous solve (not only fault-slip traction) so a missed `AssembleDirichletLoading` site is caught.

---

### [R-405] LOW [POSSIBLE] domain/elasticity_operator.hpp:MaterialCoefficient — SetTime not forwarded to the external coefficient

**Category:** QUALITY / ASSUMPTION (forward-looking)

**Description:**
`MaterialCoefficient` overrides `Eval` to delegate to `external_`, but does not override `SetTime`. If a time-dependent external `Coefficient` were ever supplied, `mfem` machinery calling `SetTime(t)` on the `MaterialCoefficient` would not propagate to `external_`, so its `Eval` would use a stale time. Material λ,μ are time-independent, so this is currently harmless — flagged because the delegation is incomplete.

**Suggested fix (defensive):**
```diff
    void SetExternal(Coefficient *c) { MFEM_VERIFY(c != nullptr, "..."); external_ = c; }
+
+   void SetTime(real_t t) override
+   {
+      Coefficient::SetTime(t);
+      if (external_) { external_->SetTime(t); }
+   }
```

---

## Summary
- Critical issues: **0** (the Phase 3c functor + test and the Stage-1a infrastructure are correct and verified; const-path is bit-for-bit by construction and 461+46 elasticity tests pass).
- Moderate issues: **2** — R-401 (heterogeneous ctor → silent NaN RHS/traction; latent behind the driver abort, CRITICAL if Stage 2 precedes Stage 1b), R-402 (Makefile stale-test risk: functor edits don't rebuild the test).
- Low issues: **3** — R-403 (GetModel() null-deref in hetero mode), R-404 (Stage-1b scope must cover all four λ,μ functions, not just traction), R-405 (MaterialCoefficient SetTime not forwarded).
- Plan compliance:
  - **Phase 3c — FULL** for its verification-only core: `Bp5AnalyticStressSource` + the 1e-10 parity test (np1/np4) close the R-007 `tau_pre_` gap; R-002 honored (no `StressSourceKind` enum extension); the config flag + `ComputeParams` dispatch correctly deferred to Phase 5 §B (the driver does not call `ComputeParams` yet).
  - **Phase 2b — PARTIAL (Stage 1a only, as intended):** operator infra + const-path bit-for-bit done; Stages 1b (λ,μ site generalization), 2 (driver wiring), 3 (slab test) pending.
- Verdict: **PASS WITH FIXES** — apply R-402 now (build correctness). Apply R-401 (guard) before any further work, and treat R-401's removal as the explicit Stage-1b deliverable so Stage 2 cannot precede Stage 1b. R-403/R-404/R-405 are hardening / scoping.

## Unreviewed Areas
- Stage 1b/2/3 code (not yet written) — out of scope; R-401/R-404 are forward constraints on them.
- The Phase 5 `bp5_analytic` driver dispatch (the functor's production caller) — out of scope; note for Phase 5: it must construct `Bp5AnalyticStressSource` from the SAME `BP5Params` the geometry used and pass `P_p=0` (bp5.sigma_n is already effective), and replicate the ±basis-constancy guard the test uses.
- Full `make test` aggregate not run; verified via the targeted elasticity/spatial tests + driver-object recompiles.
```
