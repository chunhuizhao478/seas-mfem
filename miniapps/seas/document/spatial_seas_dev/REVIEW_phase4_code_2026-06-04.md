# Code Review: Phase 4 — scalable iterative solver (CG/GMRES + BoomerAMG) — 2026-06-04

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` §Phase 4.
- Files reviewed:
  - `domain/domain_config.hpp` (ksp_* knobs), `domain/elasticity_operator.hpp` (members + ctor extraction),
  - `domain/elasticity_operator_assembly.inl` (CG_AMG/GMRES_AMG dispatch), `domain/elasticity_operator_traction.inl` (residual check),
  - `drivers/spatial_seas_driver.cpp` (plumbing), `tests/unit/test_spatial_seas_iterative_vs_direct.cpp`, `Makefile`.
- Domain context: CLAUDE.md (extreme-care domain/, no silent fallback), the plan (R-006), project memory ([[no-local-reproducer]] — no production-mesh runs locally), the implementer completion report.
- Method: three adversarial passes; build + run exercised (iterative_vs_direct np1/np4 ≈1e-11; elasticity 461/46; faultgeom-parity np4; bp5-analytic np4; het ctor/slab). Verified the reuse guard and the serial-path hardcoding directly.

## Findings

### [R-401] MODERATE [DEVIATION] tests/unit/test_spatial_seas_iterative_vs_direct.cpp — local reference is GMRES_AMG (iterative), not MUMPS (direct)

**Category:** DEVIATION (plan acceptance criterion)

**Description:**
The plan's Phase-4 acceptance criterion is literally *"`||u_CGAMG - u_MUMPS||/||u_MUMPS|| < 1e-7` and the recovered fault traction agrees `< 1e-7`."* The implemented local test compares CG_AMG against **GMRES_AMG** (another AMG-preconditioned Krylov solver), not against the MUMPS **direct** solver. The implementer documented this: MUMPS bus-errors in `dmumps_scatter_dist_rhs_` on the tiny Cartesian fixture (reproduced np=1 and np=4), and the production BP5 mesh — where MUMPS is stable — is Frontera-only ([[no-local-reproducer]]). The substitution is mathematically sound (both Krylov methods converge to the unique `K⁻¹b`; the differing preconditioners — CG has `SetElasticityOptions`, GMRES does not — make a shared-error coincidence unlikely), so it validates CG_AMG convergence. BUT it does NOT validate against an independent **direct** factorization, so a systematic Krylov/AMG-shared error (or a wrong assembled K shared by both) would pass. The CG-AMG-vs-MUMPS check is therefore still **owed**.

**Trigger:** the local test never exercises a direct solve; the direct comparison is deferred.

**Actual behavior:** local equivalence is iterative-vs-iterative.

**Expected behavior:** the plan's iterative-vs-**direct** comparison runs somewhere before Phase 4 is considered closed.

**Suggested fix:** add the CG_AMG-vs-MUMPS comparison to the **Phase-4 Frontera Stage-C job** (the perfgraph/scaling sbatch), on a real BP5 mesh where MUMPS is stable, and record the relative-error result. Keep the local GMRES_AMG cross-check as the fast unit test. Concretely, the Stage-C driver/sbatch should run the same `--config` with `[solver].type="cg_amg"` and `"mumps"` and diff the station/traction output `< 1e-7`.
```diff
  # Stage-C Frontera sbatch (to be generated): on bp5_tandem_exact.msh,
  #   run cg_amg and mumps, assert ||u_cg - u_mumps||/||u_mumps|| < 1e-7.
```

**Test case:**
```
// test_R401_cg_vs_mumps_on_bp5 (FRONTERA, real mesh):
//   build op(CG_AMG) and op(MUMPS) on bp5_tandem_exact.msh, same slip+time,
//   Solve both, assert global ||du||/||u|| < 1e-7 and ||dT||/||T|| < 1e-7.
//   (MUMPS is stable on this mesh; it crashes only on the tiny local fixture.)
```

---

### [R-402] LOW [POSSIBLE] domain/* — preconditioner-reuse ("AMG setup runs once") is structurally correct but not verified by a test

**Category:** ASSUMPTION / coverage

**Description:**
The plan's acceptance criterion "AMG setup runs **once** over an N-solve sequence (Caliper region count == 1 / a setup-time counter increments once)" is satisfied **structurally** — `Solve` guards `if (!stiffness_assembled_) { AssembleStiffness(); }` (traction.inl:494) and `AssembleStiffness` sets `stiffness_assembled_=true` at its end, so the solver/AMG dispatch runs only on the first solve. But there is **no automated test** asserting that a 2nd/Nth `Solve` does not rebuild the AMG (the criterion explicitly asks for a region-count/counter check). A future refactor that drops the call-site guard would silently regress to per-solve AMG rebuilds with no test failure.

**Suggested fix:** add a reuse test — solve twice and assert the second solve does not re-trigger assembly (e.g., add a `mutable int amg_setup_count_` incremented in the CG_AMG/GMRES_AMG branch and assert it == 1 after two solves):
```cpp
// test_R402_amg_setup_once:
//   op(CG_AMG); op.Solve(t1, slip, u1); op.Solve(t2, slip, u2);
//   assert op.NumStiffnessAssemblies() == 1 (add a small accessor/counter).
```

**Test case:** as above (requires a tiny setup counter/accessor on the operator).

---

### [R-403] LOW [DEVIATION] domain/elasticity_operator_setup.inl:SetupSolver — serial CG path ignores ksp_rtol/ksp_maxit

**Category:** DEVIATION / inconsistency

**Description:**
Phase 4 threaded `ksp_rtol_`/`ksp_atol_`/`ksp_maxit_`/`amg_print_level_` into the **parallel** dispatch only. The serial `SetupSolver` branch still hardcodes `cg->SetRelTol(1e-12); SetMaxIter(10000); SetPrintLevel(0)`, so a serial caller that sets `DomainConfig.ksp_rtol` gets it silently ignored. The residual check (which IS parameterized: `10*ksp_rtol_`) then uses a threshold unrelated to the serial solver's actual tolerance. Harmless today (the QD/SAF target is parallel; serial is a test-only path, and the serial CG converges to 1e-12 < the check threshold), but inconsistent.

**Suggested fix:** thread the same knobs into the serial branch:
```diff
       else
       {
          auto *cg = new CGSolver();
-         cg->SetRelTol(1e-12);
-         cg->SetAbsTol(0.0);
-         cg->SetMaxIter(10000);
-         cg->SetPrintLevel(0);
+         cg->SetRelTol(ksp_rtol_);
+         cg->SetAbsTol(ksp_atol_);
+         cg->SetMaxIter(ksp_maxit_);
+         cg->SetPrintLevel(amg_print_level_);
          solver_.reset(cg);
       }
```
(Defaults reproduce 1e-10 rel-tol; the serial tests do not set ksp_rtol so behavior is unchanged except the rel-tol default 1e-12 → 1e-10, still far tighter than any acceptance tol — confirm the serial elasticity tests still pass.)

---

### [R-404] LOW tests/unit/test_spatial_seas_iterative_vs_direct.cpp:GlobalRelL2 — silent partial compare on size mismatch

**Category:** EDGE_CASE

**Description:**
`GlobalRelL2` loops to `std::min(a.Size(), b.Size())`, so if the two vectors ever differed in size it would silently compare only the common prefix and report a misleadingly small error. The test does guard with separate `GetVSize`/traction-size `Check`s, so this cannot currently mislead — but the helper itself is unsafe in isolation.

**Suggested fix:**
```diff
 real_t GlobalRelL2(const Vector &a, const Vector &b, MPI_Comm comm)
 {
+   MFEM_VERIFY(a.Size() == b.Size(), "GlobalRelL2: size mismatch "
+               << a.Size() << " vs " << b.Size());
    real_t ld = 0.0, lb = 0.0;
-   const int n = std::min(a.Size(), b.Size());
+   const int n = a.Size();
```

---

### [R-405] LOW [tracked] Phase 4 Stage C (perfgraph + multi-resolution scaling) not implemented

**Category:** DEVIATION (deferred)

**Description:**
Two Phase-4 acceptance criteria — the Caliper **runtime-report comparing cg_amg/gmres_amg/mumps** (AMG-setup time, mean per-solve, iteration count) and the **AMG-iteration-count-bounded-across-two-BP5-resolutions** scaling smoke — are not implemented; they run on production meshes (Frontera, per project memory). The Caliper `MFEM_PERF_SCOPE` annotations are in place (no-op without Caliper), so the instrumentation is ready, but the runs are owed.

**Suggested fix:** generate the Frontera Stage-C sbatch (build with `MFEM_USE_CALIPER=YES`, run `cg_amg`/`gmres_amg`/`mumps` on two BP5 resolutions, emit the runtime-report + iteration counts, and fold the R-401 CG-vs-MUMPS check into the same job). Track as the remaining Phase-4 deliverable.

---

## Summary
- Critical issues: **0** — the dispatch is validated (CG_AMG converges to GMRES_AMG ≈1e-11 at np=1/np=4); preconditioner reuse is correctly guarded (`if (!stiffness_assembled_)`); no BP5 regression (BP5 = MUMPS, untouched).
- Moderate issues: **1** — R-401 (local reference is GMRES_AMG not MUMPS-direct; the direct comparison is owed on Frontera).
- Low issues: **4** — R-402 (reuse not test-verified), R-403 (serial path ignores ksp knobs), R-404 (GlobalRelL2 size-mismatch), R-405 (Stage-C perfgraph/scaling deferred).
- Plan compliance: **PARTIAL** — the local-implementable Phase-4 work (dispatch, config, residual check, driver plumbing, Caliper annotations, solver-equivalence test, preconditioner reuse) is complete and verified; 2 of 4 acceptance criteria (perfgraph + multi-resolution scaling) and the literal CG-vs-MUMPS-direct comparison are Frontera-deferred (Stage C), consistent with project memory.
- Verdict: **PASS WITH FIXES** — no blocking bug for proceeding to Phase 5. Apply R-403/R-404 cheaply now; R-401/R-402/R-405 are the Frontera Stage-C deliverable + a reuse test. The R-006 `SetElasticityOptions`-on-DG decision is empirically validated (CG converged).

## Unreviewed Areas
- The MUMPS / SuperLU / STRUMPACK dispatch branches (untouched by Phase 4; out of scope).
- Behavior on the production BP5 mesh (Frontera; the local fixture is the small Cartesian fault mesh).
- Caliper-enabled build (`MFEM_USE_CALIPER=YES`) — the macros are no-ops in the functional build tested here.
```
