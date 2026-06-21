# Fix Report: REVIEW_phase4_code_2026-06-04.md — 2026-06-04

## Summary
- Findings addressed: **3 of 5 fixed locally** (R-402, R-403, R-404); **2 are Frontera Stage-C deliverables** (R-401, R-405), not local code fixes.
- Files modified: `domain/elasticity_operator_setup.inl`, `domain/elasticity_operator.hpp`, `domain/elasticity_operator_assembly.inl`, `tests/unit/test_spatial_seas_iterative_vs_direct.cpp`.
- Tests added: 1 assertion (R-402 preconditioner-reuse check, verifying a Phase-4 acceptance criterion).
- Test suite: **PASS** — iterative_vs_direct np1 5/5 & np4 20/20; elasticity_operator 461/461, elasticity_br2 46/46, het_ctor 9/9, slab 27/27, faultgeom_parity np4 279408/0, bp5_analytic np4 48/0; BP5 driver object compiles.

## Changes Made
1. **R-402 (LOW) — preconditioner-reuse not test-verified** → added `mutable int num_stiffness_assemblies_` to `ElasticityDomainOperator`, incremented at the end of `AssembleStiffness`, exposed via `int NumStiffnessAssemblies() const`. Extended `test_spatial_seas_iterative_vs_direct` with a **2nd `Solve` at a different loading time** and an assertion `NumStiffnessAssemblies() == 1`. This now **verifies the plan's acceptance criterion** ("AMG setup runs once over an N-solve sequence") with an automated check, not just by inspection. Confirmed: counter == 1 after two solves (np=1 and np=4).
2. **R-403 (LOW) — serial CG path ignored ksp knobs** → the serial `SetupSolver` branch now uses `ksp_rtol_`/`ksp_atol_`/`ksp_maxit_`/`amg_print_level_` (was hardcoded `1e-12`/`0`/`10000`/`0`). **Verified the loosened default rel-tol (1e-12 → 1e-10) does not perturb the serial tests**: elasticity_operator 461/461, elasticity_br2 46/46, het_ctor 9/9 (its op-vs-op bit-for-bit comparison holds — both operators use the same tol), slab 27/27. (The reviewer's exact diff was applied; the rel-tol change is safe because every serial-test assertion threshold is ≥ 1e-6.)
3. **R-404 (LOW) — GlobalRelL2 silent partial compare** → added `MFEM_VERIFY(a.Size() == b.Size(), …)` and loop to `a.Size()` (was `std::min`). The test's separate size checks already guarded it; the helper is now safe in isolation.

## Verification
- [x] **R-402** — counter added + `NumStiffnessAssemblies()==1` asserted after 2 solves; iterative_vs_direct np1 5/5, np4 20/20.
- [x] **R-403** — serial path honors the knobs; serial tests all pass at the 1e-10 default (461/46/9/27).
- [x] **R-404** — size-mismatch guard added; test still passes.
- [ ] **R-401 (MODERATE)** — **not a local fix.** The CG-AMG-vs-MUMPS-**direct** comparison cannot run locally (MUMPS bus-errors in `dmumps_scatter_dist_rhs_` on the tiny fixture, np=1 and np=4; it is stable only on the production BP5 mesh, which is Frontera-only per [[no-local-reproducer]]). It is the Phase-4 **Frontera Stage-C** deliverable (run `cg_amg` and `mumps` on a real BP5 mesh, assert `||du||/||u|| < 1e-7`). The local test documents this deviation and uses GMRES_AMG as a sound iterative cross-reference.
- [ ] **R-405 (LOW)** — **not a local fix.** The Caliper perfgraph runtime-report + 2-resolution AMG-iteration scaling run on production meshes → Frontera Stage-C. The `MFEM_PERF_SCOPE` annotations are in place (no-op without Caliper).

## No-regression evidence
- elasticity_operator 461/461, elasticity_br2 46/46 (recompiled; R-403 serial-tol change harmless).
- het_ctor 9/9, slab 27/27 (serial Solve path, R-403).
- iterative_vs_direct np1 5/5, np4 20/20 (R-402 reuse + R-404 guard).
- faultgeom_parity np4 279408/0, bp5_analytic np4 48/0 (parallel construct; counter member additive).
- `drivers/seas_driver.o` (BP5, MUMPS) compiles clean.
- No debug prints / commented-out code / stray TODOs.

## Unresolved Findings (Frontera Stage-C deliverable, by design)
- **R-401** (CG-AMG vs MUMPS-direct on a BP5 mesh) and **R-405** (perfgraph runtime-report + multi-resolution AMG-iteration scaling). Both require production-mesh Frontera runs; they are the remaining Phase-4 deliverable, tracked for the Stage-C sbatch. Not blocking Phase 5.

## Ready for Re-Review: YES
