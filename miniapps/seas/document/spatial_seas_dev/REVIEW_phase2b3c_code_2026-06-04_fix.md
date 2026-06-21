# Fix Report: REVIEW_phase2b3c_code_2026-06-04.md — 2026-06-04

## Summary
- Findings addressed: **5 of 5** (R-401, R-402, R-403, R-405 fixed with code; R-404 is a forward scope constraint with no Stage-1a code change — recorded for Stage 1b).
- Files modified: `domain/elasticity_operator.hpp`, `Makefile`.
- Tests added: 0 new unit-test files (rationale per finding below). All existing tests re-run and pass.
- Test suite: **PASS** — elasticity_operator 461/461, elasticity_br2 46/46, faultgeom_parity np1 279372/0 & np4 279408/0, bp5_analytic np1 12/0 & np4 48/0; `seas_driver.o` + `spatial_seas_driver.o` recompile clean.

## Changes Made
1. **R-401 (MODERATE)** — heterogeneous ctor silent-NaN → added an unconditional `MFEM_ABORT` guard in the `Coefficient&` ctor (`domain/elasticity_operator.hpp`), after the member setup and before `InitOperator()`. It names the four functions that still read scalar `lambda_val_/mu_val_` and is marked **"REMOVE in Stage 1b"**. This makes a half-wired heterogeneous operator fail loudly instead of silently producing NaN RHS/traction, and enforces the ordering constraint that Stage 2 (driver wiring) cannot precede Stage 1b. **Root-cause note:** the true fix is Stage 1b (generalize the four functions to `coeff.Eval(T,ip)`); the guard is the correct *interim* per the review, and its removal is the explicit Stage-1b deliverable.
2. **R-402 (MODERATE)** — Makefile stale-test risk → added `spatial/code/spatial_stress.hpp`, `fault/fault_geometry_safs_templated.inl`, and `config/bp5_params.hpp` as explicit prerequisites of `$(TEST_SPATIAL_SEAS_BP5_ANALYTIC_OBJ)` (they are not in `$(SEAS_HEADERS)`). Verified: `touch spatial/code/spatial_stress.hpp` now forces a recompile (previously "up to date").
3. **R-403 (LOW)** — `GetModel()` null-deref → added `MFEM_VERIFY(model_ != nullptr, …)` before `return *model_`, with a message pointing heterogeneous-mode callers to the coefficients. Constant-ctor callers (e.g. `test_boundary_config_operator`) are unaffected (`model_` valid).
4. **R-405 (LOW)** — `MaterialCoefficient::SetTime` not forwarded → added a `SetTime(real_t)` override that calls `Coefficient::SetTime` and forwards to `external_` when set. Defensive (material λ,μ are time-independent today).

## Verification
- [x] **R-401** — guard present (`grep "Stage 1b"` → 4 hits incl. the `MFEM_ABORT`). The header + all dependents recompile clean. *No process-level unit test added:* `MFEM_ABORT` terminates the process (no `MFEM_USE_EXCEPTIONS` in this build), so an in-process "expect abort" unit test is impractical; the guard is unconditional in the ctor body (verified by inspection + compilation). The proper **functional** test — construct via the `Coefficient&` ctor, `Solve` + `ComputeTraction`, assert `isfinite` and bit-for-bit vs the constant ctor — is the Stage-1b deliverable (it requires the guard to be removed first), exactly as the review's R-401 test case states.
- [x] **R-402** — verified directly: `touch spatial/code/spatial_stress.hpp && make seas_test_spatial_seas_bp5_analytic` recompiles the object (0 "up to date" lines on the first post-touch run; up-to-date on the immediate second run). This is a Makefile-dependency behavior, verified by the rebuild check rather than a unit test.
- [x] **R-403** — guard present; `seas_test_elasticity_operator` (which calls the constant-ctor path incl. `GetShearModulus`/model access) passes 461/461.
- [ ] **R-404** — *not a Stage-1a code change.* It corrects the Stage-1b scope: the 16 `lambda_val_/mu_val_` sites live in `AssembleSlipContributionIP`, `AssembleSlipContributionIPShared`, `AssembleDirichletLoading`, and `ComputeTractionImpl` (not just traction recovery). Recorded as the Stage-1b checklist + the requirement that the Stage-3 slab test exercise a Dirichlet-loaded heterogeneous solve. No edit applicable now.
- [x] **R-405** — override present (`grep "external_->SetTime"` → 1 hit); build clean.

## No-regression evidence
- `seas_test_elasticity_operator` 461/461; `seas_test_elasticity_br2` 46/46 (const-path correctness guard, unaffected by the additive `SetTime`/`GetModel` changes).
- `seas_test_spatial_seas_faultgeom_parity` np1 279372/0, np4 279408/0.
- `seas_test_spatial_seas_bp5_analytic` np1 12/0, np4 48/0.
- BP5 driver (`drivers/seas_driver.o`) and QD driver (`drivers/spatial_seas_driver.o`) recompile clean; dynamic driver unaffected (does not include `elasticity_operator.hpp`).
- No debug prints / commented-out code / stray TODOs introduced. The R-401 `MFEM_ABORT` is a documented interim safety guard (marked for Stage-1b removal), not a leftover.

## Unresolved Findings
- **R-404** — forward scope constraint only; enforced when Stage 1b is implemented (no current code to change).

## Ready for Re-Review: YES
