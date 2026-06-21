# Fix Report: REVIEW_phase2b_stages1b3_code_2026-06-04.md — 2026-06-04

## Summary
- Findings addressed: **5 of 5** (R-501 fixed at the driver level + scope clarified; R-502/503/504/505 fixed).
- Files modified: `drivers/spatial_seas_driver.cpp`, `domain/elasticity_operator.hpp`, `tests/unit/test_elasticity_heterogeneous_ctor.cpp`, `tests/unit/test_elasticity_heterogeneous_slab.cpp`, `Makefile`.
- Files created: `tests/unit/fault_mesh_fixture.hpp` (shared mesh fixture, R-503).
- Tests added: 1 assertion (R-504 structural-equality check). R-502 verified by a contradictory-config run.
- Test suite: **PASS** — het_ctor 9/9, slab 27/27, elasticity_operator 461/461, elasticity_br2 46/46; constant + depth-profile driver dry-runs exit 0; R-502 contradictory config correctly rejected (exit 1).

## Changes Made
1. **R-501 (MODERATE)** — `drivers/spatial_seas_driver.cpp` §4.4: added `&& cfg.material.kind != spatial::MaterialKind::Constant` to `sidecar_requested`, so an EXPLICIT `kind=constant` is honored and never triggers a sidecar load (regardless of the `use_sidecar=true` default). **Scope clarification (important):** the symptom originally reproduced — a constant config *without* `use_sidecar=false` aborting — is actually emitted by the **shared parser** (`spatial_friction.cpp:925-932`: when `use_sidecar` is true it requires `dataset_root`/`override_path`), at parse time, *before* the driver's material build. That parser check is a config-completeness validation (arguably correct) and is shared no-touch code, so it is NOT changed here. The driver-level gating fix correctly handles the genuine driver bug — `kind=constant` + `use_sidecar=true` + a *valid* sidecar path, where the pre-fix driver loaded the sidecar and ignored the explicit `kind=constant`; the post-fix driver uses the constant material. (Full reproduction of that path needs an HDF5 sidecar, so it is verified by gating-logic inspection + the no-regression runs, not an automated local test.) The practical convention stands: a constant QD config sets `[velocity].use_sidecar=false` (as the smoke config does).
2. **R-502 (LOW)** — same file: added an `MFEM_VERIFY` rejecting the contradictory `[material].kind="sidecar_hdf5"` + `[velocity].use_sidecar=false` (no silent constant fallback). **Verified:** such a config now aborts with "kind=\"sidecar_hdf5\" requires [velocity].use_sidecar=true" (exit 1).
3. **R-503 (LOW)** — created `tests/unit/fault_mesh_fixture.hpp` (`mfem::seas::test::{AddFaultBoundaryElements,CreateTestMesh3D}`, header-inline) and replaced the duplicated copies in `test_elasticity_heterogeneous_ctor.cpp` and `test_elasticity_heterogeneous_slab.cpp` with `#include` + `using test::CreateTestMesh3D;`. The pre-existing copy in `test_elasticity_operator.cpp` (the 461-test file) was left untouched (out of scope — "do not refactor"; it also carries a `CreateTestMesh3DTet` variant). Both obj rules gained `tests/unit/fault_mesh_fixture.hpp` as a prerequisite (rebuild-trigger correctness, the R-402 lesson).
4. **R-504 (LOW)** — `test_elasticity_heterogeneous_ctor.cpp`: added an explicit `GetFESpace().GetVSize()` equality check + a comment documenting that the bit-for-bit Vector comparisons rely on identical FE-space size/ordering (guaranteed by the deterministic `CreateTestMesh3D`). The test now runs 9 checks (was 8).
5. **R-505 (LOW)** — `domain/elasticity_operator.hpp`: changed the heterogeneous ctor's default `DGMethod method = DGMethod::BR2` → `DGMethod::IP` (the only supported heterogeneous method), so a caller relying on the default no longer hits the BR2 abort. No behavior change for current callers (driver + tests pass `IP` explicitly).

## Verification
- [x] **R-501** — gating fix applied; constant smoke (`use_sidecar=false`) → "material: constant" exit 0; depth-profile → "material: depth_profile_1d" exit 0; no regression. Driver-level `kind=constant`-honoring logic verified by inspection (full sidecar-path repro needs an HDF5 file). Parser-level symptom documented as shared-parser config validation (not changed).
- [x] **R-502** — contradictory `sidecar_hdf5` + `use_sidecar=false` config aborts with the new message (exit 1), verified by run.
- [x] **R-503** — shared `fault_mesh_fixture.hpp`; both tests rebuild from it and pass (het_ctor 9/9, slab 27/27). Obj rules depend on the header.
- [x] **R-504** — structural-size assert present; het_ctor 9/9.
- [x] **R-505** — default method now IP; elasticity_operator 461/461, elasticity_br2 46/46 (unaffected — they use the scalar/model ctors).

## No-regression evidence
- `seas_test_elasticity_operator` 461/461; `seas_test_elasticity_br2` 46/46 (recompiled against the changed header).
- `seas_test_elasticity_heterogeneous_ctor` 9/9; `seas_test_elasticity_heterogeneous_slab` 27/27 (ratio 2.05128).
- Constant + depth-profile QD driver dry-runs: exit 0, correct material banner.
- No debug prints / commented-out code / stray TODOs introduced.

## Unresolved Findings
- None outright unresolved. **R-501 partial scope note:** the most user-visible symptom (a constant config must set `use_sidecar=false`) is enforced by the shared parser, not the driver; making constant "just work" with the default `use_sidecar=true` would require relaxing the shared `spatial_friction.cpp` velocity validation (used by the dynamic driver too) — deferred as a shared-parser decision rather than a driver fix. The driver-level gating bug R-501 targets IS fixed.

## Ready for Re-Review: YES
