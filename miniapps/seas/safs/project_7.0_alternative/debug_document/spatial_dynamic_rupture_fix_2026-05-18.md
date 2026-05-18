# Fix Report — spatial_dynamic_rupture rev-3 (2026-05-18)

Review: `spatial_dynamic_rupture_review_2026-05-18.md`

## Summary

- Findings addressed: **11 of 12** (R-001 is informational; the other 11 are fixed)
- Files modified:
  - `miniapps/seas/spatial/code/spatial_friction.cpp`
  - `miniapps/seas/spatial/code/spatial_velocity.cpp`
  - `miniapps/seas/spatial/code/spatial_stress.cpp`
  - `miniapps/seas/fault/fault_geometry_safs_templated.inl`
  - `miniapps/seas/io/field_coefficient.hpp`
  - `miniapps/seas/tests/unit/test_spatial_friction_config.cpp`
  - `miniapps/seas/tests/unit/test_spatial_friction_resolver.cpp`
  - `miniapps/seas/tests/unit/test_spatial_stress_bundle.cpp`
- Tests added: **7** (4 for the parser, 2 regressions for R-002/R-003, 1 RS eta_auto; plus rewritten S-2/S-3/S-4)
- Test suite: **PASS** — all 4 phase-0..3b binaries plus the existing `seas_test_compute_safs_params` (13/13 byte-exact regression):
  - `seas_test_spatial_friction_config`        — **24 / 24**
  - `seas_test_spatial_friction_resolver`      — **110 / 110**
  - `seas_test_spatial_velocity_bundle`        — **6 / 6**
  - `seas_test_spatial_stress_bundle`          — **6 / 6**
  - `seas_test_spatial_constant_stress_source` — **45 / 45**
  - `seas_test_compute_safs_params`            — **13 / 13** (unchanged)

## Changes Made

1. **R-002** (`ResolveForcedRupture` IP corner→centroid) — `spatial_friction.cpp` lines ~1019–1027: replaced `IntegrationPoint ip; ip.Init(0);` with `Geometries.GetCenter(mesh.GetElementBaseGeometry(e))`.  Comment updated to explain why the centroid is the correct stop-gap until Phase 5 wires `FaultGeometry::fault_dof_ip(i)`.
2. **R-003** (`resolve_rs_impl` eta_auto same bug) — `spatial_friction.cpp` lines ~895–910: same centroid replacement.
3. **R-004** (missing-block parser silence) — `spatial_friction.cpp:parse_root`: added an explicit `MFEM_VERIFY(root.contains("<block>"), ...)` for each of the seven mandatory top-level blocks before entering the per-block conditional.  Missing `[stress]`, `[time]`, `[material_constant_fallback]`, or `[numerics]` now aborts with a precise message.
4. **R-005** (fake stress-bundle tests) — `test_spatial_stress_bundle.cpp` rewritten: pulled in the same BP5-mesh `BuildFixture` pattern used by `test_compute_safs_params.cpp`; S-2 and S-3 now genuinely invoke `ApplyCsmStressSidecar(spec, geom)` inside a forked child; S-4 added — end-to-end invocation with a synthetic HDF5 sidecar that covers the BP5 mesh bbox (the real SAFS CSM sidecar is in UTM 11N coords and cannot be paired with the BP5 1000m mesh).
5. **R-006** (missing `dt_initial` validator) — `spatial_friction.cpp:parse_root` time block: added `MFEM_VERIFY(cfg.time.dt_initial > 0.0 || cfg.time.dt_initial == -1.0, ...)`.  Test added (`T-15`).
6. **R-007** (`toml_bool` no type guard) — `spatial_friction.cpp:toml_bool`: added `v.is_boolean()` check + `MFEM_ABORT` with the project-standard message format.  Test added (`T-16`).
7. **R-008** (`ApplyCsmStressSidecar` BP2-geom guard) — `spatial_stress.cpp:apply_csm_impl`: added an early `MFEM_VERIFY(geom.IsBP5(), ...)` before any potentially-uninitialised geom accessor.
8. **R-009** (misleading test name + missing eta_auto coverage) — `test_spatial_friction_resolver.cpp`: renamed `R_1_rs_defaults` → `R_1a_rs_defaults_eta_explicit` (matches what it actually tests); added a new `R_1b_rs_eta_auto_per_dof` that genuinely exercises the `eta = 0.5 * sqrt(mu*rho)` path via the serial-mesh overload + `MakeConstant`.
9. **R-010** (templated empty-N early-out drift) — `fault_geometry_safs_templated.inl`: removed `tau_pre_.SetSize(0)` from the early-out so the templated overload's post-condition matches the non-templated one exactly when `num_fault_dofs_ == 0`.
10. **R-011** (misleading L_3 comment) — `test_spatial_friction_resolver.cpp:L_3_last_match_wins`: comment changed from "rule B's mu_s = 0.9 violates the per-DOF inequality" → "rule B's mu_s = 0.9 SATISFIES the per-DOF inequality".
11. **R-012** (`AbortContainmentFailure` duplicated) — `io/field_coefficient.hpp`: moved `AbortContainmentFailure` and `AbortRangeFailure` from `private:` to `public:` with a justifying comment; `spatial_velocity.cpp:load_impl` now calls `FieldProjector::AbortContainmentFailure(...)` instead of duplicating the message text.

## Unresolved Findings

- **R-001** (Phases H, 4, 5, 6 missing) — **not fixed**: this is an informational scope-completeness report.  The plan rev-3 expects Phase H (heterogeneous `WaveOperator` + `GodunovFluxPool` + `LSW_ForcedRupture` dispatch), Phase 4 (`spatial_dyn_driver.cpp`), Phase 5 (`spatial_setup.hpp` + new `FaultGeometry` ctor + `DRIVER_TAG_V1` extension), and Phase 6 (sbatches + verify/compare scripts) to be implemented.  None of these were attempted in this fix pass — they are out of scope for "fix the bugs found in the review" and would each require their own implementation cycle.

## New / Renamed Tests

- `T-13_missing_stress_block_aborts`            — covers R-004 (parser)
- `T-14_missing_time_block_aborts`              — covers R-004 (parser)
- `T-15_dt_initial_zero_aborts`                 — covers R-006
- `T-16_use_pml_string_aborts`                  — covers R-007
- `F-3_T_uses_centroid_not_corner`              — covers R-002 with a 1-hex + Vs(z) Coefficient
- `R-9_eta_auto_uses_centroid`                  — covers R-003 with the same fixture
- `R-1a_rs_defaults_eta_explicit` (renamed) +
  `R-1b_rs_eta_auto_per_dof`                    — covers R-009
- `S-2_empty_path_aborts_through_apply`,
  `S-3_wrong_kind_aborts_through_apply`,
  `S-4_end_to_end_invocation_synthetic`         — covers R-005

## Ready for Re-Review: YES
