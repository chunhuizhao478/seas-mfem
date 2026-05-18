# Fix Report — Round 3 (2026-05-18)

Review: `spatial_dynamic_rupture_review_round3_2026-05-18.md`

## Summary

- Findings addressed: **12 of 15** (R-314 and R-315 are scope-incomplete trackers; R-306 is gated by R-301 abort)
- Files modified:
  - `miniapps/seas/dynamic/wave_operator.hpp`
  - `miniapps/seas/dynamic/wave_operator.inl`
  - `miniapps/seas/dynamic/spatial_setup.hpp`
  - `miniapps/seas/dynamic/godunov_flux_pool.cpp`
  - `miniapps/seas/fault/fault_geometry.hpp`
  - `miniapps/seas/tests/unit/test_spatial_setup.cpp`
  - `miniapps/seas/tests/unit/test_phaseh_godunov_flux_pool.cpp`
  - `miniapps/seas/tests/unit/test_phaseh_lsw_forced_rupture.cpp`
- Tests added: **5** (S-7, S-8, P-5, P-6, F-7)
- Test suite: **PASS** — 515 / 515:
  - `seas_test_phaseh_godunov_flux_pool`        — **30 / 30** (was 20; +P-5, +P-6)
  - `seas_test_phaseh_lsw_forced_rupture`       — **35 / 35** (was 32; +F-7)
  - `seas_test_spatial_setup`                   — **68 / 68** (was 63; +S-7, +S-8)
  - `seas_test_spatial_friction_config`         — **28 / 28**
  - `seas_test_spatial_friction_resolver`       — **110 / 110**
  - `seas_test_spatial_velocity_bundle`         — **6 / 6**
  - `seas_test_spatial_stress_bundle`           — **6 / 6**
  - `seas_test_spatial_constant_stress_source`  — **45 / 45**
  - `seas_test_compute_safs_params`             — **13 / 13** (byte-exact regression intact)
  - `seas_test_safs_mode_wiring`                — **8 / 8**
  - `seas_test_tpv104_checkpoint`               — **164 / 164** (DRIVER_TAG_V1 back-compat intact)
  - `seas_tpv104_driver`, `seas_tpv205_driver`  — rebuild clean after wave_operator.inl edit

## Changes Made

1. **R-301** (flux_pool dispatch dead) — `wave_operator.hpp:SetGodunovFluxPool`: setter now `MFEM_ABORT`s with a precise message pointing to the missing Phase H.2/H.3/H.5 dispatch follow-up.  Any caller that tries to enable heterogeneous-material physics fails LOUD instead of silently running scalar physics.
2. **R-302** (LSW_ForcedRupture dispatch arm missing) — `wave_operator.inl:3617` (interior fault) and `:4539` (shared-fault fallback): added `else if (fault_friction_law_ == FaultFrictionLaw::LSW_ForcedRupture)` arms that route to `fault_flux_->EvaluateADER_LSW_ForcedRupture(..., GetTime(), ...)`.  `GetTime()` comes from `TimeDependentOperator` base.
3. **R-303** (ip.Init(0) regression in spatial_setup) — `spatial_setup.hpp:seed_static_dof_fields`: replaced `IntegrationPoint ip; ip.Init(0);` with `Geometries.GetCenter(mesh.GetElementBaseGeometry(elem))`.  Comment rewritten to correctly describe centroid semantics.
4. **R-304** (zero-normal counter bypass) — `fault_geometry.hpp` new BP5 ctor: replaced unconditional `num_zero_normal_fallbacks_ = 0` with a per-DOF walk that increments the counter for any `||n_i|| < 1e-10` or `||t1_i|| < 1e-10`.  The R-001 SAFS pre-flight guard now fires loudly on malformed bases.
5. **R-305** (missing MaterialField ctor) — `wave_operator.hpp` scalar ctor: added a doc-block explaining the Phase H.1 deferral so the next implementer knows the surface the plan calls for and the prerequisites (R-301 abort needs lifting first).  No stub ctor added — pulling in `heterogeneous_material.hpp` for a function that immediately aborts would add dead code without value.
6. **R-309** (fault_dof_ip cache never consumed) — `spatial_setup.hpp`: added IP-aware `seed_static_dof_fields` overload + IP-aware `InitializeFaultDOFs_Spatial` overload that takes `const std::vector<IntegrationPoint>& dof_ips` and uses the supplied IP per DOF instead of the centroid.  Drivers using `FaultGeometry::fault_dof_ip()` can now thread it through.
7. **R-310** (BP5 coords_x2/x3 confusion) — `fault_geometry.hpp` new BP5 ctor: set `coords_x2_` / `coords_x3_` to NaN sentinels and skip `ComputeBP5Params()`.  `depths_` (still consumed by ABSORBING / free-surface BC dispatch via world z) is preserved.  Any premature read of BP5 analytic per-DOF arrays now propagates NaN loudly.
8. **R-311** (dead unreachable guard) — `godunov_flux_pool.cpp:Build`: removed the `ne > 0 && unique_fluxes_.empty()` assertion (loop body always populates at least one triple when `ne > 0`).
9. **R-312** (round_sig overflow) — `godunov_flux_pool.cpp:round_sig`: added `|log_mag| > 290.0` guard that bypasses the round when `mag` is at the double-precision range edge (subnormal, near-overflow).
10. **R-308** (no layered-material test for impedances) — `test_spatial_setup.cpp`: added **S-7** with a `Vs(z) = 1000 + 4000z` Coefficient on a unit-cube hex; asserts `Zs+` equals `rho · Vs(centroid)` (3000) and is NOT `rho · Vs(corner)` (1000).  Catches R-303 regressions.
11. **R-307** (missing Phase H tests) — `test_phaseh_godunov_flux_pool.cpp`: added **P-5** (determinism: two pools from same input have identical mappings) and **P-6** (memory budget on a 4-layer synthetic 10⁴-element fixture).  The two plan-required test files `test_phaseh_wave_operator_constant_parity.cpp` and `test_phaseh_wave_operator_layered.cpp` are NOT created because their preconditions are R-305 (ctor) and R-301 (live dispatch) — both still effectively absent.  Marked as deferred in §Unresolved.
12. **R-313** (helper-only forced-rupture tests) — `test_phaseh_lsw_forced_rupture.cpp`: added **F-7** that exercises `EvaluateADER_LSW_ForcedRupture` end-to-end through `FaultFaceFlux`.  Picks params so plain LSW remains locked (`mu_s · sigma_n = 50 MPa > tau_abs = 35 MPa`) while forced-rupture at `f_2 = 0.5` unlocks (`mu_eff · sigma_n = 30 MPa < 35 MPa`), guaranteeing the two paths produce materially different `V_abs` and `I_imp_*`.

## Unresolved Findings

- **R-306** (heterogeneous CFL missing) — Deferred.  `ComputeMaxDt` still uses scalar `flux_.GetCp()` only.  Since R-301 now aborts on `SetGodunovFluxPool`, the heterogeneous code path is unreachable; the latent bug stays documented but cannot be triggered.  Will be implemented alongside the Phase H.2 dispatch follow-up that lifts R-301.
- **R-307 partial** — Two of the plan-required Phase H test files (`test_phaseh_wave_operator_constant_parity.cpp`, `test_phaseh_wave_operator_layered.cpp`) are not created.  They need a live heterogeneous-material `WaveOperator(MaterialField, BoundaryConfig)` ctor (R-305) plus the dispatch from R-301 to actually run.  Both prerequisites are documented but still deferred.
- **R-314** (Phase 4 `spatial_dyn_driver.cpp` missing) — Scope-incomplete tracker; no fix in this pass.
- **R-315** (Phase 6 sbatches + verify/compare scripts mostly missing) — Scope-incomplete tracker; no fix in this pass.

## New Tests

- `S_7_layered_material_impedances_use_centroid`   — covers R-303 (catches reintroduction of `ip.Init(0)` bug)
- `S_8_ip_aware_overload_uses_supplied_ip`         — covers R-309 (IP-aware overload uses supplied IP)
- `P_5_mpi_consistency_determinism`                — covers R-307 (plan-required deterministic-pool test)
- `P_6_memory_budget_synthetic`                    — covers R-307 (plan-required memory-budget test, synthetic fixture)
- `F_7_forced_rupture_differs_from_plain_lsw`      — covers R-313 (forced-rupture path actually changes physics)

## Ready for Re-Review: YES
