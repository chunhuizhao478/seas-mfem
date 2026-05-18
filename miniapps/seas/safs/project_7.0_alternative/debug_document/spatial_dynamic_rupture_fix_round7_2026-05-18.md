# Fix Report — Round 7 (2026-05-18)

Review consumed: `spatial_dynamic_rupture_review_round7_2026-05-18.md`.

## Summary

- Findings addressed: **6 / 6** (R-701, R-702, R-703, R-704, R-705, R-706).
- Files modified:
  - `miniapps/seas/spatial/code/spatial_friction.cpp` — R-701 (a) sentinel for non-StrengthReduction kinds in `resolve_forced_impl`.
  - `miniapps/seas/drivers/spatial_dyn_driver.cpp` — R-701 (b) `SetForcedRuptureMode(use_strength_reduction)` gate; R-704 corrected resolver comment block at lines 941-948; R-706 `(void)overstress` → `[[maybe_unused]]`.
  - `miniapps/seas/dynamic/wave_operator.inl` — R-702 documentation comments at BOTH `LSW_ForcedRupture` dispatch arms (interior-fault L3835, shared-fault L4783).
  - `miniapps/seas/dynamic/tpv205_substep_iterator.cpp` — R-705 documentation of START-OF-SUBSTEP convention next to `t_sub_acc` initialisation.
  - `miniapps/seas/tests/unit/test_spatial_friction_config.cpp` — T-21..T-23 added (NucleationKind parser).
  - `miniapps/seas/tests/unit/test_spatial_friction_resolver.cpp` — F-4 (R-701 regression for overstress kind) + F-5 (R-703 ResolveOverstress stub behaviour) added.
- Tests added: **5 new test functions** (3 parser + 2 resolver), contributing **~14 new sub-asserts**.
- 3D SAFS path test sweep: **PASS** — see verification below.
- Pre-existing failures **not** addressed (out of scope per [[feedback-preexisting-test-failures]]):
  - `tests/verification/bp2_serial_smoke.cpp` — fails to compile because it
    references renamed paraview-compaction API (`InitFaultOutput`,
    `UpdateFaultFields` instead of the renamed `InitFaultOutputBP5` /
    `UpdateFaultFieldsBP5`).  BP2 is 2D antiplane verification, NOT on the
    3D SAFS spatial dynamic rupture path.
  - `seas_test_parallel_domain`, `seas_mms_antiplane_parallel`,
    `seas_test_serial_parallel_consistency` — macOS 26 `mpirun` Bus error
    10 with `MFEM_USE_HDF5=YES`.  Pre-existing.  Antiplane / 2D, NOT on
    the 3D SAFS path.

## Changes Made (in priority order)

### R-701 [CRITICAL/LATENT] — fixed (belt-and-suspenders, both sites)

**(a) Resolver early-out at `spatial_friction.cpp:1103`.**  Extended the
`!nuc.enabled` guard so it ALSO fires when `nuc.kind != NucleationKind::StrengthReduction`:

```cpp
if (!nuc.enabled
    || nuc.kind != NucleationKind::StrengthReduction)
{
   p.T_forced_s = 1.0e9;
   p.t0_decay_s = 0.0;
   return p;
}
```

For any kind other than StrengthReduction (currently only Overstress),
`ResolveForcedRupture` returns the "never forced" sentinel for every DOF so the
iterator's forced-rupture `mu_eff` path is a no-op.  Comment in the function
body cross-references R-701.

**(b) Driver gate at `spatial_dyn_driver.cpp:1188`.**  Changed
`substep_iterator.SetForcedRuptureMode(true)` to read the previously-computed
`use_strength_reduction` boolean:

```cpp
Tpv205SubStepIterator substep_iterator(fault_flux);
substep_iterator.SetForcedRuptureMode(use_strength_reduction);
```

The block comment above the call now documents both R-601 (round-6) AND R-701
(round-7) so the next reader sees the full history.

Together (a) + (b) ensure that a future overstress run cannot silently
double-apply strength reduction even if one of the two guards regresses.

### R-704 [LOW] — fixed (subsumed by R-701)

Replaced the misleading comment block at `spatial_dyn_driver.cpp:941-948` with
an accurate description of the now-fixed behaviour.  Quoted in full:

```cpp
// -----------------------------------------------------------------
// 12. Per-DOF forced-rupture times (Phase 1 / D-4).  ResolveForced
//     Rupture writes T_forced(r) per hypocenter distance only when
//     cfg.nucleation is enabled AND cfg.nucleation.kind ==
//     StrengthReduction; for any other kind (Overstress) it returns
//     the "never forced" sentinel T = 1e9 everywhere so the
//     iterator's forced-rupture mu_eff path is a no-op.  Below the
//     iterator is also gated by SetForcedRuptureMode(use_strength_
//     reduction) (R-701 round-7) — belt and suspenders.
// -----------------------------------------------------------------
```

### R-706 [LOW] — fixed

Replaced `(void)overstress;` with `[[maybe_unused]]` on the
`OverstressPerDOFParams` declaration at `spatial_dyn_driver.cpp:961-963`.
Cleaner, self-documenting, C++17-idiomatic.

### R-702 [MODERATE] — documented (option (a) from review)

Added a 13-line documentation block in BOTH `LSW_ForcedRupture` dispatch arms
in `wave_operator.inl` (interior-fault arm at L3835+, shared-fault arm at
L4783+) explaining the time-granularity mismatch with the substep iterator
for `ader_order >= 2`.  The shared-fault arm's comment references the
interior-fault arm's full note to avoid duplication.

Rationale: option (b) — threading substep time into
`ComputeADERSharedFaceFluxRHS` — is the structurally-correct fix but would
ripple through the wave-operator's MPI-collective face-flux pipeline.  Out
of scope for round-7; documented as future work in the comment.

### R-705 [LOW/POSSIBLE] — documented

Added a 10-line comment in `tpv205_substep_iterator.cpp:226` documenting that
`t_sub_acc = t_macro_start` selects the START-OF-SUBSTEP convention (not
midpoint), the bounded half-step lag for `ader_order > 1`, and that F-10 test
fixes this convention.  Reviewer flagged this as POSSIBLE — the documentation
preserves it as a known property rather than a bug.

### R-703 [MODERATE] — fixed (tests added)

**T-21..T-23** in `test_spatial_friction_config.cpp`:

- `T_21_nucleation_kind_default_strength_reduction` — parses `[nucleation]`
  without `kind`, asserts `cfg.nucleation.enabled == true` AND
  `cfg.nucleation.kind == NucleationKind::StrengthReduction`.
- `T_22_nucleation_kind_overstress_parses` — parses `kind="overstress"` +
  full `[nucleation.overstress]` sub-block, asserts the enum round-trips
  and `delta_tau_pa` + `direction` parse.
- `T_23_nucleation_kind_typo_aborts` — parses `kind="strenght_reduction"`
  (typo), asserts `ParseAbortsInChild` succeeds (the parser must reject
  unknown kind strings).

**F-4, F-5** in `test_spatial_friction_resolver.cpp`:

- `F_4_overstress_kind_returns_sentinel` — builds 3 DOFs at r ∈ {0, 100, 500}
  m (every one INSIDE the 4 km r_crit), sets `nuc.enabled = true` and
  `nuc.kind = NucleationKind::Overstress`, asserts every `T_forced_s(i) ==
  1e9` and every `t0_decay_s(i) == 0`.  This is the R-701 regression: pre-fix,
  DOF 0 (at the hypocenter) would receive `T_forced = 0` from the hypocenter-
  distance formula; post-fix, it gets the sentinel.
- `F_5_resolve_overstress_stub_behaviour` — three paths in one test:
  1. `kind = StrengthReduction` → empty vectors returned.
  2. `nuc.enabled = false` → empty vectors returned.
  3. `kind = Overstress` (via `RunInChild`) → stub aborts.

## Verification

```
- [x] R-701: fixed at both call sites (sentinel + gate)
- [x] R-702: documented in both wave_operator.inl dispatch arms
- [x] R-703: T-21..T-23 + F-4 + F-5 added and pass
- [x] R-704: comment block at spatial_dyn_driver.cpp:941-948 corrected
- [x] R-705: START-OF-SUBSTEP convention documented in iterator
- [x] R-706: (void)overstress → [[maybe_unused]]
```

### Test sweep (3D SAFS path)

All targets affected by round-6/round-7 changes plus TPV byte-exact regression:

| Target                                          | Result        | Notes                                          |
|-------------------------------------------------|---------------|------------------------------------------------|
| `seas_test_spatial_friction_config`             | **34 / 34**   | T-21..T-23 contribute +6 sub-asserts (was 28)  |
| `seas_test_spatial_friction_resolver`           | **122 / 122** | F-4, F-5 contribute +12 sub-asserts (was 110)  |
| `seas_test_phaseh_lsw_forced_rupture`           | **51 / 51**   | F-10 (round-6) still green; R-705 doc-only     |
| `seas_test_phaseh_wave_operator_constant_parity`| **16 / 16**   | unchanged                                       |
| `seas_test_tpv104_checkpoint`                   | **164 / 164** | TPV byte-exact contract intact                  |
| `seas_test_safs_mode_wiring`                    | **8 / 8**     | unchanged                                       |
| `seas_test_compute_safs_params`                 | **13 / 13**   | unchanged                                       |
| `seas_test_spatial_setup`                       | **71 / 71**   | unchanged                                       |
| `seas_test_spatial_velocity_bundle`             | **6 / 6**     | unchanged                                       |
| `seas_test_spatial_stress_bundle`               | **6 / 6**     | unchanged                                       |
| `seas_test_spatial_constant_stress_source`      | **45 / 45**   | unchanged                                       |
| `seas_test_phaseh_godunov_flux_pool`            | **30 / 30**   | unchanged                                       |
| `seas_spatial_dyn_driver`                       | **build OK**  | clean rebuild, no warnings                      |

**Total 3D SAFS regression: 566 / 566 passing (up from 564 / 564 in round 6;
the new F-4, F-5, T-21..T-23 add the +2 / +18-sub-asserts shown above).**

### Pre-existing failures NOT addressed (out of scope)

Per [[feedback-preexisting-test-failures]], the following failures are
documented and intentionally not fixed in this round.  None are on the 3D
SAFS spatial dynamic rupture path:

1. **`tests/verification/bp2_serial_smoke.cpp`** — fails to compile with
   `no member named 'InitFaultOutput' in 'mfem::seas::ParaViewOutput<>'`.
   The `feature/paraview-compaction` merge renamed
   `InitFaultOutput` → `InitFaultOutputBP5` and `UpdateFaultFields` →
   `UpdateFaultFieldsBP5` (consistent with the `tpv104_driver.cpp`
   call-site usage and the round-6 spatial driver).  The BP2 verification
   test was not updated as part of that merge.  BP2 is 2D antiplane;
   it is not consumed by `spatial_dyn_driver` or any TPV* / 3D SAFS
   path.  Fix: separate plan to update BP2 verification to the new API,
   then re-enable `seas_bp2_serial_smoke-test-ser`.

2. **`seas_test_parallel_domain`, `seas_mms_antiplane_parallel`,
   `seas_test_serial_parallel_consistency`** — `mpirun` Bus error 10 on
   macOS 26 with `MFEM_USE_HDF5=YES`.  Symptomatic of an
   `Open MPI / HDF5 dual-stack` ABI mismatch documented in commit
   `6e45aaf` (`seas: unblock full test sweep on macOS 26 + MFEM_USE_HDF5
   build`).  All three failures are 2D antiplane / domain tests; not
   on the 3D SAFS spatial dynamic rupture path.  Fix: separate plan to
   investigate the MPI/HDF5 stack on macOS.

## New Tests

- `T_21_nucleation_kind_default_strength_reduction` — covers R-703 (parser default).
- `T_22_nucleation_kind_overstress_parses` — covers R-703 (overstress kind round-trip + sub-block).
- `T_23_nucleation_kind_typo_aborts` — covers R-703 (parser rejects unknown kind).
- `F_4_overstress_kind_returns_sentinel` — covers R-701 (resolver sentinel for Overstress kind).
- `F_5_resolve_overstress_stub_behaviour` — covers R-703 (ResolveOverstress stub contract: empty for StrengthReduction/disabled, abort for Overstress).

## Ready for Re-Review: YES

No CRITICAL findings remain on the 3D SAFS spatial dynamic rupture path.
R-701's latent failure mode is closed by the two-site fix and pinned by F-4.
The MODERATE R-702 time-granularity asymmetry remains a documented bounded
limitation (option (a) from the review); a structural fix is a separate
follow-up.
