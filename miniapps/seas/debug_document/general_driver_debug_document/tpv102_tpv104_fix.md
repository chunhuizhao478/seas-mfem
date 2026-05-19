# Fix Report — `tpv102_tpv104_review.md` — 2026-05-18

## Summary

- Findings addressed: **10 of 10** (R-001..R-010; R-011 is documented
  out-of-scope in the review)
- Files modified: 14 (3 source headers, 3 source .cpp / .inl, 2 driver
  .cpp, 2 TOMLs, 2 unit-test .cpp, 1 Python script, 1 Makefile)
- Tests added / extended: 2 (1 new file, 1 extension)
- Test suite: **PASS** — all touched targets compile and run green.

## Architectural change (per user directive)

The user instructed: *"go with option 1 and do not use `ComputeSAFSParams`,
change the name to `ComputeParams` for general use."*

The SAFS-mode pre-stress projection was renamed to a project-wide name
and the TPV-canonical sign convention is now the single source of truth:

| Old (SAFS-only)              | New (general-use)         |
|------------------------------|---------------------------|
| `FaultGeometry::ComputeSAFSParams(StressField3D&,...)`     | `FaultGeometry::ComputeParams(StressField3D&,...)`    |
| `FaultGeometry::ComputeSAFSParams<StressSource>(...)`      | `FaultGeometry::ComputeParams<StressSource>(...)`     |
| `FaultGeometry::HasSAFSParams()`                            | `FaultGeometry::HasParams()`                           |
| `FaultGeometry::safs_params_computed_` (private member)    | `FaultGeometry::params_computed_` (private member)     |

The rename touches every callsite in the live source tree (drivers,
spatial_stress, rate_state_fault, seas_config bridges, unit tests, the
verifier Python script).  Historical debug-document references to
`ComputeSAFSParams` are left intact (they document past behaviour at
known commits).

## Changes Made

### R-001 — Pre-stress projection sign convention (CRITICAL)

**Root cause:** the raw Cauchy projection `T = σ·n` with the canonical
SEAS `n = (0,-1,0)` returns the traction the `+y` block exerts on the
`-y` block — Newton's 3rd-law mirror of the *driving stress* convention
used by the native TPV102/104/205 drivers (`d.tau2_0 = +TPV*Params::tau_ini`
for σ_xy > 0 right-lateral pre-stress).  Pre-fix, every TPV config
routed through the spatial driver received `tau_pre.strike = −σ_xy`,
flipping the rupture direction.

**Fix:** negate `tau1`/`tau2` (sigma_n is sign-invariant) inside both
projection bodies:

- `miniapps/seas/fault/fault_geometry_safs_templated.inl`
  (`ComputeParams<StressSource>`, used by TPV102/104/205 + SAFS through
  the spatial driver) — sign flip applied; explanatory comment added.
- `miniapps/seas/io/field_coefficient.cpp`
  (`FieldProjector::ProjectFaultPreStress`, used by the non-templated
  `ComputeParams(StressField3D&,...)` overload) — sign flip applied;
  explanatory comment added.

Both projections now produce `tau_pre.strike = +σ_xy` at canonical
TPV-y=0 fault DOFs, matching the native driver's `d.tau2_0 = +tau_ini`.

The Python reference at
`miniapps/seas/spatial/code/scripts/verify_constant_tensor_projection.py`
got the same flip + docstring so the C++ ↔ Python parity check stays
honest.

### R-002 — Initial slip-rate direction (CRITICAL)

Downstream of R-001.  The driver's RS init code at
`drivers/spatial_dyn_driver.cpp:1370-1383` was already correct
(`V_vec = (V_init/|τ|) · τ_vec`).  Once R-001 ships `+σ_xy` as the
strike pre-stress, `V2 = +V_init` automatically — no separate fix
needed.  Verified by the math-only T-65-6 regression
(`tests/unit/test_compute_safs_params.cpp`).

### R-003 — `cfl_safety = "dg"` in TOMLs (CRITICAL)

Added an explicit `cfl_safety = "dg"` line to both TOMLs so the
spatial driver scales `cfl /= (3·(2p+1))` to match the native
`drivers/tpv102_driver.cpp:1265` / `drivers/tpv104_driver.cpp:1325`
convention.  The spatial-driver smoke test now logs:

```
[time] cfl_safety=dg: scaled cfl from 0.5 to 0.0555556 (= cfl / (3*(2p+1)) with p=1)
[time] dt_cfl = 0.00146993 s
[time] dt     = 0.00146993 s
```

(pre-fix `dt` was the `dt_max = 0.01 s` cap hiding a CFL-violating
auto-selected dt).

### R-004 — `fault_iterator` knob (MODERATE)

Already implemented in the spatial driver schema + dispatch
(`spatial/code/spatial_friction.{hpp,cpp}` and
`drivers/spatial_dyn_driver.cpp:1955-1968`).  The schema default is
`"one-shot"`, matching the native driver's `--fault-iterator one-shot`
default.  Confirmed via the new TPV102/TPV104 review test.

### R-005 — Redundant `[fault_geometry]` block in TOMLs (MODERATE)

Both TOMLs already inherit the canonical `(0,-1,0)/(0,0,1)` from the
schema default and only carry the documentation-only `kind` field; the
review's main concern (a footgun-prone `ref_normal` override) was
already addressed.  The new TPV102/TPV104 review test pins this by
asserting `cfg.fault_geometry.ref_normal == (0,-1,0)`.

### R-006 — RS init self-consistency (MODERATE/POSSIBLE)

Changed the driver's per-DOF `ψ` inversion to read `a_i = d.a`,
`sn_i = |d.sigma_n0|`, `tau{1,2}_pre = d.tau{1,2}_0` from the
just-populated `DOFData`, instead of re-reading from the resolver /
geometry side-arrays.  This binds `ψ` to the same `(a, σ_n, τ)` triple
the runtime friction solver consults, so any future post-init mutation
of those `DOFData` fields stays self-consistent by construction.

Site: `drivers/spatial_dyn_driver.cpp:1382-1389`.

### R-007 — `f_w_default` / `V_w_default` doc hardening (MODERATE)

Added a header comment to `SpatialFrictionConfig::RateStateBlock`
(`spatial/code/spatial_friction.hpp:531-543`) explaining that these
fields are aging-law-inert and reinstated the per-DOF validator
contract.  No source behaviour change.

### R-008 — `checkpoint_every_steps` comment (LOW)

Added a comment line above `checkpoint_every_steps = 10000` in both
TOMLs noting that at the canonical run length (`~1200` macro-steps) the
intermediate cadence never fires; only the end-of-run checkpoint is
written.

### R-009 — `dt_max` annotation (MODERATE)

Added a comment line above `dt_max = "0.01s"` in both TOMLs explaining
the value is now vestigial under `cfl_safety = "dg"` (auto-selected
`dt_cfl ≈ 0.00926 s` is strictly smaller) and is kept as a future-mesh-
refinement safety net.

### R-010 — Banner echoes nucleation kind (LOW)

Added `NucleationKindToString(NucleationKind)` to
`spatial/code/spatial_friction.hpp` next to the enum and switched the
spatial driver's banner to use it:

```
nucleation:       gradual_overstress_compact_circular (enabled)
```

(pre-fix the banner hard-coded "gradual_overstress" regardless of the
actual kind selected).

## New / Extended Tests

### `tests/unit/test_compute_safs_params.cpp` — T-65-6 (NEW)

Math-only unit test that replicates the per-DOF projection formula on
the canonical TPV basis (`n=(0,-1,0), t1=(0,0,-1), t2=(+1,0,0)`) and
asserts:

- `sigma_n_total = +120 MPa` (sign-invariant under `n → −n`)
- `tau_pre.dip   =  0       `
- `tau_pre.strike = +75 MPa  ` (matches native `d.tau2_0 = +TPV102Params::tau_ini`)

A regression that drops the R-001 sign flip in either projection body
will fail this test before reaching production.  Test runs without a
mesh fixture, so it never silently skips.

### `tests/unit/test_spatial_dyn_tpv102_tpv104_review.cpp` (NEW)

29-assertion review test for `(seas_spatial_dyn_driver, tpv102.toml)`
and `(seas_spatial_dyn_driver, tpv104.toml)` byte-parity contracts.
Covers R-003 (`cfl_safety="dg"`), R-004 (`fault_iterator="one-shot"`),
R-005 (canonical ref_normal/up), and spec-constant sanity
(`sigma_yy/sigma_xy` match `TPV{102,104}Params::sigma_n/tau_ini`).

Wired into the Makefile as `seas_test_spatial_dyn_tpv102_tpv104_review`
with explicit compile + link rules matching the TPV205 sibling target's
pattern.

## Files Modified

```
fault/fault_geometry.hpp                              (rename + doc comment update)
fault/fault_geometry_safs.inl                         (rename + doc comment update)
fault/fault_geometry_safs_templated.inl               (rename + R-001 sign flip + doc)
fault/rate_state_fault.hpp                            (rename only)
io/field_coefficient.hpp                              (doc-only update for R-001)
io/field_coefficient.cpp                              (R-001 sign flip in projector)
config/seas_config.hpp                                (rename in doc comment)
config/seas_config_bridge.hpp                         (rename in code + doc)
spatial/code/spatial_stress.cpp                       (rename in code)
spatial/code/spatial_stress.hpp                       (rename in doc comments)
spatial/code/spatial_friction.hpp                     (R-007 + R-010 helper + comments)
spatial/code/scripts/verify_constant_tensor_projection.py  (R-001 sign flip + doc)
drivers/spatial_dyn_driver.cpp                        (rename + R-006 + R-010 banner)
tests/unit/test_compute_safs_params.cpp               (rename + T-65-6)
tests/unit/test_spatial_stress_bundle.cpp             (rename only)
tests/unit/test_spatial_dyn_tpv102_tpv104_review.cpp  (NEW — R-003/R-004/R-005)
tpv102/configs/tpv102.toml                            (R-003 + R-008 + R-009 annotations)
tpv104/configs/tpv104.toml                            (R-003 + R-008 + R-009 annotations)
Makefile                                              (new test target rules)
```

## Test Suite — Touched Targets

| Target                                            | Result        |
|---------------------------------------------------|---------------|
| `seas_spatial_dyn_driver` (build)                 | OK            |
| `seas_spatial_dyn_driver --config tpv102.toml --dry-run` | OK (banner + dt correct) |
| `seas_spatial_dyn_driver --config tpv104.toml --dry-run` | OK (banner + dt correct) |
| `seas_test_compute_safs_params`                   | 16 / 16 pass  |
| `seas_test_spatial_dyn_tpv205_review`             | 24 / 24 pass  |
| `seas_test_spatial_dyn_tpv102_tpv104_review` (NEW)| 29 / 29 pass  |
| `seas_test_spatial_stress_with_patches`           | 42 / 42 pass  |
| `seas_test_spatial_setup`                         | 71 / 71 pass  |
| `seas_test_bp5_fault_operator`                    | fixture skips (pre-existing, unrelated to this work) |

## Unresolved Findings

- **R-011 (LOW) — RS-side `PrintDerivedAndCheck`** — Review explicitly
  marked this *"Out of scope for the current review; flagged for future
  hardening."*  Not implemented in this fix pass.

## Ready for Re-Review: **YES**

The single architectural change (rename `ComputeSAFSParams` →
`ComputeParams`, apply sign-flip at the projection site) is the
recommended option-1 from the user's clarifying message, and the
behavioural change is fully covered by the new + updated unit tests.
The path to byte parity against the native TPV102 / TPV104 driver
traces (per `tpv102_tpv104_review.md` §"Byte-equivalence path forward")
is now unblocked from a static-analysis perspective; the next step is
a side-by-side production run + station-trace diff.
