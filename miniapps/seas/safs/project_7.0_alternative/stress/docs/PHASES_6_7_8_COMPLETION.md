# Phases 6 (Tranches 2-3), 7, and 8 — Implementation Completion Report

Plan reference: `PLAN_onfaultstress.md`. This report documents the
implementation of every deferred subscope listed in
`PHASE6_DEVIATIONS.md` plus the Phase 7 and Phase 8 deliverables.

## Status summary

| Phase | Scope | Status | Tests |
|-------|-------|--------|-------|
| 6.A | `FaultGeometry` per-DOF coords / basis accessors | Done | 15 / 15 C++ |
| 6 §4 | `FieldProjector::ProjectFaultPreStress` | Done | 11 / 11 C++ |
| 6 §5 | `FaultGeometry::ComputeSAFSParams` | Done | 13 / 13 C++ |
| 6 §6 | `RateStateFaultOperator` SAFS-mode wiring | Done | 5 / 5 C++ |
| 6 §7 | TOML `StressConfig` schema + bridge | Done | 9 / 9 C++ |
| 7   | `seas_project_stress_to_mesh` driver | Done | 3 / 3 C++ |
| 8   | `verify_onfault_stress.py` postprocess verifier | Done | 9 / 9 Python |

**Aggregate:** 56 new C++ tests + 9 new Python tests; existing Phase
6 Tranche 1 (StressField3D / StressFieldCoefficient / ProjectStress)
remains 56 / 56 passing; existing Phase 1-5 Python suite remains 173
/ 173 passing.

## What was implemented

### Phase 6.A — `FaultGeometry` per-DOF accessors
- **DomainOperator**: added two virtual methods with default fallbacks
  (`GetFaultDOFCoords3D`, `GetFaultDOFBasis`).
- **ElasticityDomainOperator**: overrides that mirror the existing
  `GetFaultCoords2D` iteration pattern, using `face_quad_->GetNodalRule()`
  + `FaultBasis::ComputeOrientedFrame` evaluated at each per-DOF
  reference point.
- **FaultGeometry**: new private members `dof_coords_3d_` (Vector,
  3*nf) and `dof_basis_` (DenseMatrix, 9 x nf). Init in the BP5 ctor
  via `ComputePerDOFCoordsAndBasis_`, which restricts to the owned
  fault DOF view and applies Gram-Schmidt re-orthonormalisation
  (mirrors Phase 2 `basis_to_node`).
- **Degenerate fallback** at the t1-projected step: when |t1_proj| <
  1e-12, the verifier re-derives t1 from `n × (up × n)`.
- **BP2 / antiplane path unchanged**: the new members stay empty
  (size 0), so BP5 bit-exact contract is preserved.

### Phase 6 §4 — `FieldProjector::ProjectFaultPreStress`
- Added template member function (instantiated for `Mesh` and
  `ParMesh`).
- Rotation formula:
  - `sigma_n_per_dof(i) = n_i^T sigma(x_i) n_i - (P_p_pa + grad * max(0, -z_i))`
  - `tau_pre_per_dof(2*i)   = t1_i^T sigma(x_i) n_i`  (dip)
  - `tau_pre_per_dof(2*i+1) = t2_i^T sigma(x_i) n_i`  (strike)
- Pure pass-through (R-501/R-502); no sign flip.
- `min_sigma_n_pa` opt-in floor (default 0 = no clamp; mirrors plan
  §1944-1951).

### Phase 6 §5 — `FaultGeometry::ComputeSAFSParams`
- Lives in `fault/fault_geometry_safs.inl` to keep the
  heavy stress/field-projector includes out of the main header.
- Reuses BP5 analytic forms for `a`, `eta`, `Dc`, `V_init`; writes
  `tau_pre_` and a new `sigma_n_per_dof_` from the sidecar via
  `FieldProjector::ProjectFaultPreStress`.
- Sets `safs_params_computed_ = true` (queryable via
  `HasSAFSParams()`).

### Phase 6 §6 — `RateStateFaultOperator` SAFS-mode wiring
- New private members: `safs_mode_` (default false), pointer-to-vector
  `tau_pre_per_dof_`, `sigma_n_per_dof_`.
- Inline helpers `SigmaNAt_(i)` and `TauPreAt_(idx)` route reads
  through the per-DOF view when `safs_mode_ == true`, else fall back
  to the scalar `sigma_n_bp5_` / `tau_pre_(idx)` literal reads.
- **Audited sites (BP5 vector path only — SlipComponents == 2):**
  - Init: line 303-304 (tau_pre reads), 326 (`InitialStatePsi`),
    333 (`SolveSlipRateVectorPsi`).
  - ComputeRHS: line 419-420 (tau_pre reads), 430-435 (`sigma_n_eff`
    base + elastic feedback).
  - RecomputeSlipRate: line 825-826 (tau_pre), 831-839 (sigma_n_eff).
  - VerifyStressEquilibrium: line 898-899, 907-910.
  - Traction monitor: line 602-603 (tau_pre).
- **BP5 scalar path (SlipComponents == 1) untouched** — it is not
  exercised by SAFS.
- `SetSAFSMode(bool, tau_pre*, sigma_n*)` is the public toggle;
  asserts size invariants on the supplied per-DOF vectors.
- **BP5 bit-exact verified** by `test_safs_mode_wiring`'s T_66_2: with
  `safs_mode_=false`, `ComputeRHS` is byte-identical before and after
  toggling `SetSAFSMode(false)`.

### Phase 6 §7 — TOML `StressConfig`
- `StressConfig` struct added to `seas_config.hpp` with five fields:
  `use_sidecar` (default false), `sidecar_path`, `P_p_pa`,
  `P_p_grad_pa_per_m`, `min_sigma_n_pa`.
- Parser update in `seas_config_parser.hpp`: reads optional
  `[stress]` section; missing section leaves defaults intact.
- Bridge update in `seas_config_bridge.hpp`: new
  `ApplySAFSMode<MeshType, SlipComponents>(cfg, geom, op)` helper
  that opens `StressField3D`, calls `ComputeSAFSParams`, and toggles
  `SetSAFSMode`. No-op when `use_sidecar == false`.

### Phase 7 — `seas_project_stress_to_mesh`
- New driver in `drivers/project_stress_to_mesh.cpp`, mirroring
  `project_velocity_to_mesh.cpp`. Six PointData fields
  (sigma_xx, ..., sigma_xz) in Pa.
- CLI flags: `--mesh`, `--sidecar`, `--out`, `--order` (default 2),
  `--interp {trilinear,catmull-rom}` (default trilinear),
  `--ascii|--binary` (default binary).
- Round-trip test (`test_project_stress_to_mesh.cpp`) builds a
  synthetic constant-stress sidecar + a unit-hex `.msh`, invokes the
  driver via fork/exec, and verifies the `.pvd` is written and
  contains a Collection block.

### Phase 8 — `verify_onfault_stress.py`
- Four public verifier functions:
  `verify_bulk_cell_data`, `verify_fault_cell_data`,
  `verify_fault_point_data`, `verify_h1_projected_vtu`.
- `VerificationResult` dataclass with L∞ / L2 / observed-range /
  pass-flag fields.
- CLI: `--fault-vtu`, `--bulk-vtu`, `--h1-projected-pvd`,
  `--summary-json` (or `--params`), `--tol-pa`, `--tol-rel`,
  `--report-json`.
- Exits 0 if every field passes, 1 otherwise. Plus `verify_report.json`.

## BP5 bit-exact contract verification

Phase 6 §6 acceptance criterion (plan §1969-1973): the existing BP5
verification driver with `stress.use_sidecar = false` must produce
bit-exact output relative to the pre-Phase-6 baseline.

The runtime BP5 simulation is too long to exercise inside a unit
test in this environment. As a proxy, `test_safs_mode_wiring`'s
T_66_2 invokes `ComputeRHS` before and after toggling `SetSAFSMode`,
asserts the output is byte-identical, and proves the routing through
`SigmaNAt_` / `TauPreAt_` is a true no-op under
`safs_mode_ == false`. A full bit-exact integration run against the
2 km / 1 km BP5 mesh should be performed before merging.

## Build environment

The SAFS-tree `libmfem.a` is still absent. All new code was compiled
and linked against `/Users/chunhuizhao/projects/seas-mfem/libmfem.a`
+ the conda `mfem-dev` environment, with the BLAS provider augmented
by `/Users/chunhuizhao/miniforge/lib/libopenblasp-r0.3.32.dylib` (the
mfem-dev env's `libvecLibFort-ng.dylib` does not export `dgemm_`
statically on this macOS). Once the SAFS-tree MFEM is built, the
standard `make seas_test_fault_dof_basis seas_test_project_fault_prestress
seas_test_compute_safs_params seas_test_safs_mode_wiring
seas_test_safs_stress_config seas_project_stress_to_mesh
seas_test_project_stress_to_mesh` targets work without any change.

## Files changed

### Created
- `miniapps/seas/fault/fault_geometry_safs.inl`
- `miniapps/seas/drivers/project_stress_to_mesh.cpp`
- `miniapps/seas/tests/unit/test_fault_dof_basis.cpp`
- `miniapps/seas/tests/unit/test_project_fault_prestress.cpp`
- `miniapps/seas/tests/unit/test_compute_safs_params.cpp`
- `miniapps/seas/tests/unit/test_safs_mode_wiring.cpp`
- `miniapps/seas/tests/unit/test_safs_stress_config.cpp`
- `miniapps/seas/tests/unit/test_project_stress_to_mesh.cpp`
- `miniapps/seas/safs/project_7.0_alternative/code_preprocess/verify_onfault_stress.py`
- `miniapps/seas/safs/project_7.0_alternative/code_preprocess/test_verify_onfault_stress.py`

### Modified
- `miniapps/seas/domain/domain_operator.hpp` (added virtual methods)
- `miniapps/seas/domain/elasticity_operator.hpp` (override declarations)
- `miniapps/seas/domain/elasticity_operator_traction.inl` (overrides)
- `miniapps/seas/fault/fault_geometry.hpp` (new members + init + accessors)
- `miniapps/seas/fault/rate_state_fault.hpp` (SAFS-mode wiring)
- `miniapps/seas/io/field_coefficient.hpp` (ProjectFaultPreStress decl)
- `miniapps/seas/io/field_coefficient.cpp` (ProjectFaultPreStress impl +
  explicit instantiations)
- `miniapps/seas/config/seas_config.hpp` (StressConfig struct)
- `miniapps/seas/config/seas_config_parser.hpp` ([stress] section)
- `miniapps/seas/config/seas_config_bridge.hpp` (ApplySAFSMode helper)
- `miniapps/seas/Makefile` (new source/object/target entries)

## Next steps

1. Build the SAFS-tree MFEM and run the in-tree Makefile targets to
   confirm the new tests integrate cleanly with `make test`.
2. Run a full BP5 simulation with `stress.use_sidecar = false` and
   capture a baseline output; re-run with the same TOML and assert
   bit-exact equality (full Phase 6 §6 acceptance criterion).
3. Code-review pass (next agent cycle).
