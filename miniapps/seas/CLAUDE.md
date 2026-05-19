# SEAS Miniapp - Claude Code Instructions

MFEM-based SIPG/BR2 Discontinuous Galerkin code for SCEC SEAS benchmark problems (BP1, BP2, BP5). Couples a DG elasticity domain solver with regularized Dieterich-Ruina rate-and-state friction on an embedded fault interface via a quasi-dynamic approximation. Time integration uses adaptive Dormand-Prince RK45.

## Before ANY Code Change

1. Read `ARCHITECTURE.md` for class hierarchy, data flow, and execution model
2. Read debug documents in `debug_document/bp5_debug_document/` for history of what has been tried
3. Run `make test` (unit tests) to confirm baseline passes
4. For BP5 changes: check Tandem reference at `/Users/chunhuizhao/projects/tandem/`

## Reference Materials

| Resource | Location |
|----------|----------|
| Tandem source code | `/Users/chunhuizhao/projects/tandem/` |
| Porting reports | `/Users/chunhuizhao/projects/mfem/miniapps/seas/document/` |
| Debug history | `debug_document/bp5_debug_document/bp5_debug_v*.md` |
| SCEC BP5 spec | https://strike.scec.org/cvws/seas/download/SEAS_BP5_QD.pdf |
| Architecture | `ARCHITECTURE.md` (this directory) |

## Canonical Coordinate System (project-wide)

**This is the single source of truth for fault geometry in mfem-seas.**
Every benchmark in this project (BP1, BP2, BP5, TPV102, TPV104, TPV205,
SAFS) is set up in the same Cartesian frame:

| Axis | Meaning | Sign convention |
|------|---------|-----------------|
| **x** | along-strike (east) | + = east |
| **y** | fault-normal | + = away from the canonical "+" side |
| **z** | vertical | **z = 0 at free surface; z < 0 below** (depth is `\|z\|`) |

The fault is the **y = 0** plane (vertical strike-slip).  The
`FaultBasis` (Tandem convention) is seeded with:

```
ref_normal = (0, -1, 0)
up         = (0,  0, 1)
```

This is **hard-coded** in `dynamic/wave_operator.inl:349-350`.  The
`SpatialFrictionConfig` schema default in
`spatial/code/spatial_friction.hpp::FaultGeometrySpec` matches.  TOML
configs (TPV102/104/205, BP5/SAFS) should NOT re-state these values;
they inherit them from the schema default.  The spatial driver
hard-fails at startup if a TOML overrides them — see REVIEW R-003 in
`debug_document/general_driver_debug_document/REVIEW.md`.

Derived fault-local frame (from `FaultBasis::ComputeOrientedFrame`,
`strike = up × n_ref`, `dip = strike × n_ref`):

```
normal   n  = (0, -1, 0)
strike   t2 = (+1, 0, 0)
dip      t1 = (0,  0, -1)   (down into earth)
```

Consequently in `DOFData`:
- `V1`, `slip1`, `tau1_0`, `tau1_corr` are the **dip** components.
- `V2`, `slip2`, `tau2_0`, `tau2_corr` are the **strike** components.

TOML coordinates also use this frame.  E.g. for TPV205's `[hypocenter]`
at spec depth 7.5 km below the surface:

```
[hypocenter]
x = 0.0
y = 0.0          # on the fault
z = -7500.0      # MESH z = negative spec depth
```

Stress patches, friction depth bounds, nucleation centres, etc. all
use the same `(x, y, z)` mesh-coord triple.  When the SCEC spec lists a
quantity at "depth 7.5 km", encode it as `center_z_m = -7500.0` in the
TOML.

### Outliers

- **TPV31** (`tpv31/configs/tpv31.toml`): the spec describes a fault on
  the **z = 0** plane with **y** as the depth axis (orthogonal
  geometry).  TPV31 currently overrides ref_normal/up to its own
  frame, but this trips the spatial driver's R-003 guard.  Rotating
  TPV31's mesh + TOML into the canonical y = 0 / z = depth frame is a
  follow-up task tracked separately.  Until then, TPV31 cannot be run
  through `seas_spatial_dyn_driver`.

## Critical Numerical Details

These lessons were learned through extensive debugging (v1-v62). Violating any of them causes simulation failure.

### Sign Conventions

- **Slip rate direction**: PARALLEL to traction tau (not antiparallel). `V_vec = (V_abs / tau_abs) * tau_vec`. Antiparallel sign creates positive feedback -> unbounded growth (debug v8).
- **Pre-stress direction**: PARALLEL to initial velocity. `tau_pre = tau0_scalar * Vi / |Vi|` (debug v8).
- **Dip direction**: (0, 0, +1) = downward into earth. Was wrongly (0, 0, -1) = upward (debug v10).
- **DG face sign**: `sign = (nor(0) > 0) ? -1.0 : 1.0` applied consistently to slip embedding and traction extraction.
- **Normal stress**: sigma_n > 0 = compression (geology convention).
- **Depth coordinate**: Z=0 at surface, Z<0 is depth. `Wf` is positive, fault extends from z=0 to z=-Wf.
- **Fault-local tangent frame (BP5 / TPV102, project-wide)**: Uses the
  `FaultBasis` (Tandem) convention — `tangent1 = dip, tangent2 = strike`.
  For TPV102's vertical y=0 fault with `ref_normal=(0,-1,0)` and `up=(0,0,1)`:
  `can_t1 = (0, 0, -1)` (down-dip), `can_t2 = (+1, 0, 0)` (along strike).
  Therefore in `DOFData`: `V1/slip1/tau1_0/tau1_corr` are the **dip** components
  and `V2/slip2/tau2_0/tau2_corr` are the **strike** components.  TPV102 is
  pure strike-slip, so `tau2_0 = tau_ini`, `V2 = V_ini`, and `tau1_0 = 0`,
  `V1 = 0` (debug v7.0.0 R-801).  Pre-R-801 the interior-fault branch used
  `GodunovFlux::BuildFrame` (t1=strike) while the shared-fault branch used
  BP5 (t1=dip); the two have been unified on the BP5 convention.

### Friction Solver

- **Must use Brent's method**, not Newton. With large psi/a, Newton fails because F(V_lo=1e-30) > 0 while true solution is V << 1e-30 (debug v1).
- **Degenerate bracket guard**: When psi << 0, check `if (Fb >= 0.0) return V_hi = tau/eta` (frictionless limit). Allows RK45 to reject gracefully (debug v7).
- **No artificial tau floor**: Use `if (tau_abs <= 0.0) return 0`, not `if (tau_abs < 1e-30)`. The latter creates discontinuity in ODE RHS (debug v13 H7).
- **Vector slip (BP5)**: Solve for |V| via Brent in log10(V) space, then decompose direction separately.

### DG Method

- **IP penalty must include elasticity tensor coupling** for 3D. Scalar-only penalty causes blowup during nucleation (debug v12-v13 H3).
- **BR2 face_int uses unnormalized normal** from `CalcOrtho`, while `TestNormal` uses unit normal. These are different by |J_F| (debug v11-v12).
- **Shared faces**: Must use `GetSharedFaceTransformations()`. Only Elem1 is local; Elem2 via `pfes->GetFaceNbrFE()`. Slip assembly: only add Elem1 RHS contributions.
- **Stiffness matrix**: Assembled once and reused across all time steps and RK stages.

### Boundary Conditions (BP5)

- **All non-free-surface boundaries must get Dirichlet loading**: `u_X = sgn(Y) * Vp * t / 2`. With only bottom-boundary Dirichlet, nucleation zone starves (6.3 vs 50 kPa/yr loading rate), causing 3.7x longer recurrence (debug v9, v14-v15).
- **BCMode::FarField** is correct: Dirichlet on attrs 1-4 (vertical faces), Natural on 5-6 (top/bottom).
- Tandem mesh tags: Physical Surface 1 = Natural, 3 = Fault, 5 = Dirichlet.

### Time Stepping

- **Initial dt**: `0.01 * L_nuc / V_nuc` (not 0.5*). Too large -> RK45 stage amplification during nucleation (debug v7).
- **Max dt**: `0.1 * seconds_per_year` for SCEC output compliance (debug v13 H4).
- **RK45 tolerances**: atol=1e-7, rtol=1e-50 (pure absolute, matches PETSc TS default).
- **Parallel error reduction**: `MPI_Allreduce(MPI_MAX)` prevents deadlock from rank divergence.

### Initialization

- **4-phase**: (1) PreInit: slip=0, psi=steady-state. (2) Solve domain -> traction. (3) Init: compute psi from equilibrium. (4) Verify equilibrium error < 1e-6.
- **Per-DOF Dc in BP5**: `PsiToTheta(psi)` must use per-DOF Dc (0.13m in nucleation, 0.14m elsewhere), not global constant (debug v13 H1).

### Tandem vs SCEC Spec

- Tandem uses V_nuc=0.01, no delta_tau. SCEC spec: V_nuc=0.03, delta_tau=eta*V_nuc.
- Made configurable: `--V-nuc` and `--delta-tau-factor` CLI flags.
- **Follow SCEC spec unless explicitly matching Tandem** for cross-verification.

## Files Requiring Extreme Care

Any change to these files requires running full verification tests:

| File | Risk | Why |
|------|------|-----|
| `friction/dieterich_ruina.hpp` | Friction solver, Brent method | Wrong bracket/tolerance -> silent wrong V |
| `domain/elasticity_operator.hpp` | DG assembly, traction, BCs | Changes affect every simulation output |
| `fault/fault_basis.hpp` | Coordinate transforms | Sign error -> positive feedback -> blowup |
| `fault/rate_state_fault.hpp` | State evolution, slip assembly | Wrong interleaving -> corrupted state |
| `solver/seas_operator.hpp` | Coupling logic, init sequence | Wrong phase order -> bad equilibrium |
| `solver/time_stepper.hpp` | RK45, error control | Parallel error reduction critical for MPI |
| `config/bp5_params.hpp` | Spatial parameter functions | Wrong a(z) or Dc -> wrong friction regime |

## Known limitation — Gmsh `.msh` format

`mfem::Mesh::ReadGmshMesh` is a **Gmsh v2.2-only ASCII parser** (the
function accepts `version >= 2.2` but unconditionally runs the v2.2
body; there is no v4 branch).  All `.msh` artifacts consumed by this
miniapp's drivers and tests MUST be emitted in Gmsh v2.2 format
(`gmsh ... -format msh22 -o file.msh`).  A v4.x file aborts at
`mesh/mesh_readers.cpp:1628` with the misleading message
`Gmsh file : vertices indices are not unique` (the indices are fine;
the parser is misaligned on the v4 `$Nodes` header).

The SAFS mesher (`safs/project_7.0_alternative/meshing/code/run_nwcut_meshing.py`)
emits v2.2 directly.  When writing a new mesh producer or
documenting a `gmsh` command for a contributor, default to `msh22`.
Full investigation: `safs/project_7.0_alternative/meshing/docs/DEBUG_msh4_mfem_incompat.md`.

## Building and Testing

### Environment

```bash
conda activate mfem-dev       # For building (mpicxx, MPI, MUMPS, etc.)
conda activate pythonenv       # For Gmsh mesh generation
```

### Build

```bash
cd miniapps/seas
make all                       # Build everything (tests + drivers)
make seas_pseas                # Build parallel BP2 driver only
make seas_bp5_full             # Build BP5 full verification driver
```

### Testing

```bash
# Unit tests (fast, ~2 min total)
make test

# Individual test targets
make test-friction             # Friction law validation
make test-elasticity-operator  # DG elasticity assembly
make test-bp5-integration      # BP5 integration test
make test-bp5-smoke            # BP5 parallel smoke test (4 ranks)

# Parallel tests
mpirun -np 4 seas_test_parallel_elasticity
mpirun -np 8 seas_test_bp5_parallel_smoke

# Full verification (long-running)
mpirun -np 8 seas_bp5_full --mesh bp5/mesh/bp5_1000m.msh --tfinal 56844000000
```

### What Constitutes a Regression

1. V_max monotonically increasing (should peak then decrease in events)
2. Traction > 1 GPa (physical limit ~25 MPa for sigma_n = 25 MPa)
3. dt going to zero or NaN
4. zeroIn bracket failure (F(a) and F(b) same sign)
5. Recurrence >> 300 years or << 100 years (expected ~240 years)
6. Event slip << 1m or >> 10m
7. Significant dip slip in BP5 (should be ~0 for pure strike-slip)

## Do NOT

- Change sign conventions without understanding the full chain (slip -> traction -> friction -> state)
- Hardcode numerical constants; derive from parameters/mesh for generalizability
- Use Newton-Raphson for the friction solver (Brent is required)
- Apply scalar-only DG penalty for 3D elasticity (must include elasticity tensor)
- Mix structural refactoring with numerical changes in one commit
- Revert a previous fix without citing the debug document and getting explicit approval
- Assume Tandem is always correct; follow SCEC spec when they differ

## ParaView output mode

(Per plan §Phase 2b.  See `document/io_dev/PLAN_paraview_compaction_2026-04-28.md`.)

| Build flag                | Default fault back end | Notes |
|---------------------------|-----------------------|-------|
| `MFEM_USE_HDF5=YES`       | `FaultOutputMode::Hdf5` (single `.vtkhdf` per run) | production target; ParaView 5.11+ |
| `MFEM_USE_HDF5=NO`        | `FaultOutputMode::Vtu` (one binary VTU per cycle, raw appended UInt64) | non-HDF5 fallback |

Driver CLI (tpv102 / tpv104 / tpv205):

- `--paraview-fault-hdf5` — explicit opt-in (rejected at parse time on non-HDF5 builds).
- `--paraview-fault-vtu` — force the binary VTU back end (overrides the HDF5 default).
- `--paraview-fault-legacy-ascii` — revert to the pre-Phase-1 per-rank ASCII writer.  Debugging only; implies `--paraview-fault-vtu` and is incompatible with `--paraview-fault-hdf5`.

`seas_driver.cpp` does not currently construct a `ParaViewOutput`; the plan-mandated BP5 default of `--no-volume-pv` is therefore moot for that driver.  When ParaView output is added to `seas_driver.cpp` in a future change, the default volume save MUST be off unless `--volume-pv-dt X` is set explicitly.

### Phase 4 deferred deviation (R-313)

Plan §Phase 4 step 1 calls for `seas_driver.cpp` to default to `--no-volume-pv` for BP5 production runs.  This is **deferred** because `seas_driver.cpp` currently uses `ProbeOutput` and `ParallelBP5BenchmarkOutput` for its output — there is no `seas::ParaViewOutput` instance in that driver to apply the gate to.  Adding ParaView output to `seas_driver.cpp` is a substantial new feature out of scope for the io_dev compaction work.

Until ParaView output is wired into `seas_driver.cpp`:
- BP5 production runs through `seas_bp5_full` write probe-based station output and benchmark CSVs only — no volume PVD/VTU and no fault PVD/VTU.
- The `--no-volume-pv` / `--volume-pv-dt X` flags are accepted by `tpv102/104/205` only.
- When PV is added to `seas_driver.cpp`, the new defaults from plan §Phase 4 step 1 must be honoured at that time.

## ZFP lossy output

(Per plan §Phase 2d.3 + §Phase 6.  Requires `MFEM_USE_HDF5=YES` AND `MFEM_USE_H5Z_ZFP=YES` AND `HDF5_PLUGIN_PATH` set to the directory containing `libh5zzfp.{so,dylib}` at run time.)

Recommended defaults:

| Driver             | `--paraview-volume-zfp-tol`  | `--paraview-bulk-zfp-tol`    | `--paraview-fault-zfp-tol` |
|--------------------|------------------------------|------------------------------|----------------------------|
| seas (BP5)         | 1e-3 (m, displacement)       | n/a (no secondary collection) | 1e-12 (slip-rate floor)    |
| tpv102 / 104 / 205 | 1e-3 (primary velocity)      | 1e-3 (secondary stresses)    | 1e-12                      |

**Fault tolerance must be smaller than bulk/volume** because slip-rate spans 1e-9..1e0 m/s during nucleation→event; using a bulk-style 1e-3 on fault would erase six decades of dynamic range.  Bulk velocity, by contrast, spans ~3 decades and 1e-3 is invisible in ParaView.

### Phase 6 — volume PV via VTKHDF + ZFP (R-310 SEMANTIC FLIP)

As of plan §Phase 6, the **primary** volume `pv_dc_` member of `seas::ParaViewOutput` defaults to `mfem::ParaViewHDFDataCollection` (single `.vtkhdf` file per run) on `MFEM_USE_HDF5=YES` builds AND `MFEM_PARALLEL_HDF5=YES` for parallel ParMesh runs.  Driver flags:

- `--paraview-volume-vtu` — force the legacy per-rank VTU back end (Phase 1 default).  Use this for backwards compatibility with downstream scripts that hard-code `<output>/ParaView/<basename>_<rank>_<cycle>.vtu` paths.
- `--paraview-volume-hdf5` — explicit opt-in (rejected at parse time on non-HDF5 builds).
- `--paraview-volume-zfp-tol <tol>` / `--paraview-volume-deflate-level <N>` — apply compression to the **primary** volume collection (e.g., velocity in TPV* / displacement in BP5).  These flags are no-ops in `--paraview-volume-vtu` mode (warn-and-ignore on rank 0).

**`--paraview-bulk-*` semantics changed in Phase 6.4**: these flags now route to the **secondary** `pv_bulk_out` collection (TPV* stress fields), NOT the primary volume collection.  Pre-Phase-6 the bulk flags were a one-time rank-0 warning on the primary collection; that wiring has been removed.  BP5 has no secondary collection, so `--paraview-bulk-*` on BP5 is a one-time warning that the flag has no effect.

Cross-driver-uniform default collection name: `"volume"` for the primary collection (per R-305).  Files land at `<output_prefix>/volume.vtkhdf` (HDF5 mode) or `<output_prefix>/ParaView/volume/...vtu` (VTU mode).

## Pre-submit output-size estimator

Phase 7 ships `miniapps/seas/scripts/estimate_output_size.py`, a standard-library Python tool that estimates `du -sb $OUTPUT_DIR` before submitting a job.  Quick check:

```
python3 miniapps/seas/scripts/estimate_output_size.py \
    --driver bp5 \
    --inline-mesh \
    --tfinal 250yr \
    --paraview \
    --paraview-fault-zfp-tol 1e-12 \
    --paraview-max-snapshots 5000 \
    --no-volume-pv \
    --np 800 \
    --scratch-quota 1TB
```

The estimator parses the same `--paraview-*` flags the C++ drivers accept, integrates the `AdaptiveSchedule` over `tfinal`, and multiplies per-write bytes by the driver's registered-field schema.  Compression ratios marked `PLACEHOLDER` in `_io_size_compression.py` produce a one-time `UserWarning` per `(field_kind, filter)` pair; suppress with `--quiet`.

Tests: `pytest miniapps/seas/scripts/test_estimate_output_size.py` (28 currently passing; 3 reference-run skips pending Frontera replication).  The drift-detection test re-grep's the four driver source files for `RegisterDomainField("name", ...)` / `RegisterField("name", ...)` and fails when the source registers a string-literal field name not in `_io_size_schemas.py`.

## Spatial driver nucleation mechanism

`seas_spatial_dyn_driver` (the SAFS dynamic-rupture driver) supports nucleation via TWO orthogonal mechanisms:

1. **Static initial-stress patches** (`[stress] kind = "constant_tensor_with_patches"`).  Background constant Cauchy tensor + one or more `[[stress.patch]]` rectangular overrides applied at t = 0.  This is the spec-exact path for TPV205 — the patches are baked into `tau_pre_` by `FaultGeometry::ComputeSAFSParams<StressSource>` and have ZERO time dependence.  Per-component overrides (NaN sentinel = "inherit background"); last-match-wins on overlap.  Use this when the benchmark's nucleation is part of the *initial condition*, not a time-dependent perturbation.

2. **Time-dependent shear-stress perturbations** (`[nucleation]` block).  Per-DOF Δτ accumulator that ramps smoothly from 0 to full amplitude over `[0, T_nuc_s]` using the SCEC smoothStep function.  Two shapes:
   - `kind = "gradual_overstress"`: Gaussian centred on `(center_*_m)` with e-fold radii `radius_dip_m`/`radius_strike_m`.
   - `kind = "square_overstress"`: indicator-function rectangular patches via `[[nucleation.square_overstress.patch]]` array (shares one `T_nuc_s`).
   - Per-sub-step apply writes `ΔS(t, Δt) · F(r) · Δτ` into `DOFData[i].tau1_nuc` (dip) and `tau2_nuc` (strike); summed over `[0, T_nuc_s]` the increments telescope to the full target.  At `t ≥ T_nuc_s` the accumulator is a no-op.
   - The temporal factor is `smoothStep(t, t0) = 0` for `t ≤ 0`, `exp(τ²/(t·(t − 2·t0)))` for `0 < t < t0`, `1` for `t ≥ t0` — `C∞` everywhere except `t = 0`.

3. **Instantaneous overstress** (`[nucleation] kind = "instantaneous_overstress_circular"`).  TPV31-style: circular cosine-tapered Δτ written ONCE into `DOFData::tau{1,2}_nuc` at init (no ramp), with per-DOF µ-scaling (`Δτ(r) · µ(point)/mu_ref_pa`).  Spec formula: full amplitude inside `r ≤ radius_inner_m`, cosine taper `0.5·(1+cos(π·(r-ri)/(ro-ri)))` between `ri` and `ro`, zero outside.  The driver reads µ via `MaterialField::EvalAt` at every fault DOF.

**TPV205 specifically** uses the static-patches path (3 patches at ±7.5 km and 0 along-strike, all at z = 7.5 km, sigma_xy = 81.6/78/62 MPa) — see `tpv205/configs/tpv205.toml`.  Encoding via `gradual_overstress` was previously documented as an "adequate approximation" but produced three documented deviations from spec (Gaussian vs square, ramp vs instantaneous, single patch vs three) and has been retired.

**TPV31 specifically** uses `instantaneous_overstress_circular` nucleation + `depth_proportional` stress + `depth_profile_1d` material + a `[[friction.slip_weakening.spatial]]` rule with the new `cohesion_grad_pa_per_m`/`cohesion_ref_depth_m`/`cohesion_floor_pa` linear-taper fields (for the spec's `C₀(y) = max(0, 0.000425 MPa/m · (2400 - y))` ramp).  See `tpv31/configs/tpv31.toml`.  Sign convention note: SEAS is compression-POSITIVE (matches the proven TPV102/TPV104/TPV205 native drivers — `config/tpv*_params.hpp::sigma_n = 120e6` Pa "(positive compression)").  The TPV31 spec uses compression-NEGATIVE; the toml flips the normal-component signs (`σ_xx, σ_zz: -60 → +60 MPa`) but keeps shear signs (right-lateral `σ_xz = +30` is positive in both conventions).

User-facing schema: `safs/project_7.0_alternative/document/spatial_friction_config_schema.md`.  Workflow runbook: `safs/project_7.0_alternative/debug_document/spatial_workflow_safs_runbook_2026-05-18.md`.

The native TPV104 (rate-state) and TPV205 (`drivers/tpv205_driver.cpp` — uses `TPV205Params::ComputeTau2_0_TPV205` directly) drivers continue to use their own nucleation paths verbatim — the spatial-driver schema does NOT touch the TPV*/BP5 byte-exact regression contract.

## Bi-material Riemann (Phase R)

Bi-material exact linearised Riemann solver per Pelties et al. 2012 / SeisSol now runs on every interior face of any `WaveOperator` constructed via the `(MaterialField, BoundaryConfig)` ctor.  Verified on TPV205 (homogeneous regression — relative parity ≤ 1.4e-12 on serial, ≤ 1.6e-10 on ParMesh np=4 via the relaxed `T-PHASEH-SCALAR-PARITY` gate) and exercised on a 1D bi-material smoke fixture (`test_phaser_dispatch_smoke` R.2.T-3: two-cell heterogeneous fixture produces per-face flux matrices that differ from the homogeneous case by > 1e-6 relative).  See `safs/project_7.0_alternative/document/PLAN_phase_R_exact_bimaterial_riemann_rev3.md` for the full plan and acceptance gates.

The scalar `WaveOperator(mesh, p, lambda, mu, rho, bc)` ctor is unchanged — every existing TPV / BP5 driver target produces BYTE-IDENTICAL output to the pre-Phase-R baseline (Commit b18bed5).  The bi-material code path is reached only when the heterogeneous ctor is used (`owned_flux_pool_ != nullptr`); the dispatch wraps the FOUR interior-face sites (RK4 local, RK4 shared, ADER local, ADER shared) while the fault-side imposed-state sites stay verbatim.

`MixedFluxMode != None` is mutually exclusive with the heterogeneous ctor (R.2 R-002) — both replace the same interior-face flux; combining them would double-modify it.  `SetMixedFluxMode` aborts with a precise message if the heterogeneous ctor was used.
