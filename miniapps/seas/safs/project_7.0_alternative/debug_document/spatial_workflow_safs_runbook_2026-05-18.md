# SAFS Dynamic-Rupture Workflow Runbook (2026-05-18)

**Status:** living runbook for the constant-material SAFS dynamic-rupture path.
**Driver:** `seas_spatial_dyn_driver` (`drivers/spatial_dyn_driver.cpp`).
**Friction:** linear slip-weakening (LSW) with per-DOF spatial heterogeneity.
**Nucleation:** `gradual_overstress` (smoothStep-ramped, Gaussian-shaped per-DOF δτ accumulator).
**Bulk material:** constant `(λ, μ, ρ)` from `[material_constant_fallback]`; CVM-based heterogeneous material is deferred to a separate branch (Phase R bi-material Riemann solver).
**Verification:** the first concrete dataset run end-to-end through this workflow is a TPV205-derived TOML; details in `spatial_workflow_tpv205_verification_2026-05-18.md` (companion debug doc, to be written).

This runbook is the single user-facing entry point for "I have a new region / dataset; how do I run a dynamic-rupture simulation on it?" Each section names the **dataset producer**, the **TOML key it populates**, and the **validation step** that closes the loop.

---

## 1. Pipeline overview

The driver is a four-layer data-driven pipeline. Each layer has a producer (Python or external) that emits an artifact; the TOML config tells the driver where to find the artifact; a parser/resolver consumes it into per-DOF arrays; the bulk physics is fixed and shared.

```
                ┌──────────────────────────────────────────────────────────┐
                │                                                          │
   raw data ───►│  producer (Python / external)  ──►  artifact on disk    │
                │                                                          │
                └──────────────────┬───────────────────────────────────────┘
                                   │  TOML [layer] block names the artifact
                                   ▼
                ┌──────────────────────────────────────────────────────────┐
                │  spatial driver loader + resolver                        │
                │   - parse TOML                                           │
                │   - load artifact                                        │
                │   - project / interpolate / rotate to per-DOF arrays     │
                └──────────────────┬───────────────────────────────────────┘
                                   │  per-DOF Vectors
                                   ▼
                ┌──────────────────────────────────────────────────────────┐
                │  FaultGeometry + FaultFaceFlux::DOFData[i]               │
                │   (shared by every nucleation kind, every friction law)  │
                └──────────────────┬───────────────────────────────────────┘
                                   │
                                   ▼
                ┌──────────────────────────────────────────────────────────┐
                │  WaveOperator + Tpv205SubStepIterator (ADER time-stepping)│
                │  → ParaView VTKHDF + V1 checkpoint                       │
                └──────────────────────────────────────────────────────────┘
```

| Layer       | TOML block          | Producer                                                                  | Artifact                                                                |
|-------------|---------------------|---------------------------------------------------------------------------|-------------------------------------------------------------------------|
| Mesh        | `[mesh]`            | `meshing/code/run_nwcut_meshing.py`                                       | `meshing/results/msh/*.msh` (Gmsh v2.2)                                  |
| Material    | `[material_constant_fallback]` (CVM sidecar deferred)                       | `velocity/code/build_velocity_<model>.py`                              | `velocity/results/<model>/velocity_safs.h5` (CVM sidecar; not yet consumed) |
| Stress      | `[stress]`          | `stress/code/build_stress_csm_safs.py` OR inline TOML scalars             | `stress/results/*.h5` OR none (constant tensor)                          |
| Friction    | `[friction.slip_weakening]` + `[[…spatial]]`                                | hand-authored TOML referencing `friction/slip-weakening/geoffrey2010.md` | TOML in-file                                                            |
| Nucleation  | `[nucleation]`      | hand-authored TOML (hypocenter location, Δτ amplitude, ramp duration)     | TOML in-file                                                            |
| Numerics    | `[numerics]`        | hand-authored                                                             | TOML in-file                                                            |
| Time        | `[time]`            | hand-authored                                                             | TOML in-file                                                            |
| Output      | `[output]`          | hand-authored                                                             | TOML in-file                                                            |

---

## 2. Step-by-step workflow

### Step 1 — Build the mesh

```bash
conda activate pythonenv  # Gmsh dependency
cd miniapps/seas/safs/project_7.0_alternative/meshing
python code/run_nwcut_meshing.py --resolution 1000 --lcfar 3000
# emits meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh
```

**Validation:**
- File exists; `wc -l` > 0.
- File begins with `$MeshFormat\n2.2 0 8` (Gmsh v2.2 required — `mfem::Mesh::ReadGmshMesh` is v2.2-only; see `miniapps/seas/CLAUDE.md` "Known limitation").
- Physical surface tags: fault = 101, top = 102, bottom = 103, sides = 104 (consumed by `BoundaryConfig` in the driver).
- For SAFS the bbox should contain the production volume in UTM 11N coordinates; sanity-check via `python -c "import meshio; m = meshio.read('…msh'); print(m.points.min(0), m.points.max(0))"`.

### Step 2 — (Optional today) Build the CVM velocity sidecar

```bash
conda activate pythonenv
cd miniapps/seas/safs/project_7.0_alternative/velocity
python code/build_velocity_cvmh.py            # produces results/cvmh/velocity_safs.h5
# OR equivalent for cvm_s4.26.m01, multiscale_statewise_cvm
```

**Validation:**
- `python code/bbox_check.py results/cvmh/velocity_safs.h5` confirms the sidecar bbox contains the mesh bbox.
- `h5ls results/cvmh/velocity_safs.h5` shows `lambda`, `mu`, `rho` datasets with consistent shapes.

**STATUS NOTE:** the sidecar is currently produced but NOT consumed by the spatial driver. `Mode::Coefficient` aborts at the `WaveOperator(MaterialField)` ctor (Phase H Stage 2 + Phase R bi-material Riemann not yet landed). The sidecar artifact is committed forward so it is ready when heterogeneous-material work merges from the separate branch. For all runs today, use `[material_constant_fallback]` and pass `--no-sidecar-material`.

### Step 3 — Build the regional stress sidecar (CSM projection)

```bash
conda activate pythonenv
cd miniapps/seas/safs/project_7.0_alternative/stress
python code/build_stress_csm_safs.py
# emits stress/results/stress_csm_safs.h5
python code/verify_onfault_stress.py results/stress_csm_safs.h5
# Prints: per-DOF τ_pre magnitude histogram, σ_n_eff histogram,
#         |τ_pre| / (μ_s · σ_n_eff) ratio histogram (must be < 1
#         everywhere outside the planned nucleation patch).
```

**Validation:**
- `verify_onfault_stress.py` exits 0.
- Sign conventions confirmed: σ_n > 0 = compression; τ_pre direction parallel to expected loading direction; tangent basis `t1 = dip, t2 = strike` (Tandem / FaultBasis convention; CLAUDE.md).
- If using `[stress] kind = "constant_tensor"` instead (testing primitive, D-1), the verification reduces to: are the six `sigma_*_pa` values consistent with the regional stress regime you intend? Use a single-rule TOML spreadsheet to map (S_Hmax magnitude, S_Hmax azimuth, S_v, S_hmin) → six `σ_*` components, and document the mapping in your run's TOML `description` field.

### Step 4 — Author the TOML config

Start from `friction/EXAMPLE_spatial_friction_slip_weakening_safs.toml` and edit. The full schema reference is `document/spatial_friction_config_schema.md`. The minimum-viable block set for a SAFS LSW run:

```toml
[meta]
schema_version = 1
law            = "slip_weakening"
description    = "<region>-<scenario> LSW, geoffrey2010 reference parameters, gradual_overstress nucleation"

[material_constant_fallback]
lambda = 32.0e9
mu     = 32.0e9
rho    = 2670.0

[pore_pressure]
P_p_pa            = 0.0
P_p_grad_pa_per_m = 0.0
min_sigma_n_pa    = 0.0

[mesh]
path  = "safs/project_7.0_alternative/meshing/results/msh/safs_fault_box_nwcut_1000m_lcfar3000.msh"
order = 1

[velocity]
# Carried for forward compatibility; ignored under --no-sidecar-material.
model         = "cvmh"
dataset_root  = "safs/project_7.0_alternative"
override_path = ""

[stress]
kind = "constant_tensor"   # or "sidecar_hdf5" + sidecar_path = "..."
sigma_xx_pa = 80.0e6
sigma_yy_pa = 80.0e6
sigma_zz_pa = 80.0e6
sigma_xy_pa = 20.0e6
sigma_yz_pa =  0.0
sigma_xz_pa =  0.0

[numerics]
ader_order = 2
mixed_flux = "none"
cfl        = 0.5
use_pml    = false

[time]
tfinal     = "12s"
t_initial  = 0.0
dt_initial = "auto"
dt_max     = "0.1s"

[output]
output_dir              = "output_<region>_<scenario>"
restart_prefix          = "cp"
paraview_volume         = "off"      # default-OFF per ParaView parity plan
paraview_bulk           = "off"
paraview_fault          = "hdf5"
paraview_volume_dt      = "0.05s"
paraview_bulk_dt        = "0.05s"
paraview_fault_dt       = "0.001s"
paraview_volume_zfp_tol = 1.0e-3
paraview_bulk_zfp_tol   = 1.0e-3
paraview_fault_zfp_tol  = 1.0e-12
max_snapshots           = 5000
checkpoint_every_steps  = 10000

[nucleation]
kind = "gradual_overstress"
[nucleation.gradual_overstress]
center_x_m       =  480000.0    # hypocenter UTM 11N x
center_y_m       = 3750000.0    # hypocenter UTM 11N y
center_z_m       =   -7000.0    # hypocenter depth (z<0 below free surface)
radius_dip_m     =     3000.0   # Gaussian e-fold radius down-dip
radius_strike_m  =     3000.0   # Gaussian e-fold radius along-strike
delta_tau_dip_pa =        0.0   # full-ramp Δτ in dip direction
delta_tau_strike_pa = 25.0e6    # full-ramp Δτ in strike direction
T_nuc_s          =        1.0   # ramp completes at t = T_nuc_s
t0_smooth_s      =        0.5   # smoothStep midpoint (typically T_nuc_s / 2)

[friction.slip_weakening]
mu_s_default     = 1.1
mu_d_default     = 0.5
d_c_default      = 0.5
cohesion_default = 0.0

[[friction.slip_weakening.spatial]]
kind    = "barrier"
z_min_m = -20000.0
z_max_m = -15000.0
```

### Step 5 — Dry-run

```bash
conda activate mfem-dev   # mpicxx, MFEM, HDF5
cd miniapps/seas
make seas_spatial_dyn_driver
./seas_spatial_dyn_driver --config <your.toml> --no-sidecar-material --dry-run --print-derived
```

**What `--dry-run --print-derived` produces (and what to check):**
- TOML parse OK, every validator gate green.
- Mesh load OK; bbox + element count printed.
- Boundary attribute mapping OK (`fault=101, top=102, bottom=103, sides=104`).
- `MaterialField::Mode::Constant` from `[material_constant_fallback]`.
- Stress source applied; per-DOF τ_pre and σ_n_eff histograms printed.
- Friction resolver: per-region histograms of `μ_s`, `μ_d`, `d_c`.
- Per-DOF `L_nuc = μ · d_c / ((μ_s − μ_d) · σ_n_eff)`: assert ≥ 10 mesh elements in the nucleation patch.
- Per-DOF `|τ_pre| / (μ_s · σ_n_eff)`: must be < 1 *everywhere outside the planned nucleation patch*. If any DOF outside the patch exceeds this, the initial conditions will spontaneously rupture and your run is ill-posed before you start.
- Per-DOF gradual-overstress amplitude `F(r) · |Δτ|`: max should equal `|Δτ|` at the hypocenter and decay smoothly to zero outside `~ 3 · radius_*`.
- CFL Δt min/max from heterogeneous-CFL pass (Phase H Stage 1) — under constant material these are identical and equal `cfl · h_min / cp`.
- Max-shear direction + σ_1 azimuth (sanity vs. input dataset documentation).
- Exit code 0.

**If `--print-derived` reports any of:**
- `|τ_pre| / (μ_s · σ_n_eff) > 1` outside the nucleation patch → adjust stress tensor / friction defaults so the fault is sub-critical at t = 0.
- `L_nuc < 10` elements → refine mesh or increase `d_c`.
- `F(r) · |Δτ|` peak inside the asperity does NOT exceed the local `(μ_s − μ_d) · σ_n_eff` overshoot needed to nucleate → increase `delta_tau_*_pa` until it does.
- CFL Δt too small to fit `[time].tfinal` in your Frontera wall budget → increase `cfl` (max 0.5 by validation rule), or coarsen the mesh.

### Step 6 — Local smoke run (8 ranks, fast tfinal)

```bash
# Edit your TOML temporarily: tfinal = "0.5s", paraview_fault_dt = "0.005s"
mpirun -np 8 ./seas_spatial_dyn_driver --config <your.toml> --no-sidecar-material
# Run completes in ~minutes on a workstation.  Output goes to <output_dir>/.
```

**Validation (post-smoke):**
- Exit code 0.
- `<output_dir>/fault.vtkhdf` exists, opens in ParaView 5.11+.
- Max slip-rate in the rupture core ≥ 1e-3 m/s by t ≈ T_nuc_s (rupture has begun).
- Visual check: rupture front propagates radially outward from the hypocenter at sub-Rayleigh speed; no spurious nucleation elsewhere.
- No NaN in any output field; σ_n > 0 everywhere; V_max peaks then decays.

### Step 7 — Production run (Frontera)

```bash
cd miniapps/seas
sbatch jobs/safs/spatial_dyn_production_normal_48hr_safs.sbatch
# (sbatch + verify_spatial_dyn_smoke_safs.py to be added in Phase 6)
```

**Resource budget targets (for a 1000 m SAFS mesh × 12 s tfinal):**
- 32 nodes × 50 ranks/node = 1600 ranks (Frontera normal queue).
- ~24-48 hr wall depending on `cfl` and per-step cost.
- V1 checkpoint every 10000 steps; ParaView fault HDF every 0.001 s; volume/bulk OFF.
- Expected output volume: 50-200 GB (fault HDF only with ZFP at tol 1e-12).

### Step 8 — Restart (long runs)

```bash
mpirun -np 1600 ./seas_spatial_dyn_driver --config <your.toml> --no-sidecar-material \
    --restart <output_dir>/cp
```

The driver verifies `driver_tag == "spatial_dyn"` in the checkpoint header (Phase 5b) and aborts on mismatch. Pre-tag checkpoints are accepted with a one-line warning.

### Step 9 — Post-run verification

```bash
# Per-station traces from fault.vtkhdf
python miniapps/seas/safs/project_7.0_alternative/spatial/code/scripts/verify_spatial_dyn_smoke_safs.py \
    <output_dir>/fault.vtkhdf  # Phase 6 of spatial_dynamic_rupture_plan.md
```

**Standard accuracy claims (the hierarchy you build to defend the result):**
1. **Sub-rung gates** — initial equilibrium (Step 5), smoke OK (Step 6), no NaN / no V_max growth.
2. **Mesh convergence** — repeat the run at 500 m mesh; slip-rate at a small set of named fault stations agrees within ~5% L∞ with the 1000 m run.
3. **Cross-code reference** — for SAFS, the closest match is to run a TPV205-derived configuration through this driver (see `spatial_workflow_tpv205_verification_2026-05-18.md`) and confirm the rupture-front arrival time at TPV205's named stations agrees with native `seas_tpv205_driver` to within `T_nuc_s/2 + 5%`. The `T_nuc_s/2` offset is the expected difference between gradual-overstress (ramped) and TPV205-native (instantaneous) nucleation; this offset must be documented in your run log.
4. **Independent code** — if you can run the same dataset through Tandem or SeisSol, the slip-rate at named stations should agree within a few percent. This is the strongest external accuracy claim.

---

## 3. Boundary conditions, sign conventions, and gotchas

### Boundary attribute mapping

The SAFS `.msh` files emit four physical surface tags:

| Tag | Surface       | Driver BC                                         |
|-----|---------------|---------------------------------------------------|
| 101 | fault         | dynamic-rupture interface (`FaultFaceFlux`)       |
| 102 | top           | free surface                                      |
| 103 | bottom        | absorbing                                          |
| 104 | sides (×4)    | absorbing                                          |

Confirm via `--dry-run` log line `"BoundaryConfig: fault=101, top=102, bottom=103, sides=104"`. If your mesh uses different tags, edit `[mesh]` accordingly (TOML key not yet wired — currently the mapping is hard-coded; see review finding R-107).

### Sign conventions (project-wide; see `miniapps/seas/CLAUDE.md`)

- σ_n > 0 = compression (geological convention).
- Depth: z = 0 at surface, z < 0 below; `min(0, -z_dof)` is the absolute depth.
- Dip direction: `(0, 0, +1)` = downward into earth.
- Fault tangent frame: `t1 = dip, t2 = strike` (Tandem FaultBasis).
- Pre-stress direction: parallel to expected loading velocity, NOT antiparallel (the parallel choice avoids positive feedback in the friction solver).
- Cauchy stress tensor in `[stress] kind = "constant_tensor"`: EAST-NORTH-UP frame, compression positive (matches `stress/code/hickman_and_zoback_*`).

### Numerical gotchas

- **CFL collapse:** under constant material today, CFL is uniform; once CVM sidecars are wired, soft basin elements will drive a much smaller `Δt`. Inspect the `--print-derived` `per_elem_h_/cp_max` histogram.
- **Brent friction solver:** required (Newton fails with large ψ). Already wired in `dieterich_ruina.hpp`; spatial driver inherits.
- **Per-DOF interleaving:** in `DOFData[i]`, slot 1 = dip, slot 2 = strike. The gradual-overstress accumulator writes both `tau1_nuc` (dip) and `tau2_nuc` (strike); legacy TPV104 wrote only `tau2_nuc` (strike-slip only). The SAFS curvilinear fault has variable strike, so both components are non-zero in general.
- **MPI rank divergence in RK error:** `MPI_Allreduce(MPI_MAX)` already in place in `time_stepper.hpp`; spatial driver inherits.

### Output sizing

Default for production runs (per ParaView parity plan):
- `paraview_volume = "off"`, `paraview_bulk = "off"`, `paraview_fault = "hdf5"`.
- This is a default-OFF posture: a SAFS 1000 m mesh × 0.001 s fault cadence × 12 s with volume+bulk ON would emit multi-TB output.
- Opt in per collection only when you have a specific question (debugging a particular waveform; producing a paper figure).
- Use `python scripts/estimate_output_size.py --driver spatial_dyn …` pre-submit to estimate `du -sb` of `<output_dir>` before launching.

---

## 4. Layer-by-layer "is this dataset ready?" checklist

Print this and tick boxes before launching a new dataset on Frontera.

### Mesh

- [ ] `.msh` exists in `meshing/results/msh/`
- [ ] File starts with `$MeshFormat\n2.2 0 8`
- [ ] Physical surface tags 101 / 102 / 103 / 104 present
- [ ] Bounding box covers the production region (UTM check)
- [ ] Element count and `h_min` match the resolution claim

### Material (constant fallback)

- [ ] `[material_constant_fallback]` lambda, mu, rho documented (which region average)
- [ ] `--no-sidecar-material` set on the CLI (or in the run script)
- [ ] (Future heterogeneous-material work: separate branch; not required here)

### Stress

- [ ] `kind = "constant_tensor"` OR `kind = "sidecar_hdf5"` chosen with documented rationale
- [ ] If sidecar: `verify_onfault_stress.py` exit 0 on the production sidecar
- [ ] If constant tensor: six `σ_*` values traceable to a regional-stress reference (Hickman & Zoback, or equivalent)
- [ ] `--dry-run` reports `|τ_pre| / (μ_s · σ_n_eff) < 1` outside the nucleation patch
- [ ] σ_1 azimuth and max-shear direction match input dataset documentation

### Friction

- [ ] `mu_s_default / mu_d_default / d_c_default` traceable to `friction/slip-weakening/geoffrey2010.md` or another cited reference
- [ ] Spatial rules document order is intentional (last-match wins per key)
- [ ] No `mu_s > 1.0e5` typed by a user (use `kind = "barrier"` for locked DOFs)
- [ ] `--dry-run` reports per-region `(μ_s, μ_d, d_c)` histograms consistent with the rules

### Nucleation

- [ ] `kind = "gradual_overstress"` selected
- [ ] `center_*_m` hypocenter coordinates in the mesh bbox
- [ ] `radius_*_m` ≥ 3 × element edge length (so the Gaussian is resolved)
- [ ] `delta_tau_*_pa` chosen so `F(r=0) · |Δτ|` exceeds the local nucleation budget `(μ_s − μ_d) · σ_n_eff`
- [ ] `T_nuc_s` ≤ `tfinal / 10` (ramp completes well before the run ends)
- [ ] `t0_smooth_s` typically `≈ T_nuc_s / 2`
- [ ] `--dry-run` reports a non-empty set of fault DOFs receiving the accumulator

### Numerics + time + output

- [ ] `cfl ≤ 0.5`, `ader_order = 2` (production default)
- [ ] `tfinal` parses with the `s` suffix
- [ ] `paraview_volume`/`paraview_bulk` default OFF; only opt in with a specific question
- [ ] `checkpoint_every_steps` chosen so a checkpoint covers ≤ 30 min wall (so a Frontera-job restart loses < 30 min progress)
- [ ] `python scripts/estimate_output_size.py` estimate fits the scratch quota

---

## 5. Common failure modes and remediation

| Symptom                                                       | Root cause                                                                                  | Remediation                                                                |
|---------------------------------------------------------------|---------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| Ctor abort `"Phase H Stage 2 not yet wired"`                  | TOML sets `velocity.model = "cvmh"` without `--no-sidecar-material`                          | Add `--no-sidecar-material` or set the model to constant-fallback path     |
| `mfem::Mesh::ReadGmshMesh: vertices indices are not unique`   | `.msh` is Gmsh v4.x; MFEM parser is v2.2-only                                                | Re-emit with `gmsh ... -format msh22 -o file.msh` (see CLAUDE.md "Known limitation") |
| Spontaneous rupture from t=0 outside the nucleation patch     | `|τ_pre|` exceeds `μ_s · σ_n_eff` somewhere outside the patch                                | Increase `mu_s_default`, or reduce stress magnitude, or refine spatial rules |
| Rupture never nucleates                                       | `F(r=0) · |Δτ|` insufficient relative to `(μ_s − μ_d) · σ_n_eff` budget                       | Increase `delta_tau_*_pa` or refine `radius_*_m`                            |
| `V_max` monotonically increasing                              | Pre-stress sign error (slip-rate antiparallel to traction)                                  | Confirm `tau_pre || V_init` per CLAUDE.md sign convention                   |
| `dt → 0` during nucleation                                    | Initial Δt too large; RK45 stage amplification                                              | Set `dt_initial = "auto"` (resolves to `0.01 · L_nuc / V_nuc`)              |
| `MFEM_VERIFY: F(a) and F(b) same sign in Brent`               | σ_n turned negative somewhere (over-effective-pressure)                                     | Add `min_sigma_n_pa = 1e5` clamp; check pore pressure block                |
| Restart aborts with `"driver_tag mismatch"`                   | Checkpoint was produced by a different driver (TPV104, BP5, etc.)                            | Use a checkpoint produced by `seas_spatial_dyn_driver`                      |
| Frontera job runs out of disk quota mid-run                   | `paraview_volume = "hdf5"` on a SAFS-scale mesh × ms fault cadence                          | Set volume/bulk OFF; run `estimate_output_size.py` before submit            |

---

## 6. What's not yet in this workflow

The following are explicitly deferred (heterogeneous-material work, separate branch):

- CVM sidecar consumption: `Mode::Coefficient` aborts at `WaveOperator(MaterialField)` ctor today. Once Phase H Stage 2 + Phase R bi-material Riemann land, drop `--no-sidecar-material` and the same TOML runs with heterogeneous bulk material.
- Per-fault-QP impedance (`InitializeImpedancesPerQP`): a no-op under constant material; needed once `(λ, μ, ρ)` varies along the fault.
- Two-CVM comparison (`cvm_s4.26.m01` vs `multiscale_statewise_cvm`): driver CLI already accepts both names; gated on Phase H Stage 2.

The following are in scope for THIS workflow but not yet implemented:

- The `gradual_overstress` resolver itself (`spatial_nucleation.{hpp,cpp}`, `NucleationKind::GradualOverstress` enum, TOML sub-block parser, per-sub-step accumulator wired into `Tpv205SubStepIterator`). This is the immediate next code change; see the plan revisions in `spatial_dynamic_rupture_plan.md` and `spatial_friction_config_schema.md` for the interface.
- ParaView parity plan Phases 1–6 (CLI flag surface, build guards, default-OFF, secondary stress collection, volume velocity, SAFS fault statics including `nuc_amplitude` and `nuc_radial_factor`).
- Phase 6 sbatch + `verify_spatial_dyn_smoke_safs.py`.

The following is a verification companion to this runbook:

- `spatial_workflow_tpv205_verification_2026-05-18.md` — the TPV205-input validation procedure that closes the loop on "does this driver compute the right answer for a known-good setup?"

---

## 7. Quick-reference commands

```bash
# Build
conda activate mfem-dev
cd miniapps/seas
make seas_spatial_dyn_driver

# Dry-run (~seconds)
./seas_spatial_dyn_driver --config <PATH.toml> --no-sidecar-material --dry-run --print-derived

# Local smoke (8 ranks, edit tfinal to "0.5s" first)
mpirun -np 8 ./seas_spatial_dyn_driver --config <PATH.toml> --no-sidecar-material

# Production (Frontera, will land in Phase 6)
sbatch jobs/safs/spatial_dyn_production_normal_48hr_safs.sbatch

# Restart
mpirun -np N ./seas_spatial_dyn_driver --config <PATH.toml> --no-sidecar-material --restart <output_dir>/cp

# Output size estimate
python scripts/estimate_output_size.py --driver spatial_dyn --inline-mesh \
    --tfinal 12s --paraview --paraview-fault-zfp-tol 1e-12 \
    --paraview-max-snapshots 5000 --no-volume-pv --np 1600 --scratch-quota 1TB
```
