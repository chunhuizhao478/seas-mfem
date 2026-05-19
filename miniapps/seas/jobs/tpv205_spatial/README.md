# TPV205 — heterogeneous Riemann sbatch jobs

These sbatch scripts drive **SCEC TPV5 / TPV205** through
`seas_spatial_dyn_driver` (config-driven driver) so the run exercises
the **heterogeneous Riemann interior-face flux** — the purpose of branch
`feature/heterogeneous_riemann_solver`.

The native `seas_tpv205_driver` (byte-parity baseline) is unchanged and
lives at `jobs/tpv205/`.

## Scripts

| Script | Queue | Nodes × ranks | Wall | Mesh   | Target tfinal |
|--------|-------|---------------|------|--------|---------------|
| `tpv205_spatial_dyn_200m_p1_O2_dev.sbatch`    | development | 8 × 400 | 2 h  | 200 m | early-rupture smoke (`~0.2 s`) |
| `tpv205_spatial_dyn_200m_p1_O2_normal.sbatch` | normal      | 16 × 800 | 24 h | 200 m | `~10 s` |
| `tpv205_spatial_dyn_100m_p1_O2_normal.sbatch` | normal      | 64 × 3200 | 48 h | 100 m | `~6 s` (resume to reach 12 s) |

## What each script verifies

The post-run summary parses the log + station files and emits a
verdict block.  Three gates must all pass for `STATUS: PASS`:

1. **Heterogeneous Riemann engaged.**  Grep for the
   `[wave_operator] BimaterialFlux precomputation:` banner.  A
   `[wave] interior_flux=scalar` line is a regression — the scalar ctor
   was selected, which means the canonical TPV205 TOML re-introduced
   the byte-parity opt-out.  See code-fix R-001.
2. **ParaView output written.**  Either `fault.vtkhdf` exists under
   `<output_dir>/` or a `*.pvd` index is present.  A
   `ParaView output: OFF` banner is a regression — `paraview_enabled`
   not set in the TOML.  See code-fix R-004.
3. **TPV205 SCEC station traces written.**  Grep for
   `[stations] TPV205 station writer active`.  Missing means
   `[problem].tag != "tpv205"` in the loaded config.  See code-fix R-007.

## Prerequisites

- MFEM built with `MFEM_USE_PETSC=YES`, `MFEM_USE_MUMPS=YES`,
  `MFEM_USE_HDF5=YES`.
- Mesh files at:
  - `miniapps/seas/tpv205/mesh/tpv2053d_200m.msh` (200 m scripts)
  - `miniapps/seas/tpv205/mesh/tpv2053d_100m.msh` (100 m script)
  generated with **Gmsh v2.2 format**:
  ```
  gmsh -format msh22 -3 tpv205/mesh/tpv2053d_200m.geo \
                  -o tpv205/mesh/tpv2053d_200m.msh
  ```
  The driver hard-aborts at preflight with this exact command if the
  `.msh` is missing.  See code-fix R-006.

## Resuming a partial run (100 m only)

The TOML's `checkpoint_every_steps = 10000` writes a checkpoint
periodically (and a final checkpoint at the end of the wall-time
window).  To continue:

```
ibrun ./seas_spatial_dyn_driver \
      --config <RESULT_DIR>/tpv205_100m.toml \
      --output-dir <NEW_RESULT_DIR> \
      --restart <PREV_RESULT_DIR>/cp_<step> \
      --verify-dispatch
```

The `--output-dir` MUST differ from the parent of `--restart`
(the driver hard-fails on overlap to prevent overwrite).

## Spec encoding (verbatim per `tpv205/configs/tpv205.toml`)

| Quantity | Value | Source |
|----------|-------|--------|
| ρ        | 2670 kg/m³ | `[material_constant_fallback]` |
| c_p      | 6000 m/s | derived from `lambda + 2μ = ρ·c_p²` |
| c_s      | 3464 m/s | derived from `μ = ρ·c_s²` |
| σ_n      | 120 MPa | `[stress].sigma_yy_pa` |
| τ_back   | 70 MPa  | `[stress].sigma_xy_pa` |
| τ_nuc    | 81.6 MPa | nucleation `[[stress.patch]]` at (0, 0, −7500) |
| τ_left   | 78.0 MPa | left `[[stress.patch]]` at (−7500, 0, −7500) |
| τ_right  | 62.0 MPa | right `[[stress.patch]]` at (+7500, 0, −7500) |
| μ_s      | 0.677 (rupture area) / 1e6 (barrier) | `mu_s_default` + 3 barrier rules |
| μ_d      | 0.525 | `mu_d_default` |
| d_c      | 0.40 m | `d_c_default` |

Rupture initiates spontaneously at t = 0+ inside the central patch
(μ_s·σ_n = 81.24 MPa < 81.6 MPa).
