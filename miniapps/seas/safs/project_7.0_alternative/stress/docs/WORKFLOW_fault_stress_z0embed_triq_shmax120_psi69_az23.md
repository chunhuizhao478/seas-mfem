# WORKFLOW: regional stress → on-fault traction VTU (triq mesh, Ψ = 69°)

**Target artefact (this doc explains exactly how it is produced):**

```
stress/results/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax120_psi69_az23/
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax120_psi69_az23_fault_stress.vtu
    └ siblings: ..._bulk_stress.vtu, ..._summary.json, verify_report.json
```

Written 2026-06-09. This is the **redo** of the `hz1671_psi37_az351` run
(`WORKFLOW_fault_stress_z0embed_hz1671_psi37_az351.md`) with three changes:

1. **Ψ = 69°** between SHmax and the local SAF strike (was 37°), i.e.
   `SHmax_az = 314 + 69 ≡ 23°` cw-from-N — the H&Z deep-SAF orientation
   (same azimuth as `demo_safod()`).
2. **New magnitudes**: `SHmax = 120`, `Shmin = 55`, `Sv = 50`, `P_p = 20` MPa
   (was the H&Z 1671 m set 113 / 49 / 45 / 16).
3. **New input mesh**: the `triq` (triangle-quality remeshed) z0embed mesh
   `experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh`
   (60 658 fault triangles / 1 164 469 tets, vs 42 358 / 3 693 060 for the
   meshing-pipeline z0embed mesh).

Read `stress/README.md` §Workflow for the general pipeline; read this for the
exact commands and the numbers that came out.

---

## TL;DR — the commands that produced the target files

```bash
conda activate pythonenv

# Stage 0 — split the .msh into fault/bulk VTUs (the triq mesh lives in
# experimental_mesh_refinement/, not meshing/results/, so the split is
# emitted next to the .msh):
cd experimental_mesh_refinement
python ../meshing/code/msh_to_vtu.py \
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh \
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq
# → ..._triq_fault.vtu (6.1 MB), ..._triq_bulk.vtu (35.1 MB)

# Stages 1–5 — projection with ALL five regional parameters overridden:
cd ../stress/code
python project_to_fault_stress.py \
    ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax120_psi69_az23 \
    --write-bulk \
    --SHmax 120 --Shmin 55 --Sv 50 --P_p 20 \
    --SHmax-az 23
```

- **arg 1** = input fault VTU (Stage 0 output).
- **arg 2** = `mesh_base` — output stem *and* output sub-directory name
  (`_extract_lc_tag()` falls back to the full base for long stems, same as
  the psi37 run).
- `--strike-hint-az` is left at the default **314°** (SAF strikes N46°W).

---

## Filename tag decode

`safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq` `_shmax120` `_psi69` `_az23`

| Segment | Meaning | Source |
|---|---|---|
| `safs_fault_box_nwcut_500m_lcfar3000_z0embed` | same box/cut/resolution/far-field/z0-cut lineage as the psi37 run | meshing pipeline |
| `triq` | triangle-quality remeshed variant of the z0embed mesh | `experimental_mesh_refinement/` |
| **`shmax120`** | magnitude set `SHmax=120, Shmin=55, Sv=50, P_p=20` MPa (full set echoed in `_summary.json["params"]`) | CLI overrides |
| **`psi69`** | **Ψ = 69°** — angle between SHmax and the SAF strike hint (`|23 − 314| mod 180 = 69`). The *deep* H&Z value (weak-fault geometry), vs the shallow Ψ = 37° of the previous run | `--SHmax-az 23` |
| **`az23`** | **SHmax azimuth = 23°** cw-from-N (= N23°E, the `demo_safod()` orientation) | `--SHmax-az 23` |

The tags are operator-supplied labels baked into `mesh_base`; the script does
**not** parse them — the actual numeric parameters are the ones echoed in
`_summary.json["params"]`.

---

## Stage 1 — regional stress parameters

From `_summary.json["params"]` of this run:

| Parameter | Value | psi37 run | Note |
|---|---|---|---|
| `SHmax_MPa` | **120.0** | 113.0 | max horizontal principal stress |
| `Shmin_MPa` | **55.0** | 49.0 | min horizontal principal stress |
| `Sv_MPa` | **50.0** | 45.0 | vertical stress |
| `P_p_MPa` | **20.0** | 16.0 | pore pressure |
| `SHmax_azimuth_deg` | **23.0** | 351.0 | Ψ = 69° against the 314° strike hint |
| `fault_strike_azimuth_hint_deg` | 314.0 | 314.0 | SAF strikes N46°W (default) |
| `depth_model` | `constant` | `constant` | σ⁰ depth-independent; all gradients = 0 |
| `rake_sense` | right-lateral | right-lateral | |

`Sv = 50 < Shmin = 55 < SHmax = 120` — still a transitional SS↔reverse regime
(Sv is the least principal stress), as in the psi37 run.

---

## Stage 2 — bulk Cauchy tensor σ⁰ at z = 0 (compression POSITIVE)

Derivation — convert the principal stress direction/magnitude [SHmax, Shmin]
to the [East, North, Up] frame:

[1] map azimuth angle to math angle:

```
alpha = (90 - azimuth) mod 360 = (90 - 23) mod 360 = 67 deg
```

(SHmax points N23E geographically = 67 deg ccw from East in math convention.)

[2] principal stress tensor (compression negative, H&Z source convention):

```
Sigma'_HZ = diag(-SHmax, -Shmin, -Sv) = diag(-120, -55, -50)  MPa
```

[3] rotation matrix: c = cos(67) = 0.39073, s = sin(67) = 0.92050

```
        [ c  -s   0 ]   [ 0.39073  -0.92050   0 ]
Rz(a) = [ s   c   0 ] = [ 0.92050   0.39073   0 ]
        [ 0   0   1 ]   [ 0         0         1 ]
```

[4] rotate into [East, North, Up] and flip sign (compression positive):

```
sigma0 = -sigma0_HZ = Rz * diag(120, 55, 50) * Rz^T

sigma_EE = SHmax*c^2 + Shmin*s^2 = 120*0.15267 + 55*0.84733 =  64.924
sigma_NN = SHmax*s^2 + Shmin*c^2 = 120*0.84733 + 55*0.15267 = 110.076
sigma_EN = (SHmax - Shmin)*c*s  = 65*0.35967                =  23.379
```

[5] stress components in the [East, North, Up] frame:

```
         [  64.924   23.379    0.0 ]
sigma0 = [  23.379  110.076    0.0 ]   MPa  (compression POSITIVE)
         [   0.0      0.0     50.0 ]
```

Matches the `sigma_global_at_z0_MPa` echo in `_summary.json` to round-off.
(The psi37/az351 run had `sigma_EN = -9.889`; the sign flip is the SHmax axis
moving from N9°W — math angle 99°, where c·s < 0 — to N23°E — math angle 67°,
where c·s > 0.)

---

## Stages 3–5 — identical machinery to the psi37 run

Per-triangle Tandem basis (`s = up × n̂`, `d = s × n̂`, strike harmonised to the
314° hint), `resolve_traction` per cell, area-weighted node averaging, then
fault/bulk VTU + summary writers. For this run:

- `fault_n_cells = 60658`, `fault_n_degenerate_basis = 0`
- `bulk_n_cells = 1164469`
- `fault_vtu_nan_vertex_count = 0`

Outputs in `results/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax120_psi69_az23/`:

- `…_fault_stress.vtu` (14 MB) — cell data `sigma_n_total_MPa_cell`,
  `sigma_n_eff_MPa_cell`, `tau_strike_MPa_cell`, `tau_dip_MPa_cell`, per-cell
  basis vectors + traction; point data (area-weighted node averages)
  `sigma_n_total_MPa`, `sigma_n_eff_MPa`, `tau_strike_MPa`, `tau_dip_MPa`,
  `tau_magnitude_MPa`, `rake_deg`, `mu_apparent`.
- `…_bulk_stress.vtu` (27 MB) — six-component `sigma_<ij>_MPa` cell data at
  every tet centroid (constant tensor for `depth_model = constant`).
- `…_summary.json` — parameter echo, σ⁰ at z = 0, counts, per-field stats.

---

## This run's results (`_summary.json["stats"]`)

| Field | min | median | max | psi37 median |
|---|---:|---:|---:|---:|
| `sigma_n_total_MPa` | 53.89 | 98.18 | 120.00 | 73.72 |
| `sigma_n_eff_MPa`   | 33.89 | 78.18 | 100.00 | 57.72 |
| `tau_strike_MPa`    | −31.96 | 3.81 | 28.52 | 26.47 |
| `tau_dip_MPa`       | −35.00 | −27.29 | 11.73 | −16.01 |
| `tau_magnitude_MPa` | 0.08 | 31.68 | 35.00 | 30.61 |
| `rake_deg`          | −179.90 | −76.15 | 179.76 | −29.54 |
| `mu_apparent`       | 0.001 | 0.411 | 0.639 | 0.562 |

Physical reading of the Ψ = 37° → 69° rotation:

- **σ_n up, τ_strike down** — SHmax now lies closer to the fault-normal
  direction (misoriented / "weak fault" geometry): median σ_n_total rises
  from 73.7 to 98.2 MPa while the right-lateral median τ_strike collapses
  from 26.5 to 3.8 MPa.
- **Median μ_apparent drops 0.56 → 0.41** with a long tail to ~0.001 on
  near-perpendicular segments. (|τ| stays high in the median only because
  `tau_dip` grows — see next point — the *strike-parallel* drive is what Ψ
  kills.)
- **τ_dip more negative** (−16 → −27 MPa median): the horizontal normal
  stress resolved onto dipping segments now contrasts more strongly with
  `Sv = 50` (fault-normal horizontal stress ≈ 112 MPa at Ψ = 69° vs ≈ 82 MPa
  at Ψ = 37°).
- Extremes: `sigma_n_total max = 120.0 = SHmax` (segments locally
  perpendicular to SHmax) and `|tau| max = 35.0` — both attained because the
  curved SAF sweeps through the critical orientations.

---

## Stage 6 — verification (run and PASSED)

```bash
cd stress
BASE=safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax120_psi69_az23
python code/verify_onfault_stress.py \
    --fault-vtu    results/$BASE/${BASE}_fault_stress.vtu \
    --bulk-vtu     results/$BASE/${BASE}_bulk_stress.vtu \
    --summary-json results/$BASE/${BASE}_summary.json \
    --report-json  results/$BASE/verify_report.json
```

**Result: 12 / 12 fields pass.** Bulk cell-data exact (L∞ = 0); fault
cell-data L∞ ≤ 4.5e-08 Pa; fault point-data L∞ ≤ 6.0e-08 Pa — machine
round-off of a ~1e8 Pa magnitude, far below the 1 Pa (cell) / 1e3 Pa (point)
tolerances. Full numbers in `results/$BASE/verify_report.json`.

---

## Conventions (carried end-to-end, identical to all prior runs)

- Frame: UTM Zone 11 N metres, `(x = east, y = north, z = up)`, `z ≤ 0` underground.
- σ⁰: **compression POSITIVE** (SEAS internal) after the single source-site flip.
- `sigma_n_eff = sigma_n_total − P_p`, `P_p > 0`.
- Tandem fault basis `s = up × n̂`, `d = s × n̂`; `tau_strike > 0` ⇒ right-lateral.
- Units: **MPa** in all VTU/JSON.

---

## References

- `stress/docs/WORKFLOW_fault_stress_z0embed_hz1671_psi37_az351.md` — the Ψ = 37° run this one redoes.
- `stress/README.md` — general (multi-variant, batch) pipeline reference.
- `stress/code/hickman_and_zoback_regional_stress_projection.py` — H&Z analytic core.
- `stress/code/project_to_fault_stress.py` — projection driver (CLI, writers).
- `stress/code/verify_onfault_stress.py` — full-population analytic verifier.
- Hickman, S. & Zoback, M. (2004), *Stress orientations and magnitudes in the SAFOD pilot hole*, GRL — `stress/reference/`.
