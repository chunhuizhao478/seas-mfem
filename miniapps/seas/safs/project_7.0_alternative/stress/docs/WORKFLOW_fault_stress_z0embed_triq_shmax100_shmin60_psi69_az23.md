# WORKFLOW: minimum-mu regional stress (triq mesh, SHmax=100, Shmin=60, Ψ=69°)

> **SUPERSEDED** by the Sv = 55 iteration (mu_max 0.510 → 0.425):
> `WORKFLOW_fault_stress_z0embed_triq_shmax100_shmin60_sv55_psi69_az23.md`.
> This doc remains the record of the decision chain (Ψ-sweep degeneracy,
> error-box grid search, corner-pinned optimum).

**Target artefact:**

```
stress/results/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax100_shmin60_psi69_az23/
    ..._fault_stress.vtu  ..._bulk_stress.vtu  ..._summary.json  verify_report.json
```

Written 2026-06-09. Final outcome of the mu-minimization study that started
from the `..._triq_shmax120_psi69_az23` run (see its workflow doc for the
pipeline mechanics — identical here; only the magnitudes changed).

## Decision chain (analysis artefacts in stress/results/)

1. **Ψ sweep at SHmax=120/Shmin=55** (`psi_sweep_triq_shmax120_pp20/`):
   field-max mu_apparent is PINNED at the orientation bound
   `mu* = (s1'-s3')/(2 sqrt(s1' s3')) = 0.639` for every Ψ in [52, 86] —
   the curved, moderately dipping SAF surface (median dip 59.8°, 22% of
   cells dip < 45°) always contains a critically oriented patch. Ψ only
   moves the critically-stressed AREA (min near Ψ ≈ 71-73°).
2. **SHmax variants** (`psi_sweep_triq_shmax113p6_pp20/`, `..._105p8_...`):
   confirmed max mu == bound at all Ψ; the bound is controlled by the
   magnitudes, not the orientation.
3. **H&Z (2004) Fig. 4a error bars at the ~2 km anchor** (user-supplied
   bounds): SHmax ∈ [100, 130], Shmin ∈ [40, 60] MPa.
4. **(SHmax, Shmin) grid search** at fixed Ψ=69, Sv=50, P_p=20
   (`grid_shmax_shmin_triq_psi69_pp20/`, script
   `code/grid_shmax_shmin_mu.py`, 13 x 11 = 143 pairs):
   **(SHmax=100, Shmin=60) minimizes EVERY metric simultaneously**
   (field-max, p99, p90, median, area fractions). Mechanism:
   - mu_max is the orientation bound with `s1' = SHmax - P_p`,
     `s3' = min(Shmin, Sv) - P_p`. Lower SHmax lowers s1'.
   - For Shmin < Sv = 50, Shmin becomes the least principal stress and the
     bound RISES steeply (at SHmax=100: Shmin=40 → bound 0.75).
   - For Shmin ≥ 50, s3' = Sv' = 30 is fixed; larger Shmin only shrinks
     the SHmax-Shmin contrast → weaker vertical-plane shear → lower
     median/tails.
   - Hence the optimum sits at the (low-SHmax, high-Shmin) corner of the
     error box. mu decreases monotonically toward that corner, so the
     result is corner-pinned, not an interior optimum.

## Final parameter set (constant with depth, anchor ~2 km)

| Parameter | Value | vs shmax120 run |
|---|---|---|
| SHmax | **100 MPa** (H&Z 2 km lower bound) | 120 |
| Shmin | **60 MPa** (H&Z 2 km upper bound) | 55 |
| Sv | 50 MPa | same |
| P_p | 20 MPa | same |
| SHmax azimuth | 23° (Ψ = 69° vs 314° strike hint) | same |

Effective principal stresses: s1' = 80, s2' = 40, s3' = 30 MPa.
Crustal frictional state mu* = 50 / (2 sqrt(80*30)) = 0.510 — comfortably
sub-Byerlee (margin 0.09 below a mu_s = 0.6 friction law everywhere).

## Reproduce

```bash
conda activate pythonenv
cd stress/code
python project_to_fault_stress.py \
    ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax100_shmin60_psi69_az23 \
    --write-bulk --SHmax 100 --Shmin 60 --Sv 50 --P_p 20 --SHmax-az 23
# grid search behind the choice:
python grid_shmax_shmin_mu.py \
    --fault-vtu ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
    --out-dir ../results/grid_shmax_shmin_triq_psi69_pp20 \
    --SHmax-range 100 130 2.5 --Shmin-range 40 60 2 --psi 69
```

## Results (`_summary.json["stats"]`, 60658 cells, verified 12/12)

| Field | min | median | max | shmax120 run |
|---|---:|---:|---:|---:|
| sigma_n_total_MPa | 52.78 | 84.76 | 100.00 | 53.89 / 98.18 / 120.00 |
| sigma_n_eff_MPa   | 32.78 | 64.76 | 80.00 | 33.89 / 78.18 / 100.00 |
| tau_strike_MPa    | -19.67 | 2.34 | 17.55 | -31.96 / 3.81 / 28.52 |
| tau_dip_MPa       | -25.00 | -19.82 | 8.48 | -35.00 / -27.29 / 11.73 |
| tau_magnitude_MPa | 0.06 | 22.14 | 25.00 | 0.08 / 31.68 / 35.00 |
| mu_apparent       | 0.001 | 0.348 | **0.510** | 0.001 / 0.411 / 0.639 |

Verifier: 12/12 fields at machine round-off (same writer/verifier path as
all prior runs).

## Caveats

- The optimum is **corner-pinned**: it sits exactly at the edge of the
  H&Z error box, so it inherits the full read-off uncertainty of Fig. 4a.
- Minimizing mu also reduces the driving shear stress everywhere
  (|tau| max drops 35 → 25 MPa, median 31.7 → 22.1): a dynamic-rupture
  setup must still verify that nucleation/propagation is possible at the
  reduced stress level.
- SHmax = 100 with Shmin = 60 halves the horizontal differential stress
  (40 vs 65 MPa); trace-parallel vertical segments keep their weak-fault
  character (mu_app ~ 0.18 at Ψ = 69°), consistent with
  Townend & Zoback (2004) Set B (see
  `REFERENCES_saf_principal_stress_for_simulation_2026-06-09.md`).
- Shmin = 60 > Sv = 50 puts the shallow regime mildly transpressional
  (SHmax > Shmin > Sv), consistent with the transitional shallow Mojave
  observations (Hickman, Zoback & Healy 1988), but a constant-with-depth
  tensor still cannot rotate to SHmax > Sv > Shmin at seismogenic depth.
