# WORKFLOW: minimum-mu regional stress, Sv = 55 (triq, SHmax=100, Shmin=60, Ψ=69°)

**Target artefact (CURRENT preferred minimum-mu stress field):**

```
stress/results/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax100_shmin60_sv55_psi69_az23/
    ..._fault_stress.vtu  ..._bulk_stress.vtu  ..._summary.json  verify_report.json
```

Written 2026-06-09. Supersedes the Sv = 50 run
(`WORKFLOW_fault_stress_z0embed_triq_shmax100_shmin60_psi69_az23.md`, which
holds the full decision chain: Ψ-sweep degeneracy, H&Z error-box grid
search, corner-pinned optimum). This iteration raises **Sv 50 → 55 MPa**
after the Sv-sensitivity analysis showed Sv is the controlling clamp on
the dipping critical planes:

```
mu_max bound = (s1' - s3') / (2 sqrt(s1' s3')),  s1' = SHmax - P_p = 80
               s3' = min(Shmin, Sv) - P_p
Sv = 45 -> 0.615   (grid_shmax_shmin_triq_psi69_pp20_sv45/)
Sv = 50 -> 0.510   (grid_shmax_shmin_triq_psi69_pp20/)
Sv = 55 -> 0.425   (grid_shmax_shmin_triq_psi69_pp20_sv55/)  <- this run
```

Sv = 55 corresponds to the lithostat at ~2.1 km (26 MPa/km) — the depth of
the deepest Hickman & Zoback (2004) Fig. 4a measurements — where
hydrostatic P_p ~ 20.8 ~ 20 MPa: the anchor is self-consistent.
Shmin has NO effect on mu_max once Shmin >= Sv (plateau verified in
`grid_shmin_extend_triq_psi69_pp20/`); it stays at its error-bar ceiling.

## Parameters

| Parameter | Value | provenance |
|---|---|---|
| SHmax | 100 MPa | H&Z 2 km error-bar floor (grid argmin) |
| Shmin | 60 MPa | H&Z 2 km error-bar ceiling (grid argmin) |
| Sv | **55 MPa** | lithostat at the 2.1 km anchor |
| P_p | 20 MPa | hydrostatic at the anchor |
| SHmax azimuth | 23° cw-from-N (Ψ = 69° vs strike hint 314°) | Townend & Zoback 2004 southern-SAF mean |

Effective principals: s1' = 80, s2' = 40, s3' = 35 MPa (Sv least).

## Reproduce

```bash
conda activate pythonenv
cd stress/code
python grid_shmax_shmin_mu.py \
    --fault-vtu ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
    --out-dir ../results/grid_shmax_shmin_triq_psi69_pp20_sv55 \
    --SHmax-range 100 130 2.5 --Shmin-range 40 60 2 --psi 69 --Sv 55
python project_to_fault_stress.py \
    ../../experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_fault.vtu \
    safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq_shmax100_shmin60_sv55_psi69_az23 \
    --write-bulk --SHmax 100 --Shmin 60 --Sv 55 --P_p 20 --SHmax-az 23
```

## Bulk Cauchy tensor σ⁰ at z = 0 — derivation

Convert the principal stress direction/magnitude [SHmax, Shmin] to the
[East, North, Up] frame:

[1] map azimuth angle to math angle:

```
alpha = (90 - azimuth) mod 360 = (90 - 23) mod 360 = 67 deg
```

(SHmax points N23E geographically = 67 deg ccw from East in math convention.)

[2] principal stress tensor (compression negative, H&Z source convention):

```
Sigma'_HZ = diag(-SHmax, -Shmin, -Sv) = diag(-100, -60, -55)  MPa
```

[3] rotation matrix: c = cos(67) = 0.39073, s = sin(67) = 0.92050

```
        [ c  -s   0 ]   [ 0.39073  -0.92050   0 ]
Rz(a) = [ s   c   0 ] = [ 0.92050   0.39073   0 ]
        [ 0   0   1 ]   [ 0         0         1 ]
```

[4] rotate into [East, North, Up] and flip sign (compression positive):

```
sigma0 = -sigma0_HZ = Rz * diag(100, 60, 55) * Rz^T

sigma_EE = SHmax*c^2 + Shmin*s^2 = 100*0.15267 + 60*0.84733 =  66.107
sigma_NN = SHmax*s^2 + Shmin*c^2 = 100*0.84733 + 60*0.15267 =  93.893
sigma_EN = (SHmax - Shmin)*c*s  = 40*0.35967                =  14.387
```

[5] stress components in the [East, North, Up] frame:

```
         [  66.107   14.387    0.0 ]
sigma0 = [  14.387   93.893    0.0 ]   MPa  (compression POSITIVE)
         [   0.0      0.0     55.0 ]
```

Matches the `sigma_global_at_z0_MPa` echo in `_summary.json` exactly.
(vs the shmax120 run: sigma_EN 23.379 -> 14.387 because the horizontal
differential stress fell 65 -> 40 MPa; sigma_NN 110.08 -> 93.89 from the
SHmax reduction; Up entry is the new Sv = 55.)

## Results (60658 cells, verifier 12/12 at machine round-off)

| Field | min | median | max | Sv=50 run | shmax120 run |
|---|---:|---:|---:|---:|---:|
| sigma_n_total_MPa | 57.50 | 86.04 | 100.00 | 84.76 max 100 | 98.18 max 120 |
| sigma_n_eff_MPa   | 37.50 | 66.04 | 80.00 | 64.76 / 80.0 | 78.18 / 100.0 |
| tau_strike_MPa    | -19.67 | 2.34 | 17.55 | identical | -31.96/3.81/28.52 |
| tau_dip_MPa       | -22.50 | -17.64 | 7.57 | -25.0/-19.8/8.5 | -35.0/-27.3/11.7 |
| tau_magnitude_MPa | 0.05 | 20.24 | 22.50 | 25.0 max | 35.0 max |
| mu_apparent       | 0.001 | 0.311 | **0.425** | 0.348 / 0.510 | 0.411 / 0.639 |

(tau_strike is Sv-independent — horizontal-plane shear only involves the
horizontal principals; the dip shear and mu shrink with the higher clamp.)

Margin below a mu_s = 0.6 friction law: 0.175 everywhere.
Grid at Sv=55: zero fault area above mu 0.55 for SHmax <= ~115; at the
optimum, A% > 0.50 = 0 (max 0.425 < 0.5) — no critically-stressed patches
at any conventional threshold.

## Caveats

- Same corner-pinning caveat as the Sv = 50 run: (100, 60) sits on the
  edge of the H&Z error box.
- Driving stress is further reduced (|tau| max 22.5 MPa, median 20.2);
  check nucleation/propagation viability downstream.
- Shmin = 60 > Sv = 55: mildly transpressional ordering at the anchor,
  consistent with shallow Mojave observations (Hickman et al. 1988).
