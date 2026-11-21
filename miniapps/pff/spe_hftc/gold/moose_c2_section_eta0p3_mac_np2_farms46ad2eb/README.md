# M-C2b: MOOSE run of the spe_hftc 2D section, eta = 0.3, damage on, to 10 us

- Run: 2026-11-21, Andy's Mac (arm64), `mpiexec -n 2`, 2001 s wall. Same
  inputs (farms 46ad2eb), mesh (seas-mfem 1d691e5) and command as M-C2a with
  `eta=0.3`. 119 steps, all CONVERGED_RELATIVE; fixed-point passes
  1-40:1 41-59:2 60-62:3 63-70:2 71-119:1.
- `elasticity_absorbing_csv_6sig.csv`: all 19 CSV columns, 6 significant figures.
- Exodus at 10 us: max d 0.1657 at the wall; no node with d >= 0.5. Lower
  thresholds (distance from the wall): d >= 0.15 to 0.183 mm, d >= 0.1 to
  0.520 mm, d >= 0.05 to 0.803 mm, d >= 0.01 to 1.107 mm.

## pff C2 (`spe_hftc/configs/verification/c2_section_eta0p3.toml`)

| t (us) | crack energy pff / MOOSE (J/m) | pff d_max | KE | SE |
| --- | --- | --- | --- | --- |
| 2.43 | 0.0372 / 0.0344 | 0.0077 | 200.7 / 200.8 | 258.9 / 258.4 |
| 3.43 | 2.083 / 2.077 | 0.089 | 193.8 / 194.0 | 287.7 / 287.7 |
| 4.43 | 4.898 / 4.906 | 0.158 | 182.3 / 182.6 | 294.3 / 294.4 |
| 5.13 | 5.335 / 5.340 | 0.1656 | 177.2 / 177.4 | 292.3 / 292.3 |
| 10.0 | 5.336 / 5.340 | 0.1656 (MOOSE 0.1657) | 202.6 / 202.9 | 223.6 / 223.7 |

Crack energy within 0.3 % from 3 us on (0.08 % at 10 us), d_max within
0.1 %, KE and SE within 0.2 %. The first damage (2.2 us) differs by 8 % in
absolute terms of 0.003 J/m.
