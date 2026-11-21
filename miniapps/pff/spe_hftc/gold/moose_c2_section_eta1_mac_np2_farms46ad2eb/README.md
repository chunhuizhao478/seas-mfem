# M-C2a: MOOSE run of the spe_hftc 2D section, eta = 1, damage on, to 10 us

- Run: 2026-11-21, Andy's Mac (arm64), `mpiexec -n 2`, 3061 s wall. Inputs
  byte-copied from farms_cdms `pulsepower/spe_hftc/code_2d/{elasticity_absorbing,fracture}.i`
  (farms 46ad2eb), mesh = seas-mfem 1d691e5 `spe_hftc/mesh/spe_section_h0p3.msh`.
- Command: `farms-opt -i elasticity_absorbing.i Np=1 eta=1.0 zeta_pulse=1e-8
  mesh_file=spe_section_h0p3.msh Executioner/end_time=1e-5
  Outputs/checkpoint/enable=false` (SDKROOT/PATH JIT env as the other Mac runs).
- 119 steps (dt as M-C1), every step CONVERGED_RELATIVE; fixed-point passes
  1-30:1 31-90:2 91-102:3 103-119:2.
- `elasticity_absorbing_csv_6sig.csv`: all 19 CSV columns, 6 significant figures.
- Exodus (last frame, 10 us): max d 0.9510, d >= 0.5 to 1.654 mm and d >= 0.9
  to 0.609 mm from the wall (= the fracture_extent postprocessors).

## pff C2 (`spe_hftc/configs/verification/c2_section_eta1.toml`, single damage pass per explicit step)

| t (us) | crack energy pff/MOOSE (J/m) | d >= 0.5 extent (mm) | d >= 0.9 (mm) | KE | SE | v_max |
| --- | --- | --- | --- | --- | --- | --- |
| 2.13 | 10.64 / 10.67 | 0 / 0 | 0 / 0 | 463.4 / 463.6 | 553.3 / 553.1 | 95.7 / 86.0 |
| 3.43 | 38.40 / 38.35 | 0.519 / 0.519 | 0 / 0 | 465.6 / 466.2 | 633.3 / 633.3 | 78.6 / 71.1 |
| 5.03 | 62.81 / 63.02 | 1.107 / 1.107 | 0 / 0 | 430.6 / 430.8 | 633.9 / 634.0 | 65.4 / 59.8 |
| 6.93 | 74.89 / 74.84 | 1.415 / 1.447 | 0 / 0 | 407.2 / 407.9 | 604.9 / 605.5 | 54.8 / 51.1 |
| 8.03 | 75.97 / 75.93 | 1.654 / 1.654 | 0 / 0 | 407.3 / 408.0 | 585.1 / 585.7 | 50.4 / 47.3 |
| 10.0 | 76.57 / 76.58 | 1.877 / 1.654 | 0.609 / 0.609 | 424.0 / 424.8 | 542.3 / 542.5 | 44.2 / 41.8 |

- Crack energy within 0.3 % at every frame (0.01 % at 10 us); KE and SE
  within 0.2 %; the d >= 0.5 front sits on the same node as MOOSE through
  9.5 us. At 9.51 us one more node (r - r_b = 1.877 mm, the tip of an
  incipient radial crack) crosses d = 0.5 in pff; MOOSE's d there stays below
  0.5 at 10 us. d_max 0.949 vs 0.951.
- v_max is 6-11 % higher in pff throughout, as in M-C1 (lumped mass at the
  front).
