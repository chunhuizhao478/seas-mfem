# M-C1: MOOSE run of the spe_hftc 2D section with damage frozen at 0

- Run: 2026-11-21, Andy's Mac (arm64), `mpiexec -n 2`, 3193 s wall (MOOSE
  3185 s, 586 MB). Inputs byte-copied from farms_cdms
  `pulsepower/spe_hftc/code_2d/{elasticity_absorbing,fracture}.i` (farms
  46ad2eb, tree dirty only outside these files), mesh = seas-mfem 1d691e5
  `miniapps/pff/spe_hftc/mesh/spe_section_h0p3.msh` (blob 100ae829).
- Command: `farms-opt -i elasticity_absorbing.i Np=1 eta=0.45 zeta_pulse=1e-8
  mesh_file=spe_section_h0p3.msh Executioner/end_time=2.5e-5 "MultiApps/active="
  "Transfers/active=" Outputs/checkpoint/enable=false` (with the SDKROOT/PATH
  JIT env of the other Mac runs). 269 steps: dt 1e-8 growing x1.25 to 3e-8
  (held to 0.6 us), then to 1e-7.
- `elasticity_absorbing_csv_6sig.csv`: 9 of the 19 CSV columns at 6
  significant figures (the full file stays in ~/pff_moose_refs/M-C1_section_d0/).
  MOOSE signs: external_work < 0; full_input_energy = -external_work -
  damping_work (dashpots only; the Rayleigh loss is not tracked).

## pff against M-C1 (pff `out/c1_section`, cfl 0.4, dt 3.43e-8, lumped mass)

| t (us) | v_max pff/MOOSE | p_wall | W_ext | dashpot work | KE | SE |
| --- | --- | --- | --- | --- | --- | --- |
| 0.29 | 113.9 / 99.5 | 2.149e8 / 2.144e8 | 106.3 / 107.3 | 0 / 0 | 49.6 / 48.5 | 52.2 / 46.5 |
| 1.03 | 90.8 / 80.5 | 1.653e8 / 1.648e8 | 410.4 / 409.5 | 0 / 0 | 191.9 / 191.2 | 200.3 / 199.7 |
| 3.43 | 59.4 / 53.7 | 6.18e7 / 6.16e7 | 686.7 / 671.0 | 0 / 0 | 257.2 / 257.5 | 384.9 / 385.0 |
| 12.0 | 29.9 / 28.4 | 1.81e6 / 1.81e6 | 652.9 / 636.3 | 65.5 / 77.7 | 250.9 / 248.4 | 242.1 / 239.5 |
| 20.0 | 11.5 / 11.1 | 6.80e4 / 6.78e4 | 651.7 / 635.1 | 458.2 / 462.6 | 46.9 / 46.9 | 45.2 / 45.2 |
| 25.0 | 6.45 / 6.43 | 8970 / 8820 | 651.7 / 635.2 | 534.6 / 538.5 | 6.22 / 6.13 | 9.28 / 9.15 |

- Kinetic and strain energy agree to 0.1-1 % from 1 us on; peak wall
  pressure 2.151e8 vs 2.144e8.
- v_max (nodal max |v|) is 13-14 % higher in pff in the first microsecond and
  converges to 0.4 % by 25 us. The sharp front at the wall is where a
  lumped-mass explicit scheme and a consistent-mass implicit Newmark scheme
  (dt 1e-8-3e-8 there) disperse differently; the energies, which are
  integrals, agree.
- External work differs by 16 J/m (2.5 %) from about 2 us on while KE + SE
  agree (642 at 3.43 us in both): pff integrates it with the trapezoidal rule
  and closes its energy budget to 3e-4 J/m (Kelvin-Voigt loss 44.6 J/m by
  3.43 us); MOOSE's FarmsExternalWork is a one-sided rule at dt up to 1e-7
  and its budget leaves a 28.5 J/m gap at 3.43 us.
