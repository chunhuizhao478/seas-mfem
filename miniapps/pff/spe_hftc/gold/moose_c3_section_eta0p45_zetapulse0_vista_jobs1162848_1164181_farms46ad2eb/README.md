# MOOSE 2D section, 20 pulses, eta 0.45, zeta_pulse = 0 (Vista)

- Inputs: farms_cdms 46ad2eb spe_hftc elasticity_absorbing.i / fracture.i, zeta_pulse=0 by CLI (zeta_quench 1e-5 unchanged), farms mesh (gmsh 4.1).
- Vista jobs 1162848 (TIMEOUT at t=0.8656 ms, 12:00:23) + 1164181 (`--recover` from out_cp/9200, COMPLETED, 5:10:13). Run dir
  /scratch/10024/zhaochun/spe2d_vista/section_homog_eta0p45_absorb_zetapulse0/cases_vista/att_homog_eta0p45_absorb_zp0/.
  The solver log of 1162848 was overwritten by the continuation; out.csv/out.e/out_cp are complete.
- `pulse_ends.csv`: out.csv rows (deduplicated by time, last kept) nearest t = k*60 us, k = 1..20.
- **Restart artefact:** dissipated_energy_total = dissipated_energy_dynamic - dissipated_energy_first_step, and farms
  FirstStepElementIntegralVariablePostprocessor keeps its first-step value in non-restartable members, so on --recover it
  re-captures at step 9201. All rows after t = 0.8656 ms are short by a constant 436.661 J/m (checked on 4 overlap rows,
  difference 436.66103-436.66104). `crack_energy_corrected` adds it back. Other columns are continuous across the restart.
- MOOSE dt 1e-7 with fixed-point cutbacks at every t_free switch (zeta 0 -> 1e-5): 184 + 67 non-converged reports, all retried.

## pff comparison (consistent mass, mass quench alpha 2e5/s, section_homog_eta0p45_absorb.toml, 1 rank, 5169 s)

| pulse | ext d>=0.5 mm MOOSE / pff | crack J/m MOOSE / pff | input J/m MOOSE / pff |
|---|---|---|---|
| 2 | 1.43 / 1.36 | 32.1 / 30.7 | 1306 / 1335 |
| 5 | 12.36 / 10.80 | 144.5 / 135.6 | 3320 / 3381 |
| 10 | 18.68 / 15.42 | 267.1 / 240.8 | 6700 / 6937 |
| 13 | 21.42 / 19.75 | 355.8 / 347.9 | 8792 / 9320 |
| 15 | 21.42 / 22.82 | 437.5 / 460.1 | 10274 / 11062 |
| 18 | 32.28 / 32.28 | 584.0 / 694.8 | 12879 / 14160 |
| 20 | 32.28 / 32.28 | 703.8 / 887.9 | 14841 / 16483 |

32.28 mm = mesh corner (max node radius - r_hole). pff budget error at 1.2 ms: -38.8 J/m.

## pff time-step convergence (cfl 0.4 vs 0.2, consistent mass, np4 cloud)

Same 20-pulse section, half the time step (dt 2.33e-8 -> 1.17e-8, 51420 -> 102820 steps, 5899 s -> 11297 s).
Damage table at cfl 0.2: pff_damage_table_cfl0p2.csv (last frame at 1140 us).

| Pulse | crack cfl0.4 | crack cfl0.2 | extent 0.4/0.2 (mm) | W_ext 0.4/0.2 | energy err 0.4/0.2 | cracks @8/@15 mm (both) |
|---|---|---|---|---|---|---|
| 5  | 135.6 | 134.8 | 10.80 / 10.80 | 3381 / 3383 | -1.26 / -0.32 | 4 / 0 |
| 10 | 240.8 | 240.4 | 15.42 / 15.42 | 6937 / 6924 | -7.99 / -1.91 | 4 / 4 |
| 15 | 460.1 | 464.8 | 22.82 / 22.91 | 11062 / 11085 | -20.2 / -5.25 | 4 / 5 |
| 19 | 786.7 | 799.4 | 32.28 / 32.28 | 15298 / 15360 | -34.4 / -9.39 | 3 / 4 |
| 20 | 887.9 | 895.6 | 32.28 / 32.28 | 16483 / 16551 | -38.8 / -10.5 | - |

Crack counts identical at every pulse; energy error drops 3.7x for 2x dt (about second order). MOOSE corrected 703.8 J/m at pulse 20 is outside the pff dt spread (0.9%).
