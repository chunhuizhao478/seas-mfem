# M-C2c: spe_hftc 2D section, eta = 1, zeta_pulse = 0, one pulse to 10 us (MOOSE)

- farms@cdms 46ad2eb, pulsepower/spe_hftc/code_2d/elasticity_absorbing.i + fracture.i (byte-identical), mesh spe_section_h0p3.msh.
- Andy's Mac, 2 ranks, 2026-11-21 22:16:45 to 23:31:47 PST (4502.6 s), rc 0; 119 steps, all CONVERGED_RELATIVE.
- `mpiexec -n 2 farms-opt -i elasticity_absorbing.i Np=1 eta=1.0 zeta_pulse=0 mesh_file=spe_section_h0p3.msh Executioner/end_time=1e-5 Outputs/checkpoint/enable=false`
- Same as M-C2a except zeta_pulse = 0. Mac folder ~/pff_moose_refs/M-C2c_section_eta1_zetapulse0/.

At 10 us (MOOSE vs pff c2_section_eta1_zetapulse0, cloud np4):

| | MOOSE | pff |
|---|---|---|
| crack energy (J/m) | 82.48 | 97.54 |
| d_max | 0.998 | 1.0 |
| d >= 0.5 extent (mm) | 2.446 | 3.238 |
| d >= 0.9 extent (mm) | 1.737 | 2.490 |
| KE + SE (J/m) | 1103.7 | 1093.7 |
| external work (J/m) | 1161.6 | 1190.3 |

Crack energy agrees within 2% to 6 us, then grows faster in pff. With zeta_pulse = 1e-8 (M-C2a) the codes agree to 0.01% at 10 us. Open: time step (MOOSE Newmark dt = 1e-7 after 1 us vs pff explicit 3.8e-8).

pff time-step check (cloud np4, 2026-11-22): cfl 0.4 -> 0.2 changes the 10 us crack energy 97.54 -> 96.07 J/m; d >= 0.5 extent 3.238 mm in both, d >= 0.9 extent 2.490 mm in both. So pff is converged in dt; the MOOSE dt = 1e-7 step is the next suspect (M-C2d, dt <= 3e-8).
