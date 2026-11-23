# M-C2f: MOOSE cmame 2D solid benchmark, first pulse (reference for pff C2/C3)

- Inputs: farms_cdms 46ad2eb `cmame_revision/2d_solid/benchmark/{elasticity,fracture}.i`, unmodified.
  Mesh `2d_mesh/2d_mesh.msh` sha256 4681910f1c82b65f789663f70cc2a27bb7a91a60f95379ee9365295e639b01b8 (14,418 nodes).
- Run: Andy's Mac, conda env `moose`, `mpiexec -n 2 farms-opt -i elasticity.i Executioner/end_time=1e-5 Outputs/checkpoint/enable=false`,
  2026-11-22 21:37-22:57 PST, 4794 s wall. Fixed dt 5e-8, 200 steps, all fixed-point iterations converged (max 7 of 10).
- Files: `energies_1us.csv` (MOOSE CSV rows every 1 us; full 6-sig CSV is on the Mac at
  ~/pff_moose_refs/M-C2f_cmame2d_solid_benchmark_pulse1/2d_solid/benchmark/csv_6g.txt),
  `damage_frames.csv` (from exodus; r measured from the origin, borehole wall at r = 2 mm; d = 0.9 at 1 us is the seeded box).

## pff comparison (seas-mfem cdms, configs/2d_solid_benchmark.toml, consistent mass, 1 rank, 1707 steps, 217 s)

pff crack energy is minus the first-step (seeded) value, matching MOOSE dissipated_energy_total; pff extents + r_hole.

| t | quantity | MOOSE | pff |
|---|---|---|---|
| 2 us | KE / SE | 1.32864 / 5.40347 | 1.32765 / 5.40698 |
| 2 us | crack | 0.259503 | 0.258388 |
| 5 us | KE / SE | 0.107226 / 3.86319 | 0.107046 / 3.86388 |
| 5 us | crack / damping | 1.85315 / 2.33453 | 1.85079 / 2.32820 |
| 5 us | r(d>=0.5) / r(d>=0.9) mm | 5.9430 / 5.6273 | 5.9430 / 5.6273 |
| 10 us | KE / SE | 0.163953 / 0.20739 | 0.159335 / 0.210448 |
| 10 us | crack / damping | 1.88331 / 3.27761 | 1.88197 / 3.26956 |
| 10 us | r(d>=0.5) / r(d>=0.9) mm | 5.9430 / 5.7480 | 5.9430 / 5.7480 |

Energies in J/m. pff budget error at 10 us: -7e-6 J/m. Crack growth stops near 5.5 us in both.
