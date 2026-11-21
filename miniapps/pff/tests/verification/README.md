# verification

| Gate | Program / script | Result |
| --- | --- | --- |
| R1 gold (RACCOON mode I) | `raccoon_mode1_gold.sh` | nodal fields to 1e-10 |
| R1 full curve vs MOOSE M-R1 | `raccoon_mode1_moose.sh` | 4.4e-12 (gate 1e-9), same pass count every step |
| C1 energy (free faces, ζ = 0) | `pff_driver -fd -c spe_hftc/configs/verification/c1_section_free_zeta0.toml -tf 3.82e-4` | 0.048 % of the input work over 10 021 steps (gate 0.1 %) |
| C1 pressurised cavity (spe_hftc PLAN A.8) | `mpirun -np 2 tests/verification/cavity_hoop` | explicit vs static FE 0.35 % (gate 1 %); static FE vs Lamé 3.5 % (gate 5 %, first order in h: 1.74 % at h/2) |
| C1 dashpot reflection | `tests/verification/dashpot_reflection -w p` / `-w s` | 0.24 % amplitude (gate 2 %) |
| C1 v_max(t) vs MOOSE M-C1 | (pending the MOOSE run) | |
| C2 one pulse, 10 µs, η = 1 | `pff_driver -c spe_hftc/configs/verification/c2_section_eta1.toml -tf 1e-5` | d_max 0.949, d ≥ 0.5 to 1.88 mm, d ≥ 0.9 to 0.61 mm, crack 76.6 J/m, energy error 0.13 J/m; vs M-C2a pending |
| C2 one pulse, 10 µs, η = 0.3 | `... c2_section_eta0p3.toml -tf 1e-5` | d_max 0.166, crack 5.34 J/m; vs M-C2b pending |
| C2 at the 2026-10-21 MOOSE setting (free, ζ = 0) | `... c2_section_eta{1,0p3}_free_zeta0.toml -tf 1e-5` | η = 1: d_max 1, 3.24 mm, 97.5 J/m (MOOSE then: 0.97, 1.76 mm, 77 J/m); η = 0.3: d_max 0.25, 6.8 J/m (MOOSE then: 0.18, 5.6 J/m). Different mesh build and old input; not a gate |

Build: `make MFEM_BUILD_DIR=... tests/verification/<name>`. The annulus mesh of
`cavity_hoop` is `meshes/annulus_rb4p15_R16p6_h0p3.msh` (`meshes/make_annulus.py`).
MOOSE comparisons cite the run (machine, job id or date, farms commit).
