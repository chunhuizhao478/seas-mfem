# verification

| Gate | Program / script | Result |
| --- | --- | --- |
| R1 gold (RACCOON mode I) | `raccoon_mode1_gold.sh` | nodal fields to 1e-10 |
| R1 full curve vs MOOSE M-R1 | `raccoon_mode1_moose.sh` | 1.6e-5, peak to 1e-11 |
| C1 energy (free faces, ζ = 0) | `pff_driver -fd -c spe_hftc/configs/verification/c1_section_free_zeta0.toml -tf 3.82e-4` | 0.048 % of the input work over 10 021 steps (gate 0.1 %) |
| C1 pressurised cavity (spe_hftc PLAN A.8) | `mpirun -np 2 tests/verification/cavity_hoop` | explicit vs static FE 0.35 % (gate 1 %); static FE vs Lamé 3.5 % (gate 5 %, first order in h: 1.74 % at h/2) |
| C1 dashpot reflection | `tests/verification/dashpot_reflection -w p` / `-w s` | 0.24 % amplitude (gate 2 %) |
| C1 v_max(t) vs MOOSE M-C1 | (pending the MOOSE run) | |

Build: `make MFEM_BUILD_DIR=... tests/verification/<name>`. The annulus mesh of
`cavity_hoop` is `meshes/annulus_rb4p15_R16p6_h0p3.msh` (`meshes/make_annulus.py`).
MOOSE comparisons cite the run (machine, job id or date, farms commit).
