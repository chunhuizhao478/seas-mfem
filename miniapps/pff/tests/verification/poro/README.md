# Poroelastic verification (phase C4)

Driver-level checks of the hydromech path of `pff_driver` (explicit u, backward
Euler p, single pass) against closed forms, with d frozen (`-fd`). No MOOSE
input: these are the standard consolidation problems.

| case | what it checks | result (2026-11-26, cloud) |
|---|---|---|
| `undrained_column.toml` | k = 0: undrained modulus H_u = K + 4G/3 + alpha^2 M and wave speed, from the first rebound of a step-loaded column | max work -0.9 %, rebound time -0.8 % (40 elements along the column) |
| `terzaghi_column.toml` | 1D consolidation, one-way drainage: pressure at the sealed end and settlement against Terzaghi's series | pressure within 0.32 %, settlement within 0.17 % to T_v = 0.5 |
| `deleeuw_cylinder.toml` | de Leeuw's plane-strain cylinder (Mandel-Cryer effect): centre pressure against the residue series of the Laplace solution | peak 0.11 %, all times within 0.45 % to T = 0.3 |

    python3 make_column_mesh.py --L 0.01 --W 0.0005 --ny 40 --nx 2 -o column_L10_h0p25.msh
    python3 make_quarter_disc_mesh.py --a 0.01 --nr 24 -o quarter_disc_a10_nr24.msh
    ./pff_driver -fd -c tests/verification/poro/undrained_column.toml -o OUT   # 1 rank, 3 s
    ./pff_driver -fd -c tests/verification/poro/terzaghi_column.toml  -o OUT   # 1 rank, 4 min
    mpirun -np 2 ./pff_driver -fd -c tests/verification/poro/deleeuw_cylinder.toml -o OUT   # 5 min
    python3 tests/verification/poro/check_poro_verification.py {undrained|terzaghi|deleeuw} OUT/<name>_c1.csv

The step load is applied from rest; the consolidation cases remove its waves
with mass-proportional damping (alpha_mass near critical for the first mode),
so the comparison starts after a few wave transits (T_v > 0.01). The cylinder
replaces Mandel's slab of PLAN.md (Andy, 2026-11-26): the slab needs a rigid
plate boundary that no production case uses.
