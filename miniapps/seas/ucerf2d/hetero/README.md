# UCERF 2.5D heterogeneity inputs (velocity and initial stress at 10 km)

The tool that builds the 2D input file of the heterogeneity runs (plan
`heterogeneity-plan-rev4`, Step A, phase H0), and its checks. The model is the
horizontal slice at z = 10 km of the UCERF 2.5D map-view cycle
(`ucerf2d/mesh/ucerf2d_saf_sjf_500m_tilted_D415km_R70.msh`). Its frame is the
azimuthal equidistant projection about (37.25 N, -119.25 E), R = 6371 km, x
east and y north, in metres.

| file | what |
|---|---|
| `hetero_lib.py` | frame (AEQD forward and inverse, grid convergence), fault points from the mesh, UCVM cvmsi queries, the so2 stress method (`CsmField` with a UTM-11 Delaunay, `magnitudes_regime`, `horizontal_tensor`, `fault_traction`) |
| `fault_tractions.py` | the so2 prestress at the 825 fault points (every 1 km): sigma_n', tau, the regime, CVM Vp, Vs, rho and L_b; writes `fault_tractions.csv/.npz` and `summary.json` (with the k sweep) |
| `check_vs_so2.py` | gate S0: our 10 km stress against the so2 deck's `stress_volume.nc` at z = -10 km (sigma_H, sigma_h, theta and the on-fault tractions) |
| `build_hetero_grid.py` | the input file: a 1 km grid over the mesh box. It holds vp, vs and rho (CVM-S4.26.M01 at 10 km, with a 50 km linear taper to the deck constants outside the model) and s_xx, s_yy, s_xy and s_zz (the so2 stress, effective, Pa, tension positive). NetCDF-4, float64, zlib |
| `qc_hetero_grid.py` | gates Q1 (the file at the fault points against `fault_tractions.py`) and Q2 (the deck constants outside the taper, constant s_zz), with maps |

## The so2 stress at 10 km

These are so2's equations from `~/Downloads/so2/so2_setup/so2_stress/docs/METHOD.md`, Steps 1, 2 and 9 to 11.

- S_v' = 171.116 MPa: -s_zz of the current so2 deck's `stress_volume.nc` (md5 24877d65).
- k = sigma_1 / sigma_3 = 1.45.
- theta and R come from the YHSM-2013 CSM tensor, linear in the UTM-11 Delaunay triangles of its nodes (as so2).
- The regime comes from A_phi, using the A_phi-bar mixture in mixed triangles.
- theta is measured from true north and rotated into the mesh frame by the grid convergence gamma: alpha = theta + gamma (plan D10).

The T_d0 edit of so2 is not applied, because the faults are vertical.

## Run (Mac, conda `pythonenv`; UCVM in `~/Downloads/ucvm/ucvm_install`)

    cd miniapps/seas
    PY=~/miniforge/envs/pythonenv/bin/python
    OUT=~/projects/seas-mfem-fully-dynamic-fault-opening-hetergeneity-results
    $PY ucerf2d/hetero/fault_tractions.py --out $OUT/so2-stress-10km                 # ~10 s
    $PY ucerf2d/hetero/check_vs_so2.py $OUT/so2-stress-10km/fault_tractions.npz \
        ~/Downloads/so2/so2_inputs/alt_c200_ctm100_visco_plast_aphi_cvmsi_phi3040_dc0p2_pph1lam045_fw035h2nw_h4c005_h6k100_h7k100_vw360n_t200_1hz/stress_volume.nc \
        $OUT/so2-stress-10km/check_vs_so2.json
    $PY ucerf2d/hetero/build_hetero_grid.py --out $OUT/grid-z10km/ucerf2d_z10km_cvms426_so2.nc   # ~1 min
    $PY ucerf2d/hetero/qc_hetero_grid.py $OUT/grid-z10km/ucerf2d_z10km_cvms426_so2.nc \
        $OUT/so2-stress-10km/fault_tractions.npz $OUT/grid-z10km/qc

The reader in seas-mfem is `io/grid_field_2d.hpp`, with the unit test
`tests/unit/test_grid_field_2d.cpp` (`make seas_test_grid_field_2d`).

## Using the file in a cycle deck

    [velocity]
    kind     = "grid"
    file     = "ucerf2d/hetero/ucerf2d_z10km_cvms426_so2.nc"
    sampling = "quadrature"      # or "cell_average" (SeisSol's rule); both regimes

    [stress]
    kind = "grid_tensor"
    file = "ucerf2d/hetero/ucerf2d_z10km_cvms426_so2.nc"

`outside = "clamp"` is optional in both sections; the default `"abort"` stops at any point off the grid.

- **Quadrature mode:** the volume integrals of K (with the 2.5D β) and the mass read the grid at every quadrature point. The fault- and interior-face terms take one value per element, the stiffest quadrature point.
- **η:** ½√(μρ) from the grid at each fault point.

Only `seas_fd_cycle_driver` reads these kinds; the other drivers refuse them. `make_test_grid.py` writes the uniform and perturbed files of the gates. The unit test is `tests/cycle/test_cycle_hetero.cpp` (`make -f cycle/cycle.mk test-cycle-hetero`).
