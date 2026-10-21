# TPV26 benchmark data

Reference fault time series for SCEC TPV26 (elastic bulk; TPV27 is the same
problem with Drucker-Prager off-fault plasticity), one subfolder per code, as
in `tpv205/benchmark_data/`. Problem description:
`../benchmark_document/TPV26_27_Description_v13.pdf`. Data supplied by
Chunhui Zhao, 2026-10-20.

## `scec_drdg3d/`

DRDG3D (Wenqiang Zhang, 2023/02/07), 200 m on the fault, order 4.
12 on-fault stations, one file each:

    tpv26_drdg3d_x2_<strike_km>_x3_<depth_km>.txt
    (x2, x3) = (-15,10) (-5,0) (-5,5) (-5,10) (-5,15) (0,10) (5,10)
               (10,10) (15,0) (15,5) (15,10) (15,15) km
    x2 = along strike, x3 = down-dip depth (x3 = 0 is the surface).

11,706 samples at dt = 0.00102518 s, t = 0 to 12.0 s (the description asks
for 13 s). Columns (m, m/s, MPa):

    t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress
    h = horizontal = along strike; v = vertical = along dip.

- `n-stress` is the effective normal stress, compression negative
  (-175.64 MPa at 10 km depth at t = 0). This repo uses compression
  positive; flip the sign when overlaying.
- At t = 0 the stresses match the description's initial-stress formulas
  at 0, 5 and 10 km depth (e.g. 13.83 MPa shear and -87.82 MPa effective
  normal stress at 5 km). At 15 km, where the depth taper begins, they are
  0.4 % (shear, 41.33 vs 41.50 MPa) and 0.03 % (normal) below the formula;
  read as DRDG3D's element interpolation smoothing the kink of the taper
  (inferred, not confirmed with the author).
- Not included: the off-fault surface stations and the rupture-time
  contour file (`cplot`).
