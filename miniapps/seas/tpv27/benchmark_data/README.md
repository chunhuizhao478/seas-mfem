# TPV27 benchmark data

Reference fault time series for SCEC TPV27 (TPV26 with Drucker-Prager
off-fault plasticity), one subfolder per code, as in
`tpv205/benchmark_data/`. Problem description:
`../../tpv26/benchmark_document/TPV26_27_Description_v13.pdf`. Data supplied
by Chunhui Zhao, 2026-10-20.

## `scec_drdg3d/`

DRDG3D (Wenqiang Zhang, 2023/02/07), 200 m on the fault, order 4.
12 on-fault stations, one file each, the same stations and format as
`../../tpv26/benchmark_data/scec_drdg3d/`:

    tpv27_drdg3d_x2_<strike_km>_x3_<depth_km>.txt
    (x2, x3) = (-15,10) (-5,0) (-5,5) (-5,10) (-5,15) (0,10) (5,10)
               (10,10) (15,0) (15,5) (15,10) (15,15) km
    x2 = along strike, x3 = down-dip depth (x3 = 0 is the surface).

11,706 samples at dt = 0.00102518 s, t = 0 to 12.0 s (the description asks
for 13 s). Columns (m, m/s, MPa):

    t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress
    h = horizontal = along strike; v = vertical = along dip.

- `n-stress` is the effective normal stress, compression negative; flip the
  sign when overlaying (this repo uses compression positive).
- The two surface stations (x2 = -5 and 15 km, x3 = 0) never slip. This is
  read as physical, not a data error: near the surface the rock yields at
  about c cos(phi) = 1.34 MPa, which keeps the fault shear stress there
  (peak 1.74 and 1.67 MPa) below the 4 MPa frictional cohesion. Only one
  code is available to confirm it.
- Not included: the off-fault surface stations and the rupture-time
  contour file (`cplot`).
