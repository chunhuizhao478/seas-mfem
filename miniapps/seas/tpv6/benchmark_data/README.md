# TPV6 / TPV7 SCEC reference data

## Authoritative reference: `scec_drdg3d/`
DRDG3D (Wenqiang Zhang, 2023; 200 m on fault, O3) — the SCEC TPV6 reference traces
this Part-C setup is validated against.  10 files = 5 ON-FAULT stations x 2 sides:

    tpv6_drdg3d_{near,far}side_x2_<strike_km>_x3_<depth_km>.txt
    stations: (x2 strike, x3 depth) = (0,0) (-12,0) (12,0) (-12,7.5) (12,7.5) km
    (x2 = along-strike, x3 = down-dip depth; "on fault", surface-breaking at x3=0.)

These are **on-fault, split-node** stations: per the SCEC TPV6/7 spec (Part I) we
report, on EACH side of the fault, the absolute **displacement** and **velocity**
(NOT slip / slip-rate — those are for the contour plots).  Columns (MKS):

    t  h-disp  h-vel  h-stress  v-disp  v-vel  v-stress  n-disp  n-vel  n-stress
    h = horizontal = ALONG-STRIKE ; v = vertical = ALONG-DIP ; n = fault-NORMAL.

### Sign / convention gotchas (carefully verified — match these in any comparison)
- **near vs far material (CONFIRMED, spec p.3 + this data):**
  - **nearside = STRONG / FAST** rock: vp2=6000, vs2=3464, rho2=2670 (high impedance).
  - **farside  = WEAK  / SLOW** rock: vp1=3750, vs1=2165, rho1=2225 (TPV6); for TPV7
    the far side is vp1=5000, vs1=2887, rho1=2670 (lower contrast).
  - Empirical check: the WEAK (far) side moves ~3-4x faster than the strong (near)
    side (peak |h-vel| far/near = 16.0/3.9 at (0,0), 18.1/7.4 at (12,0), ...).
- **normal stress sign:** drdg3d reports `n-stress = -120 MPa` at t=0, i.e.
  COMPRESSION-NEGATIVE.  THIS REPO uses COMPRESSION-POSITIVE (sigma_n = +120 MPa);
  flip the sign of the n-stress column when overlaying.
- **traction is single-valued:** h-stress / v-stress / n-stress are ~the same on the
  near and far files (the fault traction is continuous); the per-side DIFFERENCE is
  in disp / vel (the two sides move differently — that is the bi-material signature).
- **physical y-side / handedness:** spec p.3 fixes near=strong / far=weak but does
  NOT pin which is +y vs -y.  This repo's `tpv6.toml` puts near(strong) on the y<0
  side (`[material.halfspace_across_fault] n_y = -1`).  Whether that matches drdg3d's
  along-strike rupture DIRECTIVITY (the ±12 km stations are asymmetric) is a FRONTERA
  validation: if the rupture is strike-mirrored vs drdg3d, flip `n_y` to +1 (swaps the
  strong/weak y-sides).  The x2=0 (hypocenter) station is symmetric and matches either way.

## What this run produces (to compare against `scec_drdg3d/`)
- Fault output (ParaView, `--paraview-fault-*`): slip-rate contours — ALWAYS produced.
- Per-side on-fault station traces (per-side disp/vel + traction): produced by the
  per-side station writer `dynamic/tpv6_stations.hpp` — the documented-deferred Part-C
  remainder (additive; the run is functional without it).  When it lands it must emit
  the strong side as "nearside", the weak side as "farside", and flip n-stress to
  compression-positive for a direct overlay (see `tpv6/visualize_results.py`).

## TPV7
TPV7 is the SCEC ILL-POSED low-contrast case (grid-dependent BY DESIGN); compare to
the SCEC ENSEMBLE SPREAD, not a single reference (Harris/Day et al. 2007).
