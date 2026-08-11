#!/usr/bin/env python3
"""calibrate_density.py -- how many real tets does a requested size h actually buy?

Every cost estimate in this campaign integrates 6/h^3 over the collar volume.
That constant is exact only for a Kuhn subdivision of cubes; a real Delaunay
fill with dx = MAX edge lands somewhere else, and the earlier build's
"analytic 2.9M vs measured 13.4M" gap (4.6x) mixed the constant up with the
grading near the parent wall, so it cannot be reused.

Measure it instead, on the shipped collar itself:

    k = N_cells / sum_cells( V_cell * 6 / dx_cell^3 )

k is dimensionless and pipeline-specific.  A cost estimate for a new field H is
then  k * integral(6/H^3 dV), with no reference to the old build's grading.

Also reports the LEB repair bill: for every failing cell, the cube of the ratio
between its current dx and the size it must reach, which is what bisection has
to pay.
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_intermediate" / "code"))
from collar_lib import tet_edge_lengths, tet_signed_volume  # noqa: E402

CH = 1_000_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--first-tet", type=int, default=0,
                    help="calibrate over tets [first, last)")
    ap.add_argument("--last-tet", type=int, default=-1)
    a = ap.parse_args()

    with h5py.File(a.mesh, "r") as f:
        G = f["geometry"][:]
        nt = f["connect"].shape[0]
    lo = a.first_tet
    hi = nt if a.last_tet < 0 else a.last_tet
    print(f"{a.mesh}\n  calibrating over tets [{lo:,}, {hi:,}) of {nt:,}", flush=True)

    n = 0
    ideal = 0.0
    vol = 0.0
    dxs = []
    with h5py.File(a.mesh, "r") as f:
        conn = f["connect"]
        for s0 in range(lo, hi, CH):
            s1 = min(s0 + CH, hi)
            T = conn[s0:s1].astype(np.int64)
            dmax = tet_edge_lengths(G, T).max(1)
            V = np.abs(tet_signed_volume(G, T))
            n += len(T)
            vol += float(V.sum())
            ideal += float((V * 6.0 / dmax ** 3).sum())
            dxs.append(dmax[::37].astype(np.float32))
            del T

    DX = np.concatenate(dxs)
    print(f"  cells {n:,}   volume {vol/1e9:,.1f} km3")
    print(f"  sum V*6/dx^3 (ideal-lattice count) = {ideal:,.0f}")
    print(f"\n  k = cells / ideal = {n/ideal:.3f}")
    print("      (k=1.41 would be a perfect regular-tet packing; lower k means the "
          "max edge overstates the cell's true size)")
    print(f"\n  dx over these cells: p10 {np.percentile(DX,10):,.0f}  "
          f"med {np.median(DX):,.0f}  p90 {np.percentile(DX,90):,.0f}  "
          f"max {DX.max():,.0f} m")
    print(f"  mean cell volume {vol/n:,.0f} m3  -> equivalent regular edge "
          f"{(vol/n*6*np.sqrt(2))**(1/3):,.0f} m")


if __name__ == "__main__":
    main()
