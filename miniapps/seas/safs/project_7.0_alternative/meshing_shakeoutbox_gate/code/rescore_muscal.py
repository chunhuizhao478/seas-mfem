#!/usr/bin/env python3
"""rescore_muscal.py -- re-score a census failure dump against MUSCAL.

The dump already holds every cell that failed the gate on the DECK nc, with its
barycentre, max edge and the deck's Vs.  Re-sampling those same barycentres from
MUSCAL answers the only question that matters before spending another 20M cells:
how many of those failures are real, and how many are an artefact of the deck
file's 250 m vertical binning?
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_intermediate" / "code"))
from collar_lib import VsGrid            # noqa: E402
from muscal_vs import MuscalVs           # noqa: E402

CVM = ("/Users/chunhuizhao/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_"
       "CASE2_intermediate_plast_phi30_40_gradedfw_k1p40_attenuation_deep40km/safs_material_cvm.nc")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--cvm", default=CVM)
    ap.add_argument("--gate", type=float, action="append", required=True)
    a = ap.parse_args()

    S = np.load(a.dump)
    B, DX, VD, RD = S["bary"], S["dx"], S["vs"], S["ratio"]
    SURF = S["on_surface"]
    print(f"{a.dump}: {len(B):,} cells that fail the deck-nc gate "
          f"{float(S['gate'])}\n")

    mus = MuscalVs(backup=VsGrid(a.cvm))
    VM = mus.at(B.astype(np.float64))
    mus.report()
    RM = VM / DX

    print(f"\n  Vs at the same barycentres:")
    print(f"    deck nc  min {VD.min():>7,.0f}  p10 {np.percentile(VD,10):>7,.0f}  "
          f"med {np.median(VD):>7,.0f}")
    print(f"    MUSCAL   min {np.nanmin(VM):>7,.0f}  p10 {np.nanpercentile(VM,10):>7,.0f}  "
          f"med {np.nanmedian(VM):>7,.0f}")
    ratio = VM / np.maximum(VD, 1e-9)
    print(f"    MUSCAL / deck   p10 {np.nanpercentile(ratio,10):.2f}  "
          f"med {np.nanmedian(ratio):.2f}  p90 {np.nanpercentile(ratio,90):.2f}")

    for g in sorted(a.gate):
        still = np.isfinite(RM) & (RM < g)
        print(f"\n  gate {g}: {int(still.sum()):,} of {len(B):,} still fail on MUSCAL "
              f"({100*still.mean():.1f} % of the deck-nc failures)")
        if still.any():
            print(f"    worst {np.nanmin(RM[still]):.4f}   "
                  f"on the free surface {int(SURF[still].sum()):,}")
            z = B[still][:, 2]
            print(f"    depth  med {np.median(z):,.0f}  p90 {np.percentile(z,90):,.0f}")
            for lo, hi in ((0, -125), (-125, -500), (-500, -2000), (-2000, -40001)):
                m = (z <= lo) & (z > hi)
                if m.any():
                    print(f"      {lo:>7,} .. {hi:>8,} m : {int(m.sum()):>9,}")

    # how much of the deck-nc demand is pure binning: cells whose barycentre is
    # in the deck's z=0 bin but which MUSCAL resolves as competent rock
    inbin = B[:, 2] > -125.0
    print(f"\n  cells with barycentre in the deck's z=0 bin (|z| < 125 m): "
          f"{int(inbin.sum()):,} ({100*inbin.mean():.1f} %)")
    if inbin.any():
        print(f"    deck Vs med {np.median(VD[inbin]):,.0f} m/s  ->  "
              f"MUSCAL Vs med {np.nanmedian(VM[inbin]):,.0f} m/s "
              f"({np.nanmedian(VM[inbin])/np.median(VD[inbin]):.2f}x)")


if __name__ == "__main__":
    main()
