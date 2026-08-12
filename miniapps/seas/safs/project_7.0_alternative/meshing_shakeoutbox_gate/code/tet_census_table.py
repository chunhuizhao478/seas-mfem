#!/usr/bin/env python3
"""tet_census_table.py -- element counts, ALT small-domain vs ShakeOut-box.

Measures every mesh that is still on disk (never quotes a remembered number) and
normalises by footprint AREA and by VOLUME, because a raw tet count across two
different domains compares nothing: the ShakeOut box is 3.24x the plan area of
the original ALT footprint, so a mesh can triple in size while getting COARSER
per unit volume.

Footprints:
  small domain  the deployed ALT rectangle, 253.2 x 420.4 km ROTATED 30 deg to
                the SAF strike -> 106,476 km2 (half its own bbox)
  ShakeOut box  axis-aligned 714.000 x 483.000 km -> 344,862 km2
  collar        the difference, 238,386 km2
Both are 0 .. -40,000 m deep. The `small` TIER is deliberately excluded: its
ShakeOut product is only -22,101 m deep against safalt_small_deep40km's -40,000,
so the two are different lineages and pairing them reads as a collar REMOVING
cells (1,804,173 vs 2,131,818).
"""
import os

import h5py

A_SMALL = 106_476.0        # km2, rotated ALT footprint
A_BOX = 344_862.0          # km2, ShakeOut box
DEPTH = 40.0               # km

A = "/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative"

# (label, path, footprint km2, depth km) -- None path = deleted, count recorded
ROWS = [
    ("SMALL DOMAIN (deployed ALT footprint, 106,476 km2)", None, None, None),
    ("  intermediate safalt_0d5Hz_p3_deep40km",
     f"{A}/meshing_deep40km/results/safalt_0d5Hz_p3_deep40km.puml.h5", A_SMALL, DEPTH),
    ("  heavy (old)  safalt_fb200_deep40km_refine2",
     f"{A}/meshing_deep40km_refine2/results/safalt_fb200_deep40km_refine2.puml.h5", A_SMALL, DEPTH),
    ("  heavy (1Hz)  safalt_fb200_deep40km_1Hz_p5",
     f"{A}/meshing_deep40km_1Hz_p5_leb/results/safalt_fb200_deep40km_1Hz_p5.puml.h5", A_SMALL, DEPTH),
    ("SHAKEOUT BOX (344,862 km2 = 3.24x the plan area)", None, None, None),
    ("  intermediate FINAL (0.5 Hz p3 closed)",
     f"{A}/meshing_shakeoutbox_gate/results/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5",
     A_BOX, DEPTH),
    ("  heavy        FINAL (1 Hz p5)",
     f"{A}/meshing_shakeoutbox_gate/results/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5",
     A_BOX, DEPTH),
]
# superseded products, deleted from disk; counts measured before deletion
DELETED = [
    ("  intermediate v1 (uniform collar, superseded)", 29_385_401, A_BOX, DEPTH),
    ("  heavy v1 (wrong parent, superseded)", 135_567_603, A_BOX, DEPTH),
]


def count(p):
    with h5py.File(p, "r") as f:
        return f["connect"].shape[0], f["geometry"].shape[0]


def main():
    print(f"{'mesh':<52}{'tets':>15}{'verts':>14}"
          f"{'tets/km2':>11}{'tets/km3':>11}")
    print("-" * 103)
    got = {}
    for lab, p, area, dep in ROWS:
        if p is None:
            print(f"\n{lab}")
            continue
        if not os.path.exists(p):
            print(f"{lab:<52}{'(not on disk)':>15}")
            continue
        nt, nv = count(p)
        got[lab.strip().split()[0] + lab.strip().split()[1]] = nt
        print(f"{lab:<52}{nt:>15,}{nv:>14,}"
              f"{nt/area:>11,.0f}{nt/(area*dep):>11,.1f}")
        got[lab] = nt
    print("\nSUPERSEDED (deleted from disk; counts measured before deletion)")
    for lab, nt, area, dep in DELETED:
        print(f"{lab:<52}{nt:>15,}{'-':>14}{nt/area:>11,.0f}{nt/(area*dep):>11,.1f}")
        got[lab] = nt

    # the comparisons that matter
    def g(k):
        return next(v for kk, v in got.items() if isinstance(kk, str) and k in kk)

    i_small = g("safalt_0d5Hz_p3_deep40km")
    i_box = g("intermediate FINAL")
    h_old = g("safalt_fb200_deep40km_refine2")
    h_par = g("safalt_fb200_deep40km_1Hz_p5")
    h_box = g("heavy        FINAL")

    print("\n" + "=" * 103)
    print("INTERMEDIATE  small domain -> ShakeOut box")
    print(f"  tets            {i_small:>15,} -> {i_box:>15,}   "
          f"x{i_box/i_small:.3f}  (+{i_box-i_small:,})")
    print(f"  plan area                106,476 km2 ->       344,862 km2   x3.239")
    print(f"  tets per km3    {i_small/(A_SMALL*DEPTH):>15.1f} -> "
          f"{i_box/(A_BOX*DEPTH):>15.1f}   x{(i_box/(A_BOX*DEPTH))/(i_small/(A_SMALL*DEPTH)):.3f}")
    print(f"  -> the added domain is {100*(i_box-i_small)/i_box:.1f} % of the new mesh, "
          f"over {100*(A_BOX-A_SMALL)/A_BOX:.1f} % of the new area")

    print("\nHEAVY  small domain -> ShakeOut box")
    print(f"  parent used     {h_par:>15,} (1 Hz p5)   [old refine2 was {h_old:,}]")
    print(f"  tets            {h_par:>15,} -> {h_box:>15,}   "
          f"x{h_box/h_par:.3f}  (+{h_box-h_par:,})")
    print(f"  tets per km3    {h_par/(A_SMALL*DEPTH):>15.1f} -> "
          f"{h_box/(A_BOX*DEPTH):>15.1f}   x{(h_box/(A_BOX*DEPTH))/(h_par/(A_SMALL*DEPTH)):.3f}")
    print(f"  -> the added domain is {100*(h_box-h_par)/h_box:.1f} % of the new mesh, "
          f"over {100*(A_BOX-A_SMALL)/A_BOX:.1f} % of the new area")

    print("\nHEAVY vs INTERMEDIATE, on the same ShakeOut box")
    print(f"  {h_box:,} / {i_box:,} = x{h_box/i_box:.2f}")


if __name__ == "__main__":
    main()
