#!/usr/bin/env python3
"""shakeout_d_domain.py -- the CORRECT ShakeOut domain, and what ALT must cover.

WHY THE FIRST BOX WAS WRONG. The ShakeOut-box meshes were built to
`01_source_data/grid.xml` (lon -121.5..-114.0, lat 32..36), which
`SHAKEOUT_SOURCE_FACTS.md` says in as many words is the ShakeMap **product**
grid -- a north-up map footprint of ~307,000 km2, about 1.7x the largest domain
anyone actually simulated. It is where the output was gridded, not where the
simulation lived.

The only ShakeOut domain ever published with corner coordinates is ShakeOut-D's
(Olsen et al. 2009 GRL, Figure 1 caption):

    121W 34.5N ; 118.9511292W 36.621696N ; 116.032285W 31.082920N ;
    113.943965W 33.122341N

with two dropped decimal points restored (36.621696 / 31.082920). This script
re-verifies that restoration instead of trusting it: with the points corrected
the four corners must form an exact rectangle, sides 300 x 600 km, corner angles
90 deg, long axis azimuth 130/310, area 180,000 km2.

Then it reports the relationship to the deployed ALT footprint, because a
frozen-parent extension can only ADD: if the parent pokes outside ShakeOut-D,
the buildable domain is the UNION, not ShakeOut-D itself.
"""
import numpy as np
from pyproj import Transformer

# Olsen et al. 2009 GRL Fig. 1, decimals restored
SHAKEOUT_D_LONLAT = np.array([
    [-121.000000, 34.500000],
    [-118.9511292, 36.621696],
    [-116.032285, 31.082920],
    [-113.943965, 33.122341],
])
# deployed ALT footprint (rotated ~30 deg to the SAF strike), UTM 11N
ALT_CORNERS = np.array([[362120., 3996867.], [235495., 3777546.],
                        [599606., 3567327.], [726231., 3786648.]])
UTM11N = "EPSG:32611"


def order_rect(C):
    c = C.mean(0)
    return C[np.argsort(np.arctan2(C[:, 1] - c[1], C[:, 0] - c[0]))]


def rect_report(C, name):
    C = order_rect(C)
    s = [np.linalg.norm(C[(i + 1) % 4] - C[i]) for i in range(4)]
    ang = []
    for i in range(4):
        a = C[i - 1] - C[i]
        b = C[(i + 1) % 4] - C[i]
        ang.append(np.degrees(np.arccos(
            np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))))
    area = 0.5 * abs(sum(C[i, 0] * C[(i + 1) % 4, 1] - C[(i + 1) % 4, 0] * C[i, 1]
                         for i in range(4)))
    # azimuth of the long axis, clockwise from north
    i_long = int(np.argmax(s))
    e = C[(i_long + 1) % 4] - C[i_long]
    az = (np.degrees(np.arctan2(e[0], e[1]))) % 180.0
    print(f"  {name}")
    print(f"    sides  {'  '.join(f'{v/1000:,.3f}' for v in s)} km")
    print(f"    angles {'  '.join(f'{v:.2f}' for v in ang)} deg")
    print(f"    area   {area/1e6:,.0f} km2")
    print(f"    long-axis azimuth {az:.2f} / {(az+180)%360:.0f} deg")
    return C, area


def inside(P, C):
    C = order_rect(C)
    ins = np.ones(len(P), bool)
    for i in range(4):
        a, b = C[i], C[(i + 1) % 4]
        e = b - a
        ins &= (e[0] * (P[:, 1] - a[1]) - e[1] * (P[:, 0] - a[0])) >= 0.0
    return ins


def main():
    tr = Transformer.from_crs("EPSG:4326", UTM11N, always_xy=True)
    x, y = tr.transform(SHAKEOUT_D_LONLAT[:, 0], SHAKEOUT_D_LONLAT[:, 1])
    D = np.column_stack([x, y])
    print("SHAKEOUT-D, projected to UTM 11N")
    for (lo, la), (px, py) in zip(SHAKEOUT_D_LONLAT, D):
        print(f"    {lo:>12.6f} {la:>10.6f}  ->  E {px:>12,.1f}  N {py:>13,.1f}")
    Dr, Darea = rect_report(D, "ShakeOut-D rectangle")
    ok = (abs(Darea / 1e6 - 180000) / 180000 < 0.01)
    print(f"    -> matches the published 600 x 300 km / 180,000 km2: "
          f"{'YES' if ok else 'NO'}\n")

    Ar, Aarea = rect_report(ALT_CORNERS, "deployed ALT footprint")
    print()

    # containment both ways, on a dense boundary sample (corners alone lie)
    def densify(C, n=2000):
        C = order_rect(C)
        out = []
        for i in range(4):
            a, b = C[i], C[(i + 1) % 4]
            t = np.linspace(0, 1, n)[:, None]
            out.append(a + t * (b - a))
        return np.vstack(out)

    a_in_d = inside(densify(Ar), Dr)
    d_in_a = inside(densify(Dr), Ar)
    print(f"  ALT boundary inside ShakeOut-D : {100*a_in_d.mean():.1f} %")
    print(f"  ShakeOut-D boundary inside ALT : {100*d_in_a.mean():.1f} %")

    # area of ALT outside D, by sampling the ALT interior
    rng = np.random.default_rng(0)
    c = Ar.mean(0)
    u = Ar[1] - Ar[0]
    v = Ar[3] - Ar[0]
    s = rng.random((400000, 2))
    P = Ar[0] + s[:, :1] * u + s[:, 1:] * v
    frac_out = 1.0 - inside(P, Dr).mean()
    print(f"  ALT area OUTSIDE ShakeOut-D    : {100*frac_out:.1f} % "
          f"= {Aarea*frac_out/1e6:,.0f} km2")
    print("\n  => a frozen-parent extension can only ADD, so the buildable domain")
    print("     is ShakeOut-D UNION the ALT parent, not ShakeOut-D alone.")


if __name__ == "__main__":
    main()
