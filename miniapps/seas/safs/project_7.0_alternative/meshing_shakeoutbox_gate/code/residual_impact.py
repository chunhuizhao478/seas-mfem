#!/usr/bin/env python3
"""residual_impact.py -- what the remaining sub-gate cells actually cost you.

A "failure" is a cell whose resolved f = (p/4)*Vs/dx is below the target, so a
non-zero count does mean the mesh is not compliant EVERYWHERE.  The question
that matters is whether those cells sit anywhere that changes an answer, so
report, for each one:

  * the frequency it DOES resolve (not just that it missed);
  * whether it lies inside the ShakeOut v1 comparison footprint, which is where
    PGV is actually read;
  * its distance from the fault, and whether it is in the frozen parent block or
    out in the far-field collar.
"""
import argparse

import numpy as np

# ShakeOut v1 declared grid box, densely sampled and projected to UTM 11N
SHAKEOUT_E = (74758.0, 783423.2)
SHAKEOUT_N = (3540435.7, 3993324.3)
# the ALT parent footprint: rectangle rotated ~30 deg to the SAF strike
CORNERS = np.array([[362120., 3996867.], [235495., 3777546.],
                    [599606., 3567327.], [726231., 3786648.]])


def order_rect(C):
    c = C.mean(0)
    return C[np.argsort(np.arctan2(C[:, 1] - c[1], C[:, 0] - c[0]))]


def inside_rect(X, Y, C):
    C = order_rect(C)
    ins = np.ones(X.shape, bool)
    for i in range(4):
        a, b = C[i], C[(i + 1) % 4]
        e = b - a
        ins &= (e[0] * (Y - a[1]) - e[1] * (X - a[0])) >= 0.0
    return ins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--gate", type=float, required=True)
    ap.add_argument("--order", type=int, required=True, help="3 or 5")
    ap.add_argument("--label", default="")
    ap.add_argument("--fault-npz", default=None)
    a = ap.parse_args()

    S = np.load(a.dump)
    B, DX, VS, R = S["bary"], S["dx"], S["vs"], S["ratio"]
    m = R < a.gate
    B, DX, VS, R = B[m], DX[m], VS[m], R[m]
    p = a.order
    f_target = a.gate * p / 4.0
    f_res = R * p / 4.0

    print(f"=== {a.label or a.dump}")
    print(f"  target {f_target:.3f} Hz at p{p}  (gate Vs/dx >= {a.gate})")
    print(f"  cells below it: {len(R):,}\n")
    print("  frequency these cells DO resolve:")
    print(f"    worst {f_res.min():.3f} Hz   p10 {np.percentile(f_res,10):.3f}   "
          f"median {np.median(f_res):.3f}   p90 {np.percentile(f_res,90):.3f} Hz")
    for lo, hi in ((0.0, 0.25), (0.25, 0.35), (0.35, 0.45), (0.45, 0.5),
                   (0.5, 0.75), (0.75, 1.0)):
        k = (f_res >= lo * f_target / 0.5) & (f_res < hi * f_target / 0.5) if False else \
            (f_res >= lo) & (f_res < hi)
        if k.any():
            print(f"      {lo:.2f} - {hi:.2f} Hz : {int(k.sum()):>7,}")

    inside_so = ((B[:, 0] >= SHAKEOUT_E[0]) & (B[:, 0] <= SHAKEOUT_E[1])
                 & (B[:, 1] >= SHAKEOUT_N[0]) & (B[:, 1] <= SHAKEOUT_N[1]))
    in_parent = inside_rect(B[:, 0], B[:, 1], CORNERS)
    print(f"\n  inside the ShakeOut v1 comparison footprint : "
          f"{int(inside_so.sum()):,} of {len(R):,} ({100*inside_so.mean():.1f} %)")
    print(f"  inside the FROZEN PARENT footprint (near-fault region): "
          f"{int(in_parent.sum()):,} ({100*in_parent.mean():.1f} %)")

    # distance from the parent rectangle = a lower bound on distance to the fault,
    # since the fault lies well inside it
    C = order_rect(CORNERS)
    d = np.full(len(B), np.inf)
    for i in range(4):
        p0, p1 = C[i], C[(i + 1) % 4]
        e = p1 - p0
        t = np.clip(((B[:, :2] - p0) @ e) / (e @ e), 0.0, 1.0)
        proj = p0 + t[:, None] * e
        d = np.minimum(d, np.linalg.norm(B[:, :2] - proj, axis=1))
    d = np.where(in_parent, 0.0, d)
    print(f"\n  distance OUTSIDE the parent footprint (lower bound on "
          f"fault distance):")
    print(f"    min {d.min()/1000:.1f} km   p10 {np.percentile(d,10)/1000:.1f}   "
          f"median {np.median(d)/1000:.1f}   max {d.max()/1000:.1f} km")
    print(f"\n  barycentre depth: median {np.median(B[:,2]):,.0f} m   "
          f"deepest {B[:,2].min():,.0f} m")
    print(f"  bounding box of the residual: E {B[:,0].min():,.0f}..{B[:,0].max():,.0f}   "
          f"N {B[:,1].min():,.0f}..{B[:,1].max():,.0f}")


if __name__ == "__main__":
    main()
