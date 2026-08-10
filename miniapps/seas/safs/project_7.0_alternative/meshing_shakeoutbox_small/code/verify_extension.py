#!/usr/bin/env python3
"""verify_extension.py -- acceptance gates for a collar-extended ALT mesh.

The claim being tested is "the parent is untouched and the domain now covers
ShakeOut".  Each check is stated so a FAIL is unambiguous.

  P1  parent block bit-identical      geometry[:nv0], connect[:nt0], boundary[:nt0]
  P2  fault triangle multiset         identical to the parent's
  P3  fault area delta                exactly 0.000e+00 m^2
  P4  free surface flat + area        z constant; area = new footprint
  P5  BC round-trip                   tagged faces == hull + 2 x fault
  P6  no inverted tets
  P7  domain covers the ShakeOut v1 bbox
  P8  volume frequency gate           f = Vs/dx, worst + failure count
  P9  quality vs parent               eta and edge stats, regressions flagged

Usage:
    python verify_extension.py --parent <p.puml.h5> --new <n.puml.h5> \
        [--cvm <cvm.nc> --gate 0.6667]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_lib import (BC_ABSORBING, BC_DYNAMIC_RUPTURE, BC_FREE_SURFACE,
                        LOCAL_FACES, SHAKEOUT_E, SHAKEOUT_N, VsGrid, face_code,
                        faces_with_code, tet_edge_lengths, tet_eta,
                        tet_signed_volume)

EI = np.array([0, 0, 0, 1, 1, 2])
EJ = np.array([1, 2, 3, 2, 3, 3])
CH = 3_000_000
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<34} {detail}")


def read(path):
    import h5py
    with h5py.File(path, "r") as f:
        return (f["geometry"][:], f["connect"][:].astype(np.int64),
                f["boundary"][:].astype(np.int32))


def tri_area(P, tri):
    V = P[tri]
    return np.linalg.norm(np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]), axis=1) * 0.5


def canon_tris(P, tri):
    """Coordinate-based canonical form: order-independent, index-independent."""
    W = np.round(P[tri] * 1e3).astype(np.int64)
    idx = np.lexsort((W[:, :, 2], W[:, :, 1], W[:, :, 0]), axis=1)
    W = np.take_along_axis(W, idx[:, :, None], axis=1).reshape(len(W), 9)
    return W[np.lexsort([W[:, k] for k in range(8, -1, -1)])]


def hull_faces(tets):
    nt = len(tets)
    t = tets.astype(np.int32) if tets.max() < 2 ** 31 - 1 else tets
    faces = np.concatenate([t[:, LOCAL_FACES[s]] for s in range(4)])
    key = np.sort(faces, axis=1)
    order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    ks = key[order]
    dup = np.all(ks[1:] == ks[:-1], axis=1)
    uniq = np.ones(len(ks), bool)
    uniq[:-1] &= ~dup
    uniq[1:] &= ~dup
    return int(uniq.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--cvm", default=None)
    ap.add_argument("--gate", type=float, default=0.6667)
    ap.add_argument("--order", type=int, default=4, help="SeisSol ORDER (p = order-1)")
    a = ap.parse_args()

    G0, C0, B0 = read(a.parent)
    G1, C1, B1 = read(a.new)
    nv0, nt0 = len(G0), len(C0)
    print(f"parent {nt0:,} tets / {nv0:,} verts")
    print(f"new    {len(C1):,} tets / {len(G1):,} verts   "
          f"(+{len(C1)-nt0:,} tets = +{100*(len(C1)-nt0)/nt0:.2f} %)\n")

    # P1 -- the parent block must be carried over verbatim
    ok = (len(G1) >= nv0 and len(C1) >= nt0
          and np.array_equal(G1[:nv0], G0)
          and np.array_equal(C1[:nt0], C0))
    # boundary differs ONLY where the old side wall became interior
    bdiff = np.flatnonzero(B1[:nt0] != B0)
    changed_ok = True
    for t in bdiff:
        for s in range(4):
            o = int(face_code(B0[t:t + 1], s)[0])
            n = int(face_code(B1[t:t + 1], s)[0])
            if o != n and not (o == BC_ABSORBING and n == 0):
                changed_ok = False
                break
    check("P1 parent block bit-identical", ok and changed_ok,
          f"geometry+connect identical; {len(bdiff):,} tets had an absorbing "
          f"side face turned interior")

    # P2/P3 -- the fault
    f0, _, _ = faces_with_code(C0, B0, BC_DYNAMIC_RUPTURE)
    f1, _, _ = faces_with_code(C1, B1, BC_DYNAMIC_RUPTURE)
    c0, c1 = canon_tris(G0, f0), canon_tris(G1, f1)
    check("P2 fault triangle multiset", c0.shape == c1.shape and np.array_equal(c0, c1),
          f"{len(f0):,} -> {len(f1):,} facets")
    a0, a1 = tri_area(G0, f0).sum(), tri_area(G1, f1).sum()
    check("P3 fault area delta == 0", abs(a1 - a0) == 0.0,
          f"{a0/1e6:,.6f} -> {a1/1e6:,.6f} km2   delta {a1-a0:.3e} m2")

    # P4 -- free surface
    s1, _, _ = faces_with_code(C1, B1, BC_FREE_SURFACE)
    z = G1[s1][:, :, 2]
    lo, hi = G1[:, 0].min(), G1[:, 0].max()
    ylo, yhi = G1[:, 1].min(), G1[:, 1].max()
    foot = (hi - lo) * (yhi - ylo) / 1e6
    area = tri_area(G1, s1).sum() / 1e6
    check("P4 free surface flat", float(np.ptp(z)) < 1e-6,
          f"z spread {np.ptp(z):.3e} m at z = {z.flat[0]:.6f}")
    check("P4 free surface area == footprint", abs(area - foot) / foot < 1e-6,
          f"{area:,.1f} km2 vs box {foot:,.1f} km2")

    # P5 -- BC round trip
    ft = sum(int((face_code(B1, s) == BC_DYNAMIC_RUPTURE).sum()) for s in range(4))
    tagged = sum(int((face_code(B1, s) != 0).sum()) for s in range(4))
    hf = hull_faces(C1)
    check("P5 BC round-trip", tagged == hf + ft,
          f"tagged {tagged:,} == hull {hf:,} + fault x2 {ft:,}")

    # P6 -- inverted
    neg = 0
    for s0 in range(0, len(C1), CH):
        neg += int((tet_signed_volume(G1, C1[s0:s0 + CH]) < 0).sum())
    check("P6 no inverted tets", neg == 0, f"{neg} negative-volume tets")

    # P7 -- coverage
    cov = (lo <= SHAKEOUT_E[0] and hi >= SHAKEOUT_E[1]
           and ylo <= SHAKEOUT_N[0] and yhi >= SHAKEOUT_N[1])
    check("P7 covers ShakeOut v1 bbox", cov,
          f"E {lo:,.0f}..{hi:,.0f} vs {SHAKEOUT_E[0]:,.0f}..{SHAKEOUT_E[1]:,.0f} | "
          f"N {ylo:,.0f}..{yhi:,.0f} vs {SHAKEOUT_N[0]:,.0f}..{SHAKEOUT_N[1]:,.0f}")

    # P8 -- the volume gate, judged SEPARATELY on the parent and the collar.
    #
    # A single whole-mesh verdict is the wrong test.  The collar's contract is
    # "honour the parent's own rule", and the three ALT tiers do not share one:
    # the intermediate and heavy hold f >= 0.6667 volume-wide, while the small
    # mesh never did (its own worst is 0.0782).  So: if the parent is clean the
    # collar must be clean too, and in every case the collar must be no worse
    # than the mesh it extends.
    if a.cvm:
        vs = VsGrid(a.cvm)
        p = a.order - 1

        def gate_block(lo, hi):
            worst, nfail, n = np.inf, 0, 0
            for s0 in range(lo, hi, CH):
                Tc = C1[s0:min(s0 + CH, hi)]
                Pc = G1[Tc]
                dmax = np.linalg.norm(Pc[:, EI] - Pc[:, EJ], axis=2).max(1)
                r = vs.at(Pc.mean(1)) / dmax
                worst = min(worst, float(r.min()))
                nfail += int((r < a.gate).sum())
                n += len(Tc)
            return worst, nfail, n

        wp, fp, npar = gate_block(0, nt0)
        wc, fc, ncol = gate_block(nt0, len(C1))
        print(f"\n  gate f = Vs/dx (gate {a.gate}, = {a.gate*p/4:.4f} Hz at p{p}):")
        print(f"    parent block  worst {wp:.4f}  failures {fp:,} of {npar:,}")
        print(f"    collar block  worst {wc:.4f}  failures {fc:,} of {ncol:,}")
        if fp == 0:
            check(f"P8 gate Vs/dx >= {a.gate} in the collar", fc == 0,
                  f"parent was clean; collar worst {wc:.4f}, {fc:,} failures")
        else:
            check("P8 collar no worse than parent", wc >= wp,
                  f"parent worst {wp:.4f} (parent has {fp:,} failures of its own "
                  f"-- this tier carries no volume gate), collar worst {wc:.4f}")

    # P9 -- quality
    def stats(G, C):
        eta_min, eta_all, ed_min, n05, n1, n3, e50, e100 = 1e9, [], 1e9, 0, 0, 0, 0, 0
        for s0 in range(0, len(C), CH):
            Tc = C[s0:s0 + CH]
            e = tet_eta(G, Tc)
            d = tet_edge_lengths(G, Tc)
            eta_min = min(eta_min, float(e.min()))
            ed_min = min(ed_min, float(d.min()))
            n05 += int((e < 0.05).sum()); n1 += int((e < 0.1).sum()); n3 += int((e < 0.3).sum())
            e50 += int((d < 50).sum()); e100 += int((d < 100).sum())
            eta_all.append(e[::17])
        return dict(eta_min=eta_min, eta_med=float(np.median(np.concatenate(eta_all))),
                    n05=n05, n1=n1, n3=n3, ed_min=ed_min, e50=e50, e100=e100)

    # P9 -- quality.  Counting slivers over the WHOLE mesh cannot be the test:
    # the counts are absolute, so adding any cells at all reads as a
    # "regression".  The parent's own cells are bit-identical (P1), so what is
    # actually under test is the COLLAR's quality -- reported on its own block.
    s0_ = stats(G0, C0)
    s1_ = stats(G1, C1)
    sc_ = stats(G1, C1[nt0:])
    print("\n  quality (parent | collar | merged):")
    for k in ("eta_min", "eta_med", "n05", "n1", "n3", "ed_min", "e50", "e100"):
        print(f"    {k:<9} {s0_[k]:>14,.4f} | {sc_[k]:>14,.4f} | {s1_[k]:>14,.4f}")
    for k in ("n05", "n1", "n3", "e50", "e100"):
        assert s1_[k] == s0_[k] + sc_[k], f"{k} does not decompose -- parent was touched"
    print("    (every count decomposes exactly as parent + collar, so the "
          "parent's own cells are untouched)")
    check("P9 collar has no eta < 0.05 sliver", sc_["n05"] == 0,
          f"collar eta_min {sc_['eta_min']:.4f}, eta<0.1 {sc_['n1']:,}, "
          f"eta<0.3 {sc_['n3']:,} of {len(C1)-nt0:,}")
    check("P9 collar eta_min >= parent eta_min", sc_["eta_min"] >= s0_["eta_min"],
          f"parent {s0_['eta_min']:.4f} vs collar {sc_['eta_min']:.4f}")
    check("P9 collar adds no short edge", sc_["ed_min"] >= s0_["ed_min"],
          f"parent min edge {s0_['ed_min']:,.2f} m vs collar {sc_['ed_min']:,.2f} m")

    bad = [n for n, ok, _ in RESULTS if not ok]
    print(f"\n{'ALL CHECKS PASS' if not bad else 'FAILED: ' + ', '.join(bad)}"
          f"   ({len(RESULTS)-len(bad)}/{len(RESULTS)})")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
