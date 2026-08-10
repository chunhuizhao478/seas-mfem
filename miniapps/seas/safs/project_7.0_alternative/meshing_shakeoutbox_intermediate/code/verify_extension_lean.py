#!/usr/bin/env python3
"""verify_extension_lean.py -- the same acceptance gates as verify_extension.py,
but STREAMING, for the 149M-tet heavy mesh.

verify_extension.py holds both meshes' `connect` as int64 simultaneously:
3.9 GB (parent) + 4.8 GB (new) plus geometry and boundary, and it peaked at
15.9 GB with 75 MB of free RAM on this machine before being killed.  Nothing
here ever materialises a full connect array -- every pass reads chunks straight
from HDF5 and holds only its own window.

Checks (same IDs as verify_extension.py):
  P1  parent block bit-identical      geometry + connect, compared chunk by chunk
  P2  fault triangle multiset         identical to the parent's
  P3  fault area delta                exactly 0.000e+00 m^2
  P4  free surface flat + area
  P6  no inverted tets
  P7  domain covers the ShakeOut v1 bbox
  P8  volume gate, parent block vs collar block
  P9  quality, parent block vs collar block

P5 is deliberately absent: the independent hull census sorts 4x nt faces
(~12 GB here).  merge_collar.py checks the equivalent identity in O(1) and that
identity was cross-validated against the full census on the small and
intermediate meshes.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collar_lib import (BC_DYNAMIC_RUPTURE, BC_FREE_SURFACE, LOCAL_FACES,
                        SHAKEOUT_E, SHAKEOUT_N, VsGrid, face_code,
                        tet_edge_lengths, tet_eta, tet_signed_volume)

EI = np.array([0, 0, 0, 1, 1, 2])
EJ = np.array([1, 2, 3, 2, 3, 3])
CH = 2_000_000
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<34} {detail}", flush=True)


def faces_of_code(path, code, chunk=CH):
    """Stream out every face carrying `code`, as vertex-index triples."""
    out = []
    with h5py.File(path, "r") as f:
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        for s0 in range(0, nt, chunk):
            s1 = min(s0 + chunk, nt)
            Bc = B[s0:s1]
            hit = np.zeros(s1 - s0, bool)
            for s in range(4):
                hit |= face_code(Bc, s) == code
            if not hit.any():
                continue
            idx = np.nonzero(hit)[0]
            Tc = conn[s0:s1][idx].astype(np.int64)
            Bs = Bc[idx]
            for s in range(4):
                m = face_code(Bs, s) == code
                if m.any():
                    out.append(Tc[m][:, list(LOCAL_FACES[s])])
    return np.vstack(out) if out else np.zeros((0, 3), np.int64)


def canon(P, tri):
    W = np.round(P[tri] * 1e3).astype(np.int64)
    i = np.lexsort((W[:, :, 2], W[:, :, 1], W[:, :, 0]), axis=1)
    W = np.take_along_axis(W, i[:, :, None], axis=1).reshape(len(W), 9)
    return W[np.lexsort([W[:, k] for k in range(8, -1, -1)])]


def tri_area(P, tri):
    V = P[tri]
    return np.linalg.norm(np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]), axis=1) * 0.5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--cvm", default=None)
    ap.add_argument("--gate", type=float, default=0.6667)
    ap.add_argument("--order", type=int, default=4)
    a = ap.parse_args()

    with h5py.File(a.parent, "r") as f:
        G0 = f["geometry"][:]
        nt0 = f["connect"].shape[0]
    with h5py.File(a.new, "r") as f:
        G1 = f["geometry"][:]
        nt1 = f["connect"].shape[0]
    nv0 = len(G0)
    print(f"parent {nt0:,} tets / {nv0:,} verts")
    print(f"new    {nt1:,} tets / {len(G1):,} verts   "
          f"(+{nt1-nt0:,} = +{100*(nt1-nt0)/nt0:.2f} %)\n", flush=True)

    # P1 -- geometry prefix, then connect chunk by chunk
    ok = len(G1) >= nv0 and np.array_equal(G1[:nv0], G0)
    if ok:
        with h5py.File(a.parent, "r") as fp, h5py.File(a.new, "r") as fn:
            cp, cn = fp["connect"], fn["connect"]
            for s0 in range(0, nt0, CH):
                s1 = min(s0 + CH, nt0)
                if not np.array_equal(cp[s0:s1], cn[s0:s1]):
                    ok = False
                    break
    check("P1 parent block bit-identical", ok, "geometry + connect streamed")

    # P2 / P3 -- the fault
    f0 = faces_of_code(a.parent, BC_DYNAMIC_RUPTURE)
    f1 = faces_of_code(a.new, BC_DYNAMIC_RUPTURE)
    c0, c1 = canon(G0, f0), canon(G1, f1)
    check("P2 fault triangle multiset", c0.shape == c1.shape and np.array_equal(c0, c1),
          f"{len(f0):,} -> {len(f1):,} facets")
    a0, a1 = tri_area(G0, f0).sum(), tri_area(G1, f1).sum()
    check("P3 fault area delta == 0", (a1 - a0) == 0.0,
          f"{a0/1e6:,.6f} -> {a1/1e6:,.6f} km2   delta {a1-a0:.3e} m2")
    del f0, f1, c0, c1

    # P4 -- free surface
    s1_ = faces_of_code(a.new, BC_FREE_SURFACE)
    z = G1[s1_][:, :, 2]
    lo, hi = G1[:, 0].min(), G1[:, 0].max()
    ylo, yhi = G1[:, 1].min(), G1[:, 1].max()
    foot = (hi - lo) * (yhi - ylo) / 1e6
    area = tri_area(G1, s1_).sum() / 1e6
    check("P4 free surface flat", float(np.ptp(z)) < 1e-6,
          f"z spread {np.ptp(z):.3e} m at z = {z.flat[0]:.6f}")
    check("P4 free surface area == footprint", abs(area - foot) / foot < 1e-6,
          f"{area:,.1f} km2 vs box {foot:,.1f} km2")
    del s1_, z

    # P7
    check("P7 covers ShakeOut v1 bbox",
          lo <= SHAKEOUT_E[0] and hi >= SHAKEOUT_E[1]
          and ylo <= SHAKEOUT_N[0] and yhi >= SHAKEOUT_N[1],
          f"E {lo:,.0f}..{hi:,.0f} | N {ylo:,.0f}..{yhi:,.0f}")

    # P6 / P8 / P9 -- one streaming pass over the new mesh
    vs = VsGrid(a.cvm) if a.cvm else None
    p = a.order - 1
    neg = 0
    blk = {"parent": dict(worst=np.inf, nfail=0, n=0, eta=np.inf, n05=0, n1=0,
                          ed=np.inf, med=[]),
           "collar": dict(worst=np.inf, nfail=0, n=0, eta=np.inf, n05=0, n1=0,
                          ed=np.inf, med=[])}
    with h5py.File(a.new, "r") as f:
        conn = f["connect"]
        for s0 in range(0, nt1, CH):
            s1 = min(s0 + CH, nt1)
            T = conn[s0:s1].astype(np.int64)
            P = G1[T]
            neg += int((tet_signed_volume(G1, T) < 0).sum())
            et = tet_eta(G1, T)
            ed = tet_edge_lengths(G1, T)
            dmax = ed.max(1)
            r = (vs.at(P.mean(1)) / dmax) if vs else None
            for nm, m in (("parent", np.arange(s0, s1) < nt0),
                          ("collar", np.arange(s0, s1) >= nt0)):
                if not m.any():
                    continue
                d = blk[nm]
                d["n"] += int(m.sum())
                d["eta"] = min(d["eta"], float(et[m].min()))
                d["n05"] += int((et[m] < 0.05).sum())
                d["n1"] += int((et[m] < 0.1).sum())
                d["ed"] = min(d["ed"], float(ed[m].min()))
                d["med"].append(et[m][::29])
                if vs:
                    d["worst"] = min(d["worst"], float(r[m].min()))
                    d["nfail"] += int((r[m] < a.gate).sum())
            del T, P, et, ed
    check("P6 no inverted tets", neg == 0, f"{neg} negative-volume tets")

    if vs:
        print(f"\n  gate f = Vs/dx (gate {a.gate}, = {a.gate*p/4:.4f} Hz at p{p}):")
        for nm in ("parent", "collar"):
            d = blk[nm]
            print(f"    {nm:<7} block  worst {d['worst']:.4f}  "
                  f"failures {d['nfail']:,} of {d['n']:,}")
        if blk["parent"]["nfail"] == 0:
            check("P8 gate held in the collar", blk["collar"]["nfail"] == 0,
                  f"collar worst {blk['collar']['worst']:.4f}")
        else:
            check("P8 collar no worse than parent",
                  blk["collar"]["worst"] >= blk["parent"]["worst"],
                  f"parent worst {blk['parent']['worst']:.4f} "
                  f"({blk['parent']['nfail']:,} pre-existing failures), "
                  f"collar worst {blk['collar']['worst']:.4f}")

    print("\n  quality (parent block | collar block):")
    for k, lab in (("eta", "eta_min"), ("n05", "eta<0.05"), ("n1", "eta<0.1"),
                   ("ed", "min edge")):
        print(f"    {lab:<9} {blk['parent'][k]:>14,.4f} | {blk['collar'][k]:>14,.4f}")
    for nm in ("parent", "collar"):
        blk[nm]["medv"] = float(np.median(np.concatenate(blk[nm]["med"])))
    print(f"    {'eta_med':<9} {blk['parent']['medv']:>14,.4f} | "
          f"{blk['collar']['medv']:>14,.4f}")
    check("P9 collar adds no short edge", blk["collar"]["ed"] >= blk["parent"]["ed"],
          f"parent min edge {blk['parent']['ed']:,.2f} m vs "
          f"collar {blk['collar']['ed']:,.2f} m")
    check("P9 collar has no eta < 0.05 sliver", blk["collar"]["n05"] == 0,
          f"collar eta_min {blk['collar']['eta']:.4f}, "
          f"eta<0.05 {blk['collar']['n05']:,}, eta<0.1 {blk['collar']['n1']:,}")

    bad = [n for n, ok in RESULTS if not ok]
    print(f"\n{'ALL CHECKS PASS' if not bad else 'FAILED: ' + ', '.join(bad)}"
          f"   ({len(RESULTS)-len(bad)}/{len(RESULTS)})")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
