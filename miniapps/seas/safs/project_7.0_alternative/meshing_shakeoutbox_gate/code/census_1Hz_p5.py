#!/usr/bin/env python3
"""census_1Hz_p5.py -- locate every cell that fails the 1 Hz @ p5 resolution gate.

Conventions (locked; identical to check_fsurf_freq.py and the 0.5 Hz campaign):
    dx  = element MAX edge
    Vs  = sqrt(mu/rho) at the element BARYCENTER, NEAREST-GRID from the deck CVM
    resolved f at order p = (p/4) * Vs/dx    (p segments/element, 4 samples/wavelength)

    p3: resolved 0.5 Hz  <=>  Vs/dx >= 0.6667
    p5: resolved 1.0 Hz  <=>  Vs/dx >= 0.8000   <-- THIS GATE

Emits the census + the failing cell list and their required dx for the LEB refiner.

    python census_1Hz_p5.py --mesh X.puml.h5 --cvm safs_material_cvm.nc \
        --gate 0.8 --out results/census_round0
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent))
from puml_io import LOCAL_FACES, BC_DYNAMIC_RUPTURE  # noqa: E402

PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
NVMAX = 1 << 25
CH = 4_000_000


class Vs:
    """Nearest-grid Vs, matching check_fsurf_freq.py exactly."""

    def __init__(self, path):
        import netCDF4 as ncdf
        d = ncdf.Dataset(str(path))
        self.X, self.Y, self.Z = d["x"][:].data, d["y"][:].data, d["z"][:].data
        raw = d["data"][:]
        self.mu = np.asarray(raw["mu"], np.float64)
        self.rho = np.asarray(raw["rho"], np.float64)
        d.close()

    def at(self, p):
        ix = np.clip(np.rint((p[:, 0] - self.X[0]) / (self.X[1] - self.X[0])).astype(int),
                     0, len(self.X) - 1)
        iy = np.clip(np.rint((p[:, 1] - self.Y[0]) / (self.Y[1] - self.Y[0])).astype(int),
                     0, len(self.Y) - 1)
        iz = np.clip(np.rint((p[:, 2] - self.Z[0]) / (self.Z[1] - self.Z[0])).astype(int),
                     0, len(self.Z) - 1)
        return np.sqrt(np.maximum(self.mu[iz, iy, ix], 0.0) / self.rho[iz, iy, ix])


def pct(a, q):
    return float(np.percentile(a, q)) if len(a) else float("nan")


def fault_sets(conn, codes, nv, nt):
    """Unique fault edge keys (sorted) and the fault-vertex mask."""
    fkeys, fverts = [], []
    for s0 in range(0, nt, CH):
        T = conn[s0:s0 + CH].astype(np.int64)
        cch = codes[s0:s0 + CH]
        for s in range(4):
            sel = np.nonzero(cch[:, s] == BC_DYNAMIC_RUPTURE)[0]
            if not sel.size:
                continue
            tri = T[sel][:, LOCAL_FACES[s]]
            fverts.append(tri.ravel())
            for i, j in ((0, 1), (1, 2), (0, 2)):
                a = np.minimum(tri[:, i], tri[:, j])
                b = np.maximum(tri[:, i], tri[:, j])
                fkeys.append(a * NVMAX + b)
    fkeys = np.unique(np.concatenate(fkeys))
    is_fault_v = np.zeros(nv, bool)
    is_fault_v[np.unique(np.concatenate(fverts))] = True
    return fkeys, is_fault_v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--cvm", required=True)
    ap.add_argument("--gate", type=float, default=0.8, help="Vs/dx floor (0.8 = 1 Hz at p5)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    vs = Vs(args.cvm)
    f = h5py.File(args.mesh, "r")
    conn, bnd = f["connect"], f["boundary"]
    pts = f["geometry"][:]
    nt, nv = conn.shape[0], len(pts)
    print(f"[mesh] {nt:,} tets, {nv:,} verts", flush=True)

    bu = bnd[:].view(np.uint32)
    codes = np.stack([((bu >> np.uint32(8 * s)) & np.uint32(0xFF)).astype(np.uint8)
                      for s in range(4)], 1)
    del bu
    fkeys, is_fault_v = fault_sets(conn, codes, nv, nt)
    del codes
    print(f"[fault] {len(fkeys):,} unique fault edges, {int(is_fault_v.sum()):,} fault "
          f"vertices  ({time.time()-t0:.0f}s)", flush=True)

    ii = np.array([p[0] for p in PAIRS])
    jj = np.array([p[1] for p in PAIRS])
    fail_idx, fail_f, fail_dx, fail_vs, fail_z, cls_edge, cls_vert = [], [], [], [], [], [], []
    n_fail, fmin_all = 0, np.inf
    for s0 in range(0, nt, CH):
        s1 = min(s0 + CH, nt)
        T = conn[s0:s1].astype(np.int64)
        P = pts[T]
        E = np.linalg.norm(P[:, ii] - P[:, jj], axis=2)
        kmax = E.argmax(1)
        dx = E[np.arange(len(T)), kmax]
        bary = P.mean(1)
        v = vs.at(bary)
        fr = np.where(dx > 0, v / dx, 0.0)
        fmin_all = min(fmin_all, float(fr.min()))
        bad = np.nonzero(fr < args.gate)[0]
        if bad.size:
            a = np.minimum(T[bad, ii[kmax[bad]]], T[bad, jj[kmax[bad]]])
            b = np.maximum(T[bad, ii[kmax[bad]]], T[bad, jj[kmax[bad]]])
            key = a * NVMAX + b
            pos = np.clip(np.searchsorted(fkeys, key), 0, max(len(fkeys) - 1, 0))
            on_fault_edge = fkeys[pos] == key if len(fkeys) else np.zeros(len(key), bool)
            pinned = (is_fault_v[a] | is_fault_v[b]) & ~on_fault_edge
            fail_idx.append(bad + s0)
            fail_f.append(fr[bad]); fail_dx.append(dx[bad]); fail_vs.append(v[bad])
            fail_z.append(bary[bad][:, 2])
            cls_edge.append(on_fault_edge); cls_vert.append(pinned)
            n_fail += bad.size
        if (s0 // CH) % 8 == 0:
            print(f"   {s1:,}/{nt:,}  fails {n_fail:,}  ({time.time()-t0:.0f}s)", flush=True)
    f.close()

    base = dict(mesh=args.mesh, cvm=args.cvm, gate=args.gate, n_tets=int(nt),
                worst_vs_over_dx=fmin_all,
                resolved_everywhere_p3_Hz=0.75 * fmin_all,
                resolved_everywhere_p5_Hz=1.25 * fmin_all)
    if n_fail == 0:
        base.update(n_fail=0)
        json.dump(base, open(str(out) + ".json", "w"), indent=2)
        print(f"\n[GATE] PASS -- 0 cells below Vs/dx = {args.gate}; worst {fmin_all:.4f}")
        print(f"  resolved everywhere: {0.75*fmin_all:.4f} Hz at p3, "
              f"{1.25*fmin_all:.4f} Hz at p5")
        return

    idx = np.concatenate(fail_idx); fr = np.concatenate(fail_f)
    dxv = np.concatenate(fail_dx); vsv = np.concatenate(fail_vs); zv = np.concatenate(fail_z)
    ce = np.concatenate(cls_edge); cv = np.concatenate(cls_vert)
    need = vsv / args.gate
    short = dxv / need
    o = np.argsort(fr)
    base.update(
        n_fail=int(n_fail), fail_frac=n_fail / nt,
        classes=dict(fault_edge_exempt=int(ce.sum()), fault_vertex_pinned=int(cv.sum()),
                     free=int((~ce & ~cv).sum())),
        shortfall=dict(med=pct(short, 50), p90=pct(short, 90), p99=pct(short, 99),
                       max=float(short.max())),
        vs=dict(min=float(vsv.min()), med=pct(vsv, 50), max=float(vsv.max())),
        dx=dict(med=pct(dxv, 50), max=float(dxv.max())),
        need_dx=dict(med=pct(need, 50), min=float(need.min())),
        z=dict(p05=pct(zv, 5), med=pct(zv, 50), p95=pct(zv, 95),
               min=float(zv.min()), max=float(zv.max())),
        worst10=[dict(f=float(fr[k]), dx=float(dxv[k]), vs=float(vsv[k]), z=float(zv[k]),
                      fault_edge=bool(ce[k]), fault_vertex=bool(cv[k])) for k in o[:10]])
    edges = np.array([-40000, -20000, -10000, -5000, -2000, -1000, -500, -250, -125, 1.0])
    hb = np.bincount(np.clip(np.searchsorted(edges, zv) - 1, 0, 8), minlength=9)
    base["depth_bins"] = [[float(edges[i]), float(edges[i + 1]), int(hb[i])] for i in range(9)]

    np.save(str(out) + "_failcells.npy", idx.astype(np.int64))
    np.save(str(out) + "_needdx.npy", need.astype(np.float64))
    json.dump(base, open(str(out) + ".json", "w"), indent=2)

    print(f"\n[GATE] FAIL -- {n_fail:,} of {nt:,} cells ({100*n_fail/nt:.4f}%) below "
          f"Vs/dx = {args.gate}")
    print(f"  worst Vs/dx {fr.min():.4f}  ->  resolved {0.75*fmin_all:.4f} Hz at p3, "
          f"{1.25*fmin_all:.4f} Hz at p5")
    print(f"  classes: fault-edge exempt {int(ce.sum()):,} | fault-vertex pinned "
          f"{int(cv.sum()):,} | free {int((~ce&~cv).sum()):,}")
    print(f"  shortfall dx/need: med {pct(short,50):.3f}  p90 {pct(short,90):.3f}  "
          f"p99 {pct(short,99):.3f}  max {short.max():.3f}")
    print(f"  Vs  min {vsv.min():.1f}  med {pct(vsv,50):.1f} m/s")
    print(f"  dx  med {pct(dxv,50):.1f} m   need med {pct(need,50):.1f} m")
    print(f"  z   p05 {pct(zv,5):.0f}  med {pct(zv,50):.0f}  p95 {pct(zv,95):.0f} m")
    print("  depth bins (barycenter z):")
    for lo, hi, n in base["depth_bins"]:
        print(f"     [{lo:>7.0f},{hi:>7.0f})  {n:>10,}  {100*n/n_fail:5.1f}%")
    print(f"\n[out] {out}.json / _failcells.npy / _needdx.npy   ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
