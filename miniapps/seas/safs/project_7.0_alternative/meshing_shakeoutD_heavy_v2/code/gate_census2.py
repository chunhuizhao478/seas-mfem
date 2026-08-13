#!/usr/bin/env python3
"""gate_census2.py -- the resolution gate scored on BOTH velocity cubes.

Same convention as the deployed gate_census.py (raw f = Vs/dx, dx = element MAX
edge, Vs nearest-grid at the element BARYCENTRE), reported three ways:

  deck    Vs from safs_material_cvm.nc      -- what SeisSol reads
  native  Vs from MUSCAL.nc                 -- the model as published
  BOTH    Vs = min(deck, native)            -- the acceptance number

Usage: gate_census2.py <mesh.puml.h5> <deck.nc> <muscal.nc|-> <gate> [n_parent_tets]
"""
import sys
from pathlib import Path
import numpy as np, h5py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from material import Material

PAIRS = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
mesh, deck, mus = sys.argv[1], sys.argv[2], sys.argv[3]
GATE = float(sys.argv[4]); NPAR = int(sys.argv[5]) if len(sys.argv) > 5 else 0
CH = 3_000_000
BOX = (132000.0, 724000.0, 3490000.0, 4055000.0)

with h5py.File(mesh) as f:
    V = f["geometry"][:].astype(np.float64); T = f["connect"][:].astype(np.int64)
nt = len(T); print(f"elements {nt:,}  nodes {len(V):,}  (parent {NPAR:,})")
mat = Material(deck, None if mus == "-" else mus, box=BOX)

nan_tot = [0, 0, 0, 0]   # muscal_nan: all, all_n, collar, collar_n
n = dict(deck=0, native=0, both=0)
nc = dict(deck=0, native=0, both=0)
w = dict(deck=np.inf, native=np.inf, both=np.inf)
zs = []; dxs = []; vss = []; bbs = []; isp = []
for s in range(0, nt, CH):
    t = T[s:s+CH]; v = V[t]
    dx = np.max(np.stack([np.linalg.norm(v[:,b]-v[:,a],axis=1) for a,b in PAIRS],1),1)
    b = v.mean(1)
    i, j = mat._dij(b); vd = mat.Vd[mat._dk(b[:,2]), j, i]
    if mat.M is not None:
        mi, mj = mat._mij(b); vm = mat.Vm[mat._mk(b[:,2]), mj, mi]
        vn = np.where(np.isfinite(vm), vm, vd)
    else:
        vn = vd
    vb = np.minimum(vd, vn)
    coll = np.arange(s, s+len(t)) >= NPAR
    if mat.M is not None:
        nan_tot[0] += int((~np.isfinite(vm)).sum()); nan_tot[1] += len(vm)
        nan_tot[2] += int(((~np.isfinite(vm)) & coll).sum()); nan_tot[3] += int(coll.sum())
    for k, vv in (("deck", vd), ("native", vn), ("both", vb)):
        f_ = np.where(vv > 0, vv/dx, 0.0)
        n[k] += int((f_ < GATE).sum()); nc[k] += int(((f_ < GATE) & coll).sum())
        w[k] = min(w[k], float(f_.min()))
        if k == "native":
            m = f_ < GATE
            if m.any():
                zs.append(b[m,2]); dxs.append(dx[m]); vss.append(vv[m])
                bbs.append(b[m]); isp.append(~coll[m])
    del v, dx, b, vd, vn, vb
if mat.M is not None and nan_tot[1]:
    print(f"\n[muscal coverage] barycentres OUTSIDE MUSCAL (NaN -> deck value stands): "
          f"{nan_tot[0]:,} of {nan_tot[1]:,} = {100*nan_tot[0]/nan_tot[1]:.2f} %"
          + (f"   collar {nan_tot[2]:,}/{nan_tot[3]:,} = {100*nan_tot[2]/max(nan_tot[3],1):.2f} %"
             if nan_tot[3] else ""))
print(f"\n{'scored on':>8} | {'sub-gate':>12} | {'in collar':>11} | worst raw f | resolved p3 / p5")
for k in ("deck", "native", "both"):
    print(f"{k:>8} | {n[k]:12,} | {nc[k]:11,} | {w[k]:11.4f} | "
          f"{0.75*w[k]:.4f} Hz / {1.25*w[k]:.4f} Hz")
if zs:
    zs = np.concatenate(zs); dxs = np.concatenate(dxs); vss = np.concatenate(vss)
    print(f"\n min(deck,native) failures: {len(zs):,} = {100*len(zs)/nt:.4f} %")
    print(f"   depth  med {np.median(zs):9.1f}  min {zs.min():10.1f}  max {zs.max():8.1f}")
    print(f"   above z=-125 m: {int((zs>-125).sum()):,}   above z=-500 m: {int((zs>-500).sum()):,}")
    print(f"   dx  {dxs.min():.0f} .. {dxs.max():.0f} m ;  Vs  {vss.min():.0f} .. {vss.max():.0f} m/s")
    for lo, hi, lab in [(0,-125,"    0..-125"),(-125,-500," -125..-500"),
                        (-500,-3000," -500..-3k"),(-3000,-1e9,"  -3k..deep")]:
        m = (zs <= lo) & (zs > hi)
        if m.any():
            print(f"   {lab}: {int(m.sum()):9,}  dx med {np.median(dxs[m]):7.0f}  "
                  f"Vs med {np.median(vss[m]):7.0f}")
    bbs = np.vstack(bbs); isp = np.concatenate(isp)
    print(f"\n   in frozen PARENT block: {int(isp.sum()):,}   in collar: {int((~isp).sum()):,}")
    # distance to the fault -- the number the deliverable rests on
    FACE = [(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
    with h5py.File(mesh) as fh:
        B = fh["boundary"][:].astype(np.int64)
    tri = []
    for k in range(4):
        mm = ((B >> (8*k)) & 0xFF) == 3
        if mm.any(): tri.append(T[mm][:, list(FACE[k])])
    if tri:
        from scipy.spatial import cKDTree
        F = np.concatenate(tri); Pc = V[F].mean(1)
        dist, _ = cKDTree(Pc).query(bbs, k=1, workers=-1)
        print(f"   fault facets {len(F):,}")
        print("   dist to fault |  count | worst raw f")
        for lo, hi in [(0,1000),(1000,5000),(5000,10000),(10000,20000),(20000,1e9)]:
            mm = (dist>=lo) & (dist<hi)
            print(f"     {lo/1000:5.0f}-{hi/1000:<7.0f} km {int(mm.sum()):7,}")
        for r in (5000, 10000, 20000):
            mm = dist <= r
            print(f"   within {r/1000:4.0f} km: {int(mm.sum()):6,} failures"
                  + ("" if mm.sum() else "  -- CLEAN"))
