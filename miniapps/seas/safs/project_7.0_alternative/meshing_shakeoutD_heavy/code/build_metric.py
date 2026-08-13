#!/usr/bin/env python3
"""build_metric.py -- isotropic vertex metric for the mmg pass on the base fill.

Same design field as size_budget.py, evaluated at the base mesh's own vertices:
    h = min( h_fault(d) , Vs_MUSCAL / gate , hmax )   then floored at --hmin

Two deliberate choices:
  * Vs is the vertex's NEAREST-GRID MUSCAL value, not a pooled minimum.  mmg only
    has to get close; exact gate acceptance (barycentre Vs vs MAX edge) is closed
    afterwards by the chunked LEB pass, which is the tool that can actually see
    the acceptance metric.  Pooling here would double-count the safety margin.
  * The field is supplied UNGRADED.  mmg enforces gradation itself via -hgrad, and
    pre-grading would fight it.

--hmin is the cost lever: the fault triangles are frozen at their 115 m spacing
either way, so flooring the metric coarsens only the transition and far field.
Prints the predicted tet count for a sweep of floors before writing anything.
"""
import argparse, numpy as np, netCDF4 as ncdf
from scipy.spatial import cKDTree
from pyproj import Transformer

ap = argparse.ArgumentParser()
ap.add_argument("--muscal", default="/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc")
ap.add_argument("--gate", type=float, default=0.8)
ap.add_argument("--hmax", type=float, default=5000.0)
ap.add_argument("--hmin", type=float, default=0.0, help="0 = report sweep only")
ap.add_argument("--k", type=float, default=18.8992)
ap.add_argument("--out", default="build_tmp/metric.sol")
a = ap.parse_args()

d = np.load("build_tmp/fill.npz")
P, T, plcT, MARK = d["P"], d["T"], d["plcT"], d["MARK"]
print(f"[base] {len(T):,} tets  {len(P):,} verts")

Pc = plcT[MARK == 7]
Pc = P[Pc].mean(1)[::13]
tree = cKDTree(Pc); DCAP = 25000.0
print(f"[fault] {int((MARK==7).sum()):,} facets -> {len(Pc):,} decimated centroids")
dist, _ = tree.query(P, k=1, distance_upper_bound=DCAP, workers=-1)
dist = np.where(np.isfinite(dist), dist, DCAP)

PD = np.array([0, 125, 250, 500, 1000, 2000, 4000, 7500, 15000], float)
PH = np.array([115, 115, 141, 215, 262, 529, 700, 1000, 1500], float)
hf = np.interp(dist, PD, PH, right=PH[-1])
hf = np.where(dist > 15000, PH[-1] + (dist - 15000) * 0.20, hf)
hf = np.where(dist >= DCAP - 1.0, np.inf, hf)

m = ncdf.Dataset(a.muscal)
mlon = m["longitude"][:].data.astype(np.float64)
mlat = m["latitude"][:].data.astype(np.float64)
mdep = m["depth"][:].data.astype(np.float64)
VS = np.asarray(m["vs"][:], np.float32)
fwd = Transformer.from_crs("EPSG:32611", "EPSG:4326", always_xy=True)
lo, la = fwd.transform(P[:, 0], P[:, 1])
def near(ax, val):
    i = np.clip(np.searchsorted(ax, val), 0, len(ax) - 1); im = np.maximum(i - 1, 0)
    return np.where(np.abs(ax[im] - val) < np.abs(ax[i] - val), im, i).astype(np.int32)
ii, jj, kk = near(mlon, lo), near(mlat, la), near(mdep, np.maximum(-P[:, 2], 0.0))
vs = VS[kk, jj, ii].astype(np.float64)
nfin = np.isfinite(vs)
print(f"[muscal] finite at {100*nfin.mean():.2f} % of base vertices "
      f"(Vs {np.nanmin(vs):.0f}..{np.nanmax(vs):.0f} m/s)")
vs = np.where(nfin, vs, np.nanmedian(vs))
hg = vs / a.gate

h = np.minimum(np.minimum(hf, hg), a.hmax)
lim = np.where(hf <= np.minimum(hg, a.hmax), "fault",
      np.where(hg <= a.hmax, "gate", "cap"))
print(f"[field] h {h.min():.0f} .. {h.max():.0f} m   median {np.median(h):.0f}")
for nm in ("fault", "gate", "cap"):
    print(f"   limited by {nm:6s}: {100*(lim==nm).mean():5.2f} %")

# predicted count: N = k * sum(Vol / h_bary^3), h_bary = harmonic-ish vertex mean
V = P[T]
vol = np.abs(np.einsum('ij,ij->i', V[:,1]-V[:,0], np.cross(V[:,2]-V[:,0], V[:,3]-V[:,0])))/6.0
print(f"[predict] k={a.k}  base volume {vol.sum()/1e9:,.0f} km3")
for floor in (0, 150, 250, 400, 600, 900, 1400, 2000):
    hh = np.maximum(h, floor) if floor else h
    hb = hh[T].mean(1)
    N = a.k * np.sum(vol / hb**3)
    print(f"   hmin {floor:5.0f} m -> {N/1e6:9.2f} M tets"
          + ("   <== requested" if floor == a.hmin else ""))

if a.hmin > 0:
    hh = np.maximum(h, a.hmin)
    with open(a.out, "w") as f:
        f.write("MeshVersionFormatted 2\nDimension 3\n\nSolAtVertices\n")
        f.write(f"{len(hh)}\n1 1\n\n")
        np.savetxt(f, hh.reshape(-1, 1), fmt="%.6g")
        f.write("\nEnd\n")
    print(f"[out] {a.out}  ({len(hh):,} values, hmin {a.hmin:.0f} m)")
