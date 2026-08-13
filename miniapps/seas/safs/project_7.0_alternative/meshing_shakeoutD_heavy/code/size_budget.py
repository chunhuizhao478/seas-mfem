#!/usr/bin/env python3
"""size_budget.py -- what will the ShakeOut-D heavy mesh COST, before building it.

Integrates the design size field over the 600 x 300 x 80 km domain and reports the
expected tet count, the RAM the build stages need, and the run-time footprint.
Nothing is meshed here; this is the tripwire the project rule asks for -- measure
the bill and stop for the user if it trips.

Size field, in the order applied:
  1. fault zone   h_fault(d) -- the MEASURED profile of the existing small-domain
                  heavy mesh, which is the treatment being preserved
  2. gate         h_gate = Vs_MUSCAL / gate_raw   (gate 0.8 = 1 Hz at p5, which
                  also gives >= 0.6 Hz at p3, so the 0.5 Hz p3 target is implied)
  3. cap          hmax
  4. gradation    |grad h| <= g, swept to convergence -- this is what makes the
                  transition GRADUAL and it is not free; the cost is reported per g

Cell count uses N = k * integral(dV / h^3) with k CALIBRATED on the existing mesh
(k = 18.90, against 8.485 for ideal regular tets) -- so the estimate reflects this
pipeline's real behaviour, not an idealisation.
"""
import argparse
import numpy as np, h5py
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("--fault-mesh")
ap.add_argument("--fault-npz", default="build_tmp/fault_surface.npz",
                help="preferred source of fault centroids; avoids re-reading the parent h5")
ap.add_argument("--muscal", default="/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc")
ap.add_argument("--gate", type=float, default=0.8)
ap.add_argument("--hmax", type=float, default=5000.0)
ap.add_argument("--depth", type=float, default=80000.0)
ap.add_argument("--dxy", type=float, default=1500.0)
ap.add_argument("--k", type=float, default=18.8992)
ap.add_argument("--grads", default="0.10,0.15,0.20,0.30")
ap.add_argument("--no-fault", action="store_true",
                help="drop the fault-proximity term and cost the GATE ALONE -- i.e. what "
                     "the 1 Hz-at-p5 requirement actually costs, separated from what "
                     "preserving the fault-zone treatment costs")
a = ap.parse_args()

Q = np.load("build_tmp/domain_corners_utm.npy")     # W,N,E,S
W, N, E, S = Q
u = (S - W) / np.linalg.norm(S - W); v = (N - W) / np.linalg.norm(N - W)
L1, L2 = np.linalg.norm(S - W), np.linalg.norm(N - W)
print(f"[domain] {L1/1000:.1f} x {L2/1000:.1f} km, depth {a.depth/1000:.0f} km, "
      f"bearing {np.degrees(np.arctan2(u[0],u[1]))%360:.1f} deg")

# ---- grid in rotated coords ------------------------------------------------
ns, nt_ = int(L1/a.dxy)+1, int(L2/a.dxy)+1
ss = np.linspace(0, L1, ns); tt = np.linspace(0, L2, nt_)
zs = -np.concatenate([np.arange(0, 3000, 100), np.arange(3000, 12000, 500),
                      np.arange(12000, a.depth+1, 2000)])
print(f"[grid] {ns} x {nt_} x {len(zs)} = {ns*nt_*len(zs):,} nodes "
      f"({a.dxy:.0f} m lateral)")
SS, TT = np.meshgrid(ss, tt, indexing="ij")
XY = W[None,None,:] + SS[...,None]*u + TT[...,None]*v      # (ns,nt,2)
xy = XY.reshape(-1,2)

# ---- distance to the fault -------------------------------------------------
FACEV=[(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
# The h5 path re-reads the 5.5 GB parent (boundary alone is 1.07 GB) purely to
# recover facet centroids.  build_tmp/fault_surface.npz already holds exactly those
# triangles, deduplicated -- and on this shared box the other build's mmg can be
# holding 13 GB, which is the one situation that has actually crashed this project.
# So prefer the npz when it is available; the h5 path stays as the fallback.
if a.fault_npz:
    Fn=np.load(a.fault_npz); Pc=Fn["P"][Fn["T"]].mean(1); nfacet=len(Fn["T"])
else:
    with h5py.File(a.fault_mesh) as f:
        V=f["geometry"][:].astype(np.float64); C=f["connect"]; B=f["boundary"][:].astype(np.int64)
        tri=[]
        for s in range(0, C.shape[0], 8_000_000):
            c=C[s:s+8_000_000][:].astype(np.int64); b=B[s:s+8_000_000]
            for k in range(4):
                m=((b>>(8*k))&0xFF)==3
                if m.any(): tri.append(c[m][:,list(FACEV[k])])
        F=np.vstack(tri); Pc=V[F].mean(1); nfacet=len(F)//2
# DECIMATE + TIGHT BOUND.  Querying 2.76 M facet centroids with a 60 km bound is the
# documented KD trap in this project: the bound only helps if it is TIGHT, and at a
# 1500 m grid the fault's ~115 m facet spacing is absurd resolution for a distance
# field.  Beyond ~20 km the GATE sets h anyway, so the exact distance is irrelevant.
# Decimate to a TARGET COUNT, not a fixed stride: the h5 path yields every facet
# twice (once per adjacent tet) while the npz is already deduplicated, so a hardcoded
# ::13 would halve the sample density on the npz path and quietly change the field.
Pc = Pc[::max(1, len(Pc) // 425_000)]
print(f"[fault] {nfacet:,} facets -> {len(Pc):,} decimated centroids for the distance field")
tree=cKDTree(Pc)
DCAP=25000.0

# ---- MUSCAL ---------------------------------------------------------------
import netCDF4 as ncdf
from pyproj import Transformer
m=ncdf.Dataset(a.muscal)
mlon=m["longitude"][:].data.astype(np.float64); mlat=m["latitude"][:].data.astype(np.float64)
mdep=m["depth"][:].data.astype(np.float64)
fwd=Transformer.from_crs("EPSG:32611","EPSG:4326",always_xy=True)
lo,la=fwd.transform(xy[:,0],xy[:,1])
def near(ax,val):
    i=np.clip(np.searchsorted(ax,val),0,len(ax)-1); im=np.maximum(i-1,0)
    return np.where(np.abs(ax[im]-val)<np.abs(ax[i]-val),im,i).astype(np.int32)
ii,jj = near(mlon,lo), near(mlat,la)
VS=np.asarray(m["vs"][:], np.float32)
print(f"[muscal] subset lookup ready; {100*np.isfinite(VS).mean():.2f} % finite")

# fault-zone profile (measured medians of the existing mesh)
PD = np.array([0, 125, 250, 500, 1000, 2000, 4000, 7500, 15000], float)
PH = np.array([115, 115, 141, 215, 262, 529, 700, 1000, 1500], float)

nz=len(zs)
tot_naive=0.0; vol_tot=0.0
res={g:0.0 for g in [float(x) for x in a.grads.split(",")]}
H3=np.zeros((len(xy), nz), np.float32)
DIST=np.zeros((len(xy), nz), np.float32)
for kz,z in enumerate(zs):
    d,_=tree.query(np.column_stack([xy, np.full(len(xy), z)]), k=1,
                   distance_upper_bound=DCAP, workers=-1)
    d=np.where(np.isfinite(d), d, DCAP)
    hf=np.interp(d, PD, PH, right=PH[-1])
    hf=np.where(d>15000, PH[-1]+(d-15000)*0.20, hf)
    # Beyond the distance cap the fault imposes NO constraint.  Leaving the
    # extrapolation in place instead pins the whole far field at h_fault(DCAP)
    # = 3500 m and silently inflates the count ~2.9x -- the distance cap is a
    # SEARCH optimisation, not a size rule.
    hf=np.where(d>=DCAP-1.0, np.inf, hf)
    if a.no_fault:
        hf=np.full_like(hf, np.inf)   # gate-only costing: fault imposes nothing
    kdep=near(mdep, -z)
    vs=VS[kdep, jj, ii]
    hg=np.where(np.isfinite(vs), vs/a.gate, a.hmax)
    H3[:,kz]=np.minimum(np.minimum(hf,hg), a.hmax); DIST[:,kz]=d
np.save("build_tmp/H_nograd.npy", H3)
dz=np.abs(np.gradient(zs))
cellvol=(a.dxy**2)*dz[None,:]
for g in res:
    H=H3.copy().reshape(ns,nt_,nz)
    for _ in range(200):
        b=H.copy()
        np.minimum(H[1:],  H[:-1]+g*a.dxy, out=H[1:])
        np.minimum(H[:-1], H[1:] +g*a.dxy, out=H[:-1])
        np.minimum(H[:,1:],  H[:,:-1]+g*a.dxy, out=H[:,1:])
        np.minimum(H[:,:-1], H[:,1:] +g*a.dxy, out=H[:,:-1])
        np.minimum(H[:,:,1:],  H[:,:,:-1]+g*dz[None,None,1:], out=H[:,:,1:])
        np.minimum(H[:,:,:-1], H[:,:,1:] +g*dz[None,None,1:], out=H[:,:,:-1])
        if np.abs(H-b).max()<1.0: break
    DIST3=DIST.reshape(-1,nz)
    dens=(cellvol[None,:,:].reshape(1,-1,nz)[0]/H.reshape(-1,nz)**3)
    n=a.k*float(dens.sum())
    res[g]=n
    print(f"[grad {g:.2f}] h {H.min():7.1f} .. {H.max():7.1f} m   ->  N = {n/1e6:9.1f} M tets")
    if abs(g-0.15)<1e-9:
        print("   distance-band breakdown (compare vs the EXISTING mesh's measured counts):")
        REF={(0,250):46652841,(250,500):24626584,(500,1000):13668277,
             (1000,2000):11369523,(2000,4000):2734082,(4000,7500):1549036,
             (7500,15000):612690}
        for (lo,hi),ref in REF.items():
            m=(DIST3>=lo)&(DIST3<hi)
            nb=a.k*float(dens[m].sum())
            print(f"     {lo:6d}..{hi:<6d} m : model {nb/1e6:8.2f} M   existing {ref/1e6:8.2f} M   "
                  f"ratio {nb/ref:5.2f}")
        m=DIST3>=15000
        print(f"     >15 km        : model {a.k*float(dens[m].sum())/1e6:8.2f} M   "
              f"existing {9523142/1e6:8.2f} M (smaller+shallower domain)")
print(f"\n[domain volume] {(L1*L2*a.depth)/1e9:,.0f} km3   "
      f"(existing small-domain heavy: 3,979,863 km3 = {3979863/((L1*L2*a.depth)/1e9)*100:.1f} %)")
