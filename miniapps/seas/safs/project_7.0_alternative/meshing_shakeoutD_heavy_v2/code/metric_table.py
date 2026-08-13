#!/usr/bin/env python3
"""metric_table.py -- cost of PRESERVING the fault band while flooring the rest.

A global metric floor is the wrong lever here: the fault treatment being preserved
IS a 115 m triangulation, and flooring it to 900 m would ask mmg to hang 900 m tets
off frozen 115 m faces -- anisotropic slivers by construction.  The honest lever is
"keep the measured profile out to R, floor beyond it", so this reports N(R, floor).
"""
import numpy as np, netCDF4 as ncdf
from scipy.spatial import cKDTree
from pyproj import Transformer
d = np.load("build_tmp/fill.npz"); P, T, plcT, MARK = d["P"], d["T"], d["plcT"], d["MARK"]
Pc = P[plcT[MARK == 7]].mean(1)[::13]
dist, _ = cKDTree(Pc).query(P, k=1, distance_upper_bound=25000.0, workers=-1)
dist = np.where(np.isfinite(dist), dist, 25000.0)
PD = np.array([0,125,250,500,1000,2000,4000,7500,15000], float)
PH = np.array([115,115,141,215,262,529,700,1000,1500], float)
hf = np.interp(dist, PD, PH, right=PH[-1])
hf = np.where(dist > 15000, PH[-1]+(dist-15000)*0.20, hf)
hf = np.where(dist >= 24999.0, np.inf, hf)
m = ncdf.Dataset("/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc")
mlon=m["longitude"][:].data.astype(float); mlat=m["latitude"][:].data.astype(float)
mdep=m["depth"][:].data.astype(float); VS=np.asarray(m["vs"][:], np.float32)
lo,la = Transformer.from_crs("EPSG:32611","EPSG:4326",always_xy=True).transform(P[:,0],P[:,1])
def near(ax,v):
    i=np.clip(np.searchsorted(ax,v),0,len(ax)-1); im=np.maximum(i-1,0)
    return np.where(np.abs(ax[im]-v)<np.abs(ax[i]-v),im,i).astype(np.int32)
vs = VS[near(mdep,np.maximum(-P[:,2],0.)), near(mlat,la), near(mlon,lo)].astype(float)
vs = np.where(np.isfinite(vs), vs, np.nanmedian(vs))
h0 = np.minimum(np.minimum(hf, vs/0.8), 5000.0)
V=P[T]; vol=np.abs(np.einsum('ij,ij->i',V[:,1]-V[:,0],np.cross(V[:,2]-V[:,0],V[:,3]-V[:,0])))/6.
def N(h): return 18.8992*np.sum(vol/h[T].mean(1)**3)/1e6
print(f"full design field (no floor): {N(h0):8.2f} M tets")
print("\n         floor beyond R  ->  M tets")
print("  R (m) |" + "".join(f"{f:>9.0f}" for f in (250,400,600,900,1400)))
print("  ------+" + "-"*45)
for R in (0, 250, 500, 1000, 2000, 4000):
    row = []
    for f in (250,400,600,900,1400):
        row.append(N(np.where(dist <= R, h0, np.maximum(h0, f))))
    print(f"  {R:5.0f} |" + "".join(f"{x:9.2f}" for x in row))
