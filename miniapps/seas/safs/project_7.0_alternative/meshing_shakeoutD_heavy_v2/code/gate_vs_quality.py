#!/usr/bin/env python3
"""gate_vs_quality.py -- are the deep gate failures a SIZING miss or SLIVERS?

The census found 389,055 failing cells below 10 km with mean max-edge 11,307 m,
against a metric capped at 5000 m -- edges that long should not exist if the
metric were being honoured isotropically.  A stretched tet has one long edge and a
small inradius, so it fails the gate (dx = MAX edge) while its nominal size is fine.
This cross-tabulates the gate against Joe-Liu eta to tell the two apart, because
the fix differs: refine (LEB) for a sizing miss, re-optimise for slivers -- and
bisecting a sliver just makes two slivers.
"""
import numpy as np, pandas as pd, netCDF4 as ncdf
from pyproj import Transformer
import sys
F = sys.argv[1] if len(sys.argv) > 1 else "build_tmp/s1.mmg_out.mesh"
# MEASURED, not hardcoded -- see medit_hdr.py. The inherited constants
# (LV,LT = 7, 13428231 / NV,NT = 12042113, 70641218) are PREFERRED's layout and
# silently mis-parse any other mesh rather than failing.
from medit_hdr import medit_sections
_sec = medit_sections(F)
LV, NV = _sec["Vertices"]
LT, NT = _sec["Tetrahedra"]
GATE = 0.8
print(f"[hdr] Vertices {NV:,} @ {LV:,}   Tetrahedra {NT:,} @ {LT:,}")
PAIRS = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
P = pd.read_csv(F, sep=r"\s+", header=None, skiprows=LV+1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine="c").to_numpy()
m = ncdf.Dataset("/Users/chunhuizhao/Downloads/muscal_nc/MUSCAL.nc")
mlon=m["longitude"][:].data.astype(float); mlat=m["latitude"][:].data.astype(float)
mdep=m["depth"][:].data.astype(float); VS=np.asarray(m["vs"][:], np.float32)
fwd = Transformer.from_crs("EPSG:32611","EPSG:4326",always_xy=True)
def near(ax,v):
    i=np.clip(np.searchsorted(ax,v),0,len(ax)-1); im=np.maximum(i-1,0)
    return np.where(np.abs(ax[im]-v)<np.abs(ax[i]-v),im,i).astype(np.int32)
CH = 8_000_000
EB = [0.0,0.05,0.1,0.3,0.5,0.7,1.01]
tab = np.zeros((2, len(EB)-1), np.int64)   # [pass/fail][eta bin]
deep_fail_eta = []; deep_fail_dx = []
for s in range(0, NT, CH):
    n = min(CH, NT-s)
    T = pd.read_csv(F, sep=r"\s+", header=None, skiprows=LT+1+s, nrows=n,
                    usecols=range(4), dtype=np.int32, engine="c").to_numpy()-1
    V = P[T]
    e = np.stack([np.linalg.norm(V[:,a]-V[:,b],axis=1) for a,b in PAIRS], 1)
    dx = e.max(1)
    vol = np.abs(np.einsum('ij,ij->i', V[:,1]-V[:,0],
                 np.cross(V[:,2]-V[:,0], V[:,3]-V[:,0])))/6.0
    # Joe-Liu: eta = 12*(3V)^(2/3) / sum(edge^2)
    eta = 12.0*np.cbrt(3.0*vol)**2 / np.maximum((e**2).sum(1), 1e-300)
    bc = V.mean(1)
    lo,la = fwd.transform(bc[:,0], bc[:,1])
    vs = VS[near(mdep,np.maximum(-bc[:,2],0.)), near(mlat,la), near(mlon,lo)].astype(float)
    vs = np.where(np.isfinite(vs), vs, np.inf)
    bad = (vs/dx) < GATE
    ib = np.clip(np.digitize(eta, EB)-1, 0, len(EB)-2)
    for k in range(len(EB)-1):
        tab[0,k] += int(((~bad)&(ib==k)).sum()); tab[1,k] += int((bad&(ib==k)).sum())
    dm = bad & (bc[:,2] < -10000.0)
    if dm.any():
        deep_fail_eta.append(eta[dm].astype(np.float32))
        deep_fail_dx.append(dx[dm].astype(np.float32))
    del T,V,e,dx,vol,eta,bc,vs,bad,ib
    print(f"  ..{s+n:,}", flush=True)
print("\n                    " + "".join(f"{a:>8}-{b:<5}" for a,b in zip(EB[:-1],EB[1:])))
print("  gate PASS        " + "".join(f"{v:>13,}" for v in tab[0]))
print("  gate FAIL        " + "".join(f"{v:>13,}" for v in tab[1]))
tot = tab.sum(0)
print("  fail rate        " + "".join(f"{100*tab[1,k]/max(tot[k],1):>12.2f}%" for k in range(len(tot))))
de = np.concatenate(deep_fail_eta); dd = np.concatenate(deep_fail_dx)
print(f"\n[deep failures, z < -10 km]  n = {len(de):,}")
print(f"  eta   min {de.min():.5f}  median {np.median(de):.5f}  mean {de.mean():.5f}")
print(f"  share with eta < 0.10 : {100*(de<0.10).mean():.1f} %")
print(f"  share with eta < 0.30 : {100*(de<0.30).mean():.1f} %")
print(f"  max edge  median {np.median(dd):,.0f} m   max {dd.max():,.0f} m")
