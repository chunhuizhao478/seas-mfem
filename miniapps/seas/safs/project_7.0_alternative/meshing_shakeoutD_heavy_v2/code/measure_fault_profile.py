#!/usr/bin/env python3
"""measure_fault_profile.py -- the fault-zone size treatment to be PRESERVED.

Element max-edge as a function of distance to the dynamic-rupture surface, on the
existing small-domain heavy mesh.  This profile IS the "fault zone treatment": the
new large-domain mesh must reproduce it, and everything beyond it is far field
sized by the resolved-frequency gate.
"""
import sys
import numpy as np, h5py
from scipy.spatial import cKDTree
PAIRS=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
FACE=[(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
CH=4_000_000
mesh=sys.argv[1]
with h5py.File(mesh) as f:
    V=f["geometry"][:].astype(np.float64); nt=f["connect"].shape[0]
    C=f["connect"]; B=f["boundary"][:].astype(np.int64)
    tri=[]
    for s in range(0,nt,CH):
        c=C[s:s+CH][:].astype(np.int64); b=B[s:s+CH]
        for k in range(4):
            m=((b>>(8*k))&0xFF)==3
            if m.any(): tri.append(c[m][:,list(FACE[k])])
    F=np.vstack(tri)
    Pc=V[F].mean(1)
    tree=cKDTree(Pc)
    print(f"[fault] {len(F)//2:,} facets; KD on {len(Pc):,} facet centroids")
    EDGES=np.array([0,125,250,500,1000,2000,4000,7500,15000,30000,60000,1e12])
    nb=np.zeros(len(EDGES)-1,np.int64); acc=[[] for _ in range(len(EDGES)-1)]
    volb=np.zeros(len(EDGES)-1)
    for s in range(0,nt,CH):
        p=V[C[s:s+CH][:].astype(np.int64)]
        bc=p.mean(1)
        dx=np.max(np.stack([np.linalg.norm(p[:,b]-p[:,a],axis=1) for a,b in PAIRS],1),1)
        d=p[:,1:]-p[:,0:1]
        vv=np.abs(np.einsum('ij,ij->i',np.cross(d[:,0],d[:,1]),d[:,2]))/6.0
        dist,_=tree.query(bc,k=1,distance_upper_bound=60000.0,workers=-1)
        idx=np.clip(np.searchsorted(EDGES,dist,side="right")-1,0,len(EDGES)-2)
        for i in range(len(EDGES)-1):
            m=idx==i
            if m.any():
                nb[i]+=int(m.sum()); volb[i]+=float(vv[m].sum())
                acc[i].append(dx[m][::max(1,m.sum()//20000)].astype(np.float32))
        del p,bc,dx,d,vv,dist,idx
print(f"\n{'band (m)':>16} | {'cells':>12} | {'share':>7} | {'vol km3':>10} | dx med  p90   max")
for i in range(len(EDGES)-1):
    if not nb[i]: continue
    a=np.concatenate(acc[i])
    hi = "inf" if EDGES[i+1]>1e11 else f"{EDGES[i+1]:.0f}"
    print(f"{EDGES[i]:8.0f}..{hi:>6} | {nb[i]:12,} | {100*nb[i]/nt:6.2f}% | {volb[i]/1e9:10,.1f} | "
          f"{np.median(a):6.1f} {np.percentile(a,90):6.1f} {a.max():7.1f}")
print(f"\ntotal {nt:,} cells, {volb.sum()/1e9:,.1f} km3")
np.save("build_tmp/fault_profile.npy", np.array([EDGES[:-1], nb, volb], dtype=object), allow_pickle=True)
