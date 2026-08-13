#!/usr/bin/env python3
"""calibrate_density.py -- how many tets does a size field h(x) actually buy?

A regular tet of edge h has volume h^3/(6*sqrt2) = h^3/8.485, so a naive count is
N = integral dV / (h^3/8.485).  Real graded meshes do not hit that.  This measures
the CALIBRATION CONSTANT k in N = k * integral(dV / h^3) on a mesh we already have,
so the new domain's estimate is anchored on this pipeline's actual behaviour rather
than on an idealisation.

For each cell of the EXISTING mesh we know its volume and its max edge, so
    integral dV / h^3  ~=  sum_cells V_cell / dx_cell^3
and k = N_cells / that sum.
"""
import sys
import numpy as np, h5py
PAIRS=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
CH=4_000_000
mesh=sys.argv[1]
with h5py.File(mesh) as f:
    V=f["geometry"][:].astype(np.float64); nt=f["connect"].shape[0]; C=f["connect"]
    tot=0.0; vol=0.0
    for s in range(0,nt,CH):
        p=V[C[s:s+CH][:].astype(np.int64)]
        d=p[:,1:]-p[:,0:1]
        vv=np.abs(np.einsum('ij,ij->i',np.cross(d[:,0],d[:,1]),d[:,2]))/6.0
        dx=np.max(np.stack([np.linalg.norm(p[:,b]-p[:,a],axis=1) for a,b in PAIRS],1),1)
        tot+=float((vv/dx**3).sum()); vol+=float(vv.sum())
        del p,d,vv,dx
print(f"mesh {mesh.split('/')[-1]}")
print(f"  cells {nt:,}   volume {vol/1e9:,.1f} km3")
print(f"  sum V/dx^3 = {tot:,.1f}   ->  k = N / that = {nt/tot:.4f}")
print(f"  (ideal regular-tet k would be 8.485)")
