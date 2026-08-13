#!/usr/bin/env python3
"""extract_fault_trace.py -- the fault's SURFACE TRACE, cached for plotting.

Takes every edge of the dynamic-rupture triangulation whose BOTH endpoints sit at
the free surface, so branches and separate strands render as they actually are.
A fitted polyline would misrepresent the SAFS fault system, which is multi-strand;
and a footprint band is not a trace (that mistake is on record in this project).
"""
import sys
import numpy as np, h5py
FACE=[(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
mesh, out = sys.argv[1], sys.argv[2]
ZTOP = float(sys.argv[3]) if len(sys.argv)>3 else -50.0
with h5py.File(mesh) as f:
    V=f["geometry"][:].astype(np.float64); C=f["connect"]; B=f["boundary"][:].astype(np.int64)
    tri=[]
    for s in range(0,C.shape[0],8_000_000):
        c=C[s:s+8_000_000][:].astype(np.int64); b=B[s:s+8_000_000]
        for k in range(4):
            m=((b>>(8*k))&0xFF)==3
            if m.any(): tri.append(c[m][:,list(FACE[k])])
F=np.vstack(tri)
print(f"[fault] {len(F)//2:,} facets, z {V[np.unique(F),2].min():.1f}..{V[np.unique(F),2].max():.1f}")
surf = V[:,2] > ZTOP
e = np.concatenate([F[:,[0,1]], F[:,[1,2]], F[:,[0,2]]])
e = e[surf[e[:,0]] & surf[e[:,1]]]
e = np.unique(np.sort(e,axis=1), axis=0)
print(f"[trace] {len(e):,} surface edges (z > {ZTOP:.0f} m), "
      f"{len(np.unique(e)):,} vertices")
seg = V[e][:, :, :2]
L = np.linalg.norm(seg[:,1]-seg[:,0], axis=1).sum()
print(f"[trace] total edge length {L/1000:,.1f} km")
np.savez_compressed(out, seg=seg)
print(f"[out] {out}")
