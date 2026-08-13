#!/usr/bin/env python3
"""extract_fault_surface.py -- the fault triangulation to be REUSED verbatim.

Pulls the dynamic-rupture surface out of the small-domain heavy mesh as a
standalone triangulated surface (vertices + triangles, re-indexed), plus the
boundary loop of that surface -- its TOP edge is the trace where the fault
daylights, and the PLC's top surface must contain it exactly or the fault will
not reach the free surface.
"""
import sys
import numpy as np, h5py
FACE=[(0,2,1),(0,1,3),(1,2,3),(0,3,2)]
# ZTOP -- which fault-boundary vertices count as TRACE (snapped to z=0) rather
# than genuine buried edge.  MUST equal trace_polylines.py's ZTOP exactly.
#
# Measured on the deployed ALT heavy mesh: 7,785 fault-boundary vertices already
# daylight (they ARE free-surface vertices, at exactly z=0), and a 390-vertex
# fringe sits buried at -5.03 .. -142.46 m.  The fringe histogram is
#   -10..-5: 3   -25..-10: 20   -50..-25: 266   -60..-50: 63   -100..-60: 36   -200..-100: 2
# so -50 cuts THROUGH the -60..-50 group and severs the trace in 3 places
# (measured: 7 chains, 8 phantom ends, ~4.1 km unembedded).  -60 takes the whole
# -60..-50 group -> 352 snapped, 4 chain ends = 2 lateral tips + 2 T-junctions.
# Do not go past ~-65: at -70 a vertex on a descending lateral edge (neighbours
# at 0 and -134.98 m) folds a fault triangle into the z=0 plane.
ZTOP=-60.0
src, out = sys.argv[1], sys.argv[2]
with h5py.File(src) as f:
    V=f["geometry"][:].astype(np.float64); C=f["connect"]; B=f["boundary"][:].astype(np.int64)
    tri=[]
    for s in range(0,C.shape[0],8_000_000):
        c=C[s:s+8_000_000][:].astype(np.int64); b=B[s:s+8_000_000]
        for k in range(4):
            m=((b>>(8*k))&0xFF)==3
            if m.any(): tri.append(c[m][:,list(FACE[k])])
F=np.vstack(tri)
# each DR face appears twice (once per side); keep one copy, consistently oriented
key=np.sort(F,axis=1)
o=np.lexsort((key[:,2],key[:,1],key[:,0])); key=key[o]; F=F[o]
assert np.all(key[0::2]==key[1::2]), "a DR face is not shared by exactly two tets"
T=F[0::2]
vid=np.unique(T); remap=-np.ones(len(V),np.int64); remap[vid]=np.arange(len(vid))
P=V[vid]; T2=remap[T]
print(f"[fault] {len(T2):,} triangles  {len(P):,} vertices")
print(f"  x {P[:,0].min():.1f}..{P[:,0].max():.1f}  y {P[:,1].min():.1f}..{P[:,1].max():.1f}  z {P[:,2].min():.1f}..{P[:,2].max():.1f}")
# boundary loop = edges used exactly once
e=np.concatenate([T2[:,[0,1]],T2[:,[1,2]],T2[:,[0,2]]])
e=np.sort(e,axis=1)
u,c=np.unique(e,axis=0,return_counts=True)
bnd=u[c==1]
print(f"[boundary] {len(bnd):,} edges used once  ({len(np.unique(bnd)):,} vertices)")
zb=P[np.unique(bnd),2]
print(f"  boundary z {zb.min():.2f} .. {zb.max():.2f}")
top=np.unique(bnd)[zb>ZTOP]
print(f"  TOP edge (z>{ZTOP:.0f}): {len(top):,} vertices, z {P[top,2].min():.2f}..{P[top,2].max():.2f}")
# FLAT-TOP DEVIATION, deliberate and documented.  The daylighting edge is at z=0.00
# at the 1st/50th/99th percentile, but 46 nodes excurse to -49.9..+25.8 m.  A free
# surface forced through those excursions cannot be triangulated without the fault
# poking through it: tetgen -d gave 11 self-intersections, and every attempt to
# refine the surface to follow them made it WORSE (16, 21, 63) -- finer triangles
# track the excursion more faithfully and intersect more.  Snapping the edge to z=0
# moves 46 of 1,386,066 fault vertices (0.003 %) by <=50 m, all within 50 m of the
# surface, and yields an exactly FLAT top -- which the MUSCAL-native deck stack
# requires anyway (it maps depth = -z and is documented as flat-top only).
nmoved=int((np.abs(P[top,2])>1e-9).sum()); dmax=float(np.abs(P[top,2]).max())
P[top,2]=0.0
print(f"  FLATTENED: {nmoved} of {len(top):,} trace vertices snapped to z=0 (max move {dmax:.2f} m)")
# ALSO clamp any INTERIOR fault vertex above the flat top.  Exactly 3 of 1,386,066
# sit 0.83-2.49 m above z=0, in one cluster -- and they are the whole reason tetgen
# reported 11 self-intersections: a flat free surface cannot avoid a fault that
# pokes through it.  Found by reading tetgen's own skipped-face list; three earlier
# guesses (taper, local refinement, edge flattening) all made it WORSE.
na=int((P[:,2]>1e-9).sum()); da=float(P[:,2].max())
P[:,2]=np.minimum(P[:,2],0.0)
print(f"  CLAMPED: {na} interior fault vertices above z=0 pulled down (max move {da:.2f} m)")
# triangle quality of the surface
a=P[T2[:,1]]-P[T2[:,0]]; b=P[T2[:,2]]-P[T2[:,0]]
ar=0.5*np.linalg.norm(np.cross(a,b),axis=1)
el=np.stack([np.linalg.norm(P[T2[:,i]]-P[T2[:,j]],axis=1) for i,j in ((0,1),(1,2),(0,2))],1)
print(f"[quality] area {ar.sum()/1e6:,.2f} km2   edge min {el.min():.2f} med {np.median(el):.2f} max {el.max():.2f} m")
np.savez_compressed(out, P=P, T=T2, bnd=bnd, top=top)
print(f"[out] {out}")
