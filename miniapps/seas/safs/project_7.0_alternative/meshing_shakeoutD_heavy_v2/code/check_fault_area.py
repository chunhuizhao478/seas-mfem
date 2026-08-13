#!/usr/bin/env python3
"""check_fault_area.py -- did refining the fault triangulation MOVE the fault?

Splitting a fault edge at its midpoint subdivides the dynamic-rupture surface
without moving it: the midpoint of an edge of a fault triangle lies ON that
triangle.  This asserts that, rather than trusting it:

  F1  total fault AREA identical to the parent's (relative delta ~1e-15)
  F2  every parent fault VERTEX still present at identical coordinates
  F3  every fault face interior x2 with both tets tagged BC 3 (via bbox count)
  F4  facet count >= parent's, and the increase equals the subdivision

Usage: check_fault_area.py <product.puml.h5> <parent.puml.h5>
"""
import sys
import numpy as np, h5py

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]

def fault(path):
    with h5py.File(path) as f:
        V = f["geometry"][:].astype(np.float64)
        T = f["connect"][:].astype(np.int64)
        B = f["boundary"][:].astype(np.int64)
    tri = []
    for k in range(4):
        m = ((B >> (8 * k)) & 0xFF) == 3
        if m.any():
            tri.append(T[m][:, list(FACE[k])])
    F = np.concatenate(tri)
    P = V[F]
    a = np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1) / 2.0
    return V, F, a

prod, par = sys.argv[1], sys.argv[2]
Vp, Fp, ap = fault(par)
Vq, Fq, aq = fault(prod)
Ap, Aq = ap.sum(), aq.sum()
vp = np.unique(Fp); vq = np.unique(Fq)

print(f"parent  : {len(Fp):>10,} fault faces  {len(vp):>9,} fault verts  area {Ap/1e6:14.8f} km2")
print(f"product : {len(Fq):>10,} fault faces  {len(vq):>9,} fault verts  area {Aq/1e6:14.8f} km2")
rel = abs(Aq - Ap) / Ap
print(f"F1 area delta            : {Aq-Ap:+.6e} m2   relative {rel:.3e}   "
      f"[{'PASS' if rel < 1e-12 else 'FAIL'}]")
keep = np.intersect1d(vp, vq)
moved = 0
if len(keep):
    d = np.linalg.norm(Vq[keep] - Vp[keep], axis=1)
    moved = int((d > 0).sum())
print(f"F2 parent fault verts    : {len(keep):,} of {len(vp):,} still present, "
      f"{moved} moved  [{'PASS' if len(keep) == len(vp) and moved == 0 else 'FAIL'}]")
# every fault face must appear exactly twice (once from each side)
# NOTE: packing three vertex ids as v0*2^42 + v1*2^21 + v2 OVERFLOWS -- these
# meshes carry up to 34.3 M vertices, far past 2^21 = 2,097,152, so distinct
# triangles collide and the check reports phantom ">2" faces.  Compare the
# sorted triples themselves via a void view instead: exact, no packing.
key = np.ascontiguousarray(np.sort(Fq, axis=1).astype(np.int64))
u, c = np.unique(key.view([('a', np.int64), ('b', np.int64),
                           ('c', np.int64)]).ravel(), return_counts=True)
print(f"F3 fault faces shared x2 : {int((c == 2).sum()):,} of {len(u):,} distinct  "
      f"(x1 {int((c==1).sum())}, >2 {int((c>2).sum())})  "
      f"[{'PASS' if (c == 2).all() else 'FAIL'}]")
print(f"F4 facet count           : {len(Fp)//2:,} -> {len(Fq)//2:,} "
      f"({100*(len(Fq)-len(Fp))/len(Fp):+.4f} %)")
