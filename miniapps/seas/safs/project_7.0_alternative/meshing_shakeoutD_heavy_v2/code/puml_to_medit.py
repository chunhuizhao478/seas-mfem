#!/usr/bin/env python3
"""puml_to_medit.py -- PUML/HDF5 -> MEDIT .mesh directly, for an mmg -optim pass.

Direct because the alternative chain (puml -> gmsh22 ASCII -> meshio read ->
MEDIT) writes and re-reads a multi-GB intermediate for nothing; mmg consumes
MEDIT, so go there once. pandas.to_csv is used for the bulk blocks -- np.savetxt
is ~10x slower at this size and this file is ~87 M tets.

Surface handling: every boundary face is emitted as a Triangle with its SAFS ref
(101 fault / 102 top / 104 absorbing) AND listed in RequiredTriangles. With
`-nosurf` mmg already refuses to move the surface, but Required makes the fault
freeze explicit and survives if -nosurf is ever dropped.

Fault faces are deduped to ONE triangle per face: an interior crack appears on
both incident tets, and emitting it twice makes mmg see a duplicate surface.
"""
import argparse
import sys
import time

import h5py
import numpy as np
import pandas as pd

FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
BC_TO_REF = {3: 101, 1: 102, 5: 104}
CH = 6_000_000

ap = argparse.ArgumentParser()
ap.add_argument("puml")
ap.add_argument("out")
a = ap.parse_args()
t0 = time.time()

f = h5py.File(a.puml)
V = f["geometry"][:]
C = f["connect"]
B = f["boundary"]
NT = C.shape[0]
print(f"[in] {NT:,} tets  {len(V):,} verts", flush=True)

tri, ref = [], []
for s in range(0, NT, CH):
    t = C[s:s + CH][:].astype(np.int64)
    b = B[s:s + CH][:].astype(np.int64)
    for k in range(4):
        code = (b >> (8 * k)) & 0xFF
        for bc, rf in BC_TO_REF.items():
            m = code == bc
            if m.any():
                tri.append(t[m][:, list(FACE[k])])
                ref.append(np.full(int(m.sum()), rf, np.int32))
    del t, b
T = np.vstack(tri)
R = np.concatenate(ref)
del tri, ref
# dedupe: the fault is an interior crack and appears on BOTH incident tets
key = np.sort(T, axis=1)
_, first = np.unique(key, axis=0, return_index=True)
T, R = T[first], R[first]
del key, first
print(f"[surface] {len(T):,} unique triangles  " +
      "  ".join(f"{r}:{int((R == r).sum()):,}" for r in (101, 102, 104)), flush=True)

with open(a.out, "w") as fo:
    fo.write("MeshVersionFormatted 2\nDimension 3\n\nVertices\n%d\n" % len(V))
    pd.DataFrame({0: V[:, 0], 1: V[:, 1], 2: V[:, 2], 3: 0}).to_csv(
        fo, sep=" ", header=False, index=False, float_format="%.10g")
    fo.write("\nTriangles\n%d\n" % len(T))
    pd.DataFrame({0: T[:, 0] + 1, 1: T[:, 1] + 1, 2: T[:, 2] + 1, 3: R}).to_csv(
        fo, sep=" ", header=False, index=False)
    fo.write("\nRequiredTriangles\n%d\n" % len(T))
    pd.DataFrame({0: np.arange(1, len(T) + 1)}).to_csv(
        fo, sep=" ", header=False, index=False)
    fo.write("\nTetrahedra\n%d\n" % NT)
    for s in range(0, NT, CH):
        t = C[s:s + CH][:].astype(np.int64) + 1
        pd.DataFrame({0: t[:, 0], 1: t[:, 1], 2: t[:, 2], 3: t[:, 3],
                      4: np.ones(len(t), np.int32)}).to_csv(
            fo, sep=" ", header=False, index=False)
        del t
        print(f"  ..{min(s+CH, NT):,}/{NT:,}", flush=True)
    fo.write("\nEnd\n")
print(f"[out] {a.out}  ({time.time()-t0:.0f} s)")
