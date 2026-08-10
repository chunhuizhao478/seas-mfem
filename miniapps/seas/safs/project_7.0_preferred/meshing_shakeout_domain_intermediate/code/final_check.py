#!/usr/bin/env python3
"""final_check.py -- post-modification state of a PUML mesh: z range, lid
flatness, inverted tets, minimum edge.  Streams, so it is safe on 10^8 cells."""
import sys
import numpy as np
import h5py

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
CH = 2_000_000


def fc(b, s):
    return ((np.ascontiguousarray(b, np.int32).view(np.uint32) >> np.uint32(8 * s))
            & np.uint32(0xFF)).astype(np.int32)


for path in sys.argv[1:]:
    with h5py.File(path, "r") as f:
        P = f["geometry"][:]
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        lid, neg, emin = set(), 0, np.inf
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            T = conn[s0:s1].astype(np.int64)
            Bc = B[s0:s1]
            p = P[T]
            v = np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                          np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])) / 6.0
            neg += int((v < 0).sum())
            for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
                emin = min(emin, float(np.linalg.norm(p[:, a] - p[:, b], axis=1).min()))
            for s in range(4):
                m = fc(Bc, s) == 1
                if m.any():
                    lid.update(np.unique(T[m][:, list(LOCAL_FACES[s])]).tolist())
            del T, Bc, p, v
    lid = np.array(sorted(lid))
    zl = P[lid][:, 2]
    print(f"\n=== {path.split('/')[-1]}   {nt:,} tets")
    print(f"  z range    {P[:,2].min():,.4f} .. {P[:,2].max():,.4f}")
    print(f"  lid strays {int((np.abs(zl) > 1e-6).sum()):,} of {len(lid):,}"
          f"   (z {zl.min():.4f} .. {zl.max():.4f})")
    print(f"  inverted   {neg}")
    print(f"  min edge   {emin:,.3f} m")
