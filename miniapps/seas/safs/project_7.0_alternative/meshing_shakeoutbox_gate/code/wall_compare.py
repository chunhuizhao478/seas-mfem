#!/usr/bin/env python3
"""wall_compare.py -- is the vertical absorbing wall the same on two parents?

A frozen-parent collar is built against ONE wall and welded by exact coordinate
key.  The intermediate parent and `refine2` were measured to share a
bit-identical 26,704-triangle wall, which is why one collar served both.  The
1 Hz @ p5 parent was produced from `refine2` by conforming bisection, and its
report shows absorbing facets going 57,790 -> 77,127 -- so the wall almost
certainly moved and a second collar is needed.  Confirm that rather than assume
it, and report the exact triangle/vertex counts a new collar must match.

Streaming: only wall faces are retained (tens of thousands), never a full
connect array.
"""
import argparse
import hashlib
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "meshing_shakeoutbox_intermediate" / "code"))
from collar_lib import BC_ABSORBING, LOCAL_FACES, face_code  # noqa: E402

CH = 2_000_000


def wall_tris(path, vert_tol=0.5):
    """Vertical absorbing faces, as (n,3,3) coordinates."""
    out = []
    nabs = 0
    with h5py.File(path, "r") as f:
        G = f["geometry"][:]
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            Bc = B[s0:s1]
            hit = np.zeros(s1 - s0, bool)
            for s in range(4):
                hit |= face_code(Bc, s) == BC_ABSORBING
            if not hit.any():
                continue
            idx = np.nonzero(hit)[0]
            Tc = conn[s0:s1][idx].astype(np.int64)
            Bs = Bc[idx]
            for s in range(4):
                m = face_code(Bs, s) == BC_ABSORBING
                if not m.any():
                    continue
                tri = Tc[m][:, list(LOCAL_FACES[s])]
                nabs += len(tri)
                V = G[tri]
                n = np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0])
                ln = np.linalg.norm(n, axis=1, keepdims=True)
                n = n / np.where(ln == 0.0, 1.0, ln)
                out.append(V[np.abs(n[:, 2]) < vert_tol])
    W = np.vstack(out) if out else np.zeros((0, 3, 3))
    return W, nabs


def canon(W):
    K = np.round(W * 1e3).astype(np.int64)
    i = np.lexsort((K[:, :, 2], K[:, :, 1], K[:, :, 0]), axis=1)
    K = np.take_along_axis(K, i[:, :, None], axis=1).reshape(len(K), 9)
    return K[np.lexsort([K[:, k] for k in range(8, -1, -1)])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", action="append", required=True)
    a = ap.parse_args()

    info = []
    for p in a.mesh:
        W, nabs = wall_tris(p)
        C = canon(W)
        h = hashlib.sha1(C.tobytes()).hexdigest()[:16]
        nv = len(np.unique(C.reshape(-1, 3), axis=0))
        z = W[:, :, 2]
        print(f"{Path(p).name}")
        print(f"    absorbing facets total {nabs:,}   VERTICAL (the wall) {len(W):,}")
        print(f"    wall vertices {nv:,}   z {z.min():,.1f} .. {z.max():,.1f}")
        print(f"    canonical sha1 {h}", flush=True)
        info.append((Path(p).name, len(W), nv, h))

    print("\n  pairwise wall identity:")
    for i in range(len(info)):
        for j in range(i + 1, len(info)):
            same = info[i][3] == info[j][3]
            print(f"    {info[i][0][:44]:<46} vs {info[j][0][:44]:<46} "
                  f"{'IDENTICAL' if same else 'DIFFERENT'}")


if __name__ == "__main__":
    main()
