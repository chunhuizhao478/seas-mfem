#!/usr/bin/env python3
"""identify_fault_source.py -- which conditioned strand set IS the deployed fault?

Stage A cannot start until the source geometry is pinned down. The deployed
heavy mesh's fault is 13073.313420396 km2 over 3 tags (101/102/103), and the
build tree holds seven candidate conditioned strand sets (faults_corefined, _cr,
_fine, _fine2, _fine3, _knot, _off) plus the raw CFM STLs. Guessing wrong here
would mean rebuilding a DIFFERENT fault, which no downstream check would catch
as a geometry error -- the gates and quality checks would all pass on the wrong
surface.

So: compute area and bbox per tag for every candidate and match against the
deployed mesh, which is measured from the mesh itself, not quoted.
"""
import argparse
from pathlib import Path

import h5py
import numpy as np

LOCAL_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))
BC_DR = 3
CH = 1_000_000


def face_code(b, s):
    return ((np.ascontiguousarray(b, np.int32).view(np.uint32) >> np.uint32(8 * s))
            & np.uint32(0xFF)).astype(np.int32)


def tri_area(V):
    return np.linalg.norm(np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]), axis=1) * 0.5


def read_off(p):
    """OFF reader that tokenises the whole file.

    Reading line-by-line with a fixed vertex count breaks on these files: face
    lines carry a leading valence and some writers wrap or pad, so a positional
    reader silently mis-slices. Tokenising is O(file) and cannot drift.
    """
    tok = open(p).read().split()
    assert tok[0] == "OFF", f"{p} is not an OFF file"
    nv, nf = int(tok[1]), int(tok[2])
    V = np.array(tok[4:4 + 3 * nv], dtype=np.float64).reshape(nv, 3)
    F = np.empty((nf, 3), np.int64)
    i = 4 + 3 * nv
    for k in range(nf):
        m = int(tok[i])
        F[k] = [int(tok[i + 1]), int(tok[i + 2]), int(tok[i + 3])]
        i += m + 1
    return V, F


def read_stl(p):
    """ASCII or binary STL -> (n,3,3)."""
    raw = open(p, "rb").read(512)
    if raw[:5] == b"solid" and b"facet" in raw:
        V = []
        for line in open(p):
            s = line.split()
            if s and s[0] == "vertex":
                V.append([float(s[1]), float(s[2]), float(s[3])])
        return np.array(V).reshape(-1, 3, 3)
    d = np.fromfile(p, dtype=np.uint8)
    n = int(np.frombuffer(d[80:84].tobytes(), np.uint32)[0])
    rec = d[84:84 + n * 50].reshape(n, 50)
    f = np.frombuffer(rec[:, 12:48].tobytes(), np.float32).reshape(n, 3, 3)
    return f.astype(np.float64)


def clipped_area_below(T, z0):
    """Total area of the parts of triangles T lying below the plane z = z0.

    Straddling triangles are clipped exactly (Sutherland-Hodgman against one
    half-space, then fan-triangulated), because dropping or keeping them whole
    would bias the comparison by far more than the differences being resolved.
    """
    tot = 0.0
    for tri in T:
        below = tri[:, 2] <= z0
        k = int(below.sum())
        if k == 3:
            tot += float(np.linalg.norm(np.cross(tri[1] - tri[0],
                                                 tri[2] - tri[0])) * 0.5)
            continue
        if k == 0:
            continue
        poly = []
        for i in range(3):
            a, b = tri[i], tri[(i + 1) % 3]
            if a[2] <= z0:
                poly.append(a)
            if (a[2] - z0) * (b[2] - z0) < 0:
                t = (z0 - a[2]) / (b[2] - a[2])
                poly.append(a + t * (b - a))
        if len(poly) >= 3:
            P = np.array(poly)
            for i in range(1, len(P) - 1):
                tot += float(np.linalg.norm(np.cross(P[i] - P[0],
                                                     P[i + 1] - P[0])) * 0.5)
    return tot


def deployed(mesh):
    with h5py.File(mesh, "r") as f:
        G = f["geometry"][:]
        B = f["boundary"][:].astype(np.int32)
        conn = f["connect"]
        nt = conn.shape[0]
        A = 0.0
        lo = np.full(3, np.inf)
        hi = np.full(3, -np.inf)
        n = 0
        for s0 in range(0, nt, CH):
            s1 = min(s0 + CH, nt)
            Bc = B[s0:s1]
            hit = np.zeros(s1 - s0, bool)
            for s in range(4):
                hit |= face_code(Bc, s) == BC_DR
            if not hit.any():
                continue
            idx = np.nonzero(hit)[0]
            Tc = conn[s0:s1][idx].astype(np.int64)
            Bs = Bc[idx]
            for s in range(4):
                m = face_code(Bs, s) == BC_DR
                if m.any():
                    V = G[Tc[m][:, list(LOCAL_FACES[s])]]
                    A += tri_area(V).sum()
                    n += len(V)
                    lo = np.minimum(lo, V.reshape(-1, 3).min(0))
                    hi = np.maximum(hi, V.reshape(-1, 3).max(0))
    return A / 2.0, n // 2, lo, hi     # each fault face is stored twice


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--root", required=True)
    a = ap.parse_args()

    A, n, lo, hi = deployed(a.mesh)
    print(f"DEPLOYED fault: {A/1e6:,.6f} km2, {n:,} unique triangles")
    print(f"  bbox E {lo[0]:,.0f}..{hi[0]:,.0f}  N {lo[1]:,.0f}..{hi[1]:,.0f}  "
          f"z {lo[2]:,.0f}..{hi[2]:,.0f}\n")

    root = Path(a.root)
    cands = sorted([d for d in root.iterdir() if d.is_dir()
                    and list(d.glob("fault_*.off"))])
    print(f"{'candidate':<22}{'tris':>10}{'area km2':>15}{'z_min':>10}{'z_max':>8}"
          f"{'d(area)':>12}")
    for d in cands:
        tot = 0.0
        ntri = 0
        zlo, zhi = np.inf, -np.inf
        for p in sorted(d.glob("fault_*.off")):
            V, F = read_off(p)
            tot += tri_area(V[F]).sum()
            ntri += len(F)
            zlo = min(zlo, V[:, 2].min())
            zhi = max(zhi, V[:, 2].max())
        print(f"{d.name:<22}{ntri:>10,}{tot/1e6:>15,.6f}{zlo:>10,.0f}{zhi:>8,.0f}"
              f"{(tot-A)/1e6:>12,.3f}")

    # THE DECISIVE TEST: every candidate reaches z = +2,130 while the deployed
    # fault stops at z = 0, i.e. the build clips the strands at the flat lid to
    # make them daylight. Comparing RAW areas therefore cannot separate the
    # candidates (they sit within 1.2 km2 of each other); comparing the area
    # BELOW z = 0, exactly clipped, can.
    print(f"\n  area strictly BELOW z=0 (exact triangle-plane clip)")
    print(f"    {'candidate':<22}{'area km2':>15}{'delta vs deployed':>20}")
    for d in cands:
        tot = 0.0
        for p in sorted(d.glob("fault_*.off")):
            V, F = read_off(p)
            tot += clipped_area_below(V[F], 0.0)
        print(f"    {d.name:<22}{tot/1e6:>15,.6f}{(tot-A)/1e6:>20,.6f}")

    for p in sorted(root.glob("fault_*.stl")):
        V = read_stl(p)
        print(f"{p.name:<22}{len(V):>10,}{tri_area(V).sum()/1e6:>15,.6f}"
              f"{V[:,:,2].min():>10,.0f}{V[:,:,2].max():>8,.0f}")


if __name__ == "__main__":
    main()
