#!/usr/bin/env python3
"""Checks M1-M5 of plan rev5 section 10.2 on a mesh made by build_ucerf3d_mesh.py.

  M1  every fault and creep vertex lies within --m1_tol (1 m) of its source
      surface (the chain slice moved to the vertex depth).
  M2  counts and quality: tet count, minimum dihedral angle >= 10 deg,
      h_min (shortest edge), fault-face edge lengths.
  M3  (ucerf meshes) the smallest distance between vertices of different
      fault surfaces is >= the gap target, and >= 2 local fault elements.
  M4  attribute counts; every fault, creep, top, bottom and wall triangle is
      a face of the tets (2 tets for interior faces, 1 for the outer
      boundary); creep faces meet friction faces only along the friction
      bottom edge.
  M5  sha256 of the .msh equals the one in the builder's JSON.

Usage: check_ucerf3d_mesh.py MESH.msh [--json MESH.json] [--out report.json]
Exit 0 when all checks pass, 1 otherwise.  Standard library + numpy only.
"""
import argparse
import collections
import hashlib
import json
import math
import os
import sys

import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_msh22(path):
    nodes = None
    tris, tets = [], []
    with open(path) as f:
        line = f.readline()
        while line:
            if line.startswith("$MeshFormat"):
                v = f.readline().split()
                if v[0] != "2.2" or v[1] != "0":
                    sys.exit(f"{path}: not gmsh 2.2 ASCII ({v})")
            elif line.startswith("$Nodes"):
                n = int(f.readline())
                a = np.loadtxt(f, max_rows=n)
                ids = a[:, 0].astype(np.int64)
                nodes = np.zeros((ids.max() + 1, 3))
                nodes[ids] = a[:, 1:4]
            elif line.startswith("$Elements"):
                n = int(f.readline())
                for _ in range(n):
                    p = f.readline().split()
                    t, nt = int(p[1]), int(p[2])
                    tags = list(map(int, p[3:3 + nt]))
                    v = list(map(int, p[3 + nt:]))
                    if t == 2:
                        tris.append([tags[0], tags[1]] + v)
                    elif t == 4:
                        tets.append([tags[0]] + v)
            line = f.readline()
    return nodes, np.array(tris, dtype=np.int64), np.array(tets, dtype=np.int64)


def dihedral_min_deg(X):
    """Minimum dihedral angle (deg) of each tet; X (n, 4, 3)."""
    faces = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
    nrm = []
    for i, (a, b, c) in enumerate(faces):
        n = np.cross(X[:, b] - X[:, a], X[:, c] - X[:, a])
        # orient outward: away from the opposite vertex
        s = np.sign(((X[:, i] - X[:, a]) * n).sum(1))
        n = -n * s[:, None]
        nrm.append(n / np.linalg.norm(n, axis=1)[:, None])
    mins = np.full(len(X), 180.0)
    for i in range(4):
        for j in range(i + 1, 4):
            c = np.clip(-(nrm[i] * nrm[j]).sum(1), -1.0, 1.0)
            mins = np.minimum(mins, np.degrees(np.arccos(c)))
    return mins


def seg_dist(P, Q):
    a, b = Q[:-1], Q[1:]
    ab = b - a
    L2 = np.maximum((ab * ab).sum(1), 1e-30)
    d = np.full(len(P), np.inf)
    for i in range(0, len(P), 256):
        p = P[i:i + 256, None, :]
        t = np.clip(((p - a) * ab).sum(2) / L2, 0.0, 1.0)
        r = p - a - t[..., None] * ab
        d[i:i + 256] = np.sqrt((r * r).sum(2)).min(1)
    return d


def surfaces_from_json(rec):
    """(name, top trace (n, 2), horizontal shift per metre of depth)."""
    out = []
    if "chains" in rec:
        for n, c in rec["chains"].items():
            out.append((n, np.array(c["top_trace_m"]), c["dip_deg"], c["dipdir_deg"] or 0.0))
    else:
        f = rec["fault"]
        out.append(("fault", np.array(f["top_trace_m"]), f["dip_deg"], 0.0))
    res = []
    for n, tr, dip, dd in out:
        if dip < 90.0:
            u = np.array([math.sin(math.radians(dd)), math.cos(math.radians(dd))])
            rate = u / math.tan(math.radians(dip))
        else:
            rate = np.zeros(2)
        res.append((n, tr, rate))
    return res


def close_pairs(P, lab, r):
    """Smallest distance between points of P with different labels, among
    pairs closer than r (grid hashing); inf if none."""
    key = np.floor(P / r).astype(np.int64)
    cells = collections.defaultdict(list)
    for i, k in enumerate(map(tuple, key)):
        cells[k].append(i)
    best, where = np.inf, None
    offs = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for k, idx in cells.items():
        idx = np.array(idx)
        nb = []
        for o in offs:
            q = cells.get((k[0] + o[0], k[1] + o[1], k[2] + o[2]))
            if q:
                nb.extend(q)
        nb = np.array(nb)
        for i in idx:
            m = lab[nb] != lab[i]
            if m.any():
                d = np.linalg.norm(P[nb[m]] - P[i], axis=1)
                j = d.argmin()
                if d[j] < best:
                    best, where = d[j], (int(i), int(nb[m][j]))
    return best, where


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("msh")
    ap.add_argument("--json")
    ap.add_argument("--out")
    ap.add_argument("--m1_tol", type=float, default=1.0)
    ap.add_argument("--min_dihedral", type=float, default=10.0)
    args = ap.parse_args()
    js = args.json or os.path.splitext(args.msh)[0] + ".json"
    rec = json.load(open(js))
    A = rec["attributes"]
    X, tris, tets = read_msh22(args.msh)
    rep, ok = {"mesh": os.path.basename(args.msh)}, True

    def check(name, passed, **info):
        nonlocal ok
        ok &= bool(passed)
        rep[name] = dict(passed=bool(passed), **info)
        print(f"{'PASS' if passed else 'FAIL'}  {name}: " + ", ".join(f"{k} {v}" for k, v in info.items()))

    # ---- M5
    check("M5 sha256", sha256(args.msh) == rec["msh_sha256"], sha256=rec["msh_sha256"][:16] + "...")

    # ---- M4 attribute counts and conformity
    attr = tris[:, 0]
    counts = {int(a): int((attr == a).sum()) for a in np.unique(attr)}
    want = {int(k.split(":")[1]): v for k, v in rec["counts"].items() if k.startswith("2:")}
    check("M4 counts", counts == want, triangles=counts, tets=len(tets))
    T = np.sort(tets[:, 1:5], axis=1)
    fk = collections.Counter()
    for a, b, c in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        for f in map(tuple, T[:, [a, b, c]]):
            fk[f] += 1
    inner = {A["friction"], A["creep"]}
    bad = collections.Counter()
    for t in tris:
        need = 2 if t[0] in inner else 1
        if fk.get(tuple(sorted(t[2:5])), 0) != need:
            bad[int(t[0])] += 1
    check("M4 faces of the tets", not bad, wrong_by_attribute=dict(bad))
    fric = tris[attr == A["friction"]][:, 2:5]
    creep = tris[attr == A["creep"]][:, 2:5]
    if len(creep):
        shared = np.intersect1d(np.unique(fric), np.unique(creep))
        zb = -rec["depths_km"]["friction"] * 1e3 if "depths_km" in rec else -rec["fault"]["depth_km"] * 1e3
        dz = np.abs(X[shared, 2] - zb).max() if len(shared) else 0.0
        check("M4 creep meets friction at its bottom edge only", len(shared) > 0 and dz < 1.0,
              shared_vertices=len(shared), max_dz_m=round(float(dz), 6), bottom_z_m=zb)

    # ---- M1 vertices on their source surface
    surf = surfaces_from_json(rec)
    fv = np.unique(np.concatenate([fric.ravel(), creep.ravel()]) if len(creep) else fric.ravel())
    P = X[fv]
    d = np.full(len(P), np.inf)
    lab = np.zeros(len(P), dtype=np.int64)
    for k, (n, tr, rate) in enumerate(surf):
        q = P[:, :2] - rate[None, :] * (-P[:, 2:3])
        dk = seg_dist(q, tr)
        lab[dk < d] = k
        d = np.minimum(d, dk)
    # vertical distance of the top vertices: they must be at z = 0 or below
    check("M1 vertices on source surfaces", d.max() <= args.m1_tol and P[:, 2].max() <= 1e-6,
          vertices=len(P), max_dist_m=round(float(d.max()), 4), tol_m=args.m1_tol)

    # ---- M2 quality
    Xt = X[tets[:, 1:5]]
    dmin = dihedral_min_deg(Xt)
    E = np.concatenate([np.linalg.norm(Xt[:, i] - Xt[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)])
    Xf = X[fric]
    ef = np.concatenate([np.linalg.norm(Xf[:, i] - Xf[:, (i + 1) % 3], axis=1) for i in range(3)])
    check("M2 quality", dmin.min() >= args.min_dihedral, tets=len(tets),
          min_dihedral_deg=round(float(dmin.min()), 2), n_below_15deg=int((dmin < 15).sum()),
          h_min_m=round(float(E.min()), 1), h_max_m=round(float(E.max()), 1),
          fault_edge_m={"min": round(float(ef.min()), 1), "median": round(float(np.median(ef)), 1),
                        "max": round(float(ef.max()), 1)})
    rep["M2 quality"]["dihedral_hist_deg"] = np.histogram(dmin, bins=[0, 5, 10, 15, 20, 30, 45, 90])[0].tolist()

    # ---- M3 gaps between surfaces
    if "gap_target_km" in rec:
        g = rec["gap_target_km"] * 1e3
        best, where = close_pairs(P, lab, 2.0 * g)
        hloc = None
        if where is not None:
            # local fault edge length: median edge of fault/creep faces within 1.5 gaps
            allf = np.concatenate([fric, creep]) if len(creep) else fric
            cen = X[allf].mean(1)
            mid = 0.5 * (P[where[0]] + P[where[1]])
            near = np.linalg.norm(cen - mid, axis=1) < max(1.5 * best, 1e3)
            Xn = X[allf[near]]
            hloc = float(np.median(np.concatenate([np.linalg.norm(Xn[:, i] - Xn[:, (i + 1) % 3], axis=1)
                                                    for i in range(3)])))
        elems = best / hloc if hloc else np.inf
        info = dict(gap_target_m=g, min_vertex_gap_m=round(float(best), 1) if np.isfinite(best) else "> 2 gap",
                    local_fault_edge_m=round(hloc, 1) if hloc else None,
                    elements_across=round(float(elems), 2) if np.isfinite(elems) else None)
        if where is not None:
            info["between"] = [surf[lab[where[0]]][0], surf[lab[where[1]]][0]]
            info["at_depth_km"] = round(float(-0.5 * (P[where[0], 2] + P[where[1], 2])) / 1e3, 2)
        check("M3 surface gaps", (not np.isfinite(best)) or (best >= g - 1.0 and elems >= 2.0), **info)
        rep["M3 surface gaps"]["geometry_min_km"] = {k: v["min_km"] for k, v in rec["gaps"].items()}

    rep["all_passed"] = bool(ok)
    if args.out:
        json.dump(rep, open(args.out, "w"), indent=1)
    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
