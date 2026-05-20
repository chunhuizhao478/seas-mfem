#!/usr/bin/env python3
"""locate_fault_slivers.py — pinpoint low-quality fault triangles near a point.

Complements check_mesh_quality.py: that tool reports aggregate Joe-Liu eta /
triangle-q gates + histograms; this one is a pure-stdlib (no meshio/numpy)
*spatial* probe that answers "how good are the fault elements AROUND this
coordinate, vs the fault-wide distribution?".  Written for the 2026-05-20
post-event-blowup debug (spatial_dynamic_rupture_postevent_debug_2026-05-20.md),
where the spurious north-tip velocity coincided with sliver triangles at the
shallow fault-surface trace.

Gmsh v2.2 ASCII only.  Fault = physical surface tag (default 101).

Usage:
    python3 locate_fault_slivers.py MESH.msh \
        --target 448338 3798440 -249.7 --radius 3000 [--fault-tag 101]
"""
import argparse
import math


def read_fault_tris(path, fault_tag):
    nodes = {}
    tris = []
    with open(path) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln == "$Nodes":
            n = int(lines[i + 1]); i += 2
            for k in range(n):
                p = lines[i + k].split()
                nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            i += n; continue
        if ln == "$Elements":
            n = int(lines[i + 1]); i += 2
            for k in range(n):
                p = lines[i + k].split()
                etype = int(p[1]); ntags = int(p[2])
                phys = int(p[3]) if ntags >= 1 else -1
                if etype == 2 and phys == fault_tag:   # 3-node triangle
                    tris.append(tuple(map(int, p[3 + ntags:])))
            i += n; continue
        i += 1
    return nodes, tris


def tri_quality(p):
    a = math.dist(p[0], p[1]); b = math.dist(p[1], p[2]); c = math.dist(p[2], p[0])
    edges = sorted([a, b, c])
    s = (a + b + c) / 2.0
    area = math.sqrt(max(s * (s - a) * (s - b) * (s - c), 0.0))

    def angle(o, x, y):
        v1 = [x[d] - o[d] for d in range(3)]; v2 = [y[d] - o[d] for d in range(3)]
        n1 = math.sqrt(sum(t * t for t in v1)); n2 = math.sqrt(sum(t * t for t in v2))
        if n1 == 0 or n2 == 0:
            return 0.0
        cv = max(-1.0, min(1.0, sum(v1[d] * v2[d] for d in range(3)) / (n1 * n2)))
        return math.degrees(math.acos(cv))

    min_angle = min(angle(p[0], p[1], p[2]), angle(p[1], p[2], p[0]),
                    angle(p[2], p[0], p[1]))
    # triangle-q = 4*sqrt(3)*A / sum(edge^2), in [0,1], 1 = equilateral
    q = (4.0 * math.sqrt(3.0) * area) / (a * a + b * b + c * c) if area > 0 else 0.0
    return dict(min_edge=edges[0], max_edge=edges[2],
                aspect=edges[2] / edges[0] if edges[0] > 0 else float("inf"),
                min_angle=min_angle, area=area, tri_q=q)


def pct(vals, q):
    v = sorted(vals); k = (len(v) - 1) * q
    f = math.floor(k); c = math.ceil(k)
    return v[int(k)] if f == c else v[f] * (c - k) + v[c] * (k - f)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh")
    ap.add_argument("--target", nargs=3, type=float, required=True,
                    metavar=("X", "Y", "Z"))
    ap.add_argument("--radius", type=float, default=3000.0)
    ap.add_argument("--fault-tag", type=int, default=101)
    ap.add_argument("--n-closest", type=int, default=8)
    args = ap.parse_args()

    nodes, tris = read_fault_tris(args.mesh, args.fault_tag)
    print(f"nodes={len(nodes)}  fault triangles (phys {args.fault_tag})={len(tris)}")

    target = tuple(args.target)
    allq, near = [], []
    for t in tris:
        p = [nodes[n] for n in t]
        cen = tuple(sum(p[k][d] for k in range(3)) / 3.0 for d in range(3))
        qd = tri_quality(p); qd["centroid"] = cen
        allq.append(qd)
        d = math.dist(cen, target)
        if d <= args.radius:
            qd["dist"] = d; near.append(qd)

    keys = ["min_edge", "max_edge", "aspect", "min_angle", "tri_q"]

    def summary(name, qs):
        print(f"\n== {name} (n={len(qs)}) ==")
        for key in keys:
            vals = [q[key] for q in qs]
            print(f"  {key:10s} min={min(vals):.4g}  p10={pct(vals,0.10):.4g}  "
                  f"median={pct(vals,0.5):.4g}  max={max(vals):.4g}")

    summary("ALL fault triangles", allq)
    summary(f"within {args.radius:.0f} m of {target}", near)

    near.sort(key=lambda q: q["dist"])
    print(f"\n== {args.n_closest} closest fault triangles to {target} ==")
    for q in near[:args.n_closest]:
        c = q["centroid"]
        print(f"  d={q['dist']:7.1f}  min_edge={q['min_edge']:7.1f}  "
              f"aspect={q['aspect']:6.2f}  min_angle={q['min_angle']:5.2f}deg  "
              f"tri_q={q['tri_q']:.3f}  centroid=({c[0]:.0f},{c[1]:.0f},{c[2]:.0f})")

    print("\n== 5 worst min-angle fault triangles mesh-wide ==")
    for q in sorted(allq, key=lambda q: q["min_angle"])[:5]:
        c = q["centroid"]
        print(f"  min_angle={q['min_angle']:5.2f}deg  aspect={q['aspect']:6.2f}  "
              f"tri_q={q['tri_q']:.3f}  centroid=({c[0]:.0f},{c[1]:.0f},{c[2]:.0f})")


if __name__ == "__main__":
    main()
