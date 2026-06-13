"""regenerate_boundary.py — rebuild the SAFv4 boundary shell (DEM top,
flat bottom, 4 vertical side walls) at a MODERATE, well-shaped resolution,
watertight by construction.

Motivation (user feedback 2026-06-12): the graded DEM remesh left far-field
triangles up to 43 km and the raw GOCAD bottom/ribbons were 6-10 km median.
Those coarse, high-aspect boundary surfaces produce badly-shaped boundary
tets in the volume mesh.  The fix is a moderate, near-uniform boundary
(~2.5 km far field, fine only near the fault trace on the DEM) so the
adjacent tets stay well shaped.

Geometry (measured): the domain footprint is a 4-sided polygon; each of the
4 ribbons is a straight vertical wall; the bottom is exactly planar; the DEM
is a single-valued (x,y) heightfield.  ONE resampled footprint drives the
DEM rim, the bottom rim, and the ribbon tops/bottoms, so every shared edge
is identical on both incident surfaces -> watertight.

The DEM is regenerated (it is corefine-coupled to the faults downstream;
corefine re-establishes the fault trace, so the DEM triangulation here only
needs to be graded fine near the trace and moderate elsewhere).

CLI
---
    python regenerate_boundary.py --extracted-dir DIR --trace-points XYZ
        --out-dir DIR [--h-near 500] [--h-far 2500] [--d-near 1500]
        [--d-far 30000] [--boundary-size 2500] [--z-bottom -39329.4]

Writes dem.stl, bottom.stl, ribbon_<NE30|NE120|NW60|NW150>.stl + a
regen_report.json.  The fault STLs and manifest are NOT touched.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import meshio
import numpy as np

from extract_safv4_surfaces import write_ascii_stl, chains_from_edges, border_edges


def load_tri(path):
    m = meshio.read(str(path))
    return (np.asarray(m.points, dtype=np.float64),
            np.asarray(m.cells_dict["triangle"], dtype=np.int64))


def ordered_rim(pts, tris):
    """Ordered closed rim polygon (vertex indices) of a heightfield/planar
    surface from its single border loop."""
    be = border_edges(tris)
    chains = chains_from_edges(be)
    chains.sort(key=len, reverse=True)
    ch = chains[0]
    if ch[0] == ch[-1]:
        ch = ch[:-1]
    return ch


def find_corners(poly_xy, n_corners=4):
    """Indices of the `n_corners` sharpest turns in a closed xy polygon."""
    n = len(poly_xy)
    ang = np.zeros(n)
    for i in range(n):
        a = poly_xy[(i - 1) % n]
        b = poly_xy[i]
        c = poly_xy[(i + 1) % n]
        v1 = b - a
        v2 = c - b
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 == 0 or n2 == 0:
            ang[i] = 0
            continue
        cosang = np.clip(np.dot(v1, v2) / (n1 * n2), -1, 1)
        ang[i] = np.arccos(cosang)        # turning angle
    # greedily pick the n_corners largest turns, spaced apart
    order = np.argsort(-ang)
    picked = []
    for idx in order:
        if all(min(abs(idx - p), n - abs(idx - p)) > n // (3 * n_corners)
               for p in picked):
            picked.append(int(idx))
        if len(picked) == n_corners:
            break
    return sorted(picked)


def resample_polyline(P, target):
    """Resample an open polyline (M,3) to segments <= target, keeping the
    endpoints; returns (K,3) including both ends."""
    out = [P[0]]
    for i in range(len(P) - 1):
        a, b = P[i], P[i + 1]
        L = np.linalg.norm(b - a)
        n = max(1, int(np.ceil(L / target)))
        for k in range(1, n + 1):
            out.append(a + (b - a) * k / n)
    return np.asarray(out)


def _in_poly(xy, poly):
    inside = np.zeros(len(xy), dtype=bool)
    x, y = xy[:, 0], xy[:, 1]
    x0, y0 = poly[:, 0], poly[:, 1]
    x1, y1 = np.roll(x0, -1), np.roll(y0, -1)
    for i in range(len(poly)):
        dy = y1[i] - y0[i]
        if dy == 0.0:
            continue
        cond = (y0[i] > y) != (y1[i] > y)
        xin = (x1[i] - x0[i]) * (y - y0[i]) / dy + x0[i]
        inside ^= cond & (x < xin)
    return inside


def quality_polygon(boundary_uv, target, hole_uv=None):
    """Quality-meshed triangulation of a simple polygon (boundary_uv, in
    loop order) at edge ~`target` via Shewchuk's `triangle` with min-angle
    30 deg + area bound; `Y` forbids Steiner on the boundary so the
    polygon's `len(boundary_uv)` vertices keep indices 0..nb-1 (-> the
    shared rim stays identical across surfaces -> watertight).  Returns
    (verts_uv (N,2), tris (M,3)); verts[:nb] == boundary_uv."""
    import triangle as tr
    nb = len(boundary_uv)
    segs = np.column_stack([np.arange(nb), (np.arange(nb) + 1) % nb])
    A = {"vertices": np.asarray(boundary_uv, dtype=np.float64),
         "segments": segs}
    area = 0.5 * target * target          # ~equilateral target area
    out = tr.triangulate(A, f"pq30a{area:.0f}Y")
    verts = np.asarray(out["vertices"], dtype=np.float64)
    tris = np.asarray(out["triangles"], dtype=np.int64)
    # boundary vertices must be unchanged & first (Y guarantees this)
    if not np.allclose(verts[:nb], boundary_uv):
        raise RuntimeError("triangle moved boundary vertices")
    return verts, tris


def delaunay_in_polygon(boundary_xy, interior_xy, z_of=None):
    """Constrained Delaunay (Shewchuk's `triangle`) of the footprint
    polygon: the boundary segments (the first len(boundary_xy) vertices,
    in loop order) are honored exactly, so the output boundary loop IS the
    footprint -> watertight with the ribbons/bottom.  `pYY` forbids Steiner
    points on segments and in the interior, so vertices keep their input
    order (footprint nodes stay at indices 0..nb-1).  Returns
    (points (N,3), tris (M,3))."""
    import triangle as tr
    nb = len(boundary_xy)
    verts = np.vstack([boundary_xy, interior_xy])
    segs = np.column_stack([np.arange(nb), (np.arange(nb) + 1) % nb])
    # q28 = min-angle 28 deg (interior Steiner only; `Y` keeps the footprint
    # rim segments un-split so the shared rim stays identical -> watertight).
    out = tr.triangulate({"vertices": verts, "segments": segs}, "pq28Y")
    pts2 = np.asarray(out["vertices"], dtype=np.float64)
    simp = np.asarray(out["triangles"], dtype=np.int64)
    if not np.allclose(pts2[:nb], boundary_xy):
        raise RuntimeError("triangle moved the footprint rim vertices")
    z = z_of(pts2) if z_of is not None else np.zeros(len(pts2))
    return np.column_stack([pts2, z]), simp


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--extracted-dir", type=Path, required=True)
    ap.add_argument("--trace-points", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--h-near", type=float, default=500.0)
    ap.add_argument("--h-far", type=float, default=2500.0)
    ap.add_argument("--d-near", type=float, default=1500.0)
    ap.add_argument("--d-far", type=float, default=30000.0)
    ap.add_argument("--boundary-size", type=float, default=2500.0)
    ap.add_argument("--z-bottom", type=float, default=-39329.4)
    args = ap.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    dem_pts, dem_tris = load_tri(args.extracted_dir / "dem.stl")
    # heightfield interpolator z(x,y) from the original DEM
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import cKDTree
    z_interp = LinearNDInterpolator(dem_pts[:, :2], dem_pts[:, 2])
    dem_xy_tree = cKDTree(dem_pts[:, :2])

    def z_of(xy):
        z = z_interp(xy)
        bad = ~np.isfinite(z)
        if bad.any():                      # outside convex hull: nearest
            _, j = dem_xy_tree.query(xy[bad])
            z[bad] = dem_pts[j, 2]
        return z

    # --- footprint: ordered DEM rim, find 4 corners, resample sides -------
    rim = ordered_rim(dem_pts, dem_tris)
    rim_xy = dem_pts[rim][:, :2]
    corners = find_corners(rim_xy, 4)
    # 4 sides between consecutive corners
    sides_xy = []
    for s in range(4):
        i0, i1 = corners[s], corners[(s + 1) % 4]
        if i1 > i0:
            seg = rim[i0:i1 + 1]
        else:
            seg = rim[i0:] + rim[:i1 + 1]
        side_pts = dem_pts[seg]
        side_rs = resample_polyline(side_pts, args.boundary_size)
        sides_xy.append(side_rs)
    # stitch resampled sides into ONE footprint loop (drop duplicate
    # corners) while tracking each side's [start, end] node range.
    foot = []
    side_ranges = []
    cur = 0
    for s in range(4):
        seg = sides_xy[s][:-1]              # drop last (shared with next)
        foot.append(seg)
        side_ranges.append((cur, cur + len(seg)))   # [start, end) in foot
        cur += len(seg)
    foot = np.vstack(foot)
    foot_xy = foot[:, :2]
    foot_top_z = z_of(foot_xy)
    n_foot = len(foot_xy)
    print(f"footprint: {len(rim)} rim nodes -> {n_foot} resampled "
          f"(corners at {corners})", flush=True)

    # assign each of the 4 footprint sides to its original ribbon tag
    rib_names = ["NE30", "NE120", "NW60", "NW150"]
    rib_trees = {}
    for rn in rib_names:
        rp, _ = load_tri(args.extracted_dir / f"ribbon_{rn}.stl")
        rib_trees[rn] = cKDTree(rp[:, :2])
    side_rib = []
    for (a, b) in side_ranges:
        mid = foot_xy[(a + b) // 2]
        dists = {rn: rib_trees[rn].query(mid)[0] for rn in rib_names}
        side_rib.append(min(dists, key=dists.get))
    print(f"  side->ribbon: {side_rib}", flush=True)

    # --- DEM interior points: graded by distance to trace -----------------
    trace = np.loadtxt(args.trace_points)
    if trace.ndim == 1:
        trace = trace.reshape(1, -1)
    tr_tree = cKDTree(trace[:, :2])

    def hsize(xy):
        d = tr_tree.query(xy)[0]
        t = np.clip((d - args.d_near) / (args.d_far - args.d_near), 0, 1)
        return args.h_near + t * (args.h_far - args.h_near)

    # Poisson-ish: coarse background grid + fine ring near trace, rejected
    # to inside-footprint, then thinned so no two points are closer than
    # 0.7*local h (a cheap blue-noise).
    xmin, ymin = foot_xy.min(0)
    xmax, ymax = foot_xy.max(0)
    cand = []
    g = args.h_far
    xs = np.arange(xmin, xmax + g, g)
    ys = np.arange(ymin, ymax + g, g)
    GX, GY = np.meshgrid(xs, ys)
    cand.append(np.column_stack([GX.ravel(), GY.ravel()]))
    # fine band near trace
    gf = args.h_near
    for tp in trace[:, :2]:
        k = int(np.ceil(args.d_near / gf))
        lx = np.arange(tp[0] - args.d_near, tp[0] + args.d_near + gf, gf)
        ly = np.arange(tp[1] - args.d_near, tp[1] + args.d_near + gf, gf)
        MX, MY = np.meshgrid(lx, ly)
        cand.append(np.column_stack([MX.ravel(), MY.ravel()]))
    cand = np.vstack(cand)
    cand = cand[_in_poly(cand, foot_xy)]
    # blue-noise thinning via a uniform spatial-hash grid (O(N)).  Cell
    # size = min local spacing (h_near); a candidate is accepted only if no
    # already-accepted point lies within 0.7*h_local, checked over the 5x5
    # cell neighborhood (covers radius up to 2 cells = 2*h_near; for
    # coarse regions 0.7*h_far may exceed that, so cells scale with h_far
    # and the neighborhood is widened accordingly).
    order = np.argsort(hsize(cand))         # keep fine points first
    cand = cand[order]
    cell = args.h_near
    reach = int(np.ceil(0.7 * args.h_far / cell)) + 1
    grid: dict = {}
    # seed with footprint nodes
    for fx in foot_xy:
        key = (int(fx[0] // cell), int(fx[1] // cell))
        grid.setdefault(key, []).append(fx)
    accepted = []
    h_all = hsize(cand)
    for p, h in zip(cand, h_all):
        cx, cy = int(p[0] // cell), int(p[1] // cell)
        r = int(np.ceil(0.7 * h / cell)) + 1
        ok = True
        for dx in range(-r, r + 1):
            for dy in range(-r, r + 1):
                for q2 in grid.get((cx + dx, cy + dy), ()):  # noqa
                    if (q2[0] - p[0]) ** 2 + (q2[1] - p[1]) ** 2 < (0.7 * h) ** 2:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                break
        if ok:
            accepted.append(p)
            grid.setdefault((cx, cy), []).append(p)
    interior_xy = np.array(accepted) if accepted else np.zeros((0, 2))
    print(f"DEM interior points: {len(interior_xy)} "
          f"(grid {len(cand)} candidates)", flush=True)

    # --- DEM triangulation (heightfield) ----------------------------------
    dem_new_pts, dem_new_tris = delaunay_in_polygon(
        foot_xy, interior_xy, z_of=z_of)
    # orient +z
    v01 = dem_new_pts[dem_new_tris[:, 1], :2] - dem_new_pts[dem_new_tris[:, 0], :2]
    v02 = dem_new_pts[dem_new_tris[:, 2], :2] - dem_new_pts[dem_new_tris[:, 0], :2]
    cw = (v01[:, 0] * v02[:, 1] - v01[:, 1] * v02[:, 0]) < 0
    dem_new_tris[cw] = dem_new_tris[cw][:, [0, 2, 1]]

    # --- bottom (flat at z_bottom): quality-meshed footprint --------------
    bot_uv, bot_tris = quality_polygon(foot_xy, args.boundary_size)
    bot_pts = np.column_stack([bot_uv, np.full(len(bot_uv), args.z_bottom)])
    # orient -z (outward = down)
    v01 = bot_pts[bot_tris[:, 1], :2] - bot_pts[bot_tris[:, 0], :2]
    v02 = bot_pts[bot_tris[:, 2], :2] - bot_pts[bot_tris[:, 0], :2]
    ccw = (v01[:, 0] * v02[:, 1] - v01[:, 1] * v02[:, 0]) > 0
    bot_tris[ccw] = bot_tris[ccw][:, [0, 2, 1]]

    # --- ribbons (vertical walls): quality-mesh each side in its plane ----
    # Each side is straight in xy, so (u=horizontal-along-side, v=z) is a
    # faithful 2D parametrization.  The 2D boundary loop is top (footprint
    # nodes, descending corner columns subdivided at boundary_size) so
    # adjacent sides subdivide the SHARED corner column identically (same
    # top_z, z_bottom, level count) -> watertight.  `triangle` quality-meshes
    # the interior with near-equilateral triangles.
    size = args.boundary_size

    def corner_levels(top_z):
        n = max(1, int(np.ceil((top_z - args.z_bottom) / size)))
        return np.array([top_z + (args.z_bottom - top_z) * k / n
                         for k in range(n + 1)])      # top..bottom inclusive

    rib_surf = {rn: None for rn in rib_names}      # ribbon_<rn> -> (pts3, tris)
    for s in range(4):
        a0 = side_ranges[s][0]
        a1 = side_ranges[(s + 1) % 4][0]               # next corner index
        top_nodes = (list(range(a0, a1 + 1)) if a1 > a0
                     else list(range(a0, n_foot)) + list(range(0, a1 + 1)))
        side_poly = foot_xy[top_nodes]                 # actual (wavy) side
        seg_len = np.linalg.norm(np.diff(side_poly, axis=0), axis=1)
        u = np.concatenate([[0.0], np.cumsum(seg_len)])  # arc-length per node
        side_len = u[-1]
        topz = foot_top_z[top_nodes]

        def uv_to_xy(uu):
            """Map arc-length uu (M,) to xy on the actual side polyline."""
            j = np.clip(np.searchsorted(u, uu, side="right") - 1, 0,
                        len(u) - 2)
            t = np.where(seg_len[j] > 0, (uu - u[j]) / seg_len[j], 0.0)
            return side_poly[j] + (side_poly[j + 1] - side_poly[j]) * t[:, None]
        # 2D boundary loop (u,v): top L->R, right column top->bottom,
        # bottom R->L, left column bottom->top
        right_v = corner_levels(topz[-1])[1:]          # below-top..bottom
        left_v = corner_levels(topz[0])[1:-1][::-1]    # interior only
        loop = []
        for i in range(len(top_nodes)):
            loop.append((u[i], topz[i]))
        for v in right_v:
            loop.append((u[-1], v))
        for i in range(len(top_nodes) - 2, -1, -1):
            loop.append((u[i], args.z_bottom))
        for v in left_v:
            loop.append((u[0], v))
        loop = np.asarray(loop)
        uv, tris = quality_polygon(loop, size)
        # map (u,v) -> 3D along the ACTUAL side polyline (handles wavy sides
        # so the ribbon top matches the DEM rim node-for-node)
        xy = uv_to_xy(np.clip(uv[:, 0], 0.0, side_len))
        p3 = np.column_stack([xy, uv[:, 1]])
        if rib_surf[side_rib[s]] is None:
            rib_surf[side_rib[s]] = (p3, tris)
        else:                                          # 2 sides same tag
            p0, t0 = rib_surf[side_rib[s]]
            rib_surf[side_rib[s]] = (np.vstack([p0, p3]),
                                     np.vstack([t0, tris + len(p0)]))

    # --- watertight check (build a combined node table) -------------------
    surf = {}
    surf["dem"] = (dem_new_pts, dem_new_tris)
    surf["bottom"] = (bot_pts, bot_tris)
    for rn in rib_names:
        if rib_surf.get(rn) is not None:
            surf[f"ribbon_{rn}"] = rib_surf[rn]

    # global weld at 1mm
    allp = []
    offsets = {}
    for nm, (p, t) in surf.items():
        offsets[nm] = len(allp)
        allp.extend(p.tolist())
    allp = np.asarray(allp)
    q = np.round(allp / 1e-3).astype(np.int64)
    _, inv = np.unique(q, axis=0, return_inverse=True)
    glob_tris = []
    for nm, (p, t) in surf.items():
        gt = inv[t + offsets[nm]]
        glob_tris.append(gt)
    GT = np.concatenate(glob_tris)
    ge = np.sort(np.concatenate([GT[:, [0, 1]], GT[:, [1, 2]], GT[:, [2, 0]]]),
                 axis=1)
    uge, cge = np.unique(ge, axis=0, return_counts=True)
    n_open = int((cge == 1).sum())
    n_nonman = int((cge > 2).sum())
    print(f"boundary shell: {len(uge)} edges, open={n_open}, >2={n_nonman}",
          flush=True)
    if n_open or n_nonman:
        print("ABORT: regenerated boundary shell not watertight.",
              file=sys.stderr)
        return 2

    # --- write STLs --------------------------------------------------------
    for nm, (p, t) in surf.items():
        used = np.unique(t.ravel())
        remap = -np.ones(len(p), dtype=np.int64)
        remap[used] = np.arange(used.size)
        write_ascii_stl(args.out_dir / f"{nm}.stl", p[used], remap[t], name=nm)

    def edstats(p, t):
        e = np.concatenate([np.linalg.norm(p[t[:, a]] - p[t[:, b]], axis=1)
                            for a, b in [(0, 1), (1, 2), (2, 0)]])
        return float(e.min()), float(np.median(e)), float(e.max())

    report = {"footprint_nodes": n_foot, "corners": corners,
              "shell_edges": int(len(uge)), "open": n_open,
              "surfaces": {}}
    print("\nregenerated boundary surface edge stats [min/med/max m | n_tri]:")
    for nm, (p, t) in surf.items():
        mn, md, mx = edstats(p, t)
        report["surfaces"][nm] = {"n_tri": int(len(t)),
                                  "edge_min": mn, "edge_med": md,
                                  "edge_max": mx}
        print(f"  {nm:14s}: {mn:6.0f} / {md:6.0f} / {mx:7.0f} | {len(t):,}")
    with open(args.out_dir / "regen_report.json", "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"\nwrote {len(surf)} boundary STLs to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
