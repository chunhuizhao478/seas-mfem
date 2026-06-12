"""extend_fault_borders.py — Phase 3.1 of PLAN_mesh_quality_safv4_remesh.

Extrude each fault's border chains PAST their contacted entities so every
contact becomes a transversal crossing that Phase 3's corefine can resolve:

  dem (trace)     : +z by --h-ext      (default 1500 m), verified strictly
                    above the DEM (auto-raised if not);
  bottom          : -z by --h-ext-bot  (default 500 m) below the bottom plane;
  fault junction  : guest side only, along the local fault-plane outward
                    direction by --h-ext-junc (default 750 m);
  interior_tip    : no extension;
  ribbon          : none observed in SAFv4 — abort if present.

Before extrusion, EVERY border edge longer than --border-split (default
750 m) is split in half recursively (pure refinement: midpoints lie on the
existing border segments, geometry unchanged).  isotropic_remeshing
preserves borders verbatim, so without this the coarse GOCAD border
segmentation (~2 km) would survive into the 500 m remesh and sink the
fault edge-length and quality gates near tips and strip bases.

Outputs (to --out-dir):
  fault_<basename>_ext.stl    (ASCII STL, %.15g)
  trace_points.xyz            (all dem-class border vertices, for the
                               remesh_graded DEM pre-pass sizing field)
  extension_report.json       (per fault per class: counts, h used,
                               verification results)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import meshio
import numpy as np

from extract_safv4_surfaces import write_ascii_stl, chains_from_edges

Z_TOL = 1e-3


def _qkey_arr(pts: np.ndarray, tol: float = 1e-6) -> np.ndarray:
    return np.round(np.asarray(pts, dtype=np.float64) / tol).astype(np.int64)


class FaultPatch:
    """Mutable triangle mesh with border-edge class tracking."""

    def __init__(self, points: np.ndarray, tris: np.ndarray):
        self.pts = [tuple(p) for p in points]
        self.tris = [tuple(int(x) for x in t) for t in tris]

    def points_array(self) -> np.ndarray:
        return np.asarray(self.pts, dtype=np.float64)

    def tris_array(self) -> np.ndarray:
        return np.asarray(self.tris, dtype=np.int64)

    def border_edge_tri(self) -> dict:
        """border (u,v) sorted -> (tri index, opposite vertex w)."""
        cnt: dict = {}
        for ti, (a, b, c) in enumerate(self.tris):
            for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
                k = (min(u, v), max(u, v))
                cnt.setdefault(k, []).append((ti, w))
        return {k: v[0] for k, v in cnt.items() if len(v) == 1}

    def split_border_edge(self, u: int, v: int, ti: int, w: int) -> int:
        """Insert midpoint of border edge (u,v); fan-split triangle ti.
        Returns new vertex index."""
        pm = tuple((np.asarray(self.pts[u]) + np.asarray(self.pts[v])) / 2.0)
        mi = len(self.pts)
        self.pts.append(pm)
        a, b, c = self.tris[ti]
        # replace (u,v,w)-cyclic triangle with two, preserving orientation
        tri = [a, b, c]
        # find positions of u and v in tri
        iu, iv = tri.index(u), tri.index(v)
        t1 = list(tri)
        t1[iv] = mi
        t2 = list(tri)
        t2[iu] = mi
        self.tris[ti] = tuple(t1)
        self.tris.append(tuple(t2))
        return mi


def classify_local_borders(patch: FaultPatch, manifest_entry: dict,
                           tol: float = 1e-3) -> dict:
    """Map each border edge (sorted local pair) to a class using the
    manifest border_classes chains, matched to local STL vertices by
    nearest neighbor (robust to the %.15g STL round-trip, unlike
    bucket-quantized keys which can straddle)."""
    from scipy.spatial import cKDTree
    pts = patch.points_array()
    tree = cKDTree(pts)
    edge2cls = {}
    for cls, info in manifest_entry["border_classes"].items():
        for ch in info["chains"]:
            co = np.asarray(ch["coords"], dtype=np.float64)
            d, idx = tree.query(co, workers=-1)
            if float(np.max(d)) > tol:
                raise RuntimeError(
                    f"chain coord match failed for class {cls}: "
                    f"max dist {float(np.max(d)):.6f} m > tol {tol}")
            for a, b in zip(idx[:-1], idx[1:]):
                a, b = int(a), int(b)
                edge2cls[(min(a, b), max(a, b))] = cls
    out = {}
    for (u, v), (ti, w) in patch.border_edge_tri().items():
        out[(u, v)] = edge2cls.get((min(u, v), max(u, v)), "interior_tip")
    return out


def dem_z_lookup(dem_pts: np.ndarray, dem_tris: np.ndarray):
    """Return z_dem(x, y) via centroid-KDTree candidate triangles +
    barycentric test."""
    from scipy.spatial import cKDTree
    cen = dem_pts[dem_tris].mean(axis=1)[:, :2]
    tree = cKDTree(cen)
    P = dem_pts

    def z_at(xy: np.ndarray) -> np.ndarray:
        out = np.full(len(xy), np.nan)
        k = min(64, len(dem_tris))
        _, idx = tree.query(xy, k=k, workers=-1)
        for i, cand in enumerate(idx):
            x, y = xy[i]
            for t in np.atleast_1d(cand):
                a, b, c = dem_tris[t]
                x1, y1 = P[a, 0], P[a, 1]
                x2, y2 = P[b, 0], P[b, 1]
                x3, y3 = P[c, 0], P[c, 1]
                det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
                if det == 0:
                    continue
                l1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det
                l2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det
                l3 = 1.0 - l1 - l2
                if min(l1, l2, l3) < -1e-9:
                    continue
                out[i] = l1 * P[a, 2] + l2 * P[b, 2] + l3 * P[c, 2]
                break
        return out

    return z_at


def host_normal_dirs(patch: FaultPatch, edges: list,
                     host_pts: np.ndarray, host_tris: np.ndarray) -> dict:
    """Per-vertex junction extrusion direction = the HOST surface normal at
    the nearest host triangle, oriented away from the guest body.

    Rationale (measured 2026-06-12): extruding the guest IN ITS OWN PLANE
    slices a band ~h/tan(dihedral) wide of host triangles nearly parallel
    to them; at the SAF-Garnet junction (13.4 deg dihedral) that
    configuration drove CGAL::PMP::corefine into a >30-min tight loop in
    its finalize step.  A host-normal strip crosses the host at 90 deg and
    intersects it only along the contact line."""
    from scipy.spatial import cKDTree
    pts = patch.points_array()
    bet = patch.border_edge_tri()
    hcen = host_pts[host_tris].mean(axis=1)
    hn = np.cross(host_pts[host_tris[:, 1]] - host_pts[host_tris[:, 0]],
                  host_pts[host_tris[:, 2]] - host_pts[host_tris[:, 0]])
    hn /= np.linalg.norm(hn, axis=1, keepdims=True)
    tree = cKDTree(hcen)
    nodes = sorted({n for e2 in edges for n in e2})
    nearest = {n: int(tree.query(pts[n])[1]) for n in nodes}
    side = 0.0
    for (u, v) in edges:
        ti, w = bet[(min(u, v), max(u, v))]
        side += float(np.dot(hn[nearest[u]], pts[w] - pts[u]))
    sgn = -1.0 if side > 0 else 1.0   # away from the guest body side
    return {n: sgn * hn[nearest[n]] for n in nodes}


def outward_inplane_dirs(patch: FaultPatch, edges: list) -> dict:
    """Per-vertex outward in-plane direction for a set of border edges:
    average over incident edges of normalize(cross(tri_normal, edge)),
    oriented away from the adjacent triangle's third vertex."""
    pts = patch.points_array()
    bet = patch.border_edge_tri()
    acc: dict = {}
    for (u, v) in edges:
        ti, w = bet[(min(u, v), max(u, v))]
        a, b, c = patch.tris[ti]
        n_t = np.cross(pts[b] - pts[a], pts[c] - pts[a])
        nm = np.linalg.norm(n_t)
        if nm == 0:
            continue
        n_t /= nm
        e = pts[v] - pts[u]
        d = np.cross(n_t, e)
        dm = np.linalg.norm(d)
        if dm == 0:
            continue
        d /= dm
        mid = (pts[u] + pts[v]) / 2.0
        if np.dot(d, pts[w] - mid) > 0:
            d = -d
        for n in (u, v):
            acc.setdefault(n, []).append(d)
    return {n: (lambda s: s / np.linalg.norm(s))(np.sum(ds, axis=0))
            for n, ds in acc.items() if np.linalg.norm(np.sum(ds, axis=0)) > 0}


def extrude_chain(patch: FaultPatch, chain: list, offsets: dict) -> int:
    """Append a strip along `chain` with per-vertex offset vectors.

    Strip triangles are wound OPPOSITE to the base triangle's traversal of
    each border edge so the composite stays orientable (CGAL's polygon-mesh
    loader rejects non-orientable input).  Returns triangles added."""
    bet = patch.border_edge_tri()
    top = {}
    pts = patch.points_array()
    for n in chain:
        t = tuple(pts[n] + offsets[n])
        top[n] = len(patch.pts)
        patch.pts.append(t)
    n_added = 0
    for a, b in zip(chain[:-1], chain[1:]):
        ti, _w = bet[(min(a, b), max(a, b))]
        x, y, z = patch.tris[ti]
        base_traverses_ab = (x, y) == (a, b) or (y, z) == (a, b) \
            or (z, x) == (a, b)
        if base_traverses_ab:
            # base uses a->b; strip must use b->a
            patch.tris.append((b, a, top[a]))
            patch.tris.append((b, top[a], top[b]))
        else:
            patch.tris.append((a, b, top[b]))
            patch.tris.append((a, top[b], top[a]))
        n_added += 2
    return n_added


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--extracted-dir", type=Path, required=True,
                    help="Phase 2 output dir (manifest.json + STLs)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--h-ext", type=float, default=1500.0)
    ap.add_argument("--h-ext-bot", type=float, default=500.0)
    ap.add_argument("--h-ext-junc", type=float, default=750.0)
    ap.add_argument("--eps-detach", type=float, default=30.0,
                    help="pull junction/bottom border vertices back toward "
                         "the guest body by this much before extruding, so "
                         "the guest BODY no longer touches the host "
                         "bit-exactly (exact near-parallel contact drove "
                         "PMP::corefine into a tight finalize loop); the "
                         "strip is lengthened by the same amount")
    ap.add_argument("--border-split", type=float, default=750.0,
                    help="recursively split border edges longer than this "
                         "before extrusion (pure refinement)")
    ap.add_argument("--tol", type=float, default=1e-3)
    args = ap.parse_args(argv)

    man = json.loads((args.extracted_dir / "manifest.json").read_text())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    entries = {e["basename"]: e for e in man["meshes"]}
    fault_names = sorted(n for n in entries if n.startswith("fault_"))

    dem_mesh = meshio.read(str(args.extracted_dir / "dem.stl"))
    dem_pts = np.asarray(dem_mesh.points, dtype=np.float64)
    dem_tris = np.asarray(dem_mesh.cells_dict["triangle"], dtype=np.int64)
    z_at = dem_z_lookup(dem_pts, dem_tris)

    bbot = entries["bottom"]["bbox"]
    if abs(bbot["z_lo"] - bbot["z_hi"]) > Z_TOL:
        print("error: bottom surface is not planar per manifest", file=sys.stderr)
        return 1
    z_bottom = bbot["z_lo"]

    # guest set per fault (junction class extruded on the guest only)
    guests = {j["guest"]: j["host"] for j in man["junctions"].values()}

    report = {"h_ext": args.h_ext, "h_ext_bot": args.h_ext_bot,
              "h_ext_junc": args.h_ext_junc,
              "border_split": args.border_split, "z_bottom": z_bottom,
              "faults": {}}
    trace_pts_out = []

    for fn in fault_names:
        e = entries[fn]
        mesh = meshio.read(str(args.extracted_dir / e["file"]))
        patch = FaultPatch(np.asarray(mesh.points, dtype=np.float64),
                           np.asarray(mesh.cells_dict["triangle"],
                                      dtype=np.int64))
        cls_map = classify_local_borders(patch, e, args.tol)
        if any(c.startswith("ribbon") for c in cls_map.values()):
            print(f"ABORT: {fn} has a ribbon-contact border edge — not "
                  f"expected in SAFv4.", file=sys.stderr)
            return 2

        # sanity: every manifest class except interior_tip must be present
        man_classes = set(e["border_classes"])
        loc_classes = set(cls_map.values())
        missing = {c for c in man_classes if c not in loc_classes}
        if missing:
            print(f"ABORT: {fn}: manifest border classes {missing} not "
                  f"recovered from the STL (coord-match failure).",
                  file=sys.stderr)
            return 2

        # ---- recursive border split ---------------------------------
        n_split = 0
        work = True
        while work:
            work = False
            pts = patch.points_array()
            bet = patch.border_edge_tri()
            for (u, v), (ti, w) in sorted(bet.items()):
                L = float(np.linalg.norm(pts[u] - pts[v]))
                if L > args.border_split:
                    cls = cls_map.get((u, v), cls_map.get((v, u),
                                                          "interior_tip"))
                    mi = patch.split_border_edge(u, v, ti, w)
                    cls_map[(min(u, mi), max(u, mi))] = cls
                    cls_map[(min(mi, v), max(mi, v))] = cls
                    cls_map.pop((u, v), None)
                    n_split += 1
                    work = True
                    break  # topology changed; rebuild maps

        # group post-split border edges per class
        by_class: dict = {}
        for k, cls in cls_map.items():
            by_class.setdefault(cls, []).append(k)
        frep = {"n_border_split": n_split,
                "classes": {c: len(v) for c, v in sorted(by_class.items())}}

        # ---- extrusions ---------------------------------------------
        pts = patch.points_array()
        n_strip_tris = 0
        # dem (trace): +z, verified above DEM with auto-raise
        h_ext = args.h_ext
        if "dem" in by_class:
            for attempt in range(3):
                chains = chains_from_edges(np.asarray(by_class["dem"]))
                # verification BEFORE extrusion commits: sample tops
                ok = True
                worst = 0.0
                for ch in chains:
                    co = pts[ch]
                    seg = np.linalg.norm(np.diff(co, axis=0), axis=1)
                    n_samp = np.maximum(1, (seg // 100).astype(int))
                    samples = [co[0]]
                    for i2 in range(len(ch) - 1):
                        for s2 in range(1, n_samp[i2] + 1):
                            samples.append(co[i2] + (co[i2 + 1] - co[i2])
                                           * s2 / n_samp[i2])
                    S = np.asarray(samples)
                    zd = z_at(S[:, :2])
                    top = S[:, 2] + h_ext
                    bad = np.nan_to_num(zd, nan=-np.inf) >= top
                    if bad.any():
                        ok = False
                        worst = max(worst, float((zd - top)[bad].max()))
                if ok:
                    break
                h_ext = h_ext + worst + 500.0
                print(f"  {fn}: raising --h-ext to {h_ext:.0f} "
                      f"(DEM violation {worst:.0f} m)")
            for ch in chains:
                offs = {n: np.array([0.0, 0.0, h_ext]) for n in ch}
                n_strip_tris += extrude_chain(patch, ch, offs)
            frep["h_ext_used"] = h_ext
            # trace sizing sites (post-split)
            for k in by_class["dem"]:
                trace_pts_out.extend([tuple(pts[k[0]]), tuple(pts[k[1]])])
        # bottom: NO strip.  The border lies EXACTLY in the bottom plane
        # (validated below; midpoint splits stay in-plane), so the contact
        # is a tangential edge contact: corefine inserts the border line
        # into the bottom mesh conformally with only a handful of short
        # segments.  A -z strip with an eps offset instead slices ~30 m
        # slivers along the entire 100 km border (thousands of sub-floor
        # polyline edges -> cascading-collapse blow-up in polyline_cleanup,
        # measured 2026-06-12).
        if "bottom" in by_class:
            zs = np.concatenate([[pts[k[0]][2], pts[k[1]][2]]
                                 for k in by_class["bottom"]])
            if np.abs(zs - z_bottom).max() > 1.0:
                print(f"ABORT: {fn} bottom-class border nodes deviate "
                      f"{np.abs(zs - z_bottom).max():.3f} m from the bottom "
                      f"plane.", file=sys.stderr)
                return 2
            frep["bottom_contact_kept_inplane"] = len(by_class["bottom"])
        # junctions: guest only, in-plane outward
        for cls in sorted(by_class):
            if not cls.startswith("fault_"):
                continue
            if f"{fn} -> {cls}" not in man["junctions"]:
                continue  # this fault is not the guest for that pair
            hmesh = meshio.read(str(args.extracted_dir / entries[cls]["file"]))
            dirs = host_normal_dirs(
                patch, by_class[cls],
                np.asarray(hmesh.points, dtype=np.float64),
                np.asarray(hmesh.cells_dict["triangle"], dtype=np.int64))
            eps = args.eps_detach
            displaced = set()
            for ch in chains_from_edges(np.asarray(by_class[cls])):
                for n in ch:
                    d = dirs.get(n)
                    if d is not None and n not in displaced:
                        p = np.asarray(patch.pts[n])
                        patch.pts[n] = tuple(p - eps * d)
                        displaced.add(n)
                offs = {n: dirs.get(n, np.array([0.0, 0.0, 0.0]))
                        * (args.h_ext_junc + eps) for n in ch}
                n_strip_tris += extrude_chain(patch, ch, offs)
            frep[f"junction_extruded_{cls}"] = len(by_class[cls])

        frep["n_strip_tris"] = n_strip_tris
        out = args.out_dir / f"{fn}_ext.stl"
        write_ascii_stl(out, patch.points_array(), patch.tris_array(),
                        name=fn)
        report["faults"][fn] = frep
        print(f"  {fn}: split {n_split} border edges, +{n_strip_tris} "
              f"strip tris -> {out.name}")

    with open(args.out_dir / "trace_points.xyz", "w") as f:
        seen = set()
        for p in trace_pts_out:
            if p in seen:
                continue
            seen.add(p)
            f.write(f"{p[0]:.15g} {p[1]:.15g} {p[2]:.15g}\n")
    with open(args.out_dir / "extension_report.json", "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"wrote trace_points.xyz ({len(seen)} sites) and "
          f"extension_report.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
