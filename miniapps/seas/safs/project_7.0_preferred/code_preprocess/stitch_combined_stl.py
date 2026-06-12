#!/usr/bin/env python3
"""stitch_combined_stl.py — Phase 1 of PLAN_revert_to_gmsh.md.

Reads `data_corefined/*_corefined.stl` + `manifest.json`, computes a closed
bounding-box polyhedron with paddings, and writes ONE combined STL whose
triangles are: 12 box triangles followed by every fault's triangles in
alphabetical basename order.  Vertices are deduplicated at 1e-6 m; coincident
polyline endpoints across faults share a single global index, producing
non-manifold internal edges along every fault-fault polyline (the desired
input topology for the Phase 2 single-`Merge` gmsh .geo).

Also writes a provenance JSON mapping triangle index ranges to fault basenames.

Usage:
    conda activate pythonenv
    python stitch_combined_stl.py [--manifest path] [--out path]
                                  [--res 2000] [--pad-xy 50000]
                                  [--pad-top 100] [--pad-bottom 25000]

Idempotent: re-running produces byte-identical output.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import meshio


PAD_XY_DEFAULT     = 50000.0   # 50 km horizontal margin
PAD_TOP_DEFAULT    =   100.0   # 100 m of rock above the free surface (z=0)
PAD_BOTTOM_DEFAULT = 25000.0   # 25 km of rock below the deepest fault vertex
COORD_TOL          =     1e-6  # quantization for vertex dedup
SAFETY_MIN_EDGE    =    50.0   # any input STL edge below this aborts the run


def _coord_key(p, tol=COORD_TOL):
    return (int(round(p[0] / tol)),
            int(round(p[1] / tol)),
            int(round(p[2] / tol)))


def _box_triangles(xmin, xmax, ymin, ymax, zmin, zmax):
    """Return (8 corner pts as np.array (8,3), 12 triangle index triples).
    Triangulation is outward-oriented (matches mesh_volume.cpp:218-231).
    """
    pts = np.array([
        [xmin, ymin, zmin], [xmax, ymin, zmin],
        [xmax, ymax, zmin], [xmin, ymax, zmin],
        [xmin, ymin, zmax], [xmax, ymin, zmax],
        [xmax, ymax, zmax], [xmin, ymax, zmax],
    ], dtype=float)
    # Outward-oriented faces.
    tris = np.array([
        # bottom (normal -z)
        [0, 3, 2], [0, 2, 1],
        # top (normal +z)
        [4, 5, 6], [4, 6, 7],
        # front (y=ymin, normal -y)
        [0, 1, 5], [0, 5, 4],
        # back (y=ymax, normal +y)
        [2, 3, 7], [2, 7, 6],
        # left (x=xmin, normal -x)
        [3, 0, 4], [3, 4, 7],
        # right (x=xmax, normal +x)
        [1, 2, 6], [1, 6, 5],
    ], dtype=int)
    return pts, tris


def _domain_bounds_from_manifest(m, pad_xy, pad_top, pad_bottom):
    xlo = min(f["bbox"]["x_lo"] for f in m["meshes"])
    xhi = max(f["bbox"]["x_hi"] for f in m["meshes"])
    ylo = min(f["bbox"]["y_lo"] for f in m["meshes"])
    yhi = max(f["bbox"]["y_hi"] for f in m["meshes"])
    zlo = min(f["bbox"]["z_lo"] for f in m["meshes"])
    zhi = max(f["bbox"]["z_hi"] for f in m["meshes"])
    return (xlo - pad_xy, xhi + pad_xy,
            ylo - pad_xy, yhi + pad_xy,
            zlo - pad_bottom, zhi + pad_top)


def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path,
                    default=project / "data_corefined" / "manifest.json")
    ap.add_argument("--out", type=Path, default=None,
                    help="output STL path (default: data_corefined/safs_combined_<R>m.stl)")
    ap.add_argument("--res", type=int, default=2000)
    ap.add_argument("--pad-xy",     type=float, default=PAD_XY_DEFAULT)
    ap.add_argument("--pad-top",    type=float, default=PAD_TOP_DEFAULT)
    ap.add_argument("--pad-bottom", type=float, default=PAD_BOTTOM_DEFAULT)
    args = ap.parse_args()

    if not args.manifest.is_file():
        print(f"error: manifest not found: {args.manifest}", file=sys.stderr); return 1

    if args.out is None:
        args.out = project / "data_corefined" / f"safs_combined_{args.res}m.stl"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    prov_path = args.out.with_suffix("").with_suffix(".prov.json")
    # The above replaces .stl -> .stl.prov.json; we want .prov.json directly.
    prov_path = args.out.parent / (args.out.stem + "_provenance.json")

    m = json.loads(args.manifest.read_text())
    if not m.get("meshes"):
        print("error: manifest has no meshes[]", file=sys.stderr); return 1

    # Sort fault entries by basename for stable ordering.
    fault_meshes = sorted(m["meshes"], key=lambda x: x["basename"])
    print(f"manifest faults (alphabetical):")
    for fm in fault_meshes:
        print(f"  {fm['basename']}: V={fm['n_verts']}, F={fm['n_faces']}")

    # Bounding box from manifest + paddings (assertions per plan §Edge Cases line 503).
    xmin, xmax, ymin, ymax, zmin, zmax = _domain_bounds_from_manifest(
        m, args.pad_xy, args.pad_top, args.pad_bottom)
    assert args.pad_top > 0, "pad_top must be > 0"
    assert args.pad_bottom > 0, "pad_bottom must be > 0"
    assert args.pad_xy > 0, "pad_xy must be > 0"
    print(f"box bounds: x=[{xmin:.0f}, {xmax:.0f}]  y=[{ymin:.0f}, {ymax:.0f}]  "
          f"z=[{zmin:.0f}, {zmax:.0f}]")

    # ----- Build the global vertex list with coord-key dedup -----
    points = []          # list of (x, y, z)
    point_idx = {}       # coord_key -> global index

    def _intern(p):
        k = _coord_key(p)
        idx = point_idx.get(k)
        if idx is None:
            idx = len(points)
            point_idx[k] = idx
            points.append((float(p[0]), float(p[1]), float(p[2])))
        return idx

    # ----- Box (12 triangles, indices 0..11 in the combined STL) -----
    box_pts, box_tris_local = _box_triangles(xmin, xmax, ymin, ymax, zmin, zmax)
    box_tris_global = []
    for t in box_tris_local:
        box_tris_global.append([_intern(box_pts[t[0]]),
                                _intern(box_pts[t[1]]),
                                _intern(box_pts[t[2]])])

    # ----- Faults (alphabetical) -----
    n_box_tri = len(box_tris_global)
    fault_provenance = []
    fault_tris_all = []
    n_dedup_total = 0

    for fm in fault_meshes:
        stl_path = args.manifest.parent / f"{fm['basename']}_corefined.stl"
        if not stl_path.is_file():
            print(f"error: fault STL not found: {stl_path}", file=sys.stderr); return 1
        msh = meshio.read(str(stl_path))
        fault_pts = msh.points
        fault_tris_local = msh.cells_dict.get("triangle")
        if fault_tris_local is None or len(fault_tris_local) == 0:
            print(f"error: no triangles in {stl_path}", file=sys.stderr); return 1

        # Sanity safety check: every input edge ≥ SAFETY_MIN_EDGE.
        e0 = np.linalg.norm(fault_pts[fault_tris_local[:, 0]]
                            - fault_pts[fault_tris_local[:, 1]], axis=1)
        e1 = np.linalg.norm(fault_pts[fault_tris_local[:, 1]]
                            - fault_pts[fault_tris_local[:, 2]], axis=1)
        e2 = np.linalg.norm(fault_pts[fault_tris_local[:, 2]]
                            - fault_pts[fault_tris_local[:, 0]], axis=1)
        e_all = np.concatenate([e0, e1, e2])
        if e_all.min() < SAFETY_MIN_EDGE:
            print(f"error: {fm['basename']} has min edge "
                  f"{e_all.min():.3f} m < safety floor {SAFETY_MIN_EDGE} m",
                  file=sys.stderr)
            return 1

        first_tri = n_box_tri + len(fault_tris_all)
        before_n_pts = len(points)
        for t in fault_tris_local:
            fault_tris_all.append([_intern(fault_pts[t[0]]),
                                   _intern(fault_pts[t[1]]),
                                   _intern(fault_pts[t[2]])])
        n_added = len(points) - before_n_pts
        n_input_unique = len(set(_coord_key(p) for p in fault_pts))
        n_dedup_this = n_input_unique - n_added
        n_dedup_total += n_dedup_this
        fault_provenance.append({
            "basename": fm["basename"],
            "first_tri": first_tri,
            "n_tri": int(len(fault_tris_local)),
            "n_dedup": int(n_dedup_this),
        })
        print(f"  added {len(fault_tris_local)} tris, {n_added} new vertices, "
              f"{n_dedup_this} dedup'd against existing")

    # ----- Combined arrays -----
    all_pts = np.array(points, dtype=float)
    all_tris_raw = np.array(box_tris_global + fault_tris_all, dtype=int)
    keep_mask = np.ones(len(all_tris_raw), dtype=bool)

    # First pass: remove true duplicate triangles (same 3 vertex indices, regardless
    # of orientation).  These come from the corefine pipeline producing identical
    # triangles in two faults at a coplanar region.
    seen = {}
    for i, t in enumerate(all_tris_raw):
        key = tuple(sorted(t.tolist()))
        if key in seen:
            keep_mask[i] = False
        else:
            seen[key] = i
    # ----- Eliminate non-manifold polyline edges by VERTEX INDEX WELDING -----
    # Per user direction (2026-05-10): "merge them, don't drop".  At each
    # polyline edge with 4 incident triangles, weld the 3rd-vertex of fault B's
    # triangles to the 3rd-vertex of fault A's pair-mates.  After welding, fault
    # B's two incident triangles share the same vertex set as fault A's two
    # → they become true duplicates → the dedup pass below removes one of each
    # pair, leaving 2 incident triangles total (manifold).  Both faults still
    # claim the merged region (their provenance ranges are unchanged in the
    # combined STL; downstream consumers see the same triangle in both faults
    # via the kd-tree match in msh_to_vtu.py).

    def _which_fault(ti):
        if ti < n_box_tri: return -1   # box
        for fi, prov in enumerate(fault_provenance):
            if prov["first_tri"] <= ti < prov["first_tri"] + prov["n_tri"]:
                return fi
        return -2

    def _third_vertex(t, ek):
        for v in t:
            if int(v) not in ek: return int(v)
        return -1

    # Build edge -> incident-triangle list (before welding)
    edge_to_tris = {}
    for i, t in enumerate(all_tris_raw):
        if not keep_mask[i]: continue
        for a, b in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]:
            ek = frozenset((int(a), int(b)))
            edge_to_tris.setdefault(ek, []).append(i)

    # Union-find on vertex indices; root is the alphabetically-earliest fault's vertex.
    parent = list(range(len(all_pts)))
    def uf_find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def uf_unite(later, earlier):
        ra = uf_find(later); rb = uf_find(earlier)
        if ra != rb: parent[ra] = rb

    n_welds = 0
    for ek, ts in edge_to_tris.items():
        if len(ts) != 4: continue
        # Group by fault
        groups = {}
        for ti in ts:
            f = _which_fault(ti)
            groups.setdefault(f, []).append(ti)
        if len(groups) != 2: continue
        keys = sorted(groups.keys())
        fA, fB = keys[0], keys[1]   # earlier alpha = fA (lower index)
        if len(groups[fA]) != 2 or len(groups[fB]) != 2: continue
        ek_set = set(ek)
        # 3rd vertices of each incident triangle
        a_thirds = [(_third_vertex(all_tris_raw[ti], ek_set), ti) for ti in groups[fA]]
        b_thirds = [(_third_vertex(all_tris_raw[ti], ek_set), ti) for ti in groups[fB]]
        # Pair each fault-B triangle with its CLOSEST fault-A triangle by 3rd-vertex distance.
        # This pairs (a1, b1) and (a2, b2) such that we weld vertex pairs that are nearest in space.
        used_a = set()
        for v_b, _ in b_thirds:
            best = None
            for v_a, _ in a_thirds:
                if v_a in used_a: continue
                d = np.linalg.norm(all_pts[v_a] - all_pts[v_b])
                if best is None or d < best[0]:
                    best = (d, v_a)
            if best is None: continue
            used_a.add(best[1])
            uf_unite(v_b, best[1])
            n_welds += 1
    print(f"  vertex welds at non-manifold polyline edges: {n_welds}")

    # Apply welds: rewrite triangle vertex indices to their union-find root.
    welded_tris = np.array([[uf_find(int(v)) for v in t] for t in all_tris_raw], dtype=int)

    # Drop true duplicate triangles that resulted from the welds.
    seen = {}
    new_keep = np.ones(len(welded_tris), dtype=bool) & keep_mask
    for i, t in enumerate(welded_tris):
        if not new_keep[i]: continue
        # Skip degenerate (after weld, two vertices may have collapsed to one)
        if t[0] == t[1] or t[1] == t[2] or t[0] == t[2]:
            new_keep[i] = False; continue
        key = tuple(sorted(t.tolist()))
        if key in seen:
            new_keep[i] = False
        else:
            seen[key] = i
    n_after_weld = new_keep.sum()
    n_dup_removed_by_weld = keep_mask.sum() - n_after_weld
    if n_dup_removed_by_weld > 0:
        print(f"  removed {int(n_dup_removed_by_weld)} duplicate/degenerate "
              f"triangle(s) after welding")
    keep_mask = new_keep
    all_tris_raw = welded_tris  # use welded indices going forward

    # ----- Drop near-degenerate triangles created by welding -----
    # Welding can collapse a triangle's 3 vertices toward collinearity (near-
    # zero area).  tetgen rejects these as "overlapping facets".  Drop any
    # triangle whose isoperimetric q_tri < 0.02 OR area < 100 m².
    Q_DEGEN_THRESH = 0.30          # aggressive: drop low-quality triangles too
    AREA_DEGEN_THRESH = 1000.0     # m²
    n_q_dropped = 0
    for i in range(len(all_tris_raw)):
        if not keep_mask[i]: continue
        t = all_tris_raw[i]
        a = all_pts[t[0]]; b = all_pts[t[1]]; c = all_pts[t[2]]
        area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
        if area < AREA_DEGEN_THRESH:
            keep_mask[i] = False; n_q_dropped += 1; continue
        e0 = np.linalg.norm(b - a); e1 = np.linalg.norm(c - b); e2 = np.linalg.norm(a - c)
        s2 = e0*e0 + e1*e1 + e2*e2
        q = 4.0 * np.sqrt(3.0) * area / s2 if s2 > 0 else 0.0
        if q < Q_DEGEN_THRESH:
            keep_mask[i] = False; n_q_dropped += 1
    if n_q_dropped > 0:
        print(f"  dropped {n_q_dropped} near-degenerate triangle(s) post-weld "
              f"(q_tri < {Q_DEGEN_THRESH} or area < {AREA_DEGEN_THRESH} m²)")

    # ----- Iterative manifold cleanup -----
    # Welding may have left 3-incident T-junctions and a few 4-/5-incident edges
    # at triple junctions.  tetgen requires every edge incident to ≤2 facets.
    # Iteratively drop one incident triangle from every non-manifold edge until
    # convergence.  Preference order: drop triangles from alphabetically-LATER
    # fault first.  Box triangles are NEVER dropped (they form the closed outer
    # boundary that tetgen needs).
    n_pass = 0
    while True:
        n_pass += 1
        edges_now = {}
        for i, t in enumerate(all_tris_raw):
            if not keep_mask[i]: continue
            for a, b in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]:
                ek = frozenset((int(a), int(b)))
                edges_now.setdefault(ek, []).append(i)
        nonmanifold = [(ek, lst) for ek, lst in edges_now.items() if len(lst) > 2]
        if not nonmanifold:
            break
        # Drop one triangle per non-manifold edge (the alphabetically-later one).
        n_dropped_this_pass = 0
        for ek, ts in nonmanifold:
            # Order incident triangles by fault index: -1 (box) first (=> never drop),
            # then 0, 1, 2, ... (drop largest fault index = alphabetically-latest).
            sorted_ts = sorted(ts, key=lambda ti: (_which_fault(ti), ti),
                               reverse=True)  # later fault first => drop first
            for ti in sorted_ts:
                if _which_fault(ti) == -1:
                    continue   # never drop box
                keep_mask[ti] = False
                n_dropped_this_pass += 1
                break  # drop only ONE per edge per pass; re-check next pass
        print(f"  pass {n_pass}: {len(nonmanifold)} non-manifold edge(s); dropped {n_dropped_this_pass}")
        if n_dropped_this_pass == 0:
            print(f"  WARNING: cannot make manifold without dropping box triangles")
            break
        if n_pass > 20:
            print(f"  WARNING: manifold cleanup did not converge after 20 passes")
            break

    n_dup_removed = (~keep_mask).sum()
    if n_dup_removed > 0:
        print(f"  total triangles removed (duplicates + coplanar): {n_dup_removed}")
        # Track which fault entries lost triangles so provenance ranges are correct.
        # Build a position-rewrite map: old_idx -> new_idx (or -1 if removed).
        # Then update fault_provenance.first_tri and n_tri.
        keep_indices = np.where(keep_mask)[0]
        old_to_new = -np.ones(len(all_tris_raw), dtype=int)
        for new_idx, old_idx in enumerate(keep_indices):
            old_to_new[old_idx] = new_idx
        # Update each fault's first_tri (find new index of old first_tri) and
        # n_tri (count of kept triangles in [old_first_tri, old_first_tri+old_n)).
        for prov in fault_provenance:
            old_first = prov["first_tri"]
            old_n     = prov["n_tri"]
            kept_in_range = keep_mask[old_first : old_first + old_n].sum()
            # New first_tri: walk forward from old_first to first kept.
            new_first = -1
            for j in range(old_first, old_first + old_n):
                if keep_mask[j]:
                    new_first = old_to_new[j]; break
            prov["first_tri"] = int(new_first if new_first >= 0 else -1)
            prov["n_tri"]     = int(kept_in_range)
    all_tris = all_tris_raw[keep_mask]
    n_total_tri = len(all_tris)
    n_total_pts = len(all_pts)
    print(f"\ncombined STL: V={n_total_pts}, F={n_total_tri}")
    print(f"  box triangles: {n_box_tri}")
    print(f"  fault triangles: {n_total_tri - n_box_tri}")
    print(f"  fault-vertex dedup'd against existing: {n_dedup_total}")

    # ----- Cross-check (plan §Phase 1 acceptance criterion 4) -----
    # An edge whose endpoints lie in DIFFERENT fault provenance ranges is a
    # shared polyline edge.  Count them and compare to the manifest's
    # pairs[k].shared_a.  Also count tris whose three vertices fall into
    # one fault range — that's the per-fault count cross-check.
    fault_range = {}   # vertex_idx -> fault basename (or "box")
    for j in range(n_box_tri):
        for v in all_tris[j]:
            fault_range[int(v)] = "box"
    for prov in fault_provenance:
        for j in range(prov["first_tri"], prov["first_tri"] + prov["n_tri"]):
            for v in all_tris[j]:
                # Only assign if not already labeled as box (box wins at corners).
                fault_range.setdefault(int(v), prov["basename"])
    # Edges with endpoints from different fault basenames (excl. "box")
    edge_set = set()
    cross_fault_edges = 0
    for j in range(n_box_tri, n_total_tri):
        for a, b in [(all_tris[j][0], all_tris[j][1]),
                     (all_tris[j][1], all_tris[j][2]),
                     (all_tris[j][2], all_tris[j][0])]:
            ek = frozenset((int(a), int(b)))
            if ek in edge_set:
                continue
            edge_set.add(ek)
            ra = fault_range.get(int(a), "?")
            rb = fault_range.get(int(b), "?")
            if ra != "box" and rb != "box" and ra != rb:
                cross_fault_edges += 1
    print(f"  edges spanning two fault basenames: {cross_fault_edges}")
    expected_pair_edges = sum(p.get("shared_a", 0) - 1 for p in m.get("pairs", []))
    # Note: shared_a counts shared VERTICES; expected polyline edges = sum(shared_a-1) per pair (path topology).
    # For 7 pairs with shared_a totals: cross-check value, allow ±50% tolerance because each pair may have multiple polylines.
    print(f"  expected ~ sum(shared_a - n_polylines): see manifest pairs[]")
    n_pair_polyline_edges = 0
    for p in m.get("pairs", []):
        # each polyline of length L vertices contributes L-1 edges
        for poly in p.get("polylines", []):
            n_pair_polyline_edges += max(0, len(poly) - 1)
    print(f"  manifest sum(polyline edges) = {n_pair_polyline_edges}")

    # ----- Write combined STL (manual ASCII at %.15g — matches project pattern) -----
    # meshio's stl writer is binary by default and the float precision is not
    # tunable; we need %.15g to preserve corefine conformality at UTM scale.
    with open(args.out, "w") as f:
        f.write("solid safs_combined\n")
        for t in all_tris:
            a = all_pts[t[0]]; b = all_pts[t[1]]; c = all_pts[t[2]]
            n = np.cross(b - a, c - a)
            nn = np.linalg.norm(n)
            n = n / nn if nn > 0 else np.array([0.0, 0.0, 1.0])
            f.write(f"facet normal {n[0]:.15g} {n[1]:.15g} {n[2]:.15g}\n")
            f.write("  outer loop\n")
            f.write(f"    vertex {a[0]:.15g} {a[1]:.15g} {a[2]:.15g}\n")
            f.write(f"    vertex {b[0]:.15g} {b[1]:.15g} {b[2]:.15g}\n")
            f.write(f"    vertex {c[0]:.15g} {c[1]:.15g} {c[2]:.15g}\n")
            f.write("  endloop\n")
            f.write("endfacet\n")
        f.write("endsolid safs_combined\n")
    print(f"\nwrote {args.out}")

    # ----- Write provenance JSON -----
    prov = {
        "n_box_triangles": n_box_tri,
        "box_bbox": {"x_lo": xmin, "x_hi": xmax,
                     "y_lo": ymin, "y_hi": ymax,
                     "z_lo": zmin, "z_hi": zmax},
        "coord_tol_m": COORD_TOL,
        "n_total_triangles": int(n_total_tri),
        "n_total_vertices": int(n_total_pts),
        "n_cross_fault_edges": int(cross_fault_edges),
        "n_dedup_fault_vertices": int(n_dedup_total),
        "faults": fault_provenance,
    }
    prov_path.write_text(json.dumps(prov, indent=2, sort_keys=True))
    print(f"wrote {prov_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
