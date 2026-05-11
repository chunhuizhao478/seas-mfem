#!/usr/bin/env python3
"""merge_corefined_faults.py — Transitive vertex-snap merge of multiple
corefined fault STLs into ONE STL.  Use when two (or more) faults are
geologically close enough that treating them as separate causes residual
sliver tets at their interface in the downstream tetgen mesh.

Algorithm:
    1. Read each input STL into a global vertex list with coord-key dedup.
    2. Build a kd-tree over INTER-fault vertex pairs (verts from
       different source STLs).
    3. For every inter-fault pair within `--snap-tol` of each other, snap
       BOTH verts to their midpoint (union-find merge).
    4. Apply remap to all triangles.  Drop:
         - Triangles with collapsed verts (≥2 verts identical) → degenerate
         - Triangles whose vertex set duplicates another → keep one only
    5. Write merged STL with all surviving triangles in one block.
    6. Update manifest.json: replace the merged faults' entries with a
       single new entry; keep non-merged faults unchanged.

Usage:
    conda activate pythonenv
    python merge_corefined_faults.py \\
        --merge SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m \\
                SAFS-SAFZ-MULT-Southern_San_Andreas_fault_and_Banning-CFM6_2000m \\
                SAFS-SAFZ-SBMT-Garnet_Hill_fault-CFM6_2000m \\
                SAFS-SAFZ-SBMT-San_Andreas_fault-CFM6_2000m \\
        --new-name SAFS-SAFZ-BigSAF_2000m \\
        --snap-tol 500
"""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import meshio
from scipy.spatial import cKDTree


def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corefined-dir", type=Path,
                    default=project / "data_corefined")
    ap.add_argument("--manifest", type=Path,
                    default=project / "data_corefined" / "manifest.json")
    ap.add_argument("--merge", nargs="+", required=True,
                    help="space-separated list of fault basenames to merge")
    ap.add_argument("--new-name", type=str, required=True,
                    help="basename for the merged fault (no _corefined.stl suffix)")
    ap.add_argument("--snap-tol", type=float, default=500.0,
                    help="snap any inter-fault vertex pair within this "
                         "distance to their midpoint [m] (default 500)")
    ap.add_argument("--coord-tol", type=float, default=1e-3,
                    help="coord-key dedup tolerance for INTRA-fault vertex "
                         "lookup [m] (default 1e-3)")
    ap.add_argument("--out-manifest", type=Path, default=None,
                    help="path to write updated manifest.json "
                         "(default: overwrite --manifest)")
    args = ap.parse_args()

    if not args.manifest.is_file():
        print(f"error: manifest not found: {args.manifest}", file=sys.stderr); return 1
    if args.out_manifest is None:
        args.out_manifest = args.manifest

    m = json.loads(args.manifest.read_text())
    if not m.get("meshes"):
        print("error: manifest has no meshes[]", file=sys.stderr); return 1

    # Validate that every requested basename exists.
    known_basenames = {fm["basename"] for fm in m["meshes"]}
    missing = [b for b in args.merge if b not in known_basenames]
    if missing:
        print(f"error: basename(s) not in manifest: {missing}", file=sys.stderr); return 1
    if len(args.merge) < 2:
        print("error: --merge must list at least 2 basenames", file=sys.stderr); return 1

    print(f"merging {len(args.merge)} faults into '{args.new_name}':")
    for b in args.merge:
        print(f"  - {b}")
    print(f"snap-tol = {args.snap_tol} m")

    # ----- read each STL into a global vertex list with INTRA dedup -----
    points       = []
    point_idx    = {}
    fault_origin = []   # source fault index (0..N-1) per global vertex
    triangles    = []   # list of [vid0, vid1, vid2]
    tri_origin   = []   # source fault index per triangle

    def _coord_key(p):
        return (int(round(p[0] / args.coord_tol)),
                int(round(p[1] / args.coord_tol)),
                int(round(p[2] / args.coord_tol)))

    def _intern(p, src_idx):
        k = _coord_key(p)
        idx = point_idx.get(k)
        if idx is None:
            idx = len(points)
            point_idx[k] = idx
            points.append((float(p[0]), float(p[1]), float(p[2])))
            fault_origin.append(src_idx)
        return idx

    for src_idx, basename in enumerate(args.merge):
        stl_path = args.corefined_dir / f"{basename}_corefined.stl"
        if not stl_path.is_file():
            print(f"error: STL not found: {stl_path}", file=sys.stderr); return 1
        msh = meshio.read(str(stl_path))
        ftris = msh.cells_dict.get("triangle")
        if ftris is None or len(ftris) == 0:
            print(f"error: no triangles in {stl_path}", file=sys.stderr); return 1
        n_before = len(points)
        for t in ftris:
            ia = _intern(msh.points[t[0]], src_idx)
            ib = _intern(msh.points[t[1]], src_idx)
            ic = _intern(msh.points[t[2]], src_idx)
            triangles.append([ia, ib, ic])
            tri_origin.append(src_idx)
        n_added = len(points) - n_before
        print(f"  loaded {stl_path.name}: V={len(msh.points)} → +{n_added} new globals  F={len(ftris)}")

    pts_arr   = np.asarray(points, dtype=float)
    tris_arr  = np.asarray(triangles, dtype=np.int32)
    fault_origin = np.asarray(fault_origin, dtype=np.int32)
    tri_origin = np.asarray(tri_origin, dtype=np.int32)
    print(f"\nbefore inter-fault snap: V={len(pts_arr)}  F={len(tris_arr)}")

    # ----- INTER-fault snap: GREEDY pair-by-pair, NO chains -----
    # Sort pairs by distance ascending; for each pair, accept only if NEITHER
    # vert is already in a snapped pair.  This bounds every snapped cluster to
    # exactly 2 verts → cluster diameter ≤ args.snap_tol → no long edges
    # created in incident triangles.
    #
    # Earlier union-find allowed transitive chains: V1—V2—V3—V4—V5 each
    # within snap_tol, but V1↔V5 separated by 4·snap_tol.  Centroid of such
    # chain is several km from the chain endpoints; root vert (lowest index)
    # gets moved to centroid; triangles that referenced the root now have
    # multi-km edges to the original neighbors of the root → bad fault tris
    # → bad sliver tets stacked on them.  Greedy avoids all that.
    tree = cKDTree(pts_arr)
    pairs = tree.query_pairs(r=args.snap_tol, output_type='ndarray')
    inter_pairs = [(int(a), int(b)) for a, b in pairs
                   if fault_origin[a] != fault_origin[b]]
    # Compute distance for each pair; sort ascending.
    pair_dist = []
    for a, b in inter_pairs:
        d = float(np.linalg.norm(pts_arr[a] - pts_arr[b]))
        pair_dist.append((d, a, b))
    pair_dist.sort()
    print(f"  candidate pairs within {args.snap_tol}m: {len(pairs)}")
    print(f"  INTER-fault pairs (candidates): {len(inter_pairs)}")

    used = np.zeros(len(pts_arr), dtype=bool)
    remap = np.arange(len(pts_arr), dtype=np.int64)
    accepted = []
    for d, a, b in pair_dist:
        if used[a] or used[b]:
            continue
        # Snap: lower-index vert becomes the survivor at midpoint;
        # higher-index vert is remapped to it.
        lo, hi = (a, b) if a < b else (b, a)
        used[a] = True; used[b] = True
        remap[hi] = lo
        accepted.append((lo, hi))
    print(f"  INTER-fault pairs (accepted, greedy non-chain): {len(accepted)}")

    # Move surviving (lo) vert to midpoint of its pair.  Each snapped vert
    # appears in exactly ONE accepted pair (greedy), so the move is bounded
    # by snap_tol/2 — much less than tolerance, no chain blow-up.
    new_pts = pts_arr.copy()
    for lo, hi in accepted:
        new_pts[lo] = 0.5 * (pts_arr[lo] + pts_arr[hi])

    # Apply remap to triangles, drop degenerates and duplicates.
    new_tris = []
    seen = set()
    n_drop_degen = 0
    n_drop_dup   = 0
    for t in tris_arr:
        new_t = sorted(int(remap[v]) for v in t)
        if len(set(new_t)) < 3:
            n_drop_degen += 1
            continue
        key = tuple(new_t)
        if key in seen:
            n_drop_dup += 1
            continue
        seen.add(key)
        new_tris.append(new_t)
    new_tris = np.asarray(new_tris, dtype=np.int32)

    print(f"  triangles dropped (degenerate after snap): {n_drop_degen}")
    print(f"  triangles dropped (duplicate after snap)  : {n_drop_dup}")
    print(f"after merge: V={len(np.unique(new_tris.ravel()))}  F={len(new_tris)}")

    # Compact vertex array (drop any verts no longer used).
    used_v = np.unique(new_tris.ravel())
    compact_pts = new_pts[used_v]
    inv = -np.ones(len(new_pts), dtype=np.int64)
    inv[used_v] = np.arange(len(used_v))
    compact_tris = inv[new_tris].astype(np.int32)
    print(f"after compaction: V={len(compact_pts)}  F={len(compact_tris)}")

    # ----- write merged STL -----
    out_stl = args.corefined_dir / f"{args.new_name}_corefined.stl"
    with open(out_stl, "w") as f:
        f.write(f"solid {args.new_name}\n")
        for t in compact_tris:
            a = compact_pts[t[0]]; b = compact_pts[t[1]]; c = compact_pts[t[2]]
            n = np.cross(b - a, c - a)
            mag = np.linalg.norm(n)
            n = n / mag if mag > 0 else np.array([0.0, 0.0, 1.0])
            f.write(f"facet normal {n[0]:.15g} {n[1]:.15g} {n[2]:.15g}\n")
            f.write("  outer loop\n")
            for p in (a, b, c):
                f.write(f"    vertex {p[0]:.15g} {p[1]:.15g} {p[2]:.15g}\n")
            f.write("  endloop\n")
            f.write("endfacet\n")
        f.write(f"endsolid {args.new_name}\n")
    print(f"\nwrote {out_stl}")

    # ----- update manifest -----
    # Drop entries for merged faults; add one new entry.
    new_meshes = []
    for fm in m["meshes"]:
        if fm["basename"] not in args.merge:
            new_meshes.append(fm)
    # Build entry for new merged fault.
    bb_x = (compact_pts[:, 0].min(), compact_pts[:, 0].max())
    bb_y = (compact_pts[:, 1].min(), compact_pts[:, 1].max())
    bb_z = (compact_pts[:, 2].min(), compact_pts[:, 2].max())
    new_meshes.append({
        "index": len(new_meshes),
        "input": f"{args.new_name}.off",
        "output": f"{args.new_name}_corefined.off",
        "basename": args.new_name,
        "n_verts_in": len(compact_pts),
        "n_faces_in": len(compact_tris),
        "n_verts": len(compact_pts),
        "n_faces": len(compact_tris),
        "bbox": {
            "x_lo": float(bb_x[0]), "x_hi": float(bb_x[1]),
            "y_lo": float(bb_y[0]), "y_hi": float(bb_y[1]),
            "z_lo": float(bb_z[0]), "z_hi": float(bb_z[1]),
        },
        "edge_stats": {"min": 0.0, "p1": 0.0, "median": 0.0, "max": 0.0},  # not recomputed here
        "tri_q_stats": {"min": 0.0, "p1": 0.0, "median": 0.0},             # not recomputed here
        "n_edges_below_floor": 0,
        "n_tri_q_below_threshold": 0,
        "merge_provenance": {
            "merged_basenames": list(args.merge),
            "snap_tol_m": args.snap_tol,
            "n_inter_pairs_snapped": len(inter_pairs),
            "n_dropped_degenerate": n_drop_degen,
            "n_dropped_duplicate":  n_drop_dup,
        },
    })
    # Re-index meshes alphabetically (stable order for downstream tools).
    new_meshes = sorted(new_meshes, key=lambda x: x["basename"])
    for i, fm in enumerate(new_meshes):
        fm["index"] = i

    # Drop pairs[] entries that referenced merged faults (since they're
    # now intra-fault).  Keep pairs not involving any merged fault.
    old_idx_to_new = {}  # old manifest index → new index, or None if merged
    old_basename_to_new = {}
    for old_fm in m["meshes"]:
        if old_fm["basename"] in args.merge:
            old_basename_to_new[old_fm["basename"]] = "MERGED"
        else:
            new_idx = next(i for i, nfm in enumerate(new_meshes)
                           if nfm["basename"] == old_fm["basename"])
            old_basename_to_new[old_fm["basename"]] = new_idx
        old_idx_to_new[old_fm["index"]] = old_basename_to_new[old_fm["basename"]]
    new_pairs = []
    n_pairs_dropped = 0
    for p in m.get("pairs", []):
        bi = m["meshes"][p["i"]]["basename"]
        bj = m["meshes"][p["j"]]["basename"]
        ri = old_basename_to_new.get(bi)
        rj = old_basename_to_new.get(bj)
        if ri == "MERGED" and rj == "MERGED":
            n_pairs_dropped += 1
            continue
        # If exactly one is merged, the pair becomes (merged_new_idx, other)
        if ri == "MERGED":
            ri = next(i for i, nfm in enumerate(new_meshes) if nfm["basename"] == args.new_name)
        if rj == "MERGED":
            rj = next(i for i, nfm in enumerate(new_meshes) if nfm["basename"] == args.new_name)
        np_pair = dict(p)
        np_pair["i"], np_pair["j"] = (ri, rj) if ri < rj else (rj, ri)
        new_pairs.append(np_pair)
    print(f"\nmanifest pairs[] entries dropped (intra-fault now): {n_pairs_dropped}")

    new_manifest = dict(m)
    new_manifest["meshes"] = new_meshes
    new_manifest["pairs"] = new_pairs
    args.out_manifest.write_text(json.dumps(new_manifest, indent=2))
    print(f"wrote updated manifest: {args.out_manifest}")
    print(f"  faults now: {[fm['basename'] for fm in new_meshes]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
