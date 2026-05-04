"""SAFS mesh driver — wraps gmsh CLI, computes domain box, writes side-cars.

Implements PLAN.md Phase 5 + PLAN_origin.md Phase 3 + PLAN_domain.md Phase 1.

Reads:
    transform.json   — origin (validated against ts_to_stl.py's run)
    bbox.json        — per-fault local-frame bounding boxes

Writes:
    safs_includes.geo  — one Merge per included fault (consumed by safs.geo)
    domain_box.json    — resolved bulk-volume extents
    sizing.json        — resolved size-field parameters + measured tet count
    <output>.msh       — the tagged Gmsh mesh
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

from safs_origin import ORIGIN_UTM, read_transform_json


def _filter_bbox(per_fault_bbox: dict, included: list[str]) -> dict:
    missing = sorted(set(included) - set(per_fault_bbox))
    if missing:
        raise SystemExit(
            f"--include-fault asked for {missing} but bbox.json has only "
            f"{sorted(per_fault_bbox)}; rerun ts_to_stl.py with the same "
            f"--include-fault list"
        )
    return {k: v for k, v in per_fault_bbox.items() if k in included}


def _compute_box(per_fault_bbox: dict,
                 buf_x: float, buf_y: float,
                 depth: float, round_to: float) -> dict:
    """Compute bulk box extents from per-fault local-frame bboxes.

    Box snaps to round_to grid; floor/ceil ensures the resulting box always
    contains the buffered fault footprint.  See PLAN_domain.md Phase 1
    step 2 for the exact rule.
    """
    if not per_fault_bbox:
        raise ValueError("no faults to compute box from")
    x_fault_min = min(b["xmin"] for b in per_fault_bbox.values())
    x_fault_max = max(b["xmax"] for b in per_fault_bbox.values())
    y_fault_min = min(b["ymin"] for b in per_fault_bbox.values())
    y_fault_max = max(b["ymax"] for b in per_fault_bbox.values())
    z_fault_min = min(b["zmin"] for b in per_fault_bbox.values())
    z_fault_max = max(b["zmax"] for b in per_fault_bbox.values())

    if round_to <= 0:
        raise ValueError("--round-to must be > 0")
    x_min = math.floor((x_fault_min - buf_x) / round_to) * round_to
    x_max = math.ceil ((x_fault_max + buf_x) / round_to) * round_to
    y_min = math.floor((y_fault_min - buf_y) / round_to) * round_to
    y_max = math.ceil ((y_fault_max + buf_y) / round_to) * round_to
    z_top = 0.0
    z_bot = -float(depth)

    return {
        "x_min_m": x_min, "x_max_m": x_max,
        "y_min_m": y_min, "y_max_m": y_max,
        "z_top_m": z_top, "z_bot_m": z_bot,
        "buf_x_m": float(buf_x), "buf_y_m": float(buf_y),
        "depth_m": float(depth), "round_to_m": float(round_to),
        "Lx_m": x_max - x_min, "Ly_m": y_max - y_min,
        "Lz_m": z_top - z_bot,
        "fault_xmin_m": x_fault_min, "fault_xmax_m": x_fault_max,
        "fault_ymin_m": y_fault_min, "fault_ymax_m": y_fault_max,
        "fault_zmin_m": z_fault_min, "fault_zmax_m": z_fault_max,
    }


def _validate_clearance(box: dict, per_fault_bbox: dict,
                         buf_x: float, buf_y: float,
                         depth: float) -> None:
    """Hard-fail if any fault is too close to a box face.

    Per PLAN_domain.md Phase 1 step 4: clearance >= 0.1 * buffer on the
    relevant axis.  Fails with the offending fault's short-name.
    """
    margin_x = 0.1 * buf_x
    margin_y = 0.1 * buf_y
    margin_z = 0.1 * depth
    for short, b in per_fault_bbox.items():
        if b["xmin"] - box["x_min_m"] < margin_x:
            raise SystemExit(
                f"fault {short} is {b['xmin'] - box['x_min_m']:.0f} m from "
                f"x_min boundary; need >= {margin_x:.0f} m. "
                f"Increase --buf-x."
            )
        if box["x_max_m"] - b["xmax"] < margin_x:
            raise SystemExit(
                f"fault {short} is {box['x_max_m'] - b['xmax']:.0f} m from "
                f"x_max boundary; need >= {margin_x:.0f} m. "
                f"Increase --buf-x."
            )
        if b["ymin"] - box["y_min_m"] < margin_y:
            raise SystemExit(
                f"fault {short} is {b['ymin'] - box['y_min_m']:.0f} m from "
                f"y_min boundary; need >= {margin_y:.0f} m. "
                f"Increase --buf-y."
            )
        if box["y_max_m"] - b["ymax"] < margin_y:
            raise SystemExit(
                f"fault {short} is {box['y_max_m'] - b['ymax']:.0f} m from "
                f"y_max boundary; need >= {margin_y:.0f} m. "
                f"Increase --buf-y."
            )
        if b["zmin"] - box["z_bot_m"] < margin_z:
            raise SystemExit(
                f"fault {short} is {b['zmin'] - box['z_bot_m']:.0f} m from "
                f"z_bot boundary; need >= {margin_z:.0f} m. "
                f"Increase --depth."
            )


def _validate_sizing(res_f: float, res_ff: float, ramp_dist: float) -> None:
    if res_f <= 0 or res_ff <= 0 or ramp_dist <= 0:
        raise SystemExit("--res-f, --res-ff, --ramp-dist must all be > 0")
    if res_f >= res_ff:
        raise SystemExit(
            f"size-field inversion: --res-f ({res_f}) must be < --res-ff "
            f"({res_ff})"
        )
    if res_f < ramp_dist / 50.0:
        raise SystemExit(
            f"--res-f ({res_f}) is too fine relative to --ramp-dist "
            f"({ramp_dist}); rule: res_f >= ramp_dist / 50"
        )
    if res_ff > ramp_dist * 1.5:
        raise SystemExit(
            f"--res-ff ({res_ff}) is too coarse relative to --ramp-dist "
            f"({ramp_dist}); rule: res_ff <= ramp_dist * 1.5"
        )


def _read_ascii_stl(path: Path) -> tuple[list[tuple[float, float, float]],
                                            list[tuple[int, int, int]]]:
    """Read an ASCII STL.  Returns (vertices, triangles) with per-triangle
    vertex indices into the local vertex list."""
    verts: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []
    cur: list[tuple[float, float, float]] = []
    with open(path, "r") as fh:
        for line in fh:
            t = line.strip().split()
            if len(t) >= 4 and t[0] == "vertex":
                cur.append((float(t[1]), float(t[2]), float(t[3])))
            elif t and t[0] == "endloop":
                if len(cur) != 3:
                    raise ValueError(
                        f"{path}: outer loop with {len(cur)} vertices"
                    )
                a = len(verts); b = a + 1; c = a + 2
                verts.extend(cur)
                tris.append((a, b, c))
                cur = []
    return verts, tris


def _extract_polyline_endpoints(included: list[str], stl_dir: Path,
                                  snap_m: float = 0.1) -> list[tuple[float, float, float]]:
    """For each pair of fault STLs, find vertex coords shared between
    them at `snap_m` precision.  Returns a list of unique 3-D coords
    (each appearing in 2+ fault STLs) — these are the cross-fault
    polyline endpoints that delimit every shared boundary curve.

    Used by Option L (sliver-fix, 2026-04-30) to define the locations
    where local mesh refinement is imposed.

    A vertex coord is considered "shared" when, after rounding to
    `snap_m` grid, it appears in at least 2 different fault STLs.
    For a polyline of length L between faults A and B with K segments,
    K+1 endpoints (= K+1 shared vertices) will be returned.
    """
    import numpy as np
    snap = float(snap_m)
    if snap <= 0.0:
        snap = 0.1
    # vertex-key -> set of fault indices
    key_to_faults: dict[tuple[int, int, int], set[int]] = {}
    for k, short in enumerate(included):
        p = stl_dir / f"{short}.stl"
        if not p.exists():
            raise SystemExit(f"missing STL: {p}")
        verts, _ = _read_ascii_stl(p)
        v_arr = np.asarray(verts, dtype=np.float64)
        snapped = np.round(v_arr / snap).astype(np.int64)
        for row in snapped:
            key = (int(row[0]), int(row[1]), int(row[2]))
            key_to_faults.setdefault(key, set()).add(k)
    # Recover float coords from a representative vertex per shared key.
    polyline_pts: list[tuple[float, float, float]] = []
    for key, fset in key_to_faults.items():
        if len(fset) >= 2:
            polyline_pts.append((key[0] * snap, key[1] * snap, key[2] * snap))
    return polyline_pts


def _write_polyline_endpoints_geo(out_path: Path,
                                   pts: list[tuple[float, float, float]],
                                   default_size: float = 1.0e9) -> None:
    """Write a gmsh .geo fragment defining `pl_pts[]` as a list of
    Point IDs corresponding to every cross-fault polyline endpoint.
    Included by safs.geo when `local_refine = 1`.

    The default Point characteristic length is set to 1e9 m (huge)
    so that EVEN IF gmsh's mesher inserts these Points as Steiner
    vertices in the bulk volume (which it can do despite
    `Mesh.MeshSizeFromPoints = 0`), the Point-driven local mesh size
    is effectively unbounded and won't densify the mesh.  The actual
    sizing comes from Field[3]=Distance + Field[4]=Threshold, which
    impose `local_refine_size` within `local_refine_radius`.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("// Auto-generated by generate_safs_mesh.py:_write_polyline_endpoints_geo.\n")
        fh.write(f"// {len(pts)} cross-fault polyline endpoint(s).\n")
        fh.write("// DO NOT edit by hand — regenerated each run when --local-refine.\n")
        fh.write("\n")
        if not pts:
            # Empty list — write a sentinel so safs.geo's Field[3] doesn't
            # error.  Field with empty PointsList has degenerate Distance
            # = +inf everywhere, so Field[4] uses SizeMax = res_f → no
            # local refinement.  This is correct for the no-shared-vertex
            # case (single fault or fully disjoint set).
            fh.write("pl_pts[] = {};\n")
            return
        fh.write("pl_pts[] = {};\n")
        for x, y, z in pts:
            # gmsh idiom: capture newp into a temporary, then define
            # the Point at that ID and append it to pl_pts.
            fh.write(f"_pid = newp;\n")
            fh.write(f"Point(_pid) = {{{x:+.17e}, {y:+.17e}, {z:+.17e}, {default_size}}};\n")
            fh.write(f"pl_pts[] += _pid;\n")


def _combine_stls(out_path: Path, included: list[str],
                   stl_dir: Path, snap_m: float = 0.01) -> dict:
    """Combine per-fault ASCII STLs into a single ASCII STL with vertex
    deduplication at `snap_m` precision so shared CFM vertices at
    fault-fault junctions become a single point in the combined mesh.

    Without this, Gmsh's STL Merge sees two faults each carrying a copy
    of the shared vertex; HXT then "filters" the duplicates and the
    constraint topology breaks.

    Returns metadata dict (per-fault triangle ranges, total counts).
    """
    import numpy as np
    snap = snap_m
    all_verts: list[tuple[int, int, int]] = []
    all_tris: list[tuple[int, int, int]] = []
    fault_ranges: dict[str, list[int]] = {}
    # R-006: read each per-fault STL exactly once and cache.  Two reads of
    # the same file with cursor arithmetic between them is a silent-
    # corruption hazard if the read returns even one different vertex.
    cache: dict[str, tuple[list, list]] = {}
    for short in included:
        p = stl_dir / f"{short}.stl"
        if not p.exists():
            raise SystemExit(f"missing STL: {p}")
        verts, tris = _read_ascii_stl(p)
        cache[short] = (verts, tris)
        # Snap each vertex to the global snap grid; remap triangle indices.
        v_arr = np.asarray(verts, dtype=np.float64)
        snapped = np.round(v_arr / snap).astype(np.int64)
        snap_keys = [tuple(int(x) for x in row) for row in snapped]
        # local index -> global snapped key
        n_start = len(all_tris)
        # Append, then dedup at the end via a single np.unique pass.
        offset = len(all_verts)
        all_verts.extend(snap_keys)
        for a, b, c in tris:
            all_tris.append((a + offset, b + offset, c + offset))
        fault_ranges[short] = [n_start, len(all_tris)]
    arr = np.asarray(all_verts, dtype=np.int64)
    unique_keys, inverse = np.unique(arr, axis=0, return_inverse=True)
    # For each unique snapped key, recover the float coordinate from the
    # FIRST occurrence (preserves sub-snap precision of the surviving
    # representative).  Reuse the cached reads (R-006).
    raw_coords: dict[int, tuple[float, float, float]] = {}
    cursor = 0
    for short in included:
        verts, _ = cache[short]
        for v in verts:
            r = int(inverse[cursor])
            if r not in raw_coords:
                raw_coords[r] = v
            cursor += 1
    # Rewrite triangles with global unique vertex indices.
    new_tris_pre = [tuple(int(inverse[i]) for i in tri) for tri in all_tris]
    # R-007: drop degenerate triangles WHILE remapping fault_ranges so the
    # ranges still index into the post-drop triangle list.  Also drop
    # DUPLICATE triangles — these can arise when two faults' meshes
    # share a polyline and gmsh's Merge welds vertices tightly enough
    # to make two triangles (one from each fault) end up with the
    # same 3 vertex IDs.  HXT rejects such duplicates with
    # "Found two duplicated facets".
    new_tris: list[tuple[int, int, int]] = []
    range_remap: dict[str, list[int]] = {s: [0, 0] for s in fault_ranges}
    n_dropped_degenerate = 0
    n_dropped_duplicate = 0
    seen_keys: set[frozenset[int]] = set()
    for short, (start, end) in fault_ranges.items():
        out_start = len(new_tris)
        for k in range(start, end):
            a, b, c = new_tris_pre[k]
            if a == b or b == c or a == c:
                n_dropped_degenerate += 1
                continue
            key = frozenset((a, b, c))
            if key in seen_keys:
                n_dropped_duplicate += 1
                continue
            seen_keys.add(key)
            new_tris.append((a, b, c))
        out_end = len(new_tris)
        range_remap[short] = [out_start, out_end]
    fault_ranges = range_remap

    # Write combined ASCII STL.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("solid SAFS:combined\n")
        fmt = "+.17e"
        for a, b, c in new_tris:
            va = raw_coords[a]; vb = raw_coords[b]; vc = raw_coords[c]
            ax, ay, az = va; bx, by, bz = vb; cx, cy, cz = vc
            nx_ = (by - ay) * (cz - az) - (bz - az) * (cy - ay)
            ny_ = (bz - az) * (cx - ax) - (bx - ax) * (cz - az)
            nz_ = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
            nrm = (nx_ * nx_ + ny_ * ny_ + nz_ * nz_) ** 0.5
            if nrm == 0.0:
                continue
            nx_ /= nrm; ny_ /= nrm; nz_ /= nrm
            fh.write(
                f"  facet normal {nx_:{fmt}} {ny_:{fmt}} {nz_:{fmt}}\n"
                "    outer loop\n"
                f"      vertex {ax:{fmt}} {ay:{fmt}} {az:{fmt}}\n"
                f"      vertex {bx:{fmt}} {by:{fmt}} {bz:{fmt}}\n"
                f"      vertex {cx:{fmt}} {cy:{fmt}} {cz:{fmt}}\n"
                "    endloop\n  endfacet\n"
            )
        fh.write("endsolid SAFS:combined\n")
    return {
        "n_unique_vertices": int(unique_keys.shape[0]),
        "n_triangles": int(len(new_tris)),
        "n_dropped_degenerate": int(n_dropped_degenerate),
        "n_dropped_duplicate": int(n_dropped_duplicate),
        "fault_ranges": fault_ranges,
        "snap_m": snap_m,
    }


def _write_includes(path: Path, included: list[str], stl_dir: Path,
                     combined_stl: Path | None = None) -> None:
    # R-008: guard against silent override when caller mismatches args.
    if combined_stl is not None and len(included) < 2:
        raise ValueError(
            f"combined_stl was supplied but only {len(included)} fault(s) "
            f"included; the combine step is unnecessary at this scale"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    missing = []
    for short in included:
        p = stl_dir / f"{short}.stl"
        if not p.exists():
            missing.append(p)
    if missing:
        raise SystemExit(
            f"missing STL file(s): {missing}; "
            f"rerun ts_to_stl.py --include-fault ..."
        )
    with open(path, "w") as fh:
        fh.write("// Auto-generated by generate_safs_mesh.py — do not edit.\n")
        if combined_stl is not None:
            # Single Merge of the pre-deduplicated combined STL (multi-fault
            # path).  Avoids gmsh's "duplicate point filter" eating shared
            # fault-fault junction vertices.
            fh.write(f'Merge "{combined_stl.resolve()}";\n')
        else:
            for short in included:
                fh.write(f'Merge "{stl_dir.resolve() / (short + ".stl")}";\n')


def _measured_counts(msh_path: Path) -> dict:
    """Count tets/triangles per physical tag using meshio."""
    import meshio
    import numpy as np
    m = meshio.read(msh_path)
    out: dict = {"n_points": int(m.points.shape[0])}
    gp = m.cell_data_dict.get("gmsh:physical", {})
    for celltype, tags in gp.items():
        u, c = np.unique(tags, return_counts=True)
        for tt, cc in zip(u, c):
            out[f"{celltype}_tag_{int(tt)}"] = int(cc)
    return out


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate the SAFS Gmsh mesh from cleaned per-fault STLs.",
    )
    parser.add_argument("--include-fault", action="append", default=None,
                        metavar="SHORT_NAME")
    parser.add_argument("--stl-dir", required=True, type=Path)
    parser.add_argument("--bbox-json", required=True, type=Path)
    parser.add_argument("--transform-json", required=True, type=Path)
    parser.add_argument("--geo-file", default=None, type=Path,
                        help="Path to safs.geo; defaults to the script's "
                             "sibling.")
    parser.add_argument("--includes-geo", default=None, type=Path,
                        help="Where to write the auto-generated "
                             "safs_includes.geo; defaults to "
                             "<geo-file's directory>/safs_includes.geo.")
    parser.add_argument("--buf-x", type=float, default=50_000.0,
                        metavar="METRES")
    parser.add_argument("--buf-y", type=float, default=50_000.0,
                        metavar="METRES")
    parser.add_argument("--depth", type=float, default=50_000.0,
                        metavar="METRES")
    parser.add_argument("--round-to", type=float, default=10_000.0,
                        metavar="METRES")
    parser.add_argument("--res-f", type=float, default=1_000.0,
                        metavar="METRES",
                        help="On-fault target edge length (default 1000 m).")
    parser.add_argument("--res-ff", type=float, default=20_000.0,
                        metavar="METRES",
                        help="Far-field edge length (default 20 km).")
    parser.add_argument("--ramp-dist", type=float, default=30_000.0,
                        metavar="METRES",
                        help="Linear ramp width (default 30 km).")
    parser.add_argument("--tube-radius", type=float, default=5_000.0,
                        metavar="METRES",
                        help="Uniform near-fault tube radius — distance "
                             "from any fault inside which the bulk mesh is "
                             "uniform at --res-f (default 5 km).  Set to "
                             "0 to ramp immediately from --res-f to "
                             "--res-ff (legacy behaviour).")
    parser.add_argument("--algo3d", type=int, default=10, choices=(4, 10),
                        help="Mesh.Algorithm3D — 10 (default) = HXT, "
                             "enforces conformal fault interfaces.  Use "
                             "4 (Frontal-Delaunay) ONLY as a fallback for "
                             "non-conformal multi-fault input that HXT "
                             "rejects (validator check 5 will then fail).")
    parser.add_argument("--allow-non-conformal", action="store_true",
                        help="Skip the cross-fault crossing scan that "
                             "refuses non-conformal multi-fault input.  "
                             "Use this only when intentionally meshing a "
                             "raw (un-conformalized) multi-fault set with "
                             "Algorithm3D=4 as a fallback path.")
    parser.add_argument("--combine-snap-m", type=float, default=0.01,
                        metavar="METRES",
                        help="Vertex snap tolerance when combining "
                             "per-fault STLs into a single conformal "
                             "STL (multi-fault only).  Vertices closer "
                             "than this become a single shared vertex.  "
                             "Must exceed gmsh's internal STL merge "
                             "tolerance to take effect.  Default: 0.01 m.")
    parser.add_argument("--local-refine", action="store_true",
                        help="Option L (sliver-fix, 2026-04-30): write "
                             "polyline_endpoints.geo (vertex coords shared "
                             "by 2+ fault STLs) and pass `local_refine=1` "
                             "to safs.geo, which adds a Distance+Threshold "
                             "field that imposes a smaller mesh size near "
                             "every cross-fault polyline endpoint.  Helps "
                             "the worst-γ slivers concentrated at fault-"
                             "fault intersections.  Cost: tet count "
                             "typically grows 1.3-2x.")
    parser.add_argument("--local-refine-size-m", type=float, default=500.0,
                        help="Min mesh size (m) at polyline endpoints "
                             "when --local-refine is on.  Default 500.")
    parser.add_argument("--local-refine-radius-m", type=float, default=2000.0,
                        help="Radius (m) of local refinement zone around "
                             "each polyline endpoint.  Default 2000.")
    parser.add_argument("--origin-easting", type=float, default=ORIGIN_UTM[0])
    parser.add_argument("--origin-northing", type=float, default=ORIGIN_UTM[1])
    parser.add_argument("--origin-elevation", type=float, default=ORIGIN_UTM[2])
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Output .msh path.")
    args = parser.parse_args(argv)

    # 1. Load metadata
    tx = read_transform_json(args.transform_json)
    bbox = json.loads(args.bbox_json.read_text())

    # 2. Cross-check origin
    flag_origin = (args.origin_easting, args.origin_northing,
                   args.origin_elevation)
    file_origin = tuple(tx["origin_utm_m"])
    if any(abs(a - b) > 1e-9 for a, b in zip(flag_origin, file_origin)):
        raise SystemExit(
            f"origin mismatch: transform.json says {file_origin}, CLI flags "
            f"say {flag_origin}; refusing to mesh against inconsistent "
            f"metadata"
        )

    # 3. Filter to --include-fault subset (default: every fault in bbox.json)
    if args.include_fault:
        bbox = _filter_bbox(bbox, args.include_fault)
        included = list(args.include_fault)
    else:
        included = sorted(bbox)

    # 4. Validate sizing-field consistency
    _validate_sizing(args.res_f, args.res_ff, args.ramp_dist)
    # Tube-radius validation (PLAN Phase 3 §Edge cases).
    if args.tube_radius < 0.0:
        raise SystemExit(f"--tube-radius must be >= 0; got {args.tube_radius}")
    if args.tube_radius >= args.buf_x or args.tube_radius >= args.buf_y:
        raise SystemExit(
            f"--tube-radius {args.tube_radius} m engulfs the box buffer "
            f"(buf_x={args.buf_x} m, buf_y={args.buf_y} m); there would "
            f"be no far-field.  Reduce --tube-radius or increase "
            f"--buf-x/--buf-y."
        )

    # 5. Compute domain box
    box = _compute_box(bbox, args.buf_x, args.buf_y,
                       args.depth, args.round_to)
    _validate_clearance(box, bbox, args.buf_x, args.buf_y, args.depth)

    # 5b. PLAN Phase 4 — conformal-input detection.  When --stl-dir
    # contains the Phase 2 conformalize output (signaled by
    # ``triangle_to_fault.json``), skip the cross-fault crossing scan
    # — the input is already conformal by construction.  Otherwise,
    # for multi-fault input, run the scan and refuse non-conformal
    # input (unless --allow-non-conformal is set).
    has_t2f = (args.stl_dir / "triangle_to_fault.json").exists()
    is_conformal_input = has_t2f
    if len(included) > 1 and not is_conformal_input \
            and not args.allow_non_conformal:
        try:
            from fault_intersect import scan_cross_fault_crossings
        except ImportError as exc:
            raise SystemExit(
                f"cannot import fault_intersect to gate non-conformal "
                f"multi-fault input: {exc}"
            )
        clearance_m = float(tx.get("free_surface_clearance_m", 0.0))
        crossings = scan_cross_fault_crossings(
            args.stl_dir, included, clearance_m=clearance_m,
        )
        if crossings:
            n_pairs = sum(len(v) for v in crossings.values())
            details = "\n  ".join(
                f"{a} × {b}: {len(plines)} polyline(s)"
                for (a, b), plines in sorted(crossings.items())
            )
            raise SystemExit(
                f"refusing to mesh: {len(crossings)} cross-fault "
                f"intersection cluster(s), {n_pairs} polyline(s):\n  "
                f"{details}\n"
                f"Run conformalize_faults.py first to insert the "
                f"intersection curves into both triangulations, OR "
                f"split into disjoint subsets, OR pass "
                f"--allow-non-conformal to mesh anyway with "
                f"--algo3d=4 (validator check 5 will fail)."
            )

    # 6. Resolve geo + includes paths
    script_dir = Path(__file__).resolve().parent
    geo_file = args.geo_file or (script_dir / "safs.geo")
    if not geo_file.exists():
        raise SystemExit(f"safs.geo not found at {geo_file}")
    includes_path = args.includes_geo or (geo_file.parent / "safs_includes.geo")
    # When more than one fault is included, pre-merge the per-fault STLs
    # into a single combined STL with shared fault-fault junction vertices
    # unified.  Otherwise gmsh's STL merge filters them and HXT cannot
    # recover the constraint topology.
    if len(included) > 1:
        # R-005: combined STL is a pipeline artefact, not user input;
        # write it next to the .msh, not into the input STL dir.
        args.output.parent.mkdir(parents=True, exist_ok=True)
        combined_stl = args.output.parent / "safs_combined.stl"
        meta = _combine_stls(combined_stl, included, args.stl_dir,
                              snap_m=args.combine_snap_m)
        (args.output.parent / "combine_meta.json").write_text(
            json.dumps(meta, indent=2) + "\n"
        )
        print(f"  combined STL  : {combined_stl} "
              f"({meta['n_unique_vertices']} verts, "
              f"{meta['n_triangles']} tris, "
              f"snap={meta['snap_m']} m)")
        _write_includes(includes_path, included, args.stl_dir,
                         combined_stl=combined_stl)
    else:
        _write_includes(includes_path, included, args.stl_dir)

    # If the geo lives elsewhere than the includes path, copy it next to the
    # includes so the relative `Include "safs_includes.geo";` resolves.
    if includes_path.parent != geo_file.parent:
        target_geo = includes_path.parent / geo_file.name
        shutil.copy2(geo_file, target_geo)
        run_geo = target_geo
    else:
        run_geo = geo_file

    # Option L: extract cross-fault polyline endpoints and write
    # polyline_endpoints.geo next to the run_geo.  safs.geo includes
    # this file when local_refine=1.  Always write the file (with empty
    # list when no shared vertices) so the Include directive doesn't
    # fail; the gmsh-side conditional gates whether the field is used.
    polyline_geo_path = run_geo.parent / "polyline_endpoints.geo"
    if args.local_refine and len(included) >= 2:
        pl_pts = _extract_polyline_endpoints(
            included, args.stl_dir,
            snap_m=max(args.combine_snap_m, 0.1))
        _write_polyline_endpoints_geo(polyline_geo_path, pl_pts)
        print(f"  polyline endpoints: {len(pl_pts)} "
              f"(written to {polyline_geo_path})")
    else:
        # Write an empty list so the Include doesn't fail when local_refine=0.
        _write_polyline_endpoints_geo(polyline_geo_path, [])

    # 7. Side-car: domain_box.json
    out_dir = args.output.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "domain_box.json").write_text(
        json.dumps(box, indent=2) + "\n"
    )

    # 8. Invoke gmsh
    cmd = [
        "gmsh", "-3", str(run_geo),
        "-setnumber", "x_min", repr(box["x_min_m"]),
        "-setnumber", "x_max", repr(box["x_max_m"]),
        "-setnumber", "y_min", repr(box["y_min_m"]),
        "-setnumber", "y_max", repr(box["y_max_m"]),
        "-setnumber", "z_top", repr(box["z_top_m"]),
        "-setnumber", "z_bot", repr(box["z_bot_m"]),
        "-setnumber", "res_f", repr(args.res_f),
        "-setnumber", "res_ff", repr(args.res_ff),
        "-setnumber", "ramp_dist", repr(args.ramp_dist),
        "-setnumber", "tube_radius", repr(args.tube_radius),
        "-setnumber", "algo3d", repr(int(args.algo3d)),
        "-setnumber", "local_refine", "1" if args.local_refine else "0",
        "-setnumber", "local_refine_size", repr(args.local_refine_size_m),
        "-setnumber", "local_refine_radius", repr(args.local_refine_radius_m),
        "-o", str(args.output),
    ]
    print(" ".join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(f"gmsh failed with exit code {r.returncode}")

    # 9. Side-car: sizing.json (with measured counts)
    counts = _measured_counts(args.output)
    sizing = {
        "res_f_m": float(args.res_f),
        "res_ff_m": float(args.res_ff),
        "ramp_dist_m": float(args.ramp_dist),
        "tube_radius_m": float(args.tube_radius),
        "algo3d": int(args.algo3d),
        "conformal_input": bool(is_conformal_input),
        "fields": {"distance_id": 1, "threshold_id": 2,
                   "background_field_id": 2,
                   "distance_sampling_per_surface": 100},
        "measured": counts,
    }
    (out_dir / "sizing.json").write_text(json.dumps(sizing, indent=2) + "\n")

    print(f"\ngenerate_safs_mesh OK: wrote {args.output}")
    print(f"  domain_box.json: {out_dir/'domain_box.json'}")
    print(f"  sizing.json    : {out_dir/'sizing.json'}")
    print(f"  measured       : {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
