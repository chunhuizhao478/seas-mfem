"""TSurf → cleaned STL converter.

Implements PLAN.md Phase 2 + PLAN_origin.md Phase 2.

Pipeline per fault (.ts at the chosen --res):
    parse → utm_to_local → clamp z>(-Z0) → collapse 1mm duplicates →
    drop zero-area triangles → write binary STL

Side-car outputs (one per run, written into --out-dir's parent):
    transform.json      — origin / convention metadata (PLAN_origin.md schema)
    bbox.json           — per-fault local-frame bounding box, used by the
                          generate_safs_mesh.py driver
    cleanup_log_<res>m.csv  — per-fault cleanup metrics
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

from audit_ts_quality import find_overlap_pairs, parse_tsurf
from safs_origin import (
    FAULT_SHORT_NAMES,
    ORIGIN_UTM,
    short_name_from_path,
    utm_to_local,
    write_transform_json,
)


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
def _clamp_freesurface(V: np.ndarray, Z0: float,
                        clearance: float = 0.0) -> tuple[np.ndarray, int]:
    """Clamp local z above the free-surface plane.

    Threshold is `-Z0 - clearance`.  With the production origin (Z0 = 0) and
    `clearance = 1.0`, vertices with z > -1 m get pinned to z = -1 m.  A
    non-zero clearance is required when the bulk-volume free surface lies
    at z = -Z0 (i.e., z = 0 with default origin) to keep fault triangles
    strictly inside the bulk volume — Gmsh/HXT's PLC recovery rejects fault
    edges that lie exactly on the box's top face.
    """
    threshold = -Z0 - clearance
    mask = V[:, 2] > threshold
    n = int(np.sum(mask))
    if n == 0:
        return V, 0
    V = V.copy()
    V[mask, 2] = threshold
    return V, n


def _dedup_vertices(V: np.ndarray, T: np.ndarray
                    ) -> tuple[np.ndarray, np.ndarray, int]:
    """Snap to 1 mm grid and collapse coincident vertices.

    Returns (V_unique, T_remapped, n_collapsed).  Triangles whose three
    indices reduce to fewer than 3 unique values after the collapse are
    dropped — the count of those drops is reported separately by the caller.
    """
    n_in = V.shape[0]
    grid = np.round(V * 1000.0).astype(np.int64)
    # np.unique with return_inverse gives row → unique-row map.
    _, inverse, _ = np.unique(grid, axis=0, return_inverse=True, return_counts=True)
    n_out = int(inverse.max()) + 1
    # Use the *first* original vertex's float coords for each unique row to
    # preserve sub-mm precision (the snap is only for matching, not for the
    # output coordinate).
    V_out = np.empty((n_out, 3), dtype=np.float64)
    seen = np.zeros(n_out, dtype=bool)
    for i in range(n_in):
        r = int(inverse[i])
        if not seen[r]:
            V_out[r] = V[i]
            seen[r] = True
    T_out = inverse[T]
    n_collapsed = n_in - n_out
    return V_out, T_out, n_collapsed


def _repair_overlaps(V: np.ndarray, T: np.ndarray
                      ) -> tuple[np.ndarray, int, list[tuple[int, int]]]:
    """Drop the second triangle of each coplanar same-side overlap pair.

    Two triangles that share an edge AND lie in the same plane AND have
    their third vertices on the same side of that edge in the shared
    plane overlap each other on the surface — HXT rejects this with
    "exactly self-intersecting facets, dihedral=0".  The CFM TSurf
    decimator occasionally produces a handful of such pairs; we drop
    one triangle per pair.

    Returns (T_kept, n_dropped, dropped_pairs).
    The first element of each pair survives; the second is dropped.
    The function iterates the pair list and skips already-dropped
    triangles so cascades from a single triangle in multiple pairs are
    handled correctly.
    """
    pairs = find_overlap_pairs(V, T)
    if not pairs:
        return T, 0, []
    drop: set[int] = set()
    pairs_recorded: list[tuple[int, int]] = []
    for i, j in pairs:
        if i in drop or j in drop:
            continue
        drop.add(j)
        pairs_recorded.append((i, j))
    keep_mask = np.ones(T.shape[0], dtype=bool)
    for d in drop:
        keep_mask[d] = False
    return T[keep_mask], int(np.sum(~keep_mask)), pairs_recorded


def _drop_degenerate(V: np.ndarray, T: np.ndarray
                     ) -> tuple[np.ndarray, int, int]:
    """Drop triangles that are degenerate after dedup or by zero area.

    Returns (T_out, n_index_degenerate, n_zero_area).
    """
    # Index-degenerate: any two of the three indices equal.
    deg_idx = (T[:, 0] == T[:, 1]) | (T[:, 1] == T[:, 2]) | (T[:, 0] == T[:, 2])
    n_idx = int(np.sum(deg_idx))
    T = T[~deg_idx]

    # Zero area.
    a = V[T[:, 0]]
    b = V[T[:, 1]]
    c = V[T[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    keep = areas >= 1e-6
    n_zero = int(np.sum(~keep))
    T = T[keep]
    return T, n_idx, n_zero


# ---------------------------------------------------------------------------
# Binary STL writer
# ---------------------------------------------------------------------------
def _write_ascii_stl(path: Path, V: np.ndarray, T: np.ndarray,
                      short_name: str) -> None:
    """Write an ASCII STL with full double-precision coordinates.

    Binary STL stores coordinates as float32, which truncates ~1 mm at SAFS
    coordinate magnitudes (10^5 m) and triggers HXT constraint-recovery
    failures on dense fault triangulations.  ASCII STL with %.17e gives
    the implementer full double-precision round-trip.
    """
    n = T.shape[0]
    if n == 0:
        raise ValueError(f"refusing to write empty STL: {path}")
    a = V[T[:, 0]]
    b = V[T[:, 1]]
    c = V[T[:, 2]]
    normals = np.cross(b - a, c - a)
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    norms = np.where(norms > 0.0, norms, 1.0)
    normals = normals / norms

    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = "%.17e"
    with open(path, "w") as fh:
        fh.write(f"solid SAFS:{short_name}:cleaned\n")
        for i in range(n):
            nx, ny, nz = normals[i]
            ax, ay, az = a[i]
            bx, by, bz = b[i]
            cx, cy, cz = c[i]
            fh.write(
                f"  facet normal {nx:{'+.17e'}} {ny:{'+.17e'}} {nz:{'+.17e'}}\n"
                "    outer loop\n"
                f"      vertex {ax:{'+.17e'}} {ay:{'+.17e'}} {az:{'+.17e'}}\n"
                f"      vertex {bx:{'+.17e'}} {by:{'+.17e'}} {bz:{'+.17e'}}\n"
                f"      vertex {cx:{'+.17e'}} {cy:{'+.17e'}} {cz:{'+.17e'}}\n"
                "    endloop\n"
                "  endfacet\n"
            )
        fh.write(f"endsolid SAFS:{short_name}:cleaned\n")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# bbox.json
# ---------------------------------------------------------------------------
def _bbox_for(V: np.ndarray) -> dict:
    return {
        "xmin": float(V[:, 0].min()), "xmax": float(V[:, 0].max()),
        "ymin": float(V[:, 1].min()), "ymax": float(V[:, 1].max()),
        "zmin": float(V[:, 2].min()), "zmax": float(V[:, 2].max()),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
_LOG_COLUMNS = [
    "short_name",
    "input_file",
    "n_input_vrtx", "n_input_trgl",
    "n_clamped_z",
    "n_dup_collapsed",
    "n_index_degenerate_dropped",
    "n_zero_area_dropped",
    "n_overlap_pairs_dropped",
    "final_n_vrtx",
    "final_n_trgl",
    "final_xmin", "final_xmax",
    "final_ymin", "final_ymax",
    "final_zmin", "final_zmax",
    "stl_sha256",
    "stl_path",
]


def _select_files(cfm_dir: Path, res: int,
                  include_fault: list[str] | None) -> list[Path]:
    files = sorted(cfm_dir.glob(f"*_{res}m.ts"))
    if not files:
        raise FileNotFoundError(f"no *_{res}m.ts files under {cfm_dir}")
    keep: list[Path] = []
    for p in files:
        try:
            short = short_name_from_path(p)
        except KeyError:
            continue
        if include_fault and short not in include_fault:
            continue
        keep.append(p)
    if include_fault:
        found = {short_name_from_path(p) for p in keep}
        missing = sorted(set(include_fault) - found)
        if missing:
            raise FileNotFoundError(
                f"--include-fault asked for {missing} but no matching "
                f"*_{res}m.ts files were found under {cfm_dir}"
            )
    if not keep:
        raise FileNotFoundError(f"no files selected (filter={include_fault})")
    return keep


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert CFM .ts files to cleaned binary STL in the "
                    "SAFS local frame, plus transform.json + bbox.json.",
    )
    parser.add_argument("--cfm-dir", required=True, type=Path)
    parser.add_argument("--res", required=True, type=int,
                        choices=(500, 1000, 2000))
    parser.add_argument("--include-fault", action="append", default=None,
                        metavar="SHORT_NAME")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Directory to write STL files into.")
    parser.add_argument("--transform-json", default=None, type=Path,
                        help="transform.json path; defaults to "
                             "<out-dir>/../transform.json")
    parser.add_argument("--bbox-json", default=None, type=Path,
                        help="bbox.json path; defaults to "
                             "<out-dir>/../bbox.json")
    parser.add_argument("--log-csv", default=None, type=Path,
                        help="Cleanup log CSV path; defaults to "
                             "<out-dir>/../cleanup_log_<res>m.csv")
    parser.add_argument("--origin-easting", type=float, default=ORIGIN_UTM[0])
    parser.add_argument("--origin-northing", type=float, default=ORIGIN_UTM[1])
    parser.add_argument("--origin-elevation", type=float, default=ORIGIN_UTM[2])
    parser.add_argument("--free-surface-clearance", type=float, default=1.0,
                        metavar="METRES",
                        help="Clamp z above (-Z0 - clearance) to that value. "
                             "Required to be > 0 so fault edges sit strictly "
                             "below the bulk box's top face (HXT PLC needs "
                             "this clearance to avoid edge-edge intersections "
                             "between the fault trace and the box top).  "
                             "Default: 1.0 m.")
    parser.add_argument("--repair-overlaps", action="store_true",
                        help="Detect and drop coplanar same-side triangle "
                             "overlap pairs in the CFM source.  Required "
                             "for some CFM (resolution, fault) combinations "
                             "(e.g. Pinto Mountain at every resolution; "
                             "Mill Creek strand at 2000 m).  Default: off.")
    args = parser.parse_args(argv)

    origin = (args.origin_easting, args.origin_northing, args.origin_elevation)
    Z0 = args.origin_elevation
    clearance = float(args.free_surface_clearance)
    if clearance < 0.0:
        parser.error("--free-surface-clearance must be >= 0")

    # Resolve default side-car paths.
    parent = args.out_dir.parent
    transform_path = args.transform_json or (parent / "transform.json")
    bbox_path     = args.bbox_json     or (parent / "bbox.json")
    log_path      = args.log_csv       or (parent / f"cleanup_log_{args.res}m.csv")

    files = _select_files(args.cfm_dir, args.res, args.include_fault)

    log_rows: list[dict] = []
    bbox_data: dict[str, dict] = {}

    for p in files:
        short = short_name_from_path(p)
        ts = parse_tsurf(p)
        V_utm = ts.V
        T = ts.T

        n_in_v = int(V_utm.shape[0])
        n_in_t = int(T.shape[0])

        # 1. UTM → local
        V = utm_to_local(V_utm, origin=origin)

        # 2. Clamp z above (-Z0 - clearance)
        V, n_clamped = _clamp_freesurface(V, Z0, clearance=clearance)

        # 3. Collapse 1 mm duplicates (rebuild T)
        V, T, n_dup = _dedup_vertices(V, T)

        # 4. Drop degenerate triangles
        T, n_idx, n_zero = _drop_degenerate(V, T)

        # 4b. Repair coplanar same-side overlaps (HXT-fatal CFM defect).
        if args.repair_overlaps:
            T, n_overlap, _ = _repair_overlaps(V, T)
        else:
            n_overlap = 0

        if T.shape[0] == 0:
            raise RuntimeError(
                f"{short}: cleanup left zero triangles; aborting"
            )

        # 5. Drop unreferenced vertices so the STL doesn't carry dead points.
        used = np.unique(T)
        remap = -np.ones(V.shape[0], dtype=np.int64)
        remap[used] = np.arange(used.size, dtype=np.int64)
        V = V[used]
        T = remap[T]

        # 6. Write STL (ASCII, full double precision)
        stl_path = args.out_dir / f"{short}.stl"
        _write_ascii_stl(stl_path, V, T, short)

        # 7. Per-fault bbox + log
        bb = _bbox_for(V)
        bbox_data[short] = bb
        log_rows.append({
            "short_name": short,
            "input_file": p.name,
            "n_input_vrtx": n_in_v,
            "n_input_trgl": n_in_t,
            "n_clamped_z": n_clamped,
            "n_dup_collapsed": n_dup,
            "n_index_degenerate_dropped": n_idx,
            "n_zero_area_dropped": n_zero,
            "n_overlap_pairs_dropped": n_overlap,
            "final_n_vrtx": int(V.shape[0]),
            "final_n_trgl": int(T.shape[0]),
            "final_xmin": bb["xmin"], "final_xmax": bb["xmax"],
            "final_ymin": bb["ymin"], "final_ymax": bb["ymax"],
            "final_zmin": bb["zmin"], "final_zmax": bb["zmax"],
            "stl_sha256": _sha256(stl_path),
            "stl_path": str(stl_path.relative_to(parent)) if parent in stl_path.parents
                         else str(stl_path),
        })

    # Write side-cars.
    transform_path.parent.mkdir(parents=True, exist_ok=True)
    write_transform_json(
        transform_path,
        origin=origin,
        extra={
            "generated_by": "ts_to_stl.py",
            "input_resolution_m": args.res,
            "input_files": [str(p) for p in files],
            "stl_dir": str(args.out_dir),
            "fault_short_names": [short_name_from_path(p) for p in files],
            "free_surface_clearance_m": clearance,
            "repair_overlaps": bool(args.repair_overlaps),
        },
    )

    bbox_path.parent.mkdir(parents=True, exist_ok=True)
    import json
    with open(bbox_path, "w") as fh:
        json.dump(bbox_data, fh, indent=2)
        fh.write("\n")

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_LOG_COLUMNS)
        w.writeheader()
        for r in log_rows:
            w.writerow(r)

    print(f"ts_to_stl OK: {len(log_rows)} fault(s)")
    print(f"  STL dir       : {args.out_dir}")
    print(f"  transform.json: {transform_path}")
    print(f"  bbox.json     : {bbox_path}")
    print(f"  cleanup log   : {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
