#!/usr/bin/env python3
"""
clean_freesurface_mesh.py — Clip a fault surface at z = 0 (free surface)
and isotropically remesh, producing a clean ASCII STL ready for gmsh.

Two clipping methods are supported.  Both end the same way: an isotropic
explicit remesh at TARGET_LENGTH for triangle quality.

  --method clip  (default, recommended)
      Per-triangle linear interpolation at z = 0 (the algorithm in
      ts_to_stl.clip_triangle_at_z0).  For each input triangle:
        * 3 verts at/below 0  -> kept whole
        * 3 verts above       -> dropped
        * 1 vert below, 2 above -> emit 1 triangle (b, a', c') with
          a', c' = z=0 intersections on edges b->a, b->c
        * 2 verts below, 1 above -> emit 2 triangles tiling the surviving
          quadrilateral
      The new top trace lies exactly on z = 0 by construction; no
      snap-strip hack is needed.  Coincident clip-points along shared
      edges are welded by `meshing_remove_duplicate_vertices`.

  --method vdel  (legacy, brute-force)
      MeshLab-GUI-style:
        1. Select all vertices with z > 0  (cap).
        2. Select faces with any (z0,z1,z2) > 0.
        3. meshing_remove_selected_vertices_and_faces.
        4. Snap remaining "top strip" (z > -SNAP_STRIP) onto z = 0.
      Drops every triangle that touches the cap, even those with two
      perfectly good vertices below z = 0; the survivors form a jagged
      free edge that the snap-strip then flattens.  Kept for comparison
      with the `clip` output.

Default TARGET_LENGTH is inferred from the filename (e.g. ..._2000m.ts ->
2000.0).  Default SNAP_STRIP is 10 % of TARGET_LENGTH (vdel only).

Usage (single file):
    conda activate pythonenv
    python clean_freesurface_mesh.py INPUT OUTPUT
        [--method {clip,vdel}] [--target-length M]
        [--snap-strip M] [--iterations N] [--keep-reproject]

Batch:
    python clean_freesurface_mesh.py --batch [--method {clip,vdel}]
    # processes every *_unclipped.stl in ../data_preprocess/ and writes
    # cleaned STLs to ../data_cleanfreesurf/<basename>_clean_<method>.stl.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pymeshlab as ml

# Reuse the .ts parser and the per-triangle clipper from ts_to_stl.py.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ts_to_stl import parse_ts as parse_ts_file  # noqa: E402
from ts_to_stl import clip_triangle_at_z0, EPS  # noqa: E402


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------

def _load_verts_tris(path: Path):
    """Return (verts, tris) regardless of input format.

    `verts` is a dict id -> (x, y, z); `tris` is a list of (ia, ib, ic).
    For .ts files we reuse the ts_to_stl parser (1-indexed ids).
    For other formats meshio loads them and we use 0-indexed ids.
    """
    if path.suffix.lower() == ".ts":
        return parse_ts_file(path)
    import meshio
    m = meshio.read(str(path))
    verts = {i: (float(p[0]), float(p[1]), float(p[2]))
             for i, p in enumerate(m.points)}
    tris = []
    for cb in m.cells:
        if cb.type == "triangle":
            tris.extend((int(t[0]), int(t[1]), int(t[2]))
                        for t in cb.data)
    if not tris:
        raise ValueError(f"no triangle cells in {path}")
    return verts, tris


def _ts_to_meshset(path: Path) -> ml.MeshSet:
    """Load a .ts directly into a pymeshlab MeshSet (used by --method vdel
    when the input is a .ts; for .stl input vdel uses load_new_mesh)."""
    verts, tris = parse_ts_file(path)
    sorted_ids = sorted(verts.keys())
    id_map = {old: new for new, old in enumerate(sorted_ids)}
    V = np.array([verts[i] for i in sorted_ids], dtype=np.float64)
    F = []
    for a, b, c in tris:
        if a in id_map and b in id_map and c in id_map:
            F.append([id_map[a], id_map[b], id_map[c]])
    F = np.asarray(F, dtype=np.int32)
    ms = ml.MeshSet()
    ms.add_mesh(ml.Mesh(vertex_matrix=V, face_matrix=F), "raw_fault")
    return ms


def _infer_target_length(path: Path) -> float | None:
    """Pull the resolution from a filename like '..._2000m.ts' or
    '..._500m_unclipped.stl'."""
    m = re.search(r"_(\d+)m(?:_[A-Za-z0-9]+)?\.(ts|stl|obj|ply)$", path.name)
    return float(m.group(1)) if m else None


# ----------------------------------------------------------------------
# Method: clip (linear interpolation at z = 0)
# ----------------------------------------------------------------------

def clip_to_arrays(verts: dict, tris: list) -> tuple[np.ndarray,
                                                     np.ndarray,
                                                     dict]:
    """Apply ts_to_stl.clip_triangle_at_z0 to every triangle.

    Returns flat (V, F) numpy arrays with duplicated vertices along shared
    edges; the caller welds them via `meshing_remove_duplicate_vertices`.
    """
    out_pts: list[tuple[float, float, float]] = []
    out_faces: list[list[int]] = []
    n_kept = n_dropped = n_clip2 = n_clip1 = 0
    for ia, ib, ic in tris:
        if ia not in verts or ib not in verts or ic not in verts:
            continue
        result = clip_triangle_at_z0(verts[ia], verts[ib], verts[ic])
        if not result:
            n_dropped += 1
            continue
        z = [verts[i][2] for i in (ia, ib, ic)]
        n_above = sum(1 for zi in z if zi > EPS)
        if n_above == 0:
            n_kept += 1
        elif n_above == 1:
            n_clip2 += 1   # 1 above => 2 below; clipper emits 2 triangles
        elif n_above == 2:
            n_clip1 += 1   # 2 above => 1 below; clipper emits 1 triangle
        for tri in result:
            base = len(out_pts)
            out_pts.extend(tri)
            out_faces.append([base, base + 1, base + 2])
    V = np.asarray(out_pts, dtype=np.float64)
    F = np.asarray(out_faces, dtype=np.int32)
    stats = {
        "n_in": len(tris),
        "n_out": len(out_faces),
        "n_kept_whole": n_kept,
        "n_dropped": n_dropped,
        "n_clipped_2to1": n_clip1,
        "n_clipped_1to2": n_clip2,
    }
    return V, F, stats


def clean_freesurface_clip(in_path: Path, out_path: Path,
                           target_length: float,
                           eps_clip: float | None = None,
                           iterations: int = 3,
                           reproject: bool = False,
                           verbose: bool = True) -> None:
    if eps_clip is None:
        eps_clip = 0.10 * target_length
    verts, tris = _load_verts_tris(in_path)
    if verbose:
        n_above = sum(1 for v in verts.values() if v[2] > EPS)
        print(f"loaded {in_path.name}: {len(verts)} verts, {len(tris)} tris "
              f"({n_above} verts above z=0)")

    # Pre-clip snap: any input vertex within eps_clip of z = 0 is moved
    # exactly onto the plane *before* clipping.  Without this, a vertex at
    # (e.g.) z = -3 m + two neighbours above z = 0 produces a near-degenerate
    # triangle whose remaining "below" corner is essentially on z = 0; that
    # sliver propagates into the fault triangulation and forces gmsh to wrap
    # it with low-quality tets.  Snapping converts those near-zero vertices
    # into clean boundary vertices that the clipper keeps without producing
    # slivers.
    if eps_clip > 0:
        n_snapped = 0
        for k in verts:
            x, y, z = verts[k]
            if -eps_clip < z < eps_clip:
                verts[k] = (x, y, 0.0)
                n_snapped += 1
        if verbose:
            print(f"  pre-clip snap: |z| < {eps_clip} m -> z = 0  "
                  f"({n_snapped} verts)")

    V, F, stats = clip_to_arrays(verts, tris)
    if verbose:
        print(f"  clip: in={stats['n_in']}, out={stats['n_out']} "
              f"(kept_whole={stats['n_kept_whole']}, "
              f"dropped={stats['n_dropped']}, "
              f"clipped 2->1={stats['n_clipped_2to1']}, "
              f"clipped 1->2={stats['n_clipped_1to2']})")

    ms = ml.MeshSet()
    ms.add_mesh(ml.Mesh(vertex_matrix=V, face_matrix=F), "clipped_fault")

    # Weld coincident clip-points and drop any null faces.
    ms.meshing_remove_duplicate_vertices()
    ms.meshing_remove_null_faces()
    ms.meshing_remove_unreferenced_vertices()
    try:
        ms.meshing_repair_non_manifold_edges(method="Remove Faces")
    except Exception as exc:
        if verbose:
            print(f"  (repair_non_manifold_edges skipped: {exc})",
                  file=sys.stderr)

    m = ms.current_mesh()
    if verbose:
        print(f"  after weld:    {m.vertex_number()} verts, "
              f"{m.face_number()} faces")

    if verbose:
        print(f"  isotropic remesh: targetlen={target_length} m, "
              f"iters={iterations}, reproject={reproject}")
    ms.meshing_isotropic_explicit_remeshing(
        iterations=iterations,
        adaptive=False,
        selectedonly=False,
        targetlen=ml.PureValue(target_length),
        featuredeg=30.0,
        checksurfdist=False,
        splitflag=True,
        collapseflag=True,
        swapflag=True,
        smoothflag=True,
        reprojectflag=reproject,
    )
    m = ms.current_mesh()
    if verbose:
        print(f"  after remesh:  {m.vertex_number()} verts, "
              f"{m.face_number()} faces")

    ms.meshing_remove_duplicate_vertices()
    ms.meshing_remove_unreferenced_vertices()
    ms.meshing_remove_null_faces()

    # Defensive snap: reproject smoothing in the remesh can drift the top
    # trace by < 1 m.  Pull anything within 1 m of the plane back onto z = 0.
    drift_tol = 1.0
    cond = f"((z > -{drift_tol}) && (z < {drift_tol}))"
    ms.compute_selection_by_condition_per_vertex(condselect=cond)
    n_drift = ms.current_mesh().selected_vertex_number()
    if n_drift:
        ms.compute_coord_by_function(x="x", y="y", z="0", onselected=True)
        if verbose:
            print(f"  drift snap: {n_drift} verts within {drift_tol} m of "
                  f"z=0 -> z=0")

    ms.compute_selection_by_condition_per_vertex(condselect="(z > 1e-6)")
    n_above = ms.current_mesh().selected_vertex_number()
    if n_above:
        print(f"  WARNING: {n_above} verts have z > 1e-6 after remesh.",
              file=sys.stderr)
    elif verbose:
        print(f"  ok: all output verts have z <= 1e-6")
    ms.compute_selection_by_condition_per_vertex(condselect="(z > 1e30)")

    box = ms.current_mesh().bounding_box()
    if verbose:
        print(f"  output bbox: "
              f"x=[{box.min()[0]:.1f}, {box.max()[0]:.1f}], "
              f"y=[{box.min()[1]:.1f}, {box.max()[1]:.1f}], "
              f"z=[{box.min()[2]:.1f}, {box.max()[2]:.1f}]")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    ms.save_current_mesh(str(out_path), binary=False)
    if verbose:
        size_mb = out_path.stat().st_size / 1e6
        print(f"wrote {out_path} ({size_mb:.2f} MB, ascii STL)")


# ----------------------------------------------------------------------
# Method: vdel (legacy vertex-delete)
# ----------------------------------------------------------------------

def clean_freesurface_vdel(in_path: Path, out_path: Path,
                           target_length: float,
                           snap_strip: float | None = None,
                           iterations: int = 3,
                           reproject: bool = False,
                           verbose: bool = True) -> None:
    if snap_strip is None:
        snap_strip = 0.10 * target_length

    if in_path.suffix.lower() == ".ts":
        ms = _ts_to_meshset(in_path)
    else:
        ms = ml.MeshSet()
        ms.load_new_mesh(str(in_path))
    m = ms.current_mesh()
    if verbose:
        print(f"loaded {in_path.name}: {m.vertex_number()} verts, "
              f"{m.face_number()} faces")
        print(f"  snap_strip = {snap_strip} m")

    ms.compute_selection_by_condition_per_vertex(condselect="(z > 0)")
    n_v_cap = ms.current_mesh().selected_vertex_number()
    ms.compute_selection_by_condition_per_face(
        condselect="(z0 > 0) || (z1 > 0) || (z2 > 0)")
    n_f_cap = ms.current_mesh().selected_face_number()
    if verbose:
        print(f"  cap: {n_v_cap} verts above z=0, {n_f_cap} cap faces")
    ms.meshing_remove_selected_vertices_and_faces()
    if verbose:
        m = ms.current_mesh()
        print(f"  after delete:  {m.vertex_number()} verts, "
              f"{m.face_number()} faces")

    ms.meshing_remove_duplicate_vertices()
    ms.meshing_remove_unreferenced_vertices()
    ms.meshing_remove_null_faces()
    try:
        ms.meshing_repair_non_manifold_edges(method="Remove Faces")
    except Exception as exc:
        if verbose:
            print(f"  (repair_non_manifold_edges skipped: {exc})",
                  file=sys.stderr)

    cond = f"(z > {-snap_strip})"
    ms.compute_selection_by_condition_per_vertex(condselect=cond)
    n_strip = ms.current_mesh().selected_vertex_number()
    if verbose:
        print(f"  snap strip {cond}: {n_strip} verts -> z=0")
    ms.compute_coord_by_function(x="x", y="y", z="0", onselected=True)

    if verbose:
        print(f"  isotropic remesh: targetlen={target_length} m, "
              f"iters={iterations}, reproject={reproject}")
    ms.meshing_isotropic_explicit_remeshing(
        iterations=iterations,
        adaptive=False,
        selectedonly=False,
        targetlen=ml.PureValue(target_length),
        featuredeg=30.0,
        checksurfdist=False,
        splitflag=True,
        collapseflag=True,
        swapflag=True,
        smoothflag=True,
        reprojectflag=reproject,
    )
    m = ms.current_mesh()
    if verbose:
        print(f"  after remesh:  {m.vertex_number()} verts, "
              f"{m.face_number()} faces")

    ms.meshing_remove_duplicate_vertices()
    ms.meshing_remove_unreferenced_vertices()
    ms.meshing_remove_null_faces()

    ms.compute_selection_by_condition_per_vertex(condselect="(z > 1e-6)")
    n_above = ms.current_mesh().selected_vertex_number()
    if n_above:
        print(f"  WARNING: {n_above} verts have z > 1e-6 after remesh.",
              file=sys.stderr)
    elif verbose:
        print(f"  ok: all output verts have z <= 1e-6")
    ms.compute_selection_by_condition_per_vertex(condselect="(z > 1e30)")

    box = ms.current_mesh().bounding_box()
    if verbose:
        print(f"  output bbox: "
              f"x=[{box.min()[0]:.1f}, {box.max()[0]:.1f}], "
              f"y=[{box.min()[1]:.1f}, {box.max()[1]:.1f}], "
              f"z=[{box.min()[2]:.1f}, {box.max()[2]:.1f}]")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    ms.save_current_mesh(str(out_path), binary=False)
    if verbose:
        size_mb = out_path.stat().st_size / 1e6
        print(f"wrote {out_path} ({size_mb:.2f} MB, ascii STL)")


# ----------------------------------------------------------------------
# Top-level dispatcher and CLI
# ----------------------------------------------------------------------

def clean_freesurface(input_path: Path,
                      output_path: Path,
                      method: str = "clip",
                      target_length: float | None = None,
                      snap_strip: float | None = None,
                      eps_clip: float | None = None,
                      iterations: int = 3,
                      reproject: bool = False,
                      verbose: bool = True) -> None:
    in_path = Path(input_path)
    out_path = Path(output_path)
    if target_length is None:
        target_length = _infer_target_length(in_path)
        if target_length is None:
            raise ValueError(
                f"could not infer target_length from {in_path.name!r}; "
                f"pass --target-length")

    if method == "clip":
        clean_freesurface_clip(in_path, out_path,
                               target_length=target_length,
                               eps_clip=eps_clip,
                               iterations=iterations,
                               reproject=reproject,
                               verbose=verbose)
    elif method == "vdel":
        clean_freesurface_vdel(in_path, out_path,
                               target_length=target_length,
                               snap_strip=snap_strip,
                               iterations=iterations,
                               reproject=reproject,
                               verbose=verbose)
    else:
        raise ValueError(f"unknown method: {method!r}")


def _batch(here: Path, method: str, iterations: int,
           reproject: bool) -> int:
    """Process every *_unclipped.stl in ../data_preprocess/ and write
    cleaned outputs to ../data_cleanfreesurf/<base>_clean_<method>.stl."""
    project = here.parent
    in_dir = project / "data_preprocess"
    out_dir = project / "data_cleanfreesurf"
    inputs = sorted(in_dir.glob("*_unclipped.stl"))
    if not inputs:
        print(f"no *_unclipped.stl found under {in_dir}", file=sys.stderr)
        return 1
    for in_path in inputs:
        base = in_path.name.replace("_unclipped.stl", "")
        out_path = out_dir / f"{base}_clean_{method}.stl"
        print(f"\n=== {in_path.name} (method={method}) ===")
        clean_freesurface(in_path, out_path, method=method,
                          iterations=iterations, reproject=reproject)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", type=Path,
                    help="input mesh file (.ts or pymeshlab-readable format)")
    ap.add_argument("output", nargs="?", type=Path,
                    help="output ASCII STL path")
    ap.add_argument("--method", choices=("clip", "vdel"), default="clip",
                    help="clipping method: 'clip' = per-triangle linear "
                         "interpolation at z=0 (default, recommended); "
                         "'vdel' = brute-force vertex-delete + snap-strip "
                         "(legacy)")
    ap.add_argument("--target-length", type=float, default=None,
                    help="isotropic remesh target edge length in world "
                         "units (default: parsed from filename)")
    ap.add_argument("--snap-strip", type=float, default=None,
                    help="vdel only: vertices with z > -SNAP_STRIP are "
                         "snapped to z = 0 (default: 10%% of "
                         "--target-length)")
    ap.add_argument("--eps-clip", type=float, default=None,
                    help="clip only: pre-clip snap radius. Any input "
                         "vertex with |z| < EPS_CLIP is snapped to "
                         "z = 0 *before* the per-triangle linear "
                         "interpolation, eliminating near-degenerate "
                         "clip triangles when an original vertex sat "
                         "just below z = 0 (default: 10%% of "
                         "--target-length; pass 0 to disable)")
    ap.add_argument("--iterations", type=int, default=3,
                    help="isotropic remesh iterations (default: 3)")
    ap.add_argument("--keep-reproject", action="store_true",
                    help="keep reproject-to-original-surface step in the "
                         "isotropic remesher (default OFF, so the cleaned "
                         "top trace stays at z = 0)")
    ap.add_argument("--batch", action="store_true",
                    help="process every *_unclipped.stl in "
                         "../data_preprocess/ and write to "
                         "../data_cleanfreesurf/, ignoring INPUT/OUTPUT "
                         "positional args")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent

    if args.batch:
        return _batch(here, args.method, args.iterations,
                      args.keep_reproject)

    if args.input is None or args.output is None:
        ap.error("INPUT and OUTPUT are required unless --batch is set")
    if not args.input.is_file():
        print(f"error: input not found: {args.input}", file=sys.stderr)
        return 1

    try:
        clean_freesurface(args.input, args.output,
                          method=args.method,
                          target_length=args.target_length,
                          snap_strip=args.snap_strip,
                          eps_clip=args.eps_clip,
                          iterations=args.iterations,
                          reproject=args.keep_reproject)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
