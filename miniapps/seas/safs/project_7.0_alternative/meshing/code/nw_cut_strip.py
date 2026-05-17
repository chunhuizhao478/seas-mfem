#!/usr/bin/env python3
"""
nw_cut_strip.py — Hard-cut the alternative ALT6 long-strip fault STLs at
the NW-most location of the preferred MJVS SAF mesh.

The alternative SAF representation
`SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_*_clean_clip.stl` extends
~236 km further NW than the preferred SAF representation
`SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl`. This
driver:

  1. Loads the preferred mesh and finds its NW-most surface vertex
     (the "anchor"), where "NW projection" is
     `s(p) = (-x + y) / sqrt(2)`.
  2. Defines a vertical cutting plane (parallel to the Z axis) passing
     through that anchor with horizontal normal
     `n_h = (-1, +1, 0) / sqrt(2)`.
  3. For every alternative ALT6 STL in the alternative
     `meshing/results/stl_cleaned/`, clips each triangle to the half-space
     `n_h . p <= s_anchor` (linear interpolation along straddling
     edges, exactly mirroring the `clip_triangle_at_z0` idiom in
     `ts_to_stl.py`). Triangles entirely on the keep side are kept
     whole; triangles entirely on the drop (NW) side are removed.
  4. Welds duplicate vertices (pymeshlab if available, else a
     pure-NumPy weld) and writes a sibling `<base>_nwcut.stl`.

Usage:
    conda activate pythonenv
    python nw_cut_strip.py --batch                  # cut all ALT6 STLs
    python nw_cut_strip.py --print-anchor PATH      # inspect the anchor
    python nw_cut_strip.py INPUT OUTPUT             # cut a single file

The three modes (`--print-anchor`, `--batch`, positional `INPUT OUTPUT`)
are mutually exclusive.

Run unit tests with:
    cd project_7.0_alternative/meshing/code && pytest -q
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np


# ----------------------------------------------------------------------
# Module-level constants
# ----------------------------------------------------------------------

EPS: float = 1.0e-9                                   # consistent with ts_to_stl.EPS

NW_DIRECTION_XY: np.ndarray = np.array(
    [-1.0, 1.0], dtype=np.float64) / np.sqrt(2.0)     # horizontal NW unit vector

DEFAULT_REFERENCE: Path = (
    Path(__file__).resolve().parents[3]
    / "project_7.0_preferred" / "data_cleanfreesurf"
    / "SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m_clean_clip.stl"
)
DEFAULT_ALT_DIR: Path = (
    Path(__file__).resolve().parents[1] / "results" / "stl_cleaned"
)
DEFAULT_ALT_GLOB: str = (
    "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_*_clean_clip.stl"
)
DEFAULT_SUFFIX: str = "_nwcut"

# Pre-clip plane snap.  Vertices within `snap_tol` metres of the cutting
# plane are projected exactly onto it before `clip_mesh_at_plane` runs.
# This eliminates output slivers whose origin is an input vertex sitting
# at sub-LC distance from the plane (see PLAN_nw_hard_cut.md, "sliver
# floor" analysis).  None means "auto" -> 0.1 x median input edge length.
DEFAULT_SNAP_TOL: Optional[float] = None
DEFAULT_SNAP_TOL_FRACTION: float = 0.1   # used when snap_tol is None

# Post-clip close-vertex merge.  Two near-plane vertices that both get
# snapped onto the plane can land within sub-LC of each other, producing
# a sliver cut-line edge that the snap alone can't kill (the snap moves
# vertices but does not merge them).  After clip + dup-vertex cleanup,
# pymeshlab's `meshing_merge_close_vertices` welds any pair within
# `merge_tol` metres.  None means "auto" -> snap_tol.  Pass 0 to disable.
DEFAULT_MERGE_TOL: Optional[float] = None


# ----------------------------------------------------------------------
# Phase 1 — anchor extraction
# ----------------------------------------------------------------------

def load_vertices(path: Path) -> np.ndarray:
    """Return an (N, 3) float64 array of triangle-referenced vertex
    coordinates from any mesh file readable by `meshio.read`. ASCII STL,
    binary STL, PLY, OBJ, and `.ts` (via `ts_to_stl.parse_ts`) all work.

    Orphan vertices that no triangle references are stripped (pymeshlab-
    cleaned STLs and raw `.ts` files can both carry orphans that would
    otherwise pollute the NW anchor).

    Raises:
        FileNotFoundError if the path does not exist.
        ValueError       if the file contains no triangle cells.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"mesh file not found: {p}")

    if p.suffix.lower() == ".ts":
        # Local import: ts_to_stl is a sibling module.
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from ts_to_stl import parse_ts  # noqa: E402
        verts_dict, tris = parse_ts(p)
        if not tris:
            raise ValueError(f"no triangles in {p}")
        # Build the set of IDs that any triangle references; only those
        # vertices are returned.
        referenced = set()
        for ia, ib, ic in tris:
            referenced.add(ia)
            referenced.add(ib)
            referenced.add(ic)
        # Some referenced IDs may not exist in verts_dict (malformed .ts);
        # silently drop them, matching the behaviour of
        # clean_freesurface_mesh.clip_to_arrays.
        rows = [verts_dict[i] for i in sorted(referenced) if i in verts_dict]
        if not rows:
            raise ValueError(f"no triangle-referenced vertices in {p}")
        return np.ascontiguousarray(np.array(rows, dtype=np.float64))

    import meshio  # imported lazily so test fixtures can stub if needed
    m = meshio.read(str(p))
    referenced_idx = set()
    for cb in m.cells:
        if cb.type == "triangle":
            for tri in cb.data:
                referenced_idx.add(int(tri[0]))
                referenced_idx.add(int(tri[1]))
                referenced_idx.add(int(tri[2]))
    if not referenced_idx:
        raise ValueError(f"no triangles in {p}")
    idx = np.array(sorted(referenced_idx), dtype=np.int64)
    pts = np.asarray(m.points, dtype=np.float64)[idx]
    return np.ascontiguousarray(pts, dtype=np.float64)


def nw_anchor(verts: np.ndarray) -> tuple[np.ndarray, float]:
    """Return (p_anchor, s_max) where p_anchor is the (3,) float64 vertex
    of the input that maximises s(p) = (-x + y) / sqrt(2), and s_max is
    that maximum projection. Ties within EPS are broken by smallest Z
    (deterministic; matches 'top trace' intuition).
    """
    if verts.ndim != 2 or verts.shape[1] != 3:
        raise ValueError(
            f"verts must have shape (N, 3); got {verts.shape}")
    if verts.shape[0] == 0:
        raise ValueError("verts is empty")
    # Operand order matters: NW_DIRECTION_XY @ verts[:, :2] would raise
    # ValueError on any real mesh (see PLAN_nw_hard_cut.md Phase 1 step 3).
    s = verts[:, :2] @ NW_DIRECTION_XY                # shape (N,)
    s_max = float(s.max())
    mask = s >= s_max - EPS
    candidates = np.where(mask)[0]
    # Among tied rows pick the one with smallest z.
    z_among = verts[candidates, 2]
    winner = candidates[int(np.argmin(z_among))]
    p_anchor = np.ascontiguousarray(verts[winner], dtype=np.float64)
    return p_anchor, s_max


def cutting_plane(p_anchor: np.ndarray) -> tuple[np.ndarray, float]:
    """Return (n, c) such that the cutting plane is {p : n . p = c} with
    n = (-1, +1, 0)/sqrt(2) (horizontal NW normal) and c = n . p_anchor.
    The keep half-space is {p : n . p <= c}.
    """
    if p_anchor.shape != (3,):
        raise ValueError(
            f"p_anchor must have shape (3,); got {p_anchor.shape}")
    n = np.array([NW_DIRECTION_XY[0], NW_DIRECTION_XY[1], 0.0],
                 dtype=np.float64)
    c = float(n @ p_anchor)
    return n, c


# ----------------------------------------------------------------------
# Phase 2 — generic vertical-plane triangle clipper
# ----------------------------------------------------------------------

def signed_dist(p, n: np.ndarray, c: float) -> float:
    """Return n . p - c. > 0 means the point is on the drop (NW) side;
    <= 0 means the point is on the keep (SE) side.
    """
    return float(n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - c)


def lerp_to_plane(p_keep, p_drop,
                  n: np.ndarray, c: float) -> tuple[float, float, float]:
    """Linearly interpolate between p_keep (d <= 0) and p_drop (d > 0) to
    find the point where the segment crosses the plane n . p = c.

    Mirrors `ts_to_stl.lerp_to_z0` (t = -zb / (za - zb)) generalised to a
    generic plane: t = -d_keep / (d_drop - d_keep).
    """
    d_keep = signed_dist(p_keep, n, c)
    d_drop = signed_dist(p_drop, n, c)
    denom = d_drop - d_keep
    if abs(denom) < EPS:
        # Segment is parallel to / lies in the plane to within float64
        # noise. Snap p_keep onto the plane: p_keep - d_keep * n.
        return (p_keep[0] - d_keep * n[0],
                p_keep[1] - d_keep * n[1],
                p_keep[2] - d_keep * n[2])
    t = -d_keep / denom
    return (p_keep[0] + t * (p_drop[0] - p_keep[0]),
            p_keep[1] + t * (p_drop[1] - p_keep[1]),
            p_keep[2] + t * (p_drop[2] - p_keep[2]))


def clip_triangle_at_plane(p0, p1, p2,
                           n: np.ndarray, c: float) -> list:
    """Clip a triangle to the half-space n . p <= c. Returns 0, 1, or 2
    sub-triangles preserving original CCW vertex order (so outward
    normals are unchanged).

    Algorithm mirrors `ts_to_stl.clip_triangle_at_z0` with z replaced by
    signed distance from the plane.
    """
    pts = [p0, p1, p2]
    d = [signed_dist(p, n, c) for p in pts]
    keep = [di <= EPS for di in d]
    n_keep = sum(keep)

    if n_keep == 3:
        return [(p0, p1, p2)]
    if n_keep == 0:
        return []

    if n_keep == 1:
        # 1 keep + 2 drop -> single output triangle.
        bi = keep.index(True)
        ai = (bi + 1) % 3
        ci = (bi + 2) % 3
        b = pts[bi]
        a = pts[ai]
        cc = pts[ci]
        a_clip = lerp_to_plane(b, a, n, c)
        c_clip = lerp_to_plane(b, cc, n, c)
        return [(b, a_clip, c_clip)]

    # n_keep == 2: 2 keep + 1 drop -> two output triangles.
    ai = keep.index(False)                            # the drop-side vertex
    bi = (ai + 1) % 3
    ci = (ai + 2) % 3
    a = pts[ai]
    b = pts[bi]
    cc = pts[ci]
    ab = lerp_to_plane(b, a, n, c)
    ac = lerp_to_plane(cc, a, n, c)
    return [(ab, b, cc), (ab, cc, ac)]


def _median_edge_length(V: np.ndarray, F: np.ndarray) -> float:
    """Median edge length of the triangle mesh (V, F).  Used as the
    natural scale for the pre-clip snap tolerance.  Returns 0.0 for
    degenerate input."""
    if F.shape[0] == 0:
        return 0.0
    e0 = np.linalg.norm(V[F[:, 1]] - V[F[:, 0]], axis=1)
    e1 = np.linalg.norm(V[F[:, 2]] - V[F[:, 1]], axis=1)
    e2 = np.linalg.norm(V[F[:, 0]] - V[F[:, 2]], axis=1)
    return float(np.median(np.concatenate([e0, e1, e2])))


def snap_vertices_to_plane(V: np.ndarray, n: np.ndarray, c: float,
                           snap_tol: float) -> tuple[np.ndarray, int]:
    """Project every vertex with `|signed_dist(v, plane)| <= snap_tol`
    exactly onto the plane n . p = c.

    Applied BEFORE clip_mesh_at_plane to give the linear-interpolation
    clipper a positive minimum-feature tolerance: any input vertex within
    snap_tol of the plane becomes exactly on the plane (signed_dist == 0)
    and the existing `keep = [di <= EPS]` rule classifies it as kept,
    so triangles that previously produced sub-snap_tol output edges via
    lerp now either stay whole, are dropped wholesale, or are clipped
    along an edge that already lies on the plane -- no slivers possible.

    snap_tol = 0 disables the snap (legacy behaviour).  Returns
    (V_snapped, n_snapped); V is never mutated in place.
    """
    if snap_tol <= 0.0 or V.shape[0] == 0:
        return V, 0
    V = np.ascontiguousarray(V, dtype=np.float64)
    d = V @ n - c
    mask = np.abs(d) <= snap_tol
    n_snapped = int(mask.sum())
    if n_snapped == 0:
        return V, 0
    V_out = V.copy()
    V_out[mask] -= d[mask, None] * n[None, :]
    return V_out, n_snapped


def clip_mesh_at_plane(verts: np.ndarray, faces: np.ndarray,
                       n: np.ndarray, c: float
                       ) -> tuple[np.ndarray, np.ndarray, dict]:
    """Apply `clip_triangle_at_plane` to every face.

    verts: (N, 3) float64. faces: (M, 3) int32 indices into verts.
    Returns (V_out, F_out, stats) with V_out a possibly-duplicated vertex
    array (caller welds), F_out shape (M', 3), and stats:
        n_in, n_out,
        n_kept_whole, n_dropped, n_clipped_2to1, n_clipped_1to2

    Case naming follows `clean_freesurface_mesh.clip_to_arrays`:
        n_above = number of vertices on the drop side (d > EPS)
        - n_above == 0  -> kept_whole
        - n_above == 1  -> 2 keep + 1 drop -> emit 2 triangles ('1to2')
        - n_above == 2  -> 1 keep + 2 drop -> emit 1 triangle  ('2to1')
        - n_above == 3  -> dropped
    """
    if verts.dtype != np.float64:
        verts = np.ascontiguousarray(verts, dtype=np.float64)
    out_pts: list = []
    out_faces: list = []
    n_kept = n_dropped = n_clip2to1 = n_clip1to2 = 0
    for ia, ib, ic in faces:
        p0 = (float(verts[ia, 0]), float(verts[ia, 1]), float(verts[ia, 2]))
        p1 = (float(verts[ib, 0]), float(verts[ib, 1]), float(verts[ib, 2]))
        p2 = (float(verts[ic, 0]), float(verts[ic, 1]), float(verts[ic, 2]))
        d = (signed_dist(p0, n, c),
             signed_dist(p1, n, c),
             signed_dist(p2, n, c))
        n_above = sum(1 for di in d if di > EPS)
        result = clip_triangle_at_plane(p0, p1, p2, n, c)
        if n_above == 0:
            n_kept += 1
        elif n_above == 1:
            n_clip1to2 += 1
        elif n_above == 2:
            n_clip2to1 += 1
        else:  # n_above == 3
            n_dropped += 1
        for tri in result:
            base = len(out_pts)
            out_pts.extend(tri)
            out_faces.append([base, base + 1, base + 2])
    V_out = np.asarray(out_pts, dtype=np.float64).reshape(-1, 3) \
            if out_pts else np.zeros((0, 3), dtype=np.float64)
    F_out = np.asarray(out_faces, dtype=np.int32).reshape(-1, 3) \
            if out_faces else np.zeros((0, 3), dtype=np.int32)
    stats = {
        "n_in": int(len(faces)),
        "n_out": int(len(out_faces)),
        "n_kept_whole": n_kept,
        "n_dropped": n_dropped,
        "n_clipped_2to1": n_clip2to1,
        "n_clipped_1to2": n_clip1to2,
    }
    return V_out, F_out, stats


# ----------------------------------------------------------------------
# Phase 3 — driver
# ----------------------------------------------------------------------

def _load_alt_mesh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a mesh as (verts (N, 3) float64, faces (M, 3) int32)."""
    import meshio
    m = meshio.read(str(path))
    pts = np.ascontiguousarray(m.points, dtype=np.float64)
    if pts.shape[1] != 3:
        raise ValueError(
            f"expected (N, 3) points in {path}; got {pts.shape}")
    tri_blocks = [np.asarray(cb.data, dtype=np.int32)
                  for cb in m.cells if cb.type == "triangle"]
    if not tri_blocks:
        raise ValueError(f"no triangles in {path}")
    faces = np.ascontiguousarray(np.concatenate(tri_blocks, axis=0),
                                 dtype=np.int32)
    return pts, faces


def _weld_numpy(V: np.ndarray, F: np.ndarray,
                tol: float = 1.0e-6
                ) -> tuple[np.ndarray, np.ndarray]:
    """Pure-NumPy fallback weld used when pymeshlab is unavailable.

    Coordinates are quantised to `tol` metres before deduplication; faces
    are remapped to the deduplicated vertex set, and any face with a
    repeated vertex (degenerate after the weld) is dropped.
    """
    if V.shape[0] == 0:
        return V, F
    q = np.round(V / tol).astype(np.int64)
    _, inverse = np.unique(q, axis=0, return_inverse=True)
    # Reduce to canonical vertex coords by averaging within each cluster.
    n_unique = int(inverse.max()) + 1
    V_new = np.zeros((n_unique, 3), dtype=np.float64)
    counts = np.zeros(n_unique, dtype=np.int64)
    np.add.at(V_new, inverse, V)
    np.add.at(counts, inverse, 1)
    V_new /= counts[:, None]
    F_new = inverse[F]
    # Drop null faces (any two vertex indices equal).
    mask = (
        (F_new[:, 0] != F_new[:, 1])
        & (F_new[:, 1] != F_new[:, 2])
        & (F_new[:, 0] != F_new[:, 2])
    )
    F_new = F_new[mask].astype(np.int32)
    # Drop unreferenced verts.
    used = np.unique(F_new)
    remap = -np.ones(V_new.shape[0], dtype=np.int64)
    remap[used] = np.arange(used.shape[0])
    V_new = V_new[used]
    F_new = remap[F_new].astype(np.int32)
    return np.ascontiguousarray(V_new), np.ascontiguousarray(F_new)


def _drift_snap(V: np.ndarray, n: np.ndarray, c: float,
                drift_tol: float = 1.0) -> tuple[np.ndarray, int]:
    """Project any vertex with `|signed_dist| < drift_tol` exactly onto
    the plane n . p = c. Mirrors the post-remesh "drift snap" pattern in
    `clean_freesurface_mesh.clean_freesurface_clip` (where vertices near
    z = 0 are pulled to z = 0 to undo sub-metre remesh drift).

    Without this snap, pymeshlab's ASCII STL writer (which uses %.6e =
    7 significant digits) jiggles UTM-magnitude vertices by up to ~0.5 m,
    so vertices that were exactly on the cutting plane in float64
    (because of lerp cancellation) end up slightly off after the round-
    trip and the post-cut `s_max ≤ c + EPS` invariant fails.

    Returns (V_snapped, n_snapped).
    """
    if V.shape[0] == 0:
        return V, 0
    d = V @ n - c
    mask = np.abs(d) < drift_tol
    n_snapped = int(mask.sum())
    if n_snapped:
        V = V.copy()
        V[mask] -= d[mask, None] * n[None, :]
    return V, n_snapped


def _drop_degenerate_triangles(V: np.ndarray, F: np.ndarray,
                               eps_edge: float = 1.0e-4
                               ) -> tuple[np.ndarray, np.ndarray, int]:
    """Drop triangles whose shortest edge is below `eps_edge` metres.
    These are float-precision residuals of the lerp clipper (e.g. when
    the keep vertex was snapped onto the plane, the lerp parameter
    becomes 0 and the output triangle's three vertices are equal to
    machine precision but not bit-identical, so pymeshlab's
    `meshing_remove_duplicate_vertices` does not merge them and the
    null-face filter therefore does not catch the resulting tri).

    eps_edge = 1e-4 m (0.1 mm) is well below any real fault feature and
    only catches genuinely degenerate output.  Returns (V, F_clean, n_dropped).
    """
    if F.shape[0] == 0:
        return V, F, 0
    e0 = np.linalg.norm(V[F[:, 1]] - V[F[:, 0]], axis=1)
    e1 = np.linalg.norm(V[F[:, 2]] - V[F[:, 1]], axis=1)
    e2 = np.linalg.norm(V[F[:, 0]] - V[F[:, 2]], axis=1)
    keep = (e0 >= eps_edge) & (e1 >= eps_edge) & (e2 >= eps_edge)
    n_dropped = int((~keep).sum())
    return V, np.ascontiguousarray(F[keep]), n_dropped


def _cleanup_pymeshlab(V_out: np.ndarray, F_out: np.ndarray,
                       n: np.ndarray, c: float,
                       verbose: bool,
                       merge_tol: float = 0.0
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Run the pymeshlab cleanup chain (weld dups, optionally merge
    close vertices, drop nulls, drop orphans, repair non-manifold
    edges) followed by a degenerate-tri drop and a drift snap onto
    the cutting plane. Falls back to a NumPy-only weld + snap if
    pymeshlab cannot be imported.

    `merge_tol` (metres) > 0 enables `meshing_merge_close_vertices`
    after the dup-vertex pass; pairs within that absolute distance
    are welded.  This is the only step that can collapse cut-line
    sliver edges produced by two near-plane vertices both snapping
    onto the plane.  0 disables it.
    """
    try:
        import pymeshlab as ml
    except Exception as exc:
        if verbose:
            print(f"  pymeshlab unavailable ({exc!r}); using NumPy weld",
                  file=sys.stderr)
        V_clean, F_clean = _weld_numpy(V_out, F_out)
    else:
        ms = ml.MeshSet()
        ms.add_mesh(ml.Mesh(vertex_matrix=V_out, face_matrix=F_out),
                    "nwcut")
        ms.meshing_remove_duplicate_vertices()
        if merge_tol > 0.0:
            ms.meshing_merge_close_vertices(
                threshold=ml.PureValue(float(merge_tol)))
        ms.meshing_remove_null_faces()
        ms.meshing_remove_unreferenced_vertices()
        try:
            ms.meshing_repair_non_manifold_edges(method="Remove Faces")
        except Exception as exc:
            if verbose:
                print(f"  (repair_non_manifold_edges skipped: {exc})",
                      file=sys.stderr)
        m = ms.current_mesh()
        V_clean = np.ascontiguousarray(m.vertex_matrix(), dtype=np.float64)
        F_clean = np.ascontiguousarray(m.face_matrix(), dtype=np.int32)
    # Pymeshlab merges only bit-identical duplicates; the lerp clipper
    # produces "essentially identical" vertices (sub-LSB float drift)
    # when t == 0 numerically, which survive as pseudo-tris with
    # sub-micron edges.  Drop them by an absolute edge threshold.
    V_clean, F_clean, n_degen = _drop_degenerate_triangles(V_clean, F_clean)
    if verbose and n_degen:
        print(f"  dropped {n_degen} degenerate triangles "
              f"(min edge < 1e-4 m, sub-LSB lerp residuals)")
    V_clean, n_snapped = _drift_snap(V_clean, n, c)
    if verbose and n_snapped:
        print(f"  drift snap: {n_snapped} vertices within 1 m of plane "
              f"-> projected onto plane")
    return V_clean, F_clean


def _write_ascii_stl(V: np.ndarray, F: np.ndarray, path: Path,
                     name: str = "nwcut") -> None:
    """Write an ASCII STL with %.6f vertex precision (UTM-magnitude
    coords keep 13+ significant digits, well within float64). Mirrors
    `ts_to_stl.write_stl_ascii`. We do *not* use
    `pymeshlab.save_current_mesh(..., binary=False)` because pymeshlab
    formats vertices with `%.6e` (only 7 significant digits), which on
    UTM coords ~4e6 m loses ~0.5 m of precision; that round-trip drift
    moves clipped vertices off the cutting plane.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for ia, ib, ic in F:
            p0 = V[ia]; p1 = V[ib]; p2 = V[ic]
            ux, uy, uz = p1 - p0
            vx, vy, vz = p2 - p0
            nx = uy * vz - uz * vy
            ny = uz * vx - ux * vz
            nz = ux * vy - uy * vx
            mag = (nx * nx + ny * ny + nz * nz) ** 0.5
            if mag > 0.0:
                nx /= mag; ny /= mag; nz /= mag
            f.write(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n")
            f.write("    outer loop\n")
            for p in (p0, p1, p2):
                f.write(f"      vertex {p[0]:.6f} {p[1]:.6f} "
                        f"{p[2]:.6f}\n")
            f.write("    endloop\n  endfacet\n")
        f.write(f"endsolid {name}\n")


def nw_cut_file(input_path: Path, output_path: Path,
                plane: tuple,
                cleanup: bool = True, verbose: bool = True,
                snap_tol: Optional[float] = DEFAULT_SNAP_TOL,
                merge_tol: Optional[float] = DEFAULT_MERGE_TOL) -> dict:
    """Apply (n, c) = `plane` to `input_path` and write `output_path`.

    `snap_tol` (metres) is the pre-clip plane-snap tolerance: input
    vertices within `snap_tol` of the cutting plane are projected onto
    it before clipping, so the linear-interpolation clipper cannot
    produce output edges shorter than `snap_tol`.  When `snap_tol` is
    None (the default), it auto-derives from the input mesh as
    DEFAULT_SNAP_TOL_FRACTION x median input edge length.  Set
    `snap_tol = 0.0` to disable the snap (legacy behaviour, will
    re-introduce sliver triangles along the cut line).

    `merge_tol` (metres) is the post-clip close-vertex weld tolerance:
    after dup-vertex cleanup, pairs within `merge_tol` are merged.
    This is the only step that can collapse cut-line slivers caused
    by two near-plane vertices both snapping onto the plane.  None ->
    auto = `snap_tol`.  0 disables.

    Returns the stats dict from `clip_mesh_at_plane` plus:
        bbox_in, bbox_out: 2x3 ndarrays (min, max) of (x, y, z)
        s_max_in:          NW projection max of the input mesh
        s_max_out:         NW projection max of the cleaned output mesh
        snap_tol:          tolerance actually used
        n_pre_snap:        vertices snapped to plane before clipping
        merge_tol:         close-vertex merge tolerance actually used

    Raises ValueError("refusing to write empty STL: every triangle was
    clipped away") if every input triangle would be removed.
    """
    n, c = plane
    n = np.asarray(n, dtype=np.float64).reshape(3)
    c = float(c)
    in_path = Path(input_path)
    out_path = Path(output_path)

    verts_alt, faces_alt = _load_alt_mesh(in_path)
    bbox_in = np.array([verts_alt.min(axis=0), verts_alt.max(axis=0)])
    s_in = verts_alt[:, :2] @ NW_DIRECTION_XY
    s_max_in = float(s_in.max())
    if verbose:
        print(f"  loaded {in_path.name}: "
              f"{verts_alt.shape[0]} verts, {faces_alt.shape[0]} faces")
        print(f"  in  bbox: x=[{bbox_in[0, 0]:.1f}, {bbox_in[1, 0]:.1f}], "
              f"y=[{bbox_in[0, 1]:.1f}, {bbox_in[1, 1]:.1f}], "
              f"z=[{bbox_in[0, 2]:.1f}, {bbox_in[1, 2]:.1f}]")
        print(f"  in  s_max = {s_max_in:.3f}")

    # Pre-clip plane snap.  Without this, vertices within sub-LC of the
    # plane create sliver output triangles whose edges propagate through
    # the gmsh discrete-fault embedding and force bad-quality tets in
    # the bulk mesh.
    median_edge = _median_edge_length(verts_alt, faces_alt)
    if snap_tol is None:
        snap_tol_used = DEFAULT_SNAP_TOL_FRACTION * median_edge
    else:
        snap_tol_used = float(snap_tol)
    if merge_tol is None:
        merge_tol_used = snap_tol_used
    else:
        merge_tol_used = float(merge_tol)
    verts_alt, n_pre_snap = snap_vertices_to_plane(
        verts_alt, n, c, snap_tol_used)
    if verbose:
        print(f"  pre-clip snap: tol={snap_tol_used:.3f} m "
              f"(median input edge = {median_edge:.1f} m); "
              f"{n_pre_snap} vertices projected onto plane")
        print(f"  post-clip merge_tol={merge_tol_used:.3f} m")

    V_out, F_out, stats = clip_mesh_at_plane(verts_alt, faces_alt, n, c)

    if F_out.shape[0] == 0:
        raise ValueError(
            "refusing to write empty STL: every triangle was clipped away")

    if cleanup:
        V_clean, F_clean = _cleanup_pymeshlab(
            V_out, F_out, n, c, verbose, merge_tol=merge_tol_used)
    else:
        V_clean, F_clean = V_out, F_out

    _write_ascii_stl(V_clean, F_clean, out_path)

    bbox_out = np.array([V_clean.min(axis=0), V_clean.max(axis=0)])
    s_out = V_clean[:, :2] @ NW_DIRECTION_XY
    s_max_out = float(s_out.max())
    s_min_out = float(s_out.min())
    nw_end = np.ascontiguousarray(V_clean[int(s_out.argmax())],
                                  dtype=np.float64)
    se_end = np.ascontiguousarray(V_clean[int(s_out.argmin())],
                                  dtype=np.float64)

    stats.update({
        "bbox_in": bbox_in,
        "bbox_out": bbox_out,
        "s_max_in": s_max_in,
        "s_max_out": s_max_out,
        "s_min_out": s_min_out,
        "nw_endpoint": nw_end,
        "se_endpoint": se_end,
        "n_verts_in": int(verts_alt.shape[0]),
        "n_verts_out": int(V_clean.shape[0]),
        "snap_tol": float(snap_tol_used),
        "n_pre_snap": int(n_pre_snap),
        "median_input_edge": float(median_edge),
        "merge_tol": float(merge_tol_used),
    })

    if verbose:
        print(f"  clip: in={stats['n_in']}, out={stats['n_out']} "
              f"(kept_whole={stats['n_kept_whole']}, "
              f"dropped={stats['n_dropped']}, "
              f"clipped 2->1={stats['n_clipped_2to1']}, "
              f"clipped 1->2={stats['n_clipped_1to2']})")
        print(f"  out bbox: x=[{bbox_out[0, 0]:.1f}, {bbox_out[1, 0]:.1f}], "
              f"y=[{bbox_out[0, 1]:.1f}, {bbox_out[1, 1]:.1f}], "
              f"z=[{bbox_out[0, 2]:.1f}, {bbox_out[1, 2]:.1f}]")
        print(f"  out s_max = {s_max_out:.3f}  "
              f"(plane c = {c:.3f}; delta = {s_max_out - c:.2e})")
        print(f"  NW endpoint: x={nw_end[0]:.3f}  y={nw_end[1]:.3f}  "
              f"z={nw_end[2]:.3f}  s={s_max_out:.3f}")
        print(f"  SE endpoint: x={se_end[0]:.3f}  y={se_end[1]:.3f}  "
              f"z={se_end[2]:.3f}  s={s_min_out:.3f}")
        print(f"  strip length along NW-SE = {(s_max_out - s_min_out):.1f} m")
        size_mb = out_path.stat().st_size / 1e6
        print(f"  wrote {out_path.name} ({size_mb:.2f} MB, ascii STL)")

    if stats["n_dropped"] == 0 and stats["n_clipped_2to1"] == 0 \
            and stats["n_clipped_1to2"] == 0:
        print(f"WARNING: no triangles dropped from {in_path.name}; the "
              f"reference anchor may be NW of every input vertex (did "
              f"you swap --reference?)", file=sys.stderr)

    return stats


# ----------------------------------------------------------------------
# Top-level pipeline
# ----------------------------------------------------------------------

def _print_anchor(reference: Path) -> int:
    verts_ref = load_vertices(reference)
    p_anchor, s_max = nw_anchor(verts_ref)
    print(f"anchor: x={p_anchor[0]:.6e}  y={p_anchor[1]:.6e}  "
          f"z={p_anchor[2]:.6e}  s={s_max:.6e}")
    return 0


def _resolve_plane(reference: Path,
                   verbose: bool) -> tuple[tuple[np.ndarray, float], float,
                                           np.ndarray]:
    """Load the reference mesh, compute the cutting plane, return
    (plane, s_max, p_anchor)."""
    if not reference.is_file():
        print(f"reference file not found: {reference}", file=sys.stderr)
        sys.exit(2)
    verts_ref = load_vertices(reference)
    p_anchor, s_max = nw_anchor(verts_ref)
    n, c = cutting_plane(p_anchor)
    if verbose:
        print(f"reference: {reference.name}")
        print(f"  anchor: x={p_anchor[0]:.3f}  y={p_anchor[1]:.3f}  "
              f"z={p_anchor[2]:.3f}")
        print(f"  s_max  = {s_max:.6f}")
        print(f"  plane: n={tuple(round(x, 6) for x in n)}  c={c:.6f}")
    return (n, c), s_max, p_anchor


def _run_single(input_path: Path, output_path: Path,
                reference: Path, verbose: bool,
                snap_tol: Optional[float] = DEFAULT_SNAP_TOL,
                merge_tol: Optional[float] = DEFAULT_MERGE_TOL) -> int:
    plane, _, _ = _resolve_plane(reference, verbose)
    if not input_path.is_file():
        print(f"input file not found: {input_path}", file=sys.stderr)
        return 1
    try:
        nw_cut_file(input_path, output_path, plane,
                    verbose=verbose, snap_tol=snap_tol,
                    merge_tol=merge_tol)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


def _run_batch(reference: Path, alt_dir: Path, glob: str, suffix: str,
               out_dir: Optional[Path], verbose: bool,
               snap_tol: Optional[float] = DEFAULT_SNAP_TOL,
               merge_tol: Optional[float] = DEFAULT_MERGE_TOL) -> int:
    if not alt_dir.is_dir():
        print(f"alt directory not found: {alt_dir}", file=sys.stderr)
        return 1
    candidates = sorted(alt_dir.glob(glob))
    inputs = [p for p in candidates if not p.stem.endswith(suffix)]
    skipped = [p for p in candidates if p.stem.endswith(suffix)]
    for s in skipped:
        print(f"  skip already-cut file: {s.name}", file=sys.stderr)
    if not inputs:
        print(f"no files matched {alt_dir / glob}", file=sys.stderr)
        return 1

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)

    plane, _, _ = _resolve_plane(reference, verbose)

    rc = 0
    for in_path in inputs:
        out_name = f"{in_path.stem}{suffix}.stl"
        if out_dir is not None:
            out_path = out_dir / out_name
        else:
            out_path = in_path.with_name(out_name)
        if verbose:
            print(f"\n=== {in_path.name} -> {out_path} ===")
        try:
            nw_cut_file(in_path, out_path, plane,
                        verbose=verbose, snap_tol=snap_tol,
                        merge_tol=merge_tol)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            rc = 1
            continue
    return rc


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def _parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", type=Path,
                    help="single-file mode: path to the alternative STL "
                         "to cut")
    ap.add_argument("output", nargs="?", type=Path,
                    help="single-file mode: path to write the cut STL")
    ap.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE,
                    help="path to the reference mesh whose NW corner "
                         "defines the cut anchor "
                         "(default: %(default)s)")
    ap.add_argument("--alt-dir", type=Path, default=DEFAULT_ALT_DIR,
                    help="directory glob'd in --batch mode "
                         "(default: %(default)s)")
    ap.add_argument("--glob", type=str, default=DEFAULT_ALT_GLOB,
                    help="glob pattern for --batch mode "
                         "(default: %(default)s)")
    ap.add_argument("--suffix", type=str, default=DEFAULT_SUFFIX,
                    help="suffix appended to each input stem before "
                         ".stl (default: %(default)s)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="directory to write cut STLs into in --batch "
                         "mode (default: write next to each input)")
    ap.add_argument("--batch", action="store_true",
                    help="run over every file matching --alt-dir/--glob")
    ap.add_argument("--print-anchor", type=Path, default=None,
                    metavar="PATH",
                    help="print the NW anchor of PATH and exit "
                         "(no clipping performed)")
    ap.add_argument("--snap-tol", type=float, default=None,
                    metavar="METRES",
                    help="pre-clip plane-snap tolerance: input vertices "
                         "within this distance of the cutting plane are "
                         "projected onto it before clipping, preventing "
                         "sliver triangles along the cut line.  Default: "
                         f"{DEFAULT_SNAP_TOL_FRACTION} x median input "
                         "edge length (auto).  Pass 0 to disable.")
    ap.add_argument("--merge-tol", type=float, default=None,
                    metavar="METRES",
                    help="post-clip close-vertex weld tolerance: pairs "
                         "within this distance are merged after the snap "
                         "+ clip.  Catches cut-line slivers where two "
                         "near-plane vertices both snap onto the plane "
                         "and end up close.  Default: equal to snap_tol "
                         "(auto).  Pass 0 to disable.")
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="verbose progress output")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="suppress progress output (overrides -v)")
    return ap.parse_args(argv)


def main(argv: Optional[list] = None) -> int:
    args = _parse_args(argv)
    verbose = args.verbose and not args.quiet

    # Mutual exclusion: --print-anchor / --batch / (INPUT OUTPUT).
    has_print = args.print_anchor is not None
    has_batch = bool(args.batch)
    has_pos = (args.input is not None) or (args.output is not None)
    n_modes = sum((has_print, has_batch, has_pos))
    if n_modes > 1:
        print("--print-anchor / --batch / INPUT OUTPUT are mutually "
              "exclusive", file=sys.stderr)
        return 2
    if n_modes == 0:
        print("no mode given: use --print-anchor PATH, --batch, or "
              "INPUT OUTPUT", file=sys.stderr)
        return 2
    if has_pos and (args.input is None or args.output is None):
        print("INPUT and OUTPUT must both be supplied", file=sys.stderr)
        return 2

    if has_print:
        return _print_anchor(args.print_anchor)
    if has_batch:
        return _run_batch(args.reference, args.alt_dir, args.glob,
                          args.suffix, args.out_dir, verbose,
                          snap_tol=args.snap_tol,
                          merge_tol=args.merge_tol)
    return _run_single(args.input, args.output, args.reference, verbose,
                       snap_tol=args.snap_tol,
                       merge_tol=args.merge_tol)


if __name__ == "__main__":
    sys.exit(main())
