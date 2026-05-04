#!/usr/bin/env python
"""Per-fault surface remeshing via the MMG Platform's `mmgs_O3` binary.

Tier 1 of `REVIEW_interior_subdivision_and_collapse.md`: addresses the
needle-shaped fault triangles that survive `break_fault_wedges.py`'s
polyline-edge splits.  `mmgs_O3` does split + collapse + edge-flip +
vertex relocation in one pass, with `-hgrad` size-gradient control,
`-hausd` Hausdorff geometry-deviation bound, and explicit
`RequiredEdges` / `RequiredVertices` for cross-fault polyline
preservation.

Pipeline position:
    corefine_faults → break_fault_wedges → THIS → generate_safs_mesh

Cross-fault conformity invariant:
    Every vertex coord whose snap_key appears in ≥ 2 faults is marked
    as a `RequiredVertex` in the medit input to mmgs_O3.  Every edge
    whose both endpoints are polyline vertices is marked as a
    `RequiredEdge`.  After `mmgs_O3 -nr` (no auto angle-detection;
    explicit-required only), polyline geometry is bit-identical to
    the input — therefore bit-identical across faults by induction
    on the input invariant.

What mmgs CAN do (interior of each fault):
    - Insert Steiner points inside elongated triangles.
    - Collapse near-coincident interior vertices (needle removal).
    - Flip interior edges to improve aspect ratio.
    - Smooth interior vertex locations (Laplacian-like, geometry-bounded
      by `-hausd`).

What mmgs CANNOT do here (by design):
    - Modify polyline-edge geometry (RequiredEdges).
    - Move polyline vertices (RequiredVertices).
    - Cross-fault triangle changes (it sees only one fault at a time).

Usage:
    python mmgs_remesh_per_fault.py \\
        --in-stl-dir output/.../stl_wedge_broken \\
        --out-stl-dir output/.../stl_mmgs_remeshed \\
        --include-fault A --include-fault B [...] \\
        [--hmin 200] [--hmax 1000] [--hgrad 1.3] [--hausd 20] \\
        [--snap-m 0.1]

The `mmgs_O3` binary must be on PATH (e.g., `conda install -c
conda-forge mmgsuite` in the active environment).
"""
from __future__ import annotations
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

# Make the sibling module importable so we can reuse its STL I/O and
# the polyline-vertex-key helper.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import break_fault_wedges as bfw  # noqa: E402


# ---------------------------------------------------------------------------
# Medit-format I/O.  Medit is mmg's native ASCII format.  We write only
# the subset we need (Vertices / Edges / RequiredEdges / RequiredVertices /
# Triangles) and read back the same.
#
# Reference: https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf
# (Frey, 2001, "MEDIT: An interactive mesh visualization software")
# ---------------------------------------------------------------------------
def _write_medit_surface(path: Path,
                          V: np.ndarray, T: np.ndarray,
                          required_vertex_indices: list[int],
                          required_edges: list[tuple[int, int]]
                          ) -> None:
    """Write a 3-D surface triangulation to the medit `.mesh` format
    with explicit RequiredVertices and RequiredEdges blocks.

    Vertex / edge / triangle indices in the file are 1-based (medit
    convention).  Each vertex line is `x y z ref`, each edge `v1 v2 ref`,
    each triangle `v1 v2 v3 ref`.  Ref defaults to 0.
    """
    if V.ndim != 2 or V.shape[1] != 3:
        raise ValueError(
            f"V must be (N, 3); got shape {V.shape}")
    if T.ndim != 2 or T.shape[1] != 3:
        raise ValueError(
            f"T must be (M, 3); got shape {T.shape}")
    n_v = V.shape[0]
    n_t = T.shape[0]
    n_e = len(required_edges)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("MeshVersionFormatted 2\n\n")
        f.write("Dimension 3\n\n")

        f.write(f"Vertices\n{n_v}\n")
        for v in V:
            f.write(f"{v[0]:+.17e} {v[1]:+.17e} {v[2]:+.17e} 0\n")
        f.write("\n")

        f.write(f"Triangles\n{n_t}\n")
        for t in T:
            # 1-based indexing.
            f.write(f"{int(t[0])+1} {int(t[1])+1} {int(t[2])+1} 0\n")
        f.write("\n")

        if n_e > 0:
            f.write(f"Edges\n{n_e}\n")
            for (a, b) in required_edges:
                f.write(f"{a+1} {b+1} 0\n")
            f.write("\n")
            f.write(f"RequiredEdges\n{n_e}\n")
            for i in range(n_e):
                f.write(f"{i+1}\n")  # 1-based edge index
            f.write("\n")

        if required_vertex_indices:
            f.write(f"RequiredVertices\n{len(required_vertex_indices)}\n")
            for vi in required_vertex_indices:
                f.write(f"{vi+1}\n")
            f.write("\n")

        f.write("End\n")


def _read_medit_surface(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse vertices and triangles from a medit `.mesh` file.

    Skips Edges / RequiredEdges / RequiredVertices blocks — we only
    need V and T to write back to STL.  Indices in the input file are
    1-based; output arrays use 0-based indexing.
    """
    if not path.exists():
        raise FileNotFoundError(f"medit file not found: {path}")

    verts: list[tuple[float, float, float]] = []
    tris: list[tuple[int, int, int]] = []

    tokens: list[str] = []
    with path.open() as f:
        for line in f:
            # Strip comments after `#`.
            s = line.split("#", 1)[0].strip()
            if not s:
                continue
            tokens.extend(s.split())

    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok == "Vertices":
            i += 1
            n = int(tokens[i]); i += 1
            for _ in range(n):
                x = float(tokens[i]); y = float(tokens[i+1])
                z = float(tokens[i+2])
                # ref token consumed
                i += 4
                verts.append((x, y, z))
        elif tok == "Triangles":
            i += 1
            n = int(tokens[i]); i += 1
            for _ in range(n):
                a = int(tokens[i]) - 1
                b = int(tokens[i+1]) - 1
                c = int(tokens[i+2]) - 1
                i += 4  # ref consumed
                tris.append((a, b, c))
        elif tok == "End":
            break
        else:
            # Unknown / skipped block (Edges, RequiredEdges, RequiredVertices,
            # Normals, Tangents, etc.).  Most blocks of form `<NAME> <N>
            # <N rows of K tokens>` — but the row width varies per block.
            # The simplest robust skip: if next token is an int, treat as
            # block size and consume that many lines, each guessed.
            i += 1
            if i >= len(tokens):
                break
            try:
                n = int(tokens[i]); i += 1
            except ValueError:
                continue
            # Heuristic: skip n entries, each containing 1-4 tokens.
            # We can't always tell; safest is to advance by the next
            # appearance of a known section keyword.
            sections = {"Vertices", "Edges", "Triangles", "Quadrilaterals",
                         "Tetrahedra", "Hexahedra", "RequiredVertices",
                         "RequiredEdges", "RequiredTriangles", "Ridges",
                         "Normals", "Tangents", "NormalAtVertices",
                         "TangentAtVertices", "End"}
            # Move i forward token by token until we hit a section keyword
            # or End; that bounds the unknown block.
            while i < len(tokens) and tokens[i] not in sections:
                i += 1
    return (np.asarray(verts, dtype=np.float64),
            np.asarray(tris, dtype=np.int64))


# ---------------------------------------------------------------------------
# Core: per-fault remesh.
# ---------------------------------------------------------------------------
def _polyline_edges_in_fault(
        V: np.ndarray, T: np.ndarray,
        polyline_keys: set[tuple[int, int, int]],
        snap_m: float
        ) -> list[tuple[int, int]]:
    """Return the unique list of triangle edges (v_a, v_b) where BOTH
    endpoints are polyline vertices.

    Edges are returned with vertex indices into `V`, in the convention
    `(min(a, b), max(a, b))`.  This is the set we mark as
    `RequiredEdges` in the medit input.
    """
    edges: set[tuple[int, int]] = set()
    for t in T:
        for (a, b) in ((int(t[0]), int(t[1])),
                       (int(t[1]), int(t[2])),
                       (int(t[2]), int(t[0]))):
            ka = bfw._snap_key(V[a], snap_m)
            kb = bfw._snap_key(V[b], snap_m)
            if ka in polyline_keys and kb in polyline_keys:
                edges.add((a, b) if a < b else (b, a))
    return sorted(edges)


def _polyline_vertices_in_fault(
        V: np.ndarray,
        polyline_keys: set[tuple[int, int, int]],
        snap_m: float
        ) -> list[int]:
    """Return the list of vertex indices in V whose snap_key is in
    `polyline_keys` (i.e., shared with at least one other fault).
    """
    out: list[int] = []
    for vi in range(V.shape[0]):
        if bfw._snap_key(V[vi], snap_m) in polyline_keys:
            out.append(vi)
    return out


def _run_mmgs(in_mesh: Path, out_mesh: Path,
               hmin: float, hmax: float, hgrad: float, hausd: float,
               extra_args: list[str] | None = None,
               binary: str = "mmgs_O3",
               log_path: Path | None = None) -> None:
    """Invoke `mmgs_O3` on a medit input.

    Always passes `-nr` so that ridges are determined by the explicit
    `RequiredEdges` block in the input file rather than auto-detected
    from dihedral angles.  Auto-detection inside a single fault would
    incorrectly mark internal kinks as ridges and prevent legitimate
    smoothing.
    """
    cmd = [binary,
           "-in", str(in_mesh),
           "-out", str(out_mesh),
           "-hmin", repr(float(hmin)),
           "-hmax", repr(float(hmax)),
           "-hgrad", repr(float(hgrad)),
           "-hausd", repr(float(hausd)),
           "-nr"]
    if extra_args:
        cmd.extend(extra_args)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w") as fh:
            fh.write("CMD: " + " ".join(cmd) + "\n\n")
            fh.flush()
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(
            f"mmgs_O3 failed (rc={r.returncode}) on {in_mesh}; "
            f"see log at {log_path}" if log_path else
            f"mmgs_O3 failed (rc={r.returncode}) on {in_mesh}")


def _canonical_snap_polyline_vertices(
        V: np.ndarray,
        polyline_keys: set[tuple[int, int, int]],
        snap_m: float
        ) -> int:
    """Replace every polyline-vertex coord (snap_key in polyline_keys)
    with its canonical position `snap_key * snap_m`.  Modifies V
    in-place and returns the count of vertices snapped.

    Mmgs introduces machine-epsilon perturbation on RequiredVertex
    coords (observed: differences at ~1e-17 m).  This is far below
    HXT's tolerance, but to make cross-fault conformity bit-exact
    (so test assertions hold and the polyline is unambiguously
    welded by gmsh's STL Merge), we snap polyline vertices back to
    canonical after mmgs returns.  Same operation as
    break_fault_wedges.snap_polyline_vertices_to_canonical, but
    applied per-fault on the post-mmgs output.
    """
    n_snapped = 0
    for vi in range(V.shape[0]):
        key = bfw._snap_key(V[vi], snap_m)
        if key in polyline_keys:
            canon = np.array(key, dtype=np.float64) * snap_m
            if not np.array_equal(V[vi], canon):
                V[vi] = canon
                n_snapped += 1
    return n_snapped


def remesh_per_fault(
        in_stl_dir: Path, out_stl_dir: Path, include_fault: list[str],
        hmin: float = 200.0, hmax: float = 1000.0,
        hgrad: float = 1.3, hausd: float = 20.0,
        snap_m: float = 0.1,
        binary: str = "mmgs_O3",
        keep_intermediate: bool = False,
        nomove: bool = False,
        noswap: bool = False,
        noinsert: bool = False
        ) -> dict:
    """For each fault in `include_fault`, run `mmgs_O3` with:

        - polyline VERTICES marked Required (cannot move),
        - polyline EDGES marked Required (cannot split / collapse / flip),
        - everything else (interior triangles, apex vertices) free for
          mmgs to optimize.

    Returns a per-fault report dict suitable for JSON serialization.
    """
    if not in_stl_dir.is_dir():
        raise SystemExit(f"--in-stl-dir not found: {in_stl_dir}")
    if not include_fault:
        raise SystemExit("--include-fault must appear at least once")
    if hmin <= 0 or hmax <= 0:
        raise ValueError(f"hmin/hmax must be > 0; got hmin={hmin} hmax={hmax}")
    if hmin > hmax:
        raise ValueError(f"hmin must be <= hmax; got hmin={hmin} hmax={hmax}")
    if hgrad < 1.0:
        raise ValueError(f"hgrad must be >= 1.0; got {hgrad}")
    if hausd <= 0:
        raise ValueError(f"hausd must be > 0; got {hausd}")
    out_stl_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load each fault's STL into memory.
    raw: dict[str, tuple[np.ndarray, np.ndarray, str]] = {}
    for short in include_fault:
        path = in_stl_dir / f"{short}.stl"
        if not path.exists():
            raise SystemExit(f"missing input STL: {path}")
        V, T, solid = bfw._read_ascii_stl(path)
        # Dedup at snap precision so triangle-edge detection works
        # consistently with break_fault_wedges.
        V, T = bfw._dedup_vertices(V, T, tol=snap_m)
        raw[short] = (V, T, solid)
        print(f"  load {short}: V={V.shape[0]}, T={T.shape[0]}",
              file=sys.stderr)

    # 2. Identify cross-fault polyline vertex snap_keys.
    fault_meshes = {short: bfw.FaultMesh(short, V, T, solid, snap_m)
                    for short, (V, T, solid) in raw.items()}
    polyline_keys = bfw.collect_polyline_vertex_keys(fault_meshes)
    print(f"  cross-fault polyline vertex keys: {len(polyline_keys)}",
          file=sys.stderr)

    # 3. Per-fault: build medit, run mmgs, read back, write STL.
    per_fault: dict[str, dict] = {}
    intermediate_dir = out_stl_dir / "_mmgs_work"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    for short, (V, T, solid) in raw.items():
        req_v = _polyline_vertices_in_fault(V, polyline_keys, snap_m)
        req_e = _polyline_edges_in_fault(V, T, polyline_keys, snap_m)

        in_mesh = intermediate_dir / f"{short}.mesh"
        out_mesh = intermediate_dir / f"{short}_remeshed.mesh"
        log_path = intermediate_dir / f"{short}_mmgs.log"
        _write_medit_surface(in_mesh, V, T, req_v, req_e)

        extra: list[str] = []
        if nomove:
            extra.append("-nomove")
        if noswap:
            extra.append("-noswap")
        if noinsert:
            extra.append("-noinsert")
        _run_mmgs(in_mesh, out_mesh,
                   hmin=hmin, hmax=hmax, hgrad=hgrad, hausd=hausd,
                   extra_args=extra,
                   binary=binary, log_path=log_path)

        new_V, new_T = _read_medit_surface(out_mesh)
        # Restore bit-exact cross-fault conformity: mmgs may have
        # perturbed RequiredVertex coords at machine epsilon.  Snap
        # every polyline vertex back to its canonical position so
        # both faults emit identical coords for shared snap_keys.
        n_snapped = _canonical_snap_polyline_vertices(
            new_V, polyline_keys, snap_m)
        out_path = out_stl_dir / f"{short}.stl"
        bfw._write_ascii_stl(out_path, new_V, new_T,
                              solid or f"SAFS:{short}:mmgs")
        print(f"  remesh {short}: V={V.shape[0]}->{new_V.shape[0]}, "
              f"T={T.shape[0]}->{new_T.shape[0]}, "
              f"required_v={len(req_v)}, required_e={len(req_e)}, "
              f"post-mmgs canonical-snap={n_snapped}",
              file=sys.stderr)

        per_fault[short] = {
            "V_in":  int(V.shape[0]),
            "V_out": int(new_V.shape[0]),
            "T_in":  int(T.shape[0]),
            "T_out": int(new_T.shape[0]),
            "n_required_vertices": len(req_v),
            "n_required_edges":    len(req_e),
            "n_canonical_snap_post_mmgs": int(n_snapped),
        }

    if not keep_intermediate:
        shutil.rmtree(intermediate_dir, ignore_errors=True)

    return {
        "snap_m": snap_m,
        "hmin": hmin, "hmax": hmax, "hgrad": hgrad, "hausd": hausd,
        "n_polyline_vertex_keys": len(polyline_keys),
        "binary": binary,
        "per_fault": per_fault,
    }


# ---------------------------------------------------------------------------
# CLI.
# ---------------------------------------------------------------------------
def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Per-fault surface remeshing via MMG mmgs_O3.")
    p.add_argument("--in-stl-dir", required=True, type=Path)
    p.add_argument("--out-stl-dir", required=True, type=Path)
    p.add_argument("--include-fault", action="append", default=[],
                   metavar="SHORT_NAME")
    p.add_argument("--hmin", type=float, default=200.0,
                   help="mmgs -hmin (minimal mesh edge length, m).  "
                        "Default 200.")
    p.add_argument("--hmax", type=float, default=1000.0,
                   help="mmgs -hmax (maximal mesh edge length, m).  "
                        "Default 1000.")
    p.add_argument("--hgrad", type=float, default=1.3,
                   help="mmgs -hgrad (maximum edge-length ratio between "
                        "adjacent edges).  Default 1.3.")
    p.add_argument("--hausd", type=float, default=20.0,
                   help="mmgs -hausd (Hausdorff geometry-deviation "
                        "bound, m).  Default 20.")
    p.add_argument("--snap-m", type=float, default=0.1,
                   help="Coordinate snap precision for cross-fault "
                        "polyline vertex matching.  Default 0.1 m.")
    p.add_argument("--mmgs-binary", default="mmgs_O3",
                   help="Path or name of the mmgs binary.  Default "
                        "mmgs_O3 (must be on PATH).")
    p.add_argument("--keep-intermediate", action="store_true",
                   help="Retain the intermediate medit files and "
                        "mmgs logs in <out>/_mmgs_work/ for inspection.")
    p.add_argument("--nomove", action="store_true",
                   help="Pass -nomove to mmgs_O3 (disable vertex "
                        "relocation).  Required for multi-fault input "
                        "where a fault's interior vertex could be "
                        "moved into another fault's surface, creating "
                        "self-intersection in the combined mesh.")
    p.add_argument("--noswap", action="store_true",
                   help="Pass -noswap to mmgs_O3 (disable edge swap).")
    p.add_argument("--noinsert", action="store_true",
                   help="Pass -noinsert to mmgs_O3 (disable point "
                        "insertion / deletion).")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    report = remesh_per_fault(
        in_stl_dir=args.in_stl_dir,
        out_stl_dir=args.out_stl_dir,
        include_fault=args.include_fault,
        hmin=args.hmin, hmax=args.hmax,
        hgrad=args.hgrad, hausd=args.hausd,
        snap_m=args.snap_m,
        binary=args.mmgs_binary,
        keep_intermediate=args.keep_intermediate,
        nomove=args.nomove,
        noswap=args.noswap,
        noinsert=args.noinsert,
    )
    report["mmgs_flags"] = {
        "nomove": bool(args.nomove),
        "noswap": bool(args.noswap),
        "noinsert": bool(args.noinsert),
    }

    report_path = args.report_json or (
        args.out_stl_dir / "mmgs_remesh_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as fh:
        json.dump(report, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
