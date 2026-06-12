#!/usr/bin/env python3
"""corefine_faults.py — Python orchestrator for the CGAL corefine pipeline.

Reads `data_cleanfreesurf/*_<R>m_clean_clip.stl`, converts each to OFF,
invokes the C++ `corefine_set` binary, and writes
`data_corefined/*_<R>m_corefined.stl` plus the `manifest.json`.

The C++ binary is expected at
    code_preprocess/corefine_cgal/build/corefine_set
(build with `conda activate cgal-61 && cmake .. && make`).

Usage:
    conda activate pythonenv
    python corefine_faults.py [--res 2000] [--mesh-edge-size 1500.0] ...
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import meshio


def _stl_to_off(stl_path: Path, off_path: Path) -> tuple[int, int]:
    """Read STL via meshio, write OFF at full double precision.
    Returns (n_verts, n_tris).
    """
    m = meshio.read(str(stl_path))
    tris = m.cells_dict.get("triangle")
    if tris is None or len(tris) == 0:
        raise RuntimeError(f"no triangle cells in {stl_path}")
    # meshio's OFF writer precision varies by version; we write our own
    # at %.15g to match the C++ writer.
    with open(off_path, "w") as f:
        f.write("OFF\n")
        f.write(f"{len(m.points)} {len(tris)} 0\n")
        for p in m.points:
            f.write(f"{p[0]:.15g} {p[1]:.15g} {p[2]:.15g}\n")
        for t in tris:
            f.write(f"3 {int(t[0])} {int(t[1])} {int(t[2])}\n")
    return len(m.points), len(tris)


def _off_to_stl(off_path: Path, stl_path: Path) -> tuple[int, int]:
    """Read OFF, write ASCII STL with %.15g precision.

    meshio's ASCII STL writer uses Python's default str(float) format
    which yields ~7 significant digits — at UTM Y ~ 4×10^6 m that's
    ~0.1 m absolute precision, far too coarse for polyline-vertex
    conformality (gmsh dedup tolerance is 1 mm).  We write STL ourselves
    so polyline vertices that were bit-identical between meshes A and
    B in the OFF stay bit-identical in the STL.
    """
    import numpy as np
    m = meshio.read(str(off_path))
    tris = m.cells_dict.get("triangle")
    if tris is None or len(tris) == 0:
        raise RuntimeError(f"no triangle cells in {off_path}")
    pts = m.points
    with open(stl_path, "w") as f:
        f.write("solid corefined\n")
        for tri in tris:
            a, b, c = pts[tri[0]], pts[tri[1]], pts[tri[2]]
            n = np.cross(b - a, c - a)
            mag = np.linalg.norm(n)
            if mag > 0:
                n = n / mag
            f.write(f"facet normal {n[0]:.15g} {n[1]:.15g} {n[2]:.15g}\n")
            f.write("  outer loop\n")
            for p in (a, b, c):
                f.write(f"    vertex {p[0]:.15g} {p[1]:.15g} {p[2]:.15g}\n")
            f.write("  endloop\n")
            f.write("endfacet\n")
        f.write("endsolid corefined\n")
    return len(pts), len(tris)


def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", type=Path,
                    default=project / "data_cleanfreesurf",
                    help="directory of *_clean_clip.stl")
    ap.add_argument("--out-dir", type=Path,
                    default=project / "data_corefined",
                    help="output directory")
    ap.add_argument("--res", type=int, default=2000,
                    help="input fixture resolution suffix (only --res files scanned)")
    ap.add_argument("--mesh-edge-size", type=float, default=1500.0,
                    help="passed through to corefine_set")
    ap.add_argument("--min-edge", type=float, default=100.0)
    ap.add_argument("--polyline-spacing", type=float, default=None,
                    help="default: max(0.5*mesh-edge-size, min-edge)")
    ap.add_argument("--features-angle-bound", type=float, default=60.0)
    ap.add_argument("--cgal-bin", type=Path,
                    default=here / "corefine_cgal" / "build" / "corefine_set")
    ap.add_argument("--workdir", type=Path, default=None,
                    help="default: <here>/work/<R>m")
    ap.add_argument("--keep-intermediate", action="store_true",
                    help="keep workdir after run (else removed)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.workdir is None:
        args.workdir = here / "work" / f"{args.res}m"
    if args.polyline_spacing is None:
        args.polyline_spacing = max(0.5 * args.mesh_edge_size, args.min_edge)
    elif args.polyline_spacing < args.min_edge:
        args.polyline_spacing = args.min_edge

    # Locate inputs.
    inputs = sorted(args.in_dir.glob(f"*_{args.res}m_clean_clip.stl"))
    if not inputs:
        print(f"error: no inputs in {args.in_dir} for res={args.res}m",
              file=sys.stderr)
        return 1
    if args.verbose:
        print(f"found {len(inputs)} input STL(s) at {args.res}m")

    # Verify cgal-bin.
    if not args.cgal_bin.is_file():
        print(f"error: corefine_set binary not found at {args.cgal_bin}\n"
              f"  to build:  conda activate cgal-61 &&\n"
              f"             (cd {here / 'corefine_cgal'} && mkdir -p build && cd build &&\n"
              f"              cmake .. && make)\n",
              file=sys.stderr)
        return 1

    # Setup workdir/in and workdir/out (R-002: separate dirs).
    in_subdir  = args.workdir / "in"
    out_subdir = args.workdir / "out"
    if in_subdir.exists():
        shutil.rmtree(in_subdir)
    if out_subdir.exists():
        shutil.rmtree(out_subdir)
    in_subdir.mkdir(parents=True)
    out_subdir.mkdir(parents=True)

    # Convert STL -> OFF in workdir/in.
    for stl in inputs:
        base = stl.name.replace("_clean_clip.stl", "")
        off = in_subdir / f"{base}.off"
        nv, nt = _stl_to_off(stl, off)
        if args.verbose:
            print(f"  STL->OFF  {stl.name} -> {off.name}  V={nv} F={nt}")

    # Run corefine_set.
    cmd = [str(args.cgal_bin), str(in_subdir), str(out_subdir),
           "--ext", ".off",
           "--mesh-edge-size", str(args.mesh_edge_size),
           "--min-edge", str(args.min_edge),
           "--polyline-spacing", str(args.polyline_spacing),
           "--features-angle-bound", str(args.features_angle_bound),
           "--manifest", str(out_subdir / "manifest.json")]
    if args.verbose:
        cmd.append("--verbose")
        print("running:", " ".join(cmd))
    res = subprocess.run(cmd)
    if res.returncode not in (0, 7):
        # 0 = ok, 7 = some pair non-conformal (still wrote outputs).
        print(f"error: corefine_set exited {res.returncode}", file=sys.stderr)
        return res.returncode

    # Convert _corefined.off -> _corefined.stl in out_dir.
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for off in sorted(out_subdir.glob("*_corefined.off")):
        stl = args.out_dir / f"{off.stem}.stl"
        nv, nt = _off_to_stl(off, stl)
        if args.verbose:
            print(f"  OFF->STL  {off.name} -> {stl.name}  V={nv} F={nt}")

    # Copy manifest.
    manifest_src = out_subdir / "manifest.json"
    manifest_dst = args.out_dir / "manifest.json"
    shutil.copy(manifest_src, manifest_dst)

    # Verify pairwise conformality.
    with open(manifest_dst) as f:
        manifest = json.load(f)
    n_conformal = sum(1 for p in manifest["pairs"] if p["conformal"])
    n_pairs    = len(manifest["pairs"])
    if args.verbose or n_conformal != n_pairs:
        print(f"  manifest: {n_conformal} of {n_pairs} pairs conformal")
        for p in manifest["pairs"]:
            tag = "OK " if p["conformal"] else "FAIL"
            print(f"    [{tag}] ({p['i']},{p['j']}) shared={p['shared_a']}/{p['shared_b']}"
                  f" max_diff={p['max_diff_m']:.3e}")

    # Cleanup workdir.
    if not args.keep_intermediate:
        shutil.rmtree(args.workdir)

    # Summary.
    print(f"\nwrote {len(inputs)} corefined STLs to {args.out_dir}")
    print(f"manifest: {manifest_dst}")
    print(f"per-fault stats:")
    for x in manifest["meshes"]:
        print(f"  {x['basename'][:55]:<55} V={x['n_verts']:>5} F={x['n_faces']:>5}"
              f" edge_med={x['edge_stats']['median']:.0f}m"
              f" q_min={x['tri_q_stats']['min']:.3f}"
              f" q_med={x['tri_q_stats']['median']:.3f}")

    return 0 if n_conformal == n_pairs else 7


if __name__ == "__main__":
    sys.exit(main())
