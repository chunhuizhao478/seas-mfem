"""mmg_cli_cleanup.py — MMG3D volume sliver-removal via the STANDALONE mmg3d_O3
binary, using -opnbdy to preserve the embedded (open-boundary / crack) fault.

Why this instead of mmg_sliver_cleanup.py (mmgpy):  mmgpy 0.12 does NOT expose
MMG3D_IPARAM_opnbdy, so its remesh discards the embedded fault surfaces (the
faults are crack surfaces with free tips + non-manifold junctions; they do not
bound a subdomain, so without -opnbdy MMG treats them as interior-to-a-
homogeneous-domain and drops them — verified: all 529k fault tris lost).  The
standalone mmg3d_O3 (conda-forge `mmgsuite`) supports -opnbdy, which preserves
exactly these surfaces.

Flow: read .msh -> write MEDIT (all surface tris = Triangles+RequiredTriangles,
tets ref 1) -> run `mmg3d_O3 -opnbdy -optim -noinsert -nosurf ...` -> read the
MEDIT output -> VERIFY every fault triangle preserved (per-tag count) and no
fault vertex moved -> write Gmsh v2.2 ASCII.  If the fault is damaged, the
output is NOT written (exit 2).

Usage:
    conda activate pythonenv
    python mmg_cli_cleanup.py IN.msh OUT.msh --mmg-bin /path/to/mmg3d_O3 \
        [--fault-tags 101,102,103] [--hmin 50] [--hmax 200000] [--hgrad 3] \
        [--mem-mb 28000] [--no-optim] [--allow-insert] [--keep-medit] [--verbose]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import meshio
import numpy as np
import mmgpy

import check_mesh_quality as cmq
from mmg_sliver_cleanup import _write_medit, _eta_stats, _print_eta


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile", type=Path)
    ap.add_argument("outfile", type=Path)
    ap.add_argument("--mmg-bin", default="mmg3d_O3",
                    help="path to the mmg3d_O3 executable")
    ap.add_argument("--fault-tags", default="101,102,103")
    ap.add_argument("--hmin", type=float, default=50.0)
    ap.add_argument("--hmax", type=float, default=200000.0)
    ap.add_argument("--hgrad", type=float, default=3.0)
    ap.add_argument("--mem-mb", type=int, default=28000)
    ap.add_argument("--no-optim", action="store_true")
    ap.add_argument("--allow-insert", action="store_true",
                    help="drop -noinsert (lets MMG insert AND collapse/delete)")
    ap.add_argument("--flex-tags", default="",
                    help="comma surface tags MMG MAY remesh (everything else is "
                         "RequiredTriangles=immovable); empty=require all + "
                         "-nosurf. e.g. 201 lets only the DEM free surface flex "
                         "to fix fault/free-surface junction slivers, while fault "
                         "+ box boundary stay exact.")
    ap.add_argument("--hausd", type=float, default=None,
                    help="MMG -hausd: max surface deviation [m] for flex surfaces "
                         "(bounds how far the free surface may move)")
    ap.add_argument("--move-tol-m", type=float, default=1e-3)
    ap.add_argument("--keep-medit", action="store_true")
    ap.add_argument("--skip-run", action="store_true",
                    help="reuse an existing <out>.mmg_out.mesh (skip MEDIT "
                         "write + mmg3d_O3 run); only re-do readback/verify/write")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    if not args.infile.is_file():
        print(f"ERROR: input not found: {args.infile}", file=sys.stderr); return 1
    fault_tags = set(int(x) for x in args.fault_tags.split(","))

    print(f"reading {args.infile} ...")
    m = meshio.read(str(args.infile))
    pts = np.asarray(m.points, dtype=np.float64)
    phys = m.cell_data.get("gmsh:physical")
    tetb, trib, trir = [], [], []
    for bi, cb in enumerate(m.cells):
        if cb.type == "tetra":
            tetb.append(np.asarray(cb.data, dtype=np.int64))
        elif cb.type == "triangle":
            t = np.asarray(cb.data, dtype=np.int64)
            r = (np.asarray(phys[bi], dtype=np.int64) if phys is not None
                 else np.zeros(len(t), dtype=np.int64))
            trib.append(t); trir.append(r)
    if not tetb:
        print("ERROR: no tetra block", file=sys.stderr); return 1
    tets = np.concatenate(tetb) if len(tetb) > 1 else tetb[0]
    surf = np.concatenate(trib) if trib else np.zeros((0, 3), np.int64)
    surf_ref = np.concatenate(trir) if trir else np.zeros(0, np.int64)
    print(f"  nodes={len(pts):,}  tets={len(tets):,}  surf_tris={len(surf):,}")

    pre_fault_counts = {int(t): int((surf_ref == t).sum()) for t in fault_tags}
    fmask = np.isin(surf_ref, list(fault_tags))
    pre_fault_xyz = pts[np.unique(surf[fmask])]
    pre_nonfault_xyz = (pts[np.unique(surf[~fmask])] if (~fmask).any()
                        else np.zeros((0, 3)))
    print(f"  fault tris per tag (input): {pre_fault_counts}")
    print("\n=== input bulk tet quality ===")
    pre = _eta_stats(pts, tets); _print_eta("input", pre)

    medit_in = args.outfile.with_suffix(".mmg_in.mesh")
    medit_out = args.outfile.with_suffix(".mmg_out.mesh")
    if args.skip_run:
        if not medit_out.is_file():
            print(f"ERROR: --skip-run but {medit_out} missing", file=sys.stderr)
            return 3
        print(f"\n--skip-run: reusing existing {medit_out}")
    else:
        flex = set(int(x) for x in args.flex_tags.split(",") if x.strip())
        all_surf_refs = set(int(x) for x in np.unique(surf_ref).tolist())
        req_refs = None if not flex else sorted(all_surf_refs - flex)
        print(f"\nwriting MEDIT input -> {medit_in} "
              f"(flex={sorted(flex) or 'none'}; required={req_refs or 'all'}) ...")
        _write_medit(medit_in, pts, tets, surf, surf_ref, require_refs=req_refs)
        cmd = [args.mmg_bin, "-in", str(medit_in), "-out", str(medit_out),
               "-opnbdy", "-hgrad", str(args.hgrad),
               "-hmin", str(args.hmin), "-hmax", str(args.hmax),
               "-m", str(args.mem_mb), "-v", "5" if args.verbose else "1"]
        if not flex:
            cmd.append("-nosurf")
        if args.hausd is not None:
            cmd += ["-hausd", str(args.hausd)]
        if not args.no_optim:
            cmd.append("-optim")
        if not args.allow_insert:
            cmd.append("-noinsert")
        print("\nrunning:", " ".join(cmd))
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print(f"ERROR: mmg3d_O3 exited {r.returncode}", file=sys.stderr)
            return 3
        if not medit_out.is_file():
            alt = medit_in.with_suffix(".o.mesh")
            if alt.is_file():
                medit_out = alt
            else:
                print("ERROR: mmg output .mesh not found", file=sys.stderr)
                return 3

    # Read back via mmgpy (meshio's MEDIT reader rejects mmg's RequiredTriangles).
    print(f"\nreading MMG output {medit_out} (via mmgpy) ...")
    om = mmgpy.read(str(medit_out))
    new_pts = np.asarray(om.get_vertices(), dtype=np.float64)
    new_tets = np.asarray(om.get_tetrahedra_with_refs()[0], dtype=np.int64)
    tt, rt = om.get_triangles_with_refs()
    new_tris = np.asarray(tt, dtype=np.int64)
    new_triref = np.asarray(rt, dtype=np.int64)
    print(f"after MMG: nodes={len(new_pts):,}  tets={len(new_tets):,}  "
          f"surf_tris={len(new_tris):,}")

    post_fault_counts = {int(t): int((new_triref == t).sum()) for t in fault_tags}
    print(f"  fault tris per tag (output): {post_fault_counts}")
    fault_ok = (post_fault_counts == pre_fault_counts)
    move_max = 0.0
    if fmask.any():
        from scipy.spatial import cKDTree
        pfm = np.isin(new_triref, list(fault_tags))
        post_xyz = new_pts[np.unique(new_tris[pfm])] if pfm.any() else np.zeros((0, 3))
        if len(post_xyz):
            d, _ = cKDTree(pre_fault_xyz).query(post_xyz, k=1)
            move_max = float(d.max())
        print(f"  fault-vertex max displacement: {move_max:.3e} m "
              f"(tol {args.move_tol_m:g})")
        if move_max > args.move_tol_m:
            fault_ok = False

    # informational: free-surface / boundary deviation (allowed under --surf-modify)
    if len(pre_nonfault_xyz):
        from scipy.spatial import cKDTree
        nf_post = np.isin(new_triref, list(fault_tags), invert=True) & (new_triref > 0)
        post_nf_xyz = (new_pts[np.unique(new_tris[nf_post])] if nf_post.any()
                       else np.zeros((0, 3)))
        if len(post_nf_xyz):
            dd, _ = cKDTree(pre_nonfault_xyz).query(post_nf_xyz, k=1)
            print(f"  non-fault surface (DEM/boundary) deviation from input: "
                  f"max {dd.max():.2f} m  p99 {np.percentile(dd, 99):.2f} m  "
                  f"med {np.median(dd):.2f} m  (moved verts: "
                  f"{int((dd > 1.0).sum()):,}/{len(dd):,})")

    print("\n=== output bulk tet quality ===")
    post = _eta_stats(new_pts, new_tets); _print_eta("output", post)
    print(f"\n  Δ sliver(eta<0.1): {pre['n_sliver']:,} -> {post['n_sliver']:,}"
          f"   Δ eta_min: {pre['eta_min']:.4f} -> {post['eta_min']:.4f}"
          f"   Δ edge_min: {pre['edge_min']:.2f} -> {post['edge_min']:.2f} m")

    if not args.keep_medit:
        for p in (medit_in, medit_out):
            try: p.unlink()
            except OSError: pass

    if not fault_ok:
        print("\nFAULT DAMAGED by MMG (per-tag count changed or vertices moved "
              "beyond tol). NOT writing output (exit 2).", file=sys.stderr)
        return 2

    keep = new_triref > 0
    n0 = int((~keep).sum())
    if n0:
        print(f"  note: dropping {n0} MMG-generated ref-0 triangle(s)")
    cells, physical = [("tetra", new_tets)], [np.ones(len(new_tets), int)]
    for tg in sorted(set(new_triref[keep].tolist())):
        sel = new_triref == tg
        cells.append(("triangle", new_tris[sel]))
        physical.append(np.full(int(sel.sum()), int(tg), dtype=int))
    out = meshio.Mesh(points=new_pts, cells=cells,
                      cell_data={"gmsh:physical": physical,
                                 "gmsh:geometrical": [p.copy() for p in physical]})
    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(str(args.outfile), out, file_format="gmsh22", binary=False)
    print(f"\nwrote {args.outfile} (v2.2 ASCII): {len(new_pts):,} nodes, "
          f"{len(new_tets):,} tets, {int(keep.sum()):,} surface tris")
    print("fault preserved: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
