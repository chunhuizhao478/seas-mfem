#!/usr/bin/env python3
"""mmg_refine_sizemap.py -- metric-driven LOCAL REFINEMENT via standalone mmg3d_O3.

Adaptation of project_7.0_preferred/meshing/code/mmg_cli_cleanup.py for the
0.5 Hz free-surface refinement (meshing_0d5Hz):
  * passes a vertex metric (-sol metric.sol) built by build_metric_sol.py;
  * insertion ENABLED (this is a refinement, not just a cleanup);
  * fault (tag 101) triangles are RequiredTriangles -- byte-frozen, verified
    afterward (per-tag count unchanged + max vertex displacement < tol);
  * top 102 and absorbing 104 are flex (may be split; -hausd bounds deviation;
    both are planar so splits are geometry-safe);
  * -opnbdy preserves the embedded crack fault (see reference docstring).
"""
import argparse
import subprocess
import sys
from pathlib import Path

import fastmsh as meshio   # memory-lean drop-in: meshio.read OOMs at 48.9M tets
import numpy as np

REF_CODE = Path("/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/"
                "seas/safs/project_7.0_preferred/meshing/code")
sys.path.insert(0, str(REF_CODE))
from mmg_sliver_cleanup import _eta_stats, _print_eta  # noqa: E402

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def _write_medit_reqtets(path, pts, tets, surf, surf_ref, require_refs,
                         required_tets):
    """MEDIT .mesh like mmg_sliver_cleanup._write_medit, plus a
    RequiredTetrahedra block freezing every tet OUTSIDE the refinement region
    (mmg may only modify the free band -> guarantees locality / minimum element
    increase; diagnosed 2026-07-07: without it mmg re-refines the whole graded
    far field, +56% tets)."""
    with open(path, "w") as f:
        f.write("MeshVersionFormatted 2\nDimension 3\n\n")
        f.write(f"Vertices\n{len(pts)}\n")
        np.savetxt(f, np.column_stack([pts, np.zeros(len(pts))]),
                   fmt="%.15g %.15g %.15g %d")
        f.write(f"\nTriangles\n{len(surf)}\n")
        np.savetxt(f, np.column_stack([surf + 1, surf_ref]), fmt="%d")
        req = np.nonzero(np.isin(surf_ref, list(require_refs)))[0] + 1
        f.write(f"\nRequiredTriangles\n{len(req)}\n")
        np.savetxt(f, req.reshape(-1, 1), fmt="%d")
        f.write(f"\nTetrahedra\n{len(tets)}\n")
        np.savetxt(f, np.column_stack([tets + 1, np.ones(len(tets), int)]),
                   fmt="%d")
        if required_tets is not None and len(required_tets):
            f.write(f"\nRequiredTetrahedra\n{len(required_tets)}\n")
            np.savetxt(f, (required_tets + 1).reshape(-1, 1), fmt="%d")
        f.write("\nEnd\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("infile", type=Path)
    ap.add_argument("outfile", type=Path)
    ap.add_argument("--binary", action="store_true",
                    help="gmsh22 BINARY intermediate")
    ap.add_argument("sol", type=Path, help="vertex metric .sol (same point order)")
    ap.add_argument("--mmg-bin",
                    default="/Users/chunhuizhao/miniforge/envs/mmg/bin/mmg3d_O3")
    ap.add_argument("--fault-tags", default="101")
    ap.add_argument("--flex-tags", default="102,104")
    ap.add_argument("--hgrad", type=float, default=1.3)
    ap.add_argument("--hausd", type=float, default=30.0)
    ap.add_argument("--hmin", type=float, default=250.0)
    ap.add_argument("--hmax", type=float, default=25000.0)
    ap.add_argument("--mem-mb", type=int, default=20000)
    ap.add_argument("--move-tol-m", type=float, default=1e-3)
    ap.add_argument("--dilate", type=int, default=2,
                    help="free-region dilation layers beyond the demand tets")
    ap.add_argument("--free-zmax", type=float, default=None,
                    help="only tets with barycenter z < this may be freed "
                         "(volume mode: giant deep cells reach the surface; "
                         "undilated freeing wrecked the fine band, 2026-07-08)")
    ap.add_argument("--free-coarse", action="store_true",
                    help="also free tets whose sol DEMANDS COARSENING "
                         "(sol > 1.3*h_max); for the careful bulk-coarsen pass")
    ap.add_argument("--protect-tets", type=Path, default=None,
                    help="npy of tet indices to freeze IN ADDITION to the "
                         "locality mask (marginal-band protection)")
    ap.add_argument("--no-freeze", action="store_true",
                    help="skip RequiredTetrahedra (rely on the h_max metric "
                         "baseline for far-field stability; probe for the "
                         "insertion-filter blockage, 2026-07-07)")
    ap.add_argument("--seed-shift", type=float, default=0.0,
                    help="tiny uniform metric scale (1+s) to vary the draw")
    ap.add_argument("--keep-medit", action="store_true")
    ap.add_argument("--medit-only", action="store_true",
                    help="write the MEDIT input and EXIT.  Lets mmg3d then run "
                         "bare with the whole machine: the driver needs ~19 GB to "
                         "build the locality freeze and `del` does not reliably "
                         "return it, so driver+mmg concurrently reach ~35 GB of 36.")
    ap.add_argument("--reuse-mmg-out", action="store_true",
                    help="skip the mmg3d subprocess and read an existing "
                         ".mmg_out.mesh.  Lets mmg3d be run BARE from the shell "
                         "with the whole machine: `del` does not reliably return "
                         "the driver's ~12 GB to the OS, and driver+mmg together "
                         "reach ~32 GB of 36 with swap 93% full.")
    ap.add_argument("--reuse-medit", action="store_true",
                    help="reuse an existing .mmg_in.mesh instead of rewriting it "
                         "(np.savetxt over 122M tets costs ~20 min); only valid "
                         "when the mesh, sol and freeze flags are unchanged")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    for p in (args.infile, args.sol):
        if not p.is_file():
            print(f"ERROR: missing {p}", file=sys.stderr)
            return 1
    fault_tags = set(int(x) for x in args.fault_tags.split(","))
    flex = set(int(x) for x in args.flex_tags.split(","))

    print(f"reading {args.infile} ...")
    m = meshio.read(str(args.infile))
    pts = np.asarray(m.points, dtype=np.float64)
    phys = m.cell_data.get("gmsh:physical")
    tetb, trib, trir = [], [], []
    for bi, cb in enumerate(m.cells):
        if cb.type == "tetra":
            tetb.append(np.asarray(cb.data, dtype=np.int64))
        elif cb.type == "triangle":
            trib.append(np.asarray(cb.data, dtype=np.int64))
            trir.append(np.asarray(phys[bi], dtype=np.int64))
    tets = np.concatenate(tetb)
    surf = np.concatenate(trib)
    surf_ref = np.concatenate(trir)
    print(f"  nodes={len(pts):,} tets={len(tets):,} surf_tris={len(surf):,}")

    pre_fault_counts = {int(t): int((surf_ref == t).sum()) for t in fault_tags}
    fmask = np.isin(surf_ref, list(fault_tags))
    pre_fault_xyz = pts[np.unique(surf[fmask])]
    print(f"  fault tris (input): {pre_fault_counts}")
    pre = _eta_stats(pts, tets)
    _print_eta("input", pre)

    # optional draw variation: scale the metric slightly (mmg is nondeterministic
    # anyway, but a tiny scale guarantees distinct draws)
    sol_use = args.sol
    if args.seed_shift:
        txt = args.sol.read_text().splitlines()
        k = txt.index("1 1") + 1
        vals = np.array([float(v) for v in txt[k:] if v.strip() and v.strip() != "End"])
        vals *= (1.0 + args.seed_shift)
        sol_use = args.outfile.with_suffix(".draw.sol")
        with open(sol_use, "w") as f:
            f.write("MeshVersionFormatted 2\nDimension 3\n\nSolAtVertices\n"
                    f"{len(vals)}\n1 1\n")
            np.savetxt(f, vals.reshape(-1, 1), fmt="%.6g")
            f.write("\nEnd\n")

    # --- locality mask: free only tets in the refinement halos (+1 layer) ---
    txt = args.sol.read_text().split()
    k = txt.index("SolAtVertices")
    nsol = int(txt[k + 1])
    if nsol != len(pts):
        print(f"ERROR: sol has {nsol} values for {len(pts)} points", file=sys.stderr)
        return 1
    sol_vals = np.asarray(txt[k + 4:k + 4 + nsol], dtype=float)
    h_max = np.zeros(len(pts))
    for i, j in PAIRS:
        a, b = tets[:, i], tets[:, j]
        L = np.linalg.norm(pts[a] - pts[b], axis=1)
        np.maximum.at(h_max, a, L); np.maximum.at(h_max, b, L)
    free_v = sol_vals < 0.98 * h_max
    if args.free_coarse:
        free_v |= sol_vals > 1.3 * h_max
    free_t = free_v[tets].any(1)
    for _ in range(args.dilate):                     # decouple surface cells from
        free_v2 = np.zeros(len(pts), bool)           # the frozen frontier (draw-4
        free_v2[np.unique(tets[free_t])] = True      # fan-cell lesson, 2026-07-07)
        free_t = free_v2[tets].any(1)
    if args.free_zmax is not None:
        bz = pts[tets].mean(1)[:, 2]
        free_t &= bz < args.free_zmax
    if args.protect_tets is not None:
        prot = np.load(args.protect_tets)
        free_t[prot] = False
        print(f"  marginal-band protection: {len(prot):,} tets forced frozen")
    req_tets = np.flatnonzero(~free_t)
    if args.no_freeze:
        req_tets = np.zeros(0, dtype=np.int64)
        print("locality: --no-freeze (no RequiredTetrahedra)")
    else:
        print(f"locality: {int(free_t.sum()):,} free tets / {len(req_tets):,} frozen "
              f"(RequiredTetrahedra); free verts {int(free_v.sum()):,}")

    medit_in = args.outfile.with_suffix(".mmg_in.mesh")
    medit_out = args.outfile.with_suffix(".mmg_out.mesh")
    all_refs = set(int(x) for x in np.unique(surf_ref).tolist())
    req_refs = sorted(all_refs - flex)
    if args.reuse_medit and medit_in.is_file() and medit_in.stat().st_size > 0:
        print(f"reusing existing MEDIT {medit_in} "
              f"({medit_in.stat().st_size/2**30:.2f} GB)")
    else:
        print(f"writing MEDIT (required={req_refs}, flex={sorted(flex)}) ...")
        _write_medit_reqtets(medit_in, pts, tets, surf, surf_ref,
                             require_refs=req_refs, required_tets=req_tets)

    cmd = [args.mmg_bin, "-in", str(medit_in), "-sol", str(sol_use),
           "-out", str(medit_out), "-opnbdy",
           "-hgrad", str(args.hgrad), "-hausd", str(args.hausd),
           "-hmin", str(args.hmin), "-hmax", str(args.hmax),
           "-m", str(args.mem_mb), "-v", "5" if args.verbose else "1"]
    # free the driver's copy before handing the machine to mmg3d.  pts/tets/
    # req_tets are never referenced again after this point (mmg's output is read
    # back into new_pts/new_tets), but at 122M tets they pin ~6.6 GB for the whole
    # mmg run -- which, against a 26 GB mmg cap on a 36 GB box with swap already
    # 93% full, is the difference between finishing and being OOM-killed.
    n_tets_before = len(tets)
    del tets, pts, req_tets
    import gc as _gc
    _gc.collect()
    if args.medit_only:
        print(f"--medit-only: wrote {medit_in} "
              f"({medit_in.stat().st_size/2**30:.2f} GB); now run mmg3d bare:")
        print("  " + " ".join(cmd))
        return 0
    if args.reuse_mmg_out and medit_out.is_file() and medit_out.stat().st_size > 0:
        print(f"reusing existing mmg output {medit_out} "
              f"({medit_out.stat().st_size/2**30:.2f} GB) -- subprocess skipped")
        r = subprocess.CompletedProcess(cmd, 0)
    else:
        print("running:", " ".join(cmd))
        r = subprocess.run(cmd)
    if r.returncode != 0:
        print(f"ERROR: mmg3d_O3 exited {r.returncode}", file=sys.stderr)
        return 3
    if not medit_out.is_file():
        alt = medit_in.with_suffix(".o.mesh")
        if alt.is_file():
            medit_out = alt
        else:
            print("ERROR: mmg output not found", file=sys.stderr)
            return 3

    print(f"reading MMG output {medit_out} (mmgpy) ...")
    import mmgpy
    om = mmgpy.read(str(medit_out))
    new_pts = np.asarray(om.get_vertices(), dtype=np.float64)
    new_tets = np.asarray(om.get_tetrahedra_with_refs()[0], dtype=np.int64)
    tt, rt = om.get_triangles_with_refs()
    new_tris = np.asarray(tt, dtype=np.int64)
    new_triref = np.asarray(rt, dtype=np.int64)
    print(f"after MMG: nodes={len(new_pts):,} tets={len(new_tets):,} "
          f"surf_tris={len(new_tris):,}  (delta tets {len(new_tets)-n_tets_before:+,})")

    post_fault_counts = {int(t): int((new_triref == t).sum()) for t in fault_tags}
    print(f"  fault tris (output): {post_fault_counts}")
    fault_ok = post_fault_counts == pre_fault_counts
    from scipy.spatial import cKDTree
    pfm = np.isin(new_triref, list(fault_tags))
    post_xyz = new_pts[np.unique(new_tris[pfm])] if pfm.any() else np.zeros((0, 3))
    move_max = 0.0
    if len(post_xyz):
        dq, _ = cKDTree(pre_fault_xyz).query(post_xyz, k=1)
        move_max = float(dq.max())
    print(f"  fault-vertex max displacement: {move_max:.3e} m (tol {args.move_tol_m:g})")
    if move_max > args.move_tol_m:
        fault_ok = False

    post = _eta_stats(new_pts, new_tets)
    _print_eta("output", post)
    print(f"  d_sliver(eta<0.1): {pre['n_sliver']:,} -> {post['n_sliver']:,} | "
          f"eta_min {pre['eta_min']:.4f} -> {post['eta_min']:.4f} | "
          f"edge_min {pre['edge_min']:.2f} -> {post['edge_min']:.2f} m")

    if not args.keep_medit:
        for p in (medit_in, medit_out):
            try:
                p.unlink()
            except OSError:
                pass

    if not fault_ok:
        print("FAULT DAMAGED by MMG -- NOT writing output (exit 2).", file=sys.stderr)
        return 2

    keep = new_triref > 0
    cells, physical = [("tetra", new_tets)], [np.ones(len(new_tets), int)]
    for tg in sorted(set(new_triref[keep].tolist())):
        sel = new_triref == tg
        cells.append(("triangle", new_tris[sel]))
        physical.append(np.full(int(sel.sum()), int(tg), dtype=int))
    out = meshio.Mesh(points=new_pts, cells=cells,
                      cell_data={"gmsh:physical": physical,
                                 "gmsh:geometrical": [p.copy() for p in physical]})
    meshio.write(str(args.outfile), out, file_format="gmsh22", binary=args.binary)
    print(f"wrote {args.outfile}; fault preserved: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
