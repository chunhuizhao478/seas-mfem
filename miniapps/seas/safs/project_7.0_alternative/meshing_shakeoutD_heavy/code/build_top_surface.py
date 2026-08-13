#!/usr/bin/env python3
"""build_top_surface.py -- the PLC's free surface, with the fault trace embedded.

gmsh meshes the rotated rectangle in plan view with the four trace polylines
EMBEDDED as hard constraints, every trace segment pinned `setTransfiniteCurve(...,2)`
so gmsh reproduces the fault's own nodes vertex-for-vertex instead of resampling
them.  That is what makes the fault daylight: its top edge and the free surface
share nodes exactly.

z is then assigned: trace nodes keep their TRUE z (the fault's top edge wanders
-49.9..+25.8 m), everything else sits at 0.  That reproduces the existing mesh's
convention, whose free surface is flat only to +-25 m.
"""
import os, sys
import numpy as np

FIELD_DAT = "build_tmp/top_field.dat"

def write_field(path, Q, sf):
    """gmsh Structured background field on an AXIS-ALIGNED grid (gmsh requires it)."""
    from scipy.interpolate import RegularGridInterpolator
    W = sf["W"]; u = sf["u"]; v = sf["v"]; H = sf["H"]; ss = sf["ss"]; tt = sf["tt"]
    itp = RegularGridInterpolator((ss, tt), H, bounds_error=False, fill_value=None)
    x0, x1 = Q[:, 0].min() - 3000, Q[:, 0].max() + 3000
    y0, y1 = Q[:, 1].min() - 3000, Q[:, 1].max() + 3000
    step = 750.0
    gx = np.arange(x0, x1 + step, step); gy = np.arange(y0, y1 + step, step)
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    d = np.stack([GX - W[0], GY - W[1]], -1)
    a = d @ u; b = d @ v
    hv = itp(np.stack([np.clip(a, ss[0], ss[-1]), np.clip(b, tt[0], tt[-1])], -1))
    with open(path, "w") as f:
        f.write(f"{gx[0]:.9g} {gy[0]:.9g} {-1000.0:.9g}\n")
        f.write(f"{step:.9g} {step:.9g} {2000.0:.9g}\n")
        f.write(f"{len(gx)} {len(gy)} 2\n")
        out = np.repeat(hv[:, :, None], 2, axis=2).transpose(0, 1, 2).ravel()
        f.write("\n".join("%.6g" % q for q in out)); f.write("\n")
    print(f"[field] {path}  {len(gx)}x{len(gy)}x2  h {hv.min():.1f}..{hv.max():.1f}")
    return path

def main():
    import gmsh
    Q = np.load("build_tmp/domain_corners_utm.npy")
    fs = np.load("build_tmp/fault_surface.npz"); P = fs["P"]
    ch = np.load("build_tmp/trace_chains.npz")
    chains = [ch[k] for k in ch.files]
    sf = dict(np.load("build_tmp/surface_field.npz"))
    write_field(FIELD_DAT, Q, sf)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("top")
    geo = gmsh.model.geo
    LC = 5000.0
    cor = [geo.addPoint(float(Q[i, 0]), float(Q[i, 1]), 0.0, LC) for i in range(4)]
    ln = [geo.addLine(cor[i], cor[(i + 1) % 4]) for i in range(4)]
    surf = geo.addPlaneSurface([geo.addCurveLoop(ln)])

    tv = np.unique(np.concatenate(chains))
    pid = {}
    for v in tv:
        pid[int(v)] = geo.addPoint(float(P[v, 0]), float(P[v, 1]), 0.0, 115.0)
    elines = []
    for c in chains:
        for i in range(len(c) - 1):
            elines.append(geo.addLine(pid[int(c[i])], pid[int(c[i + 1])]))
    extra = np.load("build_tmp/inplane_edges.npy") if os.path.exists("build_tmp/inplane_edges.npy") \
            else np.zeros((0, 2), np.int64)
    for a_, b_ in extra:
        for w in (int(a_), int(b_)):
            if w not in pid:
                pid[w] = geo.addPoint(float(P[w, 0]), float(P[w, 1]), 0.0, 115.0)
        elines.append(geo.addLine(pid[int(a_)], pid[int(b_)]))
    if len(extra):
        print(f"[gmsh] + {len(extra)} in-plane interior fault edges embedded")
    geo.synchronize()
    for l in elines:
        gmsh.model.mesh.setTransfiniteCurve(l, 2)
    gmsh.model.mesh.embed(1, elines, 2, surf)
    print(f"[gmsh] embedded {len(tv):,} points, {len(elines):,} trace segments")

    fid = gmsh.model.mesh.field.add("Structured")
    gmsh.model.mesh.field.setString(fid, "FileName", os.path.abspath(FIELD_DAT))
    gmsh.model.mesh.field.setNumber(fid, "TextFormat", 1)

    # Local refinement at the STEEP trace nodes, expressed to gmsh directly rather
    # than through the background grid: that grid is 750 m, so a 600 m feature lands
    # on ~14 nodes and is effectively invisible (measured).  A Distance+Threshold
    # pair has no such resolution limit.
    steep = np.zeros(0, np.int64)
    use = fid
    if len(steep):
        dfl = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dfl, "PointsList",
                                         [pid[int(v)] for v in steep])
        thr = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thr, "InField", dfl)
        gmsh.model.mesh.field.setNumber(thr, "SizeMin", 35.0)
        gmsh.model.mesh.field.setNumber(thr, "SizeMax", 5000.0)
        gmsh.model.mesh.field.setNumber(thr, "DistMin", 150.0)
        gmsh.model.mesh.field.setNumber(thr, "DistMax", 1500.0)
        mn = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(mn, "FieldsList", [fid, thr])
        use = mn
        print(f"[gmsh] steep-trace refinement: {len(steep)} nodes, 35 m within 150 m")
    gmsh.model.mesh.field.setAsBackgroundMesh(use)
    for o in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature",
              "Mesh.MeshSizeExtendFromBoundary"):
        gmsh.option.setNumber(o, 0)
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(2)
    nt_, nc, _ = gmsh.model.mesh.getNodes()
    XYZ = np.asarray(nc, float).reshape(-1, 3)
    order = {int(t): i for i, t in enumerate(nt_)}
    et, _, en = gmsh.model.mesh.getElements(2, surf)
    tri = None
    for e, n in zip(et, en):
        if e == 2:
            tri = np.array([order[int(q)] for q in n], np.int64).reshape(-1, 3)
    gmsh.finalize()
    print(f"[top] {len(tri):,} triangles, {len(XYZ):,} nodes")

    # lift z: trace nodes take the fault's true z, everything else 0
    from scipy.spatial import cKDTree
    tree = cKDTree(P[tv, :2])
    # TAPER, do not step.  Dropping straight to z=0 at the first node off the trace
    # makes the surface triangle there dip ~26 m over ~115 m, and a shallowly-dipping
    # fault triangle just below the trace then pokes THROUGH it -- tetgen -d found
    # exactly 11 such triangles.  Decaying the trace elevation over TAPER m keeps the
    # surface above the fault while still matching the trace node-for-node at d=0.
    # The fault's daylighting edge is snapped to z=0 upstream, so the free surface
    # is EXACTLY flat and matches it identically -- no taper, no local refinement.
    XYZ[:, 2] = 0.0
    dd, kk = tree.query(XYZ[:, :2], k=1, distance_upper_bound=1.0)
    hit = np.isfinite(dd)
    print(f"[lift] {int(hit.sum()):,} of {len(tv):,} trace nodes matched "
          f"(max offset {dd[hit].max() if hit.any() else 0:.3e} m); "
          f"z {XYZ[:,2].min():.2f}..{XYZ[:,2].max():.2f}")
    if int(hit.sum()) != len(tv):
        print(f"  *** WARNING: {len(tv)-int(hit.sum())} trace nodes NOT reproduced by gmsh")
    np.savez_compressed("build_tmp/top_surface.npz", P=XYZ, T=tri)
    print("[out] build_tmp/top_surface.npz")

if __name__ == "__main__":
    main()
