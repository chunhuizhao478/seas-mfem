#!/usr/bin/env python3
"""Build a fully UNSTRUCTURED TPV102 mesh with equal volumes across the fault.

Goal: reduce the Eq.(18) volume ratio  r_v = max{V_A/V_B, V_B/V_A}  to 1 while
keeping a Delaunay tetrahedral (unstructured) mesh -- NOT a structured prism
strip.

Method ("half + mirror"):
  1. Mesh ONLY the y >= 0 half-space (x in [-X,X], y in [0,Ymax], z in [Zmin,0])
     with an unstructured 3D Delaunay mesh.  The fault footprint rectangle
     (x in [-FHL,FHL], z in [-FW,0]) is imprinted on the y=0 boundary face so
     its extent is exact and it can be refined.
  2. Reflect the half mesh across the fault plane y=0 to create the y<0 half.
     Nodes on y=0 are shared; every y>0 node gets a mirror partner at -y.
     Every tet gets a mirror tet (orientation flipped to keep positive volume).
  Because the -y half is the exact geometric mirror of the +y half, the two
  tets sharing any fault triangle are mirror partners with identical volume,
  so r_v = 1 to machine precision -- yet every element is an unstructured tet.

Physical groups (TPV102 driver defaults): 1 = free surface (z=0),
3 = fault (y=0 within footprint), 5 = absorbing (4 sides + bottom),
volume 1 = bulk.

Usage:
  python3 make_symmetric_unstructured_mesh.py <out.msh> [lc_fault] [lc_far]
"""
import sys
import math
from collections import defaultdict

# ---- domain / fault geometry (SCEC TPV102, metres) --------------------------
XMAX = 60e3
YMAX = 60e3
ZMIN = -60e3
FHL = 18e3          # fault half length along strike  -> x in [-18, 18] km
FW = 18e3           # fault width (depth)             -> z in [-18, 0] km

# nucleation patch (for extra refinement, matches tpv102_200m.geo)
X_NUCL = 0.0
W_NUCL = 7.5e3      # depth of hypocentre
R_NUCL = 3e3


def build_half_mesh(lc_fault, lc_far, lc_nucl):
    """Mesh the y>=0 half-space; return (coords[N,3], tets[M,4]) as 0-based."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("tpv102_half")
    occ = gmsh.model.occ

    # half box: y in [0, YMAX]
    box = occ.addBox(-XMAX, 0.0, ZMIN, 2 * XMAX, YMAX, -ZMIN)

    # fault footprint rectangle, built in xy-plane then rotated onto y=0.
    # addRectangle(x,y,z,dx,dy) lies in z=const plane; place at z=0 spanning
    # x in [-FHL,FHL], y in [-FW,0], then R_x(+pi/2): (x,y,0)->(x,0,+y) ...
    # we want z in [-FW,0], so rotate (x, y, 0) -> (x, 0, y): R_x(-pi/2).
    rect = occ.addRectangle(-FHL, -FW, 0.0, 2 * FHL, FW)
    occ.rotate([(2, rect)], 0, 0, 0, 1, 0, 0, math.pi / 2)  # (x,y,0)->(x,0,-y): z in [0,FW]?
    occ.synchronize()
    # check orientation of rect and fix if needed below via bbox query

    # imprint the rectangle onto the box (embed the fault surface in the y=0 face)
    occ.fragment([(3, box)], [(2, rect)])
    occ.synchronize()

    # locate the footprint surface on y=0 within the fault bbox
    eps = 1.0
    fault_surfs = gmsh.model.getEntitiesInBoundingBox(
        -FHL - eps, -eps, -FW - eps, FHL + eps, eps, FW + eps, 2)
    fault_tags = [t for (d, t) in fault_surfs]
    if not fault_tags:
        # rotation went the wrong way (z in [0,FW]); flip rectangle z-sign query
        fault_surfs = gmsh.model.getEntitiesInBoundingBox(
            -FHL - eps, -eps, -eps, FHL + eps, eps, FW + eps, 2)
        fault_tags = [t for (d, t) in fault_surfs]
    print(f"  footprint surface tag(s): {fault_tags}")

    # ---- size field: fine on/near fault, graded outward (matches .geo) -------
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", fault_tags)

    f_grade = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        f_grade, "F", f"0.05*F{f_dist} + (F{f_dist}/2.5e3)^2 + {lc_fault}")

    # nucleation patch refinement (ball around hypocentre on the fault)
    f_nuc_dist = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        f_nuc_dist, "F",
        f"sqrt((x-{X_NUCL})^2 + y^2 + (z+{W_NUCL})^2)")
    f_nuc = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_nuc, "InField", f_nuc_dist)
    gmsh.model.mesh.field.setNumber(f_nuc, "SizeMin", lc_nucl)
    gmsh.model.mesh.field.setNumber(f_nuc, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(f_nuc, "DistMin", R_NUCL)
    gmsh.model.mesh.field.setNumber(f_nuc, "DistMax", 3 * R_NUCL)

    f_min = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", [f_grade, f_nuc])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_min)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)     # Frontal-Delaunay 2D
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)   # Delaunay 3D
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

    gmsh.model.mesh.generate(3)

    # ---- extract nodes & tets ------------------------------------------------
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    ncoords = ncoords.reshape(-1, 3)
    tag2idx = {int(t): i for i, t in enumerate(ntags)}
    coords = ncoords.copy()

    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    tets = None
    for et, en in zip(etypes, enodes):
        if et == 4:  # 4-node tetra
            conn = en.reshape(-1, 4)
            tets = [[tag2idx[int(v)] for v in row] for row in conn]
    gmsh.finalize()
    if tets is None:
        raise RuntimeError("no tetrahedra produced")
    return coords, tets


def tet_volume(p0, p1, p2, p3):
    ax, ay, az = p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]
    bx, by, bz = p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]
    cx, cy, cz = p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]
    return (ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx)) / 6.0


def mirror_full(coords, tets, ytol=1e-6):
    """Reflect the y>=0 half mesh across y=0 -> full symmetric mesh."""
    n = len(coords)
    # mirror-node index for each original node (shared if on y=0)
    mirror_idx = [-1] * n
    full_coords = [tuple(c) for c in coords]
    for i in range(n):
        x, y, z = coords[i]
        if abs(y) <= ytol:
            mirror_idx[i] = i               # on fault plane -> shared
        else:
            mirror_idx[i] = len(full_coords)
            full_coords.append((x, -y, z))

    full_tets = []
    for t in tets:
        # ensure positive volume for the +y tet
        p = [full_coords[k] for k in t]
        if tet_volume(*p) < 0:
            t = [t[0], t[2], t[1], t[3]]
        full_tets.append(t)
        # mirror tet: reflect node ids; reflection flips orientation -> swap 2
        m = [mirror_idx[k] for k in t]
        m = [m[0], m[2], m[1], m[3]]
        full_tets.append(m)
    return full_coords, full_tets


def classify_and_write(coords, tets, out, ytol=1e-6, gtol=1.0):
    """Build boundary/fault physical groups geometrically and write msh22."""
    # face -> list of tet ids (each tet contributes its 4 triangular faces)
    face2tet = defaultdict(list)
    FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for ti, t in enumerate(tets):
        for f in FACES:
            key = tuple(sorted((t[f[0]], t[f[1]], t[f[2]])))
            face2tet[key].append(ti)

    free_tris, abs_tris, fault_tris = [], [], []
    for key, owners in face2tet.items():
        a, b, c = (coords[key[0]], coords[key[1]], coords[key[2]])
        cy = (a[1] + b[1] + c[1]) / 3.0
        cz = (a[2] + b[2] + c[2]) / 3.0
        cx = (a[0] + b[0] + c[0]) / 3.0
        if len(owners) == 1:                       # boundary face
            if abs(a[2]) <= gtol and abs(b[2]) <= gtol and abs(c[2]) <= gtol:
                free_tris.append(key)              # z = 0 free surface
            else:
                abs_tris.append(key)               # outer box: sides + bottom
        else:                                      # internal face
            on_y0 = (abs(a[1]) <= ytol and abs(b[1]) <= ytol and abs(c[1]) <= ytol)
            in_fp = (-FHL - gtol <= cx <= FHL + gtol) and (-FW - gtol <= cz <= gtol)
            if on_y0 and in_fp:
                fault_tris.append(key)             # fault footprint

    # ---- write Gmsh 2.2 ascii ------------------------------------------------
    with open(out, "w") as o:
        o.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        o.write("$PhysicalNames\n4\n")
        o.write('2 1 "free"\n2 3 "fault"\n2 5 "absorb"\n3 1 "rock"\n')
        o.write("$EndPhysicalNames\n")
        o.write("$Nodes\n")
        o.write(f"{len(coords)}\n")
        for i, (x, y, z) in enumerate(coords, 1):
            o.write(f"{i} {x:.10g} {y:.10g} {z:.10g}\n")
        o.write("$EndNodes\n")
        o.write("$Elements\n")
        ntri = len(free_tris) + len(fault_tris) + len(abs_tris)
        o.write(f"{ntri + len(tets)}\n")
        eid = 0
        for phys, group in ((1, free_tris), (3, fault_tris), (5, abs_tris)):
            for k in group:
                eid += 1
                o.write(f"{eid} 2 2 {phys} {phys} {k[0]+1} {k[1]+1} {k[2]+1}\n")
        for t in tets:
            eid += 1
            o.write(f"{eid} 4 2 1 1 {t[0]+1} {t[1]+1} {t[2]+1} {t[3]+1}\n")
        o.write("$EndElements\n")

    return dict(free=len(free_tris), fault=len(fault_tris),
                absorb=len(abs_tris), tets=len(tets), nodes=len(coords))


def report_rv(coords, tets, fault_keys=None):
    """Compute Eq.(18) r_v stats from the full mesh."""
    face2vol = defaultdict(list)
    FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    # identify fault faces: internal faces on y=0 within footprint
    face2tet = defaultdict(list)
    for ti, t in enumerate(tets):
        v = abs(tet_volume(*[coords[k] for k in t]))
        for f in FACES:
            key = tuple(sorted((t[f[0]], t[f[1]], t[f[2]])))
            face2tet[key].append(v)
    ratios = []
    for key, vols in face2tet.items():
        a, b, c = coords[key[0]], coords[key[1]], coords[key[2]]
        on_y0 = abs(a[1]) <= 1e-6 and abs(b[1]) <= 1e-6 and abs(c[1]) <= 1e-6
        cx = (a[0]+b[0]+c[0])/3.0
        cz = (a[2]+b[2]+c[2])/3.0
        in_fp = (-FHL-1 <= cx <= FHL+1) and (-FW-1 <= cz <= 1)
        if on_y0 and in_fp and len(vols) == 2:
            ratios.append(max(vols[0]/vols[1], vols[1]/vols[0]))
    ratios.sort()
    n = len(ratios)
    import statistics
    return dict(n=n, mn=ratios[0], mean=statistics.fmean(ratios),
                med=ratios[n//2], p95=ratios[int(0.95*(n-1))],
                p99=ratios[int(0.99*(n-1))], mx=ratios[-1],
                gt15=sum(1 for r in ratios if r > 1.5),
                gt20=sum(1 for r in ratios if r > 2.0))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "tpv102_sym.msh"
    lc_fault = float(sys.argv[2]) if len(sys.argv) > 2 else 200.0
    lc_far = float(sys.argv[3]) if len(sys.argv) > 3 else 5000.0
    lc_nucl = min(lc_fault, 200.0)

    print(f"[1/4] meshing y>=0 half  (lc_fault={lc_fault}, lc_far={lc_far}) ...")
    coords, tets = build_half_mesh(lc_fault, lc_far, lc_nucl)
    print(f"      half mesh: {len(coords):,} nodes, {len(tets):,} tets")

    print("[2/4] mirroring across y=0 ...")
    fcoords, ftets = mirror_full(coords, tets)
    print(f"      full mesh: {len(fcoords):,} nodes, {len(ftets):,} tets")

    print(f"[3/4] tagging + writing {out} ...")
    info = classify_and_write(fcoords, ftets, out)
    print(f"      free={info['free']:,}  fault={info['fault']:,}  "
          f"absorb={info['absorb']:,}  tets={info['tets']:,}")

    print("[4/4] Eq.(18) r_v on the new fault ...")
    s = report_rv(fcoords, ftets)
    print(f"      fault faces        : {s['n']:,}")
    print(f"      r_v min/mean/median: {s['mn']:.6f} / {s['mean']:.6f} / {s['med']:.6f}")
    print(f"      r_v 95/99/MAX      : {s['p95']:.6f} / {s['p99']:.6f} / {s['mx']:.6f}")
    print(f"      r_v > 1.5          : {s['gt15']:,} ({100.0*s['gt15']/s['n']:.2f}%)")
    print(f"      r_v > 2.0          : {s['gt20']:,} ({100.0*s['gt20']/s['n']:.2f}%)")


if __name__ == "__main__":
    main()
