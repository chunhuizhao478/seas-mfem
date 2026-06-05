#!/usr/bin/env python3
"""Export the Eq.(18) volume-ratio r_v as a spatial field on the fault surface.

Eq. (18) (Zhang et al., 2023, JGR Solid Earth, 10.1029/2022JB025817):

    r_v = max{ V_A/V_B , V_B/V_A }

where A and B are the two tetrahedra sharing a fault triangle, V_A/V_B their
volumes.  This script writes one r_v value per fault triangle into a VTU
(VTK unstructured grid) so the distribution can be visualised on the fault.

Usage:
    export_eq18_fault_vtu.py <mesh.msh> <fault_phys> <rock_phys> <out.vtu>
"""
import sys
from collections import defaultdict

TRI = 2   # Gmsh 3-node triangle
TET = 4   # Gmsh 4-node tetrahedron


def tet_volume(p0, p1, p2, p3):
    ax, ay, az = p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]
    bx, by, bz = p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]
    cx, cy, cz = p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]
    det = (ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx))
    return abs(det) / 6.0


def main(path, fault_phys, rock_phys, out):
    nodes = {}
    fault_tris = []           # ordered list of (n0,n1,n2) fault triangles
    fault_face_set = set()    # frozenset keys for O(1) membership
    tets = []

    with open(path) as f:
        line = f.readline()
        while line:
            tok = line.strip()
            if tok == "$Nodes":
                count = int(f.readline())
                for _ in range(count):
                    p = f.readline().split()
                    nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
                f.readline()
            elif tok == "$Elements":
                count = int(f.readline())
                for _ in range(count):
                    p = f.readline().split()
                    etype = int(p[1]); ntags = int(p[2]); phys = int(p[3])
                    conn = list(map(int, p[3+ntags:]))
                    if etype == TRI and phys == fault_phys:
                        fault_tris.append(tuple(conn))
                        fault_face_set.add(frozenset(conn))
                    elif etype == TET and phys == rock_phys:
                        tets.append(conn)
                f.readline()
                break
            line = f.readline()

    # Map each fault face -> list of adjacent tet volumes.
    face_to_vols = defaultdict(list)
    for conn in tets:
        v0, v1, v2, v3 = conn
        vol = tet_volume(nodes[v0], nodes[v1], nodes[v2], nodes[v3])
        for face in ((v0, v1, v2), (v0, v1, v3), (v0, v2, v3), (v1, v2, v3)):
            key = frozenset(face)
            if key in fault_face_set:
                face_to_vols[key].append(vol)

    # Per-triangle scalar fields, in the order fault_tris are listed.
    import math
    rv = []; vmin = []; vmax = []; bad = 0
    for tri in fault_tris:
        vols = face_to_vols.get(frozenset(tri), [])
        if len(vols) == 2:
            a, b = vols
            rv.append(max(a/b, b/a))
            vmin.append(min(a, b)); vmax.append(max(a, b))
        else:
            rv.append(0.0); vmin.append(0.0); vmax.append(0.0); bad += 1

    # Compact point set: only nodes used by the fault triangles, remapped 0..N-1.
    used = sorted({n for tri in fault_tris for n in tri})
    remap = {nid: i for i, nid in enumerate(used)}

    npts = len(used)
    ncells = len(fault_tris)
    with open(out, "w") as o:
        o.write('<?xml version="1.0"?>\n')
        o.write('<VTKFile type="UnstructuredGrid" version="0.1" '
                'byte_order="LittleEndian">\n')
        o.write('  <UnstructuredGrid>\n')
        o.write(f'    <Piece NumberOfPoints="{npts}" NumberOfCells="{ncells}">\n')
        # points
        o.write('      <Points>\n')
        o.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        o.write(' '.join(f'{nodes[n][0]:.6g} {nodes[n][1]:.6g} {nodes[n][2]:.6g}'
                         for n in used))
        o.write('\n        </DataArray>\n      </Points>\n')
        # cells
        o.write('      <Cells>\n')
        o.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        o.write(' '.join(str(remap[n]) for tri in fault_tris for n in tri))
        o.write('\n        </DataArray>\n')
        o.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        o.write(' '.join(str(3*(i+1)) for i in range(ncells)))
        o.write('\n        </DataArray>\n')
        o.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        o.write(' '.join('5' for _ in range(ncells)))   # 5 = VTK_TRIANGLE
        o.write('\n        </DataArray>\n      </Cells>\n')
        # cell data
        o.write('      <CellData Scalars="rv">\n')
        for name, data in (("rv", rv),
                           ("log10_rv", [math.log10(x) if x > 0 else 0.0 for x in rv]),
                           ("vol_min", vmin), ("vol_max", vmax)):
            o.write(f'        <DataArray type="Float64" Name="{name}" format="ascii">\n')
            o.write(' '.join(f'{x:.6g}' for x in data))
            o.write('\n        </DataArray>\n')
        o.write('      </CellData>\n')
        o.write('    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n')

    finite = [x for x in rv if x > 0]
    print(f"wrote {out}")
    print(f"  fault triangles : {ncells:,}  (points {npts:,})")
    print(f"  faces w/o 2 tets: {bad}")
    if finite:
        print(f"  r_v  min/mean/max : {min(finite):.4f} / "
              f"{sum(finite)/len(finite):.4f} / {max(finite):.4f}")
        # locate the worst triangle
        imax = max(range(ncells), key=lambda i: rv[i])
        c = [sum(nodes[n][k] for n in fault_tris[imax])/3.0 for k in range(3)]
        print(f"  max r_v={rv[imax]:.4f} at fault centroid "
              f"(x={c[0]:.0f}, y={c[1]:.0f}, z={c[2]:.0f})")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
