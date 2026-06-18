#!/usr/bin/env python3
"""
make_station_vtu.py — build ParaView VTU(s) that highlight the SAFS station locations.

Reads the committed station files written by make_stations.py:
  - safs_receivers.annotated.txt   12 off-fault free-surface receivers (z=0)
  - safs_pickpoints.annotated.txt   9 on-fault pickpoints (snapped to fault facets)
and writes:
  - safs_stations.vtu        point cloud of ALL stations, VTK_VERTEX cells, with point-data
                             fields {id, kind, group_id, depth_km} for coloring/labeling.
  - safs_stations_legend.csv id,label,kind,group,x,y,z   (cross-reference for the `id` labels)
  - safs_fault_trace.vtu     (optional, if safs_mesh.puml.h5 is present) the surface fault-trace
                             vertices (z~0) as a point cloud, for spatial context.

Coordinates are UTM 11N meters — the SAME frame as safs-surface.xdmf / safs-GME-surface.xdmf,
so the station VTU overlays those directly in ParaView.

In ParaView: open safs_stations.vtu -> Glyph filter (Sphere, scale ~1-2 km) to make the points
visible -> color by `kind` (0=receiver, 1=pickpoint) or `group_id` -> enable Point Labels by `id`
and consult safs_stations_legend.csv for the names.

Deps: numpy (always); h5py only for the optional fault-trace (skipped if unavailable/missing mesh).
Usage: python3 make_station_vtu.py   [--mesh safs_mesh.puml.h5]   [--no-trace]
"""
import argparse
import os
import re
import numpy as np

# group label-prefix -> (group_id, kind)   kind: 0=receiver(off-fault), 1=pickpoint(on-fault)
GROUPS = {
    "epicenter":      (0, 0),
    "fault_normal":   (1, 0),
    "along_strike":   (2, 0),
    "hypocenter":     (3, 1),
    "downdip":        (4, 1),
    "strike":         (5, 1),
}
GROUP_NAMES = {0: "epicenter", 1: "fault_normal", 2: "along_strike",
               3: "hypocenter", 4: "downdip", 5: "strike"}


def classify(label):
    for prefix, (gid, kind) in GROUPS.items():
        if label.startswith(prefix):
            return gid, kind
    return -1, -1


def parse_annotated(path):
    """Yield (label, x, y, z) from a make_stations *.annotated.txt file."""
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s.startswith("#"):
                continue
            body = s[1:].strip()
            toks = body.split()
            if len(toks) < 4:
                continue
            # last three tokens must be floats; the rest is the label
            try:
                x, y, z = float(toks[-3]), float(toks[-2]), float(toks[-1])
            except ValueError:
                continue
            label = " ".join(toks[:-3])
            if label.lower() == "label":          # header row "label x y z"
                continue
            rows.append((label, x, y, z))
    return rows


def write_points_vtu(path, pts, point_data):
    """Write a VTU of isolated points (VTK_VERTEX, type=1). point_data: name -> (np.array, vtk_type)."""
    n = len(pts)
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" '
                'byte_order="LittleEndian">\n  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{n}" NumberOfCells="{n}">\n')
        # points
        f.write('      <Points>\n        <DataArray type="Float64" NumberOfComponents="3" '
                'format="ascii">\n')
        f.write("          " + " ".join(f"{v:.4f}" for v in np.asarray(pts).ravel()) + "\n")
        f.write("        </DataArray>\n      </Points>\n")
        # cells: one vertex per point
        f.write("      <Cells>\n")
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n          '
                + " ".join(str(i) for i in range(n)) + "\n        </DataArray>\n")
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n          '
                + " ".join(str(i + 1) for i in range(n)) + "\n        </DataArray>\n")
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n          '
                + " ".join("1" for _ in range(n)) + "\n        </DataArray>\n")
        f.write("      </Cells>\n")
        # point data
        if point_data:
            f.write("      <PointData>\n")
            for name, (arr, vtktype) in point_data.items():
                arr = np.asarray(arr)
                fmt = (lambda v: f"{int(v)}") if vtktype.startswith("Int") else (lambda v: f"{v:.4f}")
                f.write(f'        <DataArray type="{vtktype}" Name="{name}" format="ascii">\n          '
                        + " ".join(fmt(v) for v in arr) + "\n        </DataArray>\n")
            f.write("      </PointData>\n")
        f.write("    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n")


def build_stations(recv_path, pick_path, out_vtu, out_legend):
    rows = []
    for path in (recv_path, pick_path):
        if not os.path.exists(path):
            print(f"WARNING: {path} not found — skipping")
            continue
        rows += parse_annotated(path)
    if not rows:
        raise SystemExit("no station rows parsed; check the annotated.txt files")

    pts, ids, kinds, gids, depth_km = [], [], [], [], []
    for i, (label, x, y, z) in enumerate(rows):
        gid, kind = classify(label)
        pts.append((x, y, z))
        ids.append(i)
        kinds.append(kind)
        gids.append(gid)
        depth_km.append(-z / 1000.0)

    write_points_vtu(out_vtu, pts, {
        "id":       (np.array(ids), "Int32"),
        "kind":     (np.array(kinds), "Int32"),       # 0=receiver, 1=pickpoint
        "group_id": (np.array(gids), "Int32"),
        "depth_km": (np.array(depth_km), "Float64"),
    })

    with open(out_legend, "w") as f:
        f.write("id,label,kind,group,x,y,z\n")
        for i, (label, x, y, z) in enumerate(rows):
            gid, kind = classify(label)
            kname = "receiver" if kind == 0 else "pickpoint"
            f.write(f"{i},{label},{kname},{GROUP_NAMES.get(gid,'?')},{x:.4f},{y:.4f},{z:.4f}\n")

    nr = sum(1 for k in kinds if k == 0)
    npk = sum(1 for k in kinds if k == 1)
    print(f"wrote {out_vtu}  ({len(rows)} stations: {nr} receivers + {npk} pickpoints)")
    print(f"wrote {out_legend}")


def build_trace(mesh_path, out_vtu, fault_bc=3, z_tol=1.0, epicenter=(606971.0, 3707270.0),
                radius=60000.0):
    try:
        import h5py
    except ImportError:
        print("h5py not available — skipping fault-trace VTU")
        return
    if not os.path.exists(mesh_path):
        print(f"{mesh_path} not found — skipping fault-trace VTU")
        return
    with h5py.File(mesh_path, "r") as f:
        geom = f["geometry"][:]
        connect = f["connect"][:].astype(np.int64)
        boundary = f["boundary"][:].astype(np.int64)
    face_sets = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]])
    tris = []
    for j in range(4):
        sel = np.where(((boundary >> (8 * j)) & 0xFF) == fault_bc)[0]
        if sel.size:
            tris.append(connect[sel][:, face_sets[j]])
    if not tris:
        print(f"no fault faces (BC={fault_bc}) — skipping trace")
        return
    fault_vids = np.unique(np.vstack(tris).ravel())
    fv = geom[fault_vids]
    on_surf = fv[np.abs(fv[:, 2]) < z_tol]
    near = on_surf[np.linalg.norm(on_surf[:, :2] - np.array(epicenter), axis=1) < radius]
    trace = np.unique(np.round(near, 3), axis=0)
    if trace.shape[0] == 0:
        print("no surface-trace vertices found — skipping trace")
        return
    write_points_vtu(out_vtu, trace, {"on_trace": (np.ones(len(trace)), "Int32")})
    print(f"wrote {out_vtu}  ({trace.shape[0]} fault-trace surface vertices)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receivers", default="safs_receivers.annotated.txt")
    ap.add_argument("--pickpoints", default="safs_pickpoints.annotated.txt")
    ap.add_argument("--mesh", default="safs_mesh.puml.h5")
    ap.add_argument("--out", default="safs_stations.vtu")
    ap.add_argument("--legend", default="safs_stations_legend.csv")
    ap.add_argument("--trace-out", default="safs_fault_trace.vtu")
    ap.add_argument("--no-trace", action="store_true")
    a = ap.parse_args()
    build_stations(a.receivers, a.pickpoints, a.out, a.legend)
    if not a.no_trace:
        build_trace(a.mesh, a.trace_out)


if __name__ == "__main__":
    main()
