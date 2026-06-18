#!/usr/bin/env python3
"""
make_stations.py — SAFS station generator for the LSW post-processing activity.

Reads the fault geometry from the PUML mesh (safs_mesh.puml.h5) and writes:
  - safs_receivers.dat   off-fault, free-surface (z=0) seismogram stations (6.1 / 6.4)
  - safs_pickpoints.dat  ON-fault pickpoints, snapped to the nearest fault facet (6.3 / 6.4)

Both .dat files are COORDS-ONLY ("x y z" per line, meters, UTM 11N) because SeisSol's
receiver/pickpoint parsers split each non-whitespace line on whitespace and float-cast every
token — a "#" comment line would crash them. Human-readable labels go to the *.annotated.txt
companions instead.

Method (see POSTPROCESS_ACTIVITY_PLAN.md Phase 5 "How to pick station points", Route 1):
  OFF-FAULT receivers may sit anywhere (SeisSol finds the containing tet), so they are placed
  RELATIVE to the surface fault trace: an epicentral station, a fault-NORMAL profile, and an
  along-strike (directivity) line offset ~1.5 km off the trace toward NW.
  ON-FAULT pickpoints MUST lie on the curved fault, so each target (along-strike s, depth d) is
  SNAPPED to the nearest fault-facet centroid via a KD-tree — never hand-typed.

PUML decoding (verified against SeisSol source):
  - boundary[i] packs the 4 face BCs as 8 bits each: bc(face j) = (boundary[i] >> (8*j)) & 0xff
    (src/Geometry/PUMLReader.h:33). Fault BC = 3 (free surface = 1, absorbing = 5).
  - the face->vertex map for the file's boundary index is the PUML Numbering
    (submodules/PUML/Numbering.h:43) {{1,0,2},{0,1,3},{1,2,3},{2,0,3}}, i.e. the vertex SETS
    face0={0,1,2} face1={0,1,3} face2={1,2,3} face3={0,2,3}. (Sets are what matter for a
    centroid; this differs from SeisSol's internal MeshTools::FACE2NODES at faces 2/3.)

Dependencies: numpy, h5py, scipy (cKDTree). Run in an env that has them (e.g. `conda run -n
pythonenv python make_stations.py`). Run ONCE before submitting and commit the generated .dat files.
"""

import sys
import numpy as np
import h5py
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Configuration (hardcoded per the plan; edit here to retarget)
# ---------------------------------------------------------------------------
MESH = "safs_mesh.puml.h5"
EPICENTER = np.array([606971.0, 3707270.0])              # (x, y) surface point above hypocenter
HYPOCENTER = np.array([606971.0, 3707270.0, -4965.62])   # on-fault source (validation anchor)

FAULT_BC = 3                  # PUML boundary tag for the dynamic-rupture fault
SURFACE_Z_TOL = 1.0          # |z| < tol [m] selects free-surface (trace) vertices (geom z==0 exact)
TRACE_RADIUS = 60000.0       # keep trace vertices within this map distance of the epicenter [m]
MARCH_MAXSTEP = 2500.0       # max gap when marching along the trace [m] (prevents strand jumps)

# off-fault receiver geometry
FAULT_NORMAL_OFFSETS = [-20000.0, -10000.0, -5000.0, -2000.0,
                         2000.0, 5000.0, 10000.0, 20000.0]   # [m], +/- across the fault
ALONG_STRIKE_DIST = [10000.0, 20000.0, 40000.0]              # NW distances along the trace [m]
ALONG_STRIKE_OFFSET = 1500.0                                 # push receivers off the trace [m]

# on-fault pickpoint targets. Because the fault DIPS, a pickpoint is chosen as the fault facet
# whose MAP position is near the target column AND whose depth is closest to the target depth
# (a plain 3D nearest-neighbour would drift off the dipping plane). Both guarantee an on-fault point.
DOWNDIP_DEPTHS = [2000.0, 5000.0, 8000.0, 12000.0]          # [m] column near the epicenter map pos.
DOWNDIP_MAP_R = 8000.0                                      # [m] map radius for the downdip column
ALONG_STRIKE_PP_DIST = [0.0, 10000.0, 20000.0, 40000.0]     # NW distances at seismogenic depth [m]
ALONG_STRIKE_PP_DEPTH = 7000.0                              # [m]
ALONG_STRIKE_MAP_R = 5000.0                                 # [m] map radius for along-strike picks

HYPO_SNAP_TOL = 800.0        # validation: nearest fault centroid to the hypocenter must be < this


def read_fault_geometry(mesh_path):
    """Return (vertices Nx3, fault_centroids Kx3, fault_vertex_ids unique array)."""
    with h5py.File(mesh_path, "r") as f:
        geom = f["geometry"][:]                 # (Nv, 3) float64
        connect = f["connect"][:].astype(np.int64)   # (Ne, 4)
        boundary = f["boundary"][:].astype(np.int64) # (Ne,)
    # PUML face -> vertex sets (Numbering.h)
    face_sets = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]])
    tris = []
    for j in range(4):
        bcj = (boundary >> (8 * j)) & 0xFF
        sel = np.where(bcj == FAULT_BC)[0]
        if sel.size:
            tris.append(connect[sel][:, face_sets[j]])
    if not tris:
        raise RuntimeError(f"No fault faces (BC={FAULT_BC}) found in {mesh_path}")
    tris = np.vstack(tris)                       # (Nf, 3) global vertex ids
    centroids = geom[tris].mean(axis=1)          # (Nf, 3)
    fault_vids = np.unique(tris.ravel())
    return geom, centroids, fault_vids


def extract_trace(geom, fault_vids):
    """Surface fault-trace vertices (z~0), restricted to a neighborhood of the epicenter."""
    fv = geom[fault_vids]
    on_surf = np.abs(fv[:, 2]) < SURFACE_Z_TOL
    trace = fv[on_surf][:, :2]                   # (Nt, 2) xy
    if trace.shape[0] < 3:
        raise RuntimeError("Fewer than 3 fault-trace (z~0) vertices found.")
    near = np.linalg.norm(trace - EPICENTER, axis=1) < TRACE_RADIUS
    trace = trace[near]
    # dedupe (mesh vertices on the trace are unique already, but be safe)
    trace = np.unique(np.round(trace, 3), axis=0)
    return trace


def local_tangent(trace, p, radius=3000.0):
    """Principal (along-strike) direction of trace points within `radius` of point p (xy)."""
    d = np.linalg.norm(trace - p, axis=1)
    nb = trace[d < radius]
    if nb.shape[0] < 2:
        nb = trace[np.argsort(d)[:5]]
    c = nb - nb.mean(axis=0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    t = vt[0]
    return t / np.linalg.norm(t)


def nw_sign(t):
    """Return +/-1 so that sign*t points NW (x decreasing, y increasing)."""
    score = (-t[0] + t[1])
    return 1.0 if score >= 0 else -1.0


def march_trace(trace, start_xy, heading, targets, maxstep=MARCH_MAXSTEP):
    """March along the trace from start_xy in `heading`, recording (xy, heading, dist) when
    each cumulative distance in `targets` is first reached. Follows the curve via a forward-aligned
    nearest-vertex step with a max-gap guard (avoids jumping to a different fault strand)."""
    targets = sorted(targets)
    out = []
    cur = start_xy.astype(float).copy()
    head = heading / np.linalg.norm(heading)
    dist = 0.0
    ti = 0
    for _ in range(10000):
        if ti >= len(targets):
            break
        d = trace - cur
        dd = np.linalg.norm(d, axis=1)
        fwd = d @ head
        cand = np.where((dd > 1.0) & (dd < maxstep) & (fwd > 0.0))[0]
        if cand.size == 0:
            break                                # strand ended before reaching remaining targets
        j = cand[np.argmax(fwd[cand] / dd[cand])]  # most forward-aligned next vertex
        step = dd[j]
        head = (trace[j] - cur) / step           # follow the local curvature
        cur = trace[j].astype(float).copy()
        dist += step
        while ti < len(targets) and dist >= targets[ti]:
            out.append((cur.copy(), head.copy(), targets[ti]))
            ti += 1
    return out


def pick_at_depth(centroids, xy, depth, map_radius):
    """Return the index of the fault centroid within `map_radius` (xy) of `xy` whose depth is
    closest to -`depth`, or None if no facet lies within the radius. Keeps the point ON the dipping
    fault at the requested depth instead of drifting off-plane like a 3D nearest-neighbour would."""
    mapd = np.linalg.norm(centroids[:, :2] - xy, axis=1)
    near = np.where(mapd < map_radius)[0]
    if near.size == 0:
        return None
    return int(near[np.argmin(np.abs(centroids[near, 2] - (-depth)))])


def main():
    mesh = sys.argv[1] if len(sys.argv) > 1 else MESH
    geom, centroids, fault_vids = read_fault_geometry(mesh)
    print(f"fault facets: {centroids.shape[0]:,}   fault vertices: {fault_vids.size:,}")

    tree = cKDTree(centroids)
    dH, iH = tree.query(HYPOCENTER)
    print(f"validation: nearest fault centroid to the hypocenter = {dH:.1f} m")
    if dH > HYPO_SNAP_TOL:
        raise RuntimeError(
            f"Hypocenter is {dH:.1f} m from the nearest fault centroid (> {HYPO_SNAP_TOL} m). "
            "The PUML fault decode is likely wrong (face map / BC tag) — do NOT ship these files.")

    trace = extract_trace(geom, fault_vids)
    print(f"trace vertices (z~0, within {TRACE_RADIUS/1000:.0f} km): {trace.shape[0]:,}")

    # epicenter station = trace vertex nearest the epicenter
    p0 = trace[np.argmin(np.linalg.norm(trace - EPICENTER, axis=1))]
    t0 = local_tangent(trace, p0)
    n0 = np.array([-t0[1], t0[0]])               # rotate tangent 90 deg -> fault normal (xy)
    head_nw = nw_sign(t0) * t0

    # ---- off-fault receivers --------------------------------------------------
    receivers = []   # (label, x, y, z)
    receivers.append(("epicenter", p0[0], p0[1], 0.0))
    for off in FAULT_NORMAL_OFFSETS:
        q = p0 + off * n0
        side = "NE" if off > 0 else "SW"
        receivers.append((f"fault_normal_{side}_{abs(off)/1000:.0f}km", q[0], q[1], 0.0))
    marched = march_trace(trace, p0, head_nw, ALONG_STRIKE_DIST)
    for (xy, head, s) in marched:
        nrm = np.array([-head[1], head[0]])      # local normal at the marched point
        q = xy + ALONG_STRIKE_OFFSET * nrm
        receivers.append((f"along_strike_NW_{s/1000:.0f}km", q[0], q[1], 0.0))
    if len(marched) < len(ALONG_STRIKE_DIST):
        print(f"WARNING: trace march reached only {len(marched)}/{len(ALONG_STRIKE_DIST)} "
              "along-strike targets (strand may be shorter than 40 km NW).")

    # ---- on-fault pickpoints (every chosen point IS a fault facet centroid -> on-fault) --------
    # (label, facet_index): hypocenter via 3D nearest; the rest via map-proximity + closest depth.
    pp_picks = []
    _, iH = tree.query(HYPOCENTER)
    pp_picks.append(("hypocenter", iH))
    for d in DOWNDIP_DEPTHS:                                  # downdip column near the epicenter
        idx = pick_at_depth(centroids, EPICENTER, d, DOWNDIP_MAP_R)
        if idx is None:
            print(f"WARNING: no fault facet within {DOWNDIP_MAP_R/1000:.0f} km of the epicenter "
                  f"at ~{d/1000:.0f} km depth (downdip target skipped).")
        else:
            pp_picks.append((f"downdip_{d/1000:.0f}km", idx))
    march_pp = {0.0: p0}                                      # along-strike at seismogenic depth
    for (xy, head, s) in march_trace(trace, p0, head_nw, [d for d in ALONG_STRIKE_PP_DIST if d > 0]):
        march_pp[s] = xy
    for s in ALONG_STRIKE_PP_DIST:
        if s not in march_pp:
            continue
        idx = pick_at_depth(centroids, march_pp[s], ALONG_STRIKE_PP_DEPTH, ALONG_STRIKE_MAP_R)
        if idx is None:
            print(f"WARNING: no fault facet within {ALONG_STRIKE_MAP_R/1000:.0f} km of the "
                  f"s={s/1000:.0f} km trace point (along-strike pickpoint skipped).")
        else:
            pp_picks.append((f"strike_NW_{s/1000:.0f}km", idx))

    pickpoints = []   # (label, x, y, z)  deduped by facet
    seen = set()
    for label, idx in pp_picks:
        c = centroids[idx]
        key = int(idx)
        if key in seen:
            print(f"  (skip duplicate pickpoint '{label}' -> same facet as a prior target)")
            continue
        seen.add(key)
        pickpoints.append((f"{label}(z={-c[2]/1000:.1f}km)", c[0], c[1], c[2]))

    # ---- write files ----------------------------------------------------------
    def write(coords_path, ann_path, rows, kind):
        with open(coords_path, "w") as fc:
            for (_, x, y, z) in rows:
                fc.write(f"{x:.4f} {y:.4f} {z:.4f}\n")
        with open(ann_path, "w") as fa:
            fa.write(f"# {kind} (annotated; SeisSol reads the coords-only {coords_path})\n")
            fa.write("# label                                  x            y            z\n")
            for (label, x, y, z) in rows:
                fa.write(f"# {label:38s} {x:12.4f} {y:12.4f} {z:12.4f}\n")
        print(f"wrote {len(rows):2d} {kind:10s} -> {coords_path}  (+ {ann_path})")

    write("safs_receivers.dat", "safs_receivers.annotated.txt", receivers, "receivers")
    write("safs_pickpoints.dat", "safs_pickpoints.annotated.txt", pickpoints, "pickpoints")


if __name__ == "__main__":
    main()
