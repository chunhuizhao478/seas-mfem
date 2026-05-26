#!/usr/bin/env python3
"""Mark the Vmax-blowup points and fault elements (job 7748818, LSW SAFS Dc2).

Reads the DIAG-ONSET argmax records from the run log, aggregates them per
physical fault location (peak |V| and number of steps it was the global
argmax), maps each location onto the fault-surface triangulation, and writes:

  * <out>/safs_fault_blowup_marked.vtu  -- fault surface with cell/point arrays
        blowup_peakV   : peak |V| (m/s) of the nearest logged hotspot (0 elsewhere)
        blowup_loghits : how many steps that hotspot was the global argmax
        is_blowup      : 1 if peakV > VBLOW (clearly unphysical), else 0
  * <out>/safs_fault_blowup_points.vtu  -- the hotspots as VTK vertices (Glyph me)
  * <out>/safs_fault_blowup_points.csv  -- x,y,z,peakV,loghits,is_blowup

Color the fault by is_blowup / log10(blowup_peakV); overlay the points file as
spheres (Glyph).  Physical rupture peaked at ~10-15 m/s, so VBLOW=50 m/s flags
the runaway DOFs only.
"""
import csv
import re
import sys
import numpy as np
import meshio

LOG   = "/Users/chunhuizhao/Documents/spatial_dyn_resDc2_7748818.log"
FAULT = ("/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/"
         "project_7.0_alternative/meshing/results/vtu/"
         "safs_fault_box_nwcut_500m_lcfar3000_z0embed_fault.vtu")
OUTDIR = ("/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/safs/"
          "project_7.0_alternative/meshing/results/vtu")

VBLOW    = 50.0     # m/s; |V| above this is a runaway (physical peak ~10-15)
MARK_R   = 500.0    # m; light up fault cells within this radius of a hotspot

# ---- 1. aggregate DIAG-ONSET hotspots from the log -----------------------
# [DIAG-ONSET] ... xyz=(X,Y,Z) V=<vabs> V1=.. V2=.. ...
pat = re.compile(r"xyz=\(([^)]*)\)\s+V=([0-9.eE+\-]+)")
peak, hits = {}, {}
with open(LOG) as fh:
    for line in fh:
        if not line.startswith("[DIAG-ONSET]"):
            continue
        m = pat.search(line)
        if not m:
            continue
        xyz = tuple(float(v) for v in m.group(1).split(","))
        v = float(m.group(2))
        hits[xyz] = hits.get(xyz, 0) + 1
        if v > peak.get(xyz, -1.0):
            peak[xyz] = v

if not peak:
    sys.exit("no [DIAG-ONSET] records found in " + LOG)

pts   = np.array(list(peak.keys()), dtype=float)          # (N,3)
pV    = np.array([peak[tuple(p)] for p in pts])
pH    = np.array([hits[tuple(p)] for p in pts], dtype=float)
isbl  = (pV > VBLOW).astype(np.int32)
print(f"hotspots: {len(pts)} total, {int(isbl.sum())} with peakV>{VBLOW:g} m/s")
print("worst by peak |V|:")
for i in np.argsort(-pV)[:8]:
    print(f"  ({pts[i,0]:.0f},{pts[i,1]:.0f},{pts[i,2]:.1f})  peakV={pV[i]:.3g}  hits={int(pH[i])}")

# ---- 2. read the fault-surface mesh --------------------------------------
fm = meshio.read(FAULT)
P  = fm.points
tri = None
for cb in fm.cells:
    if cb.type == "triangle":
        tri = cb.data
        break
if tri is None:
    sys.exit("no triangle cells in " + FAULT)
cent = P[tri].mean(axis=1)                                # (Ncell,3) centroids
print(f"fault mesh: {P.shape[0]} points, {tri.shape[0]} triangles")

# ---- 3. map hotspots onto cells (nearest + within MARK_R) ----------------
cell_peakV = np.zeros(tri.shape[0])
cell_hits  = np.zeros(tri.shape[0])
pt_peakV   = np.zeros(P.shape[0])
# only mark with blowup hotspots (runaways); change to `range(len(pts))` for all
sel = np.where(isbl == 1)[0]
for i in sel:
    d = np.linalg.norm(cent - pts[i], axis=1)
    near = np.where(d <= MARK_R)[0]
    if near.size == 0:
        near = np.array([int(d.argmin())])               # always mark nearest
    cell_peakV[near] = np.maximum(cell_peakV[near], pV[i])
    cell_hits[near]  = np.maximum(cell_hits[near],  pH[i])
    for c in near:
        pt_peakV[tri[c]] = np.maximum(pt_peakV[tri[c]], pV[i])

cell_isbl = (cell_peakV > 0).astype(np.int32)
print(f"marked fault triangles: {int(cell_isbl.sum())}")

# ---- 4. write the marked fault surface -----------------------------------
out_mesh = meshio.Mesh(
    points=P,
    cells=[("triangle", tri)],
    cell_data={"blowup_peakV":  [cell_peakV],
               "blowup_loghits": [cell_hits],
               "is_blowup":      [cell_isbl]},
    point_data={"blowup_peakV": pt_peakV},
)
p_marked = OUTDIR + "/safs_fault_blowup_marked.vtu"
meshio.write(p_marked, out_mesh)
print("wrote", p_marked)

# ---- 5. write the hotspot points (all of them, with is_blowup flag) ------
vcells = [("vertex", np.arange(len(pts)).reshape(-1, 1))]
pts_mesh = meshio.Mesh(
    points=pts,
    cells=vcells,
    point_data={"peakV": pV, "loghits": pH, "is_blowup": isbl},
)
p_pts = OUTDIR + "/safs_fault_blowup_points.vtu"
meshio.write(p_pts, pts_mesh)
print("wrote", p_pts)

p_csv = OUTDIR + "/safs_fault_blowup_points.csv"
with open(p_csv, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["x", "y", "z", "peakV_mps", "log_argmax_hits", "is_blowup"])
    for i in np.argsort(-pV):
        w.writerow([f"{pts[i,0]:.3f}", f"{pts[i,1]:.3f}", f"{pts[i,2]:.3f}",
                    f"{pV[i]:.6g}", int(pH[i]), int(isbl[i])])
print("wrote", p_csv)
