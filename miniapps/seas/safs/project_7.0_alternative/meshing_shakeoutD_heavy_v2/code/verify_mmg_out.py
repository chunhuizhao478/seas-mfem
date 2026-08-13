#!/usr/bin/env python3
"""verify_mmg_out.py -- check the mmg MEDIT output INDEPENDENTLY of the driver.

The driver refused to write because the fault facet count fell by one. That is a
strict-equality guard, not a geometry check (it also reported fault vertices moved
2e-11 m), so the mesh is judged here on its own terms: volume, quality, dt, the
free-surface plane, and exactly WHICH fault facet went missing.

Reads MEDIT by byte offset + pandas rather than meshio: the file is 3.06 GB ASCII
and 70.6 M tets, so a generic reader is both too slow and too heavy.
"""
import numpy as np, pandas as pd, sys
F = sys.argv[1] if len(sys.argv) > 1 else "build_tmp/s1.mmg_out.mesh"
# MEASURED, not hardcoded. These were PREFERRED's byte layout; on any other mesh
# they land mid-vertex-block and pandas parses whatever is there WITHOUT error.
from medit_hdr import medit_sections
_sec = medit_sections(F)
LV, nv = _sec["Vertices"]
LT, nt = _sec["Tetrahedra"]
LS, ns = _sec["Triangles"]
print(f"[hdr] verts {nv:,} @ {LV:,}   tets {nt:,} @ {LT:,}   tris {ns:,} @ {LS:,}")
def block(start_hdr, n, ncol, dt):
    return pd.read_csv(F, sep=r"\s+", header=None, skiprows=start_hdr+1,
                       nrows=n, usecols=range(ncol), dtype=dt,
                       engine="c").to_numpy()
P = block(LV, nv, 3, np.float64)
print(f"[verts] {len(P):,}   z {P[:,2].min():.4f} .. {P[:,2].max():.4f} m")
S = block(LS, ns, 4, np.int32)
tri, ref = S[:, :3] - 1, S[:, 3]
u, c = np.unique(ref, return_counts=True)
print("[surface] " + "  ".join(f"{a}:{b:,}" for a, b in zip(u, c)))
top = tri[ref == 102]
zt = P[np.unique(top), 2]
print(f"[free surface] {len(top):,} tris, {len(np.unique(top)):,} nodes, "
      f"z {zt.min():.6f} .. {zt.max():.6f} m")

# --- which fault facet was lost? match centroids against the PLC fault ---
from scipy.spatial import cKDTree
d = np.load("build_tmp/fill.npz")
plc = d["plcT"][d["MARK"] == 7]; Pb = d["P"]
Cb = Pb[plc].mean(1)
Co = P[tri[ref == 101]].mean(1)
print(f"[fault] PLC {len(Cb):,}  mmg-out {len(Co):,}  delta {len(Cb)-len(Co)}")
dd, ii = cKDTree(Co).query(Cb, k=1, workers=-1)
lost = np.flatnonzero(dd > 1.0)          # 1 m: facets are 21-124 m, vertices moved 2e-11
print(f"[fault] centroid match: max residual on matched {np.sort(dd)[-len(lost)-1]:.3e} m")
Ab = 0.5*np.linalg.norm(np.cross(Pb[plc[:,1]]-Pb[plc[:,0]], Pb[plc[:,2]]-Pb[plc[:,0]]), axis=1)
print(f"[fault] {len(lost)} PLC facets with NO match in the output:")
for i in lost:
    t = Pb[plc[i]]
    print(f"   area {Ab[i]:9.1f} m2   centroid {t.mean(0)[0]:.1f} {t.mean(0)[1]:.1f} "
          f"{t.mean(0)[2]:8.2f}   z {t[:,2].min():.2f}..{t[:,2].max():.2f}")
# DERIVED. 7.96255e9 was PREFERRED's fault area; ALT's is 6.5377e9, so the
# hardcoded value would misreport the percentage by 22 % on this tree.
FAULT_AREA = float(Ab.sum())
print(f"[fault] total area {FAULT_AREA/1e6:,.2f} km2 (derived)")
print(f"[fault] lost area {Ab[lost].sum():.1f} m2 = {100*Ab[lost].sum()/FAULT_AREA:.3e} % of fault")
del Cb, Co, plc, Pb, d, S, tri, ref, top
