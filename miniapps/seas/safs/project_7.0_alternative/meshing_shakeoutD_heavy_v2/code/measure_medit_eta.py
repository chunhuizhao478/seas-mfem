#!/usr/bin/env python3
"""measure_medit_eta.py -- eta + skin-integrity check on a (patch-scale) MEDIT mesh.

Second arg: the extractor's meta npz. Every skin vertex must still exist in the
output EXACTLY (max displacement 0) or the patch cannot be stitched back, and the
fault triangle count must be unchanged.
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, 'code')
from medit_hdr import medit_sections
from scipy.spatial import cKDTree

F = sys.argv[1]
META = sys.argv[2] if len(sys.argv) > 2 else None
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

sec = medit_sections(F)
LV, NV = sec['Vertices']
LT, NT = sec['Tetrahedra']
NS = sec['Triangles'][1] if 'Triangles' in sec else 0
P = pd.read_csv(F, sep=r'\s+', header=None, skiprows=LV + 1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine='c').to_numpy()
T = pd.read_csv(F, sep=r'\s+', header=None, skiprows=LT + 1, nrows=NT,
                usecols=range(4), dtype=np.int64, engine='c').to_numpy() - 1
p = P[T]
e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
d6 = np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                      p[:, 3] - p[:, 0])) / 6.0
ss = (e ** 2).sum(1)
eta = np.where(ss > 0, 12.0 * np.cbrt((3.0 * d6) ** 2) / ss, 0.0)
print(f'{F}')
print(f'  {NT:,} tets  {NS:,} tris   eta_min {eta.min():.4f}   '
      f'<0.05 {int((eta < 0.05).sum()):,}   <0.1 {int((eta < 0.1).sum()):,}   '
      f'med {np.median(eta):.3f}   edge_min {e.min():.2f} m')
if META:
    z = np.load(META)
    S = z['skinP']
    d, _ = cKDTree(P).query(S, k=1)
    print(f'  skin {len(S):,} verts: max displacement {d.max():.3e} m   '
          f'(fault tris expected {int(z["nfault"]):,})')
