#!/usr/bin/env python3
"""sliver_localize.py -- where exactly are v1's eta<0.05 tets, and how big is a
repair patch that contains them?

The repair-vs-rebuild decision hangs on locality: if the sliver population lives
in a thin fault-hugging band, a LOCAL repair (patch extraction with a frozen
skin, or a tiled pass) is feasible without touching the other ~86M tets.
Reports the distance-to-fault distribution of the slivers and the densest
hotspots, and dumps their barycenters for the patch extractor.
"""
import numpy as np, h5py
from scipy.spatial import cKDTree

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
CH = 4_000_000
MESH = '../meshing_shakeoutD_heavy/results/safalt_shakeoutD_heavy_1Hz_p5.puml.h5'

f = h5py.File(MESH)
V = f['geometry'][:]
C = f['connect']
B = f['boundary']
NT = C.shape[0]

ft = []
for s in range(0, NT, CH):
    t = C[s:s + CH][:].astype(np.int64)
    b = B[s:s + CH][:].astype(np.int64)
    for k in range(4):
        m = ((b >> (8 * k)) & 0xFF) == 3
        if m.any():
            ft.append(t[m][:, list(FACE[k])])
    del t, b
FT = np.vstack(ft)
Fc = V[FT].mean(1)[::10]
del ft, FT
tree = cKDTree(Fc)
print(f'[fault] {len(Fc):,} decimated centroids', flush=True)

bar, et = [], []
for s in range(0, NT, CH):
    t = C[s:s + CH][:].astype(np.int64)
    p = V[t]
    e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
    d6 = np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                          p[:, 3] - p[:, 0])) / 6.0
    ss = (e ** 2).sum(1)
    eta = np.where(ss > 0, 12.0 * np.cbrt((3.0 * d6) ** 2) / ss, 0.0)
    m = eta < 0.05
    if m.any():
        bar.append(p[m].mean(1))
        et.append(eta[m])
    del t, p, e, d6, ss, eta
Pb = np.vstack(bar)
E = np.concatenate(et)
d, _ = tree.query(Pb, k=1, distance_upper_bound=5000.0, workers=-1)
d = np.where(np.isfinite(d), d, 5000.0)

print(f'\n[slivers] eta<0.05: {len(Pb):,}')
for q in (50, 90, 99, 99.9, 100):
    print(f'  distance-to-fault pct {q:>5}: {np.percentile(d, q):>8.1f} m')
print(f'  z range {Pb[:,2].min():,.0f} .. {Pb[:,2].max():,.0f}   '
      f'z median {np.median(Pb[:,2]):,.0f}')

key = np.floor(Pb / 4000.0).astype(np.int64)
u, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
o = np.argsort(cnt)[::-1]
print('\n[top 8 hotspot cells, 4 km grid]')
for j in o[:8]:
    c = Pb[inv == j].mean(0)
    print(f'  {cnt[j]:>8,} slivers   E {c[0]:>9,.0f}  N {c[1]:>11,.0f}  z {c[2]:>9,.0f}')
np.savez('build_tmp/sliver_locs.npz', P=Pb, eta=E, d=d)
print('\n[out] build_tmp/sliver_locs.npz')
