#!/usr/bin/env python3
"""sliver_localize_medit.py -- where are the eta<0.05 tets in a MEDIT mesh?

Same job as sliver_localize.py but for the mid-pipeline MEDIT state: distance to
the fault (ref 101 triangles), z structure, hotspots, and how many survivors sit
on the fault vs free interior -- the inputs for deciding between another global
pass and targeted patch repair.
"""
import argparse
import sys
import numpy as np, pandas as pd
sys.path.insert(0, 'code')
from medit_hdr import medit_sections
from scipy.spatial import cKDTree

_ap = argparse.ArgumentParser()
_ap.add_argument('mesh')
_ap.add_argument('--thr', type=float, default=0.05,
                 help='eta threshold defining a "sliver" (round 2 of surgery uses 0.1)')
_a = _ap.parse_args()
F = _a.mesh
THR = _a.thr
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
CH = 8_000_000

sec = medit_sections(F)
LV, NV = sec['Vertices']
LT, NT = sec['Tetrahedra']
LS, NS = sec['Triangles']
P = pd.read_csv(F, sep=r'\s+', header=None, skiprows=LV + 1, nrows=NV,
                usecols=range(3), dtype=np.float64, engine='c').to_numpy()
S = pd.read_csv(F, sep=r'\s+', header=None, skiprows=LS + 1, nrows=NS,
                usecols=range(4), dtype=np.int64, engine='c').to_numpy()
ftri = S[S[:, 3] == 101][:, :3] - 1
Fc = P[ftri].mean(1)[::10]
tree = cKDTree(Fc)
fverts = np.zeros(NV, bool)
fverts[np.unique(ftri)] = True
print(f'[fault] {len(ftri):,} tris -> {len(Fc):,} KD centroids', flush=True)

bar, et, nfv = [], [], []
for s in range(0, NT, CH):
    n = min(CH, NT - s)
    T = pd.read_csv(F, sep=r'\s+', header=None, skiprows=LT + 1 + s, nrows=n,
                    usecols=range(4), dtype=np.int32, engine='c').to_numpy() - 1
    p = P[T]
    e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
    d6 = np.abs(np.einsum('ij,ij->i', np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                          p[:, 3] - p[:, 0])) / 6.0
    ss = (e ** 2).sum(1)
    eta = np.where(ss > 0, 12.0 * np.cbrt((3.0 * d6) ** 2) / ss, 0.0)
    m = eta < THR
    if m.any():
        bar.append(p[m].mean(1))
        et.append(eta[m])
        nfv.append(fverts[T[m]].sum(1))
    del T, p, e, d6, ss, eta
    print(f'  ..{s+n:,}/{NT:,}', flush=True)
Pb = np.vstack(bar)
E = np.concatenate(et)
NF = np.concatenate(nfv)
d, _ = tree.query(Pb, k=1, distance_upper_bound=5000.0, workers=-1)
d = np.where(np.isfinite(d), d, 5000.0)

print(f'\n[slivers] eta<{THR:g}: {len(Pb):,}')
for q in (50, 90, 99, 100):
    print(f'  distance-to-fault pct {q:>3}: {np.percentile(d, q):>8.1f} m')
print(f'  z median {np.median(Pb[:,2]):,.0f}   range {Pb[:,2].min():,.0f} .. {Pb[:,2].max():,.0f}')
print(f'  fault-vertex count on sliver: 0v {int((NF==0).sum()):,}  1-2v {int(((NF>0)&(NF<3)).sum()):,}  >=3v {int((NF>=3).sum()):,}')
key = np.floor(Pb / 4000.0).astype(np.int64)
u, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
o = np.argsort(cnt)[::-1]
print('\n[top 6 hotspots, 4 km grid]  (cumulative share)')
cum = 0
for j in o[:6]:
    c = Pb[inv == j].mean(0)
    cum += cnt[j]
    print(f'  {cnt[j]:>7,}  E {c[0]:>9,.0f}  N {c[1]:>11,.0f}  z {c[2]:>9,.0f}   ({100*cum/len(Pb):.1f} %)')
print(f'\n[coverage] cells needed for 50/90/99 % of slivers: '
      f'{np.searchsorted(np.cumsum(np.sort(cnt)[::-1])/len(Pb), [0.5,0.9,0.99])+1}'
      f'  of {len(u):,} occupied 4 km cells')
np.savez('build_tmp/sliver_locs_optimI.npz', P=Pb, eta=E, d=d, nf=NF)
