#!/usr/bin/env python3
"""make_surgery_boxes.py -- turn residual sliver positions into disjoint refill boxes.

Reads the localizer's npz (P = sliver barycenters), clusters by grid linkage,
inflates each cluster's bbox by a margin, enforces a minimum half-extent (a box
must give tetgen room to build good interior tets), and merges overlapping boxes
until all are disjoint -- overlapping surgeries on one global mesh would double-
drop tets.
"""
import argparse
import json
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--locs', default='build_tmp/sliver_locs_optimI.npz')
ap.add_argument('--out', default='build_tmp/surgery_boxes.json')
ap.add_argument('--link', type=float, default=1000.0, help='cluster linkage (m)')
ap.add_argument('--margin', type=float, default=400.0)
ap.add_argument('--min-half', type=float, default=650.0)
a = ap.parse_args()

P = np.load(a.locs)['P']
print(f'[in] {len(P):,} sliver positions')

# grid-linkage clustering
key = np.floor(P / a.link).astype(np.int64)
u, inv = np.unique(key, axis=0, return_inverse=True)
# union adjacent occupied cells (26-neighbourhood)
parent = np.arange(len(u))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i
cell = {tuple(c): i for i, c in enumerate(u)}
for i, c in enumerate(u):
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                j = cell.get((c[0] + dx, c[1] + dy, c[2] + dz))
                if j is not None:
                    pi, pj = find(i), find(j)
                    if pi != pj:
                        parent[pj] = pi
lab = np.array([find(i) for i in inv])
boxes = []
for L in np.unique(lab):
    Q = P[lab == L]
    lo = Q.min(0) - a.margin
    hi = Q.max(0) + a.margin
    c = 0.5 * (lo + hi)
    h = np.maximum(0.5 * (hi - lo), a.min_half)
    boxes.append([float(c[0]), float(c[1]), float(h[0]), float(h[1]),
                  float(c[2] - h[2]), float(c[2] + h[2]), int(len(Q))])

# merge overlapping boxes until disjoint
def overlap(b1, b2):
    return (abs(b1[0] - b2[0]) < b1[2] + b2[2] and
            abs(b1[1] - b2[1]) < b1[3] + b2[3] and
            b1[4] < b2[5] and b2[4] < b1[5])
merged = True
while merged:
    merged = False
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if overlap(boxes[i], boxes[j]):
                b1, b2 = boxes[i], boxes[j]
                x0 = min(b1[0] - b1[2], b2[0] - b2[2]); x1 = max(b1[0] + b1[2], b2[0] + b2[2])
                y0 = min(b1[1] - b1[3], b2[1] - b2[3]); y1 = max(b1[1] + b1[3], b2[1] + b2[3])
                z0 = min(b1[4], b2[4]); z1 = max(b1[5], b2[5])
                boxes[i] = [0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (x1 - x0),
                            0.5 * (y1 - y0), z0, z1, b1[6] + b2[6]]
                del boxes[j]
                merged = True
                break
        if merged:
            break

boxes.sort(key=lambda b: -b[6])
for k, b in enumerate(boxes):
    print(f'  box {k:02d}: {b[6]:>4} slivers   c=({b[0]:,.0f},{b[1]:,.0f})  '
          f'h=({b[2]:,.0f},{b[3]:,.0f})  z {b[4]:,.0f}..{b[5]:,.0f}')
json.dump(boxes, open(a.out, 'w'))
print(f'[out] {a.out}  ({len(boxes)} disjoint boxes, {sum(b[6] for b in boxes)} slivers)')
