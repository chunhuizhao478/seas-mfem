#!/usr/bin/env python3
"""gate_boxes.py -- surgery boxes for GATE failures, each with an edge target.

Locates every cell below the 1 Hz p5 gate (Vs/dx < 0.8 on native MUSCAL),
clusters them like make_surgery_boxes, and annotates each box with the edge
length its refill must not exceed:

    edge_target = min(Vs among the box's failures) / gate * safety

The refill turns that into a tetgen volume cap (-a), so the fresh interior
CANNOT recreate the failure -- unlike bisection, which fights the CVM bin
treadmill one split at a time (measured: accept-drive drained ~45 cells/round
against 5,564).
"""
import argparse
import json
import sys
import numpy as np, h5py
sys.path.insert(0, 'code')
from material import Material

BOX = (132000.0, 724000.0, 3490000.0, 4055000.0)
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
CH = 6_000_000

ap = argparse.ArgumentParser()
ap.add_argument('--mesh', required=True, help='PUML')
ap.add_argument('--cvm', required=True)
ap.add_argument('--muscal', required=True)
ap.add_argument('--gate', type=float, default=0.8)
ap.add_argument('--safety', type=float, default=0.85)
ap.add_argument('--link', type=float, default=1500.0)
ap.add_argument('--margin', type=float, default=500.0)
ap.add_argument('--min-half', type=float, default=700.0)
ap.add_argument('--tile', type=float, default=0.0,
                help='fixed disjoint TILE mode (m). The linkage clustering chains '
                     'scattered basin failures into 60x130 km sheet boxes, which a '
                     'blanket -a refill would over-refine absurdly. Tiles partition '
                     'by barycentre, are disjoint by construction, and stay small.')
ap.add_argument('--out', default='build_tmp/gate_boxes.json')
a = ap.parse_args()

mat = Material(a.cvm, a.muscal, box=BOX, source='muscal', verbose=False)
f = h5py.File(a.mesh)
V = f['geometry'][:]
C = f['connect']
NT = C.shape[0]
P, VS = [], []
for s in range(0, NT, CH):
    t = C[s:s + CH][:].astype(np.int64)
    p = V[t]
    dx = np.max(np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1), 1)
    vb = mat.at(p.mean(1))
    m = np.where(vb > 0, vb / dx, 0.0) < a.gate
    if m.any():
        P.append(p[m].mean(1))
        VS.append(vb[m])
    del t, p, dx, vb
P = np.vstack(P)
VS = np.concatenate(VS)
print(f'[gate] {len(P):,} failing cells   Vs {VS.min():.0f}..{VS.max():.0f} m/s')

if a.tile > 0:
    T = a.tile
    key = np.floor(P[:, :2] / T).astype(np.int64)
    u, inv = np.unique(key, axis=0, return_inverse=True)
    boxes = []
    for i in range(len(u)):
        m = inv == i
        z0 = float(P[m, 2].min() - 400.0)
        et = float(VS[m].min()) / a.gate * a.safety
        cx, cy = (u[i] + 0.5) * T
        boxes.append([float(cx), float(cy), T / 2, T / 2, z0, 0.0, int(m.sum()), round(et, 1)])
    boxes.sort(key=lambda b: -b[6])
    for k, b in enumerate(boxes[:8]):
        print(f'  tile {k:02d}: {b[6]:>5} fails  c=({b[0]:,.0f},{b[1]:,.0f})  '
              f'h={b[2]:,.0f}  z {b[4]:,.0f}..0  edge<= {b[7]:,.0f} m')
    import json as _j
    _j.dump(boxes, open(a.out, 'w'))
    print(f'[out] {a.out}  ({len(boxes)} disjoint tiles, {sum(b[6] for b in boxes)} fails)')
    raise SystemExit(0)

# cluster (same grid-linkage as make_surgery_boxes)
key = np.floor(P / a.link).astype(np.int64)
u, inv = np.unique(key, axis=0, return_inverse=True)
parent = np.arange(len(u))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i
cell = {tuple(c): i for i, c in enumerate(u)}
for i, c in enumerate(u):
    for dx_ in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                j = cell.get((c[0] + dx_, c[1] + dy, c[2] + dz))
                if j is not None:
                    pi, pj = find(i), find(j)
                    if pi != pj:
                        parent[pj] = pi
lab = np.array([find(i) for i in inv])
boxes = []
for L in np.unique(lab):
    m = lab == L
    Q = P[m]
    lo = Q.min(0) - a.margin
    hi = Q.max(0) + a.margin
    c = 0.5 * (lo + hi)
    h = np.maximum(0.5 * (hi - lo), a.min_half)
    et = float(VS[m].min()) / a.gate * a.safety
    boxes.append([float(c[0]), float(c[1]), float(h[0]), float(h[1]),
                  float(c[2] - h[2]), float(min(c[2] + h[2], 0.0)), int(m.sum()), round(et, 1)])

def overlap(b1, b2):
    return (abs(b1[0] - b2[0]) < b1[2] + b2[2] and abs(b1[1] - b2[1]) < b1[3] + b2[3]
            and b1[4] < b2[5] and b2[4] < b1[5])
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
                            0.5 * (y1 - y0), z0, z1, b1[6] + b2[6], min(b1[7], b2[7])]
                del boxes[j]
                merged = True
                break
        if merged:
            break
boxes.sort(key=lambda b: -b[6])
for k, b in enumerate(boxes[:12]):
    print(f'  box {k:02d}: {b[6]:>5} fails  c=({b[0]:,.0f},{b[1]:,.0f})  '
          f'h=({b[2]:,.0f},{b[3]:,.0f})  z {b[4]:,.0f}..{b[5]:,.0f}  edge<= {b[7]:,.0f} m')
json.dump(boxes, open(a.out, 'w'))
print(f'[out] {a.out}  ({len(boxes)} boxes, {sum(b[6] for b in boxes)} fails)')
