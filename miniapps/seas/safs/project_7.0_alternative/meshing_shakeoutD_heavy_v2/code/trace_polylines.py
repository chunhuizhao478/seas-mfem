#!/usr/bin/env python3
"""trace_polylines.py -- ordered polylines of the fault's TOP boundary edge.

The PLC's top surface has to carry the trace as a hard constraint, node for node,
or the fault will not daylight.  gmsh needs that as ORDERED curves, not a soup of
edges, so the top boundary of the fault surface is walked into chains here.
Multi-strand and branching are expected: chains simply end at a valence != 2 node.
"""
import numpy as np
from collections import defaultdict
d = np.load("build_tmp/fault_surface.npz")
P, T, bnd = d["P"], d["T"], d["bnd"]
ZTOP = -60.0   # MUST equal extract_fault_surface.py's ZTOP
m = (P[bnd[:, 0], 2] > ZTOP) & (P[bnd[:, 1], 2] > ZTOP)
E = bnd[m]
print(f"[trace] {len(E):,} top boundary edges, {len(np.unique(E)):,} vertices")
adj = defaultdict(list)
for a, b in E:
    adj[int(a)].append(int(b)); adj[int(b)].append(int(a))
val = {v: len(n) for v, n in adj.items()}
ends = [v for v, k in val.items() if k != 2]
print(f"  nodes with valence != 2 (chain ends / junctions): {len(ends)}")
used = set(); chains = []
def walk(s, n):
    ch = [s, n]; used.add((min(s, n), max(s, n)))
    while val.get(ch[-1], 0) == 2:
        a, b = adj[ch[-1]]
        nxt = a if a != ch[-2] else b
        k = (min(ch[-1], nxt), max(ch[-1], nxt))
        if k in used: break
        used.add(k); ch.append(nxt)
    return ch
for s in ends:
    for n in adj[s]:
        if (min(s, n), max(s, n)) not in used:
            chains.append(walk(s, n))
for a, b in E:                                  # any closed loops left
    k = (min(int(a), int(b)), max(int(a), int(b)))
    if k not in used:
        chains.append(walk(int(a), int(b)))
tot = sum(len(c) - 1 for c in chains)
L = sum(float(np.linalg.norm(np.diff(P[c], axis=0), axis=1).sum()) for c in chains)
print(f"[chains] {len(chains)} polylines, {tot:,} segments (of {len(E):,}), total {L/1000:,.1f} km")
print(f"  lengths: " + ", ".join(f"{float(np.linalg.norm(np.diff(P[c],axis=0),axis=1).sum())/1000:.1f}" 
                                  for c in sorted(chains, key=len, reverse=True)[:8]) + " km ...")
seg = np.array([len(c) for c in chains])
print(f"  segments per chain: min {seg.min()} med {int(np.median(seg))} max {seg.max()}")
# Fault edges that lie IN the flat-top plane but are NOT on the fault boundary.
# They join two daylighting nodes across a near-horizontal fault patch, so they sit
# exactly in the free surface and MUST be embedded too -- otherwise they cut across
# top-surface triangles and tetgen reports self-intersections (measured: 4 such
# edges, 69.9 m each, and they were the entire cause of the 11 skipped facets).
z0 = np.abs(P[:, 2]) < 1e-9
EA = np.unique(np.sort(np.concatenate([T[:, [0,1]], T[:, [1,2]], T[:, [0,2]]]), axis=1), axis=0)
ipl = EA[z0[EA[:, 0]] & z0[EA[:, 1]]]
bs = set(map(tuple, np.sort(bnd, axis=1)))
extra = np.array([e for e in ipl if tuple(e) not in bs], np.int64).reshape(-1, 2)
print(f"[in-plane] {len(extra)} interior fault edges lie in the z=0 plane -> embed them too")
np.save("build_tmp/inplane_edges.npy", extra)
np.savez_compressed("build_tmp/trace_chains.npz",
                    **{f"c{i}": np.array(c, np.int64) for i, c in enumerate(chains)})
print("[out] build_tmp/trace_chains.npz")
