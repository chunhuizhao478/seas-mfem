#!/usr/bin/env python3
"""gate_vs_eta.py -- cross-tabulate the frequency gate against element quality.

The mesh-build skill is explicit that this must be run BEFORE choosing a repair
tool, because gate failures are two populations needing OPPOSITE fixes:

  eta < 0.05  slivers        bisecting one yields two slivers -- LEB is WRONG here
  eta > 0.7   merely coarse  well-shaped and just too big     -- LEB is RIGHT here

Reported previously on a sibling mesh as a 35.6 % vs 0.61 % fail rate. Without
this split, one tool is applied to both populations and the residual never clears.

Also reports where each population sits relative to the frozen fault, because a
tet with a FACE on the fault has 3 of 4 vertices immovable and no smoothing pass
can reshape it -- that decides whether repair is even available.

Everything chunked; --muscal is REQUIRED (Material silently falls back to the
deck cube without it, which overstates gate failures ~27x).
"""
import argparse
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "code")
sys.path.insert(0, "../meshing_shakeoutD_heavy/code")
from material import Material

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
FACE = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
BOX = (132000.0, 724000.0, 3490000.0, 4055000.0)
CH = 4_000_000

ap = argparse.ArgumentParser()
ap.add_argument("--mesh", required=True)
ap.add_argument("--cvm", required=True)
ap.add_argument("--muscal", required=True)
ap.add_argument("--gate", type=float, default=0.8)
a = ap.parse_args()

mat = Material(a.cvm, a.muscal, box=BOX, source="muscal", verbose=False)
f = h5py.File(a.mesh)
V = f["geometry"][:]
C = f["connect"]
B = f["boundary"]
NT = C.shape[0]

# fault centroids for the distance field, from the mesh itself
ftri = []
for s in range(0, NT, CH):
    t = C[s:s + CH][:].astype(np.int64)
    b = B[s:s + CH][:].astype(np.int64)
    for k in range(4):
        m = ((b >> (8 * k)) & 0xFF) == 3
        if m.any():
            ftri.append(t[m][:, list(FACE[k])])
    del t, b
FT = np.vstack(ftri)
Fc = V[FT].mean(1)[::10]
del ftri, FT
tree = cKDTree(Fc)
print(f"[fault] {len(Fc):,} decimated centroids for the distance field", flush=True)

EB = [0.0, 0.05, 0.1, 0.3, 0.7, 1.01]                 # eta bands
DB = [0.0, 250.0, 1000.0, 5000.0, 1e9]                # distance-to-fault bands
ncell = np.zeros((len(EB) - 1, len(DB) - 1), np.int64)
nfail = np.zeros_like(ncell)
onfault = np.zeros(len(EB) - 1, np.int64)             # tets with a FACE on the fault

for s in range(0, NT, CH):
    t = C[s:s + CH][:].astype(np.int64)
    b = B[s:s + CH][:].astype(np.int64)
    p = V[t]
    e = np.stack([np.linalg.norm(p[:, j] - p[:, i], axis=1) for i, j in PAIRS], 1)
    dx = e.max(1)
    d6 = np.abs(np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                          p[:, 3] - p[:, 0])) / 6.0
    ssum = (e ** 2).sum(1)
    eta = np.where(ssum > 0, 12.0 * np.cbrt((3.0 * d6) ** 2) / ssum, 0.0)
    cen = p.mean(1)
    vb = mat.at(cen)
    fail = np.where(vb > 0, vb / dx, 0.0) < a.gate
    dist, _ = tree.query(cen, k=1, distance_upper_bound=6000.0, workers=-1)
    dist = np.where(np.isfinite(dist), dist, 1e9)
    hasface = np.zeros(len(t), bool)
    for k in range(4):
        hasface |= ((b >> (8 * k)) & 0xFF) == 3
    ei = np.clip(np.searchsorted(EB, eta, side="right") - 1, 0, len(EB) - 2)
    di = np.clip(np.searchsorted(DB, dist, side="right") - 1, 0, len(DB) - 2)
    np.add.at(ncell, (ei, di), 1)
    np.add.at(nfail, (ei, di), fail)
    np.add.at(onfault, ei, hasface)
    del t, b, p, e, dx, d6, eta, cen, vb, fail, dist, hasface, ei, di

tot = ncell.sum()
print(f"\n[mesh] {NT:,} tets   gate failures {nfail.sum():,} ({100*nfail.sum()/tot:.3f} %)\n")
hdr = ["0-250 m", "250-1k", "1-5 km", ">5 km"]
print(f"{'eta band':>12} {'cells':>13} {'% mesh':>8} {'fail':>11} {'fail rate':>10}   " +
      "  ".join(f"{h:>11}" for h in hdr))
for i in range(len(EB) - 1):
    c = ncell[i].sum(); fl = nfail[i].sum()
    if c == 0:
        continue
    band = f"{EB[i]:.2f}-{EB[i+1]:.2f}"
    print(f"{band:>12} {c:>13,} {100*c/tot:7.3f}% {fl:>11,} {100*fl/max(c,1):9.2f}%   " +
          "  ".join(f"{ncell[i][j]:>11,}" for j in range(len(DB) - 1)))
print()
print(f"{'eta band':>12} {'tets with a FACE ON the fault (unreshapeable)':>50}")
for i in range(len(EB) - 1):
    if ncell[i].sum() == 0:
        continue
    band = f"{EB[i]:.2f}-{EB[i+1]:.2f}"
    print(f"{band:>12} {onfault[i]:>20,}  ({100*onfault[i]/max(ncell[i].sum(),1):.2f} % of the band)")
print()
lo = ncell[0].sum() + ncell[1].sum()          # eta < 0.1
print(f"VERDICT: gate failures in eta<0.1 cells: {nfail[0].sum()+nfail[1].sum():,} "
      f"(these need RE-SHAPING, not bisection)")
print(f"         gate failures in eta>0.7 cells: {nfail[3].sum()+nfail[4].sum():,} "
      f"(these need BISECTION)")
print(f"         slivers (eta<0.1) total: {lo:,}")
