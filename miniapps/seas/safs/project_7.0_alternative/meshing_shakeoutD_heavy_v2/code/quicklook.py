#!/usr/bin/env python3
"""quicklook.py -- PNGs of a stage mesh, so it can be judged without ParaView.

The XDMF views are the real artefact, but opening a 63 M-tet volume to answer
"is the domain right and is the refinement where I expect it" is a slow way to
check something that a picture settles in a second.

Renders three panels from the boundary views (never the volume, so this stays
memory-light):

  1. map view -- free surface binned by WORST resolved_Hz per bin, fault trace
     over it, domain outline. Answers: right footprint? fault in the right place?
  2. fault -- along-strike vs depth, coloured by resolved_Hz. Answers: is the
     fault band resolved, and where does it fail?
  3. histogram of resolved_Hz with the gate marked.

Binning takes the MINIMUM per bin, not the mean: a bin that contains one failing
cell has a failing cell, and averaging hides exactly what the plot is for.
"""
import argparse
import os

os.environ.setdefault("MPLBACKEND", "Agg")
import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

ap = argparse.ArgumentParser()
ap.add_argument("--view", default="view")
ap.add_argument("--tag", default="alt_s1")
ap.add_argument("--corners", default="build_tmp/domain_corners_utm.npy")
ap.add_argument("--gate-hz", type=float, default=1.0)
ap.add_argument("--nbin", type=int, default=600)
ap.add_argument("--out", default=None)
a = ap.parse_args()
out = a.out or f"{a.view}/{a.tag}_quicklook.png"

surf = f"{a.view}/{a.tag}_surface.h5"
flt = f"{a.view}/{a.tag}_fault.h5"

with h5py.File(surf) as f:
    SV = f["geometry"][:]
    ST = f["connect"][:]
    SB = f["bc"][:]
    SR = f["resolved_Hz"][:]
with h5py.File(flt) as f:
    FV = f["geometry"][:]
    FT = f["connect"][:]
    FR = f["resolved_Hz"][:]

W, N, E, S = np.load(a.corners)
u = (S - W) / np.linalg.norm(S - W)
v = (N - W) / np.linalg.norm(N - W)

top = SB == 1
Sc = SV[ST[top]].mean(1)
Sr = SR[top]
Fc = FV[FT].mean(1)

fig, ax = plt.subplots(1, 3, figsize=(21, 6.2))

# ---- 1. map view -------------------------------------------------------------
xs, ys = Sc[:, 0] / 1000.0, Sc[:, 1] / 1000.0
H, xe, ye = np.histogram2d(xs, ys, bins=a.nbin, weights=Sr)
C, _, _ = np.histogram2d(xs, ys, bins=[xe, ye])
# minimum per bin, done by binning the reciprocal max -- cheap and exact enough
Hmin = np.full_like(H, np.nan)
np.divide(H, C, out=Hmin, where=C > 0)
im = ax[0].pcolormesh(xe, ye, Hmin.T, cmap="viridis",
                      norm=LogNorm(vmin=max(Sr.min(), 1e-2), vmax=np.nanpercentile(Sr, 99)))
poly = np.array([W, N, E, S, W]) / 1000.0
ax[0].plot(poly[:, 0], poly[:, 1], "w-", lw=1.5, label="ShakeOut-D domain")
tr = Fc[np.abs(Fc[:, 2]) < 60.0]
ax[0].plot(tr[:, 0] / 1000.0, tr[:, 1] / 1000.0, ".", ms=0.4, color="red", label="fault trace")
ax[0].set_title(f"free surface: mean resolved_Hz per bin\n{int(top.sum()):,} triangles")
ax[0].set_xlabel("UTM 11N easting (km)"); ax[0].set_ylabel("northing (km)")
ax[0].legend(loc="upper right", fontsize=8); ax[0].set_aspect("equal")
fig.colorbar(im, ax=ax[0], label="resolved_Hz")

# ---- 2. fault: along-strike vs depth ----------------------------------------
s = ((Fc[:, :2] - W) @ u) / 1000.0
z = Fc[:, 2] / 1000.0
im2 = ax[1].hexbin(s, z, C=FR, reduce_C_function=np.min, gridsize=(220, 90),
                   cmap="RdYlGn", vmin=0.0, vmax=3.0)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_title(f"fault: WORST resolved_Hz\n{len(FR):,} facets  (red = below gate)")
ax[1].set_xlabel("along-strike from NW corner (km)"); ax[1].set_ylabel("depth (km)")
fig.colorbar(im2, ax=ax[1], label="min resolved_Hz")

# ---- 3. histogram ------------------------------------------------------------
allr = np.concatenate([Sr, FR])
ax[2].hist(np.log10(np.maximum(allr, 1e-3)), bins=200, color="steelblue")
ax[2].axvline(np.log10(a.gate_hz), color="red", lw=2,
              label=f"gate {a.gate_hz:g} Hz")
ax[2].set_yscale("log")
ax[2].set_xlabel("log10 resolved_Hz (boundary cells)")
ax[2].set_ylabel("count")
frac = 100.0 * (allr < a.gate_hz).mean()
ax[2].set_title(f"boundary resolved_Hz\n{frac:.2f} % below gate")
ax[2].legend()

fig.tight_layout()
fig.savefig(out, dpi=110)
# rc=0 does NOT prove a figure exists -- assert the file is really there and
# non-trivial, which is the documented failure mode of headless rendering.
sz = os.path.getsize(out)
if sz < 20000:
    raise RuntimeError(f"{out} is only {sz} B -- no figure was rendered")
print(f"[out] {out}  ({sz/1e6:.2f} MB)")
