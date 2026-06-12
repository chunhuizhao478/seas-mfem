#!/usr/bin/env python3
"""Render the Eq.(18) fault volume-ratio r_v (baseline vs rvsmooth) as a
before/after map on the TPV102 fault (the y=0 plane -> x-z view).

Reads the two ASCII VTUs written by export_eq18_fault_vtu.py (per-triangle
`rv` cell data), draws a shared-colormap tripcolor with the Zhang (2023)
r_v = 1.5 SSO ceiling marked.

Usage:  conda activate pythonenv && python plot_rv_before_after.py
"""
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

HERE = __import__("pathlib").Path(__file__).resolve().parent
PAIRS = [("baseline tpv102_200m", HERE / "tpv102_200m_fault_rv.vtu"),
         ("rvsmooth tpv102_200m_rvsmooth", HERE / "tpv102_200m_rvsmooth_fault_rv.vtu")]
OUT = HERE / "tpv102_200m_rv_before_after.png"


def read_vtu(path):
    """Return (pts Nx3, tris Mx3, rv M) from a minimal ascii UnstructuredGrid VTU."""
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")
    def arr(parent, name=None):
        for da in parent.findall("DataArray"):
            if name is None or da.get("Name") == name:
                return np.fromstring(da.text.strip().replace("\n", " "), sep=" ")
        raise KeyError(name)
    pts = arr(piece.find("Points")).reshape(-1, 3)
    cells = piece.find("Cells")
    conn = arr(cells, "connectivity").astype(int).reshape(-1, 3)
    rv = arr(piece.find("CellData"), "rv")
    return pts, conn, rv


fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharex=True, sharey=True)
vmin, vmax = 1.0, 2.2          # shared scale across both panels
for ax, (label, path) in zip(axes, PAIRS):
    pts, tris, rv = read_vtu(path)
    # Fault is the y=0 plane: plot in (x, z), both in km.
    x, z = pts[:, 0] / 1e3, pts[:, 2] / 1e3
    triang = Triangulation(x, z, tris)
    tpc = ax.tripcolor(triang, facecolors=rv, cmap="turbo", vmin=vmin, vmax=vmax)
    n_bad = int((rv > 1.5).sum())
    ax.set_title(f"{label}\nmax r_v={rv.max():.3f}  mean={rv.mean():.3f}  "
                 f"(#faces r_v>1.5: {n_bad})", fontsize=11)
    ax.set_xlabel("along-strike x [km]")
    ax.set_aspect("equal")
    cb = fig.colorbar(tpc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("r_v = max(V_A/V_B, V_B/V_A)")
    # mark the 1.5 SSO ceiling on the colorbar
    cb.ax.axhline((1.5 - vmin) / (vmax - vmin), color="black", lw=1.5, ls="--")
axes[0].set_ylabel("depth z [km]")
fig.suptitle("TPV102 200 m fault — Eq.(18) volume ratio r_v  "
             "(Zhang 2023 SSO ceiling = 1.5, dashed)", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUT, dpi=130)
print(f"wrote {OUT}")
