#!/usr/bin/env python3
"""Figures for a UCERF 3-D mesh made by build_ucerf3d_mesh.py.

  <stem>_surfaces.png  3-D view of the friction (attribute 3) and creep
                       (attribute 6) faces, and the San Gorgonio Pass zoom.
  <stem>_gaps.png      gap between saf2 and saf1 / saf3 against depth,
                       before and after the cut-back, with the gap target.

Usage: plot_ucerf3d_mesh.py MESH.msh [--json MESH.json] [--outdir DIR]
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_ucerf3d_mesh as B          # noqa: E402
from check_ucerf3d_mesh import read_msh22   # noqa: E402

COL = {"saf1": "#d9534f", "saf2": "#e69f00", "saf3": "#b2182b", "sjf1": "#3a7bd5", "sjf2": "#2166ac"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("msh")
    ap.add_argument("--json")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    rec = json.load(open(a.json or os.path.splitext(a.msh)[0] + ".json"))
    stem = os.path.join(a.outdir, os.path.splitext(os.path.basename(a.msh))[0])
    X, tris, _ = read_msh22(a.msh)
    A = rec["attributes"]
    ent2chain = {}
    for n, t in rec["gmsh_surface_tags"].items():
        for e in t["friction"] + t["creep"]:
            ent2chain[e] = n

    # ---------------- surfaces
    fig = plt.figure(figsize=(13, 5.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    for attr, alpha in ((A["friction"], 0.9), (A["creep"], 0.35)):
        sel = tris[tris[:, 0] == attr]
        for n in COL:
            s = sel[[ent2chain.get(e) == n for e in sel[:, 1]]]
            if not len(s):
                continue
            P = X[s[:, 2:5]] / 1e3
            ax.add_collection3d(Poly3DCollection(P, facecolor=COL[n] if attr == A["friction"] else "0.6",
                                                 edgecolor="k", linewidth=0.05, alpha=alpha))
    allp = X[tris[np.isin(tris[:, 0], [A["friction"], A["creep"]])][:, 2:5]].reshape(-1, 3) / 1e3
    lo, hi = allp.min(0), allp.max(0)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], 0)
    ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1], 6 * (0 - lo[2])))
    ax.view_init(elev=28, azim=-115)
    ax.set_xlabel("x east (km)"); ax.set_ylabel("y north (km)"); ax.set_zlabel("z (km)")
    ax.set_title("friction faces (colour) and creep faces (grey), vertical x6", fontsize=9)

    # zoom at San Gorgonio Pass: map view of the slices at 0, 10, 25 km
    ax2 = fig.add_subplot(1, 2, 2)
    ch = B.name_chains(B.read_chains_2d(B.SRC_2D))
    full2 = B.Surface("saf2", ch["saf2"], B.SAF2_DIP_DEG, B.SAF2_DIPDIR_DEG)
    d_lo = rec["depths_km"]["creep_bottom"] * 1e3
    for d, ls in ((0.0, "-"), (10e3, "--"), (d_lo, ":")):
        for n in ("saf1", "saf3"):
            c = np.array(rec["chains"][n]["top_trace_m"]) / 1e3
            ax2.plot(c[:, 0], c[:, 1], color=COL[n], lw=1.2)
        s2 = B.Surface("saf2", np.array(rec["chains"]["saf2"]["top_trace_m"]), B.SAF2_DIP_DEG,
                       B.SAF2_DIPDIR_DEG, ref_depth=0.0)
        p = s2.at(d) / 1e3
        q = full2.at(d) / 1e3
        ax2.plot(q[:, 0], q[:, 1], color="0.6", lw=0.8, ls=ls)
        ax2.plot(p[:, 0], p[:, 1], color=COL["saf2"], lw=1.6, ls=ls, label=f"saf2 at {d / 1e3:.0f} km")
    c2 = full2.chain / 1e3
    ax2.set_xlim(c2[:, 0].min() - 20, c2[:, 0].max() + 20)
    ax2.set_ylim(c2[:, 1].min() - 25, c2[:, 1].max() + 20)
    ax2.set_aspect("equal")
    ax2.set_xlabel("x east (km)"); ax2.set_ylabel("y north (km)")
    cut = rec["saf2_cut_km"]
    ax2.set_title(f"San Gorgonio Pass: saf2 slices (grey: before the cut-back);\n"
                  f"cut {cut['saf3_end']:.2f} km at the saf3 end, {cut['saf1_end']:.2f} km at the saf1 end", fontsize=9)
    ax2.legend(fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(stem + "_surfaces.png", dpi=160)
    plt.close()

    # ---------------- gaps
    S = {n: B.Surface(n, c) for n, c in ch.items()}
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    cut2 = B.Surface("saf2", np.array(rec["chains"]["saf2"]["top_trace_m"]), B.SAF2_DIP_DEG,
                     B.SAF2_DIPDIR_DEG, ref_depth=0.0)
    for n in ("saf1", "saf3"):
        before = B.min_gap(full2, S[n], 0.0, d_lo)
        after = B.min_gap(cut2, S[n], 0.0, d_lo)
        ax.plot([g / 1e3 for _, g in before], [-d / 1e3 for d, _ in before], color=COL[n], lw=0.8, ls="--",
                label=f"saf2-{n} before")
        ax.plot([g / 1e3 for _, g in after], [-d / 1e3 for d, _ in after], color=COL[n], lw=1.6,
                label=f"saf2-{n} after")
    ax.axvline(rec["gap_target_km"], color="k", lw=0.8, ls=":", label=f"target {rec['gap_target_km']:.1f} km")
    ax.axhline(-rec["depths_km"]["friction"], color="0.5", lw=0.6)
    ax.set_xlabel("smallest gap (km)"); ax.set_ylabel("z (km)")
    ax.set_xlim(0, 20)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("saf2 gap to its neighbours, every 0.5 km in depth", fontsize=9)
    plt.tight_layout()
    plt.savefig(stem + "_gaps.png", dpi=160)
    plt.close()
    print(stem + "_surfaces.png", stem + "_gaps.png")


if __name__ == "__main__":
    main()
