#!/usr/bin/env python3
# The bend convergence test's comparison (Chunhui Zhao, 2026-10-30: "prepare
# the bend convergence test"): for each run of run_bends_mac.sh, how fast
# sigma_n' falls with slip at the bend, and how that rate changes with h.
#
# Usage (from miniapps/seas):
#   python3 bend_convergence/scripts/analyse_bends.py RESULTS [--out DIR]
# RESULTS holds the runs <tag>_<commit>/ (run_bends_mac.sh), tag = sharp05,
# arc05L05 or arc05L10 with _h500, _h250 or _h125.  Writes DIR/summary.txt,
# DIR/summary.json and DIR/bend_convergence.pdf (DIR default
# RESULTS/compare).
#
# Inputs per run: <prefix>_fault_snapshots.vtkhdf (one line cell per fault
# face, its two points, so each vertex appears once per face: the two sides
# of a kink keep their own sigma_n'), RUN_INFO.txt and run.log.
#
# Per fault point j (sigma_n' compression positive, MPa; D the strike slip
# from t = 0, sign chosen so that the largest |D| is positive):
#   k_j    the least-squares slope of sigma_n'_j against D_j over the
#          snapshots with D_j >= D_FIT (0.02 m), at least 3 of them [MPa/m];
#          k < 0: sigma_n' falls as the bend slips
#   the bend points: the sharp bend's vertex (0, 0), an arc's mid point
#          (s = 0) and its two ends (s = -+arc / 2)
# Reported per run: end state (exit code, the abort line of the log, the
# last snapshot's time and regime), the lowest sigma_n' of the run (where,
# the slip there, k there), the most negative k along the fault (where), k
# at the bend points, and sigma_n' at D = 0.25 m at the lowest point.
# The test (README.md): a sharp bend whose min k grows about as 1 / h is
# set by the mesh; an arc whose k profile stays the same as h halves is
# converged.
import argparse
import glob
import json
import math
import os
import re
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

D_FIT = 0.02      # m: the slip from which k is fitted
D_COMMON = 0.25   # m: sigma_n' at this slip, every run
THETA = math.radians(5.0)
SHAPES = {"sharp05": ("sharp 5 deg", 0.0), "arc05L05": ("5 deg over a 5 km arc", 5000.0),
          "arc05L10": ("5 deg over a 10 km arc", 10000.0)}
HS = (500, 250, 125)
YR = 365.0 * 86400.0
# one hue, light to dark as h falls (the order is the data's)
COL = {500: "#86b6ef", 250: "#2a78d6", 125: "#0d366b"}

ap = argparse.ArgumentParser()
ap.add_argument("results")
ap.add_argument("--out", default=None)
A = ap.parse_args()
out = A.out or os.path.join(A.results, "compare")
os.makedirs(out, exist_ok=True)


def bend_points(arc):
    """The bend's named points (x, y) in m: mesh/bend.geo's curve."""
    if arc == 0.0:
        return {"vertex": (0.0, 0.0)}
    R, s0 = arc / THETA, -arc / 2.0
    p = -s0 / R
    return {"arc start": (s0, 0.0),
            "arc mid": (s0 + R * math.sin(p), R * (1.0 - math.cos(p))),
            "arc end": (s0 + R * math.sin(THETA), R * (1.0 - math.cos(THETA)))}


def read_run(d, tag):
    fn = glob.glob(os.path.join(d, "*_fault_snapshots.vtkhdf"))
    if len(fn) != 1:
        return None, f"{len(fn)} fault snapshot files"
    f = h5py.File(fn[0], "r")
    g = f["VTKHDF"]
    n = int(np.ravel(g["NumberOfPoints"])[0])
    X = g["Points"][:n, :2]
    t = g["Steps"]["Values"][:]
    off = g["Steps"]["PointDataOffsets"]
    get = lambda k: np.array([g["PointData"][k][int(o):int(o) + n].ravel() for o in off[k][:]])
    S, D = get("normal_stress"), get("slip_strike")
    reg = f["fd_cycle"]["regime"][:] if "fd_cycle" in f else np.zeros(len(t), int)
    D = D - D[0]
    if abs(D.min()) > abs(D.max()):
        D = -D
    return dict(X=X, t=t, S=S, D=D, reg=reg), None


def slopes(S, D):
    k = np.full(S.shape[1], np.nan)
    for j in range(S.shape[1]):
        m = D[:, j] >= D_FIT
        if m.sum() >= 3:
            k[j] = np.polyfit(D[m, j], S[m, j], 1)[0]
    return k


def end_state(d):
    info = open(os.path.join(d, "RUN_INFO.txt")).read() if os.path.exists(os.path.join(d, "RUN_INFO.txt")) else ""
    log = open(os.path.join(d, "run.log"), errors="replace").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    rc = re.search(r"exit (-?\d+)", info)
    ab = re.search(r"(MFEM abort[^\n]{0,160}|Verification failed[^\n]{0,160})", log)
    return (int(rc.group(1)) if rc else None), (ab.group(1).strip() if ab else "none")


runs, lines = {}, []
for name in SHAPES:
    for h in HS:
        tag = f"{name}_h{h}"
        ds = sorted(glob.glob(os.path.join(A.results, tag + "_*/")))
        if not ds:
            lines.append(f"{tag}: no run folder")
            continue
        d = ds[-1]
        r, err = read_run(d, tag)
        if r is None:
            lines.append(f"{tag}: {err} in {d}")
            continue
        X, S, D, t = r["X"], r["S"], r["D"], r["t"]
        k = slopes(S, D)
        i, j = np.unravel_index(np.argmin(S), S.shape)
        jk = int(np.nanargmin(k)) if np.isfinite(k).any() else -1
        u = np.maximum.accumulate(D[:, j])
        s_common = float(np.interp(D_COMMON, u, S[:, j])) if u[-1] > D_COMMON else float("nan")
        k_sec = float((S[i, j] - S[0, j]) / D[i, j]) if D[i, j] > 1e-6 else float("nan")
        bp = {}
        for lab, (x, y) in bend_points(SHAPES[name][1]).items():
            dist = np.hypot(X[:, 0] - x, X[:, 1] - y)
            near = np.nonzero(dist <= dist.min() + 1.0)[0]   # the point once per face
            bp[lab] = [dict(x_km=float(X[q, 0] / 1e3), y_km=float(X[q, 1] / 1e3),
                            k=float(k[q]), D_end=float(D[-1, q]), S_min=float(S[:, q].min())) for q in near]
        rc, ab = end_state(d)
        runs[tag] = dict(dir=d, h=h, shape=name, exit=rc, abort=ab, snapshots=len(t),
                         t_end_yr=float(t[-1] / YR), t_end_s=float(t[-1]), regime_end=int(r["reg"][-1]),
                         S_min=float(S[i, j]), S_min_x_km=float(X[j, 0] / 1e3), S_min_y_km=float(X[j, 1] / 1e3),
                         S_min_t_s=float(t[i]), D_at_S_min=float(D[i, j]), k_at_S_min=float(k[j]),
                         k_sec_at_S_min=k_sec,
                         k_min=float(k[jk]) if jk >= 0 else float("nan"),
                         k_min_x_km=float(X[jk, 0] / 1e3) if jk >= 0 else float("nan"),
                         S_at_D_common=s_common, bend=bp,
                         profile=dict(x_km=(X[:, 0] / 1e3).tolist(), k=[None if not np.isfinite(v) else float(v) for v in k]),
                         curve=dict(D=D[:, j].tolist(), S=S[:, j].tolist()))
        q = runs[tag]
        lines.append(f"{tag}: exit {rc}, abort: {ab}")
        lines.append(f"  {len(t)} snapshots to t = {t[-1]:.6e} s ({t[-1] / YR:.3f} yr), last regime {q['regime_end']}")
        lines.append(f"  lowest sigma_n' {q['S_min']:.3f} MPa at ({q['S_min_x_km']:.3f}, {q['S_min_y_km']:.3f}) km, "
                     f"t = {q['S_min_t_s']:.3f} s, slip there {q['D_at_S_min']:.3f} m; there k {q['k_at_S_min']:.3f}, "
                     f"secant (sigma_n' - sigma_n'(0)) / slip {k_sec:.3f} MPa/m, sigma_n' at D = {D_COMMON} m "
                     f"{s_common:.3f} MPa")
        lines.append(f"  most negative k {q['k_min']:.3f} MPa/m at x = {q['k_min_x_km']:.3f} km")
        for lab, ps in bp.items():
            lines.append(f"  {lab}: " + "; ".join(f"({p['x_km']:.3f}, {p['y_km']:.3f}) km k {p['k']:.3f} MPa/m, "
                                                   f"D_end {p['D_end']:.3f} m, min sigma_n' {p['S_min']:.3f} MPa" for p in ps))

# the convergence tables against h, with the ratio per halving of h (about
# 2: set by the mesh; about 1: converged)
for key, what in (("k_min", "most negative fitted k along the fault"),
                  ("k_sec_at_S_min", "secant k at the run's lowest sigma_n'")):
    lines.append("")
    lines.append(f"{what} [MPa/m] by h (x ratio to the coarser h)")
    for name in SHAPES:
        row, prev = [], None
        for h in HS:
            q = runs.get(f"{name}_h{h}")
            v = q[key] if q else float("nan")
            row.append(f"h {h}: {v:.3f}" + (f" (x{v / prev:.2f})" if prev and np.isfinite(v) else ""))
            prev = v if np.isfinite(v) and v != 0 else None
        lines.append(f"  {name}: " + ", ".join(row))
txt = "\n".join(lines) + "\n"
open(os.path.join(out, "summary.txt"), "w").write(txt)
json.dump({k: {a: b for a, b in v.items() if a not in ("profile", "curve")} for k, v in runs.items()},
          open(os.path.join(out, "summary.json"), "w"), indent=1)
print(txt, end="")

plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.edgecolor": "#55646e",
                     "axes.labelcolor": "#18232b", "xtick.color": "#55646e", "ytick.color": "#55646e",
                     "axes.spines.top": False, "axes.spines.right": False})
with PdfPages(os.path.join(out, "bend_convergence.pdf")) as pdf:
    # 1: k along the fault, one panel per shape
    fig, axs = plt.subplots(3, 1, figsize=(8.5, 9.5), sharex=True, constrained_layout=True)
    for ax, name in zip(axs, SHAPES):
        for h in HS:
            q = runs.get(f"{name}_h{h}")
            if not q:
                continue
            x = np.array(q["profile"]["x_km"])
            kk = np.array([np.nan if v is None else v for v in q["profile"]["k"]])
            o = np.argsort(x, kind="stable")
            ax.plot(x[o], kk[o], color=COL[h], lw=1.6, marker="o", ms=2.5, label=f"h = {h} m")
        for lab, (bx, by) in bend_points(SHAPES[name][1]).items():
            ax.axvline(bx / 1e3, color="#9aa6ad", lw=0.8, ls=":")
        ax.axhline(0.0, color="#9aa6ad", lw=0.8)
        ax.set_title(SHAPES[name][0], loc="left", color="#18232b")
        ax.set_ylabel("k = d sigma_n' / d slip (MPa/m)")
        ax.legend(frameon=False, loc="lower left")
    axs[-1].set_xlabel("x (km); dotted: the bend points")
    fig.suptitle(f"sigma_n' against slip along the fault (fit over slip >= {D_FIT} m)", x=0.02, ha="left")
    pdf.savefig(fig); plt.close(fig)
    # 2: the most negative k against h
    fig, ax = plt.subplots(figsize=(6.5, 4.8), constrained_layout=True)
    for name, mk in zip(SHAPES, ("o", "s", "^")):
        hh = [h for h in HS if f"{name}_h{h}" in runs and np.isfinite(runs[f"{name}_h{h}"]["k_min"])]
        if hh:
            ax.loglog(hh, [-runs[f"{name}_h{h}"]["k_min"] for h in hh], marker=mk, ms=8, lw=1.6,
                      color="#18232b", mfc="white", label=SHAPES[name][0])
    q = runs.get("sharp05_h500")
    if q and np.isfinite(q["k_min"]):
        ax.loglog(HS, [-q["k_min"] * 500.0 / h for h in HS], color="#9aa6ad", lw=1.0, ls="--", label="1 / h from sharp, h = 500 m")
    ax.set_xlabel("h on the fault (m)"); ax.set_ylabel("-(most negative k) (MPa/m)")
    ax.set_xticks(HS); ax.set_xticklabels([str(h) for h in HS]); ax.minorticks_off()
    ax.invert_xaxis(); ax.legend(frameon=False)
    ax.set_title("A mesh-set response grows as 1 / h; a converged one stays flat", loc="left", color="#18232b")
    pdf.savefig(fig); plt.close(fig)
    # 3: sigma_n' against slip at each run's lowest point
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True, constrained_layout=True)
    for ax, name in zip(axs, SHAPES):
        for h in HS:
            q = runs.get(f"{name}_h{h}")
            if q:
                ax.plot(q["curve"]["D"], q["curve"]["S"], color=COL[h], lw=1.6, label=f"h = {h} m")
        ax.axhline(0.0, color="#9aa6ad", lw=0.8)
        ax.set_title(SHAPES[name][0], loc="left", color="#18232b"); ax.set_xlabel("slip (m)")
        ax.legend(frameon=False)
    axs[0].set_ylabel("sigma_n' at the run's lowest point (MPa)")
    pdf.savefig(fig); plt.close(fig)
print(f"wrote {out}/summary.txt, summary.json, bend_convergence.pdf")
