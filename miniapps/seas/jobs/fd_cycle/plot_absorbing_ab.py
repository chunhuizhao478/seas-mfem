#!/usr/bin/env python3
# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
"""plot_absorbing_ab.py — fully dynamic BP5 cycle, P5: the station pages of
two runs of the same event, reflecting ([cycle.absorbing] faces = "off",
grey, beneath) and absorbing ("outer", red, on top), as
vista_event1_absorbing.sbatch makes them.  Page 1: max slip rate, moment
rate, the growth of the difference between the runs, and per station the
times it first exceeds 1e-4 and 1e-2 against the earliest time a P wave
from the fault can come back to it off the outer boundary (mirror images
of the BP5 fault rectangle x in [-50, 50] km, z in [-40, 0] km in the box
x +-200 km, y +-100 km, z >= -100 km; Vp of the BP5 material).  The
difference of a station is the largest over its columns of |a - b| / max|a|
(log10 columns compared as values).  CONTROL_DIR, optional, is the "off"
deck on another number of ranks: its difference from OFF_DIR is round-off
growing in the event, the floor the A/B difference is read against.
Then two stations per page: on-fault slip, log10 slip rate and shear
stress along strike; off-fault velocity along strike, normal to the fault
and down.  The overview page is also written as <OUT>_overview.png.
Needs numpy and matplotlib.

  python3 plot_absorbing_ab.py OFF_DIR OUTER_DIR PREFIX OUT.pdf "title" [CONTROL_DIR]
"""
import glob, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

off_dir, out_dir, pre, pdf_path, title = sys.argv[1:6]
ctl_dir = sys.argv[6] if len(sys.argv) > 6 else None
VP = np.sqrt(3 * 32038120320.0 / 2670.0)          # lambda = mu: Vp = sqrt(3 mu / rho)
FAULT = ((-50e3, 50e3), (0.0, 0.0), (-40e3, 0.0))  # x, y, z ranges of the fault
BOUNDS = {"bottom": (2, (-100e3,)), "y walls": (1, (100e3, -100e3)),
          "x walls": (0, (200e3, -200e3))}   # axis, plane positions
TOLS = (1e-4, 1e-2)
GREY, RED, BLACK = "#9a9a9a", "#d62728", "#222222"


def read(path):
    head, rows = [], []
    for l in open(path):
        if l.startswith("#"):
            head.append(l)
        elif l.strip():
            rows.append([float(v) for v in l.split()])
    return head, np.array(rows)


def station_xyz(head):
    s = head[0]
    i = s.index("(x=")
    parts = s[i + 1:s.index(")", i)].replace(" m", "").split(", ")
    return np.array([float(p.split("=")[1]) for p in parts])


def box_dist(p, box):
    d = [max(lo - p[k], 0.0, p[k] - hi) for k, (lo, hi) in enumerate(box)]
    return float(np.sqrt(np.dot(d, d)))


def earliest(s):
    """Earliest arrival at s of a P wave from the fault off each outer boundary."""
    out = {}
    for name, (ax, planes) in BOUNDS.items():
        t = np.inf
        for P in planes:
            box = list(FAULT)
            lo, hi = FAULT[ax]
            box[ax] = tuple(sorted((2 * P - lo, 2 * P - hi)))   # mirror across x_ax = P
            t = min(t, box_dist(s, box) / VP)
        out[name] = t
    return out


def difference(a, b, logcols):
    """The station's difference on the rows both runs wrote at the same times
    (a shorter run's last row, at its stop time, is left out): the largest
    over the columns of |a - b| / max|a|, log10 columns as values."""
    n = min(len(a), len(b))
    same = a[:n, 0] == b[:n, 0]
    n = n if same.all() else int(np.argmin(same))
    assert n > 1, "the runs wrote rows at different times"
    worst = np.zeros(n)
    for c in range(1, a.shape[1]):
        x, y = a[:n, c], b[:n, c]
        if c in logcols:
            x, y = 10.0 ** x, 10.0 ** y
        sc = np.max(np.abs(x))
        worst = np.maximum(worst, np.abs(x - y) / sc if sc > 0 else np.abs(x - y))
    return a[:n, 0], worst


def crossings(t, worst):
    """The first time worst exceeds each of TOLS (None: never)."""
    out = []
    for tol in TOLS:
        k = np.nonzero(worst > tol)[0]
        out.append(t[k[0]] if len(k) else None)
    return out


def mark(ax, eb):
    ls = {"bottom": "--", "y walls": "-.", "x walls": ":"}
    for name, t in eb.items():
        ax.axvline(t, color="k", lw=0.7, ls=ls[name], alpha=0.6)


rows = []
pages = []
curves = []
for kind, cols, logcols, labels in (
        ("fltst", (1, 3, 5), {3, 4, 7}, ("slip strike (m)", "log10 V strike (m/s)", "tau strike (MPa)")),
        ("blkst", (5, 4, 6), set(), ("vel along strike (m/s)", "vel fault-normal (m/s)", "vel e3 = down (m/s)"))):
    for fo in sorted(glob.glob(f"{off_dir}/{pre}_{kind}_*.txt")):
        name = os.path.basename(fo)[len(pre) + 1:-4]
        fa = f"{out_dir}/{pre}_{name}.txt"
        if not os.path.exists(fa):
            continue
        ho, a = read(fo)
        _, b = read(fa)
        s = station_xyz(ho)
        eb = earliest(s)
        t, w = difference(a, b, logcols)
        tc = (None,) * len(TOLS)
        if ctl_dir:
            _, c = read(f"{ctl_dir}/{pre}_{name}.txt")
            tc_t, tc_w = difference(a, c, logcols)
            tc = crossings(tc_t, tc_w)
            curves.append(("control", tc_t, tc_w))
        curves.append(("ab", t, w))
        rows.append((name, s, crossings(t, w), min(eb.values()), eb, tc))
        pages.append((name, a, b, cols, labels, eb, crossings(t, w)))

with PdfPages(pdf_path) as pdf:
    # overview
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle(title, fontsize=12)
    for d, col, lab, lw in ((off_dir, GREY, "reflecting (off)", 1.0), (out_dir, RED, "absorbing (outer)", 1.4)):
        _, g = read(f"{d}/{pre}_global.txt")
        axs[0, 0].plot(g[:, 0], g[:, 1], color=col, lw=lw, label=lab)
        _, so = read(f"{d}/{pre}_source.txt")
        axs[0, 1].plot(so[:, 0], so[:, 2], color=col, lw=lw, label=lab)
    axs[0, 0].set(xlabel="t (s)", ylabel="log10 max V (m/s)", title="max slip rate on the fault")
    axs[0, 1].set(xlabel="t (s)", ylabel="moment rate (N m/s)", title="moment rate over A")
    for ax in axs[0]:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    ax = axs[1, 0]
    for kind, col, lab in (("ab", RED, "absorbing vs reflecting"),
                           ("control", BLACK, "reflecting, other rank count (round-off)")):
        cs = [(t, w) for k, t, w in curves if k == kind]
        if not cs:
            continue
        n = min(len(t) for t, _ in cs)
        ax.semilogy(cs[0][0][:n], np.max([w[:n] for _, w in cs], axis=0), color=col, lw=1.4,
                    label=lab)
    tb = min(r[3] for r in rows)
    ax.axvline(tb, color="k", lw=0.8, ls="--")
    ax.text(tb, 0.02, f" earliest reflection\n {tb:.1f} s", transform=ax.get_xaxis_transform(),
            fontsize=7, va="bottom")
    ax.set_ylim(1e-16, 10)
    ax.set(xlabel="t (s)", ylabel="largest relative difference, all stations",
           title="difference between the runs")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)
    axs[1, 1].axis("off")
    fmt = lambda v: "    -" if v is None else f"{v:5.1f}"
    txt = [f"Vp = {VP:.0f} m/s; fault x in [-50, 50] km, z in [-40, 0] km",
           "bound: earliest P reflection from the fault off the bottom",
           "(z = -100 km), y walls (+-100 km) or x walls (+-200 km);",
           "dashed / dash-dot / dotted lines on the station pages", "",
           "first time (s) the difference exceeds 1e-4 and 1e-2:",
           f"  {'station':<24s} {'A/B':>11s}" + (f" {'control':>11s}" if ctl_dir else "") + "  bound"]
    for r in rows:
        txt.append(f"  {r[0]:<24s} {fmt(r[2][0])} {fmt(r[2][1])}"
                   + (f" {fmt(r[5][0])} {fmt(r[5][1])}" if ctl_dir else "") + f"  {r[3]:5.1f}")
    axs[1, 1].text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=6.5)
    fig.tight_layout()
    pdf.savefig(fig)
    fig.savefig(pdf_path[:-4] + "_overview.png", dpi=110)
    plt.close(fig)
    # stations, two per page
    for i in range(0, len(pages), 2):
        fig, axs = plt.subplots(2, 3, figsize=(11, 8.5), squeeze=False)
        for r, (name, a, b, cols, labels, eb, td) in enumerate(pages[i:i + 2]):
            for k, (c, lab) in enumerate(zip(cols, labels)):
                ax = axs[r, k]
                ax.plot(a[:, 0], a[:, c], color=GREY, lw=1.0, label="reflecting")
                ax.plot(b[:, 0], b[:, c], color=RED, lw=1.3, label="absorbing")
                mark(ax, eb)
                ax.set(xlabel="t (s)", ylabel=lab)
                ax.grid(alpha=0.3)
            tds = "never" if td[0] is None else f"{td[0]:.2f} s"
            axs[r, 0].set_title(f"{name}: difference > 1e-4 from {tds}, earliest reflection "
                                f"{min(eb.values()):.2f} s", fontsize=8, loc="left")
            axs[r, 0].legend(fontsize=7)
        for r in range(len(pages[i:i + 2]), 2):
            for k in range(3):
                axs[r, k].axis("off")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
print(f"{pdf_path}: {len(pages)} stations")
for r in rows:
    ts = lambda v: ", ".join("never" if x is None else f"{x:.2f} s" for x in v)
    print(f"  {r[0]:<24s} > 1e-4, 1e-2 at {ts(r[2])}" + (f"; control {ts(r[5])}" if ctl_dir else "")
          + f"; bound {r[3]:.3f} s ({', '.join(f'{k} {v:.2f}' for k, v in r[4].items())})")
