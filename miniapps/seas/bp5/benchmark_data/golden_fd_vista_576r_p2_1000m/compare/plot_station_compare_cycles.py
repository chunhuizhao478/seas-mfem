#!/usr/bin/env python3
"""BP5-FD multi-cycle station comparison: our run vs BICyclE and SBEM.

Usage: python3 -I plot_station_compare_cycles.py RUN_DIR BICYCLE_DIR SBEM_DIR OUT.pdf OUT_TABLE.txt [LABEL]

Pages: the whole run (log10 V, slip, tau vs t in yr; references drawn over
the same span and beyond, not cut), then one page per event k (log10 V at
the 10 stations against t - t0_k, each code aligned at its own event start
t0_k = first time after the previous event that any station reaches
log10 V >= -1), then the event-k slip at the 10 stations.
"""
import glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

YR = 365.25 * 24 * 3600.0
STATIONS = [(-36, 0), (-16, 0), (0, 0), (16, 0), (36, 0),
            (-24, 10), (-16, 10), (0, 10), (16, 10), (0, 22)]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p or line.lstrip().startswith("#"):
                continue
            try:
                r = [float(v) for v in p[:8]]
            except ValueError:
                continue
            if len(r) == 8:
                rows.append(r)
    a = np.array(rows)
    return a[np.argsort(a[:, 0], kind="stable")]

def files(d, pat, rx):
    out = {}
    for p in glob.glob(os.path.join(d, pat)):
        m = re.search(rx, os.path.basename(p))
        if m:
            out[(int(m.group(1)), int(m.group(2)))] = p
    return out

def name(s):
    return "strk%+03ddp%+03d" % s

def style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=7)

def event_starts(data, key, tmax):
    """Times where any station first reaches log10 V >= -1, events >= 50 yr apart."""
    cand = []
    for s in STATIONS:
        a = data[s][key]
        up = a[:, 3] >= -1.0
        i = np.nonzero(up[1:] & ~up[:-1])[0] + 1
        cand += list(a[i, 0])
        if up[0]:
            cand.append(a[0, 0])
    cand = sorted(c for c in cand if c <= tmax)
    ev = []
    for c in cand:
        if not ev or c - ev[-1] > 50 * YR:
            ev.append(c)
    return ev

def main():
    run_dir, bic_dir, sbem_dir, out, table = sys.argv[1:6]
    label = sys.argv[6] if len(sys.argv) > 6 else "our run"
    series = (("bic", "#eb6834", 1.0, "BICyclE"),
              ("sbem", "#1baf7a", 1.0, "SBEM"),
              ("run", "#2a78d6", 1.6, label))
    src = {"run": files(run_dir, "*fltst_strk*dp*.txt", r"strk([+-]\d+)dp([+-]\d+)"),
           "bic": files(bic_dir, "*", r"x2_(-?\d+)_x3_(\d+)"),
           "sbem": files(sbem_dir, "*", r"x2_(-?\d+)_x3_(\d+)")}
    data = {s: {k: load(f[s]) for k, f in src.items() if s in f} for s in STATIONS}
    t_now = min(data[s]["run"][-1, 0] for s in STATIONS)
    ev = {k: event_starts(data, k, t_now if k == "run" else 1e30) for k in ("run", "bic", "sbem")}
    nev = len(ev["run"])
    xmax_yr = (t_now / YR) * 1.05

    with PdfPages(out) as pdf:
        for col, lab, ylim in ((3, "log10 V_strike (m/s)", (-14, 1.2)),
                               (1, "slip_strike (m)", None),
                               (5, "tau_strike (MPa)", None)):
            fig, axs = plt.subplots(2, 5, figsize=(16, 6.4), sharex=True)
            for ax, s in zip(axs.flat, STATIONS):
                for key, c, lw, lb in series:
                    a = data[s][key]
                    m = a[:, 0] <= xmax_yr * YR
                    ax.plot(a[m, 0] / YR, a[m, col], color=c, lw=lw, label=lb)
                ax.set_xlim(0, xmax_yr)
                if ylim: ax.set_ylim(*ylim)
                ax.set_title(name(s), fontsize=9, color=INK)
                style(ax)
            for ax in axs[1]: ax.set_xlabel("time (yr)", fontsize=8, color=INK2)
            for ax in axs[:, 0]: ax.set_ylabel(lab, fontsize=8, color=INK2)
            axs[0, 0].legend(fontsize=7, frameon=False)
            fig.suptitle("%s, 0-%.0f yr (%s to %.1f yr; references not cut)" % (lab, xmax_yr, label, t_now / YR),
                         fontsize=11, color=INK)
            fig.text(0.5, 0.005, "Event starts (yr): %s %s; BICyclE %s; SBEM %s" % (
                label, ", ".join("%.1f" % (e / YR) for e in ev["run"]),
                ", ".join("%.1f" % (e / YR) for e in ev["bic"][:nev + 1]),
                ", ".join("%.1f" % (e / YR) for e in ev["sbem"][:nev + 1])),
                ha="center", fontsize=7, color=INK2)
            fig.tight_layout(rect=(0, 0.03, 1, 0.95))
            pdf.savefig(fig); plt.close(fig)

        for k in range(nev):
            for col, lab, inc, ylim in ((3, "log10 V_strike (m/s)", False, (-12, 1.2)),
                                        (1, "event slip (m)", True, None)):
                fig, axs = plt.subplots(2, 5, figsize=(16, 6.4), sharex=True)
                for ax, s in zip(axs.flat, STATIONS):
                    for key, c, lw, lb in series:
                        if k >= len(ev[key]): continue
                        a = data[s][key]; t0 = ev[key][k]
                        t = a[:, 0] - t0
                        m = (t >= -20) & (t <= 120)
                        y = a[m, col]
                        if inc: y = y - np.interp(t0 - 20, a[:, 0], a[:, 1])
                        ax.plot(t[m], y, color=c, lw=lw, label=lb)
                    ax.set_xlim(-20, 120)
                    if ylim: ax.set_ylim(*ylim)
                    ax.set_title(name(s), fontsize=9, color=INK)
                    style(ax)
                for ax in axs[1]: ax.set_xlabel("t - t0 (s)", fontsize=8, color=INK2)
                for ax in axs[:, 0]: ax.set_ylabel(lab, fontsize=8, color=INK2)
                axs[0, 0].legend(fontsize=7, frameon=False)
                fig.suptitle("Event %d, %s, aligned at each code's own start t0 (%s %.2f yr, BICyclE %s, SBEM %s)" % (
                    k + 1, lab, label, ev["run"][k] / YR,
                    "%.2f yr" % (ev["bic"][k] / YR) if k < len(ev["bic"]) else "-",
                    "%.2f yr" % (ev["sbem"][k] / YR) if k < len(ev["sbem"]) else "-"),
                    fontsize=11, color=INK)
                fig.tight_layout(rect=(0, 0, 1, 0.95))
                pdf.savefig(fig); plt.close(fig)

    with open(table, "w") as f:
        f.write("# %s vs BICyclE vs SBEM, latest %s time %.6e s (%.2f yr)\n" % (label, label, t_now, t_now / YR))
        f.write("# event start t0 (yr, first station log10 V >= -1) and recurrence (yr)\n")
        for k in range(max(nev, 1)):
            row = []
            for key in ("run", "bic", "sbem"):
                e = ev[key]
                if k < len(e):
                    r = "%.2f" % ((e[k] - e[k - 1]) / YR) if k > 0 else "-"
                    row.append("%.2f (rec %s)" % (e[k] / YR, r))
                else:
                    row.append("-")
            f.write("event %d: %s | BICyclE %s | SBEM %s\n" % (k + 1, row[0], row[1], row[2]))
        f.write("# event slip (m), t0-20 s to t0+300 s: %s / BICyclE / SBEM\n" % label)
        for s in STATIONS:
            cells = []
            for k in range(nev):
                v = []
                for key in ("run", "bic", "sbem"):
                    a = data[s][key]
                    if k < len(ev[key]):
                        t0 = ev[key][k]
                        v.append("%.3f" % (np.interp(t0 + 300, a[:, 0], a[:, 1]) - np.interp(t0 - 20, a[:, 0], a[:, 1])))
                    else:
                        v.append("-")
                cells.append("ev%d %s" % (k + 1, "/".join(v)))
            f.write("%s  %s\n" % (name(s), "  ".join(cells)))
        f.write("# total slip (m) at %.2f yr: %s / BICyclE / SBEM\n" % (t_now / YR, label))
        for s in STATIONS:
            v = ["%.3f" % np.interp(t_now, data[s][k][:, 0], data[s][k][:, 1]) for k in ("run", "bic", "sbem")]
            f.write("%s  %s\n" % (name(s), " / ".join(v)))
    print(open(table).read())

if __name__ == "__main__":
    main()
