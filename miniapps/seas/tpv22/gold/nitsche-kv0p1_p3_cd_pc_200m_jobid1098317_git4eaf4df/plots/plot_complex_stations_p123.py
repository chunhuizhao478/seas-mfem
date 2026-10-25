#!/usr/bin/env python3
"""TPV12 / TPV15 / TPV22 station comparison: Nitsche p=1, p=2, p=3 against one SCEC reference.

The pages of the benchmarks thread's TPV31 p1+p2+p3 PDF (plot_tpv31_stations_p123.py and
plot_misfit_maps.py, /mnt/project-files/benchmarks/results/tpv31_200m_p123_vista_plots/):
  - misfit maps (arrival, final slip, peak slip rate minus the reference), one row per run,
    highest order on top, read from a compare_nitsche_fault_stations.py table;
  - an overview of slip rate at every station (rows x3, columns x2);
  - per station (a) slip rate, (b) front zoom, (c) slip, (d) shear stress, (e) normal stress
    (autoscaled), (f) the other slip-rate component (autoscaled).
Reference black, p=1 orange dashed, p=2 green, p=3 red drawn last.  Files are read through
read_stations / find_file (by name) / load of scripts/compare_nitsche_fault_stations.py.
Changes from the TPV31 scripts (inputs only): one reference from the command line; --comp picks
the main slip component (h = strike, v = dip; TPV12 slips down dip), panel (f) shows the other;
stations named with a fault tag (main_, branch_, fault1_, fault2_) get one overview page and one
row of maps per fault; a station with no reference file is drawn with the runs alone; any run
given is drawn (a run that did not finish can be left out); a run with more than --max-points
samples (LTS writes every small step) is drawn at every k-th sample, as the TPV31 p=1 plots drew
every 10th; the misfits come from the compare table at full resolution.

Usage: plot_complex_stations_p123.py --seas DIR --prob tpv12 --comp v --t-end 8 --tag TAG --out DIR
         --ref "LABEL|DIR" --run "KEY|LABEL|DIR" [--run ...] [--cmp COMPARE_TXT]
  KEY is the run's label in COMPARE_TXT (e.g. p1, p2, p3); runs are given lowest order first.
"""
import argparse, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ap = argparse.ArgumentParser()
ap.add_argument("--seas", required=True); ap.add_argument("--prob", required=True)
ap.add_argument("--comp", choices=("h", "v"), required=True); ap.add_argument("--t-end", type=float, required=True)
ap.add_argument("--tag", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--ref", required=True); ap.add_argument("--run", action="append", required=True)
ap.add_argument("--cmp", default="")
ap.add_argument("--max-points", type=int, default=40000,
                help="draw every k-th sample of a run so that at most this many remain (LTS runs write every"
                     " small step); the misfits come from COMPARE_TXT at full resolution")
a = ap.parse_args()
sys.path.insert(0, os.path.join(a.seas, "scripts"))
from compare_nitsche_fault_stations import read_stations, find_file, load, _parse_x2x3  # noqa: E402

PROB = a.prob.upper()
C1, C2 = a.comp, ("v" if a.comp == "h" else "h")
N1, N2 = ("strike", "dip") if a.comp == "h" else ("dip", "strike")
REF0, refdir = a.ref.split("|", 1)
RUNS = [r.split("|", 2) for r in a.run]                 # (key, label, dir), lowest order first
STY = [dict(color="red", lw=1.5, ls="-"), dict(color="green", lw=1.1, ls="-"),
       dict(color="orange", lw=1.0, ls="--")]           # last run, second last, third last
SOURCES = [(REF0, REF0, refdir, "scec_mpa", dict(color="black", lw=1.0, ls="-"))]
for i, (key, lab, d) in enumerate(RUNS):
    SOURCES.append((key, lab, d, "nitsche", STY[len(RUNS) - 1 - i]))
TOP_KEY, TOP_LAB = RUNS[-1][0], RUNS[-1][1]

stations = read_stations(os.path.join(a.seas, f"{a.prob}/stations/{a.prob}_fault_stations.csv"))
data = {}
stride = {}
for key, lab, d, fmt, _ in SOURCES:
    for st in stations:
        f = find_file(d, st, True)
        if f is None:
            if fmt == "nitsche":
                sys.exit(f"missing file for {lab} station {st['name']} in {d}")
            print(f"no {REF0} file for {st['name']}: runs drawn alone")
            continue
        x = load(f, fmt)
        k = -(-len(x["t"]) // a.max_points)
        if fmt == "nitsche" and k > 1:
            x = {c: (v[::k] if isinstance(v, np.ndarray) else v) for c, v in x.items()}
            stride[key] = k
        data[(key, st["name"])] = x
for key, k in stride.items():
    print(f"{key}: every {k}th sample drawn")

# misfits per (run key, station) from the compare table: (dt, slip %, peak %) or a word
mf = {}
if a.cmp:
    for l in open(a.cmp):
        w = l.split()
        if len(w) < 2 or "x2_" not in w[0]:
            continue
        m = re.search(r"dt\s+(\S+) s\s+slip\s+(\S+)%\s+peak\s+(\S+)%", l)
        if m:
            mf[(w[1], w[0])] = tuple(float(m.group(i)) for i in (1, 2, 3))
        elif "neither side passes" in l:
            mf[(w[1], w[0])] = "quiet"


def where(st):
    p = _parse_x2x3(st["name"])
    return (p[0] or "", p[1], p[2]) if p else ("", st["x"] / 1e3 + 0.0, -st["z"] / 1e3 + 0.0)


def fname(tag):
    return tag if tag.startswith("fault") else f"{tag} fault"


def arrival(d):
    v = np.hypot(d["h-slip-rate"], d["v-slip-rate"])
    i = np.nonzero(v > 1e-3)[0]
    return float(d["t"][i[0]]) if i.size else float("nan")


def draw(ax, nm, key, t0=0.0, t1=a.t_end):
    for k, lab, _, _, sty in SOURCES:
        d = data.get((k, nm))
        if d is None or key not in d:
            continue
        m = (d["t"] >= t0 - 1e-9) & (d["t"] <= t1 + 1e-9)
        ax.plot(d["t"][m], d[key][m], label=lab, **sty)
    ax.set_xlim(t0, t1)
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=7)


faults = []
for st in stations:
    if where(st)[0] not in faults:
        faults.append(where(st)[0])
os.makedirs(a.out, exist_ok=True)
pdf_path = os.path.join(a.out, f"{a.tag}_stations.pdf")
pages = 0


def save(fig, pdf, name):
    global pages
    pdf.savefig(fig)
    pages += 1
    fig.savefig(os.path.join(a.out, f"{a.tag}_p{pages - 1:02d}_{name}.png"), dpi=110)
    plt.close(fig)


with PdfPages(pdf_path) as pdf:
    # p00: misfit maps, one row per (run, fault), highest order on top
    rows = [(k, lab, f) for k, lab, _ in reversed(RUNS) for f in faults]
    fig, axs = plt.subplots(len(rows), 3, figsize=(15, 4.6 * len(rows)), squeeze=False)
    for ri, (k, lab, fault) in enumerate(rows):
        sts = [st for st in stations if where(st)[0] == fault]
        xs = sorted({where(st)[1] for st in sts}); zs = sorted({where(st)[2] for st in sts})
        for ci, (title, lim) in enumerate([("arrival misfit (s)", 0.35), ("final slip misfit (%)", 6.0),
                                           ("peak slip rate misfit (%)", 50.0)]):
            ax = axs[ri, ci]
            G = np.full((len(zs), len(xs)), np.nan)
            txt = {}
            for st in sts:
                _, x, z = where(st); zi, xi = zs.index(z), xs.index(x)
                v = mf.get((k, st["name"]))
                if isinstance(v, tuple) and np.isfinite(v[ci]):
                    G[zi, xi] = v[ci]
                    txt[(zi, xi)] = f"{v[ci]:+.2f}" if ci == 0 else f"{v[ci]:+.1f}"
                elif v == "quiet":
                    txt[(zi, xi)] = "quiet"
                elif (REF0, st["name"]) not in data:
                    txt[(zi, xi)] = "no ref"
                elif isinstance(v, tuple):
                    ra, ka = arrival(data[(REF0, st["name"])]), arrival(data[(k, st["name"])])
                    txt[(zi, xi)] = "ref quiet" if not np.isfinite(ra) else ("run quiet" if not np.isfinite(ka) else "n/a")
                else:
                    txt[(zi, xi)] = "n/a"
            im = ax.imshow(G, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
            for (zi, xi), t in txt.items():
                v = G[zi, xi]
                ax.text(xi, zi, t, ha="center", va="center", fontsize=8.5,
                        color="white" if np.isfinite(v) and abs(v) > 0.6 * lim else "black")
            ax.set_xticks(range(len(xs)), [f"{x:g}" for x in xs]); ax.set_yticks(range(len(zs)), [f"{z:g}" for z in zs])
            ax.set_xlabel("x2 along strike (km)"); ax.set_ylabel("x3 (km)")
            ax.set_title(f"{lab}{', ' + fname(fault) if fault else ''}: {title}", fontsize=10)
            fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle(f"SCEC {PROB}: {', '.join(r[1] for r in reversed(RUNS))}\nminus {REF0} at the "
                 f"{len(stations)} on-fault stations (blue = early / low)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.6 / (4.6 * len(rows))))
    save(fig, pdf, "misfit_maps")

    # overview: main-component slip rate, rows x3, columns x2, one page per fault
    order = []
    for fault in faults:
        sts = [st for st in stations if where(st)[0] == fault]
        xs = sorted({where(st)[1] for st in sts}); zs = sorted({where(st)[2] for st in sts})
        fig, axs = plt.subplots(len(zs), len(xs), figsize=(max(9, 4.7 * len(xs)), max(6, 2.4 * len(zs) + 1)),
                                sharex=True, squeeze=False)
        for ax in axs.flat:
            ax.set_visible(False)
        for zi, z in enumerate(zs):
            for xi, x in enumerate(xs):
                st = next((s for s in sts if where(s)[1] == x and where(s)[2] == z), None)
                if st is None:
                    continue
                ax = axs[zi][xi]; ax.set_visible(True)
                order.append(st["name"])
                draw(ax, st["name"], f"{C1}-slip-rate")
                ax.set_title(st["name"], fontsize=8)
                if xi == 0 or not axs[zi][xi - 1].get_visible():
                    ax.set_ylabel(f"x3 {z:g} km\nslip rate (m/s)", fontsize=8)
                ax.set_xlabel("time (s)", fontsize=8)
        h, l = next(ax for ax in axs.flat if ax.get_visible()).get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=4, fontsize=9)
        fig.suptitle(f"SCEC {PROB}{', ' + fname(fault) if fault else ''}: {N1} slip rate at the {len(sts)} "
                     f"on-fault stations (rows: x3; columns: x2 = {' / '.join(f'{x:g}' for x in xs)} km)", fontsize=12)
        fig.tight_layout(rect=(0, 0.05, 1, 0.96))
        save(fig, pdf, f"overview{'_' + fault if fault else ''}")

    for nm in order:
        st = next(s for s in stations if s["name"] == nm)
        tag, x, z = where(st)
        has_ref = (REF0, nm) in data
        ta = arrival(data[(REF0, nm)]) if has_ref else float("nan")
        zlab = REF0
        if not np.isfinite(ta):                 # no reference arrival: centre the zoom on the top run
            ta, zlab = arrival(data[(TOP_KEY, nm)]), TOP_LAB
            if not np.isfinite(ta):
                ta = 0.5
        fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
        draw(axs[0][0], nm, f"{C1}-slip-rate")
        axs[0][0].set_title(f"(a) {N1} slip rate (m/s)", fontsize=9)
        draw(axs[0][1], nm, f"{C1}-slip-rate", ta - 0.5, ta + 2.5)
        axs[0][1].set_title(f"(b) {N1} slip rate (m/s), {zlab} arrival {ta:.3f} s -0.5 / +2.5 s", fontsize=9)
        draw(axs[0][2], nm, f"{C1}-slip")
        axs[0][2].set_title(f"(c) {N1} slip (m)", fontsize=9)
        draw(axs[1][0], nm, f"{C1}-shear-stress")
        axs[1][0].set_title(f"(d) {N1} shear stress (MPa)", fontsize=9)
        draw(axs[1][1], nm, "n-stress")
        axs[1][1].set_title("(e) normal stress, compression positive (MPa), autoscaled", fontsize=9)
        draw(axs[1][2], nm, f"{C2}-slip-rate")
        axs[1][2].set_title(f"(f) {N2} slip rate (m/s), autoscaled", fontsize=9)
        for ax in list(axs[1]) + [axs[0][1]]:
            ax.set_xlabel("time (s)", fontsize=8)
        h, l = axs[0][0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=5, fontsize=8)
        v = mf.get((TOP_KEY, nm))
        if not has_ref:
            vs = f"no {REF0} file at this station"
        elif v == "quiet":
            vs = f"vs {REF0}: neither passes 1 mm/s"
        elif isinstance(v, tuple) and np.isfinite(v[0]):
            vs = f"vs {REF0}: arrival {v[0]:+.3f} s, final slip {v[1]:+.2f} %, peak slip rate {v[2]:+.1f} %"
        elif isinstance(v, tuple):
            who = REF0 if not np.isfinite(arrival(data[(REF0, nm)])) else "this run"
            vs = f"vs {REF0}: {who} does not pass 1 mm/s here"
        else:
            vs = f"vs {REF0}: n/a"
        fig.suptitle(f"SCEC {PROB} station {nm} ({fname(tag) + ', ' if tag else ''}x2 {x:g} km, x3 {z:g} km)\n"
                     f"{TOP_LAB} {vs}", fontsize=10)
        fig.tight_layout(rect=(0, 0.05, 1, 0.94))
        save(fig, pdf, nm)
print(f"wrote {pdf_path}: {pages} pages")
