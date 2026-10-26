#!/usr/bin/env python3
"""TPV205-2D station pages: Nitsche runs against a SCEC 2-D reference (PyLith).

The style of the TPV31 station plots (benchmarks thread, plot_tpv31_stations.py):
p01 = overview of the misfits per station, then one 2 x 3 page per station
(slip rate, slip rate zoomed on the reference arrival, slip, shear stress,
normal stress, shear / normal), PNG at 110 dpi and one PDF of all pages.
The reference is black; the runs follow in the order given and the LAST run
is drawn last, on top, in red.  Measures and the Nitsche reader are those of
scripts/compare_nitsche_fault_stations.py (imported): arrival = first
|V| > 1 mm/s, peak = max |V|, slip = |slip| at the end of the common window.

Usage:
  plot_tpv205_2d_stations.py SEAS_DIR REF_LABEL=REF_DIR OUT_DIR T_END LABEL=RUN_DIR [LABEL=RUN_DIR ...]
SEAS_DIR = miniapps/seas of the checkout (stations from
tpv205/stations/tpv205_2d_fault_stations.csv).  A reference file is found by
the station name or key (as the compare script does) or, failing that, by
the SCEC name faultst<strike in 0.1 km, signed 3 digits>dp<depth in 0.1 km,
3 digits> (e.g. faultst-045dp075); its columns are read from its header
line when it has one (t h-slip h-slip-rate h-shear-stress ... n-stress),
else positionally; stresses in MPa.
"""
import json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

seas = sys.argv[1]
ref_lab, ref_dir = sys.argv[2].split("=", 1)
out, t_end = sys.argv[3], float(sys.argv[4])
runs = [a.split("=", 1) for a in sys.argv[5:]]
if not runs:
    sys.exit("give at least one LABEL=RUN_DIR")
sys.path.insert(0, os.path.join(seas, "scripts"))
import compare_nitsche_fault_stations as C

os.makedirs(out, exist_ok=True)
stations = C.read_stations(os.path.join(seas, "tpv205/stations/tpv205_2d_fault_stations.csv"))
KEYS = ["t", "h-slip", "h-slip-rate", "h-shear-stress", "v-slip", "v-slip-rate",
        "v-shear-stress", "n-stress"]


def ref_file(st):
    f = C.find_file(ref_dir, st)
    if f:
        return f
    m = re.match(r"x2_(-?[\d.]+)_x3_([\d.]+)$", st["name"])
    s_km, d_km = float(m.group(1)), float(m.group(2))
    # PyLith 2-D files are named tpv205_pylith2d_x2_-12.0_x3_7.5.txt (one
    # decimal), which the exact-key lookup misses: match x2 and x3 by value,
    # on-fault files only (not the off-fault _body_ files).
    num = [os.path.join(ref_dir, fn) for fn in sorted(os.listdir(ref_dir))
           if "_body_" not in fn and (mm := re.search(r"_x2_(-?[\d.]+)_x3_(-?[\d.]+)\.txt$", fn))
           and abs(float(mm.group(1)) - s_km) < 1e-9 and abs(float(mm.group(2)) - d_km) < 1e-9]
    if len(num) == 1:
        return num[0]
    pat = re.compile(r"faultst\s*([+-]?)0*(\d+)dp0*(\d+)", re.I)
    hits = []
    for fn in sorted(os.listdir(ref_dir)):
        mm = pat.search(fn)
        if not mm:
            continue
        s = (-1 if mm.group(1) == "-" else 1) * int(mm.group(2)) / 10.0
        if abs(s - s_km) < 1e-9 and abs(int(mm.group(3)) / 10.0 - d_km) < 1e-9:
            hits.append(os.path.join(ref_dir, fn))
    if len(hits) != 1:
        sys.exit(f"{ref_dir}: {len(hits)} files for station {st['name']} "
                 f"(files: {', '.join(sorted(os.listdir(ref_dir))[:20])})")
    return hits[0]


def load_ref(path):
    """SCEC text file, stresses in MPa -> dict like C.load."""
    names, rows = None, []
    for line in open(path):
        w = line.split()
        if not w:
            continue
        if line.lstrip().startswith("#"):
            continue
        try:
            rows.append([float(v) for v in w])
        except ValueError:
            if "h-slip" in w or "h-slip-rate" in w:
                names = w
    if not rows:
        sys.exit(f"{path}: no numeric rows")
    n = min(len(r) for r in rows)
    A = np.array([r[:n] for r in rows])
    cols = names[:n] if names else KEYS[:n]
    d = {}
    for i, k in enumerate(cols):
        k = "t" if k in ("t", "time", "Time") else k
        d[k] = A[:, i]
    for k in ("h-slip", "h-slip-rate", "h-shear-stress"):
        if k not in d:
            sys.exit(f"{path}: no {k} column (columns {cols})")
    for k in ("v-slip", "v-slip-rate"):
        d.setdefault(k, np.zeros_like(d["t"]))
    if "n-stress" in d:
        d["n-stress"] = np.abs(d["n-stress"])
    d["_cols"] = cols
    return d


def ratio(d):
    return np.abs(d["h-shear-stress"]) / np.maximum(d["n-stress"], 1e-9)


srcs, ref_paths = {}, {}
for st in stations:
    s = {}
    for lab, rd in runs:
        f = C.find_file(rd, st)
        if not f:
            sys.exit(f"{rd}: no station file for {st['name']}")
        s[lab] = C.load(f, "nitsche")
    ref_paths[st["name"]] = ref_file(st)
    s["_ref"] = load_ref(ref_paths[st["name"]])
    srcs[st["name"]] = s
t_win = min([t_end] + [float(v["t"][-1]) for s in srcs.values() for v in s.values()])

COL = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
sty = {}
for i, (lab, _) in enumerate(runs):
    last = i == len(runs) - 1
    sty[lab] = dict(color="#d62728", lw=1.6) if last else dict(color=COL[i % len(COL)], lw=1.1, alpha=0.85)
REF = dict(color="k", lw=1.0)
last = runs[-1][0]

summary, table = [], []
for i, st in enumerate(stations):
    s = srcs[st["name"]]
    m = {k: C.measures(v, t_win, 0.1, 0.5) for k, v in s.items()}
    mf = {lab: C.misfit(m[lab], m["_ref"]) for lab, _ in runs}
    xkm = st["x"] / 1e3 + 0.0
    row = dict(name=st["name"], strike_km=xkm, ref_file=os.path.basename(ref_paths[st["name"]]),
               ref=dict(arrival=m["_ref"]["arrival"], peak=m["_ref"]["peak"], slip=m["_ref"]["slip"]))
    for lab, _ in runs:
        row[lab] = dict(arrival=m[lab]["arrival"], peak=m[lab]["peak"], slip=m[lab]["slip"],
                        d_arr=mf[lab]["d_arrival"], r_slip=mf[lab]["r_slip"], r_peak=mf[lab]["r_peak"])
    summary.append(row)
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    fig.suptitle(f"SCEC TPV205-2D station {st['name']} (strike {xkm:g} km)\n"
                 f"{last} vs {ref_lab}: arrival {mf[last]['d_arrival']:+.3f} s, final slip "
                 f"{100 * mf[last]['r_slip']:+.2f} %, peak slip rate {100 * mf[last]['r_peak']:+.1f} %",
                 fontsize=12)
    ta = m["_ref"]["arrival"]
    if not np.isfinite(ta):
        ta = m[last]["arrival"] if np.isfinite(m[last]["arrival"]) else 0.5
    panels = [(ax[0, 0], "h-slip-rate", "(a) slip rate (m/s)", None),
              (ax[0, 1], "h-slip-rate", f"(b) slip rate (m/s), {ref_lab} arrival {ta:.3f} s -0.5 / +2.5 s", (ta - 0.5, ta + 2.5)),
              (ax[0, 2], "h-slip", "(c) slip (m)", None),
              (ax[1, 0], "h-shear-stress", "(d) shear stress (MPa)", None),
              (ax[1, 1], "n-stress", "(e) normal stress, compression positive (MPa)", None),
              (ax[1, 2], "ratio", "(f) shear / normal stress", None)]
    for a, k, title, xl in panels:
        for key, lab, sk in [("_ref", ref_lab, REF)] + [(l, l, sty[l]) for l, _ in runs]:
            d = s[key]
            if k in ("n-stress", "ratio") and "n-stress" not in d:
                continue
            w = d["t"] <= t_win + 1e-9
            y = ratio(d) if k == "ratio" else d[k]
            a.plot(d["t"][w], y[w], label=lab, **sk)
        a.set_title(title, fontsize=10); a.grid(alpha=0.3)
        a.set_xlim(xl if xl else (0, t_win))
        if xl:
            tops = [np.max(s[key]["h-slip-rate"][(s[key]["t"] >= xl[0]) & (s[key]["t"] <= xl[1])], initial=0.0)
                    for key in ["_ref"] + [l for l, _ in runs]]
            if max(tops) > 0:
                a.set_ylim(-0.05 * max(tops), 1.08 * max(tops))
        if a in ax[1]:
            a.set_xlabel("time (s)")
    sn = [float(v["n-stress"][0]) for v in s.values() if "n-stress" in v]
    if sn:
        lo, hi = min(sn), max(sn)
        ax[1, 1].set_ylim(lo - 0.05 * hi, hi + 0.05 * hi)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=len(l), frameon=True, fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    fn = f"p{i + 2:02d}_{st['name'].replace('.', 'p')}.png"
    fig.savefig(os.path.join(out, fn), dpi=110); plt.close(fig)
    row["png"] = fn

# overview: misfits per station against the reference, one marker series per run
xs = [r["strike_km"] for r in summary]
fig, ax = plt.subplots(1, 3, figsize=(15, 5.6))
for a, key, title, scale in [(ax[0], "d_arr", "arrival misfit (s)", 1.0),
                             (ax[1], "r_slip", "final slip misfit (%)", 100.0),
                             (ax[2], "r_peak", "peak slip rate misfit (%)", 100.0)]:
    for lab, _ in runs:
        sk = sty[lab]
        a.plot(xs, [r[lab][key] * scale for r in summary], "o-", color=sk["color"], lw=sk["lw"],
               label=lab, zorder=3 if lab == last else 2)
    a.axhline(0.0, color="k", lw=0.8)
    for x, r in zip(xs, summary):
        v = r[last][key] * scale
        a.annotate(f"{v:+.3f}" if key == "d_arr" else (f"{v:+.2f}" if key == "r_slip" else f"{v:+.1f}"), (x, v),
                   textcoords="offset points",
                   xytext=(0, 7), ha="center", fontsize=8, color="#d62728")
    a.set_xticks(xs, [f"{x:g}" for x in xs]); a.set_xlabel("station strike (km)")
    a.set_title(title, fontsize=11); a.grid(alpha=0.3)
ax[0].legend(fontsize=9)
fig.suptitle(f"SCEC TPV205-2D: runs minus {ref_lab} at the 7 on-fault stations, window 0 to {t_win:g} s "
             f"(labels: {last})", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(os.path.join(out, "p01_overview.png"), dpi=110); plt.close(fig)

# text table and summary
with open(os.path.join(out, "compare_vs_ref.txt"), "w") as f:
    f.write(f"TPV205-2D, window 0 to {t_win:g} s; arrival = first |V| > 1 mm/s; slip at the window end\n")
    f.write(f"{'station':16s} {'source':14s} {'arrival s':>10s} {'peak m/s':>9s} {'slip m':>8s}"
            f" {'d_arr s':>8s} {'d_slip %':>9s} {'d_peak %':>9s}\n")
    for r in summary:
        f.write(f"{r['name']:16s} {ref_lab:14s} {r['ref']['arrival']:10.4f} {r['ref']['peak']:9.4f} "
                f"{r['ref']['slip']:8.4f}   ({r['ref_file']})\n")
        for lab, _ in runs:
            q = r[lab]
            f.write(f"{'':16s} {lab:14s} {q['arrival']:10.4f} {q['peak']:9.4f} {q['slip']:8.4f} "
                    f"{q['d_arr']:+8.4f} {100 * q['r_slip']:+9.2f} {100 * q['r_peak']:+9.1f}\n")
json.dump(dict(reference=ref_lab, runs=[l for l, _ in runs], t_window=t_win, stations=summary),
          open(os.path.join(out, "summary.json"), "w"), indent=1, default=float)

# one PDF of all pages (Pillow, 110 dpi)
from PIL import Image
pages = ["p01_overview.png"] + [r["png"] for r in summary]
ims = [Image.open(os.path.join(out, p)).convert("RGB") for p in pages]
ims[0].save(os.path.join(out, "tpv205_2d_stations.pdf"), save_all=True, append_images=ims[1:], resolution=110)
print("t_window", t_win, "pages", len(pages), "->", os.path.join(out, "tpv205_2d_stations.pdf"))
print(open(os.path.join(out, "compare_vs_ref.txt")).read())
