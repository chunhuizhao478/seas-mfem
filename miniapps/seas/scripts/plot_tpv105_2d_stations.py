#!/usr/bin/env python3
"""TPV105-2D station pages: Nitsche runs against SCEC 2-D references (MDSBI, SPECFEM2D).

The style of the TPV31 / TPV205-2D station plots: p01 = overview of the misfits
per station, then one page per station, PNG at 110 dpi and one PDF of all pages.
TPV105 adds thermal pressurization, so a station page is 3 x 3: slip rate, slip
rate zoomed on the first reference's arrival, slip; shear stress, effective
normal stress, shear / normal; temperature, pore pressure, state psi.
References are black (first), grey (second); the runs follow in the order given
and the LAST run is drawn last, on top, in red.  Misfits are against the first
reference; the second reference's own misfit against the first is drawn in grey
to show the spread between the reference codes.  Measures are those of
scripts/compare_nitsche_fault_stations.py (imported): arrival = first
|V| > 1 mm/s, peak = max |V|, slip = |slip| at the end of the common window;
plus peak temperature and pore pressure at the window end.

Usage:
  plot_tpv105_2d_stations.py SEAS_DIR OUT_DIR T_END REF1=DIR REF2=DIR -- LABEL=RUN_DIR [LABEL=RUN_DIR ...]
A reference file is found by the station name or key, or by the SCEC name
faultst<strike in 0.1 km, signed>dp<depth in 0.1 km>; its columns come from its
header line (t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate
v-shear-stress n-stress psi temperature pressure); stresses and pressure in MPa,
temperature in K.  A run file is the Nitsche station layout
t slip1 slip2 V1 V2 tau1 tau2 sigma_n psi T p (SI; sigma_n = effective normal
stress with thermal pressurization).
"""
import json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

a = sys.argv[1:]
if "--" not in a:
    sys.exit(__doc__)
k = a.index("--")
seas, out, t_end = a[0], a[1], float(a[2])
refs = [s.split("=", 1) for s in a[3:k]]
runs = [s.split("=", 1) for s in a[k + 1:]]
if not refs or not runs:
    sys.exit("give at least one REF=DIR and one LABEL=RUN_DIR")
sys.path.insert(0, os.path.join(seas, "scripts"))
import compare_nitsche_fault_stations as C

os.makedirs(out, exist_ok=True)
stations = C.read_stations(os.path.join(seas, "tpv105/stations/tpv105_2d_fault_stations.csv"))


def num(s):
    try:
        return float(s)
    except ValueError:
        return None


def ref_file(rd, st):
    f = C.find_file(rd, st)
    if f:
        return f
    m = re.match(r"x2_(-?[\d.]+)_x3_([\d.]+)$", st["name"])
    s_km, d_km = float(m.group(1)), float(m.group(2))
    hits = []
    for fn in sorted(os.listdir(rd)):
        mm = re.search(r"faultst\s*([+-]?)0*(\d+)dp0*(\d+)", fn, re.I)
        if mm:
            s = (-1 if mm.group(1) == "-" else 1) * int(mm.group(2)) / 10.0
            if abs(s - s_km) < 1e-9 and abs(int(mm.group(3)) / 10.0 - d_km) < 1e-9:
                hits.append(os.path.join(rd, fn))
            continue
        mm = re.search(r"_x2_(-?[\d.]+)_x3_(-?[\d.]+)\.txt$", fn)
        if mm and "_body_" not in fn and abs(float(mm.group(1)) - s_km) < 1e-9 \
                and abs(float(mm.group(2)) - d_km) < 1e-9:
            hits.append(os.path.join(rd, fn))
    if len(hits) != 1:
        sys.exit(f"{rd}: {len(hits)} files for station {st['name']} "
                 f"(files: {', '.join(sorted(os.listdir(rd))[:24])})")
    return hits[0]


def load_ref(path):
    names, rows = None, []
    for line in open(path):
        w = line.split()
        if not w or line.lstrip().startswith("#"):
            continue
        v = [num(x) for x in w]
        if None in v:
            if "h-slip" in w or "h-slip-rate" in w:
                names = w
            continue
        rows.append(v)
    if not rows or not names:
        sys.exit(f"{path}: {'no numeric rows' if not rows else 'no header line with h-slip'}")
    n = min(len(r) for r in rows)
    A = np.array([r[:n] for r in rows])
    d = {}
    for i, k in enumerate(names[:n]):
        d["t" if k in ("t", "time", "Time") else k] = A[:, i]
    for k in ("h-slip", "h-slip-rate", "h-shear-stress"):
        if k not in d:
            sys.exit(f"{path}: no {k} column (columns {names[:n]})")
    for k in ("v-slip", "v-slip-rate"):
        d.setdefault(k, np.zeros_like(d["t"]))
    if "n-stress" in d:
        d["n-stress"] = np.abs(d["n-stress"])
    d["_cols"] = names[:n]
    return d


def load_run(path):
    rows = []
    for line in open(path):
        w = line.split()
        if not w or line.lstrip().startswith("#"):
            continue
        v = [num(x) for x in w]
        if None not in v:
            rows.append(v)
    n = min(len(r) for r in rows)
    if n < 11:
        sys.exit(f"{path}: {n} columns, the TPV105 run layout needs 11")
    A = np.array([r[:n] for r in rows])
    return {"t": A[:, 0], "v-slip": A[:, 1], "h-slip": A[:, 2], "v-slip-rate": A[:, 3],
            "h-slip-rate": A[:, 4], "v-shear-stress": A[:, 5] / 1e6, "h-shear-stress": A[:, 6] / 1e6,
            "n-stress": np.abs(A[:, 7]) / 1e6, "psi": A[:, 8], "temperature": A[:, 9],
            "pressure": A[:, 10] / 1e6}


def ratio(d):
    return np.abs(d["h-shear-stress"]) / np.maximum(d["n-stress"], 1e-9)


def tp_measures(d, t_win):
    m = d["t"] <= t_win + 1e-9
    o = {}
    o["T_peak"] = float(np.max(d["temperature"][m])) if "temperature" in d else float("nan")
    o["p_end"] = float(np.interp(t_win, d["t"], d["pressure"])) if "pressure" in d else float("nan")
    return o


srcs, ref_paths = {}, {}
for st in stations:
    s = {}
    for lab, rd in runs:
        f = C.find_file(rd, st)
        if not f:
            sys.exit(f"{rd}: no station file for {st['name']}")
        s[lab] = load_run(f)
    for lab, rd in refs:
        p = ref_file(rd, st)
        ref_paths[(lab, st["name"])] = p
        s[lab] = load_ref(p)
    srcs[st["name"]] = s
t_win = min([t_end] + [float(v["t"][-1]) for s in srcs.values() for v in s.values()])

COL = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
RCOL = [dict(color="k", lw=1.0), dict(color="#8c8c8c", lw=1.0)]
sty = {}
for i, (lab, _) in enumerate(refs):
    sty[lab] = RCOL[min(i, 1)]
for i, (lab, _) in enumerate(runs):
    sty[lab] = dict(color="#d62728", lw=1.6) if i == len(runs) - 1 else dict(color=COL[i % len(COL)], lw=1.1, alpha=0.85)
r0, last = refs[0][0], runs[-1][0]
order = [l for l, _ in refs] + [l for l, _ in runs]
cmp_labs = [l for l, _ in refs[1:]] + [l for l, _ in runs]

summary = []
for i, st in enumerate(stations):
    s = srcs[st["name"]]
    m = {k: dict(C.measures(v, t_win, 0.1, 0.5), **tp_measures(v, t_win)) for k, v in s.items()}
    xkm = st["x"] / 1e3 + 0.0
    row = dict(name=st["name"], strike_km=xkm,
               ref_files={l: os.path.basename(ref_paths[(l, st["name"])]) for l, _ in refs})
    for lab in order:
        q = dict(arrival=m[lab]["arrival"], peak=m[lab]["peak"], slip=m[lab]["slip"],
                 T_peak=m[lab]["T_peak"], p_end=m[lab]["p_end"])
        if lab != r0:
            mf = C.misfit(m[lab], m[r0])
            q.update(d_arr=mf["d_arrival"], r_slip=mf["r_slip"], r_peak=mf["r_peak"],
                     r_T=(m[lab]["T_peak"] - m[r0]["T_peak"]) / m[r0]["T_peak"],
                     d_p=m[lab]["p_end"] - m[r0]["p_end"])
        row[lab] = q
    summary.append(row)
    L = row[last]
    fig, ax = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle(f"SCEC TPV105-2D station {st['name']} (strike {xkm:g} km)\n"
                 f"{last} vs {r0}: arrival {L['d_arr']:+.3f} s, final slip {100 * L['r_slip']:+.2f} %, "
                 f"peak slip rate {100 * L['r_peak']:+.1f} %, peak T {100 * L['r_T']:+.1f} %, "
                 f"final p {L['d_p']:+.2f} MPa", fontsize=12)
    ta = m[r0]["arrival"]
    if not np.isfinite(ta):
        ta = m[last]["arrival"] if np.isfinite(m[last]["arrival"]) else 0.5
    panels = [(ax[0, 0], "h-slip-rate", "(a) slip rate (m/s)", None),
              (ax[0, 1], "h-slip-rate", f"(b) slip rate (m/s), {r0} arrival {ta:.3f} s -0.5 / +2.5 s", (ta - 0.5, ta + 2.5)),
              (ax[0, 2], "h-slip", "(c) slip (m)", None),
              (ax[1, 0], "h-shear-stress", "(d) shear stress (MPa)", None),
              (ax[1, 1], "n-stress", "(e) effective normal stress, compression positive (MPa)", None),
              (ax[1, 2], "ratio", "(f) shear / effective normal stress", None),
              (ax[2, 0], "temperature", "(g) temperature (K)", None),
              (ax[2, 1], "pressure", "(h) pore pressure (MPa)", None),
              (ax[2, 2], "psi", "(i) state variable psi", None)]
    for a_, k, title, xl in panels:
        for lab in order:
            d = s[lab]
            need = "n-stress" if k == "ratio" else k
            if need not in d:
                continue
            w = d["t"] <= t_win + 1e-9
            y = ratio(d) if k == "ratio" else d[k]
            a_.plot(d["t"][w], y[w], label=lab, **sty[lab])
        a_.set_title(title, fontsize=10); a_.grid(alpha=0.3)
        a_.set_xlim(xl if xl else (0, t_win))
        if xl:
            tops = [np.max(s[l]["h-slip-rate"][(s[l]["t"] >= xl[0]) & (s[l]["t"] <= xl[1])], initial=0.0) for l in order]
            if max(tops) > 0:
                a_.set_ylim(-0.05 * max(tops), 1.08 * max(tops))
        if a_ in ax[2]:
            a_.set_xlabel("time (s)")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=len(l), frameon=True, fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    fn = f"p{i + 2:02d}_{st['name'].replace('.', 'p')}.png"
    fig.savefig(os.path.join(out, fn), dpi=110); plt.close(fig)
    row["png"] = fn

# overview: misfits per station against the first reference
xs = [r["strike_km"] for r in summary]
fig, ax = plt.subplots(2, 3, figsize=(15, 9))
specs = [(ax[0, 0], "d_arr", "arrival misfit (s)", 1.0, "{:+.3f}"),
         (ax[0, 1], "r_slip", "final slip misfit (%)", 100.0, "{:+.2f}"),
         (ax[0, 2], "r_peak", "peak slip rate misfit (%)", 100.0, "{:+.1f}"),
         (ax[1, 0], "r_T", "peak temperature misfit (%)", 100.0, "{:+.1f}"),
         (ax[1, 1], "d_p", "final pore pressure difference (MPa)", 1.0, "{:+.2f}")]
for a_, key, title, scale, f in specs:
    for lab in cmp_labs:
        sk = sty[lab]
        a_.plot(xs, [r[lab][key] * scale for r in summary], "o-", color=sk["color"], lw=sk["lw"],
                label=f"{lab} - {r0}", zorder=3 if lab == last else 2)
    a_.axhline(0.0, color="k", lw=0.8)
    for x, r in zip(xs, summary):
        v = r[last][key] * scale
        a_.annotate(f.format(v), (x, v), textcoords="offset points", xytext=(0, 7), ha="center",
                    fontsize=8, color="#d62728")
    a_.set_xticks(xs, [f"{x:g}" for x in xs]); a_.set_xlabel("station strike (km)")
    a_.set_title(title, fontsize=11); a_.grid(alpha=0.3)
ax[0, 0].legend(fontsize=9)
ax[1, 2].axis("off")
lines = [f"window 0 to {t_win:g} s (the shortest record)", f"misfits against {r0}", ""]
lines += [l for l in os.environ.get("PLOT_NOTE", "").split("|") if l]
ax[1, 2].text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=9)
fig.suptitle(f"SCEC TPV105-2D: runs and references minus {r0} at the 7 on-fault stations "
             f"(labels: {last})", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(out, "p01_overview.png"), dpi=110); plt.close(fig)

with open(os.path.join(out, "compare_vs_ref.txt"), "w") as fo:
    fo.write(f"TPV105-2D, window 0 to {t_win:g} s; arrival = first |V| > 1 mm/s; slip and p at the window end; "
             f"misfits against {r0}\n")
    fo.write(f"{'station':16s} {'source':14s} {'arrival s':>10s} {'peak m/s':>9s} {'slip m':>8s} {'T_peak K':>9s} "
             f"{'p_end MPa':>10s} {'d_arr s':>8s} {'d_slip %':>9s} {'d_peak %':>9s} {'d_T %':>7s} {'d_p MPa':>8s}\n")
    for r in summary:
        for lab in order:
            q = r[lab]
            line = (f"{r['name'] if lab == order[0] else '':16s} {lab:14s} {q['arrival']:10.4f} {q['peak']:9.4f} "
                    f"{q['slip']:8.4f} {q['T_peak']:9.1f} {q['p_end']:10.2f}")
            if lab != r0:
                line += (f" {q['d_arr']:+8.4f} {100 * q['r_slip']:+9.2f} {100 * q['r_peak']:+9.1f} "
                         f"{100 * q['r_T']:+7.1f} {q['d_p']:+8.2f}")
            else:
                line += f"   ({r['ref_files'][lab]})"
            fo.write(line + "\n")
json.dump(dict(references=[l for l, _ in refs], runs=[l for l, _ in runs], t_window=t_win, stations=summary),
          open(os.path.join(out, "summary.json"), "w"), indent=1, default=float)

from PIL import Image
pages = ["p01_overview.png"] + [r["png"] for r in summary]
ims = [Image.open(os.path.join(out, p)).convert("RGB") for p in pages]
ims[0].save(os.path.join(out, "tpv105_2d_stations.pdf"), save_all=True, append_images=ims[1:], resolution=110)
print("t_window", t_win, "pages", len(pages), "->", os.path.join(out, "tpv105_2d_stations.pdf"))
print(open(os.path.join(out, "compare_vs_ref.txt")).read())
