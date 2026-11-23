#!/usr/bin/env python3
"""TPV34 station pages: Nitsche runs against the SCEC submissions.

The standard station-plot PDF (as the TPV31 pages), extended by the off-fault
stations and the rupture-front page:

  page 1      on-fault overview: arrival, final slip and peak slip rate of the
              last --run minus the primary reference, on the 5 x 7 station grid
  page 2      off-fault overview: peak h, v, n velocity of the last --run minus
              the primary reference (%), at the 28 map positions, 0 and 2.4 km
  page 3      rupture front: contours of the rupture time, last --run against
              the primary reference's cplot file (when both are there)
  pages 4..   one page per on-fault station (35): slip rate, slip rate zoom,
              slip, shear stress, normal stress, shear / normal
  then        one page per off-fault station (56): h, v, n velocity and
              displacement

References are drawn first (grey, the primary one black), the runs after
them, the last --run in red on top.  Velocity axes are linear.

Readers and the arrival / slip / peak measures are those of
scripts/compare_nitsche_fault_stations.py (imported).

Usage (from miniapps/seas):
  python3 tpv34/scripts/plot_tpv34_stations.py --out OUT_DIR \\
      --run "p1 200 m=RUN_DIR_P1" --run "p2 200 m=RUN_DIR_P2" \\
      [--data tpv34/benchmark_data] [--primary seissol_50m_o4] \\
      [--refs seissol_50m_o4,eqdyna_25m,...] [--t-end 20] [--pdf NAME.pdf]

RUN_DIR holds the driver's <prefix>_station_<name>.dat, <prefix>_offfault_
<name>.dat and <prefix>_rupture_time_r*.csv.  A reference <tag> is the folder
DATA/scec_<tag> with files tpv34_<tag>_<station>.txt and tpv34_<tag>_cplot.txt.
Writes OUT_DIR/pNNN_*.png, summary.json, compare_vs_scec.txt and the PDF.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SEAS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(SEAS, "scripts"))
import compare_nitsche_fault_stations as C  # noqa: E402

BODY_KEYS = ["t", "h-disp", "h-vel", "v-disp", "v-vel", "n-disp", "n-vel"]
RUN_COLORS = ["#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]   # earlier runs
LAST = dict(color="#d62728", lw=1.6, ls="-")                # the last run, on top
DPI = 110


def read_csv_stations(path):
    out = []
    for l in open(path):
        if not l.strip() or l.startswith("#") or l.startswith("name"):
            continue
        n, x, y, z = l.strip().split(",")
        out.append(dict(name=n, x=float(x), y=float(y), z=float(z)))
    return out


def one_file(pattern):
    hits = sorted(glob.glob(pattern))
    if len(hits) > 1:
        sys.exit(f"more than one file matches {pattern}: {hits}")
    return hits[0] if hits else None


def load_fault(d, name, is_run):
    if is_run:
        f = one_file(os.path.join(d, f"*_station_{name}.dat"))
        return C.load(f, "nitsche") if f else None
    f = one_file(os.path.join(d, f"*_{name}.txt"))
    return C.load(f, "scec_mpa") if f else None


def load_body(d, name, is_run):
    f = one_file(os.path.join(d, f"*_offfault_{name}.dat" if is_run else f"*_{name}.txt"))
    if not f:
        return None
    A = C._numeric_rows(f)
    if A.shape[1] < 7:
        sys.exit(f"{f}: {A.shape[1]} columns, the off-fault layout needs 7")
    return {k: A[:, i] for i, k in enumerate(BODY_KEYS)}


def rel(a, b):
    return (a - b) / b if b != 0.0 and np.isfinite(b) else float("nan")


def styles(refs, primary, runs):
    """(key, label, style) in drawing order: references, then runs."""
    out = []
    greys = ["0.45", "0.6", "0.7", "0.5", "0.65", "0.55", "0.75"]
    dashes = ["--", ":", "-.", (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (1, 1)), (0, (6, 2, 1, 2))]
    j = 0
    for tag in refs:
        if tag == primary:
            continue
        out.append(("ref:" + tag, tag, dict(color=greys[j % 7], lw=0.9, ls=dashes[j % 7])))
        j += 1
    out.append(("ref:" + primary, primary, dict(color="k", lw=1.1, ls="-")))
    for i, (lab, _) in enumerate(runs[:-1]):
        out.append(("run:" + lab, lab, dict(color=RUN_COLORS[i % 4], lw=1.2, ls="-")))
    out.append(("run:" + runs[-1][0], runs[-1][0], LAST))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--run", action="append", required=True, help="LABEL=DIR, repeatable; the last is drawn on top")
    ap.add_argument("--data", default=os.path.join(SEAS, "tpv34/benchmark_data"))
    ap.add_argument("--primary", default="seissol_50m_o4")
    ap.add_argument("--refs", default=None, help="comma list of tags; default: every DATA/scec_* folder")
    ap.add_argument("--t-end", type=float, default=20.0)
    ap.add_argument("--pdf", default="tpv34_station_plots.pdf")
    a = ap.parse_args()

    runs = []
    for r in a.run:
        if "=" not in r:
            sys.exit(f"--run {r!r}: expected LABEL=DIR")
        lab, d = r.split("=", 1)
        if not os.path.isdir(d):
            sys.exit(f"--run {lab}: no folder {d}")
        runs.append((lab, d))
    if a.refs:
        refs = a.refs.split(",")
    else:
        refs = sorted(os.path.basename(p)[5:] for p in glob.glob(os.path.join(a.data, "scec_*")) if os.path.isdir(p))
    if a.primary not in refs:
        sys.exit(f"primary reference {a.primary} is not among the references {refs}")
    for tag in refs:
        if not os.path.isdir(os.path.join(a.data, "scec_" + tag)):
            sys.exit(f"no folder {os.path.join(a.data, 'scec_' + tag)}")
    src_dir = {"ref:" + t: os.path.join(a.data, "scec_" + t) for t in refs}
    src_dir.update({"run:" + lab: d for lab, d in runs})
    sty = styles(refs, a.primary, runs)
    last, prim = "run:" + runs[-1][0], "ref:" + a.primary
    os.makedirs(a.out, exist_ok=True)
    fst = read_csv_stations(os.path.join(SEAS, "tpv34/stations/tpv34_fault_stations.csv"))
    bst = read_csv_stations(os.path.join(SEAS, "tpv34/stations/tpv34_body_stations.csv"))
    pages, text = {}, []
    missing = []

    # ------------------------------------------------------------ on-fault
    fsum = []
    for st in fst:
        s = {}
        for key, _, _ in sty:
            d = load_fault(src_dir[key], st["name"], key.startswith("run:"))
            if d is None:
                missing.append((key, st["name"]))
            else:
                s[key] = d
        if last not in s or prim not in s:
            sys.exit(f"station {st['name']}: the last run or the primary reference has no file")
        t_win = min([a.t_end] + [float(v["t"][-1]) for v in s.values()])
        m = {k: C.measures(v, t_win, 0.1, 0.5) for k, v in s.items()}
        xkm, zkm = st["x"] / 1e3, -st["z"] / 1e3 + 0.0
        row = dict(name=st["name"], strike_km=xkm, depth_km=zkm, t_win=t_win)
        for key, lab, _ in sty:
            if key in m and key != prim:
                mf = C.misfit(m[key], m[prim])
                row[lab] = dict(d_arrival=mf["d_arrival"], r_slip=mf["r_slip"], r_peak=mf["r_peak"],
                                arrival=m[key]["arrival"], slip=m[key]["slip"], peak=m[key]["peak"])
        row[a.primary] = dict(arrival=m[prim]["arrival"], slip=m[prim]["slip"], peak=m[prim]["peak"])
        fsum.append(row)
        mf = row[runs[-1][0]]
        fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
        fig.suptitle(f"SCEC TPV34 on-fault station {st['name']} (strike {xkm:g} km, depth {zkm:g} km)\n"
                     f"{runs[-1][0]} vs {a.primary}: arrival {mf['d_arrival']:+.3f} s, final slip "
                     f"{100 * mf['r_slip']:+.2f} %, peak slip rate {100 * mf['r_peak']:+.1f} %", fontsize=12)
        ta = m[prim]["arrival"]
        if not np.isfinite(ta):
            ta = m[last]["arrival"] if np.isfinite(m[last]["arrival"]) else 0.5
        panels = [(ax[0, 0], "h-slip-rate", "(a) strike slip rate (m/s)", None),
                  (ax[0, 1], "h-slip-rate", f"(b) strike slip rate (m/s), {a.primary} arrival {ta:.3f} s -0.5 / +2.5 s",
                   (ta - 0.5, ta + 2.5)),
                  (ax[0, 2], "h-slip", "(c) strike slip (m)", None),
                  (ax[1, 0], "h-shear-stress", "(d) strike shear stress (MPa)", None),
                  (ax[1, 1], "n-stress", "(e) normal stress, compression positive (MPa)", None),
                  (ax[1, 2], "ratio", "(f) shear / normal stress", None)]
        for axx, k, title, xl in panels:
            for key, lab, stl in sty:
                if key not in s:
                    continue
                d = s[key]
                if k != "ratio" and k not in d:
                    continue
                if k == "ratio" and "n-stress" not in d:
                    continue
                w = d["t"] <= t_win + 1e-9
                y = np.abs(d["h-shear-stress"]) / np.maximum(d["n-stress"], 1e-9) if k == "ratio" else d[k]
                axx.plot(d["t"][w], y[w], label=lab, **stl)
            axx.set_title(title, fontsize=10)
            axx.grid(alpha=0.3)
            axx.set_xlim(xl if xl else (0, t_win))
            if xl:
                top = 0.0
                for key in (last, prim):
                    sel = (s[key]["t"] >= xl[0]) & (s[key]["t"] <= xl[1])
                    if sel.any():
                        top = max(top, float(np.max(s[key]["h-slip-rate"][sel])))
                if top > 0:
                    axx.set_ylim(-0.05 * top, 1.08 * top)
            if axx in ax[1]:
                axx.set_xlabel("time (s)")
        sn = [float(v["n-stress"][0]) for v in s.values() if "n-stress" in v]
        lo, hi = min(sn), max(sn)
        ax[1, 1].set_ylim(max(0.0, lo - 0.15 * hi), hi + 0.15 * hi)
        ax[1, 2].set_ylim(0.3, 0.7)
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=min(6, len(l)), frameon=True, fontsize=9)
        fig.tight_layout(rect=(0, 0.06, 1, 0.93))
        pages[("f", st["name"])] = fig
        plt.close(fig)

    # ----------------------------------------------------------- off-fault
    bsum = []
    for st in bst:
        s = {}
        for key, _, _ in sty:
            d = load_body(src_dir[key], st["name"], key.startswith("run:"))
            if d is None:
                missing.append((key, st["name"]))
            else:
                s[key] = d
        if last not in s or prim not in s:
            sys.exit(f"station {st['name']}: the last run or the primary reference has no file")
        t_win = min([a.t_end] + [float(v["t"][-1]) for v in s.values()])
        row = dict(name=st["name"], strike_km=st["x"] / 1e3, normal_km=st["y"] / 1e3, depth_km=-st["z"] / 1e3 + 0.0)
        meas = {}
        for key, d in s.items():
            w = d["t"] <= t_win + 1e-9
            meas[key] = {c: dict(peak=float(np.max(np.abs(d[c + "-vel"][w]))),
                                 final=float(np.interp(t_win, d["t"], d[c + "-disp"]))) for c in "hvn"}
        for key, lab, _ in sty:
            if key in meas and key != prim:
                row[lab] = {c: dict(r_peak=rel(meas[key][c]["peak"], meas[prim][c]["peak"]),
                                    d_final=meas[key][c]["final"] - meas[prim][c]["final"],
                                    peak=meas[key][c]["peak"], final=meas[key][c]["final"]) for c in "hvn"}
        row[a.primary] = meas[prim]
        bsum.append(row)
        r = row[runs[-1][0]]
        fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
        fig.suptitle(f"SCEC TPV34 off-fault station {st['name']} ({row['normal_km']:g} km off fault, "
                     f"{row['strike_km']:g} km along strike, {row['depth_km']:g} km depth)\n"
                     f"{runs[-1][0]} vs {a.primary}, peak velocity: h {100 * r['h']['r_peak']:+.1f} %, "
                     f"v {100 * r['v']['r_peak']:+.1f} %, n {100 * r['n']['r_peak']:+.1f} %", fontsize=12)
        names = dict(h="horizontal, along strike", v="vertical, positive down", n="fault normal")
        for j, c in enumerate("hvn"):
            for i, (q, unit) in enumerate((("vel", "m/s"), ("disp", "m"))):
                axx = ax[i, j]
                for key, lab, stl in sty:
                    if key in s:
                        d = s[key]
                        w = d["t"] <= t_win + 1e-9
                        axx.plot(d["t"][w], d[f"{c}-{q}"][w], label=lab, **stl)
                axx.set_title(f"({'abcdef'[3 * i + j]}) {names[c]} "
                              f"{'velocity' if q == 'vel' else 'displacement'} ({unit})", fontsize=10)
                axx.grid(alpha=0.3)
                axx.set_xlim(0, t_win)
                if i == 1:
                    axx.set_xlabel("time (s)")
        h, l = ax[0, 0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=min(6, len(l)), frameon=True, fontsize=9)
        fig.tight_layout(rect=(0, 0.06, 1, 0.93))
        pages[("b", st["name"])] = fig
        plt.close(fig)

    # ------------------------------------------------- overview: on-fault
    lab_last = runs[-1][0]
    xs = sorted({r["strike_km"] for r in fsum})
    zs = sorted({r["depth_km"] for r in fsum})
    fig, ax = plt.subplots(1, 3, figsize=(15, 6.2))
    for axx, key, title, scale, lim in [(ax[0], "d_arrival", "arrival misfit (s)", 1.0, 0.35),
                                        (ax[1], "r_slip", "final slip misfit (%)", 100.0, 6.0),
                                        (ax[2], "r_peak", "peak slip rate misfit (%)", 100.0, 50.0)]:
        G = np.full((len(zs), len(xs)), np.nan)
        for r in fsum:
            G[zs.index(r["depth_km"]), xs.index(r["strike_km"])] = r[lab_last][key] * scale
        im = axx.imshow(G, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        for (zi, xi), v in np.ndenumerate(G):
            axx.text(xi, zi, f"{v:+.2f}" if key == "d_arrival" else f"{v:+.1f}", ha="center", va="center",
                     fontsize=8.5, color="white" if abs(v) > 0.6 * lim else "black")
        axx.set_xticks(range(len(xs)), [f"{x:g}" for x in xs])
        axx.set_yticks(range(len(zs)), [f"{z:g}" for z in zs])
        axx.set_xlabel("strike (km)")
        axx.set_ylabel("depth (km)")
        axx.set_title(title, fontsize=11)
        fig.colorbar(im, ax=axx, shrink=0.85)
    fig.suptitle(f"SCEC TPV34: {lab_last} minus {a.primary} at the {len(fsum)} on-fault stations "
                 f"(blue = early / low)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    pages["over_f"] = fig
    plt.close(fig)

    # ------------------------------------------------ overview: off-fault
    fig, ax = plt.subplots(2, 3, figsize=(15, 9.5))
    depths = sorted({r["depth_km"] for r in bsum})
    for i, dep in enumerate(depths[:2]):
        for j, c in enumerate("hvn"):
            axx = ax[i, j]
            rr = [r for r in bsum if r["depth_km"] == dep]
            v = np.array([100 * r[lab_last][c]["r_peak"] for r in rr])
            sc = axx.scatter([r["strike_km"] for r in rr], [r["normal_km"] for r in rr], c=v, cmap="RdBu_r",
                             vmin=-30, vmax=30, s=520, edgecolors="0.3", linewidths=0.5)
            for r, val in zip(rr, v):
                axx.text(r["strike_km"], r["normal_km"], f"{val:+.0f}", ha="center", va="center", fontsize=8,
                         color="white" if abs(val) > 18 else "black")
            axx.plot([-15, 15], [0, 0], color="k", lw=2)
            axx.set_xlim(-24, 24)
            axx.set_ylim(-18.5, 18.5)
            axx.set_title(f"peak {c} velocity misfit (%), depth {dep:g} km", fontsize=10)
            axx.set_xlabel("along strike (km)")
            axx.set_ylabel("fault normal (km)")
            axx.grid(alpha=0.3)
            fig.colorbar(sc, ax=axx, shrink=0.8)
    fig.suptitle(f"SCEC TPV34: {lab_last} minus {a.primary}, peak velocity at the {len(bsum)} off-fault "
                 f"stations (black line = fault trace)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pages["over_b"] = fig
    plt.close(fig)

    # ------------------------------------------------------ rupture front
    rt_files = sorted(glob.glob(os.path.join(runs[-1][1], "*_rupture_time_r*.csv")))
    cplot = one_file(os.path.join(src_dir[prim], "*_cplot.txt"))
    if rt_files and cplot:
        P = [C._numeric_rows(f, ncol_max=4) for f in rt_files
             if any(l.strip() and not l.startswith("#") for l in open(f))]   # ranks without fault: header only
        P = np.vstack(P)
        R = C._numeric_rows(cplot, ncol_max=3)      # j (along strike, m), k (down dip, m), t (s)
        fig, axx = plt.subplots(1, 1, figsize=(15, 8))
        levels = np.arange(0.5, 12.01, 0.5)
        ok = P[:, 3] >= 0
        if ok.sum() > 3:
            axx.tricontour(P[ok, 0] / 1e3, -P[ok, 2] / 1e3, P[ok, 3], levels=levels, colors=LAST["color"],
                           linewidths=1.4)
        okr = R[:, 2] < 1e8
        cs = axx.tricontour(R[okr, 0] / 1e3, R[okr, 1] / 1e3, R[okr, 2], levels=levels, colors="k",
                            linewidths=0.9)
        axx.clabel(cs, levels[1::2], fontsize=8, fmt="%.0f s")
        axx.plot([], [], color=LAST["color"], lw=1.4, label=lab_last)
        axx.plot([], [], color="k", lw=0.9, label=a.primary)
        axx.set_xlim(-15, 15)
        axx.set_ylim(15, 0)
        axx.set_aspect("equal")
        axx.set_xlabel("along strike (km)")
        axx.set_ylabel("depth (km)")
        axx.legend(loc="lower right")
        axx.set_title(f"SCEC TPV34 rupture front every 0.5 s (first slip rate above 1 mm/s): {lab_last} "
                      f"and {a.primary}", fontsize=12)
        fig.tight_layout()
        pages["front"] = fig
        plt.close(fig)
    else:
        text.append(f"rupture-front page left out: rupture_time files {len(rt_files)}, cplot {cplot}")

    # --------------------------------------------------------------- write
    order = ["over_f", "over_b"] + (["front"] if "front" in pages else []) \
        + [("f", s["name"]) for s in fst] + [("b", s["name"]) for s in bst]
    pngs = []
    for i, k in enumerate(order):
        tag = k if isinstance(k, str) else k[1]
        fn = os.path.join(a.out, f"p{i + 1:03d}_{tag}.png")
        pages[k].savefig(fn, dpi=DPI)
        pngs.append(fn)
    from PIL import Image
    ims = [Image.open(f).convert("RGB") for f in pngs]
    ims[0].save(os.path.join(a.out, a.pdf), save_all=True, append_images=ims[1:], resolution=DPI)

    with open(os.path.join(a.out, "compare_vs_scec.txt"), "w") as f:
        f.write(f"# TPV34: misfits against {a.primary}; t_end {a.t_end} s (lowered per station to the earliest end)\n")
        f.write("# on-fault: station source arrival(s) d_arrival(s) slip(m) r_slip(%) peak(m/s) r_peak(%)\n")
        for r in fsum:
            p = r[a.primary]
            f.write(f"{r['name']:18s} {a.primary:22s} {p['arrival']:8.3f} {'':>8s} {p['slip']:8.4f} {'':>7s} "
                    f"{p['peak']:8.4f}\n")
            for key, lab, _ in sty:
                if lab in r and key != prim:
                    q = r[lab]
                    f.write(f"{r['name']:18s} {lab:22s} {q['arrival']:8.3f} {q['d_arrival']:+8.3f} {q['slip']:8.4f} "
                            f"{100 * q['r_slip']:+7.2f} {q['peak']:8.4f} {100 * q['r_peak']:+7.1f}\n")
        f.write("# off-fault: station source  peak |vel| h v n (m/s)  r_peak h v n (%)  final disp h v n (m)\n")
        for r in bsum:
            for key, lab, _ in sty:
                if lab in r and key != prim:
                    q = r[lab]
                    f.write(f"{r['name']:22s} {lab:22s} " + " ".join(f"{q[c]['peak']:8.4f}" for c in "hvn") + "  "
                            + " ".join(f"{100 * q[c]['r_peak']:+7.1f}" for c in "hvn") + "  "
                            + " ".join(f"{q[c]['final']:+8.4f}" for c in "hvn") + "\n")
        med = lambda v: float(np.nanmedian(np.abs(v)))  # noqa: E731
        for lab, _ in runs:
            f.write(f"# median |misfit| {lab}: arrival {med([r[lab]['d_arrival'] for r in fsum]):.3f} s, slip "
                    f"{100 * med([r[lab]['r_slip'] for r in fsum]):.2f} %, peak slip rate "
                    f"{100 * med([r[lab]['r_peak'] for r in fsum]):.1f} %; off-fault peak velocity h "
                    f"{100 * med([r[lab]['h']['r_peak'] for r in bsum]):.1f} %, v "
                    f"{100 * med([r[lab]['v']['r_peak'] for r in bsum]):.1f} %, n "
                    f"{100 * med([r[lab]['n']['r_peak'] for r in bsum]):.1f} %\n")
        for t in text:
            f.write("# " + t + "\n")
        for key, name in missing:
            f.write(f"# no file: {key} {name}\n")
    json.dump(dict(primary=a.primary, refs=refs, runs=[l for l, _ in runs], t_end=a.t_end, fault=fsum, body=bsum),
              open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    print(f"pages {len(pngs)}; missing files {len(missing)}; {os.path.join(a.out, a.pdf)}")
    print(open(os.path.join(a.out, "compare_vs_scec.txt")).read().split("# median", 1)[-1] if runs else "")


if __name__ == "__main__":
    main()
