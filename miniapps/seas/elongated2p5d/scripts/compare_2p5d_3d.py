#!/usr/bin/env python3
"""Compare map-view 2.5D runs with the 3-D buried reference (Level 2).

The 2.5D note's Level 2 test (Chunhui Zhao, 2026-10-30;
/mnt/project-files/2d-inplane/nonplanar-literature/
inplane-2p5d-derivation-rev8.pdf, Section 7.2): an elongated vertical
strike-slip rupture, L = 120 km by W = 15 km, run in 3-D
(elongated2p5d/configs/elongated_3d_buried.toml) and in map view with the
2.5D term (elongated_2d_2p5d_*.toml) and without it
(elongated_2d_plane_strain.toml).  The map-view slip is the depth-averaged
slip of the 3-D fault, so every 3-D quantity below is averaged over the
width first.

Sources (Nitsche driver output folders):
  --run3d DIR               the 3-D run: <prefix>_station_x2_<x>_z_<z>.dat
                            (elongated_3d_stations.csv: x every 2.5 km, z =
                            -7 ... 7 km), <prefix>_rupture_time_r*.csv and,
                            when the run had --dump-every, <prefix>_fault_r*.csv
  --run2d LABEL=DIR[:C]     a map-view run (repeatable, order kept):
                            <prefix>_station_x2_<x>.dat and the rupture times;
                            C = its c (1.178 wa2020, 1.273 lehner, 1.571 Liang)
                            for the implied-c estimate, omitted for the
                            plane-strain control
The prefix is read from the folder's *_station_conventions.txt.

Measures (all on the along-strike window |x| in [--x-min, --x-max], 15 to
57.5 km by default: past x = W the rupture length is > 2 W, the note's
"once L/W > 2"):
  final slip   at t_end (the earliest end of all runs): the 2-D station's
               slip2, against the 3-D depth average at the same x
               (trapezoid over the 15 depths with slip 0 at the fault's
               edges z = +-7.5 km; with a sqrt edge profile this reads up to
               about 1% low, so the dump average below is also printed
               when the 3-D run has dumps)
  arrival      first time the depth-averaged slip rate (2-D: V2) exceeds
               1 mm/s, and the 3-D arrival at mid-depth (z = 0) for reference
  peak rate    max of the depth-averaged slip rate
  speed        the steady rupture speed: a least-squares line of x against
               the arrival over x in [20, 50] km on each side
  implied c    for a 2.5D run, c_run x D_2D / D_3D averaged over the window
               (the steady 2.5D slip goes as 1 / c: Weng & Ampuero 2019 fit
               their 3-D runs with gamma = 0.96 / pi)
Pass (the note): final slip within 10% and arrival within 5% ("a few
percent") at every station of the window.

Output (--out DIR): summary.txt (the tables), summary.json, and
compare_2p5d_3d.pdf (final slip and arrival along strike, then the slip
rate and slip at five stations).

Example (from miniapps/seas):
  python3 elongated2p5d/scripts/compare_2p5d_3d.py --run3d R/3d \\
      --run2d wa2020=R/wa2020:1.17810 --run2d lehner=R/lehner:1.27324 \\
      --run2d "plane strain=R/plane_strain" --out R/compare
"""

import argparse
import glob
import json
import os
import sys

import numpy as np

TRAPZ = getattr(np, "trapezoid", None) or np.trapz

V_ARRIVAL = 1.0e-3      # m/s (SCEC rupture-time threshold)
W_FAULT = 15.0e3
Z_EDGE = 7.5e3


def prefix_of(d):
    f = glob.glob(os.path.join(d, "*_station_conventions.txt"))
    if len(f) != 1:
        sys.exit(f"ERROR: {d}: expected one *_station_conventions.txt, found {len(f)}")
    return os.path.basename(f[0])[:-len("_station_conventions.txt")]


def read_station(path):
    a = np.loadtxt(path, comments="#", ndmin=2)
    return {"t": a[:, 0], "slip": a[:, 2], "V": a[:, 4]}


def km(s):
    return float(s) * 1e3


def fmt_km(x):
    return ("%g" % (x / 1e3))


def stations_2d(d, p):
    out = {}
    for f in glob.glob(os.path.join(d, p + "_station_x2_*.dat")):
        name = os.path.basename(f)[len(p) + len("_station_x2_"):-4]
        if name.endswith("_copies") or "_z_" in name:
            continue
        out[km(name)] = read_station(f)
    return out


def stations_3d(d, p):
    out = {}
    for f in glob.glob(os.path.join(d, p + "_station_x2_*_z_*.dat")):
        name = os.path.basename(f)[len(p) + len("_station_x2_"):-4]
        if name.endswith("_copies"):
            continue
        xs, zs = name.split("_z_")
        out.setdefault(km(xs), {})[km(zs)] = read_station(f)
    return out


def depth_average(by_z, t):
    """Trapezoid over z of slip and V on the common times t, 0 at the edges."""
    zs = sorted(by_z)
    z = np.array([-Z_EDGE] + zs + [Z_EDGE])
    S = np.zeros((len(z), len(t)))
    V = np.zeros((len(z), len(t)))
    for i, zz in enumerate(zs):
        s = by_z[zz]
        S[i + 1] = np.interp(t, s["t"], s["slip"])
        V[i + 1] = np.interp(t, s["t"], s["V"])
    w = np.zeros(len(z))
    w[:-1] += 0.5 * np.diff(z)
    w[1:] += 0.5 * np.diff(z)
    w /= (z[-1] - z[0])
    return w @ S, w @ V


def arrival(t, V):
    i = np.nonzero(np.abs(V) > V_ARRIVAL)[0]
    if len(i) == 0:
        return float("nan")
    i = i[0]
    if i == 0:
        return float(t[0])
    # linear crossing between the two samples
    v0, v1 = abs(V[i - 1]), abs(V[i])
    return float(t[i - 1] + (t[i] - t[i - 1]) * (V_ARRIVAL - v0) / (v1 - v0))


def read_dumps(d, p):
    """Last dump of every rank: x, z, strike slip (copies averaged by point)."""
    files = sorted(glob.glob(os.path.join(d, p + "_fault_r*.csv")))
    if not files:
        return None
    rows = []
    for f in files:
        a = np.loadtxt(f, delimiter=",", comments="#", ndmin=2)
        if a.size == 0:
            continue
        rows.append(a[a[:, 0] == a[:, 0].max()])   # this rank's last dump (all dump the same steps)
    a = np.concatenate(rows)
    last = a[:, 0].max()
    a = a[a[:, 0] == last]
    t = float(a[0, 1])
    key = np.round(a[:, 2:5], 3)
    u, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    slip = np.bincount(inv, weights=a[:, 8]) / np.bincount(inv)
    return {"t": t, "step": int(last), "x": u[:, 0], "z": u[:, 2], "slip": slip}


def dump_depth_average(dm, xs, half_bin=250.0):
    out = []
    zg = np.linspace(-Z_EDGE, Z_EDGE, 301)
    for x in xs:
        m = np.abs(dm["x"] - x) <= half_bin
        if m.sum() < 4:
            out.append(float("nan"))
            continue
        z, s = dm["z"][m], dm["slip"][m]
        o = np.argsort(z)
        z, s = z[o], s[o]
        if z[0] > -Z_EDGE + 1.0:
            z, s = np.r_[-Z_EDGE, z], np.r_[0.0, s]
        if z[-1] < Z_EDGE - 1.0:
            z, s = np.r_[z, Z_EDGE], np.r_[s, 0.0]
        out.append(float(TRAPZ(np.interp(zg, z, s), zg) / (2.0 * Z_EDGE)))
    return np.array(out)


def speed(xs, ta, lo=20e3, hi=50e3):
    res = {}
    for side, sgn in (("+x", 1.0), ("-x", -1.0)):
        m = (sgn * xs >= lo) & (sgn * xs <= hi) & np.isfinite(ta)
        if m.sum() < 3:
            res[side] = float("nan")
            continue
        k = np.polyfit(ta[m], sgn * xs[m], 1)[0]
        res[side] = float(k)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run3d", required=True)
    ap.add_argument("--run2d", action="append", required=True, help="LABEL=DIR[:C]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--x-min", type=float, default=15.0, help="km")
    ap.add_argument("--x-max", type=float, default=57.5, help="km")
    ap.add_argument("--slip-tol", type=float, default=0.10)
    ap.add_argument("--arrival-tol", type=float, default=0.05)
    ap.add_argument("--plot-x", default="10,20,30,40,50", help="km, the time-series pages")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    p3 = prefix_of(a.run3d)
    s3 = stations_3d(a.run3d, p3)
    if not s3:
        sys.exit(f"ERROR: no 3-D stations {p3}_station_x2_*_z_*.dat in {a.run3d}")
    runs = []
    for spec in a.run2d:
        if "=" not in spec:
            sys.exit(f"ERROR: --run2d {spec}: expected LABEL=DIR[:C]")
        label, rest = spec.split("=", 1)
        c = None
        if ":" in rest:
            rest, cs = rest.rsplit(":", 1)
            c = float(cs)
        p = prefix_of(rest)
        st = stations_2d(rest, p)
        if not st:
            sys.exit(f"ERROR: no 2-D stations {p}_station_x2_*.dat in {rest}")
        runs.append({"label": label, "dir": rest, "prefix": p, "c": c, "st": st})

    # one window for all: t_end = the earliest end
    ends = [min(s["t"][-1] for zz in s3.values() for s in zz.values())]
    ends += [min(s["t"][-1] for s in r["st"].values()) for r in runs]
    t_end = min(ends)
    tg = np.linspace(0.0, t_end, 4001)

    xs = np.array(sorted(x for x in s3 if all(x in r["st"] for r in runs)))
    if len(xs) == 0:
        sys.exit("ERROR: no station x common to the 3-D run and every map-view run")
    D3, V3, A3, A3mid, P3 = [], [], [], [], []
    for x in xs:
        S, V = depth_average(s3[x], tg)
        D3.append(S[-1]); V3.append(V); A3.append(arrival(tg, V)); P3.append(float(np.max(np.abs(V))))
        z0 = min(s3[x], key=abs)
        A3mid.append(arrival(s3[x][z0]["t"], s3[x][z0]["V"]))
    D3, A3, A3mid, P3 = map(np.array, (D3, A3, A3mid, P3))
    dm = read_dumps(a.run3d, p3)
    D3dump = dump_depth_average(dm, xs) if dm is not None else None

    win = (np.abs(xs) >= a.x_min * 1e3 - 1.0) & (np.abs(xs) <= a.x_max * 1e3 + 1.0)
    lines = []
    lines.append(f"Level 2 (2.5D vs 3-D), t_end = {t_end:.3f} s (earliest end), window |x| in "
                 f"[{a.x_min:g}, {a.x_max:g}] km ({int(win.sum())} stations)")
    lines.append(f"3-D: {a.run3d} (prefix {p3})" + (f"; last dump step {dm['step']} at t = {dm['t']:.3f} s"
                                                     if dm is not None else "; no dumps"))
    summary = {"t_end": t_end, "window_km": [a.x_min, a.x_max], "x_km": (xs / 1e3).tolist(),
               "3d": {"final_slip": D3.tolist(), "arrival": A3.tolist(), "arrival_mid": A3mid.tolist(),
                      "peak_rate": P3.tolist(), "speed": speed(xs, A3),
                      "final_slip_dump": D3dump.tolist() if D3dump is not None else None},
               "runs": []}
    hdr = f"{'x km':>7} {'D3 m':>8}" + ("" if D3dump is None else f" {'D3dump':>8}") + f" {'t3 s':>7} {'t3mid':>7}"
    for r in runs:
        hdr += f" | {r['label'][:12]:>12} {'dD':>7} {'dt':>7}"
    rows = [hdr]
    allpass = True
    for r in runs:
        D, A, P, Vs = [], [], [], []
        for x in xs:
            s = r["st"][x]
            D.append(float(np.interp(t_end, s["t"], s["slip"])))
            A.append(arrival(s["t"], s["V"]))
            P.append(float(np.max(np.abs(s["V"][s["t"] <= t_end]))))
            Vs.append(np.interp(tg, s["t"], s["V"]))
        r["D"], r["A"], r["P"], r["V"] = np.array(D), np.array(A), np.array(P), Vs
        with np.errstate(divide="ignore", invalid="ignore"):   # stations the rupture missed: nan
            r["dD"] = (r["D"] - D3) / D3
            r["dA"] = (r["A"] - A3) / A3
        ok_s = bool(np.all(np.abs(r["dD"][win]) <= a.slip_tol))
        ok_a = bool(np.all(np.abs(r["dA"][win]) <= a.arrival_tol))
        r["pass"] = ok_s and ok_a
        if r["c"] is not None:
            allpass = allpass and r["pass"]
        r["speed"] = speed(xs, r["A"])
        r["c_implied"] = (float(r["c"] * np.mean(r["D"][win] / D3[win])) if r["c"] is not None else None)
        summary["runs"].append({"label": r["label"], "dir": r["dir"], "c": r["c"],
                                "final_slip": r["D"].tolist(), "arrival": r["A"].tolist(),
                                "peak_rate": r["P"].tolist(), "rel_slip": r["dD"].tolist(),
                                "rel_arrival": r["dA"].tolist(), "speed": r["speed"],
                                "c_implied": r["c_implied"], "pass_slip": ok_s, "pass_arrival": ok_a})
    for i, x in enumerate(xs):
        row = f"{x / 1e3:7.1f} {D3[i]:8.4f}" + ("" if D3dump is None else f" {D3dump[i]:8.4f}") \
              + f" {A3[i]:7.3f} {A3mid[i]:7.3f}"
        for r in runs:
            row += f" | {r['D'][i]:12.4f} {100 * r['dD'][i]:+6.1f}% {100 * r['dA'][i]:+6.1f}%"
        rows.append(row + ("" if win[i] else "   (outside the window)"))
    lines += rows
    lines.append("")
    sp = summary["3d"]["speed"]
    lines.append(f"3-D: steady rupture speed {sp['+x']:.0f} m/s (+x), {sp['-x']:.0f} m/s (-x); "
                 f"window mean final slip {np.mean(D3[win]):.4f} m")
    for r in runs:
        w = win
        lines.append(f"{r['label']}: final slip max |dD| {100 * np.max(np.abs(r['dD'][w])):.1f}% "
                     f"(mean {100 * np.mean(r['dD'][w]):+.1f}%), arrival max |dt| "
                     f"{100 * np.max(np.abs(r['dA'][w])):.1f}% (mean {100 * np.mean(r['dA'][w]):+.1f}%), "
                     f"speed {r['speed']['+x']:.0f} / {r['speed']['-x']:.0f} m/s"
                     + (f", implied c {r['c_implied']:.3f} (run c {r['c']:.3f})" if r["c"] is not None else "")
                     + (": PASS" if r["pass"] else ": FAIL") + ("" if r["c"] is not None else " (control)"))
    lines.append("2.5D runs " + ("all PASS" if allpass else "do not all pass")
                 + f" (slip within {100 * a.slip_tol:g}%, arrival within {100 * a.arrival_tol:g}%)")
    text = "\n".join(lines) + "\n"
    open(os.path.join(a.out, "summary.txt"), "w").write(text)
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    sys.stdout.write(text)

    if a.no_plot:
        return 0
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    colors = ["#1f5fa8", "#c4571a", "#2e8b57", "#7a4fa0", "#8a8a8a"]
    pdf = PdfPages(os.path.join(a.out, "compare_2p5d_3d.pdf"))
    fig, ax = plt.subplots(2, 1, figsize=(8.5, 8), sharex=True)
    ax[0].plot(xs / 1e3, D3, "k-o", ms=3, lw=2, label="3-D, depth average")
    if D3dump is not None:
        ax[0].plot(xs / 1e3, D3dump, "k--", lw=1, label="3-D, depth average (last dump)")
    for i, r in enumerate(runs):
        ax[0].plot(xs / 1e3, r["D"], "-s", ms=3, color=colors[i % 5], label=r["label"])
    ax[0].set_ylabel("final slip (m)")
    ax[0].legend(fontsize=8)
    ax[0].set_title(f"Final slip at t = {t_end:.2f} s")
    ax[1].plot(xs / 1e3, A3, "k-o", ms=3, lw=2, label="3-D, depth-averaged rate > 1 mm/s")
    ax[1].plot(xs / 1e3, A3mid, "k:", lw=1, label="3-D, mid-depth")
    for i, r in enumerate(runs):
        ax[1].plot(xs / 1e3, r["A"], "-s", ms=3, color=colors[i % 5], label=r["label"])
    ax[1].set_xlabel("x along strike (km)")
    ax[1].set_ylabel("arrival (s)")
    ax[1].legend(fontsize=8)
    for axx in ax:
        axx.axvspan(-a.x_min, a.x_min, color="0.92", zorder=0)
        axx.grid(alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)
    for xk in [float(v) for v in a.plot_x.split(",")]:
        i = int(np.argmin(np.abs(xs - xk * 1e3)))
        fig, ax = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)
        ax[0].plot(tg, V3[i], "k-", lw=2, label="3-D, depth average")
        S3, _ = depth_average(s3[xs[i]], tg)
        ax[1].plot(tg, S3, "k-", lw=2, label="3-D, depth average")
        for j, r in enumerate(runs):
            s = r["st"][xs[i]]
            ax[0].plot(s["t"], s["V"], color=colors[j % 5], lw=1, label=r["label"])
            ax[1].plot(s["t"], s["slip"], color=colors[j % 5], lw=1, label=r["label"])
        ax[0].set_ylabel("slip rate (m/s)")
        ax[1].set_ylabel("slip (m)")
        ax[1].set_xlabel("t (s)")
        ax[0].set_title(f"x = {xs[i] / 1e3:g} km")
        ax[0].legend(fontsize=8)
        for axx in ax:
            axx.grid(alpha=0.3)
            axx.set_xlim(0, t_end)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    pdf.close()
    print(f"wrote {a.out}/summary.txt, summary.json, compare_2p5d_3d.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
