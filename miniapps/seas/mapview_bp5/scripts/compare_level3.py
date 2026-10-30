#!/usr/bin/env python3
"""Compare map-view 2.5D BP5 cycle runs with the 3-D BP5-FD gold (Level 3).

The 2.5D note's Level 3 test (Chunhui Zhao, 2026-10-30;
/mnt/project-files/2d-inplane/nonplanar-literature/
inplane-2p5d-derivation-rev8.pdf, Section 7.3): SEAS BP5 seen from above with
the 2.5D finite-width term (config/fd_cycle/mapview_bp5_2p5d_fd.toml, buried,
W = 12 km; mapview_bp5_2p5d_surface_fd.toml, surface breaking, W = 16 km),
against the fully-dynamic-cycle gold
bp5/benchmark_data/golden_fd_vista_576r_p2_1000m (commit 3926e9c8).  The
note asks for recurrence, slip per event, moment, rupture extent and the
interseismic stressing rate within about 10-20%.

The map-view slip is the slip averaged over the width W, so the 3-D
quantities are averaged over depth 0..W first.  The gold's SCEC upload has
the whole depth profile only at x = 0 (slip_2_depth.dat, stress_2_depth.dat)
and the along-strike profile only at depth 10 km (slip_2_strike.dat,
stress_2_strike.dat), so:
  at x = 0       slip, stress drop and stressing rate are depth averages
                 over 0..W (trapezoid on the file's 1 km points), and the
                 10 km values are printed beside them;
  along strike   the gold curve is its 10 km profile (not a depth average:
                 the exact average needs the fault snapshots on Vista).

Sources:
  --gold DIR             the gold's scec_upload folder (default: the repo's)
  --run LABEL=DIR:W:C    a 2.5D run folder (repeatable, order kept), its
                         width W [m] and factor c; the folder holds the
                         driver's <prefix>_events.txt, _switch.txt,
                         _fltst_strk+00.txt, _profile_strike_slip_2.txt and
                         the run log run.log (the hand-back lines)
Events: the first event of each model starts from the nucleation patch at
t = 0, so the means use events 2 onward (recurrence: the intervals between
the starts of consecutive events from event 2 on).

Measures (each run against the gold, as run / gold - 1):
  recurrence     mean interval [yr]
  slip x = 0     event slip at the centre: the run's station, the gold's
                 depth average over 0..W
  moment         the run's mu W (integral of the event slip along strike, from
                 the profile, 1 km points); the gold's mu x catalog rupture
                 area x catalog average slip
  stress drop    shear stress at the centre before minus after each event
  stress rate    the centre's interseismic stressing rate [kPa/yr], two ways:
                 "mean" = the rise from just after an event to just before the
                 next over that time; "mid" = a linear fit over the middle
                 half of the interval (the gold loads faster early and late,
                 so mid < mean there)
  rupture length the run's _events.txt rupture_area (map view: the length
                 of the DOFs with |V| >= 1 mm/s, times 1 m) and the extent of
                 event slip > 1 m on its profile; the gold's 10 km profile
                 stops at +-38 km, so only a lower bound is given
  c to match     c x (mean rate / gold mean rate)^-1 and c x (slip / gold slip): the c
                 at this W that would match each measure, the 2.5D stressing
                 rate going as c and the slip as 1 / c
  hand-backs     every "back to QD" line of run.log: max V of the dynamic
                 phase and of the QD balance, against V_on = 0.01 m/s
Flags: within 10%, within 20%, or outside.

Output (--out DIR): summary.txt, summary.json and compare_level3.pdf (slip and
shear stress at x = 0 against time, max V against time, event slip along
strike for event 2 and the last event).

Example (from miniapps/seas):
  python3 mapview_bp5/scripts/compare_level3.py \\
      --run buried=R/mapview_bp5_2p5d_f27ca1ac:12000:1.1780972450961724 \\
      --run surface=R/mapview_bp5_2p5d_surface_f27ca1ac:16000:0.5890486225480862 \\
      --out R/compare
"""

import argparse
import glob
import json
import math
import os
import re

import numpy as np

YR = 3.15576e7       # Julian year [s], the gold's
MU = 32038120320.0   # BP5's shear modulus [Pa]
V_ON = 1.0e-2        # the BP5-FD switch [m/s]
TRAPZ = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2 renamed it
GOLD_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                            "bp5", "benchmark_data", "golden_fd_vista_576r_p2_1000m",
                            "scec_upload")


# ---------------------------------------------------------------- readers
def read_scec_rows(path):
    """Numeric rows of a SCEC file (comments and the field-name lines skipped)."""
    rows = []
    for line in open(path):
        tok = line.split()
        if not tok or line.startswith("#"):
            continue
        try:
            rows.append([float(t) for t in tok])
        except ValueError:
            continue
    return np.array(rows)


def read_profile(path, scec):
    """A profile file: first row the positions, then time, log10 Vmax, values."""
    a = read_scec_rows(path) if scec else np.loadtxt(path, comments="#", ndmin=2)
    return a[0, 2:], a[1:, 0], a[1:, 2:]


def before_after(t, vals, t0, t1):
    """Rows just before t0 and just after t1 (the file's own rows)."""
    i0 = max(np.searchsorted(t, t0, side="right") - 1, 0)
    i1 = min(np.searchsorted(t, t1, side="left"), len(t) - 1)
    return vals[i0], vals[i1]


def depth_average(z, v, W):
    """Trapezoid average of v(z) over z in [0, W] (z ascending, from 0)."""
    zz = np.append(z[z < W], W)
    vv = np.append(v[z < W], np.interp(W, z, v))
    return TRAPZ(vv, zz) / W


def interseismic_rates(t, tau, starts, ends):
    """Slope [per s] of tau over the middle half of each interval end_k .. start_k+1."""
    out = []
    for k in range(len(starts) - 1):
        a, b = ends[k], starts[k + 1]
        lo, hi = a + 0.25 * (b - a), a + 0.75 * (b - a)
        m = (t >= lo) & (t <= hi)
        if np.count_nonzero(m) >= 3:
            out.append(np.polyfit(t[m], tau[m], 1)[0])
        else:
            out.append(float("nan"))
    return np.array(out)


def interseismic_means(t, tau, starts, ends):
    """(tau just before start_k+1 - tau just after end_k) / (start_k+1 - end_k) [per s]."""
    out = []
    for k in range(len(starts) - 1):
        a, b = before_after(t, tau, starts[k + 1], ends[k])
        out.append((a - b) / (starts[k + 1] - ends[k]))
    return np.array(out)


def mean_from_2(x):
    x = np.asarray(x, dtype=float)
    return float(np.nanmean(x[1:])) if x.size > 1 else float("nan")


# ---------------------------------------------------------------- the gold
def load_gold(d):
    g = {}
    cat = read_scec_rows(os.path.join(d, "catalog.dat"))
    g["t_start"], g["t_end"] = cat[:, 1], cat[:, 2]
    g["area"], g["slip_avg"] = cat[:, 3], cat[:, 6]
    g["moment"] = MU * g["area"] * g["slip_avg"]
    g["z"], g["t_zs"], g["slip_z"] = read_profile(os.path.join(d, "slip_2_depth.dat"), True)
    _, g["t_zt"], g["tau_z"] = read_profile(os.path.join(d, "stress_2_depth.dat"), True)
    g["x"], g["t_xs"], g["slip_x"] = read_profile(os.path.join(d, "slip_2_strike.dat"), True)
    glob_rows = read_scec_rows(os.path.join(d, "global.dat"))
    g["t_glob"], g["logV_glob"] = glob_rows[:, 0], glob_rows[:, 1]
    return g


def gold_at_width(g, W):
    """The gold's centre measures, depth-averaged over 0..W and at 10 km."""
    i10 = int(np.argmin(abs(g["z"] - 10.0e3)))
    out = {"slip": [], "slip10": [], "drop": [], "drop10": []}
    for t0, t1 in zip(g["t_start"], g["t_end"]):
        s0, s1 = before_after(g["t_zs"], g["slip_z"], t0, t1)
        out["slip"].append(depth_average(g["z"], s1 - s0, W))
        out["slip10"].append(s1[i10] - s0[i10])
        a0, a1 = before_after(g["t_zt"], g["tau_z"], t0, t1)
        out["drop"].append(depth_average(g["z"], a0 - a1, W))
        out["drop10"].append(a0[i10] - a1[i10])
    tau_avg = np.array([depth_average(g["z"], row, W) for row in g["tau_z"]])
    out["rate"] = interseismic_rates(g["t_zt"], tau_avg, g["t_start"], g["t_end"]) * YR * 1e3
    out["rate10"] = interseismic_rates(g["t_zt"], g["tau_z"][:, i10], g["t_start"],
                                       g["t_end"]) * YR * 1e3
    out["mrate"] = interseismic_means(g["t_zt"], tau_avg, g["t_start"], g["t_end"]) * YR * 1e3
    out["mrate10"] = interseismic_means(g["t_zt"], g["tau_z"][:, i10], g["t_start"],
                                        g["t_end"]) * YR * 1e3
    slip_avg_t = np.array([depth_average(g["z"], row, W) for row in g["slip_z"]])
    out["series_slip"] = (g["t_zs"], slip_avg_t, g["slip_z"][:, i10])
    out["series_tau"] = (g["t_zt"], tau_avg, g["tau_z"][:, i10])
    return out


# ---------------------------------------------------------------- a run
def find_prefix(d):
    f = sorted(glob.glob(os.path.join(d, "*_events.txt")))
    if not f:
        raise SystemExit("no *_events.txt in %s" % d)
    return f[0][: -len("_events.txt")]


def load_run(label, d, W, c):
    p = find_prefix(d)
    r = {"label": label, "dir": d, "W": W, "c": c, "prefix": os.path.basename(p)}
    ev = np.loadtxt(p + "_events.txt", comments="#", ndmin=2)
    if ev.size == 0:
        ev = np.zeros((0, 7))
    r["t_start"], r["t_end"] = ev[:, 1], ev[:, 2]
    r["length_events"], r["slip_avg"] = ev[:, 3] / 1.0, ev[:, 6]
    st = np.loadtxt(p + "_fltst_strk+00.txt", comments="#", ndmin=2)
    r["t_st"], r["slip_st"], r["tau_st"] = st[:, 0], st[:, 1], st[:, 5]
    r["x"], r["t_pr"], r["slip_pr"] = read_profile(p + "_profile_strike_slip_2.txt", False)
    gl = np.loadtxt(p + "_global.txt", comments="#", ndmin=2)
    r["t_glob"], r["logV_glob"] = gl[:, 0], gl[:, 1]
    r["slip"], r["drop"], r["moment"], r["length_profile"], r["event_slip_x"] = [], [], [], [], []
    for t0, t1 in zip(r["t_start"], r["t_end"]):
        s0, s1 = before_after(r["t_st"], r["slip_st"], t0, t1)
        r["slip"].append(s1 - s0)
        a0, a1 = before_after(r["t_st"], r["tau_st"], t0, t1)
        r["drop"].append(a0 - a1)
        p0, p1 = before_after(r["t_pr"], r["slip_pr"], t0, t1)
        dD = p1 - p0
        r["event_slip_x"].append(dD)
        r["moment"].append(MU * W * TRAPZ(dD, r["x"]))
        big = r["x"][dD > 1.0]
        r["length_profile"].append(float(big.max() - big.min()) if big.size else 0.0)
    r["rate"] = interseismic_rates(r["t_st"], r["tau_st"], r["t_start"], r["t_end"]) * YR * 1e3
    r["mrate"] = interseismic_means(r["t_st"], r["tau_st"], r["t_start"], r["t_end"]) * YR * 1e3
    r["handbacks"] = []
    log = os.path.join(d, "run.log")
    if os.path.exists(log):
        pat = re.compile(r"back to QD at t = (\S+) s .*max V dynamic (\S+) m/s, QD balance (\S+) m/s")
        for line in open(log):
            m = pat.search(line)
            if m:
                r["handbacks"].append([float(m.group(i)) for i in (1, 2, 3)])
        r["stops"] = [line.strip() for line in open(log) if line.startswith("STOP")]
        r["status"] = [line.split(":", 1)[1].strip() for line in open(log)
                       if line.strip().startswith("status:")]
    else:
        r["stops"], r["status"] = [], []
    return r


# ---------------------------------------------------------------- compare
def flag(ratio):
    if not np.isfinite(ratio):
        return "n/a"
    e = abs(ratio - 1.0)
    return "within 10%" if e <= 0.10 else ("within 20%" if e <= 0.20 else "outside 20%")


def compare(g, r):
    gw = gold_at_width(g, r["W"])
    rec_r = np.diff(r["t_start"])[1:] / YR if len(r["t_start"]) > 2 else np.array([])
    rec_g = np.diff(g["t_start"])[1:] / YR
    rows = []

    def add(name, unit, run_v, gold_v, gold_note=""):
        ratio = run_v / gold_v if gold_v and np.isfinite(gold_v) else float("nan")
        rows.append({"measure": name, "unit": unit, "run": run_v, "gold": gold_v,
                     "ratio": ratio, "flag": flag(ratio), "gold_note": gold_note})

    add("recurrence", "yr", float(np.mean(rec_r)) if rec_r.size else float("nan"),
        float(np.mean(rec_g)), "events 2-7")
    add("slip x = 0", "m", mean_from_2(r["slip"]), mean_from_2(gw["slip"]),
        "depth avg 0..W (10 km: %.2f)" % mean_from_2(gw["slip10"]))
    add("moment", "N m", mean_from_2(r["moment"]), mean_from_2(g["moment"]), "mu A D (catalog)")
    add("stress drop x = 0", "MPa", mean_from_2(r["drop"]), mean_from_2(gw["drop"]),
        "depth avg 0..W (10 km: %.2f)" % mean_from_2(gw["drop10"]))
    add("stress rate x = 0, mean", "kPa/yr", mean_from_2(r["mrate"]), float(np.nanmean(gw["mrate"][1:])),
        "depth avg 0..W (10 km: %.1f)" % float(np.nanmean(gw["mrate10"][1:])))
    add("stress rate x = 0, mid", "kPa/yr", mean_from_2(r["rate"]), float(np.nanmean(gw["rate"][1:])),
        "depth avg 0..W (10 km: %.1f)" % float(np.nanmean(gw["rate10"][1:])))
    c_rate = r["c"] * rows[4]["gold"] / rows[4]["run"] if np.isfinite(rows[4]["run"]) else float("nan")
    c_slip = r["c"] * rows[1]["run"] / rows[1]["gold"] if np.isfinite(rows[1]["run"]) else float("nan")
    hb = np.array(r["handbacks"]) if r["handbacks"] else np.zeros((0, 3))
    return {
        "label": r["label"], "dir": r["dir"], "prefix": r["prefix"], "W_m": r["W"], "c": r["c"],
        "events": int(len(r["t_start"])), "event_starts_yr": (r["t_start"] / YR).tolist(),
        "intervals_yr": (np.diff(r["t_start"]) / YR).tolist(),
        "slip_x0_m": [float(v) for v in r["slip"]], "drop_x0_MPa": [float(v) for v in r["drop"]],
        "rate_x0_kPa_per_yr": [float(v) for v in r["rate"]],
        "mean_rate_x0_kPa_per_yr": [float(v) for v in r["mrate"]],
        "moment_Nm": [float(v) for v in r["moment"]],
        "length_events_m": r["length_events"].tolist(), "length_profile_m": r["length_profile"],
        "rows": rows, "c_to_match_rate": c_rate, "c_to_match_slip": c_slip,
        "handbacks": hb.tolist(), "handbacks_all_below_V_on": bool(np.all(hb[:, 2] < V_ON)) if len(hb) else None,
        "max_qd_balance_V": float(hb[:, 2].max()) if len(hb) else None,
        "stops": r["stops"], "status": r["status"], "gold_width": gw,
    }


def write_summary(out, g, res):
    lines = ["Level 3: map-view 2.5D BP5 against the 3-D BP5-FD gold (fdc 3926e9c8)",
             "gold: %d events, starts %s yr" % (len(g["t_start"]),
                                                ", ".join("%.1f" % v for v in g["t_start"] / YR)),
             "means over events 2 onward; ratio = run / gold", ""]
    for c in res:
        lines.append("== %s  (%s, W = %.0f m, c = %.6g)" % (c["label"], c["prefix"], c["W_m"], c["c"]))
        lines.append("   status %s; %d events at %s yr" % (
            "; ".join(c["status"]) or "?", c["events"],
            ", ".join("%.1f" % v for v in c["event_starts_yr"])))
        lines.append("   intervals [yr]: %s" % ", ".join("%.1f" % v for v in c["intervals_yr"]))
        lines.append("   rupture length [km]: events file %s; profile slip > 1 m %s" % (
            ", ".join("%.1f" % (v / 1e3) for v in c["length_events_m"]),
            ", ".join("%.1f" % (v / 1e3) for v in c["length_profile_m"])))
        lines.append("   %-22s %12s %12s %8s  %-12s %s" % ("measure", "run", "gold", "ratio", "flag", "gold value"))
        for row in c["rows"]:
            lines.append("   %-22s %12.4g %12.4g %8.3f  %-12s %s [%s]" % (
                row["measure"], row["run"], row["gold"], row["ratio"], row["flag"],
                row["gold_note"], row["unit"]))
        lines.append("   c to match at this W: %.4g (mean stress rate), %.4g (slip at x = 0)" % (
            c["c_to_match_rate"], c["c_to_match_slip"]))
        hb = c["handbacks"]
        lines.append("   hand-backs: %d, QD balance max V %s; all below V_on: %s" % (
            len(hb), "%.3g m/s" % c["max_qd_balance_V"] if hb else "n/a", c["handbacks_all_below_V_on"]))
        for t, vd, vq in hb:
            lines.append("      t = %.6g yr: dynamic %.3g -> QD balance %.3g m/s" % (t / YR, vd, vq))
        for s in c["stops"]:
            lines.append("   %s" % s)
        lines.append("")
    lines.append("Gold along strike: the 10 km depth profile only (to +-38 km); the gold's "
                 "rupture extends past +-38 km.")
    txt = "\n".join(lines) + "\n"
    open(os.path.join(out, "summary.txt"), "w").write(txt)
    js = [{k: v for k, v in c.items() if k != "gold_width"} for c in res]
    json.dump(js, open(os.path.join(out, "summary.json"), "w"), indent=1, default=float)
    return txt


def plot(out, g, runs, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    colors = ["#1f5f8b", "#b5462f", "#3f7d4e", "#8a7a3a"]
    with PdfPages(os.path.join(out, "compare_level3.pdf")) as pdf:
        for key, ylab, gi in (("slip", "slip at x = 0 [m]", "series_slip"),
                              ("tau", "shear stress at x = 0 [MPa]", "series_tau")):
            fig, ax = plt.subplots(figsize=(10, 5))
            for k, (r, c) in enumerate(zip(runs, res)):
                y = r["slip_st"] if key == "slip" else r["tau_st"]
                ax.plot(r["t_st"] / YR, y, color=colors[k % 4], lw=1.2,
                        label="%s (W %.0f km)" % (r["label"], r["W"] / 1e3))
                tg, ya, y10 = c["gold_width"][gi]
                ax.plot(tg / YR, ya, color=colors[k % 4], lw=0.9, ls="--",
                        label="gold, depth avg 0-%.0f km" % (r["W"] / 1e3))
            tg, ya, y10 = res[0]["gold_width"][gi] if res else (None, None, None)
            if tg is not None:
                ax.plot(tg / YR, y10, color="k", lw=0.9, label="gold, depth 10 km")
            ax.set_xlabel("time [yr]")
            ax.set_ylabel(ylab)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(g["t_glob"] / YR, g["logV_glob"], color="k", lw=0.9, label="gold")
        for k, r in enumerate(runs):
            ax.plot(r["t_glob"] / YR, r["logV_glob"], color=colors[k % 4], lw=1.0, label=r["label"])
        ax.axhline(math.log10(V_ON), color="0.5", lw=0.6, ls=":")
        ax.set_xlabel("time [yr]")
        ax.set_ylabel("log10 max V [m/s]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        pdf.savefig(fig)
        plt.close(fig)
        for which in (1, -1):
            fig, ax = plt.subplots(figsize=(10, 4.5))
            if len(g["t_start"]) > 1:
                k = which if which > 0 else len(g["t_start"]) - 1
                s0, s1 = before_after(g["t_xs"], g["slip_x"], g["t_start"][k], g["t_end"][k])
                ax.plot(g["x"] / 1e3, s1 - s0, color="k", lw=1.0,
                        label="gold event %d, depth 10 km" % (k + 1))
            for j, (r, c) in enumerate(zip(runs, res)):
                n = len(r["event_slip_x"])
                if n > 1:
                    k = which if which > 0 else n - 1
                    ax.plot(r["x"] / 1e3, r["event_slip_x"][k], color=colors[j % 4], lw=1.2,
                            label="%s event %d" % (r["label"], k + 1))
            ax.set_xlabel("x along strike [km]")
            ax.set_ylabel("event slip [m]")
            ax.set_title("event 2" if which > 0 else "last event")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gold", default=GOLD_DEFAULT)
    ap.add_argument("--run", action="append", required=True, help="LABEL=DIR:W:C")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    g = load_gold(a.gold)
    runs = []
    for spec in a.run:
        label, rest = spec.split("=", 1)
        d, W, c = rest.rsplit(":", 2)
        runs.append(load_run(label, d, float(W), float(c)))
    res = [compare(g, r) for r in runs]
    print(write_summary(a.out, g, res), end="")
    plot(a.out, g, runs, res)


if __name__ == "__main__":
    main()
