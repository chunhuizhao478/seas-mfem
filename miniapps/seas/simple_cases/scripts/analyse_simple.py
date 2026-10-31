#!/usr/bin/env python3
# The simple cases' summary (Chunhui Zhao, 2026-10-31: "go ahead with the
# simple cases, run them on local Mac"; plan "2D plane strain simple tests
# plan", "What we measure"): one table over the thirteen runs of
# run_simple_mac.sh, and the sigma_n' histories at the bend.
#
# Usage (from miniapps/seas):
#   python3 simple_cases/scripts/analyse_simple.py RESULTS [--out DIR]
# RESULTS holds the runs <case>_<commit>/ (run_simple_mac.sh).  Writes
# DIR/summary.txt, DIR/summary.json and DIR/simple_cases.pdf (DIR default
# RESULTS/compare).  Needs numpy and matplotlib; h5py for the fault
# snapshots (left out with a note when missing).
#
# Per run (the files of seas_fd_cycle_driver, prefix = case):
#   end state      RUN_INFO.txt (exit code, wall time), run.log (abort,
#                  STOP, wall limit lines; the first "[FRIC-GUARD]
#                  sigma_n_eff <= 0" line: where tension first reached the
#                  friction law; the floor's count line)
#   events         <case>_events.txt (t_start, t_end, rupture length =
#                  rupture_area / 1 m in 2-D, mean slip), <case>_catalog.txt
#                  (peak max V)
#   hand-backs     run.log's "back to QD" lines (max V dynamic, QD
#                  balance) and <case>_handback_phase<N>.txt: the fault-wide
#                  offset median(tau_qd - tau_dyn) over the DOFs, the largest
#                  |tau_qd - tau_dyn| and V_qd and where
#   the bend       the station files bend_w, bend_e (1 m from the vertex on
#                  each face), bend_s<+-km> and fltst_s<km>: sigma_n' at
#                  t = 0, its lowest value (when; in an event or between),
#                  the first time it is <= 0, sigma_n' at the end of each
#                  event, and the slip split into event slip (inside the
#                  events file's windows) and creep (the rest)
#   whole fault    <case>_fault_snapshots.vtkhdf: the lowest sigma_n' of any
#                  DOF over the snapshots (where, when)
import argparse
import glob
import json
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

YR = 365.0 * 86400.0
ORDER = ["Q0_straight", "Q1_sharp015", "Q2_sharp05", "Q3_arc05L10", "Q4_sharp015vs",
         "D0_straight", "D1_sharp015", "D2_sharp05", "D3_arc05L10", "D4_sharp015vs",
         "D1f_sharp015_floor", "D2f_sharp05_floor", "D5f_arc015L10_floor"]
BEND = ["bend_w", "bend_e"]
NEAR = ["bend_s-3.0", "bend_s-1.0", "bend_s+1.0", "bend_s+2.0", "bend_s+3.0",
        "bend_s-5.0", "bend_s+5.0"]
# geometry colours (categorical, fixed order): straight, 1.5 deg, 5 deg, arc, creeping 1.5 deg,
# 1.5 deg arc
GEO_COL = {"straight": "#4a5560", "sharp015": "#2a78d6", "sharp05": "#c0392b",
           "arc05L10": "#1e8a5a", "sharp015vs": "#b8860b", "arc015L10": "#7b4fb5"}


def geo_of(case):
    return case.split("_")[1]


def read_table(path):
    """Rows of a seas text table (comment lines '#'), as a 2-D array."""
    if not os.path.exists(path):
        return None
    rows = [ln.split() for ln in open(path, errors="replace") if ln.strip() and not ln.startswith("#")]
    if not rows:
        return np.zeros((0, 0))
    n = min(len(r) for r in rows)
    return np.array([[float(v) for v in r[:n]] for r in rows])


def scan_log(path):
    """The lines of run.log the summary needs (read once, line by line)."""
    out = dict(abort=None, stop=None, wall=None, fric_first=None, fric_count=0,
               fric_dropped=0, floor=None, backs=[], tables=[], status=None,
               final_yr=None, events=None, phases=None)
    if not os.path.exists(path):
        return out
    for ln in open(path, errors="replace"):
        s = ln.strip()
        if s.startswith("[run_simple_mac]"):
            m = re.search(r"\((\d+) in all\)", s)
            out["fric_dropped"] = int(m.group(1)) if m else 0
        elif "[FRIC-GUARD] sigma_n_eff <= 0" in s:
            out["fric_count"] += 1
            if out["fric_first"] is None:
                m = re.search(r"x2=(\S+) x3=(\S+) sn_eff=(\S+)", s)
                out["fric_first"] = dict(x_km=float(m.group(1)) / 1e3, y_km=float(m.group(2)) / 1e3,
                                         sn_eff_mpa=float(m.group(3)) / 1e6) if m else s[:160]
        elif ("MFEM abort" in s or "Verification failed" in s) and out["abort"] is None:
            out["abort"] = s[:200]
        elif s.startswith("STOP:"):
            out["stop"] = s[:200]
        elif "wall limit" in s and "reached" in s:
            out["wall"] = s[:200]
        elif s.startswith("sigma_n' floor") and out["floor"] is None:
            out["floor"] = s[:200]
        elif s.startswith("back to QD at t"):
            m = re.search(r"max V dynamic (\S+) m/s, QD balance (\S+) m/s", s)
            t = re.search(r"at t = (\S+) s", s)
            if m:
                out["backs"].append(dict(t_s=float(t.group(1)), V_dyn=float(m.group(1)),
                                         V_qd=float(m.group(2))))
        elif s.startswith("status:"):
            out["status"] = s.split(":", 1)[1].strip()
        elif s.startswith("final time:"):
            m = re.search(r"\(([\d.eE+-]+) yr\)", s)
            out["final_yr"] = float(m.group(1)) if m else None
        elif s.startswith("events:"):
            out["events"] = s.split(":", 1)[1].strip()
        elif s.startswith("dynamic phases:"):
            out["phases"] = s.split(":", 1)[1].strip()
    if out["fric_dropped"]:
        out["fric_count"] = out["fric_dropped"]
    return out


def handback_table(path):
    a = read_table(path)
    if a is None or a.size == 0:
        return None
    x, y = a[:, 0], a[:, 1]
    dtau = a[:, 10] - a[:, 9]
    dsn = a[:, 14] - a[:, 13]
    j, k = int(np.argmax(np.abs(dtau))), int(np.argmax(np.abs(a[:, 6])))
    return dict(dofs=len(a), offset_median_mpa=float(np.median(dtau)),
                dtau_max_mpa=float(dtau[j]), dtau_max_at_km=[float(x[j] / 1e3), float(y[j] / 1e3)],
                dsn_max_mpa=float(dsn[np.argmax(np.abs(dsn))]),
                Vqd_max=float(a[k, 6]), Vqd_max_at_km=[float(x[k] / 1e3), float(y[k] / 1e3)],
                sn_min_qd_mpa=float(a[:, 14].min()))


def in_events(t, ev):
    m = np.zeros(len(t), bool)
    for t0, t1 in ev:
        m |= (t >= t0) & (t <= t1)
    return m


def station(path, ev):
    """sigma_n' and slip at one station; ev = [(t_start, t_end)] in s."""
    a = read_table(path)
    if a is None or a.shape[0] < 1 or a.shape[1] < 9:
        return None
    t, slip, sn = a[:, 0], a[:, 1], a[:, 8]
    i = int(np.argmin(sn))
    ine = in_events(t, ev)
    dslip = np.diff(slip)
    # an increment belongs to an event when its mid time is inside one
    inc_ev = in_events(0.5 * (t[1:] + t[:-1]), ev)
    neg = np.nonzero(sn <= 0.0)[0]
    after = []
    for t0, t1 in ev:
        k = np.searchsorted(t, t1)
        after.append(float(sn[min(k, len(t) - 1)]))
    return dict(sn0=float(sn[0]), sn_min=float(sn[i]), sn_min_t_yr=float(t[i] / YR),
                sn_min_in_event=bool(ine[i]), sn_end=float(sn[-1]),
                first_le0_t_yr=(float(t[neg[0]] / YR) if len(neg) else None),
                first_le0_in_event=(bool(ine[neg[0]]) if len(neg) else None),
                sn_after_event=after, slip_end=float(slip[-1] - slip[0]),
                slip_events=float(dslip[inc_ev].sum()), slip_creep=float(dslip[~inc_ev].sum()),
                t=t, sn=sn, slip=slip)


def fault_snapshots(path):
    try:
        import h5py
    except ImportError:
        return dict(note="h5py missing: fault snapshots left out")
    if not os.path.exists(path):
        return dict(note="no fault snapshot file")
    with h5py.File(path, "r") as f:
        g = f["VTKHDF"]
        n = int(np.ravel(g["NumberOfPoints"])[0])
        X = g["Points"][:n, :2]
        t = g["Steps"]["Values"][:]
        off = g["Steps"]["PointDataOffsets"]["normal_stress"][:]
        S = np.array([g["PointData"]["normal_stress"][int(o):int(o) + n].ravel() for o in off])
    if S.size == 0:
        return dict(note="no snapshots")
    i, j = np.unravel_index(np.argmin(S), S.shape)
    neg = np.nonzero((S <= 0.0).any(axis=1))[0]
    return dict(snapshots=len(t), sn_min_mpa=float(S[i, j] / (1e6 if abs(S).max() > 1e4 else 1.0)),
                sn_min_at_km=[float(X[j, 0] / 1e3), float(X[j, 1] / 1e3)], sn_min_t_yr=float(t[i] / YR),
                first_le0_t_yr=(float(t[neg[0]] / YR) if len(neg) else None))


def run_info(d):
    p = os.path.join(d, "RUN_INFO.txt")
    s = open(p).read() if os.path.exists(p) else ""
    rc = re.search(r"exit (-?\d+) after (\d+) s", s)
    cm = re.search(r"^commit (\S+)", s, re.M)
    return dict(exit=int(rc.group(1)) if rc else None, wall_s=int(rc.group(2)) if rc else None,
                commit=cm.group(1) if cm else None)


def analyse(case, d):
    pre = os.path.join(d, case)
    r = dict(case=case, dir=d, **run_info(d))
    lg = scan_log(os.path.join(d, "run.log"))
    r.update({k: lg[k] for k in ("abort", "stop", "wall", "fric_first", "fric_count", "floor",
                                 "status", "final_yr", "phases")})
    evt = read_table(pre + "_events.txt")
    cat = read_table(pre + "_catalog.txt")
    ev, evrows = [], []
    if evt is not None and evt.size:
        for k, row in enumerate(evt):
            ev.append((row[1], row[2]))
            pk = float(cat[k, 3]) if cat is not None and cat.shape[0] > k and cat.shape[1] > 3 else float("nan")
            evrows.append(dict(event=int(row[0]), t_start_yr=float(row[1] / YR),
                               duration_s=float(row[2] - row[1]), length_km=float(row[3] / 1e3),
                               slip_mean_m=float(row[6]), peak_V=pk))
    r["events"] = evrows
    backs = lg["backs"]
    for n, b in enumerate(backs, start=1):
        b["table"] = handback_table(f"{pre}_handback_phase{n}.txt")
    r["handbacks"] = backs
    st = {}
    for name in BEND + NEAR + [f"fltst_s{k:+03d}" for k in range(-12, 13, 3)]:
        s = station(f"{pre}_{name}.txt", ev)
        if s is not None:
            st[name] = s
    r["stations"] = st
    r["fault"] = fault_snapshots(pre + "_fault_snapshots.vtkhdf")
    return r


def fmt(v, f="{:.3g}"):
    return "-" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f.format(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(a.results, "compare")
    os.makedirs(out, exist_ok=True)

    runs, missing = {}, []
    for case in ORDER:
        ds = sorted(glob.glob(os.path.join(a.results, case + "_*/")))
        if not ds:
            missing.append(case)
            continue
        runs[case] = analyse(case, ds[-1].rstrip("/"))

    L = ["Simple cases (simple_cases/README.md): sigma_n' in MPa, compression positive; t in yr.", ""]
    L.append("Columns: exit code; wall time (s); t_end (yr); events; lowest sigma_n' at bend_w / bend_e;")
    L.append("first sigma_n' <= 0 at bend_w or bend_e (yr; E in an event, Q between events); where the")
    L.append("friction law first saw sigma_n' <= 0 (FRIC-GUARD, km); per hand-back the QD balance max V")
    L.append("(m/s) / the fault-wide offset median(tau_qd - tau_dyn) (MPa).")
    L.append("")
    L.append(f"{'case':20s} {'exit':>4s} {'wall':>6s} {'t_end':>8s} {'ev':>3s} "
             f"{'sn_min w/e':>12s} {'bend <= 0':>11s} {'FRIC-GUARD at':>18s}  {'hand-backs':s}")
    for case, r in runs.items():
        bw, be = r["stations"].get("bend_w"), r["stations"].get("bend_e")
        snm = f"{fmt(bw and bw['sn_min'], '{:.1f}')}/{fmt(be and be['sn_min'], '{:.1f}')}"
        ff = r["fric_first"]
        fg = (f"({ff['x_km']:.2f}, {ff['y_km']:.2f})" if isinstance(ff, dict) else "-")
        firsts = sorted((s["first_le0_t_yr"], s["first_le0_in_event"]) for s in (bw, be)
                        if s and s["first_le0_t_yr"] is not None)
        b0 = f"{firsts[0][0]:.1f} {'E' if firsts[0][1] else 'Q'}" if firsts else "-"
        hb = "; ".join(f"{b['V_qd']:.2g}/{fmt(b['table'] and b['table']['offset_median_mpa'], '{:+.3f}')}"
                       for b in r["handbacks"]) or "-"
        L.append(f"{case:20s} {fmt(r['exit'], '{:d}'):>4s} {fmt(r['wall_s'], '{:d}'):>6s} "
                 f"{fmt(r['final_yr'], '{:.1f}'):>8s} {len(r['events']):>3d} {snm:>12s} {b0:>11s} {fg:>18s}  {hb}")
    L.append("")
    for case, r in runs.items():
        L.append(f"== {case}  ({r['dir']})")
        L.append(f"   commit {r['commit']}; exit {r['exit']}; wall {r['wall_s']} s; status {r['status']}; "
                 f"final {fmt(r['final_yr'], '{:.2f}')} yr; dynamic phases {r['phases']}")
        for k in ("abort", "stop", "wall", "floor"):
            if r[k]:
                L.append(f"   {k}: {r[k]}")
        if r["fric_count"]:
            L.append(f"   FRIC-GUARD sigma_n_eff <= 0: {r['fric_count']} lines; first {r['fric_first']}")
        for e in r["events"]:
            L.append(f"   event {e['event']}: t {e['t_start_yr']:.3f} yr, {e['duration_s']:.1f} s, "
                     f"length {e['length_km']:.1f} km, mean slip {e['slip_mean_m']:.3f} m, peak V {fmt(e['peak_V'])} m/s")
        for n, b in enumerate(r["handbacks"], start=1):
            tb = b["table"] or {}
            L.append(f"   hand-back {n}: max V dyn {b['V_dyn']:.3g} -> QD balance {b['V_qd']:.3g} m/s; "
                     f"offset median {fmt(tb.get('offset_median_mpa'), '{:+.3f}')} MPa, largest "
                     f"|dtau| {fmt(tb.get('dtau_max_mpa'), '{:+.3f}')} MPa at {tb.get('dtau_max_at_km')}, "
                     f"V_qd max at {tb.get('Vqd_max_at_km')}, min sn_qd {fmt(tb.get('sn_min_qd_mpa'), '{:.2f}')}")
        for name, s in r["stations"].items():
            if not name.startswith("bend"):
                continue
            L.append(f"   {name:11s} sn {s['sn0']:.2f} -> min {s['sn_min']:.2f} at {s['sn_min_t_yr']:.4f} yr "
                     f"({'event' if s['sn_min_in_event'] else 'between events'}), end {s['sn_end']:.2f}; "
                     f"first <= 0 {fmt(s['first_le0_t_yr'], '{:.4f}')} yr"
                     f"{'' if s['first_le0_in_event'] is None else (' (event)' if s['first_le0_in_event'] else ' (between events)')}"
                     f"; after events "
                     f"{', '.join(f'{v:.2f}' for v in s['sn_after_event']) or '-'}; slip {s['slip_end']:.3f} m "
                     f"(events {s['slip_events']:.3f}, creep {s['slip_creep']:.3f})")
        fs = r["fault"]
        if "note" in fs:
            L.append(f"   fault snapshots: {fs['note']}")
        else:
            L.append(f"   fault snapshots: {fs['snapshots']}; lowest sn {fs['sn_min_mpa']:.2f} at "
                     f"{fs['sn_min_at_km']} km, {fs['sn_min_t_yr']:.4f} yr; first <= 0 "
                     f"{fmt(fs['first_le0_t_yr'], '{:.4f}')} yr")
        L.append("")
    if missing:
        L.append("no run folder: " + ", ".join(missing))
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(L) + "\n")

    def strip(r):
        r = dict(r)
        r["stations"] = {k: {kk: vv for kk, vv in v.items() if kk not in ("t", "sn", "slip")}
                         for k, v in r["stations"].items()}
        return r
    json.dump({k: strip(v) for k, v in runs.items()}, open(os.path.join(out, "summary.json"), "w"), indent=1)

    with PdfPages(os.path.join(out, "simple_cases.pdf")) as pdf:
        for regime, title in (("Q", "QD only"), ("D", "fully dynamic")):
            cs = [c for c in runs if c.startswith(regime)]
            if not cs:
                continue
            fig, axs = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True)
            for c in cs:
                for ax, name in zip(axs, BEND):
                    s = runs[c]["stations"].get(name)
                    if s is None:
                        continue
                    ls = "--" if c.endswith("floor") else "-"
                    ax.plot(s["t"] / YR, s["sn"], ls, color=GEO_COL[geo_of(c)], lw=1.4, label=c)
            for ax, name in zip(axs, BEND):
                ax.axhline(0.0, color="#888", lw=0.8)
                ax.set_ylabel(f"sigma_n' at {name} (MPa)")
                ax.grid(alpha=0.25)
            axs[0].legend(fontsize=8, ncol=3)
            axs[1].set_xlabel("t (yr)")
            fig.suptitle(f"{title}: sigma_n' 1 m from the bend vertex, west (top) and east (bottom) face")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
        for c, r in runs.items():
            for ax, name in zip(axs, BEND):
                s = r["stations"].get(name)
                if s is None:
                    continue
                ls = ":" if c.startswith("Q") else ("--" if c.endswith("floor") else "-")
                ax.plot(s["slip"] - s["slip"][0], s["sn"], ls, color=GEO_COL[geo_of(c)], lw=1.3, label=c)
        for ax, name in zip(axs, BEND):
            ax.axhline(0.0, color="#888", lw=0.8)
            ax.set_xlabel(f"slip at {name} (m)")
            ax.set_ylabel("sigma_n' (MPa)")
            ax.grid(alpha=0.25)
        axs[0].legend(fontsize=7, ncol=2)
        fig.suptitle("sigma_n' against slip at the bend (QD only dotted, dynamic solid, floor dashed)")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print("\n".join(L[:len(runs) + 8]))
    print("wrote", os.path.join(out, "summary.txt"), "summary.json", "simple_cases.pdf")


if __name__ == "__main__":
    main()
