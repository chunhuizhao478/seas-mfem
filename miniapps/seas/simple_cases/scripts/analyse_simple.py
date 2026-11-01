#!/usr/bin/env python3
# The simple cases' summary (Chunhui Zhao, 2026-10-31: "go ahead with the
# simple cases, run them on local Mac"; plan "2D plane strain simple tests
# plan", "What we measure"): one table over the twenty-eight runs of
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
#   vertex jumps   the same file, per dynamic phase (2026-11-01, P0 of
#                  vertex-normal-plan-rev1): at p = 1 the snapshot points
#                  are per face, so each interior fault vertex has one DOF
#                  per face at the same place; the jump there is the
#                  largest minus the smallest sigma_n' of those DOFs.  The
#                  largest jump over the phase's dynamic snapshots, for the
#                  vertices with |x| <= 5 km (the arc and the corner) and
#                  over the whole fault, and the |x| <= 5 km one in the
#                  phase's last snapshot.  The file is opened without HDF5
#                  file locking, so a run still writing it does not stop
#                  the script (it then reads the snapshots written so far).
#   tip ends       the same file, per dynamic phase (2026-11-01, the tip
#                  zones): at the west and east free tips (the fault points
#                  of smallest and largest x), sigma_n' minus its t = 0
#                  value at the phase's last snapshot, the tip DOF's slip
#                  there and its peak |V| in the phase; and the largest
#                  |sigma_n' - sigma_n'(0)| of the DOFs within 1 km of the
#                  tip in x over the phase's snapshots (signed, where).
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
         "D1f_sharp015_floor", "D2f_sharp05_floor", "D5f_arc015L10_floor",
         "D5f_arc015L10_h250_floor", "D5f_arc015L10_h125_floor",
         "D5f_arc015L10_pf1_ev1", "D5f_arc015L10_pf2_ev1", "D5f_arc015L10_pf4_ev1",
         # vertex-normal plan P3 (event 1) and P4 (Chunhui Zhao, 2026-11-01)
         "D0_straight_ev1", "D0_straight_avg_ev1", "D5f_arc015L10_avg_ev1",
         "D5f_arc015L10_avg_floor", "D5f_arc015L10_h250_avg_floor",
         # the smoothed bends, 250 m, averaged (Chunhui Zhao, 2026-11-01: "we can
         # smooth the trace bends if needed")
         "D6f_arc03L10_h250_avg_floor", "D7f_arc08L10_h250_avg_floor",
         "D8f_arc03L2p5_h250_avg_floor", "D9f_arc08L2p5_h250_avg_floor",
         "D10f_sbend03L05_h250_avg_floor",
         # the tip zones on P4's 250 m D5f (Chunhui Zhao, 2026-11-01: "add
         # functionality about velocity strengthening patch near the tip")
         "D5f_arc015L10_h250_avg_floor_T1a08", "D5f_arc015L10_h250_avg_floor_T2w8",
         "D5f_arc015L10_h250_avg_floor_T3a08w8"]
# the vertex jumps' |x| range [km] where it is not 5 km: the +3 / -3 deg pair's
# arcs reach |x| = 6.25 km
XMAX_KM = {"D10f_sbend03L05_h250_avg_floor": 7.0}
BEND = ["bend_w", "bend_e"]
NEAR = ["bend_s-3.0", "bend_s-1.0", "bend_s+1.0", "bend_s+2.0", "bend_s+3.0",
        "bend_s-5.0", "bend_s+5.0"]
# geometry colours (categorical, fixed order): straight, 1.5 deg, 5 deg, arc, creeping 1.5 deg,
# 1.5 deg arc
GEO_COL = {"straight": "#4a5560", "sharp015": "#2a78d6", "sharp05": "#c0392b",
           "arc05L10": "#1e8a5a", "sharp015vs": "#b8860b", "arc015L10": "#7b4fb5"}
# the refined D5f runs: lighter steps of the 1.5 deg arc's purple
CASE_COL = {"D5f_arc015L10_h250_floor": "#a07fd0", "D5f_arc015L10_h125_floor": "#c6b2e6",
            # the P0 penalty scan: one orange hue, darker with the factor
            "D5f_arc015L10_pf1_ev1": "#e8a33d", "D5f_arc015L10_pf2_ev1": "#c96a24",
            "D5f_arc015L10_pf4_ev1": "#8f3a12",
            # averaged vertex normals (P3, P4): one teal hue, lighter when refined
            "D0_straight_avg_ev1": "#5f8f8c", "D5f_arc015L10_avg_ev1": "#0f7f7a",
            "D5f_arc015L10_avg_floor": "#0f7f7a", "D5f_arc015L10_h250_avg_floor": "#55b3ad",
            # the smoothed bends: blue 3 deg, red 8 deg (lighter: the 2.5 km arc),
            # green the +3 / -3 deg pair
            "D6f_arc03L10_h250_avg_floor": "#1f5fbf", "D8f_arc03L2p5_h250_avg_floor": "#6f9fe6",
            "D7f_arc08L10_h250_avg_floor": "#b3261e", "D9f_arc08L2p5_h250_avg_floor": "#e77a6f",
            "D10f_sbend03L05_h250_avg_floor": "#2f8a3e",
            # the tip zones: one brown hue, T1 / T2 / T3 dark to light
            "D5f_arc015L10_h250_avg_floor_T1a08": "#7a4a1e",
            "D5f_arc015L10_h250_avg_floor_T2w8": "#b07436",
            "D5f_arc015L10_h250_avg_floor_T3a08w8": "#d9a66b"}


def geo_of(case):
    return case.split("_")[1]


def col_of(case):
    # not CASE_COL.get(case, GEO_COL[...]): the default is looked up first, and
    # the smoothed bends' geometries have no GEO_COL entry
    return CASE_COL[case] if case in CASE_COL else GEO_COL[geo_of(case)]


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


def open_h5(path):
    """The snapshot file, read only and without HDF5 file locking where this
    h5py has the option (a run may still be writing it)."""
    import h5py
    try:
        return h5py.File(path, "r", locking=False)
    except TypeError:
        return h5py.File(path, "r")


def fault_snapshots(path):
    try:
        import h5py  # noqa: F401
    except ImportError:
        return dict(note="h5py missing: fault snapshots left out")
    if not os.path.exists(path):
        return dict(note="no fault snapshot file")
    with open_h5(path) as f:
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


def vertex_jumps(path, x_max_km=5.0, tol_m=1e-3):
    """Per dynamic phase, the sigma_n' jump at the fault vertices (file
    comment): {phase: dict(snapshots, jump_mpa, at_km, t_s, jump_all_mpa,
    all_at_km, jump_end_mpa)}, t_s counted from the phase's first snapshot;
    None without h5py or the file."""
    try:
        import h5py  # noqa: F401
    except ImportError:
        return None
    if not os.path.exists(path):
        return None
    with open_h5(path) as f:
        g = f["VTKHDF"]
        n = int(np.ravel(g["NumberOfPoints"])[0])
        X = g["Points"][:n, :2]
        t = g["Steps"]["Values"][:]
        off = g["Steps"]["PointDataOffsets"]["normal_stress"][:]
        fc = f["fd_cycle"]
        k = min(len(off), len(fc["regime"]), len(fc["phase"]))
        regime, phase = fc["regime"][:k], fc["phase"][:k]
        rows = [i for i in range(k) if regime[i] == 1]
        S = {i: g["PointData"]["normal_stress"][int(off[i]):int(off[i]) + n].ravel() for i in rows}
    if not rows:
        return {}
    scale = 1e6 if max(np.abs(s).max() for s in S.values()) > 1e4 else 1.0
    key = np.round(X / tol_m).astype(np.int64)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    order = np.argsort(inv, kind="stable")
    starts = np.r_[0, np.nonzero(np.diff(inv[order]))[0] + 1]
    shared = cnt >= 2
    # one representative point per vertex group
    gx = X[order[starts]]
    near = shared & (np.abs(gx[:, 0]) <= x_max_km * 1e3)
    out = {}
    for ph in sorted(set(int(phase[i]) for i in rows)):
        best, best_all, end = (-1.0, None, None), (-1.0, None), None
        prow = [i for i in rows if int(phase[i]) == ph]
        for i in prow:
            s = S[i][order] / scale
            jmp = np.maximum.reduceat(s, starts) - np.minimum.reduceat(s, starts)
            if near.any():
                j = int(np.argmax(np.where(near, jmp, -1.0)))
                if jmp[j] > best[0]:
                    best = (float(jmp[j]), j, float(t[i] - t[prow[0]]))
                end = float(jmp[near].max())
            ja = int(np.argmax(np.where(shared, jmp, -1.0)))
            if shared[ja] and jmp[ja] > best_all[0]:
                best_all = (float(jmp[ja]), ja)
        out[ph] = dict(snapshots=len(prow),
                       jump_mpa=(best[0] if best[1] is not None else None),
                       at_km=([float(gx[best[1], 0] / 1e3), float(gx[best[1], 1] / 1e3)]
                              if best[1] is not None else None),
                       t_s=best[2] if best[1] is not None else None,
                       jump_all_mpa=(best_all[0] if best_all[1] is not None else None),
                       all_at_km=([float(gx[best_all[1], 0] / 1e3), float(gx[best_all[1], 1] / 1e3)]
                                  if best_all[1] is not None else None),
                       jump_end_mpa=end)
    return out


def tip_ends(path, near_m=1000.0):
    """Per dynamic phase, sigma_n', slip and V at the two free tips (file
    comment): {phase: {"west"/"east": dict(...)}}; None without h5py or the
    file."""
    try:
        import h5py  # noqa: F401
    except ImportError:
        return None
    if not os.path.exists(path):
        return None
    with open_h5(path) as f:
        g = f["VTKHDF"]
        n = int(np.ravel(g["NumberOfPoints"])[0])
        X = g["Points"][:n, :2]
        t = g["Steps"]["Values"][:]
        po = g["Steps"]["PointDataOffsets"]
        fc = f["fd_cycle"]
        k = min(len(po["normal_stress"]), len(fc["regime"]), len(fc["phase"]))
        regime, phase = fc["regime"][:k], fc["phase"][:k]
        rows = [i for i in range(k) if regime[i] == 1]
        if k == 0 or not rows:
            return {}

        def field(name, i):
            o = int(po[name][i])
            return g["PointData"][name][o:o + n].ravel()
        sn0 = field("normal_stress", 0)
        scale = 1e6 if np.abs(sn0).max() > 1e4 else 1.0
        tips = {"west": int(np.argmin(X[:, 0])), "east": int(np.argmax(X[:, 0]))}
        near = {"west": X[:, 0] <= X[tips["west"], 0] + near_m,
                "east": X[:, 0] >= X[tips["east"], 0] - near_m}
        out = {}
        for ph in sorted(set(int(phase[i]) for i in rows)):
            prow = [i for i in rows if int(phase[i]) == ph]
            res = {e: dict(V_max=0.0, near_dsn=0.0, near_at_km=None) for e in tips}
            for i in prow:
                ds = (field("normal_stress", i) - sn0) / scale
                V = np.abs(field("slip_rate_strike", i))
                for e, j in tips.items():
                    r = res[e]
                    r["V_max"] = max(r["V_max"], float(V[j]))
                    m = np.nonzero(near[e])[0]
                    q = m[int(np.argmax(np.abs(ds[m])))]
                    if abs(ds[q]) > abs(r["near_dsn"]):
                        r["near_dsn"] = float(ds[q])
                        r["near_at_km"] = [float(X[q, 0] / 1e3), float(X[q, 1] / 1e3)]
            last = prow[-1]
            ds = (field("normal_stress", last) - sn0) / scale
            sl = field("slip_strike", last)
            for e, j in tips.items():
                res[e].update(tip_at_km=[float(X[j, 0] / 1e3), float(X[j, 1] / 1e3)],
                              sn0=float(sn0[j] / scale), dsn_end=float(ds[j]),
                              slip_end=float(sl[j]))
            out[ph] = res
    return out


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
    r["vertex_jumps"] = vertex_jumps(pre + "_fault_snapshots.vtkhdf", x_max_km=XMAX_KM.get(case, 5.0))
    r["tip_ends"] = tip_ends(pre + "_fault_snapshots.vtkhdf")
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
        # exactly <case>_<commit>/: a prefix match would take D0_straight_ev1_<commit>
        # for D0_straight
        ds = sorted(d for d in glob.glob(os.path.join(a.results, case + "_*/"))
                    if re.fullmatch(re.escape(case) + r"_[0-9a-f]{7,40}", os.path.basename(d.rstrip("/"))))
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
    L.append(f"{'case':36s} {'exit':>4s} {'wall':>6s} {'t_end':>8s} {'ev':>3s} "
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
        L.append(f"{case:36s} {fmt(r['exit'], '{:d}'):>4s} {fmt(r['wall_s'], '{:d}'):>6s} "
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
        vj = r["vertex_jumps"]
        for ph, v in sorted((vj or {}).items()):
            L.append(f"   vertex jump, dynamic phase {ph} ({v['snapshots']} snapshots): "
                     f"|x| <= {XMAX_KM.get(case, 5.0):g} km "
                     f"{fmt(v['jump_mpa'], '{:.2f}')} MPa at {v['at_km']} km, {fmt(v['t_s'], '{:.2f}')} s "
                     f"into the phase, {fmt(v['jump_end_mpa'], '{:.2f}')} at its end; whole fault "
                     f"{fmt(v['jump_all_mpa'], '{:.2f}')} at {v['all_at_km']} km")
        L.append("")
    L.append("Vertex jumps (MPa), the largest per dynamic phase, |x| <= 5 km unless noted (vertex jumps, "
             "file comment):")
    for case, r in runs.items():
        vj = r["vertex_jumps"]
        if vj:
            L.append(f"   {case:36s} " + " / ".join(fmt(v["jump_mpa"], "{:.2f}") for _, v in sorted(vj.items()))
                     + (f"   (|x| <= {XMAX_KM[case]:g} km)" if case in XMAX_KM else ""))
    L.append("")
    L.append("Tip ends (tip ends, file comment), per dynamic phase: sigma_n' - sigma_n'(0) at the "
             "tip DOF at the phase's end (MPa) west / east; the largest within 1 km of the tip "
             "(MPa) west / east; the tip DOF's slip (m) and peak |V| (m/s) west / east:")
    for case, r in runs.items():
        te = r.get("tip_ends")
        if not te:
            continue
        for ph, v in sorted(te.items()):
            w, e = v["west"], v["east"]
            L.append(f"   {case:36s} phase {ph}: tip {w['dsn_end']:+.2f} / {e['dsn_end']:+.2f}; "
                     f"within 1 km {w['near_dsn']:+.2f} / {e['near_dsn']:+.2f}; slip "
                     f"{w['slip_end']:.3f} / {e['slip_end']:.3f}; V {w['V_max']:.3g} / {e['V_max']:.3g}")
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
                    ls = "--" if "floor" in c else "-"
                    ax.plot(s["t"] / YR, s["sn"], ls, color=col_of(c), lw=1.4, label=c)
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
                ls = ":" if c.startswith("Q") else ("--" if "floor" in c else "-")
                ax.plot(s["slip"] - s["slip"][0], s["sn"], ls, color=col_of(c), lw=1.3, label=c)
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
