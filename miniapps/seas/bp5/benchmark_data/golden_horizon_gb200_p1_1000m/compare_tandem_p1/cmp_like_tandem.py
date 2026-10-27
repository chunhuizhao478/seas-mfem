"""Compare BP5-QD full-cycle runs the way compare_tandem_p1 (Vista gold
1074099 vs Tandem p=1) does:
  - years of 365.25 days;
  - a station's onset = first upward crossing of log10 V_strike = -3,
    interpolated linearly; an event's onset = the earliest station onset;
  - the peak log10 V_strike is searched from onset - 1 h to onset + 0.5 yr;
  - slip one year after the onset;
  - interseismic window k -> k+1 = [onset_k + 1 yr, onset_{k+1} - 1 day);
    differences interpolate one code onto the other's sample times in both
    directions and keep a point only where the interpolated code is sampled
    at least as densely there (bracketing interval of the interpolated code
    <= half the span of the sampled code's two neighbouring intervals).
Pairs: GPU-CPU, GPU-Tandem, CPU-Tandem.
Usage: cmp_like_tandem.py GPU_DIR CPU_DIR TANDEM_DIR OUT_DIR"""
import sys, os, json, csv
import numpy as np

GPU, CPU, TDM, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
YR = 365.25 * 86400.0
STATIONS = ["strk-36dp+00", "strk-16dp+00", "strk+00dp+00", "strk+16dp+00", "strk+36dp+00",
            "strk-24dp+10", "strk-16dp+10", "strk+00dp+10", "strk+16dp+10", "strk+00dp+22"]

def tname(s):
    return f"bp5qd_tandem_p1_x2_{int(s[4:7])}_x3_{int(s[9:12])}.txt"

def load(code, s):
    if code == "Tandem":
        rows = [[float(v) for v in ln.split()] for ln in open(os.path.join(TDM, tname(s)))
                if ln.strip() and ln.strip()[0] not in "#t"]
        a = np.array(rows)
    else:
        a = np.loadtxt(os.path.join(GPU if code == "GPU" else CPU, f"bp5_cycle_p1_fltst_{s}.txt"), comments="#")
    # time(s), slip_strike, log10 V_strike, tau_strike  (same column order in both formats)
    return {"t": a[:, 0], "slip": a[:, 1], "lv": a[:, 3], "tau": a[:, 5], "state": a[:, 7]}

CODES = ["GPU", "CPU", "Tandem"]
D = {c: {s: load(c, s) for s in STATIONS} for c in CODES}
T_END = {c: min(D[c][s]["t"][-1] for s in STATIONS) for c in CODES}

def crossings(t, lv, thr=-3.0):
    """Upward crossings of thr, linearly interpolated (s)."""
    below = lv[:-1] < thr
    up = np.flatnonzero(below & (lv[1:] >= thr))
    out = [t[i] + (thr - lv[i]) * (t[i + 1] - t[i]) / (lv[i + 1] - lv[i]) for i in up]
    if lv[0] >= thr:
        out = [t[0]] + out
    return np.array(out)

def event_onsets(code):
    allc = np.sort(np.concatenate([crossings(D[code][s]["t"], D[code][s]["lv"]) for s in STATIONS]))
    ev = []
    for x in allc:   # cluster crossings into events: a new event starts > 1 yr after the last onset
        if not ev or x - ev[-1] > YR:
            ev.append(x)
    return np.array(ev)

ON = {c: event_onsets(c) for c in CODES}

def station_event(code, s, k):
    d = D[code][s]; t0 = ON[code][k]
    cr = crossings(d["t"], d["lv"])
    cr = cr[(cr >= t0 - YR) & (cr < t0 + 0.5 * YR)]
    onset = cr[0] if len(cr) else None
    m = (d["t"] >= t0 - 3600.0) & (d["t"] <= t0 + 0.5 * YR)
    i = np.flatnonzero(m)[np.argmax(d["lv"][m])]
    slip1 = np.interp(t0 + YR, d["t"], d["slip"]) if t0 + YR <= d["t"][-1] else None
    return onset, d["lv"][i], d["t"][i], slip1

def dense_mask(ta, tb):
    """For points ta (sampled code A), True where code B (interpolated) is locally at least as dense."""
    ha = np.empty_like(ta)
    ha[1:-1] = 0.5 * (ta[2:] - ta[:-2]); ha[0] = ta[1] - ta[0]; ha[-1] = ta[-1] - ta[-2]
    j = np.clip(np.searchsorted(tb, ta) - 1, 0, len(tb) - 2)
    hb = tb[j + 1] - tb[j]
    inside = (ta >= tb[0]) & (ta <= tb[-1])
    return inside & (hb <= ha * (1 + 1e-12))

def compare_pair(a, b):
    nev = min(len(ON[a]), len(ON[b]))
    res = {"pair": f"{a} - {b}", "event_onsets": [], "events": [], "interseismic": []}
    for k in range(nev):
        res["event_onsets"].append({"event": k + 1, a: ON[a][k] / YR, b: ON[b][k] / YR,
                                    "diff_s": ON[a][k] - ON[b][k]})
        for s in STATIONS:
            oa, pa, tpa, sa = station_event(a, s, k)
            ob, pb, tpb, sb = station_event(b, s, k)
            res["events"].append({
                "event": k + 1, "station": s,
                "onset_diff_s": (oa - ob) if (oa is not None and ob is not None) else None,
                "peak_log10V_a": pa, "peak_log10V_b": pb, "peak_log10V_diff": pa - pb,
                "peak_time_diff_s": tpa - tpb,
                "slip_onset_plus_1yr_a_m": sa, "slip_onset_plus_1yr_b_m": sb,
                "slip_diff_m": (sa - sb) if (sa is not None and sb is not None) else None})
    t_end = min(T_END[a], T_END[b])
    for k in range(nev):
        lo = max(ON[a][k], ON[b][k]) + YR
        if k + 1 >= nev:
            continue   # as compare_tandem_p1: windows between events only
        hi = min(ON[a][k + 1], ON[b][k + 1]) - 86400.0
        if hi <= lo:
            continue
        for s in STATIONS:
            da, db = D[a][s], D[b][s]
            best = {"tau": (0.0, None), "slip": (0.0, None), "lv": (0.0, None)}
            npts = 0
            for (P, Q, sign) in ((da, db, 1.0), (db, da, -1.0)):   # P sampled, Q interpolated
                m = (P["t"] >= lo) & (P["t"] < hi)
                tp = P["t"][m]
                if len(tp) == 0:
                    continue
                keep = dense_mask(P["t"], Q["t"])[m]
                tq = tp[keep]; npts += len(tq)
                for f in ("tau", "slip", "lv"):
                    d = np.abs(np.interp(tq, Q["t"], Q[f]) - P[f][m][keep])
                    if len(d) and d.max() > best[f][0]:
                        i = d.argmax(); best[f] = (d[i], tq[i] / YR)
            res["interseismic"].append({
                "window": f"{k+1}->{k+2}" if k + 1 < nev else f"{k+1}->end", "station": s,
                "points_compared": npts,
                "max_abs_dtau_MPa": best["tau"][0], "at_yr_tau": best["tau"][1],
                "max_abs_dslip_m": best["slip"][0], "at_yr_slip": best["slip"][1],
                "max_abs_dlog10V": best["lv"][0], "at_yr_log10V": best["lv"][1]})
    return res

pairs = [("GPU", "CPU"), ("GPU", "Tandem"), ("CPU", "Tandem")]
summary = {"years": "365.25 days", "end_yr": {c: T_END[c] / YR for c in CODES},
           "n_events": {c: len(ON[c]) for c in CODES}, "pairs": {}}
for a, b in pairs:
    r = compare_pair(a, b)
    tag = f"{a}_vs_{b}".lower()
    summary["pairs"][tag] = r
    for name in ("events", "interseismic"):
        rows = r[name]
        with open(os.path.join(OUT, f"{tag}_{name}.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)

# console digest
DEPTH6 = ["strk+00dp+00", "strk+00dp+10", "strk+00dp+22", "strk-16dp+10", "strk+16dp+10", "strk-24dp+10"]
SURF4 = ["strk-36dp+00", "strk-16dp+00", "strk+16dp+00", "strk+36dp+00"]
print("events:", {c: len(ON[c]) for c in CODES}, " ends (yr):", {c: round(T_END[c] / YR, 1) for c in CODES})
print("event onsets (yr, GPU):", ", ".join(f"{x/YR:.3f}" for x in ON["GPU"]))
for tag, r in summary["pairs"].items():
    print(f"\n== {r['pair']}")
    print("  onset diff (s):", ", ".join(f"{e['diff_s']:+.2f}" for e in r["event_onsets"]))
    for grp, lab in ((DEPTH6, "6 stations (center + depth)"), (SURF4, "4 off-center surface")):
        ev = [e for e in r["events"] if e["station"] in grp]
        it = [e for e in r["interseismic"] if e["station"] in grp]
        pk = max(abs(e["peak_log10V_diff"]) for e in ev)
        sl = max(abs(e["slip_diff_m"]) for e in ev if e["slip_diff_m"] is not None)
        print(f"  {lab}: max|peak log10V diff| {pk:.2e}, max|slip(onset+1yr) diff| {sl:.2e} m, "
              f"between events max|dtau| {max(e['max_abs_dtau_MPa'] for e in it):.2e} MPa, "
              f"max|dslip| {max(e['max_abs_dslip_m'] for e in it):.2e} m, "
              f"max|dlog10V| {max(e['max_abs_dlog10V'] for e in it):.2e}")
