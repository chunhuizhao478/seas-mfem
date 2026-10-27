#!/usr/bin/env python3
"""BP3-QD (2-D, quasi-dynamic cycles): seas-mfem against Tandem p = 1, station by station.

Gate G3 of the 2-D in-plane design (Chunhui Zhao, 2026-10-27: "can we only plot tandem p = 1
with our results"; the gold is bp3/gold/qd-2d_p1_d60_thrust_25m_mac_git9503fe2).  One PDF of
PNG pages at 110 dpi: the event table, the onsets against event number, one page per fault
station (12) and per surface station (8), and a notes page.  Each station page shows 0 to the
cut, ours minus Tandem between events, the last event's window and ours minus Tandem there.

Events are matched by index.  Both series are cut at T_CUT = the earlier of the two codes' last
output times; both codes must have the same number of events, each starting at least 1 yr before
the cut.  365-day years.  Onsets of both codes come from the station data with one definition
(first time max V over the 12 fault stations > 1e-3 m/s after more than 1 yr below).
Differences: the locally denser series is interpolated onto the sparser one's times, never the
other way round (decided sample by sample).  Tandem lighter underneath, ours last, on top, red,
heavier.  Tandem's surface velocity is du/dt by finite differences (its probes write no velocity).

Usage (from anywhere):
  plot_bp3_qd_tandem_compare.py OURS_SCEC_DIR OURS_EVENTS_TXT TANDEM_DIR OUT_DIR
OURS_SCEC_DIR = the cycle/scec_bp3_export.py files of the run (fltst_dp*.dat, srfst_fn*.dat),
OURS_EVENTS_TXT = the driver's <prefix>_events.txt, TANDEM_DIR = Tandem's probe files
(fltst_dp*.csv, srfst_fn*.csv).  Needs numpy, matplotlib and Pillow.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

if len(sys.argv) != 5:
    sys.exit(__doc__)
OURS_DIR, OURS_EVENTS, TANDEM_DIR, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
YR = 365 * 86400.0
FAULT = [f"dp{d:03d}" for d in (0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350)]
SURF = ["fn-32", "fn-16", "fn-08", "fn-00", "fn+00", "fn+08", "fn+16", "fn+32"]
OURS, TP1 = "seas-mfem p=1", "Tandem p=1 (Mac)"
CODES = [(TP1, "#1f77b4"), (OURS, "#d62728")]                 # Tandem p = 1 and ours only
STY = {lab: (dict(color=c, lw=2.0, alpha=1.0, zorder=10) if lab == OURS else dict(color=c, lw=1.0, alpha=0.6, zorder=2))
       for lab, c in CODES}
W0, W1 = -100.0, 300.0                                  # event window, s from each code's own onset


def numeric(path, sep=None):
    rows = []
    for l in open(path):
        s = l.strip()
        if not s or s.startswith("#") or s[0].isalpha() or s[0] == '"':
            continue
        try:
            rows.append([float(x) for x in (s.split(",") if sep else s.split())])
        except ValueError:
            pass
    return np.array(rows)


def fault(code, st):
    """dict t, slip (m), logV (log10 m/s), tau (MPa, |.|), sn (MPa)."""
    if code == TP1:
        a = numeric(f"{TANDEM_DIR}/fltst_{st}.csv", sep=",")
        return dict(t=a[:, 0], slip=a[:, 2], logV=np.log10(np.abs(a[:, 4]) + 1e-300), tau=np.abs(a[:, 3]), sn=a[:, 5])
    a = numeric(f"{OURS_DIR}/fltst_{st}.dat")
    return dict(t=a[:, 0], slip=a[:, 1], logV=a[:, 2], tau=np.abs(a[:, 3]), sn=a[:, 4])


def surf(code, st):
    """dict t, d1, d2 (SCEC: e2 downward), v1, v2."""
    if code == OURS:
        a = numeric(f"{OURS_DIR}/srfst_{st}.dat")
        return dict(t=a[:, 0], d1=a[:, 1], d2=a[:, 2], v1=a[:, 3], v2=a[:, 4], vfd=False)
    a = numeric(f"{TANDEM_DIR}/srfst_{st}.csv", sep=",")
    t, d1, d2 = a[:, 0], a[:, 1], -a[:, 2]              # Tandem u1 is up; SCEC disp_2 is down
    return dict(t=t, d1=d1, d2=d2, v1=np.gradient(d1, t), v2=np.gradient(d2, t), vfd=True)


def cut(d, tc):
    m = d["t"] <= tc
    return {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in d.items()}


def _onto(ts, ys, ti, yi):
    """yi interpolated onto the samples ts, kept only where ti is locally at least as dense:
    the gap of ti bracketing the sample is no longer than the larger gap of ts around it."""
    j = np.clip(np.searchsorted(ti, ts), 1, ti.size - 1)
    gi = ti[j] - ti[j - 1]
    dts = np.diff(ts)
    hs = np.maximum(np.r_[dts[:1], dts], np.r_[dts, dts[-1:]])
    keep = (gi <= hs * (1 + 1e-9)) & (ts >= ti[0]) & (ts <= ti[-1])
    return ts[keep], np.interp(ts[keep], ti, yi), ys[keep]


def diff_series(ta, ya, tb, yb, lo, hi):
    """(times, ours - ref) over [lo, hi]: wherever one series is locally denser it is
    interpolated onto the other's samples (decided sample by sample, so a series is never
    interpolated across one of its own output gaps)."""
    ma, mb = (ta >= lo) & (ta <= hi), (tb >= lo) & (tb <= hi)
    if ma.sum() < 2 or mb.sum() < 2:
        return np.array([]), np.array([])
    t1, ia, yb1 = _onto(tb[mb], yb[mb], ta, ya)          # ours interpolated onto ref's samples
    t2, ib, ya2 = _onto(ta[ma], ya[ma], tb, yb)          # ref interpolated onto ours
    ts = np.r_[t1, t2]; d = np.r_[ia - yb1, ya2 - ib]
    o = np.argsort(ts, kind="stable")
    return ts[o], d[o]


def events(F):
    """Onsets (first time max V over the 12 fault stations > 1e-3 m/s, after > 1 yr below), ends
    (last time > 1e-3 within 2000 s), peak V and its station."""
    ts = np.concatenate([F[st]["t"][F[st]["logV"] > -3] for st in FAULT])
    ts.sort()
    on = [ts[0]] + [ts[i] for i in range(1, ts.size) if ts[i] - ts[i - 1] > YR] if ts.size else []
    ev = []
    for o in on:
        end, pk, pkst = -np.inf, -np.inf, ""
        for st in FAULT:
            d = F[st]; m = (d["t"] >= o - 1) & (d["t"] <= o + 2000)
            i = np.nonzero((d["logV"] > -3) & m)[0]
            if i.size:
                end = max(end, d["t"][i[-1]])
            if m.any():
                j = np.argmax(np.where(m, d["logV"], -np.inf))
                if d["logV"][j] > pk:
                    pk, pkst = d["logV"][j], st
        ev.append(dict(onset=o, end=end, peak=10 ** pk, peak_st=pkst))
    return ev


# ---------------------------------------------------------------- load, events, cut
Ffull = {code: {st: fault(code, st) for st in FAULT} for code, _ in CODES}
EVfull = {code: events(Ffull[code]) for code, _ in CODES}
ours_ev = numeric(OURS_EVENTS)                         # the driver's catalog: complete events
N = len(ours_ev)
T_END = {c: min(Ffull[c][st]["t"][-1] for st in FAULT) for c, _ in CODES}
T_CUT = min(T_END.values())
for code, _ in CODES:
    assert len(EVfull[code]) == N, (code, len(EVfull[code]), N)
    assert EVfull[code][N - 1]["onset"] + YR < T_CUT, (code, EVfull[code][N - 1]["onset"] / YR, T_CUT / YR)
t_end_ours = T_END[OURS]
F = {c: {st: cut(Ffull[c][st], T_CUT) for st in FAULT} for c, _ in CODES}
S = {c: {st: cut(surf(c, st), T_CUT) for st in SURF} for c, _ in CODES}
EV = {c: EVfull[c][:N] for c, _ in CODES}
for e_k in range(N):                                   # station onsets vs the catalog (all DOFs): within 60 s
    assert abs(EV[OURS][e_k]["onset"] - ours_ev[e_k, 1]) < 60.0, (e_k, EV[OURS][e_k]["onset"], ours_ev[e_k, 1])


def slip1(code, k, st):
    d = F[code][st]; o = EV[code][k]["onset"]; j = np.searchsorted(d["t"], o + YR)
    return d["t"][min(j, d["t"].size - 1)], d["slip"][min(j, d["t"].size - 1)]


def slip_pair(a, b, k, st):
    """slip of a and b one year after each code's own onset of event k, read at the same time after
    onset: the sparser code's first sample at or after onset + 1 yr is the point; the denser code is
    interpolated at the same relative time (afterslip makes the sample time matter)."""
    da, db = F[a][st], F[b][st]
    ea, eb = EV[a][k]["onset"], EV[b][k]["onset"]
    sa = np.median(np.diff(da["t"][(da["t"] > ea + 0.5 * YR) & (da["t"] < ea + 1.5 * YR)]))
    sb = np.median(np.diff(db["t"][(db["t"] > eb + 0.5 * YR) & (db["t"] < eb + 1.5 * YR)]))
    if sa >= sb:
        ta, ya = slip1(a, k, st); yb = np.interp(ta - ea + eb, db["t"], db["slip"])
    else:
        tb, yb = slip1(b, k, st); ya = np.interp(tb - eb + ea, da["t"], da["slip"])
    return {a: ya, b: yb}


def slip_diff(a, b, k, st):
    """ours minus ref, slip one year after each code's own onset of event k (slip_pair)."""
    p = slip_pair(a, b, k, st)
    return p[a] - p[b]


# ---------------------------------------------------------------- (a) event table
L = [f"BP3-QD dip 60 thrust, 0 to {T_CUT / YR:.2f} yr, events 1 to {N} (cut = the earlier end of the two runs;"
     " 365-day years)",
     f"onset = first time max slip rate over the 12 fault stations > 1e-3 m/s; dur = last time > 1e-3 m/s minus onset;",
     "peak = max slip rate over the fault stations in the event; slip = each code's slip one year after its own onset, both codes",
     "read at the same time after onset (m; dp = km x 10)", ""]
hdr = f"{'ev':>2s} {'code':18s} {'onset (yr)':>10s} {'dur (s)':>7s} {'peak V':>7s} {'at':>5s} " + " ".join(f"{s[2:]:>6s}" for s in FAULT)
L.append(hdr)
for k in range(N):
    for code, _ in CODES:
        e = EV[code][k]
        L.append(f"{k + 1:2d} {code:18s} {e['onset'] / YR:10.4f} {e['end'] - e['onset']:7.1f} {e['peak']:7.4f} {e['peak_st']:>5s} "
                 + " ".join(f"{slip_pair(OURS, TP1, k, s)[code]:6.3f}" for s in FAULT))
    for ref, tag in ((TP1, "ours - T p1"),):
        e, r = EV[OURS][k], EV[ref][k]
        L.append(f"{k + 1:2d} {tag:18s} {(e['onset'] - r['onset']) / YR:+10.4f} {(e['end'] - e['onset']) - (r['end'] - r['onset']):+7.1f} "
                 f"{e['peak'] - r['peak']:+7.4f} {'':5s} " + " ".join(f"{slip_diff(OURS, ref, k, s):+6.3f}" for s in FAULT))
    L.append("")
L.append("mean recurrence (yr) and onset drift")
L.append(f"{'code':18s} {'from ev 2':>10s} {'from ev 3':>10s} {'ours-code onset ev1':>20s} {'ev' + str(N):>10s} {'slope (yr/event)':>17s}")
REC = {}
for code, _ in CODES:
    o = np.array([e["onset"] for e in EV[code]]) / YR
    REC[code] = (np.mean(np.diff(o[1:])), np.mean(np.diff(o[2:])))
    d = np.array([e["onset"] for e in EV[OURS]]) / YR - o
    L.append(f"{code:18s} {REC[code][0]:10.3f} {REC[code][1]:10.3f} {d[0]:+20.4f} {d[-1]:+10.4f} {np.polyfit(np.arange(N), d, 1)[0]:+17.5f}")
L.append("(mean recurrence from ev 2 = mean of onset(k+1) - onset(k) for k >= 2, i.e. intervals ev2->ev3 ... ev{0}; from ev 3 likewise)".format(N))
open(os.path.join(OUT, "event_table.txt"), "w").write("\n".join(L) + "\n")
pages, n = [], 1
PER = 62
for i0 in range(0, len(L), PER):
    fig = plt.figure(figsize=(17, 11))
    fig.text(0.01, 0.99, "\n".join(L[i0:i0 + PER] if i0 == 0 else L[:5] + ["(continued)"] + L[i0:i0 + PER]),
             family="monospace", fontsize=7.6, va="top")
    fn = f"p{n:02d}_event_table_{i0 // PER + 1}.png"; fig.savefig(os.path.join(OUT, fn), dpi=110); plt.close(fig)
    pages.append(fn); n += 1

# ---------------------------------------------------------------- (b) summary figure
fig, ax = plt.subplots(1, 3, figsize=(18, 5.5))
k1 = np.arange(1, N + 1)
for code, _ in CODES:
    o = np.array([e["onset"] for e in EV[code]]) / YR
    ax[0].plot(k1, o, "o-", label=code, **STY[code])
    if code != OURS:
        oo = np.array([e["onset"] for e in EV[OURS]]) / YR
        ax[1].plot(k1, oo - o, "o-", label=f"ours - {code}", color=STY[code]["color"], lw=1.2)
        ax[2].plot(k1[1:], np.diff(oo) - np.diff(o), "o-", label=f"ours - {code}", color=STY[code]["color"], lw=1.2)
ax[0].set_title("onset time (yr) against event number"); ax[0].set_ylabel("yr")
ax[1].set_title("onset: ours minus each reference (yr)"); ax[1].axhline(0, color="k", lw=0.6)
ax[2].set_title("interval before the event: ours minus each reference (yr)"); ax[2].axhline(0, color="k", lw=0.6)
for a in ax:
    a.set_xlabel("event number"); a.grid(alpha=0.3); a.set_xticks(k1); a.legend(fontsize=8)
fig.suptitle(f"BP3-QD dip 60 thrust, events 1 to {N}: mean recurrence from ev 2 / ev 3 (yr): "
             + ", ".join(f"{c} {REC[c][0]:.2f}/{REC[c][1]:.2f}" for c, _ in CODES), fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fn = f"p{n:02d}_summary_onsets.png"; fig.savefig(os.path.join(OUT, fn), dpi=110); plt.close(fig); pages.append(fn); n += 1


# ---------------------------------------------------------------- largest differences
def front(d, on, thr=-3.0):
    rel = d["t"] - on; m = (rel > W0) & (rel < W1); i = np.nonzero((d["logV"] > thr) & m)[0]
    if not i.size or i[0] == 0:
        return np.nan
    k = i[0]; t0, t1, v0, v1 = rel[k - 1], rel[k], d["logV"][k - 1], d["logV"][k]
    return t0 + (thr - v0) * (t1 - t0) / (v1 - v0)


LAST = N - 1
diffs = []
for ref in (TP1,):
    for st in FAULT:
        a, b = F[OURS][st], F[ref][st]
        ea, eb = EV[OURS][LAST]["onset"], EV[ref][LAST]["onset"]
        fa, fb = front(a, ea), front(b, eb)
        for q, unit in (("slip", "m"), ("logV", "log10 m/s"), ("tau", "MPa")):
            ts, d = diff_series(a["t"] - ea, a[q], b["t"] - eb, b[q], W0, W1)
            keep = np.ones(ts.size, bool)
            for f in (fa, fb):
                if np.isfinite(f):
                    keep &= np.abs(ts - f) >= 2.0
            if keep.any():
                j = np.argmax(np.abs(d[keep])); diffs.append((ref, st, q, f"ev {N} window, off front", d[keep][j], unit, ts[keep][j], "s"))
            if ref == TP1:
                # interseismic: away from every event of either code, onset - 1 day to onset + 30 days
                ts, d = diff_series(a["t"], a[q], b["t"], b[q], 0, T_CUT)
                keep = np.ones(ts.size, bool)
                for c in (OURS, ref):
                    for e in EV[c]:
                        keep &= ~((ts >= e["onset"] - 86400) & (ts <= e["onset"] + 30 * 86400))
                if keep.any():
                    j = np.argmax(np.abs(d[keep])); diffs.append((ref, st, q, "0 to cut, interseismic", d[keep][j], unit, ts[keep][j] / YR, "yr"))
with open(os.path.join(OUT, "largest_differences.txt"), "w") as f:
    f.write("ours minus reference; the denser series interpolated onto the sparser one's times\n")
    f.write(f"ev {N} window: -100 to +300 s from each code's own onset, excluding +-2 s around either code's front "
            "(log10 V crossing -3) at that station\n")
    f.write("interseismic: 0 to the cut, excluding onset - 1 day to onset + 30 days of every event of either code\n")
    for row in sorted(diffs, key=lambda r: (r[0], r[3], r[2], -abs(r[4]))):
        ref, st, q, win, d, unit, t, tu = row
        f.write(f"{ref:18s} {win:26s} {q:5s} {st:6s} {d:+12.5g} {unit:10s} at t = {t:.6g} {tu}\n")
FR = []
for st in FAULT:
    a = front(F[OURS][st], EV[OURS][LAST]["onset"]); b = front(F[TP1][st], EV[TP1][LAST]["onset"])
    FR.append(f"{st:6s} {a:8.3f} {b:8.3f} {a - b:+8.3f}")
open(os.path.join(OUT, "front_offsets.txt"), "w").write(
    f"event {N}: front arrival (log10 V crosses -3), s after each code's own onset\n"
    f"{'st':6s} {'ours':>8s} {'T p1':>8s} {'d(p1)':>8s}\n" + "\n".join(FR) + "\n")

# ---------------------------------------------------------------- (c) fault station pages
# Rows: 0 to the cut; ours minus Tandem p = 1 there (between events: onset - 1 day to onset + 30 days
# of every event of either code left out); the event-N window; ours minus Tandem p = 1 there.
DCOL = "#d62728"


def interseismic(ts):
    keep = np.ones(ts.size, bool)
    for c in (OURS, TP1):
        for e in EV[c]:
            keep &= ~((ts >= e["onset"] - 86400) & (ts <= e["onset"] + 30 * 86400))
    return keep


def diff_panels(axl, axe, A, B, qs_long, qs_ev):
    """ours minus Tandem p = 1, the locally denser series interpolated point by point."""
    oa, ob = EV[OURS][LAST]["onset"], EV[TP1][LAST]["onset"]
    for a, q in zip(axl, qs_long):
        ts, d = diff_series(A["t"], A[q], B["t"], B[q], 0, T_CUT); k = interseismic(ts)
        a.plot(ts[k] / YR, d[k], ".", ms=1.5, color=DCOL); a.axhline(0, color="k", lw=0.5)
    for a, q in zip(axe, qs_ev):
        ts, d = diff_series(A["t"] - oa, A[q], B["t"] - ob, B[q], W0, W1)
        a.plot(ts, d, ".-", ms=2, lw=0.6, color=DCOL); a.axhline(0, color="k", lw=0.5); a.set_xlim(W0, W1)


for st in FAULT:
    fig, ax = plt.subplots(4, 3, figsize=(15, 15))
    for code, _ in CODES:
        d, on = F[code][st], EV[code][LAST]["onset"]
        for a, q in zip(ax[0], ("slip", "logV", "tau")):
            a.plot(d["t"] / YR, d[q], label=code, **STY[code])
        rel = d["t"] - on; m = (rel >= W0) & (rel <= W1)
        for a, q in zip(ax[2], ("logV", "tau", "slip")):
            a.plot(rel[m], d[q][m], label=code, **STY[code])
    diff_panels(ax[1], ax[3], F[OURS][st], F[TP1][st], ("slip", "logV", "tau"), ("logV", "tau", "slip"))
    for a, t in zip(ax[0], ("slip (m)", "log10 slip rate (m/s)", "shear stress (MPa)")):
        a.set_title(t + f", 0 to {T_CUT / YR:.1f} yr", fontsize=10); a.set_xlabel("time (yr)")
    for a, t in zip(ax[1], ("slip (m)", "log10 slip rate", "shear stress (MPa)")):
        a.set_title("ours - Tandem p=1, " + t + ", between events", fontsize=10); a.set_xlabel("time (yr)")
    for a, t in zip(ax[2], ("log10 slip rate (m/s)", "shear stress (MPa)", "slip (m)")):
        a.set_title(t + f", event {N}", fontsize=10); a.set_xlabel("time from each code's own onset (s)"); a.set_xlim(W0, W1)
    for a, t in zip(ax[3], ("log10 slip rate", "shear stress (MPa)", "slip (m)")):
        a.set_title(f"ours - Tandem p=1, " + t + f", event {N}", fontsize=10); a.set_xlabel("time from each code's own onset (s)")
    for a in ax.flat:
        a.grid(alpha=0.3)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=10)
    fig.suptitle(f"BP3-QD dip 60 thrust, fault station {st} ({int(st[2:]) / 10:g} km down-dip): seas-mfem p=1 against Tandem p=1\n"
                 f"event {N} onsets (yr): " + ", ".join(f"{c} {EV[c][LAST]['onset'] / YR:.4f}" for c, _ in CODES)
                 + "; between events = onset - 1 day to onset + 30 days left out", fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fn = f"p{n:02d}_fltst_{st}.png"; fig.savefig(os.path.join(OUT, fn), dpi=110); plt.close(fig); pages.append(fn); n += 1

# ---------------------------------------------------------------- (d) surface station pages
for st in SURF:
    fig, ax = plt.subplots(4, 4, figsize=(18, 15))
    for code, _ in CODES:
        d = S[code][st]
        on = EV[code][LAST]["onset"]; rel = d["t"] - on; m = (rel >= W0) & (rel <= W1)
        lab = code + (" (du/dt)" if d["vfd"] else "")
        for j, q in enumerate(("d1", "d2", "v1", "v2")):
            ax[0, j].plot(d["t"] / YR, d[q], label=lab, **STY[code])
            ax[2, j].plot(rel[m], d[q][m], label=lab, **STY[code])
    diff_panels(ax[1], ax[3], S[OURS][st], S[TP1][st], ("d1", "d2", "v1", "v2"), ("d1", "d2", "v1", "v2"))
    names = ("displacement 1 (m)", "displacement 2, down (m)", "velocity 1 (m/s)", "velocity 2, down (m/s)")
    for j, t in enumerate(names):
        ax[0, j].set_title(t + f", 0 to {T_CUT / YR:.1f} yr", fontsize=10); ax[0, j].set_xlabel("time (yr)")
        ax[1, j].set_title("ours - Tandem p=1, " + t + ", between events", fontsize=9); ax[1, j].set_xlabel("time (yr)")
        ax[2, j].set_title(t + f", event {N}", fontsize=10); ax[2, j].set_xlabel("time from each code's own onset (s)")
        ax[2, j].set_xlim(W0, W1)
        ax[3, j].set_title("ours - Tandem p=1, " + t + f", event {N}", fontsize=9); ax[3, j].set_xlabel("time from each code's own onset (s)")
    for a in ax.flat:
        a.grid(alpha=0.3)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=10)
    x = {"fn-32": -32, "fn-16": -16, "fn-08": -8, "fn-00": "0-", "fn+00": "0+", "fn+08": 8, "fn+16": 16, "fn+32": 32}[st]
    fig.suptitle(f"BP3-QD dip 60 thrust, surface station {st} (x1 = {x} km): seas-mfem p=1 against Tandem p=1 "
                 "(Tandem velocity = du/dt; between events = onset - 1 day to onset + 30 days left out)", fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fn = f"p{n:02d}_srfst_{st}.png"; fig.savefig(os.path.join(OUT, fn), dpi=110); plt.close(fig); pages.append(fn); n += 1

# ---------------------------------------------------------------- (e) notes page
PK = {}
for c, _ in CODES:
    v = [EV[c][k]["peak"] for k in range(1, N)]; PK[c] = (min(v), max(v))
rs = []
for k in range(1, N):
    a, b = F[OURS]["dp000"], F[TP1]["dp000"]; oa, ob = EV[OURS][k]["onset"], EV[TP1][k]["onset"]
    m = (a["t"] >= oa - 1) & (a["t"] <= oa + 2000)
    rs.append(np.max(10 ** np.interp(a["t"][m] - oa + ob, b["t"], b["logV"])))
PK["rs"] = (min(rs), max(rs))
D150 = {}
for q, key in (("tau", "tau"), ("logV", "v")):
    ts, d = diff_series(F[OURS]["dp150"]["t"], F[OURS]["dp150"][q], F[TP1]["dp150"]["t"], F[TP1]["dp150"][q], 0, T_CUT)
    k = interseismic(ts); e = ts < 100 * YR
    D150[key + "0"] = d[k & e][np.argmax(np.abs(d[k & e]))]; D150[key + "1"] = d[k & ~e][np.argmax(np.abs(d[k & ~e]))]
nb = []
for st in ("dp125", "dp175"):
    ts, d = diff_series(F[OURS][st]["t"], F[OURS][st]["tau"], F[TP1][st]["t"], F[TP1][st]["tau"], 0, 100 * YR)
    nb.append(np.max(np.abs(d[interseismic(ts)])))
D150["nb"] = max(nb)
gap = np.median(np.diff(F[OURS]["dp025"]["t"][(F[OURS]["dp025"]["t"] > EV[OURS][LAST]["end"] + 86400) &   # past the 0.1 s tail rows
                                              (F[OURS]["dp025"]["t"] < EV[OURS][LAST]["end"] + 2 * YR)]))
NOTES = [
    "Notes",
    "",
    f"1. Cut.  Ours ends at {T_END[OURS] / YR:.4f} yr and Tandem p = 1 at {T_END[TP1] / YR:.4f} yr; both are cut at the earlier,",
    f"   {T_CUT / YR:.4f} yr.  Each has {N} events, the last starting {(T_CUT - max(EV[c][LAST]['onset'] for c in EV)) / YR:.1f} yr before the cut.",
    "   Ours = the run's SCEC export (cycle/scec_bp3_export.py).",
    "   Onsets of both codes come from the station data with one definition (first time max V over the 12 fault stations",
    "   > 1e-3 m/s after more than 1 yr below).  Events are matched by index.",
    "",
    "2. The post-event output gap.  Ours writes a row every 0.1 s while max V > 1e-3 m/s, then, 10 s after it drops below,",
    f"   every dt_slow = 0.1 yr (median row spacing 1 day to 2 yr after event {N}: {gap / 86400:.1f} days).  The afterslip in the first",
    "   hours to days after each event is therefore not sampled: off-rupture slip-rate peaks that the references show",
    "   15 to 30 minutes after the event (dp025, dp175, dp200) fall into that gap.  Up to the end of our dense rows our",
    "   slip rate matches Tandem p = 1.  It is an output-schedule limit, not a solution difference.",
    "   The interseismic differences below therefore exclude onset - 1 day to onset + 30 days of every event.",
    "",
    "3. The near-trace surface pair fn-00 / fn+00.  Ours evaluates both stations exactly on the fault trace (x1 = 0), as",
    "   the limit from each side, so they differ only by the trace slip (dp000, ~1e-5 m).  Tandem's sit at x1 = -1 m and",
    "   +1 m, so they also carry the displacement gradient across 2 m of the hanging and foot walls; the pair",
    "   difference is larger there.  Neither is an error; the station definitions differ.",
    "",
    "4. Peak magnitudes differ because of OUTPUT SAMPLING, adaptive (Tandem) against fixed (ours), not the solution.",
    "   Tandem writes a row when the fault-wide max V has changed by 1 % since its last row (rtol 0.01), or after 0.3 yr;",
    "   ours writes every 0.1 s while max V > 1e-3 m/s, else every 0.1 yr.  So Tandem writes every few ms where max V",
    "   changes fast and nothing for up to about 2 s on the mid-rupture plateau.  From event 2 on the rupture breaks the",
    f"   surface and the peak is at dp000, a spike of about 0.05 s: Tandem catches {PK[TP1][0]:.1f} to {PK[TP1][1]:.1f} m/s, ours whichever 0.1 s row",
    f"   lands in it ({PK[OURS][0]:.1f} to {PK[OURS][1]:.1f} m/s); Tandem on our rows (each code's own onset as zero) gives {PK['rs'][0]:.1f} to {PK['rs'][1]:.1f} m/s.",
    "   At dp050, events 3 on, the front passes inside Tandem's mid-rupture gap: its file misses the arrival and the peak.",
    "   Use the slip and onset columns, not peak V.",
    "",
    "5. dp150 (15 km down-dip, where a starts to rise) differs from Tandem p = 1 before event 1, between about 30 and",
    f"   70 yr: shear stress up to {D150['tau0']:+.2f} MPa and log10 V up to {D150['v0']:+.2f} (ours minus Tandem, 0 to 100 yr).  After 100 yr, between",
    f"   events, the largest is {D150['tau1']:+.4f} MPa and {D150['v1']:+.4f} in log10 V.  dp125 and dp175 agree to {D150['nb']:.0e} MPa over 0 to 100 yr.",
    "   Both codes write densely there, so this is not sampling.  The probe sits on the mesh vertex at the kink of a",
    "   (x_d = 15 km), where two fault faces meet and each keeps its own value; the codes likely read different faces",
    "   (ours the up-dip one).  Not confirmed.",
    "",
    "6. Tandem p = 1 velocity at the surface stations is du/dt by finite differences (its output has no velocity).",
    "",
    "7. Differences are ours minus reference, the locally denser series interpolated onto the sparser one's times.",
    "   Wherever one series is locally denser, it is interpolated onto the other one's samples (decided sample by sample).",
]
fig = plt.figure(figsize=(15, 8.5))
fig.text(0.02, 0.97, "\n".join(NOTES), family="monospace", fontsize=10, va="top")
fn = f"p{n:02d}_notes.png"; fig.savefig(os.path.join(OUT, fn), dpi=110); plt.close(fig); pages.append(fn); n += 1

ims = [Image.open(os.path.join(OUT, p)).convert("RGB") for p in pages]
ims[0].save(os.path.join(OUT, "cycle_compare_ours_vs_tandem_p1_1500yr.pdf"), save_all=True, append_images=ims[1:], resolution=110)
open(os.path.join(OUT, "page_order.txt"), "w").write("\n".join(pages) + "\n")
print("\n".join(L))
print(f"\nT_CUT = {T_CUT / YR:.4f} yr; {len(pages)} pages -> {os.path.join(OUT, 'cycle_compare_ours_vs_tandem_p1_1500yr.pdf')}")
