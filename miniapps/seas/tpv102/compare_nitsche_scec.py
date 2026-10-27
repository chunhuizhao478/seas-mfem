#!/usr/bin/env python3
"""TPV102: compare Nitsche station output with the SCEC references and gold.

Plan N-07 / P7 acceptance (cross-method bands, information, not pass/fail):
  rupture arrival (first h-slip-rate > 1 mm/s)  2 %
  peak slip rate (max |h-slip-rate|)            10 %
  final slip at t_end (h-slip)                  3 %
  log10 theta at t_end                          0.1 decades
  sigma_n (column 8) against BOTH gold mixedflux and upwind
  sign agreement of every signed column against the SCEC files

Inputs
  --ours DIR       our station files  <prefix>_station_<name>.dat
                   (t slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta, SI;
                   1 = dip = SCEC vertical, 2 = strike = SCEC horizontal,
                   tau = SCEC shear stress -T.t, sigma_n compression +)
  --prefix NAME    station-file prefix (default tpv102)
  --bench DIR      tpv102/benchmark_data (scec_drdg3d/, scec_pylith/)
  --gold DIR ...   gold result dirs (same layout as ours)

The SCEC files carry  t h-slip h-slip-rate h-shear-stress v-slip
v-slip-rate v-shear-stress n-stress log-theta  (MPa; DRDG3D writes
n-stress = -120, so sigma_n is compared by magnitude).  The gold files were
written by the velocity-stress path; its stored slip rate is parallel to the
stored total traction (dynamic/fault_face_flux.cpp:260-262) and the stored
traction is the total tau_0 + tau_nuc + correction (:313-319) with tau2_0 =
+75 MPa, i.e. the SCEC signs: the C6 conversion is the identity on these
columns (read from the routine, confirmed on the data at t = 0 and at the
peak).  Scope: Chunhui Zhao, 2026-10-14 (P7).
"""
import argparse, glob, os, sys
import numpy as np

STATIONS = ["flt_0_3", "flt_0_7.5", "flt_0_12", "flt_9_7.5", "flt_12_3",
            "flt_12_12", "flt_n9_7.5", "flt_n12_3", "flt_n12_12"]
# SCEC file names: x2 = along strike (km), x3 = depth (km)
SCEC_NAME = {"flt_0_3": "x2_0_x3_3", "flt_0_7.5": "x2_0_x3_7.5", "flt_0_12": "x2_0_x3_12",
             "flt_9_7.5": "x2_9_x3_7.5", "flt_12_3": "x2_12_x3_3", "flt_12_12": "x2_12_x3_12",
             "flt_n9_7.5": "x2_-9_x3_7.5", "flt_n12_3": "x2_-12_x3_3",
             "flt_n12_12": "x2_-12_x3_12"}
BANDS = {"arrival": 0.02, "peak": 0.10, "slip": 0.03, "logtheta": 0.1}


def load_ours(path):
    """-> dict of SCEC-named columns in SCEC units (m, m/s, MPa)."""
    A = np.loadtxt(path, comments="#", ndmin=2)
    return {"t": A[:, 0], "h-slip": A[:, 2], "h-slip-rate": A[:, 4],
            "h-shear-stress": A[:, 6] / 1e6, "v-slip": A[:, 1], "v-slip-rate": A[:, 3],
            "v-shear-stress": A[:, 5] / 1e6, "n-stress": A[:, 7] / 1e6, "log-theta": A[:, 8]}


def load_scec(path):
    rows = []
    for line in open(path):
        w = line.split()
        if not w or line.startswith("#") or w[0][0] not in "0123456789-.":
            continue
        try:
            rows.append([float(v) for v in w[:9]])
        except ValueError:
            continue
    A = np.array(rows)
    keys = ["t", "h-slip", "h-slip-rate", "h-shear-stress", "v-slip", "v-slip-rate",
            "v-shear-stress", "n-stress", "log-theta"]
    # PyLith's header says "up-dip", DRDG3D's says "vertical"; the raw values
    # agree (e.g. v-slip at 12 s, x2=12 x3=12: 0.0258 vs 0.0264), so both
    # files share one convention and neither is negated.
    return {k: A[:, i] for i, k in enumerate(keys)}


def metrics(d, t_end):
    t = d["t"]
    v = d["h-slip-rate"]
    above = np.nonzero(v > 1e-3)[0]
    arr = t[above[0]] if above.size else np.nan
    m = t <= t_end + 1e-9
    peak = np.max(np.abs(v[m])) if m.any() else np.nan
    tpk = t[m][np.argmax(np.abs(v[m]))] if m.any() else np.nan
    slip = np.interp(t_end, t, d["h-slip"])
    lth = np.interp(t_end, t, d["log-theta"])
    return arr, peak, slip, lth, tpk


def sign_agreement(ours, ref, key, frac_floor=0.01):
    """Fraction of common times where both are significant and share a sign.
    'significant' = |value| > frac_floor * max|ref|."""
    t = ours["t"]
    r = np.interp(t, ref["t"], ref[key])
    o = ours[key]
    scale = np.max(np.abs(ref[key])) if ref[key].size else 0.0
    if scale == 0.0:
        return np.nan, 0
    sig = (np.abs(r) > frac_floor * scale) & (np.abs(o) > frac_floor * scale)
    n = int(sig.sum())
    if n == 0:
        return np.nan, 0
    return float(np.mean(np.sign(r[sig]) == np.sign(o[sig]))), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--prefix", default="tpv102")
    ap.add_argument("--bench", default=os.path.join(os.path.dirname(__file__), "benchmark_data"))
    ap.add_argument("--gold", nargs="*", default=[])
    a = ap.parse_args()

    ours = {}
    for s in STATIONS:
        p = os.path.join(a.ours, f"{a.prefix}_station_{s}.dat")
        if os.path.exists(p):
            ours[s] = load_ours(p)
    if not ours:
        sys.exit(f"no station files under {a.ours}")
    t_end = min(min(d["t"][-1] for d in ours.values()), 12.0)
    print(f"# t_end = {t_end:.4f} s (min of our last time and 12 s)")
    refs = {}
    for code in ("drdg3d", "pylith"):
        for s in STATIONS:
            f = os.path.join(a.bench, f"scec_{code}", f"tpv102_{code}_{SCEC_NAME[s]}.txt")
            if os.path.exists(f):
                refs[(code, s)] = load_scec(f)
    golds = {}
    for g in a.gold:
        for s in STATIONS:
            f = os.path.join(g, f"tpv102_station_{s}.dat")
            if os.path.exists(f):
                gname = os.path.basename(os.path.dirname(g.rstrip("/"))) \
                    if os.path.basename(g.rstrip("/")) == "results" else os.path.basename(g.rstrip("/"))
                golds[(gname.split("_")[0] + "_" + gname.split("_")[1], s)] = load_ours(f)

    hdr = (f"{'station':11s} {'src':10s} {'arrival':>9s} {'peak V':>8s} {'slip':>8s} {'log th':>7s}"
           f" {'t_peak':>7s}")
    print(hdr)
    for s in STATIONS:
        if s not in ours:
            print(f"{s:11s} (missing)")
            continue
        mo = metrics(ours[s], t_end)
        print(f"{s:11s} {'nitsche':10s} {mo[0]:9.4f} {mo[1]:8.3f} {mo[2]:8.4f} {mo[3]:7.3f} {mo[4]:7.3f}")
        for code in ("drdg3d", "pylith"):
            if (code, s) not in refs:
                continue
            mr = metrics(refs[(code, s)], t_end)
            rel = [(mo[0] - mr[0]) / mr[0], (mo[1] - mr[1]) / mr[1], (mo[2] - mr[2]) / mr[2]]
            dl = mo[3] - mr[3]
            flags = ["" if abs(rel[0]) <= BANDS["arrival"] else "*",
                     "" if abs(rel[1]) <= BANDS["peak"] else "*",
                     "" if abs(rel[2]) <= BANDS["slip"] else "*",
                     "" if abs(dl) <= BANDS["logtheta"] else "*"]
            print(f"{'':11s} {code:10s} {mr[0]:9.4f} {mr[1]:8.3f} {mr[2]:8.4f} {mr[3]:7.3f} {mr[4]:7.3f}"
                  f"   d: {100*rel[0]:+6.2f}%{flags[0]:1s} {100*rel[1]:+6.1f}%{flags[1]:1s}"
                  f" {100*rel[2]:+6.2f}%{flags[2]:1s} {dl:+6.3f}{flags[3]:1s}")
        for (gname, gs), gd in golds.items():
            if gs != s:
                continue
            mg = metrics(gd, t_end)
            print(f"{'':11s} {gname[:10]:10s} {mg[0]:9.4f} {mg[1]:8.3f} {mg[2]:8.4f} {mg[3]:7.3f} {mg[4]:7.3f}")

    print("\n# sign agreement against the SCEC files (fraction of significant samples)")
    print(f"{'station':11s} {'ref':7s} " + " ".join(f"{k:>15s}" for k in
          ("h-slip", "h-slip-rate", "h-shear-stress", "v-slip", "v-slip-rate", "v-shear-stress")))
    for s in STATIONS:
        if s not in ours:
            continue
        for code in ("drdg3d", "pylith"):
            if (code, s) not in refs:
                continue
            cells = []
            for k in ("h-slip", "h-slip-rate", "h-shear-stress", "v-slip", "v-slip-rate",
                      "v-shear-stress"):
                fr, n = sign_agreement(ours[s], refs[(code, s)], k)
                cells.append(f"{'n/a':>15s}" if n == 0 else f"{fr:9.4f} ({n:4d})")
            print(f"{s:11s} {code:7s} " + " ".join(cells))

    print("\n# sigma_n: max |sigma_n - 120 MPa| over the run (MPa)")
    for s in STATIONS:
        if s not in ours:
            continue
        nz = ours[s]["n-stress"] != 0.0   # some files write an all-zero row at t = 0
        line = f"{s:11s} nitsche {np.max(np.abs(np.abs(ours[s]['n-stress'][nz]) - 120.0)):9.4f}"
        for (gname, gs), gd in golds.items():
            if gs == s:
                m = (gd["t"] <= t_end + 1e-9) & (gd["n-stress"] != 0.0)
                line += f"  {gname[:18]} {np.max(np.abs(np.abs(gd['n-stress'][m]) - 120.0)):9.4f}"
        for code in ("drdg3d", "pylith"):
            if (code, s) in refs:
                r = refs[(code, s)]
                m = (r["t"] <= t_end + 1e-9) & (r["n-stress"] != 0.0)
                line += f"  {code} {np.max(np.abs(np.abs(r['n-stress'][m]) - 120.0)):9.4f}"
        print(line)
    print("\n# '*' = outside the band (arrival 2 %, peak 10 %, slip 3 %, log theta 0.1): "
          "a finding to explain, not an automatic failure")


if __name__ == "__main__":
    main()
