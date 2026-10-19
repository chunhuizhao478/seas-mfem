#!/usr/bin/env python3
"""Compare Nitsche fault-station output with SCEC references and golds.

Works for any problem whose stations are listed in a station file written
by scripts/make_fault_stations.py (TPV102, TPV205, TPV104, ...).  It
generalizes tpv102/compare_nitsche_scec.py: the station list, the
reference directories and the file formats are arguments, not constants.

Sources (each LABEL=DIR, repeatable; the order given is kept):
  --run  LABEL=DIR          Nitsche driver output, <prefix>_station_<name>.dat
                            columns  t slip1 slip2 V1 V2 tau1 tau2 sigma_n state
                            (SI; 1 = dip, 2 = strike; state = log10 theta,
                            psi or slip path by law).  Give runs in
                            refinement order for the convergence check.
  --ref  LABEL=DIR          SCEC reference files (MPa):  t h-slip h-slip-rate
                            h-shear-stress v-slip v-slip-rate v-shear-stress
                            [n-stress] [state].  The first --ref is the one
                            the pass test uses.
  --gold LABEL=DIR[:FMT]    golds; FMT = scec_si (SAFS writer, Pa; default)
                            or nitsche (TPV102 writer layout).  A gold's
                            station files are in its results/ subfolder.
A station's file in DIR is the one whose stem ends in "_station_<name>" or
in "_<key>", key = x2_<strike km>_x3_<depth km> from the station file.

Per-side mode (--per-side; TPV6/TPV7, plan 2026-10-19 section 9): every
source gives two files per station, one for each side of the fault, with the
DRDG3D split-node columns
    t h-disp h-vel h-stress v-disp v-vel v-stress n-disp n-vel n-stress
(FMT side_si: stresses in Pa, the driver's [output.stations] per_side files
<prefix>_station_<name>_{farside,nearside}.dat, the default for --run;
side_mpa: stresses in MPa, the DRDG3D files *_{farside,nearside}_<key>.txt,
the default for --ref and --gold).  A side's file is the one whose stem ends
in "_station_<name>_<side>" or in "_<side>_<key>".  Slip and slip rate are
far minus near (SCEC h-slip > 0 right-lateral), the stresses are the far
file's (the traction is single-valued).

Measures per station (t <= t_end):
  arrival   first time |V| > 1 mm/s (SCEC definition), |V| = |(V_h, V_v)|
  peak      max |V|
  slip      |(slip_h, slip_v)| at t_end
            (t_end = --t-end, lowered to the earliest end of any source,
            with a warning, so every source is measured on one window)
  state     state column at t_end (when both sides have it)
  noise     RMS of sigma_n, tau_h, tau_v as in the TPV102 and damping-test
            reports (plan section 8): from the station's arrival (first
            |V_h| > 1 mm/s) + --noise-after s (0.5) to t_end, the residual
            after a centred moving average of --noise-window s (0.1) in
            samples of the source's own median dt, w samples dropped at each
            end; |sigma_n| for the normal stress

Pass test (plan rev 1.3, section 8; accepted by Chunhui Zhao 2026-10-16),
against the first --ref, for the LAST --run:
  1 arrival within --arrival-tol s (0.1) at every station
  2 final slip within --slip-tol (2 %) where the reference slips > 0.1 m
  3 peak slip rate within --peak-tol (10 %) at >= --peak-min-stations
  4 the median |misfit| of 1-3 decreases along the --run order
  (5, no abort, is read from the driver log, not here)
  6 per-side mode only: each side's final h-disp within --side-tol (3 %) at
    every station where the reference side moves more than --slip-min; the
    far/near ratio of final h-disp is printed for every source
With two or more --ref, the spread between the references is printed too.

Usage (from miniapps/seas):
  python3 scripts/compare_nitsche_fault_stations.py \
      --stations tpv205/stations/tpv205_fault_stations.csv --prefix tpv205 \
      --run p1_200m=OUT1 --run p2_200m=OUT2 \
      --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4 \
      --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 \
      --gold gold_p2=tpv205/gold/mixedflux_p2_rk45_jobid7762880_git04674fc/results \
      --peak-min-stations 14 --json summary.json
  TPV6 (per side, 5 stations, peaks at 4 of 5):
  python3 scripts/compare_nitsche_fault_stations.py --per-side \
      --stations tpv6/stations/tpv6_fault_stations.csv --prefix tpv6 \
      --run p1_100m=OUT1 --run p2_100m=OUT2 --run p3_100m=OUT3 \
      --ref drdg3d=tpv6/benchmark_data/scec_drdg3d --peak-min-stations 4
  TPV7: the same with --prefix tpv7, --stations tpv7/stations/tpv7_fault_stations.csv
  and --ref drdg3d=tpv7/benchmark_data/scec_drdg3d.
  --list-columns prints the column map of every file found and stops.

Scope: Chunhui Zhao, 2026-10-16 (TPV205/TPV104 plan rev 1.3, phase B0);
per-side mode 2026-10-19 (TPV6/TPV7/TPV31 plan rev 1, item A5).
"""
import argparse
import csv
import json
import os
import re
import sys

import numpy as np

KEYS = ["t", "h-slip", "h-slip-rate", "h-shear-stress", "v-slip", "v-slip-rate",
        "v-shear-stress", "n-stress", "state"]
SIDE_KEYS = ["t", "h-disp", "h-vel", "h-stress", "v-disp", "v-vel", "v-stress",
             "n-disp", "n-vel", "n-stress"]
SIDES = ("farside", "nearside")


# ------------------------------------------------------------------ stations
def read_stations(path):
    rows = [l for l in open(path) if l.strip() and not l.lstrip().startswith("#")]
    out = []
    for r in csv.DictReader(rows):
        x, y, z = float(r["x"]), float(r["y"]), float(r["z"])
        key = f"x2_{x / 1000.0:g}_x3_{-z / 1000.0 + 0.0:g}"
        out.append({"name": r["name"].strip(), "x": x, "y": y, "z": z, "key": key})
    if not out:
        sys.exit(f"no stations in {path}")
    return out


def find_file(d, st, side=None, file_prefix=None):
    """side None: the station's one file; side 'farside'/'nearside': that
    side's file.  file_prefix, when given, must start the file name."""
    tail = "" if side is None else "_" + side
    hits = []
    for f in sorted(os.listdir(d)):
        stem, ext = os.path.splitext(f)
        if ext not in (".dat", ".txt"):
            continue
        if file_prefix and not f.startswith(file_prefix):
            continue
        if stem.endswith("_station_" + st["name"] + tail) or stem.endswith(tail + "_" + st["key"]):
            hits.append(os.path.join(d, f))
    if len(hits) > 1:
        sys.exit(f"{d}: more than one file for station {st['name']}: {hits}")
    return hits[0] if hits else None


# ------------------------------------------------------------------- loaders
def _numeric_rows(path, ncol_max=10):
    rows = []
    for line in open(path):
        w = line.split()
        if not w or line.lstrip().startswith("#") or w[0][0] not in "0123456789-.+":
            continue
        try:
            rows.append([float(v) for v in w[:ncol_max]])
        except ValueError:
            continue
    if not rows:
        sys.exit(f"{path}: no numeric rows")
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def load(path, fmt):
    """-> dict KEYS -> array (m, m/s, MPa); missing columns are absent."""
    A = _numeric_rows(path)
    nc = A.shape[1]
    d = {}
    if fmt == "nitsche":
        if nc < 8:
            sys.exit(f"{path}: {nc} columns, Nitsche layout needs >= 8")
        d = {"t": A[:, 0], "h-slip": A[:, 2], "h-slip-rate": A[:, 4],
             "h-shear-stress": A[:, 6] / 1e6, "v-slip": A[:, 1], "v-slip-rate": A[:, 3],
             "v-shear-stress": A[:, 5] / 1e6, "n-stress": A[:, 7] / 1e6}
        if nc >= 9:
            d["state"] = A[:, 8]
    elif fmt in ("scec_si", "scec_mpa"):
        if nc < 7:
            sys.exit(f"{path}: {nc} columns, SCEC layout needs >= 7")
        sc = 1e-6 if fmt == "scec_si" else 1.0
        for i, k in enumerate(KEYS[:nc]):
            d[k] = A[:, i] * (sc if k.endswith("stress") else 1.0)
    else:
        sys.exit(f"unknown format {fmt}")
    if "n-stress" in d:
        d["n-stress"] = np.abs(d["n-stress"])      # DRDG3D TPV102 writes -120
    d["_ncol"] = nc
    return d


def load_side(path, fmt):
    """One side's split-node file -> dict SIDE_KEYS -> array (m, m/s, MPa)."""
    A = _numeric_rows(path)
    if A.shape[1] < 10:
        sys.exit(f"{path}: {A.shape[1]} columns, the per-side layout needs 10")
    if fmt not in ("side_si", "side_mpa"):
        sys.exit(f"{path}: format {fmt} is not a per-side format (side_si, side_mpa)")
    sc = 1e-6 if fmt == "side_si" else 1.0
    return {k: A[:, i] * (sc if k.endswith("stress") else 1.0) for i, k in enumerate(SIDE_KEYS)}


def combine_sides(far, near, path_far):
    """Far and near files of one station -> the KEYS dict of load(): slip and
    slip rate far minus near, stresses from the far file; plus the final
    per-side h-disp series for pass item 6."""
    t = far["t"]
    if near["t"].shape != t.shape or np.max(np.abs(near["t"] - t)) > 1e-9 * max(1.0, t[-1]):
        sys.exit(f"{path_far}: the near-side file has other sample times")
    d = {"t": t,
         "h-slip": far["h-disp"] - near["h-disp"], "h-slip-rate": far["h-vel"] - near["h-vel"],
         "h-shear-stress": far["h-stress"],
         "v-slip": far["v-disp"] - near["v-disp"], "v-slip-rate": far["v-vel"] - near["v-vel"],
         "v-shear-stress": far["v-stress"], "n-stress": np.abs(far["n-stress"]),
         "_h_far": far["h-disp"], "_h_near": near["h-disp"], "_ncol": 10}
    return d


# ------------------------------------------------------------------ measures
def noise_rms(d, k, window, after):
    """Noise of column k of d (already cut at t_end): the measure of the
    TPV102 and damping-test reports (rms() in their build scripts).  The
    arrival uses |V_h| where those scripts used V_h > 1 mm/s; the two agree
    for every benchmark here, whose ruptures slip in +h."""
    t = d["t"]
    i = np.nonzero(np.abs(d["h-slip-rate"]) > 1e-3)[0]
    if i.size == 0:
        return float("nan")
    m = t >= t[i[0]] + after
    if m.sum() < 3:
        return float("nan")
    w = max(1, int(round(window / np.median(np.diff(t[m])))))
    x = np.abs(d[k][m]) if k == "n-stress" else d[k][m]
    if x.size <= 2 * w:
        return float("nan")
    r = (x - np.convolve(x, np.ones(w) / w, mode="same"))[w:-w]
    return float(np.sqrt(np.mean(r * r)))


def measures(d, t_end, window, after):
    t = d["t"]
    V = np.hypot(d["h-slip-rate"], d["v-slip-rate"])
    m = t <= t_end + 1e-12
    above = np.nonzero((V > 1e-3) & m)[0]
    out = {
        "arrival": float(t[above[0]]) if above.size else float("nan"),
        "peak": float(np.max(V[m])) if m.any() else float("nan"),
        "t_peak": float(t[m][np.argmax(V[m])]) if m.any() else float("nan"),
        "slip": float(np.hypot(np.interp(t_end, t, d["h-slip"]), np.interp(t_end, t, d["v-slip"]))),
    }
    if "state" in d:
        out["state"] = float(np.interp(t_end, t, d["state"]))
    if "_h_far" in d:
        out["h_far"] = float(np.interp(t_end, t, d["_h_far"]))
        out["h_near"] = float(np.interp(t_end, t, d["_h_near"]))
    dc = {k: v[m] for k, v in d.items() if isinstance(v, np.ndarray)}
    for k, lab in (("n-stress", "noise_sn"), ("h-shear-stress", "noise_th"), ("v-shear-stress", "noise_tv")):
        if k in d:
            out[lab] = noise_rms(dc, k, window, after)
    return out


def misfit(mo, mr):
    def rel(a, b):
        return (a - b) / b if b not in (0.0,) and np.isfinite(b) else float("nan")
    return {"d_arrival": mo["arrival"] - mr["arrival"], "r_slip": rel(mo["slip"], mr["slip"]),
            "r_peak": rel(mo["peak"], mr["peak"]), "ref_slip": mr["slip"]}


def _fmt(v, f):
    return ("{:" + f + "}").format(v) if np.isfinite(v) else "nan".rjust(len(("{:" + f + "}").format(0.0)))


# ---------------------------------------------------------------------- main
def parse_src(s, default_fmt):
    if "=" not in s:
        sys.exit(f"source {s!r}: expected LABEL=DIR")
    lab, rest = s.split("=", 1)
    fmt = default_fmt
    if ":" in rest and rest.rsplit(":", 1)[1] in ("scec_si", "scec_mpa", "nitsche", "side_si", "side_mpa"):
        rest, fmt = rest.rsplit(":", 1)
    if not os.path.isdir(rest):
        sys.exit(f"source {lab}: {rest} is not a directory")
    return lab, rest, fmt


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stations", required=True)
    ap.add_argument("--prefix", default=None, help="our file prefix (only used in messages)")
    ap.add_argument("--run", action="append", default=[])
    ap.add_argument("--ref", action="append", default=[])
    ap.add_argument("--gold", action="append", default=[])
    ap.add_argument("--t-end", type=float, default=12.0)
    ap.add_argument("--noise-window", type=float, default=0.1, help="moving-average width (s)")
    ap.add_argument("--noise-after", type=float, default=0.5, help="noise starts this long after arrival (s)")
    ap.add_argument("--arrival-tol", type=float, default=0.1)
    ap.add_argument("--slip-tol", type=float, default=0.02)
    ap.add_argument("--slip-min", type=float, default=0.1)
    ap.add_argument("--peak-tol", type=float, default=0.10)
    ap.add_argument("--peak-min-stations", type=int, default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--per-side", action="store_true",
                    help="two split-node files per station (far, near); see the header")
    ap.add_argument("--side-tol", type=float, default=0.03, help="pass item 6 tolerance")
    ap.add_argument("--ref-prefix", default=None,
                    help="per-side reference and gold files must start with this "
                         "(default: <--prefix>_, e.g. tpv7_drdg3d_...)")
    ap.add_argument("--list-columns", action="store_true")
    a = ap.parse_args()

    stations = read_stations(a.stations)
    if a.per_side:
        srcs = ([("run",) + parse_src(s, "side_si") for s in a.run]
                + [("ref",) + parse_src(s, "side_mpa") for s in a.ref]
                + [("gold",) + parse_src(s, "side_mpa") for s in a.gold])
    else:
        srcs = ([("run",) + parse_src(s, "nitsche") for s in a.run]
                + [("ref",) + parse_src(s, "scec_mpa") for s in a.ref]
                + [("gold",) + parse_src(s, "scec_si") for s in a.gold])
    if not srcs:
        sys.exit("give at least one --run, --ref or --gold")
    ref_prefix = a.ref_prefix if a.ref_prefix is not None else (a.prefix + "_" if a.prefix else None)
    data = {}
    for kind, lab, d, fmt in srcs:
        for st in stations:
            if a.per_side:
                fp = None if kind == "run" else ref_prefix
                ff, fn = find_file(d, st, "farside", fp), find_file(d, st, "nearside", fp)
                if (ff is None) != (fn is None):
                    sys.exit(f"{d}: station {st['name']} has only its {'far' if fn is None else 'near'}-side file")
                if ff is None:
                    continue
                f = ff
                data[(lab, st["name"])] = combine_sides(load_side(ff, fmt), load_side(fn, fmt), ff)
            else:
                f = find_file(d, st)
                if f is None:
                    continue
                data[(lab, st["name"])] = load(f, fmt)
            if a.list_columns:
                dd = data[(lab, st["name"])]
                print(f"{lab:14s} {st['name']:16s} {fmt:8s} {dd['_ncol']} cols "
                      f"t={dd['t'][0]:.4g}..{dd['t'][-1]:.4g}  {os.path.basename(f)}")
    if a.list_columns:
        for kind, lab, d, fmt in srcs:
            n = sum(1 for st in stations if (lab, st["name"]) in data)
            print(f"# {lab}: {n} of {len(stations)} stations found in {d}")
        return
    runs = [lab for kind, lab, d, fmt in srcs if kind == "run"]
    refs = [lab for kind, lab, d, fmt in srcs if kind == "ref"]
    labels = [lab for kind, lab, d, fmt in srcs]
    # One common end time for every source, so a run that stops early is
    # never compared with a reference measured further on.
    t_end = a.t_end
    for lab in labels:
        ends = [data[(lab, s["name"])]["t"][-1] for s in stations if (lab, s["name"]) in data]
        if ends and min(ends) < t_end - 1e-9:
            print(f"# WARNING: {lab} ends at {min(ends):.4f} s < t_end = {t_end} s; "
                  f"every source is measured to t = {min(ends):.4f} s")
            t_end = min(ends)
    M = {}
    for (lab, name), d in data.items():
        M[(lab, name)] = measures(d, t_end, a.noise_window, a.noise_after)

    print(f"# stations {a.stations}; t_end {t_end} s; noise: {a.noise_window} s centred moving average,"
          f" from arrival + {a.noise_after} s")
    print(f"{'station':16s} {'source':14s} {'arrival':>8s} {'peak V':>8s} {'slip':>8s} {'state':>9s}"
          f" {'sn rms':>7s} {'th rms':>7s} {'tv rms':>7s}   misfit vs {refs[0] if refs else '-'}")
    for st in stations:
        n = st["name"]
        for lab in labels:
            if (lab, n) not in M:
                continue
            m = M[(lab, n)]
            line = (f"{n:16s} {lab:14s} {_fmt(m['arrival'], '8.4f')} {_fmt(m['peak'], '8.3f')}"
                    f" {_fmt(m['slip'], '8.4f')} {_fmt(m.get('state', float('nan')), '9.4f')}"
                    f" {_fmt(m.get('noise_sn', float('nan')), '7.3f')} {_fmt(m.get('noise_th', float('nan')), '7.3f')}"
                    f" {_fmt(m.get('noise_tv', float('nan')), '7.3f')}")
            if refs and lab != refs[0] and (refs[0], n) in M:
                f = misfit(m, M[(refs[0], n)])
                line += (f"   dt {f['d_arrival']:+.3f} s  slip {100 * f['r_slip']:+6.2f}%"
                         f"  peak {100 * f['r_peak']:+6.1f}%")
            print(line)
    summary = {"stations": a.stations, "t_end": t_end, "noise_window": a.noise_window,
               "noise_after": a.noise_after,
               "measures": {f"{lab}|{n}": m for (lab, n), m in M.items()}}

    # noise medians
    print("\n# median noise RMS over stations (MPa)")
    for lab in labels:
        vals = {k: [M[(lab, s['name'])].get(k) for s in stations if (lab, s['name']) in M]
                for k in ("noise_sn", "noise_th", "noise_tv")}
        cells = [f"{k[-2:]} {np.nanmedian([v for v in vs if v is not None]):.3f}"
                 if any(v is not None for v in vs) else f"{k[-2:]}   n/a" for k, vs in vals.items()]
        print(f"{lab:14s} " + "  ".join(cells))

    if not refs:
        _write(a.json, summary)
        return
    r0 = refs[0]
    # reference spread
    for r in refs[1:]:
        fs = [misfit(M[(r, s["name"])], M[(r0, s["name"])]) for s in stations
              if (r, s["name"]) in M and (r0, s["name"]) in M]
        if fs:
            print(f"\n# reference spread {r} vs {r0}: max |dt| {np.nanmax([abs(f['d_arrival']) for f in fs]):.3f} s,"
                  f" max |slip| {100 * np.nanmax([abs(f['r_slip']) for f in fs]):.2f}%,"
                  f" peak within {100 * a.peak_tol:.0f}% at "
                  f"{sum(abs(f['r_peak']) <= a.peak_tol for f in fs)} of {len(fs)}")
    # convergence + pass test
    med = {}
    for lab in runs:
        fs = [misfit(M[(lab, s["name"])], M[(r0, s["name"])]) for s in stations
              if (lab, s["name"]) in M and (r0, s["name"]) in M]
        if not fs:
            continue
        med[lab] = {"arrival": float(np.nanmedian([abs(f["d_arrival"]) for f in fs])),
                    "slip": float(np.nanmedian([abs(f["r_slip"]) for f in fs if f["ref_slip"] > a.slip_min])),
                    "peak": float(np.nanmedian([abs(f["r_peak"]) for f in fs]))}
    if med:
        print(f"\n# median |misfit| vs {r0} along the run order (pass item 4)")
        for lab in runs:
            if lab in med:
                print(f"{lab:14s} arrival {med[lab]['arrival']:.4f} s  slip {100 * med[lab]['slip']:.2f}%"
                      f"  peak {100 * med[lab]['peak']:.2f}%")
    summary["median_misfit"] = med
    if runs:
        last = runs[-1]
        fs = {s["name"]: misfit(M[(last, s["name"])], M[(r0, s["name"])]) for s in stations
              if (last, s["name"]) in M and (r0, s["name"]) in M}
        missing = [s["name"] for s in stations if s["name"] not in fs]
        p1 = all(abs(f["d_arrival"]) <= a.arrival_tol for f in fs.values()) and not missing
        p2 = all(abs(f["r_slip"]) <= a.slip_tol for f in fs.values() if f["ref_slip"] > a.slip_min) and not missing
        npk = sum(abs(f["r_peak"]) <= a.peak_tol for f in fs.values())
        p3 = None if a.peak_min_stations is None else npk >= a.peak_min_stations
        order = [lab for lab in runs if lab in med]
        p4 = None
        if len(order) >= 2:
            p4 = all(med[order[i + 1]][k] < med[order[i]][k] for i in range(len(order) - 1)
                     for k in ("arrival", "slip", "peak"))
        print(f"\n# pass test for run {last} against {r0}"
              + (f" (stations missing: {', '.join(missing)})" if missing else ""))
        print(f"  1 arrival within {a.arrival_tol} s at every station: {'PASS' if p1 else 'FAIL'}")
        print(f"  2 slip within {100 * a.slip_tol:.0f}% where ref slip > {a.slip_min} m: {'PASS' if p2 else 'FAIL'}")
        print(f"  3 peak within {100 * a.peak_tol:.0f}% at {npk} of {len(fs)} stations"
              + ("" if p3 is None else f" (need {a.peak_min_stations}): {'PASS' if p3 else 'FAIL'}"))
        print("  4 medians shrink along the run order: "
              + ("n/a (one run)" if p4 is None else ("PASS" if p4 else "FAIL")))
        summary["pass"] = {"run": last, "ref": r0, "1": p1, "2": p2, "3": p3, "3_count": npk, "4": p4,
                           "missing": missing}
        if a.per_side:
            print(f"\n# per-side final h-disp (m) and far/near ratio")
            for st in stations:
                n = st["name"]
                cells = []
                for lab in labels:
                    if (lab, n) in M:
                        m = M[(lab, n)]
                        ratio = m["h_far"] / m["h_near"] if m["h_near"] != 0.0 else float("nan")
                        cells.append(f"{lab} {m['h_far']:+.3f}/{m['h_near']:+.3f} ({_fmt(ratio, '.2f')})")
                print(f"{n:16s} " + "   ".join(cells))
            bad, n6 = [], 0
            for st in stations:
                n = st["name"]
                if (last, n) not in M or (r0, n) not in M:
                    continue
                for side, k in (("far", "h_far"), ("near", "h_near")):
                    ref = M[(r0, n)][k]
                    if abs(ref) <= a.slip_min:
                        continue
                    n6 += 1
                    r = (M[(last, n)][k] - ref) / ref
                    if abs(r) > a.side_tol:
                        bad.append(f"{n} {side} {100 * r:+.1f}%")
            p6 = not bad and not missing and n6 > 0
            print(f"  6 each side's final h-disp within {100 * a.side_tol:.0f}% where the reference side "
                  f"moves > {a.slip_min} m ({n6} sides): {'PASS' if p6 else 'FAIL'}"
                  + (f" ({'; '.join(bad)})" if bad else ""))
            summary["pass"]["6"] = p6
            summary["pass"]["6_bad"] = bad
    _write(a.json, summary)


def _write(path, summary):
    if path:
        with open(path, "w") as f:
            json.dump(summary, f, indent=1, default=float)
        print(f"# wrote {path}")


if __name__ == "__main__":
    main()
