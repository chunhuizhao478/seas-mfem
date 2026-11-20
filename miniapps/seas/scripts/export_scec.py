#!/usr/bin/env python3
"""Write a Nitsche run's outputs as SCEC CVWS upload files (on-fault, off-fault, cplot).

Reads one run's results directory (the driver's --output-dir, or a gold's
results/ folder) and writes, into --out, one file per station and one
contour-plot file per fault in the formats of the SCEC benchmark
descriptions (e.g. TPV26/27 v13 Parts 8-10, TPV12/13 v6 Part 7, TPV22/23
v08 Parts 6-8).  The run is not changed; nothing is recomputed.

On-fault stations, from <prefix>_station_<name>.dat
(t slip1 slip2 V1 V2 tau1 tau2 sigma_n state; SI; 1 = dip, 2 = strike;
frame in <prefix>_station_conventions.txt):

    t        = t
    h-slip   = slip2    h-slip-rate = V2    h-shear-stress = tau2 / 1e6
    v-slip   = slip1    v-slip-rate = V1    v-shear-stress = tau1 / 1e6
    n-stress = -sigma_n / 1e6     (SCEC: MPa, positive in extension; TPV101-104:
                                   +sigma_n / 1e6, positive in compression)
    psi / log-theta = state       (rate-and-state problems only)

slip = u(far) - u(near) along shat (strike) and dhat (down-dip), the SCEC
signs: h > 0 right lateral, v > 0 far side down.  These are the columns
scripts/compare_nitsche_fault_stations.py compares with the SCEC reference
files (--signed checks the signs).

Per-side stations (TPV6/TPV7, [output.stations] per_side = true), from
<prefix>_station_<name>_{farside,nearside}.dat (t h_disp h_vel h_stress
v_disp v_vel v_stress n_disp n_vel n_stress; m, m/s, Pa, SCEC signs):
written as is with the stresses in MPa, one file per side.

Off-fault stations, from <prefix>_offfault_<name>.dat (both off-fault
writers; t h-disp h-vel v-disp v-vel n-disp n-vel, m and m/s, SCEC signs):
written as is with the SCEC header.

Times: the driver writes t with 11 significant digits, so equal steps read
back unequal (~1e-6 relative); the server refuses its digital filters on
unequal steps.  When every row lies within 1e-3 dt of t0 + i dt (all
fixed-step runs), t is rewritten as t0 + i dt with the mean dt, and the
driver's extra row at the last step (a shorter interval) is dropped.

Contour plot, from <prefix>_rupture_time_r<rank>.csv (x y z t per owned
fault DOF, t = first |V| > 1 mm/s, -1 = never): the DG copies of a node
are merged (earliest time), every node goes to the nearest fault plane of
the problem's frames (scripts/make_fault_stations.py _FRAMES; otherwise
the plane Y = 0, strike +X, down-dip -Z), and

    j = (x - o) . strike,   k = (x - o) . down-dip     (m; t = 1.0E+09 never)

so k is the distance down-dip (TPV12: along the 60 deg dip, not depth).
Problems with two faults get one file per fault: cplot_main and
cplot_branch (TPV15), cplot_1 and cplot_2 (TPV22).

The time series keep every output row up to --t-end (default: the
description's duration where the table below has it, else all rows);
--every n keeps every n-th row (for the ~10 MB upload limit).
The time column is written with 20.12e, the others with 14.6e, as the
descriptions recommend; the header names the time step (the first
interval) and the largest relative deviation of the other intervals.

Usage (from miniapps/seas):
    python3 -I scripts/export_scec.py --problem tpv205 \\
        --run-dir tpv205/gold/<gold>/results --out <dir> \\
        --element-size "200 m, p = 2" --code-version <commit> --t-end 12

Scope: Chunhui Zhao, 2026-11-19 (SCEC database upload, thread "Inventory
TPV golds for SCEC database").
"""
import argparse
import datetime
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_fault_stations as mfs  # noqa: E402  (fault frames of the complex geometries)

NEVER = 1.0e9
ON_BASE = ["t", "h-slip", "h-slip-rate", "h-shear-stress", "v-slip", "v-slip-rate",
           "v-shear-stress", "n-stress"]
ON_DESC = {"t": "Time (s)", "h-slip": "horizontal slip (m)",
           "h-slip-rate": "horizontal slip rate (m/s)",
           "h-shear-stress": "horizontal shear stress (MPa)", "v-slip": "vertical slip (m)",
           "v-slip-rate": "vertical slip rate (m/s)",
           "v-shear-stress": "vertical shear stress (MPa)", "n-stress": "normal stress (MPa)",
           "psi": "state variable psi (dimensionless)",
           "log-theta": "log10 of the state variable theta (log10 s)"}
SIDE = ["t", "h-disp", "h-vel", "h-stress", "v-disp", "v-vel", "v-stress", "n-disp", "n-vel",
        "n-stress"]
SIDE_DESC = {"t": "Time (s)", "h-disp": "horizontal displacement (m)",
             "h-vel": "horizontal velocity (m/s)", "h-stress": "horizontal shear stress (MPa)",
             "v-disp": "vertical displacement (m)", "v-vel": "vertical velocity (m/s)",
             "v-stress": "vertical shear stress (MPa)", "n-disp": "normal displacement (m)",
             "n-vel": "normal velocity (m/s)", "n-stress": "normal stress (MPa)"}
OFF = ["t", "h-disp", "h-vel", "v-disp", "v-vel", "n-disp", "n-vel"]
OFF_DESC = {"t": "Time (s)", "h-disp": "horizontal displacement (m)",
            "h-vel": "horizontal velocity (m/s)", "v-disp": "vertical displacement (m)",
            "v-vel": "vertical velocity (m/s)", "n-disp": "normal displacement (m)",
            "n-vel": "normal velocity (m/s)"}

# problem -> (SCEC name, on-fault fields, per-side, n-stress sign, duration s)
# From the SCEC descriptions and upload instructions (strike.scec.org/cvws, read
# 2026-11-19): the field list line is checked by the server; n-stress is positive in
# extension ("ext", the default) except TPV101-104, positive in compression ("comp");
# duration = "Run the model for times from 0.0 to ... seconds" (None: not read).
PROBLEMS = {
    "tpv102": ("TPV102", ON_BASE + ["log-theta"], False, "comp", 12.0),
    "tpv104": ("TPV104", ON_BASE + ["psi"], False, "comp", 12.0),
    "tpv205": ("TPV205", ON_BASE[:7], False, "ext", 12.0),
    "tpv6": ("TPV6", SIDE, True, "ext", 12.0),
    "tpv7": ("TPV7", SIDE, True, "ext", 12.0),
    "tpv31": ("TPV31", ON_BASE, False, "ext", 15.0),
    "tpv26": ("TPV26", ON_BASE, False, "ext", 13.0),
    "tpv27": ("TPV27", ON_BASE, False, "ext", 13.0),
    "tpv12": ("TPV12", ON_BASE, False, "ext", 8.0),
    "tpv15": ("TPV15", ON_BASE, False, "ext", 12.0),
    "tpv22": ("TPV22", ON_BASE, False, "ext", 15.0),
    "tpv34": ("TPV34", ON_BASE, False, "ext", 20.0),
}
UPLOAD_LIMIT = 10e6   # "about 10 MB" per file (TPV205 description)


def rows(path):
    out = []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#") or not (s[0].isdigit() or s[0] in "+-."):
                continue
            out.append([float(c) for c in s.split()])
    if not out:
        sys.exit(f"export_scec: {path}: no data rows")
    n = len(out[0])
    if any(len(r) != n for r in out):
        sys.exit(f"export_scec: {path}: rows of different lengths")
    return out


def thin(data, t_end, every):
    """Rows with t <= t_end (+ half an interval), every n-th one, regular times."""
    if t_end is not None:
        dt = data[1][0] - data[0][0] if len(data) > 1 else 0.0
        data = [r for r in data if r[0] <= t_end + 0.5 * dt]
        if data[-1][0] < t_end - 1.01 * dt:   # the last row may fall within one step
            sys.exit(f"export_scec: the run ends at t = {data[-1][0]} s < --t-end {t_end} s")
    if every > 1:
        data = data[::every]
    data, dt, jit = regular(data)
    return data


def regular(data):
    """The driver writes t with 11 significant digits, so equal steps read back unequal
    (relative jitter ~1e-6), and the server then refuses its filters.  When every row
    but the last lies within 1e-3 dt of t0 + i dt (a fixed-step run; the driver's sum
    t += dt drifts ~1e-5 dt over 2e5 steps), t is rewritten as t0 + i dt (dt the mean
    of those intervals) and a shorter last interval (the driver's extra row at the last
    step) is dropped.  -> (rows, dt or None, max jitter)."""
    if len(data) < 3:
        return data, None, 0.0
    t = [r[0] for r in data]
    last_short = abs((t[-1] - t[-2]) - (t[-2] - t[-3])) > 1e-3 * (t[-2] - t[-3])
    body = data[:-1] if last_short else data
    n = len(body)
    dt = (body[-1][0] - body[0][0]) / (n - 1)
    jit = max(abs(body[i][0] - (body[0][0] + i * dt)) for i in range(n)) / dt
    if jit > 1e-3:
        return data, None, jit
    return [[body[0][0] + i * dt] + r[1:] for i, r in enumerate(body)], dt, jit


def step_info(data):
    if len(data) < 2:
        return 0.0, 0.0
    dts = [b[0] - a[0] for a, b in zip(data, data[1:])]
    dt0 = dts[0]
    dev = max(abs(d - dt0) for d in dts[:-1]) / dt0 if len(dts) > 1 else 0.0
    return dt0, dev


def fmt_row(r):
    return f"{r[0]:20.12e}" + "".join(f" {v:14.6e}" for v in r[1:])


def header(a, kind, location, data, fields, desc, extra=()):
    dt0, dev = step_info(data)
    h = [f"# SCEC {kind} file, written by scripts/export_scec.py from {a.run_dir}",
         "# This is the file header:",
         f"# problem={a.scec_problem}", f"# author={a.author}", f"# date={a.date}",
         f"# code={a.code}"]
    if a.code_version:
        h.append(f"# code_version={a.code_version}")
    if a.element_size:
        h.append(f"# element_size={a.element_size}")
    h += [f"# time_step={dt0:.12e}", f"# num_time_steps={len(data)}",
          f"# location={location}"]
    if dev > 1e-9:
        h.append(f"# note: the time steps are not all equal (largest relative deviation "
                 f"{dev:.2e}; the last interval may be shorter)")
    h += list(extra)
    h += [f"# Column #{i + 1} = {desc[k]}" for i, k in enumerate(fields)]
    h += ["# The line below lists the names of the data fields:", " ".join(fields),
          "# Here is the time-series data."]
    return h


def write(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


_NAME = re.compile(r"^(?:([A-Za-z]+[0-9]*)_)?x2_(-?[0-9.]+)_x3_(-?[0-9.]+)$")


def scec_station(name):
    """x2_<s>_x3_<d> -> faultst<s x10>dp<d x10> (the TPV26/34 form); a fault tag stays as
    a prefix (main_, fault2_ ...).  flt_ (TPV102) names: the same.  Else the name."""
    m = _NAME.match(name)
    s = d = tag = None
    if m:
        tag, s, d = m.group(1), float(m.group(2)), float(m.group(3))
    else:
        m = re.match(r"^flt_(n?)([0-9.]+)_([0-9.]+)$", name)
        if m:
            s = -float(m.group(2)) if m.group(1) else float(m.group(2))
            d = float(m.group(3))
    if s is None:
        return name, None

    def f(v, signed=True):
        return ("-" if (signed and v < 0) else "") + f"{round(abs(v) * 10):03d}"
    loc = f"on fault, {s:g} km along strike, {d:g} km down-dip" + (f" ({tag})" if tag else "")
    return (f"{tag}_" if tag else "") + f"faultst{f(s)}dp{f(d, False)}", loc


def on_fault(a, fields, per_side, nsign):
    names = sorted({m.group(1) for f in os.listdir(a.run_dir)
                    for m in [re.match(re.escape(a.prefix) + r"_station_(.+)\.dat$", f)] if m
                    and not m.group(1).endswith("_copies")})
    if per_side:
        names = sorted({n.rsplit("_", 1)[0] for n in names
                        if n.endswith("_farside") or n.endswith("_nearside")})
    else:
        names = [n for n in names if not n.endswith(("_farside", "_nearside"))]
    if not names:
        sys.exit(f"export_scec: no {a.prefix}_station_*.dat in {a.run_dir}")
    written = []
    for n in names:
        sname, loc = scec_station(n)
        loc = loc or f"on fault, station {n}"
        if per_side:
            for side in ("farside", "nearside"):
                data = rows(os.path.join(a.run_dir, f"{a.prefix}_station_{n}_{side}.dat"))
                if len(data[0]) < 10:
                    sys.exit(f"export_scec: {n}_{side}: {len(data[0])} columns, need 10")
                data = [[r[0], r[1], r[2], r[3] / 1e6, r[4], r[5], r[6] / 1e6, r[7], r[8],
                         r[9] / 1e6] for r in data]
                data = thin(data, a.t_end, a.every)
                lines = header(a, "on-fault (split node)", f"{loc}, {side[:-4]} side", data,
                               SIDE, SIDE_DESC, (f"# our station name: {n}",))
                p = os.path.join(a.out, f"{sname}_{side}.dat")
                write(p, lines + [fmt_row(r) for r in data])
                written.append(p)
            continue
        data = rows(os.path.join(a.run_dir, f"{a.prefix}_station_{n}.dat"))
        nc = len(data[0])
        state = [k for k in fields if k in ("psi", "log-theta")]
        if nc < 8 or (state and nc < 9):
            sys.exit(f"export_scec: {n}: {nc} columns, need {9 if state else 8}")
        out = []
        for r in data:
            v = {"t": r[0], "h-slip": r[2], "h-slip-rate": r[4], "h-shear-stress": r[6] / 1e6,
                 "v-slip": r[1], "v-slip-rate": r[3], "v-shear-stress": r[5] / 1e6,
                 "n-stress": nsign * r[7] / 1e6}
            if state:
                v[state[0]] = r[8]
            out.append([v[k] for k in fields])
        out = thin(out, a.t_end, a.every)
        lines = header(a, "on-fault", loc, out, fields, ON_DESC, (f"# our station name: {n}",))
        p = os.path.join(a.out, f"{sname}.dat")
        write(p, lines + [fmt_row(r) for r in out])
        written.append(p)
    return written


def off_fault(a):
    paths = sorted(glob.glob(os.path.join(a.run_dir, f"{a.prefix}_offfault_*.dat")))
    written = []
    for src in paths:
        name = os.path.basename(src)[len(a.prefix) + len("_offfault_"):-4]
        loc = None
        with open(src) as f:
            for ln in f:
                m = re.search(r"at \(([^)]*)\) m|location= *\(([^)]*)\) m", ln)
                if m:
                    loc = "off fault, (" + (m.group(1) or m.group(2)) + ") m, global frame"
                    break
        data = rows(src)
        if len(data[0]) != 7:
            sys.exit(f"export_scec: {src}: {len(data[0])} columns, the off-fault layout has 7")
        data = thin(data, a.t_end, a.every)
        lines = header(a, "off-fault", loc or f"off fault, station {name}", data, OFF, OFF_DESC)
        p = os.path.join(a.out, f"{name}.dat")
        write(p, lines + [fmt_row(r) for r in data])
        written.append(p)
    return written


def frames(problem):
    if problem in mfs._FRAMES:
        return mfs._FRAMES[problem]
    return {"fault": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0))}


def cplot(a, problem):
    paths = sorted(glob.glob(os.path.join(a.run_dir, f"{a.prefix}_rupture_time_r*.csv")))
    if not paths:
        return [], f"no {a.prefix}_rupture_time_r*.csv: no contour-plot file"
    nodes = {}
    ndof = 0
    for p in paths:
        for r in rows(p):
            if len(r) != 4:
                sys.exit(f"export_scec: {p}: expected 'x y z t_rupture'")
            ndof += 1
            key = tuple(round(c / a.tol) for c in r[:3])
            t = r[3] if r[3] >= 0.0 else NEVER
            if key in nodes:
                nodes[key][3] = min(nodes[key][3], t)
            else:
                nodes[key] = list(r[:3]) + [t]
    fr = frames(problem)
    per = {f: [] for f in fr}
    worst = 0.0
    for x in nodes.values():
        best = None
        for f, (o, st, dp) in fr.items():
            d = [x[c] - o[c] for c in range(3)]
            nrm = [st[1] * dp[2] - st[2] * dp[1], st[2] * dp[0] - st[0] * dp[2],
                   st[0] * dp[1] - st[1] * dp[0]]
            dist = abs(sum(d[c] * nrm[c] for c in range(3)))
            j = sum(d[c] * st[c] for c in range(3))
            k = sum(d[c] * dp[c] for c in range(3))
            # a branch (TPV15) starts at the junction: nodes behind it belong to the other fault
            if f == "branch" and j < -a.tol:
                dist = math.inf
            if best is None or dist < best[0]:
                best = (dist, f, j, k + 0.0)
        worst = max(worst, best[0])
        per[best[1]].append((best[2], best[3], x[3]))
    if worst > a.plane_tol:
        sys.exit(f"export_scec: a rupture-time node lies {worst:.3g} m off every fault plane "
                 f"of {problem} (tolerance {a.plane_tol} m)")
    written = []
    for f, pts in per.items():
        if not pts:
            continue
        pts.sort(key=lambda r: (r[1], r[0]))
        ruptured = sum(1 for r in pts if r[2] < NEVER)
        h = [f"# SCEC contour-plot file, written by scripts/export_scec.py from {a.run_dir}",
             "# This is the file header:", f"# problem={a.scec_problem}",
             f"# author={a.author}", f"# date={a.date}", f"# code={a.code}"]
        if a.code_version:
            h.append(f"# code_version={a.code_version}")
        if a.element_size:
            h.append(f"# element_size={a.element_size}")
        h += [f"# fault={f}: {len(pts)} nodes ({ruptured} ruptured); a node's time is the "
              f"earliest of its DG copies; {ndof} DOFs in {len(paths)} rank files",
              "# Column #1 = Distance along strike (m)",
              "# Column #2 = Distance down-dip (m)",
              "# Column #3 = Rupture time (s)", "j k t"]
        # file names of the descriptions: cplot, cplot_main / cplot_branch (TPV15),
        # cplot_1 / cplot_2 (TPV22)
        name = "cplot" if len(fr) == 1 else "cplot_" + {"fault1": "1", "fault2": "2"}.get(f, f)
        p = os.path.join(a.out, name + ".dat")
        write(p, h + [f"{j:14.6e} {k:14.6e} {t:14.6e}" for j, k, t in pts])
        written.append(p)
    return written, ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--problem", required=True, choices=sorted(PROBLEMS))
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefix", default=None, help="file prefix (default: the problem)")
    ap.add_argument("--author", default="Chunhui Zhao")
    ap.add_argument("--date", default=datetime.date.today().strftime("%Y/%m/%d"))
    ap.add_argument("--code", default="seas-mfem (seas_nitsche_dyn_driver)")
    ap.add_argument("--code-version", default="", help="the commit of the run")
    ap.add_argument("--element-size", default="", help="e.g. '200 m, p = 2'")
    ap.add_argument("--t-end", type=float, default=None,
                    help="last time (s); default: the description's duration")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--fields", default=None, help="override the on-fault field list")
    ap.add_argument("--tol", type=float, default=1e-6, help="DG copy merge tolerance (m)")
    ap.add_argument("--plane-tol", type=float, default=1.0,
                    help="largest distance of a rupture-time node from its fault plane (m)")
    a = ap.parse_args()
    scec, fields, per_side, nconv, duration = PROBLEMS[a.problem]
    a.scec_problem = scec
    if a.t_end is None:
        a.t_end = duration
    a.prefix = a.prefix or a.problem
    if a.fields:
        fields = a.fields.split()
        bad = [k for k in fields if k not in (SIDE_DESC if per_side else ON_DESC)]
        if bad or fields[0] != "t":
            sys.exit(f"export_scec: --fields {a.fields!r}: unknown {bad} (or t not first)")
    if a.every < 1:
        sys.exit("export_scec: --every must be >= 1")
    if not os.path.isdir(a.run_dir):
        sys.exit(f"export_scec: no directory {a.run_dir}")
    os.makedirs(a.out, exist_ok=True)
    on = on_fault(a, fields, per_side, 1.0 if nconv == "comp" else -1.0)
    off = off_fault(a)
    cp, why = cplot(a, a.problem)
    big = [p for p in on + off + cp if os.path.getsize(p) > UPLOAD_LIMIT]
    print(f"{scec}: {len(on)} on-fault, {len(off)} off-fault, {len(cp)} contour-plot files "
          f"in {a.out}" + (f" ({why})" if why else ""))
    for p in big:
        print(f"  WARNING {p}: {os.path.getsize(p) / 1e6:.1f} MB > the ~10 MB upload limit "
              "(use --every)")


if __name__ == "__main__":
    main()
