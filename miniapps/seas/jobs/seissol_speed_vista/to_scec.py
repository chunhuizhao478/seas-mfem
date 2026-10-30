#!/usr/bin/env python3
"""SeisSol receiver output -> the station files our comparisons read.

Usage:
  python3 to_scec.py --problem tpv205|tpv104|tpv102 --repo <repo root> \
      (--seissol-dir <run dir>/output | --seissol-prefix <run dir>/output/<prob>) \
      --out <dir> [--commit <sha>]
  --seissol-dir D means prefix D/<prob> (the decks' OutputFile = 'output/<prob>',
  the folder run_pair_vista.sh creates and passes).

Writes, in --out, one file per MFEM station, named and laid out as the MFEM
Nitsche driver writes them (dynamic/fault_station_file.hpp:310-345,
drivers/nitsche_dyn_driver.cpp:2005-2017 and :2055-2070), which is what
tpv102/compare_nitsche_scec.py (load_ours), scripts/compare_nitsche_fault_stations.py
and jobs/vista_profiling/compare_stations.py read:

  <prob>_station_<name>.dat     '#' header lines, then rows
      time slip1 slip2 V1 V2 tau1 tau2 sigma_n <state>        (SI units)
  <prob>_station_conventions.txt the column map below
  surface/<prob>_station_<name>.dat   off-fault receivers, rows
      time vx vy vz   (the native MFEM surface writer layout,
      dynamic/tpv104_setup.hpp:255-265; kept in surface/ so the fault-
      station globs *_station_*.dat of the compare scripts do not pick
      them up)

Station names and order come from the MFEM station CSV
(<prob>/stations/<prob>_fault_stations.csv); SeisSol numbers its on-fault
receivers by line of the pickpoint file, which check_parity.py verifies is
the CSV order (faultreceiver-<k> = CSV row k, 1-based,
DynamicRupture/Output/OutputManager.cpp:331-337).  Off-fault receiver k =
line k of <prob>_receivers.dat = the MFEM native surface list
(ResultWriter/ReceiverWriter.cpp:93-98).

SeisSol files at commit a8e7f04f
--------------------------------
On-fault  <prefix>-faultreceiver-<5 digits>.dat (OutputManager.cpp:331-337):
  TITLE = "..."                         (:393-402)
  VARIABLES = "Time" ,"SRs" ,"SRd" ,... (:403-407; one name per active
                                         component, Labels in
                                         ResultWriter/FaultWriterExecutor.cpp:73-75)
  # Receiver number k / # x1 / # x2 / # x3   (:416-421, the point after
                                         projection onto the fault face)
  # P_01 / # T_s1 / # T_d1               (:447-452, initial stress)
  rows: time then the active variables, tab separated (:556-571)
  Variable order (DataTypes.h:105-127) with the deck's OutputMask
  1 0 0 1 1 1 0 1 0 0 0 0:
    SRs SRd  slip rate along strike / dip      (ReceiverBasedOutput.cpp:271-275,
             TractionsAndFailure: :438-447)
    Mud StV  friction coefficient, state       (:285-289; StV = psi for FL=103,
             theta [s] for FL=3, 0 for LSW: Output/RateAndState.h:25-27,
             ReceiverBasedOutput.h:129)
    Ts0 Td0 Pn0  TOTAL tractions (initial + nucleation + change) on the
             strike / dip / normal axes, tension positive (:308-325)
    Sls Sld  slip along strike / dip           (:346-367)
    ASl      accumulated slip magnitude, the path int |V| dt (:302-306)
  Axes (Misc.cpp:19-34, Geometry/MeshReader.cpp:139-163): with the deck's
  reference point (0, -1e5, 0) the plus side is y < 0, n = +y,
  strike = (1, 0, 0), dip = (0, 0, -1).
Off-fault <prefix>-receiver-<5 digits>.dat (ReceiverWriter.cpp:93-149):
  TITLE, VARIABLES = "Time","s_xx","s_yy","s_zz","s_xy","s_yz","s_xz","v1","v2","v3"
  (Equations/elastic/Model/Datastructures.h:36-37), # x1/# x2/# x3, rows.

Column map (MFEM frame: nhat = (0,-1,0), shat = (1,0,0), dhat = (0,0,-1),
nitsche_dyn_driver.cpp:131-137 and the run's _station_conventions.txt):
  slip1 = Sld     slip2 = Sls     (SCEC v-slip, h-slip)
  V1    = SRd     V2    = SRs     (SCEC v-slip-rate, h-slip-rate)
  tau1  = Td0     tau2  = Ts0     (SCEC v-/h-shear stress, Pa)
  sigma_n = -Pn0                  (compression positive, Pa)
  9th:  tpv205 slip_path = ASl;  tpv104 psi = StV;
        tpv102 log10_theta = log10(StV)  (StV is theta for FL=3; MFEM writes
        log10((Dc/V0) exp((psi - f0)/b)), the same theta)
  Signs: SeisSol's strike and dip axes are MFEM's shat and dhat; the shear
  tractions are n.S.strike / n.S.dip with n = -nhat, equal to -T.shat /
  -T.dhat (check_parity.py verifies this numerically at t = 0, 0.5, 2 s).
  SeisSol's slip rate is (unlimited - updated traction)/eta (same sign as the
  shear traction it relaxes), as MFEM's V is anti-parallel to X = -tau, so
  no sign flip; the script warns if the final h-slip has the opposite sign
  to the initial h-shear stress.
"""
import argparse, glob, math, os, re, sys

SURFACE_NAMES = {
    "tpv205": [],
    "tpv104": ["surf_0_9", "surf_0_n9", "surf_12_6", "surf_12_n6", "surf_n12_6", "surf_n12_n6"],
    "tpv102": ["surf_0_9", "surf_0_n9", "surf_12_6", "surf_12_n6", "surf_n12_6", "surf_n12_n6"],
}
STATE = {"tpv205": ("slip_path", "ASl", None),
         "tpv104": ("psi", "StV", None),
         "tpv102": ("log10_theta", "StV", "log10")}


def read_seissol(path):
    """-> (names, rows, header dict)."""
    names, rows, hdr = None, [], {}
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("TITLE"):
            continue
        if s.startswith("VARIABLES"):
            names = [n.strip().strip('"') for n in re.split(r"\s*,\s*", s.split("=", 1)[1].strip())]
            continue
        if s.startswith("#"):
            w = s[1:].split()
            if len(w) >= 2:
                try:
                    hdr[w[0]] = float(w[-1])
                except ValueError:
                    pass
            continue
        rows.append([float(v) for v in s.split()])
    if names is None:
        raise ValueError(f"{path}: no VARIABLES line")
    for r in rows:
        if len(r) != len(names):
            raise ValueError(f"{path}: row with {len(r)} values, header has {len(names)}")
    return names, rows, hdr


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--problem", required=True, choices=sorted(STATE))
    ap.add_argument("--repo", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--seissol-prefix", help="the deck's OutputFile, e.g. run/output/tpv104")
    g.add_argument("--seissol-dir", help="folder of OutputFile; prefix = DIR/<problem>")
    ap.add_argument("--out", required=True)
    ap.add_argument("--commit", default="", help="SeisSol commit / binary id to record in the headers")
    a = ap.parse_args()
    p = a.problem
    if a.seissol_dir:
        a.seissol_prefix = os.path.join(a.seissol_dir, p)
    csv = os.path.join(a.repo, "miniapps/seas", p, "stations", f"{p}_fault_stations.csv")
    st = []
    for line in open(csv):
        if line.startswith("#") or line.startswith("name"):
            continue
        n, x, y, z = line.strip().split(",")
        st.append((n, float(x), float(y), float(z)))
    os.makedirs(a.out, exist_ok=True)
    label, scol, fn = STATE[p]
    need = ["Time", "SRs", "SRd", "Ts0", "Td0", "Pn0", "Sls", "Sld", scol]
    nwarn = 0
    for k, (name, x, y, z) in enumerate(st, 1):
        src = f"{a.seissol_prefix}-faultreceiver-{k:05d}.dat"
        if not os.path.exists(src):
            sys.exit(f"ERROR: missing {src} (station {name})")
        names, rows, hdr = read_seissol(src)
        miss = [c for c in need if c not in names]
        if miss:
            sys.exit(f"ERROR: {src} lacks columns {miss}; the deck's &Pickpoint OutputMask must "
                     "enable SlipRate, FrictionAndState, TotalTractions, Slip and AccumulatedSlip")
        ix = {c: names.index(c) for c in need}
        px, py, pz = hdr.get("x1", math.nan), hdr.get("x2", math.nan), hdr.get("x3", math.nan)
        dist = math.dist((x, y, z), (px, py, pz))
        dst = os.path.join(a.out, f"{p}_station_{name}.dat")
        nonmono = sum(1 for i in range(1, len(rows)) if rows[i][0] <= rows[i - 1][0])
        with open(dst, "w") as f:
            f.write(f"# fault station {name}: requested ({x:.12g}, {y:.12g}, {z:.12g}) m, SeisSol "
                    f"receiver point ({px:.12g}, {py:.12g}, {pz:.12g}) m (projected onto the fault "
                    f"face; values at its nearest Gauss point), distance {dist:.12g} m\n")
            f.write(f"# Columns: time slip1 slip2 V1 V2 tau1 tau2 sigma_n {label}\n")
            f.write(f"# 1 = dip, 2 = strike; frame and signs in {p}_station_conventions.txt; "
                    f"from {os.path.basename(src)}{' (' + a.commit + ')' if a.commit else ''}\n")
            for r in rows:
                sv = r[ix[scol]]
                if fn == "log10":
                    sv = math.log10(max(sv, 1e-300))
                vals = (r[ix["Time"]], r[ix["Sld"]], r[ix["Sls"]], r[ix["SRd"]], r[ix["SRs"]],
                        r[ix["Td0"]], r[ix["Ts0"]], -r[ix["Pn0"]], sv)
                f.write(" ".join(f"{v:.10e}" for v in vals) + "\n")
        if rows and rows[-1][ix["Sls"]] * rows[0][ix["Ts0"]] < 0:
            print(f"WARNING {name}: final h-slip {rows[-1][ix['Sls']]:.4g} m has the opposite sign "
                  f"to the initial h-shear {rows[0][ix['Ts0']]:.4g} Pa: check the sign map")
            nwarn += 1
        print(f"{name:16s} <- {os.path.basename(src)}: {len(rows)} rows, point {dist:.3g} m from the "
              f"requested one{', %d non-increasing times' % nonmono if nonmono else ''}")
    with open(os.path.join(a.out, f"{p}_station_conventions.txt"), "w") as h:
        h.write("# SeisSol pickpoint output converted by jobs/seissol_speed_vista/to_scec.py\n"
                "# frame: nhat = (0, -1, 0) (MFEM pin normal), shat = (1, 0, 0), dhat = (0, 0, -1);\n"
                "#        SeisSol n = (0, 1, 0), strike = shat, dip = dhat (reference point (0, -1e5, 0))\n"
                f"# columns: t slip1 slip2 V1 V2 tau1 tau2 sigma_n {label}\n"
                "#   slip1 = Sld  slip2 = Sls  V1 = SRd  V2 = SRs  tau1 = Td0  tau2 = Ts0  sigma_n = -Pn0\n"
                f"#   {label} = " + {"slip_path": "ASl (accumulated slip path, m)",
                                     "psi": "StV (state variable psi)",
                                     "log10_theta": "log10(StV), StV = theta (s), aging law"}[label] + "\n")
    surf = SURFACE_NAMES[p]
    rec = sorted(glob.glob(f"{a.seissol_prefix}-receiver-*.dat"))
    if len(rec) != len(surf):
        print(f"WARNING: {len(rec)} off-fault receiver files, {len(surf)} MFEM surface stations")
        nwarn += 1
    if surf:
        os.makedirs(os.path.join(a.out, "surface"), exist_ok=True)
    for k, name in enumerate(surf, 1):
        src = f"{a.seissol_prefix}-receiver-{k:05d}.dat"
        if not os.path.exists(src):
            print(f"WARNING: missing {src} ({name})")
            nwarn += 1
            continue
        names, rows, hdr = read_seissol(src)
        iv = [names.index(c) for c in ("Time", "v1", "v2", "v3")]
        dst = os.path.join(a.out, "surface", f"{p}_station_{name}.dat")
        with open(dst, "w") as f:
            f.write(f"# {p.upper()} surface station: {name} (x={hdr.get('x1')}, "
                    f"fault_perp={hdr.get('x2')}, z={hdr.get('x3')}) from {os.path.basename(src)}"
                    f"{' (' + a.commit + ')' if a.commit else ''}\n")
            f.write("# Columns: time vx vy vz\n")
            for r in rows:
                f.write(" ".join(f"{r[i]:.10e}" for i in iv) + "\n")
        print(f"{name:16s} <- {os.path.basename(src)}: {len(rows)} rows")
    return 1 if nwarn else 0


if __name__ == "__main__":
    sys.exit(main())
