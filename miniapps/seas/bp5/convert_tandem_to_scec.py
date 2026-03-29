#!/usr/bin/env python3
"""
Convert Tandem raw fault probe output to SCEC BP5 benchmark format.

Tandem's fault_probe_output writes raw quantities:
  Time | state(psi) | slip0 | slip1 | traction0 | traction1 | slip-rate0 | slip-rate1 | normal-stress

SCEC BP5 benchmark format (8 columns):
  Time(s) | Slip_2(m) | Slip_3(m) | Slip_rate_2(log10 m/s) | Slip_rate_3(log10 m/s) |
  Shear_stress_2(MPa) | Shear_stress_3(MPa) | State(log10 s)

Component mapping (BP5 fault at Y=0 in Tandem coordinates):
  Tandem component 0 = dip (Z direction)     -> SCEC component 3
  Tandem component 1 = along-strike (X dir)  -> SCEC component 2

Usage:
    python convert_tandem_to_scec.py /path/to/tandem/output/fltst_ -o benchmark_data/

    # Specify polynomial degree for filename
    python convert_tandem_to_scec.py /path/to/fltst_ -o benchmark_data/ -p 1

    # Negate stress if sign convention is flipped
    python convert_tandem_to_scec.py /path/to/fltst_ -o benchmark_data/ --negate-stress
"""

import argparse
import glob
import math
import os
import re
import sys
from datetime import datetime

import numpy as np

# ---- BP5 physical parameters (from bp5.lua) ----
V0 = 1.0e-6       # reference slip rate (m/s)
F0 = 0.6          # reference friction coefficient
B = 0.03           # rate-and-state parameter b

# Nucleation zone geometry (Tandem coordinates: x=along-strike, z=depth negative)
H_S = 2.0          # depth of top of transition zone
H_T = 2.0          # transition zone width
H = 12.0           # seismogenic depth extent
W = 12.0           # nucleation zone along-strike width
L_HALF = 30.0      # half-length of VW zone (l/2 = 60/2)

L_NUC = 0.13       # critical slip distance in nucleation zone (m)
L_DEFAULT = 0.14   # critical slip distance elsewhere (m)

# Minimum slip rate for log10 (avoid -inf)
V_FLOOR = 1.0e-300

# 10 SCEC BP5 standard on-fault stations
# (name, x2_km along-strike, x3_km depth-positive-down)
# In Tandem coords: x = x2 (along-strike), z = -x3 (depth is negative)
STATIONS = [
    ("strk-36dp+00", -36,  0),
    ("strk-16dp+00", -16,  0),
    ("strk+00dp+00",   0,  0),
    ("strk+16dp+00",  16,  0),
    ("strk+36dp+00",  36,  0),
    ("strk-24dp+10", -24, 10),
    ("strk-16dp+10", -16, 10),
    ("strk+00dp+10",   0, 10),
    ("strk+16dp+10",  16, 10),
    ("strk+00dp+22",   0, 22),
]


def in_nucleation(x_strike, x3_depth, eps=1e-3):
    """Check if point is in BP5 nucleation zone.

    Args:
        x_strike: along-strike coordinate (km, Tandem X)
        x3_depth: depth positive down (km, SCEC x3)
    """
    d = x3_depth  # depth positive down
    s = x_strike  # along-strike
    return (H_S + H_T <= d + eps and d - eps <= H_S + H_T + H and
            -L_HALF <= s + eps and s - eps <= -L_HALF + W)


def get_L(x_strike, x3_depth):
    """Get critical slip distance L (Dc) at a point."""
    if in_nucleation(x_strike, x3_depth):
        return L_NUC
    return L_DEFAULT


def psi_to_log10_theta(psi, L):
    """Convert Tandem state variable psi to SCEC log10(theta).

    theta = (L / V0) * exp((psi - f0) / b)
    log10(theta) = log10(L / V0) + (psi - f0) / (b * ln(10))
    """
    return math.log10(L / V0) + (psi - F0) / (B * math.log(10.0))


def safe_log10(v):
    """Compute log10 of slip rate magnitude, with floor."""
    v_abs = abs(v)
    if v_abs < V_FLOOR:
        return math.log10(V_FLOOR)
    return math.log10(v_abs)


def parse_station_name(filename):
    """Extract station name from Tandem output filename.

    E.g., 'fltst_strk+00dp+10.csv' -> 'strk+00dp+10'
    """
    base = os.path.basename(filename)
    # Remove prefix and extension
    m = re.search(r'(strk[+-]\d+dp[+-]\d+)', base)
    if m:
        return m.group(1)
    return None


def station_to_coords(station_name):
    """Convert station name to (x2_km, x3_km).

    E.g., 'strk-36dp+00' -> (-36, 0)
          'strk+00dp+10' -> (0, 10)
    """
    m = re.match(r'strk([+-]\d+)dp([+-]\d+)', station_name)
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def read_tandem_probe(filepath):
    """Read Tandem fault probe output file.

    Returns (header_lines, column_names, data_array).
    """
    header_lines = []
    column_names = None
    data_rows = []

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                header_lines.append(stripped)
                continue

            # Check if CSV (comma-separated) or space-separated
            if ',' in stripped:
                parts = [p.strip().strip('"') for p in stripped.split(',')]
            else:
                parts = stripped.split()

            # Try parsing as numeric data
            try:
                row = [float(x) for x in parts]
                data_rows.append(row)
            except ValueError:
                # Non-numeric line = column header
                if column_names is None:
                    column_names = parts
                continue

    data = np.array(data_rows) if data_rows else None
    return header_lines, column_names, data


def convert_file(input_path, output_path, poly_degree=1, negate_stress=False,
                 element_size_m=1000):
    """Convert a single Tandem probe file to SCEC format."""

    station_name = parse_station_name(input_path)
    if station_name is None:
        print(f"  WARNING: Cannot parse station from {input_path}, skipping")
        return False

    x2, x3 = station_to_coords(station_name)
    if x2 is None:
        print(f"  WARNING: Cannot parse coordinates from {station_name}, skipping")
        return False

    L = get_L(x2, x3)

    header_lines, col_names, data = read_tandem_probe(input_path)
    if data is None or len(data) == 0:
        print(f"  WARNING: No data in {input_path}, skipping")
        return False

    ncols = data.shape[1]
    print(f"  Station {station_name}: x2={x2}km, x3={x3}km, L={L}m, "
          f"{len(data)} timesteps, {ncols} columns")

    # Identify columns
    # Expected: Time | psi | slip0 | slip1 | traction0 | traction1 |
    #           slip-rate0 | slip-rate1 | normal-stress
    if ncols == 9:
        i_time = 0
        i_psi = 1
        i_slip0, i_slip1 = 2, 3        # dip, strike
        i_tau0, i_tau1 = 4, 5          # dip, strike
        i_vr0, i_vr1 = 6, 7            # dip, strike
        # i_sn = 8                      # normal stress (not in SCEC output)
    elif ncols == 8:
        # No normal-stress column
        i_time = 0
        i_psi = 1
        i_slip0, i_slip1 = 2, 3
        i_tau0, i_tau1 = 4, 5
        i_vr0, i_vr1 = 6, 7
    else:
        print(f"  WARNING: Unexpected column count {ncols} in {input_path}")
        print(f"  Column names: {col_names}")
        return False

    # Build SCEC output array (8 columns)
    nrows = len(data)
    out = np.zeros((nrows, 8))

    sign = -1.0 if negate_stress else 1.0

    for i in range(nrows):
        t = data[i, i_time]
        psi = data[i, i_psi]
        slip_dip = data[i, i_slip0]     # component 0 = dip
        slip_strike = data[i, i_slip1]  # component 1 = strike
        tau_dip = data[i, i_tau0]       # component 0 = dip
        tau_strike = data[i, i_tau1]    # component 1 = strike
        vr_dip = data[i, i_vr0]        # component 0 = dip
        vr_strike = data[i, i_vr1]     # component 1 = strike

        out[i, 0] = t                              # Time (s)
        out[i, 1] = slip_strike                     # Slip_2 (m)
        out[i, 2] = slip_dip                        # Slip_3 (m)
        out[i, 3] = safe_log10(vr_strike)           # log10(V_2) (m/s)
        out[i, 4] = safe_log10(vr_dip)              # log10(V_3) (m/s)
        out[i, 5] = sign * tau_strike               # Shear_stress_2 (MPa)
        out[i, 6] = sign * tau_dip                  # Shear_stress_3 (MPa)
        out[i, 7] = psi_to_log10_theta(psi, L)     # log10(theta) (s)

    # Write SCEC format
    with open(output_path, 'w') as f:
        f.write("# This is the file header:\n")
        f.write("# problem=SEAS Benchmark BP5-QD\n")
        f.write("# code=tandem\n")
        f.write(f"# modeler=converted from raw Tandem output\n")
        f.write(f"# date={datetime.now().strftime('%c')}\n")
        f.write(f"# element_size={element_size_m} m\n")
        f.write(f"# polynomial_degree={poly_degree}\n")
        f.write(f"# location=({x2}, 0, {-x3})\n")
        f.write(f"# L={L}\n")
        f.write(f"# num_time_steps={nrows}\n")
        f.write("# Column #1 = Time (s)\n")
        f.write("# Column #2 = Slip_2 (m)\n")
        f.write("# Column #3 = Slip_3 (m)\n")
        f.write("# Column #4 = Slip_rate_2 (log10 m/s)\n")
        f.write("# Column #5 = Slip_rate_3 (log10 m/s)\n")
        f.write("# Column #6 = Shear_stress_2 (MPa)\n")
        f.write("# Column #7 = Shear_stress_3 (MPa)\n")
        f.write("# Column #8 = State (log10 s)\n")
        f.write("# The line below lists the names of the data fields\n")
        f.write("t slip_2 slip_3 slip_rate_2 slip_rate_3 "
                "shear_stress_2 shear_stress_3 state\n")
        f.write("# Here is the time-series data.\n")

        for i in range(nrows):
            f.write(f" {out[i,0]:.16E}  "
                    f" {out[i,1]:.6E}  "
                    f" {out[i,2]:.6E}  "
                    f" {out[i,3]:.6E}  "
                    f" {out[i,4]:.6E}  "
                    f" {out[i,5]:.6E}  "
                    f" {out[i,6]:.6E}  "
                    f" {out[i,7]:.6E}\n")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Convert Tandem fault probe output to SCEC BP5 benchmark format")
    parser.add_argument("prefix",
                        help="Tandem output prefix (e.g., path/to/fltst_)")
    parser.add_argument("-o", "--output-dir", default=".",
                        help="Output directory for SCEC-format files")
    parser.add_argument("-p", "--poly-degree", type=int, default=1,
                        help="Polynomial degree (for output filename)")
    parser.add_argument("--element-size", type=int, default=1000,
                        help="Element size in meters (for header)")
    parser.add_argument("--negate-stress", action="store_true",
                        help="Negate stress sign (if Tandem convention differs)")
    parser.add_argument("--ext", default=".csv",
                        help="Tandem output file extension (default: .csv)")
    args = parser.parse_args()

    # Find all probe files matching the prefix
    pattern = f"{args.prefix}*{args.ext}"
    files = sorted(glob.glob(pattern))

    if not files:
        # Try without extension
        pattern = f"{args.prefix}*"
        files = sorted(glob.glob(pattern))
        # Filter out directories
        files = [f for f in files if os.path.isfile(f)]

    if not files:
        print(f"ERROR: No files found matching {pattern}")
        sys.exit(1)

    print(f"Found {len(files)} probe files")
    os.makedirs(args.output_dir, exist_ok=True)

    converted = 0
    for filepath in files:
        station_name = parse_station_name(filepath)
        if station_name is None:
            continue

        x2, x3 = station_to_coords(station_name)
        if x2 is None:
            continue

        # Output filename: bp5qd_tandem_p{degree}_x2_{x2}_x3_{x3}.txt
        out_name = f"bp5qd_tandem_p{args.poly_degree}_x2_{x2}_x3_{x3}.txt"
        out_path = os.path.join(args.output_dir, out_name)

        success = convert_file(filepath, out_path,
                               poly_degree=args.poly_degree,
                               negate_stress=args.negate_stress,
                               element_size_m=args.element_size)
        if success:
            converted += 1
            print(f"  -> {out_path}")

    print(f"\nConverted {converted}/{len(files)} files to SCEC format in {args.output_dir}/")


if __name__ == "__main__":
    main()
