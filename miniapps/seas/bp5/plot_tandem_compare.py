#!/usr/bin/env python3
"""Quick comparison plot of specific Tandem BP5 output files."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from visualize_results import load_bp5_file, plot_station, plot_closeup

import re


def parse_filename(filepath):
    """Extract (polynomial degree, x2, x3) from filename."""
    base = os.path.basename(filepath)
    m = re.search(r'_p(\d+)_x2_(-?\d+)_x3_(-?\d+)', base)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    return None, None, None


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Compare specific Tandem BP5 files")
    parser.add_argument("files", nargs="+", help="BP5 SCEC-format files to compare")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--save", action="store_true", default=True)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    os.makedirs(args.output_dir, exist_ok=True)

    colors = ["#000000", "#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
    linestyles = ["--", "-", "-.", ":"]

    # Group files by station (x2, x3)
    station_files = {}
    for f in args.files:
        p, x2, x3 = parse_filename(f)
        if p is None:
            print(f"WARNING: Cannot parse {f}, skipping")
            continue
        key = (x2, x3)
        if key not in station_files:
            station_files[key] = []
        station_files[key].append((f, p))

    for (x2, x3), file_list in sorted(station_files.items()):
        station_name = f"fltst_strk{x2:+03d}dp{x3:+03d}"
        print(f"\nStation x2={x2}km, x3={x3}km:")

        datasets = []
        for i, (filepath, pdeg) in enumerate(sorted(file_list, key=lambda x: -x[1])):
            data = load_bp5_file(filepath)
            if data is None:
                print(f"  WARNING: No data in {filepath}")
                continue
            label = f"Tandem p{pdeg}"
            npts = len(data["time_s"])
            t_max = data["time_yr"][-1]
            print(f"  {label}: {npts} pts, {t_max:.2f} years")
            datasets.append((label, data, colors[i % len(colors)],
                             linestyles[i % len(linestyles)]))

        if not datasets:
            continue

        # Full time range
        fname = os.path.join(args.output_dir,
                             f"bp5_tandem_x2_{x2}_x3_{x3}.png")
        plot_station(datasets, station_name, x2, x3, save_path=fname)

        # Closeup (first year)
        fname_close = os.path.join(args.output_dir,
                                   f"bp5_tandem_x2_{x2}_x3_{x3}_closeup.png")
        plot_closeup(datasets, station_name, x2, x3, t_max_yr=1.0,
                     save_path=fname_close)

    print(f"\nPlots saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
