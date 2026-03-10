#!/usr/bin/env python3
"""
Scaling analysis postprocessing script for SEAS BP1.

Parses SCALING_DATA lines from SLURM .out files and generates
strong/weak scaling plots and CSV summaries.
"""

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCALING_RE = re.compile(
    r"SCALING_DATA:\s+"
    r"np=(\d+)\s+"
    r"elements=(\d+)\s+"
    r"steps=(\d+)\s+"
    r"elapsed=([\d.eE+\-]+)\s+"
    r"mesh=(\S+)\s+"
    r"solver=(\S+)\s+"
    r"dg=(\S+)"
)


def parse_scaling_files(directory):
    """Recursively find *.out files and extract SCALING_DATA lines."""
    records = []
    root = Path(directory)
    for out_file in sorted(root.rglob("*.out")):
        with open(out_file, "r") as f:
            for line in f:
                m = SCALING_RE.search(line)
                if m:
                    records.append(
                        {
                            "np": int(m.group(1)),
                            "elements": int(m.group(2)),
                            "steps": int(m.group(3)),
                            "elapsed": float(m.group(4)),
                            "mesh": m.group(5),
                            "solver": m.group(6),
                            "dg": m.group(7),
                            "source_file": str(out_file),
                        }
                    )
    return records


def print_table(records, title):
    """Print parsed records as a formatted table."""
    if not records:
        print(f"\n{title}: No data found.")
        return
    df = pd.DataFrame(records)
    cols = ["np", "elements", "steps", "elapsed", "mesh", "solver", "dg"]
    print(f"\n{'=' * 80}")
    print(f"  {title}")
    print(f"{'=' * 80}")
    print(df[cols].to_string(index=False))
    print()


def strong_scaling_analysis(records, output_dir):
    """Run strong scaling analysis: wall time, speedup, efficiency."""
    df = pd.DataFrame(records).sort_values("np").reset_index(drop=True)
    print_table(records, "Strong Scaling Data")

    # Base case: smallest np
    base = df.iloc[0]
    n_base = base["np"]
    t_base = base["elapsed"]

    df["speedup"] = t_base / df["elapsed"]
    df["efficiency"] = df["speedup"] / (df["np"] / n_base)

    # Save CSV
    csv_path = os.path.join(output_dir, "strong_scaling.csv")
    df[["np", "elements", "steps", "elapsed", "speedup", "efficiency", "mesh", "solver", "dg"]].to_csv(
        csv_path, index=False
    )
    print(f"CSV saved to {csv_path}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    cores = df["np"].values
    elapsed = df["elapsed"].values
    speedup = df["speedup"].values
    efficiency = df["efficiency"].values

    # Panel 1: Wall time vs cores (log-log)
    ax = axes[0]
    ax.loglog(cores, elapsed, "o-", markersize=8)
    ax.set_xlabel("Number of cores")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Wall Time vs Cores")
    ax.grid(True, which="both", ls="--", alpha=0.5)

    # Panel 2: Speedup vs cores (log-log) with ideal line
    ax = axes[1]
    ax.loglog(cores, speedup, "o-", markersize=8, label="Measured")
    ideal_speedup = cores / n_base
    ax.loglog(cores, ideal_speedup, "k--", alpha=0.6, label="Ideal")
    ax.set_xlabel("Number of cores")
    ax.set_ylabel("Speedup")
    ax.set_title("Speedup vs Cores")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)

    # Panel 3: Parallel efficiency vs cores (semilog-x)
    ax = axes[2]
    ax.semilogx(cores, efficiency * 100, "o-", markersize=8)
    ax.axhline(100, color="k", ls="--", alpha=0.6, label="100%")
    ax.set_xlabel("Number of cores")
    ax.set_ylabel("Parallel efficiency (%)")
    ax.set_title("Parallel Efficiency vs Cores")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)

    fig.tight_layout()
    fig_path = os.path.join(output_dir, "strong_scaling.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"Plot saved to {fig_path}")


def weak_scaling_analysis(records, output_dir):
    """Run weak scaling analysis: wall time and efficiency."""
    df = pd.DataFrame(records).sort_values("np").reset_index(drop=True)
    print_table(records, "Weak Scaling Data")

    # Base case: smallest np
    base = df.iloc[0]
    t_base = base["elapsed"]

    df["efficiency"] = t_base / df["elapsed"]

    # Save CSV
    csv_path = os.path.join(output_dir, "weak_scaling.csv")
    df[["np", "elements", "steps", "elapsed", "efficiency", "mesh", "solver", "dg"]].to_csv(
        csv_path, index=False
    )
    print(f"CSV saved to {csv_path}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    cores = df["np"].values
    elapsed = df["elapsed"].values
    efficiency = df["efficiency"].values

    # Panel 1: Wall time vs cores (semilog-x) with ideal flat line
    ax = axes[0]
    ax.semilogx(cores, elapsed, "o-", markersize=8, label="Measured")
    ax.axhline(t_base, color="k", ls="--", alpha=0.6, label="Ideal")
    ax.set_xlabel("Number of cores")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Wall Time vs Cores")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)

    # Panel 2: Weak efficiency vs cores (semilog-x)
    ax = axes[1]
    ax.semilogx(cores, efficiency * 100, "o-", markersize=8)
    ax.axhline(100, color="k", ls="--", alpha=0.6, label="100%")
    ax.set_xlabel("Number of cores")
    ax.set_ylabel("Weak scaling efficiency (%)")
    ax.set_title("Weak Scaling Efficiency vs Cores")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)

    fig.tight_layout()
    fig_path = os.path.join(output_dir, "weak_scaling.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"Plot saved to {fig_path}")


def validate(strong_dir, weak_dir):
    """Run validation checks on parsed data. Return True if all pass."""
    all_pass = True

    dirs = {}
    if strong_dir:
        dirs["strong"] = strong_dir
    if weak_dir:
        dirs["weak"] = weak_dir

    for label, d in dirs.items():
        records = parse_scaling_files(d)

        # Check 1: at least one record
        if len(records) > 0:
            print(f"PASS  [{label}] Found {len(records)} SCALING_DATA line(s)")
        else:
            print(f"FAIL  [{label}] No SCALING_DATA lines found in {d}")
            all_pass = False
            continue

        # Check 2: all elapsed > 0
        bad_elapsed = [r for r in records if r["elapsed"] <= 0]
        if not bad_elapsed:
            print(f"PASS  [{label}] All elapsed values > 0")
        else:
            print(f"FAIL  [{label}] {len(bad_elapsed)} record(s) with elapsed <= 0")
            all_pass = False

        # Check 3: np positive integers
        bad_np = [r for r in records if r["np"] <= 0]
        if not bad_np:
            print(f"PASS  [{label}] All np values are positive integers")
        else:
            print(f"FAIL  [{label}] {len(bad_np)} record(s) with np <= 0")
            all_pass = False

    return all_pass


def main():
    parser = argparse.ArgumentParser(
        description="Scaling analysis for SEAS BP1 SLURM output files."
    )
    parser.add_argument(
        "--strong-dir",
        type=str,
        default=None,
        help="Directory containing strong scaling .out files (searched recursively)",
    )
    parser.add_argument(
        "--weak-dir",
        type=str,
        default=None,
        help="Directory containing weak scaling .out files (searched recursively)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=".",
        help="Where to save plots and CSV (default: current dir)",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Run validation checks instead of plotting",
    )
    args = parser.parse_args()

    if args.validate:
        ok = validate(args.strong_dir, args.weak_dir)
        sys.exit(0 if ok else 1)

    os.makedirs(args.output_dir, exist_ok=True)

    if args.strong_dir:
        records = parse_scaling_files(args.strong_dir)
        if records:
            strong_scaling_analysis(records, args.output_dir)
        else:
            print(f"WARNING: No SCALING_DATA found in {args.strong_dir}")

    if args.weak_dir:
        records = parse_scaling_files(args.weak_dir)
        if records:
            weak_scaling_analysis(records, args.output_dir)
        else:
            print(f"WARNING: No SCALING_DATA found in {args.weak_dir}")

    if not args.strong_dir and not args.weak_dir:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
