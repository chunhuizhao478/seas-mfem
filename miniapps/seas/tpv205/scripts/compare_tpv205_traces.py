#!/usr/bin/env python3
"""
compare_tpv205_traces.py — A/B comparison harness for the TPV205
through-bi-material verification (plan §R.4 Detailed Req 3,
PLAN_phase_R_exact_bimaterial_riemann_rev3.md).

For each on-fault station, loads the station's .dat trace from a
NEW seas_spatial_dyn_driver run and from a REFERENCE source (either
the pre-Phase-R `tpv205_driver` gold reference or the SCEC DRDG3D
trace bank), then computes per-field:

   rms = sqrt(mean((new - ref)^2))
   peak = max(|new - ref|)

For each (station, field) pair we report rms and peak relative to
peak(ref).  Pass criteria:

  Per the plan §R.4 Detailed Req 3:
    rms  < 5e-3 * peak(ref)  AND  peak < 1e-2 * peak(ref)
    on every station per field (in the "gold" comparison).
  A second pass against the SCEC DRDG3D reference uses a wider band
  (typically 5-10% relative peak error) because DRDG3D is a different
  code.

Usage:
  python compare_tpv205_traces.py \\
      --new      tpv205/out/results            \\
      --gold     tpv205/gold/results_*/results \\
      --tol-rms  5e-3 \\
      --tol-peak 1e-2

  python compare_tpv205_traces.py \\
      --new       tpv205/out/results \\
      --reference tpv205/benchmark_data/DRDG3D_200m_O4 \\
      --tol-rms   5e-2 \\
      --tol-peak  1e-1
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import Dict, List, Tuple

import numpy as np


# Per-station fields we compare.  Names follow the existing
# tpv205_driver.cpp station .dat schema (one column per field).
FIELDS = [
    "slip-rate-1",
    "slip-rate-2",
    "slip-1",
    "slip-2",
    "traction-1",
    "traction-2",
]


def load_dat(path: str) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Load a SCEC-format on-fault station .dat file.

    The header is parsed for column names; the data rows are read
    via numpy.loadtxt.  Returns (time, {field_name: column}).
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    with open(path) as fh:
        lines = fh.readlines()
    header_line = None
    for line in lines:
        if line.startswith("# t") or line.lower().startswith("#t"):
            header_line = line.lstrip("#").strip()
            break
    if header_line is None:
        # Fallback: positional t-column-0 + FIELDS in order.
        cols = ["t"] + FIELDS
    else:
        cols = header_line.replace(",", " ").split()
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.shape[1] < len(cols):
        raise RuntimeError(
            f"{path}: expected ≥{len(cols)} columns, got {data.shape[1]}")
    out = {}
    for i, name in enumerate(cols):
        out[name] = data[:, i]
    return out["t"], out


def resample_to(t_ref: np.ndarray, t_src: np.ndarray,
                y_src: np.ndarray) -> np.ndarray:
    """Linear resample y_src(t_src) onto t_ref grid."""
    return np.interp(t_ref, t_src, y_src)


def station_files(dirpath: str) -> List[str]:
    """Return all SCEC-style station .dat files in `dirpath`."""
    return sorted(
        glob.glob(os.path.join(dirpath, "faultst*.dat"))
        + glob.glob(os.path.join(dirpath, "*station*.dat"))
        + glob.glob(os.path.join(dirpath, "x2_*x3_*.dat"))
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--new", required=True,
                    help="directory of new-run station .dat files")
    ap.add_argument("--gold",
                    help="reference directory (per-station .dat)")
    ap.add_argument("--reference",
                    help="alt reference (e.g. DRDG3D), same layout")
    ap.add_argument("--tol-rms", type=float, default=5e-3,
                    help="rms tolerance, relative to peak(ref)")
    ap.add_argument("--tol-peak", type=float, default=1e-2,
                    help="peak tolerance, relative to peak(ref)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    ref_dir = args.gold or args.reference
    if ref_dir is None:
        print("error: --gold or --reference must be provided",
              file=sys.stderr)
        return 2

    new_files = station_files(args.new)
    if not new_files:
        print(f"error: no station .dat files in {args.new}",
              file=sys.stderr)
        return 2

    n_fail = 0
    n_total = 0
    for new_path in new_files:
        base = os.path.basename(new_path)
        ref_path = os.path.join(ref_dir, base)
        if not os.path.isfile(ref_path):
            print(f"WARN: no reference for {base}; skipping")
            continue
        try:
            t_new, fields_new = load_dat(new_path)
            t_ref, fields_ref = load_dat(ref_path)
        except Exception as e:
            print(f"ERROR loading {base}: {e}")
            n_fail += 1
            continue

        print(f"\n[{base}]  (n_new={len(t_new)}, n_ref={len(t_ref)})")
        for field in FIELDS:
            if field not in fields_ref or field not in fields_new:
                continue
            n_total += 1
            ref = fields_ref[field]
            new = resample_to(t_ref, t_new, fields_new[field])
            peak_ref = float(np.max(np.abs(ref)))
            if peak_ref == 0.0:
                # Field is identically zero in reference — only pass if
                # new is also (approximately) zero.
                if np.max(np.abs(new)) < args.tol_peak:
                    status = "PASS"
                else:
                    status = "FAIL"
                    n_fail += 1
                print(f"  {field:18s}  peak_ref=0 -> {status}")
                continue
            rms_err = float(np.sqrt(np.mean((new - ref) ** 2)))
            peak_err = float(np.max(np.abs(new - ref)))
            rms_rel = rms_err / peak_ref
            peak_rel = peak_err / peak_ref
            ok = rms_rel < args.tol_rms and peak_rel < args.tol_peak
            if not ok:
                n_fail += 1
            status = "PASS" if ok else "FAIL"
            if not args.quiet or not ok:
                print(f"  {field:18s}  rms_rel={rms_rel:.2e}  "
                      f"peak_rel={peak_rel:.2e}  {status}")

    print(f"\n=== {n_total - n_fail} / {n_total} (field, station) pairs "
          f"within tol (rms<{args.tol_rms:.0e}, peak<{args.tol_peak:.0e}) ===")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
