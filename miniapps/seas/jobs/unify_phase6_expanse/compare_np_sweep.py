#!/usr/bin/env python3
"""Cross-np station-trace comparison for the unify-plan Phase 6 np sweep.

Usage:  python3 compare_np_sweep.py <ref_dir> <other_dir> [<other_dir> ...]

Each directory is one np run's --output-dir containing `*_station_*.dat`
trace files (any prefix; matched by the `_station_<name>.dat` suffix).  For
every station present in the reference directory, prints the maximum
absolute difference over the full time series against each other
directory, plus the relative-to-reference measure.

Exit status: 0 when every station in every comparison is within TOL
(default 1e-9 absolute OR 1e-9 relative — print precision of the trace
files; override with UNIFY_CMP_TOL), 1 otherwise.  Missing stations or
row-count mismatches (different dt/steps — a real red flag for the
unified path) are failures.

Standard library + numpy only.
"""

import glob
import os
import sys

import numpy as np

TOL = float(os.environ.get("UNIFY_CMP_TOL", "1e-9"))


def station_map(d):
    out = {}
    for f in glob.glob(os.path.join(d, "*_station_*.dat")):
        name = os.path.basename(f).split("_station_", 1)[1][:-len(".dat")]
        out[name] = f
    return out


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    ref_dir, other_dirs = argv[1], argv[2:]
    ref = station_map(ref_dir)
    if not ref:
        print(f"FAIL: no station files found in reference dir {ref_dir}")
        return 1

    failed = False
    print(f"reference: {ref_dir}  ({len(ref)} stations)   tol={TOL:g}")
    header = f"{'station':>24} " + " ".join(f"{os.path.basename(d.rstrip('/')):>14}"
                                            for d in other_dirs)
    print(header)
    for name in sorted(ref):
        a = np.loadtxt(ref[name])
        row = [f"{name:>24}"]
        for d in other_dirs:
            other = station_map(d)
            if name not in other:
                row.append(f"{'MISSING':>14}")
                failed = True
                continue
            b = np.loadtxt(other[name])
            if a.shape != b.shape:
                row.append(f"{'SHAPE ' + str(b.shape):>14}")
                failed = True
                continue
            dmax = float(np.abs(a - b).max())
            rel = float((np.abs(a - b)
                         / np.maximum(np.abs(a), 1.0)).max())
            ok = dmax < TOL or rel < TOL
            row.append(f"{dmax:>13.2e}{'' if ok else '*'}")
            if not ok:
                failed = True
        print(" ".join(row))

    print()
    if failed:
        print(f"FAIL: at least one station exceeds tol={TOL:g} (marked *) "
              "or is missing/misshapen.")
        print("      For the unified substep path + deterministic station "
              "tie-break this indicates a real cross-np divergence — see "
              "PLAN_unify_interior_shared_fault_substep_2026-07-09.md "
              "Phase 6 (the STOP condition).")
        return 1
    print("PASS: all stations np-independent within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
