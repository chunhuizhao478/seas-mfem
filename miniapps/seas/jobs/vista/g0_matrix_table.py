#!/usr/bin/env python3
"""Tabulate seas_bp5_full --matrix-stats / --exit-after-assembly JSON files
(GPU plan G0: matrix dimensions and nonzero counts per rank)."""
import json
import os
import sys

INT32_MAX = 2**31 - 1

rows = []
for path in sys.argv[1:]:
    try:
        d = json.load(open(path))
    except (OSError, ValueError) as e:
        print(f"unreadable: {path}: {e}")
        continue
    ln = d.get("local_nnz", {})
    rows.append((os.path.basename(path).replace("_matrix_stats.json", ""),
                 d["n_ranks"], d["global_n"], d["global_nnz"],
                 d.get("offd_fraction", float("nan")),
                 ln.get("max", -1), ln.get("mean", float("nan")),
                 d["nnz_per_row"]["max"]))
rows.sort(key=lambda r: (r[0].rsplit("_r", 1)[0], r[1]))
print(f"{'case':<22}{'ranks':>6}{'global n':>14}{'global nnz':>16}"
      f"{'offd':>7}{'max local nnz':>16}{'int32':>7}{'nnz/row max':>13}")
for name, nr, n, nnz, offd, lmax, lmean, rmax in rows:
    fits = "yes" if 0 <= lmax <= INT32_MAX else ("?" if lmax < 0 else "NO")
    print(f"{name:<22}{nr:>6}{n:>14,}{nnz:>16,}{offd:>7.3f}{lmax:>16,}"
          f"{fits:>7}{rmax:>13}")
global_ok = all(r[2] <= INT32_MAX for r in rows)
print(f"global n fits int32 in every case: {'yes' if global_ok else 'NO'}")
