"""probe_summarise.py — scan a directory of TPV104 probe output files
and print a per-probe row-count / column-count summary.

Usage:
    python3 probe_summarise.py path/to/probe-dir

Prints one line per probe file:
    <path>  probe=<name>  rows=<N>  t_range=[<tmin>, <tmax>]  qp_ids={...}
"""

from __future__ import annotations

import os
import re
import sys
from typing import List, Tuple

import numpy as np


_HEADER_RE = re.compile(r"#\s*probe=(\S+)")


def summarise_file(path: str) -> Tuple[str, int, float, float, List[int]]:
    """Return (probe_name, n_rows, tmin, tmax, qp_ids)."""
    probe_name = ""
    rows: List[List[float]] = []
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = _HEADER_RE.match(s)
                if m:
                    probe_name = m.group(1)
                continue
            parts = s.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if not rows:
        return probe_name, 0, 0.0, 0.0, []
    arr = np.array(rows)
    t   = arr[:, 0]
    qp  = arr[:, 1].astype(int)
    return probe_name, arr.shape[0], float(t.min()), float(t.max()), \
           sorted(set(qp.tolist()))


def main(argv: List[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {argv[0]} <probe-dir>", file=sys.stderr)
        return 2
    d = argv[1]
    if not os.path.isdir(d):
        print(f"{d}: not a directory", file=sys.stderr)
        return 2
    any_printed = False
    for name in sorted(os.listdir(d)):
        if not name.startswith("tpv104_probe_"):
            continue
        if not name.endswith(".txt"):
            continue
        path = os.path.join(d, name)
        probe, n, tmin, tmax, qp_ids = summarise_file(path)
        print(f"{path}  probe={probe}  rows={n}  "
              f"t_range=[{tmin:.3e}, {tmax:.3e}]  "
              f"qp_ids={qp_ids}")
        any_printed = True
    if not any_printed:
        print(f"{d}: no tpv104_probe_*.txt files found", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
