"""probe_diff.py — TPV104 Phase-3 probe-output comparison.

Reads two probe files (one from MFEM, one from the reference FVW
runtime) in the format emitted by
``dynamic/tpv104_substep_iterator.cpp`` under
``SEAS_DIAG_TPV104_STATE``:

    # probe=<name>  code=<MFEM|reference|...>  rank=<r>  nprocs=<n>
    # <optional column-name comment>
    t  qp_id  <field_1>  <field_2>  ...  <field_K>
    <numeric rows>

After parsing, applies ``tpv104_column_map.apply_reference_to_mfem``
to the reference-side fields so their channel names and signs align
with MFEM's.  Returns per-channel max / rel / absolute diffs.

Usage (CLI):

    python3 probe_diff.py --probe=trial_traction \\
        --mfem mfem_probe_trial_traction_rank0.txt \\
        --ref  ref_probe_trial_traction_rank0.txt \\
        [--max-rel 1e-3]

Exit code: 0 if every channel's max rel diff is below ``--max-rel``
(default 1.0 = off).  Non-zero on mismatch.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

import tpv104_column_map as colmap


# Probe → (list of extra-column names in emission order).
# Must stay in sync with the iterator's emission block.
PROBE_COLUMNS: Dict[str, List[str]] = {
    "trial_traction":   ["sigma_n_trial", "tau1_trial", "tau2_trial"],
    "state_evolution":  ["psi_in", "V", "L", "dt", "V_w", "a", "b",
                         "V0", "f0", "muW", "psi_out"],
    "friction_coeff":   ["V", "psi", "a", "V0", "mu"],
    "slip_rate":        ["Theta", "psi", "sigma_n", "eta_s", "a", "V_abs"],
    "corrected_imposed": ["sigma_n_corr", "tau1_corr", "tau2_corr",
                          "Q_imp_plus_vn",  "Q_imp_plus_vt1",  "Q_imp_plus_vt2",
                          "Q_imp_minus_vn", "Q_imp_minus_vt1", "Q_imp_minus_vt2"],
}


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def parse_probe_file(path: str, probe_name: str
                     ) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    """Return (t, qp_id, {channel: array}) for the probe file at ``path``.

    The probe file's first non-comment line starts with ``t qp_id ...``;
    comment lines begin with ``#``.  Rows after that are numeric.
    """
    if probe_name not in PROBE_COLUMNS:
        raise ValueError(f"unknown probe_name: {probe_name!r}")
    column_names = PROBE_COLUMNS[probe_name]
    n_cols = 2 + len(column_names)      # t + qp_id + fields

    rows: List[List[float]] = []
    with open(path, "r") as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) != n_cols:
                raise ValueError(
                    f"{path}:{ln}: expected {n_cols} columns "
                    f"(t, qp_id, {len(column_names)} fields) for probe "
                    f"{probe_name!r}; got {len(parts)} ({parts!r})"
                )
            rows.append([float(x) for x in parts])

    if not rows:
        raise ValueError(f"{path}: no numeric rows found")

    data = np.array(rows, dtype=float)
    t      = data[:, 0]
    qp_id  = data[:, 1].astype(int)
    fields = {name: data[:, 2 + i] for i, name in enumerate(column_names)}
    return t, qp_id, fields


def align_by_t_qp(t_a: np.ndarray, qp_a: np.ndarray, fa: Dict[str, np.ndarray],
                  t_b: np.ndarray, qp_b: np.ndarray, fb: Dict[str, np.ndarray],
                  t_tol: float = 1e-9
                  ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Intersect rows by (t rounded to tolerance, qp_id).  Returns two
    aligned dicts (same row order, same length).

    Missing rows are dropped silently; use ``probe_summarise.py`` for a
    report of alignment coverage.
    """
    key_a = [(round(t / t_tol), int(q)) for t, q in zip(t_a, qp_a)]
    key_b_idx = {(round(t / t_tol), int(q)): i
                 for i, (t, q) in enumerate(zip(t_b, qp_b))}

    idx_a: List[int] = []
    idx_b: List[int] = []
    for ia, k in enumerate(key_a):
        ib = key_b_idx.get(k)
        if ib is not None:
            idx_a.append(ia)
            idx_b.append(ib)

    aligned_a = {k: v[idx_a] for k, v in fa.items()}
    aligned_b = {k: v[idx_b] for k, v in fb.items()}
    return aligned_a, aligned_b


# ---------------------------------------------------------------------------
# Diff
# ---------------------------------------------------------------------------
def diff_probes(probe_name: str,
                mfem_path: str,
                ref_path: str,
                apply_map: bool = True
                ) -> Dict[str, Dict[str, float]]:
    """Diff two probe files.  Returns
        { channel_name : {max_abs, max_rel, mean_abs, count} }.
    """
    t_m, qp_m, fm = parse_probe_file(mfem_path, probe_name)
    t_r, qp_r, fr = parse_probe_file(ref_path,  probe_name)

    if apply_map:
        fr = colmap.apply_reference_to_mfem(probe_name, fr)

    fm_al, fr_al = align_by_t_qp(t_m, qp_m, fm, t_r, qp_r, fr)
    if not fm_al:
        raise ValueError("no aligned rows between MFEM and reference")

    out: Dict[str, Dict[str, float]] = {}
    for name, a in fm_al.items():
        if name not in fr_al:
            continue  # reference file lacks this channel after mapping
        b = fr_al[name]
        n = min(len(a), len(b))
        if n == 0:
            continue
        diff  = np.abs(a[:n] - b[:n])
        denom = np.maximum(np.abs(b[:n]), 1e-30)
        rel   = diff / denom
        out[name] = {
            "count":    int(n),
            "max_abs":  float(np.max(diff)),
            "max_rel":  float(np.max(rel)),
            "mean_abs": float(np.mean(diff)),
        }
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--probe", required=True,
                   choices=sorted(PROBE_COLUMNS.keys()))
    p.add_argument("--mfem",  required=True, help="MFEM probe file")
    p.add_argument("--ref",   required=True, help="Reference probe file")
    p.add_argument("--no-map", action="store_true",
                   help="skip column-name swap + sign flip "
                        "(for self-diff / identity checks)")
    p.add_argument("--max-rel", type=float, default=1.0,
                   help="fail if any channel exceeds this rel diff")
    p.add_argument("--json", action="store_true",
                   help="emit JSON instead of human-readable table")
    args = p.parse_args(argv)

    result = diff_probes(args.probe, args.mfem, args.ref,
                         apply_map=not args.no_map)

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"probe: {args.probe}")
        print(f"{'channel':<25s} {'count':>8s} {'max_abs':>12s} "
              f"{'max_rel':>12s} {'mean_abs':>12s}")
        for name, r in result.items():
            print(f"{name:<25s} {r['count']:>8d} "
                  f"{r['max_abs']:>12.3e} {r['max_rel']:>12.3e} "
                  f"{r['mean_abs']:>12.3e}")

    worst = max((r["max_rel"] for r in result.values()), default=0.0)
    if worst > args.max_rel:
        print(f"FAIL: worst max_rel = {worst:.3e} > {args.max_rel:.3e}",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
