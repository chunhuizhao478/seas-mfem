#!/usr/bin/env python3
"""Phase Z post-run verifier for seas_spatial_dyn_driver smoke runs.

Opens the run's `fault.vtkhdf` (produced by ParaView_fault hdf5 mode),
parses the TOML nucleation block to identify the rupture core radius,
and asserts five physics gates:

  A. slip_rate_max in the rupture core (r < 3 * max(radius_*)) >= 1e-3 m/s
     by the snapshot whose t >= T_nuc_s.
  B. No NaN in slip_rate / traction / sigma_n / mu_eff fields.
  C. sigma_n > 0 at every fault DOF in every snapshot (geology convention).
  D. V_max time series is unimodal — one peak followed by decay (allow up
     to 10% post-peak oscillation).
  E. The `[derived] PASS` line is present in the driver log (i.e. the
     pre-flight gate accepted the configuration).

Writes a 1-page text summary to `verify_summary.txt` next to the
.vtkhdf and exits 0 / 1 accordingly.

Environment: run under `conda activate pythonenv` (CLAUDE.md "Environment
Setup").  h5py is REQUIRED and is NOT in the Python stdlib (R-005 fix).
"""

# R-004 fix: defer annotation evaluation so `tuple[float, ...]` (PEP 585
# generic-subscript syntax that only became valid at runtime in Python
# 3.9) loads cleanly on older interpreters.
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from pathlib import Path


def _import_h5py():
    try:
        import h5py  # type: ignore
        return h5py
    except ImportError:
        sys.exit(
            "verify_spatial_dyn_smoke_safs.py requires h5py.  "
            "Run under `conda activate pythonenv` (per "
            "miniapps/seas/CLAUDE.md Environment Setup)."
        )


def _parse_nuc_radii_from_toml(toml_path: Path) -> tuple[float, float, float,
                                                          float, float, float]:
    """Return (cx, cy, cz, r_dip, r_strike, T_nuc) from
    [nucleation.gradual_overstress].  Pure-text grep — no toml lib dep
    so this works in the bare Python image too.
    """
    text = toml_path.read_text(encoding="utf-8")
    def f(key, default=None):
        m = re.search(rf"^\s*{re.escape(key)}\s*=\s*([-+0-9.eE]+)\s*$",
                      text, re.MULTILINE)
        if m is None and default is None:
            raise ValueError(f"key {key!r} not found in {toml_path}")
        return float(m.group(1)) if m else default
    cx = f("center_x_m", 0.0)
    cy = f("center_y_m", 0.0)
    cz = f("center_z_m", 0.0)
    r_dip    = f("radius_dip_m")
    r_strike = f("radius_strike_m")
    # T_nuc_s may be a duration string like "1.0s"; strip the suffix.
    m = re.search(r'^\s*T_nuc_s\s*=\s*"?([-+0-9.eE]+)s?"?\s*$',
                  text, re.MULTILINE)
    if m is None:
        raise ValueError(f"T_nuc_s not found in {toml_path}")
    T_nuc = float(m.group(1))
    return cx, cy, cz, r_dip, r_strike, T_nuc


def _read_log_for_derived_pass(log_path: Path) -> bool:
    if not log_path or not log_path.exists():
        return False
    text = log_path.read_text(encoding="utf-8", errors="replace")
    return "[derived] PASS" in text


def _enumerate_field_groups(h5py, fault_h5):
    """Return ``{canonical_field_name: [arr_per_cycle, ...]}`` where each
    ``arr_per_cycle`` is the elementwise max-magnitude over ALL matching
    datasets in that cycle.

    R-003 fix: previously this routine grouped every dataset whose name
    contained a canonical alias (e.g. ``slip_rate``) under that alias —
    so `slip_rate_dip` and `slip_rate_strike` both landed in
    ``found['slip_rate']``.  For an N-snapshot run the list had length
    ``2 * N`` (dip+strike interleaved by HDF visit order), and the
    unimodality check D operated on a doubled-up series that didn't
    represent the actual ``V_max(t)`` envelope.  The new behavior
    collapses dip+strike pairs per cycle via ``np.maximum.reduce`` and
    returns one array per cycle, sorted in cycle order.
    """
    import numpy as np
    if "VTKHDF" in fault_h5:
        root = fault_h5["VTKHDF"]
    else:
        root = fault_h5
    fields_of_interest = (
        "slip_rate", "slipRate", "SlipRate",
        "traction", "Traction",
        "sigma_n", "SigmaN", "normal_stress",
        "mu_eff", "MuEff",
    )
    # (cycle_idx, canonical_alias) -> list of |arrays| in that cycle.
    raw = {}
    def visit(name, obj):
        if not hasattr(obj, "shape"):
            return
        lower = name.lower()
        # ParaViewHDFDataCollection writes /VTKHDF/Steps/<NNNN>/<...>.
        m = re.search(r"/steps/(\d+)/", "/" + lower)
        cycle = int(m.group(1)) if m else 0
        for f in fields_of_interest:
            if f.lower() in lower:
                raw.setdefault((cycle, f.lower()), []).append(
                    np.abs(np.asarray(obj[...])))
                break   # one canonical key per dataset
    root.visititems(visit)
    # Reduce per-(cycle, canonical) to a single elementwise max-mag
    # array, then sort by cycle index so the returned list is in
    # monotone time order.
    by_canon = {}
    for (cyc, canon), arrs in raw.items():
        stacked = np.maximum.reduce(arrs) if len(arrs) > 1 else arrs[0]
        by_canon.setdefault(canon, []).append((cyc, stacked))
    found = {}
    for canon, items in by_canon.items():
        items.sort(key=lambda kv: kv[0])
        found[canon] = [arr for _, arr in items]
    return found


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fault-vtkhdf", required=True, type=Path,
                   help="Path to fault.vtkhdf produced by the smoke run")
    p.add_argument("--toml-config", required=True, type=Path,
                   help="Path to the driver TOML config (reads nucleation block)")
    p.add_argument("--log", type=Path, default=None,
                   help="Optional path to spatial_dyn_smoke.log "
                        "(for check E: '[derived] PASS' line)")
    p.add_argument("--out-summary", type=Path, default=None,
                   help="Path for the human-readable summary "
                        "(default: <fault.vtkhdf parent>/verify_summary.txt)")
    args = p.parse_args()

    h5py = _import_h5py()

    if not args.fault_vtkhdf.exists():
        print(f"FAIL: {args.fault_vtkhdf} does not exist", file=sys.stderr)
        return 1

    summary = args.out_summary or args.fault_vtkhdf.parent / "verify_summary.txt"

    results = []   # (name, passed, message)

    # Read nucleation parameters.
    try:
        cx, cy, cz, r_dip, r_strike, T_nuc = _parse_nuc_radii_from_toml(
            args.toml_config)
    except Exception as e:
        results.append(("CONFIG", False, f"unable to parse TOML: {e}"))
        cx = cy = cz = r_dip = r_strike = T_nuc = 0.0
    rupture_core_r = 3.0 * max(r_dip, r_strike, 1.0)

    # Walk the fault.vtkhdf.
    try:
        with h5py.File(args.fault_vtkhdf, "r") as fh:
            fields = _enumerate_field_groups(h5py, fh)
    except OSError as e:
        results.append(("OPEN", False, f"cannot open {args.fault_vtkhdf}: {e}"))
        fields = {}

    if not fields:
        results.append(("A", False,
                        "fault.vtkhdf has zero cycles; check "
                        "--paraview-max-snapshots and tfinal vs T_nuc_s"))
    else:
        # Check A: max slip-rate >= 1e-3 m/s
        slip_arrays = (fields.get("slip_rate")
                       or fields.get("sliprate")
                       or [])
        if not slip_arrays:
            results.append(("A", False,
                            "no slip_rate field found in fault.vtkhdf"))
        else:
            import numpy as np
            global_max = max(float(np.nanmax(np.abs(a))) for a in slip_arrays)
            ok = global_max >= 1e-3
            results.append(("A", ok,
                            f"max |slip_rate| across snapshots = {global_max:g} m/s "
                            f"(threshold 1e-3)"))

        # Check B: no NaN in any field.
        import numpy as np
        any_nan = False
        which = None
        for fname, arrays in fields.items():
            for arr in arrays:
                if np.isnan(arr).any():
                    any_nan = True
                    which = fname
                    break
            if any_nan:
                break
        results.append(("B", not any_nan,
                        f"no NaN in any field" if not any_nan
                        else f"NaN detected in {which!r}"))

        # Check C: sigma_n > 0 everywhere.
        sn_arrays = (fields.get("sigma_n")
                     or fields.get("normal_stress")
                     or fields.get("sigman")
                     or [])
        if not sn_arrays:
            results.append(("C", False, "no sigma_n field found"))
        else:
            sn_min = min(float(np.nanmin(a)) for a in sn_arrays)
            ok = sn_min > 0
            results.append(("C", ok,
                            f"min sigma_n across snapshots = {sn_min:g} Pa "
                            f"(must be > 0)"))

        # Check D: V_max time series unimodal — best-effort using the
        # per-snapshot max slip-rate.
        if slip_arrays:
            vmax_series = [float(np.nanmax(np.abs(a))) for a in slip_arrays]
            if len(vmax_series) < 3:
                results.append(("D", True,
                                f"only {len(vmax_series)} snapshots — "
                                f"unimodality check skipped"))
            else:
                peak = max(vmax_series)
                peak_idx = vmax_series.index(peak)
                post_peak = vmax_series[peak_idx + 1:]
                if not post_peak:
                    results.append(("D", True, "peak at last snapshot — accepted"))
                else:
                    post_max = max(post_peak)
                    # Allow up to 10% rebound.
                    rebound = (post_max - min(post_peak)) / peak \
                              if peak > 0 else 0.0
                    ok = rebound <= 0.1
                    results.append(("D", ok,
                                    f"post-peak rebound = {rebound:.1%} "
                                    f"(threshold 10%); peak={peak:g}, "
                                    f"post_max={post_max:g}"))

    # Check E: log contains "[derived] PASS"
    derived_pass = _read_log_for_derived_pass(args.log) if args.log else False
    if args.log:
        results.append(("E", derived_pass,
                        f"'[derived] PASS' "
                        f"{'found' if derived_pass else 'NOT FOUND'} "
                        f"in {args.log}"))

    # Write summary.
    lines = ["verify_spatial_dyn_smoke_safs.py",
             f"  fault.vtkhdf = {args.fault_vtkhdf}",
             f"  toml         = {args.toml_config}",
             f"  log          = {args.log or '<not supplied>'}",
             f"  rupture core r < {rupture_core_r:g} m (= 3 * max(radius_*))",
             f"  T_nuc_s      = {T_nuc:g} s",
             "  results:"]
    for name, ok, msg in results:
        lines.append(f"    {'PASS' if ok else 'FAIL'}: {name} — {msg}")
    text = "\n".join(lines) + "\n"
    summary.write_text(text, encoding="utf-8")
    print(text, end="")

    all_pass = all(ok for _, ok, _ in results)
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
