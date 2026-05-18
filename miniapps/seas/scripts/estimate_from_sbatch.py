#!/usr/bin/env python3
"""Estimate output size for one or more sbatch jobs.

Usage:
    python3 estimate_from_sbatch.py <sbatch> [<sbatch> ...]

The script parses the `ibrun ./seas_<driver>_driver ...` line out of each
sbatch file, extracts the relevant `--mesh` / `--tfinal` /
`--paraview-*` flags, and runs the estimator on each in turn.

Output format (default: terminal-friendly):
  * 1 sbatch         -> vertical key-value block
  * 2+ sbatch        -> aligned ASCII table
  * --markdown       -> markdown table (useful for docs and PR descriptions)
  * --json           -> JSON list of per-sbatch result dicts

The PLACEHOLDER ZFP-calibration warnings from the underlying estimator
are suppressed by default; pass --verbose to see them.

Flag overrides:
    --tfinal-override <spec>     replace `--tfinal` from the sbatch
                                 (useful when the sbatch sets a long
                                 nominal tfinal but the SLURM wall is
                                 the binding constraint).
    --markdown                   markdown table output (for docs).
    --json                       JSON output (for piping into other tools).
    --verbose                    show estimator PLACEHOLDER warnings.
    --quiet                      backwards-compatible alias for default.
"""

from __future__ import annotations

import argparse
import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

from estimate_output_size import (
    _filter_name_from_args,
    estimate,
    format_bytes,
    parse_time,
)
from _io_size_mesh import read_gmsh
from _io_size_schedule import ScheduleConfig

REPO_ROOT = THIS_DIR.parent  # miniapps/seas


# ---------------------------------------------------------------------------
# sbatch parser
# ---------------------------------------------------------------------------

# Per-driver fault-tag fallback used ONLY when the mesh has no
# `$PhysicalNames` block (i.e. when auto-detection from the mesh's
# physical-group names fails).  R-402: the previous version hard-coded
# these values per-driver and produced a 180× under-count on
# `bp5_1000m.msh` (which uses tag 100, not 3); auto-detection in
# `_io_size_mesh.read_gmsh` now picks the tag from a
# `Physical Surface(...) "fault"` declaration when present.
FAULT_TAG_FALLBACK_BY_DRIVER = {
    "bp5":    3,
    "tpv102": 3,
    "tpv104": 3,
    "tpv205": 103,
}

# Match the driver executable.  Known suffixes: `_driver` (TPV*) and
# `_full` (BP5 verification).  Capture the short name in between.
DRIVER_RE = re.compile(
    r"seas_(?P<driver>tpv\d+|bp\d+)(?:_driver|_full)\b")


_PRINT_ARITH_RE = re.compile(
    r"\$\(\s*python3?\s+-c\s+['\"]print\(([0-9eE+\-*/.\s]+)\)['\"]\s*\)")


def _stash_subshells(s: str, out: List[str]) -> str:
    """Replace each `$(...)` subshell in `s` with `__SUBSHELL_N__`,
    tracking paren nesting so nested `(` / `)` don't confuse the
    scanner.  The original subshell text is appended to `out`."""
    result = []
    i = 0
    while i < len(s):
        if s[i] == "$" and i + 1 < len(s) and s[i + 1] == "(":
            depth = 1
            j = i + 2
            while j < len(s) and depth > 0:
                if s[j] == "(":
                    depth += 1
                elif s[j] == ")":
                    depth -= 1
                j += 1
            if depth == 0:
                out.append(s[i:j])
                result.append(f"__SUBSHELL_{len(out) - 1}__")
                i = j
                continue
        result.append(s[i])
        i += 1
    return "".join(result)


def _try_eval_print_arithmetic(subshell: str) -> Optional[str]:
    """If `subshell` is `$(python3 -c 'print(<arithmetic>)')` with a
    pure-arithmetic argument, evaluate it and return the result as a
    decimal string.  Otherwise return None."""
    m = _PRINT_ARITH_RE.match(subshell.strip())
    if m is None:
        return None
    # The regex already restricts to digits, e/E, +-*/, dot, whitespace;
    # safe to eval directly.
    expr = m.group(1).strip()
    try:
        return repr(eval(expr, {"__builtins__": {}}, {}))
    except Exception:
        return None


@dataclass
class ParsedSbatch:
    sbatch_path: Path
    driver: str
    mesh: Path
    flags: Dict[str, str]   # raw flag values from the ibrun line


def _join_continuations(text: str) -> List[str]:
    """Collapse trailing-backslash continuations so each logical line
    is a single string."""
    logical = []
    buf = ""
    for line in text.splitlines():
        if line.rstrip().endswith("\\"):
            buf += line.rstrip()[:-1] + " "
        else:
            buf += line
            logical.append(buf)
            buf = ""
    if buf:
        logical.append(buf)
    return logical


def parse_sbatch(path: Path) -> ParsedSbatch:
    text = path.read_text()
    lines = _join_continuations(text)

    # Find the ibrun (or bare `./seas_*`) launch line.  Driver
    # executables are named `seas_tpv102_driver`, `seas_tpv104_driver`,
    # `seas_tpv205_driver`, `seas_bp5_full` (BP5 verification driver),
    # so match either `_driver` or `_full`.
    launch = None
    for line in lines:
        if DRIVER_RE.search(line) and (
                "ibrun" in line or "./seas_" in line):
            stripped = line.strip()
            if stripped.startswith("rm ") or stripped.startswith("make "):
                continue
            launch = line
            break
    if launch is None:
        raise ValueError(
            f"could not find an `ibrun ./seas_*_driver ...` line in {path}")

    # Extract driver name.
    m = DRIVER_RE.search(launch)
    if m is None:
        raise ValueError(f"could not derive driver name from: {launch!r}")
    driver = m.group("driver")

    # Tokenize.  Use shlex to handle quoting.  $(...) subshells can
    # contain nested parens (e.g. `$(python3 -c 'print(123)')`) so a
    # plain regex won't do — scan the string and balance parens.
    subshells: List[str] = []
    launch = _stash_subshells(launch, subshells)
    try:
        tokens = shlex.split(launch)
    except ValueError as e:
        raise ValueError(f"could not shlex-split {path}: {e}")

    # Walk tokens collecting --flag value pairs.
    flags: Dict[str, str] = {}
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t.startswith("--"):
            key = t[2:]
            if i + 1 < len(tokens) and not tokens[i + 1].startswith("--"):
                flags[key] = tokens[i + 1]
                i += 2
            else:
                flags[key] = ""  # boolean flag
                i += 1
        else:
            i += 1

    # Try to recover subshell-evaluated values.  Currently only
    # `$(python3 -c 'print(<arithmetic>)')` is supported; anything
    # else stays as `__SUBSHELL_N__` and the user must use
    # --tfinal-override (or the equivalent for that flag).
    for key, val in list(flags.items()):
        m = re.match(r"__SUBSHELL_(\d+)__", val)
        if m is not None:
            sub = subshells[int(m.group(1))]
            recovered = _try_eval_print_arithmetic(sub)
            if recovered is not None:
                flags[key] = recovered

    if "mesh" not in flags:
        raise ValueError(f"--mesh not found in {path}")
    mesh_str = flags["mesh"]
    # The sbatch path is relative to the Frontera repo root
    # (`/scratch.../miniapps/seas`); on the laptop, resolve relative to
    # the local miniapps/seas directory.
    mesh_path = Path(mesh_str)
    if not mesh_path.is_absolute():
        mesh_path = REPO_ROOT / mesh_path

    return ParsedSbatch(
        sbatch_path=path,
        driver=driver,
        mesh=mesh_path,
        flags=flags,
    )


# ---------------------------------------------------------------------------
# Estimator wrapper
# ---------------------------------------------------------------------------

def _maybe_float(v: Optional[str]) -> Optional[float]:
    if v is None or v == "" or v == "__SUBSHELL__":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _maybe_int(v: Optional[str], default: int = 0) -> int:
    if v is None or v == "":
        return default
    try:
        return int(v)
    except ValueError:
        return default


def estimate_one(p: ParsedSbatch,
                 tfinal_override: Optional[str] = None) -> Dict:
    if not p.mesh.is_file():
        return {
            "label": p.sbatch_path.name,
            "missing_mesh": str(p.mesh),
        }

    # Resolve tfinal.
    tfinal_raw = tfinal_override or p.flags.get("tfinal", "")
    if tfinal_raw == "__SUBSHELL__" or tfinal_raw == "":
        raise ValueError(
            f"tfinal in {p.sbatch_path.name} is a $(...) subshell or "
            f"missing; pass --tfinal-override <spec> to set it.")
    tfinal_seconds = parse_time(tfinal_raw)

    fault_tag_fallback = FAULT_TAG_FALLBACK_BY_DRIVER.get(p.driver, 3)
    # `read_gmsh` auto-detects the fault tag from the mesh's
    # `$PhysicalNames` block when present; the fallback applies only
    # when the mesh has no physical-name block at all.
    mesh = read_gmsh(p.mesh, fault_tag=fault_tag_fallback)

    # Schedule.
    # R-704: forward the per-regime dt overrides into the schedule
    # mirror so the estimator sees the same cadence the C++ driver
    # uses.  `_maybe_float(...) or default` keeps the C++ default
    # when the sbatch does not pass the flag.
    schedule = ScheduleConfig(
        tfinal=tfinal_seconds,
        max_total_snapshots=_maybe_int(p.flags.get("paraview-max-snapshots"), 0),
        fixed_dt=_maybe_float(p.flags.get("paraview-dt")) or 0.0,
        dt_coseismic=(
            _maybe_float(p.flags.get("paraview-coseismic-dt")) or 0.01),
        dt_nucleation=(
            _maybe_float(p.flags.get("paraview-nucleation-dt")) or 1.0),
        dt_interseismic=(
            _maybe_float(p.flags.get("paraview-interseismic-dt"))
            or 3.156e7),
        output_every_n_steps=_maybe_int(p.flags.get("paraview-every"), 0),
        n_events=0,
        driver=p.driver,
    )

    # Filters.
    fault_tol = _maybe_float(p.flags.get("paraview-fault-zfp-tol"))
    fault_def = _maybe_int(p.flags.get("paraview-fault-deflate-level"), -1)
    volume_tol = _maybe_float(p.flags.get("paraview-volume-zfp-tol"))
    volume_def = _maybe_int(p.flags.get("paraview-volume-deflate-level"), -1)
    fault_filter = _filter_name_from_args(
        fault_tol, fault_def if fault_def >= 0 else None, "hdf5")
    volume_filter = _filter_name_from_args(
        volume_tol, volume_def if volume_def >= 0 else None, "hdf5")

    paraview_enabled = (
        "paraview" in p.flags
        or "paraview-adaptive" in p.flags
        or "paraview-dt" in p.flags
        or "paraview-every" in p.flags)
    volume_pv_enabled = "no-volume-pv" not in p.flags
    volume_pv_dt = _maybe_float(p.flags.get("volume-pv-dt")) or 0.0

    report = estimate(
        driver=p.driver,
        mesh=mesh,
        order=_maybe_int(p.flags.get("order"), 1),
        schedule=schedule,
        filters={"fault": fault_filter, "volume": volume_filter},
        scratch_quota_bytes=None,
        paraview_enabled=paraview_enabled,
        volume_pv_enabled=volume_pv_enabled,
        volume_pv_dt=volume_pv_dt,
        explain=False,
        station_count=7 if p.driver == "bp5" else 0,
    )
    return {
        "label":            p.sbatch_path.name,
        "driver":           p.driver,
        "mesh_name":        p.mesh.name,
        "tfinal":           tfinal_raw,
        "n_elements":       mesh.n_elements,
        "n_fault_faces":    mesh.n_fault_faces,
        "n_writes_fault":   report.n_writes_fault,
        "n_writes_volume":  report.n_writes_volume,
        "per_write_fault":  report.per_write_bytes_fault,
        "per_write_volume": report.per_write_bytes_volume,
        "fault_bytes":      report.fault_bytes,
        "volume_bytes":     report.volume_bytes,
        "total_bytes":      report.total_bytes,
        "fault_filter":     fault_filter,
        "volume_filter":    volume_filter,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _format_tfinal(tfinal_raw: str) -> str:
    """Return tfinal with a human-readable suffix where useful."""
    try:
        seconds = float(tfinal_raw)
    except (TypeError, ValueError):
        return str(tfinal_raw)
    # >1 year -> show years; >1 day -> days; >1 hour -> hours; else seconds
    if seconds >= 365.25 * 86400:
        return f"{seconds:g} s  ({seconds / (365.25 * 86400):.2f} yr)"
    if seconds >= 86400:
        return f"{seconds:g} s  ({seconds / 86400:.2f} day)"
    if seconds >= 3600:
        return f"{seconds:g} s  ({seconds / 3600:.2f} h)"
    return f"{seconds:g} s"


def print_pretty_single(r):
    """Vertical key-value layout for a single sbatch."""
    if "missing_mesh" in r:
        print(f"  {r['label']}")
        print(f"  ERROR: {r['missing_mesh']}")
        return
    label = r["label"]
    print(f"  Sbatch       {label}")
    print(f"  Driver       {r['driver']}")
    print(f"  tfinal       {_format_tfinal(r['tfinal'])}")
    print(f"  Mesh         {r['mesh_name']}  "
          f"({r['n_elements']:,} tets, {r['n_fault_faces']:,} fault faces)")
    print(f"  Writes       {r['n_writes_fault']} fault, "
          f"{r['n_writes_volume']} volume")
    print(f"  Filter       fault={r['fault_filter']}  "
          f"volume={r['volume_filter']}")
    print()
    print(f"  Per-write    fault  = {format_bytes(r['per_write_fault']):>10}"
          f"     volume = {format_bytes(r['per_write_volume']):>10}")
    print(f"  Run totals   fault  = {format_bytes(r['fault_bytes']):>10}"
          f"     volume = {format_bytes(r['volume_bytes']):>10}")
    rule = "  " + "─" * 56
    print(rule)
    print(f"  TOTAL        {format_bytes(r['total_bytes'])}")


def _short_tfinal(tfinal_raw: str) -> str:
    """Compact tfinal for the table view: '250 yr' instead of '7.89e+09 s'."""
    try:
        seconds = float(tfinal_raw)
    except (TypeError, ValueError):
        return str(tfinal_raw)
    if seconds >= 365.25 * 86400:
        return f"{seconds / (365.25 * 86400):.0f} yr"
    if seconds >= 86400:
        return f"{seconds / 86400:.1f} day"
    if seconds >= 3600:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds:g} s"


def print_pretty_table(results):
    """Aligned ASCII table for multi-sbatch runs."""
    headers = [
        "Sbatch", "Driver", "tfinal",
        "Tets", "Fault", "n_writes",
        "fault/wr", "vol/wr",
        "fault tot", "vol tot", "TOTAL",
    ]
    rows = []
    for r in results:
        if "missing_mesh" in r:
            rows.append([r["label"], "ERROR", r["missing_mesh"][:48],
                         "", "", "", "", "", "", "", ""])
            continue
        rows.append([
            r["label"],
            r["driver"],
            _short_tfinal(r["tfinal"]),
            f"{r['n_elements']:,}",
            f"{r['n_fault_faces']:,}",
            f"{r['n_writes_fault']}/{r['n_writes_volume']}",
            format_bytes(r["per_write_fault"]),
            format_bytes(r["per_write_volume"]),
            format_bytes(r["fault_bytes"]),
            format_bytes(r["volume_bytes"]),
            format_bytes(r["total_bytes"]),
        ])
    widths = [max(len(h), max((len(row[i]) for row in rows), default=0))
              for i, h in enumerate(headers)]
    line = "  " + "  ".join(h.ljust(w) for h, w in zip(headers, widths))
    print(line)
    print("  " + "  ".join("-" * w for w in widths))
    for row in rows:
        print("  " + "  ".join(c.ljust(w) for c, w in zip(row, widths)))


def print_markdown_table(results):
    """Original markdown-row layout, useful for pasting into docs."""
    cols = (
        ("Sbatch",           "label"),
        ("Driver",           "driver"),
        ("tfinal",           "tfinal"),
        ("Tets",             "n_elements"),
        ("Fault faces",      "n_fault_faces"),
        ("n_writes",         None),
        ("Per-write fault",  "per_write_fault"),
        ("Per-write volume", "per_write_volume"),
        ("Fault total",      "fault_bytes"),
        ("Volume total",     "volume_bytes"),
        ("**TOTAL**",        "total_bytes"),
    )
    print("| " + " | ".join(c[0] for c in cols) + " |")
    print("|" + "|".join(["---"] * len(cols)) + "|")
    for r in results:
        if "missing_mesh" in r:
            print(f"| {r['label']} | _missing mesh: "
                  f"{r['missing_mesh']}_ |"
                  + "|" * (len(cols) - 2))
            continue
        n_writes = (f"{r['n_writes_fault']} fault / "
                    f"{r['n_writes_volume']} vol")
        cells = [
            r["label"],
            r["driver"],
            r["tfinal"],
            f"{r['n_elements']:,}",
            f"{r['n_fault_faces']:,}",
            n_writes,
            format_bytes(r["per_write_fault"]),
            format_bytes(r["per_write_volume"]),
            format_bytes(r["fault_bytes"]),
            format_bytes(r["volume_bytes"]),
            f"**{format_bytes(r['total_bytes'])}**",
        ]
        print("| " + " | ".join(cells) + " |")


def print_table(results, fmt: str = "pretty"):
    """Dispatcher: vertical for 1 result, table otherwise; respect fmt."""
    if fmt == "markdown":
        print_markdown_table(results)
        return
    if fmt == "json":
        import json
        # bytes/int are JSON-safe; PosixPath etc. are not in our dicts.
        print(json.dumps(results, indent=2, default=str))
        return
    # default: terminal-friendly
    if len(results) == 1:
        print_pretty_single(results[0])
    else:
        print_pretty_table(results)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="estimate_from_sbatch",
        description="Estimate output size for one or more sbatch jobs.")
    parser.add_argument(
        "sbatch", nargs="+", type=Path,
        help="Path to one or more sbatch files.")
    parser.add_argument(
        "--tfinal-override", type=str, default=None,
        help="Override the --tfinal value parsed from each sbatch.  "
             "Required when the sbatch's --tfinal is a $(...) subshell.")
    fmt_group = parser.add_mutually_exclusive_group()
    fmt_group.add_argument(
        "--markdown", action="store_true",
        help="Print a markdown table (useful for docs / PR descriptions).")
    fmt_group.add_argument(
        "--json", dest="json_out", action="store_true",
        help="Print JSON list of per-sbatch result dicts.")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "--verbose", action="store_true",
        help="Show estimator PLACEHOLDER calibration warnings.")
    verbosity.add_argument(
        "--quiet", action="store_true",
        help="Suppress PLACEHOLDER warnings (default; kept for backwards "
             "compatibility).")
    args = parser.parse_args(argv)

    # Default is quiet — only --verbose brings the warnings back.
    if not args.verbose:
        import warnings
        warnings.simplefilter("ignore")

    results = []
    for sb in args.sbatch:
        try:
            parsed = parse_sbatch(sb)
            results.append(estimate_one(parsed, args.tfinal_override))
        except Exception as e:
            results.append({"label": sb.name, "missing_mesh": str(e)})

    if args.markdown:
        fmt = "markdown"
    elif args.json_out:
        fmt = "json"
    else:
        fmt = "pretty"
    print_table(results, fmt=fmt)

    # Calibration footer.  Skip for json output so the stdout stays
    # machine-parseable.
    if fmt != "json" and any("missing_mesh" not in r for r in results):
        print()
        print("  Note: ZFP compression ratios are PLACEHOLDERs pending "
              "Phase 7.6 calibration;")
        print("        estimates may be off by >2x.  "
              "Pass --verbose to see per-field warnings.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
