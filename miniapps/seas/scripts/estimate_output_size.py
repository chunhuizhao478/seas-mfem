#!/usr/bin/env python3
# Phase 7.5 of miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md.
#
# Pre-submit output-size estimator for SEAS runs.  Standard-library
# only — runs on a laptop in <5 s, no MFEM dependency.
#
# Usage:
#   python3 estimate_output_size.py \
#       --driver bp5 \
#       --mesh bp5/mesh/bp5_1000m.msh \
#       --tfinal 250yr \
#       --paraview \
#       --paraview-fault-zfp-tol 1e-12 \
#       --paraview-max-snapshots 5000 \
#       --no-volume-pv \
#       --np 800 \
#       --scratch-quota 1TB

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional

# Sibling-module imports (this file is intended to be run as a
# script from the scripts/ directory, so use plain imports).
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _io_size_compression import bytes_per_dof, classify_field
from _io_size_mesh import MeshSummary, read_gmsh, read_inline_mesh
from _io_size_schedule import (
    ScheduleConfig,
    estimate_n_volume_writes,
    estimate_n_writes,
)
from _io_size_schemas import FieldSchema, estimate_dofs, get_driver_schema


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

_TIME_SUFFIX = {
    "":     1.0,
    "s":    1.0,
    "min":  60.0,
    "hr":   3600.0,
    "h":    3600.0,
    "day":  86_400.0,
    "d":    86_400.0,
    "wk":   86_400.0 * 7,
    "w":    86_400.0 * 7,
    "yr":   3.156e7,
    "y":    3.156e7,
}

_TIME_RE = re.compile(
    r"^\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*"
    r"(s|min|hr|h|day|d|wk|w|yr|y)?\s*$")

_BYTE_SUFFIX = {
    "":   1,
    "B":  1,
    "K":  1024, "KB": 1024, "KIB": 1024,
    "M":  1024 ** 2, "MB": 1024 ** 2, "MIB": 1024 ** 2,
    "G":  1024 ** 3, "GB": 1024 ** 3, "GIB": 1024 ** 3,
    "T":  1024 ** 4, "TB": 1024 ** 4, "TIB": 1024 ** 4,
    "P":  1024 ** 5, "PB": 1024 ** 5, "PIB": 1024 ** 5,
}


def parse_time(spec: str) -> float:
    """Parse `'250yr'`, `'1.5day'`, `'3.156e7'` (default seconds).

    Supports decimal values and whitespace per plan §Phase 7.5
    Detailed Requirements step 1.  Unknown suffixes raise
    `ValueError` listing the supported set."""
    m = _TIME_RE.match(spec)
    if m is None:
        raise ValueError(
            f"unparseable time specification {spec!r}; "
            f"examples: '250yr', '1.5day', '3.156e7'")
    val = float(m.group(1))
    suffix = (m.group(2) or "").lower()
    if suffix not in _TIME_SUFFIX:
        raise ValueError(
            f"unknown time suffix {suffix!r}; supported: "
            f"{sorted(set(_TIME_SUFFIX.keys()) - {''})}")
    return val * _TIME_SUFFIX[suffix]


def parse_bytes(spec: str) -> int:
    """Parse `'1TB'`, `'500GB'`, `'4096'` (default bytes)."""
    m = re.match(r"^\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*([A-Za-z]*)\s*$",
                 spec)
    if m is None:
        raise ValueError(f"unparseable byte size {spec!r}")
    val = float(m.group(1))
    suffix = m.group(2).upper()
    if suffix not in _BYTE_SUFFIX:
        raise ValueError(
            f"unknown byte suffix {suffix!r}; "
            f"supported: B, KB, MB, GB, TB, PB (and -iB variants)")
    return int(val * _BYTE_SUFFIX[suffix])


_FORMAT_BYTES_UNITS = ("B", "KB", "MB", "GB", "TB", "PB", "EB")


def format_bytes(n: int) -> str:
    """Return a human-readable byte string with consistent
    3-significant-digit precision (R-312).

    Auto-scales B → KB → MB → GB → TB → PB so the column never
    overflows even on very large runs.  Uses 1024-base units.
    """
    if n < 0:
        return "-" + format_bytes(-n)
    if n == 0:
        return "0 B"
    scale = 0
    val = float(n)
    while val >= 1024.0 and scale < len(_FORMAT_BYTES_UNITS) - 1:
        val /= 1024.0
        scale += 1
    if val >= 100.0:
        return f"{val:.0f} {_FORMAT_BYTES_UNITS[scale]}"
    if val >= 10.0:
        return f"{val:.1f} {_FORMAT_BYTES_UNITS[scale]}"
    return f"{val:.2f} {_FORMAT_BYTES_UNITS[scale]}"


def _filter_name_from_args(zfp_tol: Optional[float],
                            deflate_level: Optional[int],
                            volume_mode: str) -> str:
    """Map (zfp_tol, deflate_level, volume_mode) → filter name as used
    in COMPRESSION_RATIOS."""
    if volume_mode == "vtu":
        return "vtu_binary"
    if zfp_tol is not None and zfp_tol > 0:
        # Snap to the nearest table tol.
        for tol in (1e-12, 1e-9, 1e-6, 1e-3):
            if zfp_tol <= tol:
                return f"zfp_{tol:.0e}".replace("e-0", "e-")
        return "zfp_1e-3"
    level = deflate_level if deflate_level is not None else 6
    level = max(0, min(9, int(level)))
    return f"deflate_{level}"


# -----------------------------------------------------------------------------
# Estimator
# -----------------------------------------------------------------------------

@dataclass
class EstimateReport:
    fault_bytes: int
    volume_bytes: int
    station_bytes: int
    total_bytes: int
    per_write_bytes_fault: int
    per_write_bytes_volume: int
    n_writes_fault: int
    n_writes_volume: int
    quota_fraction: Optional[float] = None
    explain_lines: Optional[List[str]] = None


def estimate(driver: str,
             mesh: MeshSummary,
             order: int,
             schedule: ScheduleConfig,
             filters: Dict[str, str],
             scratch_quota_bytes: Optional[int] = None,
             paraview_enabled: bool = True,
             volume_pv_enabled: bool = True,
             volume_pv_dt: float = 0.0,
             explain: bool = False,
             station_count: int = 7,
             station_columns: int = 8,
             station_rows: Optional[int] = None) -> EstimateReport:
    """Top-level estimator (plan §Phase 7.5 Interfaces)."""
    schema = get_driver_schema(driver)
    explain_lines: List[str] = []

    # Per-write bytes — fault and volume.
    fault_filter = filters.get("fault", "deflate_6")
    volume_filter = filters.get("volume", "deflate_6")

    per_write_bytes_fault = 0
    per_write_bytes_volume = 0
    static_bytes_fault = 0
    static_bytes_volume = 0

    for fld in schema:
        n_dofs = estimate_dofs(fld, mesh, order)
        if fld.is_fault:
            bpd = bytes_per_dof(fld, fault_filter, driver)
            field_bytes = int(n_dofs * fld.n_components * bpd)
            if fld.is_static:
                static_bytes_fault += field_bytes
            else:
                per_write_bytes_fault += field_bytes
            if explain:
                explain_lines.append(
                    f"  fault.{fld.name:<24s} "
                    f"dofs={n_dofs:>10d} "
                    f"comp={fld.n_components} "
                    f"bpd={bpd:.3f} "
                    f"({classify_field(fld.name, driver)}, "
                    f"{fault_filter})  "
                    f"= {format_bytes(field_bytes)}"
                    f"{'  [static]' if fld.is_static else ''}")
        else:  # volume
            bpd = bytes_per_dof(fld, volume_filter, driver)
            field_bytes = int(n_dofs * fld.n_components * bpd)
            if fld.is_static:
                static_bytes_volume += field_bytes
            else:
                per_write_bytes_volume += field_bytes
            if explain:
                explain_lines.append(
                    f"  volume.{fld.name:<23s} "
                    f"dofs={n_dofs:>10d} "
                    f"comp={fld.n_components} "
                    f"bpd={bpd:.3f} "
                    f"({classify_field(fld.name, driver)}, "
                    f"{volume_filter})  "
                    f"= {format_bytes(field_bytes)}"
                    f"{'  [static]' if fld.is_static else ''}")

    # Schedule integration.
    n_writes_fault = estimate_n_writes(schedule)
    n_writes_volume = estimate_n_volume_writes(schedule, volume_pv_dt)

    if not paraview_enabled:
        n_writes_fault = 0
        n_writes_volume = 0

    if not volume_pv_enabled:
        n_writes_volume = 0

    # Total bytes — running cycles + static-once at first save.
    fault_bytes = (n_writes_fault * per_write_bytes_fault
                   + (static_bytes_fault if n_writes_fault > 0 else 0))
    volume_bytes = (n_writes_volume * per_write_bytes_volume
                    + (static_bytes_volume if n_writes_volume > 0 else 0))

    # Stations: SEPARATE benchmark TXT output path, not part of the
    # full PV solution data.  On-fault PV output (`fault_bytes`)
    # already covers every fault face; the per-station probe TXT
    # files are small benchmark CSVs sampled independently from the
    # PV schedule.  Reported as an informational aside; NOT summed
    # into `total_bytes`.
    if station_rows is None:
        # BP5 station TXT writer samples at fault PV cadence.
        n_station_rows = max(1, n_writes_fault) if n_writes_fault else 0
    else:
        n_station_rows = station_rows
    station_bytes = (
        station_count * n_station_rows
        * (station_columns * 25 + 1))

    # Total = full PV solution data only.  Stations excluded.
    total = fault_bytes + volume_bytes

    quota_fraction = None
    if scratch_quota_bytes is not None and scratch_quota_bytes > 0:
        quota_fraction = total / scratch_quota_bytes

    if explain:
        explain_lines.append("")
        explain_lines.append(
            f"  per_write_bytes_fault  = {format_bytes(per_write_bytes_fault)}")
        explain_lines.append(
            f"  per_write_bytes_volume = {format_bytes(per_write_bytes_volume)}")
        explain_lines.append(
            f"  n_writes_fault         = {n_writes_fault}")
        explain_lines.append(
            f"  n_writes_volume        = {n_writes_volume}")
        explain_lines.append(
            f"  static_bytes_fault     = {format_bytes(static_bytes_fault)}")
        explain_lines.append(
            f"  static_bytes_volume    = {format_bytes(static_bytes_volume)}")
        explain_lines.append(
            f"  station_rows           = {n_station_rows}")
        explain_lines.append(
            f"  station_count          = {station_count}")

    return EstimateReport(
        fault_bytes=fault_bytes,
        volume_bytes=volume_bytes,
        station_bytes=station_bytes,
        total_bytes=total,
        per_write_bytes_fault=per_write_bytes_fault,
        per_write_bytes_volume=per_write_bytes_volume,
        n_writes_fault=n_writes_fault,
        n_writes_volume=n_writes_volume,
        quota_fraction=quota_fraction,
        explain_lines=explain_lines if explain else None,
    )


# -----------------------------------------------------------------------------
# Output formatting (TTY-friendly + JSON)
# -----------------------------------------------------------------------------

def _color(stream, code: str, text: str) -> str:
    if not stream.isatty():
        return text
    return f"\033[{code}m{text}\033[0m"


def _quota_color(stream, fraction: float, text: str) -> str:
    if fraction > 1.0:
        return _color(stream, "31", text)  # red
    if fraction > 0.8:
        return _color(stream, "33", text)  # yellow
    return _color(stream, "32", text)      # green


def _print_human(report: EstimateReport,
                 driver: str,
                 mesh: MeshSummary,
                 schedule: ScheduleConfig,
                 filters: Dict[str, str],
                 paraview_enabled: bool,
                 volume_pv_enabled: bool,
                 scratch_quota_bytes: Optional[int],
                 stream=sys.stdout) -> None:
    print(f"Mesh: {Path(mesh.extra.get('source_file', '<inline>')).name} "
          f"— {mesh.n_elements:,} {mesh.element_type or 'cells'}, "
          f"{mesh.n_fault_faces:,} fault faces, "
          f"{mesh.n_boundary_faces:,} boundary faces",
          file=stream)
    schema = get_driver_schema(driver)
    n_fault = sum(1 for f in schema if f.is_fault)
    n_volume = sum(1 for f in schema if f.is_volume)
    print(f"Driver: {driver} ({n_fault} fault fields, {n_volume} volume fields)",
          file=stream)
    print(f"Schedule: cap={schedule.max_total_snapshots} "
          f"over {schedule.tfinal:.3g}s "
          f"-> fault={report.n_writes_fault} writes, "
          f"volume={report.n_writes_volume} writes",
          file=stream)
    print("", file=stream)
    print("Per-write bytes (full PV solution data):", file=stream)
    if paraview_enabled:
        print(f"  fault    ({filters.get('fault', 'deflate_6')})  "
              f"= {format_bytes(report.per_write_bytes_fault):>10s}",
              file=stream)
        if volume_pv_enabled:
            print(f"  volume   ({filters.get('volume', 'deflate_6')})  "
                  f"= {format_bytes(report.per_write_bytes_volume):>10s}",
                  file=stream)
        else:
            print("  volume   (suppressed via --no-volume-pv)", file=stream)
    else:
        print("  fault    (PV disabled — no fault output)", file=stream)
        print("  volume   (PV disabled — no volume output)", file=stream)
    print("", file=stream)
    print("Total run output (PV solution data):", file=stream)
    print(f"  fault    : {format_bytes(report.fault_bytes):>10s}", file=stream)
    print(f"  volume   : {format_bytes(report.volume_bytes):>10s}", file=stream)
    print("  -----------------", file=stream)
    print(f"  TOTAL    : {format_bytes(report.total_bytes):>10s}", file=stream)
    print("", file=stream)
    if report.station_bytes > 0:
        print("Benchmark station TXT (separate path; NOT in TOTAL):",
              file=stream)
        print(f"  stations : {format_bytes(report.station_bytes):>10s}",
              file=stream)
        print("", file=stream)

    if scratch_quota_bytes is not None:
        frac = report.quota_fraction or 0.0
        text = (f"Quota check ($SCRATCH = {format_bytes(scratch_quota_bytes)}): "
                f"{frac * 100:.1f}% of quota — "
                f"{'OVER' if frac > 1.0 else 'WARN' if frac > 0.8 else 'OK'}")
        print(_quota_color(stream, frac, text), file=stream)

    if report.explain_lines is not None:
        print("", file=stream)
        print("--- Explain ---", file=stream)
        for line in report.explain_lines:
            print(line, file=stream)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    # R-503: disable prefix-matching globally — without this argparse's
    # default `allow_abbrev=True` makes `--paraview-dt 0.5` ambiguous
    # because three options share the `--paraview-dt-` prefix.
    p = argparse.ArgumentParser(
        prog="estimate_output_size",
        description="Pre-submit output-size estimator for SEAS runs.",
        allow_abbrev=False)
    p.add_argument("--driver", required=True,
                   choices=["bp5", "tpv102", "tpv104", "tpv205"])
    p.add_argument("--mesh", type=Path, default=None,
                   help="Path to a Gmsh ASCII v2/v4 mesh.")
    p.add_argument("--inline-mesh", action="store_true",
                   help="Use BP5's inline mesh (24x24x12 hex, default).")
    p.add_argument("--inline-nx", type=int, default=24)
    p.add_argument("--inline-ny", type=int, default=24)
    p.add_argument("--inline-nz", type=int, default=12)
    p.add_argument("--inline-element-type", choices=["tet", "hex"],
                   default="tet")
    p.add_argument("--fault-tag", type=int, default=3)
    p.add_argument("--order", type=int, default=1,
                   help="Polynomial order (default 1).")
    p.add_argument("--tfinal", type=str, required=True,
                   help="Final time, e.g. '250yr', '1.5day', '3.156e7' "
                        "(seconds default).")

    # Schedule flags (mirror C++ AdaptiveSchedule).
    p.add_argument("--paraview", action="store_true",
                   help="Enable ParaView output (default off).")
    p.add_argument("--no-volume-pv", action="store_true",
                   help="Disable volume PV output (BP5 production default).")
    p.add_argument("--volume-pv-dt", type=float, default=0.0,
                   help="Fixed time gap between volume PV writes (seconds).")
    p.add_argument("--paraview-max-snapshots", type=int, default=0)
    p.add_argument("--paraview-output-every-n-steps", type=int, default=0)
    p.add_argument("--paraview-fixed-dt", type=float, default=0.0)
    p.add_argument("--paraview-v-coseismic", type=float, default=1e-3)
    p.add_argument("--paraview-v-nucleation", type=float, default=1e-7)
    # R-504: regime-FIRST word order to match the C++ driver
    # (`--paraview-coseismic-dt` etc. in tpv102_driver.cpp:585-587 and
    # the equivalent block in bp5_verification_full.cpp).
    # R-409: defaults match paraview_output.hpp:168-170 (dt_coseismic=0.01,
    # dt_nucleation=1.0, dt_interseismic=1 yr).
    p.add_argument("--paraview-coseismic-dt", type=float, default=0.01,
                   dest="paraview_dt_coseismic")
    p.add_argument("--paraview-nucleation-dt", type=float, default=1.0,
                   dest="paraview_dt_nucleation")
    p.add_argument("--paraview-interseismic-dt", type=float, default=3.156e7,
                   dest="paraview_dt_interseismic")

    # Volume-side compression (Phase 6 + R-301).
    p.add_argument("--paraview-volume-vtu", action="store_true")
    p.add_argument("--paraview-volume-hdf5", action="store_true")
    p.add_argument("--paraview-volume-zfp-tol", type=float, default=None)
    p.add_argument("--paraview-volume-deflate-level", type=int, default=None)

    # Bulk-side compression (Phase 6.4 — secondary collection).
    p.add_argument("--paraview-bulk-zfp-tol", type=float, default=None)
    p.add_argument("--paraview-bulk-deflate-level", type=int, default=None)

    # Fault-side compression (Phase 2d.3).
    p.add_argument("--paraview-fault-zfp-tol", type=float, default=None)
    p.add_argument("--paraview-fault-deflate-level", type=int, default=None)
    p.add_argument("--paraview-fault-vtu", action="store_true")
    p.add_argument("--paraview-fault-hdf5", action="store_true")
    p.add_argument("--paraview-fault-legacy-ascii", action="store_true")

    # Estimator-only flags.
    p.add_argument("--np", type=int, default=1,
                   help="Rank count (informational; per-rank bytes "
                        "computed inside the writer).")
    p.add_argument("--scratch-quota", type=str, default=None,
                   help="Available scratch quota, e.g. '1TB' or '500GB'.")
    p.add_argument("--n-events", type=int, default=0,
                   help="Override estimated coseismic event count.")
    p.add_argument("--avg-event-duration", type=float, default=30.0,
                   help="Average duration of one coseismic event (s).")
    p.add_argument("--station-count", type=int, default=None,
                   help="Probe-station count (default 7 for BP5, 0 for "
                        "TPV*).")
    p.add_argument("--station-columns", type=int, default=8)
    p.add_argument("--explain", action="store_true",
                   help="Print the formula chain that produced the "
                        "estimate.")
    p.add_argument("--json", action="store_true",
                   help="Print machine-readable JSON only.")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress PLACEHOLDER-source warnings (R-309).")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.quiet:
        import warnings
        warnings.simplefilter("ignore")

    # Resolve mesh.
    if args.mesh is not None:
        mesh = read_gmsh(args.mesh, fault_tag=args.fault_tag)
    elif args.inline_mesh:
        mesh = read_inline_mesh(
            args.inline_nx, args.inline_ny, args.inline_nz,
            element_type=args.inline_element_type)
    else:
        parser.error("either --mesh or --inline-mesh must be provided")
        return 2

    # Resolve schedule.
    schedule = ScheduleConfig(
        tfinal=parse_time(args.tfinal),
        v_coseismic=args.paraview_v_coseismic,
        v_nucleation=args.paraview_v_nucleation,
        dt_coseismic=args.paraview_dt_coseismic,
        dt_nucleation=args.paraview_dt_nucleation,
        dt_interseismic=args.paraview_dt_interseismic,
        max_total_snapshots=args.paraview_max_snapshots,
        output_every_n_steps=args.paraview_output_every_n_steps,
        fixed_dt=args.paraview_fixed_dt,
        n_events=args.n_events,
        avg_event_duration_s=args.avg_event_duration,
        driver=args.driver,
    )

    # Resolve compression filters.
    if args.paraview_volume_vtu and args.paraview_volume_hdf5:
        parser.error(
            "--paraview-volume-vtu and --paraview-volume-hdf5 are "
            "mutually exclusive")
    volume_mode = "vtu" if args.paraview_volume_vtu else "hdf5"
    primary_volume_filter = _filter_name_from_args(
        args.paraview_volume_zfp_tol, args.paraview_volume_deflate_level,
        volume_mode)
    fault_mode = "vtu" if args.paraview_fault_vtu else "hdf5"
    fault_filter = _filter_name_from_args(
        args.paraview_fault_zfp_tol, args.paraview_fault_deflate_level,
        fault_mode)

    filters = {"fault": fault_filter, "volume": primary_volume_filter}

    # Quota.
    scratch_quota_bytes = (parse_bytes(args.scratch_quota)
                           if args.scratch_quota else None)

    if args.station_count is not None:
        station_count = args.station_count
    elif args.driver == "bp5":
        station_count = 7
    else:
        station_count = 0  # TPV* drivers don't emit station CSVs.

    report = estimate(
        driver=args.driver,
        mesh=mesh,
        order=args.order,
        schedule=schedule,
        filters=filters,
        scratch_quota_bytes=scratch_quota_bytes,
        paraview_enabled=args.paraview,
        volume_pv_enabled=not args.no_volume_pv,
        volume_pv_dt=args.volume_pv_dt,
        explain=args.explain,
        station_count=station_count,
        station_columns=args.station_columns,
    )

    if args.json:
        # Drop explain_lines from the JSON unless explicitly requested
        # (it's a long list of formatted strings).
        d = asdict(report)
        if not args.explain:
            d.pop("explain_lines", None)
        json.dump(d, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        _print_human(
            report, args.driver, mesh, schedule, filters,
            paraview_enabled=args.paraview,
            volume_pv_enabled=not args.no_volume_pv,
            scratch_quota_bytes=scratch_quota_bytes)

    # Exit code: 1 if quota exceeded, 0 otherwise.
    if report.quota_fraction is not None and report.quota_fraction > 1.0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
