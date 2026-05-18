# Phase 7.4 of miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md.
#
# Compression-ratio table.  Maps `(field_kind, filter)` to a
# `(bytes_per_dof, source)` pair.  R-309: every entry whose source
# is the literal string `"PLACEHOLDER"` is an initial guess pending
# §Phase 7.6 calibration; the runtime emits a one-time warning per
# (kind, filter) pair when consulted.

from __future__ import annotations

import re
import sys
import warnings
from typing import Dict, Tuple

from _io_size_schemas import FieldSchema


# (field_kind, filter) -> (bytes_per_dof, source_or_PLACEHOLDER).
COMPRESSION_RATIOS: Dict[Tuple[str, str], Tuple[float, str]] = {
    # Uncompressed baselines.
    ("any", "none"):                    (8.0, "raw IEEE-754 double"),
    ("any", "vtu_ascii"):               (25.0, "legacy ASCII baseline"),
    ("any", "vtu_binary"):              (8.0, "raw double"),

    # Smooth floating-point fields (velocity, displacement, stress).
    ("smooth_fp", "deflate_6"):         (3.2, "PLACEHOLDER"),
    ("smooth_fp", "zfp_1e-3"):          (0.4, "PLACEHOLDER"),
    ("smooth_fp", "zfp_1e-6"):          (1.6, "PLACEHOLDER"),
    ("smooth_fp", "zfp_1e-9"):          (3.2, "PLACEHOLDER"),
    ("smooth_fp", "zfp_1e-12"):         (4.8, "PLACEHOLDER"),

    # Rough floating-point fields (slip-rate during nucleation,
    # state variable).  Less compressible than smooth fields.
    ("rough_fp", "deflate_6"):          (5.5, "PLACEHOLDER"),
    ("rough_fp", "zfp_1e-3"):           (1.0, "PLACEHOLDER"),
    ("rough_fp", "zfp_1e-12"):          (5.0, "PLACEHOLDER"),

    # Integer index / rank field — extremely compressible.
    ("integer_index", "deflate_6"):     (0.5, "PLACEHOLDER"),
    # R-406: ZFP only filters floating-point datasets in VTKHDF;
    # integer connectivity / index datasets fall back to the lossless
    # deflate path via `VTKHDF::EnsureDataset` (mesh/vtkhdf.cpp:162-170).
    # The on-disk size is therefore the same as the deflate column.
    ("integer_index", "zfp_1e-3"):      (0.5, "VTKHDF deflate fallback"),
    ("integer_index", "zfp_1e-6"):      (0.5, "VTKHDF deflate fallback"),
    ("integer_index", "zfp_1e-9"):      (0.5, "VTKHDF deflate fallback"),
    ("integer_index", "zfp_1e-12"):     (0.5, "VTKHDF deflate fallback"),
    ("integer_index", "vtu_binary"):    (8.0, "raw double"),

    # Static 2D fields (param_a, param_Dc, fault_x2/x3): emitted only
    # at first save; compress as smooth_fp but with a 1× multiplier
    # since they're written once per run.
    ("static_2d", "deflate_6"):         (3.0, "PLACEHOLDER"),
    ("static_2d", "zfp_1e-12"):         (4.0, "PLACEHOLDER"),
}


# Re-export the per-driver nucleation table here per plan §Phase 7.3
# Detailed Requirements step 1: "where `NUCLEATION_DURATION_S` is a
# per-driver constant in `_io_size_compression.py` (R-307)."
NUCLEATION_DURATION_S = {
    "bp5":    86_400.0,
    "tpv102": 1.0,
    "tpv104": 1.0,
    "tpv205": 1.0,
}


# Per-(kind, filter) one-time warning latch.
_PLACEHOLDER_WARNED: set = set()


def classify_field(field_name: str, driver: str) -> str:
    """Return the compression `kind` for `field_name`:
    `smooth_fp` / `rough_fp` / `integer_index` / `static_2d` /
    `unknown`.

    Smooth-FP: velocity, displacement, traction (smooth in space and
      slowly varying in time).
    Rough-FP: slip-rate, state-variable, normal-stress (high spatial
      contrast during nucleation; less compressible).
    Integer-index: mpi_rank.
    Static-2D: param_a / param_Dc / fault_x2 / fault_x3 (emitted once
      at first save).
    """
    if field_name in {"mpi_rank"}:
        return "integer_index"
    if field_name in {"param_a", "param_Dc", "fault_x2", "fault_x3"}:
        return "static_2d"
    if field_name in {"slip_rate_dip", "slip_rate_strike",
                      "slip_rate_dip_k4", "slip_rate_strike_k4",
                      "state_variable",
                      "normal_stress", "normal_stress_k4"}:
        return "rough_fp"
    if field_name in {"velocity", "displacement",
                      "slip_dip", "slip_strike",
                      "traction_dip", "traction_strike",
                      "traction_dip_k4", "traction_strike_k4"}:
        return "smooth_fp"
    # Stress component fields on TPV* bulk side: sigma_xx / sigma_xy / ...
    if field_name.startswith("sigma_"):
        return "smooth_fp"
    return "unknown"


_FILTER_RE = re.compile(
    r"^("
    r"none|vtu_ascii|vtu_binary"
    r"|deflate_[0-9]"
    r"|zfp_1e-(?:3|6|9|12)"
    r")$")


def _validate_filter(filter_name: str) -> None:
    if not _FILTER_RE.match(filter_name):
        raise ValueError(
            f"unrecognised compression filter {filter_name!r}; "
            f"supported: none, vtu_ascii, vtu_binary, "
            f"deflate_<N> for N=0..9, zfp_<tol> for "
            f"tol in {{1e-3, 1e-6, 1e-9, 1e-12}}")


def bytes_per_dof(field: FieldSchema,
                  filter_name: str,
                  driver: str) -> float:
    """Return bytes-per-DOF for `field` under `filter_name`.

    Lookup order:
      1. (classify_field(field.name), filter_name)
      2. ("any", filter_name)
      3. fallback: 8.0 (raw double) with stderr warning

    Emits a one-time warning per (kind, filter) pair when the source
    is `"PLACEHOLDER"` (R-309).
    """
    _validate_filter(filter_name)
    kind = classify_field(field.name, driver)
    key_specific = (kind, filter_name)
    key_any = ("any", filter_name)

    entry = COMPRESSION_RATIOS.get(key_specific)
    if entry is None:
        entry = COMPRESSION_RATIOS.get(key_any)
    if entry is None:
        # Out-of-table fallback (plan §Phase 7.4 Edge Cases).
        print(
            f"warning: no compression-ratio entry for "
            f"({kind}, {filter_name}); falling back to 8.0 bytes/DOF",
            file=sys.stderr)
        return 8.0

    bpd, source = entry
    if source == "PLACEHOLDER" and key_specific not in _PLACEHOLDER_WARNED:
        _PLACEHOLDER_WARNED.add(key_specific)
        warnings.warn(
            f"compression ratio for (field_kind={kind}, "
            f"filter={filter_name}) is a PLACEHOLDER pending "
            f"Phase 7.6 calibration; the size estimate may be off "
            f"by >2x until calibration completes.",
            UserWarning,
            stacklevel=2)
    return bpd
