"""sidecar.py — schema-v1 HDF5 writer / attribute reader.

See ../../document/features_dev/data_projection_schema_v1.md for the full
schema.  Both the writer and the runtime C++ reader independently enforce
the same value constraints (defense-in-depth).
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Any

import h5py
import numpy as np


SCHEMA_VERSION = "data_projection_v1"
CANONICAL_CRS = "EPSG:32611"
CANONICAL_UNITS = "m"
CANONICAL_Z_POSITIVE = "elevation"


def _check_axis_monotone(axis: np.ndarray, name: str) -> None:
    if axis.ndim != 1:
        raise ValueError(f"axis '{name}' must be 1-D; got shape {axis.shape}")
    if axis.size < 1:
        raise ValueError(f"axis '{name}' is empty")
    if axis.size >= 2 and not np.all(np.diff(axis) > 0.0):
        raise ValueError(
            f"axis '{name}' is not strictly monotone increasing; "
            f"first non-increasing pair near index "
            f"{int(np.argmin(np.diff(axis)))}")


def _now_iso() -> str:
    return (_dt.datetime.now(_dt.timezone.utc)
            .strftime("%Y-%m-%dT%H:%M:%SZ"))


def write_sidecar(out_path: Path,
                  x: np.ndarray, y: np.ndarray, z: np.ndarray,
                  fields: dict,
                  attrs: dict,
                  field_bounds: dict) -> None:
    """Write a schema-v1 HDF5 sidecar.

    Parameters
    ----------
    out_path : Path
        Destination ``.h5`` file. Parent directory created if missing.
    x, y, z : 1-D float arrays
        Strictly monotone increasing axes in canonical CRS (UTM 11 N, m).
    fields : dict[str -> (Nx, Ny, Nz) float array]
        Each field array's shape must match ``(len(x), len(y), len(z))``
        and contain NO NaN.
    attrs : dict
        Root-level attribute overrides. Required: ``schema_version``,
        ``crs``, ``units``, ``z_positive`` (this writer fills these with
        the canonical values automatically; passing a non-canonical
        value raises). Optional: ``created_at`` (filled with current
        UTC time if absent), ``source``, ``source_crs``, ``mesh_tag``.
    field_bounds : dict[str -> (min_value, max_value)]
        Per-field sanity bounds. Every key in ``fields`` must be
        present. Every cell must satisfy
        ``min_value <= cell <= max_value`` or this writer raises.
        Each entry may also include a third element ``units``;
        otherwise the sidecar's per-field ``units`` attribute is left
        unset (raises since it is required).

    Raises
    ------
    ValueError
        On any schema violation (NaN, out-of-bounds value, axis-shape
        mismatch, missing field bounds).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    z = np.ascontiguousarray(z, dtype=np.float64)
    _check_axis_monotone(x, "x")
    _check_axis_monotone(y, "y")
    _check_axis_monotone(z, "z")

    if not isinstance(fields, dict) or not fields:
        raise ValueError("write_sidecar: 'fields' must be a non-empty dict")

    expected_shape = (x.size, y.size, z.size)
    for name, arr in fields.items():
        if not isinstance(arr, np.ndarray):
            raise ValueError(
                f"write_sidecar: field '{name}' must be a numpy array")
        if arr.shape != expected_shape:
            raise ValueError(
                f"write_sidecar: field '{name}' shape {arr.shape} "
                f"!= expected {expected_shape}")
        if not np.issubdtype(arr.dtype, np.floating):
            raise ValueError(
                f"write_sidecar: field '{name}' dtype {arr.dtype} "
                f"is not floating-point")
        if np.isnan(arr).any():
            n_nan = int(np.isnan(arr).sum())
            raise ValueError(
                f"write_sidecar: field '{name}' contains {n_nan} NaN "
                f"cell(s); v1 schema forbids NaN. Either restrict the "
                f"output grid to the source data hull or regenerate "
                f"the source dataset.")
        if name not in field_bounds:
            raise ValueError(
                f"write_sidecar: field '{name}' is missing in "
                f"field_bounds (required by v1 schema)")
        bounds = field_bounds[name]
        if len(bounds) < 2:
            raise ValueError(
                f"write_sidecar: field_bounds['{name}'] must be "
                f"(min, max) [, units]; got {bounds}")
        v_min = float(bounds[0])
        v_max = float(bounds[1])
        if not (v_min < v_max):
            raise ValueError(
                f"write_sidecar: field_bounds['{name}'] requires "
                f"min < max; got [{v_min}, {v_max}]")
        below = arr < v_min
        above = arr > v_max
        n_below = int(below.sum())
        n_above = int(above.sum())
        if n_below or n_above:
            arr_min = float(arr.min())
            arr_max = float(arr.max())
            raise ValueError(
                f"write_sidecar: field '{name}' has {n_below} cell(s) "
                f"below min={v_min} and {n_above} cell(s) above "
                f"max={v_max}; observed range [{arr_min}, {arr_max}]. "
                f"Either widen the bounds in field_bounds or "
                f"regenerate the source dataset.")

    # Validate / fill root attrs.
    out_attrs = {
        "schema_version": SCHEMA_VERSION,
        "crs": CANONICAL_CRS,
        "units": CANONICAL_UNITS,
        "z_positive": CANONICAL_Z_POSITIVE,
        "created_at": _now_iso(),
    }
    for key, val in (attrs or {}).items():
        if key in ("schema_version", "crs", "units", "z_positive"):
            # Reject non-canonical overrides.
            if str(val) != out_attrs[key]:
                raise ValueError(
                    f"write_sidecar: attr '{key}' must equal "
                    f"'{out_attrs[key]}' in v1 schema; got '{val}'")
        out_attrs[key] = val

    with h5py.File(out_path, "w") as h5:
        for k, v in out_attrs.items():
            h5.attrs[k] = v
        grid = h5.create_group("grid")
        grid.create_dataset("x", data=x)
        grid.create_dataset("y", data=y)
        grid.create_dataset("z", data=z)
        flds = h5.create_group("fields")
        for name, arr in fields.items():
            ds = flds.create_dataset(name, data=arr)
            v_min, v_max = field_bounds[name][0], field_bounds[name][1]
            ds.attrs["min_value"] = float(v_min)
            ds.attrs["max_value"] = float(v_max)
            if len(field_bounds[name]) >= 3:
                ds.attrs["units"] = str(field_bounds[name][2])
            else:
                raise ValueError(
                    f"write_sidecar: field_bounds['{name}'] must include "
                    f"a third element 'units' (required by v1 schema)")


def read_sidecar_attrs(path: Path) -> dict:
    """Return root-level attributes of an existing sidecar as a dict.
    Useful for tests and for tooling that wants to inspect a sidecar
    without loading the field arrays.
    """
    with h5py.File(Path(path), "r") as h5:
        return {k: (v.decode() if isinstance(v, bytes) else v)
                for k, v in h5.attrs.items()}
