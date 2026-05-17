"""csm_johnson_hearn_reader.py — CSV reader for the CSM Johnson & Hearn
stressing-rate slices.

The CSM raw bundle (``stress/raw/csm_johnson_hearn/``) ships one CSV per
depth (1, 3, 5, 7, 9, 11, 13, 15 km).  Each file is a regular (lon, lat)
grid covering the same geographic footprint and carries 35 columns of
which we keep the six independent components of the stressing-rate
tensor (See, Sen, Seu, Snn, Snu, Suu — in MPa/yr, tension POSITIVE).

The CSV idiosyncrasies handled here:

- The header is a block of ``#``-prefixed comment lines (~28 of them)
  followed by a single column-name header line (``LON,LAT,DEP,See,…``)
  that is NOT comment-prefixed.
- The longitude column is double-quoted (``"-119.149"``); ``np.loadtxt``
  cannot strip those, so we line-clean before parsing.
- ``SHmax_unc`` (column 11) is allowed to be blank (the ``,,`` token).
  We map blank fields to NaN and only retain the six tensor columns,
  which are always non-blank in the SAFS bundle.

Output mirrors ``raw_readers.RawSliceStack``: ``(lon_grid, lat_grid,
depths_m, fields)`` with ``fields`` keyed by the canonical schema-v1
names ``sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz``
under the East→x, North→y, Up→z mapping.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable

import numpy as np


# Sibling import of velocity/code/raw_readers.RawSliceStack so that the
# stress driver carries the same dataclass shape as the velocity driver
# (the consumer side — _resample_field, _build_utm_grid — depends on
# that exact attribute set).
_VELOCITY_CODE_DIR = (
    Path(__file__).resolve().parents[2] / "velocity" / "code"
)
if str(_VELOCITY_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_VELOCITY_CODE_DIR))
from raw_readers import RawSliceStack  # noqa: E402


# East,North,Up → x,y,z mapping for the symmetric 2nd-order tensor.
# Keys are the canonical schema-v1 names; values are the CSV column
# indices (0-based) in the CSM Johnson-Hearn raw layout:
#   col 0  = LON              col 6  = Snn
#   col 1  = LAT              col 7  = Snu
#   col 2  = DEP              col 8  = Suu
#   col 3  = See              ...
#   col 4  = Sen
#   col 5  = Seu
_TENSOR_COLUMNS: dict[str, int] = {
    "sigma_xx": 3,   # See  (E-E)
    "sigma_xy": 4,   # Sen  (E-N)
    "sigma_xz": 5,   # Seu  (E-U)
    "sigma_yy": 6,   # Snn  (N-N)
    "sigma_yz": 7,   # Snu  (N-U)
    "sigma_zz": 8,   # Suu  (U-U)
}

# Soft tolerance for round-tripping the lon/lat values that come back as
# slightly noisy floats out of the CSV — 1e-6 degrees ≈ 11 cm at this
# latitude, far below the 0.018° native sample spacing.
_LL_ROUND_DECIMALS: int = 6


def _read_one_csv(path: Path) -> tuple[float, np.ndarray]:
    """Read a single CSM Johnson-Hearn CSV slice.

    Returns ``(depth_m, body)`` where ``depth_m`` is the constant DEP
    value in metres (positive=down, per the CSV "Depth in km relative
    to mean sea level (positive: below sea-level)") and ``body`` is a
    ``(n_rows, 35)`` float64 array with blank cells filled as NaN and
    longitude quotes stripped.
    """
    rows: list[list[float]] = []
    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            # The first non-comment line is the column header
            # `LON,LAT,DEP,See,…`; skip it.
            if line.lower().startswith("lon"):
                continue
            # Strip embedded double quotes ("-119.149" → -119.149) and
            # split.  Blank fields (",,") become NaN.
            cleaned = line.replace('"', "")
            parts = cleaned.split(",")
            row = [float(p) if p.strip() else np.nan for p in parts]
            rows.append(row)
    if not rows:
        raise ValueError(
            f"csm_johnson_hearn_reader: {path} contains no data rows"
        )
    body = np.asarray(rows, dtype=np.float64)
    if body.shape[1] < 9:
        # We need at least cols 0..8 (LON, LAT, DEP, See, Sen, Seu,
        # Snn, Snu, Suu).
        raise ValueError(
            f"csm_johnson_hearn_reader: {path} has only {body.shape[1]} "
            f"columns; expected at least 9 (LON,LAT,DEP,See,Sen,Seu,"
            f"Snn,Snu,Suu)."
        )
    depth_km_vals = np.unique(body[:, 2])
    if depth_km_vals.size != 1:
        raise ValueError(
            f"csm_johnson_hearn_reader: {path} has {depth_km_vals.size} "
            f"distinct DEP values {depth_km_vals.tolist()}; each CSV "
            f"must be a single depth slice."
        )
    depth_m = float(depth_km_vals[0]) * 1000.0
    return depth_m, body


def _grid_from_body(path: Path, body: np.ndarray
                    ) -> tuple[np.ndarray, np.ndarray,
                               np.ndarray, np.ndarray]:
    """Build the (lon, lat) axes plus per-row (lon, lat) index arrays
    for reshaping the body into a regular (N_lon, N_lat) grid.

    Returns ``(lon_axis, lat_axis, lon_idx, lat_idx)``.
    Raises ``ValueError`` if the per-slice point cloud is not the
    expected regular grid (N_lon × N_lat == n_rows).
    """
    lon_vals = body[:, 0]
    lat_vals = body[:, 1]
    lon_axis = np.unique(np.round(lon_vals, _LL_ROUND_DECIMALS))
    lat_axis = np.unique(np.round(lat_vals, _LL_ROUND_DECIMALS))
    n_lon = lon_axis.size
    n_lat = lat_axis.size
    if n_lon * n_lat != body.shape[0]:
        raise ValueError(
            f"csm_johnson_hearn_reader: {path} is not a regular grid "
            f"(n_lon*n_lat = {n_lon * n_lat} != rows = {body.shape[0]}). "
            f"Distinct lon values: {n_lon}, distinct lat values: {n_lat}."
        )
    lon_idx = np.searchsorted(
        lon_axis, np.round(lon_vals, _LL_ROUND_DECIMALS)
    )
    lat_idx = np.searchsorted(
        lat_axis, np.round(lat_vals, _LL_ROUND_DECIMALS)
    )
    if (lon_idx >= n_lon).any() or (lat_idx >= n_lat).any():
        raise ValueError(
            f"csm_johnson_hearn_reader: {path} indexing failure; "
            f"likely irregular lon/lat sampling."
        )
    return lon_axis, lat_axis, lon_idx, lat_idx


def read_csm_johnson_hearn(paths: Iterable[Path]) -> RawSliceStack:
    """Parse a list of CSM Johnson & Hearn CSV slices into a single
    ``RawSliceStack``.

    Each file is one horizontal slice at a single depth.  Depth is
    taken from the per-row ``DEP`` column (km, positive=down) and
    converted to metres.  All slices must share an identical
    ``(lon_grid, lat_grid)`` after rounding to ``_LL_ROUND_DECIMALS``
    decimal places.

    Parameters
    ----------
    paths : iterable of Path
        Paths to CSM CSV slices.  Order does not matter; the returned
        stack is sorted by ascending depth.

    Returns
    -------
    RawSliceStack
        ``fields`` is keyed by the canonical schema-v1 names
        ``sigma_xx, sigma_xy, sigma_xz, sigma_yy, sigma_yz, sigma_zz``
        (insertion order preserved for downstream determinism).
        Values are in MPa/yr, tension POSITIVE (source convention).

    Raises
    ------
    ValueError
        On malformed files, mismatched grids, duplicate depths, or any
        violation of regular-grid reshape invariants.
    FileNotFoundError
        If any path does not exist.
    """
    paths = [Path(p) for p in paths]
    if not paths:
        raise ValueError("read_csm_johnson_hearn: empty path list")

    # list[(depth_m, lon_axis, lat_axis, fields_for_slice)]
    slices: list[tuple[float, np.ndarray, np.ndarray, dict]] = []
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"CSM CSV not found: {p}")
        depth_m, body = _read_one_csv(p)
        lon_axis, lat_axis, lon_idx, lat_idx = _grid_from_body(p, body)
        n_lon = lon_axis.size
        n_lat = lat_axis.size

        fields_for_slice: dict[str, np.ndarray] = {}
        for canonical_name, col in _TENSOR_COLUMNS.items():
            arr = np.full((n_lon, n_lat), np.nan, dtype=np.float64)
            arr[lon_idx, lat_idx] = body[:, col]
            if np.isnan(arr).any():
                raise ValueError(
                    f"csm_johnson_hearn_reader: {p} reshape produced "
                    f"NaN holes for field '{canonical_name}' (col "
                    f"{col}); irregular grid in source CSV."
                )
            fields_for_slice[canonical_name] = arr
        slices.append((depth_m, lon_axis, lat_axis, fields_for_slice))

    # All slices must share the same lon/lat axes.
    ref_lon = slices[0][1]
    ref_lat = slices[0][2]
    for depth_m, lon_g, lat_g, _ in slices[1:]:
        if not np.array_equal(lon_g, ref_lon):
            raise ValueError(
                f"read_csm_johnson_hearn: slice at depth {depth_m} m "
                f"has lon_grid different from reference; all CSM "
                f"slices must share one geographic grid."
            )
        if not np.array_equal(lat_g, ref_lat):
            raise ValueError(
                f"read_csm_johnson_hearn: slice at depth {depth_m} m "
                f"has lat_grid different from reference."
            )

    # Sort by ascending depth; reject duplicates.
    slices.sort(key=lambda s: s[0])
    depths = np.asarray([s[0] for s in slices], dtype=np.float64)
    if depths.size > 1 and np.any(np.diff(depths) <= 0.0):
        dup = float(depths[int(np.argmin(np.diff(depths)))])
        raise ValueError(
            f"read_csm_johnson_hearn: duplicate or non-monotone depths "
            f"in slice set near depth = {dup} m."
        )

    n_lon = ref_lon.size
    n_lat = ref_lat.size
    n_z = depths.size
    field_names = list(slices[0][3].keys())
    fields: dict[str, np.ndarray] = {
        name: np.empty((n_lon, n_lat, n_z), dtype=np.float64)
        for name in field_names
    }
    for k, (_, _, _, sl_fields) in enumerate(slices):
        for name in field_names:
            fields[name][:, :, k] = sl_fields[name]

    return RawSliceStack(
        lon_grid=ref_lon,
        lat_grid=ref_lat,
        depths_m=depths,
        fields=fields,
    )
