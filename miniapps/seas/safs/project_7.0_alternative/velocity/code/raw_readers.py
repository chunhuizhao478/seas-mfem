"""raw_readers.py — readers for raw scalar/vector field datasets.

Each reader is generic about field count: a single reader can produce a
RawSliceStack carrying any number of named fields (Vp/Vs/density for
CVM-H, sigma_xx/.../sigma_xz for stress maps, p for pore pressure).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass
class RawSliceStack:
    """A stack of horizontal slices on a regular geographic grid.

    Attributes
    ----------
    lon_grid : (N_lon,) float64, strictly increasing
    lat_grid : (N_lat,) float64, strictly increasing
    depths_m : (N_z,)   float64, monotone increasing (positive=down)
    fields   : dict[name -> (N_lon, N_lat, N_z) float64]
    """
    lon_grid: np.ndarray
    lat_grid: np.ndarray
    depths_m: np.ndarray
    fields: dict


@dataclass
class RawGrid:
    """A single 2-D grid on a regular geographic (lon, lat) net,
    optionally with multiple named fields. Used by stress/pore-pressure
    readers that have no depth axis.

    Attributes
    ----------
    lon_grid : (N_lon,) float64, strictly increasing
    lat_grid : (N_lat,) float64, strictly increasing
    fields   : dict[name -> (N_lon, N_lat) float64]
    """
    lon_grid: np.ndarray
    lat_grid: np.ndarray
    fields: dict


# ---------------------------------------------------------------------------
# CVM-H ASCII reader
# ---------------------------------------------------------------------------

_CVMH_REQUIRED_HEADER_KEYS = (
    "Depth(m)", "Lon_pts", "Lat_pts", "Total_pts",
    "Lat1", "Lon1",
    # Spacing is required for the corner-anchor computation; Lat2/Lon2
    # are NOT required because the real CVM-H archive ships values that
    # differ from the body by up to ~3 sample widths (header bug); we
    # ignore them and derive the SE corner from (Lon1, Lat1, Spacing,
    # Lon_pts, Lat_pts) directly.
    "Spacing(degree)",
)


def _parse_cvmh_header(path: Path) -> dict:
    """Parse the `# key: value` header lines of a single CVM-H slice.

    Returns a dict of header-key -> string value plus an entry
    ``"_column_line"`` holding the column-spec line that follows the
    header. Raises ``ValueError`` if the column line is missing or any
    required key is absent.
    """
    header = {}
    column_line = None
    with open(path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            if not line.startswith("#"):
                # End of header (no more comment lines); the column line
                # must have been the last `#` line. If we never saw it,
                # the file is malformed.
                break
            stripped = line[1:].strip()
            # Detect the column-spec line: it starts with "Lon" (case
            # insensitive) and contains commas, e.g.
            # `# Lon,Lat,Vp(m/s),Vs(m/s),Density(kg/m^3)`.
            if stripped.lower().startswith("lon") and "," in stripped:
                column_line = stripped
                continue
            # Otherwise it's a `key: value` pair.
            if ":" not in stripped:
                continue
            key, _, value = stripped.partition(":")
            header[key.strip()] = value.strip()
    if column_line is None:
        raise ValueError(
            f"{path}: missing column-spec line (e.g. "
            f"`# Lon,Lat,Vp(m/s),Vs(m/s),Density(kg/m^3)`) in header.")
    header["_column_line"] = column_line
    missing = [k for k in _CVMH_REQUIRED_HEADER_KEYS if k not in header]
    if missing:
        raise ValueError(
            f"{path}: missing required header keys: {missing}. "
            f"Got keys: {sorted(k for k in header if not k.startswith('_'))}")
    return header


def _column_names_from_line(column_line: str) -> list:
    """Strip unit suffixes from the column-spec line, returning lowercase
    canonical names. E.g. `Lon,Lat,Vp(m/s),Vs(m/s),Density(kg/m^3)` ->
    ``["lon", "lat", "vp", "vs", "density"]``.
    """
    out = []
    for tok in column_line.split(","):
        tok = tok.strip()
        # Drop anything in parentheses (units).
        paren = tok.find("(")
        if paren >= 0:
            tok = tok[:paren].strip()
        out.append(tok.lower())
    return out


def read_cvmh_ascii(paths: Iterable[Path]) -> RawSliceStack:
    """Parse a list of CVM-H ASCII slice files into a single
    ``RawSliceStack``.

    Each file is one horizontal slice at a single depth. Depth is taken
    from the ``# Depth(m): N`` header line (NOT the filename — the
    existing CVM-H archive uses misleading ``X.Ykm`` filename labels
    equal to ``depth_km / 10``).

    All slices must share an identical ``(lon_grid, lat_grid)``.

    Parameters
    ----------
    paths : iterable of Path
        Paths to CVM-H slice files. Order does not matter; the result
        is sorted by depth ascending.

    Returns
    -------
    RawSliceStack

    Raises
    ------
    ValueError
        On malformed files, mismatched grids, duplicate depths, or any
        violation of the per-row corner pinning.
    """
    paths = [Path(p) for p in paths]
    if not paths:
        raise ValueError("read_cvmh_ascii: empty path list")

    slices = []  # list of (depth_m, lon_grid, lat_grid, fields_dict)
    field_names = None
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"CVM-H slice not found: {p}")
        header = _parse_cvmh_header(p)
        depth_m = float(header["Depth(m)"])
        n_lon = int(header["Lon_pts"])
        n_lat = int(header["Lat_pts"])
        n_total = int(header["Total_pts"])
        if n_lon * n_lat != n_total:
            raise ValueError(
                f"{p}: Lon_pts*Lat_pts ({n_lon}*{n_lat}={n_lon*n_lat}) "
                f"!= Total_pts ({n_total})")
        col_names = _column_names_from_line(header["_column_line"])
        if "lon" not in col_names or "lat" not in col_names:
            raise ValueError(
                f"{p}: column line lacks 'lon'/'lat': {header['_column_line']}")
        # All non-(lon, lat) columns are fields.
        fld_cols = [c for c in col_names if c not in ("lon", "lat")]
        if not fld_cols:
            raise ValueError(
                f"{p}: no field columns in {header['_column_line']}")
        if field_names is None:
            field_names = fld_cols
        elif fld_cols != field_names:
            raise ValueError(
                f"{p}: field column mismatch: got {fld_cols}, "
                f"expected {field_names}")

        # Read the body. numpy.loadtxt skips '#'-prefixed lines.
        body = np.loadtxt(str(p), delimiter=",", comments="#")
        if body.ndim != 2 or body.shape != (n_total, len(col_names)):
            raise ValueError(
                f"{p}: body shape {body.shape} does not match "
                f"({n_total}, {len(col_names)})")

        # Pin orientation against the SW corner from (Lon1, Lat1) AND
        # the NE corner derived as Lon1 + (N_lon-1)*Spacing,
        # Lat1 + (N_lat-1)*Spacing.  The CVM-H archive ships
        # `Lon2`/`Lat2` values that disagree with the body by up to
        # ~3 sample widths (header bug); we therefore derive the
        # expected NE corner from the spacing instead of trusting
        # the header values, and check the header values only as a
        # soft consistency hint.
        lon1 = float(header["Lon1"])
        lat1 = float(header["Lat1"])
        spacing = float(header["Spacing(degree)"])
        lon_ne_expected = lon1 + (n_lon - 1) * spacing
        lat_ne_expected = lat1 + (n_lat - 1) * spacing
        eps_corner = 1e-3   # 0.001 deg ~ 100 m at this latitude
        # First sample row should be the SW corner (lon1, lat1).
        if (abs(body[0, 0] - lon1) > eps_corner or
                abs(body[0, 1] - lat1) > eps_corner):
            raise ValueError(
                f"{p}: first row {tuple(body[0, :2])} does not match "
                f"header (Lon1, Lat1)=({lon1}, {lat1})")
        # Last sample row should be the NE corner derived from spacing.
        if (abs(body[-1, 0] - lon_ne_expected) > eps_corner or
                abs(body[-1, 1] - lat_ne_expected) > eps_corner):
            raise ValueError(
                f"{p}: last row {tuple(body[-1, :2])} does not match "
                f"derived NE corner ({lon_ne_expected}, "
                f"{lat_ne_expected}) from Lon1+(N_lon-1)*Spacing, "
                f"Lat1+(N_lat-1)*Spacing")

        # Build the lon/lat axes from unique values, sorted ascending.
        lon_vals_full = body[:, 0].copy()
        lat_vals_full = body[:, 1].copy()
        lon_unique = np.unique(np.round(lon_vals_full, 6))
        lat_unique = np.unique(np.round(lat_vals_full, 6))
        if lon_unique.size != n_lon:
            raise ValueError(
                f"{p}: distinct lon values ({lon_unique.size}) != "
                f"Lon_pts ({n_lon}); duplicates or noise in input")
        if lat_unique.size != n_lat:
            raise ValueError(
                f"{p}: distinct lat values ({lat_unique.size}) != "
                f"Lat_pts ({n_lat})")

        # Reshape body to (N_lon, N_lat) per field, using digitize against
        # the sorted unique axes (more robust than assuming a particular
        # order in the file).
        lon_idx = np.searchsorted(lon_unique, np.round(lon_vals_full, 6))
        lat_idx = np.searchsorted(lat_unique, np.round(lat_vals_full, 6))
        if (lon_idx >= n_lon).any() or (lat_idx >= n_lat).any():
            raise ValueError(
                f"{p}: indexing failure in reshape; possibly "
                f"non-uniform lon/lat sampling")

        fields_for_slice = {}
        for ci, name in enumerate(col_names):
            if name in ("lon", "lat"):
                continue
            arr = np.full((n_lon, n_lat), np.nan, dtype=np.float64)
            arr[lon_idx, lat_idx] = body[:, ci]
            if np.isnan(arr).any():
                raise ValueError(
                    f"{p}: reshape produced NaN holes for field "
                    f"'{name}' — irregular grid?")
            fields_for_slice[name] = arr
        slices.append((depth_m, lon_unique, lat_unique, fields_for_slice))

    # Verify all slices share an identical lon/lat grid.
    ref_lon = slices[0][1]
    ref_lat = slices[0][2]
    for depth_m, lon_g, lat_g, _ in slices[1:]:
        if not np.array_equal(lon_g, ref_lon):
            raise ValueError(
                f"slice at depth {depth_m} m has lon_grid different "
                f"from reference; the CVM-H archive must use one grid "
                f"across slices")
        if not np.array_equal(lat_g, ref_lat):
            raise ValueError(
                f"slice at depth {depth_m} m has lat_grid different "
                f"from reference")

    # Sort slices by depth ascending and check for duplicate depths.
    slices.sort(key=lambda s: s[0])
    depths = np.asarray([s[0] for s in slices], dtype=np.float64)
    if depths.size > 1 and np.any(np.diff(depths) <= 0):
        dup = depths[np.argmin(np.diff(depths))]
        raise ValueError(
            f"duplicate or non-monotone depths in slice set near "
            f"depth={dup} m")

    n_lon = ref_lon.size
    n_lat = ref_lat.size
    n_z = depths.size
    field_names = list(slices[0][3].keys())
    fields = {
        name: np.empty((n_lon, n_lat, n_z), dtype=np.float64)
        for name in field_names
    }
    for k, (_, _, _, sl_fields) in enumerate(slices):
        for name in field_names:
            fields[name][:, :, k] = sl_fields[name]

    return RawSliceStack(
        lon_grid=ref_lon, lat_grid=ref_lat, depths_m=depths, fields=fields)


# ---------------------------------------------------------------------------
# Generic CSV-grid reader (forward-looking stub for stress / pore pressure)
# ---------------------------------------------------------------------------

def read_csv_grid(path: Path, columns: list, crs: str) -> RawGrid:
    """Parse a generic CSV file with columns ``lon, lat, <field1>, ...``
    on a regular geographic grid into a ``RawGrid``.

    Parameters
    ----------
    path : Path
        CSV file. Comments starting with '#' are skipped. The first
        non-comment line must be the column header.
    columns : list of str
        Column-name list (lowercase). Must start with ``"lon", "lat"``;
        the remaining entries are field names.
    crs : str
        EPSG identifier of the (lon, lat) columns. Must equal
        ``"EPSG:4326"`` in v1 (the only currently supported source CRS).

    Returns
    -------
    RawGrid
    """
    if columns[:2] != ["lon", "lat"]:
        raise ValueError(
            f"read_csv_grid: columns must start with ['lon','lat']; "
            f"got {columns}")
    if crs != "EPSG:4326":
        raise NotImplementedError(
            f"read_csv_grid: only source_crs='EPSG:4326' is supported "
            f"in v1; got '{crs}'")
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"CSV grid not found: {p}")
    body = np.loadtxt(str(p), delimiter=",", comments="#",
                      skiprows=1)
    if body.ndim != 2 or body.shape[1] != len(columns):
        raise ValueError(
            f"{p}: body shape {body.shape} does not match "
            f"({columns}) column count {len(columns)}")
    lon_vals = body[:, 0]
    lat_vals = body[:, 1]
    lon_unique = np.unique(np.round(lon_vals, 6))
    lat_unique = np.unique(np.round(lat_vals, 6))
    n_lon = lon_unique.size
    n_lat = lat_unique.size
    if n_lon * n_lat != body.shape[0]:
        raise ValueError(
            f"{p}: not a regular grid (n_lon*n_lat = {n_lon*n_lat} "
            f"!= rows = {body.shape[0]})")
    lon_idx = np.searchsorted(lon_unique, np.round(lon_vals, 6))
    lat_idx = np.searchsorted(lat_unique, np.round(lat_vals, 6))
    fields = {}
    for ci, name in enumerate(columns[2:], start=2):
        arr = np.full((n_lon, n_lat), np.nan, dtype=np.float64)
        arr[lon_idx, lat_idx] = body[:, ci]
        if np.isnan(arr).any():
            raise ValueError(
                f"{p}: reshape produced NaN holes for field '{name}'")
        fields[name] = arr
    return RawGrid(lon_grid=lon_unique, lat_grid=lat_unique, fields=fields)
