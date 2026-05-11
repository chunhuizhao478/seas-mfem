"""crs.py — coordinate-system conversion helpers for data projection.

Canonical CRS in this project: UTM zone 11 N, WGS-84 datum (EPSG:32611).
Source CRS for the CVM-H velocity raster (and most CFM-derived rasters):
geographic WGS-84 (EPSG:4326).

Both helpers are vectorised and accept array-like input.
"""

from __future__ import annotations

import numpy as np


_TRANSFORMER_GEO_TO_UTM = None
_TRANSFORMER_UTM_TO_GEO = None


def _get_transformer_geo_to_utm():
    """Lazy singleton for pyproj.Transformer (4326 -> 32611)."""
    global _TRANSFORMER_GEO_TO_UTM
    if _TRANSFORMER_GEO_TO_UTM is None:
        try:
            from pyproj import Transformer
        except ImportError as exc:
            raise ImportError(
                "pyproj is required for geographic_to_utm11n; install via "
                "`pip install pyproj`") from exc
        _TRANSFORMER_GEO_TO_UTM = Transformer.from_crs(
            "EPSG:4326", "EPSG:32611", always_xy=True)
    return _TRANSFORMER_GEO_TO_UTM


def _get_transformer_utm_to_geo():
    """Lazy singleton for pyproj.Transformer (32611 -> 4326)."""
    global _TRANSFORMER_UTM_TO_GEO
    if _TRANSFORMER_UTM_TO_GEO is None:
        try:
            from pyproj import Transformer
        except ImportError as exc:
            raise ImportError(
                "pyproj is required for utm11n_to_geographic; install via "
                "`pip install pyproj`") from exc
        _TRANSFORMER_UTM_TO_GEO = Transformer.from_crs(
            "EPSG:32611", "EPSG:4326", always_xy=True)
    return _TRANSFORMER_UTM_TO_GEO


def geographic_to_utm11n(lon, lat):
    """Convert geographic (lon, lat) in degrees (EPSG:4326) to UTM 11 N
    (EPSG:32611) eastings/northings in metres.

    Parameters
    ----------
    lon : array_like
        Longitude in decimal degrees.
    lat : array_like
        Latitude in decimal degrees. Must broadcast against ``lon``.

    Returns
    -------
    x_utm, y_utm : numpy.ndarray
        Eastings and northings in metres, shape equal to the broadcast
        of the inputs.

    Notes
    -----
    UTM is not a pure cylindrical projection: ``x_utm`` and ``y_utm``
    each depend on BOTH ``lon`` and ``lat``. A regular geographic grid
    is therefore NOT regular in UTM; callers that need a rectilinear UTM
    grid must resample (see ``data_projection.build_velocity_cvmh``).
    """
    transformer = _get_transformer_geo_to_utm()
    lon_arr = np.asarray(lon, dtype=np.float64)
    lat_arr = np.asarray(lat, dtype=np.float64)
    if lon_arr.shape != lat_arr.shape:
        # Broadcast manually to a common shape so the transformer gets
        # 1-D flat arrays (it cannot infer broadcasting itself).
        lon_arr, lat_arr = np.broadcast_arrays(lon_arr, lat_arr)
    x_flat, y_flat = transformer.transform(lon_arr.ravel(), lat_arr.ravel())
    x = np.asarray(x_flat, dtype=np.float64).reshape(lon_arr.shape)
    y = np.asarray(y_flat, dtype=np.float64).reshape(lon_arr.shape)
    return x, y


def utm11n_to_geographic(x, y):
    """Inverse of ``geographic_to_utm11n``: UTM 11 N (m) -> (lon, lat)
    in decimal degrees (EPSG:4326).

    Parameters
    ----------
    x : array_like
        Easting in metres (UTM 11 N).
    y : array_like
        Northing in metres. Must broadcast against ``x``.

    Returns
    -------
    lon, lat : numpy.ndarray
        Longitude and latitude in decimal degrees, shape equal to the
        broadcast of the inputs.
    """
    transformer = _get_transformer_utm_to_geo()
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    if x_arr.shape != y_arr.shape:
        x_arr, y_arr = np.broadcast_arrays(x_arr, y_arr)
    lon_flat, lat_flat = transformer.transform(x_arr.ravel(), y_arr.ravel())
    lon = np.asarray(lon_flat, dtype=np.float64).reshape(x_arr.shape)
    lat = np.asarray(lat_flat, dtype=np.float64).reshape(x_arr.shape)
    return lon, lat
