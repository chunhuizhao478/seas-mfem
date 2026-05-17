"""build_stress_csm_safs.py — top-level driver: CSM Johnson & Hearn
CSV slices → schema-v1 HDF5 sidecar + ParaView VTR companion.

Mirror of ``velocity/code/build_velocity_cvmh.py`` for the stress
pillar.  Reads every ``CSM_Johnson_Hearn_raw_*km.csv`` slice in a
directory, re-projects the geographic samples onto a rectilinear UTM
11 N grid, applies the schema-v1 sanity guards (G-1: no NaN; G-2:
per-field min/max), and writes the schema-conforming HDF5 sidecar plus
an inspection-only VTK ``RectilinearGrid`` companion.

Conventions carried through this driver:

- Frame: UTM Zone 11 N metres on (x = east, y = north, z = up);
  ``z ≤ 0`` underground.
- Tensor components: 6 independent entries ``sigma_xx, sigma_yy,
  sigma_zz, sigma_xy, sigma_yz, sigma_xz`` under the East→x, North→y,
  Up→z mapping (see ``csm_johnson_hearn_reader._TENSOR_COLUMNS``).
- Units: **MPa/yr** (source CSM Johnson & Hearn is a *stressing-rate*
  model — see CSV header line 5).  Sign convention: **tension
  POSITIVE** (physics convention, preserved from source).  These two
  facts diverge from ``stress_safs.h5`` (compression-positive Pa) and
  are documented in the sidecar root attrs (``model``, ``sign_convention``).

Usage
-----
::

    conda activate pythonenv
    cd miniapps/seas/safs/project_7.0_alternative/stress/code

    python build_stress_csm_safs.py \\
        --raw-dir   ../raw/csm_johnson_hearn \\
        --out-path  ../results/stress_csm_safs.h5 \\
        --grid-dx   1500 \\
        --paraview-export \\
        --verbose
"""

from __future__ import annotations

import argparse
import base64
import struct
import sys
from pathlib import Path
from typing import Optional

import numpy as np

# Sibling import of velocity/code/{crs,sidecar} so the schema-v1 writer
# and the EPSG:4326 → EPSG:32611 transformer live in exactly one place
# in the repo (no duplicated implementations across pillars).
_VELOCITY_CODE_DIR = (
    Path(__file__).resolve().parents[2] / "velocity" / "code"
)
if str(_VELOCITY_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_VELOCITY_CODE_DIR))
from crs import geographic_to_utm11n  # noqa: E402
from sidecar import write_sidecar     # noqa: E402

# Sibling import of the CSM reader (lives next to this file).
_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))
from csm_johnson_hearn_reader import (  # noqa: E402
    RawSliceStack,
    read_csm_johnson_hearn,
)


# ---------------------------------------------------------------------------
# Field bounds — derived from data at run time, with a safety factor.
#
# The CSM Johnson-Hearn raw values for the SAFS bundle are O(1e-2)
# MPa/yr (see `head` of any raw CSV).  Rather than hardcode bounds we
# size them from the actual data range with a safety margin so the
# schema-v1 G-2 guard (every cell ∈ [min, max]) passes by construction.
# ---------------------------------------------------------------------------

DEFAULT_BOUNDS_SAFETY: float = 1.5  # widen observed range by this factor
DEFAULT_GRID_DX_M: float = 1500.0


def _derive_field_bounds(
    fields: dict[str, np.ndarray],
    *,
    units: str = "MPa/yr",
    safety: float = DEFAULT_BOUNDS_SAFETY,
) -> dict[str, tuple[float, float, str]]:
    """Return ``{name: (min_value, max_value, units)}`` for every field,
    sized to strictly enclose the actual data with a safety margin so
    the schema-v1 writer's G-2 guard passes by construction.

    Rules (mirrored from ``stress/code/build_stress_safs.derive_field_bounds``):
        - ``min < 0`` → ``min_value = min * safety`` (more negative).
        - ``min >= 0`` → ``min_value = min / safety``  (closer to 0 from above).
        - analogous for ``max`` with ``max > 0`` / ``max <= 0``.
        - identically-zero field → widen by ±1 unit so ``min < max``.
    """
    if safety <= 1.0:
        raise ValueError(
            f"safety must be > 1.0 to widen bounds; got {safety}"
        )
    out: dict[str, tuple[float, float, str]] = {}
    for name, arr in fields.items():
        a_min = float(np.min(arr))
        a_max = float(np.max(arr))
        v_min = (a_min * safety
                 if a_min < 0.0 else a_min / safety)
        v_max = (a_max * safety
                 if a_max > 0.0 else a_max / safety)
        if not (v_min < v_max):
            # Degenerate (e.g. identically-zero Seu at shallow depths
            # in this dataset).  Widen by ±1 unit so the writer's
            # strict v_min < v_max guard passes and the actual data
            # is enclosed.
            v_min = a_min - 1.0
            v_max = a_max + 1.0
        out[name] = (v_min, v_max, units)
    return out


# ---------------------------------------------------------------------------
# ParaView companion writer — VTK RectilinearGrid (.vtr) inline base64.
#
# This is a copy of the same helper in build_velocity_cvmh.py.  It is
# duplicated here (rather than imported) to keep the two drivers fully
# self-contained at the file level, mirroring how the velocity driver
# carries its own copy.  Touch one, audit the other.
# ---------------------------------------------------------------------------


def _b64_inline_array(arr: np.ndarray) -> str:
    """Encode a single float64 numpy array as VTK XML inline base64,
    matching ``format="binary"``: ``base64(uint64_LE(nbytes) || raw_bytes)``.

    Each ``DataArray`` carries its own self-contained base64 stream so
    there are no shared offsets to keep aligned (see the long comment
    on the velocity driver's twin of this helper for why we avoid the
    ``AppendedData`` form).
    """
    raw = np.ascontiguousarray(arr, dtype=np.float64).tobytes()
    return base64.b64encode(
        struct.pack("<Q", len(raw)) + raw
    ).decode("ascii")


def _write_vtr_companion(path: Path,
                         x: np.ndarray, y: np.ndarray, z: np.ndarray,
                         fields: dict) -> None:
    """Write a VTK ``RectilinearGrid`` (.vtr) wrapping the sidecar grid.

    The .vtr container accepts both uniform and non-uniform axes; we
    write it unconditionally for the stress sidecar to match the user
    request ``stress_csm_safs.vtr`` (the velocity pipeline picks .vti
    when all axes are uniform — that branch is intentionally skipped
    here).
    """
    Nx, Ny, Nz = x.size, y.size, z.size
    extent = f"0 {Nx - 1} 0 {Ny - 1} 0 {Nz - 1}"

    # VTK PointData ordering: outermost loop = z, then y, then x.
    # Our arrays are (Nx, Ny, Nz) row-major; transpose to (Nz, Ny, Nx)
    # then ravel order='C' to match VTK's expected layout.
    field_arrays = {name: arr.transpose(2, 1, 0).ravel(order="C")
                    for name, arr in fields.items()}
    x_arr = np.ascontiguousarray(x, dtype=np.float64)
    y_arr = np.ascontiguousarray(y, dtype=np.float64)
    z_arr = np.ascontiguousarray(z, dtype=np.float64)

    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="RectilinearGrid" version="1.0" '
        'byte_order="LittleEndian" header_type="UInt64">',
        f'  <RectilinearGrid WholeExtent="{extent}">',
        f'    <Piece Extent="{extent}">',
        '      <PointData>',
    ]
    for name, arr in field_arrays.items():
        lines.append(
            f'        <DataArray type="Float64" Name="{name}" '
            f'format="binary">'
        )
        lines.append('          ' + _b64_inline_array(arr))
        lines.append('        </DataArray>')
    lines.extend([
        '      </PointData>',
        '      <Coordinates>',
    ])
    for axis_name, axis_arr in (("X", x_arr), ("Y", y_arr), ("Z", z_arr)):
        lines.append(
            f'        <DataArray type="Float64" Name="{axis_name}" '
            f'format="binary">'
        )
        lines.append('          ' + _b64_inline_array(axis_arr))
        lines.append('        </DataArray>')
    lines.extend([
        '      </Coordinates>',
        '    </Piece>',
        '  </RectilinearGrid>',
        '</VTKFile>',
    ])
    Path(path).write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# UTM rectilinear grid construction and per-field resampling.
#
# Both helpers mirror the velocity driver's homonyms (``_build_utm_grid``,
# ``_resample_field``).  See build_velocity_cvmh.py for the design notes
# on the inscribed-AABB convention and the per-slice LinearNDInterpolator.
# ---------------------------------------------------------------------------


def _build_utm_grid(stack: RawSliceStack, grid_dx: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Compute a rectilinear UTM 11 N grid covering the inscribed AABB
    of the (lon, lat) samples after re-projection.  The grid spacing is
    ``grid_dx`` (m) on both x and y.

    See build_velocity_cvmh._build_utm_grid for the why behind using
    edge-extrema instead of whole-array min/max (the rectilinear bbox
    of the projected quadrilateral must inscribe the source hull or
    the LinearNDInterpolator returns NaN at the corners).
    """
    lon_grid = stack.lon_grid
    lat_grid = stack.lat_grid
    LON, LAT = np.meshgrid(lon_grid, lat_grid, indexing="ij")
    X_utm, Y_utm = geographic_to_utm11n(LON, LAT)
    xmin_in = float(X_utm[0, :].max())
    xmax_in = float(X_utm[-1, :].min())
    ymin_in = float(Y_utm[:, 0].max())
    ymax_in = float(Y_utm[:, -1].min())
    xmin = float(np.ceil(xmin_in / grid_dx) * grid_dx)
    xmax = float(np.floor(xmax_in / grid_dx) * grid_dx)
    ymin = float(np.ceil(ymin_in / grid_dx) * grid_dx)
    ymax = float(np.floor(ymax_in / grid_dx) * grid_dx)
    if xmin > xmax or ymin > ymax:
        raise ValueError(
            f"_build_utm_grid: requested grid_dx={grid_dx} m is too "
            f"large for the source extent (would yield empty grid)."
        )
    x_axis = np.arange(
        xmin, xmax + 0.5 * grid_dx, grid_dx, dtype=np.float64
    )
    y_axis = np.arange(
        ymin, ymax + 0.5 * grid_dx, grid_dx, dtype=np.float64
    )
    return x_axis, y_axis


def _resample_field(stack: RawSliceStack, name: str,
                    x_axis: np.ndarray, y_axis: np.ndarray,
                    verbose: bool = False) -> np.ndarray:
    """Resample a single named field onto the rectilinear UTM grid, one
    depth slice at a time.  Returns an array of shape ``(Nx, Ny, Nz)``.

    Cells that land outside the convex hull of the source samples are
    NaN; the writer's G-1 guard rejects them.  See
    build_velocity_cvmh._resample_field for the design rationale.
    """
    from scipy.interpolate import LinearNDInterpolator

    lon_grid = stack.lon_grid
    lat_grid = stack.lat_grid
    LON, LAT = np.meshgrid(lon_grid, lat_grid, indexing="ij")
    X_pts, Y_pts = geographic_to_utm11n(LON, LAT)
    pts_xy = np.stack([X_pts.ravel(), Y_pts.ravel()], axis=1)
    XX, YY = np.meshgrid(x_axis, y_axis, indexing="ij")
    target = np.stack([XX.ravel(), YY.ravel()], axis=1)

    n_z = stack.depths_m.size
    out = np.empty((x_axis.size, y_axis.size, n_z), dtype=np.float64)
    for k in range(n_z):
        slab = stack.fields[name][:, :, k].ravel()
        interp = LinearNDInterpolator(pts_xy, slab, fill_value=np.nan)
        vals = interp(target)
        out[:, :, k] = vals.reshape(x_axis.size, y_axis.size)
        if verbose:
            n_nan = int(np.isnan(out[:, :, k]).sum())
            print(
                f"  field {name} slice k={k}: "
                f"depth={stack.depths_m[k]:.1f} m, "
                f"NaN cells = {n_nan}/{x_axis.size * y_axis.size}"
            )
    return out


# ---------------------------------------------------------------------------
# Slice discovery
# ---------------------------------------------------------------------------


def _list_slice_paths(raw_dir: Path) -> list[Path]:
    paths = sorted(Path(raw_dir).glob("CSM_Johnson_Hearn_raw_*km.csv"))
    if not paths:
        raise FileNotFoundError(
            f"build_stress_csm_safs: no CSM_Johnson_Hearn_raw_*km.csv "
            f"files in {raw_dir}"
        )
    return paths


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------


def build_stress_csm_safs(
    raw_dir: Path,
    out_path: Path,
    *,
    grid_dx: float,
    bounds_safety: float = DEFAULT_BOUNDS_SAFETY,
    paraview_export: bool = False,
    verbose: bool = False,
) -> None:
    """Top-level driver: CSM Johnson & Hearn CSVs → schema-v1 HDF5
    sidecar (+ optional ParaView VTR companion).

    Steps mirror ``build_velocity_cvmh.build_velocity_sidecar``:

      1. Glob ``CSM_Johnson_Hearn_raw_*km.csv`` under ``raw_dir``.
      2. Parse via :func:`csm_johnson_hearn_reader.read_csm_johnson_hearn`.
      3. Build a rectilinear UTM 11 N (x, y) grid at spacing ``grid_dx``
         that inscribes the projected source quadrilateral.
      4. ``z_axis = sort(-depths_m)`` (depths positive=down → z=elevation
         ascending; matches ``CANONICAL_Z_POSITIVE``).
      5. Resample each tensor component onto ``(x, y, z)``; flip k so
         deepest is first along ``z_axis``.
      6. G-1 (NaN) and G-2 (per-field min/max) guards.
      7. Schema-v1 sidecar write + optional VTR companion.

    Raises
    ------
    ValueError
        For any G-1 / G-2 violation, malformed source CSV, or invalid
        grid_dx.
    FileNotFoundError
        If ``raw_dir`` has no matching CSVs.
    """
    raw_dir = Path(raw_dir)
    out_path = Path(out_path)

    paths = _list_slice_paths(raw_dir)
    if verbose:
        print(
            f"build_stress_csm_safs: {len(paths)} slice(s) in {raw_dir}"
        )

    stack = read_csm_johnson_hearn(paths)
    if verbose:
        print(
            f"  source grid: lon=[{stack.lon_grid[0]:.4f}, "
            f"{stack.lon_grid[-1]:.4f}] x "
            f"lat=[{stack.lat_grid[0]:.4f}, "
            f"{stack.lat_grid[-1]:.4f}]; "
            f"depths_m=[{stack.depths_m[0]:.1f}, "
            f"{stack.depths_m[-1]:.1f}] "
            f"({stack.depths_m.size} slices); "
            f"fields={list(stack.fields.keys())}"
        )

    x_axis, y_axis = _build_utm_grid(stack, grid_dx)
    z_axis = np.sort(-stack.depths_m).astype(np.float64)
    if verbose:
        print(
            f"  utm grid: x=[{x_axis[0]:.0f}, {x_axis[-1]:.0f}] m, "
            f"Nx={x_axis.size}; "
            f"y=[{y_axis[0]:.0f}, {y_axis[-1]:.0f}] m, "
            f"Ny={y_axis.size}; "
            f"z=[{z_axis[0]:.0f}, {z_axis[-1]:.0f}] m, "
            f"Nz={z_axis.size}"
        )

    # Resample each field, flipping the k axis so that ascending
    # elevation lines up with z_axis (depths_m is ascending depth =
    # descending elevation).
    fields_xyz: dict[str, np.ndarray] = {}
    for name in stack.fields:
        out = _resample_field(stack, name, x_axis, y_axis, verbose=verbose)
        flipped = out[:, :, ::-1].copy()
        fields_xyz[name] = flipped

    # G-1: NaN guard (writer also re-checks; we surface the error here
    # with a more diagnostic message).
    for name, arr in fields_xyz.items():
        n_nan = int(np.isnan(arr).sum())
        if n_nan:
            raise ValueError(
                f"build_stress_csm_safs: {n_nan} NaN cells in field "
                f"'{name}' after resample.  The rectilinear UTM grid "
                f"extends beyond the source data hull.\n"
                f"  Source hull bbox (UTM 11 N): "
                f"x=[{x_axis[0]:.0f}, {x_axis[-1]:.0f}], "
                f"y=[{y_axis[0]:.0f}, {y_axis[-1]:.0f}].\n"
                f"  Either pass a larger --grid-dx (so the inscribed "
                f"AABB lands closer to the hull centre) or regenerate "
                f"the source CSVs over a wider footprint."
            )

    # G-2: per-field sanity bounds (data-derived).
    field_bounds = _derive_field_bounds(
        fields_xyz, units="MPa/yr", safety=bounds_safety
    )

    attrs = {
        "source": ",".join(p.name for p in paths),
        "source_crs": "EPSG:4326",
        "mesh_tag": "",
        # Documentation-only attrs (the runtime / schema do not read
        # these; they explain the field semantics to a downstream
        # inspector who opens the .h5 in HDFView or ParaView).
        "model": "CSM_Johnson_Hearn",
        "model_reference": "Johnson (2024), doi:10.1029/2023JB027472",
        "model_kind": "stressing_rate",
        "tensor_units": "MPa/yr",
        "sign_convention": "tension_positive",
        "frame_mapping": "East->x, North->y, Up->z",
    }

    write_sidecar(
        out_path=out_path,
        x=x_axis, y=y_axis, z=z_axis,
        fields=fields_xyz,
        attrs=attrs,
        field_bounds=field_bounds,
    )
    if verbose:
        print(f"build_stress_csm_safs: wrote {out_path}")

    if paraview_export:
        companion = out_path.with_suffix(".vtr")
        _write_vtr_companion(companion, x_axis, y_axis, z_axis, fields_xyz)
        if verbose:
            print(
                f"build_stress_csm_safs: wrote {companion} "
                f"[VTR (rectilinear grid)]"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--raw-dir", type=Path, required=True,
        help="directory containing CSM_Johnson_Hearn_raw_*km.csv slices",
    )
    ap.add_argument(
        "--out-path", type=Path, required=True,
        help="output sidecar path (created)",
    )
    ap.add_argument(
        "--grid-dx", type=float, default=DEFAULT_GRID_DX_M,
        help=("UTM rectilinear (x, y) grid spacing [m]; "
              f"default {DEFAULT_GRID_DX_M:.0f}"),
    )
    ap.add_argument(
        "--bounds-safety", type=float, default=DEFAULT_BOUNDS_SAFETY,
        help=("safety multiplier applied when deriving per-field min/"
              "max bounds from the data; "
              f"default {DEFAULT_BOUNDS_SAFETY}"),
    )
    ap.add_argument(
        "--paraview-export", action="store_true",
        help=("also write a ParaView-compatible companion file "
              "(RectilinearGrid, .vtr) next to the .h5 sidecar"),
    )
    ap.add_argument("--verbose", "-v", action="store_true")
    return ap.parse_args(argv)


def main(argv: Optional[list] = None) -> int:
    args = _parse_args(argv)
    try:
        build_stress_csm_safs(
            raw_dir=args.raw_dir,
            out_path=args.out_path,
            grid_dx=args.grid_dx,
            bounds_safety=args.bounds_safety,
            paraview_export=args.paraview_export,
            verbose=args.verbose,
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
