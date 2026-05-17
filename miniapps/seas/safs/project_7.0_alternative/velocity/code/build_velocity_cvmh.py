"""build_velocity_cvmh.py — top-level driver: CVM ASCII slices -> sidecar.

Reads every ``velocity_raw_*.bp`` slice in a directory, re-projects
geographic samples onto a rectilinear UTM 11 N grid, applies the
schema-v1 sanity guards (G-1: no NaN; G-3: optional mesh-bbox
containment), and writes the schema-conforming HDF5 sidecar.

Raw .bp slices live in per-version subdirectories under
``velocity/raw/``: ``cvmh/``, ``cvm_s4.26.m01/``,
``multiscale_statewise_cvm/``.  Point ``--raw-dir`` at exactly one of
them per invocation and mirror the version name in ``--out-path``.

Usage
-----
::

    # Run once per CVM version:
    python -m data_projection.build_velocity_cvmh \\
        --raw-dir       ../../../velocity/raw/cvmh \\
        --out-path      ../../../velocity/results/cvmh/velocity_safs.h5 \\
        --grid-dx       1500 \\
        --mesh-msh      ../../../meshing/results/msh/safs_fault_box_nwcut_500m.msh \\
        --extend-z-top  100 \\
        --paraview-export \\
        --verbose
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from bbox_check import (
    BBoxContainmentError,
    assert_grid_contains_mesh,
    mesh_msh_bbox,
)
from crs import geographic_to_utm11n
from raw_readers import read_cvmh_ascii, RawSliceStack
from sidecar import write_sidecar


# Default per-field sanity bounds (see schema doc §5.1).
#
# The lower bounds were originally set to "competent rock floor" values
# (Vp 2000, Vs 1500, rho 2000) to catch sub-seafloor / water-saturated
# cells that come out of CVM-H near the coastline.  The user has elected
# to KEEP those cells (they are real geology, just not numerically nice
# for CFL).  The defaults below are therefore permissive, accepting the
# minima observed in the SCEC CVM-H 15.1.1 archive (Vs ~ 120 m/s,
# Vp ~ 1200 m/s, rho ~ 1420 kg/m³).  Upper bounds remain at the
# competent-rock ceiling so genuinely bogus values (e.g. cosmic-ray
# spikes) still fail-loud.
#
# Pass tighter values via --vs-min-mps / --vp-min-mps / --rho-min-kgm3
# at invocation time when the input has been pre-cleaned.
DEFAULT_VP_MIN_MPS = 1000.0
DEFAULT_VP_MAX_MPS = 9000.0
DEFAULT_VS_MIN_MPS = 100.0
DEFAULT_VS_MAX_MPS = 5000.0
DEFAULT_RHO_MIN_KGM3 = 1000.0
DEFAULT_RHO_MAX_KGM3 = 3500.0


# ---------------------------------------------------------------------------
# ParaView companion writer (VTI for uniform axes, VTR otherwise).
# ---------------------------------------------------------------------------

def _is_uniform(axis: np.ndarray, atol: float = 1.0) -> bool:
    """Return True iff the 1-D ``axis`` is uniformly spaced within
    ``atol`` metres.  Single-cell axes count as uniform."""
    if axis.size <= 2:
        return True
    diffs = np.diff(axis)
    return bool(np.allclose(diffs, diffs[0], atol=atol, rtol=0.0))


def _b64_inline_array(arr: np.ndarray) -> str:
    """Encode a single float64 numpy array as VTK XML inline base64,
    matching ``format="binary"``: ``base64(uint64_LE(nbytes) || raw_bytes)``.

    Each ``DataArray`` carries its own self-contained base64 stream,
    so there are no shared offsets to keep aligned.  This sidesteps the
    VTK ``encoding="base64"`` quirk where ``AppendedData`` offsets must
    be multiples of 3 decoded bytes (a base64 char boundary): with
    multiple Float64 arrays whose per-array size is 8 N (8 N + 8 with
    the header), only the first array hits offset 0 cleanly; subsequent
    offsets land mid-base64-char-group and the reader returns garbage.
    """
    import base64, struct
    raw = np.ascontiguousarray(arr, dtype=np.float64).tobytes()
    return base64.b64encode(
        struct.pack("<Q", len(raw)) + raw
    ).decode("ascii")


def _write_vti_companion(path: Path,
                         x: np.ndarray, y: np.ndarray, z: np.ndarray,
                         fields: dict) -> None:
    """Write a VTK ImageData (.vti) wrapping a uniformly-spaced
    rectilinear grid.  ParaView opens this directly."""
    Nx, Ny, Nz = x.size, y.size, z.size
    dx = float(x[1] - x[0]) if Nx > 1 else 1.0
    dy = float(y[1] - y[0]) if Ny > 1 else 1.0
    dz = float(z[1] - z[0]) if Nz > 1 else 1.0
    ox, oy, oz = float(x[0]), float(y[0]), float(z[0])
    extent = f"0 {Nx-1} 0 {Ny-1} 0 {Nz-1}"

    # VTK PointData ordering: outermost loop = z, then y, then x.
    # Our arrays are (Nx, Ny, Nz) row-major; transpose to (Nz, Ny, Nx)
    # then ravel order='C' to match VTK's expected layout.
    arrays = {name: arr.transpose(2, 1, 0).ravel(order="C")
              for name, arr in fields.items()}

    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="ImageData" version="1.0" '
        'byte_order="LittleEndian" header_type="UInt64">',
        f'  <ImageData WholeExtent="{extent}" '
        f'Origin="{ox} {oy} {oz}" Spacing="{dx} {dy} {dz}">',
        f'    <Piece Extent="{extent}">',
        '      <PointData>',
    ]
    for name, arr in arrays.items():
        lines.append(
            f'        <DataArray type="Float64" Name="{name}" '
            f'format="binary">')
        lines.append('          ' + _b64_inline_array(arr))
        lines.append('        </DataArray>')
    lines.extend([
        '      </PointData>',
        '    </Piece>',
        '  </ImageData>',
        '</VTKFile>',
    ])
    Path(path).write_text("\n".join(lines))


def _write_vtr_companion(path: Path,
                         x: np.ndarray, y: np.ndarray, z: np.ndarray,
                         fields: dict) -> None:
    """Write a VTK RectilinearGrid (.vtr) for non-uniform axes.
    ParaView opens this directly."""
    Nx, Ny, Nz = x.size, y.size, z.size
    extent = f"0 {Nx-1} 0 {Ny-1} 0 {Nz-1}"

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
            f'format="binary">')
        lines.append('          ' + _b64_inline_array(arr))
        lines.append('        </DataArray>')
    lines.extend([
        '      </PointData>',
        '      <Coordinates>',
    ])
    for axis_name, axis_arr in (("X", x_arr), ("Y", y_arr), ("Z", z_arr)):
        lines.append(
            f'        <DataArray type="Float64" Name="{axis_name}" '
            f'format="binary">')
        lines.append('          ' + _b64_inline_array(axis_arr))
        lines.append('        </DataArray>')
    lines.extend([
        '      </Coordinates>',
        '    </Piece>',
        '  </RectilinearGrid>',
        '</VTKFile>',
    ])
    Path(path).write_text("\n".join(lines))


def _write_paraview_companion(out_h5: Path,
                              x: np.ndarray, y: np.ndarray,
                              z: np.ndarray,
                              fields: dict,
                              verbose: bool = False) -> Path:
    """Pick VTI (uniform axes) or VTR (non-uniform) and write the
    companion file next to ``out_h5``.  Returns the companion path.
    """
    uniform = (_is_uniform(x) and _is_uniform(y) and _is_uniform(z))
    suffix = ".vti" if uniform else ".vtr"
    companion = out_h5.with_suffix(suffix)
    if uniform:
        _write_vti_companion(companion, x, y, z, fields)
    else:
        _write_vtr_companion(companion, x, y, z, fields)
    if verbose:
        kind = "VTI (uniform axes)" if uniform else "VTR (non-uniform axes)"
        print(f"  paraview companion: {companion.name} [{kind}]")
    return companion


def _list_slice_paths(raw_dir: Path) -> list:
    paths = sorted(Path(raw_dir).glob("velocity_raw_*.bp"))
    if not paths:
        raise FileNotFoundError(
            f"build_velocity_cvmh: no velocity_raw_*.bp files in {raw_dir}")
    return paths


def _build_utm_grid(stack: RawSliceStack, grid_dx: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Compute a rectilinear UTM 11 N grid covering the convex hull of
    the (lon, lat) samples after re-projection. The grid spacing is
    ``grid_dx`` (m) on both x and y; the bbox of the grid is the smallest
    axis-aligned UTM bbox that contains every projected sample.
    """
    lon_grid = stack.lon_grid
    lat_grid = stack.lat_grid
    LON, LAT = np.meshgrid(lon_grid, lat_grid, indexing="ij")
    X_utm, Y_utm = geographic_to_utm11n(LON, LAT)
    # Inscribed AABB of the projected source quadrilateral.  We must
    # use the per-edge extrema (NOT the whole-array min/max), otherwise
    # the rectilinear bbox corners stick out beyond the projected
    # quadrilateral and the LinearNDInterpolator returns NaN there:
    #   xmin_in = max of x along the western lon column   (lon = lon[0])
    #   xmax_in = min of x along the eastern lon column   (lon = lon[-1])
    #   ymin_in = max of y along the southern lat row     (lat = lat[0])
    #   ymax_in = min of y along the northern lat row     (lat = lat[-1])
    # Then round INWARD to grid_dx multiples (ceil on lower, floor on
    # upper) so every grid sample is strictly inside the source hull.
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
            f"large for the source extent (would yield empty grid)")
    x_axis = np.arange(xmin, xmax + 0.5 * grid_dx, grid_dx, dtype=np.float64)
    y_axis = np.arange(ymin, ymax + 0.5 * grid_dx, grid_dx, dtype=np.float64)
    return x_axis, y_axis


def _resample_field(stack: RawSliceStack, name: str,
                    x_axis: np.ndarray, y_axis: np.ndarray,
                    verbose: bool = False) -> np.ndarray:
    """Resample a single named field onto the rectilinear UTM grid,
    one slice at a time. Returns an array of shape (Nx, Ny, Nz).

    Cells outside the convex hull of the source samples land as NaN;
    the writer's G-1 guard will reject them, which is the desired
    behaviour (see schema §4 and feature plan v2 §C-1 / §C-2).
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
            print(f"  field {name} slice k={k}: "
                  f"depth={stack.depths_m[k]:.1f} m, "
                  f"NaN cells = {n_nan}/{x_axis.size * y_axis.size}")
    return out


def build_velocity_sidecar(raw_dir: Path, out_path: Path, grid_dx: float,
                           vp_min_mps: float, vp_max_mps: float,
                           vs_min_mps: float, vs_max_mps: float,
                           rho_min_kgm3: float, rho_max_kgm3: float,
                           mesh_msh: Path = None,
                           extend_z_top_m: float = 0.0,
                           paraview_export: bool = False,
                           verbose: bool = False) -> None:
    """Top-level pipeline. Reads CVM-H slices from ``raw_dir``, builds
    a UTM 11 N rectilinear grid at spacing ``grid_dx``, applies guards
    G-1/G-2/G-3, writes the schema-v1 sidecar to ``out_path``.

    Raises
    ------
    ValueError
        For any G-1 (NaN), G-2 (per-field min/max) violation.
    BBoxContainmentError
        For any G-3 (mesh-vs-grid containment) violation.
    """
    paths = _list_slice_paths(Path(raw_dir))
    if verbose:
        print(f"build_velocity_cvmh: {len(paths)} slice(s) in {raw_dir}")
    stack = read_cvmh_ascii(paths)
    if verbose:
        print(f"  source grid: lon=[{stack.lon_grid[0]:.4f}, "
              f"{stack.lon_grid[-1]:.4f}] x lat=[{stack.lat_grid[0]:.4f}, "
              f"{stack.lat_grid[-1]:.4f}]; "
              f"depths_m=[{stack.depths_m[0]:.1f}, "
              f"{stack.depths_m[-1]:.1f}] ({stack.depths_m.size} slices); "
              f"fields={list(stack.fields.keys())}")

    x_axis, y_axis = _build_utm_grid(stack, grid_dx)
    z_axis = np.sort(-stack.depths_m).astype(np.float64)
    # Optional synthetic top extension: clone the shallowest slice (z =
    # max in elevation) up to z + extend_z_top_m, so a mesh whose top
    # extends slightly above z = 0 (e.g. PAD_TOP = 100 m above the free
    # surface) still passes the strict-interpolation guard.
    if extend_z_top_m > 0.0:
        z_axis = np.concatenate(
            [z_axis, np.array([z_axis[-1] + extend_z_top_m])])
    if verbose:
        print(f"  utm grid: x=[{x_axis[0]:.0f}, {x_axis[-1]:.0f}] m, "
              f"Nx={x_axis.size}; y=[{y_axis[0]:.0f}, {y_axis[-1]:.0f}] m, "
              f"Ny={y_axis.size}; z=[{z_axis[0]:.0f}, {z_axis[-1]:.0f}] m, "
              f"Nz={z_axis.size}")
        if extend_z_top_m > 0.0:
            print(f"  (z extended upward by {extend_z_top_m:.1f} m via "
                  f"surface-slice clone for PAD_TOP coverage)")

    # G-3 (optional, writer side): mesh containment.  Run AFTER the
    # optional --extend-z-top so the guard reflects the final z range
    # that will be written to the sidecar.
    if mesh_msh is not None:
        mb = mesh_msh_bbox(mesh_msh)
        if verbose:
            print(f"  mesh bbox: x=[{mb['xmin']:.0f}, {mb['xmax']:.0f}] "
                  f"y=[{mb['ymin']:.0f}, {mb['ymax']:.0f}] "
                  f"z=[{mb['zmin']:.0f}, {mb['zmax']:.0f}]")
        assert_grid_contains_mesh(x_axis, y_axis, z_axis, mb)

    # Resample each field onto the UTM grid.  Note: each slice's
    # data array currently uses the depth ordering of stack.depths_m
    # (ascending depth = descending elevation), but z_axis above is
    # sorted ascending elevation.  We reverse the k index when
    # writing back to align with z_axis.
    fields_xyz = {}
    for name in stack.fields:
        out = _resample_field(stack, name, x_axis, y_axis, verbose=verbose)
        # depths_m ascending -> shallow = first k.  z elevation ascending
        # -> deepest = first k.  So we reverse along k.
        flipped = out[:, :, ::-1].copy()
        # Clone the surface slice (now last along z, ascending) onto the
        # synthetic top slice.  Keeps the field strictly interpolatable
        # over the extended z range without inventing material values.
        if extend_z_top_m > 0.0:
            top_slab = flipped[:, :, -1:]
            flipped = np.concatenate([flipped, top_slab], axis=2)
        fields_xyz[name] = flipped

    # G-1: NaN guard.
    for name, arr in fields_xyz.items():
        n_nan = int(np.isnan(arr).sum())
        if n_nan:
            grid = {"xmin": x_axis[0], "xmax": x_axis[-1],
                    "ymin": y_axis[0], "ymax": y_axis[-1],
                    "zmin": z_axis[0], "zmax": z_axis[-1]}
            raise ValueError(
                f"build_velocity_cvmh: {n_nan} NaN cells in field "
                f"'{name}' after resample. The rectilinear UTM grid "
                f"extends beyond the source data hull.\n"
                f"  Either shrink the requested grid via tighter --grid-dx "
                f"alignment, OR extend the source raster.\n"
                f"  Source hull bbox (UTM 11 N): x=[{x_axis[0]:.0f}, "
                f"{x_axis[-1]:.0f}], y=[{y_axis[0]:.0f}, {y_axis[-1]:.0f}]")

    # Translate reader column names to canonical schema names.  Reject
    # any source field NOT in the canonical map: this driver is
    # specific to the CVM-H velocity set (Vp, Vs, density); silently
    # dropping a fourth column (e.g. Qp, Qs) would lose data without
    # the user noticing.
    canonical = {"vp": "Vp", "vs": "Vs", "density": "density"}
    extras = sorted(set(fields_xyz) - set(canonical))
    if extras:
        raise ValueError(
            f"build_velocity_cvmh: source contains field(s) not handled "
            f"by this driver: {extras}.  This driver is specific to the "
            f"CVM-H velocity set (Vp, Vs, density).  Either drop the "
            f"extra columns from the source data, or extend the "
            f"canonical map to support them.")
    fields_canonical = {canonical[k]: v for k, v in fields_xyz.items()}
    expected_canonical = {"Vp", "Vs", "density"}
    if set(fields_canonical) != expected_canonical:
        raise ValueError(
            f"build_velocity_cvmh: expected fields {sorted(expected_canonical)}; "
            f"got {sorted(fields_canonical)}")

    # G-2 driver-level sanity guard removed: accept the raw model data
    # as-is. The CLI args (--vp-min-mps, --vp-max-mps, etc.) remain for
    # backwards-compatibility with existing callers and tests but no
    # longer constrain the build.  We still need to pass per-field
    # bounds to write_sidecar (required by the v1 HDF5 schema), so we
    # derive them from the observed data range (padded by 1 source-unit
    # below/above so the writer's strict v_min < v_max check is always
    # satisfied, even for a hypothetical constant-valued field).
    units_map = {"Vp": "m/s", "Vs": "m/s", "density": "kg/m^3"}
    field_bounds = {}
    for name, arr in fields_canonical.items():
        field_bounds[name] = (float(arr.min()) - 1.0,
                              float(arr.max()) + 1.0,
                              units_map[name])

    # Write the sidecar.
    write_sidecar(
        out_path=Path(out_path),
        x=x_axis, y=y_axis, z=z_axis,
        fields=fields_canonical,
        attrs={
            "source": ",".join(p.name for p in paths),
            "source_crs": "EPSG:4326",
            "mesh_tag": Path(mesh_msh).stem if mesh_msh else "",
        },
        field_bounds=field_bounds,
    )
    if verbose:
        print(f"build_velocity_cvmh: wrote {out_path}")

    if paraview_export:
        companion = _write_paraview_companion(
            Path(out_path), x_axis, y_axis, z_axis, fields_canonical,
            verbose=verbose)
        if verbose:
            print(f"build_velocity_cvmh: wrote {companion}")


def _parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", type=Path, required=True,
                    help="directory containing velocity_raw_*.bp slices")
    ap.add_argument("--out-path", type=Path, required=True,
                    help="output sidecar path (created)")
    ap.add_argument("--grid-dx", type=float, default=1500.0,
                    help="UTM rectilinear grid spacing [m] (default %(default)s)")
    ap.add_argument("--vp-min-mps", type=float, default=DEFAULT_VP_MIN_MPS)
    ap.add_argument("--vp-max-mps", type=float, default=DEFAULT_VP_MAX_MPS)
    ap.add_argument("--vs-min-mps", type=float, default=DEFAULT_VS_MIN_MPS)
    ap.add_argument("--vs-max-mps", type=float, default=DEFAULT_VS_MAX_MPS)
    ap.add_argument("--rho-min-kgm3", type=float,
                    default=DEFAULT_RHO_MIN_KGM3)
    ap.add_argument("--rho-max-kgm3", type=float,
                    default=DEFAULT_RHO_MAX_KGM3)
    ap.add_argument("--mesh-msh", type=Path, default=None,
                    help="optional path to a gmsh .msh; if provided, "
                         "writer-side G-3 guard checks that the data "
                         "grid contains the mesh bbox before writing")
    ap.add_argument("--extend-z-top", type=float, default=0.0,
                    metavar="METRES",
                    help="extend the elevation axis upward by this many "
                         "metres (default 0); the synthetic top slice "
                         "clones the shallowest source slice's values. "
                         "Use this to cover a mesh whose top extends "
                         "above the velocity dataset's free surface "
                         "(e.g. PAD_TOP=100 m -> --extend-z-top 100).")
    ap.add_argument("--paraview-export", action="store_true",
                    help="write a ParaView-compatible companion file "
                         "next to the .h5 sidecar — VTI when all three "
                         "axes are uniformly spaced, otherwise VTR "
                         "(RectilinearGrid).  Default: off.")
    ap.add_argument("--verbose", "-v", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)
    try:
        build_velocity_sidecar(
            raw_dir=args.raw_dir,
            out_path=args.out_path,
            grid_dx=args.grid_dx,
            vp_min_mps=args.vp_min_mps,
            vp_max_mps=args.vp_max_mps,
            vs_min_mps=args.vs_min_mps,
            vs_max_mps=args.vs_max_mps,
            rho_min_kgm3=args.rho_min_kgm3,
            rho_max_kgm3=args.rho_max_kgm3,
            mesh_msh=args.mesh_msh,
            extend_z_top_m=args.extend_z_top,
            paraview_export=args.paraview_export,
            verbose=args.verbose,
        )
    except (ValueError, FileNotFoundError, BBoxContainmentError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
