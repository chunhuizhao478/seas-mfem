#!/usr/bin/env python3
"""build_stress_safs.py — Phase 5 schema-v1 HDF5 sidecar writer.

Produces a single canonical ``stress_safs.h5`` schema-v1 sidecar (analog
of ``velocity/results/velocity_safs.h5``) that the C++ runtime will load
via the existing ``DataField3D`` class without any C++-side changes.

The bulk Cauchy field σ_seas is **compression-positive SEAS convention**,
in **Pa**, evaluated on a uniform UTM 11 N rectilinear grid that strictly
encloses the union bbox of the supplied bulk meshes by at least
``pad_m`` on every face.

Conventions
-----------
- σ_seas is built by Phase 3's ``bulk_stress_tensor_field`` which already
  applies the single source-site sign flip (H&Z compression-negative →
  SEAS compression-positive); this module just unit-converts MPa → Pa
  (no further sign manipulation per R-501 / R-502).
- Six fields are written: ``sigma_xx``, ``sigma_yy``, ``sigma_zz``,
  ``sigma_xy``, ``sigma_yz``, ``sigma_xz``.
- Sidecar attrs: ``schema_version``, ``crs``, ``units``, ``z_positive``
  (canonical, filled by ``data_projection.sidecar.write_sidecar``), plus
  ``source``, ``source_crs``, ``mesh_tag``.

Usage
-----
::

    python build_stress_safs.py \\
        --meshes ../../meshing/results/vtu/safs_fault_box_nwcut_500m_bulk.vtu \\
                 ../../meshing/results/vtu/safs_fault_box_nwcut_1000m_bulk.vtu \\
                 ../../meshing/results/vtu/safs_fault_box_nwcut_2000m_bulk.vtu \\
        --out    ../results/stress_safs.h5 \\
        --SHmax 113.0 --Shmin 49.0 --Sv 45.0 --SHmax-az 23.0

Tests
-----
::

    cd project_7.0_alternative/stress/code && pytest -q test_build_stress_safs.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import meshio
import numpy as np


# Sibling import of `data_projection.sidecar`, mirroring
# `data_projection/build_velocity_cvmh.py`.  We insert the
# `data_projection` directory on sys.path so the package's
# `sidecar.py` is importable directly under its module name (Phase 5
# plan §1).
_DATA_PROJ_DIR = (
    Path(__file__).resolve().parent / "data_projection"
)
if str(_DATA_PROJ_DIR) not in sys.path:
    sys.path.insert(0, str(_DATA_PROJ_DIR))
from sidecar import (                                       # noqa: E402
    CANONICAL_CRS,
    CANONICAL_UNITS,
    CANONICAL_Z_POSITIVE,
    SCHEMA_VERSION,
    write_sidecar,
)

# Sibling import of Phase 3's `bulk_stress_tensor_field`.  This module
# lives next to ``project_to_fault_stress.py``; the sys.path insert is
# defensive against pytest discovery from a parent directory.
_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))
from project_to_fault_stress import bulk_stress_tensor_field  # noqa: E402


# ----------------------------------------------------------------------
# Module-level constants — Phase 5 §2
# ----------------------------------------------------------------------

DEFAULT_GRID_DX_M: float = 1000.0
DEFAULT_GRID_PAD_M: float = 2000.0
DEFAULT_BOUNDS_FACTOR: float = 1.2

DEFAULT_OUT: Path = (
    Path(__file__).resolve().parent.parent
    / "results"
    / "stress_safs.h5"
)

# Symmetry tolerance for σ_seas on the (Nz, 3, 3) stack, in Pa.
# 1e-3 Pa = 1e-9 MPa, matching Phase 3's `_SIGMA_SYMMETRY_TOL_MPA`
# after the MPa → Pa unit conversion (plan §4 line 1432–1437 after
# R-801 reconciliation).
_SIGMA_SYMMETRY_TOL_PA: float = 1.0e-3

# Six canonical field names, in the schema-v1 order documented in
# PLAN_onfaultstress.md Phase 5 §4 and §"Required by".  The order
# is preserved on disk because Python dicts are insertion-ordered.
_FIELD_NAMES: tuple[str, ...] = (
    "sigma_xx", "sigma_yy", "sigma_zz",
    "sigma_xy", "sigma_yz", "sigma_xz",
)


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def build_uniform_utm_grid(
    mesh_bbox: tuple[float, float, float, float, float, float],
    *,
    dx_m: float,
    pad_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(x, y, z)`` 1-D arrays for a uniform UTM grid that
    strictly encloses ``mesh_bbox`` by at least ``pad_m`` on every face.

    Parameters
    ----------
    mesh_bbox : (xmin, xmax, ymin, ymax, zmin, zmax)
        Axis-aligned bounding box of the SAFS mesh union, in metres.
    dx_m : float
        Uniform grid spacing on x, y, and z (all axes use the same dx).
    pad_m : float
        Padding distance applied to every face of the bbox before
        rounding outward to ``dx_m`` multiples.  The returned grid
        satisfies ``x[0] <= xmin - pad_m`` and ``x[-1] >= xmax + pad_m``
        (analogously for y, z).

    Returns
    -------
    x, y, z : 1-D ``np.ndarray`` of float64
        Strictly monotone-increasing axes in canonical UTM Zone 11 N
        (EPSG:32611, z = elevation positive — matches
        ``CANONICAL_Z_POSITIVE``).

    Raises
    ------
    ValueError
        If the bbox is malformed (``xmin >= xmax`` etc.) or
        ``dx_m`` / ``pad_m`` are non-positive.
    """
    if dx_m <= 0:
        raise ValueError(f"dx_m must be > 0; got {dx_m}")
    if pad_m < 0:
        raise ValueError(f"pad_m must be >= 0; got {pad_m}")
    xmin, xmax, ymin, ymax, zmin, zmax = (
        float(v) for v in mesh_bbox
    )
    if not (xmin < xmax and ymin < ymax and zmin < zmax):
        raise ValueError(
            f"mesh_bbox must have min < max on every axis; got "
            f"x=[{xmin}, {xmax}], y=[{ymin}, {ymax}], "
            f"z=[{zmin}, {zmax}]"
        )

    def _floor_to(v: float) -> float:
        return float(np.floor(v / dx_m) * dx_m)

    def _ceil_to(v: float) -> float:
        return float(np.ceil(v / dx_m) * dx_m)

    x_lo = _floor_to(xmin - pad_m)
    x_hi = _ceil_to(xmax + pad_m)
    y_lo = _floor_to(ymin - pad_m)
    y_hi = _ceil_to(ymax + pad_m)
    z_lo = _floor_to(zmin - pad_m)
    z_hi = _ceil_to(zmax + pad_m)

    # `+ 0.5 * dx_m` makes np.arange include the upper endpoint
    # whenever the rounded extent is an exact multiple of dx_m.
    x = np.arange(x_lo, x_hi + 0.5 * dx_m, dx_m, dtype=np.float64)
    y = np.arange(y_lo, y_hi + 0.5 * dx_m, dx_m, dtype=np.float64)
    z = np.arange(z_lo, z_hi + 0.5 * dx_m, dx_m, dtype=np.float64)
    return x, y, z


def evaluate_stress_field_on_grid(
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    *,
    SHmax: float, Shmin: float, Sv: float,
    SHmax_az_deg: float,
    depth_model: str = "constant",
    SHmax_grad: float = 0.0,
    Shmin_grad: float = 0.0,
    Sv_grad: float = 0.0,
) -> dict[str, np.ndarray]:
    """Evaluate σ_seas on a regular ``(Nx, Ny, Nz)`` UTM grid and
    return six independent ``(Nx, Ny, Nz)`` component arrays in **Pa**,
    compression-POSITIVE SEAS convention.

    The bulk-path sign flip from H&Z to SEAS is **already applied** by
    Phase 3's ``bulk_stress_tensor_field``; this function applies a
    pure MPa → Pa unit conversion (no further sign manipulation —
    R-501 / R-502 contract).

    For the homogeneous σ⁰ case the field is constant in x and y and
    depends only on z, so the vectorised implementation builds a
    ``(Nz, 3, 3)`` σ_seas stack once and broadcasts.

    Parameters
    ----------
    x, y, z : 1-D arrays of float
        Grid axes in metres.  ``z`` is elevation (z = 0 at the free
        surface, z < 0 in the subsurface) — matches
        ``CANONICAL_Z_POSITIVE``.
    SHmax, Shmin, Sv : float
        Principal-stress magnitudes at z = 0 (POSITIVE values, MPa).
    SHmax_az_deg : float
        Geological azimuth of SHmax (cw from north, degrees).
    depth_model : {"constant", "lithostatic_sv"}
    SHmax_grad, Shmin_grad, Sv_grad : float, default 0.0
        Linear depth gradients (MPa/m), forwarded to
        ``bulk_stress_tensor_field``.

    Returns
    -------
    fields : dict[str, np.ndarray]
        Six entries (``sigma_xx`` … ``sigma_xz``), each
        ``(Nx, Ny, Nz)`` float64 in Pa.  Ordering follows
        ``_FIELD_NAMES`` (schema-v1 canonical order).

    Raises
    ------
    ValueError
        On non-1-D inputs, NaN cells, or σ_seas symmetry violation.
    """
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    z = np.ascontiguousarray(z, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or z.ndim != 1:
        raise ValueError(
            f"x/y/z must be 1-D; got shapes "
            f"{x.shape}/{y.shape}/{z.shape}"
        )
    if x.size == 0 or y.size == 0 or z.size == 0:
        raise ValueError(
            f"x/y/z must be non-empty; got sizes "
            f"{x.size}/{y.size}/{z.size}"
        )

    # Build σ_seas(z) in MPa, (Nz, 3, 3); compression POSITIVE SEAS.
    sigma_seas_MPa = bulk_stress_tensor_field(
        z,
        SHmax_az_deg=SHmax_az_deg,
        depth_model=depth_model,
        SHmax_top=SHmax,
        Shmin_top=Shmin,
        Sv_top=Sv,
        SHmax_grad=SHmax_grad,
        Shmin_grad=Shmin_grad,
        Sv_grad=Sv_grad,
    )

    # Pure unit conversion — no sign flip (R-501 / R-502).
    sigma_seas_Pa = sigma_seas_MPa * 1.0e6

    # NaN guard (defence in depth — the writer also re-checks).
    if np.isnan(sigma_seas_Pa).any():
        raise ValueError(
            "evaluate_stress_field_on_grid: σ_seas contains NaN; "
            "check the z range and depth gradients for non-finite "
            "inputs"
        )

    # Symmetry assertion (single-source-of-truth check; the H&Z
    # rotation should produce symmetric tensors).
    asym = float(
        np.max(np.abs(sigma_seas_Pa - np.swapaxes(sigma_seas_Pa, -1, -2)))
    )
    if asym > _SIGMA_SYMMETRY_TOL_PA:
        raise ValueError(
            f"evaluate_stress_field_on_grid: σ_seas is not symmetric; "
            f"max off-diagonal asymmetry {asym:.3e} Pa > tol "
            f"{_SIGMA_SYMMETRY_TOL_PA:.3e} Pa"
        )

    Nx, Ny, Nz = x.size, y.size, z.size

    def _broadcast(c1: np.ndarray) -> np.ndarray:
        """Broadcast a (Nz,) component along x and y to (Nx, Ny, Nz)
        and return a contiguous owned array."""
        out = np.broadcast_to(c1[None, None, :], (Nx, Ny, Nz))
        return np.ascontiguousarray(out, dtype=np.float64)

    return {
        "sigma_xx": _broadcast(sigma_seas_Pa[:, 0, 0]),
        "sigma_yy": _broadcast(sigma_seas_Pa[:, 1, 1]),
        "sigma_zz": _broadcast(sigma_seas_Pa[:, 2, 2]),
        "sigma_xy": _broadcast(sigma_seas_Pa[:, 0, 1]),
        "sigma_yz": _broadcast(sigma_seas_Pa[:, 1, 2]),
        "sigma_xz": _broadcast(sigma_seas_Pa[:, 0, 2]),
    }


def derive_field_bounds(
    fields: dict[str, np.ndarray],
    *,
    safety_factor: float = DEFAULT_BOUNDS_FACTOR,
) -> dict[str, tuple[float, float, str]]:
    """Return ``{name: (min_value, max_value, "Pa")}`` for every field,
    sized to strictly enclose the actual data with a safety margin so
    that schema-v1's guard G-2 (every cell in ``[min, max]``) passes.

    Rules (Phase 5 §5):
        - If ``min(field) < 0``: ``min_value = min * safety_factor``
          (more negative → bound below data).
        - Else:                  ``min_value = min / safety_factor``
          (closer to zero from above → bound below data).
        - Analogous for max with ``max > 0`` / ``max <= 0``.

    The degenerate identically-zero field would collapse the bounds to
    ``[0, 0]``; we widen by ±1 Pa in that case so ``min < max`` holds
    (the sidecar writer enforces strict inequality).
    """
    if safety_factor <= 1.0:
        raise ValueError(
            f"safety_factor must be > 1.0; got {safety_factor}"
        )
    out: dict[str, tuple[float, float, str]] = {}
    for name, arr in fields.items():
        a_min = float(np.min(arr))
        a_max = float(np.max(arr))
        v_min = (a_min * safety_factor
                 if a_min < 0.0 else a_min / safety_factor)
        v_max = (a_max * safety_factor
                 if a_max > 0.0 else a_max / safety_factor)
        if not (v_min < v_max):
            # Degenerate (e.g. identically-zero field).  Widen by ±1
            # Pa so the writer's `v_min < v_max` guard passes and the
            # actual data is enclosed.
            v_min = a_min - 1.0
            v_max = a_max + 1.0
        out[name] = (v_min, v_max, "Pa")
    return out


# ----------------------------------------------------------------------
# Internal helpers — bbox extraction / union
# ----------------------------------------------------------------------


def _mesh_bbox(
    path: Path,
) -> tuple[float, float, float, float, float, float]:
    """Return the axis-aligned bbox ``(xmin, xmax, ymin, ymax, zmin,
    zmax)`` of a bulk VTU's point cloud."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"mesh not found: {p.resolve()}")
    m = meshio.read(str(p))
    pts = np.asarray(m.points, dtype=np.float64)
    if pts.size == 0:
        raise ValueError(f"mesh {p} has no points")
    return (
        float(pts[:, 0].min()), float(pts[:, 0].max()),
        float(pts[:, 1].min()), float(pts[:, 1].max()),
        float(pts[:, 2].min()), float(pts[:, 2].max()),
    )


def _union_bbox(
    bboxes: list[tuple[float, float, float, float, float, float]],
) -> tuple[float, float, float, float, float, float]:
    """Axis-aligned union of a list of bboxes."""
    if not bboxes:
        raise ValueError("_union_bbox: empty bbox list")
    arr = np.asarray(bboxes, dtype=np.float64)
    return (
        float(arr[:, 0].min()), float(arr[:, 1].max()),
        float(arr[:, 2].min()), float(arr[:, 3].max()),
        float(arr[:, 4].min()), float(arr[:, 5].max()),
    )


# ----------------------------------------------------------------------
# Top-level driver — Phase 5 §6
# ----------------------------------------------------------------------


def build_stress_safs(
    mesh_paths: list[Path],
    out_path: Path,
    *,
    SHmax: float, Shmin: float, Sv: float,
    SHmax_az_deg: float,
    depth_model: str = "constant",
    SHmax_grad: float = 0.0,
    Shmin_grad: float = 0.0,
    Sv_grad: float = 0.0,
    dx_m: float = DEFAULT_GRID_DX_M,
    pad_m: float = DEFAULT_GRID_PAD_M,
    mesh_tag: str = "safs_fault_box_nwcut",
    safety_factor: float = DEFAULT_BOUNDS_FACTOR,
) -> None:
    """Phase 5 driver: read bulk meshes, evaluate σ_seas on a uniform
    UTM grid, and write a schema-v1 HDF5 sidecar at ``out_path``.

    Steps (Phase 5 §6):
        1. Read every mesh and compute the axis-aligned union bbox.
        2. Build the uniform UTM grid via ``build_uniform_utm_grid``.
        3. Evaluate σ_seas on the grid via
           ``evaluate_stress_field_on_grid`` (returns Pa,
           compression-positive SEAS).
        4. Derive per-field bounds via ``derive_field_bounds``.
        5. Write the schema-v1 sidecar via
           ``data_projection.sidecar.write_sidecar``.
    """
    if not mesh_paths:
        raise ValueError("build_stress_safs: mesh_paths is empty")
    if depth_model == "lithostatic_sv" and all(
        g == 0.0 for g in (SHmax_grad, Shmin_grad, Sv_grad)
    ):
        print(
            "warning: depth_model='lithostatic_sv' with all gradients "
            "= 0 is equivalent to 'constant'",
            file=sys.stderr,
        )

    # 1. Union bbox across all input meshes.
    bboxes = [_mesh_bbox(Path(p)) for p in mesh_paths]
    union = _union_bbox(bboxes)

    # 2. Uniform UTM grid.
    x, y, z = build_uniform_utm_grid(
        union, dx_m=dx_m, pad_m=pad_m,
    )

    # 3. σ_seas on the grid (Pa, compression-positive SEAS).
    fields = evaluate_stress_field_on_grid(
        x, y, z,
        SHmax=SHmax, Shmin=Shmin, Sv=Sv,
        SHmax_az_deg=SHmax_az_deg,
        depth_model=depth_model,
        SHmax_grad=SHmax_grad,
        Shmin_grad=Shmin_grad,
        Sv_grad=Sv_grad,
    )

    # 4. Per-field bounds.
    field_bounds = derive_field_bounds(
        fields, safety_factor=safety_factor,
    )

    # 5. Sidecar attrs (canonical fields are filled by write_sidecar;
    # we still pass them so the writer's reject-non-canonical check
    # confirms the writer-side conventions match this module's).
    attrs = {
        "schema_version": SCHEMA_VERSION,
        "crs": CANONICAL_CRS,
        "units": CANONICAL_UNITS,
        "z_positive": CANONICAL_Z_POSITIVE,
        "source": ",".join(Path(p).name for p in mesh_paths),
        "source_crs": "EPSG:32611",
        "mesh_tag": mesh_tag,
    }

    write_sidecar(
        out_path=Path(out_path),
        x=x, y=y, z=z,
        fields=fields,
        attrs=attrs,
        field_bounds=field_bounds,
    )


# ----------------------------------------------------------------------
# CLI — Phase 5 §7
# ----------------------------------------------------------------------


def _parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--meshes", type=Path, nargs="+", required=True,
        help="one or more bulk VTU paths whose union bbox sets the grid extent",
    )
    ap.add_argument(
        "--out", type=Path, default=DEFAULT_OUT,
        help=f"output sidecar path (default {DEFAULT_OUT})",
    )
    ap.add_argument(
        "--SHmax", type=float, required=True,
        help="SHmax magnitude at z=0 (POSITIVE, MPa)",
    )
    ap.add_argument(
        "--Shmin", type=float, required=True,
        help="Shmin magnitude at z=0 (POSITIVE, MPa)",
    )
    ap.add_argument(
        "--Sv", type=float, required=True,
        help="Sv magnitude at z=0 (POSITIVE, MPa)",
    )
    ap.add_argument(
        "--SHmax-az", type=float, required=True, dest="SHmax_az",
        help="SHmax azimuth (cw from north, degrees)",
    )
    ap.add_argument(
        "--depth-model", type=str, default="constant",
        choices=("constant", "lithostatic_sv"),
        dest="depth_model",
    )
    ap.add_argument(
        "--SHmax-grad", type=float, default=0.0, dest="SHmax_grad",
        help="SHmax depth gradient (MPa/m); active with --depth-model lithostatic_sv",
    )
    ap.add_argument(
        "--Shmin-grad", type=float, default=0.0, dest="Shmin_grad",
        help="Shmin depth gradient (MPa/m)",
    )
    ap.add_argument(
        "--Sv-grad", type=float, default=0.0, dest="Sv_grad",
        help="Sv depth gradient (MPa/m)",
    )
    ap.add_argument(
        "--dx", type=float, default=DEFAULT_GRID_DX_M,
        help=f"uniform grid spacing on x, y, z (m); default {DEFAULT_GRID_DX_M}",
    )
    ap.add_argument(
        "--pad", type=float, default=DEFAULT_GRID_PAD_M,
        help=f"grid padding outside the mesh bbox (m); default {DEFAULT_GRID_PAD_M}",
    )
    ap.add_argument(
        "--mesh-tag", type=str, default="safs_fault_box_nwcut",
        dest="mesh_tag",
        help="mesh_tag attribute in the sidecar root",
    )
    return ap.parse_args(argv)


def main(argv: Optional[list] = None) -> int:
    args = _parse_args(argv)
    try:
        build_stress_safs(
            mesh_paths=list(args.meshes),
            out_path=args.out,
            SHmax=args.SHmax, Shmin=args.Shmin, Sv=args.Sv,
            SHmax_az_deg=args.SHmax_az,
            depth_model=args.depth_model,
            SHmax_grad=args.SHmax_grad,
            Shmin_grad=args.Shmin_grad,
            Sv_grad=args.Sv_grad,
            dx_m=args.dx, pad_m=args.pad,
            mesh_tag=args.mesh_tag,
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
