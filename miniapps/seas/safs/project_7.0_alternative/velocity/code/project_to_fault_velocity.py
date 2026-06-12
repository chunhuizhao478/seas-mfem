#!/usr/bin/env python3
"""
project_to_fault_velocity.py — Sample the CVM velocity-model sidecar
(Vp, Vs, density) onto the SAFS fault surface (and optionally the bulk
volume) and write ParaView VTU artefacts.

Companion of `stress/code/project_to_fault_stress.py`: same output
layout (`<out_dir>/<mesh_base>_fault_velocity.vtu` + summary JSON,
optional `<mesh_base>_bulk_velocity.vtu`), same point-data
(continuous) + cell-data (`*_cell`) split. One deliberate numerical
difference: the stress pipeline node-averages cell values because its
field depends on the per-cell fault basis; the velocity model is a
volumetric field defined everywhere, so point data here is sampled
DIRECTLY at the fault vertices (exact trilinear) and cell data at the
triangle centroids — no cell→node averaging involved.

Conventions
-----------
- Coordinate frame: UTM Zone 11N metres; (x = east, y = north, z = up,
  z = 0 at the free surface, negative below). Matches both the mesh
  and the sidecar (`crs = EPSG:32611`, `z_positive = elevation`).
- Sidecar schema (`data_projection_v1`): rectilinear grid under
  `grid/{x,y,z}` (strictly increasing, z non-uniform allowed), fields
  under `fields/{Vp,Vs,density}` indexed `(ix, iy, iz)`.
- Trilinear sampling clamps to the grid bounding box (same behaviour
  as `plot_comparison_with_dg0.py:trilinear`); points outside the box
  are counted and reported in the summary JSON.

Usage
-----
    python project_to_fault_velocity.py <INPUT_MESH.msh | *_fault.vtu>
        [--sidecar PATH] [--write-bulk]
        [--fault-name fault] [--bulk-name rock]
        [--out-dir PATH] [--mesh-base NAME]

Default sidecar:
    ../results/multiscale_statewise_cvm/velocity_safs.h5
Default out dir:
    ../results/multiscale_statewise_cvm/<mesh_base>/

ParaView verification recipe:
1. Open the produced `<base>_fault_velocity.vtu`.
2. "Surface" representation; colour by `Vs_m_per_s` — the low-velocity
   fault-zone / basin structure should be visible near the surface.
3. Compare against `density_kg_per_m3` and the derived `mu_GPa`.

Library use:

    from project_to_fault_velocity import (
        load_velocity_sidecar, trilinear_sample,
        sample_fields_at_points, count_out_of_bounds,
        write_fault_velocity_vtu, write_bulk_velocity_vtu,
        write_summary_json,
    )

Tests
-----
    cd project_7.0_alternative/velocity/code && \
        pytest -q test_project_to_fault_velocity.py
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import h5py
import meshio
import numpy as np


# ----------------------------------------------------------------------
# Sibling import of the stress-projection geometry library
# ----------------------------------------------------------------------

_VELOCITY_CODE_DIR = Path(__file__).resolve().parent
_STRESS_CODE_DIR = _VELOCITY_CODE_DIR.parents[1] / "stress" / "code"
for _d in (str(_VELOCITY_CODE_DIR), str(_STRESS_CODE_DIR)):
    if _d not in sys.path:
        sys.path.insert(0, _d)
from project_to_fault_stress import (  # noqa: E402
    extract_cells_by_physical,
    tet_geometry,
    triangle_geometry,
)


# ----------------------------------------------------------------------
# Module-level constants
# ----------------------------------------------------------------------

DEFAULT_FAULT_NAME: str = "fault"
DEFAULT_BULK_NAME: str = "rock"

# Sidecar dataset names → output VTU field names (units explicit).
SIDECAR_FIELDS: dict = {
    "Vp": "Vp_m_per_s",
    "Vs": "Vs_m_per_s",
    "density": "density_kg_per_m3",
}

# Derived field names (computed from the sampled primaries).
DERIVED_FIELDS: tuple = ("mu_GPa", "vp_vs_ratio")

DEFAULT_SIDECAR: Path = (
    _VELOCITY_CODE_DIR.parent
    / "results" / "multiscale_statewise_cvm" / "velocity_safs.h5"
)

DEFAULT_OUT_ROOT: Path = (
    _VELOCITY_CODE_DIR.parent / "results" / "multiscale_statewise_cvm"
)

_CONVENTION_STRING: str = (
    "velocity sidecar data_projection_v1; frame (east, north, up) / "
    "UTM Zone 11N metres, z=0 free surface; Vp/Vs in m/s, density in "
    "kg/m^3; mu_GPa = density * Vs^2 * 1e-9; trilinear sampling "
    "clamped to the grid bounding box"
)


# ----------------------------------------------------------------------
# Sidecar I/O
# ----------------------------------------------------------------------


@dataclass
class VelocityGrid:
    """Rectilinear velocity-model grid loaded from the sidecar HDF5."""

    x: np.ndarray                    # (nx,) strictly increasing, m
    y: np.ndarray                    # (ny,) strictly increasing, m
    z: np.ndarray                    # (nz,) strictly increasing, m (up)
    fields: dict = field(default_factory=dict)  # name -> (nx, ny, nz)
    attrs: dict = field(default_factory=dict)   # root HDF5 attrs


def _check_axis(axis: np.ndarray, name: str) -> np.ndarray:
    """Validate a grid axis: 1-D, length ≥ 2, finite, strictly
    increasing. Returns the float64 view."""
    axis = np.asarray(axis, dtype=np.float64)
    if axis.ndim != 1:
        raise ValueError(
            f"grid/{name} must be 1-D; got shape {axis.shape}"
        )
    if axis.size < 2:
        raise ValueError(
            f"grid/{name} needs at least 2 samples for interpolation; "
            f"got {axis.size}"
        )
    if not np.all(np.isfinite(axis)):
        raise ValueError(f"grid/{name} contains non-finite values")
    if not np.all(np.diff(axis) > 0):
        raise ValueError(
            f"grid/{name} must be strictly increasing"
        )
    return axis


def load_velocity_sidecar(
    path: Path,
    field_names: tuple = tuple(SIDECAR_FIELDS),
) -> VelocityGrid:
    """Load the `data_projection_v1` velocity sidecar.

    Parameters
    ----------
    path : Path
        HDF5 sidecar (e.g. `velocity_safs.h5`).
    field_names : tuple of str
        Dataset names under `fields/` to load.

    Returns
    -------
    VelocityGrid

    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    KeyError
        If a required group / dataset is missing.
    ValueError
        If an axis is invalid or a field's shape does not match the
        grid.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"velocity sidecar not found: {p.resolve()}")

    with h5py.File(str(p), "r") as h5:
        if "grid" not in h5 or "fields" not in h5:
            raise KeyError(
                f"sidecar {p.resolve()} missing 'grid' and/or 'fields' "
                f"groups; top-level keys: {sorted(h5.keys())!r}"
            )
        x = _check_axis(h5["grid/x"][:], "x")
        y = _check_axis(h5["grid/y"][:], "y")
        z = _check_axis(h5["grid/z"][:], "z")
        expected = (x.size, y.size, z.size)

        fields_out: dict = {}
        for name in field_names:
            if name not in h5["fields"]:
                raise KeyError(
                    f"sidecar {p.resolve()} missing fields/{name}; "
                    f"available: {sorted(h5['fields'].keys())!r}"
                )
            arr = np.asarray(h5["fields"][name][:], dtype=np.float64)
            if arr.shape != expected:
                raise ValueError(
                    f"fields/{name} shape {arr.shape} != grid shape "
                    f"{expected} (nx, ny, nz)"
                )
            n_bad = int((~np.isfinite(arr)).sum())
            if n_bad:
                print(
                    f"warning: fields/{name} contains {n_bad} "
                    f"non-finite samples; they will propagate as NaN "
                    f"through the trilinear stencil",
                    file=sys.stderr,
                )
            fields_out[name] = arr

        attrs = {k: h5.attrs[k] for k in h5.attrs}

    return VelocityGrid(x=x, y=y, z=z, fields=fields_out, attrs=attrs)


# ----------------------------------------------------------------------
# Trilinear sampling
# ----------------------------------------------------------------------


def trilinear_sample(
    gx: np.ndarray,
    gy: np.ndarray,
    gz: np.ndarray,
    F: np.ndarray,
    P: np.ndarray,
) -> np.ndarray:
    """Trilinear interpolation of a rectilinear field at points.

    Same algorithm as the established
    `plot_comparison_with_dg0.py:trilinear` helper (clamped
    searchsorted cell location; handles non-uniform axes), kept local
    because the plot module pulls in matplotlib / pyproj at import
    time. Points outside the grid bounding box are clamped to the
    boundary value (count them via `count_out_of_bounds` if needed).
    NaN coordinate rows yield NaN output.

    Parameters
    ----------
    gx, gy, gz : 1-D float arrays
        Strictly increasing grid axes.
    F : (nx, ny, nz) float array
        Field values indexed `(ix, iy, iz)`.
    P : (N, 3) float array
        Query points `(x, y, z)`.

    Returns
    -------
    (N,) float64 values.
    """
    gx = np.asarray(gx, dtype=np.float64)
    gy = np.asarray(gy, dtype=np.float64)
    gz = np.asarray(gz, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)
    P = np.asarray(P, dtype=np.float64)
    if P.ndim != 2 or P.shape[1] != 3:
        raise ValueError(f"P must be (N, 3); got shape {P.shape}")
    if F.shape != (gx.size, gy.size, gz.size):
        raise ValueError(
            f"F shape {F.shape} != (nx, ny, nz) = "
            f"({gx.size}, {gy.size}, {gz.size})"
        )

    nx, ny, nz = len(gx), len(gy), len(gz)
    # NaN coordinates: searchsorted(NaN) returns n → clip keeps a
    # valid index, and the NaN fraction below propagates NaN to the
    # output (np.clip(NaN, 0, 1) is NaN), which is the wanted
    # behaviour for degenerate-triangle centroids.
    px = np.clip(np.searchsorted(gx, P[:, 0]) - 1, 0, nx - 2)
    py = np.clip(np.searchsorted(gy, P[:, 1]) - 1, 0, ny - 2)
    pz = np.clip(np.searchsorted(gz, P[:, 2]) - 1, 0, nz - 2)
    fx = np.clip((P[:, 0] - gx[px]) / (gx[px + 1] - gx[px]), 0, 1)
    fy = np.clip((P[:, 1] - gy[py]) / (gy[py + 1] - gy[py]), 0, 1)
    fz = np.clip((P[:, 2] - gz[pz]) / (gz[pz + 1] - gz[pz]), 0, 1)
    c00 = F[px, py, pz] * (1 - fx) + F[px + 1, py, pz] * fx
    c10 = F[px, py + 1, pz] * (1 - fx) + F[px + 1, py + 1, pz] * fx
    c01 = F[px, py, pz + 1] * (1 - fx) + F[px + 1, py, pz + 1] * fx
    c11 = (F[px, py + 1, pz + 1] * (1 - fx)
           + F[px + 1, py + 1, pz + 1] * fx)
    return ((c00 * (1 - fy) + c10 * fy) * (1 - fz)
            + (c01 * (1 - fy) + c11 * fy) * fz)


def count_out_of_bounds(
    gx: np.ndarray,
    gy: np.ndarray,
    gz: np.ndarray,
    P: np.ndarray,
) -> int:
    """Count finite points strictly outside the grid bounding box
    (these are CLAMPED by `trilinear_sample`, not rejected)."""
    P = np.asarray(P, dtype=np.float64)
    if P.ndim != 2 or P.shape[1] != 3:
        raise ValueError(f"P must be (N, 3); got shape {P.shape}")
    finite = np.isfinite(P).all(axis=1)
    out = (
        (P[:, 0] < gx[0]) | (P[:, 0] > gx[-1])
        | (P[:, 1] < gy[0]) | (P[:, 1] > gy[-1])
        | (P[:, 2] < gz[0]) | (P[:, 2] > gz[-1])
    )
    return int((out & finite).sum())


def sample_fields_at_points(
    grid: VelocityGrid,
    P: np.ndarray,
) -> dict:
    """Sample every sidecar field at the given points and derive the
    secondary quantities.

    Derived AFTER interpolation (interpolate-then-derive), matching
    how a solver consuming (Vp, Vs, rho) point values would compute
    moduli:

    - ``mu_GPa``      = density · Vs² · 1e-9
    - ``vp_vs_ratio`` = Vp / Vs  (NaN where Vs ≤ 0)

    Returns
    -------
    dict
        Output-field name → (N,) float64 array. Keys are the values
        of ``SIDECAR_FIELDS`` plus ``DERIVED_FIELDS``.
    """
    sampled: dict = {}
    for ds_name, out_name in SIDECAR_FIELDS.items():
        sampled[out_name] = trilinear_sample(
            grid.x, grid.y, grid.z, grid.fields[ds_name], P,
        )

    vp = sampled["Vp_m_per_s"]
    vs = sampled["Vs_m_per_s"]
    rho = sampled["density_kg_per_m3"]

    sampled["mu_GPa"] = rho * vs ** 2 * 1.0e-9
    ratio = np.full_like(vp, np.nan)
    ok = np.isfinite(vs) & (vs > 0.0)
    ratio[ok] = vp[ok] / vs[ok]
    sampled["vp_vs_ratio"] = ratio
    return sampled


# ----------------------------------------------------------------------
# Mesh loading
# ----------------------------------------------------------------------


def load_input_mesh(path: Path) -> meshio.Mesh:
    """Load the input mesh: a Gmsh v2.2 `.msh` (fault triangles +
    bulk tets + physical names) or a `*_fault.vtu` written by
    `meshing/code/msh_to_vtu.py`.

    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    ValueError
        If the file contains no triangle cells.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"input mesh not found: {p.resolve()}")
    mesh = meshio.read(str(p))
    if not any(cb.type == "triangle" for cb in mesh.cells):
        raise ValueError(
            f"input mesh has no triangle cells: {p.resolve()} "
            f"(cell types present: {[cb.type for cb in mesh.cells]})"
        )
    return mesh


def _physical_tag(mesh: meshio.Mesh, phys_name: str) -> Optional[int]:
    """Look up a physical-group tag by name; None if field_data is
    absent (e.g. a VTU that lost it on write)."""
    field_data = getattr(mesh, "field_data", None) or {}
    if phys_name in field_data:
        return int(field_data[phys_name][0])
    return None


# ----------------------------------------------------------------------
# VTU writers
# ----------------------------------------------------------------------


def write_fault_velocity_vtu(
    mesh: meshio.Mesh,
    tri_conn: np.ndarray,
    point_values: dict,
    cell_values: dict,
    out_path: Path,
    fault_phys_tag: Optional[int] = None,
    convention_str: str = _CONVENTION_STRING,
) -> None:
    """Write the fault VTU: point-data primary fields + `*_cell`
    cell-data channel, mirroring `write_fault_vtu` in the stress
    pipeline (full input point cloud carried through; orphan vertices
    — those not referenced by a fault triangle — hold NaN).

    Parameters
    ----------
    mesh : meshio.Mesh
        Input mesh (point cloud reused for the output).
    tri_conn : (N_tri, 3) ndarray
        Fault triangle connectivity into ``mesh.points``.
    point_values : dict
        Output-field name → (N_pt,) array (NaN on orphan vertices).
    cell_values : dict
        Output-field name → (N_tri,) array (centroid samples).
    out_path : Path
        Destination VTU. Parent directory created if missing.
    fault_phys_tag : int or None
        If given, a constant ``gmsh:physical`` cell array is emitted
        (every output cell belongs to the fault group by
        construction).
    convention_str : str
        Metadata string attached as ``mesh.info``.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    points = np.asarray(mesh.points, dtype=np.float64)
    tri_conn = np.asarray(tri_conn, dtype=np.int64)
    N_pt = points.shape[0]
    N_tri = tri_conn.shape[0]

    point_data: dict = {}
    for name, vals in point_values.items():
        vals = np.asarray(vals, dtype=np.float64)
        if vals.shape != (N_pt,):
            raise ValueError(
                f"point field {name!r} shape {vals.shape} != "
                f"(N_pt={N_pt},)"
            )
        point_data[name] = vals

    cell_data: dict = {}
    for name, vals in cell_values.items():
        vals = np.asarray(vals, dtype=np.float64)
        if vals.shape != (N_tri,):
            raise ValueError(
                f"cell field {name!r} shape {vals.shape} != "
                f"(N_tri={N_tri},)"
            )
        cell_data[f"{name}_cell"] = [vals]
    if fault_phys_tag is not None:
        cell_data["gmsh:physical"] = [
            np.full(N_tri, fault_phys_tag, dtype=np.int64)
        ]

    fault_out = meshio.Mesh(
        points=points,
        cells=[("triangle", tri_conn)],
        point_data=point_data,
        cell_data=cell_data,
    )
    fault_out.info = {
        "velocity_convention": convention_str,
        "frame": "(east, north, up) / UTM Zone 11N",
    }
    meshio.write(str(out_path), fault_out)


def write_bulk_velocity_vtu(
    mesh: meshio.Mesh,
    tet_conn: np.ndarray,
    cell_values: dict,
    out_path: Path,
    bulk_phys_tag: Optional[int] = None,
    convention_str: str = _CONVENTION_STRING,
) -> None:
    """Write the bulk VTU as cell-data only (tet-centroid samples),
    mirroring `write_bulk_vtu` in the stress pipeline."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    points = np.asarray(mesh.points, dtype=np.float64)
    tet_conn = np.asarray(tet_conn, dtype=np.int64)
    N_tet = tet_conn.shape[0]

    cell_data: dict = {}
    for name, vals in cell_values.items():
        vals = np.asarray(vals, dtype=np.float64)
        if vals.shape != (N_tet,):
            raise ValueError(
                f"cell field {name!r} shape {vals.shape} != "
                f"(N_tet={N_tet},)"
            )
        cell_data[f"{name}_cell"] = [vals]
    if bulk_phys_tag is not None:
        cell_data["gmsh:physical"] = [
            np.full(N_tet, bulk_phys_tag, dtype=np.int64)
        ]

    bulk_out = meshio.Mesh(
        points=points,
        cells=[("tetra", tet_conn)],
        cell_data=cell_data,
    )
    bulk_out.info = {
        "velocity_convention": convention_str,
        "frame": "(east, north, up) / UTM Zone 11N",
    }
    meshio.write(str(out_path), bulk_out)


# ----------------------------------------------------------------------
# Summary JSON
# ----------------------------------------------------------------------


def _field_stats(arr: np.ndarray) -> dict:
    """min / median / max / nan_count over the finite entries of a
    1-D field (same schema as the stress pipeline's summary)."""
    arr = np.asarray(arr, dtype=np.float64).ravel()
    finite_mask = np.isfinite(arr)
    nan_count = int(arr.size - finite_mask.sum())
    if finite_mask.any():
        finite_vals = arr[finite_mask]
        return {
            "min": float(np.min(finite_vals)),
            "median": float(np.median(finite_vals)),
            "max": float(np.max(finite_vals)),
            "nan_count": nan_count,
        }
    return {
        "min": float("nan"),
        "median": float("nan"),
        "max": float("nan"),
        "nan_count": nan_count,
    }


def _json_safe(value):
    """Convert numpy scalars / Path / bytes to JSON-serialisable types."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, Path):
        return str(value)
    return value


def write_summary_json(
    out_path: Path,
    *,
    input_mesh_path: Path,
    sidecar_path: Path,
    grid: VelocityGrid,
    params: dict,
    fault_n_cells: int,
    fault_n_degenerate: int,
    fault_n_points: int,
    n_clamped: dict,
    fault_cell_values: dict,
    fault_point_values_referenced: dict,
    bulk_n_cells: Optional[int] = None,
    bulk_cell_values: Optional[dict] = None,
) -> None:
    """Write the per-run JSON summary (schema mirrors the stress
    pipeline's: input echo + per-field min/median/max/nan_count).

    ``fault_point_values_referenced`` must hold ONLY the values at
    vertices referenced by a fault triangle — orphan-vertex NaN
    padding would otherwise dominate ``nan_count``.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "input_mesh": str(input_mesh_path),
        "sidecar": str(sidecar_path),
        "sidecar_attrs": {
            k: _json_safe(v) for k, v in grid.attrs.items()
        },
        "grid": {
            "shape": [int(grid.x.size), int(grid.y.size),
                      int(grid.z.size)],
            "x_range_m": [float(grid.x[0]), float(grid.x[-1])],
            "y_range_m": [float(grid.y[0]), float(grid.y[-1])],
            "z_range_m": [float(grid.z[0]), float(grid.z[-1])],
        },
        "params": {k: _json_safe(v) for k, v in params.items()},
        "fault_n_cells": int(fault_n_cells),
        "fault_n_degenerate_cells": int(fault_n_degenerate),
        "fault_n_points": int(fault_n_points),
        "bulk_n_cells": (
            int(bulk_n_cells) if bulk_n_cells is not None else None
        ),
        "n_clamped_to_grid_bbox": {
            k: int(v) for k, v in n_clamped.items()
        },
        "stats": {
            "fault_cell": {
                name: _field_stats(vals)
                for name, vals in fault_cell_values.items()
            },
            "fault_point": {
                name: _field_stats(vals)
                for name, vals in fault_point_values_referenced.items()
            },
            **(
                {
                    "bulk_cell": {
                        name: _field_stats(vals)
                        for name, vals in bulk_cell_values.items()
                    }
                }
                if bulk_cell_values is not None else {}
            ),
        },
        "convention": _CONVENTION_STRING,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2, allow_nan=True)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------


def _build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Sample the CVM velocity sidecar (Vp, Vs, density) onto "
            "the SAFS fault surface and write VTU artefacts."
        ),
    )
    ap.add_argument(
        "input_mesh", type=Path,
        help="Gmsh v2.2 .msh (fault + bulk + physical names) or a "
             "*_fault.vtu from msh_to_vtu.py",
    )
    ap.add_argument(
        "--sidecar", type=Path, default=DEFAULT_SIDECAR,
        help=f"velocity sidecar HDF5 (default: {DEFAULT_SIDECAR})",
    )
    ap.add_argument(
        "--fault-name", default=DEFAULT_FAULT_NAME,
        help="fault physical-group name (default: fault)",
    )
    ap.add_argument(
        "--bulk-name", default=DEFAULT_BULK_NAME,
        help="bulk physical-group name (default: rock)",
    )
    ap.add_argument(
        "--write-bulk", action="store_true",
        help="also sample tet centroids and write "
             "<base>_bulk_velocity.vtu (requires .msh input or a "
             "sibling *_bulk.vtu)",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=None,
        help="output directory (default: "
             "velocity/results/multiscale_statewise_cvm/<mesh_base>)",
    )
    ap.add_argument(
        "--mesh-base", default=None,
        help="output base name (default: input stem, with a trailing "
             "'_fault' stripped for *_fault.vtu inputs)",
    )
    return ap


def _derive_mesh_base(input_mesh: Path) -> str:
    stem = Path(input_mesh).stem
    if stem.endswith("_fault"):
        stem = stem[: -len("_fault")]
    return stem


def _run(args: argparse.Namespace) -> int:
    input_mesh = Path(args.input_mesh)
    mesh_base = args.mesh_base or _derive_mesh_base(input_mesh)
    out_dir = (
        Path(args.out_dir) if args.out_dir is not None
        else DEFAULT_OUT_ROOT / mesh_base
    )

    print(f"loading sidecar: {args.sidecar}")
    grid = load_velocity_sidecar(args.sidecar)

    print(f"loading mesh: {input_mesh}")
    mesh = load_input_mesh(input_mesh)
    points = np.asarray(mesh.points, dtype=np.float64)

    tri_conn = extract_cells_by_physical(
        mesh, "triangle", args.fault_name
    )
    fault_tag = _physical_tag(mesh, args.fault_name)
    N_tri = tri_conn.shape[0]
    N_pt = points.shape[0]
    print(f"fault: {N_tri} triangles / {N_pt} mesh points")

    centroids, _normals, areas = triangle_geometry(points, tri_conn)
    n_degen = int((areas <= 0.0).sum())
    if n_degen:
        print(
            f"warning: {n_degen} degenerate (zero-area) fault "
            f"triangles; their cell values are NaN",
            file=sys.stderr,
        )

    # Cell data: trilinear at triangle centroids.
    fault_cell_values = sample_fields_at_points(grid, centroids)

    # Point data: trilinear directly at the referenced fault
    # vertices; orphan vertices (full bulk point cloud carried
    # through, msh_to_vtu.py convention) stay NaN.
    referenced = np.unique(tri_conn.ravel())
    ref_pts = points[referenced]
    ref_values = sample_fields_at_points(grid, ref_pts)
    fault_point_values: dict = {}
    for name, vals in ref_values.items():
        full = np.full(N_pt, np.nan, dtype=np.float64)
        full[referenced] = vals
        fault_point_values[name] = full

    n_clamped = {
        "fault_points": count_out_of_bounds(
            grid.x, grid.y, grid.z, ref_pts
        ),
        "fault_cells": count_out_of_bounds(
            grid.x, grid.y, grid.z, centroids
        ),
    }

    fault_vtu = out_dir / f"{mesh_base}_fault_velocity.vtu"
    write_fault_velocity_vtu(
        mesh, tri_conn, fault_point_values, fault_cell_values,
        fault_vtu, fault_phys_tag=fault_tag,
    )
    print(f"wrote {fault_vtu}")

    # Optional bulk channel.
    bulk_n_cells: Optional[int] = None
    bulk_cell_values: Optional[dict] = None
    if args.write_bulk:
        bulk_mesh = mesh
        if not any(cb.type == "tetra" for cb in mesh.cells):
            sibling = input_mesh.with_name(
                f"{mesh_base}_bulk.vtu"
            )
            if not sibling.is_file():
                print(
                    f"error: --write-bulk needs tetra cells in the "
                    f"input mesh or a sibling {sibling.name}; neither "
                    f"found",
                    file=sys.stderr,
                )
                return 1
            print(f"loading bulk mesh: {sibling}")
            bulk_mesh = meshio.read(str(sibling))
        tet_conn = extract_cells_by_physical(
            bulk_mesh, "tetra", args.bulk_name
        )
        bulk_tag = _physical_tag(bulk_mesh, args.bulk_name)
        tet_centroids, _vols = tet_geometry(
            np.asarray(bulk_mesh.points, dtype=np.float64), tet_conn
        )
        bulk_n_cells = tet_conn.shape[0]
        print(f"bulk: {bulk_n_cells} tets")
        bulk_cell_values = sample_fields_at_points(grid, tet_centroids)
        n_clamped["bulk_cells"] = count_out_of_bounds(
            grid.x, grid.y, grid.z, tet_centroids
        )
        bulk_vtu = out_dir / f"{mesh_base}_bulk_velocity.vtu"
        write_bulk_velocity_vtu(
            bulk_mesh, tet_conn, bulk_cell_values, bulk_vtu,
            bulk_phys_tag=bulk_tag,
        )
        print(f"wrote {bulk_vtu}")

    summary_path = out_dir / f"{mesh_base}_summary.json"
    write_summary_json(
        summary_path,
        input_mesh_path=input_mesh,
        sidecar_path=Path(args.sidecar),
        grid=grid,
        params={
            "fault_name": args.fault_name,
            "bulk_name": args.bulk_name,
            "write_bulk": bool(args.write_bulk),
            "mesh_base": mesh_base,
            "out_dir": out_dir,
        },
        fault_n_cells=N_tri,
        fault_n_degenerate=n_degen,
        fault_n_points=int(referenced.size),
        n_clamped=n_clamped,
        fault_cell_values=fault_cell_values,
        fault_point_values_referenced=ref_values,
        bulk_n_cells=bulk_n_cells,
        bulk_cell_values=bulk_cell_values,
    )
    print(f"wrote {summary_path}")
    return 0


def main(argv: Optional[list] = None) -> int:
    args = _build_argparser().parse_args(argv)
    return _run(args)


if __name__ == "__main__":
    sys.exit(main())
