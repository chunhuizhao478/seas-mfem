"""bbox_check.py — mesh-vs-grid bounding box containment guards.

The runtime C++ side performs the same check via
``DataField3D::ContainsBBox``; this module is the writer-side
counterpart, called from ``build_velocity_cvmh.py`` (and by user
scripts that want to verify before producing a sidecar).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


class BBoxContainmentError(ValueError):
    """Raised when the mesh bbox is not fully contained in the grid bbox."""


def mesh_msh_bbox(msh_path: Path) -> dict:
    """Return the bounding box of the points in a gmsh ``.msh`` file
    (in whatever coordinate system the mesh was written in — for SAFS
    work this is UTM 11 N, EPSG:32611).

    Uses ``meshio`` to read the mesh and inspects ``mesh.points``.

    Parameters
    ----------
    msh_path : Path
        Path to a gmsh ``.msh`` file readable by ``meshio.read``.

    Returns
    -------
    dict
        ``{"xmin": float, "xmax": float, "ymin": float, "ymax": float,
        "zmin": float, "zmax": float}``.
    """
    import meshio
    p = Path(msh_path)
    if not p.is_file():
        raise FileNotFoundError(f"mesh file not found: {p}")
    m = meshio.read(str(p))
    P = np.asarray(m.points, dtype=np.float64)
    if P.ndim != 2 or P.shape[1] < 3:
        raise ValueError(
            f"{p}: expected (N, 3) points; got shape {P.shape}")
    return {
        "xmin": float(P[:, 0].min()), "xmax": float(P[:, 0].max()),
        "ymin": float(P[:, 1].min()), "ymax": float(P[:, 1].max()),
        "zmin": float(P[:, 2].min()), "zmax": float(P[:, 2].max()),
    }


def grid_bbox(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """Return the bounding box of a rectilinear grid (its first/last
    samples per axis).
    """
    return {
        "xmin": float(x[0]), "xmax": float(x[-1]),
        "ymin": float(y[0]), "ymax": float(y[-1]),
        "zmin": float(z[0]), "zmax": float(z[-1]),
    }


def grid_contains_bbox(grid_x: np.ndarray, grid_y: np.ndarray,
                       grid_z: np.ndarray, mesh_bbox: dict,
                       eps: float = 0.0) -> bool:
    """Return True iff the rectilinear grid (`grid_x`, `grid_y`,
    `grid_z`) fully contains the supplied mesh bbox (with optional
    tolerance ``eps`` in metres).
    """
    gb = grid_bbox(grid_x, grid_y, grid_z)
    return ((mesh_bbox["xmin"] >= gb["xmin"] - eps) and
            (mesh_bbox["xmax"] <= gb["xmax"] + eps) and
            (mesh_bbox["ymin"] >= gb["ymin"] - eps) and
            (mesh_bbox["ymax"] <= gb["ymax"] + eps) and
            (mesh_bbox["zmin"] >= gb["zmin"] - eps) and
            (mesh_bbox["zmax"] <= gb["zmax"] + eps))


def _format_bbox_table(grid: dict, mesh: dict) -> str:
    return (
        "  data grid bbox (UTM 11 N, m):\n"
        f"      x = [{grid['xmin']:.2f}, {grid['xmax']:.2f}]\n"
        f"      y = [{grid['ymin']:.2f}, {grid['ymax']:.2f}]\n"
        f"      z = [{grid['zmin']:.2f}, {grid['zmax']:.2f}]\n"
        "  mesh bbox (UTM 11 N, m):\n"
        f"      x = [{mesh['xmin']:.2f}, {mesh['xmax']:.2f}]\n"
        f"      y = [{mesh['ymin']:.2f}, {mesh['ymax']:.2f}]\n"
        f"      z = [{mesh['zmin']:.2f}, {mesh['zmax']:.2f}]")


def assert_grid_contains_mesh(x: np.ndarray, y: np.ndarray, z: np.ndarray,
                              mesh_bbox: dict, eps: float = 0.0) -> None:
    """Raise ``BBoxContainmentError`` with a formatted bbox table if the
    rectilinear grid does NOT strictly contain the mesh bbox.

    The runtime C++ side performs the same check via
    ``DataField3D::ContainsBBox`` (Phase 4 pre-flight); this writer-side
    counterpart catches the problem one stage earlier.
    """
    if grid_contains_bbox(x, y, z, mesh_bbox, eps=eps):
        return
    grid = grid_bbox(x, y, z)
    raise BBoxContainmentError(
        "data grid does NOT contain the mesh bbox; v1 schema "
        "enforces interpolation-only.\n"
        + _format_bbox_table(grid, mesh_bbox)
        + "\nResolution: either shrink the mesh padding or regenerate "
        "the source dataset over a wider region.")
