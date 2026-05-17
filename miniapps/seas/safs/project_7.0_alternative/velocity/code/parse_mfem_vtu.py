"""parse_mfem_vtu.py — VTK-backed reader for MFEM .vtu output.

Used by the velocity/ plot scripts.  MFEM emits VTU 2.2 with embedded
appended-base64 (default ``--binary``) or pure ASCII (``--ascii``);
``meshio`` chokes on the binary dialect, so the historical workflow
shells out to a ``vtkXMLUnstructuredGridReader`` instead.  This module
is the in-tree replacement for the legacy ``/tmp/parse_mfem_vtu.py``
that the plot scripts add to ``sys.path``.

Public API:

    parse_vtu(path) -> dict with keys
        points      : (n_points, 3) float64
        conn        : (sum_cell_corner_counts,) int64 — flat connectivity
        offsets     : (n_cells,) int64 — last-index-exclusive offsets
                      (matches the convention in the .vtu file itself)
        cell_types  : (n_cells,) uint8 — VTK cell-type codes (tet=10)
        point_data  : dict[name -> (n_points,) float64]

The connectivity uses VTK's flat layout (matches what MFEM writes) so
``offsets[k]`` is the END index of cell k's corner list, and
``offsets[k-1]:offsets[k]`` (with offsets[-1] interpreted as 0) gives
cell k's vertex IDs.  Callers that assume per-cell corner counts are
constant (e.g. all-tet meshes) can read it as
``conn.reshape(n_cells, n_per_cell)``.

Dependency: ``vtk`` (already installed in conda env ``pythonenv``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def parse_vtu(path: str | Path) -> dict[str, Any]:
    """Read an MFEM (.vtu) file and return its data as numpy arrays.

    See module docstring for the dict schema.

    Raises FileNotFoundError if ``path`` does not exist, RuntimeError if
    the VTK reader cannot parse the file (typically a corrupted or
    truncated VTU), and KeyError if the file has zero PointData arrays.
    """
    import vtk
    from vtk.util import numpy_support as vns

    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"parse_mfem_vtu: no such file: {p}")

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(p))
    reader.Update()
    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfPoints() == 0:
        raise RuntimeError(
            f"parse_mfem_vtu: vtkXMLUnstructuredGridReader returned an "
            f"empty grid for {p}; likely a corrupted or unsupported "
            f"VTU dialect (MFEM may have changed its writer)")

    points = vns.vtk_to_numpy(grid.GetPoints().GetData())
    points = np.ascontiguousarray(points, dtype=np.float64)

    cells = grid.GetCells()
    conn = vns.vtk_to_numpy(cells.GetConnectivityArray()).astype(
        np.int64, copy=False)
    offsets_raw = vns.vtk_to_numpy(cells.GetOffsetsArray()).astype(
        np.int64, copy=False)
    # vtkCellArray.GetOffsetsArray() returns (n_cells + 1,) values
    # starting with 0; the historical plot-script convention uses the
    # (n_cells,) "end-of-cell" offsets, so drop the leading zero.
    if offsets_raw.size >= 1 and offsets_raw[0] == 0:
        offsets = offsets_raw[1:].copy()
    else:
        offsets = offsets_raw.copy()

    cell_types = vns.vtk_to_numpy(grid.GetCellTypesArray()).astype(
        np.uint8, copy=False)

    pd = grid.GetPointData()
    n_arr = pd.GetNumberOfArrays()
    if n_arr == 0:
        raise KeyError(
            f"parse_mfem_vtu: {p} has no PointData arrays (the velocity "
            f"projector should have emitted Vp/Vs/density at minimum)")
    point_data: dict[str, np.ndarray] = {}
    for i in range(n_arr):
        a = pd.GetArray(i)
        name = a.GetName()
        arr = vns.vtk_to_numpy(a)
        # MFEM writes (n_components, n_points) for vector fields; we
        # only need scalar Vp/Vs/density here, but accept any shape.
        point_data[name] = np.ascontiguousarray(arr, dtype=np.float64)

    return {
        "points": points,
        "conn": conn,
        "offsets": offsets,
        "cell_types": cell_types,
        "point_data": point_data,
    }


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("usage: parse_mfem_vtu.py PATH", file=sys.stderr)
        sys.exit(2)
    m = parse_vtu(sys.argv[1])
    print(f"points     : {m['points'].shape}")
    print(f"conn       : {m['conn'].shape}")
    print(f"offsets    : {m['offsets'].shape}")
    print(f"cell_types : {m['cell_types'].shape} "
          f"(unique: {np.unique(m['cell_types']).tolist()})")
    for k, v in m["point_data"].items():
        print(f"point_data[{k!r}] : shape {v.shape} "
              f"range [{float(np.nanmin(v)):.3f}, {float(np.nanmax(v)):.3f}]")
