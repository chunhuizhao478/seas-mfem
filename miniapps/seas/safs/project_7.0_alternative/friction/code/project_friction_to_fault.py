#!/usr/bin/env python3
"""
project_friction_to_fault.py — Project the depth-varying rate-and-state
friction parameters a(z) and b(z) onto the SAFS fault mesh and write a VTU,
mirroring the stress projection (``stress/code/project_to_fault_stress.py``).

Inputs
------
- A fault VTU (``*_fault.vtu`` written by ``meshing/code/msh_to_vtu.py``).
- Two depth-profile CSVs in the same format the C++ driver consumes
  (``spatial_friction.cpp:LoadFrictionDepthProfileCSVs``):
  each data row is ``value, depth_km`` (value FIRST); depth ascending in km.
    * ``param_a.csv``         -> a(depth)
    * ``param_a_minus_b.csv`` -> (a-b)(depth)
  b(depth) is reconstructed as ``b = a - (a-b)``.  The interpolation is
  piecewise-linear with FLAT (constant) clamping outside the sampled depth
  range, identical to the C++ ``PiecewiseLinear1D``.  Depth at a mesh node /
  cell is ``max(0, -z)`` (z in metres, up-positive), matching the resolver.

Output
------
``<out-dir>/<base>_friction_ab.vtu`` with:
  - point_data (continuous, area-weighted cell->node average; NaN on points
    not referenced by any fault triangle, as in the stress writer):
      ``a``, ``b``, ``a_minus_b``, ``depth_km``, ``vw_indicator``
      (vw_indicator = 1.0 where a-b < 0 (velocity-weakening), else 0.0)
  - cell_data (raw per-triangle, evaluated at the cell centroid depth):
      ``a_cell``, ``b_cell``, ``a_minus_b_cell``, ``depth_km_cell``,
      ``vw_indicator_cell``

Usage
-----
    conda activate pythonenv
    python project_friction_to_fault.py \
        [--fault-vtu  ../../meshing/results/vtu/safs_fault_box_nwcut_500m_lcfar3000_z0embed_fault.vtu] \
        [--param-a-csv          ../rate-and-state/param_a_vwvs11km.csv] \
        [--param-a-minus-b-csv  ../rate-and-state/param_a_minus_b_vwvs11km.csv] \
        [--fault-phys-name fault] \
        [--out-dir ../vtu] [--out-name <auto>]

ParaView recipe: open the VTU, Surface representation, colour by ``a_minus_b``
(diverging colourmap centred at 0) to see the VW (<0) / VS (>0) split; the
``vw_indicator`` field shows the transition as a hard 0/1 boundary near 11 km.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------
# Reuse the stress-projection mesh helpers (same fault VTU contract).
# ----------------------------------------------------------------------
_STRESS_CODE_DIR = (
    Path(__file__).resolve().parent.parent.parent / "stress" / "code"
)
if str(_STRESS_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_STRESS_CODE_DIR))
try:
    from project_to_fault_stress import (  # noqa: E402
        load_fault_mesh,
        extract_cells_by_physical,
        triangle_geometry,
        cell_to_node_average,
    )
except Exception as exc:  # pragma: no cover - import guard
    raise ImportError(
        f"could not import the stress-projection helpers from "
        f"{_STRESS_CODE_DIR}; this script reuses load_fault_mesh / "
        f"extract_cells_by_physical / triangle_geometry / "
        f"cell_to_node_average. Underlying error: {exc!r}"
    ) from exc

import meshio  # noqa: E402  (after sys.path tweak, mirrors stress module)


# ----------------------------------------------------------------------
# Defaults (relative to this file's directory)
# ----------------------------------------------------------------------
_HERE = Path(__file__).resolve().parent
DEFAULT_FAULT_VTU = (
    _HERE.parent.parent
    / "meshing" / "results" / "vtu"
    / "safs_fault_box_nwcut_500m_lcfar3000_z0embed_fault.vtu"
)
DEFAULT_PARAM_A_CSV = _HERE.parent / "rate-and-state" / "param_a_vwvs11km.csv"
DEFAULT_PARAM_AMB_CSV = (
    _HERE.parent / "rate-and-state" / "param_a_minus_b_vwvs11km.csv"
)
DEFAULT_OUT_DIR = _HERE.parent / "vtu"
DEFAULT_FAULT_PHYS_NAME = "fault"


# ----------------------------------------------------------------------
# CSV profile reader (matches the C++ LoadFrictionDepthProfileCSVs contract)
# ----------------------------------------------------------------------
def read_profile_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a ``value, depth_km`` profile CSV.

    Returns ``(depth_km, value)`` arrays, sorted-by-construction ascending in
    depth.  Mirrors the C++ validators: >= 2 rows, strictly increasing depth,
    all fields finite.

    Raises
    ------
    FileNotFoundError, ValueError
        With the offending path / row in the message.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"profile CSV not found: {p.resolve()}")

    depths: list[float] = []
    values: list[float] = []
    with p.open() as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [tok.strip() for tok in line.split(",")]
            if len(parts) < 2:
                raise ValueError(
                    f"{p.resolve()}:{lineno}: expected 'value, depth_km'; "
                    f"got {line!r}"
                )
            try:
                value = float(parts[0])
                depth_km = float(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"{p.resolve()}:{lineno}: non-numeric field in {line!r}"
                ) from exc
            if not (np.isfinite(value) and np.isfinite(depth_km)):
                raise ValueError(
                    f"{p.resolve()}:{lineno}: non-finite field in {line!r}"
                )
            values.append(value)
            depths.append(depth_km)

    if len(depths) < 2:
        raise ValueError(
            f"{p.resolve()}: need >= 2 data rows; got {len(depths)}"
        )
    depth_arr = np.asarray(depths, dtype=np.float64)
    value_arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.diff(depth_arr) > 0.0):
        raise ValueError(
            f"{p.resolve()}: depth column must be strictly increasing; "
            f"got {depth_arr.tolist()}"
        )
    return depth_arr, value_arr


def eval_profile(
    depth_km_grid: np.ndarray,
    value_grid: np.ndarray,
    depth_km_query: np.ndarray,
) -> np.ndarray:
    """Flat-clamped piecewise-linear interpolation (matches C++
    ``PiecewiseLinear1D``: linear inside, constant clamp outside).

    ``np.interp`` clamps to the endpoint values outside ``[x[0], x[-1]]`` by
    default, which is exactly the C++ FLAT clamp.
    """
    return np.interp(depth_km_query, depth_km_grid, value_grid)


def transition_depth_km(
    amb_depth_km: np.ndarray, amb_value: np.ndarray
) -> float | None:
    """Depth (km) where (a-b) crosses zero (VW->VS), found by linear
    interpolation on the input profile.  Returns None if (a-b) does not
    change sign across the sampled range."""
    v = amb_value
    for i in range(len(v) - 1):
        v0, v1 = v[i], v[i + 1]
        if v0 == 0.0:
            return float(amb_depth_km[i])
        if v0 < 0.0 <= v1 or v0 > 0.0 >= v1:
            d0, d1 = amb_depth_km[i], amb_depth_km[i + 1]
            # linear root of v(d) between the two samples
            return float(d0 + (0.0 - v0) * (d1 - d0) / (v1 - v0))
    return None


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--fault-vtu", type=Path, default=DEFAULT_FAULT_VTU,
                    help="fault VTU (default: 500 m z0embed fault mesh)")
    ap.add_argument("--param-a-csv", type=Path, default=DEFAULT_PARAM_A_CSV,
                    help="a(depth) profile CSV")
    ap.add_argument("--param-a-minus-b-csv", type=Path,
                    default=DEFAULT_PARAM_AMB_CSV,
                    help="(a-b)(depth) profile CSV")
    ap.add_argument("--fault-phys-name", default=DEFAULT_FAULT_PHYS_NAME,
                    help="gmsh physical-group name of the fault triangles")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                    help="output directory for the VTU")
    ap.add_argument("--out-name", default=None,
                    help="output VTU filename (default: "
                         "<fault-vtu-base sans _fault>_friction_ab.vtu)")
    args = ap.parse_args(argv)

    # --- Load profiles ---
    a_depth, a_val = read_profile_csv(args.param_a_csv)
    amb_depth, amb_val = read_profile_csv(args.param_a_minus_b_csv)
    if not np.all(a_val > 0.0):
        raise ValueError(
            f"{Path(args.param_a_csv).resolve()}: all a values must be > 0; "
            f"got {a_val.tolist()}"
        )

    tdepth = transition_depth_km(amb_depth, amb_val)

    # --- Load fault mesh + per-triangle geometry ---
    mesh = load_fault_mesh(args.fault_vtu)
    tri_conn = extract_cells_by_physical(mesh, "triangle", args.fault_phys_name)
    if tri_conn.shape[0] == 0:
        raise ValueError(
            f"no fault triangles found in {Path(args.fault_vtu).resolve()} "
            f"for physical group {args.fault_phys_name!r}"
        )
    points = np.asarray(mesh.points, dtype=np.float64)
    centroids, _normals, areas = triangle_geometry(points, tri_conn)

    # --- Depth (km) at cell centroids and at nodes: max(0, -z) / 1000 ---
    depth_cell_km = np.maximum(0.0, -centroids[:, 2]) / 1000.0

    # --- Evaluate a, (a-b), b at cell centroids ---
    a_cell = eval_profile(a_depth, a_val, depth_cell_km)
    amb_cell = eval_profile(amb_depth, amb_val, depth_cell_km)
    b_cell = a_cell - amb_cell
    vw_cell = (amb_cell < 0.0).astype(np.float64)

    # --- Cell -> node (area-weighted; NaN on orphan/non-fault points) ---
    a_node = cell_to_node_average(points, tri_conn, a_cell, areas)
    amb_node = cell_to_node_average(points, tri_conn, amb_cell, areas)
    depth_node = cell_to_node_average(points, tri_conn, depth_cell_km, areas)
    b_node = a_node - amb_node
    # VW indicator at the node from the averaged (a-b); leave NaN where the
    # node is an orphan (a-b averaged to NaN there).
    vw_node = np.where(np.isfinite(amb_node), (amb_node < 0.0).astype(np.float64),
                       np.nan)

    # --- Write VTU ---
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.out_name:
        out_name = args.out_name
    else:
        base = Path(args.fault_vtu).stem
        if base.endswith("_fault"):
            base = base[: -len("_fault")]
        out_name = f"{base}_friction_ab.vtu"
    out_path = out_dir / out_name

    fault_out = meshio.Mesh(
        points=points,
        cells=[("triangle", tri_conn)],
        point_data={
            "a": a_node,
            "b": b_node,
            "a_minus_b": amb_node,
            "depth_km": depth_node,
            "vw_indicator": vw_node,
        },
        cell_data={
            "a_cell": [a_cell],
            "b_cell": [b_cell],
            "a_minus_b_cell": [amb_cell],
            "depth_km_cell": [depth_cell_km],
            "vw_indicator_cell": [vw_cell],
        },
    )
    fault_out.info = {
        "fields": "rate-and-state a, b = a-(a-b); depth_km; vw_indicator (a-b<0)",
        "param_a_csv": str(Path(args.param_a_csv).resolve()),
        "param_a_minus_b_csv": str(Path(args.param_a_minus_b_csv).resolve()),
        "vw_vs_transition_km": "n/a" if tdepth is None else f"{tdepth:.6g}",
        "frame": "(east, north, up) / UTM Zone 11N; depth = max(0,-z)",
    }
    meshio.write(str(out_path), fault_out)

    # --- Summary / verification ---
    vw_cells = int((amb_cell < 0.0).sum())
    vs_cells = int((amb_cell >= 0.0).sum())
    shallowest_vs = (
        float(depth_cell_km[amb_cell >= 0.0].min()) if vs_cells else float("nan")
    )
    deepest_vw = (
        float(depth_cell_km[amb_cell < 0.0].max()) if vw_cells else float("nan")
    )
    print(f"fault VTU       : {Path(args.fault_vtu).resolve()}")
    print(f"a   profile     : {Path(args.param_a_csv).resolve()}")
    print(f"a-b profile     : {Path(args.param_a_minus_b_csv).resolve()}")
    print(f"fault triangles : {tri_conn.shape[0]}")
    print(f"fault depth (km): [{depth_cell_km.min():.3f}, "
          f"{depth_cell_km.max():.3f}]")
    print(f"VW->VS transition (profile zero-crossing): "
          f"{'n/a' if tdepth is None else f'{tdepth:.4f} km'}")
    print(f"on-fault cells  : VW (a-b<0) = {vw_cells}, VS (a-b>=0) = {vs_cells}")
    print(f"  deepest VW cell  = {deepest_vw:.4f} km")
    print(f"  shallowest VS cell = {shallowest_vs:.4f} km")
    print(f"a   range (cells): [{a_cell.min():.6g}, {a_cell.max():.6g}]")
    print(f"b   range (cells): [{b_cell.min():.6g}, {b_cell.max():.6g}]")
    print(f"a-b range (cells): [{amb_cell.min():.6g}, {amb_cell.max():.6g}]")
    print(f"wrote           : {out_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
