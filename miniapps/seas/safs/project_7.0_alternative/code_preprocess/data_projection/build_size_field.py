"""build_size_field.py — produce a gmsh background size field driven by
the gradient of a sidecar field (Vs by default).

For every voxel in the sidecar, the target tet edge length is::

    LC = clip(lc_far / (1 + alpha * gmag_norm), lc_min, lc_far)

where ``gmag_norm`` is ``‖∇F‖`` normalised to its sidecar-wide max
(after optional Gaussian smoothing).  The result is written as a gmsh
ScalarPoint (``SP``) PostView ``.pos`` file that the .geo template can
``Merge`` and combine with the existing distance-to-fault size field
via ``Field[Min]``.

CLI::

    python build_size_field.py
        --sidecar      PATH                    # required
        --field        Vs                      # default Vs
        --out-pos      PATH                    # required
        --lc-near      1500                    # default
        --lc-far       10000                   # default
        --lc-min       500                     # absolute floor
        --alpha        8.0                     # gradient sensitivity
        --smooth-sigma 1.0                     # Gaussian smoothing,
                                                # voxel units
        --voxel-stride 2                       # PostView downsample
        --quiet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np


# ---------------------------------------------------------------------
# Pure numerics
# ---------------------------------------------------------------------

def compute_size_field(F: np.ndarray,
                       gx: np.ndarray, gy: np.ndarray, gz: np.ndarray,
                       *,
                       lc_near: float, lc_far: float, lc_min: float,
                       alpha: float, smooth_sigma: float = 0.0
                       ) -> np.ndarray:
    """Per-voxel target tet edge length, in metres.

    Inputs are in canonical UTM 11N units (m); ``gx, gy, gz`` are 1-D
    monotone increasing axes; ``F`` has shape ``(Nx, Ny, Nz)``.

    The ``smooth_sigma`` argument is in **voxel units** (sigma=1 means
    one-voxel Gaussian kernel along each axis); a value of 0 disables
    smoothing.

    The function only takes the maximum / clip + the formula in the
    module docstring -- no I/O, no gmsh dependency.  Suitable for
    isolated unit testing.
    """
    if not (lc_min > 0.0 and lc_near > 0.0 and lc_far > 0.0):
        raise ValueError(
            f"compute_size_field: lc_min/lc_near/lc_far must be > 0; "
            f"got ({lc_min}, {lc_near}, {lc_far})")
    if not (lc_min <= lc_far):
        raise ValueError(
            f"compute_size_field: lc_min ({lc_min}) > lc_far "
            f"({lc_far}); the floor cannot exceed the ceiling")
    if alpha < 0.0:
        raise ValueError(f"compute_size_field: alpha={alpha}; must be >=0")
    if F.ndim != 3:
        raise ValueError(
            f"compute_size_field: F must be 3-D; got {F.shape}")
    if F.shape != (gx.size, gy.size, gz.size):
        raise ValueError(
            f"compute_size_field: F.shape {F.shape} != "
            f"(gx, gy, gz) sizes ({gx.size}, {gy.size}, {gz.size})")
    for arr, name in ((gx, "gx"), (gy, "gy"), (gz, "gz")):
        if arr.ndim != 1 or arr.size < 2:
            raise ValueError(
                f"compute_size_field: axis '{name}' must be 1-D with "
                f"length >= 2; got shape {arr.shape}")
        if not np.all(np.diff(arr) > 0.0):
            raise ValueError(
                f"compute_size_field: axis '{name}' must be strictly "
                f"monotone increasing")

    dx_axis = float(np.diff(gx).mean())
    dy_axis = float(np.diff(gy).mean())
    dz_axis = float(np.diff(gz).mean())

    gFx, gFy, gFz = np.gradient(F, dx_axis, dy_axis, dz_axis,
                                edge_order=2)
    gmag = np.sqrt(gFx * gFx + gFy * gFy + gFz * gFz)
    if smooth_sigma > 0.0:
        from scipy.ndimage import gaussian_filter
        gmag = gaussian_filter(gmag, sigma=smooth_sigma)
    gmag_max = float(gmag.max())
    if gmag_max <= 0.0:
        # Constant field — return uniform lc_far.
        return np.full(F.shape, lc_far, dtype=np.float64)
    gmag_norm = gmag / gmag_max
    lc = lc_far / (1.0 + alpha * gmag_norm)
    np.clip(lc, lc_min, lc_far, out=lc)
    return lc.astype(np.float64, copy=False)


# ---------------------------------------------------------------------
# .pos writer
# ---------------------------------------------------------------------

def _stride_axis(arr: np.ndarray, stride: int, name: str) -> np.ndarray:
    """Take every ``stride``-th element of ``arr``; always include the
    last element so the downsampled grid still spans the original
    bbox."""
    if stride < 1:
        raise ValueError(f"_stride_axis: stride must be >=1; got {stride}")
    out = arr[::stride]
    if out.size == 0 or out[-1] != arr[-1]:
        out = np.concatenate([out, arr[-1:]])
    if out.size < 4:
        raise ValueError(
            f"--voxel-stride={stride} drops axis '{name}' to "
            f"{out.size} voxels; gmsh PostView interpolant requires "
            f">= 4")
    return out


def write_pos_scalar_point(out_path: Path,
                           x: np.ndarray, y: np.ndarray, z: np.ndarray,
                           values: np.ndarray,
                           *,
                           view_name: str = "size_field") -> None:
    """Write a gmsh PostView ``.pos`` (legacy ASCII) using the
    ``ScalarPoint`` (``SP``) primitive.  ``values`` is an (Nx, Ny, Nz)
    array, indexed as ``values[i, j, k] = LC at (x[i], y[j], z[k])``.
    """
    if values.shape != (x.size, y.size, z.size):
        raise ValueError(
            f"write_pos_scalar_point: values shape {values.shape} != "
            f"({x.size}, {y.size}, {z.size})")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Build the body in-memory with a generator and join — fast enough
    # for the typical ~5000-voxel downsampled SAFS sidecar.
    lines: list[str] = [f'View "{view_name}" {{']
    for i, xv in enumerate(x):
        for j, yv in enumerate(y):
            for k, zv in enumerate(z):
                lc = float(values[i, j, k])
                # SP(x,y,z){lc};   — single time step.
                lines.append(f"SP({xv:.6f},{yv:.6f},{zv:.6f}){{{lc:.6f}}};")
    lines.append("};\n")
    out_path.write_text("\n".join(lines))


# ---------------------------------------------------------------------
# Public driver
# ---------------------------------------------------------------------

def build_pos(*,
              sidecar: Path,
              out_pos: Path,
              field: str = "Vs",
              lc_near: float = 1500.0,
              lc_far: float = 10000.0,
              lc_min: float = 500.0,
              alpha: float = 8.0,
              smooth_sigma: float = 1.0,
              voxel_stride: int = 2,
              verbose: bool = True) -> dict:
    """Read the sidecar field, compute the size field, write the .pos.

    Returns a metadata dict (LC stats, grid sizes) for callers /
    integration tests.
    """
    sidecar = Path(sidecar)
    out_pos = Path(out_pos)
    if not sidecar.is_file():
        raise FileNotFoundError(f"--sidecar not found: {sidecar}")

    with h5py.File(sidecar, "r") as h5:
        if "fields" not in h5 or field not in h5["fields"]:
            available = (list(h5["fields"].keys())
                         if "fields" in h5 else [])
            raise KeyError(
                f"sidecar '{sidecar}' has no field '{field}'; "
                f"available: {available}")
        gx = np.asarray(h5["grid/x"][...], dtype=np.float64)
        gy = np.asarray(h5["grid/y"][...], dtype=np.float64)
        gz = np.asarray(h5["grid/z"][...], dtype=np.float64)
        F = np.asarray(h5[f"fields/{field}"][...], dtype=np.float64)

    lc = compute_size_field(F, gx, gy, gz,
                            lc_near=lc_near, lc_far=lc_far,
                            lc_min=lc_min, alpha=alpha,
                            smooth_sigma=smooth_sigma)

    if voxel_stride > 1:
        sx = _stride_axis(gx, voxel_stride, "x")
        sy = _stride_axis(gy, voxel_stride, "y")
        sz = _stride_axis(gz, voxel_stride, "z")
        # Find the index in (gx, gy, gz) for every entry of (sx, sy, sz)
        # via searchsorted (axes are strictly monotone).
        ix = np.searchsorted(gx, sx)
        iy = np.searchsorted(gy, sy)
        iz = np.searchsorted(gz, sz)
        lc = lc[np.ix_(ix, iy, iz)]
    else:
        sx, sy, sz = gx, gy, gz

    write_pos_scalar_point(out_pos, sx, sy, sz, lc)

    meta = {
        "out_pos":        str(out_pos),
        "n_voxels":       int(lc.size),
        "lc_min":         float(lc.min()),
        "lc_max":         float(lc.max()),
        "lc_mean":        float(lc.mean()),
        "lc_floor_count": int((lc <= lc_min + 1.0e-9).sum()),
        "voxel_stride":   int(voxel_stride),
        "alpha":          float(alpha),
        "smooth_sigma":   float(smooth_sigma),
        "field":          field,
        "axes_n":         (int(sx.size), int(sy.size), int(sz.size)),
    }
    if verbose:
        print(f"size field {field}: "
              f"lc range = [{meta['lc_min']:.0f}, {meta['lc_max']:.0f}] "
              f"m, mean = {meta['lc_mean']:.0f} m;  "
              f"{meta['lc_floor_count']} of {meta['n_voxels']} voxels "
              f"hit the lc_min={lc_min:.0f} m floor; downsampled grid "
              f"{meta['axes_n']}.")
        print(f"wrote {out_pos}")
    return meta


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sidecar", type=Path, required=True)
    ap.add_argument("--field", default="Vs")
    ap.add_argument("--out-pos", type=Path, required=True,
                    dest="out_pos")
    ap.add_argument("--lc-near", type=float, default=1500.0,
                    dest="lc_near")
    ap.add_argument("--lc-far", type=float, default=10000.0,
                    dest="lc_far")
    ap.add_argument("--lc-min", type=float, default=500.0,
                    dest="lc_min")
    ap.add_argument("--alpha", type=float, default=8.0)
    ap.add_argument("--smooth-sigma", type=float, default=1.0,
                    dest="smooth_sigma")
    ap.add_argument("--voxel-stride", type=int, default=2,
                    dest="voxel_stride")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    try:
        build_pos(sidecar=args.sidecar,
                  out_pos=args.out_pos,
                  field=args.field,
                  lc_near=args.lc_near,
                  lc_far=args.lc_far,
                  lc_min=args.lc_min,
                  alpha=args.alpha,
                  smooth_sigma=args.smooth_sigma,
                  voxel_stride=args.voxel_stride,
                  verbose=not args.quiet)
    except (FileNotFoundError, KeyError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
