"""GOCAD TSurf 1 (.ts) -> STEP (.step) converter for Coreform Cubit.

Each input file contains a single triangulated surface (TFACE) with
VRTX/PVRTX/ATOM vertex records and TRGL triangle records.  This script:

    1. parses the .ts file (reusing audit_ts_quality.parse_tsurf),
    2. drops index-degenerate and zero-area triangles,
    3. builds one OpenCASCADE planar face per surviving triangle,
    4. sews the faces into a single shell (so adjacent triangles share
       their common edges and vertices in the STEP topology), and
    5. writes the result to disk as STEP AP214 via STEPControl_Writer.

The default I/O layout matches the user's request:

    safs/CFM_data/*.ts  ->  safs/CFM_data_step/*.step

so running the script with no arguments converts every .ts file under
CFM_data/ in one go.

Coordinates are kept in their source frame (UTM easting / northing in m,
elevation in m, ZPOSITIVE Elevation), matching the .ts header — the
conversion is purely a format change so a CAD-side workflow in Cubit can
combine the faults with other geometry that uses the same frame.

Requires pythonocc-core (>= 7.9).
"""
from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

# Reuse the project's existing TSurf parser so we share the same handling
# of VRTX / PVRTX / ATOM / TRGL records and triangle index conventions.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_ts_quality import TSurf, parse_tsurf  # noqa: E402

from OCC.Core.BRep import BRep_Builder  # noqa: E402
from OCC.Core.BRepBuilderAPI import (  # noqa: E402
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_Sewing,
)
from OCC.Core.gp import gp_Pnt  # noqa: E402
from OCC.Core.IFSelect import IFSelect_RetDone  # noqa: E402
from OCC.Core.Interface import Interface_Static  # noqa: E402
from OCC.Core.STEPControl import (  # noqa: E402
    STEPControl_AsIs,
    STEPControl_Reader,
    STEPControl_Writer,
)
from OCC.Core.TopAbs import TopAbs_FACE  # noqa: E402
from OCC.Core.TopExp import TopExp_Explorer  # noqa: E402
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Shape  # noqa: E402


# ---------------------------------------------------------------------------
# Triangle mesh -> OCC shape
# ---------------------------------------------------------------------------
def _filter_valid_triangles(V: np.ndarray, T: np.ndarray,
                            min_area: float
                            ) -> tuple[np.ndarray, int, int]:
    """Drop index-degenerate and (near-)zero-area triangles.

    Returns (T_kept, n_index_degenerate, n_zero_area).  The kept array
    has dtype int64 and is safe to feed into the OCC face builder.
    """
    deg = (T[:, 0] == T[:, 1]) | (T[:, 1] == T[:, 2]) | (T[:, 0] == T[:, 2])
    n_idx = int(np.sum(deg))
    T = T[~deg]

    a = V[T[:, 0]]
    b = V[T[:, 1]]
    c = V[T[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    keep = areas >= min_area
    n_zero = int(np.sum(~keep))
    return T[keep], n_idx, n_zero


def build_compound(V: np.ndarray, T: np.ndarray) -> TopoDS_Compound:
    """Build a TopoDS_Compound containing one planar face per triangle.

    Faster than sewing but the resulting STEP file has duplicate edges
    and vertices along every shared triangle boundary — Cubit will see
    each triangle as an independent surface unless told to merge.
    """
    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)

    pts = [gp_Pnt(float(v[0]), float(v[1]), float(v[2])) for v in V]
    skipped = 0
    for tri in T:
        poly = BRepBuilderAPI_MakePolygon()
        poly.Add(pts[int(tri[0])])
        poly.Add(pts[int(tri[1])])
        poly.Add(pts[int(tri[2])])
        poly.Close()
        if not poly.IsDone():
            skipped += 1
            continue
        face_maker = BRepBuilderAPI_MakeFace(poly.Wire(), True)
        if not face_maker.IsDone():
            skipped += 1
            continue
        builder.Add(compound, face_maker.Face())
    if skipped:
        print(f"    [warn] {skipped} triangles failed BRep face construction")
    return compound


def build_sewn_shape(V: np.ndarray, T: np.ndarray,
                     sewing_tol: float) -> TopoDS_Shape:
    """Build all triangle faces and sew them together.

    Sewing identifies coincident edges/vertices between neighboring
    faces and merges them, producing a topologically connected shell
    suitable for downstream CAD operations in Cubit.

    sewing_tol is in metres (the source frame's unit).  A few mm is
    appropriate for the SAFS CFM data, whose triangles meet at vertices
    that are bit-exact in the .ts file.
    """
    pts = [gp_Pnt(float(v[0]), float(v[1]), float(v[2])) for v in V]

    sewer = BRepBuilderAPI_Sewing(sewing_tol)
    skipped = 0
    for tri in T:
        poly = BRepBuilderAPI_MakePolygon()
        poly.Add(pts[int(tri[0])])
        poly.Add(pts[int(tri[1])])
        poly.Add(pts[int(tri[2])])
        poly.Close()
        if not poly.IsDone():
            skipped += 1
            continue
        face_maker = BRepBuilderAPI_MakeFace(poly.Wire(), True)
        if not face_maker.IsDone():
            skipped += 1
            continue
        sewer.Add(face_maker.Face())
    if skipped:
        print(f"    [warn] {skipped} triangles failed BRep face construction")

    sewer.Perform()
    return sewer.SewedShape()


# ---------------------------------------------------------------------------
# STEP writer
# ---------------------------------------------------------------------------
def write_step(shape: TopoDS_Shape, out_path: Path, *, product_name: str,
               schema: str = "AP214") -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Source coordinates are in metres; declare that in the STEP header so
    # importers (Cubit included) interpret the geometry correctly.
    Interface_Static.SetCVal("write.step.unit", "M")
    Interface_Static.SetCVal("write.step.schema", schema)
    Interface_Static.SetCVal("write.step.product.name", product_name)

    writer = STEPControl_Writer()
    status = writer.Transfer(shape, STEPControl_AsIs)
    if status != IFSelect_RetDone:
        raise RuntimeError(f"STEPControl_Writer.Transfer failed for "
                           f"{out_path.name}: status={status}")
    status = writer.Write(str(out_path))
    if status != IFSelect_RetDone:
        raise RuntimeError(f"STEPControl_Writer.Write failed for "
                           f"{out_path}: status={status}")


def _verify_step(out_path: Path) -> int:
    """Read the STEP file back via OCC and count its faces.

    A successful round-trip is the cheapest sanity check that the file we
    just wrote is parseable by another OCC client (Cubit's STEP importer
    uses the same OCCT kernel underneath).
    """
    reader = STEPControl_Reader()
    status = reader.ReadFile(str(out_path))
    if status != IFSelect_RetDone:
        raise RuntimeError(f"STEPControl_Reader.ReadFile failed for "
                           f"{out_path}: status={status}")
    n_roots = reader.TransferRoots()
    if n_roots <= 0:
        raise RuntimeError(f"STEP read produced no roots for {out_path}")
    shape = reader.OneShape()
    expl = TopExp_Explorer(shape, TopAbs_FACE)
    n_faces = 0
    while expl.More():
        n_faces += 1
        expl.Next()
    return n_faces


# ---------------------------------------------------------------------------
# Per-file driver
# ---------------------------------------------------------------------------
@dataclass
class ConversionResult:
    ts_path: Path
    step_path: Path
    n_vrtx_in: int
    n_trgl_in: int
    n_trgl_out: int
    n_index_degenerate: int
    n_zero_area: int
    n_faces_roundtrip: int | None
    elapsed_s: float


def convert_one(ts_path: Path, step_path: Path, *,
                sewing_tol: float, sew: bool, min_area: float,
                step_schema: str = "AP214",
                verify: bool = True) -> ConversionResult:
    t0 = time.time()
    ts: TSurf = parse_tsurf(ts_path)
    n_v_in = int(ts.V.shape[0])
    n_t_in = int(ts.T.shape[0])

    T_clean, n_idx, n_zero = _filter_valid_triangles(ts.V, ts.T, min_area)
    if T_clean.shape[0] == 0:
        raise RuntimeError(f"{ts_path.name}: no valid triangles after "
                            f"degeneracy filtering")

    if sew:
        shape = build_sewn_shape(ts.V, T_clean, sewing_tol=sewing_tol)
    else:
        shape = build_compound(ts.V, T_clean)

    product_name = ts.name or ts_path.stem
    write_step(shape, step_path, product_name=product_name, schema=step_schema)

    n_faces_rt: int | None = None
    if verify:
        n_faces_rt = _verify_step(step_path)

    return ConversionResult(
        ts_path=ts_path,
        step_path=step_path,
        n_vrtx_in=n_v_in,
        n_trgl_in=n_t_in,
        n_trgl_out=int(T_clean.shape[0]),
        n_index_degenerate=n_idx,
        n_zero_area=n_zero,
        n_faces_roundtrip=n_faces_rt,
        elapsed_s=time.time() - t0,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
_DEFAULT_IN = Path(__file__).resolve().parents[1] / "CFM_data"
_DEFAULT_OUT = Path(__file__).resolve().parents[1] / "CFM_data_step"


def _select_inputs(in_dir: Path, res: int | None,
                   include: list[str] | None) -> list[Path]:
    pattern = "*.ts" if res is None else f"*_{res}m.ts"
    files = sorted(in_dir.glob(pattern))
    if include:
        sel = [f for f in files
               if any(token in f.name for token in include)]
        if not sel:
            raise FileNotFoundError(
                f"--include matched none of the files under {in_dir}"
            )
        files = sel
    if not files:
        raise FileNotFoundError(
            f"no {pattern} files found under {in_dir}"
        )
    return files


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--in-dir", type=Path, default=_DEFAULT_IN,
                        help=f"Directory of input .ts files "
                             f"(default: {_DEFAULT_IN})")
    parser.add_argument("--out-dir", type=Path, default=_DEFAULT_OUT,
                        help=f"Directory for output .step files "
                             f"(default: {_DEFAULT_OUT})")
    parser.add_argument("--res", type=int, choices=(500, 1000, 2000),
                        default=None,
                        help="Restrict to a single CFM resolution band "
                             "(matches *_<res>m.ts).")
    parser.add_argument("--include", action="append", default=None,
                        metavar="SUBSTRING",
                        help="Convert only files whose name contains the "
                             "given substring; may be repeated.")
    parser.add_argument("--no-sew", action="store_true",
                        help="Skip OCC sewing and emit a Compound of "
                             "independent triangle faces.  Faster to "
                             "write but produces a larger STEP file with "
                             "no shared topology.")
    parser.add_argument("--sewing-tol", type=float, default=1e-3,
                        metavar="METRES",
                        help="Tolerance for OCC vertex/edge merging "
                             "during sewing.  Default: 1e-3 m (1 mm).")
    parser.add_argument("--min-area", type=float, default=1e-6,
                        metavar="M^2",
                        help="Drop triangles whose area is below this "
                             "threshold (in m^2).  Default: 1e-6 m^2.")
    parser.add_argument("--schema", choices=("AP203", "AP214", "AP242"),
                        default="AP214",
                        help="STEP application protocol.  Default: AP214.")
    parser.add_argument("--no-verify", action="store_true",
                        help="Skip the round-trip read-back check after "
                             "writing each STEP file.")
    parser.add_argument("--dry-run", action="store_true",
                        help="List the files that would be converted "
                             "and exit without writing anything.")
    args = parser.parse_args(argv)

    in_dir: Path = args.in_dir
    out_dir: Path = args.out_dir
    if not in_dir.is_dir():
        parser.error(f"--in-dir does not exist: {in_dir}")

    files = _select_inputs(in_dir, args.res, args.include)
    print(f"ts_to_step: {len(files)} input file(s) under {in_dir}")
    print(f"            -> writing to {out_dir}")
    if args.dry_run:
        for p in files:
            print(f"  would convert: {p.name}")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[ConversionResult] = []
    failures: list[tuple[Path, Exception]] = []

    for i, p in enumerate(files, 1):
        out = out_dir / (p.stem + ".step")
        print(f"[{i:2d}/{len(files)}] {p.name}")
        try:
            res = convert_one(
                p, out,
                sewing_tol=args.sewing_tol,
                sew=not args.no_sew,
                min_area=args.min_area,
                step_schema=args.schema,
                verify=not args.no_verify,
            )
        except Exception as exc:  # noqa: BLE001 - report and continue
            print(f"    [FAIL] {exc}")
            failures.append((p, exc))
            continue
        size_mb = out.stat().st_size / (1024 * 1024)
        print(f"    in : {res.n_vrtx_in} vrtx, {res.n_trgl_in} trgl "
              f"(dropped {res.n_index_degenerate} idx-degenerate, "
              f"{res.n_zero_area} zero-area)")
        verify_str = ("" if res.n_faces_roundtrip is None
                       else f", round-trip {res.n_faces_roundtrip} faces")
        print(f"    out: {res.n_trgl_out} faces -> {out.name} "
              f"({size_mb:.1f} MB, {res.elapsed_s:.1f} s{verify_str})")
        results.append(res)

    print()
    print(f"ts_to_step done: {len(results)} succeeded, {len(failures)} failed")
    if failures:
        for p, exc in failures:
            print(f"  FAIL {p.name}: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
