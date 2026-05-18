#!/usr/bin/env python3
# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory. All rights
# reserved.  See files LICENSE and NOTICE for details.  LLNL-CODE-806117.
#
# Phase 5 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
#
# Post-run packaging tool for legacy ASCII per-rank-per-cycle fault VTU
# directories.  Reads `<dir>/fault_surface.pvd` to enumerate cycles and
# their rank pieces, then optionally:
#
#   --out-tar BASENAME    -> emit `BASENAME_chunk0001.tar.zst` etc., each
#                            containing N=100 cycles' worth of files.
#                            Idempotent: existing chunk files are
#                            preserved (re-emit by deleting the chunk).
#                            Falls back to `tar.gz` if `zstd` is missing.
#   --out-hdf5 PATH       -> emit a single VTKHDF using `pyvista` (or
#                            `vtk` Python bindings) by merging per-rank
#                            pieces and replaying timesteps.
#
# This script does NOT require any seas-mfem build artefacts.  It runs
# on the same Lustre directory the simulation produced and is the
# recommended way to pack pre-Phase-1 outputs that already exist on
# disk.  For Phase-1+ outputs (single binary VTU per cycle) the legacy
# per-rank merge step is automatically a no-op.

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# PVD enumeration
# ---------------------------------------------------------------------------

@dataclass
class CycleEntry:
    cycle: int
    time: float
    file: str         # path relative to the FaultSurface directory


def parse_pvd(pvd_path: Path) -> List[CycleEntry]:
    """Read fault_surface.pvd and return one CycleEntry per <DataSet>.

    Cycle index is parsed out of the filename ('..._c123.vtu' or
    '..._c123.pvtu') because the legacy writer's <DataSet> tag does not
    expose it directly.  Time is read verbatim from the `timestep`
    attribute.
    """
    if not pvd_path.is_file():
        raise FileNotFoundError(f"PVD not found: {pvd_path}")
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    cycle_re = re.compile(r"_c(\d+)\.(p?vtu)$")
    entries: List[CycleEntry] = []
    for ds in root.iter("DataSet"):
        ts = ds.get("timestep")
        fn = ds.get("file")
        if ts is None or fn is None:
            continue
        m = cycle_re.search(fn)
        cycle = int(m.group(1)) if m else len(entries)
        entries.append(CycleEntry(cycle=cycle, time=float(ts), file=fn))
    entries.sort(key=lambda e: e.cycle)
    return entries


def detect_layout(fault_dir: Path) -> str:
    """R-311 / plan §Phase 5 Edge Cases: classify a FaultSurface dir as
    'phase1' (single-file VTU per cycle), 'legacy' (per-rank ASCII
    pieces), or 'mixed' (both).  Recognition uses the `_r` substring
    in filenames per the plan's wording: "recognising the layout by
    presence of `_r0_` substring in filenames"."""
    has_per_rank = any(
        p.is_file() and "_r" in p.name
        for p in fault_dir.glob("fault_surface_r*_c*.vtu"))
    has_single = any(
        p.is_file() and "_r" not in p.name
        for p in fault_dir.glob("fault_surface_c*.vtu"))
    if has_per_rank and has_single:
        return "mixed"
    if has_per_rank:
        return "legacy"
    return "phase1"


def collect_cycle_files(fault_dir: Path,
                        entry: CycleEntry) -> List[Path]:
    """Return the list of on-disk VTU pieces for one cycle.

    For a Phase-1+ binary writer entry (`fault_surface_c{c}.vtu`), the
    list is just that one file.  For the legacy per-rank ASCII writer
    (`fault_surface_c{c}.pvtu` index), the list is the PVTU plus every
    referenced `fault_surface_r{r}_c{c}.vtu` piece, looked up via glob.
    """
    head = fault_dir / entry.file
    files = [head]
    if head.suffix == ".pvtu":
        # Legacy: the PVTU references per-rank pieces.  Use a glob to
        # avoid parsing the PVTU itself (parser would also pick them up
        # but the glob is simpler and resilient against truncated files).
        cycle_glob = fault_dir.glob(
            f"fault_surface_r*_c{entry.cycle}.vtu")
        files.extend(sorted(cycle_glob))
    return files


# ---------------------------------------------------------------------------
# R-301 — per-rank → per-cycle merge (plan §Phase 5 step 1).
# ---------------------------------------------------------------------------

def _parse_ascii_dataarray_floats(elem) -> List[float]:
    if elem is None or elem.text is None:
        return []
    return [float(t) for t in elem.text.split()]


def _parse_ascii_dataarray_ints(elem) -> List[int]:
    if elem is None or elem.text is None:
        return []
    return [int(t) for t in elem.text.split()]


def _is_vtu_truncated(path: Path) -> bool:
    """R-310 / plan §Phase 5 Edge Cases: detect VTU files that were
    truncated mid-write (e.g. by a tar-on-login-node CPU kill).
    Heuristic: a complete VTU ends with `</VTKFile>` within the last
    256 bytes.  Returns True for missing or truncated files.
    """
    try:
        with path.open("rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - 256))
            tail = f.read().decode("utf-8", errors="replace")
    except OSError:
        return True
    return "</VTKFile>" not in tail


def _read_legacy_piece(path: Path):
    """Parse a per-rank legacy ASCII VTU piece.  Returns
    `(points_flat, triangles_flat, cell_data_by_name)` or None.
    """
    if _is_vtu_truncated(path):
        return None
    try:
        tree = ET.parse(str(path))
    except (ET.ParseError, OSError):
        return None
    root = tree.getroot()
    piece = root.find(".//Piece")
    if piece is None:
        return None
    pts_elem = piece.find("./Points/DataArray")
    if pts_elem is None or pts_elem.get("format", "ascii") != "ascii":
        return None
    points = _parse_ascii_dataarray_floats(pts_elem)
    conn_elem = piece.find("./Cells/DataArray[@Name='connectivity']")
    if conn_elem is None or conn_elem.get("format", "ascii") != "ascii":
        return None
    triangles = _parse_ascii_dataarray_ints(conn_elem)
    cell_data: Dict[str, List[float]] = {}
    cd_block = piece.find("./CellData")
    if cd_block is not None:
        for da in cd_block.findall("./DataArray"):
            name = da.get("Name")
            if name is None or da.get("format", "ascii") != "ascii":
                continue
            cell_data[name] = _parse_ascii_dataarray_floats(da)
    return points, triangles, cell_data


def merge_cycle_to_ascii_vtu(fault_dir: Path,
                             entry: CycleEntry,
                             dest_dir: Path) -> Optional[Path]:
    """Merge `fault_surface_r{r}_c{cycle}.vtu` per-rank pieces into one
    ASCII VTU at `dest_dir/fault_surface_c{cycle}.vtu`.  Mirrors Phase
    1's `GatheredFaultPack` layout: each cell owns three unique
    consecutive vertices, triangles[c] = {3*c, 3*c+1, 3*c+2} only when
    each piece itself follows that invariant.

    For Phase-1+ already-merged single-file VTUs (head suffix .vtu),
    returns the head path unchanged.  Returns None if no readable
    pieces exist.
    """
    pieces = collect_cycle_files(fault_dir, entry)
    if not pieces:
        return None
    head = pieces[0]
    if head.suffix != ".pvtu":
        return head if head.is_file() else None

    rank_pieces = pieces[1:]
    if not rank_pieces:
        return None
    merged_points: List[float] = []
    merged_triangles: List[int] = []
    merged_cell_data: Dict[str, List[float]] = {}
    field_order: List[str] = []
    n_total_cells = 0
    n_total_pts = 0
    for piece_path in rank_pieces:
        if not piece_path.is_file():
            continue
        parsed = _read_legacy_piece(piece_path)
        if parsed is None:
            print(f"warning: malformed VTU {piece_path}; skipping",
                  file=sys.stderr)
            continue
        points, triangles, cell_data = parsed
        if len(points) % 3 != 0 or len(triangles) % 3 != 0:
            print(f"warning: piece {piece_path} has non-triangle layout; "
                  f"skipping", file=sys.stderr)
            continue
        merged_points.extend(points)
        # Re-base triangle indices: piece-local index `i` becomes
        # global index `i + base_vertex` where base_vertex is the
        # number of vertices accumulated from prior pieces.
        base_vertex = n_total_pts // 3
        merged_triangles.extend(t + base_vertex for t in triangles)
        n_total_pts += len(points)
        n_total_cells += len(triangles) // 3
        for name, vals in cell_data.items():
            if name not in merged_cell_data:
                merged_cell_data[name] = []
                field_order.append(name)
            merged_cell_data[name].extend(vals)

    if n_total_cells == 0:
        return None

    out_path = dest_dir / f"fault_surface_c{entry.cycle}.vtu"
    with out_path.open("w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
        f.write('<UnstructuredGrid>\n')
        f.write(f'<Piece NumberOfPoints="{n_total_pts // 3}" '
                f'NumberOfCells="{n_total_cells}">\n')
        f.write('<Points><DataArray type="Float64" '
                'NumberOfComponents="3" format="ascii">\n')
        for i in range(0, len(merged_points), 3):
            f.write(f"{merged_points[i]} {merged_points[i+1]} "
                    f"{merged_points[i+2]}\n")
        f.write('</DataArray></Points>\n')
        f.write('<Cells>\n'
                '<DataArray type="Int32" Name="connectivity" '
                'format="ascii">\n')
        for i in range(0, len(merged_triangles), 3):
            f.write(f"{merged_triangles[i]} {merged_triangles[i+1]} "
                    f"{merged_triangles[i+2]}\n")
        f.write('</DataArray>\n')
        f.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
        for c in range(n_total_cells):
            f.write(f"{(c + 1) * 3}\n")
        f.write('</DataArray>\n')
        f.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
        for _ in range(n_total_cells):
            f.write("5\n")
        f.write('</DataArray>\n</Cells>\n<CellData>\n')
        for name in field_order:
            vals = merged_cell_data[name]
            if len(vals) != n_total_cells:
                if len(vals) < n_total_cells:
                    vals = vals + [0.0] * (n_total_cells - len(vals))
                else:
                    vals = vals[:n_total_cells]
            f.write(f'<DataArray type="Float64" Name="{name}" '
                    f'format="ascii">\n')
            for v in vals:
                f.write(f"{v}\n")
            f.write('</DataArray>\n')
        f.write('</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n')
    return out_path


# ---------------------------------------------------------------------------
# tar.zst chunk mode
# ---------------------------------------------------------------------------

def _detect_zstd() -> Optional[str]:
    return shutil.which("zstd") or None


def _have_zstandard_module() -> bool:
    try:
        import zstandard  # noqa: F401
        return True
    except ImportError:
        return False


def write_chunk_tar(out_path: Path,
                    fault_dir: Path,
                    entries: Sequence[CycleEntry]) -> int:
    """Write a single chunk tar.zst (or tar.gz fallback) for the given cycle range.

    Per plan §Phase 5 step 1, the chunk contains the MERGED per-cycle
    VTUs (one per cycle), NOT the raw per-rank pieces.  Phase-1+
    already-merged single-file VTUs are passed through unchanged.

    Returns the number of cycles successfully written into the chunk
    (the PVD is added once per chunk regardless).
    """
    if out_path.exists():
        return 0   # idempotent: skip already-present chunk

    use_zstd = _detect_zstd() is not None or _have_zstandard_module()
    suffix = ".tar.zst" if use_zstd else ".tar.gz"
    out_str = str(out_path)
    if not out_str.endswith(suffix):
        out_path = Path(out_str + suffix)
        if out_path.exists():
            return 0

    with tempfile.TemporaryDirectory(
            prefix=".pack_fault.", dir=str(out_path.parent)) as tmpdir:
        merge_dir = Path(tmpdir) / "merged"
        merge_dir.mkdir()

        # R-301: merge per-rank pieces into single per-cycle VTUs in
        # the temp dir, then tar only the merged outputs.
        files_to_add: List[Tuple[Path, str]] = []
        for entry in entries:
            merged = merge_cycle_to_ascii_vtu(fault_dir, entry, merge_dir)
            if merged is None:
                print(f"warning: cycle {entry.cycle}: no readable pieces; "
                      f"skipping", file=sys.stderr)
                continue
            # arcname uses the canonical Phase-1 single-file naming.
            files_to_add.append(
                (merged, f"fault_surface_c{entry.cycle}.vtu"))
        # Also include the PVD so ParaView can scrub the chunk's
        # contents directly.  Re-write a chunk-local PVD that points
        # at the merged filenames (independent of the source PVD).
        chunk_pvd = merge_dir / "fault_surface.pvd"
        with chunk_pvd.open("w") as f:
            f.write('<?xml version="1.0"?>\n'
                    '<VTKFile type="Collection" version="0.1">\n'
                    '<Collection>\n')
            for entry in entries:
                f.write(
                    f'<DataSet timestep="{entry.time:.17g}" '
                    f'file="fault_surface_c{entry.cycle}.vtu"/>\n')
            f.write('</Collection>\n</VTKFile>\n')
        files_to_add.append((chunk_pvd, "fault_surface.pvd"))

        if not files_to_add or len(files_to_add) <= 1:
            # R-208 follow-on: nothing meaningful to pack; skip
            # writing the chunk so re-runs can detect the absence
            # rather than see an empty tar that blocks regeneration.
            print(f"warning: chunk for cycles {entries[0].cycle}.."
                  f"{entries[-1].cycle} has no merged VTUs; not "
                  f"writing {out_path}", file=sys.stderr)
            return 0

        final_tmp = Path(tmpdir) / out_path.name
        if suffix == ".tar.zst" and _have_zstandard_module():
            import zstandard
            cctx = zstandard.ZstdCompressor(level=10, threads=-1)
            with open(final_tmp, "wb") as raw, \
                 cctx.stream_writer(raw) as zfh, \
                 tarfile.open(fileobj=zfh, mode="w|") as tar:
                for src, arcname in files_to_add:
                    tar.add(str(src), arcname=arcname)
        elif suffix == ".tar.zst":
            tar_tmp = Path(tmpdir) / (out_path.stem.replace(".tar", "")
                                      + ".tar")
            with tarfile.open(str(tar_tmp), "w") as tar:
                for src, arcname in files_to_add:
                    tar.add(str(src), arcname=arcname)
            subprocess.run(
                ["zstd", "-q", "-10", "-T0", "-f", "--rm",
                 "-o", str(final_tmp), str(tar_tmp)],
                check=True)
        else:
            with tarfile.open(str(final_tmp), "w:gz",
                              compresslevel=6) as tar:
                for src, arcname in files_to_add:
                    tar.add(str(src), arcname=arcname)
        os.replace(final_tmp, out_path)

    return len(files_to_add)


# ---------------------------------------------------------------------------
# VTKHDF mode
# ---------------------------------------------------------------------------

def write_hdf5(out_path: Path,
               fault_dir: Path,
               entries: Sequence[CycleEntry]) -> int:
    """Merge all per-rank pieces into one VTKHDF time-series.

    Requires `pyvista` (preferred) or `vtk` Python bindings.  Returns
    the number of timesteps written.

    Implementation strategy: read each cycle's piece(s) with pyvista
    `read(...)`, MergeMeshes if multiple, append to an HDF5 file using
    a transient ParaView-compatible layout under `/VTKHDF/Steps/`.
    Where pyvista's API does not expose a multi-step VTKHDF writer,
    we fall back to one VTKHDF per cycle in a sibling directory.

    Plan §Phase 5 acceptance: file < 2 GB, opens in ParaView identical
    to scrubbing the original PVD.  This implementation prefers a
    single file but degrades gracefully when the runtime lacks support.
    """
    try:
        import pyvista as pv
    except ImportError:
        try:
            import vtk  # noqa: F401
            print("warning: pyvista not installed; falling back to "
                  "raw vtk Python bindings (limited).  install pyvista "
                  "for full VTKHDF support.", file=sys.stderr)
            return _write_hdf5_via_vtk(out_path, fault_dir, entries)
        except ImportError as exc:
            raise RuntimeError(
                "VTKHDF mode requires pyvista or vtk Python bindings.  "
                "Install with `pip install pyvista` or "
                "`conda install -c conda-forge pyvista`.") from exc

    if out_path.exists():
        out_path.unlink()

    # pyvista 0.45+ supports HDFWriter; older versions write per-step
    # VTU and rely on a sibling .pvd.  We try the new path first.
    has_hdf_writer = hasattr(pv, "HDFWriter") or hasattr(pv, "save_meshio")
    n_steps = 0
    if has_hdf_writer and hasattr(pv, "HDFWriter"):
        writer = pv.HDFWriter(str(out_path))
        for entry in entries:
            mesh = _read_cycle(fault_dir, entry, pv)
            if mesh is None:
                continue
            mesh.field_data["time"] = [entry.time]
            writer.write(mesh, time=entry.time)
            n_steps += 1
        writer.close()
        return n_steps

    # Fallback: write per-cycle VTKHDF into a sibling directory and a
    # PVD index next to it (still smaller than the original per-rank
    # mess by a factor of nranks, even without a single file).
    print("warning: pyvista has no HDFWriter; emitting per-cycle "
          ".vtkhdf files in a sibling directory.", file=sys.stderr)
    sib_dir = out_path.parent / (out_path.stem + "_per_cycle")
    sib_dir.mkdir(exist_ok=True)
    pvd_root = ET.Element("VTKFile", {"type": "Collection",
                                      "version": "0.1"})
    coll = ET.SubElement(pvd_root, "Collection")
    for entry in entries:
        mesh = _read_cycle(fault_dir, entry, pv)
        if mesh is None:
            continue
        cycle_path = sib_dir / f"fault_surface_c{entry.cycle}.vtkhdf"
        try:
            mesh.save(str(cycle_path))
        except (RuntimeError, ValueError):
            cycle_path = cycle_path.with_suffix(".vtu")
            mesh.save(str(cycle_path))
        ET.SubElement(coll, "DataSet",
                      {"timestep": f"{entry.time:.17g}",
                       "file": cycle_path.name})
        n_steps += 1
    pvd_path = sib_dir / (out_path.stem + ".pvd")
    ET.ElementTree(pvd_root).write(str(pvd_path), encoding="utf-8",
                                    xml_declaration=True)
    print(f"wrote {n_steps} per-cycle VTKHDFs and PVD index at "
          f"{pvd_path}", file=sys.stderr)
    return n_steps


def _read_cycle(fault_dir: Path,
                entry: CycleEntry,
                pv) -> Optional["pv.UnstructuredGrid"]:
    pieces = collect_cycle_files(fault_dir, entry)
    if not pieces:
        return None
    head = pieces[0]
    if head.suffix == ".pvtu":
        # pyvista handles PVTU natively.
        try:
            return pv.read(str(head))
        except (RuntimeError, FileNotFoundError) as exc:
            print(f"warning: could not read {head} ({exc}); skipping "
                  f"cycle {entry.cycle}", file=sys.stderr)
            return None
    # Single-file VTU (Phase-1+ binary).
    try:
        return pv.read(str(head))
    except (RuntimeError, FileNotFoundError) as exc:
        print(f"warning: could not read {head} ({exc}); skipping "
              f"cycle {entry.cycle}", file=sys.stderr)
        return None


def _write_hdf5_via_vtk(out_path: Path,
                        fault_dir: Path,
                        entries: Sequence[CycleEntry]) -> int:
    # Bare-vtk fallback: per-cycle VTU → per-cycle VTKHDF using
    # vtkHDFWriter (VTK 9.3+).  Older VTK lacks this writer; in that
    # case we error out with an actionable message.
    import vtk
    if not hasattr(vtk, "vtkHDFWriter"):
        raise RuntimeError(
            "vtk version too old for vtkHDFWriter; install pyvista or "
            "upgrade vtk.")
    sib_dir = out_path.parent / (out_path.stem + "_per_cycle")
    sib_dir.mkdir(exist_ok=True)
    n_steps = 0
    for entry in entries:
        pieces = collect_cycle_files(fault_dir, entry)
        if not pieces:
            continue
        head = pieces[0]
        reader = (vtk.vtkXMLPUnstructuredGridReader()
                  if head.suffix == ".pvtu"
                  else vtk.vtkXMLUnstructuredGridReader())
        reader.SetFileName(str(head))
        reader.Update()
        writer = vtk.vtkHDFWriter()
        writer.SetFileName(
            str(sib_dir / f"fault_surface_c{entry.cycle}.vtkhdf"))
        writer.SetInputData(reader.GetOutput())
        writer.Write()
        n_steps += 1
    return n_steps


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Pack a legacy fault-surface VTU directory.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "fault_dir",
        type=Path,
        help="path to FaultSurface directory containing fault_surface.pvd "
             "(e.g. /scratch/.../job_xxx/FaultSurface)")
    out_group = parser.add_mutually_exclusive_group(required=True)
    out_group.add_argument(
        "--out-tar",
        type=Path,
        metavar="BASENAME",
        help="emit chunked tar.zst (or tar.gz fallback) under BASENAME.  "
             "Idempotent: skips existing chunks on re-run.")
    out_group.add_argument(
        "--out-hdf5",
        type=Path,
        metavar="PATH.vtkhdf",
        help="emit single VTKHDF file (requires pyvista).")
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100,
        help="cycles per tar chunk (default: 100)")
    parser.add_argument(
        "--max-cycles",
        type=int,
        default=0,
        help="stop after N cycles (0 = no limit; useful for testing).")
    args = parser.parse_args(argv)

    fault_dir = args.fault_dir.resolve()
    pvd_path = fault_dir / "fault_surface.pvd"
    entries = parse_pvd(pvd_path)
    if args.max_cycles > 0:
        entries = entries[: args.max_cycles]
    if not entries:
        print(f"error: no <DataSet> entries in {pvd_path}",
              file=sys.stderr)
        return 1
    print(f"found {len(entries)} cycles in {pvd_path}", file=sys.stderr)
    # R-311: log the on-disk layout so the user can tell whether the
    # merge step (R-301) had work to do.
    layout = detect_layout(fault_dir)
    print(f"detected layout: {layout}", file=sys.stderr)
    # R-204: ensure the output parent exists; otherwise the temp-dir
    # creation inside write_chunk_tar / write_hdf5 raises
    # FileNotFoundError mid-write with no actionable message.
    out_target = args.out_tar if args.out_tar is not None else args.out_hdf5
    out_target.parent.mkdir(parents=True, exist_ok=True)

    if args.out_tar is not None:
        chunk_size = max(1, args.chunk_size)
        n_chunks = (len(entries) + chunk_size - 1) // chunk_size
        total_files = 0
        for c in range(n_chunks):
            chunk_entries = entries[c * chunk_size:(c + 1) * chunk_size]
            chunk_path = args.out_tar.with_name(
                f"{args.out_tar.name}_chunk{c + 1:04d}")
            n = write_chunk_tar(chunk_path, fault_dir, chunk_entries)
            if n == 0:
                print(f"chunk {c + 1}/{n_chunks}: skipped (already "
                      f"exists)", file=sys.stderr)
            else:
                total_files += n
                print(f"chunk {c + 1}/{n_chunks}: {n} files",
                      file=sys.stderr)
        print(f"done: {total_files} files across {n_chunks} chunks",
              file=sys.stderr)
        return 0

    # out_hdf5
    n = write_hdf5(args.out_hdf5, fault_dir, entries)
    print(f"done: wrote {n} timesteps to {args.out_hdf5}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
