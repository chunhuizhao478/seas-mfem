#!/usr/bin/env python3
# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory. All rights
# reserved.  See files LICENSE and NOTICE for details.  LLNL-CODE-806117.
#
# Regression tests for miniapps/seas/scripts/pack_fault_output.py.
# Covers:
#   - R-301: per-rank pieces are MERGED into one VTU per cycle.
#   - R-310: truncated VTU pieces are detected and skipped with a
#            warning rather than packed verbatim.
#   - R-311: legacy / phase1 / mixed layouts are detected and reported.
#   - R-204: missing output parent dirs are created before packing.
#
# Run:    python3 test_pack_fault_output.py
# Exits 0 on success, 1 on failure.

from __future__ import annotations

import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "pack_fault_output.py"


def _write_legacy_piece(path: Path, value: float) -> None:
    path.write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" version="0.1">\n'
        '<UnstructuredGrid><Piece NumberOfPoints="3" NumberOfCells="1">\n'
        '<Points><DataArray type="Float64" NumberOfComponents="3" '
        'format="ascii">\n'
        '0 0 0 1 0 0 0 1 0\n'
        '</DataArray></Points>\n'
        '<Cells>\n'
        '<DataArray type="Int32" Name="connectivity" format="ascii">'
        '0 1 2</DataArray>\n'
        '<DataArray type="Int32" Name="offsets" format="ascii">3</DataArray>\n'
        '<DataArray type="UInt8" Name="types" format="ascii">5</DataArray>\n'
        '</Cells>\n<CellData>\n'
        f'<DataArray type="Float64" Name="slip" format="ascii">{value}'
        '</DataArray>\n'
        '</CellData></Piece></UnstructuredGrid></VTKFile>\n'
    )


def _make_legacy_dir(root: Path, n_cycles: int, n_ranks: int) -> Path:
    fault = root / "FaultSurface"
    fault.mkdir(parents=True, exist_ok=True)
    for c in range(n_cycles):
        for r in range(n_ranks):
            _write_legacy_piece(
                fault / f"fault_surface_r{r}_c{c}.vtu",
                value=float(f"{c}.{r}"))
        (fault / f"fault_surface_c{c}.pvtu").write_text(
            '<?xml version="1.0"?>\n'
            '<VTKFile type="PUnstructuredGrid" version="0.1">\n'
            '<PUnstructuredGrid GhostLevel="0">\n'
            '<PPoints><PDataArray type="Float64" '
            'NumberOfComponents="3"/></PPoints>\n'
            '<PCellData><PDataArray type="Float64" Name="slip"/>'
            '</PCellData>\n'
            + ''.join(
                f'<Piece Source="fault_surface_r{r}_c{c}.vtu"/>\n'
                for r in range(n_ranks))
            + '</PUnstructuredGrid></VTKFile>\n')
    pvd_lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1">',
        '<Collection>',
    ]
    for c in range(n_cycles):
        pvd_lines.append(
            f'<DataSet timestep="{float(c)}" '
            f'file="fault_surface_c{c}.pvtu"/>')
    pvd_lines.append('</Collection>')
    pvd_lines.append('</VTKFile>')
    (fault / "fault_surface.pvd").write_text("\n".join(pvd_lines) + "\n")
    return fault


def _list_chunk_members(chunk: Path) -> list:
    if chunk.suffix == ".zst":
        # Use the system zstd CLI; it's always available on Frontera
        # and on the local mfem-dev conda env.
        proc = subprocess.run(
            ["zstd", "-dc", str(chunk)],
            check=True, capture_output=True)
        with tarfile.open(fileobj=io.BytesIO(proc.stdout)) as tar:
            return sorted(m.name for m in tar.getmembers())
    with tarfile.open(str(chunk)) as tar:
        return sorted(m.name for m in tar.getmembers())


def _read_chunk_member(chunk: Path, member: str) -> bytes:
    if chunk.suffix == ".zst":
        proc = subprocess.run(
            ["zstd", "-dc", str(chunk)],
            check=True, capture_output=True)
        with tarfile.open(fileobj=io.BytesIO(proc.stdout)) as tar:
            return tar.extractfile(member).read()
    with tarfile.open(str(chunk)) as tar:
        return tar.extractfile(member).read()


class PackFaultOutputTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="pack_fault_test."))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_R301_chunk_holds_one_merged_vtu_per_cycle(self):
        """Plan §Phase 5 step 1 / R-301: per-rank pieces are MERGED."""
        fault = _make_legacy_dir(self.tmp / "src", n_cycles=2, n_ranks=4)
        out_base = self.tmp / "out" / "job"
        rc = subprocess.call(
            [sys.executable, str(SCRIPT), str(fault),
             "--out-tar", str(out_base), "--chunk-size", "10"])
        self.assertEqual(rc, 0)
        chunk = list(out_base.parent.glob("job_chunk0001*"))
        self.assertEqual(len(chunk), 1, "exactly one chunk expected")
        members = _list_chunk_members(chunk[0])
        # Expect: one merged VTU per cycle + one chunk-local PVD.
        self.assertEqual(
            members,
            ["fault_surface.pvd",
             "fault_surface_c0.vtu", "fault_surface_c1.vtu"],
            "R-301: chunk must contain merged per-cycle VTUs, NOT raw "
            "per-rank pieces")
        # The merged VTU has 4 cells (one per rank).
        merged = _read_chunk_member(chunk[0], "fault_surface_c0.vtu")
        root = ET.fromstring(merged)
        piece = root.find(".//Piece")
        self.assertIsNotNone(piece)
        self.assertEqual(int(piece.get("NumberOfCells")), 4,
                         "R-301: merged cycle has nranks cells")
        self.assertEqual(int(piece.get("NumberOfPoints")), 12,
                         "R-301: merged cycle has 3 * nranks vertices")

    def test_R310_truncated_piece_skipped_with_warning(self):
        """Plan §Phase 5 Edge Cases / R-310: malformed VTU is skipped."""
        fault = _make_legacy_dir(self.tmp / "src", n_cycles=1, n_ranks=4)
        # Truncate one rank's VTU mid-tag.
        (fault / "fault_surface_r2_c0.vtu").write_text("<?xml version=")
        out_base = self.tmp / "out" / "job"
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), str(fault),
             "--out-tar", str(out_base)],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("malformed VTU", proc.stderr,
                      "R-310: truncated VTU emits 'malformed VTU' warning")
        chunk = list(out_base.parent.glob("job_chunk0001*"))
        merged = _read_chunk_member(chunk[0], "fault_surface_c0.vtu")
        root = ET.fromstring(merged)
        piece = root.find(".//Piece")
        # Three good ranks remain; the truncated rank is dropped.
        self.assertEqual(int(piece.get("NumberOfCells")), 3,
                         "R-310: truncated rank dropped from merge")

    def test_R311_layout_detection_legacy(self):
        """Plan §Phase 5 Edge Cases / R-311: layout reporting."""
        fault = _make_legacy_dir(self.tmp / "src", n_cycles=1, n_ranks=2)
        out_base = self.tmp / "out" / "job"
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), str(fault),
             "--out-tar", str(out_base)],
            capture_output=True, text=True)
        self.assertIn("detected layout: legacy", proc.stderr,
                      "R-311: 'legacy' layout label printed for per-rank dir")

    def test_R204_creates_missing_output_parent(self):
        """R-204: missing parent directory is created before packing."""
        fault = _make_legacy_dir(self.tmp / "src", n_cycles=1, n_ranks=2)
        # Three levels deep — must be auto-created.
        target = self.tmp / "deep" / "nested" / "missing" / "job"
        rc = subprocess.call(
            [sys.executable, str(SCRIPT), str(fault),
             "--out-tar", str(target)])
        self.assertEqual(rc, 0)
        self.assertTrue(target.parent.is_dir(),
                        "R-204: parent dir auto-created")
        chunks = list(target.parent.glob("job_chunk0001*"))
        self.assertEqual(len(chunks), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
