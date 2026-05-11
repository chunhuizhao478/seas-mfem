# Phase 7.1 of miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md.
#
# Minimal Gmsh `.msh` v2 / v4 reader for the size estimator.  Returns
# a small `MeshSummary` dataclass; no MFEM dependency.
#
# Lazy parsing — only `$Nodes` and `$Elements` blocks are read.

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional


@dataclass
class MeshSummary:
    n_elements: int = 0          # total volume elements (tets/hexes/...)
    n_vertices: int = 0
    n_fault_faces: int = 0       # surface elements with the fault tag
    n_boundary_faces: int = 0    # surface elements with a non-fault tag
    element_type: str = ""       # "tet" / "hex" / "prism" / "pyramid"
    element_order: int = 1       # geometric order (1 for linear)
    fault_tag: int = 3
    # Diagnostics — populated for `--explain` mode.
    extra: Dict[str, int] = field(default_factory=dict)


# Gmsh element type IDs we care about.  Per the Gmsh manual:
#   1 = 2-node line, 2 = 3-node triangle, 3 = 4-node quad,
#   4 = 4-node tet, 5 = 8-node hex, 6 = 6-node prism, 7 = 5-node pyramid,
#   8..14 = quadratic / cubic variants (treat as the same element_type
#                                       but bump element_order).
# Surface elements (fault / boundary): triangle (2) and quad (3).
_VOLUME_TYPES = {
    4: ("tet", 1),
    5: ("hex", 1),
    6: ("prism", 1),
    7: ("pyramid", 1),
    11: ("tet", 2),  # 10-node tet
    17: ("hex", 2),  # 27-node hex (less common)
}
_SURFACE_TYPES = {2, 3, 9, 10, 16}  # triangle/quad + their higher-order variants


def read_inline_mesh(nx: int, ny: int, nz: int,
                     element_type: str = "tet",
                     fault_axis: str = "y") -> MeshSummary:
    """Inline mesh summary for drivers that pass `--inline-mesh` (BP5
    has this flag).  No file I/O; matches MFEM's `Mesh::MakeCartesian3D`
    refinement.  For TET the cube is split into 6 tets; for HEX it's
    left as 1 hex per cube.

    `fault_axis` selects the cube face the fault lies on: 'x' / 'y' /
    'z'.  BP5 inline uses y=0 (default).  Fault faces are estimated as
    the count of surface elements on that single face plane:
      tet:  2 × dim_a × dim_b   (each quad face → 2 triangles)
      hex:  dim_a × dim_b       (one quad per cube face)
    Other 5 cube faces become `n_boundary_faces`."""
    if element_type == "tet":
        n_elements = nx * ny * nz * 6
    elif element_type == "hex":
        n_elements = nx * ny * nz
    else:
        raise ValueError(
            f"unsupported inline element type {element_type!r}; "
            f"choose 'tet' or 'hex'")
    n_vertices = (nx + 1) * (ny + 1) * (nz + 1)

    if fault_axis == "x":
        fault_a, fault_b = ny, nz
    elif fault_axis == "z":
        fault_a, fault_b = nx, ny
    else:  # default y
        fault_a, fault_b = nx, nz

    if element_type == "tet":
        n_fault_faces = 2 * fault_a * fault_b
        # 6 cube faces = 1 fault + 5 other; sum of all surface faces.
        n_total_surface = 2 * 2 * (nx * ny + ny * nz + nx * nz)
    else:  # hex
        n_fault_faces = fault_a * fault_b
        n_total_surface = 2 * (nx * ny + ny * nz + nx * nz)
    n_boundary_faces = n_total_surface - n_fault_faces

    return MeshSummary(
        n_elements=n_elements,
        n_vertices=n_vertices,
        n_fault_faces=n_fault_faces,
        n_boundary_faces=n_boundary_faces,
        element_type=element_type,
        element_order=1,
        fault_tag=0,
    )


def _detect_format(first_line: str) -> str:
    """Return 'gmsh-ascii' or raise."""
    if first_line.strip() == "$MeshFormat":
        return "gmsh-ascii"
    raise ValueError(
        f"unsupported mesh format: first line is {first_line!r}.  "
        f"This estimator supports Gmsh ASCII v2 and v4 only "
        f"(produced by `gmsh -3 ... -o foo.msh -format msh`).")


def _read_format_block(it):
    """Consume the $MeshFormat ... $EndMeshFormat block; return
    (version_major, is_binary)."""
    line = next(it).strip()
    parts = line.split()
    if len(parts) < 3:
        raise ValueError(f"malformed $MeshFormat header: {line!r}")
    version = parts[0]
    is_binary = parts[1] == "1"
    if is_binary:
        raise ValueError(
            "binary Gmsh meshes are not supported by this estimator.  "
            "Re-export with `gmsh -3 -format msh22` or `-format msh4` "
            "(both ASCII).")
    # Skip until $EndMeshFormat.
    for line in it:
        if line.strip() == "$EndMeshFormat":
            break
    return int(version.split(".")[0]), is_binary


def _read_v2_nodes(it) -> int:
    """v2 $Nodes block: first line is the node count, then one line per
    node `id x y z`.  Returns node count."""
    n = int(next(it).strip())
    for _ in range(n):
        next(it)
    end = next(it).strip()
    if end != "$EndNodes":
        raise ValueError(
            f"expected $EndNodes after {n} nodes; got {end!r}")
    return n


def _read_v4_nodes(it) -> int:
    """v4 $Nodes block.  Header: numEntityBlocks numNodes minTag maxTag.
    Then per-block: entityDim entityTag parametric numNodesInBlock,
    followed by numNodesInBlock IDs and numNodesInBlock coordinate
    triples (each on its own line).  Returns total node count."""
    header = next(it).split()
    n_total = int(header[1])
    consumed = 0
    while consumed < n_total:
        block = next(it).split()
        # entityDim entityTag parametric numNodesInBlock
        n_block = int(block[3])
        # Tags (n_block lines).
        for _ in range(n_block):
            next(it)
        # Coordinates (n_block lines).
        for _ in range(n_block):
            next(it)
        consumed += n_block
    end = next(it).strip()
    if end != "$EndNodes":
        raise ValueError(
            f"expected $EndNodes after {n_total} nodes; got {end!r}")
    return n_total


def _read_v2_elements(it, summary: MeshSummary) -> None:
    """v2 $Elements block.  Per-line:
    `id type n_tags tag1 tag2 [...] node1 node2 ...`.  The first tag is
    the physical-region tag we care about."""
    n = int(next(it).strip())
    for _ in range(n):
        parts = next(it).split()
        if len(parts) < 4:
            continue
        etype = int(parts[1])
        n_tags = int(parts[2])
        tag = int(parts[3]) if n_tags >= 1 else 0
        if etype in _VOLUME_TYPES:
            kind, order = _VOLUME_TYPES[etype]
            summary.n_elements += 1
            if not summary.element_type:
                summary.element_type = kind
                summary.element_order = order
        elif etype in _SURFACE_TYPES:
            summary.n_boundary_faces += 1
            if tag == summary.fault_tag:
                summary.n_fault_faces += 1
    end = next(it).strip()
    if end != "$EndElements":
        raise ValueError(
            f"expected $EndElements after {n} elements; got {end!r}")


def _read_v4_elements(it, summary: MeshSummary) -> None:
    """v4 $Elements block.  Header: numEntityBlocks numElements minTag
    maxTag.  Per block: entityDim entityTag elementType numElementsInBlock,
    followed by numElementsInBlock lines of `id node1 node2 ...`.

    The entity tag in the block header IS the physical-region tag in v4
    (when elementary == physical, which is the Gmsh default for
    geo-generated meshes).  More precisely the entity tag identifies
    the elementary entity; mapping to the physical group requires
    parsing $Entities, which is more involved.  For the BP5 mesh
    (`gmsh -3 bp5/mesh/bp5.geo`) the physical-surface tags equal the
    elementary-surface tags by convention, so this approximation is
    correct for the in-scope BP5 / TPV mesh families."""
    header = next(it).split()
    n_total = int(header[1])
    consumed = 0
    while consumed < n_total:
        block = next(it).split()
        # entityDim entityTag elementType numElementsInBlock
        if len(block) < 4:
            raise ValueError(f"malformed v4 element block header: {block}")
        ent_tag = int(block[1])
        etype = int(block[2])
        n_block = int(block[3])
        for _ in range(n_block):
            next(it)
        if etype in _VOLUME_TYPES:
            kind, order = _VOLUME_TYPES[etype]
            summary.n_elements += n_block
            if not summary.element_type:
                summary.element_type = kind
                summary.element_order = order
        elif etype in _SURFACE_TYPES:
            summary.n_boundary_faces += n_block
            if ent_tag == summary.fault_tag:
                summary.n_fault_faces += n_block
        consumed += n_block
    end = next(it).strip()
    if end != "$EndElements":
        raise ValueError(
            f"expected $EndElements after {n_total} elements; got {end!r}")


def _read_physical_names(it) -> Dict[int, str]:
    """Parse a `$PhysicalNames ... $EndPhysicalNames` block (Gmsh v2
    and v4 share this block format).  Returns the *surface* entries
    (dim == 2) as a `{tag: name}` map; volumes and other dims are
    skipped.

    The block format is:

        $PhysicalNames
        N
        2  100  "fault"
        2    5  "absorbing"
        3   10  "bulk"
        $EndPhysicalNames
    """
    surface_names: Dict[int, str] = {}
    n = int(next(it).strip())
    for _ in range(n):
        parts = next(it).split(None, 2)
        if len(parts) < 3:
            continue
        dim = int(parts[0])
        tag = int(parts[1])
        # Strip surrounding double-quotes from the name.
        name = parts[2].strip().strip('"')
        if dim == 2:
            surface_names[tag] = name
    end = next(it).strip()
    if end != "$EndPhysicalNames":
        raise ValueError(
            f"expected $EndPhysicalNames; got {end!r}")
    return surface_names


_FAULT_PREFIX_RE = re.compile(r"^fault(?:[^A-Za-z]|$)", re.IGNORECASE)


def _resolve_fault_tag(surface_names: Dict[int, str],
                       fallback: int) -> int:
    """Pick the surface tag whose name matches the fault.  Returns
    `fallback` if no match is found.

    Match rule (R-501): two-pass case-insensitive match.
      1. Exact equality with "fault" — preferred.
      2. Name starts with "fault" followed by a non-letter or
         end-of-string (matches "fault_main", "Fault-Zone", "FAULT 1";
         does NOT match "default", "fault_main"-as-a-suffix, "faulty",
         or "non_fault_zone").

    Plain substring matching was rejected because "default" contains
    the substring "fault" (positions 2-6 are f-a-u-l-t) and would
    silently mis-route catch-all Gmsh names into the fault count.
    Python's `\\b` boundary treats underscore as a word character, so
    `\\bfault\\b` would also reject "fault_main", which is the most
    common multi-fault naming convention; the start-anchored match
    above accepts it while still rejecting "default" and "non_fault_*".
    """
    # Pass 1: exact match.
    for tag, name in surface_names.items():
        if name.lower() == "fault":
            return tag
    # Pass 2: starts-with-fault followed by non-letter or EOS.
    for tag, name in surface_names.items():
        if _FAULT_PREFIX_RE.match(name):
            return tag
    return fallback


def read_gmsh(path: Path, fault_tag: int = 3) -> MeshSummary:
    """Read a Gmsh ASCII v2 or v4 mesh and return a `MeshSummary`.

    Lazy: only `$Nodes`, `$Elements`, and `$PhysicalNames` are
    consumed; everything else is skipped.  Aborts on binary or
    non-Gmsh formats.

    `fault_tag` is a *fallback* used only when the mesh has no
    `$PhysicalNames` block or none of its surface tags is named
    "fault".  When the mesh has a `Physical Surface(...) "fault"`
    declaration, that tag wins regardless of the `fault_tag` arg
    (R-402: prevents 180× under-counts on meshes that don't follow
    the Tandem tag-3 convention, e.g. `bp5_1000m.msh` uses tag 100).
    The resolved tag is stored in `summary.fault_tag`.
    """
    path = Path(path).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"mesh not found: {path}")
    # Pass 1: locate $PhysicalNames if present, to auto-detect the
    # fault tag before we count surface elements.
    physical_names: Dict[int, str] = {}
    with path.open("r") as fh:
        first = next(fh)
        _detect_format(first)
        _read_format_block(fh)
        for line in fh:
            line = line.rstrip("\n")
            if line == "$PhysicalNames":
                physical_names = _read_physical_names(fh)
                break
            if line == "$Nodes" or line == "$Elements":
                break  # PhysicalNames must come before these blocks.
    resolved_fault_tag = _resolve_fault_tag(physical_names, fault_tag)

    summary = MeshSummary(fault_tag=resolved_fault_tag)
    summary.extra["source_file"] = str(path)
    if physical_names:
        summary.extra["physical_names"] = ",".join(
            f"{t}={n}" for t, n in physical_names.items())
    with path.open("r") as fh:
        first = next(fh)
        _detect_format(first)
        # Consume $MeshFormat header.
        version_major, _binary = _read_format_block(fh)
        # Walk the file block by block.
        for line in fh:
            line = line.rstrip("\n")
            if line == "$Nodes":
                if version_major >= 4:
                    summary.n_vertices = _read_v4_nodes(fh)
                else:
                    summary.n_vertices = _read_v2_nodes(fh)
            elif line == "$Elements":
                if version_major >= 4:
                    _read_v4_elements(fh, summary)
                else:
                    _read_v2_elements(fh, summary)
            elif line.startswith("$") and not line.startswith("$End"):
                # Unknown block — skip to its $End<Name>.
                tag = line[1:]
                end_tag = "$End" + tag
                for skip in fh:
                    if skip.strip() == end_tag:
                        break
    if summary.n_boundary_faces == 0 and summary.n_fault_faces == 0:
        # No Physical Surface tags / no surface elements at all.
        # Document via a stderr warning per plan §Phase 7.1 Edge Cases.
        print(
            f"warning: mesh {path} has no surface elements; "
            f"n_fault_faces == n_boundary_faces == 0",
            file=sys.stderr)
    return summary
