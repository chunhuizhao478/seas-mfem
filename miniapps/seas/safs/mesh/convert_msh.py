"""Convert a Gmsh .msh to viewable XML formats.

Two outputs are produced for the same input mesh:

  <stem>.vtu   — VTK XML (UnstructuredGrid).  Open in ParaView directly.
                  Carries the gmsh:physical tag as a cell-data array
                  named "physical" so you can colour by tag.
  <stem>.xml   — Legacy Dolfin XML (FEniCS).  Volume-only — boundary tags
                  are written to a side-car <stem>_facet_region.xml.

The split is required because the .msh contains BOTH triangles (boundary
+ fault) AND tets (volume).  ParaView is happy to ingest both in a single
.vtu; FEniCS-legacy expects volume cells in <stem>.xml and boundary cells
in <stem>_facet_region.xml.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import numpy as np


def _split_msh(msh_path: Path):
    """Read the .msh and return per-cell-type arrays + tag arrays."""
    import meshio
    m = meshio.read(msh_path)
    points = np.asarray(m.points, dtype=np.float64)
    gp = m.cell_data_dict.get("gmsh:physical", {})

    tris_blocks = [cb for cb in m.cells if cb.type == "triangle"]
    tets_blocks = [cb for cb in m.cells if cb.type == "tetra"]

    if not tets_blocks:
        raise SystemExit(f"{msh_path}: no tetrahedra (volume cells)")

    tris = (np.concatenate([cb.data for cb in tris_blocks], axis=0)
            if tris_blocks else np.empty((0, 3), dtype=np.int64))
    tets = np.concatenate([cb.data for cb in tets_blocks], axis=0)

    if "triangle" in gp and tris_blocks:
        tri_tags = np.asarray(gp["triangle"], dtype=np.int64)
    else:
        tri_tags = np.full(tris.shape[0], 0, dtype=np.int64)
    if "tetra" in gp and tets_blocks:
        tet_tags = np.asarray(gp["tetra"], dtype=np.int64)
    else:
        tet_tags = np.full(tets.shape[0], 0, dtype=np.int64)

    return points, tris, tri_tags, tets, tet_tags


def _write_vtu(out: Path, points, tris, tri_tags, tets, tet_tags) -> None:
    """ParaView-friendly VTU: tris + tets in one mesh, tag in 'physical'."""
    import meshio
    cells = []
    cell_data = {"physical": []}
    if tris.shape[0]:
        cells.append(("triangle", tris))
        cell_data["physical"].append(tri_tags)
    if tets.shape[0]:
        cells.append(("tetra", tets))
        cell_data["physical"].append(tet_tags)
    mesh = meshio.Mesh(points=points, cells=cells, cell_data=cell_data)
    out.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(out, mesh)


def _write_dolfin_xml(out: Path, points, tris, tri_tags,
                        tets, tet_tags) -> None:
    """Volume-only Dolfin XML + companion facet region XML.

    Both files are valid against the Dolfin XML 1.0 spec for a 3-D
    tetrahedral mesh.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    facet_path = out.with_name(out.stem + "_facet_region.xml")

    # ---------- volume mesh ----------
    n_vertices = points.shape[0]
    n_cells = tets.shape[0]
    with open(out, "w") as fh:
        fh.write('<?xml version="1.0"?>\n')
        fh.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
        fh.write('  <mesh celltype="tetrahedron" dim="3">\n')
        fh.write(f'    <vertices size="{n_vertices}">\n')
        for i, (x, y, z) in enumerate(points):
            fh.write(f'      <vertex index="{i}" x="{x:.17g}" '
                     f'y="{y:.17g}" z="{z:.17g}"/>\n')
        fh.write('    </vertices>\n')
        fh.write(f'    <cells size="{n_cells}">\n')
        for i, tet in enumerate(tets):
            fh.write(
                f'      <tetrahedron index="{i}" '
                f'v0="{int(tet[0])}" v1="{int(tet[1])}" '
                f'v2="{int(tet[2])}" v3="{int(tet[3])}"/>\n'
            )
        fh.write('    </cells>\n')
        # Volume domain markers.
        fh.write('    <domains>\n')
        fh.write(f'      <mesh_value_collection name="cell_domains" '
                 f'type="uint" dim="3" size="{n_cells}">\n')
        for i, tag in enumerate(tet_tags):
            fh.write(
                f'        <value cell_index="{i}" local_entity="0" '
                f'value="{int(tag)}"/>\n'
            )
        fh.write('      </mesh_value_collection>\n')
        fh.write('    </domains>\n')
        fh.write('  </mesh>\n')
        fh.write('</dolfin>\n')

    # ---------- facet (triangle) tags ----------
    # Build face → cell+local_face_index map so we can name each tagged
    # triangle by its containing tet + local-facet index (0..3).
    # In Dolfin tetrahedra: local facet i is the face opposite vertex i.
    tet_face_map: dict[frozenset, tuple[int, int]] = {}
    for ci, tet in enumerate(tets):
        v = [int(x) for x in tet]
        # facet local index = vertex index opposite to it
        for li in range(4):
            verts = tuple(v[j] for j in range(4) if j != li)
            tet_face_map[frozenset(verts)] = (ci, li)

    with open(facet_path, "w") as fh:
        fh.write('<?xml version="1.0"?>\n')
        fh.write('<dolfin xmlns:dolfin="http://fenicsproject.org">\n')
        fh.write('  <mesh_value_collection name="facet_region" '
                 'type="uint" dim="2" size="{}">\n'.format(int(tris.shape[0])))
        skipped = 0
        for tri, tag in zip(tris, tri_tags):
            key = frozenset(int(x) for x in tri)
            hit = tet_face_map.get(key)
            if hit is None:
                skipped += 1
                continue
            ci, li = hit
            fh.write(
                f'    <value cell_index="{ci}" local_entity="{li}" '
                f'value="{int(tag)}"/>\n'
            )
        fh.write('  </mesh_value_collection>\n')
        fh.write('</dolfin>\n')
        if skipped:
            sys.stderr.write(
                f"warning: {skipped} triangle(s) not adjacent to any tet "
                f"were dropped from {facet_path.name} (typical for "
                f"orphan boundary tris); volume mesh is unaffected\n"
            )


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert SAFS .msh to viewable XML formats "
                    "(.vtu for ParaView, .xml for Dolfin/FEniCS-legacy).",
    )
    parser.add_argument("--msh", required=True, type=Path)
    parser.add_argument("--vtu", default=None, type=Path,
                        help="Output .vtu path; defaults to "
                             "<msh stem>.vtu next to the .msh.")
    parser.add_argument("--xml", default=None, type=Path,
                        help="Output Dolfin .xml path; defaults to "
                             "<msh stem>.xml next to the .msh.")
    parser.add_argument("--no-vtu", action="store_true",
                        help="Skip the VTU output.")
    parser.add_argument("--no-xml", action="store_true",
                        help="Skip the Dolfin XML output.")
    args = parser.parse_args(argv)

    if not args.msh.exists():
        raise SystemExit(f"missing input mesh: {args.msh}")

    points, tris, tri_tags, tets, tet_tags = _split_msh(args.msh)

    vtu_path = args.vtu or args.msh.with_suffix(".vtu")
    xml_path = args.xml or args.msh.with_suffix(".xml")

    written: list[Path] = []
    if not args.no_vtu:
        _write_vtu(vtu_path, points, tris, tri_tags, tets, tet_tags)
        written.append(vtu_path)
    if not args.no_xml:
        _write_dolfin_xml(xml_path, points, tris, tri_tags,
                            tets, tet_tags)
        written.append(xml_path)
        written.append(xml_path.with_name(
            xml_path.stem + "_facet_region.xml"))

    print(f"convert_msh OK: {args.msh}")
    print(f"  n_points   : {points.shape[0]}")
    print(f"  n_triangles: {tris.shape[0]}")
    print(f"  n_tetra    : {tets.shape[0]}")
    for w in written:
        print(f"  wrote      : {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
