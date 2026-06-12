#!/usr/bin/env python3
"""generate_multifault_geo.py - Phase 3 generator for the SAFS preferred set.

Reads `data_corefined/manifest.json` and emits `code_meshing/
safs_multifault_box_<R>m.geo` that embeds all corefined faults inside a
single bounding-box volume for `gmsh -3`.

Implements PLAN_cgal_corefine_multifault.md Phase 3 (Phase 3.5 sub-phases
are NOT yet implemented; see plan for what to add later).

Usage:
    conda activate pythonenv
    python generate_multifault_geo.py --res 2000

Idempotency: re-running with the same manifest produces a byte-identical
.geo (no embedded timestamps; alphabetical fault ordering; fixed-precision
float formatting).  Verifiable via `md5 code_meshing/*.geo` before/after.
"""

import argparse
import json
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants from the plan template (lines 385-413).
PAD_XY_DEFAULT     = 50000.0   # 50 km horizontal margin
PAD_TOP_DEFAULT    =   100.0   # 100 m above z=0 (free-surface convention)
PAD_BOTTOM_DEFAULT = 25000.0   # 25 km below the deepest fault vertex
LC_FAR_MULTIPLIER  =    10.0   # default: --lc-far = 10 * lc-near
DIST_INNER_DEFAULT =  3000.0
DIST_OUTER_DEFAULT = 40000.0
GEOMETRY_TOLERANCE =     1e-3  # 1 mm; tighter than STL writer's LSB at UTM scale


# ---------------------------------------------------------------------------
# Manifest helpers.
REQUIRED_TOP_FIELDS = ("mesh_edge_size", "min_edge", "meshes")
REQUIRED_MESH_FIELDS = ("basename", "bbox", "n_faces")
REQUIRED_BBOX_FIELDS = ("x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi")


def load_manifest(path: Path) -> dict:
    """Load and validate a corefine manifest.  Exits 1 on schema violation."""
    if not path.is_file():
        print(f"error: manifest not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        m = json.load(f)
    for k in REQUIRED_TOP_FIELDS:
        if k not in m:
            print(f"error: manifest missing required field '{k}': {path}",
                  file=sys.stderr)
            sys.exit(1)
    if not isinstance(m["meshes"], list) or len(m["meshes"]) == 0:
        print(f"error: manifest 'meshes' must be a non-empty list: {path}",
              file=sys.stderr)
        sys.exit(1)
    for i, mesh in enumerate(m["meshes"]):
        for k in REQUIRED_MESH_FIELDS:
            if k not in mesh:
                print(f"error: manifest meshes[{i}] missing field '{k}'",
                      file=sys.stderr)
                sys.exit(1)
        for k in REQUIRED_BBOX_FIELDS:
            if k not in mesh["bbox"]:
                print(f"error: manifest meshes[{i}].bbox missing field '{k}'",
                      file=sys.stderr)
                sys.exit(1)
    return m


def union_bbox(meshes: list[dict]) -> dict:
    """Union of per-fault bboxes from manifest.meshes[]."""
    return {
        "x_lo": min(m["bbox"]["x_lo"] for m in meshes),
        "x_hi": max(m["bbox"]["x_hi"] for m in meshes),
        "y_lo": min(m["bbox"]["y_lo"] for m in meshes),
        "y_hi": max(m["bbox"]["y_hi"] for m in meshes),
        "z_lo": min(m["bbox"]["z_lo"] for m in meshes),
        "z_hi": max(m["bbox"]["z_hi"] for m in meshes),
    }


# ---------------------------------------------------------------------------
# Physical-surface name derivation (plan lines 449-456, 459-465).
#
# Rule: from a basename like
#     SAFS-SAFZ-COAV-Mission_Creek_fault_strand-CFM4_2000m
# strip the leading `SAFS-SAFZ-` and the trailing `_<R>m`, then replace
# every `-` with `_` to produce a gmsh-legal identifier:
#     COAV_Mission_Creek_fault_strand_CFM4
# and prefix with `fault_`:
#     fault_COAV_Mission_Creek_fault_strand_CFM4
def physical_surface_name(basename: str, res: int) -> str:
    name = basename
    prefix = "SAFS-SAFZ-"
    if name.startswith(prefix):
        name = name[len(prefix):]
    suffix = f"_{res}m"
    if name.endswith(suffix):
        name = name[: -len(suffix)]
    name = name.replace("-", "_")
    return f"fault_{name}"


# ---------------------------------------------------------------------------
# .geo emission.
def emit_geo(*,
             manifest: dict,
             res: int,
             pad_xy: float,
             pad_top: float,
             pad_bottom: float,
             lc_min: float,
             lc_near: float,
             lc_far: float,
             dist_inner: float,
             dist_outer: float) -> str:
    """Render the .geo file as a string (deterministic; no timestamps)."""
    meshes = sorted(manifest["meshes"], key=lambda m: m["basename"])
    n_faults = len(meshes)
    ub = union_bbox(meshes)

    # Sanity checks (plan §Edge Cases line 503).
    assert pad_top > 0, "pad_top must be > 0 (top is the free surface)"
    assert ub["z_hi"] + pad_top > ub["z_hi"], "pad_top did not raise zmax"
    assert pad_bottom > 0, "pad_bottom must be > 0"
    assert pad_xy > 0,     "pad_xy must be > 0"

    # Fault-tag layout (plan lines 419-429): box volume = 1; box face tags 1..6;
    # discrete-surface tags 7..6+n_faults in order of Merge.
    fault_tag_start = 7
    fault_tags = list(range(fault_tag_start, fault_tag_start + n_faults))

    out = []
    add = out.append

    add( "// =============================================================")
    add(f"// safs_multifault_box_{res}m.geo  (generated; do not edit by hand)")
    add( "//")
    add(f"// Embed the {n_faults} corefined SAFS faults inside a 3-D bounding box.")
    add(f"// Fault sources: ../data_corefined/*_{res}m_corefined.stl, produced by")
    add( "// code_preprocess/corefine_faults.py (CGAL corefine + cleanup;")
    add( "// see document/PLAN_cgal_corefine_multifault.md).")
    add( "//")
    add(f"// Generator: code_preprocess/generate_multifault_geo.py")
    add( "// Manifest:  ../data_corefined/manifest.json")
    add( "//")
    add( "// Generate mesh:")
    add(f"//   gmsh -3 safs_multifault_box_{res}m.geo -o safs_multifault_box_{res}m.msh")
    add( "//")
    add( "// All units: meters.  Coordinates: UTM.")
    add( "// =============================================================")
    add( "SetFactory(\"OpenCASCADE\");")
    add( "")
    add( "// 1. Domain bounds: union of corefined-fault bboxes + paddings.")
    add( "//    Filled in from manifest by the generator.")
    add(f"xmin_fault = {ub['x_lo']:.2f};  xmax_fault = {ub['x_hi']:.2f};")
    add(f"ymin_fault = {ub['y_lo']:.2f};  ymax_fault = {ub['y_hi']:.2f};")
    add(f"zmin_fault = {ub['z_lo']:.2f};  zmax_fault = {ub['z_hi']:.2f};   // top trace clipped at z=0")
    add( "")
    add(f"PAD_XY     = {pad_xy:.1f};   // {pad_xy/1000:.0f} km horizontal margin")
    add(f"PAD_TOP    = {pad_top:.1f};   // {pad_top:.0f} m of rock above the free surface")
    add(f"PAD_BOTTOM = {pad_bottom:.1f};   // {pad_bottom/1000:.0f} km of rock below deepest fault vertex")
    add( "")
    add( "xmin = xmin_fault - PAD_XY;")
    add( "xmax = xmax_fault + PAD_XY;")
    add( "ymin = ymin_fault - PAD_XY;")
    add( "ymax = ymax_fault + PAD_XY;")
    add( "zmin = zmin_fault - PAD_BOTTOM;")
    add( "zmax = zmax_fault + PAD_TOP;")
    add( "dx = xmax - xmin;  dy = ymax - ymin;  dz = zmax - zmin;")
    add( "")
    add( "// 2. Mesh size controls.  HARD floor: LC_MIN must equal Phase 1 --min-edge.")
    add( "//    LC_MIN, LC_NEAR, LC_FAR sourced from manifest:")
    add( "//      LC_MIN  = manifest.min_edge       (= --min-edge passed to corefine)")
    add( "//      LC_NEAR = manifest.mesh_edge_size (= --mesh-edge-size passed to corefine)")
    add(f"//      LC_FAR  = --lc-far  (default {LC_FAR_MULTIPLIER:.0f} * LC_NEAR)")
    add(f"LC_MIN     = {lc_min:.2f};")
    add(f"LC_NEAR    = {lc_near:.2f};")
    add(f"LC_FAR     = {lc_far:.2f};")
    add(f"DIST_INNER = {dist_inner:.2f};")
    add(f"DIST_OUTER = {dist_outer:.2f};")
    add( "")
    add( "Mesh.MeshSizeMin                = LC_MIN;")
    add( "Mesh.MeshSizeMax                = LC_FAR;")
    add( "Mesh.MeshSizeExtendFromBoundary = 0;")
    add( "Mesh.MeshSizeFromPoints         = 0;")
    add( "Mesh.MeshSizeFromCurvature      = 0;")
    add(f"Geometry.Tolerance              = {GEOMETRY_TOLERANCE};          // 1 mm; tighter than STL writer LSB at UTM scale")
    add( "")
    add( "// Disable gmsh's facet-overlap detector.  Default 0.1° flags any pair of")
    add( "// non-manifold incident triangles whose dihedral angle is below threshold.")
    add( "// Our combined STL has 54+ polyline edges with dihedral < 0.1° (fault pairs")
    add( "// meeting at near-coplanar angles — geological reality at the SBMT-MC × SBMT-")
    add( "// SAF intersection, etc.).  Setting to 0 lets the 3-D mesher proceed; the")
    add( "// resulting tets near these edges are necessarily thin (sliver-like) but")
    add( "// Mesh.OptimizeNetgen below should clean them up.")
    add( "Mesh.AngleToleranceFacetOverlap = 0;")
    add( "")
    add( "// 3. Box volume.  OCC tags 1..6 for box faces; volume tag 1.")
    add( "Box(1) = {xmin, ymin, zmin, dx, dy, dz};")
    add( "")
    add(f"// 4. Merge ONE combined STL containing the box (12 tris) + all {n_faults}")
    add( "//    fault triangulations as a single non-manifold discrete surface.")
    add( "//    Built by code_preprocess/stitch_combined_stl.py from data_corefined/.")
    add( "//    The single-Merge approach (vs per-fault Merges) is required because")
    add( "//    tetgen rejects multi-surface PLCs at shared polyline edges; with one")
    add( "//    discrete surface containing those edges as INTERNAL non-manifold")
    add( "//    edges, tetgen handles them via its standard non-manifold logic.")
    add(f"Merge \"../data_corefined/safs_combined_{res}m.stl\";")
    add( "")
    add(f"// gmsh assigns the merged STL discrete-surface tag {fault_tag_start}.  Embed it.")
    add(f"combined_surf = {fault_tag_start};")
    add( "Surface{combined_surf} In Volume{1};")
    add( "")
    add( "// 5. Distance-based size field on the polyline POINTS (PLAN_revert_to_gmsh.md")
    add( "//    §Phase 2 Plan B; original gmsh-plan §Variant 3.5a').  Auxiliary Point()")
    add( "//    entities are emitted at every corefined fault-fault polyline vertex; the")
    add( "//    Distance field measures distance to the nearest such point, then a")
    add( "//    Threshold ramps SizeMin=LC_NEAR at d=DIST_INNER to SizeMax=LC_FAR at")
    add( "//    d=DIST_OUTER.  Geometry.Tolerance = 1e-3 m guarantees these Points")
    add( "//    dedup against the discrete-STL polyline vertices (which are at %.15g")
    add( "//    precision in the input STL).")
    add( "POLY_PT_ID_START   = 100000;")
    add( "POLY_LINE_ID_START = 200000;")
    # Emit polyline points + lines from manifest.pairs[k].polylines.
    poly_pt_id   = 100000
    poly_line_id = 200000
    poly_pt_ids  = []
    poly_line_ids = []
    add( "// Polyline vertices and lines:")
    for pair in manifest.get("pairs", []):
        for poly in pair.get("polylines", []):
            prev_id = None
            for v in poly:
                add(f"Point({poly_pt_id}) = {{{v[0]:.6f}, {v[1]:.6f}, {v[2]:.6f}, LC_FAR}};")
                poly_pt_ids.append(poly_pt_id)
                if prev_id is not None:
                    add(f"Line({poly_line_id}) = {{{prev_id}, {poly_pt_id}}};")
                    poly_line_ids.append(poly_line_id)
                    poly_line_id += 1
                prev_id = poly_pt_id
                poly_pt_id += 1
    add( "")
    add(f"poly_pt_ids[]   = {{{', '.join(str(i) for i in poly_pt_ids)}}};")
    if poly_line_ids:
        add(f"poly_curve_ids[] = {{{', '.join(str(i) for i in poly_line_ids)}}};")
    else:
        add( "poly_curve_ids[] = {};")
    add( "")
    add( "// Distance field: distance to the nearest polyline vertex (PointsList).")
    add( "// Plan §Variant 3.5a' fallback; bypasses ambiguity of Distance over the")
    add( "// merged STL surface (which would otherwise pick up the box faces too).")
    add( "Field[1] = Distance;")
    add( "Field[1].PointsList = {poly_pt_ids[]};")
    add( "Field[1].Sampling   = 50;")
    add( "")
    add( "Field[2] = Threshold;")
    add( "Field[2].InField  = 1;")
    add( "Field[2].SizeMin  = LC_NEAR;")
    add( "Field[2].SizeMax  = LC_FAR;")
    add( "Field[2].DistMin  = DIST_INNER;")
    add( "Field[2].DistMax  = DIST_OUTER;")
    add( "Background Field = 2;")
    add( "")
    add( "// 6. Physical groups.  All 6 faults share a single 'fault_all' tag (100);")
    add( "//    downstream Python (msh_to_vtu.py) splits per-fault using the")
    add( "//    safs_combined_<R>m_provenance.json file.")
    add( "//    Box boundaries: top=200, bottom=201, sides=202.")
    add( "Physical Volume(\"rock\", 1)         = {1};")
    add(f"Physical Surface(\"fault_all\", 100) = {{combined_surf}};")
    add( "Physical Surface(\"top\",      200) = {6};")
    add( "Physical Surface(\"bottom\",   201) = {5};")
    add( "Physical Surface(\"sides\",    202) = {1, 2, 3, 4};")
    add( "")
    add( "// 7. Mesh algorithms (proven on single-fault freesurface_clip variant).")
    add( "Mesh.Algorithm   = 6;     // 2-D: Frontal-Delaunay")
    add( "Mesh.Algorithm3D = 1;     // 3-D: Delaunay (HXT/Algorithm3D=10 fails on")
    add( "                          // embedded discrete STLs)")
    add( "")
    add( "// 8. Phase 3.5b: post-mesh tet optimization (mandatory per original plan).")
    add( "//    Mesh.OptimizeNetgen typically raises min q_tet from ~0.05 to ~0.15.")
    add( "Mesh.Optimize          = 1;")
    add( "Mesh.OptimizeNetgen    = 1;")
    add( "Mesh.OptimizeThreshold = 0.3;")
    add( "Mesh.HighOrderOptimize = 0;")
    add( "")  # trailing newline
    return "\n".join(out)


# ---------------------------------------------------------------------------
# CLI.
def main() -> int:
    here = Path(__file__).resolve().parent
    project = here.parent
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path,
                    default=project / "data_corefined" / "manifest.json",
                    help="path to corefine manifest.json")
    ap.add_argument("--res", type=int, default=2000,
                    help="resolution suffix used in STL filenames and .geo name")
    ap.add_argument("--out-geo", type=Path, default=None,
                    help="output .geo path (default: "
                         "<project>/code_meshing/safs_multifault_box_<R>m.geo)")
    ap.add_argument("--pad-xy", type=float, default=PAD_XY_DEFAULT,
                    help=f"horizontal padding [m] (default: {PAD_XY_DEFAULT})")
    ap.add_argument("--pad-top", type=float, default=PAD_TOP_DEFAULT,
                    help=f"top padding above z=0 [m] (default: {PAD_TOP_DEFAULT})")
    ap.add_argument("--pad-bottom", type=float, default=PAD_BOTTOM_DEFAULT,
                    help=f"bottom padding below deepest fault [m] (default: {PAD_BOTTOM_DEFAULT})")
    ap.add_argument("--lc-near", type=float, default=None,
                    help="near-fault target edge size [m] (default: manifest.mesh_edge_size)")
    ap.add_argument("--lc-min", type=float, default=None,
                    help="hard min-edge floor [m] (default: manifest.min_edge)")
    ap.add_argument("--lc-far", type=float, default=None,
                    help=f"far-field target edge size [m] (default: {LC_FAR_MULTIPLIER:.0f} * lc-near)")
    ap.add_argument("--dist-inner", type=float, default=DIST_INNER_DEFAULT)
    ap.add_argument("--dist-outer", type=float, default=DIST_OUTER_DEFAULT)
    args = ap.parse_args()

    manifest = load_manifest(args.manifest)

    # Resolve defaults from manifest (R-008).
    lc_near = args.lc_near if args.lc_near is not None else float(manifest["mesh_edge_size"])
    lc_min  = args.lc_min  if args.lc_min  is not None else float(manifest["min_edge"])
    lc_far  = args.lc_far  if args.lc_far  is not None else (LC_FAR_MULTIPLIER * lc_near)

    if lc_min <= 0 or lc_near <= 0 or lc_far <= 0:
        print(f"error: lc_min={lc_min}, lc_near={lc_near}, lc_far={lc_far} "
              f"must all be positive", file=sys.stderr)
        return 1
    if not (lc_min <= lc_near <= lc_far):
        print(f"error: require lc_min <= lc_near <= lc_far; got "
              f"{lc_min} {lc_near} {lc_far}", file=sys.stderr)
        return 1

    out_geo = args.out_geo
    if out_geo is None:
        out_geo = project / "code_meshing" / f"safs_multifault_box_{args.res}m.geo"
    out_geo.parent.mkdir(parents=True, exist_ok=True)

    text = emit_geo(
        manifest=manifest,
        res=args.res,
        pad_xy=args.pad_xy,
        pad_top=args.pad_top,
        pad_bottom=args.pad_bottom,
        lc_min=lc_min,
        lc_near=lc_near,
        lc_far=lc_far,
        dist_inner=args.dist_inner,
        dist_outer=args.dist_outer,
    )
    out_geo.write_text(text)
    n_faults = len(manifest["meshes"])
    print(f"wrote {out_geo}  ({n_faults} fault(s), "
          f"LC_MIN={lc_min}, LC_NEAR={lc_near}, LC_FAR={lc_far})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
