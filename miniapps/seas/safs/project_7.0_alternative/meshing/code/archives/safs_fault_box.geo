// =============================================================
// safs_fault_box.geo
//
// Embed the San Andreas Fault System trace
//   (SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m)
// inside a rectangular 3-D volume.
//
// This file does NOT model a free surface.  The fault is buried
// strictly inside the box on all sides; the box top sits a comfortable
// distance ABOVE the highest fault vertex.  The user plans to revisit
// the raw .ts data and produce a true free-surface variant later.
//
// Pipeline:
//   1. The TSurf .ts file is clipped at z = 0 and (optionally) snapped
//      to z = -OFFSET_TOP by tools/ts_to_stl.py to give a clean
//      horizontal top edge.  Resulting STL sits next to this file.
//   2. Gmsh creates an OCC bounding box that fully encloses the fault.
//   3. The fault STL is merged as a discrete surface (gmsh tag 7,
//      i.e. one past the six box faces) and embedded as a 2-D
//      constraint inside the box volume via `Surface{} In Volume{}`.
//
// Note on `ClassifySurfaces`: the gmsh "discrete -> geometric" workflow
// (`ClassifySurfaces` + `CreateGeometry`) is intentionally NOT used.
// On this fault it produces multiple disconnected discrete patches whose
// shared edges are not stitched into a single closed manifold, and the
// 3-D mesher (HXT or Delaunay) reports "No closed volume".  Embedding
// the un-classified discrete surface directly works.
//
// Generate mesh:
//   gmsh -3 safs_fault_box.geo -o safs_fault_box.msh
//   gmsh -3 safs_fault_box.geo -clscale 4.0 -o safs_fault_box_coarse.msh
//
// All units: meters.  Coordinates: UTM (matches the .ts source).
// =============================================================

SetFactory("OpenCASCADE");

// -------------------------------------------------------------
// Fault bounding box (after z<=0 clipping & top snap; printed by
// tools/ts_to_stl.py --top-offset 50):
//   X  in [179308.97, 622329.47]      ( 443.02 km wide )
//   Y  in [3692278.31, 3988004.64]    ( 295.73 km tall )
//   Z  in [-16607.45, -50.00]         (  16.56 km deep, top at -50 m )
// -------------------------------------------------------------

xmin_fault =  179308.97;
xmax_fault =  622329.47;
ymin_fault = 3692278.31;
ymax_fault = 3988004.64;
zmin_fault = -16607.45;
zmax_fault =    -50.00;

// Domain padding -- edit to taste.
// The fault is buried inside the box on all six sides; PAD_TOP > 0
// places the box top strictly above the fault so the fault is
// guaranteed to be a fully interior 2-D constraint (gmsh requires
// embedded surfaces to be in the volume's interior, not on its
// boundary).  PAD_BOTTOM extends the box well below the deepest fault
// vertex.
PAD_XY     = 50000.0;   // 50 km horizontal margin
PAD_TOP    =  5000.0;   //  5 km of rock above the highest fault vertex
PAD_BOTTOM = 25000.0;   // 25 km of rock below the deepest fault vertex

xmin = xmin_fault - PAD_XY;
xmax = xmax_fault + PAD_XY;
ymin = ymin_fault - PAD_XY;
ymax = ymax_fault + PAD_XY;
zmin = zmin_fault - PAD_BOTTOM;
zmax = zmax_fault + PAD_TOP;

dx = xmax - xmin;
dy = ymax - ymin;
dz = zmax - zmin;

// -------------------------------------------------------------
// Mesh size controls -- edit to taste.
// -------------------------------------------------------------
LC_NEAR     =  1500.0;   // near-fault element edge
LC_FAR      = 10000.0;   // far-field element edge
DIST_INNER  =  3000.0;   // distance at which size starts to coarsen
DIST_OUTER  = 40000.0;   // distance beyond which size = LC_FAR

Mesh.CharacteristicLengthMin = LC_NEAR;
Mesh.CharacteristicLengthMax = LC_FAR;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;

// -------------------------------------------------------------
// Step 1: bounding box volume.
// OCC tags after Box(1) (deterministic in gmsh OCC):
//   volume  : 1
//   surfaces: 1 (xmin), 2 (xmax), 3 (ymin), 4 (ymax),
//             5 (zmin = bottom),  6 (zmax = top cap; NOT a free
//                                surface in this round -- the fault
//                                is buried entirely inside the box)
// -------------------------------------------------------------
Box(1) = {xmin, ymin, zmin, dx, dy, dz};

// -------------------------------------------------------------
// Step 2: merge the fault surface and embed it.
// Merge of an STL produces a single discrete surface whose tag is the
// next free dim-2 tag.  After Box(1) takes tags 1..6, the fault lands
// on tag 7.
// -------------------------------------------------------------
Merge "SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m.stl";

fault_tag = 7;
fault_surfs[] = {fault_tag};

// `Surface{} In Volume{}` forces the volume mesher to preserve the
// fault triangulation as interior facets.  We do NOT call
// `ClassifySurfaces`; see the file header for the reason.
Surface{fault_tag} In Volume{1};

// -------------------------------------------------------------
// Distance-based size field: refine near the fault.
// -------------------------------------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfs[]};
Field[1].Sampling     = 100;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = LC_NEAR;
Field[2].SizeMax  = LC_FAR;
Field[2].DistMin  = DIST_INNER;
Field[2].DistMax  = DIST_OUTER;

Background Field = 2;

// -------------------------------------------------------------
// Physical groups (so that downstream solvers can find tags).
//   - Volume  "rock"     : the box volume
//   - Surface "fault"    : the embedded discrete fault
//   - Surface "top"      : top face of the box (z = zmax) -- buried, NOT a free surface in this round
//   - Surface "bottom"   : bottom face (z = zmin)
//   - Surface "sides"    : the four lateral faces
// Box face tags are deterministic for a single OCC Box(1).
// -------------------------------------------------------------
Physical Volume("rock", 1)      = {1};
Physical Surface("fault", 101)  = fault_surfs[];
Physical Surface("top", 102)    = {6};
Physical Surface("bottom", 103) = {5};
Physical Surface("sides", 104)  = {1, 2, 3, 4};

// -------------------------------------------------------------
// Mesh algorithms.  Algorithm3D = 1 (Delaunay) handles embedded
// discrete surfaces robustly; HXT (Algorithm3D = 10) reports
// "No closed volume" on this geometry.
// -------------------------------------------------------------
Mesh.Algorithm   = 6;    // 2-D: Frontal-Delaunay
Mesh.Algorithm3D = 1;    // 3-D: Delaunay
