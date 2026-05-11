// =============================================================
// safs_fault_box_buried.geo
//
// Full-space variant: bury the entire (un-clipped, un-snapped) fault
// triangulation inside a 3-D box so the box top is well above the
// fault's highest vertex (z = +2058 m) and the box bottom is well
// below its deepest vertex (z = -16607 m).  The fault is therefore a
// strictly interior 2-D constraint on all six sides of the box.
//
// This is the "does the snap cause the slivers?" test: by skipping the
// z = 0 clip and the top-edge snap to z = -50 m, the fault's source
// triangulation enters the model intact (median triangle quality
// ~0.999, no triangles with q < 0.3, shortest edge 1262 m) and any
// remaining bad tets in the volume mesh are purely a meshing artefact,
// not an input artefact.
//
// Generate mesh:
//   gmsh -3 safs_fault_box_buried.geo -o safs_fault_box_buried.msh
//
// All units: meters.  Coordinates: UTM (matches the .ts source).
// =============================================================

SetFactory("OpenCASCADE");

// -------------------------------------------------------------
// Un-clipped fault bounding box (from ../preprocess_data/...
// _unclipped.stl):
//   X  in [179308.97, 622329.47]      ( 443.02 km wide )
//   Y  in [3692278.31, 3988061.92]    ( 295.78 km tall )
//   Z  in [-16607.45,    2058.15]     (  18.67 km deep, top above z=0 )
// -------------------------------------------------------------

xmin_fault =  179308.97;
xmax_fault =  622329.47;
ymin_fault = 3692278.31;
ymax_fault = 3988061.92;
zmin_fault = -16607.45;
zmax_fault =   2058.15;

// Domain padding -- edit to taste.
PAD_XY     = 50000.0;   // 50 km horizontal margin
PAD_TOP    =  5000.0;   //  5 km of rock above the fault's highest vertex
PAD_BOTTOM = 25000.0;   // 25 km of rock below the fault's deepest vertex

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
// Mesh size controls -- match the clipped variant for direct comparison.
// -------------------------------------------------------------
LC_NEAR     =  1500.0;
LC_FAR      = 10000.0;
DIST_INNER  =  3000.0;
DIST_OUTER  = 40000.0;

Mesh.CharacteristicLengthMin = LC_NEAR;
Mesh.CharacteristicLengthMax = LC_FAR;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;

// -------------------------------------------------------------
// Box volume.  OCC tags 1..6 for box faces; volume tag 1.
// -------------------------------------------------------------
Box(1) = {xmin, ymin, zmin, dx, dy, dz};

// -------------------------------------------------------------
// Merge the un-clipped fault and embed it.  Discrete fault surface
// arrives at tag 7 (one past the six box faces).  Do NOT call
// ClassifySurfaces (it splits the fault into multiple discrete patches
// whose shared edges are not stitched into a closed manifold, so the
// 3-D mesher reports "No closed volume").
// -------------------------------------------------------------
Merge "../preprocess_data/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m_unclipped.stl";

fault_tag = 7;
fault_surfs[] = {fault_tag};

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
// Physical groups.
// -------------------------------------------------------------
Physical Volume("rock", 1)      = {1};
Physical Surface("fault", 101)  = fault_surfs[];
Physical Surface("top", 102)    = {6};
Physical Surface("bottom", 103) = {5};
Physical Surface("sides", 104)  = {1, 2, 3, 4};

// -------------------------------------------------------------
// Mesh algorithms (Algorithm3D = 1 / Delaunay; HXT fails on this
// embedded discrete surface).
// -------------------------------------------------------------
Mesh.Algorithm   = 6;
Mesh.Algorithm3D = 1;
