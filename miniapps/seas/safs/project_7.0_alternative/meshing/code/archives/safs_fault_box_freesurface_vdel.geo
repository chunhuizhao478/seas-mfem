// =============================================================
// safs_fault_box_freesurface.geo
//
// Free-surface variant: embed the SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6
// fault inside a 3-D box whose top face approximates a free surface.
//
// Source fault: ../data_cleanfreesurf/<filename>.stl, produced by
// code_preprocess/clean_freesurface_mesh.py.  That pipeline drops the
// z > 0 cap of the raw fault, snaps the jagged top trace to z = 0,
// and isotropically remeshes; the surviving fault top edge sits exactly
// at z = 0.
//
// Box top is placed PAD_TOP above z = 0 (default 100 m) so the fault
// remains a strictly interior 2-D constraint and the proven
// `Surface{} In Volume{}` embedding works without ClassifySurfaces.
// PAD_TOP = 100 is a thin overburden — close enough to a free surface
// for SEAS post-processing while keeping gmsh's topology valid.  Set
// PAD_TOP = 0 to make the fault top edge coincident with the box top
// face boundary (gmsh may report "embedded surface touches volume
// boundary"; only use this if you have followed up with a Boolean
// fragment to stitch the fault trace into the top face).
//
// Generate mesh:
//   gmsh -3 safs_fault_box_freesurface.geo -o safs_fault_box_freesurface.msh
//
// All units: meters.  Coordinates: UTM (matches the .ts source).
// =============================================================

SetFactory("OpenCASCADE");

// -------------------------------------------------------------
// Cleaned-fault bounding box (from
// ../data_cleanfreesurf/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m_clean.stl):
//   X  in [179309.00, 622329.47]      ( 443.02 km wide )
//   Y  in [3692278.19, 3987901.78]    ( 295.62 km tall )
//   Z  in [-16607.45,        0.00]    (  16.61 km deep, top at z = 0 )
// -------------------------------------------------------------

xmin_fault =  179309.00;
xmax_fault =  622329.47;
ymin_fault = 3692278.19;
ymax_fault = 3987901.78;
zmin_fault = -16607.45;
zmax_fault =       0.00;

// Domain padding -- edit to taste.
PAD_XY     = 50000.0;   // 50 km horizontal margin
PAD_TOP    =   100.0;   // 100 m of rock above the free surface
                        // (thin overburden -- keeps fault strictly
                        //  interior so `Surface In Volume` is valid)
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
// Mesh size controls -- match the existing buried/clipped variants
// for direct comparability of mesh statistics.
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
// Merge the cleaned fault and embed it.  Discrete fault surface
// arrives at tag 7 (one past the six box faces).  Do NOT call
// ClassifySurfaces (it splits the fault into multiple discrete
// patches whose shared edges are not stitched into a closed
// manifold, so the 3-D mesher reports "No closed volume").
// -------------------------------------------------------------
Merge "../data_cleanfreesurf/SAFS-SAFZ-MULT-San_Andreas_fault_Fuis-ALT6_2000m_clean_vdel.stl";

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
//   - Volume  "rock"     : the box volume
//   - Surface "fault"    : the embedded discrete fault
//   - Surface "top"      : top face of the box (z = zmax) -- the
//                          approximate free surface
//   - Surface "bottom"   : bottom face (z = zmin)
//   - Surface "sides"    : the four lateral faces
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
