// SCEC TPV6/TPV7 mesh — half-space, vertical strike-slip fault reaching the
// free surface.  100 m near-fault resolution (spec).  STANDARD (asymmetric)
// mesh for Arm 2 (matrix + MIXED flux + RK).  Adapted from tpv31/mesh/
// tpv31_50m.geo (same canonical frame + surface-reaching-fault embed trick);
// the ONLY changes are lc_fault (100 m), the Physical-Surface tags (to match
// tpv6/configs/tpv6.toml [boundary]: free=101, fault=103, absorbing=105), and
// this header.
//
// The mesh carries NO material tags: TPV6/7's bi-material split is assigned by
// the config's [material.halfspace_across_fault] per element centroid sign of y
// (B3), so the SAME mesh serves both TPV6 (high contrast) and TPV7 (low).
//
// Canonical SEAS frame:  x = along-strike, y = fault-normal, z = depth (z = 0
// free surface, z < 0 below).  Fault: y = 0 plane, x in [-15,15] km, z in
// [-15,0] km (reaches the surface).  ref_normal = (0,-1,0), up = (0,0,1).
//
// Reflection-free window: cp_max = 6000 m/s, tfinal ~ 8 s -> d >= cp*t/2 = 24 km.
// The 50 km half-box clears this with margin (matches the proven TPV31 box).
//
// Generate (Frontera; .msh is gitignored, do NOT commit the .msh):
//   gmsh tpv6_100m.geo -3 -format msh22 -o tpv6_100m.msh
// (msh22 REQUIRED — the MFEM reader is Gmsh v2.2-only; see CLAUDE.md.)

// ----------------------------------------------------------------------
// Mesh-size parameters
// ----------------------------------------------------------------------
lc_fault = 100.0;           // 100 m near the fault (TPV6/7 spec)
lc_far   = 5000.0;          // coarse far field at the absorbing boundary

// ----------------------------------------------------------------------
// Domain box (x, y, z) in metres — canonical SEAS frame.
// ----------------------------------------------------------------------
xMin = -50000.0;  xMax =  50000.0;
yMin = -50000.0;  yMax =  50000.0;
zMin = -50000.0;  zMax =      0.0;   // z = 0 free surface, z < 0 below

// ----------------------------------------------------------------------
// Fault rectangle (y = 0 plane), x in [-15,15] km, z in [-15,0] km.
// ----------------------------------------------------------------------
fx0 = -15000.0;  fx1 =  15000.0;
fz0 = -15000.0;  fz1 =      0.0;     // depth 15 km up to the free surface

// ----------------------------------------------------------------------
// Vertices: 8 box corners (looking from +y):
//   z = zMax (= 0, free surface)   1(-x,-y) 2(+x,-y) / 4(-x,+y) 3(+x,+y)
//   z = zMin (= -50, bottom)       5(-x,-y) 6(+x,-y) / 8(-x,+y) 7(+x,+y)
// ----------------------------------------------------------------------
Point(1) = {xMin, yMin, zMax, lc_far};
Point(2) = {xMax, yMin, zMax, lc_far};
Point(3) = {xMax, yMax, zMax, lc_far};
Point(4) = {xMin, yMax, zMax, lc_far};
Point(5) = {xMin, yMin, zMin, lc_far};
Point(6) = {xMax, yMin, zMin, lc_far};
Point(7) = {xMax, yMax, zMin, lc_far};
Point(8) = {xMin, yMax, zMin, lc_far};

// Fault rectangle corners on the y = 0 plane.
Point(11) = {fx0, 0.0, fz0, lc_fault};
Point(12) = {fx1, 0.0, fz0, lc_fault};
Point(13) = {fx1, 0.0, fz1, lc_fault};
Point(14) = {fx0, 0.0, fz1, lc_fault};

// ----------------------------------------------------------------------
// Box edges.
// ----------------------------------------------------------------------
Line(1) = {1, 2};  Line(2) = {2, 3};  Line(3) = {3, 4};  Line(4) = {4, 1};
Line(5) = {5, 6};  Line(6) = {6, 7};  Line(7) = {7, 8};  Line(8) = {8, 5};
Line(9)  = {1, 5}; Line(10) = {2, 6}; Line(11) = {3, 7}; Line(12) = {4, 8};

// Fault rectangle edges (y = 0 plane).
Line(21) = {11, 12};  Line(22) = {12, 13};  Line(23) = {13, 14};  Line(24) = {14, 11};

Curve Loop(101) = {21, 22, 23, 24};
Plane Surface(1) = {101};      // Fault rectangle

// ----------------------------------------------------------------------
// Bounding-box faces.
// ----------------------------------------------------------------------
Curve Loop(102) = {1, 2, 3, 4};        Plane Surface(2) = {102};  // z=0 free surface
Curve Loop(103) = {5, 6, 7, 8};        Plane Surface(3) = {103};  // z=zMin bottom
Curve Loop(104) = {4, 9, -8, -12};     Plane Surface(4) = {104};  // x=xMin wall
Curve Loop(105) = {2, 11, -6, -10};    Plane Surface(5) = {105};  // x=xMax wall
Curve Loop(106) = {1, 10, -5, -9};     Plane Surface(6) = {106};  // y=yMin wall
Curve Loop(107) = {3, 12, -7, -11};    Plane Surface(7) = {107};  // y=yMax wall

// ----------------------------------------------------------------------
// Volume — embed the fault surface so the mesher splits tets along it.
// ----------------------------------------------------------------------
Surface Loop(1) = {2, 3, 4, 5, 6, 7};
Volume(1) = {1};

// Fault reaches the free surface: embed its top edge (Line 23, y=0,z=0) in the
// free-surface plane (Surface 2) so the 2D mesh conforms to the fault trace.
Curve{23} In Surface{2};
Surface{1} In Volume{1};

// ----------------------------------------------------------------------
// Mesh-size field — fast coarsening away from the fault (tpv205/tpv31 style).
// ----------------------------------------------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {1};

Field[2] = MathEval;
Field[2].F = Sprintf("0.3*F1 + (F1/2.5e3)^2 + %g", lc_fault);

Field[6] = Threshold;
Field[6].IField  = 1;
Field[6].LcMin   = lc_fault;
Field[6].LcMax   = lc_far;
Field[6].DistMin = 2*lc_fault;
Field[6].DistMax = 2*lc_fault + 0.001;

Field[7] = Min;
Field[7].FieldsList = {2, 6};
Background Field = 7;

Characteristic Length{ PointsOf{ Surface{1}; } } = lc_fault;

// ----------------------------------------------------------------------
// Physical groups — MUST match tpv6/configs/tpv6.toml [boundary]:
//   free=101, fault=103, absorbing=105.
// ----------------------------------------------------------------------
Physical Surface(101) = {2};               // free surface (z = 0)
Physical Surface(103) = {1};               // fault (y = 0)
Physical Surface(105) = {3, 4, 5, 6, 7};   // absorbing (bottom + 4 walls)
Physical Volume(1)   = {1};

// ----------------------------------------------------------------------
// Generator hints (msh22 required — see CLAUDE.md).
// ----------------------------------------------------------------------
Mesh.MshFileVersion = 2.2;
Mesh.Algorithm3D    = 1;     // Delaunay
// Do NOT enable Mesh.OptimizeNetgen (SIGABRT on surface-reaching faults; see
// tpv31_50m.geo).  The standard optimizer is sufficient.
