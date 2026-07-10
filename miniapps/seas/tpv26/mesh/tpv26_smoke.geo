// TPV26/27 Phase-0 SMOKE mesh — a TINY surface-breaking vertical
// strike-slip fault, for plumbing/smoke tests ONLY (NOT a production or
// convergence mesh; the real tpv26_{200,100,50}m.geo arrive in a later
// phase).  Structure mirrors tpv31/mesh/tpv31_50m.geo (a proven
// surface-reaching-fault template) but is shrunk to ~a few-km box so it
// meshes and runs locally in seconds.
//
// CANONICAL SEAS COORDINATE FRAME (see miniapps/seas/CLAUDE.md):
//   x = along-strike, y = fault-normal, z = vertical (z=0 surface, z<0 down).
//
// Geometry:
//   * Box:   x,y in [-3, 3] km,  z in [-3, 0] km.
//   * Fault: y = 0 plane,  x in [-1.5, 1.5] km,  z in [-1.5, 0] km
//            (surface-breaking: top edge at z = 0).
//
// Boundary attributes (match tpv26/configs/tpv26_spatial_smoke.toml
// [boundary]; tpv205 convention 103/101/105):
//   Physical Surface(101) = free surface (z = 0)            -> natural
//   Physical Surface(103) = fault (y = 0)
//   Physical Surface(105) = absorbing (four walls + bottom)
//
// Build:  gmsh -format msh22 -3 tpv26_smoke.geo -o tpv26_smoke.msh
// (CLAUDE.md "Known limitation — Gmsh .msh format": msh22 required.)

// ----------------------------------------------------------------------
// Mesh-size parameters (coarse — smoke only).
// ----------------------------------------------------------------------
lc_fault = 400.0;
lc_far   = 1200.0;

// ----------------------------------------------------------------------
// Domain box (m) — canonical SEAS frame.
// ----------------------------------------------------------------------
xMin = -3000.0;  xMax =  3000.0;
yMin = -3000.0;  yMax =  3000.0;
zMin = -3000.0;  zMax =     0.0;

// Fault rectangle (y = 0 plane), x in [-1.5,1.5] km, z in [-1.5,0] km.
fx0 = -1500.0;  fx1 =  1500.0;
fz0 = -1500.0;  fz1 =     0.0;

// ----------------------------------------------------------------------
// Box corners (looking from +y): 1-4 top (z=0), 5-8 bottom (z=zMin).
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

// Fault rectangle edges.
Line(21) = {11, 12};  Line(22) = {12, 13};
Line(23) = {13, 14};  Line(24) = {14, 11};

Curve Loop(101) = {21, 22, 23, 24};
Plane Surface(1) = {101};                 // Fault rectangle

// Box faces.
Curve Loop(102) = {1, 2, 3, 4};      Plane Surface(2) = {102};  // z=0 free surface
Curve Loop(103) = {5, 6, 7, 8};      Plane Surface(3) = {103};  // z=zMin bottom
Curve Loop(104) = {4, 9, -8, -12};   Plane Surface(4) = {104};  // x=xMin wall
Curve Loop(105) = {2, 11, -6, -10};  Plane Surface(5) = {105};  // x=xMax wall
Curve Loop(106) = {1, 10, -5, -9};   Plane Surface(6) = {106};  // y=yMin wall
Curve Loop(107) = {3, 12, -7, -11};  Plane Surface(7) = {107};  // y=yMax wall

Surface Loop(1) = {2, 3, 4, 5, 6, 7};
Volume(1) = {1};

// Embed the fault's surface-breaking top edge (Line 23, y=0,z=0) in the
// free surface so the 3D mesher conforms (mirrors tpv31_50m.geo:173).
Curve{23} In Surface{2};
Surface{1} In Volume{1};

// ----------------------------------------------------------------------
// Mesh-size field (tpv31-style fast coarsening away from the fault).
// ----------------------------------------------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {1};

Field[2] = MathEval;
Field[2].F = Sprintf("0.3*F1 + (F1/1.0e3)^2 + %g", lc_fault);

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
// Physical groups (match [boundary] in the smoke config).
// ----------------------------------------------------------------------
Physical Surface(101) = {2};             // free surface (z = 0)
Physical Surface(103) = {1};             // fault (y = 0)
Physical Surface(105) = {3, 4, 5, 6, 7}; // absorbing (bottom + four walls)
Physical Volume(1)   = {1};

// ----------------------------------------------------------------------
// Generator hints (msh22 required — see CLAUDE.md).
// ----------------------------------------------------------------------
Mesh.MshFileVersion  = 2.2;
Mesh.Algorithm3D     = 1;     // Delaunay
// Do NOT enable Mesh.OptimizeNetgen (SIGABRT on surface-reaching faults;
// see tpv31_50m.geo).
