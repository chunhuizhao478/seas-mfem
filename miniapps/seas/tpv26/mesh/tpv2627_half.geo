// SCEC TPV26 / TPV27 mesh, y >= 0 half, for a mesh that is mirror-symmetric about the fault
// plane (rheology plan rev 1, phase R1; Chunhui Zhao, 2026-10-20).  scripts/mirror_mesh.py
// reflects this half across Y = 0 (tpv26/mesh/make_tpv2627_meshes.sh does both):
//
//   gmsh -3 tpv26/mesh/tpv2627_half.geo -setnumber lc_fault 200 -format msh22 -o half.msh
//   python3 scripts/mirror_mesh.py half.msh tpv26/mesh/tpv2627_mirror_200m.msh --fault-attr 103
//
// Why the mirror: on a gmsh mesh the tets on the two sides of a fault face differ and the
// discrete fault normal stress then changes with slip (TPV205 100 m, Vista 1076419: 38-231 MPa
// against 120 MPa exact); the mirror mesh keeps the change zero for an elastic bulk
// (scripts/mirror_mesh.py, commit 41879e8).  TPV26 is elastic; TPV27 runs on the same mesh.
//
// Geometry (TPV26/27 description v13, Parts 2-3), code frame X along strike, Y fault normal
// (SCEC z, the far side at Y > 0), Z up (SCEC -y): fault X in [-20, 20] km, Z in [-20, 0] km on
// Y = 0, reaching the free surface; hypocenter (-5, 0, -10) km.  Half-space box X in
// [-Xmax, Xmax], Y in [0, Ymax] (mirrored to [-Ymax, Ymax]), Z in [Zmin, 0]; the margins
// (45 km along strike, 40 km below the fault, 60 km across) follow tpv2053d_half.geo.  The
// rupture stops at the fault's buried edges as the benchmark prescribes (no fault faces past
// them).
//
// Size field: min(lc, 0.05 d + (d / 2.5 km)^2 + lc_fault) with d the distance to the fault
// (tpv2053d_half.geo), and inside the box |Y| <= box_halfwidth, X within box_halfwidth of the
// fault ends, depth up to fault depth + box_halfwidth, at most box_factor * lc_fault: the
// off-fault plastic zone of TPV27.  box_halfwidth = 0 drops the box.
//
// Physical groups (TPV26/27 decks): surface 101 free surface (Z = 0), 103 fault, 105 absorbing,
// 99 the rest of the Y = 0 plane (interior after mirroring; mirror_mesh.py drops it);
// volume 1.
//
// Built with gmsh's built-in geometry kernel only (no OpenCASCADE: the Vista gmsh module has
// none; Vista build job 1085864 got an empty mesh from the first, OpenCASCADE version, as the
// TPV6/TPV31 .geo files did in d910b2f).  The Y = 0 face of the half-box is split into the fault
// and the rest of the plane, which share the fault's two sides and bottom edge; the fault's top
// edge lies in the free surface.

DefineConstant[ lc_fault = 200 ];        // element size on the fault (m)
DefineConstant[ lc = 5e3 ];              // largest element size (m)
DefineConstant[ box_halfwidth = 2e3 ];   // plastic-zone box half-width (m); 0 = no box
DefineConstant[ box_factor = 2 ];        // element size in the box, in units of lc_fault
Fault_length = 40e3;
Fault_width  = 20e3;
Xmax = 65e3;  Ymax = 60e3;  Zmin = -60e3;
xf = 0.5 * Fault_length;
e = 1;

// The Y = 0 face: its corners and the fault corners (the fault's top edge is in Z = 0).
Point(1) = {-Xmax, 0, 0, lc};
Point(2) = {-xf, 0, 0, lc};
Point(3) = {xf, 0, 0, lc};
Point(4) = {Xmax, 0, 0, lc};
Point(5) = {Xmax, 0, Zmin, lc};
Point(6) = {-Xmax, 0, Zmin, lc};
Point(7) = {-xf, 0, -Fault_width, lc};
Point(8) = {xf, 0, -Fault_width, lc};
Point(111) = {-5e3, 0, 0, lc};   // the two surface stations (below)
Point(112) = {15e3, 0, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 111};              // the fault's top edge, in three pieces
Line(31) = {111, 112};
Line(32) = {112, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 7};
Line(8) = {7, 8};
Line(9) = {8, 3};
Curve Loop(1) = {1, 7, 8, 9, 3, 4, 5, 6};
Curve Loop(2) = {2, 31, 32, -9, -8, -7};
Plane Surface(1) = {1};       // the Y = 0 plane outside the fault (99)
Plane Surface(2) = {2};       // the fault (103)

// The 12 on-fault stations (description v13, Part 8; tpv26/stations/tpv26_fault_stations.csv) are
// mesh vertices, so every station file holds the solution at its SCEC point.  Without them the
// nearest fault node was up to 50 m deeper and 46 m along strike (job 1094998, 200 m: sigma_n 0.3 to
// 1 % above DRDG3D at five stations).  The two surface stations, points 111 and 112, end pieces of
// the fault's top edge (the built-in kernel cannot embed a point in a curve).
Point(101) = {-15e3, 0, -10e3, lc};
Point(102) = {-5e3, 0, -5e3, lc};
Point(103) = {-5e3, 0, -10e3, lc};
Point(104) = {-5e3, 0, -15e3, lc};
Point(105) = {0, 0, -10e3, lc};
Point(106) = {5e3, 0, -10e3, lc};
Point(107) = {10e3, 0, -10e3, lc};
Point(108) = {15e3, 0, -5e3, lc};
Point(109) = {15e3, 0, -10e3, lc};
Point(110) = {15e3, 0, -15e3, lc};
Point{101:110} In Surface{2};

// The far side Y = Ymax and the faces that join it to the Y = 0 face.
Point(11) = {-Xmax, Ymax, 0, lc};
Point(14) = {Xmax, Ymax, 0, lc};
Point(15) = {Xmax, Ymax, Zmin, lc};
Point(16) = {-Xmax, Ymax, Zmin, lc};
Line(11) = {1, 11};
Line(14) = {4, 14};
Line(15) = {5, 15};
Line(16) = {6, 16};
Line(21) = {11, 14};
Line(22) = {14, 15};
Line(23) = {15, 16};
Line(24) = {16, 11};
Curve Loop(11) = {21, 22, 23, 24};          // the far side Y = Ymax
Curve Loop(12) = {1, 2, 31, 32, 3, 14, -21, -11};   // the free surface Z = 0
Curve Loop(13) = {-5, 15, 23, -16};         // the bottom Z = Zmin
Curve Loop(14) = {4, 15, -22, -14};         // the side X = Xmax
Curve Loop(15) = {6, 11, -24, -16};         // the side X = -Xmax
For k In {11:15}
   Plane Surface(k) = {k};
EndFor
Surface Loop(1) = {1, 2, 11, 12, 13, 14, 15};
Volume(1) = {1};

Field[1] = Distance;
Field[1].SurfacesList = {2};
Field[1].Sampling = 300;
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);
Field[3] = Threshold;
Field[3].InField = 1;
Field[3].SizeMin = lc_fault;
Field[3].SizeMax = lc;
Field[3].DistMin = 2 * lc_fault;
Field[3].DistMax = 2 * lc_fault + 0.001;
If (box_halfwidth > 0)
   Field[5] = Box;
   Field[5].VIn = box_factor * lc_fault;
   Field[5].VOut = lc;
   Field[5].XMin = -0.5 * Fault_length - box_halfwidth;
   Field[5].XMax = 0.5 * Fault_length + box_halfwidth;
   Field[5].YMin = -box_halfwidth;
   Field[5].YMax = box_halfwidth;
   Field[5].ZMin = -Fault_width - box_halfwidth;
   Field[5].ZMax = e;
   Field[4] = Min;
   Field[4].FieldsList = {2, 3, 5};
Else
   Field[4] = Min;
   Field[4].FieldsList = {2, 3};
EndIf
Background Field = 4;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromCurvature = 0;

Physical Surface(101) = {12};
Physical Surface(103) = {2};
Physical Surface(105) = {11, 13, 14, 15};
Physical Surface(99) = {1};
Physical Volume(1) = {1};
Mesh.MshFileVersion = 2.2;
