// Map-view rate-and-state test mesh (2D in-plane, map view; Chunhui Zhao,
// 2026-10-29: "there is no plane strain map view benchmark, so it can
// assume the same layout as tpv205-2d-like (but with rate and state) and
// run FD").
//
// The TPV205-2D layout (tpv205/mesh/tpv205_2d.geo): a straight fault in a
// 2-D whole space, plane strain, the slip in the plane (map view, the
// decks' up = [0, 0, 1]), no free surface.  Rate and state instead of slip
// weakening, and the outer sides loaded as in the UCERF 2D model.
//
// Geometry (m):
//   fault: y = 0, x in [-15, 15] km, embedded in the box with two free tips
//          (outside it the material is continuous); split at x = +-10 and
//          +-12 km, where the decks' a changes (velocity weakening for
//          |x| <= 10 km, a linear transition to velocity strengthening at
//          12 km), so the friction steps fall on mesh vertices;
//   box:   x, y in [-60, 60] km, as TPV205-2D.
//
// Physical groups (the fd_cycle layout, cycle_config.hpp: fault 3,
// Dirichlet 5, natural 1):
//   Physical Curve 3   = fault
//   Physical Curve 5   = the four sides: the far-field plate loading in the
//                        QD regime ([cycle.loading]), absorbing walls in the
//                        dynamic one ([cycle.absorbing] faces = "all")
//   Physical Surface 1 = bulk
// No free surface, so no curve carries tag 1.
//
// Size: lc_fault on the fault, growing away from it as TPV205-2D's
// (0.05 d + (d / 2.5 km)^2 + lc_fault, d the distance to the fault), capped
// at lc.
//
// Units: meters (the decks' mesh_scale = 1).
//
// Usage:
//   gmsh -2 mapview_rs.geo -format msh22 -o mapview_rs_500m.msh
//   gmsh -2 mapview_rs.geo -format msh22 -setnumber lc_fault 250 -o mapview_rs_250m.msh

DefineConstant[ lc_fault = {500, Min 1, Max 5000, Name "Fault resolution [m]" } ];
DefineConstant[ lc = {5e3, Min 100, Max 50e3, Name "Far-field resolution [m]" } ];

Fault_half_length = 15e3;
X_vw = 10e3;       // velocity weakening for |x| <= 10 km
X_vs = 12e3;       // velocity strengthening for |x| >= 12 km
Xmax = 60e3;
Ymax = 60e3;

// box
Point(1) = {-Xmax, -Ymax, 0, lc};
Point(2) = { Xmax, -Ymax, 0, lc};
Point(3) = { Xmax,  Ymax, 0, lc};
Point(4) = {-Xmax,  Ymax, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// fault, left to right, split where a changes
xs[] = {-Fault_half_length, -X_vs, -X_vw, X_vw, X_vs, Fault_half_length};
For i In {0 : #xs[] - 1}
   Point(100 + i) = {xs[i], 0, 0, lc_fault};
EndFor
For i In {0 : #xs[] - 2}
   Line(100 + i) = {100 + i, 101 + i};
EndFor
fault[] = {100 : 100 + #xs[] - 2};
Curve{fault[]} In Surface{1};

Field[1] = Distance;
Field[1].CurvesList = {fault[]};
Field[1].Sampling = 2000;
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);
Field[3] = MathEval;
Field[3].F = Sprintf("%g", lc);
Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Physical Curve(3) = {fault[]};
Physical Curve(5) = {1, 2, 3, 4};
Physical Surface(1) = {1};
Mesh.MshFileVersion = 2.2;
