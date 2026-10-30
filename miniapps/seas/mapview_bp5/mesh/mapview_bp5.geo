// Map-view BP5 mesh: the 2.5D note's Level 3 test (Chunhui Zhao,
// 2026-10-30; /mnt/project-files/2d-inplane/nonplanar-literature/
// inplane-2p5d-derivation-rev8.pdf, Section 7.3: "a map-view 2.5D version
// of SEAS BP5 ... with loading through the substrate at the BP5 plate
// rate", against the 3-D BP5-FD gold).
//
// SEAS BP5's fault seen from above: a straight vertical strike-slip fault
// y = 0, the rate-and-state segment |x| <= l_f / 2 = 50 km, and beyond it
// the plane slips at the plate rate (BP5's prescribed slip outside the
// rate-and-state region) out to the box sides.  The fault crosses the whole
// box, so no tip sits in the medium: with the 2.5D substrate moving
// +-(Vp t / 2) x on the two sides of y = 0 everywhere, a free tip would load
// the material ahead of it without bound.
//
// Geometry (m):
//   fault:  y = 0, x in [-Xmax, Xmax]; rate and state for |x| <= 50 km,
//           split at x = -30, -18 (the nucleation patch), +-30, +-32 (the
//           velocity-weakening core and the end of a's taper), so the decks'
//           steps fall on mesh vertices; creep (prescribed slip) for
//           50 < |x| <= Xmax;
//   box:    x in [-100, 100] km, y in [-60, 60] km, built as two halves
//           that share the fault line.
//
// Physical groups (the fd_cycle layout, cycle_config.hpp):
//   Physical Curve 3   = the rate-and-state segment
//   Physical Curve 5   = the two creep segments (interior Dirichlet faces:
//                        the prescribed jump Vp t, [cycle.loading] band) and
//                        the four sides (the plate motion +-(Vp t / 2) x in
//                        the QD regime, absorbing walls in the dynamic one)
//   Physical Surface 1 = bulk (both halves)
// No free surface, so no curve carries tag 1.
//
// Size: lc_fault on the rate-and-state segment, growing away from it as
// mapview_rs.geo (0.05 d + (d / 2.5 km)^2 + lc_fault, d the distance to
// that segment), capped at lc; the creep segments take the same field.
//
// Units: meters (the decks' mesh_scale = 1).
//
// Usage (from miniapps/seas, conda env pythonenv):
//   gmsh -2 mapview_bp5/mesh/mapview_bp5.geo -format msh22 \
//        -o mapview_bp5/mesh/mapview_bp5_500m.msh

DefineConstant[ lc_fault = {500, Min 1, Max 5000, Name "Fault resolution [m]" } ];
DefineConstant[ lc = {5e3, Min 100, Max 50e3, Name "Far-field resolution [m]" } ];

Xmax = 100e3;
Ymax = 60e3;

// the fault, left to right: creep, rate and state (split where the decks'
// fields change), creep
xs[] = {-Xmax, -50e3, -32e3, -30e3, -18e3, 30e3, 32e3, 50e3, Xmax};
For i In {0 : #xs[] - 1}
   Point(100 + i) = {xs[i], 0, 0, lc_fault};
EndFor
For i In {0 : #xs[] - 2}
   Line(100 + i) = {100 + i, 101 + i};
EndFor
nf = #xs[] - 1;                        // 8 segments
rs[] = {101 : 100 + nf - 2};           // the rate-and-state segments
creep[] = {100, 100 + nf - 1};         // the two creep segments

// the box corners and sides; points 100 and 100 + nf are on the sides
Point(1) = {-Xmax, -Ymax, 0, lc};
Point(2) = { Xmax, -Ymax, 0, lc};
Point(3) = { Xmax,  Ymax, 0, lc};
Point(4) = {-Xmax,  Ymax, 0, lc};
Line(1) = {1, 2};                // bottom
Line(2) = {2, 100 + nf};         // right, below the fault
Line(3) = {100 + nf, 3};         // right, above
Line(4) = {3, 4};                // top
Line(5) = {4, 100};              // left, above
Line(6) = {100, 1};              // left, below

// lower half: bottom, right, the fault from right to left, left
lower[] = {1, 2};
For i In {0 : nf - 1}
   lower[] += {-(100 + nf - 1 - i)};
EndFor
lower[] += {6};
Curve Loop(1) = {lower[]};
Plane Surface(1) = {1};

// upper half: the fault from left to right, right, top, left
upper[] = {100 : 100 + nf - 1};
upper[] += {3, 4, 5};
Curve Loop(2) = {upper[]};
Plane Surface(2) = {2};

Field[1] = Distance;
Field[1].CurvesList = {rs[]};
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

Physical Curve(3) = {rs[]};
Physical Curve(5) = {1, 2, 3, 4, 5, 6, creep[]};
Physical Surface(1) = {1, 2};
Mesh.MshFileVersion = 2.2;
