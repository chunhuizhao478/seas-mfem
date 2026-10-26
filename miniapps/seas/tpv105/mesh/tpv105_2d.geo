// SCEC TPV105-2D mesh (2D in-plane design rev 3, phase B, gate G2).
//
// TPV105 is TPV103 (rate and state, slip law with strong rate weakening, a
// whole space) plus thermal pressurization; its 2-D version is the line
// y = 7.5 km down-dip of the 3-D fault (the seven 2-D stations at x = 0,
// +-9, +-15, +-18 km, uploadTPV105_2D.pdf), a straight fault in a 2-D whole
// space, plane strain, mode II.
//
// Geometry (m):
//   fault: y = 0 across the whole box.  The description makes the fault
//          velocity strengthening outside the transition region |x| > W + w =
//          18 km and gives it no end, so the mesh has no fault tip inside the
//          box.  Split at x = +-15 km and +-18 km (the edges of the
//          transition layer of the box function B(x; W, w)) and at +-24 km,
//          the end of the finely meshed part;
//   box:   x in [-70, 70] km, y in [-50, 50] km: a P wave (6 km/s) from the
//          outermost stations (+-18 km) needs more than the 15 s of the run
//          to reach a side and come back.
//
// Physical groups (the tags of the 3-D TPV104 meshes, tpv104/mesh):
//   Physical Curve 3   = fault
//   Physical Curve 5   = absorbing, all four sides (no free surface in 2-D)
//   Physical Surface 1 = bulk (both sides of the fault)
//
// Size: lc_fault on the fault for |x| <= 24 km (6 km past the last
// station), growing away from that part as the TPV205 meshes do
// (0.05 d + (d / 2.5 km)^2 + lc_fault), capped at lc; the fault beyond
// 24 km is meshed by the same growth.
//
// Units: meters (Nitsche driver --mesh-scale 1.0).
//
// Usage:
//   gmsh -2 tpv105_2d.geo -format msh22 -o tpv105_2d_100m.msh
//   gmsh -2 tpv105_2d.geo -format msh22 -setnumber lc_fault 50 -o tpv105_2d_50m.msh

DefineConstant[ lc_fault = {100, Min 1, Max 5000, Name "Fault resolution [m]" } ];
DefineConstant[ lc = {5e3, Min 100, Max 50e3, Name "Far-field resolution [m]" } ];

W = 15e3;               // velocity-weakening half length
w = 3e3;                // transition width
Fine_half_length = 24e3;
Xmax = 70e3;
Ymax = 50e3;

// fault points, left to right: the box sides, the fine part, the transition
xs[] = {-Xmax, -Fine_half_length, -W - w, -W, W, W + w, Fine_half_length, Xmax};
For i In {0 : #xs[] - 1}
   Point(100 + i) = {xs[i], 0, 0, lc};
EndFor
Point(1) = {-Xmax, -Ymax, 0, lc};
Point(2) = { Xmax, -Ymax, 0, lc};
Point(3) = { Xmax,  Ymax, 0, lc};
Point(4) = {-Xmax,  Ymax, 0, lc};
nf = #xs[] - 1;          // fault pieces
For i In {0 : nf - 1}
   Line(100 + i) = {100 + i, 101 + i};
EndFor
fault[] = {100 : 100 + nf - 1};
fine[] = {101 : 100 + nf - 2};
// sides: bottom, right (lower, upper), top, left (upper, lower)
Line(1) = {1, 2};
Line(2) = {2, 100 + nf};
Line(3) = {100 + nf, 3};
Line(4) = {3, 4};
Line(5) = {4, 100};
Line(6) = {100, 1};
Curve Loop(1) = {1, 2, -fault[], 6};
Curve Loop(2) = {fault[], 3, 4, 5};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Field[1] = Distance;
Field[1].CurvesList = {fine[]};
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
Physical Curve(5) = {1, 2, 3, 4, 5, 6};
Physical Surface(1) = {1, 2};
Mesh.MshFileVersion = 2.2;
