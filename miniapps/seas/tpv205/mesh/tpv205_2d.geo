// SCEC TPV205-2D mesh (2D in-plane design rev 3, phase B, gate G1).
//
// TPV205-2D is the horizontal slice of TPV205 (= TPV5) at 7.5 km depth
// (TPV205/210 multi-resolution description): a straight fault in a 2-D
// whole space, plane strain, mode II.  The 3-D meshes of this folder
// (tpv2053d_*.geo, SeisSol's TPV5 geometry) set the fault length, the box
// and the size field; this file is their 2-D slice.
//
// Geometry (m):
//   fault: y = 0, x in [-15, 15] km, embedded in the box (outside it the
//          material is continuous: the TPV5 strength barrier is the fault's
//          end); split at x = +-1.5 km (nucleation patch) and x = +-6, +-9 km
//          (the two TPV5 stress patches, 7.5 +- 1.5 km) so the initial stress
//          steps fall on mesh vertices;
//   box:   x, y in [-60, 60] km, as the 3-D meshes; the fault ends are 45 km
//          and the fault 60 km from the box, more than the 36 km a P wave
//          (6 km/s) needs to come back to the fault within the 12 s run.
//
// Physical groups (the SeisSol tags of the 3-D TPV205 meshes):
//   Physical Curve 103   = fault
//   Physical Curve 105   = absorbing, all four sides (no free surface in 2-D)
//   Physical Surface 1   = bulk
//
// Size: lc_fault on the fault, growing away from it as the 3-D meshes do
// (0.05 d + (d / 2.5 km)^2 + lc_fault, d the distance to the fault), capped
// at lc.
//
// Units: meters (Nitsche driver --mesh-scale 1.0).
//
// Usage:
//   gmsh -2 tpv205_2d.geo -format msh22 -o tpv205_2d_100m.msh
//   gmsh -2 tpv205_2d.geo -format msh22 -setnumber lc_fault 200 -o tpv205_2d_200m.msh

DefineConstant[ lc_fault = {100, Min 1, Max 5000, Name "Fault resolution [m]" } ];
DefineConstant[ lc = {5e3, Min 100, Max 50e3, Name "Far-field resolution [m]" } ];

Fault_half_length = 15e3;
X_nucl = 0;
R_nucl = 1.5e3;
X_patch = 7.5e3;   // the left and right patches at x = -+7.5 km
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

// fault, left to right, split at the patch edges
xs[] = {-Fault_half_length, -X_patch - R_nucl, -X_patch + R_nucl, X_nucl - R_nucl,
        X_nucl + R_nucl, X_patch - R_nucl, X_patch + R_nucl, Fault_half_length};
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

Physical Curve(103) = {fault[]};
Physical Curve(105) = {1, 2, 3, 4};
Physical Surface(1) = {1};
Mesh.MshFileVersion = 2.2;
