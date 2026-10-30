// The 2.5D note's Level 2 test, the map-view side (Chunhui Zhao,
// 2026-10-30; /mnt/project-files/2d-inplane/nonplanar-literature/
// inplane-2p5d-derivation-rev8.pdf, Section 7.2): an elongated vertical
// strike-slip fault, L = 120 km (L / W = 8 for W = 15 km), in map view.
// tpv205/mesh/tpv205_2d.geo with the fault 4 x longer and the box grown to
// keep 45 km between the fault ends and the walls.
//
// Geometry (m):
//   fault: y = 0, x in [-60, 60] km, split at x = +-1.5 km (the nucleation
//          patch edges fall on mesh vertices); outside it the material is
//          continuous (the fault's end is the strength barrier);
//   box:   x in [-105, 105] km, y in [-60, 60] km, all four sides absorbing.
//
// Physical groups (as tpv205_2d.geo): 103 fault, 105 absorbing, 1 bulk.
// Size: lc_fault on the fault, 0.05 d + (d / 2.5 km)^2 + lc_fault away
// from it (d the distance to the fault), capped at lc.
//
// Usage (from miniapps/seas, conda env pythonenv):
//   gmsh -2 elongated2p5d/mesh/elongated_2d.geo -format msh22 \
//        -setnumber lc_fault 200 -o elongated2p5d/mesh/elongated_2d_200m.msh

DefineConstant[ lc_fault = {200, Min 1, Max 5000, Name "Fault resolution [m]" } ];
DefineConstant[ lc = {5e3, Min 100, Max 50e3, Name "Far-field resolution [m]" } ];

Fault_half_length = 60e3;
R_nucl = 1.5e3;
Xmax = 105e3;
Ymax = 60e3;

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

xs[] = {-Fault_half_length, -R_nucl, R_nucl, Fault_half_length};
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
Field[1].Sampling = 4000;
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
