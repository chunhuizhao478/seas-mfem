// The 2.5D note's Level 2 test, the 3-D reference (Chunhui Zhao,
// 2026-10-30; /mnt/project-files/2d-inplane/nonplanar-literature/
// inplane-2p5d-derivation-rev8.pdf, Section 7.2): a buried vertical
// strike-slip fault, L = 120 km by W = 15 km (L / W = 8), in a full space
// (no free surface: the note's buried case, whose c values are the ones the
// 2.5D decks use).
//
// Geometry (m):
//   fault: y = 0, x in [-60, 60] km, z in [-7.5, 7.5] km, split along the
//          3 km x 3 km nucleation patch |x|, |z| <= 1.5 km (its edges are
//          mesh edges); outside it the material is continuous;
//   box:   x in [-90, 90] km, y in [-40, 40] km, z in [-40, 40] km, all six
//          sides absorbing (30 km from the fault ends, 32.5 km above and
//          below the fault).
//
// Built-in kernel only (TACC's gmsh module has no OpenCASCADE), as
// tpv205/mesh/tpv2053d_200m.geo: the box is an extruded rectangle and the
// fault's two plane surfaces are embedded in it.
//
// Physical groups (as the TPV205 meshes): 103 fault, 105 absorbing, 1 bulk.
// Size: lc_fault on the fault, 0.05 d + (d / 2.5 km)^2 + lc_fault away
// from it (d the distance to the fault), capped at lc (tpv2053d_*.geo's
// field).
//
// Usage (from miniapps/seas, conda env pythonenv):
//   gmsh -3 elongated2p5d/mesh/elongated_3d.geo -format msh22 \
//        -setnumber lc_fault 400 -o elongated2p5d/mesh/elongated_3d_400m.msh

DefineConstant[ lc_fault = {400, Min 1, Max 5000, Name "Fault resolution [m]" } ];
DefineConstant[ lc = {5e3, Min 100, Max 50e3, Name "Far-field resolution [m]" } ];

L = 120e3;  W = 15e3;  R_nucl = 1.5e3;
Xmax = 90e3;  Ymax = 40e3;  Zmax = 40e3;

// the box: the rectangle z = -Zmax extruded to z = +Zmax
Point(1) = {-Xmax, -Ymax, -Zmax, lc};
Point(2) = { Xmax, -Ymax, -Zmax, lc};
Point(3) = { Xmax,  Ymax, -Zmax, lc};
Point(4) = {-Xmax,  Ymax, -Zmax, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
box[] = Extrude {0, 0, 2 * Zmax} { Surface{1}; };
// box[0] the top, box[1] the volume, box[2 .. 5] the four sides

// the fault, y = 0, with the nucleation patch cut out of it
Point(101) = {-L / 2, 0, -W / 2, lc_fault};
Point(102) = { L / 2, 0, -W / 2, lc_fault};
Point(103) = { L / 2, 0,  W / 2, lc_fault};
Point(104) = {-L / 2, 0,  W / 2, lc_fault};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 101};
Curve Loop(101) = {101, 102, 103, 104};

Point(201) = {-R_nucl, 0, -R_nucl, lc_fault};
Point(202) = { R_nucl, 0, -R_nucl, lc_fault};
Point(203) = { R_nucl, 0,  R_nucl, lc_fault};
Point(204) = {-R_nucl, 0,  R_nucl, lc_fault};
Line(201) = {201, 202};
Line(202) = {202, 203};
Line(203) = {203, 204};
Line(204) = {204, 201};
Curve Loop(201) = {201, 202, 203, 204};

Plane Surface(101) = {101, 201};   // the fault outside the patch
Plane Surface(201) = {201};        // the nucleation patch
Surface{101, 201} In Volume{box[1]};

Field[1] = Distance;
Field[1].SurfacesList = {101, 201};
Field[1].Sampling = 200;
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

Physical Surface(103) = {101, 201};
Physical Surface(105) = {1, box[0], box[2], box[3], box[4], box[5]};
Physical Volume(1) = {box[1]};
Mesh.MshFileVersion = 2.2;
