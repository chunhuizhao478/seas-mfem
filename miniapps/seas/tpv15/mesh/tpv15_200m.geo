// SCEC TPV15 mesh, 200 m fault resolution (complex-geometries plan rev 2,
// Section 6; Chunhui Zhao 2026-10-20).  Built from
// tpv205/mesh/tpv2053d_200m.geo (box, size field, tags) with the TPV14/15
// faults (TPV14_15_Description_v08.pdf): two vertical faults, 15 km deep,
// both reaching the free surface,
//   main fault: Y = 0, X in [-16, 12] km, continuous through the junction
//               at the origin;
//   branch:     along (cos 30, -sin 30, 0) from the junction (toward the
//               SCEC near side, to the right), from Branch_start to 12 km.
// The branch does not touch the main fault (Q1 (a), the spec's geometry):
// its first faces start Branch_start = 100 m along the branch, 50 m from
// the main fault, and the rupture must jump the gap.  The junction field
// (Field[9]) grades the cells down to that 50 m gap along the branch's
// start edge.  Q1 (c), one element away, is Branch_start = lc_fault.
// The 3 km x 3 km nucleation square (centred 8 km left of the junction,
// 7.5 km deep) is embedded in the main fault, so each face is inside or
// outside the [[stress_region]] box that carries its traction.
// Built-in kernel only (no OpenCASCADE).  Tags: 101 free surface, 103 both
// faults, 105 absorbing.
//   gmsh -3 tpv15_200m.geo -format msh22 -o tpv15_200m.msh

lc = 5e3;
lc_fault = 200;

Fault_depth = 15e3;
Main_xmin = -16e3;
Main_xmax =  12e3;
Branch_angle = 30*Pi/180.;
Branch_start = 100;       // Q1 (a); Q1 (c): lc_fault
Branch_end = 12e3;
Gap = Branch_start*Sin(Branch_angle);   // perpendicular gap to the main fault
X_nucl = -8e3;
Z_nucl = -7.5e3;
R_nucl = 1.5e3;

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax;
Ymax =  Xmax;
Zmin = -Xmax;

// Create the Volume
Point(1) = {Xmin, Ymin, 0, lc};
Point(2) = {Xmin, Ymax, 0, lc};
Point(3) = {Xmax, Ymax, 0, lc};
Point(4) = {Xmax, Ymin, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};
Extrude {0,0, Zmin} { Surface{1}; }

// Main fault: left side (100-101), trace (101-102), right side (102-103), bottom (103-100)
Point(100) = {Main_xmin, 0, -Fault_depth, lc};
Point(101) = {Main_xmin, 0, 0, lc};
Point(102) = {Main_xmax, 0, 0, lc};
Point(103) = {Main_xmax, 0, -Fault_depth, lc};
Line(100) = {100, 101};
Line(101) = {101, 102};
Line{101} In Surface{1};
Line(102) = {102, 103};
Line(103) = {103, 100};

// Nucleation square on the main fault
Point(201) = {X_nucl + R_nucl, 0, Z_nucl + R_nucl, lc};
Point(202) = {X_nucl + R_nucl, 0, Z_nucl - R_nucl, lc};
Point(203) = {X_nucl - R_nucl, 0, Z_nucl - R_nucl, lc};
Point(204) = {X_nucl - R_nucl, 0, Z_nucl + R_nucl, lc};
Line(200) = {201, 202};
Line(201) = {202, 203};
Line(202) = {203, 204};
Line(203) = {204, 201};
Curve Loop(204) = {200,201,202,203};
Plane Surface(200) = {204};

Curve Loop(105) = {100,101,102,103};
Plane Surface(100) = {105, 204};

// Branch: start edge (300-301), trace (301-302), end edge (302-303), bottom (303-300)
Bx = Cos(Branch_angle);
By = -Sin(Branch_angle);
Point(300) = {Branch_start*Bx, Branch_start*By, -Fault_depth, lc};
Point(301) = {Branch_start*Bx, Branch_start*By, 0, lc};
Point(302) = {Branch_end*Bx, Branch_end*By, 0, lc};
Point(303) = {Branch_end*Bx, Branch_end*By, -Fault_depth, lc};
Line(300) = {300, 301};
Line(301) = {301, 302};
Line{301} In Surface{1};
Line(302) = {302, 303};
Line(303) = {303, 300};
Curve Loop(305) = {300,301,302,303};
Plane Surface(300) = {305};

// There is a bug in "Attractor", we need to define a Ruled surface in FaceList
Line Loop(106) = {100,101,102,103};
Ruled Surface(101) = {106};
Line Loop(306) = {300,301,302,303};
Ruled Surface(301) = {306};

Surface{100,200,300} In Volume{1};

// Managing coarsening away from the faults (as tpv2053d_200m.geo)
Field[1] = Distance;
Field[1].FacesList = {101, 301};

Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 +(F1/2.5e3)^2 + %g", lc_fault);

// Equivalent of propagation size on element
Field[6] = Threshold;
Field[6].IField = 1;
Field[6].LcMin = lc_fault;
Field[6].LcMax = lc;
Field[6].DistMin = 2*lc_fault;
Field[6].DistMax = 2*lc_fault+0.001;

// Junction: cells of the gap size along the branch's start edge, growing
// at 0.25 m per m (LcMax = lc, so the field never caps the far mesh; the
// fault fields take over within a few km)
Field[8] = Distance;
Field[8].EdgesList = {300};
Field[8].NNodesByEdge = 2000;

Field[9] = Threshold;
Field[9].IField = 8;
Field[9].LcMin = Gap;
Field[9].LcMax = lc;
Field[9].DistMin = Gap;
Field[9].DistMax = Gap + (lc - Gap)/0.25;

Field[7] = Min;
Field[7].FieldsList = {2,6,9};

Background Field = 7;

Physical Surface(101) = {1};
Physical Surface(103) = {100,200,300};
// the extruded box faces (sides and bottom), numbered as in tpv2053d_200m.geo
Physical Surface(105) = {14,18,22,26,27};

Physical Volume(1) = {1};
Mesh.MshFileVersion = 2.2;
