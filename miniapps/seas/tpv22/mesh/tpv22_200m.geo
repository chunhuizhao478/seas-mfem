// SCEC TPV22 mesh, 200 m fault resolution (complex-geometries plan rev 2,
// Section 6; Chunhui Zhao 2026-10-20).  Built from
// tpv205/mesh/tpv2053d_200m.geo (box, size field, tags) with the TPV22
// faults (TPV22_23_Description_v08.pdf): two vertical right-lateral
// strike-slip faults, 20 km deep, both reaching the free surface,
//   fault 1: Y = 0,          X in [-25, 5] km;
//   fault 2: Y = Stepover_y, X in [-5, 25] km
// (SCEC (x, y down, z) -> (X, Y, Z) = (x, z, -y); TPV22 Stepover_y = -1600,
// a right step, extensional; TPV23 = +1000).  Both traces are embedded in
// the free surface.  Each fault surface ends at its rectangle: the border
// is locked with [fault_tip_lock] (option O5), except the trace.
// Built-in kernel only (no OpenCASCADE).  Tags: 101 free surface, 103 both
// faults, 105 absorbing.
//   gmsh -3 tpv22_200m.geo -format msh22 -o tpv22_200m.msh

lc = 5e3;
lc_fault = 200;

Fault_depth = 20e3;
F1_xmin = -25e3;
F1_xmax =   5e3;
F2_xmin =  -5e3;
F2_xmax =  25e3;
Stepover_y = -1600;

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

// Fault 1: trace (101-102), right side (102-103), bottom (103-100), left side (100-101)
Point(100) = {F1_xmin, 0, -Fault_depth, lc};
Point(101) = {F1_xmin, 0, 0, lc};
Point(102) = {F1_xmax, 0, 0, lc};
Point(103) = {F1_xmax, 0, -Fault_depth, lc};
Line(100) = {100, 101};
Line(101) = {101, 102};
Line{101} In Surface{1};
Line(102) = {102, 103};
Line(103) = {103, 100};
Curve Loop(105) = {100,101,102,103};
Plane Surface(100) = {105};

// Fault 2, the same at Y = Stepover_y
Point(200) = {F2_xmin, Stepover_y, -Fault_depth, lc};
Point(201) = {F2_xmin, Stepover_y, 0, lc};
Point(202) = {F2_xmax, Stepover_y, 0, lc};
Point(203) = {F2_xmax, Stepover_y, -Fault_depth, lc};
Line(200) = {200, 201};
Line(201) = {201, 202};
Line{201} In Surface{1};
Line(202) = {202, 203};
Line(203) = {203, 200};
Curve Loop(205) = {200,201,202,203};
Plane Surface(200) = {205};

// There is a bug in "Attractor", we need to define a Ruled surface in FaceList
Line Loop(106) = {100,101,102,103};
Ruled Surface(101) = {106};
Line Loop(206) = {200,201,202,203};
Ruled Surface(201) = {206};

Surface{100,200} In Volume{1};

// Managing coarsening away from the faults (as tpv2053d_200m.geo)
Field[1] = Distance;
Field[1].FacesList = {101, 201};

Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 +(F1/2.5e3)^2 + %g", lc_fault);

// Equivalent of propagation size on element
Field[6] = Threshold;
Field[6].IField = 1;
Field[6].LcMin = lc_fault;
Field[6].LcMax = lc;
Field[6].DistMin = 2*lc_fault;
Field[6].DistMax = 2*lc_fault+0.001;

Field[7] = Min;
Field[7].FieldsList = {2,6};

Background Field = 7;

Physical Surface(101) = {1};
Physical Surface(103) = {100,200};
// the extruded box faces (sides and bottom), numbered as in tpv2053d_200m.geo
Physical Surface(105) = {14,18,22,26,27};

Physical Volume(1) = {1};
Mesh.MshFileVersion = 2.2;
