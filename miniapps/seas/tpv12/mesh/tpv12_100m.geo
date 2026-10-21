// SCEC TPV12 mesh, 100 m fault resolution: a copy of tpv12_200m.geo with
// lc_fault = 100 (for p = 1, plan Q2).
// The rest of this header is tpv12_200m.geo's:
// (complex-geometries plan rev 2, Section 6; Chunhui Zhao 2026-10-20).  Built from
// tpv205/mesh/tpv2053d_200m.geo (box, size field, tags) with the TPV12
// fault: planar, 60 deg dip toward +Y, 30 km along strike (X in
// [-15, 15] km) by 15 km down-dip, trace on the X axis at Z = 0
// (TPV12_13_Description_v6.pdf).  Two lines are embedded in the fault:
//   - the stress jump at 13.8 km down-dip (depth 11951.15 m), so that
//     [prestress_depth_table] sees no fault face across it (option O1);
//   - the border of the 3 km x 3 km nucleation square centred 12 km
//     down-dip (mu_s = 0.54 inside), so each face is inside or outside.
// Outside the 30 x 15 km rectangle there is no fault: the medium is the
// strength barrier.  Built-in kernel only (no OpenCASCADE).  Tags: 101 free
// surface, 103 fault, 105 absorbing.
//   gmsh -3 tpv12_100m.geo -format msh22 -o tpv12_100m.msh

lc = 5e3;
lc_fault = 100;

Fault_length = 30e3;
Fault_width = 15e3;
Fault_dip = 60*Pi/180.;
W_jump = 13.8e3;          // down-dip distance of the stress jump
W_nucl = 12e3;            // down-dip centre of the nucleation square
R_nucl = 1.5e3;           // half side of the nucleation square

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax +  0.5 * Fault_width  *Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width  *Cos(Fault_dip);
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

// Create the fault: trace (101-102), jump line (104-105), bottom (100-103)
Point(100) = {-0.5*Fault_length, Fault_width*Cos(Fault_dip), -Fault_width*Sin(Fault_dip), lc};
Point(101) = {-0.5*Fault_length, 0, 0, lc};
Point(102) = { 0.5*Fault_length, 0, 0, lc};
Point(103) = { 0.5*Fault_length, Fault_width*Cos(Fault_dip), -Fault_width*Sin(Fault_dip), lc};
Point(104) = {-0.5*Fault_length, W_jump*Cos(Fault_dip), -W_jump*Sin(Fault_dip), lc};
Point(105) = { 0.5*Fault_length, W_jump*Cos(Fault_dip), -W_jump*Sin(Fault_dip), lc};
Line(101) = {101, 102};
Line{101} In Surface{1};
Line(102) = {102, 105};
Line(104) = {105, 104};
Line(105) = {104, 101};
Line(106) = {105, 103};
Line(103) = {103, 100};
Line(107) = {100, 104};

// Nucleation square: X in [-R, R], down-dip in [W_nucl - R, W_nucl + R]
Point(201) = { R_nucl, (W_nucl + R_nucl)*Cos(Fault_dip), -(W_nucl + R_nucl)*Sin(Fault_dip), lc};
Point(202) = { R_nucl, (W_nucl - R_nucl)*Cos(Fault_dip), -(W_nucl - R_nucl)*Sin(Fault_dip), lc};
Point(203) = {-R_nucl, (W_nucl - R_nucl)*Cos(Fault_dip), -(W_nucl - R_nucl)*Sin(Fault_dip), lc};
Point(204) = {-R_nucl, (W_nucl + R_nucl)*Cos(Fault_dip), -(W_nucl + R_nucl)*Sin(Fault_dip), lc};
Line(200) = {201, 202};
Line(201) = {202, 203};
Line(202) = {203, 204};
Line(203) = {204, 201};
Curve Loop(204) = {200,201,202,203};
Plane Surface(200) = {204};

// upper part (trace to the jump, with the nucleation square cut out) and
// the strip below the jump
Curve Loop(105) = {101,102,104,105};
Plane Surface(100) = {105, 204};
Curve Loop(115) = {-104,106,103,107};
Plane Surface(110) = {115};

// There is a bug in "Attractor", we need to define a Ruled surface in FaceList
Line Loop(106) = {101,102,104,105};
Ruled Surface(101) = {106};
Line Loop(116) = {-104,106,103,107};
Ruled Surface(111) = {116};

Surface{100,110,200} In Volume{1};

// Managing coarsening away from the fault (as tpv2053d_200m.geo)
Field[1] = Distance;
Field[1].FacesList = {101, 111};

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
Physical Surface(103) = {100,110,200};
// the extruded box faces (sides and bottom), numbered as in tpv2053d_200m.geo
Physical Surface(105) = {14,18,22,26,27};

Physical Volume(1) = {1};
Mesh.MshFileVersion = 2.2;
