// SCEC TPV102 benchmark mesh — structured Gmsh (built-in kernel), 1000 m fault resolution
// *** COARSE DEBUG-ONLY MESH ***  Λ_dyn ≈ 160 m, so Λ/h ≈ 0.16 here — process
// zone is NOT resolved and spontaneous rupture propagation is NOT expected.
// This mesh exists only for cheap local tests of cross-rank wave propagation
// (debug doc tpv102_debug_v1.md H1/H2) via the --debug-qnorm diagnostic.
// Never use for SCEC benchmark comparison.
// Adapted from Thomas Ulrich / SeisSol TPV example, updated to TPV102 spec.
// Process zone Lambda_dynamic ~ 160 m: h=200 m gives ~0.8 elements per Lambda.
//
// Geometry (SCEC TPV102, half-space):
//   Fault plane: Y=0 (code), vertical strike-slip
//   Free surface: Z=0 (code)  — SCEC y=0
//   Half-space: Z < 0         — SCEC y > 0 (depth)
//
// SCEC-to-code coordinate map:
//   SCEC x (strike)       -> code X
//   SCEC y (depth, >=0)   -> code -Z
//   SCEC z (fault-normal) -> code Y
//
// Fault sliding area (VW 30 km + 3 km transition each side):
//   X in [-18, 18] km along strike, Z in [-18, 0] km depth
// Hypocenter: (X, Z) = (0, -7.5 km)
//
// Physical groups (match current TPV102 driver defaults):
//   Physical Surface 1 = Z=0 free surface
//   Physical Surface 3 = fault (Y=0, over 36x18 km sliding area)
//   Physical Surface 5 = absorbing outer boundaries (4 sides + bottom)
//   Physical Volume  1 = bulk
//
// Units: meters. Use driver --mesh-scale 1.0 (NOT 1000).
//
// Usage:
//   gmsh -3 tpv102_200m.geo -o tpv102_200m.msh


lc = 10e3;
lc_fault = 1000;

// Fault sliding area (SCEC VW + transition, Eq. 4-5)
VW_half_length  = 15e3;
VW_depth        = 15e3;
Transition_width = 3e3;
Fault_half_length = VW_half_length + Transition_width;
Fault_length = 2*Fault_half_length;   // along-strike: 30 km VW + 3 km trans each side
Fault_width  = VW_depth + Transition_width;   // depth: 15 km VW + 3 km trans at bottom
Fault_dip    = 90*Pi/180.;

// Hypocenter / nucleation patch (TPV102 spec: (x0,y0) = (0, 7.5 km), R = 3 km)
X_nucl     = 0e3;
Width_nucl = 7.5e3;   // depth of hypocenter below free surface (SCEC y0)
R_nucl     = 3e3;     // nucleation radius (mesh refinement box half-width)
lc_nucl    = 1000;

// Domain (half-space), driver expects domain_half = 60 km
Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Zmin = -Xmax;

// Create the Volume (top surface at Z=0, extruded down to Z=Zmin)
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

// Create the fault (embedded rectangle at Y=0 inside the volume).
// Partition the fault surface along the SCEC VW/VS boundaries so the mesh conforms
// to the central 30x15 km VW patch and the surrounding 3 km transition layer.
Point(100) = {-Fault_half_length, Fault_width * Cos(Fault_dip), -Fault_width * Sin(Fault_dip), lc};
Point(101) = {-Fault_half_length, 0, 0, lc};
Point(102) = { Fault_half_length, 0, 0, lc};
Point(103) = { Fault_half_length, Fault_width * Cos(Fault_dip), -Fault_width * Sin(Fault_dip), lc};

Point(110) = {-VW_half_length, 0, 0, lc_fault};
Point(111) = { VW_half_length, 0, 0, lc_fault};
Point(112) = {-Fault_half_length, VW_depth * Cos(Fault_dip), -VW_depth * Sin(Fault_dip), lc_fault};
Point(113) = {-VW_half_length,  VW_depth * Cos(Fault_dip), -VW_depth * Sin(Fault_dip), lc_fault};
Point(114) = { VW_half_length,  VW_depth * Cos(Fault_dip), -VW_depth * Sin(Fault_dip), lc_fault};
Point(115) = { Fault_half_length, VW_depth * Cos(Fault_dip), -VW_depth * Sin(Fault_dip), lc_fault};
Point(116) = {-VW_half_length, Fault_width * Cos(Fault_dip), -Fault_width * Sin(Fault_dip), lc_fault};
Point(117) = { VW_half_length, Fault_width * Cos(Fault_dip), -Fault_width * Sin(Fault_dip), lc_fault};

Line(100) = {100, 112};
Line(101) = {112, 101};
Line(102) = {101, 110};
Line(103) = {110, 111};
Line(104) = {111, 102};
Line(105) = {102, 115};
Line(106) = {115, 103};
Line(107) = {103, 117};
Line(108) = {117, 116};
Line(109) = {116, 100};
Line(110) = {110, 113};
Line(111) = {111, 114};
Line(112) = {112, 113};
Line(113) = {113, 114};
Line(114) = {114, 115};
Line(115) = {113, 116};
Line(116) = {114, 117};

Line{102, 103, 104} In Surface{1};     // fault-surface intersection line, embedded in free surface

// Nucleation refinement patch on the fault: 2R x 2R box centered at (X_nucl, depth=Width_nucl)
Point(201) = {X_nucl + R_nucl, (Width_nucl + R_nucl) * Cos(Fault_dip), -(Width_nucl + R_nucl) * Sin(Fault_dip), lc_nucl};
Point(202) = {X_nucl + R_nucl, (Width_nucl - R_nucl) * Cos(Fault_dip), -(Width_nucl - R_nucl) * Sin(Fault_dip), lc_nucl};
Point(203) = {X_nucl - R_nucl, (Width_nucl - R_nucl) * Cos(Fault_dip), -(Width_nucl - R_nucl) * Sin(Fault_dip), lc_nucl};
Point(204) = {X_nucl - R_nucl, (Width_nucl + R_nucl) * Cos(Fault_dip), -(Width_nucl + R_nucl) * Sin(Fault_dip), lc_nucl};
Line(200) = {201, 202};
Line(201) = {202, 203};
Line(202) = {203, 204};
Line(203) = {204, 201};
Curve Loop(204) = {200,201,202,203};
Plane Surface(200) = {204};

Curve Loop(120) = {102, 110, -112, 101};
Plane Surface(120) = {120};
Curve Loop(121) = {103, 111, -113, -110};
Plane Surface(121) = {121, 204};
Curve Loop(122) = {104, 105, -114, -111};
Plane Surface(122) = {122};
Curve Loop(123) = {112, 115, 109, 100};
Plane Surface(123) = {123};
Curve Loop(124) = {113, 116, 108, -115};
Plane Surface(124) = {124};
Curve Loop(125) = {114, 106, 107, -116};
Plane Surface(125) = {125};

// Distance-field targets: reference the segmented plane surfaces directly.
// (Ruled Surface only supports 3- or 4-edge loops, which no longer matches
// the 10-edge outer boundary introduced by VW/VS partitioning.)
Ruled Surface(201) = {204};

Surface{120,121,122,123,124,125,200} In Volume{1};

// Mesh size control: refine toward fault, more refinement near nucleation patch
Field[1] = Distance;
Field[1].FacesList = {120, 121, 122, 123, 124, 125};

Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);

Field[3] = Distance;
Field[3].FacesList = {201};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc_nucl;
Field[4].LcMax = lc_fault;
Field[4].DistMin = R_nucl;
Field[4].DistMax = 3*R_nucl;

Field[5] = Restrict;
Field[5].IField = 4;
Field[5].FacesList = {120,121,122,123,124,125,200};

Field[6] = Threshold;
Field[6].IField = 1;
Field[6].LcMin = lc_fault;
Field[6].LcMax = lc;
Field[6].DistMin = 2*lc_fault;
Field[6].DistMax = 2*lc_fault+0.001;

Field[7] = Min;
Field[7].FieldsList = {2,5,6};

Background Field = 7;

// Physical groups — MUST match TPV102 driver defaults (bc-free=1, bc-fault=3, bc-absorb=5)
Physical Surface(1) = {1};                  // free surface (Z=0 top)
Physical Surface(3) = {120, 121, 122, 123, 124, 125, 200}; // fault (segmented VW/VS regions + nucleation patch)
// Outer side walls + bottom from Extrude (IDs read from Gmsh GUI, deterministic for this topology)
Physical Surface(5) = {14, 18, 22, 26, 27};

Physical Volume(1) = {1};
Mesh.MshFileVersion = 2.2;
