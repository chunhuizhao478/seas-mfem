// SCEC TPV104 mesh, y >= 0 half, for a mesh that is mirror-symmetric about the fault plane
// (Chunhui Zhao, 2026-10-17).  scripts/mirror_mesh.py reflects this half across Y = 0:
//
//   gmsh -3 tpv104/mesh/tpv104_half.geo -setnumber lc_fault 200 -format msh22 -o half.msh
//   python3 scripts/mirror_mesh.py half.msh tpv104/mesh/tpv104_mirror_ext22_200m.msh
//
// Why: on a gmsh mesh the tets on the two sides of a fault face differ, and the discrete
// normal stress on this vertical fault in a uniform half-space then changes with slip
// (1000 m test: 72-161 MPa from 120 MPa; the 100 m run 1076418 went tensile at the fault
// end).  The exact change is zero by mirror symmetry; on a mirror-symmetric mesh the
// discrete one is zero too (1000 m test: 120 MPa to every printed digit, 12 s).
//
// Fault: X in [-Fault_half_length, Fault_half_length], Z in [-Fault_width, 0].  The SCEC
// description makes the fault velocity-strengthening outside the 3 km transition; the
// default 22 km half-length and width leave a 4 km border beyond it (the deck's
// boxcar_taper gives a = 0.02, V_w = 1 m/s there), as DRDG3D's 44 km x 22 km TPV105-3D
// fault.  Fault_half_length = Fault_width = 18e3 gives the old fault extent.  The fault is
// split along the velocity-weakening edge (15 km), the transition edge (18 km) and the
// 6 km x 6 km nucleation square, as tpv104_200m.geo.
//
// Size field as tpv104_200m.geo: min(lc, 0.05 d + (d / 2.5 km)^2 + lc_fault), d = distance to
// the fault.  (tpv104_200m.geo's nucleation threshold is a no-op there, lc_nucl = lc_fault.)
//
// Coordinates (SCEC -> code): x -> X (strike), y (depth) -> -Z, z (fault normal) -> Y.
// Physical groups (TPV104 deck): surface 1 free surface (Z = 0), 3 fault, 5 absorbing,
// 99 the rest of the Y = 0 plane (interior after mirroring; mirror_mesh.py drops it);
// volume 1.

SetFactory("OpenCASCADE");

DefineConstant[ lc_fault = 200 ];
DefineConstant[ lc = 5e3 ];
DefineConstant[ Fault_half_length = 22e3 ];
DefineConstant[ Fault_width = 22e3 ];
VW_half_length   = 15e3;
VW_depth         = 15e3;
Transition_width = 3e3;
X_nucl = 0;  Z_nucl = -7.5e3;  R_nucl = 3e3;
Xmax = 60e3;  Ymax = 60e3;  Zmin = -60e3;

Box(1) = {-Xmax, 0, Zmin, 2 * Xmax, Ymax, -Zmin};

// fault pieces, made in the XY plane and turned into the Y = 0 plane: (x, y, z) -> (x, -z, y)
Rectangle(101) = {-Fault_half_length, -Fault_width, 0, 2 * Fault_half_length, Fault_width};
Rectangle(102) = {-(VW_half_length + Transition_width), -(VW_depth + Transition_width), 0,
                  2 * (VW_half_length + Transition_width), VW_depth + Transition_width};
Rectangle(103) = {-VW_half_length, -VW_depth, 0, 2 * VW_half_length, VW_depth};
Rectangle(104) = {X_nucl - R_nucl, Z_nucl - R_nucl, 0, 2 * R_nucl, 2 * R_nucl};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi / 2} { Surface{101:104}; }
BooleanFragments{ Volume{1}; Delete; }{ Surface{101:104}; Delete; }

e = 1;
fault() = Surface In BoundingBox{-Fault_half_length - e, -e, -Fault_width - e,
                                  Fault_half_length + e, e, e};
plane() = Surface In BoundingBox{-Xmax - e, -e, Zmin - e, Xmax + e, e, e};
plane() -= fault();
top()   = Surface In BoundingBox{-Xmax - e, -e, -e, Xmax + e, Ymax + e, e};
all()   = Abs(Boundary{ Volume{:}; });
absorb() = all();
absorb() -= fault();
absorb() -= plane();
absorb() -= top();

Field[1] = Distance;
Field[1].SurfacesList = {fault()};
Field[1].Sampling = 200;
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);
Field[3] = Threshold;
Field[3].InField = 1;
Field[3].SizeMin = lc_fault;
Field[3].SizeMax = lc;
Field[3].DistMin = 2 * lc_fault;
Field[3].DistMax = 2 * lc_fault + 0.001;
Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromCurvature = 0;

Physical Surface(1) = {top()};
Physical Surface(3) = {fault()};
Physical Surface(5) = {absorb()};
Physical Surface(99) = {plane()};
Physical Volume(1) = {Volume{:}};
Mesh.MshFileVersion = 2.2;
