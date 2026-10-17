// SCEC TPV205 mesh, y >= 0 half, for a mesh that is mirror-symmetric about the fault plane
// (Chunhui Zhao, 2026-10-17).  scripts/mirror_mesh.py reflects this half across Y = 0:
//
//   gmsh -3 tpv205/mesh/tpv2053d_half.geo -setnumber lc_fault 200 -format msh22 -o half.msh
//   python3 scripts/mirror_mesh.py half.msh tpv205/mesh/tpv2053d_mirror_200m.msh
//
// Why: see tpv104/mesh/tpv104_half.geo and scripts/mirror_mesh.py (the fault normal-stress
// change on this vertical fault in a uniform medium is zero; a gmsh mesh makes it nonzero,
// 38-231 MPa in the 100 m run 1076419; a mirror-symmetric mesh keeps it zero).
//
// Geometry as tpv2053d_200m.geo (SeisSol TPV5 setup, Thomas Ulrich): 30 km x 15 km fault on
// Y = 0, split along the three 3 km x 3 km prestress squares at 7.5 km depth (x = -7.5, 0,
// +7.5 km), box X, Y in [-60, 60] km (here Y in [0, 60] km), Z in [-60, 0] km.  The rupture
// stops at the fault edges as the benchmark prescribes.  Size field as tpv2053d_200m.geo:
// min(lc, 0.05 d + (d / 2.5 km)^2 + lc_fault), d = distance to the fault (its patch
// thresholds are no-ops there, lc_nucl = lc_fault).
//
// Physical groups (TPV205 deck): surface 101 free surface (Z = 0), 103 fault, 105 absorbing,
// 99 the rest of the Y = 0 plane (interior after mirroring; mirror_mesh.py drops it);
// volume 1.

SetFactory("OpenCASCADE");

DefineConstant[ lc_fault = 200 ];
DefineConstant[ lc = 5e3 ];
Fault_length = 30e3;
Fault_width  = 15e3;
Z_patch = -7.5e3;  R_patch = 1.5e3;
Xmax = 60e3;  Ymax = 60e3;  Zmin = -60e3;

Box(1) = {-Xmax, 0, Zmin, 2 * Xmax, Ymax, -Zmin};

// fault pieces, made in the XY plane and turned into the Y = 0 plane: (x, y, z) -> (x, -z, y)
Rectangle(101) = {-0.5 * Fault_length, -Fault_width, 0, Fault_length, Fault_width};
Rectangle(102) = {-R_patch, Z_patch - R_patch, 0, 2 * R_patch, 2 * R_patch};
Rectangle(103) = {7.5e3 - R_patch, Z_patch - R_patch, 0, 2 * R_patch, 2 * R_patch};
Rectangle(104) = {-7.5e3 - R_patch, Z_patch - R_patch, 0, 2 * R_patch, 2 * R_patch};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi / 2} { Surface{101:104}; }
BooleanFragments{ Volume{1}; Delete; }{ Surface{101:104}; Delete; }

e = 1;
fault() = Surface In BoundingBox{-0.5 * Fault_length - e, -e, -Fault_width - e,
                                  0.5 * Fault_length + e, e, e};
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

Physical Surface(101) = {top()};
Physical Surface(103) = {fault()};
Physical Surface(105) = {absorb()};
Physical Surface(99) = {plane()};
Physical Volume(1) = {Volume{:}};
Mesh.MshFileVersion = 2.2;
