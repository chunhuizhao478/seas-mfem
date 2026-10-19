// SCEC TPV31 mesh for the Nitsche dynamic-rupture driver (seas_nitsche_dyn_driver), 200 m on the
// fault.  Tetrahedra only.  The layered material is set by the deck's [material] depth_profile_1d per
// element centroid; the mesh carries no material tags.
//
// Frame (as TPV205 and TPV104): X along strike, Y fault normal, Z up; free surface Z = 0.  (The TPV31
// description puts the fault on z = 0 with y as depth; this is the same rotation as tpv31_50m.geo.)
//   Box     X, Y in [-60, 60] km, Z in [-60, 0] km (the TPV205 box; P waves at 6500 m/s need 18.5 s to
//           reach a side and come back, against the 15 s run).
//   Fault   Y = 0, X in [-15, 15] km, Z in [-15, 0] km; outside it the Y = 0 plane is welded, so slip
//           is zero at the fault border as the description requires.
//   Layers  horizontal planes at Z = -2400, -5000 and -10000 m (the depths where the velocity
//           structure jumps) are mesh surfaces across the whole box, so every element lies inside one
//           layer and the jumps fall on element faces.  They carry no physical tag.
//
// Size field: the TPV205 recipe (tpv205/mesh/tpv2053d_200m.geo) with the fault distance sampled every
// 100 m, so the fault itself is meshed at lc_fault.
//
// Physical groups (deck [boundary]): free surface 101, fault 103, absorbing 105 (bottom and the four
// sides); volume 1.
//
//   gmsh -3 tpv31_nitsche_200m.geo -format msh22 -o tpv31_nitsche_200m.msh
//
// The .msh is not committed (gmsh output depends on the platform; the Vista build job writes it and
// checks its sha256).  Chunhui Zhao, TPV6/TPV7/TPV31 plan, 2026-10-19.

SetFactory("OpenCASCADE");

lc_fault = 200;
lc = 5e3;

Fault_length = 30e3;
Fault_width = 15e3;
Xmax = 60e3;
Zmin = -60e3;

Box(1) = {-Xmax, -Xmax, Zmin, 2 * Xmax, 2 * Xmax, -Zmin};

// The fault, drawn in the Z = 0 plane and turned a quarter turn about the X axis, which takes
// (x, y, 0) to (x, 0, y).
Rectangle(11) = {-0.5 * Fault_length, -Fault_width, 0, Fault_length, Fault_width};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi / 2} { Surface{11}; }

// Layer interfaces.
Rectangle(20) = {-Xmax, -Xmax, -2400, 2 * Xmax, 2 * Xmax};
Rectangle(21) = {-Xmax, -Xmax, -5000, 2 * Xmax, 2 * Xmax};
Rectangle(22) = {-Xmax, -Xmax, -10000, 2 * Xmax, 2 * Xmax};

BooleanFragments{ Volume{1}; Delete; }{ Surface{11, 20, 21, 22}; Delete; }

eps = 1.0;
fault() = Surface In BoundingBox{-0.5 * Fault_length - eps, -eps, -Fault_width - eps,
                                  0.5 * Fault_length + eps, eps, eps};
top() = Surface In BoundingBox{-Xmax - eps, -Xmax - eps, -eps, Xmax + eps, Xmax + eps, eps};
layers() = Surface In BoundingBox{-Xmax - eps, -Xmax - eps, -2400 - eps, Xmax + eps, Xmax + eps, -2400 + eps};
layers() += Surface In BoundingBox{-Xmax - eps, -Xmax - eps, -5000 - eps, Xmax + eps, Xmax + eps, -5000 + eps};
layers() += Surface In BoundingBox{-Xmax - eps, -Xmax - eps, -10000 - eps, Xmax + eps, Xmax + eps, -10000 + eps};
absorbing() = Surface{:};
absorbing() -= fault();
absorbing() -= top();
absorbing() -= layers();

Field[1] = Distance;
Field[1].SurfacesList = {fault()};
Field[1].Sampling = 301;

Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 +(F1/2.5e3)^2 + %g", lc_fault);

Field[6] = Threshold;
Field[6].InField = 1;
Field[6].SizeMin = lc_fault;
Field[6].SizeMax = lc;
Field[6].DistMin = 2 * lc_fault;
Field[6].DistMax = 2 * lc_fault + 0.001;

Field[7] = Min;
Field[7].FieldsList = {2, 6};
Background Field = 7;

Mesh.MeshSizeMax = lc;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

Physical Surface(101) = {top()};
Physical Surface(103) = {fault()};
Physical Surface(105) = {absorbing()};
Physical Volume(1) = {Volume{:}};

Mesh.MshFileVersion = 2.2;
