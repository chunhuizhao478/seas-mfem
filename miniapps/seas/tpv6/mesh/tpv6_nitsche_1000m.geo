// SCEC TPV6 / TPV7 mesh for the Nitsche dynamic-rupture driver (seas_nitsche_dyn_driver), 1000 m on
// the fault (smoke runs only: it does not resolve the cohesive zone).  Tetrahedra only.  TPV6 and TPV7
// share this mesh: their bimaterial split is set by the deck's [material.halfspace_across_fault] per
// element centroid, so the mesh carries no material tags.
//
// Frame (as TPV205 and TPV104): X along strike, Y fault normal, Z up; free surface Z = 0.
//   Box     X, Y in [-60, 60] km, Z in [-60, 0] km (the TPV205 box; P waves at 6000 m/s need 20 s to
//           reach a side and come back, against the 12 s run).
//   Fault   Y = 0, X in [-15, 15] km, Z in [-15, 0] km (the 30 km x 15 km rupture area of the spec;
//           outside it the Y = 0 plane is welded, which is the spec's unbreakable barrier).
//   Patch   the 3 km x 3 km nucleation square (X in [-1.5, 1.5] km, Z in [-9, -6] km) is a sub-surface
//           of the fault, so its edges are element edges and each fault face gets one prestress value
//           (as on the TPV205 meshes).
//   Plane   the whole Y = 0 plane is a mesh surface, so the material interface between the two
//           half-spaces lies on element faces everywhere, not only on the fault.  Outside the fault
//           its faces carry no physical tag and are ordinary interior faces.
//
// Size field: the TPV205 recipe (tpv205/mesh/tpv2053d_200m.geo) with the fault distance sampled every
// 100 m, so the fault itself is meshed at lc_fault.
//
// Physical groups (deck [boundary]): free surface 101, fault 103, absorbing 105 (bottom and the four
// sides); volume 1.
//
//   gmsh -3 tpv6_nitsche_1000m.geo -format msh22 -o tpv6_nitsche_1000m.msh
//
// The .msh is not committed (gmsh output depends on the platform; the Vista build job writes it and
// checks its sha256).  Chunhui Zhao, TPV6/TPV7/TPV31 plan, 2026-10-19.

SetFactory("OpenCASCADE");

lc_fault = 1000;
lc = 5e3;

Fault_length = 30e3;
Fault_width = 15e3;
Xmax = 60e3;
Zmin = -60e3;

Box(1) = {-Xmax, -Xmax, Zmin, 2 * Xmax, 2 * Xmax, -Zmin};

// The Y = 0 plane, the fault and the nucleation patch, drawn in the Z = 0 plane and turned a quarter
// turn about the X axis, which takes (x, y, 0) to (x, 0, y).
Rectangle(10) = {-Xmax, Zmin, 0, 2 * Xmax, -Zmin};
Rectangle(11) = {-0.5 * Fault_length, -Fault_width, 0, Fault_length, Fault_width};
Rectangle(12) = {-1500, -9000, 0, 3000, 3000};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi / 2} { Surface{10, 11, 12}; }

BooleanFragments{ Volume{1}; Delete; }{ Surface{10, 11, 12}; Delete; }

eps = 1.0;
fault() = Surface In BoundingBox{-0.5 * Fault_length - eps, -eps, -Fault_width - eps,
                                  0.5 * Fault_length + eps, eps, eps};
plane() = Surface In BoundingBox{-Xmax - eps, -eps, Zmin - eps, Xmax + eps, eps, eps};
top() = Surface In BoundingBox{-Xmax - eps, -Xmax - eps, -eps, Xmax + eps, Xmax + eps, eps};
absorbing() = Surface{:};
absorbing() -= plane();
absorbing() -= top();

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
