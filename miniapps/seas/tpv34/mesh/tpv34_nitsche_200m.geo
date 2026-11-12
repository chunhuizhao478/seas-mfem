// SCEC TPV34 mesh for the Nitsche dynamic-rupture driver (seas_nitsche_dyn_driver), 200 m on the
// fault.  Tetrahedra only.  The material comes from the deck's 3-D table ([material] kind
// sidecar_hdf5) per element centroid; it is smooth, so the mesh has no interfaces and no material tags.
//
// Frame (as TPV205, TPV104, TPV31): X along strike, Y fault normal, Z up; free surface Z = 0.  (The
// TPV34 description puts the fault on z = 0 with y as depth: x = X, y = -Z, z = Y.)
//   Box     X, Y in [-80, 80] km, Z in [-80, 0] km: the run lasts 20 s and Vp reaches 7.86 km/s below
//           25 km depth in the CVM-H table, so a P wave from the hypocenter needs about 20 s to reach a
//           side or the bottom and come back.  The material table covers the same box.
//   Fault   Y = 0, X in [-15, 15] km, Z in [-15, 0] km; outside it the Y = 0 plane is welded, so slip
//           is zero at the fault border as the description requires.
//
// Size field: the TPV205 recipe (tpv205/mesh/tpv2053d_200m.geo), as tpv31/mesh/tpv31_nitsche_200m.geo.
//
// Physical groups (deck [boundary]): free surface 101, fault 103, absorbing 105 (bottom and the four
// sides); volume 1.
//
//   gmsh -3 tpv34_nitsche_200m.geo -format msh22 -o tpv34_nitsche_200m.msh
//
// Built with gmsh's built-in geometry kernel only (no OpenCASCADE: the Vista gmsh module has none).
//
// The .msh is not committed (gmsh output depends on the platform; the build job writes it and checks
// its sha256).  Chunhui Zhao, TPV34 plan rev 2, 2026-11-12.

lc_fault = 200;
lc = 5e3;

Fault_length = 30e3;
Fault_width = 15e3;
Xmax = 80e3;
Zmin = -80e3;
xf = 0.5 * Fault_length;

// The free surface (tag 1) and the bottom (tag 101).
zi[] = {0, Zmin};
For k In {0:1}
  Point(100 * k + 1) = {-Xmax, -Xmax, zi[k], lc};
  Point(100 * k + 2) = {Xmax, -Xmax, zi[k], lc};
  Point(100 * k + 3) = {Xmax, Xmax, zi[k], lc};
  Point(100 * k + 4) = {-Xmax, Xmax, zi[k], lc};
  For j In {1:4}
    Line(100 * k + j) = {100 * k + j, 100 * k + Modulo(j, 4) + 1};
  EndFor
  Curve Loop(100 * k + 1) = {100 * k + 1, 100 * k + 2, 100 * k + 3, 100 * k + 4};
  Plane Surface(100 * k + 1) = {100 * k + 1};
EndFor

// The box (volume 1): four vertical edges and four sides (tags 1001 to 1004).
For j In {1:4}
  Line(1000 + j) = {j, 100 + j};
EndFor
For j In {1:4}
  jn = Modulo(j, 4) + 1;
  Curve Loop(1000 + j) = {j, 1000 + jn, -(100 + j), -(1000 + j)};
  Plane Surface(1000 + j) = {1000 + j};
EndFor
Surface Loop(1) = {1, 101, 1001, 1002, 1003, 1004};
Volume(1) = {1};

// The fault (surface 2000), Y = 0, embedded in the box; its top edge is embedded in the free surface.
Point(2001) = {-xf, 0, 0, lc};
Point(2002) = {xf, 0, 0, lc};
Point(2011) = {-xf, 0, -Fault_width, lc};
Point(2012) = {xf, 0, -Fault_width, lc};
Line(2000) = {2001, 2002};
Line(2010) = {2011, 2012};
Line(2001) = {2001, 2011};
Line(2002) = {2002, 2012};
Curve Loop(2000) = {2000, 2002, -2010, -2001};
Plane Surface(2000) = {2000};
Curve{2000} In Surface{1};
Surface{2000} In Volume{1};

Field[1] = Distance;
Field[1].SurfacesList = {2000};
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

Physical Surface(101) = {1};
Physical Surface(103) = {2000};
Physical Surface(105) = {101, 1001:1004};
Physical Volume(1) = {1};

Mesh.MshFileVersion = 2.2;
