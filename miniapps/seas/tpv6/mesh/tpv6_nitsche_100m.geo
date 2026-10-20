// SCEC TPV6 / TPV7 mesh for the Nitsche dynamic-rupture driver (seas_nitsche_dyn_driver), 100 m on
// the fault (the node spacing the description asks for).  Tetrahedra only.  TPV6 and TPV7 share this
// mesh: their bimaterial split is set by the deck's [material.halfspace_across_fault] per element
// centroid, so the mesh carries no material tags.  Same geometry as the smoke mesh
// tpv6_nitsche_1000m.geo.
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
// Size field: the TPV205 recipe (tpv205/mesh/tpv2053d_100m.geo) with the fault distance sampled every
// 50 m along strike, so the fault itself is meshed at lc_fault.
//
// Physical groups (deck [boundary]): free surface 101, fault 103, absorbing 105 (bottom and the four
// sides); volume 1.
//
//   gmsh -3 tpv6_nitsche_100m.geo -format msh22 -o tpv6_nitsche_100m.msh
//
// Built with gmsh's built-in geometry kernel only (no OpenCASCADE: the Vista gmsh module has none),
// as tpv6_nitsche_1000m.geo.
//
// The .msh is not committed (gmsh output depends on the platform; the Vista build job writes it and
// checks its sha256).  Chunhui Zhao, TPV6/TPV7/TPV31 plan, 2026-10-19.

lc_fault = 100;
lc = 5e3;

Fault_length = 30e3;
Fault_width = 15e3;
Xmax = 60e3;
Zmin = -60e3;
xf = 0.5 * Fault_length;

// The Y = 0 plane: its corners, the fault corners and the nucleation patch corners.
Point(1) = {-Xmax, 0, 0, lc};
Point(2) = {-xf, 0, 0, lc};
Point(3) = {xf, 0, 0, lc};
Point(4) = {Xmax, 0, 0, lc};
Point(5) = {Xmax, 0, Zmin, lc};
Point(6) = {-Xmax, 0, Zmin, lc};
Point(7) = {-xf, 0, -Fault_width, lc};
Point(8) = {xf, 0, -Fault_width, lc};
Point(11) = {-1500, 0, -9000, lc};
Point(12) = {1500, 0, -9000, lc};
Point(13) = {1500, 0, -6000, lc};
Point(14) = {-1500, 0, -6000, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 7};
Line(8) = {7, 8};
Line(9) = {8, 3};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};
Curve Loop(1) = {1, 7, 8, 9, 3, 4, 5, 6};
Curve Loop(2) = {2, -9, -8, -7};
Curve Loop(3) = {11, 12, 13, 14};
Plane Surface(1) = {1};       // the Y = 0 plane outside the fault
Plane Surface(2) = {2, 3};    // the fault outside the patch
Plane Surface(3) = {3};       // the nucleation patch

// The two half-boxes, Y < 0 (i = 0, tags 1xx) and Y > 0 (i = 1, tags 2xx).
For i In {0:1}
  s = 2 * i - 1;
  o = 100 * (i + 1);
  Point(o + 1) = {-Xmax, s * Xmax, 0, lc};
  Point(o + 4) = {Xmax, s * Xmax, 0, lc};
  Point(o + 5) = {Xmax, s * Xmax, Zmin, lc};
  Point(o + 6) = {-Xmax, s * Xmax, Zmin, lc};
  Line(o + 1) = {1, o + 1};
  Line(o + 4) = {4, o + 4};
  Line(o + 5) = {5, o + 5};
  Line(o + 6) = {6, o + 6};
  Line(o + 11) = {o + 1, o + 4};
  Line(o + 12) = {o + 4, o + 5};
  Line(o + 13) = {o + 5, o + 6};
  Line(o + 14) = {o + 6, o + 1};
  Curve Loop(o + 1) = {o + 11, o + 12, o + 13, o + 14};      // the far side Y = s Xmax
  Curve Loop(o + 2) = {1, 2, 3, o + 4, -(o + 11), -(o + 1)};  // the free surface Z = 0
  Curve Loop(o + 3) = {-5, o + 5, o + 13, -(o + 6)};          // the bottom Z = Zmin
  Curve Loop(o + 4) = {4, o + 5, -(o + 12), -(o + 4)};        // the side X = Xmax
  Curve Loop(o + 5) = {6, o + 1, -(o + 14), -(o + 6)};        // the side X = -Xmax
  For k In {1:5}
    Plane Surface(o + k) = {o + k};
  EndFor
  Surface Loop(i + 1) = {1, 2, 3, o + 1, o + 2, o + 3, o + 4, o + 5};
  Volume(i + 1) = {i + 1};
EndFor

Field[1] = Distance;
Field[1].SurfacesList = {2, 3};
Field[1].Sampling = 601;

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

Physical Surface(101) = {102, 202};
Physical Surface(103) = {2, 3};
Physical Surface(105) = {101, 103, 104, 105, 201, 203, 204, 205};
Physical Volume(1) = {1, 2};

Mesh.MshFileVersion = 2.2;
