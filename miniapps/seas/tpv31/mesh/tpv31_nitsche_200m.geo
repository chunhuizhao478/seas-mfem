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
// Built with gmsh's built-in geometry kernel only (no OpenCASCADE: the Vista gmsh module has none),
// as tpv31_nitsche_1000m.geo.
//
// The .msh is not committed (gmsh output depends on the platform; the Vista build job writes it and
// checks its sha256).  Chunhui Zhao, TPV6/TPV7/TPV31 plan, 2026-10-19.

lc_fault = 200;
lc = 5e3;

Fault_length = 30e3;
Fault_width = 15e3;
Xmax = 60e3;
Zmin = -60e3;
xf = 0.5 * Fault_length;

// Depths of the free surface, the three layer interfaces and the bottom (zi), and of the fault's
// pieces: its top, the interfaces it crosses and its bottom (zf).
zi[] = {0, -2400, -5000, -10000, Zmin};
zf[] = {0, -2400, -5000, -10000, -Fault_width};

// The horizontal surface at depth zi[k] (tag 100 k + 1): the free surface, an interface, the bottom.
For k In {0:4}
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

// Layer k (volume k + 1) between zi[k] and zi[k + 1]: four vertical edges and four sides.
For k In {0:3}
  For j In {1:4}
    Line(1000 + 10 * k + j) = {100 * k + j, 100 * (k + 1) + j};
  EndFor
  For j In {1:4}
    jn = Modulo(j, 4) + 1;
    Curve Loop(1000 + 10 * k + j) = {100 * k + j, 1000 + 10 * k + jn,
                                     -(100 * (k + 1) + j), -(1000 + 10 * k + j)};
    Plane Surface(1000 + 10 * k + j) = {1000 + 10 * k + j};
  EndFor
  Surface Loop(k + 1) = {100 * k + 1, 100 * (k + 1) + 1, 1000 + 10 * k + 1, 1000 + 10 * k + 2,
                         1000 + 10 * k + 3, 1000 + 10 * k + 4};
  Volume(k + 1) = {k + 1};
EndFor

// The fault, Y = 0, in four pieces split at the interfaces; piece k (tag 2000 + k) lies in layer k
// and is embedded in it.  Its top edge lies in the free surface or an interface and is embedded there.
For k In {0:4}
  Point(2000 + 10 * k + 1) = {-xf, 0, zf[k], lc};
  Point(2000 + 10 * k + 2) = {xf, 0, zf[k], lc};
  Line(2000 + 10 * k) = {2000 + 10 * k + 1, 2000 + 10 * k + 2};
EndFor
For k In {0:3}
  Line(2000 + 10 * k + 1) = {2000 + 10 * k + 1, 2000 + 10 * (k + 1) + 1};
  Line(2000 + 10 * k + 2) = {2000 + 10 * k + 2, 2000 + 10 * (k + 1) + 2};
  Curve Loop(2000 + k) = {2000 + 10 * k, 2000 + 10 * k + 2, -(2000 + 10 * (k + 1)), -(2000 + 10 * k + 1)};
  Plane Surface(2000 + k) = {2000 + k};
  Curve{2000 + 10 * k} In Surface{100 * k + 1};
  Surface{2000 + k} In Volume{k + 1};
EndFor

Field[1] = Distance;
Field[1].SurfacesList = {2000:2003};
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

// Away from the fault the free surface, the interfaces and the bottom are identical squares with the
// same size field, so the same 2D algorithm gives them stacked copies of one triangulation, and the
// 2400 m and 2600 m layers between two copies fill with flat slivers (one 34 m thick in a first try)
// that set the time step.  Alternate the 2D algorithm from one horizontal surface to the next:
// Delaunay (5) on the -2400 m and -10000 m interfaces, the default Frontal-Delaunay (6) on the rest.
MeshAlgorithm Surface{101, 301} = 5;

Mesh.MeshSizeMax = lc;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

// The interfaces (101, 201, 301) carry no physical tag.
Physical Surface(101) = {1};
Physical Surface(103) = {2000:2003};
Physical Surface(105) = {401, 1001:1004, 1011:1014, 1021:1024, 1031:1034};
Physical Volume(1) = {1:4};

Mesh.MshFileVersion = 2.2;
