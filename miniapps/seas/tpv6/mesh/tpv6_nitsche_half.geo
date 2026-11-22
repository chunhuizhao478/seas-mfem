// SCEC TPV6 / TPV7 mesh, Y >= 0 half, for a mesh that is mirror-symmetric about the fault plane
// (Chunhui Zhao, 2026-11-22).  scripts/mirror_mesh.py reflects this half across Y = 0:
//
//   gmsh -3 tpv6/mesh/tpv6_nitsche_half.geo -setnumber lc_fault 100 -format msh22 -o half.msh
//   python3 scripts/mirror_mesh.py half.msh tpv6/mesh/tpv6_nitsche_mirror_100m.msh --fault-attr 103
//
// Why: TPV7 100 m p = 3 (Vista 1159765) aborted at t = 5.99 s with sigma_n_eff = -55 kPa at
// (14264, 0, -49) m, the fault end at the free surface, on tpv6_nitsche_100m.msh, whose two halves
// gmsh meshes separately (as TPV104 1076418 and 1156772, and TPV6 1085891 and 1098967).  TPV6 and
// TPV7 are bimaterial, so their sigma_n changes physically; a mirror mesh removes the part that
// comes from the two sides' tets differing (cloud TPV6 1000 m test 2026-10-25: one rock gmsh 108-128
// MPa, mirrored 120.00 flat; two rocks mirrored 88.1-140.2, DRDG3D 89.9-135.4).  The material split
// is set per element centroid by [material.halfspace_across_fault], so the mirrored mesh needs no
// material tags.
//
// Geometry and size field as tpv6_nitsche_100m.geo (the Y > 0 half-box of its loop), lc_fault set
// with -setnumber.  Physical groups (TPV6/TPV7 decks): surface 101 free surface, 103 fault, 105
// absorbing, 99 the rest of the Y = 0 plane (interior after mirroring; mirror_mesh.py drops it);
// volume 1.  Built with gmsh's built-in geometry kernel only (the Vista gmsh module has no
// OpenCASCADE).

DefineConstant[ lc_fault = 100 ];
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

// The half-box Y > 0 (tags 2xx, as the i = 1 pass of tpv6_nitsche_100m.geo).
o = 200;
Point(o + 1) = {-Xmax, Xmax, 0, lc};
Point(o + 4) = {Xmax, Xmax, 0, lc};
Point(o + 5) = {Xmax, Xmax, Zmin, lc};
Point(o + 6) = {-Xmax, Xmax, Zmin, lc};
Line(o + 1) = {1, o + 1};
Line(o + 4) = {4, o + 4};
Line(o + 5) = {5, o + 5};
Line(o + 6) = {6, o + 6};
Line(o + 11) = {o + 1, o + 4};
Line(o + 12) = {o + 4, o + 5};
Line(o + 13) = {o + 5, o + 6};
Line(o + 14) = {o + 6, o + 1};
Curve Loop(o + 1) = {o + 11, o + 12, o + 13, o + 14};      // the far side Y = Xmax
Curve Loop(o + 2) = {1, 2, 3, o + 4, -(o + 11), -(o + 1)};  // the free surface Z = 0
Curve Loop(o + 3) = {-5, o + 5, o + 13, -(o + 6)};          // the bottom Z = Zmin
Curve Loop(o + 4) = {4, o + 5, -(o + 12), -(o + 4)};        // the side X = Xmax
Curve Loop(o + 5) = {6, o + 1, -(o + 14), -(o + 6)};        // the side X = -Xmax
For k In {1:5}
  Plane Surface(o + k) = {o + k};
EndFor
Surface Loop(1) = {1, 2, 3, o + 1, o + 2, o + 3, o + 4, o + 5};
Volume(1) = {1};

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

Physical Surface(101) = {202};
Physical Surface(103) = {2, 3};
Physical Surface(105) = {201, 203, 204, 205};
Physical Surface(99) = {1};
Physical Volume(1) = {1};

Mesh.MshFileVersion = 2.2;
