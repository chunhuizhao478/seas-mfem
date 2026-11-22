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
//
// Built with gmsh's built-in geometry kernel only, as tpv26/mesh/tpv2627_half.geo (2f4ede17):
// the Vista gmsh/4.15.0 module has no OpenCASCADE, and job 1159764 got an empty mesh from the
// first, OpenCASCADE version of this file.  The Y = 0 face of the half-box is cut into five plane
// surfaces: the nucleation square, the velocity-weakening patch around it, the 3 km transition
// ring, the 4 km velocity-strengthening border (these four are the fault) and the rest of the
// plane.  The pieces share their edges, as the BooleanFragments of the old version did.


DefineConstant[ lc_fault = 200 ];
DefineConstant[ lc = 5e3 ];
DefineConstant[ Fault_half_length = 22e3 ];
DefineConstant[ Fault_width = 22e3 ];
VW_half_length   = 15e3;
VW_depth         = 15e3;
Transition_width = 3e3;
X_nucl = 0;  Z_nucl = -7.5e3;  R_nucl = 3e3;
Xmax = 60e3;  Ymax = 60e3;  Zmin = -60e3;

Xf = Fault_half_length;  Zf = -Fault_width;
Xt = VW_half_length + Transition_width;  Zt = -(VW_depth + Transition_width);
Xv = VW_half_length;  Zv = -VW_depth;

// points in the Y = 0 plane (x, 0, z); the mesh size comes from the fields below
Point(1)  = {-Xmax, 0, 0};    Point(2)  = {-Xf, 0, 0};    Point(3)  = {-Xt, 0, 0};
Point(4)  = {-Xv, 0, 0};      Point(5)  = {Xv, 0, 0};     Point(6)  = {Xt, 0, 0};
Point(7)  = {Xf, 0, 0};       Point(8)  = {Xmax, 0, 0};
Point(9)  = {Xmax, 0, Zmin};  Point(10) = {-Xmax, 0, Zmin};
Point(11) = {-Xf, 0, Zf};     Point(12) = {Xf, 0, Zf};
Point(13) = {-Xt, 0, Zt};     Point(14) = {Xt, 0, Zt};
Point(15) = {-Xv, 0, Zv};     Point(16) = {Xv, 0, Zv};
Point(17) = {X_nucl - R_nucl, 0, Z_nucl + R_nucl};  Point(18) = {X_nucl + R_nucl, 0, Z_nucl + R_nucl};
Point(19) = {X_nucl + R_nucl, 0, Z_nucl - R_nucl};  Point(20) = {X_nucl - R_nucl, 0, Z_nucl - R_nucl};

// the top edge Z = 0, Y = 0, in seven pieces
Line(1) = {1, 2};  Line(2) = {2, 3};  Line(3) = {3, 4};  Line(4) = {4, 5};
Line(5) = {5, 6};  Line(6) = {6, 7};  Line(7) = {7, 8};
Line(8) = {8, 9};  Line(9) = {9, 10};  Line(10) = {10, 1};       // rest of the Y = 0 boundary
Line(11) = {2, 11};  Line(12) = {11, 12};  Line(13) = {12, 7};   // fault outline (22 km)
Line(14) = {3, 13};  Line(15) = {13, 14};  Line(16) = {14, 6};   // transition outer edge (18 km)
Line(17) = {4, 15};  Line(18) = {15, 16};  Line(19) = {16, 5};   // velocity-weakening edge (15 km)
Line(20) = {17, 18};  Line(21) = {18, 19};  Line(22) = {19, 20};  Line(23) = {20, 17};  // nucleation

Curve Loop(1) = {1, 11, 12, 13, 7, 8, 9, 10};                  // rest of the plane (99)
Curve Loop(2) = {2, 14, 15, 16, 6, -13, -12, -11};             // 18-22 km border
Curve Loop(3) = {3, 17, 18, 19, 5, -16, -15, -14};             // transition ring
Curve Loop(4) = {4, -19, -18, -17};                            // velocity-weakening patch
Curve Loop(5) = {20, 21, 22, 23};                              // nucleation square
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4, 5};
Plane Surface(5) = {5};

// the rest of the half-box, Y in [0, Ymax]
Point(21) = {-Xmax, Ymax, 0};  Point(22) = {Xmax, Ymax, 0};
Point(23) = {Xmax, Ymax, Zmin};  Point(24) = {-Xmax, Ymax, Zmin};
Line(31) = {1, 21};  Line(32) = {8, 22};  Line(33) = {9, 23};  Line(34) = {10, 24};
Line(35) = {21, 22};  Line(36) = {22, 23};  Line(37) = {23, 24};  Line(38) = {24, 21};
Curve Loop(11) = {35, 36, 37, 38};                             // far side Y = Ymax
Curve Loop(12) = {1, 2, 3, 4, 5, 6, 7, 32, -35, -31};          // free surface Z = 0
Curve Loop(13) = {-9, 33, 37, -34};                            // bottom Z = Zmin
Curve Loop(14) = {8, 33, -36, -32};                            // side X = Xmax
Curve Loop(15) = {10, 31, -38, -34};                           // side X = -Xmax
For k In {11:15}
   Plane Surface(k) = {k};
EndFor
Surface Loop(1) = {1, 2, 3, 4, 5, 11, 12, 13, 14, 15};
Volume(1) = {1};

fault() = {2, 3, 4, 5};
plane() = {1};
top()   = {12};
absorb() = {11, 13, 14, 15};

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
