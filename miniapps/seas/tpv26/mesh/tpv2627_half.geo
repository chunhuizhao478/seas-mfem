// SCEC TPV26 / TPV27 mesh, y >= 0 half, for a mesh that is mirror-symmetric about the fault
// plane (rheology plan rev 1, phase R1; Chunhui Zhao, 2026-10-20).  scripts/mirror_mesh.py
// reflects this half across Y = 0 (tpv26/mesh/make_tpv2627_meshes.sh does both):
//
//   gmsh -3 tpv26/mesh/tpv2627_half.geo -setnumber lc_fault 200 -format msh22 -o half.msh
//   python3 scripts/mirror_mesh.py half.msh tpv26/mesh/tpv2627_mirror_200m.msh --fault-attr 103
//
// Why the mirror: on a gmsh mesh the tets on the two sides of a fault face differ and the
// discrete fault normal stress then changes with slip (TPV205 100 m, Vista 1076419: 38-231 MPa
// against 120 MPa exact); the mirror mesh keeps the change zero for an elastic bulk
// (scripts/mirror_mesh.py, commit 41879e8).  TPV26 is elastic; TPV27 runs on the same mesh.
//
// Geometry (TPV26/27 description v13, Parts 2-3), code frame X along strike, Y fault normal
// (SCEC z, the far side at Y > 0), Z up (SCEC -y): fault X in [-20, 20] km, Z in [-20, 0] km on
// Y = 0, reaching the free surface; hypocenter (-5, 0, -10) km.  Half-space box X in
// [-Xmax, Xmax], Y in [0, Ymax] (mirrored to [-Ymax, Ymax]), Z in [Zmin, 0]; the margins
// (45 km along strike, 40 km below the fault, 60 km across) follow tpv2053d_half.geo.  The
// rupture stops at the fault's buried edges as the benchmark prescribes (no fault faces past
// them).
//
// Size field: min(lc, 0.05 d + (d / 2.5 km)^2 + lc_fault) with d the distance to the fault
// (tpv2053d_half.geo), and inside the box |Y| <= box_halfwidth, X within box_halfwidth of the
// fault ends, depth up to fault depth + box_halfwidth, at most box_factor * lc_fault: the
// off-fault plastic zone of TPV27.  box_halfwidth = 0 drops the box.
//
// Physical groups (TPV26/27 decks): surface 101 free surface (Z = 0), 103 fault, 105 absorbing,
// 99 the rest of the Y = 0 plane (interior after mirroring; mirror_mesh.py drops it);
// volume 1.

SetFactory("OpenCASCADE");

DefineConstant[ lc_fault = 200 ];        // element size on the fault (m)
DefineConstant[ lc = 5e3 ];              // largest element size (m)
DefineConstant[ box_halfwidth = 2e3 ];   // plastic-zone box half-width (m); 0 = no box
DefineConstant[ box_factor = 2 ];        // element size in the box, in units of lc_fault
Fault_length = 40e3;
Fault_width  = 20e3;
Xmax = 65e3;  Ymax = 60e3;  Zmin = -60e3;

Box(1) = {-Xmax, 0, Zmin, 2 * Xmax, Ymax, -Zmin};

// the fault, made in the XY plane and turned into the Y = 0 plane: (x, y, z) -> (x, -z, y)
Rectangle(101) = {-0.5 * Fault_length, -Fault_width, 0, Fault_length, Fault_width};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi / 2} { Surface{101}; }
BooleanFragments{ Volume{1}; Delete; }{ Surface{101}; Delete; }

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
Field[1].Sampling = 300;
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);
Field[3] = Threshold;
Field[3].InField = 1;
Field[3].SizeMin = lc_fault;
Field[3].SizeMax = lc;
Field[3].DistMin = 2 * lc_fault;
Field[3].DistMax = 2 * lc_fault + 0.001;
If (box_halfwidth > 0)
   Field[5] = Box;
   Field[5].VIn = box_factor * lc_fault;
   Field[5].VOut = lc;
   Field[5].XMin = -0.5 * Fault_length - box_halfwidth;
   Field[5].XMax = 0.5 * Fault_length + box_halfwidth;
   Field[5].YMin = -box_halfwidth;
   Field[5].YMax = box_halfwidth;
   Field[5].ZMin = -Fault_width - box_halfwidth;
   Field[5].ZMax = e;
   Field[4] = Min;
   Field[4].FieldsList = {2, 3, 5};
Else
   Field[4] = Min;
   Field[4].FieldsList = {2, 3};
EndIf
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
