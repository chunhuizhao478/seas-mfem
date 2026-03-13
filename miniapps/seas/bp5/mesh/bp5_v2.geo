// BP5-QD 3D benchmark mesh for MFEM SEAS miniapp (v2)
// Adapted from Tandem: examples/tandem/3d/bp5.geo
//
// Changes from bp5.geo:
//   - H10 fix: all-boundary Dirichlet (attrs 1-6), no free surface
//   - Region-conformal fault sub-surfaces with Physical Surface tags:
//       101 = shallow VS zone      (z in [0, h_s])
//       102 = shallow transition   (z in [h_s, h_s+h_t], y in [-l/2-h_t, l/2+h_t])
//       103 = VW core (non-nuc)    (z in [h_s+h_t, h_s+h_t+H], y in [-l/2, l/2])
//       104 = deep transition      (z in [h_s+h_t+H, h_s+2*h_t+H], y in [-l/2-h_t, l/2+h_t])
//       105 = deep VS zone         (z in [h_s+2*h_t+H, W_f])
//       106 = nucleation patch     (z in [h_s+h_t, h_s+h_t+H], y in [-l/2, -l/2+w])
//
// MFEM coordinate convention:
//   x = fault-normal       (Tandem Y)
//   y = along-strike       (Tandem X)
//   z = depth, positive down (Tandem -Z)
//
// Domain: x in [-Lx, Lx], y in [-Ly, Ly], z in [0, Lz]
// Fault plane at x = 0, y in [-lf/2, lf/2], z in [0, Wf]
//
// Boundary attributes (MFEM):
//   1 = x = -Lx  (Dirichlet: u_y = -Vp*t/2)
//   2 = x = +Lx  (Dirichlet: u_y = +Vp*t/2)
//   3 = y = +Ly  (Dirichlet: u_y = sgn(x)*Vp*t/2)
//   4 = y = -Ly  (Dirichlet: u_y = sgn(x)*Vp*t/2)
//   5 = z = 0    (Dirichlet: u_y = sgn(x)*Vp*t/2)  [H10 fix: was free surface]
//   6 = z = Lz   (Dirichlet: u_y = sgn(x)*Vp*t/2)
//
// Usage:
//   gmsh -3 bp5_v2.geo -o bp5_v2_coarse.msh               (default: coarse)
//   gmsh -3 bp5_v2.geo -setnumber res_f 5 -o bp5_v2_med.msh   (medium)
//   gmsh -3 bp5_v2.geo -setnumber res_f 1 -o bp5_v2_fine.msh  (benchmark)

// --- Resolution parameters ---
DefineConstant[ res   = {40, Min 0, Max 1000, Name "Far-field resolution (km)"} ];
DefineConstant[ res_f = {10, Min 0, Max 1000, Name "Fault resolution (km)"} ];

// --- Geometric parameters (km, matching SCEC BP5 / Tandem) ---
// Domain half-sizes
Lx = 100;    // fault-normal half-width
Ly = 200;    // along-strike half-width
Lz = 100;    // depth

// Fault dimensions
l_f = 100;   // fault length along-strike
W_f = 40;    // fault width (depth extent)

// Nucleation / transition zone (from Tandem bp5.lua)
h_s = 2;     // shallow transition
h_t = 2;     // transition half-width
H   = 12;    // seismogenic zone height
l   = 60;    // seismogenic zone along-strike extent
w   = 12;    // nucleation patch width

eps = 1e-3;  // tolerance for bounding box queries

// --- Geometry ---
SetFactory("OpenCASCADE");

// Two half-volumes split at fault plane x = 0
Box(1) = {-Lx, -Ly, 0, Lx, 2*Ly, Lz};   // x in [-Lx, 0]
Box(2) = {  0, -Ly, 0, Lx, 2*Ly, Lz};    // x in [0, Lx]

// Fault surface at x = 0, y in [-lf/2, lf/2], z in [0, Wf]
// Create in XY plane then rotate -90 deg around Y to get YZ plane
fault = news;
Rectangle(fault) = {0, -l_f/2, 0, W_f, l_f};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{fault}; }

// Nucleation refinement zones (control mesh transition)
// nuc1: outer transition — y in [-(l/2+h_t), l/2+h_t], z in [h_s, h_s+2*h_t+H]
nuc1 = news;
Rectangle(nuc1) = {h_s, -(l/2+h_t), 0, 2*h_t+H, l+2*h_t};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{nuc1}; }

// nuc2: intermediate — y in [-l/2, l/2], z in [h_s+h_t, h_s+h_t+H]
nuc2 = news;
Rectangle(nuc2) = {h_s+h_t, -l/2, 0, H, l};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{nuc2}; }

// nuc3: nucleation patch — y in [-l/2, -l/2+w], z in [h_s+h_t, h_s+h_t+H]
nuc3 = news;
Rectangle(nuc3) = {h_s+h_t, -l/2, 0, H, w};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{nuc3}; }

// VS sub-region surfaces (velocity-strengthening zones)
// vs_shallow: z in [0, h_s], y in [-l_f/2, l_f/2]
vs_shallow = news;
Rectangle(vs_shallow) = {0, -l_f/2, 0, h_s, l_f};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{vs_shallow}; }

// vs_deep: z in [h_s+2*h_t+H, W_f], y in [-l_f/2, l_f/2]
vs_deep = news;
Rectangle(vs_deep) = {h_s+2*h_t+H, -l_f/2, 0, W_f-(h_s+2*h_t+H), l_f};
Rotate{ {0, 1, 0}, {0, 0, 0}, -Pi/2} { Surface{vs_deep}; }

// Boolean fragment to embed fault and all refinement/region zones into volumes
BooleanFragments{ Volume{1,2}; Delete; }{ Surface{fault,nuc1,nuc2,nuc3,vs_shallow,vs_deep}; Delete; }

// --- Identify boundary surfaces ---
// x = -Lx (attr 1)
xm() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -eps, -Lx+eps, Ly+eps, Lz+eps};
// x = +Lx (attr 2)
xp() = Surface In BoundingBox{Lx-eps, -Ly-eps, -eps, Lx+eps, Ly+eps, Lz+eps};
// y = +Ly (attr 3)
yp() = Surface In BoundingBox{-Lx-eps, Ly-eps, -eps, Lx+eps, Ly+eps, Lz+eps};
// y = -Ly (attr 4)
ym() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -eps, Lx+eps, -Ly+eps, Lz+eps};
// z = 0 (attr 5)
ztop() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -eps, Lx+eps, Ly+eps, eps};
// z = Lz (attr 6)
zbot() = Surface In BoundingBox{-Lx-eps, -Ly-eps, Lz-eps, Lx+eps, Ly+eps, Lz+eps};

// --- Identify fault sub-region surfaces (all on x = 0 plane) ---
// All fault surfaces (tag 100)
fault_surfs() = Surface In BoundingBox{-eps, -l_f/2-eps, -eps, eps, l_f/2+eps, W_f+eps};

// Tag 101: Shallow VS — y in [-l_f/2, l_f/2], z in [0, h_s]
vs_shallow_surfs() = Surface In BoundingBox{-eps, -l_f/2-eps, -eps, eps, l_f/2+eps, h_s+eps};

// Tag 102: Shallow transition — y in [-(l/2+h_t), l/2+h_t], z in [h_s, h_s+h_t]
shallow_trans_surfs() = Surface In BoundingBox{-eps, -(l/2+h_t)-eps, h_s-eps, eps, (l/2+h_t)+eps, h_s+h_t+eps};

// Tag 103: VW core (non-nuc) — y in [-l/2, l/2], z in [h_s+h_t, h_s+h_t+H]
// This captures the full VW zone including the nucleation patch
vw_core_surfs() = Surface In BoundingBox{-eps, -l/2-eps, h_s+h_t-eps, eps, l/2+eps, h_s+h_t+H+eps};

// Tag 104: Deep transition — y in [-(l/2+h_t), l/2+h_t], z in [h_s+h_t+H, h_s+2*h_t+H]
deep_trans_surfs() = Surface In BoundingBox{-eps, -(l/2+h_t)-eps, h_s+h_t+H-eps, eps, (l/2+h_t)+eps, h_s+2*h_t+H+eps};

// Tag 105: Deep VS — y in [-l_f/2, l_f/2], z in [h_s+2*h_t+H, W_f]
vs_deep_surfs() = Surface In BoundingBox{-eps, -l_f/2-eps, h_s+2*h_t+H-eps, eps, l_f/2+eps, W_f+eps};

// Tag 106: Nucleation zone — y in [-l/2, -l/2+w], z in [h_s+h_t, h_s+h_t+H]
nuc_surfs() = Surface In BoundingBox{-eps, -l/2-eps, h_s+h_t-eps, eps, -l/2+w+eps, h_s+h_t+H+eps};

// --- Mesh sizing ---
MeshSize{ PointsOf{Volume{:};} } = res;

// Fault surface: benchmark resolution
MeshSize{ PointsOf{Surface{fault_surfs()};} } = res_f;

// Nucleation zones: intermediate refinement
nuc1_surfs() = Surface In BoundingBox{-eps, -(l/2+h_t)-eps, h_s-eps, eps, (l/2+h_t)+eps, h_s+2*h_t+H+eps};
MeshSize{ PointsOf{Surface{nuc1_surfs()};} } = res_f * 2;

nuc2_surfs() = Surface In BoundingBox{-eps, -l/2-eps, h_s+h_t-eps, eps, l/2+eps, h_s+h_t+H+eps};
MeshSize{ PointsOf{Surface{nuc2_surfs()};} } = res_f * 1.5;

nuc3_surfs() = Surface In BoundingBox{-eps, -l/2-eps, h_s+h_t-eps, eps, -l/2+w+eps, h_s+h_t+H+eps};
MeshSize{ PointsOf{Surface{nuc3_surfs()};} } = res_f;

// Smooth mesh grading from fault surface into the volume
// Prevents sharp element size jumps that cause DG traction noise
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfs()};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = res_f;      // At fault: res_f (1 km)
Field[2].SizeMax = res;         // Far from fault: res (40 km)
Field[2].DistMin = 0;           // Start grading at fault
Field[2].DistMax = 80;          // Reach far-field size at 80 km distance

Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;

// --- Physical groups ---
// Boundary surfaces (MFEM boundary attributes 1-6)
Physical Surface("xm",   1) = {xm()};           // x = -Lx  (Dirichlet)
Physical Surface("xp",   2) = {xp()};           // x = +Lx  (Dirichlet)
Physical Surface("yp",   3) = {yp()};           // y = +Ly  (Dirichlet)
Physical Surface("ym",   4) = {ym()};           // y = -Ly  (Dirichlet)
Physical Surface("ztop", 5) = {ztop()};         // z = 0    (Dirichlet, H10 fix)
Physical Surface("zbot", 6) = {zbot()};         // z = Lz   (Dirichlet)

// Fault surface — all sub-surfaces (tag 100)
Physical Surface("fault", 100) = {fault_surfs()};

// Fault sub-region surfaces
Physical Surface("vs_shallow",     101) = {vs_shallow_surfs()};     // shallow VS (amax)
Physical Surface("shallow_trans",  102) = {shallow_trans_surfs()};   // shallow transition
Physical Surface("vw_core",        103) = {vw_core_surfs()};        // VW core (a0)
Physical Surface("deep_trans",     104) = {deep_trans_surfs()};      // deep transition
Physical Surface("vs_deep",        105) = {vs_deep_surfs()};        // deep VS (amax)
Physical Surface("nucleation",     106) = {nuc_surfs()};            // nucleation patch

// Volume (tag 10 — distinct from boundary tags to avoid VTK ambiguity)
Physical Volume("domain", 10) = {Volume{:}};

Mesh.MshFileVersion = 2.2;
