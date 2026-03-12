// BP5-QD 3D benchmark mesh for MFEM SEAS miniapp
// Adapted from Tandem: examples/tandem/3d/bp5.geo
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
//   5 = z = 0    (free surface, natural BC)
//   6 = z = Lz   (Dirichlet: u_y = sgn(x)*Vp*t/2)
//
// Usage:
//   gmsh -3 bp5.geo -o bp5_coarse.msh               (default: coarse)
//   gmsh -3 bp5.geo -setnumber res_f 5 -o bp5_med.msh   (medium)
//   gmsh -3 bp5.geo -setnumber res_f 1 -o bp5_fine.msh  (benchmark)

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

// Boolean fragment to embed fault and refinement zones into volumes
BooleanFragments{ Volume{1,2}; Delete; }{ Surface{fault,nuc1,nuc2,nuc3}; Delete; }

// --- Identify boundary surfaces ---
// x = -Lx (attr 1)
xm() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -eps, -Lx+eps, Ly+eps, Lz+eps};
// x = +Lx (attr 2)
xp() = Surface In BoundingBox{Lx-eps, -Ly-eps, -eps, Lx+eps, Ly+eps, Lz+eps};
// y = +Ly (attr 3)
yp() = Surface In BoundingBox{-Lx-eps, Ly-eps, -eps, Lx+eps, Ly+eps, Lz+eps};
// y = -Ly (attr 4)
ym() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -eps, Lx+eps, -Ly+eps, Lz+eps};
// z = 0 (attr 5, free surface)
ztop() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -eps, Lx+eps, Ly+eps, eps};
// z = Lz (attr 6, deep boundary)
zbot() = Surface In BoundingBox{-Lx-eps, -Ly-eps, Lz-eps, Lx+eps, Ly+eps, Lz+eps};

// Fault surface (internal, for reference — not a physical boundary)
fault_surfs() = Surface In BoundingBox{-eps, -l_f/2-eps, -eps, eps, l_f/2+eps, W_f+eps};

// --- Mesh sizing ---
MeshSize{ PointsOf{Volume{:};} } = res;
MeshSize{ PointsOf{Surface{fault_surfs()};} } = res_f;

// --- Physical groups ---
// Boundary surfaces (MFEM boundary attributes 1-6, 100)
Physical Surface("xm",   1) = {xm()};           // x = -Lx  (far-field, natural BC)
Physical Surface("xp",   2) = {xp()};           // x = +Lx  (far-field, natural BC)
Physical Surface("yp",   3) = {yp()};           // y = +Ly  (Dirichlet: +Vp/2)
Physical Surface("ym",   4) = {ym()};           // y = -Ly  (Dirichlet: -Vp/2)
Physical Surface("ztop", 5) = {ztop()};         // z = 0    (free surface)
Physical Surface("zbot", 6) = {zbot()};         // z = Lz   (deep boundary)
Physical Surface("fault", 100) = {fault_surfs()};  // fault plane (internal)
// Volume (tag 10 — distinct from boundary tags to avoid VTK ambiguity)
Physical Volume("domain", 10) = {Volume{:}};

Mesh.MshFileVersion = 2.2;
