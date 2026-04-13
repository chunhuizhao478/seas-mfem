// SCEC TPV102 benchmark mesh — Coarse (h=2km, ~5k elements)
// Fault: 30km along-strike x 15km depth, centered at origin
// Domain: 60km x 60km x 30km (half-space)
//
// Coordinate convention (matching Tandem/SEAS-MFEM):
//   X = along-strike [-30, 30] km
//   Y = fault-normal [-30, 30] km
//   Z = depth [-30, 0] km (surface at Z=0)
//
// Boundary tags:
//   Physical Surface 1 = top (free surface, Z=0)
//   Physical Surface 3 = fault (Y=0, within fault extent)
//   Physical Surface 5 = far-field sides + bottom (absorbing)
//
// Usage:
//   gmsh -3 tpv102_coarse.geo -o tpv102_coarse.msh

DefineConstant[ res = {5, Min 0, Max 100, Name "Far-field resolution (km)"} ];
DefineConstant[ res_f = {2, Min 0, Max 100, Name "Fault resolution (km)"} ];

// Domain half-sizes (km)
Lx = 30;   // along-strike
Ly = 30;   // fault-normal
Lz = 30;   // depth

// Fault dimensions (km)
l_f = 30;  // along-strike extent
W_f = 15;  // depth extent

eps = 1e-3;

SetFactory("OpenCASCADE");

// Two half-volumes split at fault plane Y=0
Box(1) = {-Lx, -Ly, -Lz, 2*Lx, Ly, Lz};     // Y in [-Ly, 0]
Box(2) = {-Lx,   0, -Lz, 2*Lx, Ly, Lz};      // Y in [0, Ly]

// Fault surface at Y=0
fault = news;
Rectangle(fault) = {-l_f/2, 0, -W_f, l_f, W_f};
Rotate{ {1, 0, 0}, {0, 0, 0}, -Pi/2} { Surface{fault}; }

v() = BooleanFragments{ Volume{1,2}; Delete; }{ Surface{fault}; Delete; };

// Identify surfaces
fault() = Surface In BoundingBox{-l_f/2-eps, -eps, -W_f-eps, l_f/2+eps, eps, eps};
top() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -eps, Lx+eps, Ly+eps, eps};
bottom() = Surface In BoundingBox{-Lx-eps, -Ly-eps, -Lz-eps, Lx+eps, Ly+eps, -Lz+eps};
absorbing() = Surface{:};
absorbing() -= top();
absorbing() -= fault();

// Mesh sizing
MeshSize{ PointsOf{Volume{:};} } = res;
MeshSize{ PointsOf{Surface{fault()};} } = res_f;

Mesh.MeshSizeExtendFromBoundary = 0;

Field[1] = Distance;
Field[1].SurfacesList = {fault()};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = res_f;
Field[2].SizeMax = res;
Field[2].DistMin = 0;
Field[2].DistMax = 20;  // grade over 20 km
Background Field = 2;

// Physical groups (Tandem convention)
Physical Surface(1) = {top()};              // free surface
Physical Surface(3) = {fault()};            // fault
Physical Surface(5) = {absorbing()};        // absorbing (sides + bottom)
Physical Volume(1) = {v()};

Mesh.MshFileVersion = 2.2;
