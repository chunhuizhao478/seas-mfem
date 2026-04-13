// BP5-QD parametric study: 2x domain size, 1000m near-fault resolution
// Based on Tandem bp5.geo (coordinate convention: X=along-strike, Y=fault-normal, Z=depth)
//
// Domain: 800 x 400 x 200 km (2x baseline in each dimension)
//   X (along-strike):  -400 to 400 km  (was -200 to 200)
//   Y (fault-normal):  -200 to 200 km  (was -100 to 100)
//   Z (depth):         -200 to 0 km    (was -100 to 0)
//
// Fault: 100 km along-strike, 40 km depth, at Y=0 (UNCHANGED)
// Near-fault resolution: 1000m (res_f=1, same as baseline)
//
// Boundary tags (Tandem convention):
//   Physical Surface 1 = top (Z=0) + bottom (Z=-200km) → Natural BC
//   Physical Surface 3 = fault → Rate-and-state friction
//   Physical Surface 5 = far-field sides → Dirichlet loading
//
// Uses distance field grading to prevent excessive far-field elements.
//
// Usage:
//   conda activate pythonenv
//   gmsh -3 bp5_2x_1000m.geo -o bp5_2x_1000m.msh

DefineConstant[ res =   {60, Min 0, Max 2000, Name "Far-field resolution (km)"} ];
DefineConstant[ res_f = {1,  Min 0, Max 1000, Name "Fault resolution (km)"} ];

h_s = 2;
h_t = 2;
H = 12;
l = 60;
W_f = 40;
l_f = 100;
w = 12;

// 2x domain dimensions
X0 = -400;
X1 = -X0;
Y0 = -200;
Y1 = -Y0;
Z0 = -200;

eps = 1e-3;

SetFactory("OpenCASCADE");
Box(1) = {X0, Y0, Z0, X1-X0, -Y0, -Z0};
Box(2) = {X0,  0, Z0, X1-X0, -Y0, -Z0};

fault = news;
Rectangle(fault) = {-l_f/2, 0, 0, l_f, W_f};
Rotate{ {1, 0, 0}, {0, 0, 0}, -Pi/2} { Surface{fault}; }

nuc1 = news;
Rectangle(nuc1) = {-l/2-h_t, 0, -h_s, l+2*h_t, 2*h_t + H};
Rotate{ {1, 0, 0}, {0, 0, -h_s}, -Pi/2} { Surface{nuc1}; }

nuc2 = news;
Rectangle(nuc2) = {-l/2, 0, -h_s-h_t, l, H};
Rotate{ {1, 0, 0}, {0, 0, -h_s-h_t}, -Pi/2} { Surface{nuc2}; }

nuc3 = news;
Rectangle(nuc3) = {-l/2, 0, -h_s-h_t, w, H};
Rotate{ {1, 0, 0}, {0, 0, -h_s-h_t}, -Pi/2} { Surface{nuc3}; }

v() = BooleanFragments{ Volume{1,2}; Delete; }{ Surface{fault,nuc1,nuc2,nuc3}; Delete; };

fault() = Surface In BoundingBox{-l_f/2-eps, -eps, -W_f-eps, l_f/2+eps, eps, W_f+eps};
top() = Surface In BoundingBox{X0-eps, Y0-eps, -eps, X1+eps, Y1+eps, eps};
bottom() = Surface In BoundingBox{X0-eps, Y0-eps, Z0-eps, X1+eps, Y1+eps, Z0+eps};
diri() = Surface{:};
diri() -= top();
diri() -= bottom();
diri() -= fault();

// Mesh sizing: uniform on fault, distance-based grading into the volume
MeshSize{ PointsOf{Volume{:};} } = res;
MeshSize{ PointsOf{Surface{fault()};} } = res_f;

// Disable boundary extension to prevent implicit grading competing with Background Field
Mesh.MeshSizeExtendFromBoundary = 0;

// Distance field grading from fault surface
// Prevents excessive elements in the enlarged far-field region
Field[1] = Distance;
Field[1].SurfacesList = {fault()};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = res_f;      // At fault: 1 km
Field[2].SizeMax = res;         // Far from fault: 60 km
Field[2].DistMin = 0;           // Start grading at fault
Field[2].DistMax = 160;         // Reach far-field size at 160 km (2x baseline 80 km)

Background Field = 2;

Physical Surface(1) = {bottom(),top()};
Physical Surface(3) = {fault()};
Physical Surface(5) = {diri()};
Physical Volume(1) = {v()};

Mesh.MshFileVersion = 2.2;
