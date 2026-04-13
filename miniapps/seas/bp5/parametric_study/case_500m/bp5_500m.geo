// BP5-QD parametric study: 500m near-fault resolution
// Based on Tandem bp5.geo (coordinate convention: X=along-strike, Y=fault-normal, Z=depth)
//
// Domain: 400 x 200 x 100 km (same as baseline 1000m case)
// Fault: 100 km along-strike, 40 km depth, at Y=0
// Near-fault resolution: 500m (res_f=0.5)
//
// Boundary tags (Tandem convention):
//   Physical Surface 1 = top (Z=0) + bottom (Z=-100km) → Natural BC
//   Physical Surface 3 = fault → Rate-and-state friction
//   Physical Surface 5 = far-field sides → Dirichlet loading
//
// Usage:
//   conda activate pythonenv
//   gmsh -3 bp5_500m.geo -o bp5_500m.msh
//   # Or override resolution: gmsh -3 bp5_500m.geo -setnumber res_f 0.4 -o bp5_400m.msh

DefineConstant[ res =   {40, Min 0, Max 1000, Name "Far-field resolution (km)"} ];
DefineConstant[ res_f = {0.5, Min 0, Max 1000, Name "Fault resolution (km)"} ];

h_s = 2;
h_t = 2;
H = 12;
l = 60;
W_f = 40;
l_f = 100;
w = 12;

X0 = -200;
X1 = -X0;
Y0 = -100;
Y1 = -Y0;
Z0 = -100;

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

MeshSize{ PointsOf{Volume{:};} } = res;
MeshSize{ PointsOf{Surface{fault()};} } = res_f;

// Disable boundary extension to prevent implicit grading competing with Background Field
Mesh.MeshSizeExtendFromBoundary = 0;

// Distance field grading from fault surface for smooth 80:1 transition
Field[1] = Distance;
Field[1].SurfacesList = {fault()};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = res_f;      // At fault: 0.5 km
Field[2].SizeMax = res;         // Far from fault: 40 km
Field[2].DistMin = 0;           // Start grading at fault
Field[2].DistMax = 80;          // Reach far-field size at 80 km (~1.5x growth, matches bp5_v2.geo)

Background Field = 2;

Physical Surface(1) = {bottom(),top()};
Physical Surface(3) = {fault()};
Physical Surface(5) = {diri()};
Physical Volume(1) = {v()};

Mesh.MshFileVersion = 2.2;
