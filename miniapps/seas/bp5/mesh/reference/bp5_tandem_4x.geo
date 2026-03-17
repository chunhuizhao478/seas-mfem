// BP5 mesh with 4x domain size for boundary truncation study.
// Domain: X in [-800, 800], Y in [-400, 400], Z in [-400, 0] km
// (4x the standard 200x100x100 km half-extents)
// Fault geometry unchanged: 100 km along-strike, 40 km down-dip.
//
// Usage:
//   conda activate pythonenv
//   gmsh -3 bp5_tandem_4x.geo -setnumber res_f 1 -o bp5_tandem_4x_1000m.msh

DefineConstant[ res =   {60, Min 0, Max 1000, Name "Resolution" } ];
DefineConstant[ res_f = {1, Min 0, Max 1000, Name "Fault resolution" } ];

h_s = 2;
h_t = 2;
H = 12;
l = 60;
W_f = 40;
l_f = 100;
w = 12;

// 4x domain extents (was X0=-200, Y0=-100, Z0=-100)
X0 = -800;
X1 = -X0;
Y0 = -400;
Y1 = -Y0;
Z0 = -400;

eps = 1e-3;

SetFactory("OpenCASCADE");
Box(1) = {X0, Y0, Z0, X1-X0, -Y0, -Z0};
Box(2) = {X0,  0, Z0, X1-X0, -Y0, -Z0};

fault = news;
Rectangle(fault) = {-l_f/2, 0, 0, l_f, W_f};
Rotate{ {1, 0, 0}, {0, 0, 0}, -Pi/2} { Surface{fault}; }
pf = PointsOf{Surface{fault};};

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

Physical Surface(1) = {bottom(),top()};
Physical Surface(3) = {fault()};
Physical Surface(5) = {diri()};
Physical Volume(1) = {v()};

Mesh.MshFileVersion = 2.2;
