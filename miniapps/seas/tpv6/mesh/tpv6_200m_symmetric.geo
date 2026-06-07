// SCEC TPV6/TPV7 SYMMETRIC mesh: y-mirror prism strip near fault + unstructured
// outer.  Arm 1 (matrix + UPWIND + ADER) mesh. 200 m FIRST-PASS (100 m for production).  Adapted from tpv31/mesh/
// tpv31_50m_hybrid.geo: a structured y-mirror prism strip on each side of the
// fault so the +y and -y near-fault halves are EXACT mirror partners by
// construction — removing the near-fault upwind-dissipation asymmetry.  The ONLY
// changes vs the TPV31 hybrid are lc_fault (100 m), the Physical-Surface tags (to
// match tpv6/configs/tpv6.toml [boundary]: free=101, fault=103, absorbing=105),
// and this header.
//
// The mesh carries NO material tags: the bi-material split is assigned by the
// config's [material.halfspace_across_fault] per element centroid (B3), so the
// SAME mesh serves TPV6 and TPV7.
//
// Geometry (canonical SEAS frame: x=strike, y=fault-normal, z=depth, z=0 free):
//   Fault:  y=0 plane, x in [-15,15] km, z in [-15,0] km (SURFACE-BREAKING at z=0)
//   Domain: x,y in [-50,50] km, z in [-50,0] km
//
// Build (Frontera; .msh gitignored — do NOT commit):
//   gmsh -format msh22 -3 tpv6_200m_symmetric.geo -o tpv6_200m_symmetric.msh
//   (msh22 required — CLAUDE.md.)  100 m on the 30x15 km fault is a LARGE mesh.

SetFactory("OpenCASCADE");

// ---- Parameters --------------------------------------------------------------
lc       = 5000.0;    // far-field size at the absorbing boundary
lc_fault = 200.0;     // near-fault size (TPV6/7 spec)

Fault_half_length = 15000.0;   // x in [-15, 15] km
Fault_width       = 15000.0;   // z in [-15, 0] km (depth, surface-breaking)

// Strip thickness in y; N_strip_y layers per side at lc_fault each.
N_strip_y = 5;
h_strip   = N_strip_y * lc_fault;     // 1000 m (5 layers of 200 m)

// Outer domain (half-space)
Xmax = 50000.0;  Xmin = -Xmax;
Ymax = 50000.0;  Ymin = -Ymax;
Zmin = -50000.0;

// ---- Outer box (full domain) -------------------------------------------------
Box(1) = {Xmin, Ymin, Zmin, Xmax - Xmin, Ymax - Ymin, -Zmin};

// ---- Fault rectangle at y=0 --------------------------------------------------
Rectangle(101) = {-Fault_half_length, -Fault_width, 0,
                  2*Fault_half_length, Fault_width};
Rotate { {1, 0, 0}, {0, 0, 0}, Pi/2 } { Surface{101}; }

// ---- Extrude fault into y-mirror prism stacks --------------------------------
plus_out[]  = Extrude { 0,  h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };
minus_out[] = Extrude { 0, -h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };

// ---- Boolean fragment: embed strips into the outer box -----------------------
v_all[] = BooleanFragments
   { Volume{1}; Delete; }
   { Volume{plus_out[1], minus_out[1]}; Delete; };

// ---- Identify volumes and surfaces by bounding box ---------------------------
eps = 1.0;     // bounding-box slack [m]

v_strip_plus()  = Volume In BoundingBox{
    -Fault_half_length - eps,    0 - eps,           -Fault_width - eps,
     Fault_half_length + eps,    h_strip + eps,      0 + eps };

v_strip_minus() = Volume In BoundingBox{
    -Fault_half_length - eps,   -h_strip - eps,     -Fault_width - eps,
     Fault_half_length + eps,    0 + eps,            0 + eps };

fault_surfaces() = Surface In BoundingBox{
    -Fault_half_length - eps, -eps, -Fault_width - eps,
     Fault_half_length + eps,  eps,  0 + eps };

free_surfaces() = Surface In BoundingBox{
    Xmin - eps, Ymin - eps, -eps,
    Xmax + eps, Ymax + eps,  eps };

abs_xmin() = Surface In BoundingBox{ Xmin - eps, Ymin - eps, Zmin - eps,
                                     Xmin + eps, Ymax + eps, 0 + eps };
abs_xmax() = Surface In BoundingBox{ Xmax - eps, Ymin - eps, Zmin - eps,
                                     Xmax + eps, Ymax + eps, 0 + eps };
abs_ymin() = Surface In BoundingBox{ Xmin - eps, Ymin - eps, Zmin - eps,
                                     Xmax + eps, Ymin + eps, 0 + eps };
abs_ymax() = Surface In BoundingBox{ Xmin - eps, Ymax - eps, Zmin - eps,
                                     Xmax + eps, Ymax + eps, 0 + eps };
abs_zmin() = Surface In BoundingBox{ Xmin - eps, Ymin - eps, Zmin - eps,
                                     Xmax + eps, Ymax + eps, Zmin + eps };

all_v() = Volume "*";
v_outer() = {};
For i In {0:#all_v()-1}
    is_strip = 0;
    For j In {0:#v_strip_plus()-1}
        If (all_v(i) == v_strip_plus(j))
            is_strip = 1;
        EndIf
    EndFor
    For j In {0:#v_strip_minus()-1}
        If (all_v(i) == v_strip_minus(j))
            is_strip = 1;
        EndIf
    EndFor
    If (is_strip == 0)
        v_outer() += {all_v(i)};
    EndIf
EndFor

// ---- Mesh size control (fast-coarsening recipe) ------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfaces()};

Field[2] = MathEval;
Field[2].F = Sprintf("0.3*F1 + (F1/2.5e3)^2 + %g", lc_fault);

Field[3] = Threshold;
Field[3].IField   = 1;
Field[3].LcMin    = lc_fault;
Field[3].LcMax    = lc;
Field[3].DistMin  = 2 * lc_fault;
Field[3].DistMax  = 2 * lc_fault + 0.001;

Field[5] = Min;
Field[5].FieldsList = {2, 3};
Background Field = 5;

// ---- Physical groups — MUST match tpv6/configs/tpv6.toml [boundary]:
//        free=101, fault=103, absorbing=105.
Physical Surface(101) = {free_surfaces()};
Physical Surface(103) = {fault_surfaces()};
Physical Surface(105) = {abs_xmin(), abs_xmax(), abs_ymin(), abs_ymax(), abs_zmin()};
Physical Volume (1)   = {v_outer(), v_strip_plus(), v_strip_minus()};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm   = 6;       // Frontal-Delaunay 2D
Mesh.Algorithm3D = 1;       // Delaunay 3D
// Do NOT enable Mesh.OptimizeNetgen (SIGABRT on surface-reaching faults).
