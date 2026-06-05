// SCEC TPV31 hybrid mesh: y-mirror prism strip near fault + unstructured outer
// =============================================================================
// Follows tpv102_200m_hybrid.geo (Pelties et al. 2012 Fig. 4a "symmetric mesh"):
// a structured y-mirror prism strip on each side of the fault, so the +y and -y
// near-fault halves are EXACT mirror partners by construction — removing the
// near-fault mesh asymmetry that seeds the on-fault sigma_n speckle (the
// upwind-flux/Zhang-2023 channel).  Mesh-control variant of tpv31_50m.geo
// (unstructured): same geometry + same physical-group IDs, only the near-fault
// connectivity is made symmetric.
//
// Geometry (canonical SEAS frame: x=strike, y=fault-normal, z=depth, z=0 free):
//   Fault:  y=0 plane, x in [-15,15] km, z in [-15,0] km (SURFACE-BREAKING at z=0)
//   Domain: x,y in [-50,50] km, z in [-50,0] km
//   Free surface: z=0; absorbing: 4 vertical walls + bottom.
//
// Physical groups (MUST match tpv31/configs/tpv31.toml [boundary]):
//   Physical Surface 101 = fault (y=0)
//   Physical Surface 102 = free surface (z=0)                 [natural_attrs]
//   Physical Surface 103 = absorbing — four vertical walls    [absorbing_attrs]
//   Physical Surface 104 = absorbing — bottom (z=-50 km)      [absorbing_attrs]
//   Physical Volume  1   = bulk
//
// Build:  gmsh -format msh22 -3 tpv31_50m_hybrid.geo -o tpv31_50m_hybrid.msh
//   (msh22 required — CLAUDE.md "Known limitation — Gmsh .msh format".)
//   NB: 50 m on a 30x15 km fault makes this a LARGE mesh (structured strip is
//   ~10 M tets); expect a multi-hundred-MB to ~GB .msh and minutes to generate.
// =============================================================================

SetFactory("OpenCASCADE");

// ---- Parameters --------------------------------------------------------------
lc       = 5000.0;    // far-field size at the absorbing boundary (= tpv31_50m lc_far)
lc_fault = 50.0;      // near-fault size (= tpv31_50m lc_fault, spec p.11)

Fault_half_length = 15000.0;   // x in [-15, 15] km
Fault_width       = 15000.0;   // z in [-15, 0] km (depth, surface-breaking)

// Strip thickness in y; N_strip_y layers per side at lc_fault each.
N_strip_y = 5;
h_strip   = N_strip_y * lc_fault;     // 250 m (5 layers of 50 m)

// Outer domain (half-space)
Xmax = 50000.0;  Xmin = -Xmax;
Ymax = 50000.0;  Ymin = -Ymax;
Zmin = -50000.0;

// ---- Outer box (full domain) -------------------------------------------------
Box(1) = {Xmin, Ymin, Zmin,
          Xmax - Xmin, Ymax - Ymin, -Zmin};

// ---- Fault rectangle at y=0 --------------------------------------------------
// OCC Rectangle defaults to xy-plane; build it spanning y in [-FW, 0], then
// rotate +90 deg about x: R_x(+pi/2): (x, y, 0) -> (x, 0, +y).
// Result: Surface{101} at y=0, x in [-FHL, FHL], z in [-FW, 0] (depth).
Rectangle(101) = {-Fault_half_length, -Fault_width, 0,
                  2*Fault_half_length, Fault_width};
Rotate { {1, 0, 0}, {0, 0, 0}, Pi/2 } { Surface{101}; }

// ---- Extrude fault into y-mirror prism stacks --------------------------------
// Layers{N} -> geometry-level structured extrusion (N layers in y, lc_fault each).
// No Recombine -> triangles map to 3-tet-per-prism layers.  Both Extrudes use
// Surface{101} as their source -> identical 2D base mesh, translated by
// +/- k*lc_fault in y, so every node at (x_i, +k*lc, z_i) has an exact partner
// at (x_i, -k*lc, z_i).
plus_out[]  = Extrude { 0,  h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };
minus_out[] = Extrude { 0, -h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };

// ---- Boolean fragment: embed strips into the outer box -----------------------
v_all[] = BooleanFragments
   { Volume{1}; Delete; }
   { Volume{plus_out[1], minus_out[1]}; Delete; };

// ---- Identify volumes and surfaces by bounding box ---------------------------
// (OCC reassigns tags after Boolean ops; query by geometry.)
eps = 1.0;     // bounding-box slack [m]

v_strip_plus()  = Volume In BoundingBox{
    -Fault_half_length - eps,    0 - eps,           -Fault_width - eps,
     Fault_half_length + eps,    h_strip + eps,      0 + eps };

v_strip_minus() = Volume In BoundingBox{
    -Fault_half_length - eps,   -h_strip - eps,     -Fault_width - eps,
     Fault_half_length + eps,    0 + eps,            0 + eps };

// Fault surface(s) within the strip (y=0 plane over the fault footprint)
fault_surfaces() = Surface In BoundingBox{
    -Fault_half_length - eps, -eps, -Fault_width - eps,
     Fault_half_length + eps,  eps,  0 + eps };

// Free surface = top z=0 over the entire outer box
free_surfaces() = Surface In BoundingBox{
    Xmin - eps, Ymin - eps, -eps,
    Xmax + eps, Ymax + eps,  eps };

// Absorbing surfaces = 4 vertical sides + bottom
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

// Outer volume(s) = everything not in the strips.
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

// ---- Mesh size control (tpv31_50m fast-coarsening recipe) ---------------------
// Field 1: Distance from the fault surface(s).
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfaces()};

// Field 2: graded size, steep (0.3) slope as in tpv31_50m.geo (a shallow slope
// over-refines on the 50 m base):  F2(d) = 0.3*d + (d/2500)^2 + lc_fault.
Field[2] = MathEval;
Field[2].F = Sprintf("0.3*F1 + (F1/2.5e3)^2 + %g", lc_fault);

// Field 3: pin lc_fault within 2*lc_fault of the fault, then release to lc.
Field[3] = Threshold;
Field[3].IField   = 1;
Field[3].LcMin    = lc_fault;
Field[3].LcMax    = lc;
Field[3].DistMin  = 2 * lc_fault;
Field[3].DistMax  = 2 * lc_fault + 0.001;

Field[5] = Min;
Field[5].FieldsList = {2, 3};
Background Field = 5;

// ---- Physical groups (TPV31 convention — match tpv31/configs/tpv31.toml) ------
Physical Surface(101) = {fault_surfaces()};
Physical Surface(102) = {free_surfaces()};
Physical Surface(103) = {abs_xmin(), abs_xmax(), abs_ymin(), abs_ymax()};
Physical Surface(104) = {abs_zmin()};
Physical Volume (1)   = {v_outer(), v_strip_plus(), v_strip_minus()};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm  = 6;       // Frontal-Delaunay 2D
Mesh.Algorithm3D = 1;      // Delaunay 3D (most reliable for mixed structured/unstructured)
// Do NOT enable Mesh.OptimizeNetgen — it SIGABRTs on this surface-reaching-fault
// geometry (see tpv31_50m.geo note); the standard optimizer suffices.
