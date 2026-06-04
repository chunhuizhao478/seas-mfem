// SCEC TPV102 hybrid mesh: y-mirror prism strip near fault + unstructured outer
// =============================================================================
//
// Geometry-identical to tpv104_200m_hybrid.geo (TPV102 and TPV104 share the SCEC
// vertical strike-slip geometry; only the friction law differs, and that lives in
// the driver config, not the mesh).  This is the "symmetric mesh" of Pelties et
// al. 2012 Fig. 4(a): the near-fault strip is a structured y-mirror prism stack,
// so the +y and -y halves are exact mirror partners by construction — removing
// the near-fault mesh asymmetry that seeds the on-fault sigma_n speckle (plan
// §3.4).  Use it as a MESH control vs the unstructured tpv102_200m.geo: pure
// upwind, no over-integration, no resample.
//
// Construction:
//   1. Mesh the fault rectangle as a 2D triangulation in the y=0 plane.
//   2. Extrude {0, +h_strip, 0} { fault } Layers{N} -> structured prism stack +y.
//   3. Extrude {0, -h_strip, 0} { fault } Layers{N} -> structured prism stack -y.
//      Both extrusions reference the SAME source surface, so the +y and -y
//      halves are exact y-mirror partners by geometric construction.
//   4. BooleanFragments outer half-space with the two strip volumes -> outer
//      is carved with a strip-shaped cavity, strip volumes embedded conformally.
//
// The Layers{} directive is a geometry-level structured-extrusion instruction
// (independent of the 3D mesh algorithm) — gmsh meshes the strip by replicating
// the 2D mesh of Surface{101} at five evenly spaced y-positions per side, then
// connecting layers into prisms (split into 3 tets each since no Recombine).
//
// Geometry (SCEC TPV102, half-space):
//   Fault plane: y = 0, vertical strike-slip, x in [-18,18] km, z in [-18,0] km
//   Free surface: z = 0
//   Half-space:   z < 0 (depth)
//   Hypocenter:   (x, z) = (0, -7.5 km)
//
// Physical groups (must match TPV102 driver defaults):
//   Physical Surface 1 = z=0 free surface
//   Physical Surface 3 = fault (y=0 over the 36 x 18 km sliding area)
//   Physical Surface 5 = absorbing outer boundaries (4 sides + bottom)
//   Physical Volume  1 = bulk
//
// Build:
//   gmsh -3 tpv102_200m_hybrid.geo -o tpv102_200m_hybrid.msh -format msh22
//
// =============================================================================

SetFactory("OpenCASCADE");

// ---- Parameters --------------------------------------------------------------
lc       = 5e3;
lc_fault = 200;

VW_half_length    = 15e3;
VW_depth          = 15e3;
Transition_width  = 3e3;
Fault_half_length = VW_half_length + Transition_width;     // 18 km
Fault_width       = VW_depth + Transition_width;           // 18 km

// Strip thickness in y; 5 layers per side at 200 m each (= lc_fault).
h_strip   = 1e3;
N_strip_y = 5;

// Outer domain (full half-space)
Xmax = 60e3;  Xmin = -Xmax;
Ymax = 60e3;  Ymin = -Ymax;
Zmin = -60e3;

// ---- Outer box (full domain) -------------------------------------------------
Box(1) = {Xmin, Ymin, Zmin,
          Xmax - Xmin, Ymax - Ymin, -Zmin};

// ---- Fault rectangle at y=0 --------------------------------------------------
// OCC Rectangle defaults to xy-plane; build it spanning y in [-FW, 0], then
// rotate +90 deg around x-axis: R_x(+pi/2): (x, y, 0) -> (x, 0, +y).
// Result: Surface{101} at y=0, x in [-FHL, FHL], z in [-FW, 0] (depth).
Rectangle(101) = {-Fault_half_length, -Fault_width, 0,
                  2*Fault_half_length, Fault_width};
Rotate { {1, 0, 0}, {0, 0, 0}, Pi/2 } { Surface{101}; }

// ---- Extrude fault into y-mirror prism stacks --------------------------------
// Layers{N} -> geometry-level structured extrusion (5 layers in y, 200 m each).
// No Recombine -> triangles map to 3-tet-per-prism layers.
// Both Extrudes use Surface{101} as their source -> identical 2D base mesh,
// translated by +/- 200, 400, 600, 800, 1000 m in y.  Every node at
// (x_i, +k*200, z_i) has an exact partner at (x_i, -k*200, z_i).
plus_out[]  = Extrude { 0,  h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };
minus_out[] = Extrude { 0, -h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };
// plus_out[1]  = +y strip volume (structured)
// minus_out[1] = -y strip volume (structured)

// ---- Boolean fragment: embed strip into outer box ----------------------------
// BooleanFragments merges the topology so the strip volumes are conformally
// embedded in the outer box (shared faces become single shared faces, not
// disjoint copies).  In OCC, fragments that fully contain or are fully
// contained by another input remain topologically intact, so the structured
// extrusion attribute on the strip volumes is preserved.
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

// ---- Mesh size control (matches tpv102_200m.geo grading scheme) --------------
// Field 1: Distance from the fault surface(s).
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfaces()};

// Field 2: graded size (lc_fault near fault, quadratic+linear growth far away).
//   F2(d) = 0.05*d + (d/2500)^2 + lc_fault
//   d=1 km -> 250 m; d=10 km -> 716 m; d=20 km -> 1264 m; d=50 km -> 3100 m.
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);

// Field 3: sharp threshold cap at 2*lc_fault from the fault.
//   Inside 400 m of fault -> lc_fault (200 m); beyond -> lc (5000 m).
Field[3] = Threshold;
Field[3].IField   = 1;
Field[3].LcMin    = lc_fault;
Field[3].LcMax    = lc;
Field[3].DistMin  = 2 * lc_fault;
Field[3].DistMax  = 2 * lc_fault + 0.001;

Field[5] = Min;
Field[5].FieldsList = {2, 3};
Background Field = 5;

// ---- Physical groups ---------------------------------------------------------
Physical Surface(1) = {free_surfaces()};
Physical Surface(3) = {fault_surfaces()};
Physical Surface(5) = {abs_xmin(), abs_xmax(), abs_ymin(), abs_ymax(),
                       abs_zmin()};
Physical Volume (1) = {v_outer(), v_strip_plus(), v_strip_minus()};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm  = 6;       // Frontal-Delaunay 2D
Mesh.Algorithm3D = 1;      // Delaunay 3D (most reliable for mixed structured/unstructured)
