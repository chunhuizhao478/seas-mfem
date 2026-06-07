// SCEC TPV102 hybrid mesh, STRIKE-SYMMETRIC variant
// =============================================================================
//
// Variant of tpv102_200m_hybrid.geo that fixes the on-fault DIP-SLIP DRIFT.
//
// Root cause (debug_document/dip_drift_local_reproduction_2026-06-06.md):
//   The dip drift is ANTI-symmetric in along-strike x; every solver operator is
//   proven exactly x->-x equivariant (compiled probes), so the ONLY strike (x)
//   asymmetry is the MESH.  tpv102_200m_hybrid.geo is symmetric across the FAULT
//   PLANE (the +/-y prism strips are mirror partners — this tamed the sigma_n
//   leak) but its fault triangulation (Frontal-Delaunay) is STRIKE-asymmetric:
//   only 1.78% of fault nodes have an exact -x mirror partner.  An x-asymmetric
//   fault mesh seeds an x-antisymmetric spurious [[v_z]] that the (clean) friction
//   rectifies into secular dip slip.
//
// Fix: make the fault-plane triangulation MIRROR-SYMMETRIC under x->-x about the
// hypocenter (x=0), while KEEPING the y-mirror prism-strip extrusion that fixes
// sigma_n.  This is the strike analogue of the Pelties 2012 y-mirror strip.
//
// Construction (vs the parent .geo, which used ONE Frontal-Delaunay fault Rectangle):
//   1. Split the fault rectangle at x=0 into a RIGHT half (x in [0,FHL]) and a
//      LEFT half (x in [-FHL,0]) sharing the x=0 edge.  Each half is an 18x18 km
//      SQUARE, so a uniform transfinite (structured) grid is natural.
//   2. Make each half a TRANSFINITE surface with MIRRORED diagonals:
//      right half "Right" diagonals, left half "Left" diagonals.  Under x->-x a
//      "Right" (/) diagonal maps to a "Left" (\) diagonal, so the two halves are
//      EXACT strike mirrors — node positions AND triangle connectivity.
//   3. Extrude BOTH halves +/- h_strip in y -> structured y-mirror prism stacks
//      (same Layers{N_strip_y} as the parent; +y/-y stay exact partners).  The
//      extrusion carries the strike symmetry into the near-fault strip.
//   4. BooleanFragments the outer half-space with the four strip volumes.
//   5. RE-ASSERT the transfinite constraint on the POST-boolean fault surfaces.
//      (BooleanFragments drops the transfinite attribute set on the pre-boolean
//      source surfaces; without this re-assertion the 3D mesher re-triangulates
//      the fault and the strike symmetry drops to ~28%.  With it: 100%.)
//
// Verified at coarse resolution (gmsh 4.15, MFEM load): manifold, MFEM-loadable,
// Physical groups {1,3,5}/Volume 1 intact, fault footprint x in [-18,18] km,
// fault triangle x-mirror = 100.00% (was 1.78%), strip y-mirror = 100.00%.
// No $Periodic section is written (the parent's reader is a brittle v2.2 parser).
// Solver code is UNCHANGED — this is a drop-in mesh control.
//
// Geometry (SCEC TPV102, half-space):
//   Fault plane: y = 0, vertical strike-slip, x in [-18,18] km, z in [-18,0] km
//   Free surface: z = 0 ;  Half-space z < 0 ;  Hypocenter (x,z) = (0, -7.5 km)
//
// Physical groups (must match TPV102 driver defaults):
//   Physical Surface 1 = z=0 free surface
//   Physical Surface 3 = fault (y=0 over the 36 x 18 km sliding area; two halves)
//   Physical Surface 5 = absorbing outer boundaries (4 sides + bottom)
//   Physical Volume  1 = bulk
//
// Build (large at 200 m; pre-generate on a login node like the parent):
//   gmsh -3 tpv102_200m_hybrid_xsym.geo -o tpv102_200m_hybrid_xsym.msh -format msh22
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

// Transfinite fault resolution: cells per square fault half-edge.
// Both halves are Fault_half_length x Fault_width = 18 x 18 km squares, and
// lc_fault must divide Fault_width (18000 / 200 = 90 -> 200 m fault cells).
N_fault   = Fault_width / lc_fault;

// Outer domain (full half-space)
Xmax = 60e3;  Xmin = -Xmax;
Ymax = 60e3;  Ymin = -Ymax;
Zmin = -60e3;

// ---- Outer box (full domain) -------------------------------------------------
Box(1) = {Xmin, Ymin, Zmin,
          Xmax - Xmin, Ymax - Ymin, -Zmin};

// ---- Fault rectangle at y=0, SPLIT at x=0 into right/left halves -------------
// OCC Rectangle defaults to the xy-plane; build each half spanning y in [-FW,0],
// then rotate +90 deg about x: R_x(+pi/2): (x, y, 0) -> (x, 0, +y).
// Result: Surface{101}=RIGHT half (x in [0,FHL], z in [-FW,0]),
//         Surface{102}=LEFT  half (x in [-FHL,0], z in [-FW,0]); shared x=0 edge.
Rectangle(101) = { 0,                 -Fault_width, 0, Fault_half_length, Fault_width};
Rectangle(102) = {-Fault_half_length, -Fault_width, 0, Fault_half_length, Fault_width};
Rotate { {1, 0, 0}, {0, 0, 0}, Pi/2 } { Surface{101}; Surface{102}; }

// ---- Structured fault halves with MIRRORED diagonals (strike symmetry) -------
// Uniform N_fault cells per side; right "Right" (/) reflects to left "Left" (\).
Transfinite Curve { Boundary{ Surface{101}; } } = N_fault + 1;
Transfinite Curve { Boundary{ Surface{102}; } } = N_fault + 1;
Transfinite Surface { 101 } Right;
Transfinite Surface { 102 } Left;

// ---- Extrude each half into +/- y prism stacks (y-mirror, structured) --------
// Layers{N} -> geometry-level structured extrusion (5 layers in y, 200 m each).
// No Recombine -> triangles map to 3-tet-per-prism layers.  Each half is
// extruded +y and -y from the SAME source half -> +y/-y are exact y-mirrors,
// and the strike-symmetric source makes them strike-symmetric too.
rp[] = Extrude { 0,  h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };  // right +y
rm[] = Extrude { 0, -h_strip, 0 } { Surface{101}; Layers{N_strip_y}; };  // right -y
lp[] = Extrude { 0,  h_strip, 0 } { Surface{102}; Layers{N_strip_y}; };  // left  +y
lm[] = Extrude { 0, -h_strip, 0 } { Surface{102}; Layers{N_strip_y}; };  // left  -y
// rp[1], rm[1], lp[1], lm[1] are the four structured strip volumes.

// ---- Boolean fragment: embed the four strips into the outer box --------------
v_all[] = BooleanFragments
   { Volume{1}; Delete; }
   { Volume{rp[1], rm[1], lp[1], lm[1]}; Delete; };

// ---- Identify surfaces by bounding box (OCC reassigns tags after Boolean) ----
eps = 1.0;     // bounding-box slack [m]

// Fault halves at y=0 (each entirely within its x-range -> the two stay separate).
fault_right() = Surface In BoundingBox{
    -eps,                  -eps, -Fault_width - eps,
     Fault_half_length+eps,  eps,  0 + eps };
fault_left()  = Surface In BoundingBox{
    -Fault_half_length-eps, -eps, -Fault_width - eps,
     eps,                    eps,  0 + eps };

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

// ---- Re-assert transfinite on the POST-boolean fault halves (critical) -------
// BooleanFragments drops the transfinite attribute from the pre-boolean source;
// re-assert it on the final fault surfaces so the strike symmetry is preserved
// (without this the 3D mesher re-triangulates the fault -> ~28% strike symmetry).
Transfinite Curve { Boundary{ Surface{fault_right()}; } } = N_fault + 1;
Transfinite Curve { Boundary{ Surface{fault_left()};  } } = N_fault + 1;
Transfinite Surface { fault_right() } Right;
Transfinite Surface { fault_left()  } Left;

// ---- Mesh size control (matches parent: graded by distance from the fault) ----
// Field 1: Distance from BOTH fault halves (controls the unstructured outer bulk).
Field[1] = Distance;
Field[1].SurfacesList = {fault_right(), fault_left()};

// Field 2: graded size (lc_fault near fault, quadratic+linear growth far away).
//   F2(d) = 0.05*d + (d/2500)^2 + lc_fault
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);

// Field 3: sharp threshold cap at 2*lc_fault from the fault.
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
Physical Surface(3) = {fault_right(), fault_left()};
Physical Surface(5) = {abs_xmin(), abs_xmax(), abs_ymin(), abs_ymax(),
                       abs_zmin()};
Physical Volume (1) = {Volume{:}};

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm  = 6;       // Frontal-Delaunay 2D (unstructured outer bulk)
Mesh.Algorithm3D = 1;      // Delaunay 3D (reliable for mixed structured/unstructured)
