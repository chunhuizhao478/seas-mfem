// SCEC TPV31 mesh — half-space model with planar strike-slip fault.
// Reference: tpv31/benchmark_document/TPV31_32_Description_v03.pdf
//             plan §R.5 step 2 (PLAN_phase_R_exact_bimaterial_riemann_rev3).
//
// CANONICAL SEAS COORDINATE FRAME (see miniapps/seas/CLAUDE.md
// "Canonical Coordinate System"; matches TPV102, TPV104, TPV205):
//   x = along-strike (east)
//   y = fault-normal
//   z = vertical (z = 0 at free surface; z < 0 below)
//
// The TPV31 spec uses a different frame (fault on z = 0 with y = depth);
// this mesh applies the rotation (x_s, y_s, z_s) → (x_s, z_s, -y_s) so
// that TPV31 runs through `seas_spatial_dyn_driver` with the same
// FaultBasis (ref_normal = (0, -1, 0), up = (0, 0, 1)) as the other
// vertical-y=0 benchmarks.  See
// `debug_document/general_driver_debug_document/REVIEW.md` R-003 and
// the TPV31 rotation follow-up.
//
// Geometry (canonical):
//   * Half-space domain:  x ∈ [-50, 50] km,
//                          y ∈ [-50, 50] km   (fault-normal extent),
//                          z ∈ [-50, 0]  km   (depth: 0 = surface).
//   * Fault:              y = 0 plane,  x ∈ [-15, 15] km,
//                          z ∈ [-15, 0]  km   (spec depth ∈ [0, 15] km).
//   * Free surface:       z = 0 plane.
//   * Absorbing walls:    the four vertical sides at x = ±50 km,
//                          y = ±50 km, and the bottom at z = -50 km.
//
// Reflection-time budget (plan §R.5 step 2 REVIEW R-001):
//   tfinal = 15 s, c_p_deepest = 6500 m/s ⇒ minimum reflection-free
//   distance d ≥ c_p · tfinal / 2 = 48.75 km.  The 50 km box clears
//   this by 1.25 km (margin ≈ 2.6%).  A future tighter run could use
//   a 30 km box + 5 km PML.
//
// Mesh refinement (per plan §R.5 step 2):
//   * 50 m near the fault (lc_fault).
//   * Ramps to 500 m at the absorbing boundary (lc_far).
//
// Boundary attribute conventions (REVIEW R-014; must match
// tpv31/configs/tpv31.toml [boundary] block):
//   Physical Surface(101) = fault plane (y = 0, x ∈ [-15,15], z ∈ [-15,0])
//   Physical Surface(102) = free surface (z = 0)
//   Physical Surface(103) = absorbing — four vertical walls (x = ±d, y = ±d)
//   Physical Surface(104) = absorbing — bottom (z = -box_extent)
//
// Build:  gmsh -format msh22 -3 tpv31_50m.geo -o tpv31_50m.msh
// (CLAUDE.md "Known limitation — Gmsh .msh format" — must use msh22,
// not the default v4 format.)

// ----------------------------------------------------------------------
// Mesh-size parameters
// ----------------------------------------------------------------------
lc_fault = 50.0;            // 50 m near the fault (spec p. 11)
lc_far   = 500.0;           // 500 m at the absorbing boundary

// ----------------------------------------------------------------------
// Domain box  (x, y, z) in metres — canonical SEAS frame.
// ----------------------------------------------------------------------
xMin = -50000.0;  xMax =  50000.0;
yMin = -50000.0;  yMax =  50000.0;
zMin = -50000.0;  zMax =      0.0;   // z = 0 free surface, z < 0 below

// ----------------------------------------------------------------------
// Fault rectangle (y = 0 plane), x ∈ [-15, 15] km, z ∈ [-15, 0] km.
// ----------------------------------------------------------------------
fx0 = -15000.0;  fx1 =  15000.0;
fz0 = -15000.0;  fz1 =      0.0;     // depth 15 km down to surface

// ----------------------------------------------------------------------
// Vertices: 8 corners of the box.  Layout (looking from +y):
//   z = zMax (= 0, free surface)
//     1 (-x,-y) -- 2 (+x,-y)
//     |              |
//     4 (-x,+y) -- 3 (+x,+y)
//   z = zMin (= -50, bottom)
//     5 (-x,-y) -- 6 (+x,-y)
//     |              |
//     8 (-x,+y) -- 7 (+x,+y)
// ----------------------------------------------------------------------
Point(1) = {xMin, yMin, zMax, lc_far};
Point(2) = {xMax, yMin, zMax, lc_far};
Point(3) = {xMax, yMax, zMax, lc_far};
Point(4) = {xMin, yMax, zMax, lc_far};
Point(5) = {xMin, yMin, zMin, lc_far};
Point(6) = {xMax, yMin, zMin, lc_far};
Point(7) = {xMax, yMax, zMin, lc_far};
Point(8) = {xMin, yMax, zMin, lc_far};

// Fault rectangle corners on the y = 0 plane.
Point(11) = {fx0, 0.0, fz0, lc_fault};
Point(12) = {fx1, 0.0, fz0, lc_fault};
Point(13) = {fx1, 0.0, fz1, lc_fault};
Point(14) = {fx0, 0.0, fz1, lc_fault};

// ----------------------------------------------------------------------
// Edges of the bounding box.
// ----------------------------------------------------------------------
// Top face (z = 0 = free surface):
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Bottom face (z = zMin):
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
// Vertical edges (top → bottom, +z → -z):
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// ----------------------------------------------------------------------
// Fault rectangle edges (y = 0 plane).
// ----------------------------------------------------------------------
Line(21) = {11, 12};
Line(22) = {12, 13};
Line(23) = {13, 14};
Line(24) = {14, 11};

Curve Loop(101) = {21, 22, 23, 24};
Plane Surface(1) = {101};      // Fault rectangle

// ----------------------------------------------------------------------
// Bounding-box faces.
// ----------------------------------------------------------------------
// z = zMax = 0 (free surface).
Curve Loop(102) = {1, 2, 3, 4};
Plane Surface(2) = {102};

// z = zMin (absorbing bottom).
Curve Loop(103) = {5, 6, 7, 8};
Plane Surface(3) = {103};

// x = xMin (absorbing wall):  edges 4 (top), 9 (down), -8 (bottom rev), -12 (up rev)
Curve Loop(104) = {4, 9, -8, -12};
Plane Surface(4) = {104};

// x = xMax:  edges 2 (top), 11 (down), -6 (bottom rev), -10 (up rev)
Curve Loop(105) = {2, 11, -6, -10};
Plane Surface(5) = {105};

// y = yMin:  edges 1 (top), 10 (down), -5 (bottom rev), -9 (up rev)
Curve Loop(106) = {1, 10, -5, -9};
Plane Surface(6) = {106};

// y = yMax:  edges 3 (top), 12 (down), -7 (bottom rev), -11 (up rev)
Curve Loop(107) = {3, 12, -7, -11};
Plane Surface(7) = {107};

// ----------------------------------------------------------------------
// Volume — embed the fault surface so the mesher splits tets along it.
// ----------------------------------------------------------------------
Surface Loop(1) = {2, 3, 4, 5, 6, 7};
Volume(1) = {1};

// The fault reaches the free surface: its top edge (Line 23, the
// y = 0, z = 0 segment x ∈ [-15,15] km) lies IN the free-surface plane
// (Surface 2 = z = 0 box top).  Embed that edge in Surface 2 so the
// free-surface 2D mesh conforms to the fault trace.  Without this the
// 3D mesher hits a segment-facet intersection along z = 0 and produces
// "No elements in volume 1".  Mirrors the `Line{101} In Surface{1}`
// trick in tpv205/mesh/tpv2053d_200m.geo (also a surface-reaching fault).
Curve{23} In Surface{2};

Surface{1} In Volume{1};

// ----------------------------------------------------------------------
// Mesh-size field: refined near the fault rectangle, coarsening
// outward.
// ----------------------------------------------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {1};

Field[2] = MathEval;
Field[2].F = Sprintf("%g + %g*F1", lc_fault, (lc_far - lc_fault) / 5000.0);

Background Field = 2;

// Force the surface MeshSize to lc_fault at the fault patch.
Characteristic Length{ PointsOf{ Surface{1}; } } = lc_fault;

// ----------------------------------------------------------------------
// Physical groups (REVIEW R-014).  Attribute IDs must match
// tpv31/configs/tpv31.toml [boundary] block.
// ----------------------------------------------------------------------
Physical Surface(101) = {1};            // fault (y = 0)
Physical Surface(102) = {2};            // free surface (z = 0)
Physical Surface(103) = {4, 5, 6, 7};   // absorbing vertical walls
Physical Surface(104) = {3};            // absorbing bottom (z = zMin)
Physical Volume(1)   = {1};

// ----------------------------------------------------------------------
// Generator hints (msh22 required — see CLAUDE.md).
// ----------------------------------------------------------------------
Mesh.MshFileVersion  = 2.2;
Mesh.Algorithm3D     = 1;     // Delaunay
// NOTE: do NOT enable Mesh.OptimizeNetgen here.  The standard gmsh
// optimizer already drives this mesh to "No ill-shaped tets" quality;
// the extra Netgen pass (OptimizeNetgen=1) aborts with SIGABRT on this
// surface-reaching-fault geometry ("illegal tets / badmax = 1e+24" in
// SwapImprove).  The working tpv205/tpv102/tpv104 meshes likewise do
// not set it.  Mesh.Optimize (standard) stays on by default.
