// =============================================================
// safs_fault_box_nwcut.geo
//
// Free-surface variant for the *NW-cut* alternative ALT6 fault.
//
// Source fault: ../data_cutnwfault/<filename>.stl, produced by
// code_preprocess/nw_cut_strip.py.  The pipeline drops every triangle
// that lies NW of the preferred MJVS SAF mesh's NW-most surface vertex,
// so the surviving fault is the SE ~283 km of the original ALT6 strip.
//
// The bounding box is no longer fitted to the *uncut* alt strip's
// bbox (443 km wide x 296 km tall) but to the *cut* strip's bbox
// (~258 km x ~147 km), padded by PAD_XY laterally and PAD_BOTTOM
// vertically.  All bbox values are injected from
// run_nwcut_meshing.py via gmsh -setnumber so the same template
// handles 500 / 1000 / 2000 m resolutions.
//
// Free-surface treatment: mirror of safs_fault_box_freesurface_clip.geo.
// PAD_TOP = 100 m of overburden keeps the fault top edge strictly
// interior to the box top face; the proven `Surface{} In Volume{}`
// embedding therefore works without ClassifySurfaces.
//
// Generate (typical invocation; see run_nwcut_meshing.py for the actual
// command-line plumbing):
//   gmsh -3 safs_fault_box_nwcut.geo \
//        -setstring fault_stl ../data_cutnwfault/<file>.stl \
//        -setnumber xmin_fault <X_MIN> -setnumber xmax_fault <X_MAX> \
//        -setnumber ymin_fault <Y_MIN> -setnumber ymax_fault <Y_MAX> \
//        -setnumber zmin_fault <Z_MIN> -setnumber zmax_fault <Z_MAX> \
//        -o safs_fault_box_nwcut_<RES>m.msh
//
// All units: meters.  Coordinates: UTM 11 N (matches the .ts source).
// =============================================================

SetFactory("OpenCASCADE");

// -------------------------------------------------------------
// Required parameters (must be provided via -setnumber / -setstring).
// Values are computed by run_nwcut_meshing.py from the cut STL's bbox
// at runtime, so the .geo template stays bbox-agnostic.
// -------------------------------------------------------------
If (!Exists(xmin_fault) || !Exists(xmax_fault) ||
    !Exists(ymin_fault) || !Exists(ymax_fault) ||
    !Exists(zmin_fault) || !Exists(zmax_fault))
    Error("safs_fault_box_nwcut.geo requires {x,y,z}{min,max}_fault to be set via -setnumber");
EndIf
If (!Exists(fault_stl))
    Error("safs_fault_box_nwcut.geo requires fault_stl to be set via -setstring");
EndIf

// Domain padding -- matched to safs_fault_box_freesurface_clip.geo.
// Three layers of overrides:
//   pad_xy     = legacy single horizontal pad applied to all 4 sides
//   pad_x / pad_y = symmetric per-axis pads (override pad_xy on that axis)
//   pad_x_lo, pad_x_hi, pad_y_lo, pad_y_hi = asymmetric per-side pads
//                  (override pad_x / pad_y on that side)
// pad_top default = 100 m (back-compat); pad_bottom default = 25 km.
PAD_XY     = 50000.0;   // 50 km horizontal margin (legacy default)
PAD_TOP_DEFAULT    = 100.0;       // 100 m of rock above the free surface
                                  // (thin overburden -- keeps fault strictly
                                  //  interior so `Surface In Volume` is valid)
PAD_BOTTOM_DEFAULT = 25000.0;     // 25 km of rock below deepest fault vert
If (!Exists(pad_x))       pad_x       = PAD_XY;              EndIf
If (!Exists(pad_y))       pad_y       = PAD_XY;              EndIf
If (!Exists(pad_x_lo))    pad_x_lo    = pad_x;               EndIf
If (!Exists(pad_x_hi))    pad_x_hi    = pad_x;               EndIf
If (!Exists(pad_y_lo))    pad_y_lo    = pad_y;               EndIf
If (!Exists(pad_y_hi))    pad_y_hi    = pad_y;               EndIf
If (!Exists(pad_bottom))  pad_bottom  = PAD_BOTTOM_DEFAULT;   EndIf
If (!Exists(pad_top))     pad_top     = PAD_TOP_DEFAULT;      EndIf

xmin = xmin_fault - pad_x_lo;
xmax = xmax_fault + pad_x_hi;
ymin = ymin_fault - pad_y_lo;
ymax = ymax_fault + pad_y_hi;
zmin = zmin_fault - pad_bottom;
zmax = zmax_fault + pad_top;

dx = xmax - xmin;
dy = ymax - ymin;
dz = zmax - zmin;

Printf("nwcut bbox: x=[%g,%g] y=[%g,%g] z=[%g,%g]  (%g x %g x %g km)",
       xmin, xmax, ymin, ymax, zmin, zmax,
       dx/1000, dy/1000, dz/1000);

// -------------------------------------------------------------
// Mesh size controls -- match safs_fault_box_freesurface_clip.geo.
// Allow override per-resolution via -setnumber lc_near, lc_far,
// dist_inner, dist_outer; otherwise use the existing project defaults.
// -------------------------------------------------------------
If (!Exists(lc_near))     lc_near     =  1500.0; EndIf
If (!Exists(lc_far))      lc_far      = 10000.0; EndIf
If (!Exists(dist_inner))  dist_inner  =  3000.0; EndIf
If (!Exists(dist_outer))  dist_outer  = 40000.0; EndIf

LC_NEAR    = lc_near;
LC_FAR     = lc_far;
DIST_INNER = dist_inner;
DIST_OUTER = dist_outer;

Mesh.CharacteristicLengthMin = LC_NEAR;
Mesh.CharacteristicLengthMax = LC_FAR;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;

// -------------------------------------------------------------
// Box volume.  OCC tags 1..6 for box faces; volume tag 1.
// -------------------------------------------------------------
Box(1) = {xmin, ymin, zmin, dx, dy, dz};

// -------------------------------------------------------------
// Merge the cut fault and embed it.  Discrete fault surface
// arrives at tag 7 (one past the six box faces).  Do NOT call
// ClassifySurfaces (it splits the fault into multiple discrete
// patches whose shared edges are not stitched into a closed
// manifold, so the 3-D mesher reports "No closed volume").
// -------------------------------------------------------------
Merge Str(fault_stl);

fault_tag = 7;
fault_surfs[] = {fault_tag};

Surface{fault_tag} In Volume{1};

// -------------------------------------------------------------
// Distance-based size field: refine near the fault.
// -------------------------------------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {fault_surfs[]};
Field[1].Sampling     = 100;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = LC_NEAR;
Field[2].SizeMax  = LC_FAR;
Field[2].DistMin  = DIST_INNER;
Field[2].DistMax  = DIST_OUTER;

// -------------------------------------------------------------
// Optional sidecar-driven background size field (G-1 / PLAN.md
// Phase 2).  When USE_SIZE_FIELD == 1, the .pos referenced by
// `size_field_pos` is merged in as PostView 0 and combined with
// the distance-to-fault Field[2] via Field[Min].
// -------------------------------------------------------------
If (!Exists(USE_SIZE_FIELD)) USE_SIZE_FIELD = 0; EndIf

// -------------------------------------------------------------
// Optional z-graded background size field.  When USE_Z_GRADED == 1,
// LC is also capped from above by a depth-dependent piecewise-
// CONSTANT step function (NOT piecewise-linear; see REVIEW R-010).
// Implementation: four `Box` fields combined with `Field[Min]`.  Each
// Box reports `VIn` inside its z-slab and `VOut = 1e22` outside, so
// LC is constant within each slab and JUMPS discontinuously at the
// slab boundaries.  Sharp jumps in target LC can produce sliver
// tets at the boundaries; if that becomes an issue, replace the
// `Box` fields with a smooth `MathEval` formulation (gmsh's
// MathEval rejects C-style ternary, so use a `sign()`/`min()`
// combination).
//
// "推荐" defaults (override via -setnumber if needed):
//   z >  -500  m:  LC = 500   (basin sediment / surface)
//   z in [-3 km, -500 m]: LC = 1000  (basin body)
//   z in [-10 km, -3 km]: LC = 2000  (transition to basement)
//   z < -10 km:    LC = 3000  (basement)
// Boundary behaviour: at z = -500, -3000, -10000 exactly, two
// adjacent Boxes both report `VIn` (their `Min`/`Max` are inclusive);
// `Field[Min]` picks the smaller, so the boundary cell adopts the
// shallower-layer LC.  Combined with Field[2] (distance-to-fault)
// via `Min` so the fault corridor's lc_near always wins where smaller.
// -------------------------------------------------------------
If (!Exists(USE_Z_GRADED)) USE_Z_GRADED = 0; EndIf
If (!Exists(z_lc_top))      z_lc_top      =  500.0; EndIf
If (!Exists(z_lc_basin))    z_lc_basin    = 1000.0; EndIf
If (!Exists(z_lc_trans))    z_lc_trans    = 2000.0; EndIf
If (!Exists(z_lc_deep))     z_lc_deep     = 3000.0; EndIf

If (USE_SIZE_FIELD == 1)
    If (!Exists(size_field_pos))
        Error("USE_SIZE_FIELD=1 requires -setstring size_field_pos PATH");
    EndIf
    Merge Str(size_field_pos);
    Field[3] = PostView;
    Field[3].ViewIndex = 0;

    Field[4] = Min;
    Field[4].FieldsList = {2, 3};
    bg_field_id = 4;
Else
    bg_field_id = 2;
EndIf

If (USE_Z_GRADED == 1)
    // gmsh MathEval lacks C-style ternary, so build the piecewise
    // depth-dependent LC out of 4 Box fields and take the per-point
    // minimum.  Each Box reports VIn inside its slab, VOut=1e22
    // outside.  Min() naturally picks the slab that contains the
    // query point.
    //
    // We need horizontal extents that cover the whole mesh; padding
    // them with 1e6 m gives a safe margin without inflating the
    // bounding-box check (Box fields only test scalar inclusion).
    box_xlo = xmin - 1.0e6;  box_xhi = xmax + 1.0e6;
    box_ylo = ymin - 1.0e6;  box_yhi = ymax + 1.0e6;
    Field[5] = Box;
    Field[5].VIn  = z_lc_top;        Field[5].VOut = 1.0e22;
    Field[5].XMin = box_xlo;         Field[5].XMax = box_xhi;
    Field[5].YMin = box_ylo;         Field[5].YMax = box_yhi;
    Field[5].ZMin =  -500.0;         Field[5].ZMax = 1.0e6;
    Field[6] = Box;
    Field[6].VIn  = z_lc_basin;      Field[6].VOut = 1.0e22;
    Field[6].XMin = box_xlo;         Field[6].XMax = box_xhi;
    Field[6].YMin = box_ylo;         Field[6].YMax = box_yhi;
    Field[6].ZMin = -3000.0;         Field[6].ZMax = -500.0;
    Field[7] = Box;
    Field[7].VIn  = z_lc_trans;      Field[7].VOut = 1.0e22;
    Field[7].XMin = box_xlo;         Field[7].XMax = box_xhi;
    Field[7].YMin = box_ylo;         Field[7].YMax = box_yhi;
    Field[7].ZMin = -10000.0;        Field[7].ZMax = -3000.0;
    Field[8] = Box;
    Field[8].VIn  = z_lc_deep;       Field[8].VOut = 1.0e22;
    Field[8].XMin = box_xlo;         Field[8].XMax = box_xhi;
    Field[8].YMin = box_ylo;         Field[8].YMax = box_yhi;
    Field[8].ZMin = -1.0e6;          Field[8].ZMax = -10000.0;
    Field[9] = Min;
    Field[9].FieldsList = {bg_field_id, 5, 6, 7, 8};
    Background Field = 9;
Else
    Background Field = bg_field_id;
EndIf

// -------------------------------------------------------------
// Physical groups.
//   - Volume  "rock"     : the box volume
//   - Surface "fault"    : the embedded discrete fault
//   - Surface "top"      : top face of the box (z = zmax) -- the
//                          approximate free surface
//   - Surface "bottom"   : bottom face (z = zmin)
//   - Surface "sides"    : the four lateral faces
// -------------------------------------------------------------
Physical Volume("rock", 1)      = {1};
Physical Surface("fault", 101)  = fault_surfs[];
Physical Surface("top", 102)    = {6};
Physical Surface("bottom", 103) = {5};
Physical Surface("sides", 104)  = {1, 2, 3, 4};

// -------------------------------------------------------------
// Mesh algorithms.  Algorithm3D = 1 (Delaunay) handles embedded
// discrete surfaces robustly; HXT (Algorithm3D = 10) reports
// "No closed volume" on this geometry.
// -------------------------------------------------------------
Mesh.Algorithm   = 6;    // 2-D: Frontal-Delaunay
Mesh.Algorithm3D = 1;    // 3-D: Delaunay
