// The bend convergence test's meshes (Chunhui Zhao, 2026-10-30: "prepare
// the bend convergence test"): one 30 km fault with a 5 deg bend in a
// +-120 km box, plane strain, map view.  Is the drop of sigma_n' per metre
// of slip at a bend set by the elasticity or by the mesh size?
//
// mapview_rs/mesh/mapview_rs.geo (the earlier bend runs' base) with:
//   box   x, y in [-120, 120] km (was +-60 km);
//   fault the arc length s in [-15, 15] km, the point s = 0 at (0, 0):
//           s < -arc / 2         on y = 0 (x = s);
//           |s| <= arc / 2       a circular arc of radius R = arc / theta,
//                                turning left by theta in all;
//           s > arc / 2          straight again, at theta to the x axis;
//         arc = 0 is the sharp bend at (0, 0).  Fault points every lc_fault
//         in s, ON the curve, so the facet angle theta lc_fault / arc falls
//         as the mesh is refined (the earlier arc meshes had 500 m facets
//         at every h).  The arc ends (s = +-2.5 and +-5 km) and the vertex
//         fall on points at lc_fault = 500, 250, 125 m.
// No split of the fault where a changes: in the decks a is a linear taper
// in x between 10 and 12 km, continuous.
//
// Physical groups (as mapview_rs.geo): 3 the fault, 5 the four sides, 1 the
// bulk.  Size: lc_fault on the fault, 0.05 d + (d / 2.5 km)^2 + lc_fault
// away from it (d the distance to the fault), capped at lc.  Meters.
//
// Usage (from this folder; the committed meshes, gmsh 4.15.2):
//   gmsh -2 bend.geo -format msh22 -setnumber theta_deg 5 -setnumber arc 0     -setnumber lc_fault 500 -o sharp05_h500.msh
//   gmsh -2 bend.geo -format msh22 -setnumber theta_deg 5 -setnumber arc 5000  -setnumber lc_fault 250 -o arc05L05_h250.msh
//   gmsh -2 bend.geo -format msh22 -setnumber theta_deg 5 -setnumber arc 10000 -setnumber lc_fault 125 -o arc05L10_h125.msh
// (the nine: sharp05, arc05L05, arc05L10 at lc_fault 500, 250, 125).
//
// Two more inputs for the simple cases (Chunhui Zhao, 2026-10-31: "go ahead
// with the simple cases"; simple_cases/README.md), defaults = the nine above:
//   box     the box half-width [m] (120 km);
//   s_bend  the arc length s of the bend centre [m] (0): the fault keeps
//           its x in [-15, 15] km on the straight part, and the vertex (or
//           arc centre) moves to x = s_bend on y = 0, e.g. s_bend = 12500
//           puts a sharp bend in the creeping part (a = 0.04 from 12 km).
//   gmsh -2 bend.geo -format msh22 -setnumber box 240e3 -setnumber theta_deg 1.5 -setnumber arc 0 -setnumber s_bend 12500 -setnumber lc_fault 500 -o sharp015vs_h500.msh

DefineConstant[ lc_fault = {500, Min 1, Max 5000, Name "Fault resolution [m]" } ];
DefineConstant[ lc = {5e3, Min 100, Max 50e3, Name "Far-field resolution [m]" } ];
DefineConstant[ theta_deg = {5, Min 0, Max 45, Name "Bend angle [deg]" } ];
DefineConstant[ arc = {0, Min 0, Max 30e3, Name "Arc length of the bend [m], 0 = sharp" } ];
DefineConstant[ box = {120e3, Min 30e3, Max 1000e3, Name "Box half-width [m]" } ];
DefineConstant[ s_bend = {0, Min -15e3, Max 15e3, Name "Arc-length position of the bend centre [m]" } ];

Fault_half_length = 15e3;
Xmax = box;
Ymax = box;
th = theta_deg * Pi / 180;

// box
Point(1) = {-Xmax, -Ymax, 0, lc};
Point(2) = { Xmax, -Ymax, 0, lc};
Point(3) = { Xmax,  Ymax, 0, lc};
Point(4) = {-Xmax,  Ymax, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// fault, s from -15 to 15 km in steps of lc_fault
n = Round(2 * Fault_half_length / lc_fault);
s0 = -arc / 2;
For i In {0 : n}
   s = -Fault_half_length + i * lc_fault - s_bend;
   If (arc == 0)
      If (s <= 0)
         x = s; y = 0;
      Else
         x = s * Cos(th); y = s * Sin(th);
      EndIf
   Else
      R = arc / th;
      If (s <= s0)
         x = s; y = 0;
      Else
         If (s <= -s0)
            p = (s - s0) / R;
            x = s0 + R * Sin(p); y = R * (1 - Cos(p));
         Else
            x = s0 + R * Sin(th) + (s + s0) * Cos(th);
            y = R * (1 - Cos(th)) + (s + s0) * Sin(th);
         EndIf
      EndIf
   EndIf
   x = x + s_bend;
   Point(1000 + i) = {x, y, 0, lc_fault};
EndFor
For i In {0 : n - 1}
   Line(1000 + i) = {1000 + i, 1001 + i};
EndFor
fault[] = {1000 : 1000 + n - 1};
Curve{fault[]} In Surface{1};

Field[1] = Distance;
Field[1].CurvesList = {fault[]};
Field[1].Sampling = 2000;
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 + (F1/2.5e3)^2 + %g", lc_fault);
Field[3] = MathEval;
Field[3].F = Sprintf("%g", lc);
Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Physical Curve(3) = {fault[]};
Physical Curve(5) = {1, 2, 3, 4};
Physical Surface(1) = {1};
Mesh.MshFileVersion = 2.2;
