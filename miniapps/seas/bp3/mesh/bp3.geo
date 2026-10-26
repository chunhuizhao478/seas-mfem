// SEAS BP3 mesh: 2-D plane strain half-space with a dipping fault
// (2D in-plane design rev 3, phase B).
//
// Ported from Tandem, examples/tandem/2d/bp3.geo (Tandem a26753d6; the file's
// last commit c0ad90f9, 2020-12-16).  Geometry, tags and size fields are
// Tandem's; the one change is the default resolution: Lf = Ls = 0.025 km, the
// 25 m the BP3 description suggests (Tandem's default is 0.020).
//
// Geometry (km, y up, the free surface y = 0):
//   fault: from (0, 0) along the dip, x = -y cot(dip); rate and state for
//          down-dip distance 0..d4 = 40 km (Wf), split at d1, d2, d3 = 15, 16,
//          18 km so the mesh conforms to the VW/VS transition;
//   creep: the fault's continuation from 40 km down to the bottom y = -d;
//   box:   x in [-w, w] at the surface, w = max(d cot(dip), d), d = 400 km;
//          the sides are parallel to the fault.
//
// Physical groups (seas-mfem tags, rule R3; the same numbers as Tandem's):
//   Physical Curve 3   = fault, 0 to 40 km down-dip
//   Physical Curve 1   = natural: the free surface y = 0 and the bottom y = -d
//   Physical Curve 5   = Dirichlet: the two sides and the creep segment
//   Physical Surface 1 = hanging wall (x > -y cot(dip)), 2 = footwall
//
// Units: kilometers (deck mesh_scale = 1000).
//
// Usage:
//   gmsh -2 bp3.geo -format msh22 -setnumber dip 60 -o bp3_d60_25m.msh
//   gmsh -2 bp3.geo -format msh22 -setnumber dip 30 -o bp3_d30_25m.msh
//   gmsh -2 bp3.geo -format msh22 -setnumber dip 90 -o bp3_d90_25m.msh
//   gmsh -2 bp3.geo -format msh22 -setnumber Lf 0.05 -setnumber Ls 0.05 -o bp3_d60_50m.msh

DefineConstant[ Lf = {0.025, Min 0, Max 10, Name "Fault resolution" } ];
DefineConstant[ Ls = {0.025, Min 0, Max 10, Name "Resolution near free surface" } ];
DefineConstant[ dip = {60, Min 0, Max 90, Name "Dipping angle" } ];

dip_rad = dip * Pi / 180.0;

h = 40;
d = 400;
w = d * Cos(dip_rad) / Sin(dip_rad);
w = w < d ? d : w;
d1 = 15;
d2 = 16;
d3 = 18;
d4 = 40;
Point(1) = {-w, 0, 0, h};
Point(2) = {w, 0, 0, h};
Point(3) = {w + d * Cos(dip_rad) / Sin(dip_rad), -d, 0, h};
Point(4) = {-w + d * Cos(dip_rad) / Sin(dip_rad), -d, 0, h};
Point(5) = {d * Cos(dip_rad) / Sin(dip_rad), -d, 0, h};
Point(6) = {0, 0, 0, h};
Point(7) = {d4 * Cos(dip_rad), -d4 * Sin(dip_rad), 0, h};
Point(8) = {d3 * Cos(dip_rad), -d3 * Sin(dip_rad), 0, h};
Point(9) = {d2 * Cos(dip_rad), -d2 * Sin(dip_rad), 0, h};
Point(10) = {d1 * Cos(dip_rad), -d1 * Sin(dip_rad), 0, h};
Line(1) = {1, 6};
Line(2) = {6, 2};
Line(3) = {2, 3};
Line(4) = {3, 5};
Line(5) = {5, 4};
Line(6) = {4, 1};
Line(7) = {7, 5};
Line(8) = {8, 7};
Line(9) = {9, 8};
Line(10) = {10, 9};
Line(11) = {6, 10};
Curve Loop(1) = {11, 10, 9, 8, 7, -4, -3, -2};
Curve Loop(2) = {11, 10, 9, 8, 7, 5, 6, 1};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
// Tandem's "bottom: free surface" variant (its other two variants, bottom
// Dirichlet and all-fault, are not ported)
Physical Curve(3) = {8, 9, 10, 11};
Physical Curve(1) = {1, 2, 4, 5};
Physical Curve(5) = {3, 6, 7};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
// Size: Lf on the fault, growing with the squared distance from the fault
// line and from the 40 km depth; Ls at the surface trace, growing linearly.
Field[1] = MathEval;
Field[1].F = Sprintf("%g + 3e-2*(x + y*%g)^2 + 2e-3*(min(0, y/%g+40))^2",
                        Lf, Cos(dip_rad)/Sin(dip_rad), Sin(dip_rad));
Field[2] = MathEval;
Field[2].F = Sprintf("%g + 0.1 * sqrt(x^2 + y^2)", Ls);
Field[3] = Min;
Field[3].FieldsList = {1, 2};
Background Field = 3;
Mesh.MshFileVersion = 2.2;
