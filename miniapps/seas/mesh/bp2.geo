// BP2 mesh: full domain, vertical antiplane fault at x=0
//
// Domain: x in [-D, +D] km, z in [-D, 0] km
// Fault:  x=0, z=0 to z=-Wf (interior interface between two surfaces)
//
// Units: kilometers (scale by 1000 when loading into MFEM)
//
// Usage:
//   gmsh -2 bp2.geo -o bp2.msh
//   gmsh -2 bp2.geo -o bp2.msh -setnumber hf 0.100   # 100m resolution
//   gmsh -2 bp2.geo -o bp2.msh -setnumber hf 0.050   # 50m resolution
//   gmsh -2 bp2.geo -o bp2.msh -setnumber h 100       # 100km far-field
//
// Physical Curve IDs match BP2BoundaryAttributes in C++ code:
//   1 = FARFIELD_LEFT  (x = -D)
//   2 = FARFIELD_RIGHT (x = +D)
//   3 = FREE_SURFACE   (z = 0)
//   4 = BOTTOM         (z = -D)

DefineConstant[ hf = {0.200, Min 0.01, Max 10, Name "Fault resolution [km]"} ];
DefineConstant[ h  = {50.0,  Min 1,    Max 200, Name "Far boundary resolution [km]"} ];
DefineConstant[ D  = {400,   Min 10,   Max 10000, Name "Domain half-size [km]"} ];

// BP2 fault zone depths [km]
d1 = 15;   // VW/VS transition (H)
d2 = 16;   // Mesh control point
d3 = 18;   // H + h (end of transition zone)
d4 = 40;   // Wf (bottom of rate-state fault)

// ==== Points ====

// Fault line (x=0): fine resolution
Point(1)  = {0, 0, 0, hf};      // Surface/fault intersection
Point(2)  = {0, -d1, 0, hf};    // VW/VS transition
Point(3)  = {0, -d2, 0, hf};    // Control point
Point(4)  = {0, -d3, 0, hf};    // End of transition
Point(5)  = {0, -d4, 0, hf};    // Bottom of rate-state fault (Wf)
Point(6)  = {0, -D, 0, h};      // Bottom of domain at x=0

// Right boundary (x = +D): coarse
Point(7)  = {D, 0, 0, h};       // Top-right
Point(8)  = {D, -D, 0, h};      // Bottom-right

// Left boundary (x = -D): coarse
Point(9)  = {-D, 0, 0, h};      // Top-left
Point(10) = {-D, -D, 0, h};     // Bottom-left

// ==== Lines ====

// Top boundary (z = 0)
Line(1) = {9, 1};     // Top-left to fault
Line(2) = {1, 7};     // Fault to top-right

// Right boundary (x = +D)
Line(3) = {7, 8};

// Bottom boundary (z = -D)
Line(4) = {8, 6};     // Bottom-right to fault-bottom
Line(5) = {6, 10};    // Fault-bottom to bottom-left

// Left boundary (x = -D)
Line(6) = {10, 9};

// Fault segments (x = 0, top to bottom)
Line(7)  = {1, 2};    // Surface to d1
Line(8)  = {2, 3};    // d1 to d2
Line(9)  = {3, 4};    // d2 to d3
Line(10) = {4, 5};    // d3 to Wf
Line(11) = {5, 6};    // Wf to bottom

// Force exact hf spacing on the fault (lines 7-10)
Transfinite Curve{7}  = d1 / hf + 1;
Transfinite Curve{8}  = (d2 - d1) / hf + 1;
Transfinite Curve{9}  = (d3 - d2) / hf + 1;
Transfinite Curve{10} = (d4 - d3) / hf + 1;

// ==== Surfaces ====

// Right half (x > 0)
Curve Loop(1) = {2, 3, 4, -11, -10, -9, -8, -7};
Plane Surface(1) = {1};

// Left half (x < 0)
Curve Loop(2) = {7, 8, 9, 10, 11, 5, 6, 1};
Plane Surface(2) = {2};

// ==== Physical groups ====

// Boundary attributes (must match BP2BoundaryAttributes in C++ code)
Physical Curve(1) = {6};        // FARFIELD_LEFT  (x = -D)
Physical Curve(2) = {3};        // FARFIELD_RIGHT (x = +D)
Physical Curve(3) = {1, 2};     // FREE_SURFACE   (z = 0)
Physical Curve(4) = {4, 5};     // BOTTOM         (z = -D)

// Two physical surfaces
Physical Surface(1) = {1};      // Right half
Physical Surface(2) = {2};      // Left half

Mesh.MshFileVersion = 2.2;
