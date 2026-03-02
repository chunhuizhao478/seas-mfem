// BP1 self-similar mesh: uniform near fault, geometric coarsening in far-field
//
// Domain: x in [-D, +D] km, z in [-D, 0] km
// Fault:  x=0 (interior interface between two surfaces)
// Uniform zone: x in [-W, +W], z in [-Zf, 0]
// Coarsening zone: geometric progression with ratio r elsewhere
//
// Units: kilometers (scale by 1000 when loading into MFEM)
//
// Usage:
//   gmsh -2 bp1_selfsimilar.geo -o bp1_ss_100m.msh -setnumber h 0.1
//   gmsh -2 bp1_selfsimilar.geo -o bp1_ss_50m.msh  -setnumber h 0.05
//   gmsh -2 bp1_selfsimilar.geo -o bp1_ss_25m.msh  -setnumber h 0.025
//
// Physical Curve IDs match SEASBoundaryTags in C++ code:
// 5 = FAULT (x = 0, interior interface)
//   1 = FARFIELD_LEFT  (x = -D)
//   2 = FARFIELD_RIGHT (x = +D)
//   3 = FREE_SURFACE   (z = 0)
//   4 = BOTTOM         (z = -D)
//
// Block layout (8 surfaces, split at x=0, x=+-W, z=-Zf):
//
//  z=0    ___________________________________
//        |        |      |      |        |
//        |  S1    |  S2  |  S3  |   S4   |
//        | coarse | unif | unif | coarse |
//        |  in x  |      |      |  in x  |
//  z=-Zf |________|______|______|________|
//        |        |      |      |        |
//        |  S5    |  S6  |  S7  |   S8   |
//        | coarse | coar | coar | coarse |
//        | x & z  | in z | in z | x & z  |
//  z=-D  |________|______|______|________|
//       x=-D    x=-W   x=0   x=W      x=D
//
// Estimated element counts:
//   h=0.1  (100m): ~94k quads
//   h=0.05  (50m): ~289k quads
//   h=0.025 (25m): ~963k quads

DefineConstant[ h  = {0.1,  Min 0.01, Max 10,  Name "Uniform resolution [km]"} ];
DefineConstant[ D  = {80,   Min 10,   Max 1000, Name "Domain half-size [km]"} ];
DefineConstant[ W  = {5,    Min 1,    Max 50,  Name "Uniform zone half-width [km]"} ];
DefineConstant[ Zf = {45,   Min 10,   Max 100, Name "Uniform zone depth [km]"} ];
DefineConstant[ r  = {1.1,  Min 1.01, Max 2,   Name "Coarsening ratio"} ];

// Computed node counts
nx_uniform = W / h;
nz_uniform = Zf / h;
nx_coarsen = Ceil(Log(1 + (D - W) * (r - 1) / h) / Log(r));
nz_coarsen = Ceil(Log(1 + (D - Zf) * (r - 1) / h) / Log(r));

// ==== Points (15 corners on a 5x3 grid) ====

// Row z = 0 (top)
Point(1)  = {-D,  0, 0};
Point(2)  = {-W,  0, 0};
Point(3)  = { 0,  0, 0};
Point(4)  = { W,  0, 0};
Point(5)  = { D,  0, 0};

// Row z = -Zf (middle)
Point(6)  = {-D, -Zf, 0};
Point(7)  = {-W, -Zf, 0};
Point(8)  = { 0, -Zf, 0};
Point(9)  = { W, -Zf, 0};
Point(10) = { D, -Zf, 0};

// Row z = -D (bottom)
Point(11) = {-D, -D, 0};
Point(12) = {-W, -D, 0};
Point(13) = { 0, -D, 0};
Point(14) = { W, -D, 0};
Point(15) = { D, -D, 0};

// ==== Lines ====

// --- Horizontal lines (left to right) ---

// Row z = 0
Line(1)  = {1, 2};    // -D to -W (coarsening)
Line(2)  = {2, 3};    // -W to 0  (uniform)
Line(3)  = {3, 4};    //  0 to W  (uniform)
Line(4)  = {4, 5};    //  W to D  (coarsening)

// Row z = -Zf
Line(5)  = {6, 7};
Line(6)  = {7, 8};
Line(7)  = {8, 9};
Line(8)  = {9, 10};

// Row z = -D
Line(9)  = {11, 12};
Line(10) = {12, 13};
Line(11) = {13, 14};
Line(12) = {14, 15};

// --- Vertical lines (top to bottom) ---

// Column x = -D
Line(13) = {1, 6};    // 0 to -Zf (uniform in z)
Line(14) = {6, 11};   // -Zf to -D (coarsening in z)

// Column x = -W
Line(15) = {2, 7};
Line(16) = {7, 12};

// Column x = 0 (fault)
Line(17) = {3, 8};
Line(18) = {8, 13};

// Column x = W
Line(19) = {4, 9};
Line(20) = {9, 14};

// Column x = D
Line(21) = {5, 10};
Line(22) = {10, 15};

// ==== Transfinite Curves ====

// Horizontal coarsening: -D to -W (coarse → fine, shrinking elements)
Transfinite Curve{1, 5, 9}   = nx_coarsen + 1 Using Progression 1/r;

// Horizontal uniform: -W to 0 and 0 to W
Transfinite Curve{2, 3, 6, 7, 10, 11} = nx_uniform + 1;

// Horizontal coarsening: W to D (fine → coarse, growing elements)
Transfinite Curve{4, 8, 12}  = nx_coarsen + 1 Using Progression r;

// Vertical uniform: 0 to -Zf
Transfinite Curve{13, 15, 17, 19, 21} = nz_uniform + 1;

// Vertical coarsening: -Zf to -D (fine → coarse, growing elements)
Transfinite Curve{14, 16, 18, 20, 22} = nz_coarsen + 1 Using Progression r;

// ==== Surfaces (8 blocks) ====

// S1: left-far, upper (coarse in x, uniform in z)
Curve Loop(1) = {1, 15, -5, -13};
Plane Surface(1) = {1};
Transfinite Surface{1} = {1, 2, 7, 6};
Recombine Surface{1};

// S2: left-fault, upper (uniform)
Curve Loop(2) = {2, 17, -6, -15};
Plane Surface(2) = {2};
Transfinite Surface{2} = {2, 3, 8, 7};
Recombine Surface{2};

// S3: right-fault, upper (uniform)
Curve Loop(3) = {3, 19, -7, -17};
Plane Surface(3) = {3};
Transfinite Surface{3} = {3, 4, 9, 8};
Recombine Surface{3};

// S4: right-far, upper (coarse in x, uniform in z)
Curve Loop(4) = {4, 21, -8, -19};
Plane Surface(4) = {4};
Transfinite Surface{4} = {4, 5, 10, 9};
Recombine Surface{4};

// S5: left-far, lower (coarse in x and z)
Curve Loop(5) = {5, 16, -9, -14};
Plane Surface(5) = {5};
Transfinite Surface{5} = {6, 7, 12, 11};
Recombine Surface{5};

// S6: left-fault, lower (coarse in z)
Curve Loop(6) = {6, 18, -10, -16};
Plane Surface(6) = {6};
Transfinite Surface{6} = {7, 8, 13, 12};
Recombine Surface{6};

// S7: right-fault, lower (coarse in z)
Curve Loop(7) = {7, 20, -11, -18};
Plane Surface(7) = {7};
Transfinite Surface{7} = {8, 9, 14, 13};
Recombine Surface{7};

// S8: right-far, lower (coarse in x and z)
Curve Loop(8) = {8, 22, -12, -20};
Plane Surface(8) = {8};
Transfinite Surface{8} = {9, 10, 15, 14};
Recombine Surface{8};

// ==== Physical groups ====

// Boundary attributes (must match SEASBoundaryTags in C++ code)
Physical Curve(1) = {13, 14};           // FARFIELD_LEFT  (x = -D)
Physical Curve(2) = {21, 22};           // FARFIELD_RIGHT (x = +D)
Physical Curve(3) = {1, 2, 3, 4};       // FREE_SURFACE   (z = 0)
Physical Curve(4) = {9, 10, 11, 12};    // BOTTOM         (z = -D)

// Fault interior interface at x=0
Physical Curve(5) = {17, 18};           // FAULT (x = 0)

// Two physical surfaces (split at x = 0)
Physical Surface(1) = {3, 4, 7, 8};     // Right half (x > 0)
Physical Surface(2) = {1, 2, 5, 6};     // Left half  (x < 0)

Mesh.MshFileVersion = 2.2;
