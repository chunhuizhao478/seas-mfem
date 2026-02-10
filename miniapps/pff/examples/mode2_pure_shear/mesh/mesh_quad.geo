// Mode 2 Pure Shear Fracture Mesh
// Domain: [0,1] x [-0.5, 0.5] matching pff.cpp
// Physical notch at y=0, from x=0 to x=0.5 (crack tip)
// Refinement in expected crack propagation region

// Use MSH format 2.2 for MFEM compatibility
Mesh.MshFileVersion = 2.2;

// Parameters - mesh sizes matching pff.cpp after local refinement
// Base mesh h = 1/30 = 0.0333, refined h = 0.0333/4 = 0.0083
lc = 0.0333;              // Base mesh size
lc_refined = 0.0083;      // Refined mesh size near crack path
lc_refined_2 = 0.0083;    // Extra refinement near crack tip

// Domain dimensions (matching pff.cpp)
Lx = 1.0;
Ly_half = 0.5;
crack_tip = 0.5;
notch_half_width = 0.001;  // Small notch opening (half-width)

// Points - domain corners and notch
Point(1) = {0, -Ly_half, 0, lc};                      // Bottom-left
Point(2) = {Lx, -Ly_half, 0, lc};                     // Bottom-right
Point(3) = {Lx, Ly_half, 0, lc};                      // Top-right
Point(4) = {0, Ly_half, 0, lc};                       // Top-left
Point(5) = {0, notch_half_width, 0, lc_refined};     // Left edge, above notch
Point(6) = {crack_tip, 0, 0, lc_refined};            // Crack tip
Point(7) = {0, -notch_half_width, 0, lc_refined};    // Left edge, below notch

// Lines
Line(1) = {1, 2};    // Bottom edge
Line(2) = {2, 3};    // Right edge
Line(3) = {3, 4};    // Top edge
Line(4) = {4, 5};    // Left edge (above notch)
Line(5) = {5, 6};    // Upper notch face
Line(6) = {6, 7};    // Lower notch face
Line(7) = {7, 1};    // Left edge (below notch)

// Line loop and surface
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

// Refine mesh along expected crack path
// From pff.cpp: oriented box with center=(0.65, -0.25), length=0.8, width=0.2
// length_direction=(1, -1.5), so crack goes down-right from notch tip
Field[1] = Box;
Field[1].VIn = lc_refined;
Field[1].VOut = lc;
Field[1].XMin = 0.45;      // Start near notch tip
Field[1].XMax = 1.0;       // To right edge
Field[1].YMin = -0.5;      // Down to bottom
Field[1].YMax = 0.1;       // Slightly above y=0
Field[1].Thickness = 0.05; // Transition thickness

// Extra refinement near crack tip
Field[2] = Box;
Field[2].VIn = lc_refined_2;
Field[2].VOut = lc;
Field[2].XMin = 0.45;      // Around crack tip
Field[2].XMax = 0.65;
Field[2].YMin = -0.1;
Field[2].YMax = 0.1;
Field[2].Thickness = 0.02;

// Use the minimum of the two fields
Field[3] = Min;
Field[3].FieldsList = {1, 2};

// Set the background mesh field
Background Field = 3;

Mesh 2; // Generate 2D mesh

// Physical groups for boundary conditions (matching pff.cpp convention)
// Attribute 1 = Bottom, 2 = Right, 3 = Top, 4 = Left
Physical Line("Bottom", 1) = {1};
Physical Line("Right", 2) = {2};
Physical Line("Top", 3) = {3};
Physical Line("Left", 4) = {4, 7};
Physical Line("Notch", 5) = {5, 6};
Physical Surface("Domain", 1) = {1};
