// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Unit tests for FaultBasis: per-face orthonormal coordinate frames on fault.
// Tests axis-aligned 3D basis, ProjectTraction, EmbedSlip, NormalStress,
// round-trip consistency, 2D degeneracy, and orientation flipping.

#include "test_macros.hpp"
#include "../../fault/fault_basis.hpp"
#include "../../integrator/dg_elasticity_ip_combined_integrator.hpp"

#include <cstdlib>

using namespace mfem;
using namespace mfem::seas;

// =============================================================================
// Helper: find all interior faces of a mesh
// =============================================================================
static Array<int> FindInteriorFaces(Mesh &mesh)
{
   Array<int> faces;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      int e1, e2;
      mesh.GetFaceElements(f, &e1, &e2);
      if (e1 >= 0 && e2 >= 0)
      {
         faces.Append(f);
      }
   }
   return faces;
}

// =============================================================================
// Helper: dot product of 3D arrays
// =============================================================================
static real_t Dot3(const real_t a[3], const real_t b[3], int dim)
{
   real_t sum = 0.0;
   for (int d = 0; d < dim; d++) { sum += a[d] * b[d]; }
   return sum;
}

// =============================================================================
// Helper: norm of 3D array
// =============================================================================
static real_t Norm3(const real_t a[3], int dim)
{
   return std::sqrt(Dot3(a, a, dim));
}

// =============================================================================
// Test 1: Axis-aligned vertical fault in 3D
// =============================================================================
void TestAxisAligned3D()
{
   std::cout << "\n=== Axis-aligned 3D Vertical Fault ===\n";

   // Create 2-element hex mesh: [0,2] x [0,1] x [0,1]
   // Interior face at x=1 with normal along x-axis
   Mesh mesh = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);

   Array<int> fault_faces = FindInteriorFaces(mesh);
   TEST_ASSERT(fault_faces.Size() == 1,
               "2-element 3D mesh has exactly 1 interior face");

   // ref_normal = (1,0,0), up = (0,0,-1) (SCEC: x3 positive downward)
   Vector ref_normal(3);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   TEST_ASSERT(fb.NumFaces() == 1, "NumFaces = 1");
   TEST_ASSERT(fb.Dimension() == 3, "Dimension = 3");
   TEST_ASSERT(fb.NumTangentComponents() == 2, "NumTangentComponents = 2");

   const FaultBasisData &b = fb.GetBasis(0);

   // Expected: n = (1,0,0)
   TEST_NEAR(b.normal[0], 1.0, 1e-12, "n[0] = 1");
   TEST_NEAR(b.normal[1], 0.0, 1e-12, "n[1] = 0");
   TEST_NEAR(b.normal[2], 0.0, 1e-12, "n[2] = 0");

   // Expected: tangent1 = dip = strike x n = (0,0,1)
   // No negation: matches Tandem's Curvilinear.cpp:facetBasis()
   TEST_NEAR(b.tangent1[0], 0.0, 1e-12, "t1[0] = 0 (dip)");
   TEST_NEAR(b.tangent1[1], 0.0, 1e-12, "t1[1] = 0 (dip)");
   TEST_NEAR(b.tangent1[2], 1.0, 1e-12, "t1[2] = 1 (dip)");

   // Expected: tangent2 = strike = (0,-1,0)
   TEST_NEAR(b.tangent2[0], 0.0, 1e-12, "t2[0] = 0 (strike)");
   TEST_NEAR(b.tangent2[1], -1.0, 1e-12, "t2[1] = -1 (strike)");
   TEST_NEAR(b.tangent2[2], 0.0, 1e-12, "t2[2] = 0 (strike)");

   // Orthonormality checks
   real_t n_len = Norm3(b.normal, 3);
   real_t t1_len = Norm3(b.tangent1, 3);
   real_t t2_len = Norm3(b.tangent2, 3);
   TEST_NEAR(n_len, 1.0, 1e-12, "|n| = 1");
   TEST_NEAR(t1_len, 1.0, 1e-12, "|t1| = 1");
   TEST_NEAR(t2_len, 1.0, 1e-12, "|t2| = 1");

   real_t n_t1 = Dot3(b.normal, b.tangent1, 3);
   real_t n_t2 = Dot3(b.normal, b.tangent2, 3);
   real_t t1_t2 = Dot3(b.tangent1, b.tangent2, 3);
   TEST_NEAR(n_t1, 0.0, 1e-12, "n . t1 = 0");
   TEST_NEAR(n_t2, 0.0, 1e-12, "n . t2 = 0");
   TEST_NEAR(t1_t2, 0.0, 1e-12, "t1 . t2 = 0");
}

// =============================================================================
// Test 2: ProjectTraction correctness
// =============================================================================
void TestProjectTraction()
{
   std::cout << "\n=== ProjectTraction ===\n";

   // Same setup as Test 1
   Mesh mesh = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);
   Array<int> fault_faces = FindInteriorFaces(mesh);

   Vector ref_normal(3);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   // With n=(1,0,0), t1=dip=(0,0,1), t2=strike=(0,-1,0):
   // Traction (10, 20, 30) MPa
   real_t traction[3] = {10.0, 20.0, 30.0};
   real_t tau_local[2];
   fb.ProjectTraction(0, traction, tau_local);

   // tau_local[0] = traction . t1 = (10,20,30).(0,0,1) = 30
   TEST_NEAR(tau_local[0], 30.0, 1e-12,
             "ProjectTraction dip = 30 for (10,20,30)");
   // tau_local[1] = traction . t2 = (10,20,30).(0,-1,0) = -20
   TEST_NEAR(tau_local[1], -20.0, 1e-12,
             "ProjectTraction strike = -20 for (10,20,30)");

   // Pure normal traction: (sigma, 0, 0) → tangential should be zero
   real_t normal_traction[3] = {50.0, 0.0, 0.0};
   real_t tau_pure_n[2];
   fb.ProjectTraction(0, normal_traction, tau_pure_n);
   TEST_NEAR(tau_pure_n[0], 0.0, 1e-12,
             "Pure normal traction: dip component = 0");
   TEST_NEAR(tau_pure_n[1], 0.0, 1e-12,
             "Pure normal traction: strike component = 0");
}

// =============================================================================
// Test 3: EmbedSlip correctness
// =============================================================================
void TestEmbedSlip()
{
   std::cout << "\n=== EmbedSlip ===\n";

   Mesh mesh = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);
   Array<int> fault_faces = FindInteriorFaces(mesh);

   Vector ref_normal(3);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   // t1=dip=(0,0,1), t2=strike=(0,-1,0)

   // slip_local = (1, 0) → delta_u should be along t1 = (0,0,1)
   real_t slip_dip[2] = {1.0, 0.0};
   real_t du[3];
   fb.EmbedSlip(0, slip_dip, du);
   TEST_NEAR(du[0], 0.0, 1e-12, "EmbedSlip(1,0): du[0] = 0");
   TEST_NEAR(du[1], 0.0, 1e-12, "EmbedSlip(1,0): du[1] = 0");
   TEST_NEAR(du[2], 1.0, 1e-12, "EmbedSlip(1,0): du[2] = 1 (dip)");

   // slip_local = (0, 1) → delta_u should be along t2 = (0,-1,0)
   real_t slip_strike[2] = {0.0, 1.0};
   fb.EmbedSlip(0, slip_strike, du);
   TEST_NEAR(du[0], 0.0, 1e-12, "EmbedSlip(0,1): du[0] = 0");
   TEST_NEAR(du[1], -1.0, 1e-12, "EmbedSlip(0,1): du[1] = -1 (strike)");
   TEST_NEAR(du[2], 0.0, 1e-12, "EmbedSlip(0,1): du[2] = 0");

   // Mixed slip: (2, 3) → delta_u = 2*(0,0,1) + 3*(0,-1,0) = (0,-3,2)
   real_t slip_mixed[2] = {2.0, 3.0};
   fb.EmbedSlip(0, slip_mixed, du);
   TEST_NEAR(du[0], 0.0, 1e-12, "EmbedSlip(2,3): du[0] = 0");
   TEST_NEAR(du[1], -3.0, 1e-12, "EmbedSlip(2,3): du[1] = -3");
   TEST_NEAR(du[2], 2.0, 1e-12, "EmbedSlip(2,3): du[2] = 2");
}

// =============================================================================
// Test 4: Round-trip consistency
// =============================================================================
void TestRoundTrip()
{
   std::cout << "\n=== Round-Trip Consistency ===\n";

   Mesh mesh = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);
   Array<int> fault_faces = FindInteriorFaces(mesh);

   Vector ref_normal(3);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   // Given a purely tangential traction (no normal component):
   // traction = 5*t1 + 7*t2 = 5*(0,0,1) + 7*(0,-1,0) = (0,-7,5)
   real_t traction[3] = {0.0, -7.0, 5.0};
   real_t tau_local[2];
   fb.ProjectTraction(0, traction, tau_local);

   // Should recover (5, 7)
   TEST_NEAR(tau_local[0], 5.0, 1e-12, "Round-trip: dip component = 5");
   TEST_NEAR(tau_local[1], 7.0, 1e-12, "Round-trip: strike component = 7");

   // Now EmbedSlip with those local components and verify tangential match
   real_t du[3];
   fb.EmbedSlip(0, tau_local, du);
   TEST_NEAR(du[0], traction[0], 1e-12,
             "Round-trip embed: du[0] matches tangential traction[0]");
   TEST_NEAR(du[1], traction[1], 1e-12,
             "Round-trip embed: du[1] matches tangential traction[1]");
   TEST_NEAR(du[2], traction[2], 1e-12,
             "Round-trip embed: du[2] matches tangential traction[2]");

   // With a traction that has a normal component:
   // traction2 = (10, -7, 5) → normal part = 10, tangential part = (0,-7,5)
   real_t traction2[3] = {10.0, -7.0, 5.0};
   fb.ProjectTraction(0, traction2, tau_local);
   TEST_NEAR(tau_local[0], 5.0, 1e-12,
             "Round-trip with normal: dip component = 5");
   TEST_NEAR(tau_local[1], 7.0, 1e-12,
             "Round-trip with normal: strike component = 7");

   // EmbedSlip only captures the tangential part
   fb.EmbedSlip(0, tau_local, du);
   TEST_NEAR(du[0], 0.0, 1e-12,
             "Embed discards normal: du[0] = 0 (not 10)");
}

// =============================================================================
// Test 5: Normal stress extraction
// =============================================================================
void TestNormalStress()
{
   std::cout << "\n=== Normal Stress ===\n";

   Mesh mesh = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);
   Array<int> fault_faces = FindInteriorFaces(mesh);

   Vector ref_normal(3);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   // Traction (10, 0, 0) with n=(1,0,0):
   // sigma_n = -traction . n = -10 (tensile, since convention is positive
   // compression, here traction points away from fault → tension)
   real_t traction[3] = {10.0, 0.0, 0.0};
   real_t sigma_n = fb.NormalStress(0, traction);
   TEST_NEAR(sigma_n, -10.0, 1e-12,
             "NormalStress: (10,0,0) with n=(1,0,0) → -10 (tension)");

   // Compressive traction: (-25e6, 0, 0) → sigma_n = 25e6
   real_t compressive[3] = {-25.0e6, 0.0, 0.0};
   sigma_n = fb.NormalStress(0, compressive);
   TEST_NEAR(sigma_n, 25.0e6, 1e-6,
             "NormalStress: compressive 25 MPa → +25e6");

   // Purely tangential traction: (0, 5, 10) → sigma_n = 0
   real_t tangential[3] = {0.0, 5.0, 10.0};
   sigma_n = fb.NormalStress(0, tangential);
   TEST_NEAR(sigma_n, 0.0, 1e-12,
             "NormalStress: purely tangential → 0");
}

// =============================================================================
// Test 6: 2D degeneracy
// =============================================================================
void Test2DDegeneracy()
{
   std::cout << "\n=== 2D Degeneracy ===\n";

   // Create 2-element quad mesh: [0,2] x [0,1]
   Mesh mesh = Mesh::MakeCartesian2D(
      2, 1, Element::QUADRILATERAL, false, 2.0, 1.0);

   Array<int> fault_faces = FindInteriorFaces(mesh);
   TEST_ASSERT(fault_faces.Size() == 1,
               "2-element 2D mesh has exactly 1 interior face");

   // ref_normal = (1,0), up = (0,-1) (depth positive downward)
   Vector ref_normal(2);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0;
   Vector up(2);
   up(0) = 0.0; up(1) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   TEST_ASSERT(fb.NumFaces() == 1, "NumFaces = 1 (2D)");
   TEST_ASSERT(fb.Dimension() == 2, "Dimension = 2");
   TEST_ASSERT(fb.NumTangentComponents() == 1, "NumTangentComponents = 1");

   const FaultBasisData &b = fb.GetBasis(0);

   // Expected: n = (1,0), tangent1 = (0,1)
   // cross = up[0]*n[1] - up[1]*n[0] = 0*0 - (-1)*1 = 1 > 0 → sign = 1
   // tangent1 = (-sign*n[1], sign*n[0]) = (0, 1)
   TEST_NEAR(b.normal[0], 1.0, 1e-12, "2D n[0] = 1");
   TEST_NEAR(b.normal[1], 0.0, 1e-12, "2D n[1] = 0");
   TEST_NEAR(b.tangent1[0], 0.0, 1e-12, "2D t1[0] = 0");
   TEST_NEAR(b.tangent1[1], 1.0, 1e-12, "2D t1[1] = 1 (depth direction)");

   // Orthonormality
   real_t n_len = Norm3(b.normal, 2);
   real_t t1_len = Norm3(b.tangent1, 2);
   TEST_NEAR(n_len, 1.0, 1e-12, "2D |n| = 1");
   TEST_NEAR(t1_len, 1.0, 1e-12, "2D |t1| = 1");
   real_t n_t1 = Dot3(b.normal, b.tangent1, 2);
   TEST_NEAR(n_t1, 0.0, 1e-12, "2D n . t1 = 0");

   // ProjectTraction: scalar output
   real_t traction[2] = {3.0, 7.0};
   real_t tau_local[1];
   fb.ProjectTraction(0, traction, tau_local);
   // tau = traction . t1 = (3,7).(0,1) = 7
   TEST_NEAR(tau_local[0], 7.0, 1e-12,
             "2D ProjectTraction: scalar traction = 7");

   // EmbedSlip
   real_t slip[1] = {5.0};
   real_t du[2];
   fb.EmbedSlip(0, slip, du);
   // du = 5*(0,1) = (0,5)
   TEST_NEAR(du[0], 0.0, 1e-12, "2D EmbedSlip: du[0] = 0");
   TEST_NEAR(du[1], 5.0, 1e-12, "2D EmbedSlip: du[1] = 5");

   // NormalStress
   real_t sigma_n = fb.NormalStress(0, traction);
   // sigma_n = -(3,7).(1,0) = -3
   TEST_NEAR(sigma_n, -3.0, 1e-12, "2D NormalStress = -3");
}

// =============================================================================
// Test 7: Orientation flip
// =============================================================================
void TestOrientationFlip()
{
   std::cout << "\n=== Orientation Flip ===\n";

   // Create a mesh with a single interior face
   Mesh mesh = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);
   Array<int> fault_faces = FindInteriorFaces(mesh);

   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   // Compute with ref_normal = (1,0,0)
   Vector ref_pos(3);
   ref_pos(0) = 1.0; ref_pos(1) = 0.0; ref_pos(2) = 0.0;
   FaultBasis fb_pos;
   fb_pos.Compute(mesh, fault_faces, ref_pos, up);

   // Compute with ref_normal = (-1,0,0) → should flip the raw normal
   Vector ref_neg(3);
   ref_neg(0) = -1.0; ref_neg(1) = 0.0; ref_neg(2) = 0.0;
   FaultBasis fb_neg;
   fb_neg.Compute(mesh, fault_faces, ref_neg, up);

   const FaultBasisData &bp = fb_pos.GetBasis(0);
   const FaultBasisData &bn = fb_neg.GetBasis(0);

   // Under Tandem's "sign baked in" convention (fault_basis.hpp:442-450),
   // Step-5 negation cancels Step-3 flip, so the STORED basis equals the
   // raw `CalcOrtho` direction regardless of `ref_normal`'s sign.  The
   // observable effect of flipping ref_normal is ONLY in `sign_flipped`
   // (diagnostic) — the stored basis vectors themselves are invariant.
   TEST_NEAR(bp.normal[0], bn.normal[0], 1e-12,
             "Basis invariant under ref_normal sign: normal[0]");
   TEST_NEAR(bp.normal[1], bn.normal[1], 1e-12,
             "Basis invariant under ref_normal sign: normal[1]");
   TEST_NEAR(bp.normal[2], bn.normal[2], 1e-12,
             "Basis invariant under ref_normal sign: normal[2]");
   TEST_NEAR(bp.tangent1[0], bn.tangent1[0], 1e-12,
             "Basis invariant: tangent1[0]");
   TEST_NEAR(bp.tangent2[0], bn.tangent2[0], 1e-12,
             "Basis invariant: tangent2[0]");

   // sign_flipped is the observable of ref_normal sign change: exactly
   // one of the two FaultBasis objects has sign_flipped set.
   TEST_ASSERT(bp.sign_flipped != bn.sign_flipped,
               "sign_flipped toggles when ref_normal sign flips");

   // Both should have unit normals
   real_t len_pos = Norm3(bp.normal, 3);
   real_t len_neg = Norm3(bn.normal, 3);
   TEST_NEAR(len_pos, 1.0, 1e-12, "Flipped pos: |n| = 1");
   TEST_NEAR(len_neg, 1.0, 1e-12, "Flipped neg: |n| = 1");

   // Both should maintain orthonormality
   real_t t1_len = Norm3(bn.tangent1, 3);
   real_t t2_len = Norm3(bn.tangent2, 3);
   TEST_NEAR(t1_len, 1.0, 1e-12, "Flipped: |t1| = 1");
   TEST_NEAR(t2_len, 1.0, 1e-12, "Flipped: |t2| = 1");
   real_t n_t1 = Dot3(bn.normal, bn.tangent1, 3);
   real_t n_t2 = Dot3(bn.normal, bn.tangent2, 3);
   real_t t1_t2 = Dot3(bn.tangent1, bn.tangent2, 3);
   TEST_NEAR(n_t1, 0.0, 1e-12, "Flipped: n . t1 = 0");
   TEST_NEAR(n_t2, 0.0, 1e-12, "Flipped: n . t2 = 0");
   TEST_NEAR(t1_t2, 0.0, 1e-12, "Flipped: t1 . t2 = 0");
}

// =============================================================================
// Test 8: Multi-face mesh (4 elements, 3 interior faces in 3D)
// =============================================================================
void TestMultiFace3D()
{
   std::cout << "\n=== Multi-Face 3D Mesh ===\n";

   // Create 4x1x1 hex mesh: [0,4] x [0,1] x [0,1]
   // This has 3 interior faces (at x=1, x=2, x=3)
   Mesh mesh = Mesh::MakeCartesian3D(
      4, 1, 1, Element::HEXAHEDRON, 4.0, 1.0, 1.0);

   Array<int> fault_faces = FindInteriorFaces(mesh);
   TEST_ASSERT(fault_faces.Size() == 3,
               "4-element 3D mesh has 3 interior faces");

   Vector ref_normal(3);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   TEST_ASSERT(fb.NumFaces() == 3, "NumFaces = 3");

   // All faces should have the same basis (all aligned with x-axis)
   for (int i = 0; i < fb.NumFaces(); i++)
   {
      const FaultBasisData &b = fb.GetBasis(i);

      // Normal should be (1,0,0) for all faces
      TEST_NEAR(b.normal[0], 1.0, 1e-12, "multi-face: n[0] = 1");

      // Orthonormality
      real_t n_len = Norm3(b.normal, 3);
      real_t t1_len = Norm3(b.tangent1, 3);
      real_t t2_len = Norm3(b.tangent2, 3);
      TEST_NEAR(n_len, 1.0, 1e-12, "multi-face: |n| = 1");
      TEST_NEAR(t1_len, 1.0, 1e-12, "multi-face: |t1| = 1");
      TEST_NEAR(t2_len, 1.0, 1e-12, "multi-face: |t2| = 1");

      real_t n_t1 = Dot3(b.normal, b.tangent1, 3);
      real_t n_t2 = Dot3(b.normal, b.tangent2, 3);
      real_t t1_t2 = Dot3(b.tangent1, b.tangent2, 3);
      TEST_NEAR(n_t1, 0.0, 1e-12, "multi-face: n . t1 = 0");
      TEST_NEAR(n_t2, 0.0, 1e-12, "multi-face: n . t2 = 0");
      TEST_NEAR(t1_t2, 0.0, 1e-12, "multi-face: t1 . t2 = 0");
   }

   // Consistency: all faces should produce the same basis
   const FaultBasisData &b0 = fb.GetBasis(0);
   for (int i = 1; i < fb.NumFaces(); i++)
   {
      const FaultBasisData &bi = fb.GetBasis(i);
      for (int d = 0; d < 3; d++)
      {
         TEST_NEAR(bi.normal[d], b0.normal[d], 1e-12,
                   "multi-face: normals consistent across faces");
         TEST_NEAR(bi.tangent1[d], b0.tangent1[d], 1e-12,
                   "multi-face: tangent1 consistent across faces");
         TEST_NEAR(bi.tangent2[d], b0.tangent2[d], 1e-12,
                   "multi-face: tangent2 consistent across faces");
      }
   }

   // ProjectTraction/EmbedSlip should work for all face indices
   real_t traction[3] = {5.0, 10.0, 15.0};
   real_t tau0[2], tau2[2];
   fb.ProjectTraction(0, traction, tau0);
   fb.ProjectTraction(2, traction, tau2);
   TEST_NEAR(tau0[0], tau2[0], 1e-12,
             "multi-face: ProjectTraction consistent (dip)");
   TEST_NEAR(tau0[1], tau2[1], 1e-12,
             "multi-face: ProjectTraction consistent (strike)");
}

// =============================================================================
// Test 9: Empty fault_faces array
// =============================================================================
void TestEmptyFaultFaces()
{
   std::cout << "\n=== Empty Fault Faces ===\n";

   Mesh mesh = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);

   Array<int> empty_faces;  // no fault faces

   Vector ref_normal(3);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;
   Vector up(3);
   up(0) = 0.0; up(1) = 0.0; up(2) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, empty_faces, ref_normal, up);

   TEST_ASSERT(fb.NumFaces() == 0, "Empty: NumFaces = 0");
   TEST_ASSERT(fb.Dimension() == 3, "Empty: Dimension = 3");
   TEST_ASSERT(fb.NumTangentComponents() == 2, "Empty: NumTangentComponents = 2");
}

// =============================================================================
// Test 10: 2D multi-face mesh
// =============================================================================
void TestMultiFace2D()
{
   std::cout << "\n=== Multi-Face 2D Mesh ===\n";

   // Create 4x1 quad mesh: [0,4] x [0,1]
   // This has 3 interior faces (at x=1, x=2, x=3)
   Mesh mesh = Mesh::MakeCartesian2D(
      4, 1, Element::QUADRILATERAL, false, 4.0, 1.0);

   Array<int> fault_faces = FindInteriorFaces(mesh);
   TEST_ASSERT(fault_faces.Size() == 3,
               "4-element 2D mesh has 3 interior faces");

   // ref_normal = (1,0), up = (0,-1)
   Vector ref_normal(2);
   ref_normal(0) = 1.0; ref_normal(1) = 0.0;
   Vector up(2);
   up(0) = 0.0; up(1) = -1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   TEST_ASSERT(fb.NumFaces() == 3, "2D multi-face: NumFaces = 3");
   TEST_ASSERT(fb.Dimension() == 2, "2D multi-face: Dimension = 2");
   TEST_ASSERT(fb.NumTangentComponents() == 1, "2D multi-face: NumTangentComponents = 1");

   // All faces should have the same basis (all aligned with x-axis)
   for (int i = 0; i < fb.NumFaces(); i++)
   {
      const FaultBasisData &b = fb.GetBasis(i);

      // Normal should be (1,0) for all faces
      TEST_NEAR(b.normal[0], 1.0, 1e-12, "2D multi-face: n[0] = 1");
      TEST_NEAR(b.normal[1], 0.0, 1e-12, "2D multi-face: n[1] = 0");

      // Tangent should be (0,1) for all faces
      TEST_NEAR(b.tangent1[0], 0.0, 1e-12, "2D multi-face: t1[0] = 0");
      TEST_NEAR(b.tangent1[1], 1.0, 1e-12, "2D multi-face: t1[1] = 1");

      // Orthonormality
      real_t n_len = Norm3(b.normal, 2);
      real_t t1_len = Norm3(b.tangent1, 2);
      TEST_NEAR(n_len, 1.0, 1e-12, "2D multi-face: |n| = 1");
      TEST_NEAR(t1_len, 1.0, 1e-12, "2D multi-face: |t1| = 1");
      real_t n_t1 = Dot3(b.normal, b.tangent1, 2);
      TEST_NEAR(n_t1, 0.0, 1e-12, "2D multi-face: n . t1 = 0");
   }

   // Consistency across faces
   const FaultBasisData &b0 = fb.GetBasis(0);
   for (int i = 1; i < fb.NumFaces(); i++)
   {
      const FaultBasisData &bi = fb.GetBasis(i);
      for (int d = 0; d < 2; d++)
      {
         TEST_NEAR(bi.normal[d], b0.normal[d], 1e-12,
                   "2D multi-face: normals consistent");
         TEST_NEAR(bi.tangent1[d], b0.tangent1[d], 1e-12,
                   "2D multi-face: tangent1 consistent");
      }
   }

   // ProjectTraction/EmbedSlip for different face indices
   real_t traction[2] = {3.0, 7.0};
   real_t tau0[1], tau2[1];
   fb.ProjectTraction(0, traction, tau0);
   fb.ProjectTraction(2, traction, tau2);
   TEST_NEAR(tau0[0], tau2[0], 1e-12,
             "2D multi-face: ProjectTraction consistent across faces");
}

// =============================================================================
// Test 11: 2D orientation flip
// =============================================================================
void TestOrientationFlip2D()
{
   std::cout << "\n=== 2D Orientation Flip ===\n";

   Mesh mesh = Mesh::MakeCartesian2D(
      2, 1, Element::QUADRILATERAL, false, 2.0, 1.0);
   Array<int> fault_faces = FindInteriorFaces(mesh);

   Vector up(2);
   up(0) = 0.0; up(1) = -1.0;

   // Compute with ref_normal = (1,0)
   Vector ref_pos(2);
   ref_pos(0) = 1.0; ref_pos(1) = 0.0;
   FaultBasis fb_pos;
   fb_pos.Compute(mesh, fault_faces, ref_pos, up);

   // Compute with ref_normal = (-1,0) → should flip the raw normal
   Vector ref_neg(2);
   ref_neg(0) = -1.0; ref_neg(1) = 0.0;
   FaultBasis fb_neg;
   fb_neg.Compute(mesh, fault_faces, ref_neg, up);

   const FaultBasisData &bp = fb_pos.GetBasis(0);
   const FaultBasisData &bn = fb_neg.GetBasis(0);

   // Under Tandem's "sign baked in" convention, Step-5 negation cancels
   // Step-3 flip, so the stored basis equals the raw CalcOrtho direction
   // regardless of ref_normal's sign.  Only `sign_flipped` observes the
   // ref_normal sign toggle.
   TEST_NEAR(bp.normal[0], bn.normal[0], 1e-12,
             "2D basis invariant under ref_normal sign: normal[0]");
   TEST_NEAR(bp.normal[1], bn.normal[1], 1e-12,
             "2D basis invariant under ref_normal sign: normal[1]");
   TEST_NEAR(bp.tangent1[0], bn.tangent1[0], 1e-12,
             "2D basis invariant: tangent1[0]");

   // sign_flipped toggles between the two ref_normal signs.
   TEST_ASSERT(bp.sign_flipped != bn.sign_flipped,
               "2D sign_flipped toggles when ref_normal sign flips");

   // Both should have unit normals
   real_t len_pos = Norm3(bp.normal, 2);
   real_t len_neg = Norm3(bn.normal, 2);
   TEST_NEAR(len_pos, 1.0, 1e-12, "2D flipped pos: |n| = 1");
   TEST_NEAR(len_neg, 1.0, 1e-12, "2D flipped neg: |n| = 1");

   // Both should maintain orthonormality
   real_t t1_len_pos = Norm3(bp.tangent1, 2);
   real_t t1_len_neg = Norm3(bn.tangent1, 2);
   TEST_NEAR(t1_len_pos, 1.0, 1e-12, "2D flipped pos: |t1| = 1");
   TEST_NEAR(t1_len_neg, 1.0, 1e-12, "2D flipped neg: |t1| = 1");

   real_t n_t1_pos = Dot3(bp.normal, bp.tangent1, 2);
   real_t n_t1_neg = Dot3(bn.normal, bn.tangent1, 2);
   TEST_NEAR(n_t1_pos, 0.0, 1e-12, "2D flipped pos: n . t1 = 0");
   TEST_NEAR(n_t1_neg, 0.0, 1e-12, "2D flipped neg: n . t1 = 0");
}

// =============================================================================
// Test 12: Non-unit ref_normal and up vectors
// =============================================================================
void TestNonUnitInputVectors()
{
   std::cout << "\n=== Non-Unit Input Vectors ===\n";

   // 3D test: ref_normal and up are non-unit, should still produce correct basis
   Mesh mesh3 = Mesh::MakeCartesian3D(
      2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);
   Array<int> faces3 = FindInteriorFaces(mesh3);

   // Unit vectors (baseline)
   Vector ref_unit(3);
   ref_unit(0) = 1.0; ref_unit(1) = 0.0; ref_unit(2) = 0.0;
   Vector up_unit(3);
   up_unit(0) = 0.0; up_unit(1) = 0.0; up_unit(2) = -1.0;

   FaultBasis fb_unit;
   fb_unit.Compute(mesh3, faces3, ref_unit, up_unit);

   // Non-unit vectors (scaled by arbitrary factors)
   Vector ref_scaled(3);
   ref_scaled(0) = 7.5; ref_scaled(1) = 0.0; ref_scaled(2) = 0.0;
   Vector up_scaled(3);
   up_scaled(0) = 0.0; up_scaled(1) = 0.0; up_scaled(2) = -42.0;

   FaultBasis fb_scaled;
   fb_scaled.Compute(mesh3, faces3, ref_scaled, up_scaled);

   const FaultBasisData &bu = fb_unit.GetBasis(0);
   const FaultBasisData &bs = fb_scaled.GetBasis(0);

   // Results should be identical
   for (int d = 0; d < 3; d++)
   {
      TEST_NEAR(bs.normal[d], bu.normal[d], 1e-12,
                "3D non-unit: normal matches unit-input result");
      TEST_NEAR(bs.tangent1[d], bu.tangent1[d], 1e-12,
                "3D non-unit: tangent1 matches unit-input result");
      TEST_NEAR(bs.tangent2[d], bu.tangent2[d], 1e-12,
                "3D non-unit: tangent2 matches unit-input result");
   }

   // Orthonormality of non-unit result
   real_t n_len = Norm3(bs.normal, 3);
   real_t t1_len = Norm3(bs.tangent1, 3);
   real_t t2_len = Norm3(bs.tangent2, 3);
   TEST_NEAR(n_len, 1.0, 1e-12, "3D non-unit: |n| = 1");
   TEST_NEAR(t1_len, 1.0, 1e-12, "3D non-unit: |t1| = 1");
   TEST_NEAR(t2_len, 1.0, 1e-12, "3D non-unit: |t2| = 1");

   // 2D test: same idea
   Mesh mesh2 = Mesh::MakeCartesian2D(
      2, 1, Element::QUADRILATERAL, false, 2.0, 1.0);
   Array<int> faces2 = FindInteriorFaces(mesh2);

   Vector ref2_unit(2);
   ref2_unit(0) = 1.0; ref2_unit(1) = 0.0;
   Vector up2_unit(2);
   up2_unit(0) = 0.0; up2_unit(1) = -1.0;

   FaultBasis fb2_unit;
   fb2_unit.Compute(mesh2, faces2, ref2_unit, up2_unit);

   Vector ref2_scaled(2);
   ref2_scaled(0) = 3.14; ref2_scaled(1) = 0.0;
   Vector up2_scaled(2);
   up2_scaled(0) = 0.0; up2_scaled(1) = -100.0;

   FaultBasis fb2_scaled;
   fb2_scaled.Compute(mesh2, faces2, ref2_scaled, up2_scaled);

   const FaultBasisData &bu2 = fb2_unit.GetBasis(0);
   const FaultBasisData &bs2 = fb2_scaled.GetBasis(0);

   for (int d = 0; d < 2; d++)
   {
      TEST_NEAR(bs2.normal[d], bu2.normal[d], 1e-12,
                "2D non-unit: normal matches unit-input result");
      TEST_NEAR(bs2.tangent1[d], bu2.tangent1[d], 1e-12,
                "2D non-unit: tangent1 matches unit-input result");
   }

   real_t n2_len = Norm3(bs2.normal, 2);
   real_t t12_len = Norm3(bs2.tangent1, 2);
   TEST_NEAR(n2_len, 1.0, 1e-12, "2D non-unit: |n| = 1");
   TEST_NEAR(t12_len, 1.0, 1e-12, "2D non-unit: |t1| = 1");
}

// =============================================================================
// Test 13: Near-collinear up and normal (boundary of MFEM_VERIFY check)
// =============================================================================
void TestNearCollinearUpNormal()
{
   std::cout << "\n=== Near-Collinear Up and Normal ===\n";

   // Note: MFEM_VERIFY calls abort() when MFEM_USE_EXCEPTIONS is not defined,
   // so we cannot test the exact collinear case with try/catch. Instead we:
   // 1. Test near-collinear (just above tolerance) → should succeed
   // 2. Verify the result is still orthonormal even with nearly-degenerate input

   // 3D: up nearly parallel to normal but with a tiny y-component
   {
      Mesh mesh = Mesh::MakeCartesian3D(
         2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0);
      Array<int> faces = FindInteriorFaces(mesh);

      Vector ref_normal(3);
      ref_normal(0) = 1.0; ref_normal(1) = 0.0; ref_normal(2) = 0.0;

      // up = (1, 1e-6, 0): nearly collinear with n=(1,0,0) but cross product
      // |up x n| = |(0, 0, -1e-6)| = 1e-6 > 1e-12 → passes check
      Vector up_near(3);
      up_near(0) = 1.0; up_near(1) = 1e-6; up_near(2) = 0.0;

      FaultBasis fb;
      fb.Compute(mesh, faces, ref_normal, up_near);

      const FaultBasisData &b = fb.GetBasis(0);

      // Should still produce orthonormal basis
      real_t n_len = Norm3(b.normal, 3);
      real_t t1_len = Norm3(b.tangent1, 3);
      real_t t2_len = Norm3(b.tangent2, 3);
      TEST_NEAR(n_len, 1.0, 1e-12, "3D near-collinear: |n| = 1");
      TEST_NEAR(t1_len, 1.0, 1e-12, "3D near-collinear: |t1| = 1");
      TEST_NEAR(t2_len, 1.0, 1e-12, "3D near-collinear: |t2| = 1");

      real_t n_t1 = Dot3(b.normal, b.tangent1, 3);
      real_t n_t2 = Dot3(b.normal, b.tangent2, 3);
      real_t t1_t2 = Dot3(b.tangent1, b.tangent2, 3);
      TEST_NEAR(n_t1, 0.0, 1e-10, "3D near-collinear: n . t1 = 0");
      TEST_NEAR(n_t2, 0.0, 1e-10, "3D near-collinear: n . t2 = 0");
      TEST_NEAR(t1_t2, 0.0, 1e-10, "3D near-collinear: t1 . t2 = 0");
   }

   // 2D: up nearly parallel to normal
   {
      Mesh mesh = Mesh::MakeCartesian2D(
         2, 1, Element::QUADRILATERAL, false, 2.0, 1.0);
      Array<int> faces = FindInteriorFaces(mesh);

      Vector ref_normal(2);
      ref_normal(0) = 1.0; ref_normal(1) = 0.0;

      // up = (1, 1e-6): cross = 1*0 - 1e-6*1 = -1e-6, |cross|=1e-6 > 1e-12
      Vector up_near(2);
      up_near(0) = 1.0; up_near(1) = 1e-6;

      FaultBasis fb;
      fb.Compute(mesh, faces, ref_normal, up_near);

      const FaultBasisData &b = fb.GetBasis(0);

      // Should still produce orthonormal basis
      real_t n_len = Norm3(b.normal, 2);
      real_t t1_len = Norm3(b.tangent1, 2);
      TEST_NEAR(n_len, 1.0, 1e-12, "2D near-collinear: |n| = 1");
      TEST_NEAR(t1_len, 1.0, 1e-12, "2D near-collinear: |t1| = 1");

      real_t n_t1 = Dot3(b.normal, b.tangent1, 2);
      TEST_NEAR(n_t1, 0.0, 1e-10, "2D near-collinear: n . t1 = 0");
   }
}

// =============================================================================
// Test unified 3-component projection (normal + dip + strike)
// Verifies that ProjectTractionToFaultDOFs with ncomp=3 gives the same
// tangential result as ncomp=2, and that the normal component satisfies
// T_n = T . n_hat (with sign_flipped handled correctly).
// =============================================================================
static void TestUnifiedNormalProjection()
{
   std::cout << "  TestUnifiedNormalProjection...\n";

   // Construct a known orthonormal basis on a planar fault at y=0
   // normal = (0, -1, 0), tangent1(dip) = (0, 0, 1), tangent2(strike) = (-1, 0, 0)
   real_t normal[3]   = {0.0, -1.0, 0.0};
   real_t tangent1[3] = {0.0,  0.0, 1.0};
   real_t tangent2[3] = {-1.0, 0.0, 0.0};

   // Quadrature: 1 point (centroid) on a reference triangle
   IntegrationRule ir(1);
   ir.IntPoint(0).x = 1.0/3.0;
   ir.IntPoint(0).y = 1.0/3.0;
   ir.IntPoint(0).weight = 0.5;  // reference triangle area

   int nbf = 1;  // 1 basis function (p=0, constant)
   int dim = 3;
   int nq = 1;

   // Basis function = 1 at the single quad point
   DenseMatrix e_q(nbf, nq);
   e_q(0, 0) = 1.0;

   // Normal length (face Jacobian determinant)
   Vector nl_q(nq);
   nl_q(0) = 1.0;

   // Known 3D traction at the quad point:
   // T = (T_x, T_y, T_z) = (3.0, -5.0, 7.0)
   // Expected projections:
   //   T . normal   = (3)(0) + (-5)(-1) + (7)(0)  = 5.0
   //   T . tangent1 = (3)(0) + (-5)(0)  + (7)(1)  = 7.0
   //   T . tangent2 = (3)(-1)+ (-5)(0)  + (7)(0)  = -3.0
   Vector traction_q(dim * nq);
   traction_q(0) = 3.0;   // T_x
   traction_q(1) = -5.0;  // T_y
   traction_q(2) = 7.0;   // T_z

   // --- Test 1: 2-component projection (tangential only, sign_flipped=false) ---
   real_t tang_vecs[2][3] = {
      {tangent1[0], tangent1[1], tangent1[2]},
      {tangent2[0], tangent2[1], tangent2[2]}
   };
   Vector trac_2comp;
   DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
      dim, 2, traction_q, nl_q, ir, nbf, e_q,
      tang_vecs, trac_2comp);

   TEST_NEAR(trac_2comp(0), 7.0, 1e-12, "2-comp dip (sf=false)");
   TEST_NEAR(trac_2comp(1), -3.0, 1e-12, "2-comp strike (sf=false)");

   // --- Test 2: 3-component projection (normal + dip + strike, sf=false) ---
   real_t full_basis[3][3] = {
      {normal[0],   normal[1],   normal[2]},
      {tangent1[0], tangent1[1], tangent1[2]},
      {tangent2[0], tangent2[1], tangent2[2]}
   };
   Vector trac_3comp;
   DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
      dim, 3, traction_q, nl_q, ir, nbf, e_q,
      full_basis, trac_3comp);

   TEST_NEAR(trac_3comp(0), 5.0, 1e-12, "3-comp normal (sf=false)");
   TEST_NEAR(trac_3comp(1), 7.0, 1e-12, "3-comp dip (sf=false)");
   TEST_NEAR(trac_3comp(2), -3.0, 1e-12, "3-comp strike (sf=false)");

   // Tangential components must match between 2-comp and 3-comp
   TEST_NEAR(trac_3comp(1), trac_2comp(0), 1e-14, "dip: 3-comp == 2-comp");
   TEST_NEAR(trac_3comp(2), trac_2comp(1), 1e-14, "strike: 3-comp == 2-comp");

   // --- Test 3: sign_flipped = true ---
   // Tandem convention: sign is baked INTO the basis vectors by the caller.
   // To simulate sf=true, pass a pre-negated basis.  The current
   // ProjectTractionToFaultDOFs API has no standalone sign_flipped flag;
   // the only sign-flip channel is via qp_data[q].sign_flipped (tested in
   // Test 5 and TestUnifiedNormalProjectionQPData below).
   real_t full_basis_sf[3][3] = {
      {-normal[0],   -normal[1],   -normal[2]},
      {-tangent1[0], -tangent1[1], -tangent1[2]},
      {-tangent2[0], -tangent2[1], -tangent2[2]}
   };
   Vector trac_3comp_sf;
   DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
      dim, 3, traction_q, nl_q, ir, nbf, e_q,
      full_basis_sf, trac_3comp_sf);

   TEST_NEAR(trac_3comp_sf(0), -5.0, 1e-12, "3-comp normal (sf=true)");
   TEST_NEAR(trac_3comp_sf(1), -7.0, 1e-12, "3-comp dip (sf=true)");
   TEST_NEAR(trac_3comp_sf(2), 3.0, 1e-12, "3-comp strike (sf=true)");

   // Sign-flipped should negate all components uniformly
   TEST_NEAR(trac_3comp_sf(0), -trac_3comp(0), 1e-14, "sf negates normal");
   TEST_NEAR(trac_3comp_sf(1), -trac_3comp(1), 1e-14, "sf negates dip");
   TEST_NEAR(trac_3comp_sf(2), -trac_3comp(2), 1e-14, "sf negates strike");

   // --- Test 4: Verify normal_traction sign convention ---
   // In the elasticity operator: normal_traction = -projected_normal
   // (positive in compression). With sf=false:
   //   projected_normal = T . n = 5.0
   //   normal_traction = -5.0 (tensile, since T_y = -5 pushes inward on y=-1 face)
   // This matches: sigma_n_eff = sn0 + normal_traction = 25 + (-5) = 20 MPa
   real_t normal_traction_val = -trac_3comp(0);  // same as elasticity_operator line 4466
   TEST_NEAR(normal_traction_val, -5.0, 1e-12,
             "normal_traction = -T.n (positive in compression)");

   // --- Test 5: per-QP data path (sign baked in, sign_flipped=true) ---
   // This exercises the qp_data branch in ProjectTractionToFaultDOFs,
   // which is the actual BP5 DG code path.  Per Tandem convention, the
   // caller bakes the sign into the basis vectors — `sign_flipped` in
   // `FaultBasisQPData` is a diagnostic field only.  So to test the
   // "sf=true negates output" semantic we negate the basis vectors
   // explicitly here (which is exactly what `FaultBasis::Compute`
   // does at line 442-450 of fault/fault_basis.hpp when the raw
   // normal opposes ref_normal).
   {
      FaultBasisQPData qpd;
      for (int d = 0; d < 3; d++)
      {
         qpd.normal[d]   = -normal[d];
         qpd.tangent1[d] = -tangent1[d];
         qpd.tangent2[d] = -tangent2[d];
      }
      qpd.nl = 1.0;
      qpd.sign_flipped = true;
      std::vector<FaultBasisQPData> qp_vec = {qpd};

      // 3-comp with qp_data, sign_flipped=true in the per-QP struct.
      // qp_data takes precedence over the centroid `tangents` array.
      Vector trac_3comp_qp;
      DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
         dim, 3, traction_q, nl_q, ir, nbf, e_q,
         full_basis, trac_3comp_qp, &qp_vec);

      // Per-QP sign_flipped=true should negate all components
      TEST_NEAR(trac_3comp_qp(0), -5.0, 1e-12, "qp_data 3-comp normal (sf=true)");
      TEST_NEAR(trac_3comp_qp(1), -7.0, 1e-12, "qp_data 3-comp dip (sf=true)");
      TEST_NEAR(trac_3comp_qp(2),  3.0, 1e-12, "qp_data 3-comp strike (sf=true)");

      // Must match the centroid sf=true result exactly
      TEST_NEAR(trac_3comp_qp(0), trac_3comp_sf(0), 1e-14, "qp vs centroid: normal");
      TEST_NEAR(trac_3comp_qp(1), trac_3comp_sf(1), 1e-14, "qp vs centroid: dip");
      TEST_NEAR(trac_3comp_qp(2), trac_3comp_sf(2), 1e-14, "qp vs centroid: strike");

      // 2-comp with qp_data for backward compat (sign_flipped=false case).
      // Rebuild qpd_nosf with the ORIGINAL (non-negated) basis vectors —
      // the `qpd` in this scope is negated (Test 5's sf=true setup), so
      // copying from it would bake in the flip here too.
      Vector trac_2comp_qp;
      FaultBasisQPData qpd_nosf;
      for (int d = 0; d < 3; d++)
      {
         qpd_nosf.normal[d]   = normal[d];
         qpd_nosf.tangent1[d] = tangent1[d];
         qpd_nosf.tangent2[d] = tangent2[d];
      }
      qpd_nosf.nl = 1.0;
      qpd_nosf.sign_flipped = false;
      std::vector<FaultBasisQPData> qp_vec_nosf = {qpd_nosf};
      DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
         dim, 2, traction_q, nl_q, ir, nbf, e_q,
         tang_vecs, trac_2comp_qp, &qp_vec_nosf);

      TEST_NEAR(trac_2comp_qp(0), 7.0, 1e-12, "qp_data 2-comp dip (sf=false)");
      TEST_NEAR(trac_2comp_qp(1), -3.0, 1e-12, "qp_data 2-comp strike (sf=false)");
   }

   std::cout << "  TestUnifiedNormalProjection... DONE\n";
}

// =============================================================================
// Test qp_data branch of unified 3-component projection
// Verifies that per-quad-point basis data overrides the centroid basis arrays,
// and that sign_flipped from qp_data negates all projected components uniformly.
// =============================================================================
static void TestUnifiedNormalProjectionQPData()
{
   std::cout << "  TestUnifiedNormalProjectionQPData...\n";

   IntegrationRule ir(1);
   ir.IntPoint(0).x = 1.0/3.0;
   ir.IntPoint(0).y = 1.0/3.0;
   ir.IntPoint(0).weight = 0.5;

   const int nbf = 1;
   const int dim = 3;
   const int nq = 1;

   DenseMatrix e_q(nbf, nq);
   e_q(0, 0) = 1.0;

   Vector nl_q(nq);
   nl_q(0) = 1.0;

   // Known 3D traction at the single quad point.
   Vector traction_q(dim * nq);
   traction_q(0) = 3.0;   // T_x
   traction_q(1) = -5.0;  // T_y
   traction_q(2) = 7.0;   // T_z

   // Deliberately wrong centroid basis. If qp_data is used correctly,
   // these values should be ignored.
   real_t wrong_basis[3][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}
   };

   std::vector<FaultBasisQPData> qp_data(1);
   auto &qd = qp_data[0];
   qd.normal[0] = 0.0;  qd.normal[1] = -1.0; qd.normal[2] = 0.0;
   qd.tangent1[0] = 0.0; qd.tangent1[1] = 0.0; qd.tangent1[2] = 1.0;
   qd.tangent2[0] = -1.0; qd.tangent2[1] = 0.0; qd.tangent2[2] = 0.0;
   qd.nl = 1.0;

   // sf=false: should use qp_data basis, not wrong_basis.
   qd.sign_flipped = false;
   Vector trac_qpd;
   DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
      dim, 3, traction_q, nl_q, ir, nbf, e_q,
      wrong_basis, trac_qpd, &qp_data);

   TEST_NEAR(trac_qpd(0), 5.0, 1e-12, "qp_data normal (sf=false)");
   TEST_NEAR(trac_qpd(1), 7.0, 1e-12, "qp_data dip (sf=false)");
   TEST_NEAR(trac_qpd(2), -3.0, 1e-12, "qp_data strike (sf=false)");

   // sf=true in qp_data: per Tandem convention, the caller bakes the sign
   // into the basis vectors; qp_data[q].sign_flipped is diagnostic only.
   // Simulate sf=true by negating qp_data's basis vectors directly.
   qd.sign_flipped = true;
   for (int d = 0; d < 3; d++) {
      qd.normal[d]   = -qd.normal[d];
      qd.tangent1[d] = -qd.tangent1[d];
      qd.tangent2[d] = -qd.tangent2[d];
   }
   Vector trac_qpd_sf;
   DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
      dim, 3, traction_q, nl_q, ir, nbf, e_q,
      wrong_basis, trac_qpd_sf, &qp_data);

   TEST_NEAR(trac_qpd_sf(0), -5.0, 1e-12, "qp_data normal (sf=true)");
   TEST_NEAR(trac_qpd_sf(1), -7.0, 1e-12, "qp_data dip (sf=true)");
   TEST_NEAR(trac_qpd_sf(2), 3.0, 1e-12, "qp_data strike (sf=true)");

   TEST_NEAR(trac_qpd_sf(0), -trac_qpd(0), 1e-14, "qp_data sf negates normal");
   TEST_NEAR(trac_qpd_sf(1), -trac_qpd(1), 1e-14, "qp_data sf negates dip");
   TEST_NEAR(trac_qpd_sf(2), -trac_qpd(2), 1e-14, "qp_data sf negates strike");

   std::cout << "  TestUnifiedNormalProjectionQPData... DONE\n";
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  FaultBasis Unit Tests\n";
   std::cout << "========================================\n";

   TestAxisAligned3D();
   TestProjectTraction();
   TestEmbedSlip();
   TestRoundTrip();
   TestNormalStress();
   Test2DDegeneracy();
   TestOrientationFlip();
   TestMultiFace3D();
   TestEmptyFaultFaces();
   TestMultiFace2D();
   TestOrientationFlip2D();
   TestNonUnitInputVectors();
   TestNearCollinearUpNormal();
   TestUnifiedNormalProjection();
   TestUnifiedNormalProjectionQPData();

   TEST_PRINT_RESULTS();
   return (num_failed > 0) ? 1 : 0;
}
