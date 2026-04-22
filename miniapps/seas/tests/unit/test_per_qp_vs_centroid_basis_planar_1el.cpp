// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.4 (regression gate): per-QP fault basis
// coincides with the centroid basis on a hand-constructed planar
// fixture.
//
// Fixture: a 2-tetrahedron mesh whose shared interior face is **exactly
// planar** with outward normal = (0, -1, 0) (the BP5/TPV102 canonical
// ref_normal).  The mesh Jacobian on a planar face is constant, so
// CalcOrtho returns the same normal at every QP, and therefore:
//
//   qp_data[q].normal   == centroid.normal   (bit-exact)
//   qp_data[q].tangent1 == centroid.tangent1 (bit-exact)
//   qp_data[q].tangent2 == centroid.tangent2 (bit-exact)
//
// This test is the MFEM-Jacobian regression gate: if FaultBasis or
// MFEM's CalcOrtho ever introduces a subtle QP-dependence on a planar
// face (e.g. an integration-rule-normalized Jacobian, or a reference-
// space orientation flip), this test fails first.
//
// The fixture is small enough (2 tets, planar face) that the
// assertions can be **bit-exact**.  No ULP tolerance.

#include "mfem.hpp"
#include "../../fault/fault_basis.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_EQ_EXACT(val, exp, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp); \
   if (v_ == e_) { num_passed++; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " \
          << msg << " (got " << std::setprecision(17) << v_ \
          << ", expected " << e_ << ", |diff|=" << std::abs(v_ - e_) \
          << ")\n"; } \
} while (0)

// --------------------------------------------------------------------------
// Build a 2-tet mesh whose shared face lies in the plane y = 0.
//
// Vertices:
//   0: (0, -1, 0)   -- below the fault (y < 0)
//   1: (0,  0, 0)   -- fault vertex
//   2: (1,  0, 0)   -- fault vertex
//   3: (0,  0, 1)   -- fault vertex  (TRIANGULAR fault face = {1,2,3})
//   4: (0,  1, 0)   -- above the fault (y > 0)
//
// Tet 0: verts (0, 1, 2, 3)   -- y <= 0 side (elem1)
// Tet 1: verts (4, 1, 2, 3)   -- y >= 0 side (elem2)
//
// Shared face: triangle with vertices (1, 2, 3) -- in the x-z plane.
// --------------------------------------------------------------------------
static Mesh BuildTwoTetPlanarFaultMesh(int fault_attr)
{
   Mesh mesh(3, 5, 2, 1, 3);  // dim=3, 5 verts, 2 elems, 1 bdr elem, sdim=3

   // Vertices
   real_t coords[5][3] = {
      { 0.0, -1.0, 0.0 },
      { 0.0,  0.0, 0.0 },
      { 1.0,  0.0, 0.0 },
      { 0.0,  0.0, 1.0 },
      { 0.0,  1.0, 0.0 }
   };
   for (int i = 0; i < 5; i++) { mesh.AddVertex(coords[i]); }

   // Tetrahedra (attribute 1 for both elements)
   int tet0[4] = {0, 1, 2, 3};
   int tet1[4] = {4, 1, 2, 3};   // Elem2 -- orientation chosen so that the
                                 // (1,2,3) face is shared with tet0.
   mesh.AddTet(tet0, 1);
   mesh.AddTet(tet1, 1);

   // Boundary element on the fault face (attribute = fault_attr).
   int fault_face[3] = {1, 2, 3};
   mesh.AddBdrTriangle(fault_face, fault_attr);

   mesh.FinalizeTopology();
   mesh.Finalize();
   return mesh;
}

int main(int argc, char **argv)
{
   std::cout << "\n=== TPV102 v9.2.0 §4.4 Regression: "
             << "Per-QP vs Centroid Basis on Planar 1-Tet Fixture ===\n";

   const int fault_attr = 3;
   Mesh mesh = BuildTwoTetPlanarFaultMesh(fault_attr);

   // Locate the fault face — the single boundary element carrying
   // fault_attr (mirrors the pattern in wave_operator.inl).
   Array<int> fault_faces;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) != fault_attr) { continue; }
      int face_idx = mesh.GetBdrElementFaceIndex(b);
      auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
      if (ftr == nullptr) { continue; }
      fault_faces.Append(face_idx);
   }

   std::cout << "  fault faces : " << fault_faces.Size() << "\n";
   num_tests++;
   if (fault_faces.Size() != 1)
   {
      num_failed++;
      std::cout << "  FAILED: expected exactly 1 fault face, got "
                << fault_faces.Size() << "\n";
      std::cout << "========================================\n";
      std::cout << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n";
      std::cout << "========================================\n";
      return 1;
   }
   num_passed++;

   // Build FaultBasis with BP5/TPV102 conventions.
   Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
   Vector up(3);         up = 0.0;         up(2) = 1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   // Face geometry is TRIANGLE; use an order=3 rule (several QPs so the
   // "per-QP == centroid" claim is exercised at multiple points).
   Geometry::Type face_geom = Geometry::TRIANGLE;
   {
      auto *ftr0 = mesh.GetInteriorFaceTransformations(fault_faces[0]);
      if (ftr0) { face_geom = ftr0->GetGeometryType(); }
   }
   const int ir_order = 3;
   const IntegrationRule &face_ir = IntRules.Get(face_geom, ir_order);
   fb.ComputeQPBasis(mesh, fault_faces, ref_normal, up, face_ir);

   const FaultBasisData &bd = fb.GetBasis(0);
   std::cout << "  QPs/face    : " << face_ir.GetNPoints() << "\n";
   std::cout << "  centroid n  = (" << bd.normal[0]   << ", "
             << bd.normal[1]   << ", " << bd.normal[2]   << ")\n";
   std::cout << "  centroid t1 = (" << bd.tangent1[0] << ", "
             << bd.tangent1[1] << ", " << bd.tangent1[2] << ")\n";
   std::cout << "  centroid t2 = (" << bd.tangent2[0] << ", "
             << bd.tangent2[1] << ", " << bd.tangent2[2] << ")\n";

   // --------- Verify centroid basis equals expected BP5 canonical -------
   // Expected (ref_normal=(0,-1,0), up=(0,0,1)):
   //   strike = normalize(up x n) = normalize((0,0,1) x (0,-1,0))
   //          = normalize((0*0 - 1*(-1), 1*0 - 0*0, 0*(-1) - 0*0))
   //          = (1, 0, 0)
   //   dip    = strike x n = (1,0,0) x (0,-1,0) = (0, 0, -1)
   //   tangent1 = dip = (0, 0, -1), tangent2 = strike = (1, 0, 0)
   TEST_EQ_EXACT(bd.normal[0],    0.0, "centroid n[0]");
   TEST_EQ_EXACT(bd.normal[1],   -1.0, "centroid n[1]");
   TEST_EQ_EXACT(bd.normal[2],    0.0, "centroid n[2]");

   TEST_EQ_EXACT(bd.tangent1[0],  0.0, "centroid t1[0] (dip)");
   TEST_EQ_EXACT(bd.tangent1[1],  0.0, "centroid t1[1] (dip)");
   TEST_EQ_EXACT(bd.tangent1[2], -1.0, "centroid t1[2] (dip)");

   TEST_EQ_EXACT(bd.tangent2[0],  1.0, "centroid t2[0] (strike)");
   TEST_EQ_EXACT(bd.tangent2[1],  0.0, "centroid t2[1] (strike)");
   TEST_EQ_EXACT(bd.tangent2[2],  0.0, "centroid t2[2] (strike)");

   // --------- Verify per-QP basis == centroid basis (bit-exact) ---------
   int qp_mismatches = 0;
   for (int q = 0; q < face_ir.GetNPoints(); q++)
   {
      const FaultBasisQPData &qpd = bd.qp_data[q];
      for (int d = 0; d < 3; d++)
      {
         if (qpd.normal[d]   != bd.normal[d])   { qp_mismatches++; }
         if (qpd.tangent1[d] != bd.tangent1[d]) { qp_mismatches++; }
         if (qpd.tangent2[d] != bd.tangent2[d]) { qp_mismatches++; }
      }
   }
   num_tests++;
   if (qp_mismatches == 0)
   {
      num_passed++;
      std::cout << "  PASSED: per-QP basis == centroid basis (bit-exact) "
                << "on all " << face_ir.GetNPoints() << " QPs\n";
   }
   else
   {
      num_failed++;
      std::cout << "  FAILED: " << qp_mismatches
                << " QP basis components differ from centroid basis "
                << "on planar fixture (expected bit-exact)\n";
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
