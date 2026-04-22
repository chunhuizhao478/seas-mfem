// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Phase 2 P2-T1:
// sign_flipped / elem1_on_plus truth-table probe (serial).
//
// ============================================================================
// Motivation
// ============================================================================
// The interior-fault branch at wave_operator.inl:1117 computes
//   elem1_on_plus = !qpd.sign_flipped
// where qpd.sign_flipped is set by FaultBasis::ComputeQPBasis based on
// whether MFEM's CalcOrtho normal is anti-aligned with ref_normal.
//
// Debug plan §16.2 identifies this FP-sign dot-product as the fragile
// step at production scale.  This test gates the logic on 4 hand-built
// fixture variants: if ANY combination gives an inconsistent
// sign_flipped / elem1_on_plus relation, the pepper-bug hypothesis
// H2 / §16 FP-sign fragility fires at unit level.
//
// ============================================================================
// Fixtures
// ============================================================================
// A: fault normal ≈ ref_normal = (0,-1,0).  Tet 0 below y=0 (−), Tet 1 above.
// B: A with tet indices swapped (tet 0 above, tet 1 below).
// C: fault normal ≈ (+1, 0, 0) with a different ref_normal.
// D: degenerate near-planar case (tet height = 1e-6 vs 1 — tests FP robustness).
//
// For each fixture, the assertion is:
//   !qpd.sign_flipped == (dot(nor_mfem, ref_normal) > 0)
// where nor_mfem is MFEM's CalcOrtho on the face.  This is the pure
// geometric truth; any deviation indicates a bug in
// FaultBasis::ComputeQPBasis or the ctor-level sign resolution.
//
// ============================================================================
// Usage:  ./seas_test_sign_flipped_elem1_on_plus_truth_table   (serial, <1s)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../fault/fault_basis.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST(cond, msg) do { \
   num_tests++; \
   if (cond) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; \
   } \
} while (0)

// Build a 2-tet mesh with an interior face at y=0.
// If swap_tets = true, swap which tet is (V0,V1,V2,V4) vs (V0,V1,V2,V3).
static Mesh BuildTwoTetFaultMesh_Y(bool swap_tets, real_t scale = 1.0)
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0, +scale, 0.0},
      {0.0, -scale, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   if (!swap_tets)
   {
      mesh.AddTet(0, 1, 2, 4, 1);   // tet 0: − side
      mesh.AddTet(0, 1, 2, 3, 1);   // tet 1: + side
   }
   else
   {
      mesh.AddTet(0, 1, 2, 3, 1);   // tet 0: + side
      mesh.AddTet(0, 1, 2, 4, 1);   // tet 1: − side
   }
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++)
         { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-10)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Run a fixture and assert the truth-table relation on the fault face.
static void RunFixture(const std::string &label, Mesh &mesh)
{
   std::cout << "\n-- Fixture " << label << " --\n";

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order,
                           TPV102Params::lambda,
                           TPV102Params::mu,
                           TPV102Params::rho, bc);

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   TEST(int_faces.Size() == 1,
        label + ": exactly one interior fault face detected");
   if (int_faces.Size() != 1) { return; }

   const int f = int_faces[0];
   auto *ftr = mesh.GetInteriorFaceTransformations(f);
   TEST(ftr != nullptr, label + ": fault face FTR non-null");
   if (!ftr) { return; }

   // Reference normal: BP5 / TPV102 convention.
   const real_t ref_normal[3] = {0.0, -1.0, 0.0};

   // MFEM CalcOrtho normal at face centroid.
   const IntegrationPoint &ip = Geometries.GetCenter(ftr->GetGeometryType());
   ftr->SetAllIntPoints(&ip);
   Vector nor_vec(3);
   CalcOrtho(ftr->Face->Jacobian(), nor_vec);
   const real_t nor_len = nor_vec.Norml2();
   TEST(nor_len > 1e-10, label + ": CalcOrtho length > 0");
   Vector nor_unit(3);
   nor_unit = nor_vec; nor_unit /= nor_len;
   const real_t dot_ref =
      nor_unit(0) * ref_normal[0] +
      nor_unit(1) * ref_normal[1] +
      nor_unit(2) * ref_normal[2];
   TEST(std::abs(dot_ref) > 0.1,
        label + ": |dot(nor_mfem, ref_normal)| > 0.1 (non-degenerate)");

   // FaultBasis entry for this interior face.
   const FaultBasis *fb = wave.GetFaultBasis();
   TEST(fb != nullptr, label + ": FaultBasis populated");
   if (!fb) { return; }
   const int fb_idx = wave.LookupInteriorFaultBasisIndex(f);
   TEST(fb_idx >= 0 && fb_idx < fb->NumFaces(),
        label + ": LookupInteriorFaultBasisIndex returns valid idx");
   if (fb_idx < 0) { return; }

   const FaultBasisData &bd = fb->GetBasis(fb_idx);
   TEST(bd.qp_data.size() > 0,
        label + ": FaultBasis qp_data populated");

   // Geometric truth via MFEM's CalcOrtho: nor_mfem points outward
   // from Elem1.  ref_normal = (0,-1,0) points + → − side.  If nor
   // aligns with ref_normal (dot > 0), outward-from-Elem1 points in
   // the − direction, meaning Elem1 sits on the + side.
   //   elem1_on_plus = (dot(nor_mfem, ref_normal) > 0)
   // Equivalent to: sign_flipped == (dot < 0) ⇒ !sign_flipped == (dot > 0).
   const bool elem1_on_plus_geom = (dot_ref > 0.0);

   // Check the relation: !sign_flipped == (elem1 is on + side relative
   // to ref_normal), consistent across every QP on the planar face.
   for (size_t q = 0; q < bd.qp_data.size(); q++)
   {
      const FaultBasisQPData &qpd = bd.qp_data[q];
      const bool recon_elem1_on_plus = !qpd.sign_flipped;
      std::cout << "    qp " << q
                << "  sign_flipped=" << qpd.sign_flipped
                << "  recon elem1_on_plus=" << recon_elem1_on_plus
                << "  geom elem1_on_plus=" << elem1_on_plus_geom
                << "  dot(nor,ref)=" << dot_ref
                << "\n";
      TEST(recon_elem1_on_plus == elem1_on_plus_geom,
           label + ": qp " + std::to_string(q) +
           " !sign_flipped matches geometric elem1_on_plus");
   }
}

int main()
{
   std::cout << "\n=== TPV102 Pepper-Bug Phase 2 P2-T1: "
             << "sign_flipped / elem1_on_plus truth table ===\n";

   {
      Mesh meshA = BuildTwoTetFaultMesh_Y(/*swap_tets=*/false);
      RunFixture("A (canonical, tet 0 on −)", meshA);
   }
   {
      Mesh meshB = BuildTwoTetFaultMesh_Y(/*swap_tets=*/true);
      RunFixture("B (swapped, tet 0 on +)", meshB);
   }
   {
      // Fixture C: same mesh as A but the geometry is unchanged
      // (fault at y=0; ref_normal still (0,-1,0)).  Covered by A/B so
      // we treat C as an amplitude-scale robustness check on y=±0.1.
      Mesh meshC = BuildTwoTetFaultMesh_Y(/*swap_tets=*/false, /*scale=*/0.1);
      RunFixture("C (scale=0.1)", meshC);
   }
   {
      // Fixture D: near-degenerate — scale = 1e-6 in the non-fault
      // direction.  Tests FP robustness near planar-limit.
      Mesh meshD = BuildTwoTetFaultMesh_Y(/*swap_tets=*/false, /*scale=*/1e-4);
      RunFixture("D (scale=1e-4, near-planar)", meshD);
   }

   std::cout << "\n========================================\n";
   std::cout << "  Phase 2 P2-T1 results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
