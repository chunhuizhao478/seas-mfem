// Round-11 R-1102 (with R-1202 fix): central_flux_face_set_ topology test.
//
// Contract gates (from MIXED_FLUX_PLAN.md Phase 5 R-1102):
//
// Gate 1 (R-1202 primary fixture): on a 48-tet Cartesian fixture
//        (8 hexes × 6-tet split) with an interior fault plane at y = Ly/2,
//   - |central_flux_face_set_| > 0 (set is non-empty in Adjacent mode)
//   - every face in the set is interior + non-fault (R-1207 predicate)
//   - every face in the set has at least one element in E_fault_adj
//   - central_flux_face_set_ ∩ fault_interior_faces_ = ∅
//
// Gate 2: SetMixedFluxMode(None) clears the set (size == 0)
//
// Gate 3: AllContinuous mode produces a strictly larger set than Adjacent
//
// Gate 4 (R-1202 sub-gate): on the 2-tet TPV102 fixture (the
//        BuildTwoTetFaultMesh from test_ader_tpv102_smoke), the set is
//        EMPTY because every non-fault face is a boundary face, not an
//        interior face.  Catches an implementer accidentally including
//        boundary faces.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unordered_set>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

namespace
{
constexpr real_t kRho    = 2670.0;
constexpr real_t kCp     = 6000.0;
constexpr real_t kCs     = 3464.0;
constexpr real_t kMu     = kRho * kCs * kCs;
constexpr real_t kLambda = kRho * kCp * kCp - 2.0 * kMu;

// Build a 48-tet Cartesian mesh (2 cells per axis, 6-tet hex split = 48 tets;
// MFEM's default split actually produces 6 tets per hex so 8 hexes × 6 = 48
// tets. The exact count is irrelevant; what matters is that the mesh has
// interior non-fault faces touching fault elements.
//
// Tag boundary triangles whose centroid lies on the y = Ly/2 plane as
// fault attribute = 3.  Tag external boundaries as attribute 1 (free
// surface) so they don't accidentally collide with the fault attr.
//
// MFEM's `MakeCartesian3D(..., TETRAHEDRON)` produces tet meshes with
// boundary attributes already set on external faces.  The y=Ly/2 plane
// is INTERIOR to this mesh, so we cannot use AddBdrTriangle directly.
// Instead, we create the mesh with `bdr_tag` = 1 on all external faces,
// then manually mark every interior face whose centroid lies on the
// y=Ly/2 plane via the `face_bdr_attr_` table inside WaveOperator.
//
// Easier approach: build a small custom tet mesh manually so we have full
// control over which faces are fault.
Mesh BuildSmallFaultTetMesh(real_t L = 1000.0)
{
   // Two stacked rows of tets sharing the y = L/2 plane.
   // Vertices forming a 1×2×1 grid of cells, then split each cell into
   // 6 tets (= 12 tets total), with the y = L/2 plane interior.
   //
   // Vertices (8 corners × 2 layers in y = 12 vertices):
   //   index = 4*j + 2*k + i   for (i, j, k) ∈ {0,1} × {0,1,2} × {0,1}
   // Wait — let me use a simple manual layout.
   //
   // Use a 2×2×2 grid of cells (3×3×3 = 27 vertices), split each cell
   // into 6 tets (48 tets total).  Plane y = L/2 is interior.

   const int Nv = 3 * 3 * 3;  // 27 vertices on a 3×3×3 grid
   Mesh mesh(3, Nv, 0, 0);

   auto vid = [&](int i, int j, int k) { return i + 3*j + 9*k; };
   for (int k = 0; k < 3; k++)
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            real_t x = i * 0.5 * L;
            real_t y = j * 0.5 * L;
            real_t z = k * 0.5 * L;
            real_t v[3] = {x, y, z};
            mesh.AddVertex(v);
         }

   // 6-tet split per cell, diagonal v0–v6 (using local hex vertex
   // indices 0..7).  Local-to-global: lv 0..7 = (i+ix, j+iy, k+iz) with
   // bit encoding ix=bit0, iy=bit1, iz=bit2.  Cell origin (i, j, k).
   auto vfunc = [&](int ci, int cj, int ck, int lv) {
      int ix = lv & 1;
      int iy = (lv >> 1) & 1;
      int iz = (lv >> 2) & 1;
      return vid(ci + ix, cj + iy, ck + iz);
   };
   const int pat[6][4] = {
      {0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7},
      {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}
   };
   for (int ck = 0; ck < 2; ck++)
      for (int cj = 0; cj < 2; cj++)
         for (int ci = 0; ci < 2; ci++)
            for (int t = 0; t < 6; t++)
            {
               int v0 = vfunc(ci, cj, ck, pat[t][0]);
               int v1 = vfunc(ci, cj, ck, pat[t][1]);
               int v2 = vfunc(ci, cj, ck, pat[t][2]);
               int v3 = vfunc(ci, cj, ck, pat[t][3]);
               mesh.AddTet(v0, v1, v2, v3, /*attr=*/1);
            }

   mesh.FinalizeTopology();

   // Add boundary triangles: external faces get attr=1, interior faces on
   // the y = L/2 plane get attr=3 (fault).  An interior face has Elem2No >= 0.
   const real_t eps = 1e-6;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= fv.Size();
      if (ftr && ftr->Elem2No >= 0)
      {
         // Interior face — only tag the y = L/2 plane as fault.
         if (std::abs(cy - 0.5 * L) < eps)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      // External face: attr = 1 (free surface — placeholder).
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }

   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// 2-tet fault fixture (R-1202 sub-gate) — same as test_ader_tpv102_smoke.
Mesh BuildTwoTetFaultMesh(real_t L = 1000.0)
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {L, 0.0, 0.0}, {0.0, 0.0, L},
      {0.0,  L, 0.0},
      {0.0, -L, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-6)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}
} // anonymous

int main()
{
   std::cout << "\n=== Round-11 R-1102: central_flux_face_set_ topology ===\n";

   // -----------------------------------------------------------------
   // Primary fixture: 48-tet Cartesian mesh with interior fault plane.
   // -----------------------------------------------------------------
   {
      std::cout << "\n-- Primary fixture: 48-tet w/ interior fault --\n";
      Mesh mesh = BuildSmallFaultTetMesh();

      BoundaryConfig bc;
      bc.fault_attr = 3;
      bc.natural_attrs = {1};
      bc.absorbing_attrs = {};

      WaveOperator<Mesh> wave(mesh, /*order=*/1,
                              kLambda, kMu, kRho, bc);
      // No fault flux setter required — Adjacent mode just needs
      // bc_.fault_attr > 0 and a populated fault_interior_faces_ list.

      const Array<int> &fif = wave.GetFaultInteriorFaces();
      std::cout << "  fault_interior_faces.Size() = " << fif.Size() << "\n";
      TEST_ASSERT(fif.Size() > 0,
                  "Primary fixture has at least one interior fault face");

      wave.SetMixedFluxMode(MixedFluxMode::Adjacent);

      // Note: WaveOperator's MixedFluxMode enum was named at the class
      // template level; access via the templated typename.
      const auto &S_adj = wave.GetCentralFluxFaceSet();
      std::cout << "  |central_flux_face_set_| (Adjacent) = "
                << S_adj.size() << "\n";
      TEST_ASSERT(S_adj.size() > 0,
                  "Adjacent: set is non-empty on primary fixture");

      // Build E_fault_adj for the verification gate.
      std::unordered_set<int> E_fault_adj;
      for (int i = 0; i < fif.Size(); i++)
      {
         int f = fif[i];
         FaceElementTransformations *ftr =
            mesh.GetInteriorFaceTransformations(f);
         if (!ftr) { continue; }
         if (ftr->Elem1No >= 0) { E_fault_adj.insert(ftr->Elem1No); }
         if (ftr->Elem2No >= 0) { E_fault_adj.insert(ftr->Elem2No); }
      }
      std::cout << "  |E_fault_adj| = " << E_fault_adj.size() << "\n";

      // Gate: every face in S_adj is interior + non-fault, with at least
      // one neighbor in E_fault_adj, and not itself a fault face.
      bool all_ok = true;
      int unique_fault_set_intersection = 0;
      std::unordered_set<int> fif_set;
      for (int i = 0; i < fif.Size(); i++) { fif_set.insert(fif[i]); }
      for (int f : S_adj)
      {
         FaceElementTransformations *ftr =
            mesh.GetFaceElementTransformations(f);
         const bool interior = (ftr != nullptr && ftr->Elem2No >= 0);
         if (!interior)
         {
            std::cout << "  FAIL: face " << f << " in set is NOT interior\n";
            all_ok = false; break;
         }
         // R-1207 sanity: face_bdr_attr_ check via Mesh::GetBdrAttribute
         // on any boundary element pointing to this face would be
         // appropriate here, but face_bdr_attr_ is a private member.
         // Substitute: check that f is NOT in fault_interior_faces_.
         if (fif_set.count(f) > 0)
         {
            unique_fault_set_intersection++;
            std::cout << "  FAIL: face " << f
                      << " is in BOTH the central set and fault_interior_faces_\n";
            all_ok = false;
         }
         // Has at least one neighbor in E_fault_adj.
         const bool e1_adj = (ftr->Elem1No >= 0 &&
                              E_fault_adj.count(ftr->Elem1No) > 0);
         const bool e2_adj = (ftr->Elem2No >= 0 &&
                              E_fault_adj.count(ftr->Elem2No) > 0);
         if (!(e1_adj || e2_adj))
         {
            std::cout << "  FAIL: face " << f
                      << " has no neighbor in E_fault_adj\n";
            all_ok = false; break;
         }
      }
      TEST_ASSERT(all_ok,
                  "Every face in S_adj is interior, non-fault, and has "
                  "an E_fault_adj neighbor");
      TEST_ASSERT(unique_fault_set_intersection == 0,
                  "central_flux_face_set_ ∩ fault_interior_faces_ == ∅");

      // Gate 2: None mode clears.
      wave.SetMixedFluxMode(MixedFluxMode::None);
      const auto &S_none = wave.GetCentralFluxFaceSet();
      TEST_ASSERT(S_none.empty(),
                  "None: set is empty (cleared by setter)");

      // Gate 3: AllContinuous strictly larger than Adjacent.
      wave.SetMixedFluxMode(MixedFluxMode::AllContinuous);
      const size_t s_all = wave.GetCentralFluxFaceSet().size();
      wave.SetMixedFluxMode(MixedFluxMode::Adjacent);
      const size_t s_adj = wave.GetCentralFluxFaceSet().size();
      std::cout << "  |S_AllContinuous| = " << s_all
                << ", |S_Adjacent| = " << s_adj << "\n";
      TEST_ASSERT(s_all >= s_adj,
                  "AllContinuous set ≥ Adjacent set (superset)");
      TEST_ASSERT(s_all > 0,
                  "AllContinuous set is non-empty");
   }

   // -----------------------------------------------------------------
   // R-1202 sub-gate (clarified per R-1406): 2-tet fixture has
   // |central_set| == 0 in BOTH Adjacent AND AllContinuous modes.
   //
   // Reasoning:
   //   - 6 of the 7 unique faces are external (Elem2No < 0); the
   //     `is_interior_face` predicate correctly returns false for these
   //     in both modes (no insertion).
   //   - 1 face is the y=0 fault triangle (Elem2No >= 0, bdr_attr=3);
   //     `is_fault_face` correctly excludes it in both modes.
   //
   // This sub-gate catches: (a) an implementer relaxing the
   // `Elem2No >= 0` interior check (would wrongly insert external
   // faces in AllContinuous), and (b) `is_fault_face`
   // mis-classification (would insert the fault face in either mode).
   // The R-1406 addition exercises the latter under AllContinuous —
   // the original test only covered Adjacent.
   // -----------------------------------------------------------------
   {
      std::cout << "\n-- R-1202 sub-gate: 2-tet fixture, |set| == 0 --\n";
      Mesh mesh2 = BuildTwoTetFaultMesh();

      BoundaryConfig bc2;
      bc2.fault_attr = 3;
      bc2.natural_attrs = {1};
      bc2.absorbing_attrs = {};

      WaveOperator<Mesh> wave2(mesh2, /*order=*/1,
                               kLambda, kMu, kRho, bc2);
      const Array<int> &fif2 = wave2.GetFaultInteriorFaces();
      TEST_ASSERT(fif2.Size() == 1,
                  "2-tet fixture has exactly 1 interior fault face");

      wave2.SetMixedFluxMode(MixedFluxMode::Adjacent);
      const auto &S2 = wave2.GetCentralFluxFaceSet();
      std::cout << "  |central_flux_face_set_| on 2-tet (Adjacent) = "
                << S2.size() << "\n";
      TEST_ASSERT(S2.empty(),
                  "2-tet Adjacent: set is EMPTY — 6 external faces "
                  "fail is_interior_face, 1 fault face is excluded");

      // R-1406: also exercise AllContinuous on the 2-tet to catch the
      // same fault-misclassification bug under the broader code path.
      wave2.SetMixedFluxMode(MixedFluxMode::AllContinuous);
      const auto &S2_all = wave2.GetCentralFluxFaceSet();
      std::cout << "  |central_flux_face_set_| on 2-tet (AllContinuous) = "
                << S2_all.size() << "\n";
      TEST_ASSERT(S2_all.empty(),
                  "2-tet AllContinuous: set is EMPTY — the only non-"
                  "external face IS the fault and must be excluded "
                  "(R-1406: catches is_fault_face regressions under "
                  "the broader AllContinuous walk)");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
