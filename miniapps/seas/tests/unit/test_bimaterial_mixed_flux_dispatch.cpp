// Phase 3 (PLAN_mixed_flux_hetero_riemann.md): per-face mixed-flux DISPATCH on
// the bi-material (matrix) operator.
//
//   Test 3.1 — per-face kernel selection: a central-set face dispatches the
//              CENTRAL flux (single-valued: F_h_e1 == F_h_e2), while a
//              non-central interior face dispatches the bi-material GODUNOV
//              flux (per-side: F_h_e1 != F_h_e2 under heterogeneity).
//   Test 3.2 — BUG-3 LOAD-BEARING guard: on a fault-adjacent central face with
//              a GENUINE material contrast (mu_e1 != mu_e2), InteriorFaceFlux_
//              gives F_h_e1 == F_h_e2 componentwise to <= 1e-12*scale.  Fails on
//              a non-single-valued (per-side-swapped) storage/dispatch; passes
//              with the Phase-2 side-symmetric store + Phase-3 copy.
//   Test 3.3 — mixed_flux=none is byte-exact Godunov: the SAME central-set face
//              dispatches the bi-material Godunov path on a None operator
//              (central branch never taken), matching the per_face_bimaterial_
//              flux_ apply bit-for-bit; enabling Adjacent CHANGES the dispatch.
//
// Fixture: a row of 5 unit hexes along x (elems 0..4).  The fault is the
// interior face at x=3 (between elems 2 and 3).  Adjacent mode therefore puts
// the fault-adjacent non-fault interior faces {(1,2),(3,4)} in the central set;
// the non-fault interior face (0,1) is NOT adjacent to the fault, so it stays a
// Godunov face.  A mu(x) gradient makes every element's material distinct, so
// every interior face here is heterogeneous (the contrast Test 3.2 needs).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/godunov_flux_bimaterial.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do {                                          \
   num_tests++;                                                              \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; }     \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: "       \
                                  << msg << "\n"; }                          \
} while (0)

namespace
{
// Expose the protected InteriorFaceFlux_ for direct per-face testing.
template <typename MeshType>
class TestableBimat : public BimaterialWaveOperator<MeshType>
{
public:
   using BimaterialWaveOperator<MeshType>::BimaterialWaveOperator;  // inherit ctor
   using BimaterialWaveOperator<MeshType>::InteriorFaceFlux_;       // protected -> public
};

// Row of `nx` unit hexes along x; the interior face at x == fault_x is tagged
// FAULT (attr 3), all external faces FREE (attr 1).
Mesh BuildHexRowFaultMesh(int nx, int fault_x)
{
   Mesh mesh(3, (nx + 1) * 4, nx, 0);
   auto vid = [](int ix, int iy, int iz) { return ix * 4 + iy * 2 + iz; };
   for (int ix = 0; ix <= nx; ++ix)
      for (int iy = 0; iy < 2; ++iy)
         for (int iz = 0; iz < 2; ++iz)
         {
            real_t v[3] = { static_cast<real_t>(ix), static_cast<real_t>(iy),
                            static_cast<real_t>(iz) };
            mesh.AddVertex(v);
         }
   for (int e = 0; e < nx; ++e)
   {
      const int v[8] = {
         vid(e, 0, 0), vid(e + 1, 0, 0), vid(e + 1, 1, 0), vid(e, 1, 0),
         vid(e, 0, 1), vid(e + 1, 0, 1), vid(e + 1, 1, 1), vid(e, 1, 1)
      };
      mesh.AddHex(v, 1);
   }
   mesh.FinalizeTopology();

   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cx += mesh.GetVertex(fv[i])[0]; }
      cx /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cx - fault_x) < 1e-9)
         {
            mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3);   // FAULT
         }
         continue;                                            // other interior: untagged
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);          // FREE (external)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

real_t WorstRel(const real_t *a, const real_t *b)
{
   real_t scale = 0.0;
   for (int i = 0; i < NUM_STATE; ++i)
   {
      scale = std::max(scale, std::max(std::abs(a[i]), std::abs(b[i])));
   }
   if (scale == 0.0) { return 0.0; }
   real_t m = 0.0;
   for (int i = 0; i < NUM_STATE; ++i)
   {
      m = std::max(m, std::abs(a[i] - b[i]) / scale);
   }
   return m;
}

// A small fixed set of (Q_self, Q_nbr) pairs exercising stress / velocity /
// mixed channels (deterministic, no RNG portability concerns).
struct StatePair { real_t self[NUM_STATE]; real_t nbr[NUM_STATE]; };

std::vector<StatePair> MakeStatePairs()
{
   std::vector<StatePair> ps(6);
   for (auto &p : ps)
   {
      for (int c = 0; c < NUM_STATE; ++c) { p.self[c] = 0.0; p.nbr[c] = 0.0; }
   }
   // Stress-dominated jump.
   for (int c = 0; c < SXY; ++c) { ps[0].self[c] = 1.0e7; ps[0].nbr[c] = 0.5e7; }
   // Shear-stress jump.
   ps[1].self[SXY] = 2.0e6; ps[1].nbr[SXY] = -1.0e6;
   ps[1].self[SXZ] = 1.0e6; ps[1].nbr[SXZ] =  3.0e6;
   // Velocity-dominated jump.
   for (int c = VX; c < NUM_STATE; ++c) { ps[2].self[c] = 1.0; ps[2].nbr[c] = -0.5; }
   // Mixed.
   ps[3].self[SXX] = 7.0e6; ps[3].self[VX] = 0.8; ps[3].self[SYZ] = -2.0e6;
   ps[3].nbr[SXX]  = -3.0e6; ps[3].nbr[VY]  = 0.4; ps[3].nbr[SXZ]  =  1.5e6;
   // (P3-3) Two more deterministic pairs to widen the per-side Godunov vs
   // single-valued central discriminator beyond the four hand-picked above.
   // Anti-symmetric stress jump across the face (every stress channel flips).
   for (int c = 0; c < SXY; ++c) { ps[4].self[c] =  4.0e6; ps[4].nbr[c] = -4.0e6; }
   ps[4].self[SXY] = -2.5e6; ps[4].nbr[SXY] =  2.5e6;
   ps[4].self[VZ]  =  0.6;   ps[4].nbr[VZ]  =  0.6;   // common-mode velocity
   // Mixed stress + multi-velocity with no symmetry, large dynamic range.
   ps[5].self[SYY] = 9.0e6; ps[5].self[SXZ] =  1.2e6; ps[5].self[VY] = -0.9;
   ps[5].self[VX]  = 0.3;
   ps[5].nbr[SYY]  = 1.0e6; ps[5].nbr[SYZ]  = -2.2e6; ps[5].nbr[VZ] =  0.7;
   ps[5].nbr[VX]   = -0.4;
   return ps;
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== Phase 3: BimaterialWaveOperator mixed-flux dispatch ===\n";

   // 5-hex row; fault interior face at x=3 (between elems 2 and 3).
   Mesh mesh = BuildHexRowFaultMesh(/*nx=*/5, /*fault_x=*/3);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;

   // Heterogeneous (Mode::Coefficient) material: mu(x) and lambda(x) grow with x
   // so every element centroid (x = 0.5, 1.5, ..., 4.5) has a distinct material
   // => every interior face is heterogeneous.  rho constant.
   auto lam_fn = [](const Vector &x) -> real_t { return 30.0e9 * (1.0 + 0.2 * x(0)); };
   auto mu_fn  = [](const Vector &x) -> real_t { return 30.0e9 * (1.0 + 0.2 * x(0)); };
   FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn);
   ConstantCoefficient  rho_c(2670.0);

   // Adjacent operator (mixed flux ON).
   TestableBimat<Mesh> op_adj(
      mesh, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);
   op_adj.SetMixedFluxMode(MixedFluxMode::Adjacent);   // (Phase 3) no longer aborts

   // None operator (mixed flux OFF; ctor default).
   TestableBimat<Mesh> op_none(
      mesh, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);

   const auto &cset = op_adj.GetCentralFluxFaceSet();
   TEST_ASSERT(!cset.empty(),
               "Adjacent mode populated a non-empty central_flux_face_set_ on the "
               "fault-bearing hex row (no silent no-op)");

   // Pick a central face f_c (a 2-sided interior face in the central set) and a
   // non-central interior non-fault face f_n.
   const Array<int> &fault_faces = op_adj.GetFaultInteriorFaces();
   auto is_fault = [&](int f) {
      for (int i = 0; i < fault_faces.Size(); ++i) { if (fault_faces[i] == f) return true; }
      return false;
   };
   int f_c = -1, f_n = -1;
   for (int f : cset)
   {
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0) { f_c = f; break; }
   }
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr || ftr->Elem2No < 0) { continue; }   // boundary / 1-sided
      if (cset.count(f) > 0) { continue; }            // central
      if (is_fault(f)) { continue; }                  // the fault itself
      f_n = f;
      break;
   }
   TEST_ASSERT(f_c >= 0, "found a 2-sided interior central-set face f_c");
   TEST_ASSERT(f_n >= 0 && cset.count(f_n) == 0,
               "found a 2-sided interior NON-central, non-fault face f_n");
   std::cout << "  central face f_c=" << f_c << "  non-central face f_n=" << f_n
             << "  (central set size " << cset.size() << ")\n";

   const real_t dummy_nor[3] = {1.0, 0.0, 0.0};   // ignored by InteriorFaceFlux_
   const std::vector<StatePair> pairs = MakeStatePairs();

   // -----------------------------------------------------------------------
   // Test 3.1 + 3.2: central face is single-valued; non-central face is per-
   // side (Godunov) under heterogeneity.
   // -----------------------------------------------------------------------
   std::cout << "\n-- Test 3.1/3.2: per-face kernel selection + single-valuedness --\n";
   real_t worst_central_diff = 0.0;   // max |F_e1 - F_e2| (central) -> must be ~0
   real_t max_godunov_diff    = 0.0;  // max |F_e1 - F_e2| (Godunov het) -> must be > 0
   for (const auto &p : pairs)
   {
      real_t Fc1[NUM_STATE], Fc2[NUM_STATE];
      op_adj.InteriorFaceFlux_(f_c, p.self, p.nbr, dummy_nor, Fc1, Fc2);
      worst_central_diff = std::max(worst_central_diff, WorstRel(Fc1, Fc2));

      real_t Fn1[NUM_STATE], Fn2[NUM_STATE];
      op_adj.InteriorFaceFlux_(f_n, p.self, p.nbr, dummy_nor, Fn1, Fn2);
      max_godunov_diff = std::max(max_godunov_diff, WorstRel(Fn1, Fn2));
   }
   TEST_ASSERT(worst_central_diff <= 1.0e-12,
               "Test 3.2: central-set face is SINGLE-VALUED under heterogeneity "
               "(F_h_e1 == F_h_e2 to <= 1e-12 rel; BUG-3 load-bearing guard)");
   TEST_ASSERT(max_godunov_diff > 1.0e-3,
               "Test 3.1: non-central interior face dispatches the bi-material "
               "GODUNOV flux (per-side: F_h_e1 != F_h_e2 under heterogeneity), so "
               "the per-face kernel selection is real");

   // Independent value check (Test 3.1): the central F* must equal the Phase-1
   // primitive applied to the SAME stored matrices (sanity that the dispatch
   // routes to per_face_central_flux_, not something else).
   {
      const auto &c = op_adj.GetPerFaceCentralFlux().at(f_c);
      real_t worst = 0.0;
      for (const auto &p : pairs)
      {
         real_t Fdisp1[NUM_STATE], Fdisp2[NUM_STATE], Fprim[NUM_STATE];
         op_adj.InteriorFaceFlux_(f_c, p.self, p.nbr, dummy_nor, Fdisp1, Fdisp2);
         BimaterialFlux::ApplyPerFaceFlux(c[0], c[1], p.self, p.nbr, Fprim);
         worst = std::max(worst, WorstRel(Fdisp1, Fprim));
      }
      TEST_ASSERT(worst == 0.0,
                  "Test 3.1: central dispatch == ApplyPerFaceFlux(per_face_central_"
                  "flux_[f_c], ...) bit-for-bit (routes to the central store)");
   }

   // phaser_dispatch_count_ delta sanity: one InteriorFaceFlux_ on an interior
   // face increments by 2 (both sides), on either kernel.
   {
      op_adj.ResetPhaserDispatchCount();
      real_t a[NUM_STATE], b[NUM_STATE];
      op_adj.InteriorFaceFlux_(f_c, pairs[0].self, pairs[0].nbr, dummy_nor, a, b);
      TEST_ASSERT(op_adj.GetPhaserDispatchCount() == 2,
                  "Test 3.1: central InteriorFaceFlux_ increments phaser_dispatch_"
                  "count_ by 2");
   }

   // -----------------------------------------------------------------------
   // Test 3.3: mixed_flux=none -> the SAME f_c dispatches Godunov (central
   // branch never taken), byte-identical to the per_face_bimaterial_flux_ apply.
   // -----------------------------------------------------------------------
   std::cout << "\n-- Test 3.3: mixed_flux=none is byte-exact Godunov on f_c --\n";
   TEST_ASSERT(op_none.GetCentralFluxFaceSet().empty(),
               "Test 3.3: None operator has an empty central set");
   {
      const auto &bm = op_none.GetPerFaceBimaterialFlux();   // [face][side][0/1]
      real_t worst = 0.0, max_diff_sides = 0.0;
      for (const auto &p : pairs)
      {
         real_t Fe1[NUM_STATE], Fe2[NUM_STATE];
         op_none.InteriorFaceFlux_(f_c, p.self, p.nbr, dummy_nor, Fe1, Fe2);
         // Reference: the Godunov per-side apply from the stored bi-material
         // matrices (this is exactly what the None path runs).
         real_t Re1[NUM_STATE], Re2[NUM_STATE];
         BimaterialFlux::ApplyPerFaceFlux(bm[f_c][0][0], bm[f_c][0][1],
                                          p.self, p.nbr, Re1);
         BimaterialFlux::ApplyPerFaceFlux(bm[f_c][1][0], bm[f_c][1][1],
                                          p.self, p.nbr, Re2);
         worst = std::max(worst, std::max(WorstRel(Fe1, Re1), WorstRel(Fe2, Re2)));
         max_diff_sides = std::max(max_diff_sides, WorstRel(Fe1, Fe2));
      }
      TEST_ASSERT(worst == 0.0,
                  "Test 3.3: None operator's InteriorFaceFlux_(f_c) == the stored "
                  "bi-material Godunov per-side apply, bit-for-bit (central branch "
                  "NOT taken)");
      TEST_ASSERT(max_diff_sides > 1.0e-3,
                  "Test 3.3: on the None path f_c is per-side Godunov (F_e1 != F_e2 "
                  "under het) — enabling Adjacent CHANGES the dispatch to "
                  "single-valued central");
   }

   std::cout << "\n==============================================\n";
   std::cout << "Summary: " << num_passed << " / " << num_tests
             << " passed (" << num_failed << " failed)\n";
   std::cout << "==============================================\n";
   return (num_failed > 0) ? 1 : 0;
}
