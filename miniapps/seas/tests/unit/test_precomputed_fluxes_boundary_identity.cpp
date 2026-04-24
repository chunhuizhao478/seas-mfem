// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
// Phase 1 Step B acceptance:
// tests/unit/test_precomputed_fluxes_boundary_identity.cpp
//
// Covers the Step B gates (plan §5.2.4):
//   P_BFACE_GODUNOV     : probe matches FreeSurfaceGodunovTotal to 1e-12
//   P_BFACE_GAMMA       : probe matches FreeSurfaceTotal         to 1e-12
//   P_BFACE_ABSORBING   : probe matches AbsorbingTotal           to 1e-12
//   test_R009_probe_linearity  : F_h(2·I) = 2·F_h(I) at bulk_bg=0
//   TestBoundaryProbeIsDeterministic  (formerly misnamed
//       test_R4_003_boundary_probe_uses_topology_normal_only — R-003):
//     repeated calls with identical `n` return bit-equal matrices.  The
//     real tangent-independence invariant is enforced at the SIGNATURE
//     level (no tangent parameters to pass) — there is nothing to perturb
//     at runtime — and the probe-vs-runtime equivalence is covered
//     separately by P_BFACE_GODUNOV / GAMMA / ABSORBING above.  This
//     test remains useful as a determinism smoke check.
//
// P_BFACE_SEISSOL_FREE_SURFACE is documented but not exercised here:
// the SeisSol-convention S = -matR21·matR11⁻¹ hand-port is deferred to
// Phase 1 Step B sign-off per plan §5.2.5.  P_BFACE_GODUNOV is the
// blocking gate; if the probe matches FreeSurfaceGodunovTotal to ULP
// and MFEM's FreeSurfaceGodunovTotal is already validated, the SeisSol
// parity receipt reduces to documenting the A^+/A^- basis convention
// used by MFEM.
//
// P_BFACE_NONZERO_BG_ABORT is also documented but not exercised at
// runtime: MFEM_VERIFY aborts the process in the default build (no
// MFEM_USE_EXCEPTIONS), so it cannot be observed without either a
// subprocess test harness or an exceptions-enabled build.

#include "mfem.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
          std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static void FillSaneRandom(real_t *Q, std::mt19937 &rng)
{
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   for (int c = 0; c < 6; c++) { Q[c] = U(rng) * 2.0e6; }
   for (int c = 6; c < 9; c++) { Q[c] = U(rng) * 2.0; }
}

static Mesh BuildCartesianTetFixture(int nx, int ny, int nz)
{
   Mesh mesh = Mesh::MakeCartesian3D(
      nx, ny, nz, Element::TETRAHEDRON, 1.0, 1.0, 1.0, false);
   mesh.FinalizeTopology();
   mesh.Finalize();
   return mesh;
}

// --------------------------------------------------------------------------
// BuildProbeFixture + Init PrecomputedFaceFluxes with one of three BC modes.
// --------------------------------------------------------------------------
enum class BcMode { Godunov, Gamma, Absorbing };

struct BdrFixture
{
   Mesh mesh;
   DG_FECollection fec;
   FiniteElementSpace fes;
   GodunovFlux flux;
   BoundaryConfig bc;
   std::vector<int> face_bdr_attr;
   PrecomputedFaceFluxes pcf;

   BdrFixture(BcMode mode)
      : mesh(BuildCartesianTetFixture(2, 2, 2)),
        fec(1, 3),
        fes(&mesh, &fec, 1),
        flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho)
   {
      bc.fault_attr = -1;
      for (int a = 1; a <= 6; a++)
      {
         if (mode == BcMode::Absorbing) { bc.absorbing_attrs.insert(a); }
         else                           { bc.natural_attrs.insert(a); }
      }

      face_bdr_attr.assign(mesh.GetNumFaces(), 0);
      for (int b = 0; b < mesh.GetNBE(); b++)
      {
         int f, o; mesh.GetBdrElementFace(b, &f, &o);
         face_bdr_attr[f] = mesh.GetBdrAttribute(b);
      }

      FreeSurfaceBCMode fs_mode = (mode == BcMode::Gamma)
                                  ? FreeSurfaceBCMode::Gamma
                                  : FreeSurfaceBCMode::Godunov;
      std::set<int> fault_set, shared_set;
      pcf.Init(mesh, fes, bc, flux, fs_mode,
               fault_set, face_bdr_attr, shared_set);
   }
};

// --------------------------------------------------------------------------
// P_BFACE_{GODUNOV,GAMMA,ABSORBING}: probe matches the runtime BC flux.
// --------------------------------------------------------------------------
static void TestBoundaryFluxMatchesRuntime(BcMode mode, const char *label)
{
   std::cout << "\n[P_BFACE_" << label << "]\n";
   BdrFixture fx(mode);

   std::mt19937 rng(2026 + (int)mode);
   real_t bg_zero[NUM_STATE] = {0};
   int n_sampled = 0, bad = 0;
   real_t worst = 0.0;

   for (int f = 0; f < fx.mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr =
         fx.mesh.GetFaceElementTransformations(f);
      if (!ftr || ftr->Elem2No >= 0) { continue; }
      int e1 = ftr->Elem1No;
      if (!fx.pcf.HasEntry(f, e1)) { continue; }
      const auto &fe = fx.pcf.GetEntry(f, e1);

      for (int trial = 0; trial < 20; trial++)
      {
         real_t Iself[NUM_STATE];
         FillSaneRandom(Iself, rng);

         real_t F_pc[NUM_STATE] = {0.0};
         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int k = 0; k < NUM_STATE; k++)
            {
               F_pc[c] += fe.nApNm1(c, k) * Iself[k];
            }
         }

         real_t F_runtime[NUM_STATE];
         switch (mode)
         {
            case BcMode::Godunov:
               fx.flux.FreeSurfaceGodunovTotal(fe.normal, Iself, bg_zero,
                                               F_runtime);
               break;
            case BcMode::Gamma:
               fx.flux.FreeSurfaceTotal(fe.normal, Iself, bg_zero,
                                        F_runtime);
               break;
            case BcMode::Absorbing:
               fx.flux.AbsorbingTotal(fe.normal, Iself, bg_zero,
                                      F_runtime);
               break;
         }

         real_t me = 0.0;
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t scale = std::max({real_t(1.0), std::abs(F_pc[c]),
                                     std::abs(F_runtime[c])});
            real_t rel = std::abs(F_pc[c] - F_runtime[c]) / scale;
            me = std::max(me, rel);
         }
         worst = std::max(worst, me);
         if (me > 1e-12) { bad++; }
         n_sampled++;
      }
   }
   std::cout << "  samples=" << n_sampled << " worst_rel=" << worst << "\n";
   TEST_ASSERT(bad == 0 && n_sampled > 0,
               std::string("P_BFACE_") + label + ": probe matches runtime "
               "to 1e-12 on every boundary face QP sample");
}

// --------------------------------------------------------------------------
// test_R009_probe_linearity: F_h(2·I) == 2·F_h(I) for each BC at bulk_bg=0.
// --------------------------------------------------------------------------
static void TestR009ProbeLinearity()
{
   std::cout << "\n[test_R009_probe_linearity]\n";
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu,
                    TPV102Params::rho);
   real_t bg_zero[NUM_STATE] = {0};
   std::mt19937 rng(42);

   const real_t test_normals[3][3] = {
      {1.0, 0.0, 0.0},
      {0.0, -1.0, 0.0},
      {0.0, 0.0, 1.0}
   };

   struct Probe {
      const char *label;
      void (GodunovFlux::*fn)(const real_t*, const real_t*, const real_t*,
                              real_t*) const;
   };
   Probe probes[] = {
      {"FreeSurfaceGodunovTotal", &GodunovFlux::FreeSurfaceGodunovTotal},
      {"FreeSurfaceTotal",        &GodunovFlux::FreeSurfaceTotal},
      {"AbsorbingTotal",          &GodunovFlux::AbsorbingTotal}
   };

   for (const auto &p : probes)
   {
      real_t worst = 0.0;
      int samples = 0, bad = 0;
      for (int n = 0; n < 3; n++)
      {
         for (int t = 0; t < 16; t++)
         {
            real_t I[NUM_STATE], F_s[NUM_STATE], F_d[NUM_STATE];
            FillSaneRandom(I, rng);
            real_t I2[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++) { I2[c] = 2.0 * I[c]; }
            (flux.*(p.fn))(test_normals[n], I,  bg_zero, F_s);
            (flux.*(p.fn))(test_normals[n], I2, bg_zero, F_d);
            real_t me = 0.0;
            for (int c = 0; c < NUM_STATE; c++)
            {
               real_t scale = std::max({real_t(1.0),
                                        std::abs(F_d[c]),
                                        std::abs(F_s[c])});
               real_t rel = std::abs(F_d[c] - 2.0 * F_s[c]) / scale;
               me = std::max(me, rel);
            }
            worst = std::max(worst, me);
            if (me > 1e-14) { bad++; }
            samples++;
         }
      }
      std::cout << "  " << p.label << "  samples=" << samples
                << " worst=" << worst << "\n";
      TEST_ASSERT(bad == 0,
                  std::string(p.label) +
                  " is linear in I_self at bulk_bg=0 (to 1e-14)");
   }
}

// --------------------------------------------------------------------------
// TestBoundaryProbeIsDeterministic  (R-003 rename):
// repeated calls with identical `n` must produce bit-equal matrices.
// The tangent-independence invariant originally documented as R4-003
// is a compile-time property of the signature (no `(t1, t2)` parameters
// to perturb); the runtime equivalence to MFEM's BC flux is covered by
// the P_BFACE_* tests above.  This remains a smoke check for
// determinism / statelessness of the probe.
// --------------------------------------------------------------------------
static void TestBoundaryProbeIsDeterministic()
{
   std::cout << "\n[TestBoundaryProbeIsDeterministic]\n";
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu,
                    TPV102Params::rho);
   real_t n[3] = {0.0, 1.0, 0.0};

   // Signature takes only (n, flux, nApNm1) -- there's literally no way
   // to pass a perturbed (t1, t2).  This is a compile-time receipt.
   // Build the matrix twice and assert equality as a smoke check.
   DenseMatrix M1(NUM_STATE), M2(NUM_STATE);
   PrecomputedFaceFluxes::BuildBoundaryMatricesGodunov(n, flux, M1);
   PrecomputedFaceFluxes::BuildBoundaryMatricesGodunov(n, flux, M2);
   real_t max_err = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      for (int j = 0; j < NUM_STATE; j++)
      {
         max_err = std::max(max_err, std::abs(M1(i,j) - M2(i,j)));
      }
   }
   TEST_ASSERT(max_err == 0.0,
               "BuildBoundaryMatricesGodunov is a pure function of "
               "(n, flux) -- deterministic across repeated calls");

   // Analogous for Gamma and Absorbing.
   DenseMatrix G1(NUM_STATE), G2(NUM_STATE);
   PrecomputedFaceFluxes::BuildBoundaryMatricesGamma(n, flux, G1);
   PrecomputedFaceFluxes::BuildBoundaryMatricesGamma(n, flux, G2);
   real_t gerr = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
      for (int j = 0; j < NUM_STATE; j++)
         gerr = std::max(gerr, std::abs(G1(i,j) - G2(i,j)));
   TEST_ASSERT(gerr == 0.0,
               "BuildBoundaryMatricesGamma deterministic");

   DenseMatrix A1(NUM_STATE), A2(NUM_STATE);
   PrecomputedFaceFluxes::BuildBoundaryMatricesAbsorbing(n, flux, A1);
   PrecomputedFaceFluxes::BuildBoundaryMatricesAbsorbing(n, flux, A2);
   real_t aerr = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
      for (int j = 0; j < NUM_STATE; j++)
         aerr = std::max(aerr, std::abs(A1(i,j) - A2(i,j)));
   TEST_ASSERT(aerr == 0.0,
               "BuildBoundaryMatricesAbsorbing deterministic");
}

int main()
{
   std::cout << "===================================================\n";
   std::cout << "Phase 1 Step B acceptance (plan 2026-04-23 §5.2.4)\n";
   std::cout << "===================================================\n";

   TestBoundaryFluxMatchesRuntime(BcMode::Godunov,   "GODUNOV");
   TestBoundaryFluxMatchesRuntime(BcMode::Gamma,     "GAMMA");
   TestBoundaryFluxMatchesRuntime(BcMode::Absorbing, "ABSORBING");
   TestR009ProbeLinearity();
   TestBoundaryProbeIsDeterministic();

   std::cout << "\n===================================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "===================================================\n";
   return (num_failed > 0) ? 1 : 0;
}
