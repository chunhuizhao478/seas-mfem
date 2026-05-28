// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_godunov_flux_bimaterial.cpp — Phase 9 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.md (AC #2:
// the 2-block bimaterial unit test).  Ported VERBATIM from the hrs-ref
// BimaterialFlux test (test_phaser_bimaterial_flux.cpp, Phase R.1 of
// PLAN_phase_R_exact_bimaterial_riemann_rev3.md) — only this banner is
// changed.  R.1.T-1 doubles as the plan's "near-equal material reduces
// the bimaterial flux to the scalar Godunov flux to round-off" edge case.
//
// Acceptance criteria covered (one test per criterion):
//   R.1.T-1 — homogeneous-limit byte-exact vs GodunovFlux::Interior
//             (100 random Q triples × 10 random material triples).
//   R.1.T-2 — P-wave analytic VELOCITY transmission
//             v_n*   = Tv · v_n^L,   Tv = 2 Z_p^L / (Z_p^L + Z_p^R).
//             (Plan rev-3 §R.1.T-2 uses the stress formula Tp = 2 Z_R /
//             (Z_L + Z_R), which is wrong for the velocity channel —
//             see REVIEW round-1 R-001 and round-2 R-101.)
//   R.1.T-3 — P-wave analytic STRESS at interface
//             Q*[SXX] = -2 Z_p^L Z_p^R / (Z_p^L + Z_p^R),
//             plus the algebraic identity Tp = 1 + R_p where
//             Tp = 2 Z_p^R / (Z_p^L + Z_p^R) (stress transmission) and
//             R_p = (Z_p^R - Z_p^L) / (Z_p^L + Z_p^R).
//   R.1.T-4 — S-wave analytic VELOCITY transmission
//             v_{t1}* = Tsv · v_{t1}^L,
//             Tsv = 2 Z_s^L / (Z_s^L + Z_s^R).  Same caveat as R.1.T-2.
//   R.1.T-5 — rotation invariance under generic unit normal.
//   R.1.T-6 — role-swap symmetry: swap (self,nbr,nor=+x) ↔
//             (nbr,self,nor=-x) gives equivalent flux up to the
//             outward-normal sign flip.
//   R.1.T-7 — qGodLocal + qGodNeighbor = I to round-off.
//   R.1.T-8 — acoustic abort (mu_self = 0 or mu_nbr = 0 → MFEM_VERIFY
//             abort).  Tested via fork/_exit.

#include "mfem.hpp"

#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/godunov_flux_bimaterial.hpp"
#include "../../dynamic/wave_state.hpp"

#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <random>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e \
   << ", tol " << _t << ")\n"; num_failed++; } \
   else { std::cout << "  PASSED: " << m << "\n"; num_passed++; } \
   } while (0)

namespace
{

bool RunInChild(const std::function<void()>& body)
{
   ::fflush(stdout); ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { return false; }
   if (pid == 0)
   {
      ::freopen("/dev/null", "w", stderr);
      try { body(); } catch (...) { ::_exit(1); }
      ::_exit(0);
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}

// Build a random unit vector from a uniformly sampled point on the
// sphere (Marsaglia 1972 reject method).
void RandomUnitVector(std::mt19937& rng, real_t out[3])
{
   std::uniform_real_distribution<real_t> u(-1.0, 1.0);
   while (true)
   {
      real_t a = u(rng), b = u(rng);
      real_t r2 = a*a + b*b;
      if (r2 >= 1.0 || r2 < 1.0e-8) { continue; }
      real_t s = std::sqrt(1.0 - r2);
      out[0] = 2.0 * a * s;
      out[1] = 2.0 * b * s;
      out[2] = 1.0 - 2.0 * r2;
      return;
   }
}

}  // namespace

// ===========================================================================
// R.1.T-1: homogeneous-limit byte-exact vs GodunovFlux::Interior.
// 100 random (Q_self, Q_nbr, nor) × 10 random (lam, mu, rho).
// ===========================================================================
static void R1_T_1_homogeneous_limit_byte_exact()
{
   std::cout << "\n[R.1.T-1] homogeneous-limit vs GodunovFlux::Interior\n";
   std::mt19937 rng(0x12345678);
   std::uniform_real_distribution<real_t> u_state(-1.0, 1.0);
   std::uniform_real_distribution<real_t> u_lam(1e9, 1e11);
   std::uniform_real_distribution<real_t> u_mu (1e9, 1e11);
   std::uniform_real_distribution<real_t> u_rho(1.5e3, 3.5e3);

   // Tolerance: in principle equal via algebraic identity, but the two
   // paths invert different 9×9 matrices via mfem::DenseMatrixInverse:
   // GodunovFlux inverts its R once at construction (cached
   // Ax_plus_/Ax_minus_), BimaterialFlux inverts matR per call.  The
   // FP rounding sequences differ enough that bit-exact agreement is
   // not achievable; 1e-10 rel / 1e-9 abs is the empirical band that
   // holds across the 1000-case sweep below.  See .hpp file header
   // (REVIEW round-1 R-005, round-2 R-102).
   const real_t tol_rel = 1e-10;
   const real_t tol_abs = 1e-9;

   for (int trip = 0; trip < 10; ++trip)
   {
      const real_t lam = u_lam(rng);
      const real_t mu  = u_mu(rng);
      const real_t rho = u_rho(rng);
      GodunovFlux flux(lam, mu, rho);
      DenseMatrix fluxLocal, fluxNeighbor;
      for (int trial = 0; trial < 100; ++trial)
      {
         real_t nor[3];
         RandomUnitVector(rng, nor);
         BimaterialFlux::BuildPerFaceFluxMatricesGlobal(
            nor, flux, flux, fluxLocal, fluxNeighbor);
         real_t Q_self[NUM_STATE], Q_nbr[NUM_STATE];
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Q_self[c] = u_state(rng);
            Q_nbr[c]  = u_state(rng);
         }
         real_t F_bi[NUM_STATE], F_god[NUM_STATE];
         BimaterialFlux::ApplyPerFaceFlux(fluxLocal, fluxNeighbor,
                                          Q_self, Q_nbr, F_bi);
         flux.Interior(nor, Q_self, Q_nbr, F_god);
         for (int c = 0; c < NUM_STATE; ++c)
         {
            real_t diff = std::abs(F_bi[c] - F_god[c]);
            real_t ref  = std::max(std::abs(F_god[c]),
                                   std::max(std::abs(F_bi[c]), 1.0));
            if (diff > tol_abs + tol_rel * ref)
            {
               std::cerr << "  trip=" << trip << " trial=" << trial
                         << " c=" << c << " F_bi=" << F_bi[c]
                         << " F_god=" << F_god[c] << " diff=" << diff
                         << " lam=" << lam << " mu=" << mu
                         << " rho=" << rho << "\n";
               TEST_ASSERT(false, "homogeneous-limit byte-exact");
               return;
            }
         }
      }
   }
   TEST_ASSERT(true, "R.1.T-1: 1000 (trip × trial) cases agree to "
               "1e-10 rel / 1e-9 abs");
}

// ===========================================================================
// R.1.T-2/T-3: P-wave analytic VELOCITY transmission + interface stress
// at n̂ = +x̂.
//   Q_self pure right-going P-wave with v_n^L = 1, sigma_nn = -Z_p^L.
//   Q_nbr = 0.
//   Expected:
//     Q*[VX]  = Tv · v_n^L = 2 Z_p^L / (Z_p^L + Z_p^R)   (velocity)
//     Q*[SXX] = -Z_p^R · Q*[VX]
//             = -2 Z_p^L Z_p^R / (Z_p^L + Z_p^R)         (stress)
//     R_p = (Z_p^R - Z_p^L) / (Z_p^L + Z_p^R)            (reflection)
//   Note: the *stress* transmission coefficient Tp = 2 Z_p^R /
//   (Z_p^L + Z_p^R) is what governs Q*[SXX]; Q*[VX] uses the
//   *velocity* coefficient Tv = 2 Z_p^L / (Z_p^L + Z_p^R) — the L
//   and R impedances swap places between the two channels.  Extract
//   Q* by:
//     Q*_FL = qGodLocal_FL · Q_self_FL + qGodNeighbor_FL · Q_nbr_FL.
// ===========================================================================
static void R1_T_2_3_p_wave_transmission_and_reflection()
{
   std::cout << "\n[R.1.T-2/3] P-wave analytic T_p and R_p\n";

   // Fix MFEM's flux/eigenvector sign convention: with col 0 having
   // SXX = -lp and VX = +cp, a pure +cp_L eigenmode of amplitude
   // α_L has v_n = α_L · cp_L and sigma_nn = -α_L · lp_L =
   // -α_L · Z_p^L · cp_L.  Setting v_n = 1 gives α_L = 1/cp_L and
   // sigma_nn = -Z_p^L.  Plan §R.1.T-2 specifies exactly this
   // construction.
   struct Case
   {
      real_t lam_L, mu_L, rho_L;
      real_t lam_R, mu_R, rho_R;
   };
   std::vector<Case> cases = {
      // (Z_p^R / Z_p^L) ratios = ~2, ~1/3, ~5 per plan.
      { 32e9, 32e9, 2670.0,
        128e9, 128e9, 2670.0 },                       // Z_R / Z_L ≈ 2
      { 64e9, 64e9, 2670.0,
        16e9, 16e9, 1500.0 },                         // Z_R / Z_L ≈ 1/3 (rough)
      { 16e9, 16e9, 2670.0,
        400e9, 400e9, 3000.0 },                       // Z_R / Z_L ≈ 5
   };
   const real_t tol_rel = 1e-10;

   for (const auto& c : cases)
   {
      const real_t cp_L  = std::sqrt((c.lam_L + 2.0*c.mu_L) / c.rho_L);
      const real_t cp_R  = std::sqrt((c.lam_R + 2.0*c.mu_R) / c.rho_R);
      const real_t Zp_L  = c.rho_L * cp_L;
      const real_t Zp_R  = c.rho_R * cp_R;
      // Tp is the STRESS / pressure transmission coefficient
      // (LeVeque §22.4 / standard acoustics).  Q*[SXX] = -Zp_L · Tp.
      const real_t Tp    = 2.0 * Zp_R / (Zp_L + Zp_R);
      // Tv is the VELOCITY / particle-motion transmission coefficient.
      // The L and R impedances swap places vs Tp.  Q*[VX] = Tv · v_n^L
      // — the physical interface velocity.  Confirmed by limiting
      // cases: rigid wall Zp_R → ∞ gives Tv → 0 (no motion at wall);
      // free surface Zp_R → 0 gives Tv → 2 (velocity doubles).
      const real_t Tv    = 2.0 * Zp_L / (Zp_L + Zp_R);
      const real_t Rp    = (Zp_R - Zp_L) / (Zp_L + Zp_R);

      // Q_L state per plan §R.1.T-2: SXX = -Z_p^L, VX = 1, others
      // zero.  Note: this is NOT a pure +cp_L eigenmode of A^L (those
      // also need SYY = SZZ = -lam_L/cp_L); the plan's setup
      // deliberately tests the case where Q_L has content in the
      // zero-mode columns of matR too, exercising the full bi-material
      // projection.
      real_t Q_L[NUM_STATE] = {0};
      Q_L[SXX] = -Zp_L;
      Q_L[VX]  = 1.0;
      real_t Q_R[NUM_STATE] = {0};

      // Build qGodLocal/qGodNeighbor in face-local frame.  Since n̂ =
      // (1,0,0), face-local frame == global frame; qGodLocal_FL acts
      // directly on Q_L (no rotation needed for the interface state).
      DenseMatrix qGodL_FL, qGodN_FL;
      BimaterialFlux::BuildGodunovStateFaceLocal(
         c.lam_L, c.mu_L, c.rho_L,
         c.lam_R, c.mu_R, c.rho_R,
         qGodL_FL, qGodN_FL);

      // Q*_FL = qGodLocal_FL · Q_L + qGodNeighbor_FL · Q_R
      real_t Q_star[NUM_STATE] = {0};
      for (int i = 0; i < NUM_STATE; ++i)
      {
         real_t s = 0.0;
         for (int j = 0; j < NUM_STATE; ++j)
         {
            s += qGodL_FL(i, j) * Q_L[j] + qGodN_FL(i, j) * Q_R[j];
         }
         Q_star[i] = s;
      }

      // R.1.T-2: physical interface velocity v_n* = Tv · v_n^L
      // (velocity transmission, NOT stress transmission).
      TEST_NEAR(Q_star[VX], Tv * Q_L[VX],
                std::max<real_t>(tol_rel * std::abs(Tv), tol_rel),
                "R.1.T-2 P-wave transmission v_n*");

      // R.1.T-3: reflected P-wave amplitude on L side.
      // The reflected wave is a LEFT-going P (col 8 of matR_L, which
      // uses neighbor material in the bi-material matR; but for the
      // reflected wave on the LEFT side we should reconstruct from
      // physics: total left-side state at the interface = incident
      // + reflected.  The interface velocity v_n^* should equal
      // (1 - R_p) · v_n^L_incident + R_p · v_n^L_incident · (...)
      // Actually the clean check is via Q*[SXX]:
      //   Q*[SXX] = -Zp_L · (1 - R_p) − (compression continuity)
      //           = -Zp_R · v_n^* (impedance ratio gives the
      //                            transmitted side's relation)
      //           = -Zp_R · Tp = -2 Zp_L Zp_R / (Zp_L + Zp_R)
      // Verify Q*[SXX] matches this analytic value.
      const real_t Q_star_SXX_expected = -2.0 * Zp_L * Zp_R / (Zp_L + Zp_R);
      TEST_NEAR(Q_star[SXX], Q_star_SXX_expected,
                std::max<real_t>(tol_rel * std::abs(Q_star_SXX_expected),
                                 tol_rel),
                "R.1.T-3 P-wave reflection via Q*[SXX] continuity");

      // Cross-check: the reflection coefficient relation
      // v_n^* = (1 + R_p) · v_n^L_incident gives Tp = 1 + R_p
      // (since the L-side state at the interface is the sum of
      // incident + reflected velocities).
      TEST_NEAR(Tp, 1.0 + Rp, 1e-12,
                "R.1.T-3 identity T_p = 1 + R_p (impedance match)");
   }
}

// ===========================================================================
// R.1.T-4: S-wave analytic transmission for the v_{t1} channel.
// ===========================================================================
static void R1_T_4_s_wave_transmission()
{
   std::cout << "\n[R.1.T-4] S-wave analytic transmission T_s\n";
   struct Case
   {
      real_t lam_L, mu_L, rho_L;
      real_t lam_R, mu_R, rho_R;
   };
   std::vector<Case> cases = {
      { 32e9, 16e9, 2670.0,
        50e9, 64e9, 2700.0 },           // mu_L != mu_R; lam_L != lam_R
      { 100e9, 50e9, 2900.0,
        20e9, 20e9, 2500.0 },           // lower mu_R
   };
   const real_t tol_rel = 1e-10;

   for (const auto& c : cases)
   {
      const real_t cs_L = std::sqrt(c.mu_L / c.rho_L);
      const real_t cs_R = std::sqrt(c.mu_R / c.rho_R);
      const real_t Zs_L = c.rho_L * cs_L;
      const real_t Zs_R = c.rho_R * cs_R;
      // Velocity (particle-motion) transmission coefficient for the
      // S-wave v_{t1} channel.  See R.1.T-2 for the Tv vs Ts (stress
      // transmission) distinction — Q*[VY] is the interface velocity
      // and follows Tsv = 2 Zs_L / (Zs_L + Zs_R), NOT Ts.
      const real_t Tsv  = 2.0 * Zs_L / (Zs_L + Zs_R);

      // Pure +cs_L S-wave on L side, y-polarised: from matR col 1
      // with α_L = 1/cs_L → v_{t1} = 1, sigma_{nt1} = -mu_L / cs_L =
      // -rho_L · cs_L = -Z_s^L.
      real_t Q_L[NUM_STATE] = {0};
      Q_L[SXY] = -Zs_L;
      Q_L[VY]  = 1.0;
      real_t Q_R[NUM_STATE] = {0};

      DenseMatrix qGodL_FL, qGodN_FL;
      BimaterialFlux::BuildGodunovStateFaceLocal(
         c.lam_L, c.mu_L, c.rho_L,
         c.lam_R, c.mu_R, c.rho_R,
         qGodL_FL, qGodN_FL);

      real_t Q_star[NUM_STATE] = {0};
      for (int i = 0; i < NUM_STATE; ++i)
      {
         real_t s = 0.0;
         for (int j = 0; j < NUM_STATE; ++j)
         {
            s += qGodL_FL(i, j) * Q_L[j] + qGodN_FL(i, j) * Q_R[j];
         }
         Q_star[i] = s;
      }

      // Interface v_{t1} = Tsv · v_{t1}^L (velocity transmission).
      TEST_NEAR(Q_star[VY], Tsv * Q_L[VY],
                std::max<real_t>(tol_rel * std::abs(Tsv), tol_rel),
                "R.1.T-4 S-wave transmission v_{t1}*");
   }
}

// ===========================================================================
// R.1.T-5: rotation invariance.  Rotate (nor, Q_self, Q_nbr) by R_ext
// and check the resulting F differs by the same rotation.
// Implementation: pick nor1 = +x̂ axis-aligned, nor2 = generic unit;
// build a face-local frame change and verify the bi-material F transforms
// covariantly.  Simpler equivalent: compare the per-face flux matrix
// to its global-frame composition for two normals related by a known
// rotation; algebraic identity gives bit-equal up to FP.
// ===========================================================================
static void R1_T_5_rotation_invariance()
{
   std::cout << "\n[R.1.T-5] rotation invariance\n";
   GodunovFlux flux_L(20e9, 18e9, 2600.0);
   GodunovFlux flux_R(50e9, 40e9, 2900.0);

   // Generic Q in some global frame.
   real_t Q_self_g[NUM_STATE] = {  1.0, -2.0,  3.0,  0.5, -0.4,
                                   0.8,  1.5, -0.7,  2.1 };
   real_t Q_nbr_g [NUM_STATE] = { -1.5,  2.4, -0.3,  1.1,  0.6,
                                  -0.9, -2.0,  1.8,  0.4 };

   // Reference: axis-aligned normal, raw Q.
   real_t nor1[3] = {1.0, 0.0, 0.0};
   DenseMatrix fL1, fN1;
   BimaterialFlux::BuildPerFaceFluxMatricesGlobal(
      nor1, flux_L, flux_R, fL1, fN1);
   real_t F1[NUM_STATE];
   BimaterialFlux::ApplyPerFaceFlux(fL1, fN1, Q_self_g, Q_nbr_g, F1);

   // Build the rotation R = T2·Tinv1 that maps the face-local
   // frame at nor1 to the face-local frame at nor2.  Applying R to
   // Q_self / Q_nbr re-expresses them as if the original face-local
   // states were observed at the nor2 viewpoint.  The bi-material
   // flux is covariant under R: F(nor2, R·Q_self, R·Q_nbr) un-
   // rotated by R^{-1} = T1·Tinv2 must equal F(nor1, Q_self, Q_nbr).
   const real_t inv_sqrt3 = 1.0 / std::sqrt(3.0);
   real_t nor2[3] = {inv_sqrt3, inv_sqrt3, inv_sqrt3};
   real_t t1_1[3], t2_1[3], t1_2[3], t2_2[3];
   GodunovFlux::BuildFrame(nor1, t1_1, t2_1);
   GodunovFlux::BuildFrame(nor2, t1_2, t2_2);

   DenseMatrix T1(NUM_STATE, NUM_STATE), Tinv1(NUM_STATE, NUM_STATE);
   DenseMatrix T2(NUM_STATE, NUM_STATE), Tinv2(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotation       (nor1, t1_1, t2_1, T1);
   GodunovFlux::BuildRotationInverse(nor1, t1_1, t2_1, Tinv1);
   GodunovFlux::BuildRotation       (nor2, t1_2, t2_2, T2);
   GodunovFlux::BuildRotationInverse(nor2, t1_2, t2_2, Tinv2);

   DenseMatrix R(NUM_STATE, NUM_STATE);
   mfem::Mult(T2, Tinv1, R);

   real_t Q_self_rot[NUM_STATE], Q_nbr_rot[NUM_STATE];
   R.Mult(Q_self_g, Q_self_rot);
   R.Mult(Q_nbr_g,  Q_nbr_rot);

   DenseMatrix fL2, fN2;
   BimaterialFlux::BuildPerFaceFluxMatricesGlobal(
      nor2, flux_L, flux_R, fL2, fN2);
   real_t F2[NUM_STATE];
   BimaterialFlux::ApplyPerFaceFlux(fL2, fN2, Q_self_rot, Q_nbr_rot, F2);

   DenseMatrix Rinv(NUM_STATE, NUM_STATE);
   mfem::Mult(T1, Tinv2, Rinv);
   real_t F2_unrot[NUM_STATE];
   Rinv.Mult(F2, F2_unrot);

   // Tolerance is relative to the GLOBAL max(|F|), not per-component.
   // A_FL contains entries O(λ, μ, 1/ρ) so |F| can be O(1e10) while
   // some components (e.g. F[SYZ] = 0 exactly here) are ~0.  The
   // rounding error in those zero components is bounded by
   // ε · ‖matrix‖ · ‖Q*_FL‖ ≈ ε · |max(F)|, not ε · |F[c]|.
   real_t F1_max = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      F1_max = std::max(F1_max, std::abs(F1[c]));
   }
   const real_t tol = 1e-12 * std::max<real_t>(F1_max, 1.0);
   bool ok = true;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      const real_t d = std::abs(F1[c] - F2_unrot[c]);
      if (d > tol) { ok = false; }
   }
   TEST_ASSERT(ok, "R.1.T-5 rotation covariance: "
                   "F(nor2, R·Q) un-rotated == F(nor1, Q)");
}

// ===========================================================================
// R.1.T-6: role-swap symmetry.  Compute F with (self=A, nbr=B, nor=+x̂).
// Then compute F' with (self=B, nbr=A, nor=-x̂) and matching swap of Q
// arguments.  F + F' should be related by the bi-material continuity
// (the same physical Q* at the interface, but each side applies its
// own A^self).  Specifically: the interface Q* is symmetric under
// (L,R,+n) ↔ (R,L,-n) — swapping the role of L and R is the same
// physical Riemann problem.
// ===========================================================================
static void R1_T_6_role_swap_symmetry()
{
   std::cout << "\n[R.1.T-6] role-swap symmetry\n";
   real_t lam_A = 30e9, mu_A = 20e9, rho_A = 2600.0;
   real_t lam_B = 60e9, mu_B = 50e9, rho_B = 2900.0;

   real_t Q_A[NUM_STATE] = {  1.0, -2.0,  3.0,  0.5, -0.4,
                              0.8,  1.5, -0.7,  2.1 };
   real_t Q_B[NUM_STATE] = { -1.5,  2.4, -0.3,  1.1,  0.6,
                             -0.9, -2.0,  1.8,  0.4 };

   // mirror_x: under face-local x̂ → -x̂ (i.e. swapping which side
   // is "self"), state-vector components transform as
   //   σ_xx, σ_yy, σ_zz, σ_yz : invariant (no x-tensor index, or
   //                            even number of them)
   //   σ_xy, σ_xz, v_x        : sign-flip (single x-tensor index)
   //   v_y, v_z               : invariant
   const real_t mirror[NUM_STATE] = { +1, +1, +1, -1, +1, -1, -1, +1, +1 };

   // The BA picture's face-local +x̂ is OPPOSITE the AB picture's
   // (each side's face normal points INTO the neighbour).  So the
   // same physical states must be mirrored before being fed to the
   // BA call.  Without this mirroring the test compares apples to
   // oranges and the assertion below cannot hold for generic Q.
   real_t Q_A_mir[NUM_STATE], Q_B_mir[NUM_STATE];
   for (int c = 0; c < NUM_STATE; ++c)
   {
      Q_A_mir[c] = mirror[c] * Q_A[c];
      Q_B_mir[c] = mirror[c] * Q_B[c];
   }

   // Build qGodLocal/Neighbor with the two role assignments in the
   // face-local frame.
   DenseMatrix qGodL_AB, qGodN_AB, qGodL_BA, qGodN_BA;
   BimaterialFlux::BuildGodunovStateFaceLocal(
      lam_A, mu_A, rho_A, lam_B, mu_B, rho_B,
      qGodL_AB, qGodN_AB);
   BimaterialFlux::BuildGodunovStateFaceLocal(
      lam_B, mu_B, rho_B, lam_A, mu_A, rho_A,
      qGodL_BA, qGodN_BA);

   // Physical invariant: at a bi-material interface, the INTERFACE-
   // CONTINUOUS quantities (normal traction σ_nn = σ_xx, shear
   // traction σ_nt1 = σ_xy and σ_nt2 = σ_xz, all three velocity
   // components) must agree between A's POV and B's POV up to the
   // face-local frame mirror.  The TANGENTIAL STRESSES σ_yy, σ_zz,
   // σ_yz are NOT interface-continuous in elastodynamics — each side
   // carries its own material's tangential stress state at the
   // interface — so Q*_AB and Q*_BA legitimately disagree on those
   // three components.  The test therefore checks the 6 continuous
   // components only.
   //   Q*_AB = qGodL_AB · Q_A         + qGodN_AB · Q_B         (A's POV)
   //   Q*_BA = qGodL_BA · mirror(Q_B) + qGodN_BA · mirror(Q_A) (B's POV)
   // For c ∈ {SXX, SXY, SXZ, VX, VY, VZ}:
   //   Q*_BA[c] = mirror[c] · Q*_AB[c]
   real_t Q_star_AB[NUM_STATE], Q_star_BA[NUM_STATE];
   for (int i = 0; i < NUM_STATE; ++i)
   {
      real_t s_ab = 0.0, s_ba = 0.0;
      for (int j = 0; j < NUM_STATE; ++j)
      {
         s_ab += qGodL_AB(i, j) * Q_A[j]     + qGodN_AB(i, j) * Q_B[j];
         s_ba += qGodL_BA(i, j) * Q_B_mir[j] + qGodN_BA(i, j) * Q_A_mir[j];
      }
      Q_star_AB[i] = s_ab;
      Q_star_BA[i] = s_ba;
   }

   // 1 = continuous (assert mirror relation), 0 = tangential stress
   // (legitimately differs between A's and B's Q*).
   const int interface_continuous[NUM_STATE] =
      { 1, 0, 0, 1, 0, 1, 1, 1, 1 };
   //   SXX SYY SZZ SXY SYZ SXZ VX VY VZ

   // Tolerance: floor at max |Q*| component (matrix products amplify
   // roundoff by the largest entry, ~1e7 here from the (lam, mu, rho)
   // scales).
   real_t Q_max = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      Q_max = std::max(Q_max, std::abs(Q_star_AB[c]));
      Q_max = std::max(Q_max, std::abs(Q_star_BA[c]));
   }
   const real_t tol = 1e-12 * std::max<real_t>(Q_max, 1.0);
   bool ok = true;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      if (!interface_continuous[c]) { continue; }
      const real_t expected = mirror[c] * Q_star_AB[c];
      const real_t diff = std::abs(Q_star_BA[c] - expected);
      if (diff > tol) { ok = false; }
   }
   TEST_ASSERT(ok, "R.1.T-6 role-swap symmetry: Q*_BA = mirror_x(Q*_AB) "
                   "on interface-continuous components "
                   "{σ_xx, σ_xy, σ_xz, v_x, v_y, v_z}");
}

// ===========================================================================
// R.1.T-7: qGodLocal + qGodNeighbor = I to round-off.
// ===========================================================================
static void R1_T_7_partition_of_identity()
{
   std::cout << "\n[R.1.T-7] qGodLocal + qGodNeighbor = I invariant\n";
   std::mt19937 rng(0xCAFEBABE);
   std::uniform_real_distribution<real_t> u_lam(1e9, 1e11);
   std::uniform_real_distribution<real_t> u_mu (1e9, 1e11);
   std::uniform_real_distribution<real_t> u_rho(1.5e3, 3.5e3);

   const real_t tol = 1e-10;
   for (int trial = 0; trial < 20; ++trial)
   {
      const real_t lam_L = u_lam(rng), mu_L = u_mu(rng), rho_L = u_rho(rng);
      const real_t lam_R = u_lam(rng), mu_R = u_mu(rng), rho_R = u_rho(rng);
      DenseMatrix qGodL, qGodN;
      BimaterialFlux::BuildGodunovStateFaceLocal(
         lam_L, mu_L, rho_L, lam_R, mu_R, rho_R, qGodL, qGodN);
      bool ok = true;
      for (int i = 0; i < NUM_STATE && ok; ++i)
      {
         for (int j = 0; j < NUM_STATE && ok; ++j)
         {
            real_t sum = qGodL(i, j) + qGodN(i, j);
            real_t target = (i == j ? 1.0 : 0.0);
            if (std::abs(sum - target) > tol) { ok = false; }
         }
      }
      if (!ok) { TEST_ASSERT(false, "qGodLocal + qGodNeighbor = I"); return; }
   }
   TEST_ASSERT(true, "R.1.T-7: 20 random bi-material triples all "
               "satisfy partition-of-identity");
}

// ===========================================================================
// R.1.T-8: acoustic abort fires.
// ===========================================================================
static void R1_T_8_acoustic_abort()
{
   std::cout << "\n[R.1.T-8] acoustic-input abort (mu = 0)\n";
   // Self acoustic.
   const bool aborted_self = RunInChild([](){
      DenseMatrix qL, qN;
      BimaterialFlux::BuildGodunovStateFaceLocal(
         32e9, 0.0, 2670.0,            // mu_self = 0
         32e9, 32e9, 2670.0,
         qL, qN);
   });
   TEST_ASSERT(aborted_self,
               "R.1.T-8 mu_self=0 aborts via MFEM_VERIFY");
   // Nbr acoustic.
   const bool aborted_nbr = RunInChild([](){
      DenseMatrix qL, qN;
      BimaterialFlux::BuildGodunovStateFaceLocal(
         32e9, 32e9, 2670.0,
         32e9, 0.0, 2670.0,            // mu_nbr = 0
         qL, qN);
   });
   TEST_ASSERT(aborted_nbr,
               "R.1.T-8 mu_nbr=0 aborts via MFEM_VERIFY");
}

// ===========================================================================
// R.1.T-extra: ApplyPerFaceFlux MFEM_VERIFY size guard fires on
// wrong-sized inputs (release-build safety; covers REVIEW R-004).
// ===========================================================================
static void R1_T_extra_apply_size_guard()
{
   std::cout << "\n[R.1.T-extra] ApplyPerFaceFlux size-guard\n";
   const bool aborted_local = RunInChild([](){
      DenseMatrix bad(3, 3);                       bad = 0.0;
      DenseMatrix good(NUM_STATE, NUM_STATE);      good = 0.0;
      real_t Q_self[NUM_STATE] = {0}, Q_nbr[NUM_STATE] = {0};
      real_t F[NUM_STATE];
      BimaterialFlux::ApplyPerFaceFlux(bad, good, Q_self, Q_nbr, F);
   });
   TEST_ASSERT(aborted_local,
               "ApplyPerFaceFlux aborts on wrong-size fluxLocal");

   const bool aborted_nbr = RunInChild([](){
      DenseMatrix good(NUM_STATE, NUM_STATE);      good = 0.0;
      DenseMatrix bad(3, 3);                       bad = 0.0;
      real_t Q_self[NUM_STATE] = {0}, Q_nbr[NUM_STATE] = {0};
      real_t F[NUM_STATE];
      BimaterialFlux::ApplyPerFaceFlux(good, bad, Q_self, Q_nbr, F);
   });
   TEST_ASSERT(aborted_nbr,
               "ApplyPerFaceFlux aborts on wrong-size fluxNeighbor");
}

// ===========================================================================
// R.1.T-extra2: BuildPerFaceFluxMatricesGlobal MFEM_VERIFY rejects a
// non-unit normal vector (covers REVIEW round-2 R-106).
// ===========================================================================
static void R1_T_extra_nor_unit_guard()
{
   std::cout << "\n[R.1.T-extra] BuildPerFaceFluxMatricesGlobal nor-unit "
                "guard\n";
   const bool aborted_long = RunInChild([](){
      GodunovFlux flux_L(20e9, 18e9, 2600.0);
      GodunovFlux flux_R(50e9, 40e9, 2900.0);
      real_t nor_nonunit[3] = {2.0, 0.0, 0.0};   // |nor| = 2
      DenseMatrix fL, fN;
      BimaterialFlux::BuildPerFaceFluxMatricesGlobal(
         nor_nonunit, flux_L, flux_R, fL, fN);
   });
   TEST_ASSERT(aborted_long,
               "BuildPerFaceFluxMatricesGlobal aborts on |nor|=2");

   const bool aborted_short = RunInChild([](){
      GodunovFlux flux_L(20e9, 18e9, 2600.0);
      GodunovFlux flux_R(50e9, 40e9, 2900.0);
      real_t nor_nonunit[3] = {0.5, 0.0, 0.0};   // |nor| = 0.5
      DenseMatrix fL, fN;
      BimaterialFlux::BuildPerFaceFluxMatricesGlobal(
         nor_nonunit, flux_L, flux_R, fL, fN);
   });
   TEST_ASSERT(aborted_short,
               "BuildPerFaceFluxMatricesGlobal aborts on |nor|=0.5");
}

// ===========================================================================
int main(int argc, char* argv[])
{
   std::cout << "test_godunov_flux_bimaterial: Phase 9 BimaterialFlux unit tests\n";
   R1_T_1_homogeneous_limit_byte_exact();
   R1_T_2_3_p_wave_transmission_and_reflection();
   R1_T_4_s_wave_transmission();
   R1_T_5_rotation_invariance();
   R1_T_6_role_swap_symmetry();
   R1_T_7_partition_of_identity();
   R1_T_8_acoustic_abort();
   R1_T_extra_apply_size_guard();
   R1_T_extra_nor_unit_guard();

   std::cout << "\n========================================\n"
             << "Results: " << num_passed << "/" << num_tests << " passed";
   if (num_failed > 0)
   {
      std::cout << " (" << num_failed << " FAILED)\n";
      return 1;
   }
   std::cout << "\n========================================\n";
   return 0;
}
