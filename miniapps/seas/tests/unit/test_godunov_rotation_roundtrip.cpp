// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug deep dive: GODUNOV-FLUX ROTATION ROUND-TRIP test.
//
// =============================================================================
// Motivation
// =============================================================================
// The pepper-bug investigation in test_adjacent_triangle_fault_uniformity
// localized a 1e11x per-step amplification of the dip-channel SXZ
// component in bulk Q.  The wave operator's fault-face dispatch composes:
//
//   Tinv_can = BuildRotationInverse(can_n, can_t1, can_t2)  // FaultBasis frame
//   T_can    = BuildRotation       (can_n, can_t1, can_t2)
//   Q_self_can = Tinv_can * Q_self_global
//   ... friction solve in canonical fault-local frame ...
//   Q_imp_global = T_can * Q_imp_can
//   flux_.Interior(can_n, Q_imp_global, ..., F_h)            // BuildFrame frame
//
// flux_.Interior internally calls BuildFrame(can_n) → t1', t2' (a DIFFERENT
// tangent pair than FaultBasis), then computes:
//
//   Tinv' = BuildRotationInverse(can_n, t1', t2')
//   T'    = BuildRotation       (can_n, t1', t2')
//   Q_local = Tinv' * Q_imp_global
//   F_local = A_x * Q_local
//   F_h     = T'   * F_local
//
// For this to give the physically correct F_h = A_n · Q_imp_global
// REGARDLESS of which valid (t1, t2) orthonormal tangent pair was used,
// the rotation operators must satisfy:
//
//   (R1) T(Tinv(σ)) ≡ σ                  for every Voigt-symmetric σ
//   (R2) F_h is identical for two valid (t1, t2) pairs  — the
//        rotation of A_x must be a true tensor transformation
//
// If (R1) fails, flux_.Interior's Tinv'-then-T' composition leaks σ into
// spurious tensor components.  Specifically, a uniform global Q with
// SXZ = 0 might roundtrip to non-zero SXZ — which would cause exactly
// the observed dip-channel pepper.
//
// This test directly exercises (R1) and (R2) for the FaultBasis tangent
// pair used by the wave operator on the TPV102-like 8-fault-triangle
// fixture (can_n = +y, FaultBasis can_t1 = -z, can_t2 = -x).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/godunov_flux.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

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

namespace {

const char *StateName(int c)
{
   static const char *kNames[NUM_STATE] =
   {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ", "VX", "VY", "VZ"};
   return (c >= 0 && c < NUM_STATE) ? kNames[c] : "?";
}

void PrintQ(const char *label, const real_t *Q)
{
   std::cout << "  " << std::setw(20) << label << ":";
   for (int c = 0; c < NUM_STATE; c++)
   {
      std::cout << "  " << StateName(c) << "="
                << std::scientific << std::setprecision(3) << std::showpos
                << Q[c];
   }
   std::cout << std::noshowpos << "\n";
}

void Roundtrip(const char *label, const real_t *n, const real_t *t1,
               const real_t *t2, const real_t *Q_global, real_t *Q_back,
               real_t *Q_local_out = nullptr)
{
   DenseMatrix Tinv(NUM_STATE, NUM_STATE), T(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   GodunovFlux::BuildRotation       (n, t1, t2, T);
   Vector Q_g(const_cast<real_t*>(Q_global), NUM_STATE);
   Vector Q_l(NUM_STATE), Q_b(NUM_STATE);
   Tinv.Mult(Q_g, Q_l);
   T.Mult(Q_l, Q_b);
   for (int c = 0; c < NUM_STATE; c++) { Q_back[c] = Q_b(c); }
   if (Q_local_out)
   { for (int c = 0; c < NUM_STATE; c++) { Q_local_out[c] = Q_l(c); } }
   PrintQ((std::string("Q_global ") + label).c_str(), Q_global);
   real_t Q_local_local[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++) { Q_local_local[c] = Q_l(c); }
   PrintQ((std::string("Q_local  ") + label).c_str(), Q_local_local);
   PrintQ((std::string("Q_back   ") + label).c_str(), Q_back);
}

real_t MaxAbsDiff(const real_t *a, const real_t *b)
{
   real_t m = 0;
   for (int c = 0; c < NUM_STATE; c++)
   {
      m = std::max(m, std::abs(a[c] - b[c]));
   }
   return m;
}

} // anonymous

int main()
{
   std::cout << "\n=== Godunov-Flux rotation T(Tinv(σ)) = σ test ===\n";

   // --- Frame A: FaultBasis canonical for the 8-tet fixture ---
   // can_n = (0, 1, 0)   (y+ outward)
   // can_t1 = (0, 0, -1) (dip = -z)
   // can_t2 = (-1, 0, 0) (strike = -x)
   const real_t can_n[3]  = { 0,  1,  0};
   const real_t can_t1[3] = { 0,  0, -1};
   const real_t can_t2[3] = {-1,  0,  0};

   // --- Frame B: BuildFrame canonical for the same n ---
   // Inside GodunovFlux::Interior, BuildFrame(can_n) gives:
   real_t bf_t1[3], bf_t2[3];
   GodunovFlux::BuildFrame(can_n, bf_t1, bf_t2);
   std::cout << "  Frame A (FaultBasis): can_t1=(" << can_t1[0] << ","
             << can_t1[1] << "," << can_t1[2] << "), can_t2=("
             << can_t2[0] << "," << can_t2[1] << "," << can_t2[2] << ")\n";
   std::cout << "  Frame B (BuildFrame): bf_t1 =(" << bf_t1[0] << ","
             << bf_t1[1] << "," << bf_t1[2] << "), bf_t2 =("
             << bf_t2[0] << "," << bf_t2[1] << "," << bf_t2[2] << ")\n";

   // --- Test A: TPV102 background Q.  SYY = sigma_n, SXY = -tau_ini, others 0
   const real_t sigma_n = 1.20e8;
   const real_t tau_ini = 7.5e7;
   real_t Q_bg[NUM_STATE] = {0};
   Q_bg[SYY] =  sigma_n;
   Q_bg[SXY] = -tau_ini;

   std::cout << "\n-- TPV102 background round-trip in Frame A (FaultBasis) --\n";
   real_t Q_back_A[NUM_STATE], Q_local_A[NUM_STATE];
   Roundtrip("A", can_n, can_t1, can_t2, Q_bg, Q_back_A, Q_local_A);
   const real_t diff_A = MaxAbsDiff(Q_bg, Q_back_A);
   std::cout << "    max |Q_back - Q| = " << std::scientific
             << std::setprecision(6) << diff_A << "\n";
   TEST_LE(diff_A / std::max(sigma_n, tau_ini), 1.0e-13,
           "Frame A: T(Tinv(Q_bg)) == Q_bg (1e-13 relative)");

   std::cout << "\n-- TPV102 background round-trip in Frame B (BuildFrame) --\n";
   real_t Q_back_B[NUM_STATE], Q_local_B[NUM_STATE];
   Roundtrip("B", can_n, bf_t1, bf_t2, Q_bg, Q_back_B, Q_local_B);
   const real_t diff_B = MaxAbsDiff(Q_bg, Q_back_B);
   std::cout << "    max |Q_back - Q| = " << std::scientific
             << std::setprecision(6) << diff_B << "\n";
   TEST_LE(diff_B / std::max(sigma_n, tau_ini), 1.0e-13,
           "Frame B: T(Tinv(Q_bg)) == Q_bg (1e-13 relative)");

   // --- Test B: dirty input — SXZ should stay 0 in round-trip
   real_t Q_dirty[NUM_STATE] = {0};
   Q_dirty[SXX] = 5.0e6;  Q_dirty[SYY] = 1.2e8;  Q_dirty[SZZ] = 1.0e6;
   Q_dirty[SXY] = -7.5e7; Q_dirty[SYZ] = 2.0e6;  Q_dirty[SXZ] = 0.0;
   Q_dirty[VX]  = 0.1;    Q_dirty[VY]  = -0.2;   Q_dirty[VZ]  = 0.3;

   std::cout << "\n-- Dirty input round-trip in Frame A --\n";
   real_t Q_back_C[NUM_STATE], Q_local_C[NUM_STATE];
   Roundtrip("C", can_n, can_t1, can_t2, Q_dirty, Q_back_C, Q_local_C);
   const real_t diff_C = MaxAbsDiff(Q_dirty, Q_back_C);
   const real_t scale = 1.2e8;
   std::cout << "    max |Q_back - Q| = " << std::scientific
             << std::setprecision(6) << diff_C << "\n";
   TEST_LE(diff_C / scale, 1.0e-13,
           "Dirty Q round-trip is identity in Frame A");

   // --- Test C: cross-frame consistency.  Apply A_n_global to Q_bg via BOTH
   //     frames; results MUST match because A_n_global is frame-independent
   //     once n is fixed.  Mimics flux_.Interior(can_n, Q, Q, F_h).
   std::cout << "\n-- Cross-frame F_h consistency: flux_.Interior with both --\n";

   // We need a GodunovFlux instance.  TPV102 material.
   // CTOR signature: GodunovFlux(real_t lambda, real_t mu, real_t rho)
   const real_t rho = 2670.0, cp = 6000.0, cs = 3464.0;
   const real_t mu = rho * cs * cs;
   const real_t lambda = rho * cp * cp - 2.0 * mu;
   GodunovFlux flux(lambda, mu, rho);

   real_t F_h_A[NUM_STATE], F_h_B[NUM_STATE];
   // Frame A is implicit in flux.Interior(can_n, Q, Q, F_h) — Interior calls
   // BuildFrame(n) internally which is FRAME B.  So both calls actually
   // use BuildFrame.  To get a "Frame A" comparison I'd need to reach
   // inside GodunovFlux.  Instead, manually replicate flux_.Interior with
   // each of the two frames using BuildRotation/BuildRotationInverse +
   // ApplySplitFlux equivalent.
   //
   // For the purpose of this unit test, the key check is:  does
   // flux.Interior(can_n, Q_bg, Q_bg, F_h) produce F_h with F_h[SXZ] = 0?
   flux.Interior(can_n, Q_bg, Q_bg, F_h_A);
   PrintQ("F_h Interior(can_n, Q_bg)", F_h_A);
   const real_t F_h_SXZ_abs = std::abs(F_h_A[SXZ]);
   const real_t F_h_scale = 0;  // for this Q_bg, A_n·Q has only specific entries
   // For elastodynamic A_y · Q with Q=(SYY=σ_n, SXY=-τ, others 0),
   // expected non-zero F_h entries:
   //   F_h[VX] = -SXY / ρ = +τ_ini / ρ  ≈ +28083 m/s²
   //   F_h[VY] = -SYY / ρ = -σ_n / ρ    ≈ -44944 m/s²
   //   all others = 0
   const real_t expected_F_VX = -Q_bg[SXY] / rho;
   const real_t expected_F_VY = -Q_bg[SYY] / rho;
   std::cout << "    expected F_h[VX] = " << std::scientific
             << std::setprecision(6) << expected_F_VX
             << "  actual = " << F_h_A[VX] << "\n"
             << "    expected F_h[VY] = " << expected_F_VY
             << "  actual = " << F_h_A[VY] << "\n";
   TEST_LE(std::abs(F_h_A[VX] - expected_F_VX), 1.0e-3,
           "F_h[VX] = -SXY/ρ");
   TEST_LE(std::abs(F_h_A[VY] - expected_F_VY), 1.0e-3,
           "F_h[VY] = -SYY/ρ");
   TEST_LE(F_h_SXZ_abs, 1.0e-12,
           "F_h[SXZ] == 0 (no y-direction flux of σ_xz in elastodynamics)");

   // Also check other expected zeros:
   for (int c : {SXX, SXY, SXZ, SYY, SZZ, VZ})
   {
      if (c == SYY || c == SXY) { continue; }  // these have nonzero entries
      // Actually expected_F_h[SYY] = -(λ+2μ)·VY = 0 since VY=0. Same for SXY etc.
      const real_t v = std::abs(F_h_A[c]);
      std::cout << "    F_h[" << StateName(c) << "] = " << std::scientific
                << std::setprecision(3) << F_h_A[c] << "\n";
      TEST_LE(v, 1.0e-12,
              std::string("F_h[") + StateName(c) + "] == 0 for Q_bg input");
   }

   // -----------------------------------------------------------------
   // CRITICAL TEST: Does the spectral decomposition recover A_x?
   // I.e., does A_x_plus + A_x_minus == A_x_full (BuildJacobian)?
   //
   // The flux.Interior pipeline assumes:
   //   F_rot = ApplySplitFlux(Q_self, Q_nbr) computes
   //           Ax_plus_ * Q_self + Ax_minus_ * Q_nbr.
   //   When Q_self == Q_nbr, this = (Ax_plus_ + Ax_minus_) * Q.
   //   For physical correctness, this MUST equal Ax_full * Q.
   //
   // If the F_h[VX] discrepancy above is due to Ax_plus_ + Ax_minus_ !=
   // Ax_full, this test will reveal it.
   // -----------------------------------------------------------------
   std::cout << "\n-- Ax decomposition consistency test --\n";
   const DenseMatrix &Ax_full = flux.GetReferenceStarMatrix(0);
   // Apply Ax_full to the LOCAL Q in BuildFrame frame.
   // Step 1: rotate Q_bg into BuildFrame-local.
   DenseMatrix Tinv_bf(NUM_STATE), T_bf(NUM_STATE);
   GodunovFlux::BuildRotationInverse(can_n, bf_t1, bf_t2, Tinv_bf);
   GodunovFlux::BuildRotation       (can_n, bf_t1, bf_t2, T_bf);
   Vector Q_g_vec(const_cast<real_t*>(Q_bg), NUM_STATE);
   Vector Q_local_bf(NUM_STATE);
   Tinv_bf.Mult(Q_g_vec, Q_local_bf);
   // Step 2: Apply Ax_full directly.
   Vector F_local_full(NUM_STATE);
   Ax_full.Mult(Q_local_bf, F_local_full);
   // Step 3: Rotate back to global.
   Vector F_global_full(NUM_STATE);
   T_bf.Mult(F_local_full, F_global_full);

   std::cout << "  Ax_full (BuildJacobian) gives:\n";
   real_t F_full_arr[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++) { F_full_arr[c] = F_global_full(c); }
   PrintQ("F_h via Ax_full", F_full_arr);
   std::cout << "  flux.Interior (Ax_plus + Ax_minus) gives:\n";
   PrintQ("F_h via Interior", F_h_A);

   // Intermediate diagnostic: print F_local before rotating back, to
   // localize whether Ax_full is wrong OR T is wrong.
   real_t F_loc_arr[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++) { F_loc_arr[c] = F_local_full(c); }
   std::cout << "  F_local (Ax_full * Q_local, Frame B):\n";
   PrintQ("F_local", F_loc_arr);

   // Hand-expected per elastodynamics:
   //   For dir=0 (local x), V=0, Q_local = (SXX=σ_n, SXY=+τ, ...):
   //     F_local[VX] = -SXX/ρ = -σ_n/ρ
   //     F_local[VY] = -SXY/ρ = -τ/ρ
   //     all others = 0.
   const real_t F_loc_VX_expect = -1.20e8 / rho;
   const real_t F_loc_VY_expect = -7.50e7 / rho;
   std::cout << "  expected F_local[VX] = " << std::scientific
             << std::setprecision(6) << F_loc_VX_expect
             << "  actual = " << F_loc_arr[VX] << "\n";
   std::cout << "  expected F_local[VY] = " << F_loc_VY_expect
             << "  actual = " << F_loc_arr[VY] << "\n";
   TEST_LE(std::abs(F_loc_arr[VX] - F_loc_VX_expect) /
           std::abs(F_loc_VX_expect), 1.0e-12,
           "Ax_full * Q_local: F_local[VX] = -SXX/ρ");
   TEST_LE(std::abs(F_loc_arr[VY] - F_loc_VY_expect) /
           std::abs(F_loc_VY_expect), 1.0e-12,
           "Ax_full * Q_local: F_local[VY] = -SXY/ρ");

   // Direct dump of Ax_full entries (with full precision):
   std::cout << std::scientific << std::setprecision(15);
   std::cout << "  Ax_full row VX (the -1/rho row, expected entry at SXX): \n";
   for (int c = 0; c < NUM_STATE; c++)
   {
      if (std::abs(Ax_full(VX, c)) > 1e-30)
      {
         std::cout << "    Ax_full(VX, " << StateName(c) << ") = "
                   << Ax_full(VX, c) << "\n";
      }
   }
   std::cout << "  Ax_full row VY (the -1/rho row, expected entry at SXY): \n";
   for (int c = 0; c < NUM_STATE; c++)
   {
      if (std::abs(Ax_full(VY, c)) > 1e-30)
      {
         std::cout << "    Ax_full(VY, " << StateName(c) << ") = "
                   << Ax_full(VY, c) << "\n";
      }
   }
   std::cout << "  Q_local SXX = " << Q_local_bf(SXX)
             << "  SXY = " << Q_local_bf(SXY) << "\n";
   std::cout << "  Expected: 1/rho = " << (1.0/rho) << "\n";
   std::cout << std::setprecision(6);

   // Now print key entries of T to see if T(VX, VY) = -1 as I expect.
   std::cout << "  T_bf velocity block (should be Q^T of "
             << "[n; t1; t2]):\n";
   for (int gi = 0; gi < 3; gi++)
   {
      std::cout << "    T(V_global[" << gi << "], V_local[*]) = (";
      for (int lj = 0; lj < 3; lj++)
      {
         std::cout << T_bf(VX + gi, VX + lj) << (lj < 2 ? ", " : ")");
      }
      std::cout << "\n";
   }

   real_t max_diff = 0;
   for (int c = 0; c < NUM_STATE; c++)
   {
      max_diff = std::max(max_diff, std::abs(F_full_arr[c] - F_h_A[c]));
   }
   std::cout << "  max |F_h_full - F_h_Interior| = " << std::scientific
             << std::setprecision(6) << max_diff << "\n";
   TEST_LE(max_diff / std::max(std::abs(F_full_arr[VX]), 1.0),
           1.0e-10,
           "Ax_plus + Ax_minus == Ax_full (spectral decomposition is "
           "correct decomposition of A_x)");

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
