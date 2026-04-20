// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.0.0 Phase 1 §3.1e — primitive identity for the Pelties-9 fix.
//
// Verifies:
//     flux_.Interior(n, Q, Q)  ==  T . A_x . T^{-1} . Q
//                              ==  A_{n} . Q
// for arbitrary unit normals n and arbitrary physically-sane Q.
//
// This is the algebraic identity that the v9.0.0 per-side fix (plan
// §14.2 and §14.3 in
// debug_document/tpv102_debug_document/tpv102_debug_v9.0.0_seissol_flux_comparison.md)
// depends on: calling Interior with Q_L = Q_R collapses
//   T * (A_x^+ * Q_rot + A_x^- * Q_rot)  ==  T * A_x * Q_rot
// because the characteristic decomposition gives A_x^+ + A_x^- = A_x
// exactly (see godunov_flux.cpp:94-105: Lambda_plus + Lambda_minus ==
// full Lambda eigenvalue matrix by construction).
//
// Any failure here means the identity is broken and the v9.0.0 fix is
// not a drop-in for Pelties 2012 eq. (9) -- reject the fix until
// understood.

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define T_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ \
      << ", |diff| " << std::abs(v_-e_) << ", tol " << t_ << ")\n"; } \
} while(0)

// Build the full x-direction Jacobian externally and evaluate
// F_expected = T * A_x * T^{-1} * Q in the global frame.  This
// replicates Interior's algorithm but uses the FULL (not split) A_x,
// so agreement with Interior(n, Q, Q) confirms A_x^+ + A_x^- = A_x.
static void ExpectedFullFlux(const GodunovFlux &flux,
                             const real_t *n, const real_t *Q, real_t *F)
{
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(n, t1, t2);
   DenseMatrix T(NUM_STATE), Tinv(NUM_STATE), Ax(NUM_STATE);
   GodunovFlux::BuildRotation(n, t1, t2, T);
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   flux.BuildJacobian(0, Ax);  // full A_x, not A_x^+/A_x^-

   real_t Q_rot[NUM_STATE], F_rot[NUM_STATE];
   Tinv.Mult(Q, Q_rot);
   Ax.Mult(Q_rot, F_rot);
   T.Mult(F_rot, F);
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.0.0 Phase 1 §3.1e — identity test\n";
   std::cout << "  flux_.Interior(n, Q, Q) == T . A_x . T^{-1} . Q\n";
   std::cout << "========================================\n\n";

   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu,
                    TPV102Params::rho);

   // 5 normals — covers three coordinate axes (incl. z-axis which
   // exercises BuildFrame's up=(1,0,0) fallback), one oblique, one
   // tilted in the xy-plane.
   real_t normals[5][3] = {
      { 1.0,  0.0,  0.0},
      { 0.0,  1.0,  0.0},
      { 0.0,  0.0,  1.0},
      { 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)},
      { 0.6, -0.8,  0.0}
   };

   // 3 physically-sane Q samples:
   //   (a) stress-only — probes A[v, sigma] rows
   //   (b) velocity-only — probes A[sigma, v] rows
   //   (c) mixed physical magnitudes
   real_t Q_samples[3][NUM_STATE] = {
      {1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.5, 0.3},
      {1.0e6, -2.0e6, 0.5e6, 3.0e5, -4.0e5, 6.0e5, 0.1, -0.2, 0.15}
   };

   int case_id = 0;
   for (int ni = 0; ni < 5; ni++)
   {
      for (int qi = 0; qi < 3; qi++)
      {
         real_t F_interior[NUM_STATE], F_expected[NUM_STATE];
         flux.Interior(normals[ni], Q_samples[qi], Q_samples[qi],
                       F_interior);
         ExpectedFullFlux(flux, normals[ni], Q_samples[qi], F_expected);

         // Tolerance per case: the absolute FP residue on a component
         // that should be exactly 0 is bounded by
         //     eps * max|Q| * max|A_x entry|
         // where max|A_x| ~ lambda+2*mu ~ O(1e10) for typical rock.
         // With Q at stress scale O(1e7), the residue floor is
         // O(eps * 1e7 * 1e10) ~ O(1e-1), observed at O(1e-6).  Use a
         // per-case scaling that also includes a |Q|_inf factor to
         // track input scale.  Relative tolerance 1e-10 * F_mag
         // catches drift on non-zero components.
         real_t F_mag = 0.0, Q_mag = 0.0;
         for (int c = 0; c < NUM_STATE; c++)
         {
            F_mag = std::max(F_mag, std::abs(F_expected[c]));
            Q_mag = std::max(Q_mag, std::abs(Q_samples[qi][c]));
         }
         // Scale floor: 1e-10 * Q_mag for stress-scale (Q~1e7 -> 1e-3),
         // with an absolute bottom of 1e-10 for velocity-only cases
         // (Q~1 -> floor ~ 1e-10).  Overall:
         const real_t case_tol = std::max(
            std::max(1.0e-10, 1.0e-10 * Q_mag),
            1.0e-10 * F_mag);

         for (int c = 0; c < NUM_STATE; c++)
         {
            T_NEAR(F_interior[c], F_expected[c], case_tol,
                   "case " + std::to_string(case_id)
                   + " ni=" + std::to_string(ni)
                   + " qi=" + std::to_string(qi)
                   + " c="  + std::to_string(c));
         }
         case_id++;
      }
   }

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";
   return (num_failed > 0) ? 1 : 0;
}
