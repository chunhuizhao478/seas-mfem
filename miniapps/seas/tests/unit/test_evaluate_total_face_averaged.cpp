// Unit test for FaultFaceFlux::EvaluateTotalFaceAveraged (option 1
// of the 2026-04-22 pepper fix).
//
// Validates:
//   1. With UNIFORM per-QP DOFData inputs, the face-averaged variant
//      gives the SAME outputs as per-QP EvaluateTotal called on any
//      one QP.  (Sanity: averaging uniform values is a no-op.)
//
//   2. With NON-UNIFORM per-QP DOFData (different psi or tau*_nuc per
//      QP), the face-averaged variant writes the SAME slip_rate,
//      tau*_corr to every QP.  Per-QP psi is preserved.
//
//   3. EvaluateTotal contract is preserved: the friction call inside
//      sees zero static prestress (R-I06-007).

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_=(val), e_=(exp), t_=(tol); \
   if (std::abs(v_-e_) < t_) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << " (got " << std::scientific << std::setprecision(6) \
                << v_ << ", expected " << e_ << ", |diff| " \
                << std::abs(v_-e_) << ", tol " << t_ << ")\n"; } \
} while (0)

namespace {
DOFData MakeLockedDOF(real_t a, real_t psi)
{
   DOFData d;
   d.Zp_plus = TPV102Params::Zp; d.Zp_minus = TPV102Params::Zp;
   d.Zs_plus = TPV102Params::Zs; d.Zs_minus = TPV102Params::Zs;
   d.eta_p = TPV102Params::Zp / 2.0;
   d.eta_s = TPV102Params::eta_s;
   d.sigma_n0 = 0;  d.tau1_0 = 0;  d.tau2_0 = 0;
   d.sigma_n_nuc = 0;  d.tau1_nuc = 0;  d.tau2_nuc = 0;
   d.a = a;  d.Dc = TPV102Params::Dc;
   d.psi = psi;
   d.slip_rate = TPV102Params::V_ini;
   d.V1 = 0; d.V2 = TPV102Params::V_ini;
   return d;
}
} // anonymous

int main()
{
   std::cout << "\n=== EvaluateTotalFaceAveraged unit test ===\n";

   FaultFaceFlux flux(TPV102Params::lambda, TPV102Params::mu,
                       TPV102Params::rho);

   const real_t a = TPV102Params::a_vw;
   const real_t psi_eq = ComputeInitialPsi(a);

   // Symmetric Q± at TPV102 equilibrium in fault-local frame
   // (SXX = sigma_n, SXZ = +tau_ini, others = 0).
   real_t Q[NUM_STATE] = {0};
   Q[SXX] = TPV102Params::sigma_n;
   Q[SXZ] = TPV102Params::tau_ini;

   const int nqp = 3;

   // --- Test 1: uniform per-QP DOFData → face-averaged matches per-QP ---
   {
      std::cout << "\n-- Test 1: uniform per-QP inputs --\n";
      DOFData per_qp[3] = {MakeLockedDOF(a, psi_eq),
                            MakeLockedDOF(a, psi_eq),
                            MakeLockedDOF(a, psi_eq)};
      // Apply nucleation drive
      for (int q = 0; q < nqp; q++)
      { per_qp[q].tau2_nuc = TPV102Params::nuc_dtau; }

      // Reference: per-QP EvaluateTotal on q=0
      DOFData ref = per_qp[0];
      real_t Qip_ref[NUM_STATE], Qim_ref[NUM_STATE];
      flux.EvaluateTotal(ref, Q, Q, Qip_ref, Qim_ref);

      // Face-averaged variant
      DOFData fa[3] = {per_qp[0], per_qp[1], per_qp[2]};
      real_t Qip_fa[NUM_STATE], Qim_fa[NUM_STATE];
      flux.EvaluateTotalFaceAveraged(fa, nqp, Q, Q, Qip_fa, Qim_fa);

      // Per-QP outputs should match the reference exactly
      for (int q = 0; q < nqp; q++)
      {
         TEST_NEAR(fa[q].slip_rate, ref.slip_rate, 1e-12,
                   std::string("T1 q=") + std::to_string(q) +
                   ": slip_rate matches per-QP reference");
         TEST_NEAR(fa[q].tau2_corr, ref.tau2_corr, 1e-6,
                   std::string("T1 q=") + std::to_string(q) +
                   ": tau2_corr matches per-QP reference");
      }
      // Q_imp matches reference
      for (int c = 0; c < NUM_STATE; c++)
      {
         TEST_NEAR(Qip_fa[c], Qip_ref[c], 1e-6,
                   std::string("T1 Q_imp_plus[") +
                   std::to_string(c) + "] matches reference");
      }
   }

   // --- Test 2: NON-UNIFORM per-QP psi → outputs are uniform ---
   {
      std::cout << "\n-- Test 2: non-uniform per-QP psi --\n";
      DOFData per_qp[3] = {MakeLockedDOF(a, psi_eq),
                            MakeLockedDOF(a, psi_eq + 0.001),
                            MakeLockedDOF(a, psi_eq - 0.001)};
      for (int q = 0; q < nqp; q++)
      { per_qp[q].tau2_nuc = TPV102Params::nuc_dtau; }
      DOFData fa[3] = {per_qp[0], per_qp[1], per_qp[2]};
      real_t Qip[NUM_STATE], Qim[NUM_STATE];
      flux.EvaluateTotalFaceAveraged(fa, nqp, Q, Q, Qip, Qim);

      // After call: per-QP psi is preserved (we don't write it),
      // but slip_rate, tau*_corr, V1, V2 are uniform across QPs.
      TEST_NEAR(fa[0].psi, psi_eq, 1e-15,
                "T2 q=0: psi preserved (= psi_eq)");
      TEST_NEAR(fa[1].psi, psi_eq + 0.001, 1e-15,
                "T2 q=1: psi preserved (= psi_eq + 0.001)");
      TEST_NEAR(fa[2].psi, psi_eq - 0.001, 1e-15,
                "T2 q=2: psi preserved (= psi_eq - 0.001)");
      // slip_rate uniform
      TEST_NEAR(fa[1].slip_rate, fa[0].slip_rate, 1e-15,
                "T2: slip_rate q=1 == q=0 (face-uniform)");
      TEST_NEAR(fa[2].slip_rate, fa[0].slip_rate, 1e-15,
                "T2: slip_rate q=2 == q=0 (face-uniform)");
      TEST_NEAR(fa[1].tau2_corr, fa[0].tau2_corr, 1e-15,
                "T2: tau2_corr q=1 == q=0 (face-uniform)");
      TEST_NEAR(fa[2].tau2_corr, fa[0].tau2_corr, 1e-15,
                "T2: tau2_corr q=2 == q=0 (face-uniform)");
   }

   // --- Test 3: NON-UNIFORM per-QP tau2_nuc → outputs use the average ---
   {
      std::cout << "\n-- Test 3: non-uniform per-QP tau2_nuc --\n";
      DOFData per_qp[3] = {MakeLockedDOF(a, psi_eq),
                            MakeLockedDOF(a, psi_eq),
                            MakeLockedDOF(a, psi_eq)};
      // tau2_nuc differs per QP
      per_qp[0].tau2_nuc = 1.0e6;
      per_qp[1].tau2_nuc = 2.0e6;
      per_qp[2].tau2_nuc = 3.0e6;
      const real_t avg_tau2_nuc = 2.0e6;  // (1+2+3)/3

      // Reference: per-QP EvaluateTotal with tau2_nuc = avg
      DOFData ref = per_qp[0];
      ref.tau2_nuc = avg_tau2_nuc;
      real_t Qip_ref[NUM_STATE], Qim_ref[NUM_STATE];
      flux.EvaluateTotal(ref, Q, Q, Qip_ref, Qim_ref);

      DOFData fa[3] = {per_qp[0], per_qp[1], per_qp[2]};
      real_t Qip[NUM_STATE], Qim[NUM_STATE];
      flux.EvaluateTotalFaceAveraged(fa, nqp, Q, Q, Qip, Qim);
      TEST_NEAR(fa[0].slip_rate, ref.slip_rate, 1e-12,
                "T3: face-avg slip_rate matches per-QP-with-avg-nuc");
      TEST_NEAR(fa[0].tau2_corr, ref.tau2_corr, 1e-6,
                "T3: face-avg tau2_corr matches per-QP-with-avg-nuc");
   }

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
