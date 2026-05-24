// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// PLAN_predictor_vs_macro_diag_2026-05-23 §Phase 1 AC#3.
//
// Guards the ARITHMETIC of the [MACRO] decomposition (in
// ComputeADERSharedFaceFluxRHS, wave_operator.inl): the trace splits the macro
// solve's trial normal traction as
//      sn_vjump = eta_p·(Q_avg_minus[VX] − Q_avg_plus[VX])
//      sn_sterm = eta_p·(Q_avg_plus[SXX]/Zp⁺ + Q_avg_minus[SXX]/Zp⁻)
// (Q_avg = I/dt on the shared path) and claims sn_vjump + sn_sterm ==
// sigma_n_trial (ComputeTrialTraction Eq. 7a).  This test calls the REAL
// production FaultFaceFlux::ComputeTrialTraction and asserts the identity to
// <= 16 ULP, incl. Zp⁺ != Zp⁻, so a wrong index/sign in the trace's
// decomposition (e.g. plus/minus swap, SXX vs VX) fails here.
//
// SCOPE: this guards the FORMULA ONLY.  The call-site placement (that [MACRO]
// is wired into ComputeADERSharedFaceFluxRHS and fires at the shared seed QP)
// is NOT exercised here — it is verified by the R-001 integration check on
// Frontera (grep '^[MACRO] .* c=(6075' with is_shared dof indices).

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

static constexpr double kEps = 2.220446049250313e-16;

static void CheckDecomposition(const std::string &label,
                               real_t Zp_p, real_t Zp_m,
                               const real_t *Qp, const real_t *Qm)
{
   DOFData data;
   data.Zp_plus  = Zp_p;
   data.Zp_minus = Zp_m;
   data.Zs_plus  = 0.6 * Zp_p;   // arbitrary positive S-impedances (for tau*)
   data.Zs_minus = 0.6 * Zp_m;
   // eta_p is consumed directly by ComputeTrialTraction; the trace uses the
   // SAME data.eta_p, so the identity is index/sign-sensitive but eta_p-value
   // independent.  Use the harmonic mean (the physical value).
   data.eta_p = (Zp_p * Zp_m) / (Zp_p + Zp_m);
   data.eta_s = (data.Zs_plus * data.Zs_minus)
                / (data.Zs_plus + data.Zs_minus);

   real_t sigma_n_trial = 0, tau1_trial = 0, tau2_trial = 0;
   FaultFaceFlux::ComputeTrialTraction(data, Qp, Qm,
                                       sigma_n_trial, tau1_trial, tau2_trial);

   // The trace's decomposition (verbatim from wave_operator.inl [MACRO]):
   const real_t sn_vjump = data.eta_p * (Qm[VX] - Qp[VX]);
   const real_t sn_sterm = data.eta_p * (Qp[SXX] / data.Zp_plus
                                         + Qm[SXX] / data.Zp_minus);
   const real_t sum = sn_vjump + sn_sterm;

   const double scale = std::max<double>(1.0, std::abs(sigma_n_trial));
   const double rel = std::abs(sum - sigma_n_trial) / scale;
   std::cout << "    [" << label << "] sn_vjump=" << sn_vjump
             << " sn_sterm=" << sn_sterm << " sum=" << sum
             << " sigma_n_trial=" << sigma_n_trial << " rel=" << rel << "\n";
   TEST_ASSERT(rel <= 16 * kEps,
               label + ": sn_vjump+sn_sterm == sigma_n_trial (rel "
               + std::to_string(rel) + " <= " + std::to_string(16*kEps) + ")");
}

int main()
{
   std::cout << "=== test_macro_diag_decomposition_identity ===\n";

   // Distinct stresses AND velocities on each side; a plus/minus swap or a
   // SXX<->VX index error would break the identity.
   real_t Qp[NUM_STATE] = { 4.0e7, -2.0, 3.0, 0.5, -0.7, 0.9, 5.2, -5.0, 6.0 };
   real_t Qm[NUM_STATE] = { 5.0e7,  2.5, -3.5, 0.2, 0.3, -0.4, -1.3, 7.0, -8.0 };

   // Matched impedances (SAFS single-material case).
   const real_t Zp = 2670.0 * 5996.25;   // rho * cp
   CheckDecomposition("matched-Zp", Zp, Zp, Qp, Qm);

   // Mismatched impedances (general bimaterial — exercises both 1/Zp± terms).
   CheckDecomposition("bimaterial-Zp", Zp, 1.7 * Zp, Qp, Qm);
   CheckDecomposition("bimaterial-Zp-rev", 1.7 * Zp, Zp, Qp, Qm);

   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed ===\n";
   return (num_failed == 0) ? 0 : 1;
}
