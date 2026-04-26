// Round-7 R-602/R-603: Tpv104SubStepIterator parity contract at O=1.
//
// Contract:
//   At sub-step quadrature O = 1 with deltaT = {dt_macro},
//   time_weights = {1.0}, the new `AdvanceWithSubStepStates` method
//   (per-sub-step pointwise Q) and the legacy `Advance` method
//   (time-integrated I) must produce IDENTICAL DOFData updates and
//   I_imp_flat outputs when:
//      Q_pointwise[0]  ==  I_flat / dt_macro
//   on every QP.  This is the SeisSol-equivalent reduction at O = 1
//   — both APIs reduce to a single ComputeStageState call on Q̄.
//
//   We verify this by:
//     1. Picking a representative DOFData state (TPV102 locked-fault
//        snapshot at rest + some perturbation in V_ini).
//     2. Running Advance on a fixed I_flat.
//     3. Running AdvanceWithSubStepStates with Q_pointwise[0] = I_flat / dt
//        on a FRESH dof_data copy.
//     4. Asserting both DOFData snapshots and both I_imp outputs match.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv104_substep_iterator.hpp"
#include "../../friction/slip_law_srw_psi.hpp"
#include "../../config/tpv104_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

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

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace
{

// Configure a representative TPV104-style DOFData at rest with a small
// perturbation so the friction solve is exercised non-trivially.
DOFData MakeTPV104DOF()
{
   DOFData d;
   d.Zp_plus  = TPV104Params::Zp;
   d.Zp_minus = TPV104Params::Zp;
   d.Zs_plus  = TPV104Params::Zs;
   d.Zs_minus = TPV104Params::Zs;
   d.eta_p    = TPV104Params::eta_p;
   d.eta_s    = TPV104Params::eta_s;
   d.sigma_n0 = TPV104Params::sigma_n;
   d.tau1_0   = 0.0;                      // pure strike-slip
   d.tau2_0   = 29.38e6;                  // representative pre-stress
   d.sigma_n_nuc = 0.0;
   d.tau1_nuc    = 0.0;
   d.tau2_nuc    = 0.0;                   // nucleation off for parity test
   d.a   = TPV104Params::a_in;
   d.Dc  = TPV104Params::L;
   d.psi = 0.564;                          // representative steady-state ψ
   d.slip_rate = 0.0;
   d.V1 = 0.0; d.V2 = 0.0;
   d.slip1 = 0.0; d.slip2 = 0.0;
   d.tau1_corr = 0.0; d.tau2_corr = 0.0;
   d.sigma_n_corr = 0.0;
   return d;
}

real_t MaxAbsDiff(const std::vector<real_t> &a,
                  const std::vector<real_t> &b)
{
   const size_t n = std::min(a.size(), b.size());
   real_t m = 0.0;
   for (size_t i = 0; i < n; i++)
   {
      real_t d = std::abs(a[i] - b[i]);
      if (d > m) { m = d; }
   }
   return m;
}

real_t MaxDOFFieldDiff(const DOFData &a, const DOFData &b)
{
   real_t m = 0.0;
   m = std::max(m, std::abs(a.psi          - b.psi));
   m = std::max(m, std::abs(a.slip_rate    - b.slip_rate));
   m = std::max(m, std::abs(a.V1           - b.V1));
   m = std::max(m, std::abs(a.V2           - b.V2));
   m = std::max(m, std::abs(a.slip1        - b.slip1));
   m = std::max(m, std::abs(a.slip2        - b.slip2));
   m = std::max(m, std::abs(a.tau1_corr    - b.tau1_corr));
   m = std::max(m, std::abs(a.tau2_corr    - b.tau2_corr));
   m = std::max(m, std::abs(a.sigma_n_corr - b.sigma_n_corr));
   return m;
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-7 R-602/R-603: SubStep iterator parity at O=1 ===\n";

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);

   // SlipLawSRWPsi for global friction scalars (b, V0, f0, muW); a_scalar
   // and V_w_default are not consulted in production-mode SRW free-function
   // dispatch (per-DOF a / V_w shadow them), so any nonzero placeholder
   // value is fine.
   SlipLawSRWPsi state_evo(
      TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
      TPV104Params::f0, TPV104Params::muW, TPV104Params::V_w_out);
   state_evo.SetProductionMode();

   Tpv104SubStepIterator iter(flux, state_evo);
   const real_t dt_macro = 1e-4;
   iter.SetSubSteps({dt_macro}, {1.0});

   // Single-DOF fixture.
   const int n = 1;
   std::vector<DOFData> dof_data_legacy = { MakeTPV104DOF() };
   std::vector<DOFData> dof_data_substep = { MakeTPV104DOF() };
   std::vector<Vector> fault_coords(n, Vector(3));
   fault_coords[0] = 0.0;
   fault_coords[0](2) = -7500.0;            // hypocenter z
   std::vector<real_t> V_w(n, TPV104Params::V_w_in);

   // Construct a representative fault-local-frame Q̄ at the fault QP.
   // Pure strike-slip rest state with a small VZ kick.
   std::vector<real_t> Q_avg_plus(NUM_STATE, 0.0);
   std::vector<real_t> Q_avg_minus(NUM_STATE, 0.0);
   Q_avg_plus[VZ]  = +0.001;     // small + side strike velocity (m/s)
   Q_avg_minus[VZ] = -0.001;     // small − side strike velocity (m/s)

   // Legacy I_flat = dt * Q̄.  AdvanceWithSubStepStates expects the
   // pointwise Q̄ directly (already at sub-step time, no integration scaling).
   std::vector<real_t> I_plus_flat(NUM_STATE * n, 0.0);
   std::vector<real_t> I_minus_flat(NUM_STATE * n, 0.0);
   for (int c = 0; c < NUM_STATE; c++)
   {
      I_plus_flat [c]  = dt_macro * Q_avg_plus[c];
      I_minus_flat[c]  = dt_macro * Q_avg_minus[c];
   }

   // Outputs.
   std::vector<real_t> I_imp_plus_legacy (NUM_STATE * n, 0.0);
   std::vector<real_t> I_imp_minus_legacy(NUM_STATE * n, 0.0);
   std::vector<real_t> I_imp_plus_substep (NUM_STATE * n, 0.0);
   std::vector<real_t> I_imp_minus_substep(NUM_STATE * n, 0.0);

   // (1) Legacy Advance — operates on I (time-integrated) and divides
   //     internally by dt_macro to get Q̄.
   iter.Advance(dof_data_legacy, fault_coords, V_w,
                I_plus_flat.data(), I_minus_flat.data(),
                dt_macro, /*t_macro_start=*/0.0,
                I_imp_plus_legacy.data(), I_imp_minus_legacy.data(),
                FrictionSolver::Method::NewtonRaphsonStable);

   // (2) AdvanceWithSubStepStates — operates on per-sub-step pointwise Q̄
   //     directly.  At O = 1, single sub-step with the SAME Q̄ as the
   //     legacy path's I/dt_macro reduction.
   std::vector<std::vector<real_t>> Q_plus_per({Q_avg_plus});
   std::vector<std::vector<real_t>> Q_minus_per({Q_avg_minus});
   iter.AdvanceWithSubStepStates(dof_data_substep, fault_coords, V_w,
                                 Q_plus_per, Q_minus_per,
                                 dt_macro, /*t_macro_start=*/0.0,
                                 I_imp_plus_substep.data(),
                                 I_imp_minus_substep.data(),
                                 FrictionSolver::Method::NewtonRaphsonStable);

   // ------------------------------------------------------------------
   // Gate 1: I_imp outputs match.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 1: I_imp outputs at O=1 --\n";
   const real_t err_I_plus  =
      MaxAbsDiff(I_imp_plus_legacy,  I_imp_plus_substep);
   const real_t err_I_minus =
      MaxAbsDiff(I_imp_minus_legacy, I_imp_minus_substep);
   std::cout << "  max |I_imp_plus_legacy − I_imp_plus_substep|   = "
             << std::scientific << err_I_plus  << "\n";
   std::cout << "  max |I_imp_minus_legacy − I_imp_minus_substep| = "
             << std::scientific << err_I_minus << "\n";
   TEST_LE(err_I_plus,  0.0,
           "I_imp_plus parity at O=1 (bit-equal)");
   TEST_LE(err_I_minus, 0.0,
           "I_imp_minus parity at O=1 (bit-equal)");

   // ------------------------------------------------------------------
   // Gate 2: DOFData updates match.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 2: DOFData updates at O=1 --\n";
   const real_t err_dof =
      MaxDOFFieldDiff(dof_data_legacy[0], dof_data_substep[0]);
   std::cout << "  max DOFData field diff = " << std::scientific
             << err_dof << "\n";
   std::cout << "    psi(legacy/substep)        = "
             << dof_data_legacy[0].psi << " / "
             << dof_data_substep[0].psi << "\n";
   std::cout << "    slip_rate(legacy/substep)  = "
             << dof_data_legacy[0].slip_rate << " / "
             << dof_data_substep[0].slip_rate << "\n";
   std::cout << "    V2(legacy/substep)         = "
             << dof_data_legacy[0].V2 << " / "
             << dof_data_substep[0].V2 << "\n";
   TEST_LE(err_dof, 0.0,
           "DOFData fields parity at O=1 (bit-equal)");

   // ------------------------------------------------------------------
   // Gate 3: Both APIs called with reasonable-magnitude inputs do not
   //         produce NaN/Inf in any output.  Defensive smoke.
   // ------------------------------------------------------------------
   std::cout << "\n-- Gate 3: outputs finite --\n";
   bool all_finite = true;
   for (real_t v : I_imp_plus_substep)  { if (!std::isfinite(v)) { all_finite = false; } }
   for (real_t v : I_imp_minus_substep) { if (!std::isfinite(v)) { all_finite = false; } }
   if (!std::isfinite(dof_data_substep[0].slip_rate) ||
       !std::isfinite(dof_data_substep[0].V1) ||
       !std::isfinite(dof_data_substep[0].V2) ||
       !std::isfinite(dof_data_substep[0].psi))
   {
      all_finite = false;
   }
   TEST_ASSERT(all_finite,
               "AdvanceWithSubStepStates produces finite outputs");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
