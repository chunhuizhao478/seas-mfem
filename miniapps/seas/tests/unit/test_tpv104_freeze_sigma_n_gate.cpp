// Regression test for the SEAS_TPV104_FREEZE_SIGMA_N env-var gate in
// FaultFaceFlux::ComputeStageState.  Catches the fluctuation-mode
// doubling bug observed on Frontera job 7677661 (sigma_n_corr = 240 MPa
// instead of 120 MPa) by exercising the gate in BOTH modes:
//   - fluctuation-Q (TPV104 production): data.sigma_n0 = 120 MPa,
//     bulk Q[SXX] = 0 → naive override of s.sigma_n_trial would give
//     sigma_n_total = 240 MPa.  The corrected gate must produce 120.
//   - total-Q (TPV102 v9.3.0): data.sigma_n0 = 0, bulk Q[SXX] = 120e6
//     → the same env value must still produce sigma_n_total = 120 MPa.
//
// The gate is exercised end-to-end via ComputeStageState (the chokepoint
// reached by Evaluate / EvaluateADER on the live TPV104 dispatch path).
//
// Built locally; no MPI, no mesh, no time loop.  Per CLAUDE.md
// "feedback_no_local_reproducer.md", this is a single-QP fixture.

#include "mfem.hpp"
#include "../../config/tpv104_params.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) <= t_) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg \
                << " (got " << v_ << ", expected " << e_ << ")\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << " (got " << v_ << ", expected " << e_ \
                << ", diff " << std::abs(v_ - e_) \
                << ", tol " << t_ << ")\n"; \
   } \
} while(0)

static DOFData MakeFluctuationModeDOF()
{
   // TPV104 production: pre-stress lives in DOFData, bulk Q is fluctuation.
   DOFData d;
   d.Zp_plus  = TPV104Params::Zp;  d.Zp_minus = TPV104Params::Zp;
   d.Zs_plus  = TPV104Params::Zs;  d.Zs_minus = TPV104Params::Zs;
   d.eta_p    = TPV104Params::eta_p;
   d.eta_s    = TPV104Params::eta_s;
   d.sigma_n0 = TPV104Params::sigma_n;   // 120 MPa
   d.tau1_0   = 0.0;
   d.tau2_0   = TPV104Params::tau_ini;
   d.sigma_n_nuc = 0.0;
   d.tau1_nuc = 0.0;
   d.tau2_nuc = 0.0;
   d.psi      = 0.5;
   d.a        = 0.008;
   d.Dc       = TPV104Params::L;
   return d;
}

static DOFData MakeTotalModeDOF()
{
   // TPV102 v9.3.0 total-Q: DOFData pre-stress is zeroed; bulk Q carries
   // the pre-stress tensor.
   DOFData d = MakeFluctuationModeDOF();
   d.sigma_n0 = 0.0;
   d.tau1_0   = 0.0;
   d.tau2_0   = 0.0;
   return d;
}

int main()
{
   std::cout << "\n=== TPV104 σ_n freeze-gate regression "
             << "(fluctuation + total mode) ===\n";

   FaultFaceFlux ff(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);

   // -----------------------------------------------------------------
   // Test 1 — fluctuation mode at REST (TPV104 init).  Bulk Q = 0.
   // The naive override would give sigma_n_total = 240 MPa.  The fix
   // must give 120 MPa.
   // -----------------------------------------------------------------
   {
      DOFData d = MakeFluctuationModeDOF();
      real_t Qp[NUM_STATE] = {0}, Qm[NUM_STATE] = {0};

      setenv("SEAS_TPV104_FREEZE_SIGMA_N", "120e6", 1);
      EvalStageState s;
      ff.ComputeStageState(d, Qp, Qm, s);
      unsetenv("SEAS_TPV104_FREEZE_SIGMA_N");

      TEST_NEAR(s.sigma_n_total, 120.0e6, 1.0,
                "fluctuation-Q: =120e6 pins sigma_n_total to 120 MPa "
                "(was 240 MPa pre-fix on job 7677661)");
   }

   // -----------------------------------------------------------------
   // Test 2 — total mode at REST.  Bulk Q[SXX] carries the pre-stress.
   // -----------------------------------------------------------------
   {
      DOFData d = MakeTotalModeDOF();
      real_t Qp[NUM_STATE] = {0}, Qm[NUM_STATE] = {0};
      Qp[SXX] = TPV104Params::sigma_n;
      Qm[SXX] = TPV104Params::sigma_n;

      setenv("SEAS_TPV104_FREEZE_SIGMA_N", "120e6", 1);
      EvalStageState s;
      ff.ComputeStageState(d, Qp, Qm, s);
      unsetenv("SEAS_TPV104_FREEZE_SIGMA_N");

      TEST_NEAR(s.sigma_n_total, 120.0e6, 1.0,
                "total-Q: =120e6 pins sigma_n_total to 120 MPa "
                "(no double-count from DOFData pre-stress)");
   }

   // -----------------------------------------------------------------
   // Test 3 — fluctuation mode with a perturbed bulk Q.  The override
   // must still pin sigma_n_total exactly at 120 MPa regardless of the
   // bulk Q[SXX] swing (which would otherwise drive sigma_n_trial).
   // -----------------------------------------------------------------
   {
      DOFData d = MakeFluctuationModeDOF();
      real_t Qp[NUM_STATE] = {0}, Qm[NUM_STATE] = {0};
      Qp[SXX] = -3.0e6;   // simulate a 3 MPa fluctuation swing
      Qm[SXX] = +3.0e6;
      Qp[VX]  = -1.0e-3;
      Qm[VX]  = +1.0e-3;

      setenv("SEAS_TPV104_FREEZE_SIGMA_N", "120e6", 1);
      EvalStageState s;
      ff.ComputeStageState(d, Qp, Qm, s);
      unsetenv("SEAS_TPV104_FREEZE_SIGMA_N");

      TEST_NEAR(s.sigma_n_total, 120.0e6, 1.0,
                "fluctuation-Q with bulk perturbation: override pins "
                "sigma_n_total at 120 MPa regardless of Q-swing");
   }

   // -----------------------------------------------------------------
   // Test 4 — env unset: byte-identical to baseline ComputeTrialTraction.
   // sigma_n_total at TPV104 rest must equal sigma_n0 = 120 MPa.
   // -----------------------------------------------------------------
   {
      DOFData d = MakeFluctuationModeDOF();
      real_t Qp[NUM_STATE] = {0}, Qm[NUM_STATE] = {0};

      unsetenv("SEAS_TPV104_FREEZE_SIGMA_N");
      EvalStageState s;
      ff.ComputeStageState(d, Qp, Qm, s);

      TEST_NEAR(s.sigma_n_total, 120.0e6, 1e-6,
                "env unset baseline: sigma_n_total = sigma_n0 + 0 + 0 = 120 MPa");
      TEST_NEAR(s.sigma_n_trial, 0.0, 1e-6,
                "env unset baseline: sigma_n_trial = 0 (Q at rest)");
   }

   // -----------------------------------------------------------------
   // Test 5 — env = 50e6: pin to 50 MPa in fluctuation mode (proves
   // the override picks up arbitrary positive values, not hard-coded 120).
   // -----------------------------------------------------------------
   {
      DOFData d = MakeFluctuationModeDOF();
      real_t Qp[NUM_STATE] = {0}, Qm[NUM_STATE] = {0};

      setenv("SEAS_TPV104_FREEZE_SIGMA_N", "50e6", 1);
      EvalStageState s;
      ff.ComputeStageState(d, Qp, Qm, s);
      unsetenv("SEAS_TPV104_FREEZE_SIGMA_N");

      TEST_NEAR(s.sigma_n_total, 50.0e6, 1.0,
                "fluctuation-Q: =50e6 pins sigma_n_total to 50 MPa");
   }

   std::cout << "\n=== Summary: " << num_passed << " passed, "
             << num_failed << " failed of " << num_tests << " ===\n";
   return num_failed == 0 ? 0 : 1;
}
