// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 3 (I-06 part A): R-003 regression — EvaluateTotal
// MUST leave data.psi pristine between entry and exit.
//
// Rationale (review-incorporation plan, R-003): the driver's coupled
// RK4-on-psi integrator in drivers/tpv102_driver.cpp passes data.psi to
// Evaluate / EvaluateTotal as the stage-local psi and relies on it being
// UNCHANGED on return.  Any write to data.psi inside the Riemann solver
// double-integrates psi across RK4 stages and silently corrupts the
// dynamics (test_rk4_psi_integration would fail first, but that test
// targets the fluctuation path only — this test closes the gap for the
// new total-stress path).
//
// A bit-exact match is the right check: no arithmetic on psi should
// happen between entry and exit (the solver READS psi via
// FrictionSolver::Solve but does not WRITE through DOFData).

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../dynamic/wave_state.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <random>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

static const real_t RHO = 2670.0;
static const real_t CP  = 6000.0;
static const real_t CS  = 3464.0;
static const real_t ZP  = RHO * CP;
static const real_t ZS  = RHO * CS;

static DOFData MakeTotalDOF(real_t psi, real_t a)
{
   DOFData d;
   d.Zp_plus = ZP; d.Zp_minus = ZP;
   d.Zs_plus = ZS; d.Zs_minus = ZS;
   d.eta_p = ZP / 2.0;
   d.eta_s = ZS / 2.0;
   d.sigma_n0 = 0;  // pre-stress baked into bulk Q under total mode
   d.tau1_0 = 0;
   d.tau2_0 = 0;
   d.a = a;
   d.Dc = 0.14;
   d.psi = psi;
   d.slip_rate = 0.0;
   return d;
}

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 3 (I-06): "
             << "R-003 EvaluateTotal psi-invariance regression ===\n";

   FaultFaceFlux flux(RHO, CP, CS);

   // 200 random total-Q fixtures covering the full (Theta > 0, Theta = 0)
   // dichotomy + extremes of psi and a.
   std::mt19937_64 rng(0xBEEF);
   std::uniform_real_distribution<double> stress_d(-2e8, 2e8);
   std::uniform_real_distribution<double> vel_d   (-5.0, 5.0);
   std::uniform_real_distribution<double> psi_d   (0.05, 3.0);
   std::uniform_real_distribution<double> a_d     (0.008, 0.025);

   int bit_exact_count = 0;
   const int N = 200;
   for (int i = 0; i < N; i++)
   {
      real_t Q_plus[NUM_STATE], Q_minus[NUM_STATE];
      for (int c = 0; c < VX; c++)
      {
         Q_plus[c]  = stress_d(rng);
         Q_minus[c] = stress_d(rng);
      }
      for (int c = VX; c < NUM_STATE; c++)
      {
         Q_plus[c]  = vel_d(rng);
         Q_minus[c] = vel_d(rng);
      }

      const real_t psi_in = psi_d(rng);
      const real_t a      = a_d(rng);
      DOFData data = MakeTotalDOF(psi_in, a);

      real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
      flux.EvaluateTotal(data, Q_plus, Q_minus,
                         Q_imp_plus, Q_imp_minus);

      if (data.psi == psi_in) { bit_exact_count++; }
      else
      {
         std::cout << "  FAILED fixture " << i
                   << ": psi drifted " << std::scientific
                   << std::setprecision(3)
                   << (data.psi - psi_in) << " from " << psi_in << "\n";
      }
   }

   std::cout << "  bit-exact psi-invariance count: " << bit_exact_count
             << " / " << N << "\n";

   TEST_ASSERT(bit_exact_count == N,
               "R-003: EvaluateTotal leaves data.psi bit-exact "
               "across " + std::to_string(N) + " random fixtures");

   // Also test the Theta == 0 shortcut (both sides identical + no
   // velocity jump + no shear → trial traction shear is zero, solver
   // short-circuits to V_abs = 0 and cannot touch psi).
   {
      real_t Q[NUM_STATE] = {0};
      Q[SXX] = 120e6;   // pure normal compression, no shear
      const real_t psi_in = 1.25;
      DOFData data = MakeTotalDOF(psi_in, 0.012);
      real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
      flux.EvaluateTotal(data, Q, Q, Q_imp_plus, Q_imp_minus);
      TEST_ASSERT(data.psi == psi_in,
                  "R-003: Theta == 0 short-circuit still leaves psi pristine");
      TEST_ASSERT(data.slip_rate == 0.0,
                  "R-003: Theta == 0 gives V_abs == 0 (no friction solve)");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
