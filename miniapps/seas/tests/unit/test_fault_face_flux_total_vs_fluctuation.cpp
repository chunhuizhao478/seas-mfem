// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 3 (I-06 part A): equivalence between Evaluate (on
// fluctuation Q + pre-stress in DOFData) and EvaluateTotal (on total Q +
// zero pre-stress in DOFData).
//
// Plan §3 Acceptance Criterion:
// "For 1000 random (Q_plus_fluc, Q_minus_fluc, pre-stress, psi, a)
//  fixtures, compute Q_imp_fluc = Evaluate(Q_fluc) and
//  Q_imp_tot  = EvaluateTotal(Q_fluc + pre).  Assert
//  |Q_imp_fluc + pre - Q_imp_tot| <= 1e-10 * |pre| per component.
//  Same for data.V1, V2, slip_rate (bit-identical modulo solver
//  tolerance ~1e-8).  Same for data.{sigma_n_corr, tau_i_corr} — both
//  paths store TOTAL, so they must match to solver tolerance."
//
// Fault-local frame (matches FaultFaceFlux::ComputeTrialTraction):
//   SXX = sigma_nn, SXY = tau_1 (nt1), SXZ = tau_2 (nt2)
//   VX  = v_n,     VY  = v_t1,          VZ  = v_t2
//
// This test is documented as a migration-reference gate only
// (plan Phase 3 AC last bullet): it does not imply that fluctuation mode
// remains a supported TPV102 runtime option after Phase 4.

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../dynamic/wave_state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
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

// TPV102 material
static const real_t RHO    = 2670.0;
static const real_t CP     = 6000.0;
static const real_t CS     = 3464.0;
static const real_t ZP     = RHO * CP;
static const real_t ZS     = RHO * CS;

// ---------------------------------------------------------------------------
// Build a DOFData for the FLUCTUATION path: pre-stress lives in sigma_n0,
// tau_{1,2}_0.  psi, a set per fixture.
// ---------------------------------------------------------------------------
static DOFData MakeFluctuationDOF(real_t sigma_n0, real_t tau1_0, real_t tau2_0,
                                   real_t psi, real_t a)
{
   DOFData d;
   d.Zp_plus = ZP; d.Zp_minus = ZP;
   d.Zs_plus = ZS; d.Zs_minus = ZS;
   d.eta_p = ZP / 2.0;
   d.eta_s = ZS / 2.0;
   d.sigma_n0 = sigma_n0;
   d.tau1_0   = tau1_0;
   d.tau2_0   = tau2_0;
   d.a  = a;
   d.Dc = 0.14;
   d.psi = psi;
   d.slip_rate = 0.0;
   d.V1 = 0; d.V2 = 0;
   d.slip1 = 0; d.slip2 = 0;
   d.tau1_corr = 0; d.tau2_corr = 0; d.sigma_n_corr = 0;
   return d;
}

// ---------------------------------------------------------------------------
// Build a DOFData for the TOTAL path: pre-stress zeroed (driver baked
// pre-stress into bulk Q under the v9.3.0 migration).
// ---------------------------------------------------------------------------
static DOFData MakeTotalDOF(real_t psi, real_t a)
{
   return MakeFluctuationDOF(0.0, 0.0, 0.0, psi, a);
}

// ---------------------------------------------------------------------------
// Run both paths on one random fixture.  Report worst |imp_fluc + pre -
// imp_tot| across all 9 components and both sides, along with the
// absolute differences in DOFData.{V1, V2, slip_rate, sigma_n_corr,
// tau1_corr, tau2_corr}.
// ---------------------------------------------------------------------------
struct OnePointResult
{
   real_t worst_imp_diff = 0.0;   // plus and minus combined
   real_t pre_scale      = 0.0;
   real_t V1_diff = 0, V2_diff = 0, V_diff = 0;
   real_t sigma_n_corr_diff = 0, tau1_corr_diff = 0, tau2_corr_diff = 0;
};

static OnePointResult RunOneFixture(const real_t Q_plus_fluc[NUM_STATE],
                                    const real_t Q_minus_fluc[NUM_STATE],
                                    real_t sigma_n0, real_t tau1_0, real_t tau2_0,
                                    real_t psi, real_t a,
                                    const FaultFaceFlux &flux)
{
   // Pre-stress state in the fault-local frame (symmetric on both sides).
   // Only the traction components (SXX, SXY, SXZ) carry pre-stress.
   real_t pre[NUM_STATE] = {0};
   pre[SXX] = sigma_n0;
   pre[SXY] = tau1_0;
   pre[SXZ] = tau2_0;

   real_t Q_plus_tot[NUM_STATE], Q_minus_tot[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_plus_tot[c]  = Q_plus_fluc[c]  + pre[c];
      Q_minus_tot[c] = Q_minus_fluc[c] + pre[c];
   }

   // Fluctuation path
   DOFData data_fluc = MakeFluctuationDOF(sigma_n0, tau1_0, tau2_0, psi, a);
   real_t Q_imp_plus_fluc[NUM_STATE], Q_imp_minus_fluc[NUM_STATE];
   flux.Evaluate(data_fluc, Q_plus_fluc, Q_minus_fluc,
                 Q_imp_plus_fluc, Q_imp_minus_fluc);

   // Total path (pre-stress zeroed in DOFData)
   DOFData data_tot = MakeTotalDOF(psi, a);
   real_t Q_imp_plus_tot[NUM_STATE], Q_imp_minus_tot[NUM_STATE];
   flux.EvaluateTotal(data_tot, Q_plus_tot, Q_minus_tot,
                      Q_imp_plus_tot, Q_imp_minus_tot);

   OnePointResult r;

   // For stress rows, |(Q_imp_fluc + pre) - Q_imp_tot| should be O(FP noise).
   // For velocity rows, pre-stress has zero velocity so pre contributes 0
   // and |Q_imp_fluc[VX] - Q_imp_tot[VX]| should also be O(FP noise).
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t e_p = std::abs((Q_imp_plus_fluc[c]  + pre[c]) - Q_imp_plus_tot[c]);
      real_t e_m = std::abs((Q_imp_minus_fluc[c] + pre[c]) - Q_imp_minus_tot[c]);
      r.worst_imp_diff = std::max({r.worst_imp_diff, e_p, e_m});
   }
   r.pre_scale = std::max({std::abs(sigma_n0), std::abs(tau1_0),
                            std::abs(tau2_0), real_t(1.0)});

   r.V1_diff = std::abs(data_fluc.V1 - data_tot.V1);
   r.V2_diff = std::abs(data_fluc.V2 - data_tot.V2);
   r.V_diff  = std::abs(data_fluc.slip_rate - data_tot.slip_rate);
   // Both paths store TOTAL in DOFData.{sigma_n_corr, tau_i_corr}, so they
   // should match to solver tolerance.
   r.sigma_n_corr_diff = std::abs(data_fluc.sigma_n_corr - data_tot.sigma_n_corr);
   r.tau1_corr_diff    = std::abs(data_fluc.tau1_corr    - data_tot.tau1_corr);
   r.tau2_corr_diff    = std::abs(data_fluc.tau2_corr    - data_tot.tau2_corr);

   return r;
}

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 3 (I-06): EvaluateTotal vs Evaluate "
                "equivalence (migration reference) ===\n";

   FaultFaceFlux flux(RHO, CP, CS);

   // 1000 random fixtures.  Ranges:
   //   sigma_n0 : 50-150 MPa (positive, compression)
   //   tau_{1,2}_0 : -50 .. +50 MPa
   //   fluctuation stresses : +/- 10 MPa
   //   fluctuation velocities : +/- 1 m/s
   //   psi : 0.1 .. 2.0 (physical range under Dieterich-Ruina)
   //   a   : 0.008 .. 0.025
   std::mt19937_64 rng(0xC0FFEE);
   std::uniform_real_distribution<double> sn_d(50e6, 150e6);
   std::uniform_real_distribution<double> tau_d(-50e6, 50e6);
   std::uniform_real_distribution<double> stress_fluc_d(-1.0e7, 1.0e7);
   std::uniform_real_distribution<double> vel_fluc_d(-1.0, 1.0);
   std::uniform_real_distribution<double> psi_d(0.1, 2.0);
   std::uniform_real_distribution<double> a_d (0.008, 0.025);

   const int N = 1000;
   // Per-component relative tolerance: 1e-10 * |pre| per plan Phase 3 AC.
   // Absolute floor 1.0 Pa to absorb FP noise when a component of pre is 0.
   const real_t kImpRelTol = 1.0e-10;
   const real_t kImpAbsFloor = 1.0;     // Pa, Pa*m/s (velocity imp rows)

   // Solver tolerance: the Brent solver inside FrictionSolver::Solve
   // has default tolerance ~1e-8 in V (relative).
   const real_t kSolverVTol = 1.0e-6;     // m/s absolute
   const real_t kSolverTauTol = 1.0e-3;   // Pa absolute

   real_t worst_imp_ratio = 0.0;
   real_t worst_V        = 0.0;
   real_t worst_V1       = 0.0;
   real_t worst_V2       = 0.0;
   real_t worst_sigma_n_corr = 0.0;
   real_t worst_tau1_corr    = 0.0;
   real_t worst_tau2_corr    = 0.0;

   int passed_per_fixture = 0;

   for (int i = 0; i < N; i++)
   {
      real_t Q_plus_fluc[NUM_STATE] = {0};
      real_t Q_minus_fluc[NUM_STATE] = {0};
      for (int c = 0; c < VX; c++)
      {
         Q_plus_fluc[c]  = stress_fluc_d(rng);
         Q_minus_fluc[c] = stress_fluc_d(rng);
      }
      for (int c = VX; c < NUM_STATE; c++)
      {
         Q_plus_fluc[c]  = vel_fluc_d(rng);
         Q_minus_fluc[c] = vel_fluc_d(rng);
      }
      const real_t sigma_n0 = sn_d(rng);
      const real_t tau1_0   = tau_d(rng);
      const real_t tau2_0   = tau_d(rng);
      const real_t psi      = psi_d(rng);
      const real_t a        = a_d(rng);

      OnePointResult r = RunOneFixture(Q_plus_fluc, Q_minus_fluc,
                                       sigma_n0, tau1_0, tau2_0, psi, a,
                                       flux);
      real_t denom = std::max(r.pre_scale, kImpAbsFloor);
      real_t ratio = r.worst_imp_diff / denom;
      worst_imp_ratio       = std::max(worst_imp_ratio,       ratio);
      worst_V               = std::max(worst_V,               r.V_diff);
      worst_V1              = std::max(worst_V1,              r.V1_diff);
      worst_V2              = std::max(worst_V2,              r.V2_diff);
      worst_sigma_n_corr    = std::max(worst_sigma_n_corr,    r.sigma_n_corr_diff);
      worst_tau1_corr       = std::max(worst_tau1_corr,       r.tau1_corr_diff);
      worst_tau2_corr       = std::max(worst_tau2_corr,       r.tau2_corr_diff);
      if (ratio <= kImpRelTol) { passed_per_fixture++; }
   }

   std::cout << "\n-- summary over " << N << " random fixtures --\n"
             << std::scientific << std::setprecision(3)
             << "  worst |imp_fluc + pre - imp_tot| / max(|pre|, 1) = "
             << worst_imp_ratio << "  (budget " << kImpRelTol << ")\n"
             << "  worst |V_fluc - V_tot|       = " << worst_V       << " m/s\n"
             << "  worst |V1_fluc - V1_tot|     = " << worst_V1      << " m/s\n"
             << "  worst |V2_fluc - V2_tot|     = " << worst_V2      << " m/s\n"
             << "  worst |sigma_n_corr_diff|    = " << worst_sigma_n_corr << " Pa\n"
             << "  worst |tau1_corr_diff|       = " << worst_tau1_corr    << " Pa\n"
             << "  worst |tau2_corr_diff|       = " << worst_tau2_corr    << " Pa\n"
             << "  fixtures meeting imp budget  = " << passed_per_fixture
             << " / " << N << "\n";

   TEST_LE(worst_imp_ratio, kImpRelTol,
           "Phase 3 AC: |imp_fluc + pre - imp_tot| <= 1e-10 * |pre|");
   TEST_LE(worst_V,  kSolverVTol,
           "Phase 3 AC: |V_fluc - V_tot| matches modulo solver tol");
   TEST_LE(worst_V1, kSolverVTol,
           "Phase 3 AC: |V1_fluc - V1_tot| matches modulo solver tol");
   TEST_LE(worst_V2, kSolverVTol,
           "Phase 3 AC: |V2_fluc - V2_tot| matches modulo solver tol");
   TEST_LE(worst_sigma_n_corr, kSolverTauTol,
           "Phase 3 AC: data.sigma_n_corr matches TOTAL on both paths");
   TEST_LE(worst_tau1_corr,    kSolverTauTol,
           "Phase 3 AC: data.tau1_corr matches TOTAL on both paths");
   TEST_LE(worst_tau2_corr,    kSolverTauTol,
           "Phase 3 AC: data.tau2_corr matches TOTAL on both paths");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
