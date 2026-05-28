// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_friction_substep_iterator_parity.cpp — Phase 5 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// Acceptance gate: each unified iterator reproduces its standalone oracle
// bit-for-bit (DOFData + I_imp_±) on a fixed multi-QP fixture, for ADER
// orders O ∈ {1, 2, 3}.
//
//   LSW  : LinearSlipWeakeningIterator                vs Tpv205SubStepIterator
//   aging: RateStateSubStepIterator<RateStateAgingPolicy>      vs Tpv102SubStepIterator
//   SRW  : RateStateSubStepIterator<RateStateSlipLawSrwPolicy> vs Tpv104SubStepIterator
//
// Both sides run the per-sub-step pointwise-Q (`AdvanceWithSubStepStates` /
// `Advance`) path on FRESH copies of the same fixture, then compare every
// mutated DOFData field and both I_imp flat arrays with an EXACT (== 0) tol.
//
// Nucleation handling (so the two sides see identical nucleation):
//   - LSW / aging: the standalone callback overload + the unified both get a
//     no-op callback (no nucleation on either side).
//   - SRW: the standalone Tpv104 `AdvanceWithSubStepStates` hard-codes
//     `ApplyNucleationIncremental_TPV104`; the unified is given a callback
//     that calls the same function on its own dof copy, so both apply an
//     identical per-sub-step nucleation increment.

#include "mfem.hpp"

#include "../../dynamic/wave_state.hpp"           // QIndex (VZ)
#include "../../dynamic/fault_face_flux.hpp"      // DOFData, FaultFaceFlux, NUM_STATE
#include "../../dynamic/friction_solver.hpp"      // FrictionSolver::Method
#include "../../dynamic/tpv205_substep_iterator.hpp"
#include "../../dynamic/tpv102_substep_iterator.hpp"
#include "../../dynamic/tpv104_substep_iterator.hpp"
#include "../../dynamic/tpv104_nucleation.hpp"    // ApplyNucleationIncremental_TPV104
#include "../../dynamic/friction_substep_iterator.hpp"
#include "../../friction/state_evolution.hpp"     // AgingLawPsi
#include "../../friction/slip_law_srw_psi.hpp"    // SlipLawSRWPsi

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_EQ0(v, msg) do {                                                  \
   num_tests++;                                                                \
   const double vv = static_cast<double>(v);                                  \
   if (vv == 0.0) { num_passed++;                                             \
      std::cout << "  PASSED: " << msg << "\n"; }                             \
   else { num_failed++;                                                        \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg                    \
                << "  (max|diff| = " << std::scientific                        \
                << std::setprecision(3) << vv << ", expected 0)\n"; }          \
} while (0)

// R-002: guard against vacuous bit-parity (all-zero / unchanged output).
#define TEST_ASSERT(cond, msg) do {                                            \
   num_tests++;                                                                \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; }       \
   else { num_failed++;                                                        \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; }         \
} while (0)

namespace
{

// Representative homogeneous material (granite-ish); the exact values are
// irrelevant to parity — only that both sides see the identical fixture.
constexpr real_t kRho = 2670.0;
constexpr real_t kCp  = 6000.0;
constexpr real_t kCs  = 3464.0;

// Base impedances derived from the material (Z = ρ·c; η = Z/2).
const real_t kZp    = kRho * kCp;
const real_t kZs    = kRho * kCs;
const real_t kEta_p = 0.5 * kZp;
const real_t kEta_s = 0.5 * kZs;

// Fill the impedance/geometry fields shared by every fixture.
void SetImpedances(DOFData &d)
{
   d.Zp_plus  = kZp;  d.Zp_minus = kZp;
   d.Zs_plus  = kZs;  d.Zs_minus = kZs;
   d.eta_p    = kEta_p;
   d.eta_s    = kEta_s;
}

// --- aging (TPV102-style) DOFData; `pert` perturbs the 2nd/3rd QP. ---
DOFData MakeAgingDOF(real_t pert)
{
   DOFData d;
   SetImpedances(d);
   d.sigma_n0 = 50.0e6;
   d.tau1_0   = 0.0;
   d.tau2_0   = 29.2e6 + pert * 1.0e6;
   d.a   = 0.010;
   d.b   = 0.015;                 // Phase 11a per-DOF b (aging policy reads d.b)
   d.Dc  = 2.0;
   d.psi = 0.60 + pert * 0.01;
   return d;
}

// --- SRW (TPV104-style) DOFData. ---
DOFData MakeSrwDOF(real_t pert)
{
   DOFData d;
   SetImpedances(d);
   d.sigma_n0 = 50.0e6;
   d.tau1_0   = 0.0;
   d.tau2_0   = 29.38e6 + pert * 1.0e6;
   d.a   = 0.008 + pert * 0.0001;
   d.Dc  = 0.40;
   d.psi = 0.564 + pert * 0.01;
   return d;
}

// --- LSW (TPV205-style) DOFData. ---
DOFData MakeLswDOF(real_t pert)
{
   DOFData d;
   SetImpedances(d);
   d.sigma_n0  = 120.0e6;
   d.tau1_0    = 0.0;
   // Super-critical: tau2_0 must EXCEED the static strength mu_s*sigma_n =
   // 0.677*120e6 = 81.2 MPa, else the fault stays locked (V2=slip2=0) and the
   // bit-parity is vacuous (R-002).  90 MPa > 81.2 MPa -> ruptures.
   d.tau2_0    = 90.0e6 + pert * 1.0e6;
   d.lsw_mu_s  = 0.677;
   d.lsw_mu_d  = 0.525;
   d.lsw_d_c   = 0.40;
   d.slip1     = 0.0;
   d.slip2     = 0.0;
   // Reset the per-macro-step diag accumulators (the driver does this before
   // each Advance); both sides start identical so the max/min compare equal.
   d.slip_rate_substep_max = 0.0;
   d.sigma_n_substep_min   = std::numeric_limits<real_t>::max();
   return d;
}

// Exact (bit-for-bit) max field diff over every DOFData field the iterators
// mutate, including the LSW-only diag accumulators and the nucleation
// channels (SRW writes them).
real_t MaxDOFFieldDiff(const DOFData &a, const DOFData &b)
{
   real_t m = 0.0;
   m = std::max(m, std::abs(a.psi                   - b.psi));
   m = std::max(m, std::abs(a.slip1                 - b.slip1));
   m = std::max(m, std::abs(a.slip2                 - b.slip2));
   m = std::max(m, std::abs(a.V1                    - b.V1));
   m = std::max(m, std::abs(a.V2                    - b.V2));
   m = std::max(m, std::abs(a.slip_rate             - b.slip_rate));
   m = std::max(m, std::abs(a.tau1_corr             - b.tau1_corr));
   m = std::max(m, std::abs(a.tau2_corr             - b.tau2_corr));
   m = std::max(m, std::abs(a.sigma_n_corr          - b.sigma_n_corr));
   m = std::max(m, std::abs(a.tau1_nuc              - b.tau1_nuc));
   m = std::max(m, std::abs(a.tau2_nuc              - b.tau2_nuc));
   m = std::max(m, std::abs(a.sigma_n_nuc           - b.sigma_n_nuc));
   m = std::max(m, std::abs(a.slip_rate_substep_max - b.slip_rate_substep_max));
   m = std::max(m, std::abs(a.sigma_n_substep_min   - b.sigma_n_substep_min));
   return m;
}

real_t MaxVecDiff(const std::vector<real_t> &a, const std::vector<real_t> &b)
{
   const size_t n = std::min(a.size(), b.size());
   real_t m = (a.size() == b.size()) ? 0.0
              : std::numeric_limits<real_t>::infinity();
   for (size_t i = 0; i < n; ++i)
   { m = std::max(m, std::abs(a[i] - b[i])); }
   return m;
}

real_t MaxDofVecDiff(const std::vector<DOFData> &a,
                     const std::vector<DOFData> &b)
{
   if (a.size() != b.size())
   { return std::numeric_limits<real_t>::infinity(); }
   real_t m = 0.0;
   for (size_t i = 0; i < a.size(); ++i)
   { m = std::max(m, MaxDOFFieldDiff(a[i], b[i])); }
   return m;
}

// R-001: NON-uniform sub-step sizes + weights (still Σδt = dt_macro, Σw = 1)
// so a [0]-vs-[o] mis-index of deltaT_/time_weights_ in RunSubSteps_ is
// observable in the bit-parity comparison.
void MakeQuadrature(int O, real_t dt_macro,
                    std::vector<real_t> &deltaT, std::vector<real_t> &weights)
{
   deltaT.assign(O, 0.0);
   weights.assign(O, 0.0);
   real_t dsum = 0.0, wsum = 0.0;
   for (int o = 0; o < O; ++o)
   {
      deltaT[o]  = 1.0 + 0.5  * o;   // 1.0, 1.5, 2.0, ...
      weights[o] = 1.0 + 0.25 * o;   // 1.0, 1.25, 1.5, ...
      dsum += deltaT[o];
      wsum += weights[o];
   }
   for (int o = 0; o < O; ++o)
   {
      deltaT[o]  = deltaT[o] / dsum * dt_macro;   // Σ = dt_macro
      weights[o] = weights[o] / wsum;             // Σ = 1
   }
}

// R-001: pointwise-Q field that VARIES per sub-step AND per QP, so a
// [0]-vs-[o] mis-index of Q_pointwise in RunSubSteps_ changes the result and
// is caught by the bit-parity comparison (a constant-Q field could not).
void MakeQField(int n, int O, std::vector<std::vector<real_t>> &Qp,
                std::vector<std::vector<real_t>> &Qm)
{
   Qp.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   Qm.assign(O, std::vector<real_t>(static_cast<size_t>(NUM_STATE) * n, 0.0));
   for (int o = 0; o < O; ++o)
   {
      for (int i = 0; i < n; ++i)
      {
         // Along-strike (VZ) velocity jump driving slip; varies with (o, i).
         const real_t s = 1.0e-3 * (1.0 + 0.1 * i) * (1.0 + 0.3 * o);
         Qp[o][static_cast<size_t>(i) * NUM_STATE + VZ] = +s;
         Qm[o][static_cast<size_t>(i) * NUM_STATE + VZ] = -s;
      }
   }
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== Phase 5: unified vs standalone sub-step iterator parity ===\n";

   const int    n        = 3;        // multi-QP fixture
   const real_t dt_macro = 1.0e-4;
   const real_t t0       = 3.7;      // R-001: non-zero so t_macro_start
                                     // propagation through RunSubSteps_ is
                                     // observable via the t-dependent callback.
   const std::vector<int> orders = {1, 2, 3};

   FaultFaceFlux flux(kRho, kCp, kCs);

   std::vector<Vector> coords(n, Vector(3));
   for (int i = 0; i < n; ++i)
   {
      coords[i] = 0.0;
      coords[i](2) = -7500.0;        // near the TPV hypocenter depth
   }

   // R-001: a t-dependent nucleation increment (writes a per-sub-step,
   // time-dependent value into tau2_nuc).  Passed to BOTH sides of the
   // LSW/aging parity so that, if the unified iterator mishandled
   // t_macro_start, the two would diverge.  Each side gets a callback bound
   // to its OWN dof copy (identical operation -> identical result when the
   // t_sub_end values match).
   const auto make_nuc_t = [](std::vector<DOFData> &dof) {
      std::vector<DOFData> *pd = &dof;   // capture the pointer (no dangling ref)
      return [pd](real_t t_sub_end, real_t dt_sub) {
         for (DOFData &d : *pd) { d.tau2_nuc += 1.0e6 * t_sub_end * dt_sub; }
      };
   };

   for (int O : orders)
   {
      std::cout << "\n-- ADER order O = " << O << " --\n";
      std::vector<real_t> deltaT, weights;
      MakeQuadrature(O, dt_macro, deltaT, weights);
      std::vector<std::vector<real_t>> Qp, Qm;
      MakeQField(n, O, Qp, Qm);

      const size_t flat = static_cast<size_t>(NUM_STATE) * n;

      // ============================ LSW ============================
      {
         std::vector<DOFData> dof_std(n), dof_uni(n);
         for (int i = 0; i < n; ++i)
         { dof_std[i] = MakeLswDOF(i); dof_uni[i] = MakeLswDOF(i); }

         std::vector<real_t> Ip_std(flat, 0.0), Im_std(flat, 0.0);
         std::vector<real_t> Ip_uni(flat, 0.0), Im_uni(flat, 0.0);

         Tpv205SubStepIterator it_std(flux);
         it_std.SetSubSteps(deltaT, weights);
         it_std.AdvanceWithSubStepStates(dof_std, coords, Qp, Qm, dt_macro, t0,
                                         Ip_std.data(), Im_std.data(),
                                         make_nuc_t(dof_std));

         LinearSlipWeakeningIterator it_uni(flux);
         it_uni.SetSubSteps(deltaT, weights);
         it_uni.Advance(dof_uni, coords, Qp, Qm, dt_macro, t0,
                        Ip_uni.data(), Im_uni.data(), make_nuc_t(dof_uni));

         TEST_EQ0(MaxDofVecDiff(dof_std, dof_uni), "LSW DOFData bit-parity");
         TEST_EQ0(MaxVecDiff(Ip_std, Ip_uni),      "LSW I_imp_plus bit-parity");
         TEST_EQ0(MaxVecDiff(Im_std, Im_uni),      "LSW I_imp_minus bit-parity");
         // R-002: confirm the friction solve actually moved the state (so the
         // bit-parity above is not vacuously true on an all-zero output).
         TEST_ASSERT(std::abs(dof_std[0].V2) > 0.0
                     || std::abs(dof_std[0].slip2) > 0.0,
                     "LSW fixture exercises a non-trivial friction solve");
      }

      // =========================== aging ===========================
      {
         std::vector<DOFData> dof_std(n), dof_uni(n);
         for (int i = 0; i < n; ++i)
         { dof_std[i] = MakeAgingDOF(i); dof_uni[i] = MakeAgingDOF(i); }

         std::vector<real_t> Ip_std(flat, 0.0), Im_std(flat, 0.0);
         std::vector<real_t> Ip_uni(flat, 0.0), Im_uni(flat, 0.0);

         const real_t b = 0.015, V0 = 1.0e-6, f0 = 0.6;
         const FrictionSolver::Method method = FrictionSolver::Method::Brent;

         AgingLawPsi law_std(b, V0, f0);
         Tpv102SubStepIterator it_std(flux, law_std);
         it_std.SetSubSteps(deltaT, weights);
         it_std.AdvanceWithSubStepStates(dof_std, coords, Qp, Qm, dt_macro, t0,
                                         Ip_std.data(), Im_std.data(),
                                         method, make_nuc_t(dof_std));

         RateStateAgingIterator it_uni(flux, AgingLawPsi(b, V0, f0), method);
         it_uni.SetSubSteps(deltaT, weights);
         it_uni.Advance(dof_uni, coords, Qp, Qm, dt_macro, t0,
                        Ip_uni.data(), Im_uni.data(), make_nuc_t(dof_uni));

         TEST_EQ0(MaxDofVecDiff(dof_std, dof_uni), "aging DOFData bit-parity");
         TEST_EQ0(MaxVecDiff(Ip_std, Ip_uni),      "aging I_imp_plus bit-parity");
         TEST_EQ0(MaxVecDiff(Im_std, Im_uni),      "aging I_imp_minus bit-parity");
         // R-002: non-triviality guard (psi must have moved off the seed, or
         // slip accumulated).
         TEST_ASSERT(dof_std[0].psi != MakeAgingDOF(0).psi
                     || std::abs(dof_std[0].slip2) > 0.0,
                     "aging fixture exercises a non-trivial friction/state solve");
      }

      // ============================ SRW ============================
      {
         std::vector<DOFData> dof_std(n), dof_uni(n);
         for (int i = 0; i < n; ++i)
         { dof_std[i] = MakeSrwDOF(i); dof_uni[i] = MakeSrwDOF(i); }

         std::vector<real_t> Ip_std(flat, 0.0), Im_std(flat, 0.0);
         std::vector<real_t> Ip_uni(flat, 0.0), Im_uni(flat, 0.0);
         // Standalone Tpv104 takes std::vector<real_t>; the unified SRW iterator
         // takes an mfem::Vector (matches the resolver's rs.V_w).  Same values.
         std::vector<real_t> V_w(n, 0.1);
         mfem::Vector V_w_mfem(n); V_w_mfem = 0.1;

         const real_t a0 = 0.008, b = 0.012, V0 = 1.0e-6, f0 = 0.6, muW = 0.1;
         const real_t Vw_default = 0.1;
         const FrictionSolver::Method method =
            FrictionSolver::Method::NewtonRaphsonStable;

         SlipLawSRWPsi law_std(a0, b, V0, f0, muW, Vw_default);
         law_std.SetProductionMode();
         Tpv104SubStepIterator it_std(flux, law_std);
         it_std.SetSubSteps(deltaT, weights);
         it_std.AdvanceWithSubStepStates(dof_std, coords, V_w, Qp, Qm,
                                         dt_macro, t0,
                                         Ip_std.data(), Im_std.data(), method);

         // The unified iterator owns its law by value; build an equivalent one.
         SlipLawSRWPsi law_uni(a0, b, V0, f0, muW, Vw_default);
         law_uni.SetProductionMode();
         RateStateSlipLawSrwIterator it_uni(flux, law_uni, method, &V_w_mfem);
         it_uni.SetSubSteps(deltaT, weights);
         // Replicate Tpv104's hard-coded per-sub-step nucleation on dof_uni.
         const auto nuc_tpv104 = [&](real_t t_sub_end, real_t dt_sub) {
            ApplyNucleationIncremental_TPV104(dof_uni, coords, t_sub_end, dt_sub);
         };
         it_uni.Advance(dof_uni, coords, Qp, Qm, dt_macro, t0,
                        Ip_uni.data(), Im_uni.data(), nuc_tpv104);

         TEST_EQ0(MaxDofVecDiff(dof_std, dof_uni), "SRW DOFData bit-parity");
         TEST_EQ0(MaxVecDiff(Ip_std, Ip_uni),      "SRW I_imp_plus bit-parity");
         TEST_EQ0(MaxVecDiff(Im_std, Im_uni),      "SRW I_imp_minus bit-parity");
         // R-002: non-triviality guard.
         TEST_ASSERT(std::abs(dof_std[0].V2) > 0.0
                     || std::abs(dof_std[0].slip2) > 0.0,
                     "SRW fixture exercises a non-trivial friction/state solve");
      }
   }

   std::cout << "\n========================================\n";
   std::cout << "  Phase 5 substep-iterator parity: " << num_passed
             << " / " << num_tests << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
