// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_tpv102_nuc_callback_parity.cpp — Phase 1 (RS primitives) unit test
// for the new nucleation-callback overload of
// `Tpv102SubStepIterator::AdvanceWithSubStepStates`
// (dynamic/tpv102_substep_iterator.{hpp,cpp}).
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, §9.1
// "Acceptance Criteria" (R-008 — pin one method).
//
// 2-QP, O=2 fixture.  Cases:
//   P1a  the new overload with FrictionSolver::Method::NewtonRaphsonStable
//        (the SAME method the plain overload defaults to) and a callback
//        wrapping ApplyNucleationIncremental_TPV102 reproduces the plain
//        overload BIT-FOR-BIT (psi, slip1/2, V1/2, tau*_corr,
//        sigma_n_corr, tau2_nuc, slip_rate, I_imp±).  Both sides use the
//        same method — a Brent-vs-default mismatch would masquerade as a
//        regression.
//   P1b  the new overload with an EMPTY callback matches a reference run
//        with nucleation suppressed (the plain overload in a t<=0 regime
//        where ApplyNucleationIncremental_TPV102 is a no-op).  A
//        meaningfulness check confirms active nucleation actually changes
//        the solve, so the parity in P1b is not vacuous.
//   P2   the new overload with FrictionSolver::Method::Brent (the
//        production RS method, CLAUDE.md) produces finite, sign-consistent
//        tau*_corr and monotone psi over several macro steps — a smoke
//        that Brent is wired, NOT a bit-for-bit equality.

#include "mfem.hpp"

#include "../../dynamic/tpv102_substep_iterator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_nucleation.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../config/tpv102_params.hpp"        // hypocenter / patch geometry
#include "../../friction/state_evolution.hpp"    // AgingLawPsi

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

// ---------------------------------------------------------------------
// Shared fixture
// ---------------------------------------------------------------------
namespace
{
constexpr int    kNQP      = 2;
constexpr int    kO        = 2;
constexpr real_t kMu       = 32.0e9;
constexpr real_t kRho      = 2670.0;
constexpr real_t kLambda   = 32.0e9;
constexpr real_t kSigmaN0  = 50.0e6;     // effective normal stress [Pa]
constexpr real_t kA        = 0.010;      // RS direct-effect parameter
constexpr real_t kDc       = 0.40;       // critical slip distance [m]
constexpr real_t kB        = 0.015;
constexpr real_t kV0       = 1.0e-6;
constexpr real_t kF0       = 0.6;
constexpr real_t kStrikeShear = 40.0e6;  // Q SXZ -> tau2_trial = kStrikeShear [Pa]

inline real_t Cs() { return std::sqrt(kMu / kRho); }
inline real_t Cp() { return std::sqrt((kLambda + 2.0 * kMu) / kRho); }

// Per-QP material impedances written into DOFData (homogeneous material).
void SetImpedances(DOFData& d)
{
   const real_t cp = Cp(), cs = Cs();
   d.Zp_plus = d.Zp_minus = kRho * cp;
   d.Zs_plus = d.Zs_minus = kRho * cs;
   d.eta_p   = 0.5 * kRho * cp;
   d.eta_s   = 0.5 * kRho * cs;
}

// Build the 2-QP fixture: DOFData, 3D fault coords, and the per-sub-step
// Q± vectors.  QP0 sits at the TPV102 hypocenter (Gaussian factor F = 1);
// QP1 sits 1 km along strike (0 < F < 1) — both inside the nucleation
// patch so active nucleation perturbs them differently.
//
// The Q SXZ component on both sides yields a strike-channel trial
// traction tau2_trial = eta_s · 2·SXZ / Zs = kStrikeShear (eta_s =
// 0.5·rho·cs, Zs = rho·cs); all other Q components are zero, so the
// normal and dip trial tractions are zero (pure strike-slip).
void BuildFixture(std::vector<DOFData>& dof_data,
                  std::vector<Vector>& coords,
                  std::vector<std::vector<real_t>>& Qp,
                  std::vector<std::vector<real_t>>& Qm)
{
   dof_data.assign(kNQP, DOFData{});
   for (auto& d : dof_data)
   {
      SetImpedances(d);
      d.sigma_n0 = kSigmaN0;
      d.tau1_0   = 0.0;
      d.tau2_0   = 0.0;
      d.tau1_nuc = d.tau2_nuc = d.sigma_n_nuc = 0.0;
      d.a   = kA;
      d.b   = kB;   // Phase 11a: iterator reads per-DOF d.b (== AgingLawPsi kB)
      d.Dc  = kDc;
      d.psi = kF0;   // off-equilibrium starting state -> psi relaxes monotonically
   }

   coords.assign(kNQP, Vector(3));
   // QP0: hypocenter (along_strike = 0, down_dip = |z| = hypo_down_dip).
   coords[0](0) = TPV102Params::hypo_along_strike;
   coords[0](1) = 0.0;
   coords[0](2) = -TPV102Params::hypo_down_dip;
   // QP1: 1 km along strike from the hypocenter, same depth.
   coords[1](0) = TPV102Params::hypo_along_strike + 1.0e3;
   coords[1](1) = 0.0;
   coords[1](2) = -TPV102Params::hypo_down_dip;

   const std::size_t words = static_cast<std::size_t>(NUM_STATE) * kNQP;
   Qp.assign(kO, std::vector<real_t>(words, 0.0));
   Qm.assign(kO, std::vector<real_t>(words, 0.0));
   for (int o = 0; o < kO; ++o)
   {
      for (int i = 0; i < kNQP; ++i)
      {
         const std::size_t base = static_cast<std::size_t>(i) * NUM_STATE;
         Qp[o][base + SXZ] = kStrikeShear;   // strike shear on + side
         Qm[o][base + SXZ] = kStrikeShear;   // strike shear on - side
      }
   }
}

// Exact (bit-for-bit) DOFData comparison over the acceptance fields.
bool DofBitEqual(const DOFData& a, const DOFData& b)
{
   return a.psi          == b.psi
       && a.slip1        == b.slip1
       && a.slip2        == b.slip2
       && a.V1           == b.V1
       && a.V2           == b.V2
       && a.tau1_corr    == b.tau1_corr
       && a.tau2_corr    == b.tau2_corr
       && a.sigma_n_corr == b.sigma_n_corr
       && a.slip_rate    == b.slip_rate
       && a.tau2_nuc     == b.tau2_nuc;
}

bool ArraysBitEqual(const std::vector<real_t>& a, const std::vector<real_t>& b)
{
   if (a.size() != b.size()) { return false; }
   for (std::size_t k = 0; k < a.size(); ++k)
   {
      if (a[k] != b[k]) { return false; }
   }
   return true;
}
}  // namespace

// =====================================================================
// P1a: new overload (NewtonRaphsonStable + nuc wrapper) == plain overload
//      bit-for-bit, nucleation ACTIVE.
// =====================================================================
static void P1a_parity_active_nucleation()
{
   std::cout << "\n[P1a] nuc-callback parity vs plain overload (active nucleation)\n";
   FaultFaceFlux flux(kRho, Cp(), Cs());
   AgingLawPsi aging(kB, kV0, kF0);
   Tpv102SubStepIterator iter(flux, aging);

   const real_t dt_macro = 0.2;
   std::vector<real_t> deltaT(kO, dt_macro / kO);
   std::vector<real_t> weights(kO, 1.0 / kO);
   iter.SetSubSteps(deltaT, weights);

   std::vector<DOFData> dof_ref, dof_new;
   std::vector<Vector> coords;
   std::vector<std::vector<real_t>> Qp, Qm;
   BuildFixture(dof_ref, coords, Qp, Qm);
   dof_new = dof_ref;   // identical starting state

   const std::size_t words = static_cast<std::size_t>(NUM_STATE) * kNQP;
   std::vector<real_t> Iip_ref(words), Iim_ref(words);
   std::vector<real_t> Iip_new(words), Iim_new(words);

   // t_macro_start = 0 -> sub-step ends 0.1, 0.2 are inside (0, nuc_T)
   // so the smoothStep increment is positive and nucleation fires.
   const real_t t0 = 0.0;

   // Reference: the existing (plain) overload, which calls
   // ApplyNucleationIncremental_TPV102 internally.
   iter.AdvanceWithSubStepStates(dof_ref, coords, Qp, Qm, dt_macro, t0,
                                 Iip_ref.data(), Iim_ref.data(),
                                 FrictionSolver::Method::NewtonRaphsonStable);

   // New overload: same method, with a callback wrapping the identical
   // nucleation function on the new run's own dof_data.
   auto nuc = [&](real_t t_sub_end, real_t dt_sub)
   {
      ApplyNucleationIncremental_TPV102(dof_new, coords, t_sub_end, dt_sub);
   };
   iter.AdvanceWithSubStepStates(dof_new, coords, Qp, Qm, dt_macro, t0,
                                 Iip_new.data(), Iim_new.data(),
                                 FrictionSolver::Method::NewtonRaphsonStable,
                                 nuc);

   bool dof_eq = true;
   for (int i = 0; i < kNQP; ++i)
   {
      dof_eq = dof_eq && DofBitEqual(dof_ref[i], dof_new[i]);
   }
   TEST_ASSERT(dof_eq, "DOFData bit-identical (plain vs new+wrapper)");
   TEST_ASSERT(ArraysBitEqual(Iip_ref, Iip_new), "I_imp_plus bit-identical");
   TEST_ASSERT(ArraysBitEqual(Iim_ref, Iim_new), "I_imp_minus bit-identical");
   // Confirm nucleation was actually active (so the parity is meaningful).
   TEST_ASSERT(dof_ref[0].tau2_nuc > 0.0,
               "nucleation fired at hypocenter QP (tau2_nuc > 0)");
}

// =====================================================================
// P1b: new overload with EMPTY callback == nucleation-suppressed reference.
// =====================================================================
static void P1b_empty_callback_suppressed()
{
   std::cout << "\n[P1b] empty callback == nucleation-suppressed reference\n";
   FaultFaceFlux flux(kRho, Cp(), Cs());
   AgingLawPsi aging(kB, kV0, kF0);
   Tpv102SubStepIterator iter(flux, aging);

   const real_t dt_macro = 0.2;
   std::vector<real_t> deltaT(kO, dt_macro / kO);
   std::vector<real_t> weights(kO, 1.0 / kO);
   iter.SetSubSteps(deltaT, weights);

   std::vector<DOFData> dof_empty, dof_suppressed, dof_active;
   std::vector<Vector> coords;
   std::vector<std::vector<real_t>> Qp, Qm;
   BuildFixture(dof_empty, coords, Qp, Qm);
   dof_suppressed = dof_empty;
   dof_active     = dof_empty;

   const std::size_t words = static_cast<std::size_t>(NUM_STATE) * kNQP;
   std::vector<real_t> Iip_e(words), Iim_e(words);
   std::vector<real_t> Iip_s(words), Iim_s(words);
   std::vector<real_t> Iip_a(words), Iim_a(words);

   // Empty-callback run: nucleation suppressed by the no-op callback.
   // (The per-QP friction pipeline does not depend on absolute time, so
   //  t_macro_start is immaterial when no nucleation is applied.)
   auto noop = [](real_t, real_t) {};
   iter.AdvanceWithSubStepStates(dof_empty, coords, Qp, Qm, dt_macro,
                                 /*t0=*/0.0, Iip_e.data(), Iim_e.data(),
                                 FrictionSolver::Method::NewtonRaphsonStable,
                                 noop);

   // Reference with nucleation suppressed: the plain overload run in a
   // t <= 0 regime where ApplyNucleationIncremental_TPV102 returns early
   // (smoothStep increment is zero for all sub-steps ending at t <= 0).
   const real_t t_neg = -(dt_macro + 1.0);   // last sub-step ends at < 0
   iter.AdvanceWithSubStepStates(dof_suppressed, coords, Qp, Qm, dt_macro,
                                 t_neg, Iip_s.data(), Iim_s.data(),
                                 FrictionSolver::Method::NewtonRaphsonStable);

   bool dof_eq = true;
   for (int i = 0; i < kNQP; ++i)
   {
      dof_eq = dof_eq && DofBitEqual(dof_empty[i], dof_suppressed[i]);
   }
   TEST_ASSERT(dof_eq, "empty callback == suppressed-nucleation reference (DOFData)");
   TEST_ASSERT(ArraysBitEqual(Iip_e, Iip_s), "I_imp_plus bit-identical");
   TEST_ASSERT(ArraysBitEqual(Iim_e, Iim_s), "I_imp_minus bit-identical");
   TEST_ASSERT(dof_empty[0].tau2_nuc == 0.0,
               "empty callback leaves tau2_nuc at zero");

   // Meaningfulness: the plain overload at t0 = 0 (nucleation ACTIVE)
   // must differ from the suppressed run — otherwise P1b is vacuous.
   iter.AdvanceWithSubStepStates(dof_active, coords, Qp, Qm, dt_macro,
                                 /*t0=*/0.0, Iip_a.data(), Iim_a.data(),
                                 FrictionSolver::Method::NewtonRaphsonStable);
   TEST_ASSERT(dof_active[0].tau2_nuc > 0.0,
               "active run accumulates tau2_nuc (nucleation has an effect)");
   TEST_ASSERT(dof_active[0].tau2_corr != dof_empty[0].tau2_corr,
               "active nucleation changes the solve vs the empty-callback run");
}

// =====================================================================
// P2: Brent smoke — finite, sign-consistent tau*_corr + monotone psi.
// =====================================================================
static void P2_brent_smoke()
{
   std::cout << "\n[P2] Brent smoke (finite, sign-consistent, monotone psi)\n";
   FaultFaceFlux flux(kRho, Cp(), Cs());
   AgingLawPsi aging(kB, kV0, kF0);
   Tpv102SubStepIterator iter(flux, aging);

   const real_t dt_macro = 0.05;
   std::vector<real_t> deltaT(kO, dt_macro / kO);
   std::vector<real_t> weights(kO, 1.0 / kO);
   iter.SetSubSteps(deltaT, weights);

   std::vector<DOFData> dof;
   std::vector<Vector> coords;
   std::vector<std::vector<real_t>> Qp, Qm;
   BuildFixture(dof, coords, Qp, Qm);

   const std::size_t words = static_cast<std::size_t>(NUM_STATE) * kNQP;
   std::vector<real_t> Iip(words), Iim(words);

   const int n_steps = 5;
   std::vector<real_t> psi_seq;
   psi_seq.reserve(n_steps);
   auto noop = [](real_t, real_t) {};   // no nucleation -> constant driving

   bool all_finite = true;
   bool sign_consistent = true;
   real_t t = 0.0;
   for (int k = 0; k < n_steps; ++k)
   {
      iter.AdvanceWithSubStepStates(dof, coords, Qp, Qm, dt_macro, t,
                                    Iip.data(), Iim.data(),
                                    FrictionSolver::Method::Brent, noop);
      t += dt_macro;
      for (int i = 0; i < kNQP; ++i)
      {
         const DOFData& d = dof[i];
         all_finite = all_finite
                      && std::isfinite(d.psi)
                      && std::isfinite(d.tau1_corr)
                      && std::isfinite(d.tau2_corr)
                      && std::isfinite(d.sigma_n_corr)
                      && std::isfinite(d.slip_rate);
         // Driving is positive strike shear -> corrected strike traction
         // stays positive and the fault slips (V_abs >= 0).
         sign_consistent = sign_consistent
                           && (d.tau2_corr > 0.0)
                           && (d.slip_rate >= 0.0);
      }
      psi_seq.push_back(dof[0].psi);
   }

   TEST_ASSERT(all_finite, "all tau*_corr / psi / slip_rate finite under Brent");
   TEST_ASSERT(sign_consistent,
               "tau2_corr > 0 (sign-consistent with positive driving), V_abs >= 0");

   // Monotone psi: the sequence must be non-increasing OR non-decreasing
   // (Brent is wired and the aging-law state relaxes without oscillating).
   bool nonincreasing = true, nondecreasing = true, moved = false;
   for (std::size_t k = 1; k < psi_seq.size(); ++k)
   {
      if (psi_seq[k] > psi_seq[k - 1]) { nonincreasing = false; }
      if (psi_seq[k] < psi_seq[k - 1]) { nondecreasing = false; }
      if (psi_seq[k] != psi_seq[k - 1]) { moved = true; }
   }
   TEST_ASSERT(nonincreasing || nondecreasing, "psi is monotone over the steps");
   TEST_ASSERT(moved, "psi actually evolves (Brent solve drives state)");
}

// =====================================================================
// P3: R-013 — a default-constructed (empty) std::function callback is
//     rejected with std::runtime_error, NOT std::bad_function_call.
// =====================================================================
static void P3_empty_callback_rejected()
{
   std::cout << "\n[P3] R-013 empty callback rejected with runtime_error\n";
   FaultFaceFlux flux(kRho, Cp(), Cs());
   AgingLawPsi aging(kB, kV0, kF0);
   Tpv102SubStepIterator iter(flux, aging);

   const real_t dt_macro = 0.2;
   std::vector<real_t> deltaT(kO, dt_macro / kO);
   std::vector<real_t> weights(kO, 1.0 / kO);
   iter.SetSubSteps(deltaT, weights);

   std::vector<DOFData> dof;
   std::vector<Vector> coords;
   std::vector<std::vector<real_t>> Qp, Qm;
   BuildFixture(dof, coords, Qp, Qm);
   const std::size_t words = static_cast<std::size_t>(NUM_STATE) * kNQP;
   std::vector<real_t> Iip(words), Iim(words);

   std::function<void(real_t, real_t)> empty;   // default-constructed -> empty
   bool threw_runtime = false, threw_bad_call = false;
   try
   {
      iter.AdvanceWithSubStepStates(dof, coords, Qp, Qm, dt_macro, /*t0=*/0.0,
                                    Iip.data(), Iim.data(),
                                    FrictionSolver::Method::Brent, empty);
   }
   catch (const std::bad_function_call&)
   {
      threw_bad_call = true;
   }
   catch (const std::runtime_error&)
   {
      threw_runtime = true;
   }
   TEST_ASSERT(threw_runtime && !threw_bad_call,
               "empty callback rejected with std::runtime_error (not bad_function_call)");
}

// =====================================================================
// P4: R-018 — broadened parity fixture (guards the verbatim-copy
//     divergence risk).  Unlike P1a (single macro-step, O=2, pure
//     strike-slip, zero pre-stress), this drives O=3 sub-steps over
//     SEVERAL macro-steps with BOTH dip and strike pre-stress + trial
//     traction, asserting the new overload (NewtonRaphsonStable + nuc
//     wrapper) stays bit-identical to the plain overload at every step.
// =====================================================================
static void P4_parity_broadened_fixture()
{
   std::cout << "\n[P4] nuc-callback parity, broadened fixture (R-018)\n";
   FaultFaceFlux flux(kRho, Cp(), Cs());
   AgingLawPsi aging(kB, kV0, kF0);
   Tpv102SubStepIterator iter(flux, aging);

   const int    O        = 3;
   const real_t dt_macro = 0.05;
   std::vector<real_t> deltaT(O, dt_macro / O);
   std::vector<real_t> weights(O, 1.0 / O);
   iter.SetSubSteps(deltaT, weights);

   const real_t cs = Cs(), cp = Cp();
   const real_t eta_s = 0.5 * kRho * cs;
   const real_t Zs    = kRho * cs;
   const real_t Sdip    = 8.0e6;    // -> tau1_trial = Sdip   (dip)
   const real_t Sstrike = 30.0e6;   // -> tau2_trial = Sstrike (strike)

   // Two QPs with BOTH dip and strike pre-stress; QP0 at the hypocentre
   // (nucleation fires), QP1 off-centre.
   std::vector<DOFData> dof_ref(kNQP), dof_new(kNQP);
   std::vector<Vector> coords(kNQP, Vector(3));
   for (int i = 0; i < kNQP; ++i)
   {
      DOFData &d = dof_ref[i];
      d.Zp_plus = d.Zp_minus = kRho * cp;
      d.Zs_plus = d.Zs_minus = Zs;
      d.eta_p   = 0.5 * kRho * cp;
      d.eta_s   = eta_s;
      d.sigma_n0 = kSigmaN0;
      d.tau1_0   = 5.0e6;    // dip pre-stress (nonzero)
      d.tau2_0   = 10.0e6;   // strike pre-stress (nonzero)
      d.tau1_nuc = d.tau2_nuc = d.sigma_n_nuc = 0.0;
      d.a   = kA;
      d.b   = kB;   // Phase 11a: iterator reads per-DOF d.b (== AgingLawPsi kB)
      d.Dc  = kDc;
      d.psi = kF0;
   }
   dof_new = dof_ref;
   coords[0](0) = TPV102Params::hypo_along_strike;
   coords[0](1) = 0.0;
   coords[0](2) = -TPV102Params::hypo_down_dip;
   coords[1](0) = TPV102Params::hypo_along_strike + 1.0e3;
   coords[1](1) = 0.0;
   coords[1](2) = -TPV102Params::hypo_down_dip;

   const std::size_t words = static_cast<std::size_t>(NUM_STATE) * kNQP;
   std::vector<std::vector<real_t>> Qp(O, std::vector<real_t>(words, 0.0));
   std::vector<std::vector<real_t>> Qm(O, std::vector<real_t>(words, 0.0));
   for (int o = 0; o < O; ++o)
   {
      for (int i = 0; i < kNQP; ++i)
      {
         const std::size_t base = static_cast<std::size_t>(i) * NUM_STATE;
         Qp[o][base + SXY] = Sdip;     // dip shear
         Qm[o][base + SXY] = Sdip;
         Qp[o][base + SXZ] = Sstrike;  // strike shear
         Qm[o][base + SXZ] = Sstrike;
      }
   }

   std::vector<real_t> Iip_ref(words), Iim_ref(words);
   std::vector<real_t> Iip_new(words), Iim_new(words);
   auto nuc = [&](real_t t_sub_end, real_t dt_sub)
   {
      ApplyNucleationIncremental_TPV102(dof_new, coords, t_sub_end, dt_sub);
   };

   const int n_steps = 3;
   bool all_eq = true;
   for (int k = 0; k < n_steps; ++k)
   {
      const real_t t0 = k * dt_macro;   // inside (0, nuc_T): nucleation active
      iter.AdvanceWithSubStepStates(dof_ref, coords, Qp, Qm, dt_macro, t0,
                                    Iip_ref.data(), Iim_ref.data(),
                                    FrictionSolver::Method::NewtonRaphsonStable);
      iter.AdvanceWithSubStepStates(dof_new, coords, Qp, Qm, dt_macro, t0,
                                    Iip_new.data(), Iim_new.data(),
                                    FrictionSolver::Method::NewtonRaphsonStable,
                                    nuc);
      for (int i = 0; i < kNQP; ++i)
      {
         all_eq = all_eq && DofBitEqual(dof_ref[i], dof_new[i]);
      }
      all_eq = all_eq && ArraysBitEqual(Iip_ref, Iip_new)
                      && ArraysBitEqual(Iim_ref, Iim_new);
   }
   TEST_ASSERT(all_eq,
               "broadened fixture: plain == new+wrapper bit-identical over "
               "3 macro-steps (O=3, dip+strike, nonzero pre-stress)");
   // Meaningful: dip slip accumulated AND nucleation fired.
   TEST_ASSERT(dof_ref[0].slip1 != 0.0,
               "broadened fixture exercises the dip channel (slip1 != 0)");
   TEST_ASSERT(dof_ref[0].tau2_nuc > 0.0,
               "broadened fixture exercises active nucleation (tau2_nuc > 0)");
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 1 test_tpv102_nuc_callback_parity\n";
   P1a_parity_active_nucleation();
   P1b_empty_callback_suppressed();
   P2_brent_smoke();
   P3_empty_callback_rejected();
   P4_parity_broadened_fixture();

   std::cout << "\n========================================\n";
   std::cout << "Phase 1 test_tpv102_nuc_callback_parity: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
