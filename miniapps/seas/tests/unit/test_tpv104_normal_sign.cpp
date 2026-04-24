// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for the TPV104 normal-stress sign convention audit
// (§4.10 Step 8 gates T_TPV104_SIGN_1..4).
//
// Decision (plan §3.11): MFEM uses σ_n > 0 = compression (geology
// convention).  BP5, TPV102, and TPV104 all follow this.  The reference
// FVW runtime uses σ_n < 0 = compression; Step 13 probe-diff tooling
// normalises the sign before comparing.
//
// The tests below guard the invariant at three levels:
//   SIGN_1 — BP5 and TPV104 configuration constants keep σ_n > 0.
//   SIGN_2 — TPV104 DOFData populated from the config keeps σ_n0 > 0.
//   SIGN_3 — the friction-strength formula uses `std::abs(σ_n)`.
//   SIGN_4 — probe-diff tool sign-flip (Step 13 deferral).

#include "test_macros.hpp"
#include "../../config/bp5_params.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../config/tpv104_params.hpp"
#include "../../dynamic/tpv104_friction_solver.hpp"
#include "../../friction/friction_coeff_stable.hpp"
#include "../../dynamic/fault_face_flux.hpp"

#include <cstdlib>
#include <cmath>

using namespace mfem;
using namespace mfem::seas;

// ----------------------------------------------------------------------------
// T_TPV104_SIGN_1 — BP5 convention guard.
// Plan text proposes static_assert but BP5Params::sigma_n is a
// non-constexpr default-initialised member, so only a runtime test is
// possible.  TPV102Params and TPV104Params use constexpr, so we can add
// static_assert for those.
// ----------------------------------------------------------------------------
void TestBP5ConventionGuard()
{
   std::cout << "\n[T_TPV104_SIGN_1] BP5 / TPV102 / TPV104 σ_n > 0 guard\n";

   // Compile-time guards for TPV102 and TPV104 (constexpr).
   static_assert(TPV102Params::sigma_n > 0,
                 "TPV102 sigma_n must be positive (geology convention)");
   static_assert(TPV104Params::sigma_n > 0,
                 "TPV104 sigma_n must be positive (geology convention)");

   // Runtime guard for BP5 (non-constexpr default-initialised member).
   const BP5Params bp5;
   TEST_ASSERT(bp5.sigma_n > 0.0,
               "BP5 sigma_n > 0 (positive compression)");
   TEST_ASSERT(TPV102Params::sigma_n > 0,
               "TPV102 sigma_n > 0");
   TEST_ASSERT(TPV104Params::sigma_n > 0,
               "TPV104 sigma_n > 0");

   // Numerical anchors — the exact values are the SCEC benchmark spec
   // and should not drift.  A test failure here means the config file
   // was edited; the maintainer must consult the plan before landing.
   TEST_NEAR(bp5.sigma_n, 25e6, 1e-3, "BP5 σ_n = 25 MPa");
   TEST_NEAR(TPV102Params::sigma_n, 120e6, 1e-3, "TPV102 σ_n = 120 MPa");
   TEST_NEAR(TPV104Params::sigma_n, 120e6, 1e-3, "TPV104 σ_n = 120 MPa");
}

// ----------------------------------------------------------------------------
// T_TPV104_SIGN_2_storage — DOFData storage-invariant gate (this test).
//
// R3-003 (review round 3): the plan's T_TPV104_SIGN_2 gate, as written,
// requires end-to-end coverage via `InitializeFaultDOFs_TPV104` (Step 3
// / Step 9).  That production-init function has not shipped yet, so we
// can only test the STORAGE invariant here — that a manually populated
// `DOFData.sigma_n0 = TPV104Params::sigma_n` survives as a positive
// value.  This is necessary but NOT sufficient; a Step-9 init bug that
// flipped the sign would not be caught by this test.
//
// The plan's T_TPV104_SIGN_2 gate is hereby split:
//   - T_TPV104_SIGN_2_storage  — closed now (this function).
//   - T_TPV104_SIGN_2_init     — DEFERRED to Step 9.  The end-to-end
//     test lives in `tests/unit/test_tpv104_setup.cpp` (to be created
//     when Step 3 / Step 9 land) and calls `InitializeFaultDOFs_TPV104`
//     then asserts `dof_data[i].sigma_n0 > 0` at every QP.
// ----------------------------------------------------------------------------
void TestSigmaN0PropagationStorage()
{
   std::cout << "\n[T_TPV104_SIGN_2_storage] σ_n0 storage invariant "
             << "(production-init _init gate deferred to Step 9)\n";

   DOFData d;
   d.sigma_n0 = TPV104Params::sigma_n;
   d.tau1_0   = 0.0;
   d.tau2_0   = TPV104Params::tau_ini;

   TEST_ASSERT(d.sigma_n0 > 0.0,
               "DOFData.sigma_n0 positive after TPV104 config");
   TEST_NEAR(d.sigma_n0, 120e6, 1e-3,
             "DOFData.sigma_n0 = 120 MPa");

   // Default nucleation channels: sigma_n_nuc starts at 0.
   TEST_NEAR(d.sigma_n_nuc, 0.0, 1e-30,
             "DOFData.sigma_n_nuc starts at 0 (TPV104 no normal-stress nuc)");

   // Total normal stress at rest = sigma_n0 + sigma_n_nuc = sigma_n0.
   const real_t sigma_n_total = d.sigma_n0 + d.sigma_n_nuc;
   TEST_ASSERT(sigma_n_total > 0.0,
               "sigma_n_total positive throughout nucleation");
}

// ----------------------------------------------------------------------------
// T_TPV104_SIGN_3 — friction strength formula uses |σ_n|.
// The Newton solver's residual is g = -(|σ_n|·μ - τ)/η - V.  Direct
// invocation with a sign-flipped σ_n (physically a tensile fault) must
// NOT produce negative friction — it takes the early-return path.  For
// positive σ_n, we verify the residual matches a hand computation that
// explicitly inserts std::abs().
// ----------------------------------------------------------------------------
void TestStrengthUsesMagnitude()
{
   std::cout << "\n[T_TPV104_SIGN_3] friction strength uses |σ_n|\n";

   const real_t tau      = 40e6;
   const real_t psi      = 0.5636;
   const real_t sigma_n  = 120e6;
   const real_t eta_s    = 4.625e6;
   const real_t a        = 0.01;
   const real_t V0       = 1e-6;

   // Positive σ_n: solver runs Newton and returns a small V.
   int iter = 0; bool conv = false;
   const real_t V_pos = SolveSlipRateNewtonStable(tau, psi, sigma_n, eta_s,
                                                  a, V0, /*V_prev=*/1e-16,
                                                  60, 1e-10, &iter, &conv);
   TEST_ASSERT(conv, "Newton converges at σ_n > 0");
   TEST_ASSERT(V_pos >= 0.0, "V ≥ 0 at positive σ_n");

   // Hand-compute residual at V_pos to confirm sign convention:
   //   g = -(|σ_n|·μ - τ)/η - V_pos
   // |σ_n| should equal σ_n for σ_n > 0, so:
   const real_t mu_at_V = friction_stable::FrictionCoefficientStable(
                             V_pos, psi, a, V0);
   const real_t g_expected = -(std::abs(sigma_n) * mu_at_V - tau) / eta_s
                             - V_pos;
   TEST_ASSERT(std::abs(g_expected) < 1e-9,
               "residual with |σ_n|·μ lands within Newton tol");

   // Tensile σ_n: solver early-returns V = τ/η without invoking μ.
   int iter_ten = 0; bool conv_ten = false;
   const real_t V_ten = SolveSlipRateNewtonStable(tau, psi, -sigma_n, eta_s,
                                                  a, V0, /*V_prev=*/1e-16,
                                                  60, 1e-10,
                                                  &iter_ten, &conv_ten);
   TEST_ASSERT(conv_ten,
               "tensile σ_n reports converged (early return)");
   TEST_NEAR(V_ten, tau / eta_s, 1e-12,
             "tensile σ_n → V = τ/η (frictionless limit)");
   TEST_ASSERT(iter_ten == 0,
               "tensile early-return uses 0 Newton iterations");

   // Sanity: a zero normal stress also takes the tensile branch (σ_n ≤ 0).
   int iter_zero = 0; bool conv_zero = false;
   const real_t V_zero = SolveSlipRateNewtonStable(tau, psi, /*σ_n=*/0.0,
                                                   eta_s, a, V0, 1e-16,
                                                   60, 1e-10,
                                                   &iter_zero, &conv_zero);
   TEST_NEAR(V_zero, tau / eta_s, 1e-12,
             "σ_n = 0 → V = τ/η (frictionless limit)");
}

// ----------------------------------------------------------------------------
// T_TPV104_SIGN_4 — probe-diff normalisation.  Step 13 deferral.
// This test exists to make the gate visible; the actual assertion
// moves to pytest in tpv104/scripts/tests/ once Step 13 lands.
// ----------------------------------------------------------------------------
void TestProbeDiffNormalizationDeferred()
{
   std::cout << "\n[T_TPV104_SIGN_4] probe-diff sign-flip — "
             << "deferred to Step 13 (pytest tpv104/scripts/tests/)\n";
   // No TEST_ASSERT — see R-012 rationale: deferred tests must not
   // inflate the pass count.
}

int main(int argc, char *argv[])
{
   TestBP5ConventionGuard();
   TestSigmaN0PropagationStorage();
   TestStrengthUsesMagnitude();

   // R3-003 deferred placeholder — _init gate is a Step-9 item:
   std::cout << "\n[T_TPV104_SIGN_2_init] production-init sign "
             << "propagation — DEFERRED to Step 9 "
             << "(test_tpv104_setup.cpp)\n";
   TestProbeDiffNormalizationDeferred();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
