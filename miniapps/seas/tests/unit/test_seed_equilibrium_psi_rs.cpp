// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_seed_equilibrium_psi_rs.cpp — Phase 1 (RS primitives) unit test
// for `spatial::SeedEquilibriumPsi_RS` (dynamic/spatial_setup.hpp).
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, §9.1
// "Files to Create".  Covers R-001 (the equilibrium seed + the d.eta_s
// vs the nonexistent rs.eta_s member name + the seed<->solve η
// consistency) and R-009 (the V_0 single-source-of-truth guard).
//
// Cases:
//   S1  steady residual after seeding < 1e-6·tau0 (independent
//       cross-check of the friction balance, NOT a reuse of
//       DieterichRuinaFriction's own asinh path).
//   S2  tau0 == 0 -> psi == f0 + b·ln(V0/V_init) (locked branch).
//   S3  V_init <= 0 -> MFEM_VERIFY aborts (fork).
//   S4  R-009: blk.V_0_default = 2e-6 (!= FrictionSolver::V0) aborts;
//       control 1e-6 passes, the seeded fault solves back to V_init,
//       and a 2e-6 seed solved against the solver's 1e-6 is O(1) off
//       (documents the silent-wrong behaviour the guard prevents).

#include "mfem.hpp"

#include "../../dynamic/spatial_setup.hpp"
#include "../../dynamic/friction_solver.hpp"   // FrictionSolver::V0
#include "../../friction/dieterich_ruina.hpp"  // residual cross-check + V solve

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e \
   << ", tol " << _t << ")\n"; num_failed++; } \
   else { std::cout << "  PASSED: " << m << "\n"; num_passed++; } \
   } while (0)

// Run `f` in a forked child; return true iff the child aborted (signal)
// or exited non-zero.  Used to exercise MFEM_VERIFY abort paths without
// killing the test process.  (Mirrors test_spatial_nucleation.cpp.)
static bool RunInChild_(const std::function<void()>& f)
{
   std::fflush(stdout);      // empty the buffer so the child does not re-dup it
   pid_t pid = fork();
   if (pid == 0)
   {
      std::fclose(stderr);   // silence the abort message
      f();
      std::_Exit(0);
   }
   int status = 0;
   waitpid(pid, &status, 0);
   return WIFSIGNALED(status) || (WIFEXITED(status) && WEXITSTATUS(status) != 0);
}

// ---------------------------------------------------------------------
// Shared synthetic fixture (plan §9.1 Acceptance):
//   tau0=29.2e6, σ_n0=49.27e6, V_init=1e-9, Dc=2.0,
//   eta=0.5·sqrt(μρ), b=0.015, V0=1e-6, f0=0.6.
// `a` is not pinned by the acceptance list; use the RateStateBlock
// default a_default = 0.010.
// ---------------------------------------------------------------------
namespace
{
constexpr real_t kMu      = 32.0e9;          // SAFS shear modulus [Pa]
constexpr real_t kRho     = 2670.0;          // SAFS density [kg/m^3]
constexpr real_t kTau0    = 29.2e6;          // strike pre-stress [Pa]
constexpr real_t kSigmaN0 = 49.27e6;         // effective normal stress [Pa]
constexpr real_t kVInit   = 1.0e-9;          // quasi-static initial slip rate [m/s] (production; R-022)
constexpr real_t kDc      = 2.0;             // critical slip distance [m]
constexpr real_t kA       = 0.010;           // direct-effect parameter
constexpr real_t kB       = 0.015;
constexpr real_t kF0      = 0.6;
constexpr real_t kV0      = 1.0e-6;          // == FrictionSolver::V0

// 0.5·sqrt(μρ) — exactly what seed_static_dof_fields writes to d.eta_s
// and FaultFaceFlux uses in the solve (homogeneous SAFS material).
inline real_t EtaS() { return 0.5 * std::sqrt(kMu * kRho); }

// Build a 1-DOF RS param/block pair and the matching DOFData.  `tau2_0`
// and `V_init`/`V0` are knobs the individual cases override.
void MakeFixture(std::vector<DOFData>& dof_data,
                 RateStatePerDOFParams& rs,
                 RateStateBlock& blk,
                 real_t tau2_0, real_t v_init, real_t v0_default)
{
   dof_data.assign(1, DOFData{});
   DOFData& d = dof_data[0];
   d.eta_s     = EtaS();
   d.sigma_n0  = kSigmaN0;
   d.tau1_0    = 0.0;          // pure strike-slip
   d.tau2_0    = tau2_0;
   d.tau1_nuc  = d.tau2_nuc = 0.0;
   d.b         = kB;           // R-028: post-init d.b (= rs.b(0)); seed guards d.b>0
   d.psi       = 0.0;          // the stub the seed overwrites

   rs.a.SetSize(1);      rs.a(0)      = kA;
   rs.b.SetSize(1);      rs.b(0)      = kB;       // unused by the seed
   rs.Dc.SetSize(1);     rs.Dc(0)     = kDc;
   rs.V_init.SetSize(1); rs.V_init(0) = v_init;
   rs.f_0.SetSize(1);    rs.f_0(0)    = kF0;      // unused by the seed
   rs.V_0.SetSize(1);    rs.V_0(0)    = v0_default; // unused by the seed
   rs.eta.SetSize(1);    rs.eta(0)    = EtaS();   // unused by the seed
   rs.sigma_n_eff.SetSize(1); rs.sigma_n_eff(0) = kSigmaN0;

   blk.f_0_default = kF0;
   blk.b_default   = kB;
   blk.V_0_default = v0_default;
}

// Independent friction balance (raw asinh; does NOT call
// DieterichRuinaFriction): σ_n·a·asinh[V/(2V0)·exp(ψ/a)] + η·V.
real_t SteadyTraction(real_t psi, real_t v, real_t sigma_n,
                      real_t eta, real_t a, real_t v0)
{
   const real_t arg = (v / (2.0 * v0)) * std::exp(psi / a);
   return sigma_n * a * std::asinh(arg) + eta * v;
}
}  // namespace

// =====================================================================
// S1: steady residual after seeding < 1e-6·tau0
// =====================================================================
static void S1_steady_residual()
{
   std::cout << "\n[S1] SeedEquilibriumPsi_RS steady residual\n";
   std::vector<DOFData> dof_data;
   RateStatePerDOFParams rs;
   RateStateBlock blk;
   MakeFixture(dof_data, rs, blk, /*tau2_0=*/kTau0, kVInit, kV0);

   SeedEquilibriumPsi_RS(dof_data, rs, blk);

   const real_t psi = dof_data[0].psi;
   TEST_ASSERT(std::isfinite(psi), "seeded psi is finite");

   const real_t pred = SteadyTraction(psi, kVInit, kSigmaN0,
                                       EtaS(), kA, kV0);
   const real_t residual = std::abs(kTau0 - pred);
   TEST_ASSERT(residual < 1e-6 * kTau0,
               "steady residual |tau0 - (sigma_n*f + eta*V_init)| < 1e-6*tau0");
}

// =====================================================================
// S2: tau0 == 0 -> psi == f0 + b·ln(V0/V_init) (locked branch)
// =====================================================================
static void S2_locked_branch()
{
   std::cout << "\n[S2] SeedEquilibriumPsi_RS tau0==0 locked branch\n";
   std::vector<DOFData> dof_data;
   RateStatePerDOFParams rs;
   RateStateBlock blk;
   MakeFixture(dof_data, rs, blk, /*tau2_0=*/0.0, kVInit, kV0);

   SeedEquilibriumPsi_RS(dof_data, rs, blk);

   const real_t expected = kF0 + kB * std::log(kV0 / kVInit);
   TEST_NEAR(dof_data[0].psi, expected, 1e-12,
             "tau0==0 seeds locked steady state f0 + b·ln(V0/V_init)");
}

// =====================================================================
// S3: V_init <= 0 -> MFEM_VERIFY aborts
// =====================================================================
static void S3_vinit_nonpositive_aborts()
{
   std::cout << "\n[S3] SeedEquilibriumPsi_RS V_init <= 0 aborts\n";
   const bool zero_aborts = RunInChild_([]() {
      std::vector<DOFData> dof_data;
      RateStatePerDOFParams rs;
      RateStateBlock blk;
      MakeFixture(dof_data, rs, blk, kTau0, /*v_init=*/0.0, kV0);
      SeedEquilibriumPsi_RS(dof_data, rs, blk);
   });
   TEST_ASSERT(zero_aborts, "V_init == 0 aborts (ill-posed steady state)");

   const bool neg_aborts = RunInChild_([]() {
      std::vector<DOFData> dof_data;
      RateStatePerDOFParams rs;
      RateStateBlock blk;
      MakeFixture(dof_data, rs, blk, kTau0, /*v_init=*/-1.0, kV0);
      SeedEquilibriumPsi_RS(dof_data, rs, blk);
   });
   TEST_ASSERT(neg_aborts, "V_init < 0 aborts (ill-posed steady state)");
}

// =====================================================================
// S4: R-009 V_0 single source of truth
// =====================================================================
static void S4_v0_guard()
{
   std::cout << "\n[S4] R-009 V_0 single-source-of-truth guard\n";

   // (a) blk.V_0_default = 2e-6 != FrictionSolver::V0 (1e-6) -> abort.
   const bool mismatch_aborts = RunInChild_([]() {
      std::vector<DOFData> dof_data;
      RateStatePerDOFParams rs;
      RateStateBlock blk;
      MakeFixture(dof_data, rs, blk, kTau0, kVInit, /*v0_default=*/2.0e-6);
      SeedEquilibriumPsi_RS(dof_data, rs, blk);
   });
   TEST_ASSERT(mismatch_aborts,
               "blk.V_0_default = 2e-6 (!= FrictionSolver::V0) aborts");

   // Precondition: the deliverable's default V_0 equals the solver's.
   TEST_NEAR(kV0, FrictionSolver::V0, 0.0,
             "fixture V_0 == FrictionSolver::V0 (control passes the guard)");

   // (b) Control 1e-6 passes; the seeded fault solves back to V_init.
   std::vector<DOFData> dof_data;
   RateStatePerDOFParams rs;
   RateStateBlock blk;
   MakeFixture(dof_data, rs, blk, kTau0, kVInit, /*v0_default=*/kV0);
   SeedEquilibriumPsi_RS(dof_data, rs, blk);
   const real_t psi_good = dof_data[0].psi;
   TEST_ASSERT(std::isfinite(psi_good), "control seed produces finite psi");

   // The force solve hard-codes FrictionSolver::V0 == 1e-6 (== kV0).
   const DieterichRuinaFriction solver(
      DieterichRuinaFriction::Constants{FrictionSolver::V0, kF0, kB, kDc});
   const real_t v_good = solver.SolveSlipRatePsi(kTau0, psi_good, kSigmaN0,
                                                 EtaS(), kA);
   TEST_ASSERT(std::abs(v_good - kVInit) <= 1e-6 * kVInit,
               "control: t=0 solve sits at V_init (|V - V_init| <= 1e-6*V_init)");

   // (c) Document the silent-wrong path the guard prevents: a ψ SEEDED
   // with V0 = 2e-6 but SOLVED against the solver's hard-coded 1e-6 is
   // O(1) off in V — never reachable in production because the guard
   // aborts, but proven here so the guard's value is explicit.
   const DieterichRuinaFriction seed_2e6(
      DieterichRuinaFriction::Constants{2.0e-6, kF0, kB, kDc});
   const real_t psi_bad = seed_2e6.InitialStatePsi(kTau0, kVInit, kSigmaN0,
                                                   EtaS(), kA);
   const real_t v_bad = solver.SolveSlipRatePsi(kTau0, psi_bad, kSigmaN0,
                                                EtaS(), kA);
   TEST_ASSERT(std::abs(v_bad - kVInit) > 0.1 * kVInit,
               "mismatched V_0 (2e-6 seed vs 1e-6 solve) is O(1) off in V");
}

// =====================================================================
// S5: R-012 — rs vectors shorter than dof_data abort (size validation)
// =====================================================================
static void S5_rs_size_mismatch_aborts()
{
   std::cout << "\n[S5] R-012 rs-vector size validation\n";

   // dof_data.size() = 2 but rs vectors size 1 -> abort before the loop.
   const bool aborts = RunInChild_([]() {
      std::vector<DOFData> dof_data(2, DOFData{});
      for (auto& d : dof_data)
      {
         d.eta_s = EtaS();
         d.sigma_n0 = kSigmaN0;
         d.tau2_0 = kTau0;
      }
      RateStatePerDOFParams rs;
      rs.a.SetSize(1);      rs.a(0)      = kA;
      rs.b.SetSize(1);      rs.b(0)      = kB;       // Phase 11a: seed reads rs.b(ii)
      rs.Dc.SetSize(1);     rs.Dc(0)     = kDc;
      rs.V_init.SetSize(1); rs.V_init(0) = kVInit;
      RateStateBlock blk;
      blk.f_0_default = kF0;
      blk.b_default   = kB;
      blk.V_0_default = kV0;
      SeedEquilibriumPsi_RS(dof_data, rs, blk);
   });
   TEST_ASSERT(aborts, "rs vectors shorter than dof_data aborts (R-012)");

   // Control: matched sizes (2 DOFs, rs size 2) -> no abort, finite psi.
   std::vector<DOFData> dof_data(2, DOFData{});
   for (auto& d : dof_data)
   {
      d.eta_s = EtaS();
      d.sigma_n0 = kSigmaN0;
      d.tau2_0 = kTau0;
      d.b = kB;   // R-028: post-init d.b; seed guards d.b>0
   }
   RateStatePerDOFParams rs;
   rs.a.SetSize(2);
   rs.b.SetSize(2);
   rs.Dc.SetSize(2);
   rs.V_init.SetSize(2);
   for (int i = 0; i < 2; ++i)
   {
      rs.a(i)      = kA;
      rs.b(i)      = kB;
      rs.Dc(i)     = kDc;
      rs.V_init(i) = kVInit;
   }
   RateStateBlock blk;
   blk.f_0_default = kF0;
   blk.b_default   = kB;
   blk.V_0_default = kV0;
   SeedEquilibriumPsi_RS(dof_data, rs, blk);
   TEST_ASSERT(std::isfinite(dof_data[0].psi) && std::isfinite(dof_data[1].psi),
               "matched sizes -> finite psi for all DOFs (R-012 control)");
}

// =====================================================================
// S6: R-014 — radiation-damping-dominated DOF (tau_eff <= 0) seeds the
//     locked steady state, not NaN.
// =====================================================================
static void S6_radiation_damping_dominated()
{
   std::cout << "\n[S6] R-014 radiation-damping-dominated seed (no NaN)\n";
   std::vector<DOFData> dof_data;
   RateStatePerDOFParams rs;
   RateStateBlock blk;
   // tau0 = 1 Pa, V_init = 1 m/s -> eta*V_init (~4.6e6 Pa) >> tau0,
   // i.e. tau_eff < 0; pre-fix InitialStatePsi would log(negative) -> NaN.
   MakeFixture(dof_data, rs, blk, /*tau2_0=*/1.0, /*v_init=*/1.0, kV0);

   SeedEquilibriumPsi_RS(dof_data, rs, blk);

   TEST_ASSERT(std::isfinite(dof_data[0].psi),
               "tau_eff <= 0 DOF seeds finite psi (R-014, no NaN)");
   const real_t expected = kF0 + kB * std::log(kV0 / 1.0);
   TEST_NEAR(dof_data[0].psi, expected, 1e-12,
             "tau_eff <= 0 falls back to locked steady state (R-014)");
}

// =====================================================================
// S7: R-017 — non-positive effective normal stress aborts (InitialStatePsi
//     would otherwise seed inf/NaN psi).
// =====================================================================
static void S7_sigma_n0_nonpositive_aborts()
{
   std::cout << "\n[S7] R-017 sigma_n0 <= 0 aborts\n";

   auto run_with_sigma_n0 = [](real_t sigma_n0) {
      std::vector<DOFData> dof_data(1, DOFData{});
      DOFData& d = dof_data[0];
      d.eta_s    = EtaS();
      d.sigma_n0 = sigma_n0;
      d.tau2_0   = kTau0;
      d.b        = kB;   // R-028: post-init d.b so the sigma_n0 guard is what fires
      RateStatePerDOFParams rs;
      rs.a.SetSize(1);      rs.a(0)      = kA;
      rs.b.SetSize(1);      rs.b(0)      = kB;       // Phase 11a: seed reads rs.b(ii)
      rs.Dc.SetSize(1);     rs.Dc(0)     = kDc;
      rs.V_init.SetSize(1); rs.V_init(0) = kVInit;
      RateStateBlock blk;
      blk.f_0_default = kF0;
      blk.b_default   = kB;
      blk.V_0_default = kV0;
      SeedEquilibriumPsi_RS(dof_data, rs, blk);
   };

   const bool zero_aborts = RunInChild_([&]() { run_with_sigma_n0(0.0); });
   TEST_ASSERT(zero_aborts, "sigma_n0 == 0 aborts (ill-posed; inf psi)");

   const bool neg_aborts = RunInChild_([&]() { run_with_sigma_n0(-1.0e6); });
   TEST_ASSERT(neg_aborts, "sigma_n0 < 0 aborts (ill-posed; NaN psi)");
}

// =====================================================================
// S8: Phase 11a — the seed reads per-DOF rs.b(ii), NOT scalar blk.b_default.
//     Two locked DOFs with distinct b must seed psi = f0 + b_i*ln(V0/V_init)
//     with their OWN b_i; blk.b_default is a decoy that would blow up psi
//     identically on both DOFs if the seed wrongly used it.
// =====================================================================
static void S8_per_dof_b_used_by_seed()
{
   std::cout << "\n[S8] Phase 11a seed uses per-DOF rs.b(ii)\n";
   std::vector<DOFData> dof_data(2, DOFData{});
   for (auto& d : dof_data)
   {
      d.eta_s    = EtaS();
      d.sigma_n0 = kSigmaN0;
      d.tau2_0   = 0.0;   // locked branch -> psi = f0 + b*ln(V0/V_init)
   }
   RateStatePerDOFParams rs;
   rs.a.SetSize(2); rs.b.SetSize(2); rs.Dc.SetSize(2); rs.V_init.SetSize(2);
   const real_t b0 = 0.012, b1 = 0.020;   // distinct per-DOF b
   rs.a(0) = kA;  rs.b(0) = b0;  rs.Dc(0) = kDc;  rs.V_init(0) = kVInit;
   rs.a(1) = kA;  rs.b(1) = b1;  rs.Dc(1) = kDc;  rs.V_init(1) = kVInit;
   dof_data[0].b = b0;  dof_data[1].b = b1;   // R-028: post-init d.b (= rs.b(i)); seed guards d.b>0
   RateStateBlock blk;
   blk.f_0_default = kF0;
   blk.b_default   = 999.0;   // DECOY: used => both psi blow up identically
   blk.V_0_default = kV0;     // R-009 guard: must equal FrictionSolver::V0
   SeedEquilibriumPsi_RS(dof_data, rs, blk);

   const real_t exp0 = kF0 + b0 * std::log(kV0 / kVInit);
   const real_t exp1 = kF0 + b1 * std::log(kV0 / kVInit);
   TEST_NEAR(dof_data[0].psi, exp0, 1e-12,
             "DOF0 seeded with rs.b(0), not blk.b_default (Phase 11a)");
   TEST_NEAR(dof_data[1].psi, exp1, 1e-12,
             "DOF1 seeded with rs.b(1), not blk.b_default (Phase 11a)");
   TEST_ASSERT(std::abs(dof_data[0].psi - dof_data[1].psi) > 1e-6,
               "distinct per-DOF b -> distinct seeded psi");
}

// =====================================================================
// S9: R-028 — SeedEquilibriumPsi_RS aborts when DOFData.b is left at the
//     NaN default, so a forgotten per-DOF b is LOUD, not silent-wrong
//     (b=0 would snap psi to f0 for psi<f0 with no NaN).
// =====================================================================
static void S9_unset_b_aborts()
{
   std::cout << "\n[S9] R-028 unset DOFData.b (NaN default) aborts in the seed\n";

   // DOFData.b default-constructs as NaN (NOT 0.0).
   DOFData fresh{};
   TEST_ASSERT(std::isnan(fresh.b),
               "DOFData.b default-constructs as NaN (R-028, not 0.0)");

   const bool aborts = RunInChild_([]() {
      std::vector<DOFData> dof_data(1, DOFData{});   // d.b left at NaN default
      DOFData& d = dof_data[0];
      d.eta_s    = EtaS();
      d.sigma_n0 = kSigmaN0;
      d.tau2_0   = kTau0;
      // intentionally do NOT set d.b
      RateStatePerDOFParams rs;
      rs.a.SetSize(1);      rs.a(0)      = kA;
      rs.b.SetSize(1);      rs.b(0)      = kB;
      rs.Dc.SetSize(1);     rs.Dc(0)     = kDc;
      rs.V_init.SetSize(1); rs.V_init(0) = kVInit;
      RateStateBlock blk;
      blk.f_0_default = kF0;
      blk.b_default   = kB;
      blk.V_0_default = kV0;
      SeedEquilibriumPsi_RS(dof_data, rs, blk);
   });
   TEST_ASSERT(aborts, "unset d.b (NaN) aborts the seed (R-028)");
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 1 test_seed_equilibrium_psi_rs\n";
   S1_steady_residual();
   S2_locked_branch();
   S3_vinit_nonpositive_aborts();
   S4_v0_guard();
   S5_rs_size_mismatch_aborts();
   S6_radiation_damping_dominated();
   S7_sigma_n0_nonpositive_aborts();
   S8_per_dof_b_used_by_seed();
   S9_unset_b_aborts();

   std::cout << "\n========================================\n";
   std::cout << "Phase 1 test_seed_equilibrium_psi_rs: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
