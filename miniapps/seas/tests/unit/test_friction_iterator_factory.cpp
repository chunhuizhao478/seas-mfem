// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_friction_iterator_factory.cpp — Phase 2 unit test for
// MakeFrictionIterator (dynamic/friction_iterator_factory.{hpp,cpp}).
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, §9.2
// "Acceptance Criteria".
//
//   F1  cfg.law == SlipWeakening -> LswFrictionIterator (WaveOpLaw()==LSW)
//   F2  cfg.law == RateState (aging) -> RateStateAgingFrictionIterator
//       (WaveOpLaw()==RateAndState)
//   F3  cfg.law == RateState with no [friction.rate_state] block -> abort
//   F4  R-009: cfg.rate_state.V_0_default != FrictionSolver::V0 -> abort;
//       control V_0 == FrictionSolver::V0 passes.
//
// Note: the slip-law-SRW abort path is unreachable in Phase 2 — there is
// no state_evolution field in the config until Phase 6, so RateState
// always selects the aging adapter.  The SRW guard placement is marked in
// friction_iterator_factory.cpp; its test arrives with the field.

#include "mfem.hpp"

#include "../../dynamic/friction_iterator_factory.hpp"
#include "../../dynamic/friction_iterator.hpp"     // FaultFrictionLaw, adapters
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_solver.hpp"       // FrictionSolver::V0
#include "../../spatial/code/spatial_friction.hpp"

#include <cstdio>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <type_traits>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

// R-015: the RS adapter binds it_.state_evo_ to its own law_ member, so a
// copy/move would dangle.  It must be non-copyable and non-movable (held
// only via unique_ptr).  Compile-time guard.
static_assert(!std::is_copy_constructible<RateStateAgingFrictionIterator>::value,
              "R-015: RateStateAgingFrictionIterator must not be copy-constructible");
static_assert(!std::is_move_constructible<RateStateAgingFrictionIterator>::value,
              "R-015: RateStateAgingFrictionIterator must not be move-constructible");
static_assert(!std::is_copy_assignable<RateStateAgingFrictionIterator>::value,
              "R-015: RateStateAgingFrictionIterator must not be copy-assignable");
static_assert(!std::is_move_assignable<RateStateAgingFrictionIterator>::value,
              "R-015: RateStateAgingFrictionIterator must not be move-assignable");

// The SRW adapter binds it_.state_evo_ to its own state_evo_ member; same
// dangling-ref footgun as the aging adapter -> must be non-copyable/movable.
static_assert(!std::is_copy_constructible<SlipLawSRWFrictionIterator>::value,
              "SlipLawSRWFrictionIterator must not be copy-constructible");
static_assert(!std::is_move_constructible<SlipLawSRWFrictionIterator>::value,
              "SlipLawSRWFrictionIterator must not be move-constructible");
static_assert(!std::is_copy_assignable<SlipLawSRWFrictionIterator>::value,
              "SlipLawSRWFrictionIterator must not be copy-assignable");
static_assert(!std::is_move_assignable<SlipLawSRWFrictionIterator>::value,
              "SlipLawSRWFrictionIterator must not be move-assignable");

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

static bool RunInChild_(const std::function<void()>& f)
{
   std::fflush(stdout);
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

namespace
{
constexpr real_t kRho = 2670.0, kCp = 6000.0, kCs = 3500.0;

spatial::SpatialFrictionConfig MakeLswConfig()
{
   spatial::SpatialFrictionConfig cfg;
   cfg.law = spatial::FrictionLawKind::SlipWeakening;
   cfg.slip_weakening = spatial::SlipWeakeningBlock{};   // defaults
   return cfg;
}

spatial::SpatialFrictionConfig MakeRateStateConfig(real_t v0)
{
   spatial::SpatialFrictionConfig cfg;
   cfg.law = spatial::FrictionLawKind::RateState;
   spatial::RateStateBlock blk;   // defaults: V_0_default = 1e-6
   blk.V_0_default = v0;
   cfg.rate_state = blk;
   return cfg;
}

// slip-law-SRW RateState config (state_evolution = SlipSRW + f_w/V_w globals).
spatial::SpatialFrictionConfig MakeSlipSRWConfig(real_t v0)
{
   spatial::SpatialFrictionConfig cfg;
   cfg.law = spatial::FrictionLawKind::RateState;
   spatial::RateStateBlock blk;
   blk.V_0_default     = v0;
   blk.state_evolution = spatial::StateEvolutionKind::SlipSRW;
   blk.f_w_default     = 0.2;
   blk.V_w_default     = 0.1;
   cfg.rate_state = blk;
   return cfg;
}

// Minimal resolved per-DOF params for the SRW factory branch.  The factory
// requires rs != nullptr and rs->V_w.Size() == rs->a.Size() > 0; `fill_vw=false`
// leaves V_w empty to exercise the size-mismatch guard.
spatial::RateStatePerDOFParams MakeRsWithVw(int n, bool fill_vw)
{
   spatial::RateStatePerDOFParams rs;
   rs.a.SetSize(n);           rs.a = 0.0127;
   rs.b.SetSize(n);           rs.b = 0.0261;
   rs.Dc.SetSize(n);          rs.Dc = 0.10;
   rs.V_init.SetSize(n);      rs.V_init = 1.0e-12;
   rs.f_0.SetSize(n);         rs.f_0 = 0.6;
   rs.V_0.SetSize(n);         rs.V_0 = 1.0e-6;
   rs.eta.SetSize(n);         rs.eta = 1.0e6;
   rs.sigma_n_eff.SetSize(n); rs.sigma_n_eff = 48.9e6;
   if (fill_vw) { rs.V_w.SetSize(n); rs.V_w = 0.1; }
   return rs;
}
}  // namespace

// =====================================================================
// F1: SlipWeakening -> LswFrictionIterator (WaveOpLaw()==LSW)
// =====================================================================
static void F1_slip_weakening()
{
   std::cout << "\n[F1] SlipWeakening -> LswFrictionIterator\n";
   FaultFaceFlux flux(kRho, kCp, kCs);
   const spatial::SpatialFrictionConfig cfg = MakeLswConfig();
   std::unique_ptr<IFrictionIterator> fr =
      MakeFrictionIterator(cfg, flux, /*rs=*/nullptr);
   TEST_ASSERT(fr != nullptr, "factory returns a non-null iterator");
   TEST_ASSERT(fr->WaveOpLaw() == FaultFrictionLaw::LSW,
               "SlipWeakening iterator WaveOpLaw() == LSW");
   TEST_ASSERT(dynamic_cast<LswFrictionIterator*>(fr.get()) != nullptr,
               "SlipWeakening iterator is a LswFrictionIterator");
}

// =====================================================================
// F2: RateState (aging) -> RateStateAgingFrictionIterator
// =====================================================================
static void F2_rate_state_aging()
{
   std::cout << "\n[F2] RateState -> RateStateAgingFrictionIterator\n";
   FaultFaceFlux flux(kRho, kCp, kCs);
   const spatial::SpatialFrictionConfig cfg =
      MakeRateStateConfig(FrictionSolver::V0);   // consistent V_0
   std::unique_ptr<IFrictionIterator> fr =
      MakeFrictionIterator(cfg, flux, /*rs=*/nullptr);
   TEST_ASSERT(fr != nullptr, "factory returns a non-null iterator");
   TEST_ASSERT(fr->WaveOpLaw() == FaultFrictionLaw::RateAndState,
               "RateState iterator WaveOpLaw() == RateAndState");
   TEST_ASSERT(dynamic_cast<RateStateAgingFrictionIterator*>(fr.get()) != nullptr,
               "RateState iterator is a RateStateAgingFrictionIterator");
}

// =====================================================================
// F3: RateState with no [friction.rate_state] block -> abort
// =====================================================================
static void F3_rate_state_missing_block_aborts()
{
   std::cout << "\n[F3] RateState without rate_state block aborts\n";
   const bool aborts = RunInChild_([]() {
      FaultFaceFlux flux(kRho, kCp, kCs);
      spatial::SpatialFrictionConfig cfg;
      cfg.law = spatial::FrictionLawKind::RateState;   // cfg.rate_state == nullopt
      auto fr = MakeFrictionIterator(cfg, flux, nullptr);
      (void) fr;
   });
   TEST_ASSERT(aborts, "RateState with no [friction.rate_state] block aborts");
}

// =====================================================================
// F4: R-009 — V_0 mismatch aborts; control (V_0 == FrictionSolver::V0) passes
// =====================================================================
static void F4_v0_guard()
{
   std::cout << "\n[F4] R-009 V_0 guard in the factory\n";
   const bool mismatch_aborts = RunInChild_([]() {
      FaultFaceFlux flux(kRho, kCp, kCs);
      const spatial::SpatialFrictionConfig cfg =
         MakeRateStateConfig(2.0e-6);   // != FrictionSolver::V0 (1e-6)
      auto fr = MakeFrictionIterator(cfg, flux, nullptr);
      (void) fr;
   });
   TEST_ASSERT(mismatch_aborts,
               "RateState V_0 != FrictionSolver::V0 aborts (R-009)");

   // Control: matching V_0 builds the iterator without aborting.
   FaultFaceFlux flux(kRho, kCp, kCs);
   const spatial::SpatialFrictionConfig cfg =
      MakeRateStateConfig(FrictionSolver::V0);
   auto fr = MakeFrictionIterator(cfg, flux, nullptr);
   TEST_ASSERT(fr != nullptr && fr->WaveOpLaw() == FaultFrictionLaw::RateAndState,
               "control: matching V_0 builds the RateState iterator (R-009)");
}

// =====================================================================
// F5: R-016 — LswFrictionIterator::Advance rejects an empty callback with
//     std::runtime_error (NOT std::bad_function_call) on the WIRED path.
// =====================================================================
static void F5_lsw_empty_callback_rejected()
{
   std::cout << "\n[F5] R-016 LswFrictionIterator empty callback rejected\n";
   FaultFaceFlux flux(kRho, kCp, kCs);
   const spatial::SpatialFrictionConfig cfg = MakeLswConfig();
   std::unique_ptr<IFrictionIterator> fr =
      MakeFrictionIterator(cfg, flux, /*rs=*/nullptr);

   std::vector<DOFData> dof;
   std::vector<Vector> coords;
   std::vector<std::vector<real_t>> Qp, Qm;
   std::function<void(real_t, real_t)> empty;   // default-constructed -> empty

   bool threw_runtime = false, threw_bad_call = false;
   try
   {
      // The adapter guard fires before forwarding, so empty fixtures are fine.
      fr->Advance(dof, coords, Qp, Qm, /*dt=*/1.0, /*t0=*/0.0,
                  /*Iimp_p=*/nullptr, /*Iimp_m=*/nullptr, empty);
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
               "LSW Advance rejects empty callback with runtime_error (R-016)");
}

// =====================================================================
// F6: RateState slip_srw -> SlipLawSRWFrictionIterator (WaveOpLaw==RateAndState)
// =====================================================================
static void F6_slip_srw_builds()
{
   std::cout << "\n[F6] RateState slip_srw -> SlipLawSRWFrictionIterator\n";
   FaultFaceFlux flux(kRho, kCp, kCs);
   const spatial::SpatialFrictionConfig cfg = MakeSlipSRWConfig(FrictionSolver::V0);
   const spatial::RateStatePerDOFParams rs = MakeRsWithVw(2, /*fill_vw=*/true);
   std::unique_ptr<IFrictionIterator> fr = MakeFrictionIterator(cfg, flux, &rs);
   TEST_ASSERT(fr != nullptr, "factory returns a non-null iterator");
   TEST_ASSERT(fr->WaveOpLaw() == FaultFrictionLaw::RateAndState,
               "slip_srw iterator WaveOpLaw() == RateAndState");
   TEST_ASSERT(dynamic_cast<SlipLawSRWFrictionIterator*>(fr.get()) != nullptr,
               "slip_srw iterator is a SlipLawSRWFrictionIterator");
}

// =====================================================================
// F7: slip_srw with rs == nullptr aborts (the SRW branch needs per-DOF V_w)
// =====================================================================
static void F7_slip_srw_null_rs_aborts()
{
   std::cout << "\n[F7] slip_srw with rs==nullptr aborts\n";
   const bool aborts = RunInChild_([]() {
      FaultFaceFlux flux(kRho, kCp, kCs);
      const spatial::SpatialFrictionConfig cfg =
         MakeSlipSRWConfig(FrictionSolver::V0);
      auto fr = MakeFrictionIterator(cfg, flux, /*rs=*/nullptr);
      (void) fr;
   });
   TEST_ASSERT(aborts, "slip_srw with rs==nullptr aborts");
}

// =====================================================================
// F8: slip_srw with V_w.Size() != a.Size() aborts (unresolved V_w guard)
// =====================================================================
static void F8_slip_srw_vw_mismatch_aborts()
{
   std::cout << "\n[F8] slip_srw with unfilled V_w aborts\n";
   const bool aborts = RunInChild_([]() {
      FaultFaceFlux flux(kRho, kCp, kCs);
      const spatial::SpatialFrictionConfig cfg =
         MakeSlipSRWConfig(FrictionSolver::V0);
      const spatial::RateStatePerDOFParams rs =
         MakeRsWithVw(2, /*fill_vw=*/false);   // V_w left empty (Size()==0)
      auto fr = MakeFrictionIterator(cfg, flux, &rs);
      (void) fr;
   });
   TEST_ASSERT(aborts, "slip_srw with V_w.Size() != a.Size() aborts");
}

// =====================================================================
// F9: slip_srw with N==0 fault DOFs builds (no abort) — R-001 regression.
//     A rank that owns no fault DOFs has rs.a.Size()==rs.V_w.Size()==0; the
//     factory is called on EVERY rank (the iterator is needed for SetSubSteps),
//     so it must build rather than MFEM_VERIFY-abort the whole MPI job.
// =====================================================================
static void F9_slip_srw_zero_fault_dofs_builds()
{
   std::cout << "\n[F9] slip_srw with N==0 fault DOFs builds (R-001)\n";
   FaultFaceFlux flux(kRho, kCp, kCs);
   const spatial::SpatialFrictionConfig cfg = MakeSlipSRWConfig(FrictionSolver::V0);
   const spatial::RateStatePerDOFParams rs = MakeRsWithVw(/*n=*/0, /*fill_vw=*/true);
   std::unique_ptr<IFrictionIterator> fr = MakeFrictionIterator(cfg, flux, &rs);
   TEST_ASSERT(fr != nullptr, "factory returns a non-null iterator at N==0");
   TEST_ASSERT(fr->WaveOpLaw() == FaultFrictionLaw::RateAndState,
               "N==0 slip_srw iterator WaveOpLaw() == RateAndState");
   TEST_ASSERT(dynamic_cast<SlipLawSRWFrictionIterator*>(fr.get()) != nullptr,
               "N==0 slip_srw iterator is a SlipLawSRWFrictionIterator (no abort)");
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 2 test_friction_iterator_factory\n";
   F1_slip_weakening();
   F2_rate_state_aging();
   F3_rate_state_missing_block_aborts();
   F4_v0_guard();
   F5_lsw_empty_callback_rejected();
   F6_slip_srw_builds();
   F7_slip_srw_null_rs_aborts();
   F8_slip_srw_vw_mismatch_aborts();
   F9_slip_srw_zero_fault_dofs_builds();

   std::cout << "\n========================================\n";
   std::cout << "Phase 2 test_friction_iterator_factory: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}
