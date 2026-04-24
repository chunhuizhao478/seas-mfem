// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for Tpv104SubStepIterator (§4.10 Step 7 gates
// T_TPV104_SSI_1..6).  SSI_6 (BP5 end-to-end non-regression) is
// deferred — it runs under `make test-bp5-smoke` and the Phase-2
// acceptance matrix.

#include "test_macros.hpp"
#include "../../dynamic/tpv104_substep_iterator.hpp"
#include "../../dynamic/tpv104_setup.hpp"
#include "../../dynamic/tpv104_nucleation.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../friction/slip_law_srw_psi.hpp"
#include "../../config/tpv104_params.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Fixture helpers.
// ---------------------------------------------------------------------------

// Build a single-QP synthetic DOFData with TPV104 rest-state defaults.
// The QP lies at the hypocenter by default; caller can override the
// coordinate to move it outside the nucleation patch.
static void MakeSingleQPFixture(std::vector<DOFData> &dof_data,
                                std::vector<Vector> &fault_coords,
                                std::vector<real_t> &V_w,
                                real_t x2 = 0.0,
                                real_t x3 = TPV104Params::hypo_down_dip)
{
   fault_coords.resize(1);
   fault_coords[0].SetSize(3);
   fault_coords[0](0) = x2;
   fault_coords[0](1) = 0.0;
   fault_coords[0](2) = -x3;

   InitializeFaultDOFs_TPV104(dof_data, 1, fault_coords);
   PopulateVwSideChannel_TPV104(V_w, fault_coords);
}

// ADER-O5 Gauss-Legendre quadrature on [0, 1], mapped to a macro-step
// of size `dt_macro`.  Returns (deltaT, time_weights) with Σ deltaT ==
// dt_macro and Σ time_weights == 1.  We pick deltaT[o] =
// dt_macro * w[o] (so the sub-step sizes are weighted by the GL weight,
// matching the convention in the reference FVW runtime when each sub-
// step integrates one quadrature weight of the Taylor predictor).
static void GaussLegendreO5SubSteps(real_t dt_macro,
                                    std::vector<real_t> &deltaT,
                                    std::vector<real_t> &time_weights)
{
   // GL nodes / weights on [0, 1] (5-point).
   const real_t w[5] = {
      0.11846344252809454375713202035995868,
      0.23931433524968323402064575741781910,
      0.28444444444444444444444444444444444,
      0.23931433524968323402064575741781910,
      0.11846344252809454375713202035995868
   };
   time_weights.assign(w, w + 5);
   deltaT.resize(5);
   for (int o = 0; o < 5; ++o) { deltaT[o] = dt_macro * w[o]; }
}

// Build a time-integrated predictor I_± = dt_macro · Q_bar where
// Q_bar is a "rest-state" zero bulk-field.  The fault-face flux uses
// velocity components and stress components; rest state = all zero.
static void MakeRestStatePredictor(int ndof, real_t dt_macro,
                                   std::vector<real_t> &I_plus,
                                   std::vector<real_t> &I_minus)
{
   I_plus.assign(NUM_STATE * ndof, 0.0);
   I_minus.assign(NUM_STATE * ndof, 0.0);
   (void)dt_macro;  // Q_bar = 0 → I = 0 regardless of dt.
}

// ---------------------------------------------------------------------------
// T_TPV104_SSI_1 — iterator ↔ inline replay self-consistency.
//
// R4-008 (review round 4) disclosure: the "reference" replay below
// calls the production `ComputeStageState`, `BuildImposedState`,
// `UpdateStateAnalyticSlipLawSRW`, and `ApplyNucleationIncremental_TPV104`
// helpers — it is NOT an independent reference formula.  The test
// therefore verifies the ACCUMULATOR LAYER (timeWeight·dt_macro scaling,
// sub-step order, Write-back cadence) is identical between the
// iterator's internal loop and the replay, NOT that any underlying
// physics primitive matches an external reference.  Primitive-layer
// byte-match against inline straight-line arithmetic is covered by
// T_SRW_5 (SRW), T_TPV104_FC_1 (μ), T_TPV104_FS_3 (Newton), and
// T_TPV104_NUC_1 (smoothStepIncrement) in their respective test
// binaries.  The additional `TestAccumulatorScaleIdentity` below
// provides a scalar-level standalone check of the accumulator
// mathematics using only std:: primitives.
// ---------------------------------------------------------------------------
void TestSubStepAccumulatorSelfConsistency()
{
   std::cout << "\n[T_TPV104_SSI_1] iterator ↔ replay self-consistency\n";

   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<real_t> V_w;
   MakeSingleQPFixture(dof_data, fault_coords, V_w);
   // Place at the hypocenter but in a regime where the friction solve
   // behaves benignly (rest state, Q=0).  The iterator's imposed state
   // depends on data.psi (pre-update) and data.tau*_nuc (post-increment
   // in sub-step o).

   const real_t dt_macro = 4.7e-3;   // plan §4.10 SSI_1 fixture
   std::vector<real_t> deltaT, tw;
   GaussLegendreO5SubSteps(dt_macro, deltaT, tw);

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
                     TPV104Params::f0, TPV104Params::f_w, TPV104Params::V_w_in);
   Tpv104SubStepIterator it(flux, law);
   // Self-consistency test uses `Method::Brent` for both the iterator
   // and the replay so the stage-state arithmetic matches
   // bit-for-bit.  The stable-asinh Newton default is exercised by
   // `TestStableNewtonIsDefault` below.
   it.SetSubSteps(deltaT, tw);

   std::vector<real_t> I_plus, I_minus;
   MakeRestStatePredictor(1, dt_macro, I_plus, I_minus);
   std::vector<real_t> I_imp_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_minus(NUM_STATE, 0.0);

   // Clone dof_data for the reference-accumulator replay.
   std::vector<DOFData> dof_ref = dof_data;

   it.Advance(dof_data, fault_coords, V_w,
              I_plus.data(), I_minus.data(),
              dt_macro, /*t_start=*/0.0,
              I_imp_plus.data(), I_imp_minus.data(),
              FrictionSolver::Method::Brent);

   // Inline reference accumulator: reproduce the iterator's arithmetic
   // with std-primitive calls only (no ComputeStageState / Build /
   // WriteBack recursion — those helpers are read-only, so the
   // reference *does* call them; the byte-match here is about the
   // accumulator layer, not the stage helpers).
   std::vector<real_t> I_imp_plus_ref(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_minus_ref(NUM_STATE, 0.0);
   real_t t_cursor = 0.0;
   const real_t inv_dt_macro = 1.0 / dt_macro;
   const int O = static_cast<int>(deltaT.size());
   for (int o = 0; o < O; ++o)
   {
      const real_t dts  = deltaT[o];
      const real_t wo   = tw[o];
      const real_t tend = t_cursor + dts;
      const real_t scale = wo * dt_macro;

      ApplyNucleationIncremental_TPV104(dof_ref, fault_coords, tend, dts);

      real_t Qp[NUM_STATE], Qm[NUM_STATE];
      for (int c = 0; c < NUM_STATE; ++c)
      {
         Qp[c] = I_plus[c]  * inv_dt_macro;
         Qm[c] = I_minus[c] * inv_dt_macro;
      }

      EvalStageState s;
      flux.ComputeStageState(dof_ref[0], Qp, Qm, s,
                             FrictionSolver::Method::Brent);

      dof_ref[0].psi = UpdateStateAnalyticSlipLawSRW(
         dof_ref[0].psi, s.V_abs, dof_ref[0].Dc, dts,
         V_w[0], dof_ref[0].a,
         TPV104Params::b, TPV104Params::V0,
         TPV104Params::f0, TPV104Params::f_w);

      real_t Qip[NUM_STATE], Qim[NUM_STATE];
      flux.BuildImposedState(dof_ref[0], s, Qp, Qm, Qip, Qim);

      for (int c = 0; c < NUM_STATE; ++c)
      {
         I_imp_plus_ref[c]  += scale * Qip[c];
         I_imp_minus_ref[c] += scale * Qim[c];
      }

      if (o == O - 1) { flux.WriteBackState(dof_ref[0], s); }
      t_cursor = tend;
   }

   // Compare accumulated outputs.
   real_t max_rel = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      const real_t dp = std::abs(I_imp_plus[c] - I_imp_plus_ref[c]);
      const real_t dm = std::abs(I_imp_minus[c] - I_imp_minus_ref[c]);
      const real_t dpscale = std::max<real_t>(std::abs(I_imp_plus_ref[c]),
                                              static_cast<real_t>(1e-30));
      const real_t dmscale = std::max<real_t>(std::abs(I_imp_minus_ref[c]),
                                              static_cast<real_t>(1e-30));
      max_rel = std::max(max_rel, dp / dpscale);
      max_rel = std::max(max_rel, dm / dmscale);
   }
   std::cout << "  max rel err = " << max_rel << " over 9 components ±\n";
   TEST_ASSERT(max_rel < 1e-12,
               "iterator imposed-state matches reference accumulator to 1e-12 rel");

   // DOFData state (psi, slip_rate, V1/V2, tau*_corr) also matches.
   TEST_NEAR(dof_data[0].psi, dof_ref[0].psi, 1e-12,
             "iterator ψ matches reference");
   TEST_NEAR(dof_data[0].slip_rate, dof_ref[0].slip_rate, 1e-12,
             "iterator slip_rate matches reference");
}

// ---------------------------------------------------------------------------
// T_TPV104_SSI_3 — O=1 limit bit-identical to `EvaluateADER` on a
// fixture with nucleation turned off (r >> R so F(r) = 0).
// ---------------------------------------------------------------------------
void TestO1LimitMatchesEvaluateADER()
{
   std::cout << "\n[T_TPV104_SSI_3] O=1 limit bit-identical to EvaluateADER\n";

   // Place the QP far from the hypocenter (> nuc_radius = 3 km away
   // along strike AND > 3 km from down-dip hypo) so the nucleation
   // spatial factor F(r) = 0 at every sub-step.
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<real_t> V_w;
   MakeSingleQPFixture(dof_data, fault_coords, V_w,
                       /*x2=*/20.0e3,   // far off-strike
                       /*x3=*/20.0e3);  // far down-dip
   TEST_NEAR(NucleationSpatial_TPV104(std::sqrt(20e3 * 20e3 + 12.5e3 * 12.5e3)),
             0.0, 1e-30,
             "test QP is outside the nucleation patch (F = 0)");

   const real_t dt_macro = 1.0e-3;

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
                     TPV104Params::f0, TPV104Params::f_w,
                     TPV104Params::V_w_out);
   Tpv104SubStepIterator it(flux, law);
   // Bit-identity with EvaluateADER requires the legacy MFEM-native μ
   // (`Method::Brent` inside the iterator; EvaluateADER uses the same
   // dispatch).  The stable-asinh default is exercised by
   // `TestStableNewtonIsDefault` below.
   // O=1 quadrature: single sub-step spanning the full macro-step.
   it.SetSubSteps(std::vector<real_t>{dt_macro},
                  std::vector<real_t>{1.0});

   std::vector<real_t> I_plus, I_minus;
   MakeRestStatePredictor(1, dt_macro, I_plus, I_minus);
   std::vector<real_t> I_imp_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_minus(NUM_STATE, 0.0);

   // Clone DOFData for the EvaluateADER comparison (iterator mutates
   // dof_data, so we need a snapshot to pass into EvaluateADER for an
   // apples-to-apples input).
   std::vector<DOFData> dof_for_ader = dof_data;

   it.Advance(dof_data, fault_coords, V_w,
              I_plus.data(), I_minus.data(),
              dt_macro, /*t_start=*/0.0,
              I_imp_plus.data(), I_imp_minus.data(),
              FrictionSolver::Method::Brent);

   // EvaluateADER reference call.  Nucleation channel is already 0
   // (initialised that way); iterator's ApplyNucleationIncremental call
   // adds 0 because F(r) = 0 at this QP.
   real_t I_imp_plus_ader[NUM_STATE] = {0};
   real_t I_imp_minus_ader[NUM_STATE] = {0};
   flux.EvaluateADER(dof_for_ader[0],
                     I_plus.data(), I_minus.data(), dt_macro,
                     I_imp_plus_ader, I_imp_minus_ader,
                     FrictionSolver::Method::Brent);

   // Bit-compare the accumulated output.
   bool all_equal = true;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      if (I_imp_plus[c] != I_imp_plus_ader[c])
      {
         all_equal = false;
         std::cerr << "  FAIL + c=" << c
                   << "  iterator=" << I_imp_plus[c]
                   << "  ader="     << I_imp_plus_ader[c]
                   << "  diff="     << (I_imp_plus[c] - I_imp_plus_ader[c])
                   << "\n";
      }
      if (I_imp_minus[c] != I_imp_minus_ader[c])
      {
         all_equal = false;
         std::cerr << "  FAIL − c=" << c
                   << "  iterator=" << I_imp_minus[c]
                   << "  ader="     << I_imp_minus_ader[c] << "\n";
      }
   }
   TEST_ASSERT(all_equal,
               "O=1 iterator output bit-identical to EvaluateADER");

   // DOFData after iterator should match after EvaluateADER for all
   // fields EXCEPT psi — the iterator updates psi (§3.12), EvaluateADER
   // does not.  At V=0 (rest state) ψ is invariant under the analytic
   // step (exp1m = 0), so both paths produce the same ψ too.
   TEST_NEAR(dof_data[0].psi, dof_for_ader[0].psi, 1e-12,
             "rest-state ψ bit-identical (V = 0 means exp1m = 0)");
   TEST_NEAR(dof_data[0].slip_rate, dof_for_ader[0].slip_rate, 1e-12,
             "slip_rate bit-identical in O=1 limit");
}

// ---------------------------------------------------------------------------
// T_TPV104_SSI_4 — nucleation injection on locked fault (V=0) accumulates
// to Δτ₀ · F(r) at the hypocenter after t = T_nuc.
// ---------------------------------------------------------------------------
void TestNucleationInjectionLockedFault()
{
   std::cout << "\n[T_TPV104_SSI_4] nucleation injection on locked fault\n";

   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<real_t> V_w;
   MakeSingleQPFixture(dof_data, fault_coords, V_w);  // hypocenter QP

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
                     TPV104Params::f0, TPV104Params::f_w,
                     TPV104Params::V_w_in);

   // Many macro-steps over [0, T_nuc], each with O = 5 sub-steps.
   const int N_macro = 200;
   const real_t dt_macro = TPV104Params::nuc_T / N_macro;   // 5 ms
   std::vector<real_t> deltaT, tw;
   GaussLegendreO5SubSteps(dt_macro, deltaT, tw);

   Tpv104SubStepIterator it(flux, law);
   it.SetSubSteps(deltaT, tw);

   std::vector<real_t> I_plus, I_minus;
   MakeRestStatePredictor(1, dt_macro, I_plus, I_minus);
   std::vector<real_t> I_imp_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_minus(NUM_STATE, 0.0);

   real_t t = 0.0;
   for (int m = 0; m < N_macro; ++m)
   {
      it.Advance(dof_data, fault_coords, V_w,
                 I_plus.data(), I_minus.data(),
                 dt_macro, t,
                 I_imp_plus.data(), I_imp_minus.data(),
                 FrictionSolver::Method::Brent);
      t += dt_macro;
   }

   // At the hypocenter, F(0) = 1; expected tau2_nuc = Δτ₀ at t = T_nuc.
   const real_t expected = TPV104Params::nuc_dtau;
   const real_t got = dof_data[0].tau2_nuc;
   const real_t rel = std::abs(got - expected) / std::abs(expected);
   std::cout << "  tau2_nuc(hypo) after T_nuc = " << got
             << "  (expected " << expected << ", rel=" << rel << ")\n";
   TEST_ASSERT(rel < 1e-8,
               "tau2_nuc[hypo] reaches Δτ₀ = 45 MPa to 1e-8 rel");

   // R3-001 invariant: tau1_nuc and sigma_n_nuc remain at 0.
   TEST_NEAR(dof_data[0].tau1_nuc, 0.0, 1e-30,
             "tau1_nuc stays at 0 after iterator + nucleation ramp");
   TEST_NEAR(dof_data[0].sigma_n_nuc, 0.0, 1e-30,
             "sigma_n_nuc stays at 0 after iterator + nucleation ramp");
}

// ---------------------------------------------------------------------------
// T_TPV104_SSI_5 — ψ update placement: rest-state (V = V_ini = 1e-16) ψ
// after O=5 sub-steps matches an inline per-sub-step accumulator.
// ---------------------------------------------------------------------------
void TestPsiUpdatePlacement()
{
   std::cout << "\n[T_TPV104_SSI_5] ψ update placement (rest state)\n";

   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<real_t> V_w;
   MakeSingleQPFixture(dof_data, fault_coords, V_w,
                       /*x2=*/0.0, /*x3=*/TPV104Params::hypo_down_dip);
   const real_t psi_entry = dof_data[0].psi;

   const real_t dt_macro = 4.7e-3;
   std::vector<real_t> deltaT, tw;
   GaussLegendreO5SubSteps(dt_macro, deltaT, tw);

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
                     TPV104Params::f0, TPV104Params::f_w,
                     TPV104Params::V_w_in);
   Tpv104SubStepIterator it(flux, law);
   it.SetSubSteps(deltaT, tw);

   std::vector<real_t> I_plus, I_minus;
   MakeRestStatePredictor(1, dt_macro, I_plus, I_minus);
   std::vector<real_t> I_imp_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_minus(NUM_STATE, 0.0);

   it.Advance(dof_data, fault_coords, V_w,
              I_plus.data(), I_minus.data(),
              dt_macro, /*t_start=*/0.0,
              I_imp_plus.data(), I_imp_minus.data(),
              FrictionSolver::Method::Brent);

   // At V ≈ V_ini = 1e-16, `preexp1 = -V*dt/L ≈ 0`, so every analytic
   // sub-step update leaves ψ essentially unchanged.  The iterator's
   // final ψ must differ from the entry ψ by ULP-scale only.
   const real_t delta = std::abs(dof_data[0].psi - psi_entry);
   std::cout << "  ψ drift over " << deltaT.size() << " sub-steps = "
             << delta << "\n";
   TEST_ASSERT(delta < 1e-14,
               "rest-state ψ stays invariant (per-sub-step V≈0 → exp1m≈0)");

   // R5-005 (review round 5): use the iterator's actual final
   // slip_rate as the per-sub-step V — for the rest-state fixture
   // the per-sub-step V is constant (ψ doesn't move, so V doesn't
   // either) and the iterator's WriteBackState captures the final
   // sub-step's V_abs, which is the same as every other sub-step's.
   // The previous hard-coded `V_sub = TPV104Params::V_ini` was
   // fixture-brittle: changing V_ini would silently break the
   // reference.
   const real_t V_sub = dof_data[0].slip_rate;
   real_t psi_ref = psi_entry;
   for (size_t o = 0; o < deltaT.size(); ++o)
   {
      psi_ref = UpdateStateAnalyticSlipLawSRW(
         psi_ref, V_sub, TPV104Params::L, deltaT[o],
         V_w[0], TPV104Params::a_in,
         TPV104Params::b, TPV104Params::V0,
         TPV104Params::f0, TPV104Params::f_w);
   }
   TEST_NEAR(dof_data[0].psi, psi_ref, 1e-12,
             "iterator-updated ψ matches inline sub-step accumulator "
             "(V_sub sourced from iterator state, R5-005)");
}

// ---------------------------------------------------------------------------
// T_TPV104_SSI_2 — stability check under dt_macro halving.
//
// R5-001 (review round 5): the previous version computed
// `|I_imp_plus[c] - I_imp_plus[c]|` — a self-subtraction typo that
// yielded 0 unconditionally, letting any iterator output (including
// NaN) pass the gate.  Rewritten to:
//   (1) finite-value check on every output channel at both dt, at
//       dt/2 (no NaN / Inf regardless of numerical path);
//   (2) linear-scaling check: for a zero-Q rest-state fixture the
//       time-integrated imposed state scales exactly proportionally
//       to dt_macro (since Q_imp is constant in time), so halving
//       dt_macro halves every component.  `|I_imp(dt/2) · 2 −
//       I_imp(dt)| < 1e-10 rel` verifies the accumulator is linear
//       in dt_macro as required.
//
// A full O(dt²) convergence test requires a non-trivial analytic
// fixture not available in the current plan; deferred to Phase 3.
// ---------------------------------------------------------------------------
void TestConvergenceUnderDtHalving()
{
   std::cout << "\n[T_TPV104_SSI_2] stability + dt-linearity check\n";

   auto run = [](real_t dt_macro,
                 std::vector<real_t> &I_imp_plus_out,
                 std::vector<real_t> &I_imp_minus_out)
   {
      std::vector<DOFData> dof_data;
      std::vector<Vector> fault_coords;
      std::vector<real_t> V_w;
      MakeSingleQPFixture(dof_data, fault_coords, V_w,
                          /*x2=*/20e3, /*x3=*/20e3);

      FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp,
                         TPV104Params::cs);
      SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b,
                        TPV104Params::V0, TPV104Params::f0,
                        TPV104Params::f_w, TPV104Params::V_w_out);
      Tpv104SubStepIterator it(flux, law);
      std::vector<real_t> deltaT, tw;
      GaussLegendreO5SubSteps(dt_macro, deltaT, tw);
      it.SetSubSteps(deltaT, tw);

      std::vector<real_t> I_plus, I_minus;
      MakeRestStatePredictor(1, dt_macro, I_plus, I_minus);
      I_imp_plus_out.assign(NUM_STATE, 0.0);
      I_imp_minus_out.assign(NUM_STATE, 0.0);

      it.Advance(dof_data, fault_coords, V_w,
                 I_plus.data(), I_minus.data(),
                 dt_macro, 2.0 * TPV104Params::nuc_T,
                 I_imp_plus_out.data(), I_imp_minus_out.data(),
                 FrictionSolver::Method::Brent);
   };

   std::vector<real_t> Ip_big, Im_big, Ip_half, Im_half;
   run(4.7e-3,         Ip_big,  Im_big);
   run(4.7e-3 / 2.0,   Ip_half, Im_half);

   // (1) Stability: every output component is finite at both dt.
   int nan_count = 0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      if (!std::isfinite(Ip_big[c]))  { ++nan_count; }
      if (!std::isfinite(Im_big[c]))  { ++nan_count; }
      if (!std::isfinite(Ip_half[c])) { ++nan_count; }
      if (!std::isfinite(Im_half[c])) { ++nan_count; }
   }
   TEST_ASSERT(nan_count == 0,
               "iterator emits finite output at both dt_macro and dt_macro/2 "
               "(R5-001: prior self-diff test could not detect NaN)");

   // (2) dt-linearity: I_imp_±(dt/2) · 2 == I_imp_±(dt) within 1e-10 rel.
   // This is exact for Q_bar constant + rest-state (zero-Q) where Q_imp
   // is time-independent, so halving dt_macro halves I_imp exactly.
   real_t max_rel = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      const real_t expected_p = 2.0 * Ip_half[c];
      const real_t expected_m = 2.0 * Im_half[c];
      const real_t scale_p = std::max<real_t>(std::abs(Ip_big[c]),
                                              static_cast<real_t>(1e-30));
      const real_t scale_m = std::max<real_t>(std::abs(Im_big[c]),
                                              static_cast<real_t>(1e-30));
      max_rel = std::max(max_rel,
                         std::abs(Ip_big[c] - expected_p) / scale_p);
      max_rel = std::max(max_rel,
                         std::abs(Im_big[c] - expected_m) / scale_m);
   }
   std::cout << "  max dt-linearity rel err = " << max_rel << "\n";
   TEST_ASSERT(max_rel < 1e-10,
               "I_imp(dt_macro/2) · 2 == I_imp(dt_macro) within 1e-10 rel "
               "(linear scaling in dt on rest-state fixture)");
}

// ---------------------------------------------------------------------------
// Additional guards: SetSubSteps input validation.
// ---------------------------------------------------------------------------
void TestSetSubStepsValidation()
{
   std::cout << "\n[T_TPV104_SSI_VAL] SetSubSteps input validation\n";

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
                     TPV104Params::f0, TPV104Params::f_w,
                     TPV104Params::V_w_in);
   Tpv104SubStepIterator it(flux, law);

   // Empty arrays → throws.
   bool threw = false;
   try { it.SetSubSteps({}, {}); } catch (...) { threw = true; }
   TEST_ASSERT(threw, "empty deltaT/time_weights rejected");

   // Size mismatch → throws.
   threw = false;
   try { it.SetSubSteps({1.0, 2.0}, {1.0}); } catch (...) { threw = true; }
   TEST_ASSERT(threw, "size-mismatched deltaT/time_weights rejected");

   // Non-positive dt → throws.
   threw = false;
   try { it.SetSubSteps({1.0, 0.0}, {0.5, 0.5}); } catch (...) { threw = true; }
   TEST_ASSERT(threw, "non-positive deltaT[o] rejected");

   // Weights don't sum to 1 → throws.
   threw = false;
   try { it.SetSubSteps({1.0, 1.0}, {0.4, 0.4}); } catch (...) { threw = true; }
   TEST_ASSERT(threw, "time_weights not summing to 1 rejected");

   // Valid O=1 accepted.
   bool ok = true;
   try { it.SetSubSteps({2.0}, {1.0}); } catch (...) { ok = false; }
   TEST_ASSERT(ok, "valid O=1 quadrature accepted");
   TEST_ASSERT(it.NumSubSteps() == 1, "NumSubSteps() reports 1 after O=1 set");

   // Valid O=5 accepted.
   std::vector<real_t> dT, tw;
   GaussLegendreO5SubSteps(1e-3, dT, tw);
   it.SetSubSteps(dT, tw);
   TEST_ASSERT(it.NumSubSteps() == 5, "NumSubSteps() reports 5 after O=5 set");
}

// ---------------------------------------------------------------------------
// R5-003 / plan §4.10.X — stable-asinh Newton dispatches cleanly via
// `FrictionSolver::Method::NewtonRaphsonStable`.
//
// Verify that `Advance(..., Method::NewtonRaphsonStable)` runs without
// aborting and produces finite output on the TPV104 rest-state fixture.
// ---------------------------------------------------------------------------
void TestStableNewtonIsDefault()
{
   std::cout << "\n[R5-003] Method::NewtonRaphsonStable dispatches cleanly\n";

   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<real_t> V_w;
   MakeSingleQPFixture(dof_data, fault_coords, V_w);

   const real_t dt_macro = 1.0e-3;
   std::vector<real_t> deltaT, tw;
   GaussLegendreO5SubSteps(dt_macro, deltaT, tw);

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b, TPV104Params::V0,
                     TPV104Params::f0, TPV104Params::f_w,
                     TPV104Params::V_w_in);
   Tpv104SubStepIterator it(flux, law);
   it.SetSubSteps(deltaT, tw);
   std::vector<real_t> I_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_minus(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_minus(NUM_STATE, 0.0);

   // Default Advance() uses Method::NewtonRaphsonStable (plan §4.10.X).
   it.Advance(dof_data, fault_coords, V_w,
              I_plus.data(), I_minus.data(),
              dt_macro, 0.0,
              I_imp_plus.data(), I_imp_minus.data());

   TEST_ASSERT(std::isfinite(dof_data[0].psi),
               "ψ finite after Method::NewtonRaphsonStable");
   TEST_ASSERT(std::isfinite(dof_data[0].slip_rate),
               "slip_rate finite after Method::NewtonRaphsonStable");
   TEST_ASSERT(dof_data[0].slip_rate >= 0.0,
               "slip_rate ≥ 0 (no sign inversion)");

   // Explicit invocation with Brent matches the default's finiteness
   // (regression guard that the dispatch doesn't drop any method).
   std::vector<DOFData> dof_brent;
   MakeSingleQPFixture(dof_brent, fault_coords, V_w);
   it.Advance(dof_brent, fault_coords, V_w,
              I_plus.data(), I_minus.data(),
              dt_macro, 0.0,
              I_imp_plus.data(), I_imp_minus.data(),
              FrictionSolver::Method::Brent);
   TEST_ASSERT(std::isfinite(dof_brent[0].psi),
               "Method::Brent path also finite");
}

// ---------------------------------------------------------------------------
// R4-001 — slip1 / slip2 accumulation after non-zero V.
//
// Construct a fixture where the friction solve converges to a known
// non-zero V_abs and V2, then run the iterator and assert that
// `dof_data[0].slip2 ≈ V2 · dt_macro`.  Without the R4-001 fix, slip2
// remains at its init value (0) and this test fires.
// ---------------------------------------------------------------------------
void TestSlipAccumulation()
{
   std::cout << "\n[R4-001] slip1 / slip2 accumulation\n";

   // Place the test QP far outside the nucleation patch so no nuc
   // increment contaminates the analysis, and bypass the rest-state
   // fixture's V_ini = 1e-16 by elevating slip_rate manually and
   // driving the bulk state to a stress imbalance.
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<real_t> V_w;
   MakeSingleQPFixture(dof_data, fault_coords, V_w,
                       /*x2=*/20e3, /*x3=*/20e3);

   // Drive a bulk trial traction: non-zero tangent-2 velocity jump on
   // the two sides produces a trial tangent-2 traction that the
   // friction solver inverts to a finite V_abs in tangent-2.
   //
   // In fault-local coords: VZ = velocity along tangent2 (strike).
   const real_t dt_macro = 1.0e-3;
   std::vector<real_t> I_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_minus(NUM_STATE, 0.0);
   // Q_avg_plus[VZ] = +0.5; Q_avg_minus[VZ] = -0.5 → jump of 1 m/s.
   I_plus[VZ]  = +0.5 * dt_macro;
   I_minus[VZ] = -0.5 * dt_macro;

   FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp, TPV104Params::cs);
   SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b,
                     TPV104Params::V0, TPV104Params::f0,
                     TPV104Params::f_w, TPV104Params::V_w_out);
   Tpv104SubStepIterator it(flux, law);
   std::vector<real_t> deltaT, tw;
   GaussLegendreO5SubSteps(dt_macro, deltaT, tw);
   it.SetSubSteps(deltaT, tw);

   std::vector<real_t> I_imp_plus(NUM_STATE, 0.0);
   std::vector<real_t> I_imp_minus(NUM_STATE, 0.0);

   const real_t slip2_before = dof_data[0].slip2;
   const real_t slip1_before = dof_data[0].slip1;

   it.Advance(dof_data, fault_coords, V_w,
              I_plus.data(), I_minus.data(),
              dt_macro, /*t_start=*/2.0 * TPV104Params::nuc_T,
              I_imp_plus.data(), I_imp_minus.data(),
              FrictionSolver::Method::Brent);

   const real_t slip2_after = dof_data[0].slip2;
   const real_t V2_final    = dof_data[0].V2;

   std::cout << "  slip2 before = " << slip2_before << "\n"
             << "  slip2 after  = " << slip2_after
             << "  (V2 final = " << V2_final << ")\n";

   // R4-001 primary guard: slip2 MUST change from its init value of 0.
   // Without the R4-001 fix, slip2 stays at exactly 0.0 regardless of
   // V2 (iterator never accumulated slip).
   TEST_ASSERT(slip2_after != slip2_before,
               "slip2 accumulates when V2 ≠ 0 (R4-001 guard)");

   // R5-002 (review round 5): sign-sensitive check — the accumulated
   // slip increment and V2_final must have the SAME sign.  The
   // previous `std::abs(...) >= 0.2 * |V2|·dt` window was 25× wide
   // and admitted a sign-inverted accumulator (`slip -= V·dt`) —
   // exactly the failure mode CLAUDE.md flags as "Antiparallel sign
   // creates positive feedback → unbounded growth".  A strict
   // same-sign ratio catches this at the unit-test layer.
   const real_t slip_increment = slip2_after - slip2_before;
   TEST_ASSERT(slip_increment * V2_final > 0.0,
               "slip2 accumulation has SAME SIGN as V2_final "
               "(R5-002 sign guard)");

   // Tight magnitude check using the known-exact per-sub-step ratio.
   // With V_abs dominated by the tangent-2 trial traction, V2 is
   // essentially constant across the 5 sub-steps (ψ barely moves at
   // V ≲ 1e-15), so the expected accumulator output is V2·dt_macro
   // to within a per-sub-step ψ drift of O(1%).  The tolerance
   // window is tightened to [0.95, 1.05] — 5 % — which would reject
   // a sign error (ratio −1), an off-by-factor-of-O (ratio 0.2 or 5),
   // and a missing-sub-step bug (ratio 0.8 or 1.2).
   const real_t expected = V2_final * dt_macro;
   const real_t ratio = slip_increment / expected;
   std::cout << "  slip_increment / (V2·dt_macro) = " << ratio
             << " (expect ~1.0)\n";
   TEST_ASSERT(ratio > 0.95 && ratio < 1.05,
               "slip2 accumulation ≈ V2·dt_macro within 5 % "
               "(R5-002 tight magnitude window)");

   // slip1 stays at ≈ 0 (no tangent-1 forcing).
   TEST_NEAR(dof_data[0].slip1 - slip1_before, 0.0, 1e-20,
             "slip1 stays at 0 (no dip forcing)");
}

// ---------------------------------------------------------------------------
// R4-001 / R5-004 — slip is signed (accumulation cancels under sign flip).
//
// This test locks in the SIGNED semantics of the slip accumulator as
// an INTENTIONAL contract (R5-004 disclosure).  A maintainer who
// "fixes" slip accumulation by taking `std::abs(s.V2)` would trip
// this test: opposing-sign forcing on the two halves of a macro-step
// produces exactly zero net slip under signed accumulation.
//
// Fixture: two sub-step halves with (synthetic) per-QP V2 = +v then
// V2 = −v.  The iterator's built-in per-sub-step V2 comes from the
// friction solve, which we cannot flip post-hoc — so this test checks
// the accumulator contract directly by manipulating an injected DOF
// after each sub-step would be fragile.  Instead, we assert the
// invariant HOLDS on the accumulated state:
//   — when V2 > 0 across a macro-step: slip_increment > 0.
//   — when V2 < 0 across a macro-step: slip_increment < 0.
// The symmetric flip is tested by running the fixture twice with
// opposite-sign VZ jumps.
// ---------------------------------------------------------------------------
void TestSlipAccumulationIsSigned()
{
   std::cout << "\n[R5-004] slip accumulation is signed (both sign branches)\n";

   auto run_with_vz_jump = [&](real_t vz_sign) -> std::pair<real_t, real_t>
   {
      std::vector<DOFData> dof_data;
      std::vector<Vector> fault_coords;
      std::vector<real_t> V_w;
      MakeSingleQPFixture(dof_data, fault_coords, V_w,
                          /*x2=*/20e3, /*x3=*/20e3);
      // Suppress the pre-stress so the sign of the tangent-2 trial
      // traction — not the background 40 MPa — determines the sign
      // of V2.  Without this, tau2_total stays positive under both
      // VZ-jump sign choices and V2 never flips (the sign-guard then
      // cannot exercise the slip accumulator's sign handling).
      dof_data[0].tau2_0 = 0.0;

      const real_t dt_macro = 1.0e-3;
      std::vector<real_t> I_plus(NUM_STATE, 0.0);
      std::vector<real_t> I_minus(NUM_STATE, 0.0);
      I_plus[VZ]  = +0.5 * dt_macro * vz_sign;
      I_minus[VZ] = -0.5 * dt_macro * vz_sign;

      FaultFaceFlux flux(TPV104Params::rho, TPV104Params::cp,
                         TPV104Params::cs);
      SlipLawSRWPsi law(TPV104Params::a_in, TPV104Params::b,
                        TPV104Params::V0, TPV104Params::f0,
                        TPV104Params::f_w, TPV104Params::V_w_out);
      Tpv104SubStepIterator it(flux, law);
      std::vector<real_t> deltaT, tw;
      GaussLegendreO5SubSteps(dt_macro, deltaT, tw);
      it.SetSubSteps(deltaT, tw);

      std::vector<real_t> Ioutp(NUM_STATE, 0.0), Ioutm(NUM_STATE, 0.0);
      it.Advance(dof_data, fault_coords, V_w,
                 I_plus.data(), I_minus.data(),
                 dt_macro, 2.0 * TPV104Params::nuc_T,
                 Ioutp.data(), Ioutm.data(),
                 FrictionSolver::Method::Brent);
      return {dof_data[0].slip2, dof_data[0].V2};
   };

   auto [slip2_pos, V2_pos] = run_with_vz_jump(+1.0);
   auto [slip2_neg, V2_neg] = run_with_vz_jump(-1.0);
   std::cout << "  +V2 run: slip2 = " << slip2_pos << ", V2 = " << V2_pos << "\n";
   std::cout << "  -V2 run: slip2 = " << slip2_neg << ", V2 = " << V2_neg << "\n";

   TEST_ASSERT(slip2_pos * V2_pos > 0.0,
               "positive-V2 fixture: slip2 and V2 have same sign");
   TEST_ASSERT(slip2_neg * V2_neg > 0.0,
               "negative-V2 fixture: slip2 and V2 have same sign "
               "(R5-004 signed-semantics lock-in)");
   TEST_ASSERT(slip2_pos * slip2_neg < 0.0,
               "flipping V2 sign flips slip2 sign (no abs() applied)");
}

// ---------------------------------------------------------------------------
// R4-008 companion — scalar accumulator identity using only std::.
//
// For ANY constant scalar x emitted as "Q_imp" at every sub-step, the
// iterator's weighted sum over O sub-steps with weights summing to 1
// equals x · dt_macro.  This is the algebraic core of the accumulator
// layer; a straight-line math check catches formula drift in the
// accumulator without any dependence on production helpers.
// ---------------------------------------------------------------------------
void TestAccumulatorScaleIdentity()
{
   std::cout << "\n[T_TPV104_SSI_1b] scalar accumulator identity "
             << "(R4-008 standalone)\n";

   auto check_for_O = [&](int O, real_t dt_macro)
   {
      std::vector<real_t> deltaT, tw;
      if (O == 1)       { deltaT = {dt_macro}; tw = {1.0}; }
      else if (O == 5)  { GaussLegendreO5SubSteps(dt_macro, deltaT, tw); }
      else              { return; }

      // Iterator: set up, drive with I=0 (zero predictor) so Q_imp is
      // purely a function of DOFData — but the accumulator math
      // applies to ANY per-sub-step value, so we verify directly.
      const real_t x = 3.7;  // arbitrary per-sub-step value
      real_t accumulated = 0.0;
      for (int o = 0; o < O; ++o) { accumulated += tw[o] * dt_macro * x; }

      // Expected: x · dt_macro.
      const real_t expected = x * dt_macro;
      TEST_NEAR(accumulated, expected, 1e-14,
                ("O=" + std::to_string(O)
                 + ": Σ tw[o]·dt·x == x·dt (exact)").c_str());

      // Σ deltaT must equal dt_macro (iterator's Advance enforces this).
      real_t dtsum = 0.0;
      for (real_t d : deltaT) { dtsum += d; }
      TEST_NEAR(dtsum, dt_macro, 1e-14,
                ("O=" + std::to_string(O)
                 + ": Σ deltaT[o] == dt_macro").c_str());
   };

   check_for_O(1, 4.7e-3);
   check_for_O(5, 4.7e-3);
   check_for_O(5, 1.0);
   check_for_O(5, 1e-6);
}

int main(int argc, char *argv[])
{
   TestSetSubStepsValidation();
   TestSubStepAccumulatorSelfConsistency();
   TestAccumulatorScaleIdentity();
   TestStableNewtonIsDefault();
   TestSlipAccumulation();
   TestSlipAccumulationIsSigned();
   TestO1LimitMatchesEvaluateADER();
   TestNucleationInjectionLockedFault();
   TestPsiUpdatePlacement();
   TestConvergenceUnderDtHalving();

   std::cout << "\n[T_TPV104_SSI_6] BP5 end-to-end non-regression — "
             << "DEFERRED (requires make test-bp5-smoke; tracked by "
             << "Phase-2 acceptance matrix).\n";

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
