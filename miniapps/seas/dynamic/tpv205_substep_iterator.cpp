// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Implementation of TPV205 sub-step fault iterator (linear slip-weakening).

#include "tpv205_substep_iterator.hpp"
#include "tpv205_friction.hpp"

#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>
#include <limits>

namespace mfem
{
namespace seas
{

Tpv205SubStepIterator::Tpv205SubStepIterator(FaultFaceFlux &flux)
   : flux_(flux)
{
}

void Tpv205SubStepIterator::SetSubSteps(std::vector<real_t> deltaT,
                                        std::vector<real_t> time_weights)
{
   if (deltaT.empty())
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::SetSubSteps: deltaT must have at least "
         "one entry.");
   }
   if (deltaT.size() != time_weights.size())
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::SetSubSteps: deltaT and time_weights "
         "must have equal size; got deltaT.size() = "
         + std::to_string(deltaT.size())
         + ", time_weights.size() = "
         + std::to_string(time_weights.size()));
   }
   for (size_t o = 0; o < deltaT.size(); ++o)
   {
      if (!std::isfinite(deltaT[o]) || deltaT[o] <= 0.0)
      {
         throw std::runtime_error(
            "Tpv205SubStepIterator::SetSubSteps: deltaT[" + std::to_string(o)
            + "] must be finite and positive; got " + std::to_string(deltaT[o]));
      }
      if (!std::isfinite(time_weights[o]))
      {
         throw std::runtime_error(
            "Tpv205SubStepIterator::SetSubSteps: time_weights["
            + std::to_string(o) + "] must be finite; got "
            + std::to_string(time_weights[o]));
      }
   }
   const real_t wsum = std::accumulate(time_weights.begin(),
                                       time_weights.end(),
                                       static_cast<real_t>(0));
   if (std::abs(wsum - 1.0) > 1e-12)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::SetSubSteps: time_weights must sum to "
         "1 within 1e-12; got sum = " + std::to_string(wsum));
   }

   deltaT_       = std::move(deltaT);
   time_weights_ = std::move(time_weights);
}

void Tpv205SubStepIterator::StepOneQP_(DOFData &d,
                                       const real_t *Q_plus,
                                       const real_t *Q_minus,
                                       real_t dt_sub,
                                       bool last_sub_step,
                                       EvalStageState &s,
                                       real_t *Q_imp_plus,
                                       real_t *Q_imp_minus)
{
   // Step 1: trial traction (pure helper — no friction solve, no Newton).
   FaultFaceFlux::ComputeTrialTraction(d, Q_plus, Q_minus,
                                       s.sigma_n_trial,
                                       s.tau1_trial,
                                       s.tau2_trial);

   // Step 2: total traction (TPV205 has zero nucleation channels).
   s.sigma_n_total = d.sigma_n0 + d.sigma_n_nuc + s.sigma_n_trial;
   s.tau1_total    = d.tau1_0   + d.tau1_nuc    + s.tau1_trial;
   s.tau2_total    = d.tau2_0   + d.tau2_nuc    + s.tau2_trial;
   // R-008: keep `s.Theta = |τ_total|` populated so downstream
   // diagnostic prints (which assume the rate-and-state invariant)
   // see the correct magnitude rather than the default 0.
   s.Theta         = std::sqrt(s.tau1_total * s.tau1_total
                               + s.tau2_total * s.tau2_total);

   // Step 3: LSW μ(δ) at current slip magnitude δ = sqrt(slip1² + slip2²).
   // R-016: read LSW-native fields directly — no field repurposing.
   const real_t delta = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
   const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
                                                       d.lsw_mu_s,
                                                       d.lsw_mu_d,
                                                       d.lsw_d_c);

   // Step 4: closed-form LSW solve.  Sets s.V_abs, s.V1, s.V2 and the
   // trial-scale corrected traction s.tau1_corr, s.tau2_corr.
   SolveLSW_TPV205(s.tau1_trial, s.tau2_trial,
                   s.tau1_total, s.tau2_total,
                   s.sigma_n_total, d.eta_s,
                   mu_eff,
                   s.V_abs, s.V1, s.V2,
                   s.tau1_corr, s.tau2_corr);

   // sigma_n is unaffected by friction; sigma_n_corr in EvalStageState
   // is the TRIAL-scale value (matches the rate-and-state path's
   // convention — `BuildImposedState` reads s.sigma_n_corr and writes
   // it into Q_imp[SXX] as the imposed normal traction).
   s.sigma_n_corr = s.sigma_n_trial;

   // Step 5: per-sub-step slip accumulation.  Mirrors TPV104's
   // R4-001 fix — slip evolves on the same grid as the friction solve.
   d.slip1 += s.V1 * dt_sub;
   d.slip2 += s.V2 * dt_sub;

   // Step 6: imposed state (Eq. 11-12 of the FaultFaceFlux pipeline).
   flux_.BuildImposedState(d, s, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // Step 7: on the final sub-step, write back V/slip_rate/tau*_corr/
   // sigma_n_corr to DOFData (with TOTAL physical traction = pre + nuc +
   // trial-scale corrected, per the FaultFaceFlux convention).
   if (last_sub_step)
   {
      flux_.WriteBackState(d, s);
   }
}

void Tpv205SubStepIterator::Advance(std::vector<DOFData> &dof_data,
                                    const std::vector<Vector> &fault_coords,
                                    const real_t *I_plus_flat,
                                    const real_t *I_minus_flat,
                                    real_t dt_macro,
                                    real_t /*t_macro_start*/,
                                    real_t *I_imp_plus_flat,
                                    real_t *I_imp_minus_flat)
{
   if (deltaT_.empty())
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::Advance: SetSubSteps has not been "
         "called; cannot iterate without a configured quadrature.");
   }
   if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::Advance: dt_macro must be finite and "
         "positive; got " + std::to_string(dt_macro));
   }
   if (I_plus_flat == nullptr || I_minus_flat == nullptr
       || I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::Advance: all four I_* pointers must be "
         "non-null.");
   }
   const int n = static_cast<int>(dof_data.size());
   if (static_cast<int>(fault_coords.size()) != n)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::Advance: dof_data and fault_coords "
         "must have equal size; got "
         + std::to_string(dof_data.size()) + ", "
         + std::to_string(fault_coords.size()));
   }

   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   const real_t rel = std::abs(dtsum - dt_macro)
                      / std::max(dt_macro, static_cast<real_t>(1e-300));
   const int O = static_cast<int>(deltaT_.size());
   const real_t sum_tol = std::max<real_t>(
      1e-12,
      static_cast<real_t>(10.0) * O
         * std::numeric_limits<real_t>::epsilon());
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::Advance: Σ deltaT[o] must equal "
         "dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro)
         + ", O = " + std::to_string(O));
   }

   const size_t nwords = static_cast<size_t>(NUM_STATE) * n;
   std::memset(I_imp_plus_flat,  0, nwords * sizeof(real_t));
   std::memset(I_imp_minus_flat, 0, nwords * sizeof(real_t));

   const real_t inv_dt_macro = 1.0 / dt_macro;

   real_t Q_avg_plus[NUM_STATE];
   real_t Q_avg_minus[NUM_STATE];
   real_t Q_imp_plus[NUM_STATE];
   real_t Q_imp_minus[NUM_STATE];

   for (int o = 0; o < O; ++o)
   {
      const real_t dt_sub      = deltaT_[o];
      const real_t weight      = time_weights_[o];
      const real_t accum_scale = weight * dt_macro;
      const bool last_sub_step = (o == O - 1);

      for (int i = 0; i < n; ++i)
      {
         DOFData      &d  = dof_data[i];
         const real_t *Ip = I_plus_flat  + static_cast<ptrdiff_t>(i) * NUM_STATE;
         const real_t *Im = I_minus_flat + static_cast<ptrdiff_t>(i) * NUM_STATE;

         for (int c = 0; c < NUM_STATE; ++c)
         {
            Q_avg_plus[c]  = Ip[c] * inv_dt_macro;
            Q_avg_minus[c] = Im[c] * inv_dt_macro;
         }

         EvalStageState s;
         StepOneQP_(d, Q_avg_plus, Q_avg_minus, dt_sub, last_sub_step,
                    s, Q_imp_plus, Q_imp_minus);

         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }
      }
   }
}

void Tpv205SubStepIterator::AdvanceWithSubStepStates(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
   real_t dt_macro,
   real_t t_macro_start,
   real_t *I_imp_plus_flat,
   real_t *I_imp_minus_flat)
{
   // Phase N (R-N-005): route through the callback overload with a
   // no-op callback so the original API stays bit-identical for any
   // test still calling this signature.
   AdvanceWithSubStepStates(dof_data, fault_coords,
                            Q_pointwise_plus_per_substep,
                            Q_pointwise_minus_per_substep,
                            dt_macro, t_macro_start,
                            I_imp_plus_flat, I_imp_minus_flat,
                            [](real_t, real_t){});
}

void Tpv205SubStepIterator::AdvanceWithSubStepStates(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
   real_t dt_macro,
   real_t t_macro_start,
   real_t *I_imp_plus_flat,
   real_t *I_imp_minus_flat,
   const std::function<void(real_t, real_t)> &nuc_callback)
{
   if (deltaT_.empty())
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::AdvanceWithSubStepStates: SetSubSteps "
         "has not been called.");
   }
   if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::AdvanceWithSubStepStates: dt_macro "
         "must be finite and positive; got " + std::to_string(dt_macro));
   }
   if (I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::AdvanceWithSubStepStates: I_imp_*_flat "
         "must both be non-null.");
   }

   const int O = static_cast<int>(deltaT_.size());
   if (static_cast<int>(Q_pointwise_plus_per_substep.size()) != O ||
       static_cast<int>(Q_pointwise_minus_per_substep.size()) != O)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::AdvanceWithSubStepStates: Q_pointwise_* "
         "must each have size O = " + std::to_string(O));
   }

   const int n = static_cast<int>(dof_data.size());
   if (static_cast<int>(fault_coords.size()) != n)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::AdvanceWithSubStepStates: dof_data and "
         "fault_coords must have equal size; got "
         + std::to_string(dof_data.size()) + ", "
         + std::to_string(fault_coords.size()));
   }

   const size_t expected_words =
      static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n);
   for (int o = 0; o < O; ++o)
   {
      if (Q_pointwise_plus_per_substep[o].size() != expected_words ||
          Q_pointwise_minus_per_substep[o].size() != expected_words)
      {
         throw std::runtime_error(
            "Tpv205SubStepIterator::AdvanceWithSubStepStates: "
            "Q_pointwise_*[" + std::to_string(o) + "] must have size "
            "NUM_STATE * n.");
      }
   }

   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   const real_t rel = std::abs(dtsum - dt_macro)
                      / std::max(dt_macro, static_cast<real_t>(1e-300));
   const real_t sum_tol = std::max<real_t>(
      1e-12,
      static_cast<real_t>(10.0) * O
         * std::numeric_limits<real_t>::epsilon());
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv205SubStepIterator::AdvanceWithSubStepStates: Σ deltaT "
         "must equal dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro));
   }

   const size_t nwords = expected_words;
   std::memset(I_imp_plus_flat,  0, nwords * sizeof(real_t));
   std::memset(I_imp_minus_flat, 0, nwords * sizeof(real_t));

   real_t Q_imp_plus[NUM_STATE];
   real_t Q_imp_minus[NUM_STATE];

   // Track running sub-step end time so the nucleation callback gets
   // the absolute simulation time `t_macro_start + Σ_{o'<=o} dt_sub`.
   real_t t_substep_end = t_macro_start;

   for (int o = 0; o < O; ++o)
   {
      const real_t dt_sub      = deltaT_[o];
      const real_t weight      = time_weights_[o];
      const real_t accum_scale = weight * dt_macro;
      const bool last_sub_step = (o == O - 1);

      // Phase N hook: nucleation callback fires BEFORE the per-QP
      // friction pipeline (mirrors tpv104_substep_iterator.cpp:295).
      // Sub-step end time is `t_macro_start + Σ_{k≤o} dt_k`.
      t_substep_end += dt_sub;
      nuc_callback(t_substep_end, dt_sub);

      const real_t *Qp_o = Q_pointwise_plus_per_substep[o].data();
      const real_t *Qm_o = Q_pointwise_minus_per_substep[o].data();

      for (int i = 0; i < n; ++i)
      {
         DOFData &d = dof_data[i];

         const real_t *Q_tilde_plus  =
            Qp_o + static_cast<ptrdiff_t>(i) * NUM_STATE;
         const real_t *Q_tilde_minus =
            Qm_o + static_cast<ptrdiff_t>(i) * NUM_STATE;

         EvalStageState s;
         StepOneQP_(d, Q_tilde_plus, Q_tilde_minus, dt_sub,
                    last_sub_step, s, Q_imp_plus, Q_imp_minus);

         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }
      }
   }
}

} // namespace seas
} // namespace mfem
