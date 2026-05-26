// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Implementation of TPV102 sub-step fault iterator.  Mirrors
// `tpv104_substep_iterator.cpp` with the AgingLawPsi / no-V_w / aging-
// law analytic ψ substitutions documented in the header.

#include "tpv102_substep_iterator.hpp"
#include "../friction/friction_coeff_stable.hpp"
#include "../friction/state_evolution.hpp"

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

Tpv102SubStepIterator::Tpv102SubStepIterator(FaultFaceFlux &flux,
                                             const AgingLawPsi &state_evo)
   : flux_(flux), state_evo_(state_evo)
{
}

void Tpv102SubStepIterator::SetSubSteps(std::vector<real_t> deltaT,
                                        std::vector<real_t> time_weights)
{
   if (deltaT.empty())
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::SetSubSteps: deltaT must have at least "
         "one entry.");
   }
   if (deltaT.size() != time_weights.size())
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::SetSubSteps: deltaT and time_weights "
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
            "Tpv102SubStepIterator::SetSubSteps: deltaT[" + std::to_string(o)
            + "] must be finite and positive; got " + std::to_string(deltaT[o]));
      }
      if (!std::isfinite(time_weights[o]))
      {
         throw std::runtime_error(
            "Tpv102SubStepIterator::SetSubSteps: time_weights["
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
         "Tpv102SubStepIterator::SetSubSteps: time_weights must sum to "
         "1 within 1e-12; got sum = " + std::to_string(wsum));
   }

   deltaT_       = std::move(deltaT);
   time_weights_ = std::move(time_weights);
}

void Tpv102SubStepIterator::Advance(std::vector<DOFData> &dof_data,
                                    const std::vector<Vector> &fault_coords,
                                    const real_t *I_plus_flat,
                                    const real_t *I_minus_flat,
                                    real_t dt_macro,
                                    real_t t_macro_start,
                                    real_t *I_imp_plus_flat,
                                    real_t *I_imp_minus_flat,
                                    FrictionSolver::Method method)
{
   if (deltaT_.empty())
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::Advance: SetSubSteps has not been "
         "called; cannot iterate without a configured quadrature.");
   }
   if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::Advance: dt_macro must be finite and "
         "positive; got " + std::to_string(dt_macro));
   }
   if (!std::isfinite(t_macro_start))
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::Advance: t_macro_start must be finite; "
         "got " + std::to_string(t_macro_start));
   }
   if (I_plus_flat == nullptr || I_minus_flat == nullptr
       || I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::Advance: all four I_* pointers must be "
         "non-null.");
   }
   const int n = static_cast<int>(dof_data.size());
   if (static_cast<int>(fault_coords.size()) != n)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::Advance: dof_data and fault_coords must "
         "have equal size; got "
         + std::to_string(dof_data.size()) + ", "
         + std::to_string(fault_coords.size()));
   }

   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   const real_t rel = std::abs(dtsum - dt_macro) / std::max(dt_macro,
                                                            1e-300);
   const int    O_size = static_cast<int>(deltaT_.size());
   const real_t sum_tol = std::max<real_t>(1e-12,
                                           static_cast<real_t>(10.0) * O_size
                                             * std::numeric_limits<real_t>::epsilon());
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::Advance: Σ deltaT[o] must equal "
         "dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro)
         + ", O = " + std::to_string(O_size));
   }

   const size_t nwords = static_cast<size_t>(NUM_STATE) * n;
   std::memset(I_imp_plus_flat,  0, nwords * sizeof(real_t));
   std::memset(I_imp_minus_flat, 0, nwords * sizeof(real_t));

   const real_t inv_dt_macro = 1.0 / dt_macro;

   real_t Q_avg_plus[NUM_STATE];
   real_t Q_avg_minus[NUM_STATE];
   real_t Q_imp_plus[NUM_STATE];
   real_t Q_imp_minus[NUM_STATE];

   real_t t_sub_cursor = t_macro_start;

   const int O = static_cast<int>(deltaT_.size());
   for (int o = 0; o < O; ++o)
   {
      const real_t dt_sub      = deltaT_[o];
      const real_t weight      = time_weights_[o];
      const real_t t_sub_end   = t_sub_cursor + dt_sub;
      const real_t accum_scale = weight * dt_macro;

      // 1. Inject one nucleation increment at the sub-step endpoint.
      ApplyNucleationIncremental_TPV102(dof_data, fault_coords,
                                        t_sub_end, dt_sub);

      // 2. Per-QP friction pipeline + ψ update + imposed-state accumulator.
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
         flux_.ComputeStageState(d, Q_avg_plus, Q_avg_minus, s, method);

         // Per-sub-step slip accumulation in both fault-tangent components.
         d.slip1 += s.V1 * dt_sub;
         d.slip2 += s.V2 * dt_sub;

         // Aging-law analytic ψ update — exact for constant V over dt_sub.
         d.psi = UpdateStateAnalytic(d.psi, s.V_abs,
                                     d.Dc, dt_sub,
                                     state_evo_.GetF0(),
                                     d.b,  // Phase 11a: per-DOF b (was state_evo_.GetB())
                                     state_evo_.GetV0());

         flux_.BuildImposedState(d, s, Q_avg_plus, Q_avg_minus,
                                 Q_imp_plus, Q_imp_minus);

         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }

         if (last_sub_step)
         {
            flux_.WriteBackState(d, s);
         }
      }

      t_sub_cursor = t_sub_end;
   }
}

void Tpv102SubStepIterator::AdvanceWithSubStepStates(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
   real_t dt_macro,
   real_t t_macro_start,
   real_t *I_imp_plus_flat,
   real_t *I_imp_minus_flat,
   FrictionSolver::Method method)
{
   if (deltaT_.empty())
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: SetSubSteps "
         "has not been called; cannot iterate without a configured "
         "quadrature.");
   }
   if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: dt_macro "
         "must be finite and positive; got " + std::to_string(dt_macro));
   }
   if (!std::isfinite(t_macro_start))
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: t_macro_start "
         "must be finite; got " + std::to_string(t_macro_start));
   }
   if (I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: I_imp_*_flat "
         "must both be non-null.");
   }

   const int O = static_cast<int>(deltaT_.size());
   if (static_cast<int>(Q_pointwise_plus_per_substep.size()) != O ||
       static_cast<int>(Q_pointwise_minus_per_substep.size()) != O)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: Q_pointwise_*"
         "_per_substep must each have size O = " + std::to_string(O)
         + "; got plus = "
         + std::to_string(Q_pointwise_plus_per_substep.size())
         + ", minus = "
         + std::to_string(Q_pointwise_minus_per_substep.size()));
   }

   const int n = static_cast<int>(dof_data.size());
   if (static_cast<int>(fault_coords.size()) != n)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: dof_data and "
         "fault_coords must have equal size; got "
         + std::to_string(dof_data.size()) + ", "
         + std::to_string(fault_coords.size()));
   }

   const size_t expected_words =
      static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n);
   for (int o = 0; o < O; o++)
   {
      if (Q_pointwise_plus_per_substep[o].size() != expected_words ||
          Q_pointwise_minus_per_substep[o].size() != expected_words)
      {
         throw std::runtime_error(
            "Tpv102SubStepIterator::AdvanceWithSubStepStates: "
            "Q_pointwise_*[" + std::to_string(o) + "] must have size "
            "NUM_STATE * n = " + std::to_string(expected_words)
            + "; got plus = "
            + std::to_string(Q_pointwise_plus_per_substep[o].size())
            + ", minus = "
            + std::to_string(Q_pointwise_minus_per_substep[o].size()));
      }
   }

   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   const real_t rel = std::abs(dtsum - dt_macro) / std::max(dt_macro,
                                                            1e-300);
   const real_t sum_tol = std::max<real_t>(
      1e-12,
      static_cast<real_t>(10.0) * O
         * std::numeric_limits<real_t>::epsilon());
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: Σ deltaT "
         "must equal dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro));
   }

   const size_t nwords = expected_words;
   std::memset(I_imp_plus_flat,  0, nwords * sizeof(real_t));
   std::memset(I_imp_minus_flat, 0, nwords * sizeof(real_t));

   real_t Q_imp_plus[NUM_STATE];
   real_t Q_imp_minus[NUM_STATE];

   real_t t_sub_cursor = t_macro_start;

   for (int o = 0; o < O; ++o)
   {
      const real_t dt_sub      = deltaT_[o];
      const real_t weight      = time_weights_[o];
      const real_t t_sub_end   = t_sub_cursor + dt_sub;
      const real_t accum_scale = weight * dt_macro;

      ApplyNucleationIncremental_TPV102(dof_data, fault_coords,
                                        t_sub_end, dt_sub);

      const bool last_sub_step = (o == O - 1);
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
         flux_.ComputeStageState(d, Q_tilde_plus, Q_tilde_minus, s, method);

         d.slip1 += s.V1 * dt_sub;
         d.slip2 += s.V2 * dt_sub;

         d.psi = UpdateStateAnalytic(d.psi, s.V_abs,
                                     d.Dc, dt_sub,
                                     state_evo_.GetF0(),
                                     d.b,  // Phase 11a: per-DOF b (was state_evo_.GetB())
                                     state_evo_.GetV0());

         flux_.BuildImposedState(d, s, Q_tilde_plus, Q_tilde_minus,
                                 Q_imp_plus, Q_imp_minus);

         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }

         if (last_sub_step)
         {
            flux_.WriteBackState(d, s);
         }
      }

      t_sub_cursor = t_sub_end;
   }
}

// Nucleation-callback variant: byte-for-byte the plain
// AdvanceWithSubStepStates body above, EXCEPT the hard-coded
// ApplyNucleationIncremental_TPV102(...) call (once per sub-step, before
// the per-QP loop) is replaced by nuc_callback(t_sub_end, dt_sub).  Used
// by the SAFS spatial driver to route Gaussian gradual-overstress
// nucleation through the aging-law iterator.
void Tpv102SubStepIterator::AdvanceWithSubStepStates(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
   real_t dt_macro,
   real_t t_macro_start,
   real_t *I_imp_plus_flat,
   real_t *I_imp_minus_flat,
   FrictionSolver::Method method,
   const std::function<void(real_t, real_t)> &nuc_callback)
{
   // R-013: reject an empty callback up front — a default-constructed
   // std::function would throw std::bad_function_call inside the sub-step
   // loop.  Pass a no-op [](real_t,real_t){} to opt out of nucleation.
   if (!nuc_callback)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: nuc_callback is "
         "empty; pass a no-op [](real_t,real_t){} to opt out of nucleation.");
   }
   if (deltaT_.empty())
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: SetSubSteps "
         "has not been called; cannot iterate without a configured "
         "quadrature.");
   }
   if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: dt_macro "
         "must be finite and positive; got " + std::to_string(dt_macro));
   }
   if (!std::isfinite(t_macro_start))
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: t_macro_start "
         "must be finite; got " + std::to_string(t_macro_start));
   }
   if (I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: I_imp_*_flat "
         "must both be non-null.");
   }

   const int O = static_cast<int>(deltaT_.size());
   if (static_cast<int>(Q_pointwise_plus_per_substep.size()) != O ||
       static_cast<int>(Q_pointwise_minus_per_substep.size()) != O)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: Q_pointwise_*"
         "_per_substep must each have size O = " + std::to_string(O)
         + "; got plus = "
         + std::to_string(Q_pointwise_plus_per_substep.size())
         + ", minus = "
         + std::to_string(Q_pointwise_minus_per_substep.size()));
   }

   const int n = static_cast<int>(dof_data.size());
   if (static_cast<int>(fault_coords.size()) != n)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: dof_data and "
         "fault_coords must have equal size; got "
         + std::to_string(dof_data.size()) + ", "
         + std::to_string(fault_coords.size()));
   }

   const size_t expected_words =
      static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n);
   for (int o = 0; o < O; o++)
   {
      if (Q_pointwise_plus_per_substep[o].size() != expected_words ||
          Q_pointwise_minus_per_substep[o].size() != expected_words)
      {
         throw std::runtime_error(
            "Tpv102SubStepIterator::AdvanceWithSubStepStates: "
            "Q_pointwise_*[" + std::to_string(o) + "] must have size "
            "NUM_STATE * n = " + std::to_string(expected_words)
            + "; got plus = "
            + std::to_string(Q_pointwise_plus_per_substep[o].size())
            + ", minus = "
            + std::to_string(Q_pointwise_minus_per_substep[o].size()));
      }
   }

   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   const real_t rel = std::abs(dtsum - dt_macro) / std::max(dt_macro,
                                                            1e-300);
   const real_t sum_tol = std::max<real_t>(
      1e-12,
      static_cast<real_t>(10.0) * O
         * std::numeric_limits<real_t>::epsilon());
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv102SubStepIterator::AdvanceWithSubStepStates: Σ deltaT "
         "must equal dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro));
   }

   const size_t nwords = expected_words;
   std::memset(I_imp_plus_flat,  0, nwords * sizeof(real_t));
   std::memset(I_imp_minus_flat, 0, nwords * sizeof(real_t));

   real_t Q_imp_plus[NUM_STATE];
   real_t Q_imp_minus[NUM_STATE];

   real_t t_sub_cursor = t_macro_start;

   for (int o = 0; o < O; ++o)
   {
      const real_t dt_sub      = deltaT_[o];
      const real_t weight      = time_weights_[o];
      const real_t t_sub_end   = t_sub_cursor + dt_sub;
      const real_t accum_scale = weight * dt_macro;

      // Only line that differs from the plain overload: route nucleation
      // through the caller-supplied callback instead of the hard-coded
      // ApplyNucleationIncremental_TPV102(dof_data, fault_coords,
      // t_sub_end, dt_sub) at the identical point (once per sub-step,
      // before the per-QP loop).
      nuc_callback(t_sub_end, dt_sub);

      const bool last_sub_step = (o == O - 1);
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
         flux_.ComputeStageState(d, Q_tilde_plus, Q_tilde_minus, s, method);

         d.slip1 += s.V1 * dt_sub;
         d.slip2 += s.V2 * dt_sub;

         d.psi = UpdateStateAnalytic(d.psi, s.V_abs,
                                     d.Dc, dt_sub,
                                     state_evo_.GetF0(),
                                     d.b,  // Phase 11a: per-DOF b (was state_evo_.GetB())
                                     state_evo_.GetV0());

         flux_.BuildImposedState(d, s, Q_tilde_plus, Q_tilde_minus,
                                 Q_imp_plus, Q_imp_minus);

         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }

         if (last_sub_step)
         {
            flux_.WriteBackState(d, s);
         }
      }

      t_sub_cursor = t_sub_end;
   }
}

} // namespace seas
} // namespace mfem
