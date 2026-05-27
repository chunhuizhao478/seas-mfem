// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/friction_substep_iterator.hpp — Phase 5 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// The method-oriented unification of the three benchmark-named sub-step
// iterators (tpv205 LSW, tpv102 aging RS, tpv104 slip-law-SRW RS).  They
// collapse into:
//
//   SubStepIteratorBase                         <- common ADER envelope
//     RateStateSubStepIterator<RateStateAgingPolicy>      (= TPV102 aging)
//     RateStateSubStepIterator<RateStateSlipLawSrwPolicy> (= TPV104 SRW)
//     LinearSlipWeakeningIterator                          (= TPV205 LSW)
//
// All three implement the driver-facing `IFrictionIterator` strategy.  The
// standalone `tpv{205,102,104}_substep_iterator.*` are LEFT UNTOUCHED as
// the byte-exact oracle; `test_friction_substep_iterator_parity` proves the
// unified classes reproduce them bit-for-bit (DOFData + I_imp_±).
//
// === 2026-05-27 PLAN DEVIATION (user-approved; see §Phase 5 of the plan) ===
// The plan's §5.1 claimed the per-sub-step envelope is "byte-for-byte
// identical across all three" and that `RunSubSteps_` performs the slip
// accumulation.  The live oracle disagrees: LSW accumulates slip *inside*
// `StepOneQP_` (tpv205:129-130) and writes `slip_rate_substep_max` /
// `sigma_n_substep_min` unconditionally (tpv205:413/417), neither of which
// the RS path does.  So `RunSubSteps_` here owns ONLY the truly-common
// envelope (validation, memset, the sub-step loop + nuc_callback + running
// t_sub_end, the Q± slice, the I_imp += accum_scale·Q_imp accumulation),
// and a richer per-QP `step_fn` does ALL per-QP work — slip accumulation,
// the friction/state solve, the LSW-only diag writes, and WriteBackState.

#ifndef MFEM_SEAS_FRICTION_SUBSTEP_ITERATOR_HPP
#define MFEM_SEAS_FRICTION_SUBSTEP_ITERATOR_HPP

#include "mfem.hpp"

#include "fault_face_flux.hpp"               // DOFData, FaultFaceFlux, EvalStageState, NUM_STATE
#include "friction_solver.hpp"               // FrictionSolver::Method
#include "friction_iterator.hpp"             // IFrictionIterator, FaultFrictionLaw
#include "../friction/state_policies.hpp"    // RateStateAgingPolicy, RateStateSlipLawSrwPolicy

#include <cstddef>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Common ADER sub-step envelope shared by the rate-and-state and
/// linear-slip-weakening iterators.  Implements the `IFrictionIterator`
/// quadrature + diagnostic accessors; the per-QP `Advance` and the
/// wave-op-law tag remain pure for the method-specific derived classes.
class SubStepIteratorBase : public IFrictionIterator
{
public:
   /// Configure the ADER sub-step quadrature.  Byte-identical validator to
   /// the standalone iterators (`tpv{205,102,104}SubStepIterator::SetSubSteps`).
   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights) override;

   const std::vector<real_t> &GetDeltaT()      const override { return deltaT_; }
   const std::vector<real_t> &GetTimeWeights() const override { return time_weights_; }

   void SetDiagNumLocalFaultQPs(int n) override { diag_num_local_fault_qps_ = n; }

protected:
   /// The common per-macro-step envelope.  `step_fn` does ALL per-QP work
   /// (slip accumulation, friction/state solve, diag writes, WriteBackState)
   /// and FILLS `Q_imp_p` / `Q_imp_m` for the current QP; this routine owns
   /// the validation, the `memset`, the sub-step loop with `nuc_callback`
   /// and the running `t_sub_end`, the per-QP `Q±` slice, and the
   /// `I_imp += accum_scale·Q_imp` accumulation.
   ///
   /// `step_fn` signature:
   ///   void step_fn(int o, int O, int i, DOFData& d,
   ///                const real_t* Qp, const real_t* Qm,
   ///                real_t dt_sub, real_t t_sub_end, bool last_sub_step,
   ///                real_t* Q_imp_p, real_t* Q_imp_m);
   ///
   /// Header-defined (template); instantiated per call site in the .cpp.
   template <class StepFn>
   void RunSubSteps_(
      const char *who,
      std::vector<DOFData> &dof_data,
      const std::vector<Vector> &fault_coords,
      const std::vector<std::vector<real_t>> &Qp_per_substep,
      const std::vector<std::vector<real_t>> &Qm_per_substep,
      real_t dt_macro,
      real_t t_macro_start,
      real_t *I_imp_plus_flat,
      real_t *I_imp_minus_flat,
      const std::function<void(real_t, real_t)> &nuc_callback,
      StepFn &&step_fn)
   {
      const std::string who_s(who);

      // R-013/R-016: reject an empty callback up front — a default-constructed
      // std::function throws std::bad_function_call inside the loop.  Pass a
      // no-op [](real_t,real_t){} to opt out of nucleation.
      if (!nuc_callback)
      {
         throw std::runtime_error(
            who_s + ": nuc_callback is empty; pass a no-op "
            "[](real_t,real_t){} to opt out of nucleation.");
      }
      if (deltaT_.empty())
      {
         throw std::runtime_error(
            who_s + ": SetSubSteps has not been called; cannot iterate "
            "without a configured quadrature.");
      }
      if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
      {
         throw std::runtime_error(
            who_s + ": dt_macro must be finite and positive; got "
            + std::to_string(dt_macro));
      }
      if (!std::isfinite(t_macro_start))
      {
         throw std::runtime_error(
            who_s + ": t_macro_start must be finite; got "
            + std::to_string(t_macro_start));
      }
      if (I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
      {
         throw std::runtime_error(
            who_s + ": I_imp_*_flat must both be non-null.");
      }

      const int O = static_cast<int>(deltaT_.size());
      if (static_cast<int>(Qp_per_substep.size()) != O ||
          static_cast<int>(Qm_per_substep.size()) != O)
      {
         throw std::runtime_error(
            who_s + ": Q_pointwise_* must each have size O = "
            + std::to_string(O));
      }

      const int n = static_cast<int>(dof_data.size());
      if (static_cast<int>(fault_coords.size()) != n)
      {
         throw std::runtime_error(
            who_s + ": dof_data and fault_coords must have equal size; got "
            + std::to_string(dof_data.size()) + ", "
            + std::to_string(fault_coords.size()));
      }

      const size_t expected_words =
         static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n);
      for (int o = 0; o < O; ++o)
      {
         if (Qp_per_substep[o].size() != expected_words ||
             Qm_per_substep[o].size() != expected_words)
         {
            throw std::runtime_error(
               who_s + ": Q_pointwise_*[" + std::to_string(o)
               + "] must have size NUM_STATE * n = "
               + std::to_string(expected_words));
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
            who_s + ": Σ deltaT must equal dt_macro within "
            + std::to_string(sum_tol) + " rel; got Σ = "
            + std::to_string(dtsum) + ", dt_macro = "
            + std::to_string(dt_macro));
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
         const bool   last_sub_step = (o == O - 1);

         // Nucleation fires once per sub-step BEFORE the per-QP loop, at the
         // absolute sub-step end time (matches tpv102:326 / tpv205:390).
         nuc_callback(t_sub_end, dt_sub);

         const real_t *Qp_o = Qp_per_substep[o].data();
         const real_t *Qm_o = Qm_per_substep[o].data();

         for (int i = 0; i < n; ++i)
         {
            DOFData &d = dof_data[i];
            const real_t *Qp_i = Qp_o + static_cast<ptrdiff_t>(i) * NUM_STATE;
            const real_t *Qm_i = Qm_o + static_cast<ptrdiff_t>(i) * NUM_STATE;

            step_fn(o, O, i, d, Qp_i, Qm_i, dt_sub, t_sub_end, last_sub_step,
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
         }

         t_sub_cursor = t_sub_end;
      }
   }

   std::vector<real_t> deltaT_;
   std::vector<real_t> time_weights_;
   int                 diag_num_local_fault_qps_ = -1;
};

/// @brief Unified rate-and-state sub-step iterator, parameterized by a
/// state-evolution policy (`RateStateAgingPolicy` = TPV102/SAFS aging;
/// `RateStateSlipLawSrwPolicy` = TPV104 slip-law-SRW).  The policy's
/// `UpdatePsi` is the single point of variation; everything else is shared.
///
/// Owns the state-evolution law BY VALUE (a documented deviation from the
/// plan's `const Law&`): self-contained, no dangling-reference footgun,
/// matches the Phase-2 adapter's ownership, and bit-equivalent (the law is a
/// scalar holder whose getters feed the analytic ψ-update free functions).
template <class StatePolicy>
class RateStateSubStepIterator : public SubStepIteratorBase
{
public:
   using Law   = typename StatePolicy::Law;
   using Extra = typename StatePolicy::Extra;

   RateStateSubStepIterator(
      FaultFaceFlux &flux, Law law,
      FrictionSolver::Method method = FrictionSolver::Method::Brent,
      Extra extra = Extra{})
      : flux_(flux), law_(std::move(law)), method_(method), extra_(extra) {}

   // Held only via unique_ptr (flux_ is a reference); forbid copy/move so a
   // member-wise copy cannot dangle (mirrors R-015 on the Phase-2 adapter).
   RateStateSubStepIterator(const RateStateSubStepIterator &) = delete;
   RateStateSubStepIterator &operator=(const RateStateSubStepIterator &) = delete;
   RateStateSubStepIterator(RateStateSubStepIterator &&) = delete;
   RateStateSubStepIterator &operator=(RateStateSubStepIterator &&) = delete;

   void Advance(std::vector<DOFData> &dof_data,
                const std::vector<Vector> &fault_coords,
                const std::vector<std::vector<real_t>> &Q_pointwise_plus,
                const std::vector<std::vector<real_t>> &Q_pointwise_minus,
                real_t dt_macro,
                real_t t_macro_start,
                real_t *I_imp_plus_flat,
                real_t *I_imp_minus_flat,
                const std::function<void(real_t, real_t)> &nuc_callback)
      override;

   FaultFrictionLaw WaveOpLaw() const override
   { return FaultFrictionLaw::RateAndState; }

private:
   FaultFaceFlux          &flux_;
   Law                     law_;     ///< owned by value (see class doc)
   FrictionSolver::Method  method_;
   Extra                   extra_;   ///< V_w side-channel (SRW); unused for aging
};

/// Method-oriented aliases (driver/factory-facing names).
using RateStateAgingIterator      = RateStateSubStepIterator<RateStateAgingPolicy>;
using RateStateSlipLawSrwIterator = RateStateSubStepIterator<RateStateSlipLawSrwPolicy>;

/// @brief Unified linear-slip-weakening sub-step iterator (= TPV205).  A
/// separate class (no friction-solver method, no ComputeStageState, no state
/// variable); its per-QP `step_fn` lifts `Tpv205SubStepIterator::StepOneQP_`
/// verbatim and adds the LSW-only `slip_rate_substep_max` /
/// `sigma_n_substep_min` diag writes + the env-gated `[SLIP]` trace.
class LinearSlipWeakeningIterator : public SubStepIteratorBase
{
public:
   explicit LinearSlipWeakeningIterator(FaultFaceFlux &flux) : flux_(flux) {}

   LinearSlipWeakeningIterator(const LinearSlipWeakeningIterator &) = delete;
   LinearSlipWeakeningIterator &operator=(const LinearSlipWeakeningIterator &) = delete;
   LinearSlipWeakeningIterator(LinearSlipWeakeningIterator &&) = delete;
   LinearSlipWeakeningIterator &operator=(LinearSlipWeakeningIterator &&) = delete;

   void Advance(std::vector<DOFData> &dof_data,
                const std::vector<Vector> &fault_coords,
                const std::vector<std::vector<real_t>> &Q_pointwise_plus,
                const std::vector<std::vector<real_t>> &Q_pointwise_minus,
                real_t dt_macro,
                real_t t_macro_start,
                real_t *I_imp_plus_flat,
                real_t *I_imp_minus_flat,
                const std::function<void(real_t, real_t)> &nuc_callback)
      override;

   FaultFrictionLaw WaveOpLaw() const override
   { return FaultFrictionLaw::LSW; }

private:
   /// Per-QP closed-form LSW pipeline — lifted verbatim from
   /// `Tpv205SubStepIterator::StepOneQP_` (tpv205_substep_iterator.cpp:75-142):
   /// trial traction -> total traction -> LSW μ(δ) -> SolveLSW_TPV205 ->
   /// sigma_n_corr -> slip accumulation -> BuildImposedState ->
   /// WriteBackState on the last sub-step.
   void StepOneQP_(DOFData &d,
                   const real_t *Q_plus,
                   const real_t *Q_minus,
                   real_t dt_sub,
                   bool last_sub_step,
                   EvalStageState &s,
                   real_t *Q_imp_plus,
                   real_t *Q_imp_minus);

   FaultFaceFlux &flux_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FRICTION_SUBSTEP_ITERATOR_HPP
