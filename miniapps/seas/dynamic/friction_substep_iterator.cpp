// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Implementation of the unified method-oriented sub-step iterators
// (Phase 5).  See friction_substep_iterator.hpp for the design + the
// 2026-05-27 RunSubSteps_/step_fn deviation note.
//
// The per-QP step bodies are verbatim lifts of the standalone oracle
// iterators' per-QP work:
//   - RS  step_fn  <- tpv102_substep_iterator.cpp:333-371 (aging) and
//                     tpv104_substep_iterator.cpp (SRW), modulo the policy
//                     ψ-update line and the #ifdef SEAS_DIAG_TPV104_STATE
//                     probes (pure file I/O — omitted; they mutate no
//                     DOFData and so do not affect parity).
//   - LSW StepOneQP_ <- tpv205_substep_iterator.cpp:75-142, plus the LSW
//                     envelope's unconditional diag writes (tpv205:413/417)
//                     and the env-gated [SLIP] trace (tpv205:419-456).

#include "friction_substep_iterator.hpp"

#include "tpv205_friction.hpp"   // LSWFrictionCoefficient_TPV205, SolveLSW_TPV205
#include "wave_state.hpp"        // QIndex (VX, SXX) for the [SLIP] decomposition
#include "fault_resample.hpp"    // Phase 3: ApplyFaultResample (rate-state Δψ resample)

#include <algorithm>
#include <cmath>
#include <cstdio>     // [SLIP] diagnostic trace (env-gated)
#include <cstdlib>    // std::getenv / std::atof for the trace gates
#include <numeric>
#include <stdexcept>
#include <string>
#include <limits>

namespace mfem
{
namespace seas
{

// ---------------------------------------------------------------------------
// SubStepIteratorBase::SetSubSteps — byte-identical validator to the
// standalone iterators (tpv102_substep_iterator.cpp:28-74, etc.).
// ---------------------------------------------------------------------------
void SubStepIteratorBase::SetSubSteps(std::vector<real_t> deltaT,
                                      std::vector<real_t> time_weights)
{
   if (deltaT.empty())
   {
      throw std::runtime_error(
         "SubStepIteratorBase::SetSubSteps: deltaT must have at least "
         "one entry.");
   }
   if (deltaT.size() != time_weights.size())
   {
      throw std::runtime_error(
         "SubStepIteratorBase::SetSubSteps: deltaT and time_weights "
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
            "SubStepIteratorBase::SetSubSteps: deltaT[" + std::to_string(o)
            + "] must be finite and positive; got " + std::to_string(deltaT[o]));
      }
      if (!std::isfinite(time_weights[o]))
      {
         throw std::runtime_error(
            "SubStepIteratorBase::SetSubSteps: time_weights["
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
         "SubStepIteratorBase::SetSubSteps: time_weights must sum to "
         "1 within 1e-12; got sum = " + std::to_string(wsum));
   }

   deltaT_       = std::move(deltaT);
   time_weights_ = std::move(time_weights);
}

// ---------------------------------------------------------------------------
// RateStateSubStepIterator<StatePolicy>::Advance
//
// Per-QP step (reproduces tpv102/tpv104 AdvanceWithSubStepStates body):
//   ComputeStageState(method) -> slip accumulation -> ψ-update (policy)
//   -> BuildImposedState -> WriteBackState on the last sub-step.
// ---------------------------------------------------------------------------
template <class StatePolicy>
void RateStateSubStepIterator<StatePolicy>::Advance(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus,
   real_t dt_macro,
   real_t t_macro_start,
   real_t *I_imp_plus_flat,
   real_t *I_imp_minus_flat,
   const std::function<void(real_t, real_t)> &nuc_callback)
{
   // Validate the per-QP side-channel (V_w for SRW; no-op for aging) before
   // the loop dereferences extra_[i].
   StatePolicy::ValidateExtra(extra_, static_cast<int>(dof_data.size()));

   // Phase 3 (fault-dealiasing §6): rate-state resamples the per-MACRO-step
   // state-variable increment Δψ onto the degree-N space.  Snapshot ψ before
   // the sub-steps; the sub-step friction solve below runs EXACTLY as today
   // (it produces I_imp / τ_corr from the un-resampled rate — no resampled
   // quantity feeds the flux), then the net Δψ over the macro step is projected
   // per face and re-applied.  Inactive (resample disabled, or R == I because
   // there is no over-integration) ⇒ byte-exact with today.
   const bool do_resample = resample_enabled_ && resample_R_ != nullptr
                            && resample_nbf_per_face_ > 0 && !dof_data.empty();
   std::vector<real_t> psi_before;
   if (do_resample)
   {
      psi_before.resize(dof_data.size());
      for (size_t i = 0; i < dof_data.size(); ++i)
      { psi_before[i] = dof_data[i].psi; }
   }

   RunSubSteps_(
      "RateStateSubStepIterator::Advance",
      dof_data, fault_coords, Q_pointwise_plus, Q_pointwise_minus,
      dt_macro, t_macro_start, I_imp_plus_flat, I_imp_minus_flat, nuc_callback,
      [this](int /*o*/, int /*O*/, int i, DOFData &d,
             const real_t *Qp_i, const real_t *Qm_i,
             real_t dt_sub, real_t /*t_sub_end*/, bool last_sub_step,
             real_t *Q_imp_plus, real_t *Q_imp_minus)
      {
         EvalStageState s;
         flux_.ComputeStageState(d, Qp_i, Qm_i, s, method_);

         // Per-sub-step slip accumulation in both fault-tangent components.
         d.slip1 += s.V1 * dt_sub;
         d.slip2 += s.V2 * dt_sub;

         // The ONLY point of variation between aging and SRW.
         d.psi = StatePolicy::UpdatePsi(law_, d, s.V_abs, dt_sub, extra_, i);

         flux_.BuildImposedState(d, s, Qp_i, Qm_i, Q_imp_plus, Q_imp_minus);

         if (last_sub_step)
         {
            flux_.WriteBackState(d, s);
         }
      });

   // Phase 3: project the net Δψ over this macro step onto the degree-N space,
   // per fault face, and re-apply it.  R = I (no over-integration) makes this a
   // no-op; do_resample == false skips it entirely — both byte-exact.  The
   // per-face QP blocks are contiguous (interior faces then shared faces, each
   // nbf QPs in the rule's order — matching the R built from that same rule).
   if (do_resample)
   {
      const DenseMatrix &R = *resample_R_;
      const int nbf  = resample_nbf_per_face_;
      const int ndof = static_cast<int>(dof_data.size());
      MFEM_VERIFY(R.Height() == nbf && R.Width() == nbf,
                  "RateStateSubStepIterator::Advance: resample R is "
                  << R.Height() << "x" << R.Width() << ", expected "
                  << nbf << "x" << nbf << " (nbf_per_face).");
      MFEM_VERIFY(ndof % nbf == 0,
                  "RateStateSubStepIterator::Advance: fault QP count " << ndof
                  << " is not a multiple of nbf_per_face " << nbf
                  << " — per-face QP blocks are not contiguous.");
      const int nfaces = ndof / nbf;
      std::vector<real_t> dpsi(nbf), proj(nbf);
      for (int f = 0; f < nfaces; ++f)
      {
         const int b = f * nbf;
         for (int q = 0; q < nbf; ++q)
         { dpsi[q] = dof_data[b + q].psi - psi_before[b + q]; }
         ApplyFaultResample(R, dpsi.data(), proj.data());   // proj = R · Δψ
         for (int q = 0; q < nbf; ++q)
         { dof_data[b + q].psi = psi_before[b + q] + proj[q]; }
      }
   }
}

// Explicit instantiations — both rate-and-state policies.
template class RateStateSubStepIterator<RateStateAgingPolicy>;
template class RateStateSubStepIterator<RateStateSlipLawSrwPolicy>;

// ---------------------------------------------------------------------------
// LinearSlipWeakeningIterator::StepOneQP_ — verbatim lift of
// Tpv205SubStepIterator::StepOneQP_ (tpv205_substep_iterator.cpp:75-142).
// ---------------------------------------------------------------------------
void LinearSlipWeakeningIterator::StepOneQP_(DOFData &d,
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
   s.Theta         = std::sqrt(s.tau1_total * s.tau1_total
                               + s.tau2_total * s.tau2_total);

   // Step 3: LSW μ(δ) at current slip magnitude δ = sqrt(slip1² + slip2²).
   const real_t delta = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
   const real_t mu_eff = LSWFrictionCoefficient_TPV205(delta,
                                                       d.lsw_mu_s,
                                                       d.lsw_mu_d,
                                                       d.lsw_d_c);

   // Step 4: closed-form LSW solve.
   SolveLSW_TPV205(s.tau1_trial, s.tau2_trial,
                   s.tau1_total, s.tau2_total,
                   s.sigma_n_total, d.eta_s,
                   mu_eff,
                   s.V_abs, s.V1, s.V2,
                   s.tau1_corr, s.tau2_corr,
                   flux_.SigmaNStrengthFloorForLSW(),
                   d.lsw_cohesion);   // Phase 10 (TPV31) additive C0 (0 ⇒ byte-exact)

   s.sigma_n_corr = s.sigma_n_trial;

   // Step 5: per-sub-step slip accumulation.
   d.slip1 += s.V1 * dt_sub;
   d.slip2 += s.V2 * dt_sub;

   // Step 6: imposed state.
   flux_.BuildImposedState(d, s, Q_plus, Q_minus, Q_imp_plus, Q_imp_minus);

   // Step 7: on the final sub-step, write back V/slip_rate/tau*_corr/sigma_n_corr.
   if (last_sub_step)
   {
      flux_.WriteBackState(d, s);
   }
}

// ---------------------------------------------------------------------------
// LinearSlipWeakeningIterator::Advance
// ---------------------------------------------------------------------------
void LinearSlipWeakeningIterator::Advance(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus,
   real_t dt_macro,
   real_t t_macro_start,
   real_t *I_imp_plus_flat,
   real_t *I_imp_minus_flat,
   const std::function<void(real_t, real_t)> &nuc_callback)
{
   // PLAN_speckle_slip_runaway Phase 1: per-sub-step diagnostic trace at
   // spiking nodes (runtime env gate; zero overhead and byte-exact when
   // SEAS_DIAG_SLIP is unset).  Parsed once — lifted from
   // tpv205_substep_iterator.cpp:363-377.
   static const bool slip_diag = [] {
      const char *e = std::getenv("SEAS_DIAG_SLIP");
      return e && e[0] && e[0] != '0'; }();
   static const double slip_v_thr = [] {
      const char *e = std::getenv("SEAS_DIAG_SLIP_VTHR");
      return e ? std::atof(e) : 10.0; }();
   static const int s_mpi_rank = [] {
      int r = 0;
#ifdef MFEM_USE_MPI
      int inited = 0; MPI_Initialized(&inited);
      if (inited) { MPI_Comm_rank(MPI_COMM_WORLD, &r); }
#endif
      return r; }();

   RunSubSteps_(
      "LinearSlipWeakeningIterator::Advance",
      dof_data, fault_coords, Q_pointwise_plus, Q_pointwise_minus,
      dt_macro, t_macro_start, I_imp_plus_flat, I_imp_minus_flat, nuc_callback,
      [this, &fault_coords](int o, int O, int i, DOFData &d,
                            const real_t *Qp_i, const real_t *Qm_i,
                            real_t dt_sub, real_t t_sub_end, bool last_sub_step,
                            real_t *Q_imp_plus, real_t *Q_imp_minus)
      {
         EvalStageState s;
         StepOneQP_(d, Qp_i, Qm_i, dt_sub, last_sub_step, s,
                    Q_imp_plus, Q_imp_minus);

         // LSW-only unconditional diag writes (tpv205:413/417): honest
         // per-macro-step sub-step |V| max and most-tensile σ_n.  Reset to
         // 0 / +inf per macro step by the driver before Advance.
         d.slip_rate_substep_max = std::max(d.slip_rate_substep_max, s.V_abs);
         d.sigma_n_substep_min   = std::min(d.sigma_n_substep_min, s.sigma_n_total);

         // [SLIP] trace — env-gated, mutates no state (tpv205:419-456).
         if (slip_diag && s.V_abs > slip_v_thr)
         {
            const real_t delta  = std::sqrt(d.slip1 * d.slip1
                                            + d.slip2 * d.slip2);
            const real_t mu_eff = LSWFrictionCoefficient_TPV205(
                                     delta, d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
            const real_t sn_pos = std::max<real_t>(s.sigma_n_total, 0.0);
            const real_t sn_vjump = d.eta_p * (Qm_i[VX] - Qp_i[VX]);
            // Precondition: d.Zp_plus / d.Zp_minus are P-wave impedances
            // (ρ·c_p > 0 for any physical material), so these divisions are
            // safe.  (R-004) This block is env-gated diagnostic-only and
            // mutates no DOFData; it never runs on the byte-exact path.
            const real_t sn_sterm = d.eta_p
                                    * (Qp_i[SXX] / d.Zp_plus
                                       + Qm_i[SXX] / d.Zp_minus);
            const int is_shared = (diag_num_local_fault_qps_ < 0) ? -1
                                  : (i >= diag_num_local_fault_qps_ ? 1 : 0);
            const Vector &xyz = fault_coords[i];
            std::fprintf(stderr,
               "[SLIP] t=%.6e o=%d/%d qp=%d last=%d V_abs=%+.6e tau_abs=%+.6e "
               "tau_str=%+.6e sigma_n_tot=%+.6e sigma_n_pos=%+.6e delta=%+.6e "
               "mu_eff=%+.6e d_c=%+.6e rank=%d is_shared=%d sn_vjump=%+.6e "
               "sn_sterm=%+.6e c=(%.1f,%.1f,%.1f)\n",
               t_sub_end, o, O, i, last_sub_step ? 1 : 0,
               s.V_abs, s.Theta, mu_eff * sn_pos, s.sigma_n_total, sn_pos,
               delta, mu_eff, d.lsw_d_c,
               s_mpi_rank, is_shared, sn_vjump, sn_sterm,
               xyz(0), xyz(1), xyz(2));
            std::fflush(stderr);
         }
      });
}

} // namespace seas
} // namespace mfem
