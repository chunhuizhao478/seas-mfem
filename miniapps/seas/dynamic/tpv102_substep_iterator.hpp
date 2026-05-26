// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// TPV102 sub-step fault iterator — per-sub-step friction + ψ + nucleation
// orchestration for the SCEC TPV102 (rate-and-state, ageing law) benchmark.
//
// Mirrors `dynamic/tpv104_substep_iterator.hpp` with these substitutions:
//   - bound state-evolution law is `AgingLawPsi` (not `SlipLawSRWPsi`);
//   - the API drops the `V_w` per-QP side-channel (TPV102 has no
//     weakening-velocity field; the SRW slip-law-with-strong-rate-weakening
//     is replaced by the regularized rate-and-state ageing law);
//   - the analytic ψ update inside the per-sub-step loop calls
//     `UpdateStateAnalytic(psi, V, Dc, dt, f0, b, V0)` (the exact aging-law
//     ψ update from `friction/state_evolution.hpp`), not
//     `UpdateStateAnalyticSlipLawSRW(...)`;
//   - the nucleation accumulator routes through
//     `ApplyNucleationIncremental_TPV102(...)`.
//
// All other contracts match the TPV104 iterator literally:
//   1. Caller supplies a macro-step ADER predictor as time-integrated
//      I_± = ∫_{t_start}^{t_start + dt_macro} Q_±(τ) dτ.
//   2. For every ADER sub-step o in [0, O):
//        - apply one nucleation increment at the sub-step endpoint,
//        - run the friction pipeline `ComputeStageState` on Q̄ = I/dt,
//        - analytic ψ update for size deltaT[o],
//        - build per-sub-step imposed state,
//        - accumulate into I_imp_±_out weighted by timeWeights[o]*dt_macro.
//   3. After the final sub-step, `WriteBackState` is called once.
//
// In the O = 1 limit (`deltaT = {dt_macro}, time_weights = {1.0}`) this
// reduces bit-identically to the legacy TPV102 `EvaluateADER` path.
//
// Files NOT touched:
//   - `dynamic/fault_face_flux.{hpp,cpp}` — shared with BP5 [C2].
//   - `dynamic/wave_operator.{hpp,inl}` — already TPV102/TPV104-shared.
//   - `friction/state_evolution.hpp` — only an additive `GetB/V0/F0`
//     accessor patch was applied to `AgingLawPsi`, behaviour unchanged.

#ifndef MFEM_SEAS_TPV102_SUBSTEP_ITERATOR_HPP
#define MFEM_SEAS_TPV102_SUBSTEP_ITERATOR_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"
#include "friction_solver.hpp"
#include "tpv102_friction_solver.hpp"
#include "wave_state.hpp"
#include "tpv102_nucleation.hpp"
#include "../friction/state_evolution.hpp"

#include <functional>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Per-sub-step fault iterator for TPV102.
///
/// Holds non-owning references to a `FaultFaceFlux` (for the friction
/// pipeline helpers) and an `AgingLawPsi` (for global friction scalars
/// `b, V0, f0` and the steady-state / rate routines).  Per-call inputs are:
///   - `dof_data`         (mutable, per-QP state)
///   - `fault_coords`     (per-QP physical coords)
///   - `I_plus / I_minus` (time-integrated predictor, flat NUM_STATE*ndof)
///   - macro-step metadata (dt_macro, t_start)
///   - `I_imp_plus_out / I_imp_minus_out`  (time-integrated output)
///
/// Sub-step quadrature `(deltaT, timeWeights)` is configured once via
/// `SetSubSteps` and reused per macro-step.
class Tpv102SubStepIterator
{
public:
   /// Bind the iterator to a flux solver and a state-evolution law.  The
   /// references MUST outlive the iterator.
   Tpv102SubStepIterator(FaultFaceFlux &flux,
                         const AgingLawPsi &state_evo);

   /// Configure the ADER sub-step quadrature.  Requirements:
   ///   - `deltaT.size() == time_weights.size()`
   ///   - `deltaT.size() >= 1`
   ///   - every `deltaT[o] > 0`
   ///   - `Σ deltaT[o] == dt_macro` (verified per-call)
   ///   - `Σ time_weights[o] == 1` (verified here)
   /// Throws `std::runtime_error` on violation.
   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights);

   /// Number of configured sub-steps.
   int NumSubSteps() const
   { return static_cast<int>(deltaT_.size()); }

   /// @brief Advance every fault QP through one macro-step.
   ///
   /// Mutates `dof_data` (psi, slip_rate, V1/V2, slip1/slip2, tau*_corr,
   /// sigma_n_corr, tau2_nuc).  Writes the accumulated time-integrated
   /// imposed state into `I_imp_plus_out / I_imp_minus_out`.
   ///
   /// Friction-solve dispatch is identical to TPV104 — `Method` enum
   /// exposes Brent, NewtonRaphson, NewtonRaphsonStable, HybridNRBisection.
   /// Default is `NewtonRaphsonStable` (stable-asinh Newton; the same
   /// Newton solver `tpv104_friction_solver.hpp` provides — the TPV102
   /// regularised rate-and-state friction COEFFICIENT has the same
   /// functional form `μ = a · arcsinh[V/(2V0)·exp(ψ/a)]` as TPV104,
   /// so the per-QP solver is law-agnostic).
   void Advance(std::vector<DOFData> &dof_data,
                const std::vector<Vector> &fault_coords,
                const real_t *I_plus_flat,
                const real_t *I_minus_flat,
                real_t dt_macro,
                real_t t_macro_start,
                real_t *I_imp_plus_flat,
                real_t *I_imp_minus_flat,
                FrictionSolver::Method method
                   = FrictionSolver::Method::NewtonRaphsonStable);

   /// SeisSol-equivalent variant: takes per-sub-step pointwise Q at fault
   /// QPs in the canonical fault-local frame.  At each sub-step `o`, the
   /// friction pipeline runs on `Q_pointwise_*_per_substep[o]` directly.
   void AdvanceWithSubStepStates(
      std::vector<DOFData> &dof_data,
      const std::vector<Vector> &fault_coords,
      const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
      const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
      real_t dt_macro,
      real_t t_macro_start,
      real_t *I_imp_plus_flat,
      real_t *I_imp_minus_flat,
      FrictionSolver::Method method
         = FrictionSolver::Method::NewtonRaphsonStable);

   /// @brief Nucleation-callback variant of `AdvanceWithSubStepStates`.
   ///
   /// Identical to the plain per-sub-step overload above EXCEPT the
   /// hard-coded `ApplyNucleationIncremental_TPV102(...)` call (once per
   /// ADER sub-step, before the per-QP loop) is replaced by
   /// `nuc_callback(t_substep_end, dt_substep)`.  This lets the SAFS
   /// spatial driver route its own (Gaussian gradual-overstress)
   /// nucleation through the aging-law iterator without hard-coding the
   /// TPV102 patch.  Mirrors `Tpv205SubStepIterator`'s callback overload.
   ///
   /// The callback MUST mutate only `DOFData::tau1_nuc / tau2_nuc /
   /// sigma_n_nuc`.  Pass `[](real_t, real_t){}` to opt out (no-op) —
   /// then the per-QP friction pipeline is byte-identical to a run with
   /// nucleation suppressed.  `method` has NO default here (the trailing
   /// `nuc_callback` is required), so callers must pass both; the SAFS
   /// driver passes `FrictionSolver::Method::Brent` (CLAUDE.md).
   ///
   /// The plain overload above is left untouched (byte-identical) for
   /// the standalone TPV102 driver; this overload is strictly additive.
   void AdvanceWithSubStepStates(
      std::vector<DOFData> &dof_data,
      const std::vector<Vector> &fault_coords,
      const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
      const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
      real_t dt_macro,
      real_t t_macro_start,
      real_t *I_imp_plus_flat,
      real_t *I_imp_minus_flat,
      FrictionSolver::Method method,
      const std::function<void(real_t /*t_substep_end*/,
                               real_t /*dt_substep*/)> &nuc_callback);

   /// Accessor for the configured sub-step sizes (test hook).
   const std::vector<real_t> &GetDeltaT() const { return deltaT_; }

   /// Accessor for the configured quadrature weights (test hook).
   const std::vector<real_t> &GetTimeWeights() const { return time_weights_; }

   /// Diagnostic-only accessor (mirrors `Tpv205SubStepIterator:146`).
   /// TPV102 has no `[SLIP]` is_shared trace, so this is a no-op store —
   /// it exists purely so the Phase-2 adapter can satisfy the
   /// `IFrictionIterator` interface uniformly across LSW and RS laws.
   void SetDiagNumLocalFaultQPs(int n) { diag_num_local_fault_qps_ = n; }

private:
   FaultFaceFlux        &flux_;
   const AgingLawPsi    &state_evo_;
   std::vector<real_t>   deltaT_;
   std::vector<real_t>   time_weights_;
   int                   diag_num_local_fault_qps_ = -1;  // interface symmetry (diag only)
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV102_SUBSTEP_ITERATOR_HPP
