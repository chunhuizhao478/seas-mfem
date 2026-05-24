// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// TPV205 sub-step fault iterator — per-sub-step LSW friction solve +
// per-component slip accumulation orchestration around `FaultFaceFlux`
// (mirrors `dynamic/tpv104_substep_iterator.hpp`'s API minus state
// evolution and nucleation).
//
// LSW (linear slip-weakening) doesn't have a state variable, so the
// per-sub-step pipeline collapses to:
//
//   for o in 0..O:
//      1. Trial traction via FaultFaceFlux::ComputeTrialTraction (pure).
//      2. Compute current slip magnitude δ = sqrt(slip1² + slip2²).
//      3. μ_eff(δ) = μ_s − (μ_s − μ_d) · min(δ/d_c, 1).
//      4. Closed-form V_abs = max(0, (|τ_total| − μ_eff·σ_n) / η_s).
//      5. Trial-scale corrected traction: τ*_corr = τ*_trial − η_s · V*.
//      6. Accumulate slip: slip{1,2} += V{1,2} · dt_sub.
//      7. Imposed state via FaultFaceFlux::BuildImposedState (pure).
//      8. Accumulate I_imp_± weighted by timeWeights[o]·dt_macro.
//   After the last sub-step, FaultFaceFlux::WriteBackState records
//   slip_rate / V1 / V2 / tau*_corr (TOTAL) / sigma_n_corr on DOFData.
//
// Nucleation is STATIC in TPV5 (patch pre-stress set at t=0; no
// time-varying perturbation), so the iterator never touches
// `tau*_nuc` / `sigma_n_nuc` — they remain at their init-time zeros.
//
// At O = 1 with `deltaT = {dt_macro}`, `time_weights = {1.0}` the
// iterator collapses to a single LSW solve + BuildImposedState.

#ifndef MFEM_SEAS_TPV205_SUBSTEP_ITERATOR_HPP
#define MFEM_SEAS_TPV205_SUBSTEP_ITERATOR_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"
#include "wave_state.hpp"

#include <functional>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Per-sub-step fault iterator for TPV205 (linear slip-weakening).
///
/// Holds a non-owning reference to a `FaultFaceFlux` (used only for the
/// pure helpers `ComputeTrialTraction`, `BuildImposedState`, and
/// `WriteBackState`).  Per-call inputs:
///   - `dof_data`         (mutable, per-QP state — V, slip, tau*_corr)
///   - `fault_coords`     (per-QP physical coords; passed for symmetry
///                         with TPV104's API even though LSW doesn't
///                         currently consult them)
///   - `I_plus / I_minus` (time-integrated predictor, NUM_STATE * ndof)
///   - macro-step metadata (`dt_macro`, `t_macro_start`)
///   - `I_imp_*_out`      (accumulated imposed states; zeroed on entry)
///
/// Configure the sub-step quadrature once via `SetSubSteps`; the
/// per-call `Σ deltaT == dt_macro` check is RELATIVE so a fixed
/// configuration handles every macro-step under varying dt.
class Tpv205SubStepIterator
{
public:
   /// Bind the iterator to a flux solver.  The reference MUST outlive
   /// the iterator.
   explicit Tpv205SubStepIterator(FaultFaceFlux &flux);

   /// Configure the ADER sub-step quadrature.  See TPV104 iterator for
   /// the constraints (deltaT > 0, equal sizes, weights sum to 1).
   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights);

   /// Number of configured sub-steps.
   int NumSubSteps() const
   { return static_cast<int>(deltaT_.size()); }

   /// @brief Advance every fault QP through one macro-step using a
   /// single time-integrated predictor (Q̄ = I/dt_macro for friction).
   ///
   /// Mutates `dof_data` (slip_rate, V1/V2, slip1/slip2, tau*_corr,
   /// sigma_n_corr).  Writes the accumulated time-integrated imposed
   /// state into `I_imp_*_flat`.
   ///
   /// R-011: the production driver routes through
   /// `AdvanceADERWithSubStep` → `AdvanceWithSubStepStates`, NOT this
   /// method.  `Advance` is kept for API symmetry with the TPV104
   /// iterator and as the entry point for any future test that wants
   /// the time-averaged Q̄ semantics.  Refactors of `StepOneQP_` must
   /// keep this code path consistent with the per-sub-step variant —
   /// at O = 1 with `deltaT = {dt_macro}, weights = {1.0}` the two
   /// methods MUST produce bit-identical DOFData when `Q_pointwise[0]
   /// == I/dt_macro`.
   void Advance(std::vector<DOFData> &dof_data,
                const std::vector<Vector> &fault_coords,
                const real_t *I_plus_flat,
                const real_t *I_minus_flat,
                real_t dt_macro,
                real_t t_macro_start,
                real_t *I_imp_plus_flat,
                real_t *I_imp_minus_flat);

   /// @brief SeisSol-equivalent variant: per-sub-step pointwise Q
   /// supplied by the predictor.
   void AdvanceWithSubStepStates(
      std::vector<DOFData> &dof_data,
      const std::vector<Vector> &fault_coords,
      const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
      const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
      real_t dt_macro,
      real_t t_macro_start,
      real_t *I_imp_plus_flat,
      real_t *I_imp_minus_flat);

   /// @brief Phase N callback-aware overload of `AdvanceWithSubStepStates`.
   ///
   /// `nuc_callback(t_substep_end, dt_substep)` is invoked ONCE per
   /// ADER sub-step BEFORE the per-QP friction pipeline (mirrors the
   /// internal call to `ApplyNucleationIncremental_TPV104` at
   /// `tpv104_substep_iterator.cpp:295`).  The callback MUST mutate
   /// ONLY `DOFData::tau1_nuc / tau2_nuc / sigma_n_nuc` channels.
   /// Callers may pass `[](real_t, real_t){}` to opt out (no-op).
   ///
   /// Bit-equivalent to the no-callback overload when `nuc_callback`
   /// is a no-op lambda; the no-callback overload now routes through
   /// this method with a no-op callback (R-N-005 — keeps the original
   /// API working byte-for-byte while exposing the new hook).
   void AdvanceWithSubStepStates(
      std::vector<DOFData> &dof_data,
      const std::vector<Vector> &fault_coords,
      const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
      const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
      real_t dt_macro,
      real_t t_macro_start,
      real_t *I_imp_plus_flat,
      real_t *I_imp_minus_flat,
      const std::function<void(real_t, real_t)> &nuc_callback);

   /// Accessors for the configured quadrature (test hook).
   const std::vector<real_t> &GetDeltaT() const { return deltaT_; }
   const std::vector<real_t> &GetTimeWeights() const { return time_weights_; }

   /// Diagnostic-only (PLAN_speckle_slip_runaway): number of LOCAL/interior
   /// fault QPs.  The dof_data layout is interior QPs [0, n_local) followed by
   /// shared QPs [n_local, n_total), so the env-gated [SLIP] trace tags each QP
   /// is_shared = (i >= n_local).  Default -1 => unknown => trace prints
   /// is_shared=-1.  Set once by the driver from GetNumLocalFaultQPs().
   void SetDiagNumLocalFaultQPs(int n) { diag_num_local_fault_qps_ = n; }

private:
   /// Per-QP LSW solve + imposed-state construction shared between
   /// `Advance` and `AdvanceWithSubStepStates`.  The caller supplies the
   /// per-sub-step Q± (either Q̄ = I/dt_macro or Q̃[o] from the predictor)
   /// and the sub-step size.
   void StepOneQP_(DOFData &d,
                   const real_t *Q_plus, const real_t *Q_minus,
                   real_t dt_sub,
                   bool last_sub_step,
                   EvalStageState &s,
                   real_t *Q_imp_plus, real_t *Q_imp_minus);

   FaultFaceFlux        &flux_;
   std::vector<real_t>   deltaT_;
   std::vector<real_t>   time_weights_;
   int                   diag_num_local_fault_qps_ = -1;  // [SLIP] is_shared tag (diag only)
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV205_SUBSTEP_ITERATOR_HPP
