// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// TPV104 sub-step fault iterator — per-sub-step friction + ψ + nucleation
// orchestration around `FaultFaceFlux` (§4.10 Step 7 of tpv104 debug plan).
//
// R-007 (2026-04-24) Step 7 of TPV104 implementation.  Wraps (does NOT
// edit) `dynamic/fault_face_flux.{hpp,cpp}` and the ADER predictor in
// `dynamic/wave_operator.inl` — both are on the [C2] do-not-touch list
// because they are shared with TPV102 and BP5.  The per-sub-step loop
// lives here; `FaultFaceFlux::ComputeStageState`,
// `FaultFaceFlux::BuildImposedState`, and `FaultFaceFlux::WriteBackState`
// are called through as read-only public helpers.
//
// Contract summary (plan §4.10 Step 7):
//   1. Caller supplies a macro-step ADER predictor as time-integrated
//      I_± = ∫_{t_start}^{t_start + dt_macro} Q_±(τ) dτ.
//   2. For every ADER sub-step o in [0, O):
//        - apply one nucleation increment at the sub-step endpoint,
//        - run the friction pipeline `ComputeStageState` on Q̄ = I/dt,
//        - analytic ψ update for size deltaT[o],
//        - build per-sub-step imposed state,
//        - accumulate into I_imp_±_out weighted by timeWeights[o]*dt_macro.
//   3. After the final sub-step, `WriteBackState` is called once so
//      TPV102-shaped DOFData semantics (slip_rate, V1/V2, tau*_corr,
//      sigma_n_corr) are preserved for the probe format.
//
// Imposed-state accumulation (plan §3.8):
//   I_imp_± = Σ_o timeWeights[o] · dt_macro · Q_imp_±^{(o)}.
// In the O = 1 limit this reduces to the TPV102 `EvaluateADER` path
// (T_TPV104_SSI_3 bit-identity).

#ifndef MFEM_SEAS_TPV104_SUBSTEP_ITERATOR_HPP
#define MFEM_SEAS_TPV104_SUBSTEP_ITERATOR_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"
#include "friction_solver.hpp"
#include "tpv104_friction_solver.hpp"
#include "wave_state.hpp"
#include "tpv104_nucleation.hpp"
#include "../friction/slip_law_srw_psi.hpp"

#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Per-sub-step fault iterator for TPV104.
///
/// Holds non-owning references to a `FaultFaceFlux` (for the friction
/// pipeline helpers) and a `SlipLawSRWPsi` (for global friction scalars
/// `b, V0, f0, muW`).  Per-call inputs are:
///   - `dof_data`       (mutable, per-QP state)
///   - `fault_coords`   (per-QP physical coords)
///   - `V_w`            (driver-owned per-QP side-channel)
///   - `I_plus / I_minus`  (time-integrated predictor, flat NUM_STATE*ndof)
///   - macro-step metadata (dt_macro, t_start)
///   - `I_imp_plus_out / I_imp_minus_out`  (time-integrated output)
///
/// The sub-step quadrature `(deltaT, timeWeights)` is configured once
/// via `SetSubSteps` and reused for every `Advance` call — every ADER
/// macro-step uses the same ADER-O quadrature on this cluster.
class Tpv104SubStepIterator
{
public:
   /// Bind the iterator to a flux solver (for stage helpers) and a
   /// state-evolution law (for global friction scalars).  The references
   /// MUST outlive the iterator.
   Tpv104SubStepIterator(FaultFaceFlux &flux,
                         const SlipLawSRWPsi &state_evo);

   /// Configure the ADER sub-step quadrature.  Requirements:
   ///   - `deltaT.size() == time_weights.size()`
   ///   - `deltaT.size() >= 1`
   ///   - every `deltaT[o] > 0`
   ///   - `Σ deltaT[o] == dt_macro` (verified per-call)
   ///   - `Σ time_weights[o] == 1` (verified here)
   /// Throws `std::runtime_error` on violation.  The single-sub-step
   /// case `deltaT = {dt_macro}, time_weights = {1.0}` is the TPV102
   /// one-shot limit (T_TPV104_SSI_3).
   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights);

   /// Number of configured sub-steps.
   int NumSubSteps() const
   { return static_cast<int>(deltaT_.size()); }

   /// @brief Advance every fault QP through one macro-step.
   ///
   /// Mutates `dof_data` (psi, slip_rate, V1/V2, slip1/slip2, tau*_corr,
   /// sigma_n_corr, tau2_nuc).  `slip1` and `slip2` are accumulated
   /// per sub-step (R4-001 review fix) using the sub-step V1/V2; the
   /// final sub-step's V1/V2 is written back via WriteBackState.
   /// Writes the accumulated time-integrated imposed state into
   /// `I_imp_plus_out / I_imp_minus_out` (flat `NUM_STATE * ndof`).
   ///
   /// Friction-solve dispatch: the iterator calls
   /// `FaultFaceFlux::ComputeStageState(..., method)` which routes
   /// through `FrictionSolver::Solve(..., method)` for the per-QP
   /// friction solve.  The `Method` enum exposes FOUR choices:
   ///
   ///   - `Method::Brent`                — MFEM μ + log10-V Brent (legacy).
   ///   - `Method::NewtonRaphson`        — MFEM μ + Newton (legacy).
   ///   - `Method::NewtonRaphsonStable`  — stable-asinh μ + Newton
   ///                                       (Plan §4.10 Step 5, DEFAULT).
   ///   - `Method::HybridNRBisection`    — MFEM μ + NR-with-bisection.
   ///
   /// R5-003 (plan §4.10.X) resolution: stable-asinh Newton now lives
   /// inside `FrictionSolver::Solve` (as `SolveNRStable`), so the
   /// iterator does NOT need to bypass `ComputeStageState` — the same
   /// dispatch routes all four methods uniformly.  `fault_face_flux.cpp`
   /// remains untouched ([C2]), BP5/TPV102 legacy paths are preserved
   /// (they default to `Method::Brent` via `FrictionSolver::Solve`'s
   /// argument default).
   ///
   /// Preconditions validated:
   ///   - `SetSubSteps` called at least once.
   ///   - `dof_data.size() == fault_coords.size() == V_w.size()`.
   ///   - `dt_macro > 0`.
   ///   - All sub-step sizes positive.
   ///
   /// @param[in,out] dof_data           Per-DOF state (resized by caller).
   /// @param[in]     fault_coords       Per-DOF physical coords.
   /// @param[in]     V_w                Per-DOF weakening velocity
   ///                                   (side-channel from Step 3).
   /// @param[in]     I_plus_flat        Time-integrated + predictor
   ///                                   (`NUM_STATE * ndof`, per-DOF
   ///                                   stride `NUM_STATE`).
   /// @param[in]     I_minus_flat       Time-integrated − predictor.
   /// @param[in]     dt_macro           Macro-step size (> 0).
   /// @param[in]     t_macro_start      Absolute time at the macro-step
   ///                                   start (nucleation endpoints use
   ///                                   this + cumulative sub-step
   ///                                   offsets).
   /// @param[out]    I_imp_plus_flat    Accumulated time-integrated +
   ///                                   imposed state (zeroed on entry).
   /// @param[out]    I_imp_minus_flat   Accumulated time-integrated −
   ///                                   imposed state.
   /// @param[in]     method             Friction solver dispatch.
   ///                                   Default `NewtonRaphsonStable`
   ///                                   (Plan §4.10 Step 5 / R5-003).
   void Advance(std::vector<DOFData> &dof_data,
                const std::vector<Vector> &fault_coords,
                const std::vector<real_t> &V_w,
                const real_t *I_plus_flat,
                const real_t *I_minus_flat,
                real_t dt_macro,
                real_t t_macro_start,
                real_t *I_imp_plus_flat,
                real_t *I_imp_minus_flat,
                FrictionSolver::Method method
                   = FrictionSolver::Method::NewtonRaphsonStable);

   /// SeisSol-equivalent variant of `Advance`: takes per-sub-step pointwise
   /// Q at fault QPs in the canonical fault-local frame, instead of a single
   /// macro-step time-integrated I.  At each sub-step `o`, the friction
   /// pipeline runs on `Q_pointwise_*_per_substep[o]` directly — matching
   /// SeisSol's `qInterpolated[o]` → `precomputeStressFromQInterpolated[o]`
   /// pattern (closes the R4-004 cadence-deviation note).
   ///
   /// At O = 1 with `deltaT = {dt_macro}`, `time_weights = {1.0}`, and
   /// `Q_pointwise_*[0]` equal to `I_plus_flat / dt_macro`, this method
   /// produces bit-identical DOFData to `Advance` (T_TPV104_SSI_3 contract
   /// extended to the per-sub-step Q path).
   ///
   /// Layout: `Q_pointwise_plus_per_substep[o]` is a flat vector of size
   /// `NUM_STATE * num_fault_qps`; entry `[i * NUM_STATE + c]` is component
   /// c of Q at fault QP i on the canonical-+ side at sub-step time
   /// `t_macro_start + Σ_{o'<o} deltaT[o'] + tau_local[o]` (the predictor
   /// supplies whichever node convention the driver chose; the iterator
   /// does not enforce a specific node distribution).
   ///
   /// @param[in,out] dof_data           Per-DOF state.
   /// @param[in]     fault_coords       Per-DOF physical coords.
   /// @param[in]     V_w                Per-DOF weakening velocity.
   /// @param[in]     Q_pointwise_plus_per_substep   O Vectors, each
   ///                                   `NUM_STATE * n` in canonical frame.
   /// @param[in]     Q_pointwise_minus_per_substep  Same.
   /// @param[in]     dt_macro           Macro-step size (> 0).
   /// @param[in]     t_macro_start      Absolute start time.
   /// @param[out]    I_imp_plus_flat    Accumulated time-integrated +
   ///                                   imposed state (zeroed on entry).
   /// @param[out]    I_imp_minus_flat   Accumulated time-integrated −
   ///                                   imposed state.
   /// @param[in]     method             Friction solver dispatch.
   void AdvanceWithSubStepStates(
      std::vector<DOFData> &dof_data,
      const std::vector<Vector> &fault_coords,
      const std::vector<real_t> &V_w,
      const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
      const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
      real_t dt_macro,
      real_t t_macro_start,
      real_t *I_imp_plus_flat,
      real_t *I_imp_minus_flat,
      FrictionSolver::Method method
         = FrictionSolver::Method::NewtonRaphsonStable);

   /// Accessor for the configured sub-step sizes (test hook).
   const std::vector<real_t> &GetDeltaT() const { return deltaT_; }

   /// R5-006 (review round 5): explicitly close all probe files owned
   /// by the SEAS_DIAG_TPV104_STATE probe writer.  The driver SHOULD
   /// call this just before MPI_Finalize so the file destructors do
   /// not race the MPI teardown on implementations that hook stdio.
   /// A no-op under the default (SEAS_DIAG_TPV104_STATE undefined)
   /// build.  Safe to call from non-iterator code.
   static void CloseAllProbeFiles();

   /// Accessor for the configured quadrature weights (test hook).
   const std::vector<real_t> &GetTimeWeights() const { return time_weights_; }

private:
   FaultFaceFlux          &flux_;
   const SlipLawSRWPsi    &state_evo_;
   std::vector<real_t>     deltaT_;
   std::vector<real_t>     time_weights_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV104_SUBSTEP_ITERATOR_HPP
