// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV104 nucleation accumulator — cumulative per-sub-step injection.
//
// R-005 (2026-04-24) Step 6 of TPV104 implementation.  Implements the
// smoothStep / smoothStepIncrement pair (SCEC "Gaussian" nucleation ramp,
// TPV104 docs eq. 14) and a per-DOF accumulator that mirrors the
// reference FVW runtime's `adjustInitialStress` pattern: instead of
// overwriting the nucleation channel with the full ramp value at each
// macro-step, the accumulator ADDS the increment
//   Δτ · [ smoothStep(t_end) − smoothStep(t_end − dt) ]
// at every ADER sub-step.  Over the [0, T_nuc] interval the sum
// telescopes to the full perturbation Δτ·F(r).
//
// Do NOT modify `dynamic/tpv102_setup_total.hpp::ApplyNucleationPrestress`
// — that path is shared with TPV102 and on the extreme-care list.  This
// file is a TPV104-only sibling.  Reason: TPV102 overwrites the channel
// each macro-step; TPV104 accumulates per-sub-step.  Duplication is
// intentional per `feedback_tpv102_bp5_no_shared_edit`.

#ifndef MFEM_SEAS_TPV104_NUCLEATION_HPP
#define MFEM_SEAS_TPV104_NUCLEATION_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"
#include "../config/tpv104_params.hpp"

#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

/// SCEC "Gaussian" smoothStep ramp.
///   smoothStep(t, t0) = 0                                    t ≤ 0
///                     = exp((t-t0)² / (t·(t-2·t0)))          0 < t < t0
///                     = 1                                    t ≥ t0
/// Matches the reference runtime's `smoothStep` and SCEC TPV104 eq. 14.
inline real_t SmoothStep_TPV104(real_t current_time, real_t t0)
{
   if (current_time <= 0.0) { return 0.0; }
   if (current_time < t0)
   {
      const real_t tau = current_time - t0;
      return std::exp(tau * tau / (current_time * (current_time - 2.0 * t0)));
   }
   return 1.0;
}

/// Per-sub-step smoothStep increment:
///   ΔS(t, dt, t0) = smoothStep(t, t0) − smoothStep(t − dt, t0).
/// At every sub-step the accumulator adds ΔS · spatial-factor · amplitude.
/// Summed over [0, T_nuc] the increments telescope to the full ramp.
///
/// R3-007 (review round 3): `dt` must be finite and strictly positive.
/// A NaN `dt` silently routes `current_time - dt = NaN` through the
/// `NaN < t0` branch (which is false), producing `smoothStep(NaN, t0) =
/// 1.0` and feeding a spurious negative increment into the accumulator.
/// A +∞ `dt` similarly inflates a single sub-step to the full ramp.
inline real_t SmoothStepIncrement_TPV104(real_t current_time, real_t dt,
                                         real_t t0)
{
   MFEM_ASSERT(std::isfinite(dt) && dt > 0.0,
               "SmoothStepIncrement_TPV104: dt must be finite and "
               "positive; got dt = " << dt);
   return SmoothStep_TPV104(current_time, t0)
          - SmoothStep_TPV104(current_time - dt, t0);
}

/// Reset the per-DOF nucleation accumulator channels to zero.
/// Call once at driver init before the time loop.  After this call, the
/// per-DOF fields `tau1_nuc`, `tau2_nuc`, `sigma_n_nuc` are all 0.
inline void ResetNucleationAccumulator_TPV104(std::vector<DOFData> &dof_data)
{
   for (auto &d : dof_data)
   {
      d.tau1_nuc    = 0.0;
      d.tau2_nuc    = 0.0;
      d.sigma_n_nuc = 0.0;
   }
}

/// Apply ONE sub-step's worth of nucleation increment to the per-DOF
/// persistent nucleation channel.
///
/// Pattern: at sub-step o with endpoint time `t_substep_end` and size
/// `dt_substep`:
///   Δ = SmoothStepIncrement_TPV104(t_substep_end, dt_substep, T_nuc)
///   for every QP i:
///     tau2_nuc[i] += Δ · NucleationSpatial_TPV104(r_i) · Δτ₀
///   tau1_nuc / sigma_n_nuc are left at 0 — TPV104 is pure strike-slip
///   with no normal-stress nucleation.
///
/// The caller (Step-7 sub-step iterator) must invoke this once per ADER
/// sub-step.  Calling it once per macro-step would underestimate the
/// perturbation because smoothStep is nonlinear in t.
///
/// @param[in,out] dof_data      Fault DOFData (only `tau2_nuc` mutated).
/// @param[in]     fault_coords  Per-QP physical coords;
///                              `x` = along-strike, `|z|` = down-dip.
/// @param[in]     t_substep_end Sub-step endpoint time [s].
/// @param[in]     dt_substep    Sub-step size [s].
inline void ApplyNucleationIncremental_TPV104(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   real_t t_substep_end,
   real_t dt_substep)
{
   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= n,
               "ApplyNucleationIncremental_TPV104: fault_coords size "
               << fault_coords.size() << " < dof_data size " << n);

   const real_t dS = SmoothStepIncrement_TPV104(t_substep_end, dt_substep,
                                                TPV104Params::nuc_T);
   // R3-004 (review round 3): clamp non-positive increments to 0.
   // `smoothStep` is monotonically non-decreasing so mathematically
   // dS ≥ 0, but floating-point round-off at the ramp boundaries can
   // produce a sub-ULP negative dS (two close-to-1 values subtracted
   // via different code paths).  `dS == 0.0` exact-equality would miss
   // this; `dS <= 0.0` both fast-paths the before/after-ramp branches
   // and guards against sub-ULP negative drift into the accumulator.
   if (dS <= 0.0) { return; }

   for (int i = 0; i < n; ++i)
   {
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));
      const real_t dx = along_strike - TPV104Params::hypo_along_strike;
      const real_t dz = down_dip     - TPV104Params::hypo_down_dip;
      const real_t r  = std::sqrt(dx * dx + dz * dz);
      const real_t F  = NucleationSpatial_TPV104(r);
      dof_data[i].tau2_nuc += dS * F * TPV104Params::nuc_dtau;
      // tau1_nuc and sigma_n_nuc are intentionally not updated
      // (pure strike-slip; no normal-stress nucleation in TPV104).
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV104_NUCLEATION_HPP
