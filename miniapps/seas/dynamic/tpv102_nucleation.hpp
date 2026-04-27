// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV102 nucleation accumulator — cumulative per-sub-step injection.
//
// Mirrors `dynamic/tpv104_nucleation.hpp` with TPV102Params substitutions:
//   - amplitude Δτ0    = 25 MPa     (vs TPV104's 45 MPa)
//   - rise time T_nuc  = 1 s        (same as TPV104)
//   - radius R         = 3 km       (same as TPV104)
//   - hypocenter       = (0, 7.5 km) (same as TPV104)
//
// The accumulator pattern is identical to TPV104's: instead of overwriting
// `tau2_0` with `tau_ini + Δτ` once per macro-step (the legacy
// `tpv102_setup.hpp::ApplyNucleation` overwrite-pattern), this header
// ADDS the increment
//     Δτ0 · F(r) · [ smoothStep(t_end) − smoothStep(t_end − dt) ]
// at every ADER sub-step into the per-DOF persistent nucleation channel
// `tau2_nuc`.  Over the [0, T_nuc] interval the increments telescope to
// the full perturbation Δτ₀·F(r), matching SCEC TPV101/102 §"Nucleation
// Method" Eqs. (9)–(11).
//
// `tpv102_setup_total.hpp::ApplyNucleationPrestress` (the legacy total-Q
// path, [C2]) is NOT modified — the new TPV102 driver uses the
// fluctuation-Q + accumulator path here, leaving the legacy header
// available for the older driver.

#ifndef MFEM_SEAS_TPV102_NUCLEATION_HPP
#define MFEM_SEAS_TPV102_NUCLEATION_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"
#include "../config/tpv102_params.hpp"

#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

/// SCEC "Gaussian" smoothStep ramp (TPV101/102 Eq. (11)).
///   smoothStep(t, t0) = 0                                    t ≤ 0
///                     = exp((t-t0)² / (t·(t-2·t0)))          0 < t < t0
///                     = 1                                    t ≥ t0
inline real_t SmoothStep_TPV102(real_t current_time, real_t t0)
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
/// Summed over [0, T_nuc] the increments telescope to the full ramp.
inline real_t SmoothStepIncrement_TPV102(real_t current_time, real_t dt,
                                         real_t t0)
{
   MFEM_ASSERT(std::isfinite(dt) && dt > 0.0,
               "SmoothStepIncrement_TPV102: dt must be finite and "
               "positive; got dt = " << dt);
   return SmoothStep_TPV102(current_time, t0)
          - SmoothStep_TPV102(current_time - dt, t0);
}

/// Reset the per-DOF nucleation accumulator channels to zero.
/// Call once at driver init before the time loop.
inline void ResetNucleationAccumulator_TPV102(std::vector<DOFData> &dof_data)
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
///   Δ = SmoothStepIncrement_TPV102(t_substep_end, dt_substep, T_nuc)
///   for every QP i:
///     tau2_nuc[i] += Δ · F(r_i) · Δτ₀
///   tau1_nuc / sigma_n_nuc are left at 0 — TPV102 is pure strike-slip
///   with no normal-stress nucleation.
///
/// Calling once per ADER sub-step is required for high-order accuracy;
/// calling once per macro-step would underestimate the perturbation
/// because smoothStep is nonlinear in t.
///
/// @param[in,out] dof_data      Fault DOFData (only `tau2_nuc` mutated).
/// @param[in]     fault_coords  Per-QP physical coords;
///                              `x` = along-strike, `|z|` = down-dip.
/// @param[in]     t_substep_end Sub-step endpoint time [s].
/// @param[in]     dt_substep    Sub-step size [s].
inline void ApplyNucleationIncremental_TPV102(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   real_t t_substep_end,
   real_t dt_substep)
{
   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= n,
               "ApplyNucleationIncremental_TPV102: fault_coords size "
               << fault_coords.size() << " < dof_data size " << n);

   const real_t dS = SmoothStepIncrement_TPV102(t_substep_end, dt_substep,
                                                TPV102Params::nuc_T);
   // Mirror the TPV104 R3-004 guard: clamp non-positive dS to 0 so
   // sub-ULP negative drift at the ramp boundaries (and the t < 0 /
   // t > T fast-path) does not leak into the accumulator.
   if (dS <= 0.0) { return; }

   for (int i = 0; i < n; ++i)
   {
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));
      const real_t dx = along_strike - TPV102Params::hypo_along_strike;
      const real_t dz = down_dip     - TPV102Params::hypo_down_dip;
      const real_t r  = std::sqrt(dx * dx + dz * dz);
      const real_t F  = NucleationSpatial(r);
      dof_data[i].tau2_nuc += dS * F * TPV102Params::nuc_dtau;
   }
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV102_NUCLEATION_HPP
