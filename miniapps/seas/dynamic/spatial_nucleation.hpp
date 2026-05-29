// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/spatial_nucleation.hpp — Phase N of
// safs/project_7.0_alternative/document/05_18_2026/PLAN_first_safs_run.md.
//
// The single nucleation mechanism for `seas_spatial_dyn_driver`:
// **gradual_overstress**.  Per-DOF shear-traction accumulator that ramps
// smoothly from 0 to a full-amplitude target `Δτ · F(r)` over the
// interval `[0, T_nuc_s]`, then is a no-op at every later sub-step.
//
//   - Spatial factor F(r): Gaussian centred on `(center_*_m)` with
//     e-fold radii `radius_dip_m` (down-dip) / `radius_strike_m`
//     (along-strike).  F = 1 at the centre; numerically zero outside
//     ~3 · radius_*.
//   - Temporal factor: SCEC smoothStep(t, t0) — C∞ everywhere except
//     t = 0; no step discontinuity like one-shot overstress.
//   - Per-sub-step: ΔS(t, Δt, t0) = SmoothStep(t, t0) - SmoothStep(t - Δt, t0)
//     telescopes exactly over [0, T_nuc_s] to the full target.
//
// The TPV104 driver continues to use `dynamic/tpv104_nucleation.hpp`
// verbatim (pure strike-slip, TPV104Params::nuc_T, etc.).  This file is
// the spatial-driver-specific sibling: generic smoothStep, generic
// Gaussian, per-DOF (dip, strike) projection, writes BOTH `tau1_nuc`
// AND `tau2_nuc`.

#ifndef MFEM_SEAS_SPATIAL_NUCLEATION_HPP
#define MFEM_SEAS_SPATIAL_NUCLEATION_HPP

#include "mfem.hpp"

#include "fault_face_flux.hpp"   // DOFData

#include <functional>
#include <vector>

namespace mfem
{
namespace seas
{
namespace spatial
{

// =====================================================================
// Temporal helpers — SCEC smoothStep ramp (generic).
// =====================================================================

/// @brief SCEC "Gaussian" smoothStep ramp (generic — no TPV104 ties).
///
///   smoothStep(t, t0) = 0                                    t ≤ 0
///                     = exp((t - t0)² / (t · (t - 2·t0)))    0 < t < t0
///                     = 1                                    t ≥ t0
///
/// C∞ everywhere except `t = 0`.  Mirror of `SmoothStep_TPV104`.
real_t SmoothStep(real_t current_time, real_t t0);

/// @brief Per-sub-step smoothStep increment.
///
/// `dt` MUST be finite and strictly positive (`MFEM_ASSERT`).  A NaN or
/// non-positive `dt` would silently route through `current_time - dt`
/// and corrupt the accumulator.
///
/// Summed over an arbitrary partition of `[0, t0]` the increments
/// telescope to `SmoothStep(t0, t0) - SmoothStep(0, t0) = 1`.
real_t SmoothStepIncrement(real_t current_time, real_t dt, real_t t0);

// =====================================================================
// Spatial helper — Gaussian factor in the fault-local frame.
// =====================================================================

/// @brief Gaussian-shaped spatial factor F(r) in the fault-local
/// `(dip, strike)` frame.
///
/// Computes `dx = (dof - center) · dip` and `ds = (dof - center) · strike`
/// (signed projections), then returns
/// `exp(− ((dx / radius_dip_m)² + (ds / radius_strike_m)²))`.
///
/// Both radii MUST be > 0 (`MFEM_ASSERT`).  Returns 0 immediately if
/// the exponent overflows below -700 to avoid denormal underflow noise.
///
/// @param[in] dof_xyz          Physical coord (3 components) at the DOF.
/// @param[in] dof_basis_dip    Per-DOF dip-direction unit vector
///                             (3 components, contiguous).
/// @param[in] dof_basis_strike Per-DOF strike-direction unit vector
///                             (3 components, contiguous).
/// @param[in] center_xyz       Hypocenter physical coords (3).
/// @param[in] radius_dip_m     e-fold radius along the dip direction (m).
/// @param[in] radius_strike_m  e-fold radius along strike (m).
/// @returns F(r) in [0, 1].
real_t GaussianFactorFaceLocal(const real_t* dof_xyz,
                               const real_t* dof_basis_dip,
                               const real_t* dof_basis_strike,
                               const real_t* center_xyz,
                               real_t radius_dip_m,
                               real_t radius_strike_m);

// =====================================================================
// Spec + resolved per-DOF parameters.
// =====================================================================

/// @brief `[nucleation.gradual_overstress]` TOML sub-block.
///
/// `radius_dip_m`, `radius_strike_m`, `T_nuc_s` MUST be > 0 when the
/// `[nucleation]` block is enabled (validated by the parser).  The
/// schema's `t0_smooth_s` field was DROPPED (R-007): `T_nuc_s` plays
/// the smoothStep `t0` role (matches TPV104's `nuc_T` convention).
struct GradualOverstressSpec
{
   real_t center_x_m          = 0.0;
   real_t center_y_m          = 0.0;
   real_t center_z_m          = 0.0;
   real_t radius_dip_m        = 0.0;   ///< > 0 required when enabled
   real_t radius_strike_m     = 0.0;   ///< > 0 required when enabled
   real_t delta_tau_dip_pa    = 0.0;   ///< may be 0
   real_t delta_tau_strike_pa = 0.0;   ///< may be 0
   real_t T_nuc_s             = 0.0;   ///< > 0 required when enabled
};

/// @brief Per-DOF resolved gradual_overstress targets.
///
///   amplitude_dip(i)    = F(r_i) · delta_tau_dip_pa
///   amplitude_strike(i) = F(r_i) · delta_tau_strike_pa
///   radial(i)           = F(r_i) ∈ [0, 1]
///
/// All three Vectors are zero-sized when nucleation is disabled.  The
/// `radial` field is published as the static `nuc_radial_factor` fault
/// field by Parity Phase 6.
struct GradualOverstressPerDOFParams
{
   Vector amplitude_dip;
   Vector amplitude_strike;
   Vector radial;
};

/// @brief Build per-DOF amplitudes from the spec + per-DOF coords +
/// per-DOF basis.
///
/// The `dof_basis` layout is `(9, N)` column-major per
/// `drivers/spatial_dyn_driver.cpp:228, 257, 275-277`: rows 0..2 =
/// normal, rows 3..5 = tangent1 = dip, rows 6..8 = tangent2 = strike.
/// MFEM `DenseMatrix` is column-major so `&dof_basis(3, i)` is a
/// contiguous 3-vector pointer (dip) and `&dof_basis(6, i)` is
/// contiguous strike (R-001 fix).
///
/// When `!enabled`, returns three zero-sized Vectors immediately.
GradualOverstressPerDOFParams ResolveGradualOverstress(
   const GradualOverstressSpec& spec,
   bool                         enabled,
   const Vector&                dof_coords_3d,
   const DenseMatrix&           dof_basis);

// =====================================================================
// Per-sub-step accumulator.
// =====================================================================

/// @brief Apply ONE ADER sub-step's worth of gradual_overstress
/// increment to the per-DOF `DOFData::tau1_nuc / tau2_nuc` fields.
///
/// Mirror of `ApplyNucleationIncremental_TPV104`, generalised:
///   * writes BOTH `tau1_nuc` (dip) and `tau2_nuc` (strike) — SAFS
///     curvilinear fault has variable strike, so amplitudes are
///     already projected into the local `(dip, strike)` frame by the
///     resolver;
///   * `sigma_n_nuc` is left at 0 (this mechanism does NOT perturb σ_n);
///   * early-return when `params.amplitude_dip.Size() == 0` (nucleation
///     disabled) or when `dS <= 0` (round-off guard).
///
/// `dof_data.size()` MUST equal `params.amplitude_*.Size()` (matches
/// the local fault-DOF count).
void ApplyGradualOverstressIncrement(
   std::vector<DOFData>&                  dof_data,
   const GradualOverstressPerDOFParams&   params,
   real_t                                 T_nuc_s,
   real_t                                 t_substep_end,
   real_t                                 dt_substep);

/// @brief Phase 14.3 — ABSOLUTE (set, not accumulate) gradual-overstress
/// forcing at stage time `t`: writes
///   dof_data[i].tau1_nuc = SmoothStep(t, T_nuc_s) · amplitude_dip(i)
///   dof_data[i].tau2_nuc = SmoothStep(t, T_nuc_s) · amplitude_strike(i)
///
/// This is the RK-stage analogue of `ApplyGradualOverstressIncrement`: the
/// per-sub-step telescoped increment `SmoothStep(t) − SmoothStep(t−dt)` would
/// DOUBLE-APPLY across RK stages that revisit a sub-interval (e.g. RK4 stages 2
/// and 3 both at t+dt/2), so the RK stepper sets the absolute target each stage
/// instead.  Summed/sampled at any partition of `[0, T_nuc_s]` the absolute
/// value equals the telescoped total (no double-apply).  `sigma_n_nuc` is left
/// at 0 (this mechanism does not perturb σ_n).  Early-returns when nucleation
/// is disabled (zero-sized params).
void ApplyGradualOverstressAbsolute(
   std::vector<DOFData>&                  dof_data,
   const GradualOverstressPerDOFParams&   params,
   real_t                                 T_nuc_s,
   real_t                                 t);

// =====================================================================
// Phase 7 — compact-circular gradual overstress (TPV102/104).
// =====================================================================

/// @brief `[nucleation.gradual_overstress_compact_circular]` TOML sub-block
/// (TPV102/104).  Moved here from spatial_friction.hpp in Phase 7 (req 1).
/// `radius_m`, `T_nuc_s` MUST be > 0 when the `[nucleation]` block is enabled
/// (validated by the parser).  Pure strike-slip: a single `delta_tau_pa` seeds
/// the strike component only.
struct GradualOverstressCompactCircularSpec
{
   real_t center_x_m   = 0.0;
   real_t center_y_m   = 0.0;
   real_t center_z_m   = 0.0;
   real_t radius_m     = 0.0;   ///< R; > 0 required when enabled
   real_t delta_tau_pa = 0.0;   ///< strike overstress amplitude [Pa]
   real_t T_nuc_s      = 0.0;   ///< > 0 required when enabled
};

/// @brief SCEC compact-bell radial factor F(r) — Eq. (13) of the
/// TPV101/102/104 spec: `F = exp(r²/(r²−R²))` for `r < R`, `0` for `r ≥ R`.
/// `R` MUST be > 0 (`MFEM_ASSERT`).  Mirror of `NucleationSpatial_TPV104`
/// (config/tpv104_params.hpp), kept here so the spatial driver does not depend
/// on the standalone TPV params header.
real_t CompactBellFactor(real_t r, real_t R);

/// @brief Per-DOF resolved compact-circular targets (strike component only).
///   amplitude_strike(i) = F(r_i) · delta_tau_pa   (seeds DOFData.tau2_nuc)
///   radial(i)           = F(r_i) ∈ [0, 1]
/// Both Vectors are zero-sized when nucleation is disabled.
struct CompactCircularPerDOFParams
{
   Vector amplitude_strike;
   Vector radial;
};

/// @brief Build per-DOF compact-circular amplitudes.  `r_i` is the in-fault-
/// plane distance from the centre, measured via the per-DOF (dip, strike)
/// basis: `r = sqrt(((dof−c)·dip)² + ((dof−c)·strike)²)` — equals the SCEC
/// `sqrt(Δalong_strike² + Δdown_dip²)` for the planar y=0 fault.  `dof_basis`
/// is the `(9, N)` column-major layout (rows 3..5 = dip, 6..8 = strike), as in
/// `ResolveGradualOverstress`.  Returns zero-sized Vectors when `!enabled`.
CompactCircularPerDOFParams ResolveGradualOverstressCompactCircular(
   const GradualOverstressCompactCircularSpec& spec,
   bool                                        enabled,
   const Vector&                               dof_coords_3d,
   const DenseMatrix&                          dof_basis);

/// @brief Apply ONE ADER sub-step's worth of compact-circular increment to
/// `DOFData::tau2_nuc` (strike).  Same smoothStep telescoping + `dS <= 0`
/// guard + post-ramp fast path as `ApplyGradualOverstressIncrement`; `tau1_nuc`
/// / `sigma_n_nuc` are left at 0 (pure strike-slip).
void ApplyGradualOverstressCompactCircularIncrement(
   std::vector<DOFData>&                  dof_data,
   const CompactCircularPerDOFParams&     params,
   real_t                                 T_nuc_s,
   real_t                                 t_substep_end,
   real_t                                 dt_substep);

/// @brief Phase 14.3 — ABSOLUTE (set, not accumulate) compact-circular forcing
/// at stage time `t`: writes
///   dof_data[i].tau2_nuc = SmoothStep(t, T_nuc_s) · amplitude_strike(i)
/// (pure strike-slip; tau1_nuc / sigma_n_nuc left at 0).  Same RK-stage
/// rationale as `ApplyGradualOverstressAbsolute`.  Early-returns when disabled.
void ApplyGradualOverstressCompactCircularAbsolute(
   std::vector<DOFData>&                  dof_data,
   const CompactCircularPerDOFParams&     params,
   real_t                                 T_nuc_s,
   real_t                                 t);

// =====================================================================
// Phase 7 — instantaneous circular overstress (TPV31, one-shot).
// =====================================================================

/// @brief `[nucleation.instantaneous_overstress_circular]` TOML sub-block
/// (TPV31).  Moved here from spatial_friction.hpp in Phase 7 (req 1).
/// `radius_m` MUST be > 0 and `taper_m` >= 0 when enabled (parser-validated).
struct InstantaneousOverstressCircularSpec
{
   real_t center_x_m   = 0.0;
   real_t center_y_m   = 0.0;
   real_t center_z_m   = 0.0;
   real_t radius_m     = 0.0;   ///< R; > 0 required when enabled
   real_t taper_m      = 0.0;   ///< cosine taper width [m] (>= 0)
   real_t delta_tau_pa = 0.0;   ///< overstress amplitude [Pa]
   /// Phase 10 (TPV31, spec p. 7): per-DOF amplitude scaling by
   /// `mu(point)/mu_ref` (the same depth-dependent shear modulus that scales
   /// the depth-proportional background stress).  `<= 0` (the default, and the
   /// value for every non-TPV31 config that omits the key) DISABLES the
   /// scaling — `amplitude_strike` is then the uniform `delta_tau_pa`, exactly
   /// as before.  The `mu(point)` lookup is supplied by the resolver's optional
   /// `mu_at_xyz` callback (the driver passes the depth-profile evaluator).
   real_t mu_ref_pa    = 0.0;
};

/// @brief Cosine-tapered radial factor for the one-shot patch (TPV31):
///   1                              r ≤ R
///   0.5·(1 + cos(π·(r−R)/taper))   R < r < R + taper
///   0                              r ≥ R + taper
/// `taper <= 0` degenerates to a hard cutoff (1 for r ≤ R, 0 beyond).  `R`
/// MUST be > 0 (`MFEM_ASSERT`).  Returns a factor in [0, 1].
real_t CosineTaperFactor(real_t r, real_t R, real_t taper);

/// @brief Per-DOF resolved instantaneous-circular amplitudes (strike only).
///   amplitude_strike(i) = CosineTaperFactor(r_i, R, taper) · delta_tau_pa
/// Zero-sized when nucleation is disabled.
struct InstantaneousOverstressPerDOFParams
{
   Vector amplitude_strike;
};

/// @brief Build per-DOF instantaneous-circular amplitudes (same in-fault-plane
/// `r` measure as `ResolveGradualOverstressCompactCircular`).  The patch is
/// seeded ONCE at t=0 by the caller (see `INucleationMethod::ApplyOnce`), which
/// adds `amplitude_strike(i)` to `DOFData::tau2_nuc`.  Returns zero-sized when
/// `!enabled`.
///
/// `mu_at_xyz` (optional): when non-empty AND `spec.mu_ref_pa > 0`, each DOF's
/// amplitude is scaled by `mu_at_xyz(x,y,z)/spec.mu_ref_pa` (TPV31 spec p. 7
/// per-DOF mu(depth)/mu_0 scaling).  Empty callback or `mu_ref_pa <= 0` ⇒ no
/// scaling (uniform `delta_tau_pa`), byte-identical to the pre-Phase-10 path.
InstantaneousOverstressPerDOFParams ResolveInstantaneousOverstressCircular(
   const InstantaneousOverstressCircularSpec& spec,
   bool                                       enabled,
   const Vector&                              dof_coords_3d,
   const DenseMatrix&                         dof_basis,
   const std::function<real_t(real_t, real_t, real_t)>& mu_at_xyz = {});

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_SPATIAL_NUCLEATION_HPP
