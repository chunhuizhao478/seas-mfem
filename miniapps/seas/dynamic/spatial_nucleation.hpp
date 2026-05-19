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

#include <cmath>
#include <functional>
#include <limits>
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

// =====================================================================
// GradualOverstressCompactCircular: SCEC TPV102/TPV104-family
// nucleation.  Same temporal smoothStep ramp as `GradualOverstress` but
// a circular *compact-support* spatial bell instead of the Gaussian:
//
//     F(r) = exp(r² / (r² − R²))   for r < R
//          = 0                     for r ≥ R
//
// where `r = sqrt(dx² + ds²)` with `(dx, ds) = (dof − center) · (dip,
// strike)` (fault-local projection).  F(0) = exp(0) = 1, F(R⁻) → 0
// (C∞ at r = R because every derivative vanishes).  No benchmark tag in
// the name; the formula is verbatim SCEC TPV102/104 spec Eq. (13).
// =====================================================================

/// @brief `[nucleation.gradual_overstress_compact_circular]` TOML sub-block.
struct GradualOverstressCompactCircularSpec
{
   real_t center_x_m          = 0.0;
   real_t center_y_m          = 0.0;
   real_t center_z_m          = 0.0;
   real_t radius_m            = 0.0;   ///< support radius R; > 0 when enabled
   real_t delta_tau_dip_pa    = 0.0;   ///< may be 0
   real_t delta_tau_strike_pa = 0.0;   ///< may be 0
   real_t T_nuc_s             = 0.0;   ///< > 0 when enabled (smoothStep t0)
};

/// @brief Per-DOF resolved compact-circular overstress targets.
struct GradualOverstressCompactCircularPerDOFParams
{
   Vector amplitude_dip;
   Vector amplitude_strike;
   Vector radial;
};

/// @brief Build per-DOF amplitudes from the spec + per-DOF coords +
/// per-DOF basis.  Same `dof_basis` layout as
/// `ResolveGradualOverstress` (rows 3..5 = dip, rows 6..8 = strike).
///
/// When `!enabled`, returns three zero-sized Vectors immediately.
GradualOverstressCompactCircularPerDOFParams
ResolveGradualOverstressCompactCircular(
   const GradualOverstressCompactCircularSpec& spec,
   bool                                        enabled,
   const Vector&                               dof_coords_3d,
   const DenseMatrix&                          dof_basis);

/// @brief Apply ONE ADER sub-step's worth of compact-circular overstress
/// increment to `DOFData::tau1_nuc / tau2_nuc`.  Same telescoping
/// semantics as `ApplyGradualOverstressIncrement`.
void ApplyGradualOverstressCompactCircularIncrement(
   std::vector<DOFData>&                                  dof_data,
   const GradualOverstressCompactCircularPerDOFParams&    params,
   real_t                                                 T_nuc_s,
   real_t                                                 t_substep_end,
   real_t                                                 dt_substep);

// =====================================================================
// SquareOverstress: rectangular indicator-function shape (sibling of
// GradualOverstress).  The temporal ramp is the same SCEC smoothStep
// (`T_nuc_s`); the spatial shape is one or more disjoint rectangular
// patches with per-patch (delta_tau_dip_pa, delta_tau_strike_pa)
// amplitudes.  At each DOF, the LAST patch containing the DOF wins
// (matches StressPatch convention).  A DOF outside every patch sees
// zero amplitude.
// =====================================================================

/// @brief One rectangular patch in the `[[nucleation.square_overstress.patch]]`
/// array.  Patch indicator is `|x - center_x_m| <= half_x_m` AND likewise
/// for y, z; missing/+inf `half_*_m` means "no constraint along that axis".
struct SquareOverstressPatch
{
   real_t center_x_m         = 0.0;
   real_t center_y_m         = 0.0;
   real_t center_z_m         = 0.0;
   real_t half_x_m           = std::numeric_limits<real_t>::infinity();
   real_t half_y_m           = std::numeric_limits<real_t>::infinity();
   real_t half_z_m           = std::numeric_limits<real_t>::infinity();
   real_t delta_tau_dip_pa    = 0.0;
   real_t delta_tau_strike_pa = 0.0;

   bool contains(real_t x, real_t y, real_t z) const
   {
      return std::abs(x - center_x_m) <= half_x_m
          && std::abs(y - center_y_m) <= half_y_m
          && std::abs(z - center_z_m) <= half_z_m;
   }
};

/// @brief `[nucleation.square_overstress]` TOML sub-block: list of
/// rectangular patches sharing one smoothStep `T_nuc_s`.
struct SquareOverstressSpec
{
   std::vector<SquareOverstressPatch> patches;
   real_t T_nuc_s = 0.0;   ///< > 0 required when enabled
};

/// @brief Per-DOF resolved square_overstress targets.  Same Vector
/// layout as GradualOverstressPerDOFParams — `radial(i)` is 1 inside
/// any patch, 0 outside (the indicator function), so downstream
/// ParaView publication (`nuc_radial_factor`) still works.
struct SquareOverstressPerDOFParams
{
   Vector amplitude_dip;
   Vector amplitude_strike;
   Vector radial;
};

/// @brief Build per-DOF amplitudes from spec + per-DOF coords.  No
/// `dof_basis` needed — patches are tested in physical coords and the
/// per-patch `delta_tau_*` is taken as already-fault-local (caller is
/// responsible for matching the dip/strike convention used elsewhere).
SquareOverstressPerDOFParams ResolveSquareOverstress(
   const SquareOverstressSpec& spec,
   bool                        enabled,
   const Vector&               dof_coords_3d);

/// @brief Apply ONE ADER sub-step's worth of square_overstress
/// increment to the per-DOF `DOFData::tau1_nuc / tau2_nuc` fields.
/// Same telescoping semantics as ApplyGradualOverstressIncrement.
void ApplySquareOverstressIncrement(
   std::vector<DOFData>&                 dof_data,
   const SquareOverstressPerDOFParams&   params,
   real_t                                T_nuc_s,
   real_t                                t_substep_end,
   real_t                                dt_substep);

// =====================================================================
// InstantaneousOverstressCircular: TPV31-style nucleation.  Applied
// ONCE at init time (no time-ramp, no per-sub-step accumulation), with
// a circular shape, cosine taper from r=radius_inner_m to
// r=radius_outer_m, and per-DOF µ-scaling.  Pure strike-direction
// perturbation by default (TPV31 is right-lateral).
//
// Spec formula (SCEC TPV31 §7):
//   r(x, y, z) = sqrt((x - cx)² + (y - cy)² + (z - cz)²)  but only the
//                in-fault-plane components contribute non-trivially —
//                for a planar fault the off-plane projection is
//                identically zero at every fault DOF, so this matches
//                the spec's 2D distance formula exactly.
//   tau_nuke(r) = delta_tau_peak_pa · (µ(point) / µ_ref_pa),  r ≤ ri
//               = (delta_tau_peak_pa / 2) · (1 + cos(π(r-ri)/(ro-ri)))
//                                          · (µ(point) / µ_ref_pa),  ri ≤ r ≤ ro
//               = 0                                       otherwise
//
// Notation: ri = radius_inner_m, ro = radius_outer_m.
// =====================================================================

/// @brief `[nucleation.instantaneous_overstress_circular]` TOML sub-block.
struct InstantaneousOverstressCircularSpec
{
   real_t center_x_m         = 0.0;
   real_t center_y_m         = 0.0;
   real_t center_z_m         = 0.0;
   real_t radius_inner_m     = 0.0;   ///< full-amplitude radius (1400 m for TPV31)
   real_t radius_outer_m     = 0.0;   ///< taper-end radius     (2000 m for TPV31)
   real_t delta_tau_peak_pa  = 0.0;   ///< peak Δτ at the reference shear modulus
   real_t mu_ref_pa          = 32.03812032e9;  ///< µ_0 in the spec (TPV5/TPV31)
   /// Apportionment of the spec's pure shear into the fault-local
   /// (dip, strike) frame.  TPV31 is pure right-lateral so dip = 0,
   /// strike = 1 (the resolver multiplies these factors by tau_nuke
   /// before writing tau{1,2}_nuc).  Negative values flip polarity.
   real_t dip_fraction       = 0.0;
   real_t strike_fraction    = 1.0;
};

/// @brief Per-DOF resolved instantaneous_overstress_circular targets.
/// Same Vector layout as GradualOverstressPerDOFParams — `radial(i)`
/// records the spatial factor in [0, 1] (excluding µ-scaling) so the
/// downstream ParaView publication (`nuc_radial_factor`) still works.
struct InstantaneousOverstressCircularPerDOFParams
{
   Vector amplitude_dip;
   Vector amplitude_strike;
   Vector radial;
};

/// @brief Build per-DOF amplitudes from the spec + per-DOF coords +
/// per-DOF µ values.
///
/// @param mu_per_dof  Vector of size N — shear modulus µ at each fault
///                    DOF (in Pa).  Driver populates this by evaluating
///                    its MaterialField at every fault-DOF physical
///                    coord.
///
/// When `!enabled`, returns three zero-sized Vectors immediately.
InstantaneousOverstressCircularPerDOFParams
ResolveInstantaneousOverstressCircular(
   const InstantaneousOverstressCircularSpec& spec,
   bool                                       enabled,
   const Vector&                              dof_coords_3d,
   const Vector&                              mu_per_dof);

/// @brief Write the full per-DOF amplitude into DOFData::tau{1,2}_nuc
/// at init time (NOT incremental).  Mirrors a one-shot t=0 overstress.
///
/// Called ONCE during driver setup, BEFORE the time loop.  After this
/// call, the per-sub-step nucleation accumulator should NOT also write
/// to tau{1,2}_nuc for this kind (the spec is instantaneous; no ramp).
void ApplyInstantaneousOverstressCircular(
   std::vector<DOFData>&                              dof_data,
   const InstantaneousOverstressCircularPerDOFParams& params);

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

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_SPATIAL_NUCLEATION_HPP
