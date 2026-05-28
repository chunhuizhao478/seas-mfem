// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_friction.hpp — Phase 1 of spatial_dynamic_rupture_plan.md.
//
// Consolidated header (rev-3 §"What is genuinely new") for everything in
// Phase 1 of the spatial SAFS dynamic-rupture plan:
//
//   - TOML schema types (SpatialFrictionConfig and friends)
//   - LoadSpatialFrictionConfig(...) — the parser entry point
//   - SpatialFrictionResolver — per-DOF rule evaluation
//   - LSWFrictionCoefficient_ForcedRupture(...) — TPV26/27 §Part 4 helper
//     (consumed by Phase H's EvaluateADER_LSW_ForcedRupture)
//   - SpatialTimeParseSeconds(...) — small "0.01s"/"12s" parser
//
// Schema reference:
//   safs/project_7.0_alternative/document/spatial_friction_config_schema.md
//
// Implementation notes:
//   - All TOML I/O is gated by SEAS_USE_TOML (defined in the Makefile
//     when extern/toml11/toml.hpp is present).  When SEAS_USE_TOML is
//     not defined, LoadSpatialFrictionConfig aborts at runtime; the
//     resolver, time parser, and friction-coefficient helper compile
//     and work unconditionally.
//   - StressSource3D concept documentation lives in spatial_stress.hpp
//     (Phase 3/3b).  The resolver in this file does not depend on the
//     stress source.

#ifndef MFEM_SEAS_SPATIAL_FRICTION_HPP
#define MFEM_SEAS_SPATIAL_FRICTION_HPP

#include "mfem.hpp"

#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/spatial_nucleation.hpp"   // GradualOverstressSpec

#include <array>
#include <cmath>     // std::abs/std::tanh in inline SCECBoxcar (req 5)
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace spatial
{

// =====================================================================
//  Time parser
// =====================================================================

/// @brief Parse a time-string like "12s", "0.001s", "1.5e-2s" into
/// seconds.  Returns the parsed value (always >= 0 per validator
/// upstream).
///
/// The dynamic-rupture plan accepts only the `s` suffix (no yr / ka /
/// Ma); the quasi-dynamic sibling plan parser accepts more.
///
/// Special value: "auto" → returns -1.0 (sentinel meaning "compute from
/// CFL"; the driver chooses an actual Δt later).
///
/// MFEM_ABORTs with the offending string on parse failure.
real_t SpatialTimeParseSeconds(const std::string& s);

// =====================================================================
//  Configuration types (mirror §Phase 0 TOML grammar)
// =====================================================================

enum class FrictionLawKind { SlipWeakening, RateState };

struct PorePressureSpec
{
   real_t P_p_pa            = 0.0;
   real_t P_p_grad_pa_per_m = 0.0;
   real_t min_sigma_n_pa    = 0.0;
};

struct MaterialConstantFallback
{
   real_t lambda = 32.0e9;
   real_t mu     = 32.0e9;
   real_t rho    = 2670.0;
};

struct MeshSpec
{
   std::string path;
   int         order = 1;
};

/// Phase 6 req 3 selectors.
/// cfl_safety: Raw = use cfl as-is; Dg = apply the DG /(2p+1)·(1/3) safety
///   factor (D2 decision).  Default Raw.
/// fault_iterator: OneShot = single EvaluateADER per macro-step; Substep =
///   the ADER sub-step iterator.  Default OneShot.
/// interior_flux: Scalar = homogeneous Godunov (mixed_flux allowed); Matrix =
///   heterogeneous/bimaterial Riemann (mixed_flux forbidden).  Default Scalar.
enum class CflSafety        { Raw, Dg };
enum class FaultIteratorKind { OneShot, Substep };
enum class InteriorFlux     { Scalar, Matrix };

struct NumericsSpec
{
   int               ader_order     = 2;
   std::string       mixed_flux     = "none";
   real_t            cfl            = 0.5;
   bool              use_pml        = false;
   CflSafety         cfl_safety     = CflSafety::Raw;            // Phase 6 req 3
   FaultIteratorKind fault_iterator = FaultIteratorKind::OneShot;
   InteriorFlux      interior_flux  = InteriorFlux::Scalar;
};

struct TimeSpec
{
   real_t tfinal     = 12.0;
   real_t t_initial  = 0.0;
   real_t dt_initial = -1.0;     // -1 sentinel ⇒ "auto"
   real_t dt_max     = 0.1;
};

struct OutputSpec
{
   std::string output_dir;
   std::string restart_prefix         = "cp";
   // Parity Phase 2 default flip: SAFS production turns volume + bulk
   // OFF by default (multi-TB output otherwise).  Fault stays "hdf5".
   std::string paraview_volume        = "off";
   std::string paraview_bulk          = "off";
   std::string paraview_fault         = "hdf5";
   real_t      paraview_volume_dt     = 0.05;
   real_t      paraview_bulk_dt       = 0.05;
   real_t      paraview_fault_dt      = 0.001;
   real_t      paraview_volume_zfp_tol = 1e-3;
   real_t      paraview_bulk_zfp_tol   = 1e-3;
   real_t      paraview_fault_zfp_tol  = 1e-12;
   int         max_snapshots          = 5000;
   int         checkpoint_every_steps = 10000;

   // Parity Phase 2 additions (9 new fields).
   bool        paraview_enabled              = false;  ///< master gate
   int         paraview_every_steps          = 0;      ///< 0 = use *_dt
   bool        paraview_fault_legacy_ascii   = false;  ///< debug only
   int         paraview_volume_deflate_level = -1;     ///< -1 = none, [0..9]
   int         paraview_bulk_deflate_level   = -1;
   int         paraview_fault_deflate_level  = -1;
   real_t      paraview_coseismic_dt         = -1.0;   ///< -1 = unset
   real_t      paraview_nucleation_dt        = -1.0;
   real_t      paraview_interseismic_dt      = -1.0;
};

enum class VelocityModel { CVMH, CVMS_4_26_M01, MultiscaleStatewise };

struct VelocitySpec
{
   VelocityModel model = VelocityModel::CVMH;
   std::string   dataset_root;
   std::string   override_path;       // "" ⇒ resolve from model
   // When false, the driver skips the velocity-sidecar load entirely
   // and uses [material_constant_fallback] for (lambda, mu, rho).  In
   // that mode `model` / `dataset_root` / `override_path` are ignored
   // and may be empty.  Equivalent to passing --no-sidecar-material on
   // the CLI; the CLI flag still wins as an override.
   bool          use_sidecar = true;
};

enum class StressSourceKind { ConstantTensor, SidecarHDF5, FaultLocalPrestress };

/// Rectangular `tau_strike` patch for the FaultLocalPrestress source (D3.2,
/// Phase 6 req 2).  Inside the box `|coord - center| <= half` (per axis), the
/// uniform background `tau_strike` is overridden by `tau_strike_pa`.  A NaN
/// `center_*` or +inf `half_*` makes that axis unconstrained (always inside),
/// so an all-default patch covers the whole fault.  Patches are applied
/// last-match-wins (later patches in the list override earlier ones).  Used by
/// TPV205 (background +70 MPa; central/left/right patches).
struct FaultLocalPatch
{
   real_t center_x_m = std::numeric_limits<real_t>::quiet_NaN();
   real_t center_y_m = std::numeric_limits<real_t>::quiet_NaN();
   real_t center_z_m = std::numeric_limits<real_t>::quiet_NaN();
   real_t half_x_m   = std::numeric_limits<real_t>::infinity();
   real_t half_y_m   = std::numeric_limits<real_t>::infinity();
   real_t half_z_m   = std::numeric_limits<real_t>::infinity();
   real_t tau_strike_pa = 0.0;   ///< right-lateral POSITIVE strike pre-stress

   /// True iff (x,y,z) lies inside this patch box.  An unconstrained axis
   /// (NaN center or non-finite half) always matches.
   bool inside(real_t x, real_t y, real_t z) const
   {
      auto axis = [](real_t c, real_t cen, real_t half) {
         return std::isnan(cen) || !std::isfinite(half)
                || std::abs(c - cen) <= half;
      };
      return axis(x, center_x_m, half_x_m)
             && axis(y, center_y_m, half_y_m)
             && axis(z, center_z_m, half_z_m);
   }
};

struct StressSpec
{
   StressSourceKind kind = StressSourceKind::ConstantTensor;
   // ConstantTensor (Phase 3b):
   real_t sigma_xx_pa = 0.0;
   real_t sigma_yy_pa = 0.0;
   real_t sigma_zz_pa = 0.0;
   real_t sigma_xy_pa = 0.0;
   real_t sigma_yz_pa = 0.0;
   real_t sigma_xz_pa = 0.0;
   // SidecarHDF5 (Phase 3):
   std::string sidecar_path;
   // FaultLocalPrestress (Phase 6 / D3.2): constant background pre-stress in
   // the canonical fault-local frame, right-lateral / compression POSITIVE.
   // Seeded directly via FaultGeometry::ComputeParamsFaultLocal — NO Cauchy
   // projection.  tau2_0 = tau_strike, tau1_0 = tau_dip, sigma_n0 = sigma_n - P_p.
   real_t tau_strike_pa = 0.0;
   real_t tau_dip_pa    = 0.0;
   real_t sigma_n_pa    = 0.0;
   // Optional rectangular tau_strike patches (last-match-wins), layered on
   // the uniform background by FaultGeometry::ApplyFaultLocalStrikePatches.
   std::vector<FaultLocalPatch> fault_local_patches;
   // Common:
   PorePressureSpec pore_pressure;
};

/// Phase 6 req 5: SCEC boxcar B(offset, half, trans) — Eq.(5) of the
/// SCEC TPV101/102/104 spec.  B = 1 for |offset| <= half, a C∞ tanh
/// transition for half < |offset| < half+trans, and 0 beyond.  `trans <= 0`
/// degenerates to a hard boxcar (1 inside, 0 outside, no transition).
/// This is a config-frame copy of `Boxcar_TPV104` (config/tpv104_params.hpp);
/// the spatial driver keeps its own so it does not pull in the standalone
/// TPV params header.
inline real_t SCECBoxcar(real_t offset, real_t half, real_t trans)
{
   const real_t ax = std::abs(offset);
   if (ax <= half)              { return 1.0; }
   if (trans <= 0.0 || ax >= half + trans) { return 0.0; }
   return 0.5 * (1.0 + std::tanh(trans / (ax - half - trans)
                                 + trans / (ax - half)));
}

/// Spatial rule (D-3 relaxed validator; R-014 explicit defaults).
/// D-4 (rev-3): the `NucleationBox` kind is REMOVED.  Spatial rules
/// never override `tau_pre_*` or `sigma_n` — pre-stress comes entirely
/// from the stress source.
struct SpatialRule
{
   // Phase 6 req 5: `BoxcarTaper` is a smooth 3-D taper region (SCEC boxcar
   // product); unlike the hard Box/Depth kinds it returns a factor in [0,1]
   // (see BoxcarTaperFactor), enabling cohesion / parameter tapers.
   enum class Kind { Depth, Box, RegionAttribute, Barrier, BoxcarTaper };
   Kind kind = Kind::Depth;

   real_t x_min_m = -std::numeric_limits<real_t>::infinity();
   real_t x_max_m =  std::numeric_limits<real_t>::infinity();
   real_t y_min_m = -std::numeric_limits<real_t>::infinity();
   real_t y_max_m =  std::numeric_limits<real_t>::infinity();
   real_t z_min_m = -std::numeric_limits<real_t>::infinity();
   real_t z_max_m =  std::numeric_limits<real_t>::infinity();
   int    region_attr = -1;

   // Per-key overrides.  NaN sentinel means "do not override".
   real_t mu_s    = std::numeric_limits<real_t>::quiet_NaN();
   real_t mu_d    = std::numeric_limits<real_t>::quiet_NaN();
   real_t d_c     = std::numeric_limits<real_t>::quiet_NaN();
   real_t cohesion= std::numeric_limits<real_t>::quiet_NaN();

   real_t a       = std::numeric_limits<real_t>::quiet_NaN();
   real_t b       = std::numeric_limits<real_t>::quiet_NaN();
   real_t Dc      = std::numeric_limits<real_t>::quiet_NaN();
   real_t V_init  = std::numeric_limits<real_t>::quiet_NaN();
   real_t f_0     = std::numeric_limits<real_t>::quiet_NaN();
   real_t V_0     = std::numeric_limits<real_t>::quiet_NaN();
   real_t sigma_n = std::numeric_limits<real_t>::quiet_NaN();
   real_t eta     = std::numeric_limits<real_t>::quiet_NaN();
   real_t V_w     = std::numeric_limits<real_t>::quiet_NaN();  // Phase 6 req 4 (SRW per-QP V_w)

   // Phase 6 req 5: BoxcarTaper geometry.  Per-axis plateau half-width
   // boxcar_half_* and tanh transition boxcar_trans_*, centred at
   // boxcar_center_*.  A non-finite half makes that axis untapered (the
   // axis contributes factor 1).  Consumed by BoxcarTaperFactor.
   real_t boxcar_center_x_m = 0.0;
   real_t boxcar_center_y_m = 0.0;
   real_t boxcar_center_z_m = 0.0;
   real_t boxcar_half_x_m   =  std::numeric_limits<real_t>::infinity();
   real_t boxcar_half_y_m   =  std::numeric_limits<real_t>::infinity();
   real_t boxcar_half_z_m   =  std::numeric_limits<real_t>::infinity();
   real_t boxcar_trans_x_m  = 0.0;
   real_t boxcar_trans_y_m  = 0.0;
   real_t boxcar_trans_z_m  = 0.0;
   // Cohesion-taper endpoints (LSW BoxcarTaper rules): cohesion ramps from
   // cohesion_inner (plateau, factor 1) to cohesion_outer (factor 0) via
   // BoxcarTaperFactor.  NaN = unset (Phase 7 resolver falls back to the
   // block default uniformly).
   real_t cohesion_inner = std::numeric_limits<real_t>::quiet_NaN();
   real_t cohesion_outer = std::numeric_limits<real_t>::quiet_NaN();

   /// `matches` semantics per §Phase 1 Detailed Req. 4:
   ///   Depth            : only z bounds checked.
   ///   Box              : all 6 coord bounds checked.
   ///   RegionAttribute  : only region_attr checked.
   ///   Barrier          : same as Depth, plus any non-infinite
   ///                      x/y bounds also apply.
   ///   BoxcarTaper      : matches wherever the boxcar factor is non-zero
   ///                      (BoxcarTaperFactor(x,y,z) > 0); coord bounds are
   ///                      ignored (the taper geometry governs the region).
   bool matches(real_t x, real_t y, real_t z, int attr) const;

   /// Phase 6 req 5: smooth taper factor in [0,1] for a BoxcarTaper rule,
   /// F = ∏_axis SCECBoxcar(coord - center, half, trans).  An axis with a
   /// non-finite half contributes factor 1 (untapered).  Returns 1.0 for
   /// any non-BoxcarTaper kind (those have no taper).
   real_t BoxcarTaperFactor(real_t x, real_t y, real_t z) const;
};

/// @brief Phase N — the single nucleation mechanism for the spatial
/// driver is `gradual_overstress` (per-DOF Δτ accumulator,
/// smoothStep-ramped, Gaussian-shaped).  The obsolete `StrengthReduction`
/// (TPV26/27 §Part 4) and `Overstress` (TPV205-style one-shot) kinds were
/// removed from the spatial code path.  Native TPV* drivers keep their
/// own paths via `dynamic/tpv104_nucleation.hpp` + the still-live
/// `FaultFrictionLaw::LSW_ForcedRupture` flux dispatch.
enum class NucleationKind
{
   GradualOverstress,                  ///< Gaussian gradual overstress (SAFS)
   GradualOverstressCompactCircular,   ///< SCEC compact bell (TPV102/104)
   InstantaneousOverstressCircular     ///< one-shot cosine-tapered patch (TPV31)
};

/// SCEC compact-circular gradual overstress (TPV102/104): radial bell
/// F(r) = exp(r^2/(r^2 - R^2)) for r<R (0 outside) times the smoothStep
/// temporal ramp over [0, T_nuc_s].  Config-only here; the resolver/applicator
/// is Phase 7 (nucleation_methods module).
struct GradualOverstressCompactCircularSpec
{
   real_t center_x_m   = 0.0;
   real_t center_y_m   = 0.0;
   real_t center_z_m   = 0.0;
   real_t radius_m     = 0.0;   ///< R; > 0 required when enabled
   real_t delta_tau_pa = 0.0;   ///< strike overstress amplitude [Pa]
   real_t T_nuc_s      = 0.0;   ///< > 0 required when enabled
};

/// Instantaneous circular overstress (TPV31): a one-shot cosine-tapered
/// strike-overstress patch.  Config-only here; the applicator is Phase 7.
struct InstantaneousOverstressCircularSpec
{
   real_t center_x_m   = 0.0;
   real_t center_y_m   = 0.0;
   real_t center_z_m   = 0.0;
   real_t radius_m     = 0.0;   ///< > 0 required when enabled
   real_t taper_m      = 0.0;   ///< cosine taper width [m] (>= 0)
   real_t delta_tau_pa = 0.0;   ///< overstress amplitude [Pa]
};

/// @brief Top-level `[nucleation]` TOML block.  When absent in the
/// TOML, `enabled == false` and the driver runs with no nucleation
/// perturbation (tau{1,2}_nuc stay at the init-time zero seeded by
/// `InitializeFaultDOFs_Spatial`).  TPV205 uses no [nucleation] block (static).
struct NucleationSpec
{
   NucleationKind        kind     = NucleationKind::GradualOverstress;
   bool                  enabled  = false;
   GradualOverstressSpec gradual_overstress;
   GradualOverstressCompactCircularSpec compact_circular;        // Phase 6 req 6
   InstantaneousOverstressCircularSpec  instantaneous_circular;  // Phase 6 req 6
};

struct SlipWeakeningBlock
{
   real_t mu_s_default     = 1.1;   // geoffrey2010.md
   real_t mu_d_default     = 0.5;
   real_t d_c_default      = 0.5;
   real_t cohesion_default = 0.0;
   std::vector<SpatialRule> spatial;
};

/// Generic flat-clamped piecewise-linear 1-D interpolant (Phase 11b).  Used for
/// the depth-profile a(z) and (a-b)(z) curves.
struct PiecewiseLinear1D
{
   std::vector<real_t> x;   ///< ascending, strictly increasing; size >= 2
   std::vector<real_t> y;   ///< same length as x

   /// Linear interpolation; FLAT (constant) clamp outside [x.front(), x.back()].
   real_t operator()(real_t xq) const;

   /// MFEM_VERIFY: size>=2, x.size()==y.size(), x strictly increasing, all finite.
   void Validate() const;
};

/// Depth profile for rate-and-state a(z) and b(z), built from two CSV files
/// (Phase 11b).  a(depth) and (a-b)(depth) are interpolated INDEPENDENTLY (the
/// two files may use different depth grids) and combined as b = a - (a-b).
/// Named to avoid collision with the MATERIAL DepthProfile1D (Phase 10,
/// dynamic/heterogeneous_material.hpp).
struct FrictionDepthProfile1D
{
   PiecewiseLinear1D a_of_depth;     ///< a(depth_m)        (from param_a.csv)
   PiecewiseLinear1D amb_of_depth;   ///< (a-b)(depth_m)    (from param_a_minus_b.csv)
   real_t a(real_t depth_m) const { return a_of_depth(depth_m); }
   real_t b(real_t depth_m) const { return a_of_depth(depth_m) - amb_of_depth(depth_m); }
};

/// `[friction.rate_state.depth_profile]` block (Phase 11b).  When `enabled`,
/// the resolver seeds per-DOF a/b from `profile` (evaluated at depth = max(0,-z)
/// in metres) instead of the scalar a_default/b_default.
struct FrictionDepthProfileSpec
{
   bool        enabled             = false;
   std::string param_a_csv;          ///< path to param_a.csv (CWD-relative, like [mesh].path)
   std::string param_a_minus_b_csv;  ///< path to param_a_minus_b.csv
   real_t      depth_to_m          = 1000.0;  ///< depth_units="km" -> 1000; "m" -> 1
   FrictionDepthProfile1D profile;            ///< built at parse time
};

/// Load + build a FrictionDepthProfile1D from the two CSVs named in `spec`
/// (Phase 11b).  Each CSV data row is `value, depth_km` (value FIRST); depth is
/// scaled by `spec.depth_to_m` to metres.  MFEM_ABORTs on a missing file, <2
/// rows, a duplicate depth, a non-finite field, or a non-positive `a` value.
FrictionDepthProfile1D LoadFrictionDepthProfileCSVs(const FrictionDepthProfileSpec& spec);

/// Rate-and-state evolution-law selector (Phase 6 req 4).  Default AgingLaw
/// (TPV102 / SAFS).  SlipLawStrongRateWeakening is TPV104 (FVW): per-QP V_w +
/// the weakening friction f_w (= muW) feed the slip-law-SRW analytic step.
enum class StateEvolutionKind { AgingLaw, SlipLawStrongRateWeakening };

struct RateStateBlock
{
   real_t f_0_default        = 0.6;
   real_t V_0_default        = 1.0e-6;
   bool   eta_auto           = true;
   real_t eta_default        = 0.0;
   real_t a_default          = 0.010;
   real_t b_default          = 0.015;
   real_t Dc_default         = 0.004;
   real_t V_init_default     = 1.0e-9;
   real_t sigma_n_default    = 50.0e6;
   // Phase 6 req 4: SRW state-evolution.  Defaults keep existing RS configs
   // byte-identical (AgingLaw ignores f_w_default / V_w_default).
   StateEvolutionKind state_evolution = StateEvolutionKind::AgingLaw;
   real_t f_w_default        = 0.1;   ///< SRW weakening friction (TPV104 muW); SRW only
   real_t V_w_default        = 0.1;   ///< SRW weakening velocity (TPV104 V_w_in); SRW only
   std::vector<SpatialRule>  spatial;
   FrictionDepthProfileSpec  depth_profile;   // Phase 11b: depth-varying a/b (optional)
};

// =====================================================================
//  Phase 6 req 1: TPV problem / boundary / fault-frame / hypocenter /
//  material config blocks.  All optional (defaults preserve existing SAFS
//  behaviour); req-7 guards validate them when present.
// =====================================================================

/// Descriptive problem tag (informational only; no code branches on it).
struct ProblemSpec { std::string tag; };

/// Boundary-attribute assignment.  `fault_attr` is the mesh attribute of the
/// fault interface; natural/absorbing list the outer-boundary attributes
/// treated as free-surface / absorbing.  The three must be disjoint and
/// fault_attr > 0 (req-7 guards).
struct BoundarySpec
{
   int              fault_attr = -1;
   std::vector<int> natural_attrs;
   std::vector<int> absorbing_attrs;
};

/// Fault-local frame: `ref_normal` / `up` define the canonical
/// (n, t1=dip, t2=strike) basis.  Defaults match the SAFS/TPV y=0 vertical
/// strike-slip fault.  Both must be unit-norm and non-parallel (req-7 guards).
struct FaultGeometrySpec
{
   std::array<real_t, 3> ref_normal = {{0.0, -1.0, 0.0}};
   std::array<real_t, 3> up         = {{0.0,  0.0, 1.0}};
   std::string           kind;   ///< informational (no code branches on it)
};

/// Hypocenter / nucleation centre.  R-008 guard: if up[2] > 0 (z increases
/// upward) the hypocenter z must be <= 0 (at/below the free surface).
struct HypocenterSpec
{
   real_t x_m = 0.0, y_m = 0.0, z_m = 0.0;
   real_t nucleation_radius_m = 0.0;
   real_t nucleation_taper_m  = 0.0;
};

/// Bulk-material model selector (req 1).  `matrix` interior_flux requires a
/// non-Constant material (the deferred req-3 guard, completed here).
///   Constant       -> [material_constant_fallback] (lambda, mu, rho)
///   DepthProfile1D  -> a 1-D depth profile (TPV31; built in Phase 10)
///   SidecarHDF5     -> CVM-H/CVM-S velocity sidecar
enum class MaterialKind { Constant, DepthProfile1D, SidecarHDF5 };

struct MaterialSpec
{
   MaterialKind kind = MaterialKind::Constant;
   std::string  profile_csv;    ///< DepthProfile1D path (built in Phase 10)
   std::string  sidecar_path;   ///< SidecarHDF5 velocity-model path
};

struct SpatialFrictionConfig
{
   int                                 schema_version = 0;
   FrictionLawKind                     law = FrictionLawKind::SlipWeakening;
   std::string                         description;
   /// Compressive normal-stress floor [Pa] applied to σ_n in the shear
   /// strength only (σ_n cap on strength; see the sliver-blowup debug
   /// plan 2026-05-26).  `< 0` is the disabled sentinel ⇒ each friction
   /// law keeps its exact current expression (LSW `max(σ_n,0)`, RS
   /// `|σ_n|`) so the TPV/BP5 byte-exact regressions are untouched.
   /// `>= 0` floors the strength's σ_n at this value, breaking the
   /// σ_n→strength→radiation feedback that drives the tensile free-slip
   /// runaway.  Parsed from `[friction].sigma_n_strength_floor_pa`.
   real_t                              sigma_n_strength_floor_pa = -1.0;
   MaterialConstantFallback            material_fallback;
   MeshSpec                            mesh;
   VelocitySpec                        velocity;
   StressSpec                          stress;
   NumericsSpec                        numerics;
   TimeSpec                            time;
   OutputSpec                          output;
   NucleationSpec                      nucleation;       // D-4
   std::optional<SlipWeakeningBlock>   slip_weakening;
   std::optional<RateStateBlock>       rate_state;
   // Phase 6 req 1: optional TPV config blocks (defaults preserve SAFS).
   ProblemSpec                         problem;
   BoundarySpec                        boundary;
   FaultGeometrySpec                   fault_geometry;
   HypocenterSpec                      hypocenter;
   MaterialSpec                        material;
};

/// Parse + validate a TOML config.  Aborts on any schema violation with
/// a precise message; never returns a half-validated SpatialFrictionConfig.
///
/// Requires SEAS_USE_TOML at compile time (extern/toml11/toml.hpp must
/// be available).  Without it, this function aborts at runtime with a
/// "rebuild with toml11" message.
SpatialFrictionConfig LoadSpatialFrictionConfig(const std::string& toml_path);

/// Same parsing logic but reads from an in-memory TOML string.  Used by
/// the unit tests so they do not depend on any on-disk file path.
SpatialFrictionConfig ParseSpatialFrictionConfigString(const std::string& toml_text);

// =====================================================================
//  Resolved per-DOF parameter vectors
// =====================================================================

struct SlipWeakeningPerDOFParams
{
   Vector mu_s, mu_d, d_c, cohesion;
};

struct RateStatePerDOFParams
{
   Vector a, b, Dc, V_init, f_0, V_0, eta, sigma_n_eff;
   Vector V_w;   ///< Phase 6 req 4: per-DOF SRW weakening velocity (SRW only;
                 ///< filled to V_w_default / per-rule V_w like `a`).
};

class SpatialFrictionResolver
{
public:
   SlipWeakeningPerDOFParams ResolveSlipWeakening(
      const SlipWeakeningBlock& cfg,
      const Vector&             dof_coords_3d,
      const Array<int>&         dof_to_attr) const;

   RateStatePerDOFParams ResolveRateState(
      const RateStateBlock&     cfg,
      const Vector&             dof_coords_3d,
      const Array<int>&         dof_to_elem,
      const Array<int>&         dof_to_attr,
      const MaterialField&      material,
      ParMesh&                  pmesh,
      const PorePressureSpec&   pp,
      const Vector&             sigma_n_total_per_dof) const;

   /// Serial-mesh overload (unit tests).  Identical semantics to the
   /// ParMesh path.
   RateStatePerDOFParams ResolveRateState(
      const RateStateBlock&     cfg,
      const Vector&             dof_coords_3d,
      const Array<int>&         dof_to_elem,
      const Array<int>&         dof_to_attr,
      const MaterialField&      material,
      mfem::Mesh&               mesh,
      const PorePressureSpec&   pp,
      const Vector&             sigma_n_total_per_dof) const;
};

// =====================================================================
//  TPV26/27 §Part 4 friction-coefficient helper (Phase H consumer)
// =====================================================================

/// Effective LSW friction coefficient at slip magnitude `delta` and
/// simulation time `t_now`, with TPV26/27 gradual forced rupture.
///
///   μ = μ_s + (μ_d − μ_s) · max(f_1(δ), f_2(t, T_forced, t_0))
///   f_1(δ) = min(δ / d_c, 1)               (standard slip-weakening)
///   f_2(t) = 0                            if t < T_forced
///          = (t - T_forced) / t_0         if T_forced <= t < T_forced + t_0
///          = 1                            if t >= T_forced + t_0
///
/// When `T_forced >= 1e8` (the "never forced" sentinel set on TPV205
/// DOFData defaults), `f_2 = 0` for every reachable simulation time,
/// so `μ` reduces to the plain TPV205 LSW formula `μ_s - (μ_s - μ_d) ·
/// min(δ/d_c, 1)` byte-identically.  The TPV205 byte-exact contract
/// (§TPV Benchmark Isolation Phase H) relies on this.
///
/// The strength-barrier short-circuit (R-002 in dynamic/tpv205_friction.hpp)
/// is mirrored here: when `mu_s >= 0.5 * mu_s_barrier`, the function
/// returns `mu_s` immediately.  `mu_s_barrier` defaults to `1.0e6`
/// (the existing sentinel used by `LSWFrictionCoefficient_TPV205`).
inline real_t LSWFrictionCoefficient_ForcedRupture(real_t delta,
                                                   real_t mu_s,
                                                   real_t mu_d,
                                                   real_t d_c,
                                                   real_t t_now,
                                                   real_t T_forced,
                                                   real_t t0_decay,
                                                   real_t mu_s_barrier = 1.0e6)
{
   // Barrier short-circuit: locked DOF, no slip-weakening, no forced rupture.
   if (mu_s >= 0.5 * mu_s_barrier) { return mu_s; }

   // f_1: standard slip-weakening factor.
   real_t f1;
   if (delta <= 0.0)
   {
      f1 = 0.0;
   }
   else if (delta >= d_c)
   {
      f1 = 1.0;
   }
   else
   {
      f1 = delta / d_c;
   }

   // f_2: forced-rupture factor.  First branch protects against
   // t0_decay == 0 (TPV205 default): when T_forced = 1e9 and t0_decay
   // = 0, the "t < T_forced" branch always fires for any reachable
   // simulation time, so f_2 = 0 with no division by zero.
   real_t f2;
   if (t_now < T_forced)
   {
      f2 = 0.0;
   }
   else if (t_now < T_forced + t0_decay)
   {
      // t0_decay > 0 is guaranteed by this branch (the "t < T_forced"
      // branch above already covered the case t0_decay == 0 because
      // T_forced + 0 == T_forced and t_now >= T_forced ⇒ t_now >=
      // T_forced + t0_decay ⇒ this branch is not entered).
      f2 = (t_now - T_forced) / t0_decay;
   }
   else
   {
      f2 = 1.0;
   }

   const real_t factor = (f1 > f2) ? f1 : f2;
   return mu_s + (mu_d - mu_s) * factor;
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif // MFEM_SEAS_SPATIAL_FRICTION_HPP
