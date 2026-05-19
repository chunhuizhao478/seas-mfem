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
#include "../../dynamic/spatial_nucleation.hpp"   // GradualOverstressSpec, SquareOverstressSpec

#include <array>
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

// =====================================================================
//  Phase R.3 new top-level sections (plan §R.3 Detailed Requirements 1
//  + plan §"Constraints" Convention constraints).
// =====================================================================

/// `[problem]` — free-form benchmark tag written into checkpoints'
/// `driver_tag` field and printed at startup.  Human-use only; no
/// code branches on this value.  Default "safs" preserves the
/// existing SAFS hardcoded behaviour.
struct ProblemSpec
{
   std::string tag = "safs";
};

/// `[boundary]` — mesh-attribute → BC-class mapping.  Single fault
/// attribute; lists for natural (free-surface / traction-free) and
/// absorbing.  Defaults match the existing SAFS hardcoded values:
/// fault=101, natural=[102], absorbing=[103, 104].
struct BoundarySpec
{
   int              fault_attr      = 101;
   std::vector<int> natural_attrs   = {102};
   std::vector<int> absorbing_attrs = {103, 104};
};

/// `[fault_geometry]` — FaultBasis seed orientation + FaultGeometry
/// ctor / friction-law dispatch selector.
///
/// CANONICAL SEAS COORDINATE CONVENTION (REVIEW R-003 fix):
/// Defaults match the Tandem-convention vertical y=0 fault used by
/// TPV102 / TPV104 / TPV205 and the wave operator's INTERNAL FaultBasis
/// at `dynamic/wave_operator.inl:349-350`:
///     ref_normal = (0, -1, 0)
///     up         = (0,  0, 1)
/// These values are project-wide canonical for any planar y=0 fault.
/// The spatial driver hard-fails at startup if a TOML overrides
/// ref_normal/up to a different vector, because the wave operator's
/// internal FaultBasis is hard-coded; an external override would put
/// pre-stress (`tau_pre`) in a different fault-local frame from the
/// runtime trial traction (REVIEW R-003).  See
/// `miniapps/seas/CLAUDE.md` "Canonical Coordinate System".
///
/// `kind` values:
///   * "bp5_safs"  — existing SAFS Phase 5a ctor with BP5Params seed.
///   * "tpv205_lsw" — TPV205-style: planar fault, LSW, hypocenter-
///                    centred strength-reduction nucleation.
///   * "tpv31_lsw"  — TPV31-style: planar fault, 1D bi-material,
///                    LSW with depth-dependent cohesion taper +
///                    tau_nuke circular zone (overstress nucleation).
///                    NOTE: TPV31's spec uses a different fault
///                    geometry (fault on z=0 with y as depth axis).
///                    Until TPV31 is rotated into the canonical
///                    y=0/depth=-z frame (follow-up task), TPV31
///                    configs explicitly override ref_normal/up and
///                    will trip the R-003 guard in the spatial driver.
struct FaultGeometrySpec
{
   std::array<real_t, 3> ref_normal = { 0.0, -1.0, 0.0 };
   std::array<real_t, 3> up         = { 0.0,  0.0, 1.0 };
   std::string           kind       = "bp5_safs";
};

/// `[hypocenter]` — nucleation patch geometry.  All units SI (metres).
/// `nucleation_taper_m` is the linear cosine taper width beyond
/// `nucleation_radius_m` over which the nucleation perturbation
/// (StrengthReduction for TPV205, Overstress for TPV31) ramps to 0.
struct HypocenterSpec
{
   real_t x                   = 0.0;
   real_t y                   = 7500.0;
   real_t z                   = 0.0;
   real_t nucleation_radius_m = 1400.0;
   real_t nucleation_taper_m  = 600.0;   // TPV31 uses 1400→2000 (taper=600)
};

/// `[material]` kind selector — picks how (lambda, mu, rho) varies
/// across the mesh.  Defaults to "constant" (uses
/// `MaterialConstantFallback` from the existing
/// `[material_constant_fallback]` block).
enum class MaterialKind
{
   Constant,           ///< MaterialField::MakeConstant from fallback
   DepthProfile1D,     ///< Piecewise-1D profile parsed below
   SidecarHDF5         ///< Existing sidecar bundle path
};

struct MaterialSpec
{
   MaterialKind                  kind         = MaterialKind::Constant;
   std::vector<DepthProfileLayer> profile_layers;   // only when kind=DepthProfile1D
   char                          depth_axis   = 'y'; // TPV31 spec p. 3
};

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

struct NumericsSpec
{
   int         ader_order = 2;
   std::string mixed_flux = "none";
   real_t      cfl        = 0.5;
   bool        use_pml    = false;

   /// CFL safety-factor convention (REVIEW R-004).
   ///   "raw"  — pass `cfl` directly to `WaveOperator::ComputeMaxDt`
   ///            (SAFS legacy behaviour; dt = cfl_mixed_flux_factor *
   ///             cfl * h_min / cp).
   ///   "dg"   — pre-scale `cfl` by the DG safety factor
   ///            1 / (3*(2*order+1)) before calling ComputeMaxDt; this
   ///            matches the native TPV102/TPV104/TPV205 drivers
   ///            (cfl_native = cfl / 9 for order=1).  Required for
   ///            byte parity with the native driver outputs.
   /// Default = "raw" preserves existing SAFS configs.
   std::string cfl_safety = "raw";

   /// Sub-step iterator quadrature (REVIEW R-005).
   ///   "one-shot" — O = 1 sub-step per macro-step (matches native
   ///                TPV205 default `--fault-iterator one-shot`).
   ///   "substep"  — O = ader_order sub-steps per macro-step (per-
   ///                sub-step ADER quadrature for tighter friction
   ///                resolution at higher cost).
   /// Default = "one-shot" for native TPV205 parity.
   std::string fault_iterator = "one-shot";

   /// Interior-face flux dispatch (REVIEW R-006).
   ///   "bimaterial" — use the heterogeneous WaveOperator ctor;
   ///                  every interior face runs the exact bi-material
   ///                  Riemann solver (Phase R).  Even Mode::Constant
   ///                  input routes through this path, with ~1e-12
   ///                  relative drift from the scalar ctor.  Default;
   ///                  required for any sidecar / depth-profile / etc.
   ///                  spatially-varying material.
   ///   "scalar"     — use the scalar WaveOperator(λ, μ, ρ) ctor;
   ///                  byte-identical to native TPV102/TPV104/TPV205
   ///                  for homogeneous material.  Aborts if
   ///                  cfg.material.kind != Constant.
   std::string interior_flux = "bimaterial";
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
};

enum class StressSourceKind
{
   ConstantTensor,
   SidecarHDF5,
   DepthProportionalToShearModulus,   // Phase R.5 (TPV31)
   ConstantTensorWithPatches          // TPV205-style static square patches
};

/// @brief Per-patch override of the background constant tensor for
/// `StressSourceKind::ConstantTensorWithPatches`.
///
/// A DOF at physical coord (x, y, z) is INSIDE the patch when
/// `|x - center_x_m| <= half_x_m` AND likewise for y, z; missing/NaN
/// `half_*_m` means "no constraint along that axis" (use `+inf`
/// internally), so a 2D square patch on a planar y=0 fault is encoded
/// by setting `half_x_m`/`half_z_m` to 1500 m and leaving `half_y_m`
/// unset.
///
/// Each `sigma_*_pa` is a NaN-sentinel override: NaN ⇒ inherit the
/// component from the background tensor; finite value ⇒ replace that
/// component inside the patch.  TPV205 patches override `sigma_xy_pa`
/// (along-strike shear) only; the normal-stress component stays at
/// 120 MPa from the background.
///
/// Multiple patches may cover the same DOF; the LAST patch in the
/// `[[stress.patch]]` array wins on per-component overrides (matches
/// the spatial-friction-rule convention: "document order; last-match
/// wins per key").
struct StressPatch
{
   real_t center_x_m = std::numeric_limits<real_t>::quiet_NaN();
   real_t center_y_m = std::numeric_limits<real_t>::quiet_NaN();
   real_t center_z_m = std::numeric_limits<real_t>::quiet_NaN();
   real_t half_x_m   = std::numeric_limits<real_t>::infinity();
   real_t half_y_m   = std::numeric_limits<real_t>::infinity();
   real_t half_z_m   = std::numeric_limits<real_t>::infinity();

   real_t sigma_xx_pa = std::numeric_limits<real_t>::quiet_NaN();
   real_t sigma_yy_pa = std::numeric_limits<real_t>::quiet_NaN();
   real_t sigma_zz_pa = std::numeric_limits<real_t>::quiet_NaN();
   real_t sigma_xy_pa = std::numeric_limits<real_t>::quiet_NaN();
   real_t sigma_yz_pa = std::numeric_limits<real_t>::quiet_NaN();
   real_t sigma_xz_pa = std::numeric_limits<real_t>::quiet_NaN();

   /// True iff (x, y, z) lies inside all 3 half-width bounds.
   bool contains(real_t x, real_t y, real_t z) const
   {
      return std::abs(x - center_x_m) <= half_x_m
          && std::abs(y - center_y_m) <= half_y_m
          && std::abs(z - center_z_m) <= half_z_m;
   }
};

/// Phase R.5: TPV31 depth-proportional pre-stress.  Per spec p. 6, the
/// six stress components are stored at the reference shear modulus
/// `mu_ref_pa` (= 32.038e9 Pa for TPV5/TPV31 layered profile) and
/// scaled at each fault DOF by the LOCAL `mu(y) / mu_ref`.  The
/// driver applies this scaling at fault-DOF init time via the
/// `DepthProportionalToShearModulusStressSource` in spatial_stress.cpp.
struct DepthProportionalStressSpec
{
   real_t sigma_xx_per_mu = 0.0;   // MPa (TOML), converted to Pa internally
   real_t sigma_yy_per_mu = 0.0;
   real_t sigma_zz_per_mu = 0.0;
   real_t sigma_xy_per_mu = 0.0;
   real_t sigma_yz_per_mu = 0.0;
   real_t sigma_xz_per_mu = 0.0;
   real_t mu_ref_pa       = 32.03812032e9;  // TPV5/TPV31 spec
};

struct StressSpec
{
   StressSourceKind kind = StressSourceKind::ConstantTensor;
   // ConstantTensor (Phase 3b) — also reused as the BACKGROUND tensor
   // for ConstantTensorWithPatches.
   real_t sigma_xx_pa = 0.0;
   real_t sigma_yy_pa = 0.0;
   real_t sigma_zz_pa = 0.0;
   real_t sigma_xy_pa = 0.0;
   real_t sigma_yz_pa = 0.0;
   real_t sigma_xz_pa = 0.0;
   // SidecarHDF5 (Phase 3):
   std::string sidecar_path;
   // DepthProportionalToShearModulus (Phase R.5 / TPV31):
   DepthProportionalStressSpec depth_proportional;
   // ConstantTensorWithPatches — zero or more static square patches
   // applied on top of the background sigma_*_pa tensor.  See
   // StressPatch.  TPV205 uses 3 patches.
   std::vector<StressPatch> patches;
   // Common:
   PorePressureSpec pore_pressure;
};

/// Spatial rule (D-3 relaxed validator; R-014 explicit defaults).
/// D-4 (rev-3): the `NucleationBox` kind is REMOVED.  Spatial rules
/// never override `tau_pre_*` or `sigma_n` — pre-stress comes entirely
/// from the stress source.
struct SpatialRule
{
   enum class Kind { Depth, Box, RegionAttribute, Barrier, BoxcarTaper };
   Kind kind = Kind::Depth;

   real_t x_min_m = -std::numeric_limits<real_t>::infinity();
   real_t x_max_m =  std::numeric_limits<real_t>::infinity();
   real_t y_min_m = -std::numeric_limits<real_t>::infinity();
   real_t y_max_m =  std::numeric_limits<real_t>::infinity();
   real_t z_min_m = -std::numeric_limits<real_t>::infinity();
   real_t z_max_m =  std::numeric_limits<real_t>::infinity();
   int    region_attr = -1;

   // BoxcarTaper geometry (SCEC TPV101/102/104 Eq. (5) C∞ tanh boxcar).
   //   B(x; W, w) = 1                       for |x − c| ≤ W
   //              = 0.5(1 + tanh(w/(s−W−w) + w/(s−W)))
   //                                          where s = |x − c|
   //                                        for W < |x − c| < W + w
   //              = 0                       for |x − c| ≥ W + w
   // The rule applies a product `B(x) · B(y) · B(z)` along the
   // constrained axes (any axis with `boxcar_half_*_m == NaN` is
   // treated as B = 1, i.e., unconstrained on that axis).  Per
   // parameter (a, b, Dc, V_init, V_w, mu_s, ...) the override blends
   // the prior baseline into the override using this product:
   //   value = baseline + (override − baseline) · B_total
   // so B_total = 1 ⇒ value = override (inside) and B_total = 0 ⇒
   // value = baseline (outside).  When `kind != BoxcarTaper`, every
   // `boxcar_*_m` is NaN/+inf and ignored.
   real_t boxcar_center_x_m = 0.0;
   real_t boxcar_center_y_m = 0.0;
   real_t boxcar_center_z_m = 0.0;
   real_t boxcar_half_x_m   = std::numeric_limits<real_t>::quiet_NaN();
   real_t boxcar_half_y_m   = std::numeric_limits<real_t>::quiet_NaN();
   real_t boxcar_half_z_m   = std::numeric_limits<real_t>::quiet_NaN();
   real_t boxcar_trans_x_m  = std::numeric_limits<real_t>::quiet_NaN();
   real_t boxcar_trans_y_m  = std::numeric_limits<real_t>::quiet_NaN();
   real_t boxcar_trans_z_m  = std::numeric_limits<real_t>::quiet_NaN();

   // Per-key overrides.  NaN sentinel means "do not override".
   real_t mu_s    = std::numeric_limits<real_t>::quiet_NaN();
   real_t mu_d    = std::numeric_limits<real_t>::quiet_NaN();
   real_t d_c     = std::numeric_limits<real_t>::quiet_NaN();
   real_t cohesion= std::numeric_limits<real_t>::quiet_NaN();

   // Depth-linear cohesion taper (TPV31-style: C0(y) = max(0, floor +
   // grad * (ref_depth - y))).  When `cohesion_grad_pa_per_m` is NaN
   // the linear taper is disabled and the `cohesion` override above
   // applies as a constant.  When set, the linear taper takes
   // precedence over `cohesion` (which is then ignored).
   //
   // Spec convention: positive `grad` means cohesion DECREASES with
   // depth (`y` increasing).  For TPV31:
   //   ref_depth = 2400 m, floor = 0 Pa, grad = 425 Pa/m
   //   → C0(y) = max(0, 0 + 425 · (2400 - y))
   //   → C0(0) = 1.02 MPa; C0(2400) = 0; C0(y>2400) = 0
   real_t cohesion_grad_pa_per_m  = std::numeric_limits<real_t>::quiet_NaN();
   real_t cohesion_ref_depth_m    = std::numeric_limits<real_t>::quiet_NaN();
   real_t cohesion_floor_pa       = 0.0;
   char   cohesion_taper_axis     = 'y';   // 'x' | 'y' | 'z'

   real_t a       = std::numeric_limits<real_t>::quiet_NaN();
   real_t b       = std::numeric_limits<real_t>::quiet_NaN();
   real_t Dc      = std::numeric_limits<real_t>::quiet_NaN();
   real_t V_init  = std::numeric_limits<real_t>::quiet_NaN();
   real_t f_0     = std::numeric_limits<real_t>::quiet_NaN();
   real_t V_0     = std::numeric_limits<real_t>::quiet_NaN();
   real_t sigma_n = std::numeric_limits<real_t>::quiet_NaN();
   real_t eta     = std::numeric_limits<real_t>::quiet_NaN();
   real_t V_w     = std::numeric_limits<real_t>::quiet_NaN();  ///< SRW only

   /// `matches` semantics per §Phase 1 Detailed Req. 4:
   ///   Depth            : only z bounds checked.
   ///   Box              : all 6 coord bounds checked.
   ///   RegionAttribute  : only region_attr checked.
   ///   Barrier          : same as Depth, plus any non-infinite
   ///                      x/y bounds also apply.
   ///   BoxcarTaper      : always matches; the C∞ boxcar product
   ///                      (computed by `BoxcarTaperFactor`) is the
   ///                      blend weight applied by the resolver.
   bool matches(real_t x, real_t y, real_t z, int attr) const;

   /// @brief SCEC C∞ Boxcar product B(x)·B(y)·B(z) for this rule.
   /// Axes with `boxcar_half_*_m == NaN` are treated as B = 1.  Each
   /// constrained axis requires both `boxcar_half_*_m` and
   /// `boxcar_trans_*_m` to be finite and ≥ 0 (validated by the
   /// parser).  Returns a value in `[0, 1]`.
   real_t BoxcarTaperFactor(real_t x, real_t y, real_t z) const;
};

/// @brief Stand-alone SCEC C∞ Boxcar function B(|s − center|; W, w).
/// Generic helper exposed for unit tests and re-use; the SpatialRule
/// `BoxcarTaperFactor` is a 3-axis product of this scalar.
real_t SCECBoxcar(real_t signed_offset, real_t half_width, real_t transition);

/// @brief Phase N — nucleation mechanisms for the spatial driver:
///   * `GradualOverstress` (Phase N): per-DOF Δτ accumulator,
///     smoothStep-ramped, Gaussian-shaped.
///   * `SquareOverstress` (TPV5/TPV205-family extension): per-DOF Δτ
///     accumulator, smoothStep-ramped, indicator-function (sharp
///     rectangular) shape with one or more disjoint patches.
///   * `InstantaneousOverstressCircular` (TPV31): per-DOF Δτ written
///     ONCE at init (no ramp), circular cosine-tapered shape, per-DOF
///     µ-scaling.
///
/// The obsolete `StrengthReduction` (TPV26/27 §Part 4) and `Overstress`
/// (TPV205-style one-shot) kinds were removed from the spatial code
/// path.  Native TPV* drivers keep their own paths via
/// `dynamic/tpv104_nucleation.hpp` + the still-live
/// `FaultFrictionLaw::LSW_ForcedRupture` flux dispatch.
enum class NucleationKind
{
   GradualOverstress,
   SquareOverstress,
   InstantaneousOverstressCircular,
   GradualOverstressCompactCircular   ///< SCEC TPV101/102/104 spec Eq. (13)
};

/// @brief Stringify a NucleationKind for banners / logs.  Returns
/// pointer to a static C string; safe to embed in std::string.
/// (R-010 of tpv102_tpv104_review.md: the spatial driver's startup
/// banner must echo the actual kind, not a hard-coded "gradual_overstress".)
inline const char* NucleationKindToString(NucleationKind k)
{
   switch (k)
   {
   case NucleationKind::GradualOverstress:
      return "gradual_overstress";
   case NucleationKind::SquareOverstress:
      return "square_overstress";
   case NucleationKind::InstantaneousOverstressCircular:
      return "instantaneous_overstress_circular";
   case NucleationKind::GradualOverstressCompactCircular:
      return "gradual_overstress_compact_circular";
   }
   return "unknown";
}

/// @brief Top-level `[nucleation]` TOML block.  When absent in the
/// TOML, `enabled == false` and the driver runs with no nucleation
/// perturbation (tau{1,2}_nuc stay at the init-time zero seeded by
/// `InitializeFaultDOFs_Spatial`).
struct NucleationSpec
{
   NucleationKind                              kind = NucleationKind::GradualOverstress;
   bool                                        enabled = false;
   GradualOverstressSpec                       gradual_overstress;
   SquareOverstressSpec                        square_overstress;
   InstantaneousOverstressCircularSpec         instantaneous_overstress_circular;
   GradualOverstressCompactCircularSpec        gradual_overstress_compact_circular;
};

struct SlipWeakeningBlock
{
   real_t mu_s_default     = 1.1;   // geoffrey2010.md
   real_t mu_d_default     = 0.5;
   real_t d_c_default      = 0.5;
   real_t cohesion_default = 0.0;
   std::vector<SpatialRule> spatial;
};

/// State evolution choice for rate-and-state friction.  `AgingLaw`
/// (Dieterich–Ruina aging) is the default and what BP5 / TPV102 use.
/// `SlipLawStrongRateWeakening` (FVW / SCEC FL=103) adds two parameters
/// `f_w` and `V_w` consumed by `SlipLawSRWPsi` — see `dynamic/`
/// `tpv104_substep_iterator.{hpp,cpp}` and `friction/slip_law_srw_psi.hpp`.
enum class StateEvolutionKind { AgingLaw, SlipLawStrongRateWeakening };

struct RateStateBlock
{
   StateEvolutionKind        state_evolution = StateEvolutionKind::AgingLaw;
   real_t f_0_default        = 0.6;
   real_t V_0_default        = 1.0e-6;
   bool   eta_auto           = true;
   real_t eta_default        = 0.0;
   real_t a_default          = 0.010;
   real_t b_default          = 0.015;
   real_t Dc_default         = 0.004;
   real_t V_init_default     = 1.0e-9;
   real_t sigma_n_default    = 50.0e6;
   // Slip-law-strong-rate-weakening parameters (only consumed when
   // state_evolution == SlipLawStrongRateWeakening).
   //
   // R-007 (tpv102_tpv104_review.md): under the AgingLaw branch
   // (TPV102), these fields are read from TOML but are NEVER consumed
   // by the friction solver — do not assume they are honoured at
   // aging-law DOFs.  A future refactor that exposes `f_w` to the
   // aging-law path must reinstate the per-DOF validator currently
   // gated on `state_evolution == SlipLawStrongRateWeakening` (see
   // `spatial_friction.cpp::resolve_rs_impl`).  TPV102 toml may omit
   // these knobs; TPV104 toml must set them to the FL=103 spec
   // (`f_w = 0.2`, `V_w_default = 1.0` outside the VW core).
   real_t f_w_default        = 0.2;     ///< weakening-friction coefficient
   real_t V_w_default        = 1.0;     ///< weakening velocity [m/s]
   std::vector<SpatialRule>  spatial;
};

struct SpatialFrictionConfig
{
   int                                 schema_version = 0;
   FrictionLawKind                     law = FrictionLawKind::SlipWeakening;
   std::string                         description;
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

   // Phase R.3 additions — every section defaults to current SAFS
   // hardcoded behaviour for backwards compatibility.
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
   /// V_w(i) is populated only when state_evolution ==
   /// SlipLawStrongRateWeakening.  Zero-sized otherwise.
   Vector V_w;
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
