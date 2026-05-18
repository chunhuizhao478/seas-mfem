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

struct NumericsSpec
{
   int         ader_order = 2;
   std::string mixed_flux = "none";
   real_t      cfl        = 0.5;
   bool        use_pml    = false;
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

enum class StressSourceKind { ConstantTensor, SidecarHDF5 };

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
   // Common:
   PorePressureSpec pore_pressure;
};

/// Spatial rule (D-3 relaxed validator; R-014 explicit defaults).
/// D-4 (rev-3): the `NucleationBox` kind is REMOVED.  Spatial rules
/// never override `tau_pre_*` or `sigma_n` — pre-stress comes entirely
/// from the stress source.
struct SpatialRule
{
   enum class Kind { Depth, Box, RegionAttribute, Barrier };
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

   /// `matches` semantics per §Phase 1 Detailed Req. 4:
   ///   Depth            : only z bounds checked.
   ///   Box              : all 6 coord bounds checked.
   ///   RegionAttribute  : only region_attr checked.
   ///   Barrier          : same as Depth, plus any non-infinite
   ///                      x/y bounds also apply.
   bool matches(real_t x, real_t y, real_t z, int attr) const;
};

/// @brief Phase N — the single nucleation mechanism for the spatial
/// driver is `gradual_overstress` (per-DOF Δτ accumulator,
/// smoothStep-ramped, Gaussian-shaped).  The obsolete `StrengthReduction`
/// (TPV26/27 §Part 4) and `Overstress` (TPV205-style one-shot) kinds were
/// removed from the spatial code path.  Native TPV* drivers keep their
/// own paths via `dynamic/tpv104_nucleation.hpp` + the still-live
/// `FaultFrictionLaw::LSW_ForcedRupture` flux dispatch.
enum class NucleationKind { GradualOverstress };

/// @brief Top-level `[nucleation]` TOML block.  When absent in the
/// TOML, `enabled == false` and the driver runs with no nucleation
/// perturbation (tau{1,2}_nuc stay at the init-time zero seeded by
/// `InitializeFaultDOFs_Spatial`).
struct NucleationSpec
{
   NucleationKind        kind     = NucleationKind::GradualOverstress;
   bool                  enabled  = false;
   GradualOverstressSpec gradual_overstress;
};

struct SlipWeakeningBlock
{
   real_t mu_s_default     = 1.1;   // geoffrey2010.md
   real_t mu_d_default     = 0.5;
   real_t d_c_default      = 0.5;
   real_t cohesion_default = 0.0;
   std::vector<SpatialRule> spatial;
};

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
