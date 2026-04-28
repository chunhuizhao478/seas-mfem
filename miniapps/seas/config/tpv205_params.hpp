// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV5 / TPV205 benchmark parameters (linear slip-weakening friction
// with multiple initial-shear-stress patches).
//
// Reference: SCEC TPV5 description (TPV5_forwebsite.pdf, Aug 4 2005)
// stored in miniapps/seas/tpv205/benchmark_document/.
//
// TPV205 differs from TPV104 in two places:
//   1. Friction law — TPV205 uses LINEAR SLIP-WEAKENING (LSW):
//        μ(δ) = μ_s − (μ_s − μ_d) · min(δ/d_c, 1)
//        τ_strength = μ(δ) · σ_n
//      No state variable, no Newton iteration; the slip-rate magnitude is
//      a CLOSED-FORM result of the radiation-damping balance:
//        V_abs = max(0, (|τ_total| − τ_strength) / η_s)
//   2. Nucleation — TPV205 nucleates by SPATIALLY-VARYING INITIAL SHEAR
//      STRESS (4 patches of 3 km × 3 km).  No time-varying perturbation;
//      the τ2_0 field carries patch-dependent values from t = 0.
//
// All other infrastructure (mesh handling, ADER predictor/corrector,
// mixed-flux dispatch, ParaView, station probes) is mode-independent and
// reused as-is from TPV104's code flow.

#ifndef MFEM_SEAS_TPV205_PARAMS_HPP
#define MFEM_SEAS_TPV205_PARAMS_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// SCEC TPV205 material, friction (LSW), and pre-stress patch parameters.
struct TPV205Params
{
   // ====== Material (homogeneous, isotropic) =============================
   static constexpr real_t rho = 2670.0;        ///< Density [kg/m³]
   static constexpr real_t cs  = 3464.0;        ///< S-wave speed [m/s]
   static constexpr real_t cp  = 6000.0;        ///< P-wave speed [m/s]
   static constexpr real_t mu  = rho * cs * cs; ///< Shear modulus [Pa]
   static constexpr real_t lambda = rho * cp * cp - 2.0 * mu;

   static constexpr real_t Zp    = rho * cp;
   static constexpr real_t Zs    = rho * cs;
   static constexpr real_t eta_s = Zs / 2.0;
   static constexpr real_t eta_p = Zp / 2.0;

   // ====== Linear slip-weakening (LSW) friction ==========================
   // SCEC TPV5 Eq. set in §10 of TPV5_forwebsite.pdf.
   //   Inside the 30 km × 15 km rupture area:
   //     μ_s = 0.677, μ_d = 0.525, d_c = 0.40 m
   //   Outside (strength barrier, §11):
   //     μ_s = 10000  (so the rupture cannot escape the rupture area),
   //     μ_d = 0.525, d_c = 0.40 m  (same dynamic friction & d_c)
   static constexpr real_t mu_s        = 0.677;     ///< Static friction (inside)
   static constexpr real_t mu_d        = 0.525;     ///< Dynamic friction
   static constexpr real_t d_c         = 0.40;      ///< Slip-weakening crit. distance [m]
   static constexpr real_t mu_s_barrier = 10000.0;  ///< Static friction (outside / barrier)

   // ====== Rupture area + barrier geometry ===============================
   // §3 of the spec: 30000 m × 15000 m rupture area; bottom + left + right
   // ends bounded by strength barrier.  Top boundary is the free surface
   // (z = 0).  In the codebase's coordinate system: x = along-strike,
   // |z| = down-dip depth; the rupture area extends
   //   |x| < 15 km and 0 ≤ depth ≤ 15 km.
   static constexpr real_t rupture_along_strike_half = 15.0e3;  ///< |x| < 15 km
   static constexpr real_t rupture_depth             = 15.0e3;  ///< 0 ≤ depth ≤ 15 km

   // ====== Pre-stress and patches ========================================
   // §10 of the spec (background) + §7-9 (three 3 km × 3 km patches).
   // The patches are CENTRED at:
   //   nucleation patch (§7) : (0, 7.5 km)            — high-stress trigger
   //   left  patch     (§9) : (-7.5 km, 7.5 km)       — also high stress (78 MPa)
   //   right patch     (§8) : (+7.5 km, 7.5 km)       — low stress (62 MPa)
   // Each patch has side-length 3 km, so it spans ±1.5 km around its centre
   // in both along-strike and along-dip.
   static constexpr real_t sigma_n   = 120.0e6;     ///< Initial normal stress [Pa]
   static constexpr real_t tau_back  = 70.0e6;      ///< Background τ_strike [Pa]
   static constexpr real_t tau_nuc   = 81.6e6;      ///< Nucleation patch τ_strike [Pa]
   static constexpr real_t tau_left  = 78.0e6;      ///< Left patch τ_strike [Pa]
   static constexpr real_t tau_right = 62.0e6;      ///< Right patch τ_strike [Pa]

   // Patch geometry (centres + half-side).
   static constexpr real_t patch_half        = 1.5e3;   ///< Half side-length [m]
   static constexpr real_t hypo_along_strike = 0.0;     ///< Nucleation centre x [m]
   static constexpr real_t hypo_down_dip     = 7.5e3;   ///< Nucleation centre depth [m]
   static constexpr real_t left_patch_x      = -7.5e3;  ///< Left  patch centre x [m]
   static constexpr real_t right_patch_x     = 7.5e3;   ///< Right patch centre x [m]
   static constexpr real_t patch_z           = 7.5e3;   ///< All three patches at depth = 7.5 km

   // Initial slip rate.  TPV5 doesn't specify an initial V; the medium
   // is at rest until the (static) initial shear stress in the
   // nucleation patch exceeds the static yield stress and rupture
   // spontaneously initiates.  V_ini = 0 is safe — `SolveLSW_TPV205`
   // guards against tau_abs == 0 before normalising V along τ_total,
   // so no divide-by-zero at t = 0.
   static constexpr real_t V_ini = 0.0;             ///< Initially at rest

   // ====== Domain & simulation ===========================================
   // The mesh in tpv205/mesh/tpv2053d_*.geo extends ±60 km horizontally and
   // 60 km below the surface; the fault (Physical Surface 103) is the
   // y = 0 plane.  This driver does not enforce these — values are kept
   // here for the banner / diagnostics only.
   static constexpr real_t domain_half = 60.0e3;
   static constexpr real_t t_final     = 12.0;      ///< §III request: 0–12 s

   // ====== Boundary attribute defaults (to match the SeisSol-style geo) =
   // The TPV205 .geo files (tpv2053d_200m.geo, _100m.geo) emit:
   //   Physical Surface 101 = free surface (top)
   //   Physical Surface 103 = fault
   //   Physical Surface 105 = absorbing (sides + bottom)
   // The driver's default --bc-free / --bc-fault / --bc-absorb match.
   static constexpr int bc_free_default    = 101;
   static constexpr int bc_fault_default   = 103;
   static constexpr int bc_absorb_default  = 105;
};

// =========================================================================
// Spatial helpers — patch lookup and strength-barrier mask.
// =========================================================================

/// Return TRUE if (along_strike, down_dip) is inside one of the three
/// 3 km × 3 km square patches (nucleation, left, or right).
inline bool InAnyPatch_TPV205(real_t along_strike, real_t down_dip,
                              int *which_patch_out = nullptr)
{
   const real_t half = TPV205Params::patch_half;
   const real_t z    = down_dip;
   // Nucleation patch (centre 0, 7.5).
   if (std::abs(along_strike - TPV205Params::hypo_along_strike) <= half &&
       std::abs(z - TPV205Params::patch_z) <= half)
   {
      if (which_patch_out) { *which_patch_out = 0; }  // 0 = nucleation
      return true;
   }
   // Right patch (centre +7.5, 7.5).
   if (std::abs(along_strike - TPV205Params::right_patch_x) <= half &&
       std::abs(z - TPV205Params::patch_z) <= half)
   {
      if (which_patch_out) { *which_patch_out = 1; }  // 1 = right
      return true;
   }
   // Left patch (centre -7.5, 7.5).
   if (std::abs(along_strike - TPV205Params::left_patch_x) <= half &&
       std::abs(z - TPV205Params::patch_z) <= half)
   {
      if (which_patch_out) { *which_patch_out = 2; }  // 2 = left
      return true;
   }
   if (which_patch_out) { *which_patch_out = -1; }
   return false;
}

/// Return TRUE if (along_strike, down_dip) is inside the 30 km × 15 km
/// rupture area where rupture is allowed (§3).  Outside this area the
/// strength barrier (μ_s = 10000) is applied (§11).
inline bool InRuptureArea_TPV205(real_t along_strike, real_t down_dip)
{
   return std::abs(along_strike) <= TPV205Params::rupture_along_strike_half &&
          down_dip                <= TPV205Params::rupture_depth &&
          down_dip                >= 0.0;
}

/// Initial along-strike shear pre-stress τ_strike at this fault QP.
///   Inside nucleation patch:                 81.6 MPa
///   Inside right patch:                      62.0 MPa
///   Inside left patch:                       78.0 MPa
///   Otherwise (background, incl. barrier):   70.0 MPa
inline real_t ComputeTau2_0_TPV205(real_t along_strike, real_t down_dip)
{
   int which = -1;
   if (InAnyPatch_TPV205(along_strike, down_dip, &which))
   {
      switch (which)
      {
         case 0: return TPV205Params::tau_nuc;    // nucleation
         case 1: return TPV205Params::tau_right;  // right (62 MPa)
         case 2: return TPV205Params::tau_left;   // left  (78 MPa)
         default: break;
      }
   }
   return TPV205Params::tau_back;                  // background (70 MPa)
}

/// Static friction coefficient μ_s at this fault QP.
///   Inside the 30 km × 15 km rupture area: 0.677
///   Outside (strength barrier):            10000
inline real_t ComputeMuS_TPV205(real_t along_strike, real_t down_dip)
{
   if (InRuptureArea_TPV205(along_strike, down_dip))
   {
      return TPV205Params::mu_s;
   }
   return TPV205Params::mu_s_barrier;
}

/// Dynamic friction coefficient μ_d at this fault QP (§7-11: uniform 0.525).
inline real_t ComputeMuD_TPV205(real_t /*along_strike*/, real_t /*down_dip*/)
{
   return TPV205Params::mu_d;
}

/// Slip-weakening critical distance d_c at this fault QP (uniform 0.40 m).
inline real_t ComputeDc_TPV205(real_t /*along_strike*/, real_t /*down_dip*/)
{
   return TPV205Params::d_c;
}

// =========================================================================
// SCEC TPV205 on-fault stations (16 entries, matching the canonical
// benchmark-trace filename suffixes used by the DRDG3D reference data
// in tpv205/benchmark_data/DRDG3D_*).
// =========================================================================

struct StationTPV205
{
   real_t x2;           ///< Along-strike coordinate [m]
   real_t x3;           ///< Down-dip coordinate [m] (depth, positive)
   const char *label;   ///< SCEC trace-filename suffix (e.g. "x2_0_x3_7.5")
};

inline constexpr StationTPV205 kStationsTPV205[16] = {
   { -12.0e3,  0.0,    "x2_-12_x3_0"   },
   { -12.0e3,  7.5e3,  "x2_-12_x3_7.5" },
   {  -7.5e3,  0.0,    "x2_-7.5_x3_0"  },
   {  -7.5e3,  7.5e3,  "x2_-7.5_x3_7.5"},
   {  -4.5e3,  0.0,    "x2_-4.5_x3_0"  },
   {  -4.5e3,  7.5e3,  "x2_-4.5_x3_7.5"},
   {     0.0,  0.0,    "x2_0_x3_0"     },
   {     0.0,  3.0e3,  "x2_0_x3_3"     },
   {     0.0,  7.5e3,  "x2_0_x3_7.5"   },
   {     0.0, 12.0e3,  "x2_0_x3_12"    },
   {   4.5e3,  0.0,    "x2_4.5_x3_0"   },
   {   4.5e3,  7.5e3,  "x2_4.5_x3_7.5" },
   {   7.5e3,  0.0,    "x2_7.5_x3_0"   },
   {   7.5e3,  7.5e3,  "x2_7.5_x3_7.5" },
   {  12.0e3,  0.0,    "x2_12_x3_0"    },
   {  12.0e3,  7.5e3,  "x2_12_x3_7.5"  },
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV205_PARAMS_HPP
