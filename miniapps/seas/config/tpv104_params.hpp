// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV104 benchmark parameters.
// References: SCEC benchmark description (strike.scec.org/cvws/tpv104docs)
//             Noda & Lapusta (2013) — slip-law with strong rate weakening
//
// R-001 (2026-04-24) Step 1 of TPV104 implementation.  This file is a
// standalone duplicate of `tpv102_params.hpp` — per
// `feedback_tpv102_bp5_no_shared_edit`, no shared helper between the two
// benchmarks so that TPV102 remains bit-identical as TPV104 evolves.
//
// §4.1 of tpv104_debug_plan_2026-04-24.md is the single source of truth
// for the numeric values below.  Phase 3 cross-verification probes
// surface any differences with the reference benchmark trace; they are
// not silently papered over.

#ifndef MFEM_SEAS_TPV104_PARAMS_HPP
#define MFEM_SEAS_TPV104_PARAMS_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// SCEC TPV104 material, friction, and loading parameters.
///
/// TPV104 differs from TPV102 in the friction law (slip law with strong
/// rate weakening — SCEC FL=103, also called "Fast Velocity
/// Weakening") and in numeric inputs:
///   - Lower critical slip distance L = 0.4 m (vs TPV102 Dc = 0.02 m).
///   - Added weakening velocity V_w (per-QP, VW core vs strengthening).
///   - Added weakening friction coefficient f_w = 0.1.
///   - Much lower initial slip rate V_ini = 1e-16 m/s (vs 1e-12 m/s).
///   - Different shear pre-stress tau_ini = 40 MPa (vs 75 MPa).
///   - Larger nucleation amplitude Δτ0 = 45 MPa (vs 25 MPa).
struct TPV104Params
{
   // ====== Material (homogeneous, isotropic) =============================
   static constexpr real_t rho = 2670.0;        ///< Density [kg/m³]
   static constexpr real_t cs  = 3464.0;        ///< S-wave speed [m/s]
   static constexpr real_t cp  = 6000.0;        ///< P-wave speed [m/s]
   static constexpr real_t mu  = rho * cs * cs; ///< Shear modulus [Pa]
   static constexpr real_t lambda = rho * cp * cp - 2.0 * mu; ///< First Lamé [Pa]

   static constexpr real_t Zp    = rho * cp;        ///< P-impedance
   static constexpr real_t Zs    = rho * cs;        ///< S-impedance
   static constexpr real_t eta_s = Zs / 2.0;        ///< Radiation damping [Pa·s/m]
   static constexpr real_t eta_p = Zp / 2.0;

   // ====== Rate-and-state friction (FVW / TPV104) =========================
   // FVW requires f0, b, V0, L, plus per-QP (a, V_w), plus global f_w.
   static constexpr real_t f0     = 0.6;            ///< Reference friction coefficient
   static constexpr real_t V0     = 1.0e-6;         ///< Reference slip rate [m/s]
   static constexpr real_t b      = 0.014;          ///< Slowness effect parameter
   static constexpr real_t L      = 0.4;            ///< Critical slip distance [m]
   static constexpr real_t f_w    = 0.1;            ///< Weakening friction coefficient
   static constexpr real_t muW    = f_w;            ///< alias used by friction coeff port

   // Direct effect `a` — spatially varying.
   static constexpr real_t a_in   = 0.01;           ///< a inside VW core
   static constexpr real_t da     = 0.01;           ///< Δa transition
   static constexpr real_t a_out  = a_in + da;      ///< a outside VW = 0.02

   // Weakening velocity `V_w` — spatially varying.
   static constexpr real_t V_w_in  = 0.1;           ///< V_w inside VW core [m/s]
   static constexpr real_t dV_w    = 0.9;           ///< ΔV_w transition [m/s]
   static constexpr real_t V_w_out = V_w_in + dV_w; ///< V_w outside VW = 1.0 [m/s]

   // ====== VW zone geometry (shared with TPV102 but values repeated here) ===
   // R-008 (review 2026-04-24): `W` is the FULL down-dip depth of the VW
   // zone; `ComputeA_TPV104` and `ComputeVw_TPV104` pass W/2 as the boxcar
   // half-width centred on the hypocenter.
   static constexpr real_t Ls = 15.0e3;             ///< Along-strike half-width [m]
   static constexpr real_t W  = 15.0e3;             ///< Down-dip FULL VW depth [m] (boxcar half-width = W/2)
   static constexpr real_t w  = 3.0e3;              ///< Transition half-width [m]
   static constexpr real_t ws = 3.0e3;              ///< Along-strike transition width [m]

   // Fault geometry (VW zone + 3 km transition on each side).
   static constexpr real_t fault_length = 36.0e3;
   static constexpr real_t fault_depth  = 18.0e3;

   // ====== Stress and loading ============================================
   static constexpr real_t sigma_n = 120.0e6;       ///< Normal stress [Pa] (positive compression)
   static constexpr real_t tau_ini = 40.0e6;        ///< Initial shear [Pa]
   static constexpr real_t V_ini   = 1.0e-16;       ///< Initial slip rate [m/s]

   // ====== Nucleation ====================================================
   static constexpr real_t hypo_along_strike = 0.0;
   static constexpr real_t hypo_down_dip     = 7.5e3;
   static constexpr real_t nuc_radius        = 3.0e3;
   static constexpr real_t nuc_dtau          = 45.0e6; ///< Δτ0 [Pa] (vs 25 MPa for TPV102)
   static constexpr real_t nuc_T             = 1.0;    ///< Rise time [s]

   // ====== Domain & simulation ===========================================
   static constexpr real_t domain_half = 60.0e3;
   static constexpr real_t t_final     = 12.0;
};

// R-001 note: the §4.1 "derived-scalar" invariants (a_out = a_in + da,
// V_w_out = V_w_in + dV_w, lambda = rho*cp^2 - 2*mu, eta_s = rho*cs/2)
// are enforced by the `= expr` definitions above — the constants are
// LITERALLY `a_in + da` etc.  Previous rounds (R-001) carried
// static_asserts that re-evaluated the same expressions and compared
// to the defined constant, which is tautological.  Intel 19.1.1 on
// Frontera trips those redundant asserts under its constexpr FP-
// contraction rules (the sum gets contracted differently at the
// definition site vs the assert site).  The asserts are intentionally
// absent here; if the derived expressions need a tolerance-based
// runtime guard, add it in ComputeA_TPV104 / ComputeVw_TPV104 or in
// TPV104 unit tests, not in constexpr-only static_assert.

// =========================================================================
// Spatial distributions (duplicated from tpv102_params.hpp per
// feedback_tpv102_bp5_no_shared_edit — both headers are standalone).
// =========================================================================

/// SCEC Boxcar B(x, W, w) — Eq. (5) of SCEC TPV101/102/104 spec.
/// B = 1 for |x| <= W, C∞ tanh transition for W < |x| < W+w, 0 beyond.
inline real_t Boxcar_TPV104(real_t x, real_t W_half, real_t w_trans)
{
   real_t ax = std::abs(x);
   if (ax <= W_half) { return 1.0; }
   if (ax >= W_half + w_trans) { return 0.0; }
   return 0.5 * (1.0 + std::tanh(w_trans / (ax - W_half - w_trans)
                                 + w_trans / (ax - W_half)));
}

/// Compute direct-effect parameter a at (along_strike, down_dip).
///
/// SCEC TPV104 Eq. (4):  a = a_in + (a_out - a_in) * (1 - B_strike * B_dip).
/// - VW core (|x| ≤ L_s AND 0 ≤ z ≤ 2W centred at z = W): a = a_in.
/// - Strengthening outside: a = a_out.
/// - Smooth C∞ transition of width w = 3 km.
///
/// The dip boxcar is centred at z_hypo = 7.5 km (= W) with half-width W/2
/// so the VW zone spans down-dip [0, 15 km].  Transition width is `w`.
inline real_t ComputeA_TPV104(real_t along_strike, real_t down_dip)
{
   real_t B_strike = Boxcar_TPV104(along_strike, TPV104Params::Ls,
                                   TPV104Params::ws);
   real_t B_dip    = Boxcar_TPV104(down_dip - TPV104Params::hypo_down_dip,
                                   TPV104Params::W / 2.0, TPV104Params::w);
   real_t B = B_strike * B_dip;
   return TPV104Params::a_out + (TPV104Params::a_in - TPV104Params::a_out) * B;
}

/// Compute weakening-velocity V_w at (along_strike, down_dip).
/// Same boxcar geometry as `ComputeA_TPV104` but interpolating between
/// V_w_in (VW core) and V_w_out (strengthening).
inline real_t ComputeVw_TPV104(real_t along_strike, real_t down_dip)
{
   real_t B_strike = Boxcar_TPV104(along_strike, TPV104Params::Ls,
                                   TPV104Params::ws);
   real_t B_dip    = Boxcar_TPV104(down_dip - TPV104Params::hypo_down_dip,
                                   TPV104Params::W / 2.0, TPV104Params::w);
   real_t B = B_strike * B_dip;
   return TPV104Params::V_w_out + (TPV104Params::V_w_in
                                   - TPV104Params::V_w_out) * B;
}

/// Nucleation spatial factor F(r) — SCEC Eq. (13).
/// F = exp(r²/(r²−R²)) for r < R, 0 for r ≥ R.
inline real_t NucleationSpatial_TPV104(real_t r)
{
   real_t R = TPV104Params::nuc_radius;
   if (r >= R) { return 0.0; }
   real_t r2 = r * r;
   real_t R2 = R * R;
   return std::exp(r2 / (r2 - R2));
}

/// Nucleation temporal factor G(t) — SCEC Eq. (14).
/// G = exp((t−T)²/(t(t−2T))) for 0 < t < T, 1 for t ≥ T, 0 for t ≤ 0.
inline real_t NucleationTemporal_TPV104(real_t t)
{
   real_t T = TPV104Params::nuc_T;
   if (t <= 0.0) { return 0.0; }
   if (t >= T)   { return 1.0; }
   return std::exp((t - T) * (t - T) / (t * (t - 2.0 * T)));
}

/// Full nucleation perturbation δτ(x, z, t) = Δτ₀ · F(r) · G(t).
inline real_t NucleationPerturbation_TPV104(real_t along_strike,
                                            real_t down_dip, real_t t)
{
   real_t dx = along_strike - TPV104Params::hypo_along_strike;
   real_t dz = down_dip     - TPV104Params::hypo_down_dip;
   real_t r  = std::sqrt(dx * dx + dz * dz);
   return TPV104Params::nuc_dtau * NucleationSpatial_TPV104(r)
          * NucleationTemporal_TPV104(t);
}

/// Compute initial ψ from the steady-state condition
///   tau_ini = sigma_n * a * asinh((V_ini/(2V0)) * exp(ψ/a))
/// Solved:
///   sinh(tau_ini/(sigma_n*a)) = (V_ini/(2V0)) * exp(ψ/a)
///   ψ = a * ln((2V0/V_ini) * sinh(tau_ini/(sigma_n*a)))
///
/// §4.1 anchor: ψ_ini(a_in=0.01) ≈ 5.6359184e-01 — matches the reference
/// TPV104 benchmark trace column 9 row 1 to 1e-8 (T_TPV104_P_1 unit test).
///
/// R-009 (review 2026-04-24): evaluated via the numerically stable
///   log(x · sinh(c)) = |c| + log((x/2) · -sign(c) · expm1(-2|c|))
/// identity, inlined here (instead of including slip_law_srw_psi.hpp
/// to avoid a layering cycle: params → friction).  Stable for large `c`
/// — e.g. `a = 0.001` would give `c ≈ 333`, where `std::sinh(c)`
/// overflows double; the logsinh form remains finite.
inline real_t ComputeInitialPsiTPV104(real_t a)
{
   const real_t arg = TPV104Params::tau_ini / (TPV104Params::sigma_n * a);
   const real_t x   = 2.0 * TPV104Params::V0 / TPV104Params::V_ini;
   const real_t sign_c = (arg >= 0.0) ? 1.0 : -1.0;
   const real_t absC   = std::abs(arg);
   // log(x · sinh(c)) via expm1 — stable for large |c|.
   return a * (absC + std::log(x / 2.0 * -sign_c
                               * std::expm1(-2.0 * absC)));
}

// =========================================================================
// Nine canonical output stations.
// Source: SCEC TPV104 benchmark-trace filename suffixes (x2, x3 in km).
// Labels match the SCEC benchmark trace naming convention and the
// probe-diff tooling (Step 13).
// =========================================================================

/// TPV104 station definition.
struct StationTPV104
{
   real_t x2;           ///< Along-strike coordinate [m]
   real_t x3;           ///< Down-dip coordinate [m]  (depth, positive)
   const char *label;   ///< SCEC benchmark filename suffix (e.g. "x2_0_x3_7.5")
};

/// Nine SCEC TPV104 fault stations (anchored to benchmark trace naming).
inline constexpr StationTPV104 kStationsTPV104[9] = {
   {    0.0,  3.0e3,  "x2_0_x3_3"     },
   {    0.0,  7.5e3,  "x2_0_x3_7.5"   },
   {    0.0, 12.0e3,  "x2_0_x3_12"    },
   {  9.0e3,  7.5e3,  "x2_9_x3_7.5"   },
   { 12.0e3,  3.0e3,  "x2_12_x3_3"    },
   { 12.0e3, 12.0e3,  "x2_12_x3_12"   },
   { -9.0e3,  7.5e3,  "x2_-9_x3_7.5"  },
   {-12.0e3,  3.0e3,  "x2_-12_x3_3"   },
   {-12.0e3, 12.0e3,  "x2_-12_x3_12"  },
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV104_PARAMS_HPP
