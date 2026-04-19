// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV102 benchmark parameters.
// Reference: https://strike.scec.org/cvws/tpv101_102docs.html

#ifndef MFEM_SEAS_TPV102_PARAMS_HPP
#define MFEM_SEAS_TPV102_PARAMS_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{
namespace seas
{

/// SCEC TPV102 material and friction parameters.
struct TPV102Params
{
   // Material (homogeneous, isotropic)
   static constexpr real_t rho = 2670.0;       ///< Density [kg/m³]
   static constexpr real_t cs = 3464.0;         ///< S-wave speed [m/s]
   static constexpr real_t cp = 6000.0;         ///< P-wave speed [m/s]
   static constexpr real_t mu = rho * cs * cs;  ///< Shear modulus [Pa]
   static constexpr real_t lambda = rho * cp * cp - 2.0 * mu;  ///< First Lame [Pa]

   // Impedances
   static constexpr real_t Zp = rho * cp;
   static constexpr real_t Zs = rho * cs;
   static constexpr real_t eta_s = Zs / 2.0;   ///< S-impedance harmonic mean (homogeneous)

   // Friction law
   static constexpr real_t f0 = 0.6;            ///< Reference friction coefficient
   static constexpr real_t V0 = 1e-6;           ///< Reference slip rate [m/s]
   static constexpr real_t b = 0.012;           ///< State evolution parameter
   static constexpr real_t Dc = 0.02;           ///< Critical slip distance [m] (= L in spec)

   // Spatial distribution of direct effect parameter a
   static constexpr real_t a_vw = 0.008;        ///< VW zone (velocity-weakening)
   static constexpr real_t a_vs = 0.016;        ///< VS zone (velocity-strengthening)
   static constexpr real_t W = 15e3;            ///< VW half-width [m] along-dip
   static constexpr real_t w = 3e3;             ///< VW-VS transition width [m]

   // Along-strike VW zone (fault is 30 km long, VW in central part)
   static constexpr real_t Ls = 15e3;           ///< Along-strike VW half-width [m]
   static constexpr real_t ws = 3e3;            ///< Along-strike transition width [m]

   // Fault geometry (VW zone + 3 km transition on each side)
   static constexpr real_t fault_length = 36e3; ///< Along-strike [m] (VW 30 + 3 km each side)
   static constexpr real_t fault_depth = 18e3;  ///< Down-dip [m] (VW 15 + 3 km at bottom)

   // Loading and initial conditions
   static constexpr real_t sigma_n = 120e6;     ///< Normal stress [Pa] (compression)
   static constexpr real_t tau_ini = 75e6;       ///< Initial shear traction [Pa]
   static constexpr real_t V_ini = 1e-12;       ///< Initial slip rate [m/s]

   // Nucleation
   static constexpr real_t hypo_along_strike = 0.0;  ///< Hypocenter along-strike [m]
   static constexpr real_t hypo_down_dip = 7.5e3;    ///< Hypocenter down-dip [m]
   static constexpr real_t nuc_radius = 3e3;          ///< Nucleation radius [m]
   static constexpr real_t nuc_dtau = 25e6;           ///< Nucleation stress perturbation [Pa]
   static constexpr real_t nuc_T = 1.0;               ///< Nucleation rise time [s]

   // Domain
   static constexpr real_t domain_half = 60e3;  ///< Domain half-size [m] (each direction)

   // Simulation
   static constexpr real_t t_final = 12.0;      ///< Final time [s]
};

/// SCEC Boxcar function B(x, W, w) — Eq. (5) of SCEC TPV101/102 spec.
/// B = 1 for |x| <= W, tanh transition for W < |x| < W+w, 0 for |x| >= W+w.
/// Uses C-infinity tanh taper (required for SCEC cross-code comparison).
inline real_t Boxcar(real_t x, real_t W, real_t w_trans)
{
   real_t ax = std::abs(x);
   if (ax <= W) { return 1.0; }
   if (ax >= W + w_trans) { return 0.0; }
   // SCEC Eq. (5): tanh transition
   return 0.5 * (1.0 + std::tanh(w_trans / (ax - W - w_trans)
                                + w_trans / (ax - W)));
}

/// Compute direct effect parameter a at position (along_strike, down_dip).
/// SCEC Eq. (4): delta_a = delta_a0 * [1 - B(x; W, w) * B(y - y0; W/2, w)]
/// VW covers the entire fault extent, transition starts OUTSIDE the fault.
inline real_t ComputeA(real_t along_strike, real_t down_dip)
{
   real_t B_strike = Boxcar(along_strike, TPV102Params::Ls, TPV102Params::ws);
   // Dip: centered at hypocenter depth y0=7.5km, half-width W/2=7.5km
   real_t B_dip = Boxcar(down_dip - TPV102Params::hypo_down_dip,
                         TPV102Params::W / 2.0, TPV102Params::w);
   real_t B = B_strike * B_dip;
   return TPV102Params::a_vs + (TPV102Params::a_vw - TPV102Params::a_vs) * B;
}

/// Nucleation perturbation spatial factor F(r).
/// F = exp(r²/(r²-R²)) for r < R, 0 for r >= R.
inline real_t NucleationSpatial(real_t r)
{
   real_t R = TPV102Params::nuc_radius;
   if (r >= R) { return 0.0; }
   real_t r2 = r * r;
   real_t R2 = R * R;
   return std::exp(r2 / (r2 - R2));
}

/// Nucleation perturbation temporal factor G(t).
/// G = exp((t-T)²/(t(t-2T))) for 0 < t < T, 1 for t >= T.
inline real_t NucleationTemporal(real_t t)
{
   real_t T = TPV102Params::nuc_T;
   if (t <= 0.0) { return 0.0; }
   if (t >= T) { return 1.0; }
   return std::exp((t - T) * (t - T) / (t * (t - 2.0 * T)));
}

/// Full nucleation perturbation delta_tau(x, z, t).
inline real_t NucleationPerturbation(real_t along_strike, real_t down_dip, real_t t)
{
   real_t dx = along_strike - TPV102Params::hypo_along_strike;
   real_t dz = down_dip - TPV102Params::hypo_down_dip;
   real_t r = std::sqrt(dx*dx + dz*dz);
   return TPV102Params::nuc_dtau * NucleationSpatial(r) * NucleationTemporal(t);
}

/// Compute initial state variable theta_ini from equilibrium.
/// From: tau_ini = sigma_n * f(V_ini, theta_ini)
/// Invert: theta_ini = (Dc/V0) * exp((tau_ini/(sigma_n*a) - f0 - a*ln(V_ini/(2*V0))) / b)
/// Using psi = f0 + b * ln(V0 * theta / Dc):
///   psi_ini such that tau_ini = sigma_n * a * asinh(V_ini/(2*V0) * exp(psi_ini/a))
inline real_t ComputeInitialPsi(real_t a)
{
   // From tau = sigma_n * a * asinh(V/(2*V0) * exp(psi/a)):
   //   sinh(tau/(sigma_n*a)) = V/(2*V0) * exp(psi/a)
   //   psi = a * ln(2*V0/V * sinh(tau/(sigma_n*a)))
   real_t arg = TPV102Params::tau_ini / (TPV102Params::sigma_n * a);
   real_t psi = a * std::log(2.0 * TPV102Params::V0 / TPV102Params::V_ini * std::sinh(arg));
   return psi;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV102_PARAMS_HPP
