// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// prestress.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.3 c):
// BuildPrestress, the single source of the fault prestress of both regimes.
//
// Steady state at the initial slip rate on every owned fault DOF (SCEC BP5
// eqs. 18, 20, 22, 23):
//   psi0  = f0 + b ln(V0 / V_bg)                 theta = Dc / V_bg (eq. 18)
//   |V|   = |V_init|                             two kernel components
//   tau   = sigma_n a asinh( |V| / (2 V0) exp(psi0 / a) )
//   tau  += eta |V|                              only if include_eta
//   tau0  = -tau V_init / |V|                    anti-parallel, kernel frame
//   sigma0 = sigma_n
// The kernel frame and its anti-parallel sign are those of the QD kernel
// (fault/rate_state_fault.hpp, config/bp5_params.hpp tau0_vec).  The
// arithmetic follows BP5Params::tau0_vec term by term.  QD receives tau0,
// V_init and sigma0 through hook H1b; the dynamic operator receives the same
// values as global vectors through hook H2 (P2).
//
// The IC set fixes include_eta: "gold" (Tandem convention, acceptance A) adds
// eta |V| (eq. 20, and eq. 23 with dtau = eta V_i); "scec_fd" does not
// (eq. 22, eq. 23 with dtau = 0).

#ifndef MFEM_SEAS_CYCLE_PRESTRESS_HPP
#define MFEM_SEAS_CYCLE_PRESTRESS_HPP

#include "mfem.hpp"
#include "fault_fields.hpp"

#include <cmath>

namespace mfem
{
namespace seas
{
namespace cycle
{

struct PrestressConfig
{
   bool include_eta = true;    ///< "gold": true; "scec_fd": false
   real_t V_background = 0.0;  ///< V_bg of psi0 (BP5: V_init = 1e-9 m/s)
};

struct Prestress
{
   Vector tau0;    ///< [2N] (dip, strike) per owned DOF, kernel frame [Pa]
   Vector sigma0;  ///< [N] normal stress, compression positive [Pa]
   Vector psi0;    ///< [N] steady-state psi (eq. 18), a check value
};

inline Prestress BuildPrestress(const FaultParams &par, const PrestressConfig &ic)
{
   MFEM_VERIFY(ic.V_background > 0.0, "BuildPrestress: V_background must be > 0");
   const spatial::RateStatePerDOFParams &rs = par.rs;
   const int N = rs.a.Size();
   MFEM_VERIFY(par.V_init.Size() == 2 * N && rs.sigma_n_eff.Size() == N &&
               rs.eta.Size() == N && rs.f_0.Size() == N && rs.b.Size() == N &&
               rs.V_0.Size() == N, "BuildPrestress: inconsistent parameter sizes");

   Prestress p;
   p.tau0.SetSize(2 * N);
   p.sigma0.SetSize(N);
   p.psi0.SetSize(N);
   for (int i = 0; i < N; i++)
   {
      const real_t Vi[2] = {par.V_init(2 * i), par.V_init(2 * i + 1)};
      const real_t Vi_abs = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);
      MFEM_VERIFY(Vi_abs > 0.0, "BuildPrestress: |V_init| = 0 at DOF " << i);
      const real_t a = rs.a(i);
      const real_t sigma_n = rs.sigma_n_eff(i);
      const real_t V0 = rs.V_0(i);

      const real_t psi0 = rs.f_0(i) + rs.b(i) * std::log(V0 / ic.V_background);
      const real_t e = std::exp(psi0 / a);
      real_t tau = sigma_n * a * std::asinh((Vi_abs / (2.0 * V0)) * e);
      if (ic.include_eta) { tau += rs.eta(i) * Vi_abs; }

      p.tau0(2 * i) = -tau * Vi[0] / Vi_abs;
      p.tau0(2 * i + 1) = -tau * Vi[1] / Vi_abs;
      p.sigma0(i) = sigma_n;
      p.psi0(i) = psi0;
   }
   return p;
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_PRESTRESS_HPP
