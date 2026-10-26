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
// arithmetic follows BP5Params::tau0_vec term by term, down to tau being one
// expression, so that a compiler that fuses a*b + c (clang on arm64 does by
// default) rounds it the same way in both.  QD receives tau0,
// V_init and sigma0 through hook H1b; the dynamic operator receives the same
// values as global vectors through hook H2 (P2).
//
// include_eta is [cycle.prestress].eta_V_init.  "gold" (Tandem convention,
// acceptance A) always adds eta |V| (eq. 20, and eq. 23 with dtau = eta V_i);
// "scec_fd" does not by default (eq. 22, eq. 23 with dtau = 0), and with
// eta_V_init = true it does, so that the dynamic regime starts at V_init
// from u = u' = 0.
//
// [cycle.prestress].a (2D in-plane design rev 3, decision D5, derivation
// 10.4) chooses the a of the asinh.  "local" (the default) is the above, the
// DOF's own a: BP5's convention (eqs 18-23, Tandem bp5.lua), steady state at
// V_init everywhere.  "amax" is SCEC BP3's (eqs 22-26, Tandem bp3.lua):
//   tau   = sigma_n a_max asinh( |V| / (2 V0) exp(psi0 / a_max) ) + eta |V|
//   psi0  = a ln[ (2 V0 / |V|) sinh( (tau - eta |V|) / (a sigma_n) ) ]
// with a_max the global maximum of the per-DOF a (GlobalMaxA, one
// allreduce) and the DOF's own a in psi0: a uniform tau0, equilibrium at
// V_init everywhere, steady state only where a = a_max.  psi0 is then the
// kernel's inversion (DieterichRuinaFriction::InitialStatePsi, called here),
// which the kernel runs with ic = "gold"; the config refuses "amax" with
// "scec_fd" (whose uniform psi0 would leave the weakening region out of
// equilibrium at t = 0) and BuildPrestress refuses it without eta |V| (the
// inversion subtracts eta |V|; BP3-FD's eq. 22 is phase F's).  The "local"
// arithmetic is unchanged to the character.

#ifndef MFEM_SEAS_CYCLE_PRESTRESS_HPP
#define MFEM_SEAS_CYCLE_PRESTRESS_HPP

#include "mfem.hpp"
#include "fault_fields.hpp"
#include "../friction/dieterich_ruina.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace mfem
{
namespace seas
{
namespace cycle
{

struct PrestressConfig
{
   bool include_eta = true;    ///< [cycle.prestress].eta_V_init ("gold": true)
   real_t V_background = 0.0;  ///< V_bg of psi0 (BP5: V_init = 1e-9 m/s)
   std::string a = "local";    ///< [cycle.prestress].a: "local" (BP5) or "amax" (BP3)
   /// a_max of a = "amax": not a deck key, the caller sets it with
   /// GlobalMaxA before BuildPrestress.
   real_t a_max = 0.0;
   bool AMax() const { return a == "amax"; }
};

struct Prestress
{
   Vector tau0;    ///< [2N] (dip, strike) per owned DOF, kernel frame [Pa]
   Vector sigma0;  ///< [N] normal stress, compression positive [Pa]
   Vector psi0;    ///< [N] psi of the IC (eq. 18; a = "amax": the inversion), a check value
};

/// a_max of [cycle.prestress].a = "amax": the maximum of the per-DOF a over
/// the owned fault DOFs of every rank (collective).
inline real_t GlobalMaxA(const FaultParams &par, MPI_Comm comm)
{
   real_t loc = 0.0, glob = 0.0;
   for (int i = 0; i < par.rs.a.Size(); i++) { loc = std::max(loc, par.rs.a(i)); }
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MFEM_VERIFY(glob > 0.0, "GlobalMaxA: no fault DOF with a > 0 on any rank");
   return glob;
}

inline Prestress BuildPrestress(const FaultParams &par, const PrestressConfig &ic)
{
   MFEM_VERIFY(ic.V_background > 0.0, "BuildPrestress: V_background must be > 0");
   MFEM_VERIFY(ic.a == "local" || ic.AMax(),
               "BuildPrestress: a must be \"local\" or \"amax\"; got \"" << ic.a << "\"");
   MFEM_VERIFY(!ic.AMax() || ic.a_max > 0.0,
               "BuildPrestress: a = \"amax\" needs a_max > 0 (GlobalMaxA)");
   MFEM_VERIFY(!ic.AMax() || ic.include_eta,
               "BuildPrestress: a = \"amax\" needs eta_V_init = true: its psi0 is the "
               "kernel's inversion, which subtracts eta |V_init| from tau0");
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
      if (ic.AMax())
      {
         // BP3 (see the file comment): the expression below with a_max, and
         // psi0 by the kernel's own inversion with the DOF's a.
         const real_t am = ic.a_max;
         const real_t em = std::exp(psi0 / am);
         const real_t tau =
            sigma_n * am * std::asinh((Vi_abs / (2.0 * V0)) * em) + rs.eta(i) * Vi_abs;
         DieterichRuinaFriction::Constants fc;
         fc.V0 = V0;
         fc.f0 = rs.f_0(i);
         fc.b = rs.b(i);
         const DieterichRuinaFriction law(fc);
         p.tau0(2 * i) = -tau * Vi[0] / Vi_abs;
         p.tau0(2 * i + 1) = -tau * Vi[1] / Vi_abs;
         p.sigma0(i) = sigma_n;
         p.psi0(i) = law.InitialStatePsi(tau, Vi_abs, sigma_n, rs.eta(i), a);
         continue;
      }
      const real_t e = std::exp(psi0 / a);
      // One expression, as in BP5Params::tau0_vec.  Adding eta |V| in a
      // second statement left tau 1 ulp off under FMA contraction (Mac,
      // clang 22), which changed the first adaptive step against
      // seas_bp5_full.
      const real_t tau =
         ic.include_eta
         ? sigma_n * a * std::asinh((Vi_abs / (2.0 * V0)) * e) + rs.eta(i) * Vi_abs
         : sigma_n * a * std::asinh((Vi_abs / (2.0 * V0)) * e);

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
