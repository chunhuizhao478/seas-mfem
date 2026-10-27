// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// slip_weakening_closure.hpp — Nitsche dynamic rupture, linear slip
// weakening (benchmark plan rev 1.3, section 4.3; SCEC TPV5 / TPV205).
//
// Linear slip weakening as the projection identity of IFaultClosure:
//
//    mu(D)  = mu_s - (mu_s - mu_d) min(D / Dc, 1),   D = |d^n| = st.slip_abs
//    tau_s  = c + mu(D) sigma_n
//    |X| <= tau_s:  V = 0                                   (stick)
//    |X| >  tau_s:  V = -(|X| - tau_s) / gamma  X / |X|     (slip)
//
// so |X - gamma V| = tau_s whenever the point slips, V anti-parallel to X
// (the convention of DieterichRuinaFriction::SolveSlipRateVectorPsi).  The
// strength uses the slip at the start of the step, as our ADER code does
// for TPV205 (dynamic/tpv205_substep_iterator.cpp: |slip| before the
// solve, then LSWFrictionCoefficient_TPV205 + SolveLSW_TPV205 with the
// radiation impedance eta_s in place of gamma); the step is therefore
// first order in dt (plan, Appendix A.1).  The arithmetic is written in
// the same order as those two functions, so the two agree to round-off
// with the sign mapped (test L-04).
//
// The state psi is the slip path int |V| dt (output only): StateRate = |V|,
// and the operator's Heun update gives the trapezoid rule.  A large mu_s
// (e.g. 10000) never yields, so no barrier sentinel is needed; the file
// carries no benchmark constant.
//
// sigma_n is the compression-positive normal stress the operator hands the
// law (the Prakash-Clifton value when that is on); the operator's range
// check keeps it >= sigma_min >= 0.
//
// Scope (Chunhui Zhao, 2026-10-16).

#ifndef MFEM_SEAS_SLIP_WEAKENING_CLOSURE_HPP
#define MFEM_SEAS_SLIP_WEAKENING_CLOSURE_HPP

#include "mfem.hpp"
#include "fault_closure.hpp"

#include <cmath>

namespace mfem
{
namespace seas
{

class SlipWeakeningClosure : public IFaultClosure
{
public:
   /// mu(D): mu_s at D <= 0, mu_d at D >= Dc, linear in between (the
   /// branches of LSWFrictionCoefficient_TPV205, without its barrier
   /// sentinel).
   static real_t FrictionCoefficient(real_t D, real_t mu_s, real_t mu_d, real_t Dc)
   {
      if (D <= 0.0) { return mu_s; }
      if (D >= Dc) { return mu_d; }
      return mu_s - (mu_s - mu_d) * (D / Dc);
   }

   real_t SolveSlipRate(const real_t X_tan[2], real_t gamma, real_t sigma_n,
                        const FaultPointState &s, real_t V_tan[2]) const override
   {
      const real_t mu = FrictionCoefficient(s.slip_abs, s.mu_s, s.mu_d, s.Dc);
      const real_t tau_s = mu * sigma_n + s.cohesion;
      const real_t X_abs = std::sqrt(X_tan[0] * X_tan[0] + X_tan[1] * X_tan[1]);
      real_t V_abs = 0.0;
      if (X_abs > tau_s) { V_abs = (X_abs - tau_s) / gamma; }
      if (V_abs > 0.0 && X_abs > 0.0)
      {
         const real_t inv_X = 1.0 / X_abs;
         V_tan[0] = -(V_abs * X_tan[0] * inv_X);
         V_tan[1] = -(V_abs * X_tan[1] * inv_X);
      }
      else
      {
         V_tan[0] = 0.0;
         V_tan[1] = 0.0;
      }
      return V_abs;
   }

   /// d(slip path)/dt = |V|.
   real_t StateRate(real_t V, const FaultPointState &) const override { return V; }

   /// No slip at t = 0: the slip path starts at zero.
   void SeedEquilibrium(real_t, real_t, real_t, FaultPointState &s) const override
   {
      s.psi = 0.0;
   }

   /// mu_s >= mu_d >= 0, Dc > 0, cohesion >= 0.
   bool CheckParameters(const FaultPointState &s) const override
   {
      return s.mu_s >= s.mu_d && s.mu_d >= 0.0 && s.Dc > 0.0 && s.cohesion >= 0.0;
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_SLIP_WEAKENING_CLOSURE_HPP
