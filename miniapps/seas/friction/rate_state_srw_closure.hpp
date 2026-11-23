// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// rate_state_srw_closure.hpp — Nitsche dynamic rupture, rate and state with
// the slip law and strong rate weakening (benchmark plan rev 1.3, section
// 4.4; SCEC TPV104).
//
// A pure adapter, like RateStateClosure (rate_state_closure.hpp): it adds
// no friction arithmetic of its own, so its doubles equal the existing
// calls' by construction (test S-01):
//
//   SolveSlipRate   -> DieterichRuinaFriction::SolveSlipRateVectorPsi with
//                      gamma in the eta slot (the friction law is the same
//                      regularized form as the aging closure's).
//   StateRate       -> SlipLawSRWPsi::Rate_SRW(V, psi, L = s.Dc, s.V_w, s.a)
//                      with b, V0, f0 from the friction constants and f_w
//                      given to the constructor.
//   SeedEquilibrium -> DieterichRuinaFriction::InitialStatePsi with eta = 0
//                      (TPV104 eq. 11).
//   CheckParameters -> a > 0, Dc > 0, V_w > 0.
//
// The SlipLawSRWPsi instance is put in production mode, so its scalar-
// default virtuals abort and only the per-point Rate_SRW runs (the R-001
// guard in slip_law_srw_psi.hpp).  dieterich_ruina.hpp, state_evolution.hpp
// and slip_law_srw_psi.hpp are not edited.
//
// The returned |V| is the Euclidean norm of the returned slip-rate vector.
//
// Scope (Chunhui Zhao, 2026-10-16).

#ifndef MFEM_SEAS_RATE_STATE_SRW_CLOSURE_HPP
#define MFEM_SEAS_RATE_STATE_SRW_CLOSURE_HPP

#include "mfem.hpp"
#include "fault_closure.hpp"
#include "dieterich_ruina.hpp"
#include "slip_law_srw_psi.hpp"

#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

class RateStateSRWClosure : public IFaultClosure
{
public:
   /// `friction` must outlive this closure.  f_w is the weakening friction
   /// coefficient (TPV104: 0.2).  The SRW law's scalar a and V_w defaults
   /// are NaN: production mode never reads them.
   RateStateSRWClosure(const DieterichRuinaFriction &friction, real_t f_w)
      : fr_(friction),
        srw_(std::numeric_limits<real_t>::quiet_NaN(), friction.GetB(), friction.GetV0(),
             friction.GetF0(), f_w, std::numeric_limits<real_t>::quiet_NaN())
   {
      srw_.SetProductionMode();
   }

   real_t SolveSlipRate(const real_t X_tan[2], real_t gamma, real_t sigma_n,
                        const FaultPointState &s, real_t V_tan[2]) const override
   {
      fr_.SolveSlipRateVectorPsi(X_tan, s.psi, sigma_n, /*eta=*/gamma, s.a, V_tan);
      return std::sqrt(V_tan[0] * V_tan[0] + V_tan[1] * V_tan[1]);
   }

   real_t StateRate(real_t V, const FaultPointState &s) const override
   {
      return srw_.Rate_SRW(V, s.psi, s.Dc, s.V_w, s.a);
   }

   void SeedEquilibrium(real_t tau0, real_t V_ini, real_t sigma_n,
                        FaultPointState &s) const override
   {
      s.psi = fr_.InitialStatePsi(tau0, V_ini, sigma_n, /*eta=*/0.0, s.a);
   }

   bool CheckParameters(const FaultPointState &s) const override
   {
      return s.a > 0.0 && s.Dc > 0.0 && s.V_w > 0.0;
   }

   const SlipLawSRWPsi &StateLaw() const { return srw_; }
   /// The friction law; the device closure of friction/device_closures.hpp
   /// reads its constants.
   const DieterichRuinaFriction &Friction() const { return fr_; }
   real_t WeakeningFriction() const { return srw_.GetMuW(); }

private:
   const DieterichRuinaFriction &fr_;
   SlipLawSRWPsi srw_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RATE_STATE_SRW_CLOSURE_HPP
