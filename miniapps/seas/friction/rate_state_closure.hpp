// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// rate_state_closure.hpp — Nitsche dynamic rupture, plan item N-04.
//
// Rate-and-state IFaultClosure: a pure adapter over the existing
// DieterichRuinaFriction (psi form) and the psi-form state laws.  It adds no
// friction arithmetic of its own, so its doubles equal the existing calls'
// by construction:
//
//   SolveSlipRate   -> DieterichRuinaFriction::SolveSlipRateVectorPsi with
//                      gamma in the eta (radiation-damping) slot
//                      (anti-parallel, friction/dieterich_ruina.hpp).
//   StateRate       -> AgingLawPsi::Rate or SlipLawPsi::Rate (Dc from s).
//   SeedEquilibrium -> DieterichRuinaFriction::InitialStatePsi with eta = 0.
//
// The only arithmetic here is the returned |V|, the Euclidean norm of the
// returned slip-rate vector.
//
// Scope (Chunhui Zhao, 2026-10-14).

#ifndef MFEM_SEAS_RATE_STATE_CLOSURE_HPP
#define MFEM_SEAS_RATE_STATE_CLOSURE_HPP

#include "mfem.hpp"
#include "fault_closure.hpp"
#include "dieterich_ruina.hpp"
#include "state_evolution.hpp"

#include <cmath>
#include <memory>

namespace mfem
{
namespace seas
{

class RateStateClosure : public IFaultClosure
{
public:
   enum class StateLaw { Aging, Slip };

   /// `friction` must outlive this closure.  The state law is built from the
   /// same (b, V0, f0) constants.
   RateStateClosure(const DieterichRuinaFriction &friction, StateLaw law)
      : fr_(friction), law_(law)
   {
      if (law == StateLaw::Aging)
      {
         aging_ = std::make_unique<AgingLawPsi>(fr_.GetB(), fr_.GetV0(), fr_.GetF0());
      }
      else
      {
         slip_ = std::make_unique<SlipLawPsi>(fr_.GetB(), fr_.GetV0(), fr_.GetF0());
      }
   }

   real_t SolveSlipRate(const real_t X_tan[2], real_t gamma, real_t sigma_n,
                        const FaultPointState &s, real_t V_tan[2]) const override
   {
      fr_.SolveSlipRateVectorPsi(X_tan, s.psi, sigma_n, /*eta=*/gamma, s.a, V_tan);
      return std::sqrt(V_tan[0] * V_tan[0] + V_tan[1] * V_tan[1]);
   }

   real_t StateRate(real_t V, const FaultPointState &s) const override
   {
      return (law_ == StateLaw::Aging) ? aging_->Rate(V, s.psi, s.Dc)
                                       : slip_->Rate(V, s.psi, s.Dc);
   }

   void SeedEquilibrium(real_t tau0, real_t V_ini, real_t sigma_n,
                        FaultPointState &s) const override
   {
      s.psi = fr_.InitialStatePsi(tau0, V_ini, sigma_n, /*eta=*/0.0, s.a);
   }

   StateLaw Law() const { return law_; }
   /// The friction law (its constants also build the state law); the device
   /// closure of friction/device_closures.hpp reads them.
   const DieterichRuinaFriction &Friction() const { return fr_; }

private:
   const DieterichRuinaFriction &fr_;
   StateLaw law_;
   std::unique_ptr<AgingLawPsi> aging_;
   std::unique_ptr<SlipLawPsi> slip_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RATE_STATE_CLOSURE_HPP
