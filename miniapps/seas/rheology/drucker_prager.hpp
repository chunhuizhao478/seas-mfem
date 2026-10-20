// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// drucker_prager.hpp — non-associative Drucker-Prager viscoplasticity with
// yielding in shear, Duvaut-Lions form (SCEC TPV27; description v13, Part 6;
// rheology plan rev 1, section 6.2 and Appendix B).
//
// With the trial stress sigma_tr = sigma_0' + C_0 : (eps^{n+1} - eps_p^n)
// (+ the correction of earlier models in the stack):
//
//    sm  = tr(sigma_tr)/3,  s = sigma_tr - sm I,  tau = sqrt(1/2 s:s)
//    Y   = max(0, c cos(phi) - sm sin(phi)),        phi = atan(nu)
//    tau <= Y:  elastic, state unchanged
//    tau >  Y:  r = exp(-dt/Tv) + (1 - exp(-dt/Tv)) Y / tau
//               d eps_p = (1 - r) s / (2 mu),  eps_p += d eps_p
//               eta += sqrt(1/2 d eps_p : d eps_p)
//    Delta sigma = -C_0 : eps_p
//
// so sigma^{n+1} = sm I + r s, the sample code's return map (p. 16-17).
// sm is the EFFECTIVE mean stress because sigma_0' already carries the
// fluid pressure (the sample code's sm + fluidpresh).  The update is the
// exact solution of the relaxation ODE at fixed strain, so it adds no step
// limit (Appendix B); Tv = 0 gives the rate-independent projection r = Y/tau.
// The state is the plastic strain, not the stress, so nothing drifts
// between steps (sigma is rebuilt from the total strain every step).
//
// eta is the accumulated plastic strain (output only); its definition is to
// be checked against SeisSol's before it is compared with other codes.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_DRUCKER_PRAGER_HPP
#define MFEM_SEAS_DRUCKER_PRAGER_HPP

#include "rheology_model.hpp"

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace seas
{

struct DruckerPragerParams
{
   real_t cohesion = 0.0;          ///< c (Pa); not the frictional cohesion
   real_t bulk_friction = 0.0;     ///< nu = tan(phi)
   real_t relaxation_time = 0.0;   ///< Tv (s); 0 = rate-independent
};

class DruckerPrager : public IRheologyModel
{
public:
   /// State layout per node: eps_p[6] (tensor shear), eta.
   static constexpr int kNumState = 7;
   static constexpr int kEta = 6;

   explicit DruckerPrager(const DruckerPragerParams &p)
      : p_(p), cos_phi_(std::cos(std::atan(p.bulk_friction))),
        sin_phi_(std::sin(std::atan(p.bulk_friction)))
   {
      MFEM_VERIFY(p.cohesion >= 0.0, "drucker_prager: cohesion_pa must be >= 0, got "
                  << p.cohesion);
      MFEM_VERIFY(p.bulk_friction >= 0.0, "drucker_prager: bulk_friction must be >= 0, got "
                  << p.bulk_friction);
      MFEM_VERIFY(p.relaxation_time >= 0.0,
                  "drucker_prager: relaxation_time_s must be >= 0, got " << p.relaxation_time);
   }

   const DruckerPragerParams &Params() const { return p_; }

   /// Y(sm) = max(0, c cos(phi) - sm sin(phi)), sm the effective mean stress.
   real_t YieldStress(real_t sm) const
   {
      return std::max(real_t(0.0), p_.cohesion * cos_phi_ - sm * sin_phi_);
   }

   /// Mean stress, deviator and tau = sqrt(J2) of a stress.
   static void Invariants(const real_t sig[6], real_t &sm, real_t s[6], real_t &tau)
   {
      sm = sym6::Trace(sig) / 3.0;
      s[0] = sig[0] - sm;
      s[1] = sig[1] - sm;
      s[2] = sig[2] - sm;
      s[3] = sig[3];
      s[4] = sig[4];
      s[5] = sig[5];
      tau = std::sqrt(0.5 * sym6::DoubleDot(s, s));
   }

   const char *Name() const override { return "drucker_prager"; }
   int NumStateVars() const override { return kNumState; }
   bool IsLinearInStrain() const override { return false; }

   void InitState(const RheologyPointContext &, real_t *q) const override
   {
      for (int i = 0; i < kNumState; i++) { q[i] = 0.0; }
   }

   bool Update(const RheologyPointContext &ctx, const real_t eps[6], real_t dt,
               real_t *q, real_t dsig[6]) const override
   {
      const real_t lam = ctx.lambda0, mu = ctx.mu0;
      real_t c_ep[6], c_eps[6], sig[6];
      sym6::ApplyIsotropic(lam, mu, q, c_ep);
      sym6::ApplyIsotropic(lam, mu, eps, c_eps);
      for (int i = 0; i < 6; i++)
      {
         sig[i] = ctx.sigma0_eff[i] + c_eps[i] - c_ep[i] + dsig[i];
      }
      real_t sm, s[6], tau;
      Invariants(sig, sm, s, tau);
      const real_t Y = YieldStress(sm);
      if (tau <= Y)
      {
         for (int i = 0; i < 6; i++) { dsig[i] -= c_ep[i]; }
         return false;
      }
      // tau > Y >= 0, so tau > 0
      const real_t r = PlasticFactor(dt, Y, tau);
      real_t dep[6];
      for (int i = 0; i < 6; i++)
      {
         dep[i] = (1.0 - r) * s[i] / (2.0 * mu);
         q[i] += dep[i];
      }
      q[kEta] += std::sqrt(0.5 * sym6::DoubleDot(dep, dep));
      sym6::ApplyIsotropic(lam, mu, q, c_ep);
      for (int i = 0; i < 6; i++) { dsig[i] -= c_ep[i]; }
      return true;
   }

   /// Plastic flow only softens.
   real_t StiffnessBound(const RheologyPointContext &, const real_t *) const override
   {
      return 1.0;
   }

   int NumDiagnostics() const override { return kNumState; }
   const char *DiagnosticName(int i) const override
   {
      static const char *names[kNumState] = {"eps_p_xx", "eps_p_yy", "eps_p_zz", "eps_p_xy",
                                             "eps_p_yz", "eps_p_xz", "eta"
                                            };
      return names[i];
   }
   void Diagnostics(const RheologyPointContext &, const real_t *q, real_t *out) const override
   {
      for (int i = 0; i < kNumState; i++) { out[i] = q[i]; }
   }

private:
   /// r = exp(-dt/Tv) + (1 - exp(-dt/Tv)) Y/tau, with -expm1 for 1 - exp.
   real_t PlasticFactor(real_t dt, real_t Y, real_t tau) const
   {
      if (p_.relaxation_time == 0.0) { return Y / tau; }
      const real_t a = -dt / p_.relaxation_time;
      const real_t one_minus_decay = -std::expm1(a);
      return std::exp(a) + one_minus_decay * (Y / tau);
   }

   DruckerPragerParams p_;
   real_t cos_phi_, sin_phi_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DRUCKER_PRAGER_HPP
