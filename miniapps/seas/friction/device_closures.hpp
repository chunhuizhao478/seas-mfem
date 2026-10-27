// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// device_closures.hpp — GPU dynamic-rupture plan D3 (plan rev2 §4.4): the
// Nitsche fault closures as MFEM_HOST_DEVICE functions, for the device
// fault step (NitscheWaveOperator with SetDeviceFault).
//
// Each function repeats the host code step for step, as
// device_rate_state.hpp does for the quasi-dynamic friction solve (GPU plan
// G4, gpu-qd-bp5), and the rate-and-state solve IS that file's:
//
//   DCSolveSlipRate  IFaultClosure::SolveSlipRate of
//                      SlipWeakeningClosure   (FrictionCoefficient and its
//                                              stick / slip projection)
//                      RateStateClosure,      DieterichRuinaFriction::
//                      RateStateSRWClosure    SolveSlipRateVectorPsi with
//                                             eta = gamma, i.e.
//                                             RSSolveSlipRateVectorPsi
//   DCStateRate      IFaultClosure::StateRate of the same laws:
//                      slip weakening         |V|
//                      aging (psi form)       AgingLawPsi::Rate = RSAgingRatePsi
//                      slip law (psi form)    SlipLawPsi::Rate
//                      strong rate weakening  SlipLawSRWPsi::Rate_SRW
//                                             (PsiSS_SRW, LogSinhStable,
//                                             IntegerPow8)
//
// so the debug device gives the host's bits; on a GPU exp, log, pow, asinh
// and expm1 may differ in the last bit (test_device_closures bounds V).  A
// kernel cannot print or abort: the rate-and-state solve sets the RSStatus
// bits of device_rate_state.hpp where the host prints or stops, and the
// operator reruns its host loop on such a step.
//
// The per-point parameters of a fault DOF are the FaultPointState fields
// the laws read, kDCParams values per DOF in the order of DCParam.

#ifndef MFEM_SEAS_DEVICE_CLOSURES_HPP
#define MFEM_SEAS_DEVICE_CLOSURES_HPP

#include "mfem.hpp"
#include "device_rate_state.hpp"
#include "fault_closure.hpp"
#include "rate_state_closure.hpp"
#include "rate_state_srw_closure.hpp"
#include "slip_weakening_closure.hpp"

#include <cmath>

namespace mfem
{
namespace seas
{

/// The closures with a device form.
enum class DeviceClosureLaw : int
{
   None = 0,
   SlipWeakening,   ///< SlipWeakeningClosure
   Aging,           ///< RateStateClosure, StateLaw::Aging
   Slip,            ///< RateStateClosure, StateLaw::Slip
   SRW              ///< RateStateSRWClosure
};

/// Per-point parameters, kDCParams per fault DOF in this order.
enum DCParam : int { kDCa = 0, kDCDc, kDCmu_s, kDCmu_d, kDCcohesion, kDCV_w, kDCParams };

/// Law constants of one closure (the per-point values are in DCParam).
struct DeviceClosure
{
   DeviceClosureLaw law = DeviceClosureLaw::None;
   RSDeviceParams rs{};   ///< friction V0 and the aging law's b, V0, f0
   real_t b = 0.0;        ///< slip law / SRW state law constants
   real_t V0 = 0.0;
   real_t f0 = 0.0;
   real_t muW = 0.0;      ///< SRW weakening friction coefficient
};

/// The per-point parameters of `s` in DCParam order.
inline void DCPackParams(const FaultPointState &s, real_t q[kDCParams])
{
   q[kDCa] = s.a;
   q[kDCDc] = s.Dc;
   q[kDCmu_s] = s.mu_s;
   q[kDCmu_d] = s.mu_d;
   q[kDCcohesion] = s.cohesion;
   q[kDCV_w] = s.V_w;
}

/// The device form of `cl`; false (and dc.law = None) when it has none.
inline bool MakeDeviceClosure(const IFaultClosure &cl, DeviceClosure &dc)
{
   dc = DeviceClosure{};
   if (dynamic_cast<const SlipWeakeningClosure *>(&cl))
   {
      dc.law = DeviceClosureLaw::SlipWeakening;
      return true;
   }
   if (auto *rs = dynamic_cast<const RateStateClosure *>(&cl))
   {
      // the state law is built from the friction law's b, V0, f0
      const DieterichRuinaFriction &fr = rs->Friction();
      dc.law = (rs->Law() == RateStateClosure::StateLaw::Aging) ? DeviceClosureLaw::Aging
                                                                : DeviceClosureLaw::Slip;
      dc.rs = RSDeviceParams{fr.GetV0(), fr.GetB(), fr.GetV0(), fr.GetF0()};
      dc.b = fr.GetB();
      dc.V0 = fr.GetV0();
      dc.f0 = fr.GetF0();
      return true;
   }
   if (auto *srw = dynamic_cast<const RateStateSRWClosure *>(&cl))
   {
      const SlipLawSRWPsi &law = srw->StateLaw();
      dc.law = DeviceClosureLaw::SRW;
      dc.rs.V0 = srw->Friction().GetV0();
      dc.b = law.GetB();
      dc.V0 = law.GetV0();
      dc.f0 = law.GetF0();
      dc.muW = law.GetMuW();
      return true;
   }
   return false;
}

/// SlipWeakeningClosure::FrictionCoefficient.
MFEM_HOST_DEVICE inline real_t DCSlipWeakeningMu(real_t D, real_t mu_s, real_t mu_d,
                                                 real_t Dc)
{
   if (D <= 0.0) { return mu_s; }
   if (D >= Dc) { return mu_d; }
   return mu_s - (mu_s - mu_d) * (D / Dc);
}

/// LogSinhStable of slip_law_srw_psi.hpp.
MFEM_HOST_DEVICE inline real_t DCLogSinhStable(real_t x, real_t c)
{
   const real_t sign_c = (c >= 0.0) ? 1.0 : -1.0;
   const real_t absC = std::fabs(c);
   return absC + std::log(x / 2.0 * -sign_c * std::expm1(-2.0 * absC));
}

/// SlipLawSRWPsi::PsiSS_SRW (IntegerPow8 written out).
MFEM_HOST_DEVICE inline real_t DCPsiSSSRW(real_t V, real_t V_w, real_t a, real_t b,
                                          real_t V0, real_t f0, real_t muW)
{
   const real_t f_LV0 = f0 - (b - a) * std::log(V / V0);
   const real_t f_LV = (static_cast<real_t>(0) < f_LV0) ? f_LV0 : static_cast<real_t>(0);
   const real_t V_over_Vw = V / V_w;
   const real_t x2 = V_over_Vw * V_over_Vw;
   const real_t x4 = x2 * x2;
   const real_t V_over_Vw_8 = x4 * x4;
   const real_t denom = std::pow(1.0 + V_over_Vw_8, static_cast<real_t>(1.0 / 8.0));
   const real_t f_ss = muW + (f_LV - muW) / denom;
   return a * DCLogSinhStable(2.0 * V0 / V, f_ss / a);
}

/// IFaultClosure::StateRate(V, st) with st.psi = psi and the per-point
/// parameters q (DCParam order).
MFEM_HOST_DEVICE inline real_t DCStateRate(const DeviceClosure &c, const real_t *q, real_t V,
                                           real_t psi)
{
   switch (c.law)
   {
      case DeviceClosureLaw::SlipWeakening: return V;
      case DeviceClosureLaw::Aging: return RSAgingRatePsi(V, psi, q[kDCDc], c.rs);
      case DeviceClosureLaw::Slip:
      {
         // SlipLawPsi::Rate
         V = (V < 1e-50) ? 1e-50 : V;   // std::max(V, 1e-50)
         const real_t psi_ss = c.f0 + c.b * std::log(c.V0 / V);
         return -(V / q[kDCDc]) * (psi - psi_ss);
      }
      case DeviceClosureLaw::SRW:
      {
         // SlipLawSRWPsi::Rate_SRW(V, psi, L = Dc, V_w, a)
         const real_t psi_ss = DCPsiSSSRW(V, q[kDCV_w], q[kDCa], c.b, c.V0, c.f0, c.muW);
         return -(V / q[kDCDc]) * (psi - psi_ss);
      }
      default: return kRSNaN;
   }
}

/// IFaultClosure::SolveSlipRate(X_tan, gamma, sigma_n, st, V_tan) with
/// st.psi = psi, st.slip_abs = slip_abs and the per-point parameters q.
/// The rate-and-state laws OR their RSStatus bits into `status`.
MFEM_HOST_DEVICE inline void DCSolveSlipRate(const DeviceClosure &c, const real_t *q,
                                             const real_t X_tan[2], real_t gamma,
                                             real_t sigma_n, real_t psi, real_t slip_abs,
                                             real_t V_tan[2], int &status)
{
   if (c.law == DeviceClosureLaw::SlipWeakening)
   {
      const real_t mu = DCSlipWeakeningMu(slip_abs, q[kDCmu_s], q[kDCmu_d], q[kDCDc]);
      const real_t tau_s = mu * sigma_n + q[kDCcohesion];
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
      return;
   }
   RSSolveSlipRateVectorPsi(X_tan[0], X_tan[1], psi, sigma_n, gamma, q[kDCa], c.rs.V0,
                            V_tan[0], V_tan[1], status);
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DEVICE_CLOSURES_HPP
