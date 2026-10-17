// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_DEVICE_RATE_STATE_HPP
#define MFEM_SEAS_DEVICE_RATE_STATE_HPP

#include "mfem.hpp"
#include "general/reducers.hpp"  // mfem::reduce, DevicePair (not in mfem.hpp)
#include <cmath>
#include <limits>

namespace mfem
{
namespace seas
{

/// @file
/// GPU plan G4: the BP5 friction solve (Dieterich-Ruina, psi form) and the
/// aging law in psi form as MFEM_HOST_DEVICE functions, for mfem::forall
/// kernels over fault nodes.
///
/// Each function repeats the host code step for step, with the same
/// constants and bracket fallbacks:
///   RSFrictionCoefficientPsi  DieterichRuinaFriction::FrictionCoefficientPsi
///   RSZeroIn                  DieterichRuinaFriction::zeroIn
///   RSSolveSlipRatePsi        DieterichRuinaFriction::SolveSlipRatePsi
///   RSSolveSlipRateVectorPsi  DieterichRuinaFriction::SolveSlipRateVectorPsi
///   RSAgingRatePsi            AgingLawPsi::Rate
/// so a CPU run gives the same bits as the host code.  On a GPU, exp, pow
/// and log10 may differ in the last bit; test_device_friction bounds the
/// difference in V.
///
/// A kernel cannot print or abort.  Where the host code prints or aborts, a
/// node sets a bit in its status word instead (RSStatus), and the caller
/// reruns the host code, which prints today's lines and aborts as before.

/// Status bits of one node's friction solve.
enum RSStatus : int
{
   kRSBadInput          = 1,  ///< non-finite psi, tau, sigma_n, a, eta or Dc
   kRSSigmaNNonPositive = 2,  ///< sigma_n_eff <= 0 (the host warns)
   kRSResidualNonFinite = 4,  ///< a residual was not finite ([FRIC-RESID])
   kRSZeroInFailed      = 8,  ///< Brent's method met a non-finite residual
   kRSNoBracket         = 16, ///< both brackets failed; V is NaN
   kRSNonFiniteV        = 32  ///< non-finite slip rate
};

/// Status bits after which the host code stops the run: kRSAbortSolve in
/// the bare solve (RecomputeSlipRate), kRSAbortRHS in ComputeRHS, which also
/// checks its inputs and the slip rate.  A failed Brent step stops the run
/// only without exceptions: with MFEM_USE_EXCEPTIONS the host catches it and
/// tries the next bracket, as the kernel does.
#ifdef MFEM_USE_EXCEPTIONS
constexpr int kRSAbortSolve = 0;
#else
constexpr int kRSAbortSolve = kRSZeroInFailed;
#endif
constexpr int kRSAbortRHS = kRSBadInput | kRSNoBracket | kRSNonFiniteV |
                            kRSAbortSolve;

/// Constants of the friction law and of the aging law (psi form).  The two
/// laws keep their own copies on the host (DieterichRuinaFriction::Constants
/// and AgingLawPsi), so they are kept apart here too.
struct RSDeviceParams
{
   real_t V0;        ///< friction law reference slip rate
   real_t aging_b;   ///< AgingLawPsi b
   real_t aging_V0;  ///< AgingLawPsi V0
   real_t aging_f0;  ///< AgingLawPsi f0
};

constexpr real_t kRSEps = std::numeric_limits<real_t>::epsilon();
constexpr real_t kRSNaN = std::numeric_limits<real_t>::quiet_NaN();

/// f(V, psi) = a * asinh[(V / 2V0) * exp(psi / a)], f = 0 at V <= 0, and the
/// large-argument form for psi / a > 700.
MFEM_HOST_DEVICE inline real_t RSFrictionCoefficientPsi(real_t V, real_t psi,
                                                        real_t a, real_t V0)
{
   if (V <= 0.0) { return 0.0; }

   const real_t psi_over_a = psi / a;
   if (psi_over_a > 700.0)
   {
      const real_t log_term = std::log(V / V0);
      const real_t f = a * log_term + psi;
      return (0.0 < f) ? f : 0.0;  // std::max(0.0, f)
   }

   const real_t sinh_arg = (V / (2.0 * V0)) * std::exp(psi_over_a);
   return a * std::asinh(sinh_arg);
}

/// Brent's method (Tandem's zeroIn) on [a, b], F(a) and F(b) of opposite
/// signs.  Where the host version stops the run (MFEM_VERIFY), this one sets
/// kRSZeroInFailed in @a status and returns NaN.
template <typename F>
MFEM_HOST_DEVICE inline real_t RSZeroIn(real_t a, real_t b, F &&Fn,
                                        int &status, real_t tol = 0.0)
{
   real_t Fa = Fn(a);
   if (Fa == 0.0) { return a; }
   real_t Fb = Fn(b);

   if (!std::isfinite(Fa) || !std::isfinite(Fb) ||
       !(std::copysign(Fa, Fb) != Fa || Fb == 0.0))
   {
      status |= kRSZeroInFailed;
      return kRSNaN;
   }

   real_t c = a, Fc = Fa;
   real_t d = b - a, e = d;

   while (Fb != 0.0)
   {
      if (std::copysign(Fb, Fc) == Fb)
      {
         c = a; Fc = Fa; d = b - a; e = d;
      }
      if (std::fabs(Fc) < std::fabs(Fb))
      {
         a = b; b = c; c = a;
         Fa = Fb; Fb = Fc; Fc = Fa;
      }
      const real_t xm = 0.5 * (c - b);
      const real_t tol1 = 2.0 * kRSEps * std::fabs(b) + 0.5 * tol;
      if (std::fabs(xm) <= tol1 || Fb == 0.0) { break; }
      if (std::fabs(e) < tol1 || std::fabs(Fa) <= std::fabs(Fb))
      {
         d = xm; e = d;  // bisection
      }
      else
      {
         const real_t s = Fb / Fa;
         real_t p, q;
         if (a != c)
         {
            // inverse quadratic interpolation
            const real_t qq = Fa / Fc;
            const real_t r = Fb / Fc;
            p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
            q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
         }
         else
         {
            // linear interpolation
            p = 2.0 * xm * s;
            q = 1.0 - s;
         }
         if (p > 0) { q = -q; } else { p = -p; }
         if (2.0 * p < 3.0 * xm * q - std::fabs(tol1 * q) &&
             p < std::fabs(0.5 * e * q))
         {
            e = d; d = p / q;
         }
         else
         {
            d = xm; e = d;  // bisection
         }
      }
      a = b; Fa = Fb;
      if (std::fabs(d) > tol1) { b += d; }
      else { b += std::copysign(tol1, xm); }
      Fb = Fn(b);
      if (!std::isfinite(Fb))
      {
         status |= kRSZeroInFailed;
         return kRSNaN;
      }
   }
   return b;
}

/// Solve tau = sigma_n * f(V, psi) + eta * V for V >= 0: Brent's method in
/// log10(V) on [-32, log10(tau/eta)], then on [-300, log10(tau/eta)]
/// (Tandem DieterichRuinaBase.h:90-132).  Both brackets failing sets
/// kRSNoBracket and returns NaN, as the host code does.
MFEM_HOST_DEVICE inline real_t RSSolveSlipRatePsi(real_t tau, real_t psi,
                                                  real_t sigma_n, real_t eta,
                                                  real_t a, real_t V0,
                                                  int &status)
{
   if (sigma_n <= 0.0)
   {
      if (eta > 0.0) { return tau / eta; }
      else { return 0.0; }
   }

   if (tau <= 0.0) { return 0.0; }

   if (eta == 0.0)
   {
      // No radiation damping: direct inversion (Tandem: Finv path)
      const real_t r = tau / sigma_n;
      const real_t V = V0 * (std::exp((r - psi) / a)
                             - std::exp(-(r + psi) / a));
      return (V < 0.0) ? 0.0 : V;  // std::max(V, 0.0)
   }

   auto fF = [&](real_t Ve) -> real_t
   {
      const real_t V = std::pow(10.0, Ve);
      const real_t f_val = RSFrictionCoefficientPsi(V, psi, a, V0);
      const real_t result = tau - sigma_n * f_val - eta * V;
      if (!std::isfinite(result)) { status |= kRSResidualNonFinite; }
      return result;
   };

   const real_t Va = -32.0;
   const real_t Vb = std::log10(tau / eta);
   const real_t Va_min = -300.0;
   real_t Flo = kRSNaN, Fhi = kRSNaN, lo = kRSNaN, hi = kRSNaN;

   auto try_bracket = [&](real_t a_try, real_t b_try, real_t &V_try) -> bool
   {
      lo = (b_try < a_try) ? b_try : a_try;  // std::min(a_try, b_try)
      hi = (a_try < b_try) ? b_try : a_try;  // std::max(a_try, b_try)
      Flo = fF(lo);
      Fhi = fF(hi);

      if (!std::isfinite(Flo) || !std::isfinite(Fhi)) { return false; }

      if (Flo == 0.0)
      {
         V_try = std::pow(10.0, lo);
         return true;
      }
      if (Fhi == 0.0)
      {
         V_try = std::pow(10.0, hi);
         return true;
      }

      if (std::copysign(Flo, Fhi) == Flo) { return false; }

      const real_t Ve = RSZeroIn(lo, hi, fF, status);
      V_try = std::pow(10.0, Ve);
      return std::isfinite(V_try);
   };

   real_t V_try = kRSNaN;
   if (try_bracket(Va, Vb, V_try)) { return V_try; }
   // Fallback bracket (Tandem: lines 121-129)
   if (try_bracket(Va_min, Vb, V_try)) { return V_try; }

   status |= kRSNoBracket;
   return kRSNaN;
}

/// Slip rate vector from the traction vector: |V| from RSSolveSlipRatePsi at
/// |tau|, direction anti-parallel to tau (Tandem DieterichRuinaBase.h:174).
MFEM_HOST_DEVICE inline void RSSolveSlipRateVectorPsi(
   real_t tau0, real_t tau1, real_t psi, real_t sigma_n, real_t eta,
   real_t a, real_t V0, real_t &V_0, real_t &V_1, int &status)
{
   const real_t tau_abs = std::sqrt(tau0 * tau0 + tau1 * tau1);
   if (tau_abs <= 0.0)
   {
      V_0 = 0.0;
      V_1 = 0.0;
      return;
   }
   const real_t V_abs = RSSolveSlipRatePsi(tau_abs, psi, sigma_n, eta, a, V0,
                                           status);
   const real_t inv_tau = 1.0 / tau_abs;
   V_0 = -(tau0 * inv_tau) * V_abs;
   V_1 = -(tau1 * inv_tau) * V_abs;
}

/// Aging law in psi form: dpsi/dt = (b*V0/Dc) * [exp((f0-psi)/b) - V/V0].
MFEM_HOST_DEVICE inline real_t RSAgingRatePsi(real_t V, real_t psi, real_t Dc,
                                              const RSDeviceParams &prm)
{
   const real_t exp_arg = (prm.aging_f0 - psi) / prm.aging_b;
   return (prm.aging_b * prm.aging_V0 / Dc) *
          (std::exp(exp_arg) - V / prm.aging_V0);
}

/// One fault node of RateStateFaultOperator::ComputeRHS (BP5 path): the
/// input checks, the slip rate and the state rate.  Writes V and dpsi/dt and
/// returns |V| (0 when the inputs are not finite: the host stops the run
/// before solving there).
MFEM_HOST_DEVICE inline real_t RSNodeRHS(
   real_t psi, real_t tau0, real_t tau1, real_t sigma_n_eff, real_t a,
   real_t eta, real_t Dc, const RSDeviceParams &prm,
   real_t &V_0, real_t &V_1, real_t &dpsi, int &status)
{
   if (!std::isfinite(psi) || !std::isfinite(tau0) || !std::isfinite(tau1) ||
       !std::isfinite(sigma_n_eff) || !std::isfinite(a) ||
       !std::isfinite(eta) || !std::isfinite(Dc))
   {
      status |= kRSBadInput;
      V_0 = kRSNaN;
      V_1 = kRSNaN;
      dpsi = kRSNaN;
      return 0.0;
   }
   if (sigma_n_eff <= 0.0) { status |= kRSSigmaNNonPositive; }

   RSSolveSlipRateVectorPsi(tau0, tau1, psi, sigma_n_eff, eta, a, prm.V0,
                            V_0, V_1, status);
   const real_t V_abs = std::sqrt(V_0 * V_0 + V_1 * V_1);
   if (!std::isfinite(V_0) || !std::isfinite(V_1) || !std::isfinite(V_abs))
   {
      status |= kRSNonFiniteV;
   }
   dpsi = RSAgingRatePsi(V_abs, psi, Dc, prm);
   return V_abs;
}

/// Largest |V| over the nodes, and the OR of their status words.
using RSReduceValue = DevicePair<real_t, int>;

/// Reducer for the fault kernels: largest |V| and the OR of the nodes'
/// status words.  fmax skips a NaN, as std::max(V_max, V) did on the host.
struct RSMaxVStatusReducer
{
   using value_type = RSReduceValue;
   static MFEM_HOST_DEVICE void Join(value_type &a, const value_type &b)
   {
      a.first = fmax(a.first, b.first);
      a.second |= b.second;
   }
   static MFEM_HOST_DEVICE void SetInitialValue(value_type &a)
   {
      a.first = 0.0;
      a.second = 0;
   }
};

/// Where a kernel finds the fault's per-node inputs (BP5 layout: state
/// [slip_dip, slip_strike, psi] per node, traction, pre-stress and slip rate
/// [dip, strike] per node).  normal and sigma_n_dof may be null: then
/// sigma_n_eff is sigma_n (+ the normal traction when given), as on the host.
struct RSFaultInputs
{
   const Vector *traction;
   const Vector *state;
   const Vector *normal;       ///< elastic normal traction, or nullptr
   const Vector *tau_pre;
   const Vector *sigma_n_dof;  ///< per-node pre-stress sigma_n, or nullptr
   real_t sigma_n;             ///< pre-stress sigma_n without sigma_n_dof
   const Vector *a;
   const Vector *eta;
};

/// RateStateFaultOperator::ComputeRHS (BP5 path) over N nodes, in one
/// reduction kernel: writes rate [V_dip, V_strike, dpsi/dt] and slip_rate
/// [V_dip, V_strike] per node and returns {max |V|, OR of status}.
/// use_dev picks the device (as mfem::reduce does); on the host the nodes
/// run in order, as the host loop did.
inline RSReduceValue RSFaultRHS(int N, bool use_dev, const RSDeviceParams &prm,
                                const RSFaultInputs &in, const Vector &Dc,
                                Vector &rate, Vector &slip_rate,
                                Array<RSReduceValue> &workspace)
{
   const real_t *T = in.traction->Read(use_dev);
   const real_t *X = in.state->Read(use_dev);
   const real_t *TN = in.normal ? in.normal->Read(use_dev) : nullptr;
   const real_t *TP = in.tau_pre->Read(use_dev);
   const real_t *SN = in.sigma_n_dof ? in.sigma_n_dof->Read(use_dev) : nullptr;
   const real_t sn0 = in.sigma_n;
   const real_t *A = in.a->Read(use_dev);
   const real_t *ETA = in.eta->Read(use_dev);
   const real_t *DC = Dc.Read(use_dev);
   real_t *R = rate.Write(use_dev);
   real_t *VS = slip_rate.Write(use_dev);
   const RSDeviceParams p = prm;

   RSReduceValue res;
   res.first = 0.0;
   res.second = 0;
   reduce(N, res, [=] MFEM_HOST_DEVICE (int i, RSReduceValue &r)
   {
      const real_t psi = X[3*i + 2];
      const real_t tau0 = TP[2*i] + T[2*i];
      const real_t tau1 = TP[2*i + 1] + T[2*i + 1];
      const real_t sigma_n_pre = SN ? SN[i] : sn0;
      const real_t sigma_n_eff = TN ? sigma_n_pre + TN[i] : sigma_n_pre;
      int status = 0;
      real_t V_0, V_1, dpsi;
      const real_t V_abs = RSNodeRHS(psi, tau0, tau1, sigma_n_eff, A[i],
                                     ETA[i], DC[i], p, V_0, V_1, dpsi, status);
      R[3*i] = V_0;
      R[3*i + 1] = V_1;
      R[3*i + 2] = dpsi;
      VS[2*i] = V_0;
      VS[2*i + 1] = V_1;
      r.first = fmax(r.first, V_abs);
      r.second |= status;
   }, RSMaxVStatusReducer{}, use_dev, workspace);
   return res;
}

/// RateStateFaultOperator::RecomputeSlipRate (BP5 path) over N nodes: the
/// slip rate alone, without the input checks of the RHS (as on the host).
/// Writes V [V_dip, V_strike] per node and returns the OR of the status
/// words.
inline int RSFaultSlipRate(int N, bool use_dev, real_t V0,
                           const RSFaultInputs &in, Vector &V,
                           Array<RSReduceValue> &workspace)
{
   const real_t *T = in.traction->Read(use_dev);
   const real_t *X = in.state->Read(use_dev);
   const real_t *TN = in.normal ? in.normal->Read(use_dev) : nullptr;
   const real_t *TP = in.tau_pre->Read(use_dev);
   const real_t *SN = in.sigma_n_dof ? in.sigma_n_dof->Read(use_dev) : nullptr;
   const real_t sn0 = in.sigma_n;
   const real_t *A = in.a->Read(use_dev);
   const real_t *ETA = in.eta->Read(use_dev);
   real_t *VO = V.Write(use_dev);

   RSReduceValue res;
   res.first = 0.0;
   res.second = 0;
   reduce(N, res, [=] MFEM_HOST_DEVICE (int i, RSReduceValue &r)
   {
      const real_t psi = X[3*i + 2];
      const real_t tau0 = TP[2*i] + T[2*i];
      const real_t tau1 = TP[2*i + 1] + T[2*i + 1];
      const real_t sigma_n_pre = SN ? SN[i] : sn0;
      const real_t sigma_n_eff = TN ? sigma_n_pre + TN[i] : sigma_n_pre;
      int status = 0;
      real_t V_0, V_1;
      RSSolveSlipRateVectorPsi(tau0, tau1, psi, sigma_n_eff, ETA[i], A[i], V0,
                               V_0, V_1, status);
      VO[2*i] = V_0;
      VO[2*i + 1] = V_1;
      r.second |= status;
   }, RSMaxVStatusReducer{}, use_dev, workspace);
   return res.second;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DEVICE_RATE_STATE_HPP
