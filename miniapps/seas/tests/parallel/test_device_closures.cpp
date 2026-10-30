// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_device_closures.cpp — GPU dynamic-rupture plan D3 acceptance: the
// device closures of friction/device_closures.hpp against the host closures
// of the Nitsche fault (IFaultClosure).
//
//   ./seas_test_device_closures [--device cuda|debug] [--n 1000000] [--seed 1]
//                               [--tol 1e-12]
//
// For each law with a device form (linear slip weakening, rate and state
// with the aging law, the slip law and strong rate weakening) --n fault
// points with random inputs go through one kernel, on the device when
// --device names one: DCSolveSlipRate (box 9) and DCStateRate at the
// returned |V| and at a second, independent (V, psi) pair (boxes 4, 10).
// The reference is the host closure, point by point: SolveSlipRate and
// StateRate with a FaultPointState holding the same values.  Inputs:
//   slip weakening  X_tan in [-60, 60] MPa per component, sigma_n in
//                   [0, 120] MPa, |d^n| in [0, 2 Dc], mu_s in [0.5, 0.8],
//                   mu_d in [0.1, mu_s], Dc in [0.1, 1] m, cohesion in
//                   [0, 2] MPa, gamma in [1e5, 1e9] Pa s / m;
//   rate and state  |X_tan| log-uniform in [1e3, 1e8] Pa, sigma_n in [1,
//                   120] MPa, psi in [0, 1.2], a in [0.005, 0.03], Dc in
//                   [0.01, 0.5] m, V_w in [0.05, 0.2] m/s, gamma as above,
//                   and the second rate at V log-uniform in [1e-20, 10] m/s.
// Then benchmark-sampled sets: TPV102 (aging law; b = 0.012, Dc = 0.02 m,
// a in [0.008, 0.016]) and TPV104 (slip law and strong rate weakening; b =
// 0.014, Dc = 0.4 m, a in [0.01, 0.02], V_w in [0.1, 1] m/s, f_w = 0.2),
// sigma_n in 120 +- 20 MPa, psi in [0.3, 1], gamma log-uniform in [1e6,
// 1e8] Pa s / m, and |X_tan| = sigma_n f(V, psi) + gamma V for a slip rate V
// log-uniform in [1e-20, 10] m/s, so every slip rate of a rupture is met;
// the second rate at V in the same range.
// Then an edge set per law, a thousand points of each branch the Nitsche
// fault can reach (gamma > 0 there, so not the rate-and-state eta = 0
// inversion):
//   slip weakening  X_tan = 0; |d^n| = 0 (mu_s); |d^n| = Dc (mu_d);
//                   sigma_n = 0 (cohesion only); |X_tan| = tau_s exactly
//                   (stick at the threshold);
//   rate and state  X_tan = 0; sigma_n <= 0 (V = |X_tan| / gamma);
//                   psi / a > 700 with the root in the second bracket (the
//                   log form of the friction coefficient); the second rate
//                   at V = 0 and at V = 1e-60 (the slip law's 1e-50 floor).
// Checks, per law and set:
//   1. |V_device - V_host| <= tol |V_host| per component and point (NaN
//      where the host gives NaN);
//   2. |r_device - r_host| <= tol (|r_host| + s + T), T the size of the
//      two terms whose difference is r at that point (StateRateTerms) and s
//      a fixed floor per law;
//   3. no RSStatus bit after which the host stops the run.
// On the CPU (and the debug device) the kernels repeat the host arithmetic
// and every difference must be exactly 0; on a GPU exp, log, pow, asinh and
// expm1 may differ in the last bit.  Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../friction/device_closures.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int failures = 0;

void Check(bool ok, const std::string &what)
{
   std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
   if (!ok) { failures++; }
}

/// Inputs of n points: X (2), sigma, psi, slip_abs, gamma, the parameters
/// (kDCParams) and the second rate's (V, psi).
struct Inputs
{
   int n = 0;
   Vector X, sig, psi, sabs, gam, q, V2, psi2;
};

/// The kernel under test: V (2 per point), r(|V|, psi), r(V2, psi2) and
/// the OR of the status words.
int RunKernel(const DeviceClosure &cl, const Inputs &in, Vector &V, Vector &r1, Vector &r2)
{
   const int n = in.n;
   V.SetSize(2 * n);
   r1.SetSize(n);
   r2.SetSize(n);
   Vector st(n);
   const real_t *X = in.X.Read(), *sig = in.sig.Read(), *psi = in.psi.Read(),
                 *sabs = in.sabs.Read(), *gam = in.gam.Read(), *q = in.q.Read(),
                 *V2 = in.V2.Read(), *psi2 = in.psi2.Read();
   real_t *dV = V.Write(), *d1 = r1.Write(), *d2 = r2.Write(), *ds = st.Write();
   mfem::forall(n, [=] MFEM_HOST_DEVICE (int i)
   {
      int status = 0;
      real_t Vl[2];
      const real_t *qi = q + kDCParams * i;
      DCSolveSlipRate(cl, qi, X + 2 * i, gam[i], sig[i], psi[i], sabs[i], Vl, status);
      dV[2 * i] = Vl[0];
      dV[2 * i + 1] = Vl[1];
      d1[i] = DCStateRate(cl, qi, std::sqrt(Vl[0] * Vl[0] + Vl[1] * Vl[1]), psi[i]);
      d2[i] = DCStateRate(cl, qi, V2[i], psi2[i]);
      ds[i] = static_cast<real_t>(status);
   });
   const real_t *hs = st.HostRead();
   int all = 0;
   for (int i = 0; i < n; i++) { all |= static_cast<int>(hs[i]); }
   return all;
}

bool Close(real_t a, real_t b, real_t tol, real_t scale)
{
   if (std::isnan(b)) { return std::isnan(a); }
   // equal, infinities included (the state rate at a |V| that underflows)
   if (a == b) { return true; }
   return std::abs(a - b) <= tol * (std::abs(b) + scale);
}

/// The same value, NaN where the other is NaN (the strong-rate-weakening
/// state rate at |V| = 0 is 0 * inf on the host too).
bool Same(real_t a, real_t b) { return a == b || (std::isnan(a) && std::isnan(b)); }

/// The size of the two terms whose difference is the state rate at (V, psi):
/// (b V0 / Dc)(exp((f0 - psi) / b) + V / V0) for the aging law and
/// (V / Dc)(|psi| + |psi_ss(V)|) for the slip law and strong rate
/// weakening; 0 for slip weakening (or a term that is not finite).  Near
/// steady state the rate is far smaller than its terms, so a last-bit
/// difference of exp or log on a GPU is bounded by the terms, not the rate.
real_t StateRateTerms(const DeviceClosure &c, const real_t *q, real_t V, real_t psi)
{
   real_t t = 0.0;
   switch (c.law)
   {
      case DeviceClosureLaw::Aging:
         t = (c.rs.aging_b * c.rs.aging_V0 / q[kDCDc])
             * (std::exp((c.rs.aging_f0 - psi) / c.rs.aging_b) + V / c.rs.aging_V0);
         break;
      case DeviceClosureLaw::Slip:
      {
         V = (V < 1e-50) ? 1e-50 : V;
         t = (V / q[kDCDc]) * (std::abs(psi) + std::abs(c.f0 + c.b * std::log(c.V0 / V)));
         break;
      }
      case DeviceClosureLaw::SRW:
         t = (V / q[kDCDc])
             * (std::abs(psi)
                + std::abs(DCPsiSSSRW(V, q[kDCV_w], q[kDCa], c.b, c.V0, c.f0, c.muW)));
         break;
      default: break;
   }
   return std::isfinite(t) ? t : 0.0;
}

void RunLaw(const std::string &name, const IFaultClosure &host, const Inputs &in, real_t tol,
            bool use_dev, real_t rate_scale)
{
   DeviceClosure cl;
   Check(MakeDeviceClosure(host, cl), name + ": the closure has a device form");
   const auto t0 = std::chrono::steady_clock::now();
   Vector V, r1, r2;
   const int status = RunKernel(cl, in, V, r1, r2);
   const double kt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
   const real_t *hV = V.HostRead(), *h1 = r1.HostRead(), *h2 = r2.HostRead();
   const real_t *X = in.X.HostRead(), *sig = in.sig.HostRead(), *psi = in.psi.HostRead(),
                 *sabs = in.sabs.HostRead(), *gam = in.gam.HostRead(), *q = in.q.HostRead(),
                 *V2 = in.V2.HostRead(), *psi2 = in.psi2.HostRead();
   long badV = 0, bad1 = 0, bad2 = 0, same = 0, nonfinite = 0;
   real_t errV = 0.0, err1 = 0.0;
   const auto t1 = std::chrono::steady_clock::now();
   for (int i = 0; i < in.n; i++)
   {
      FaultPointState s;
      const real_t *qi = q + kDCParams * i;
      s.a = qi[kDCa];
      s.Dc = qi[kDCDc];
      s.mu_s = qi[kDCmu_s];
      s.mu_d = qi[kDCmu_d];
      s.cohesion = qi[kDCcohesion];
      s.V_w = qi[kDCV_w];
      s.psi = psi[i];
      s.slip_abs = sabs[i];
      real_t Vh[2];
      host.SolveSlipRate(X + 2 * i, gam[i], sig[i], s, Vh);
      const real_t rh1 = host.StateRate(std::sqrt(Vh[0] * Vh[0] + Vh[1] * Vh[1]), s);
      FaultPointState s2 = s;
      s2.psi = psi2[i];
      const real_t rh2 = host.StateRate(V2[i], s2);
      const real_t vs = std::abs(Vh[0]) + std::abs(Vh[1]);
      for (int c = 0; c < 2; c++)
      {
         if (!Close(hV[2 * i + c], Vh[c], tol, 0.0) && !(vs == 0.0 && hV[2 * i + c] == 0.0))
         {
            badV++;
         }
         if (vs > 0.0) { errV = std::max(errV, std::abs(hV[2 * i + c] - Vh[c]) / vs); }
      }
      if (!std::isfinite(rh1) || !std::isfinite(rh2)) { nonfinite++; }
      const real_t sc1 = rate_scale
                         + StateRateTerms(cl, qi, std::sqrt(Vh[0] * Vh[0] + Vh[1] * Vh[1]),
                                          psi[i]);
      const real_t sc2 = rate_scale + StateRateTerms(cl, qi, V2[i], psi2[i]);
      if (!Close(h1[i], rh1, tol, sc1)) { bad1++; }
      if (!Close(h2[i], rh2, tol, sc2)) { bad2++; }
      if (std::isfinite(h1[i]) && std::isfinite(rh1))
      {
         err1 = std::max(err1, std::abs(h1[i] - rh1) / (std::abs(rh1) + sc1));
      }
      if (Same(hV[2 * i], Vh[0]) && Same(hV[2 * i + 1], Vh[1]) && Same(h1[i], rh1)
          && Same(h2[i], rh2)) { same++; }
   }
   const double ht = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
   std::cout << name << ": " << in.n << " points; largest relative V difference " << errV
             << ", state rate " << err1 << "; bit-identical points " << same << " / " << in.n
             << "; host state rate not finite at " << nonfinite << " points; kernel " << kt
             << " s, host loop " << ht << " s\n";
   Check(badV == 0, name + ": slip rate within tol (" + std::to_string(badV) + " outside)");
   Check(bad1 == 0, name + ": state rate at |V| within tol (" + std::to_string(bad1)
         + " outside)");
   Check(bad2 == 0, name + ": state rate at (V2, psi2) within tol (" + std::to_string(bad2)
         + " outside)");
   Check(!(status & kRSAbortSolve) && !(status & kRSNoBracket),
         name + ": no status bit that stops the run (status " + std::to_string(status) + ")");
   if (!use_dev)
   {
      Check(same == in.n, name + ": every point bit-identical on the CPU");
   }
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   std::string device_str = "cpu";
   int n = 1000000, seed = 1;
   real_t tol = 1e-12;
   for (int i = 1; i < argc; i++)
   {
      const std::string arg = argv[i];
      if (arg == "--device" && i + 1 < argc) { device_str = argv[++i]; }
      else if (arg == "--n" && i + 1 < argc) { n = std::atoi(argv[++i]); }
      else if (arg == "--seed" && i + 1 < argc) { seed = std::atoi(argv[++i]); }
      else if (arg == "--tol" && i + 1 < argc) { tol = std::atof(argv[++i]); }
      else
      {
         std::cerr << "unknown argument " << arg << "\n";
         return 2;
      }
   }
   Device device(device_str);
   if (Mpi::Root()) { device.Print(); }
   // a GPU may differ in the last bit; the CPU and the debug device may not
   const bool use_dev = Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK);
   std::cout << "test_device_closures: n " << n << ", seed " << seed << ", tol " << tol
             << ", kernels on the " << (Device::IsEnabled() ? "device" : "host") << "\n";

   std::mt19937_64 gen(static_cast<std::uint64_t>(seed));
   auto U = [&](real_t a, real_t b) { return std::uniform_real_distribution<real_t>(a, b)(gen); };
   auto LogU = [&](real_t a, real_t b) { return std::pow(10.0, U(std::log10(a), std::log10(b))); };
   auto Fill = [&](bool rs, Inputs &in, int n)
   {
      in.n = n;
      for (Vector *v : {&in.sig, &in.psi, &in.sabs, &in.gam, &in.V2, &in.psi2}) { v->SetSize(n); }
      in.X.SetSize(2 * n);
      in.q.SetSize(kDCParams * n);
      for (int i = 0; i < n; i++)
      {
         real_t *q = &in.q(kDCParams * i);
         q[kDCa] = U(0.005, 0.03);
         q[kDCDc] = rs ? U(0.01, 0.5) : U(0.1, 1.0);
         q[kDCmu_s] = U(0.5, 0.8);
         q[kDCmu_d] = U(0.1, q[kDCmu_s]);
         q[kDCcohesion] = U(0.0, 2e6);
         q[kDCV_w] = U(0.05, 0.2);
         if (rs)
         {
            const real_t th = U(0.0, 2.0 * M_PI), m = LogU(1e3, 1e8);
            in.X(2 * i) = m * std::cos(th);
            in.X(2 * i + 1) = m * std::sin(th);
            in.sig(i) = U(1e6, 120e6);
         }
         else
         {
            in.X(2 * i) = U(-60e6, 60e6);
            in.X(2 * i + 1) = U(-60e6, 60e6);
            in.sig(i) = U(0.0, 120e6);
         }
         in.psi(i) = U(0.0, 1.2);
         in.sabs(i) = U(0.0, 2.0 * q[kDCDc]);
         in.gam(i) = LogU(1e5, 1e9);
         in.V2(i) = LogU(1e-20, 10.0);
         in.psi2(i) = U(0.0, 1.2);
      }
   };

   // TPV102 constants (V0, f0, b; Dc is per point); the laws read only these
   const DieterichRuinaFriction::Constants fc{1e-6, 0.6, 0.012, 0.02};
   const DieterichRuinaFriction fr(fc);
   Inputs sw, rs;
   Fill(false, sw, n);
   Fill(true, rs, n);
   // the edge sets: random points, then one branch per thousand
   const int per = 1000;
   Inputs swe, rse;
   Fill(false, swe, 5 * per);
   Fill(true, rse, 5 * per);
   for (int i = 0; i < 5 * per; i++)
   {
      const real_t *q = &swe.q(kDCParams * i);
      switch (i / per)
      {
         case 0: swe.X(2 * i) = 0.0; swe.X(2 * i + 1) = 0.0; break;
         case 1: swe.sabs(i) = 0.0; break;
         case 2: swe.sabs(i) = q[kDCDc]; break;
         case 3: swe.sig(i) = 0.0; break;
         default:
            swe.X(2 * i) = DCSlipWeakeningMu(swe.sabs(i), q[kDCmu_s], q[kDCmu_d], q[kDCDc])
                           * swe.sig(i) + q[kDCcohesion];
            swe.X(2 * i + 1) = 0.0;
            break;
      }
      real_t *qr = &rse.q(kDCParams * i);
      switch (i / per)
      {
         case 0: rse.X(2 * i) = 0.0; rse.X(2 * i + 1) = 0.0; break;
         case 1: rse.sig(i) = (i % 2) ? 0.0 : -U(0.0, 5e6); break;
         case 2:
         {
            // tau above sigma_n f(1e-300 m/s, psi), so the second bracket
            // holds the root (as the edges of test_device_friction)
            qr[kDCa] = 0.005;
            rse.psi(i) = qr[kDCa] * U(700.5, 750.5);
            rse.sig(i) = U(1e6, 20e6);
            const real_t th = U(0.0, 2.0 * M_PI), m = U(10e6, 30e6);
            rse.X(2 * i) = m * std::cos(th);
            rse.X(2 * i + 1) = m * std::sin(th);
            break;
         }
         case 3: rse.V2(i) = 0.0; break;
         default: rse.V2(i) = 1e-60; break;
      }
   }
   // the benchmark-sampled sets
   auto FillBench = [&](const DieterichRuinaFriction &f, real_t a_lo, real_t a_hi, real_t Dc,
                        Inputs &in)
   {
      Fill(true, in, n);
      for (int i = 0; i < n; i++)
      {
         real_t *q = &in.q(kDCParams * i);
         q[kDCa] = U(a_lo, a_hi);
         q[kDCDc] = Dc;
         q[kDCV_w] = U(0.1, 1.0);
         in.sig(i) = 120e6 + U(-20e6, 20e6);
         in.psi(i) = U(0.3, 1.0);
         in.gam(i) = LogU(1e6, 1e8);
         const real_t V = LogU(1e-20, 10.0), th = U(0.0, 2.0 * M_PI);
         const real_t m = in.sig(i) * f.FrictionCoefficientPsi(V, in.psi(i), q[kDCa])
                          + in.gam(i) * V;
         in.X(2 * i) = m * std::cos(th);
         in.X(2 * i + 1) = m * std::sin(th);
         in.V2(i) = LogU(1e-20, 10.0);
         in.psi2(i) = U(0.3, 1.0);
      }
   };
   const DieterichRuinaFriction fr104(DieterichRuinaFriction::Constants{1e-6, 0.6, 0.014, 0.4});
   Inputs b102, b104;
   FillBench(fr, 0.008, 0.016, 0.02, b102);
   FillBench(fr104, 0.01, 0.02, 0.4, b104);

   RunLaw("slip weakening", SlipWeakeningClosure(), sw, tol, use_dev, 1.0);
   RunLaw("slip weakening edges", SlipWeakeningClosure(), swe, tol, use_dev, 1.0);
   const RateStateClosure aging(fr, RateStateClosure::StateLaw::Aging),
         slip(fr, RateStateClosure::StateLaw::Slip);
   const RateStateSRWClosure srw(fr, 0.2);
   RunLaw("aging law", aging, rs, tol, use_dev, fc.b * fc.V0 / 0.01);
   RunLaw("aging law edges", aging, rse, tol, use_dev, fc.b * fc.V0 / 0.01);
   RunLaw("aging law, TPV102 sample", aging, b102, tol, use_dev, fc.b * fc.V0 / 0.02);
   RunLaw("slip law", slip, rs, tol, use_dev, 1e-30);
   RunLaw("slip law edges", slip, rse, tol, use_dev, 1e-30);
   const RateStateClosure slip104(fr104, RateStateClosure::StateLaw::Slip);
   RunLaw("slip law, TPV104 sample", slip104, b104, tol, use_dev, 1e-30);
   RunLaw("strong rate weakening", srw, rs, tol, use_dev, 1e-30);
   RunLaw("strong rate weakening edges", srw, rse, tol, use_dev, 1e-30);
   const RateStateSRWClosure srw104(fr104, 0.2);
   RunLaw("strong rate weakening, TPV104 sample", srw104, b104, tol, use_dev, 1e-30);

   std::cout << (failures == 0 ? "PASSED" : "FAILED") << " (" << failures
             << " failed checks)\n";
   return failures == 0 ? 0 : 1;
}
