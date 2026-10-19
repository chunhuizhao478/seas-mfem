// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// test_device_friction.cpp — GPU plan G4 acceptance: the fault kernels of
// friction/device_rate_state.hpp against the host friction code.
//
//   ./seas_test_device_friction [--device cuda] [--n 1000000] [--seed 1]
//                               [--tol 1e-12]
//
// Three input sets, --n fault nodes each (edges: 4000), go
// through RSFaultRHS (the kernel of RateStateFaultOperator::ComputeRHS) and
// RSFaultSlipRate (the kernel of RecomputeSlipRate), on the device when
// --device names one:
//   random  psi, tau, sigma_n, the elastic normal traction, a, eta and Dc
//           drawn uniformly over wide ranges (per-node sigma_n);
//   bp5     BP5 constants, a between a0 and amax, psi at steady state for a
//           slip rate drawn log-uniformly in [1e-12, 1] m/s, and tau giving
//           a slip rate drawn log-uniformly in [1e-20, 10] m/s (sigma_n
//           25 MPa plus an elastic normal traction);
//   edges   tau = 0, eta = 0, sigma_n_eff <= 0 and psi / a > 700, a
//           thousand nodes each.
// The reference is the host code ComputeRHS called before G4, node by node:
// DieterichRuinaFriction::SolveSlipRateVectorPsi and AgingLawPsi::Rate.
// Checks, per set:
//   1. |V_kernel - V_host| / |V_host| <= tol at every node, for both kernels
//      (NaN where the host gives NaN);
//   2. |dpsi_kernel - dpsi_host| <= tol times the size of the two terms of
//      the aging law, (b V0 / Dc) (exp((f0 - psi) / b) + |V| / V0);
//   3. the kernel's V_max equals the host's max |V| within tol;
//   4. no status bit after which ComputeRHS stops the run, and the
//      sigma_n_eff <= 0 bit exactly when such a node is in the set.
// A last check feeds one non-finite psi and expects the bad-input bit.
// On the CPU the kernels repeat the host arithmetic and every difference is
// 0; on a GPU, exp, pow and log10 may differ in the last bit.  Printed: the
// largest differences, the share of bit-identical slip rates and the time of
// each kernel against the host loop.  Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../friction/device_rate_state.hpp"
#include "../../config/bp5_params.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
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

void DeviceSync()
{
   if (Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK))
   {
      MFEM_DEVICE_SYNC;
   }
}

double Seconds(std::chrono::steady_clock::time_point t0)
{
   return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        t0).count();
}

/// One input set: the fault-node vectors a kernel reads.
struct InputSet
{
   std::string name;
   int n = 0;
   Vector traction, state, normal, tau_pre, sigma_n_dof, a, eta, Dc;
   bool use_normal = false;
   bool use_sigma_n_dof = false;
   real_t sigma_n = 0.0;
   bool expect_sigma_n_nonpositive = false;

   void Resize(int n_)
   {
      n = n_;
      for (Vector *v : {&traction, &tau_pre}) { v->SetSize(2 * n); }
      state.SetSize(3 * n);
      for (Vector *v : {&normal, &sigma_n_dof, &a, &eta, &Dc})
      {
         v->SetSize(n);
      }
      for (Vector *v : {&traction, &state, &normal, &tau_pre, &sigma_n_dof,
                        &a, &eta, &Dc})
      {
         v->UseDevice(true);
         v->HostWrite();
      }
   }

   RSFaultInputs Inputs() const
   {
      RSFaultInputs in;
      in.traction = &traction;
      in.state = &state;
      in.normal = use_normal ? &normal : nullptr;
      in.tau_pre = &tau_pre;
      in.sigma_n_dof = use_sigma_n_dof ? &sigma_n_dof : nullptr;
      in.sigma_n = sigma_n;
      in.a = &a;
      in.eta = &eta;
      return in;
   }

   /// sigma_n_eff at node i, as ComputeRHS forms it.
   real_t SigmaNEff(int i) const
   {
      const real_t pre = use_sigma_n_dof ? sigma_n_dof(i) : sigma_n;
      return use_normal ? pre + normal(i) : pre;
   }
};

/// Shear traction split between the dip and strike components at a random
/// angle, as tau_pre + traction (tau_pre takes a random share).
void SetShear(InputSet &s, int i, real_t tau_abs, std::mt19937_64 &gen)
{
   std::uniform_real_distribution<real_t> angle(0.0, 2.0 * M_PI);
   std::uniform_real_distribution<real_t> share(0.0, 1.0);
   const real_t th = angle(gen);
   const real_t t0 = tau_abs * std::cos(th);
   const real_t t1 = tau_abs * std::sin(th);
   const real_t w = share(gen);
   s.tau_pre(2*i) = w * t0;
   s.tau_pre(2*i + 1) = w * t1;
   s.traction(2*i) = t0 - s.tau_pre(2*i);
   s.traction(2*i + 1) = t1 - s.tau_pre(2*i + 1);
}

void FillRandom(InputSet &s, int n, std::mt19937_64 &gen)
{
   s.name = "random";
   s.Resize(n);
   s.use_normal = true;
   s.use_sigma_n_dof = true;
   std::uniform_real_distribution<real_t> u(0.0, 1.0);
   for (int i = 0; i < n; i++)
   {
      s.state(3*i) = 20.0 * (u(gen) - 0.5);
      s.state(3*i + 1) = 20.0 * (u(gen) - 0.5);
      s.state(3*i + 2) = 0.2 + 0.8 * u(gen);            // psi
      s.a(i) = 0.004 + 0.036 * u(gen);
      s.eta(i) = 1.0e6 + 9.0e6 * u(gen);
      s.Dc(i) = 0.01 + 0.2 * u(gen);
      s.sigma_n_dof(i) = 10.0e6 + 40.0e6 * u(gen);
      s.normal(i) = 10.0e6 * (u(gen) - 0.5);            // keeps sigma_n > 0
      SetShear(s, i, 60.0e6 * u(gen), gen);
   }
}

void FillBP5(InputSet &s, int n, std::mt19937_64 &gen)
{
   const BP5Params p;
   s.name = "bp5";
   s.Resize(n);
   s.use_normal = true;
   s.use_sigma_n_dof = false;
   s.sigma_n = p.sigma_n;
   std::uniform_real_distribution<real_t> u(0.0, 1.0);
   std::normal_distribution<real_t> dsn(0.0, 2.0e6);
   const DieterichRuinaFriction::Constants fc{p.V0, p.f0, p.b, p.L0};
   DieterichRuinaFriction fr(fc);
   for (int i = 0; i < n; i++)
   {
      const real_t a = p.a0 + (p.amax - p.a0) * u(gen);
      const real_t V_ss = std::pow(10.0, -12.0 + 12.0 * u(gen));
      const real_t psi = p.f0 + p.b * std::log(p.V0 / V_ss);
      const real_t V = std::pow(10.0, -20.0 + 21.0 * u(gen));
      real_t sn = p.sigma_n + dsn(gen);
      sn = std::max(sn, 1.0e6);
      s.a(i) = a;
      s.eta(i) = p.eta();
      s.Dc(i) = (u(gen) < 0.5) ? p.L_nuc : p.L0;
      s.normal(i) = sn - p.sigma_n;
      s.state(3*i) = 5.0 * (u(gen) - 0.5);
      s.state(3*i + 1) = 5.0 * (u(gen) - 0.5);
      s.state(3*i + 2) = psi;
      const real_t tau = sn * fr.FrictionCoefficientPsi(V, psi, a) +
                         p.eta() * V;
      SetShear(s, i, tau, gen);
   }
}

void FillEdges(InputSet &s, std::mt19937_64 &gen)
{
   const BP5Params p;
   const int per = 1000;
   s.name = "edges";
   s.Resize(4 * per);
   s.use_normal = false;
   s.use_sigma_n_dof = true;
   s.expect_sigma_n_nonpositive = true;
   std::uniform_real_distribution<real_t> u(0.0, 1.0);
   for (int i = 0; i < s.n; i++)
   {
      const int kind = i / per;
      s.state(3*i) = 0.0;
      s.state(3*i + 1) = 0.0;
      s.state(3*i + 2) = 0.4 + 0.5 * u(gen);
      s.a(i) = p.a0 + (p.amax - p.a0) * u(gen);
      s.eta(i) = p.eta();
      s.Dc(i) = p.L0;
      s.sigma_n_dof(i) = p.sigma_n;
      real_t tau = 30.0e6 * u(gen);
      switch (kind)
      {
         case 0: tau = 0.0; break;                         // tau = 0
         case 1: s.eta(i) = 0.0; break;                    // eta = 0
         case 2: s.sigma_n_dof(i) = -5.0e6 * u(gen); break; // sigma_n <= 0
         default:                                          // psi / a > 700
            // tau above sigma_n f(1e-300 m/s, psi), so the second bracket
            // holds the root
            s.a(i) = p.a0;
            s.state(3*i + 2) = p.a0 * (700.5 + 50.0 * u(gen));
            tau = 10.0e6 + 20.0e6 * u(gen);
            break;
      }
      SetShear(s, i, tau, gen);
   }
}

/// The host code of ComputeRHS before G4, node by node.
struct HostResult
{
   Vector V, dpsi, scale;
   real_t V_max = 0.0;
   bool any_sigma_n_nonpositive = false;
};

HostResult RunHost(const InputSet &s, const DieterichRuinaFriction &fr,
                   const AgingLawPsi &aging, double &seconds)
{
   for (const Vector *v : {&s.traction, &s.state, &s.normal, &s.tau_pre,
                           &s.sigma_n_dof, &s.a, &s.eta, &s.Dc})
   {
      v->HostRead();
   }
   HostResult h;
   h.V.SetSize(2 * s.n);
   h.dpsi.SetSize(s.n);
   h.scale.SetSize(s.n);
   const auto t0 = std::chrono::steady_clock::now();
   for (int i = 0; i < s.n; i++)
   {
      const real_t psi = s.state(3*i + 2);
      const real_t tau_vec[2] = {s.tau_pre(2*i) + s.traction(2*i),
                                 s.tau_pre(2*i + 1) + s.traction(2*i + 1)};
      const real_t sn = s.SigmaNEff(i);
      if (sn <= 0.0) { h.any_sigma_n_nonpositive = true; }
      real_t V_vec[2];
      fr.SolveSlipRateVectorPsi(tau_vec, psi, sn, s.eta(i), s.a(i), V_vec);
      const real_t V_abs = std::sqrt(V_vec[0]*V_vec[0] + V_vec[1]*V_vec[1]);
      h.V(2*i) = V_vec[0];
      h.V(2*i + 1) = V_vec[1];
      h.dpsi(i) = aging.Rate(V_abs, psi, s.Dc(i));
      h.V_max = std::max(h.V_max, V_abs);
   }
   seconds = Seconds(t0);
   for (int i = 0; i < s.n; i++)
   {
      const real_t psi = s.state(3*i + 2);
      const real_t V_abs = std::hypot(h.V(2*i), h.V(2*i + 1));
      h.scale(i) = (aging.GetB() * aging.GetV0() / s.Dc(i)) *
                   (std::exp((aging.GetF0() - psi) / aging.GetB()) +
                    V_abs / aging.GetV0());
   }
   return h;
}

/// Largest per-node relative difference of slip-rate vectors, and the count
/// of bit-identical nodes.  NaN on both sides counts as equal.
real_t MaxRelV(const Vector &Vk, const Vector &Vh, int n, int &identical)
{
   Vk.HostRead();
   Vh.HostRead();
   real_t worst = 0.0;
   identical = 0;
   for (int i = 0; i < n; i++)
   {
      const real_t k0 = Vk(2*i), k1 = Vk(2*i + 1);
      const real_t h0 = Vh(2*i), h1 = Vh(2*i + 1);
      const bool hn = std::isnan(h0) || std::isnan(h1);
      const bool kn = std::isnan(k0) || std::isnan(k1);
      if (hn || kn)
      {
         if (hn && kn) { identical++; }
         else { worst = std::numeric_limits<real_t>::infinity(); }
         continue;
      }
      if (k0 == h0 && k1 == h1) { identical++; continue; }
      const real_t d = std::hypot(k0 - h0, k1 - h1);
      const real_t m = std::hypot(h0, h1);
      worst = std::max(worst, m > 0.0 ? d / m
                       : std::numeric_limits<real_t>::infinity());
   }
   return worst;
}

void RunSet(InputSet &s, bool use_dev, real_t tol,
            const DieterichRuinaFriction &fr, const AgingLawPsi &aging,
            const RSDeviceParams &prm)
{
   std::cout << "\n[" << s.name << "] " << s.n << " nodes\n";
   double t_host = 0.0;
   const HostResult h = RunHost(s, fr, aging, t_host);

   Vector rate(3 * s.n), slip_rate(2 * s.n), V_re(2 * s.n);
   for (Vector *v : {&rate, &slip_rate, &V_re}) { v->UseDevice(true); }
   Array<RSReduceValue> ws;
   const RSFaultInputs in = s.Inputs();

   // Once to move the inputs to the device, once timed.
   RSReduceValue res = RSFaultRHS(s.n, use_dev, prm, in, s.Dc, rate,
                                  slip_rate, ws);
   DeviceSync();
   auto t0 = std::chrono::steady_clock::now();
   res = RSFaultRHS(s.n, use_dev, prm, in, s.Dc, rate, slip_rate, ws);
   DeviceSync();
   const double t_rhs = Seconds(t0);

   int st_re = RSFaultSlipRate(s.n, use_dev, prm.V0, in, V_re, ws);
   DeviceSync();
   t0 = std::chrono::steady_clock::now();
   st_re = RSFaultSlipRate(s.n, use_dev, prm.V0, in, V_re, ws);
   DeviceSync();
   const double t_re = Seconds(t0);

   // Slip rate from the RHS kernel: rate holds [V_dip, V_strike, dpsi].
   rate.HostRead();
   Vector V_rhs(2 * s.n), dpsi(s.n);
   for (int i = 0; i < s.n; i++)
   {
      V_rhs(2*i) = rate(3*i);
      V_rhs(2*i + 1) = rate(3*i + 1);
      dpsi(i) = rate(3*i + 2);
   }
   int same_rhs = 0, same_sr = 0, same_re = 0;
   const real_t dV_rhs = MaxRelV(V_rhs, h.V, s.n, same_rhs);
   const real_t dV_sr = MaxRelV(slip_rate, h.V, s.n, same_sr);
   const real_t dV_re = MaxRelV(V_re, h.V, s.n, same_re);

   real_t d_dpsi = 0.0;
   for (int i = 0; i < s.n; i++)
   {
      if (std::isnan(h.dpsi(i)) && std::isnan(dpsi(i))) { continue; }
      const real_t d = std::abs(dpsi(i) - h.dpsi(i)) / h.scale(i);
      d_dpsi = std::max(d_dpsi, std::isnan(d)
                        ? std::numeric_limits<real_t>::infinity() : d);
   }
   const real_t d_vmax = h.V_max > 0.0
                         ? std::abs(res.first - h.V_max) / h.V_max
                         : std::abs(res.first);

   std::printf("  max rel |dV|: RHS %.3e, cached slip rate %.3e, "
               "RecomputeSlipRate %.3e\n", dV_rhs, dV_sr, dV_re);
   std::printf("  bit-identical V: %d / %d (RHS), %d / %d (RecomputeSlipRate)\n",
               same_rhs, s.n, same_re, s.n);
   std::printf("  max |d dpsi/dt| / scale %.3e, rel |d V_max| %.3e "
               "(V_max %.6e)\n", d_dpsi, d_vmax, h.V_max);
   std::printf("  status: RHS %d, RecomputeSlipRate %d\n", res.second, st_re);
   std::printf("  time: host loop %.3f s, RHS kernel %.3f s, "
               "RecomputeSlipRate kernel %.3f s\n", t_host, t_rhs, t_re);

   Check(dV_rhs <= tol && dV_sr <= tol && dV_re <= tol,
         s.name + ": slip rate within tol of the host at every node");
   Check(d_dpsi <= tol, s.name + ": state rate within tol");
   Check(d_vmax <= tol, s.name + ": V_max within tol");
   Check((res.second & kRSAbortRHS) == 0 && (st_re & kRSAbortSolve) == 0,
         s.name + ": no status bit that stops the run");
   const bool sn_bit = (res.second & kRSSigmaNNonPositive) != 0;
   Check(sn_bit == h.any_sigma_n_nonpositive &&
         sn_bit == s.expect_sigma_n_nonpositive,
         s.name + ": sigma_n_eff <= 0 bit as expected");
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);

   std::string device_str = "cpu";
   int n = 1000000;
   int seed = 1;
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
   const bool use_dev = Device::Allows(Backend::DEVICE_MASK);

   const BP5Params p;
   const DieterichRuinaFriction::Constants fc{p.V0, p.f0, p.b, p.L0};
   const DieterichRuinaFriction fr(fc);
   const AgingLawPsi aging(p.b, p.V0, p.f0);
   RSDeviceParams prm;
   prm.V0 = fc.V0;
   prm.aging_b = aging.GetB();
   prm.aging_V0 = aging.GetV0();
   prm.aging_f0 = aging.GetF0();

   std::cout << "test_device_friction: n " << n << ", seed " << seed
             << ", tol " << tol << ", kernels on the "
             << (use_dev ? "device" : "host") << "\n";

   std::mt19937_64 gen(static_cast<std::uint64_t>(seed));
   {
      InputSet s;
      FillRandom(s, n, gen);
      RunSet(s, use_dev, tol, fr, aging, prm);
   }
   {
      InputSet s;
      FillBP5(s, n, gen);
      RunSet(s, use_dev, tol, fr, aging, prm);
   }
   {
      InputSet s;
      FillEdges(s, gen);
      RunSet(s, use_dev, tol, fr, aging, prm);
   }

   // A non-finite input sets the bad-input bit (ComputeRHS stops the run).
   {
      std::cout << "\n[bad input] 4 nodes, one psi NaN\n";
      InputSet s;
      FillBP5(s, 4, gen);
      s.state(3*2 + 2) = std::numeric_limits<real_t>::quiet_NaN();
      Vector rate(12), slip_rate(8);
      rate.UseDevice(true);
      slip_rate.UseDevice(true);
      Array<RSReduceValue> ws;
      const RSReduceValue res = RSFaultRHS(4, use_dev, prm, s.Inputs(), s.Dc,
                                           rate, slip_rate, ws);
      rate.HostRead();
      std::printf("  status %d, V at that node (%g, %g)\n", res.second,
                  rate(6), rate(7));
      Check((res.second & kRSBadInput) != 0 && std::isnan(rate(6)),
            "bad input: status bit set and V is NaN");
   }

   std::cout << "\n" << (failures ? "FAILED" : "PASSED") << " ("
             << failures << " failed checks)\n";
   return failures ? 1 : 0;
}
