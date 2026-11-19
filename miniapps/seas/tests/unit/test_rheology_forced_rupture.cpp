// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Rheology plan rev 1, unit test U11 (no mesh): the TPV26/27 forced rupture
// of friction/forced_rupture.hpp and its branch in SlipWeakeningClosure.
//
//   T(r) at the hypocenter, r_c/2 and r_c against the description's
//   formula (v13, Part 5), the front speed 0.7 Vs at the hypocenter, f2 at
//   the ends and the middle of the ramp, max(f1, f2); the coefficient
//   bitwise equal to SAFS's LSWFrictionCoefficient_ForcedRupture; and with
//   forced_t0 = 0 (no [friction.forced_rupture]) the closure bitwise equal
//   to the slip-weakening arithmetic before forced rupture (TPV205 contract).
//
// Run: ./seas_test_rheology_forced_rupture

#include "mfem.hpp"
#include "../../friction/forced_rupture.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../spatial/code/spatial_friction.hpp"   // the SAFS helper (oracle)
#include "test_macros.hpp"

#include <cmath>
#include <iostream>
#include <random>

using namespace mfem;
using namespace mfem::seas;

namespace
{

// TPV26/27 (description v13, Parts 3-5; code frame X strike, Z up)
ForcedRuptureSpec Tpv26()
{
   ForcedRuptureSpec s;
   s.hypocenter[0] = -5000.0;
   s.hypocenter[1] = 0.0;
   s.hypocenter[2] = -10000.0;
   s.r_crit = 4000.0;
   s.vs = 3464.0;
   s.t0 = 0.5;
   return s;
}

/// The closure's SolveSlipRate as it was before forced rupture (copied).
real_t OldSolveSlipRate(const real_t X_tan[2], real_t gamma, real_t sigma_n,
                        const FaultPointState &s, real_t V_tan[2])
{
   const real_t mu = SlipWeakeningClosure::FrictionCoefficient(s.slip_abs, s.mu_s, s.mu_d, s.Dc);
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

} // namespace

void TestTime()
{
   std::cout << "\n=== U11: T(r) ===\n";
   const ForcedRuptureSpec s = Tpv26();
   const real_t v = 0.7 * 3464.0;
   const real_t at_hypo[3] = {-5000.0, 0.0, -10000.0};
   TEST_ASSERT(ForcedRupture::Time(s, at_hypo) == 0.0, "U11 T = 0 at the hypocenter");
   // r_c / 2 along strike and along dip: 2000/v + 0.081*4000/v*(1/(1-1/4) - 1)
   const real_t ref = 2000.0 / v + 0.081 * 4000.0 / v * (4.0 / 3.0 - 1.0);
   const real_t xs[3] = {-3000.0, 0.0, -10000.0}, xd[3] = {-5000.0, 0.0, -8000.0};
   TEST_ASSERT(std::abs(ForcedRupture::Time(s, xs) - ref) <= 1e-14 * ref,
               "U11 T(r_c/2) along strike = formula (" << ref << " s)");
   TEST_ASSERT(std::abs(ForcedRupture::Time(s, xd) - ref) <= 1e-14 * ref,
               "U11 T(r_c/2) along dip = formula");
   const real_t xc[3] = {-1000.0, 0.0, -10000.0}, xo[3] = {5000.0, 0.0, -10000.0};
   TEST_ASSERT(ForcedRupture::Time(s, xc) == 1.0e9 && ForcedRupture::Time(s, xo) == 1.0e9,
               "U11 T = 1.0e9 at r = r_c and beyond");
   // front speed dr/dT: 0.7 Vs at the hypocenter, falling to 0 at r_c
   const real_t h = 1e-3;
   const real_t v0 = h / ForcedRupture::TimeAtDistance(s, h);
   TEST_ASSERT(std::abs(v0 / 3464.0 - 0.7) < 1e-6, "U11 front speed 0.7 Vs at the hypocenter");
   bool mono = true;
   real_t prev = 0.0;
   for (int k = 1; k < 4000; k++)
   {
      const real_t T = ForcedRupture::TimeAtDistance(s, real_t(k));
      mono = mono && T > prev;
      prev = T;
   }
   TEST_ASSERT(mono, "U11 T increases with r inside r_c");
   const real_t r1 = 3999.0;
   const real_t vr = 1.0 / (ForcedRupture::TimeAtDistance(s, r1 + 0.5)
                            - ForcedRupture::TimeAtDistance(s, r1 - 0.5));
   TEST_ASSERT(vr / 3464.0 < 1e-3, "U11 front speed near r_c below 0.001 Vs (" << vr << " m/s)");

   ForcedRuptureSpec bad = s;
   bad.t0 = 0.0;
   TEST_ASSERT(ForcedRupture::Validate(bad) == "t0_s must be > 0", "U11 t0 = 0 refused");
   bad = s;
   bad.r_crit = -1.0;
   TEST_ASSERT(ForcedRupture::Validate(bad) == "r_crit_m must be > 0", "U11 r_crit < 0 refused");
   TEST_ASSERT(ForcedRupture::Validate(s).empty(), "U11 TPV26 spec accepted");
}

void TestCoefficient()
{
   std::cout << "\n=== U11: f2 and max(f1, f2) ===\n";
   const real_t ms = 0.18, md = 0.12, Dc = 0.30, T = 1.2, t0 = 0.5;
   auto mu = [&](real_t D, real_t t) { return ForcedRupture::FrictionCoefficient(D, ms, md, Dc, t, T, t0); };
   TEST_ASSERT(mu(0.0, 0.0) == ms && mu(0.0, T - 1e-9) == ms, "U11 mu_s before T");
   TEST_ASSERT(mu(0.0, T) == ms, "U11 f2 = 0 at t = T");
   TEST_ASSERT(std::abs(mu(0.0, T + 0.25) - 0.15) < 1e-15, "U11 f2 = 1/2 at T + t0/2");
   TEST_ASSERT(mu(0.0, T + t0) == md && mu(0.0, 100.0) == md, "U11 mu_d from T + t0 on");
   TEST_ASSERT(mu(Dc, 0.0) == md && mu(2.0 * Dc, 0.0) == md, "U11 slip alone weakens (f1 = 1)");
   TEST_ASSERT(std::abs(mu(0.6 * Dc, T + 0.25) - (ms + (md - ms) * 0.6)) < 1e-15,
               "U11 max(f1, f2) takes f1 when it is larger");
   TEST_ASSERT(std::abs(mu(0.2 * Dc, T + 0.25) - 0.15) < 1e-15,
               "U11 max(f1, f2) takes f2 when it is larger");
   TEST_ASSERT(ForcedRupture::FrictionCoefficient(0.0, ms, md, Dc, 1.0e8,
                                                  ForcedRupture::kNeverForced, t0) == ms,
               "U11 never forced (T = 1.0e9) at t = 1e8");
   const real_t ref = 0.18 + (0.12 - 0.18) * 0.5;
   TEST_ASSERT(ForcedRupture::FrictionCoefficient(0.0, ms, md, Dc, T + 0.25, T, t0) == ref,
               "U11 written as mu_s + (mu_d - mu_s) f (description Part 4)");
}

void TestSafsParity()
{
   std::cout << "\n=== U11: bitwise equal to SAFS LSWFrictionCoefficient_ForcedRupture ===\n";
   std::mt19937_64 g(20261020);
   std::uniform_real_distribution<real_t> u(0.0, 1.0);
   long n = 0, same = 0;
   for (int k = 0; k < 200000; k++)
   {
      const real_t ms = 0.1 + 0.8 * u(g), md = ms * u(g), Dc = 0.05 + u(g);
      const real_t D = 1.5 * Dc * u(g) * (k % 7 == 0 ? 0.0 : 1.0);
      const real_t T = (k % 5 == 0) ? 1.0e9 : 3.0 * u(g);
      const real_t t0 = 0.1 + u(g), t = 4.0 * u(g);
      const real_t a = ForcedRupture::FrictionCoefficient(D, ms, md, Dc, t, T, t0);
      const real_t b = spatial::LSWFrictionCoefficient_ForcedRupture(D, ms, md, Dc, t, T, t0);
      n++;
      same += (a == b);
   }
   TEST_ASSERT(same == n, "U11 " << same << " of " << n << " samples bitwise equal");
}

void TestClosure()
{
   std::cout << "\n=== U11: SlipWeakeningClosure ===\n";
   SlipWeakeningClosure cl;
   // no forced rupture: bitwise the arithmetic before it, whatever T and t hold
   std::mt19937_64 g(7);
   std::uniform_real_distribution<real_t> u(0.0, 1.0);
   long n = 0, same = 0;
   for (int k = 0; k < 100000; k++)
   {
      FaultPointState s;
      s.mu_s = 0.5 + u(g);
      s.mu_d = s.mu_s * u(g);
      s.Dc = 0.1 + u(g);
      s.cohesion = 1e6 * u(g);
      s.slip_abs = 1.5 * s.Dc * u(g);
      s.forced_time = (k % 2) ? 0.0 : 3.0 * u(g);
      s.forced_t0 = 0.0;
      s.t = 5.0 * u(g);
      const real_t X[2] = {(u(g) - 0.5) * 2e8, (u(g) - 0.5) * 2e8};
      const real_t gamma = 1e6 + 1e9 * u(g), sn = 1.2e8 * u(g);
      real_t Va[2], Vb[2];
      const real_t a = cl.SolveSlipRate(X, gamma, sn, s, Va);
      const real_t b = OldSolveSlipRate(X, gamma, sn, s, Vb);
      n++;
      same += (a == b && Va[0] == Vb[0] && Va[1] == Vb[1]);
   }
   TEST_ASSERT(same == n, "U11 forced_t0 = 0: " << same << " of " << n
               << " closure calls bitwise equal to the pre-forced-rupture closure");

   // forced: a stuck point is pushed to slip at T + t0
   FaultPointState s;
   s.mu_s = 0.18;
   s.mu_d = 0.12;
   s.Dc = 0.30;
   s.cohesion = 0.4e6;
   s.forced_time = 1.0;
   s.forced_t0 = 0.5;
   const real_t sn = 100e6, X[2] = {0.0, 15.0e6}, gamma = 5e6;   // between the two strengths
   real_t V[2];
   s.t = 0.9;
   TEST_ASSERT(cl.SolveSlipRate(X, gamma, sn, s, V) == 0.0, "U11 stuck before T");
   s.t = 1.5;
   const real_t Vabs = cl.SolveSlipRate(X, gamma, sn, s, V);
   TEST_ASSERT(std::abs(Vabs - (15.0e6 - 12.4e6) / gamma) < 1e-12 && V[1] < 0.0,
               "U11 slips at T + t0 with tau_s = C0 + mu_d sigma_n, anti-parallel");
   s.forced_t0 = -0.1;
   TEST_ASSERT(!cl.CheckParameters(s), "U11 forced_t0 < 0 refused by CheckParameters");
   s.forced_t0 = 0.5;
   TEST_ASSERT(cl.CheckParameters(s), "U11 TPV26 parameters accepted");
}

int main()
{
   std::cout << "Forced-rupture tests (rheology plan rev 1, U11)\n";
   TestTime();
   TestCoefficient();
   TestSafsParity();
   TestClosure();
   std::cout << "\n" << num_passed << " / " << num_tests << " passed, " << num_failed
             << " failed\n";
   return num_failed == 0 ? 0 : 1;
}
