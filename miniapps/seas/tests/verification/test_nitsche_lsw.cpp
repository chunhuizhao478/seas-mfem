// Nitsche dynamic rupture: the linear slip-weakening closure
// (friction/slip_weakening_closure.hpp; benchmark plan rev 1.3, section
// 4.3, phase B2).
//
//   L-01  stick branch: |X| <= tau_s gives V = 0 exactly, including
//         |X| = tau_s exactly.
//   L-02  slip branch: the end-of-step traction X + gamma V has magnitude
//         tau_s to 1e-14 relative and V . X < 0; mu(D) at D = 0, Dc/2, Dc,
//         2 Dc; cohesion adds to the strength; mu_s = 10000 never yields
//         while D = 0.  Also StateRate = |V|, SeedEquilibrium psi = 0 and
//         CheckParameters.
//   L-03  closure derivative: central finite differences of V(X) match the
//         closed form dV/dX in the slip branch (the T-04 analogue).
//   L-04  same rule as the ADER code: at 10^4 random inputs V equals
//         -(V of LSWFrictionCoefficient_TPV205 + SolveLSW_TPV205) with
//         eta_s = gamma and the same cohesion, to round-off (the ADER code
//         returns V parallel to the traction, the closure anti-parallel).
//
// Pointwise tests, no mesh.  Usage: ./seas_test_nitsche_lsw
//
// Scope (Chunhui Zhao, 2026-10-16).

#include "mfem.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../dynamic/tpv205_friction.hpp"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static FaultPointState TPV205Point(real_t D)
{
   FaultPointState s;
   s.mu_s = 0.677;
   s.mu_d = 0.525;
   s.Dc = 0.40;
   s.cohesion = 0.0;
   s.slip_abs = D;
   return s;
}

static real_t TauS(const FaultPointState &s, real_t sigma)
{
   return SlipWeakeningClosure::FrictionCoefficient(s.slip_abs, s.mu_s, s.mu_d, s.Dc)
          * sigma + s.cohesion;
}

// ------------------------------------------------------------------ L-01
static void TestStick()
{
   std::cout << "\n=== L-01 stick branch ===\n";
   SlipWeakeningClosure cl;
   std::mt19937_64 rng(20261016);
   std::uniform_real_distribution<real_t> U(0.0, 1.0);
   const real_t gamma = 1.7e7, sigma = 120e6;
   int bad = 0, n = 0;
   for (int k = 0; k < 10000; k++)
   {
      FaultPointState s = TPV205Point(0.8 * U(rng));
      s.cohesion = (k % 3 == 0) ? 0.0 : 2e6 * U(rng);
      const real_t ts = TauS(s, sigma);
      const real_t th = 2 * M_PI * U(rng), r = ts * U(rng);
      const real_t X[2] = {r * std::cos(th), r * std::sin(th)};
      real_t V[2] = {1.0, 1.0};
      const real_t Va = cl.SolveSlipRate(X, gamma, sigma, s, V);
      if (!(Va == 0.0 && V[0] == 0.0 && V[1] == 0.0)) { bad++; }
      n++;
   }
   // |X| = tau_s exactly: sqrt(fl(x^2)) == |x| in binary floating point
   int bad_eq = 0;
   const real_t Ds[4] = {0.0, 0.1, 0.4, 1.0};
   for (real_t D : Ds)
   {
      FaultPointState s = TPV205Point(D);
      const real_t ts = TauS(s, sigma);
      const real_t Xs[4][2] = {{ts, 0.0}, {0.0, ts}, {-ts, 0.0}, {0.0, -ts}};
      for (const auto &X : Xs)
      {
         real_t V[2] = {1.0, 1.0};
         const real_t Va = cl.SolveSlipRate(X, gamma, sigma, s, V);
         if (!(Va == 0.0 && V[0] == 0.0 && V[1] == 0.0)) { bad_eq++; }
      }
   }
   std::cout << "  " << n << " random points inside the disc: " << bad << " with V != 0; "
             << "16 points on the circle: " << bad_eq << " with V != 0\n";
   TEST_ASSERT(bad == 0, "L-01 |X| < tau_s gives V = 0 exactly (10^4 points)");
   TEST_ASSERT(bad_eq == 0, "L-01 |X| = tau_s gives V = 0 exactly");
}

// ------------------------------------------------------------------ L-02
static void TestSlip()
{
   std::cout << "\n=== L-02 slip branch, mu(D), cohesion, barrier ===\n";
   SlipWeakeningClosure cl;
   std::mt19937_64 rng(7);
   std::uniform_real_distribution<real_t> U(0.0, 1.0);
   const real_t sigma = 120e6;
   real_t worst = 0.0;
   int anti_bad = 0, n = 0;
   for (int k = 0; k < 10000; k++)
   {
      FaultPointState s = TPV205Point(0.8 * U(rng));
      s.cohesion = (k % 3 == 0) ? 0.0 : 2e6 * U(rng);
      const real_t gamma = std::pow(10.0, 6.0 + 2.0 * U(rng));
      const real_t ts = TauS(s, sigma);
      const real_t th = 2 * M_PI * U(rng), r = ts * (1.0 + 2.0 * U(rng)) + 1.0;
      const real_t X[2] = {r * std::cos(th), r * std::sin(th)};
      real_t V[2];
      const real_t Va = cl.SolveSlipRate(X, gamma, sigma, s, V);
      const real_t E0 = X[0] + gamma * V[0], E1 = X[1] + gamma * V[1];
      worst = std::max(worst, std::abs(std::hypot(E0, E1) - ts) / ts);
      if (!(Va > 0.0 && X[0] * V[0] + X[1] * V[1] < 0.0
            && std::abs(std::hypot(V[0], V[1]) - Va) <= 1e-15 * Va)) { anti_bad++; }
      n++;
   }
   std::cout << "  " << n << " random points outside the disc: max | |X + gamma V| - tau_s | / tau_s = "
             << worst << "; " << anti_bad << " not anti-parallel or |V| mismatch\n";
   TEST_ASSERT(worst < 1e-14, "L-02 |X + gamma V| = tau_s to 1e-14 relative");
   TEST_ASSERT(anti_bad == 0, "L-02 V anti-parallel to X (V . X < 0), returned |V| = |V_tan|");

   // mu(D) at D = 0, Dc/2, Dc, 2 Dc
   const real_t ms = 0.677, md = 0.525, Dc = 0.40;
   const real_t m0 = SlipWeakeningClosure::FrictionCoefficient(0.0, ms, md, Dc);
   const real_t mh = SlipWeakeningClosure::FrictionCoefficient(0.5 * Dc, ms, md, Dc);
   const real_t m1 = SlipWeakeningClosure::FrictionCoefficient(Dc, ms, md, Dc);
   const real_t m2 = SlipWeakeningClosure::FrictionCoefficient(2.0 * Dc, ms, md, Dc);
   std::cout << std::setprecision(17) << "  mu(0) = " << m0 << ", mu(Dc/2) = " << mh
             << ", mu(Dc) = " << m1 << ", mu(2 Dc) = " << m2 << "\n";
   TEST_ASSERT(m0 == ms && m1 == md && m2 == md
               && std::abs(mh - 0.5 * (ms + md)) <= 1e-15 * ms,
               "L-02 mu(0) = mu_s, mu(Dc/2) = (mu_s + mu_d)/2, mu(Dc) = mu(2 Dc) = mu_d");

   // cohesion adds: same X, c = 0 and c = 1 MPa, both slipping
   {
      const real_t gamma = 1.7e7, c = 1e6;
      FaultPointState s0 = TPV205Point(0.1), s1 = TPV205Point(0.1);
      s1.cohesion = c;
      const real_t X[2] = {-90e6, 20e6};
      real_t V0[2], V1[2];
      const real_t Va0 = cl.SolveSlipRate(X, gamma, sigma, s0, V0);
      const real_t Va1 = cl.SolveSlipRate(X, gamma, sigma, s1, V1);
      const real_t rel = std::abs((Va0 - Va1) - c / gamma) / (c / gamma);
      std::cout << std::setprecision(6) << "  cohesion 1 MPa: |V| " << Va0 << " -> " << Va1
                << " m/s, (|V0| - |V1|) / (c / gamma) - 1 = " << rel << "\n";
      TEST_ASSERT(Va1 > 0.0 && rel < 1e-9, "L-02 cohesion c lowers |V| by c / gamma");
   }

   // mu_s = 10000 never yields while D = 0
   {
      int yields = 0;
      for (int k = 0; k < 1000; k++)
      {
         FaultPointState s = TPV205Point(0.0);
         s.mu_s = 10000.0;
         const real_t th = 2 * M_PI * U(rng), r = 1e9 * U(rng);   // up to 1 GPa
         const real_t X[2] = {r * std::cos(th), r * std::sin(th)};
         real_t V[2];
         if (cl.SolveSlipRate(X, 1e6, 1e6 + 200e6 * U(rng), s, V) != 0.0) { yields++; }
      }
      TEST_ASSERT(yields == 0, "L-02 mu_s = 10000, D = 0: no yield up to |X| = 1 GPa");
   }

   // StateRate, SeedEquilibrium, CheckParameters
   {
      FaultPointState s = TPV205Point(0.3);
      s.psi = 5.0;
      cl.SeedEquilibrium(70e6, 1e-12, 120e6, s);
      const bool ok_rate = cl.StateRate(2.5, s) == 2.5 && cl.StateRate(0.0, s) == 0.0;
      FaultPointState g = TPV205Point(0.0), b1 = g, b2 = g, b3 = g, b4 = g;
      b1.mu_d = 0.7;      // mu_d > mu_s
      b2.mu_d = -0.1;     // mu_d < 0
      b3.Dc = 0.0;        // Dc = 0
      b4.cohesion = -1.0; // c < 0
      const bool ok_chk = cl.CheckParameters(g) && !cl.CheckParameters(b1)
                          && !cl.CheckParameters(b2) && !cl.CheckParameters(b3)
                          && !cl.CheckParameters(b4);
      TEST_ASSERT(s.psi == 0.0, "L-02 SeedEquilibrium: slip path starts at 0");
      TEST_ASSERT(ok_rate, "L-02 StateRate = |V| (slip path)");
      TEST_ASSERT(ok_chk, "L-02 CheckParameters: mu_s >= mu_d >= 0, Dc > 0, c >= 0");
   }
}

// ------------------------------------------------------------------ L-03
static void TestDerivative()
{
   std::cout << "\n=== L-03 dV/dX, slip branch ===\n";
   SlipWeakeningClosure cl;
   std::mt19937_64 rng(11);
   std::uniform_real_distribution<real_t> U(0.0, 1.0);
   const real_t sigma = 120e6, gamma = 1.7e7;
   real_t worst = 0.0;
   for (int k = 0; k < 200; k++)
   {
      FaultPointState s = TPV205Point(0.8 * U(rng));
      const real_t ts = TauS(s, sigma);
      const real_t th = 2 * M_PI * U(rng), r = ts * (1.05 + U(rng));
      const real_t X[2] = {r * std::cos(th), r * std::sin(th)};
      const real_t Xa = std::hypot(X[0], X[1]), xh[2] = {X[0] / Xa, X[1] / Xa};
      const real_t q = (Xa - ts) / Xa;
      real_t J[2][2], Jfd[2][2];
      for (int i = 0; i < 2; i++)
         for (int j = 0; j < 2; j++)
         {
            J[i][j] = -(xh[i] * xh[j] + q * ((i == j ? 1.0 : 0.0) - xh[i] * xh[j])) / gamma;
         }
      const real_t h = 1e-5 * Xa;
      for (int j = 0; j < 2; j++)
      {
         real_t Xp[2] = {X[0], X[1]}, Xm[2] = {X[0], X[1]}, Vp[2], Vm[2];
         Xp[j] += h;
         Xm[j] -= h;
         cl.SolveSlipRate(Xp, gamma, sigma, s, Vp);
         cl.SolveSlipRate(Xm, gamma, sigma, s, Vm);
         for (int i = 0; i < 2; i++) { Jfd[i][j] = (Vp[i] - Vm[i]) / (2 * h); }
      }
      for (int i = 0; i < 2; i++)
         for (int j = 0; j < 2; j++)
         {
            worst = std::max(worst, std::abs(Jfd[i][j] - J[i][j]) * gamma);
         }
   }
   std::cout << "  200 points: max |J_fd - J| gamma = " << worst << "\n";
   TEST_ASSERT(worst < 1e-6, "L-03 finite-difference dV/dX matches the closed form (1e-6 of 1/gamma)");
}

// ------------------------------------------------------------------ L-04
static void TestADER()
{
   std::cout << "\n=== L-04 same rule as the ADER code (TPV205) ===\n";
   SlipWeakeningClosure cl;
   std::mt19937_64 rng(205);
   std::uniform_real_distribution<real_t> U(0.0, 1.0);
   int n = 0, n_slip = 0, n_bits = 0;
   real_t worst = 0.0;
   for (int k = 0; k < 10000; k++)
   {
      FaultPointState s;
      s.mu_s = 0.5 + 0.5 * U(rng);
      s.mu_d = 0.3 + (s.mu_s - 0.3) * U(rng);
      s.Dc = 0.1 + 0.9 * U(rng);
      s.cohesion = (k % 2 == 0) ? 0.0 : 3e6 * U(rng);
      const int dk = k % 5;   // D = 0, inside, Dc, beyond, random
      s.slip_abs = (dk == 0) ? 0.0 : (dk == 1) ? 0.5 * s.Dc * U(rng)
                   : (dk == 2) ? s.Dc : (dk == 3) ? s.Dc * (1.0 + U(rng)) : 2.0 * s.Dc * U(rng);
      const real_t sigma = 200e6 * U(rng);
      const real_t gamma = std::pow(10.0, 6.0 + 2.0 * U(rng));
      const real_t X[2] = {150e6 * (2 * U(rng) - 1), 150e6 * (2 * U(rng) - 1)};
      real_t V[2];
      const real_t Va = cl.SolveSlipRate(X, gamma, sigma, s, V);
      const real_t mu = LSWFrictionCoefficient_TPV205(s.slip_abs, s.mu_s, s.mu_d, s.Dc);
      real_t Va_r, V1_r, V2_r, t1c, t2c;
      SolveLSW_TPV205(0.0, 0.0, X[0], X[1], sigma, gamma, mu, Va_r, V1_r, V2_r, t1c, t2c,
                      0.0, s.cohesion);
      const real_t Vr[2] = {-V1_r, -V2_r};
      const real_t sc = std::max(Va_r, 1e-300);
      worst = std::max({worst, std::abs(Va - Va_r) / sc, std::abs(V[0] - Vr[0]) / sc,
                        std::abs(V[1] - Vr[1]) / sc});
      if (std::memcmp(&Va, &Va_r, sizeof(real_t)) == 0 && V[0] == Vr[0] && V[1] == Vr[1])
      {
         n_bits++;
      }
      if (Va_r > 0.0) { n_slip++; }
      n++;
   }
   std::cout << "  " << n << " inputs (" << n_slip << " slipping): max relative difference "
             << worst << "; " << n_bits << " bit-identical\n";
   TEST_ASSERT(n_slip > 1000 && n - n_slip > 1000, "L-04 inputs cover both branches");
   TEST_ASSERT(worst <= 1e-15, "L-04 V == -V_ADER (eta_s = gamma) to round-off");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   std::cout << "=== Linear slip-weakening closure (L-01 to L-04) ===\n";
   TestStick();
   TestSlip();
   TestDerivative();
   TestADER();
   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, "
             << num_failed << " failed ===\n";
   return num_failed == 0 ? 0 : 1;
}
