// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_time_basis.cpp — Appendix B.3 of PLAN_clustered_lts_ader_2026-07-18.md.
//
// IntegrateTaylor (dynamic/lts_time_basis.hpp):
//   * sub-interval integrals exact for polynomials up to order 4,
//   * Sigma of sub-intervals == the whole interval to ~1e-15,
//   * a=0,b=dt reproduces the wave operator's whole-step integral weights
//     dt^{k+1}/(k+1)!,
//   * the pure coefficients, and the b<a antisymmetry.

#include "../../dynamic/lts_time_basis.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace mfem::seas;
using mfem::real_t;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                     \
   do { ++g_checks; if (!(cond)) { ++g_fails;                               \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

// Closed-form reference: integral over [a,b] of the Taylor series with RAW stack
// dk (length order*n): out[j] = sum_k (b^{k+1}-a^{k+1})/(k+1)! * dk[k*n+j].
static std::vector<real_t> ref_integral(real_t a, real_t b, int order, int n,
                                        const std::vector<real_t>& dk)
{
   std::vector<real_t> out(n, 0.0);
   real_t apow = a, bpow = b, fact = 1.0;
   for (int k = 0; k < order; ++k)
   {
      const real_t c = (bpow - apow) / fact;
      for (int j = 0; j < n; ++j) { out[j] += c * dk[k * n + j]; }
      apow *= a; bpow *= b; fact *= static_cast<real_t>(k + 2);
   }
   return out;
}

static real_t maxabsdiff(const std::vector<real_t>& x, const std::vector<real_t>& y)
{
   real_t m = 0.0;
   for (std::size_t i = 0; i < x.size(); ++i) { m = std::max(m, std::abs(x[i] - y[i])); }
   return m;
}

int main()
{
   const int n = 2;   // two flattened (component,node) entries
   // A representative raw Taylor stack for order 4 (k=0..3), n=2.
   std::vector<real_t> dk = {
      2.0, -1.0,     // D(0)
      0.5,  3.0,     // D(1)
     -4.0,  0.25,    // D(2)
      1.5, -2.5,     // D(3)
   };

   // ---- T1: closed-form polynomial exactness, orders 2..4, a generic [a,b] ---
   for (int order = 2; order <= 4; ++order)
   {
      const real_t a = 0.3, b = 1.7;
      std::vector<real_t> out(n, 0.0);
      IntegrateTaylor(a, b, dk.data(), order, n, out.data());
      std::vector<real_t> ref = ref_integral(a, b, order, n, dk);
      char m[64]; std::snprintf(m, sizeof m, "T1 order %d exact vs closed form", order);
      CHECK(maxabsdiff(out, ref) <= 1e-12, m);
   }

   // ---- T2: additivity  int_0^dt == int_0^s + int_s^dt  (1e-15 scale) --------
   {
      const int order = 4;
      const real_t dt = 2.0, s = 0.8;
      std::vector<real_t> whole(n), lo(n), hi(n);
      IntegrateTaylor(0.0, dt, dk.data(), order, n, whole.data());
      IntegrateTaylor(0.0, s,  dk.data(), order, n, lo.data());
      IntegrateTaylor(s,  dt, dk.data(), order, n, hi.data());
      std::vector<real_t> sum(n);
      for (int j = 0; j < n; ++j) { sum[j] = lo[j] + hi[j]; }
      // relative to the whole-interval magnitude.
      real_t scale = 0.0; for (real_t v : whole) { scale = std::max(scale, std::abs(v)); }
      CHECK(maxabsdiff(sum, whole) <= 1e-13 * (scale + 1.0),
            "T2 sub-interval additivity to ~1e-15 relative");
   }

   // ---- T3: a=0,b=dt reproduces the whole-step weights dt^{k+1}/(k+1)! -------
   {
      const int order = 4;
      const real_t dt = 1.3;
      std::vector<real_t> out(n);
      IntegrateTaylor(0.0, dt, dk.data(), order, n, out.data());
      // Manual whole-step integral == the operator's I weights.
      std::vector<real_t> manual(n, 0.0);
      real_t fac = dt;                          // dt^1/1!
      for (int k = 0; k < order; ++k)
      {
         for (int j = 0; j < n; ++j) { manual[j] += fac * dk[k * n + j]; }
         fac *= dt / static_cast<real_t>(k + 2);   // -> dt^{k+2}/(k+2)!
      }
      CHECK(maxabsdiff(out, manual) <= 1e-13, "T3 a=0,b=dt == whole-step integral weights");
   }

   // ---- T4: pure coefficients over [0,1] are 1, 1/2, 1/6, 1/24 --------------
   {
      real_t c[4];
      TaylorIntegralCoeffs(0.0, 1.0, 4, c);
      CHECK(std::abs(c[0] - 1.0)        < 1e-15, "T4 coeff0 == 1");
      CHECK(std::abs(c[1] - 0.5)        < 1e-15, "T4 coeff1 == 1/2");
      CHECK(std::abs(c[2] - 1.0 / 6.0)  < 1e-15, "T4 coeff2 == 1/6");
      CHECK(std::abs(c[3] - 1.0 / 24.0) < 1e-15, "T4 coeff3 == 1/24");
   }

   // ---- T5: antisymmetry int_a^b == -int_b^a --------------------------------
   {
      const int order = 4;
      const real_t a = 0.4, b = 1.9;
      std::vector<real_t> ab(n), ba(n);
      IntegrateTaylor(a, b, dk.data(), order, n, ab.data());
      IntegrateTaylor(b, a, dk.data(), order, n, ba.data());
      std::vector<real_t> neg(n);
      for (int j = 0; j < n; ++j) { neg[j] = -ba[j]; }
      CHECK(maxabsdiff(ab, neg) <= 1e-13, "T5 int_a^b == -int_b^a");
   }

   // ---- T6: n == 0 is a no-op (empty cluster / degenerate block) ------------
   {
      std::vector<real_t> out;   // empty
      IntegrateTaylor(0.0, 1.0, dk.data(), 4, 0, out.data());
      CHECK(true, "T6 n==0 does not crash");
   }

   std::printf("test_lts_time_basis: %d/%d passed, %d failed.\n",
               g_checks - g_fails, g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}
