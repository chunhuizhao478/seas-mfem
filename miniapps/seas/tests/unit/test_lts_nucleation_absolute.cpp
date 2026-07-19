// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_nucleation_absolute.cpp — Appendix B.10 of
// PLAN_clustered_lts_ader_2026-07-18.md (LTS Phase 3, D-3).
//
// The LTS fault sweep applies each cluster's nucleation forcing over its own
// cluster-contiguous global-QP range, at its own stage times, via the ABSOLUTE
// appliers' new (qp_begin, qp_end) range overloads.  This test pins the three
// D-3 properties for BOTH gradual kinds:
//
//   (1) telescoping: the absolute value at any stage time t equals the
//       telescoped incremental sum over any partition of [0, t] to 1e-15·amp;
//   (2) partition independence: two different time partitions produce identical
//       tau_nuc at a common stage time (the absolute form is time-only);
//   (3) range restriction: applying disjoint ranges [0,k) then [k,N) equals a
//       single whole-vector apply, BIT-FOR-BIT, and indices outside a range are
//       left byte-untouched.

#include "../../dynamic/spatial_nucleation.hpp"
#include "../../dynamic/nucleation_method.hpp"   // INucleationMethod range routing (D-3)

#include <cmath>
#include <cstdio>
#include <vector>

using namespace mfem::seas;
using namespace mfem::seas::spatial;
using mfem::real_t;
using mfem::Vector;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                     \
   do { ++g_checks; if (!(cond)) { ++g_fails;                               \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

namespace
{
// A deterministic, non-trivial per-DOF amplitude field.
GradualOverstressPerDOFParams make_gradual(int n)
{
   GradualOverstressPerDOFParams p;
   p.amplitude_dip.SetSize(n);
   p.amplitude_strike.SetSize(n);
   p.radial.SetSize(n);
   for (int i = 0; i < n; ++i)
   {
      p.amplitude_dip(i)    = 1.0e6 * (0.3 + 0.11 * i);
      p.amplitude_strike(i) = 2.0e6 * (0.7 - 0.05 * i);
      p.radial(i)           = 1.0;
   }
   return p;
}

CompactCircularPerDOFParams make_compact(int n)
{
   CompactCircularPerDOFParams p;
   p.amplitude_strike.SetSize(n);
   p.radial.SetSize(n);
   for (int i = 0; i < n; ++i)
   {
      p.amplitude_strike(i) = 1.3e6 * (0.4 + 0.09 * i);
      p.radial(i)           = 1.0;
   }
   return p;
}

std::vector<DOFData> zeros(int n) { return std::vector<DOFData>(n); }

real_t max_amp(const GradualOverstressPerDOFParams& p)
{
   real_t m = 0.0;
   for (int i = 0; i < p.amplitude_dip.Size(); ++i)
   {
      m = std::max(m, std::abs(p.amplitude_dip(i)));
      m = std::max(m, std::abs(p.amplitude_strike(i)));
   }
   return m;
}
}  // namespace

int main()
{
   const int n = 7;
   const real_t T_nuc = 1.5;

   // ---------------------------------------------------------------------
   // T1 — gradual: absolute at t == telescoped increment sum over [0, t].
   // ---------------------------------------------------------------------
   {
      const GradualOverstressPerDOFParams p = make_gradual(n);
      const real_t amp = max_amp(p);

      // Telescoped path: a fine, IRREGULAR partition of [0, T_nuc + margin].
      const std::vector<real_t> dts = {0.13, 0.27, 0.05, 0.4, 0.31, 0.22, 0.19,
                                       0.3, 0.5};
      std::vector<DOFData> tel = zeros(n);
      real_t t = 0.0;
      for (real_t dt : dts)
      {
         t += dt;
         ApplyGradualOverstressIncrement(tel, p, T_nuc, t, dt);

         // Absolute path evaluated at the SAME stage time.
         std::vector<DOFData> abs = zeros(n);
         ApplyGradualOverstressAbsolute(abs, p, T_nuc, t);

         real_t md = 0.0;
         for (int i = 0; i < n; ++i)
         {
            md = std::max(md, std::abs(tel[i].tau1_nuc - abs[i].tau1_nuc));
            md = std::max(md, std::abs(tel[i].tau2_nuc - abs[i].tau2_nuc));
         }
         CHECK(md <= 1e-15 * amp,
               "gradual: absolute != telescoped increment sum at stage time");
      }
   }

   // ---------------------------------------------------------------------
   // T2 — compact-circular: same telescoping property.
   // ---------------------------------------------------------------------
   {
      const CompactCircularPerDOFParams p = make_compact(n);
      real_t amp = 0.0;
      for (int i = 0; i < n; ++i)
      { amp = std::max(amp, std::abs(p.amplitude_strike(i))); }

      const std::vector<real_t> dts = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
      std::vector<DOFData> tel = zeros(n);
      real_t t = 0.0;
      for (real_t dt : dts)
      {
         t += dt;
         ApplyGradualOverstressCompactCircularIncrement(tel, p, T_nuc, t, dt);
         std::vector<DOFData> abs = zeros(n);
         ApplyGradualOverstressCompactCircularAbsolute(abs, p, T_nuc, t);
         real_t md = 0.0;
         for (int i = 0; i < n; ++i)
         { md = std::max(md, std::abs(tel[i].tau2_nuc - abs[i].tau2_nuc)); }
         CHECK(md <= 1e-15 * amp,
               "compact: absolute != telescoped increment sum at stage time");
      }
   }

   // ---------------------------------------------------------------------
   // T3 — partition independence: absolute form is time-only, so evaluating
   // at a common time t* from two different partitions gives identical tau_nuc.
   // ---------------------------------------------------------------------
   {
      const GradualOverstressPerDOFParams p = make_gradual(n);
      const real_t tstar = 0.9;  // mid-ramp, so S in (0,1)
      std::vector<DOFData> a = zeros(n), b = zeros(n);
      // Partition A reached tstar via one apply; partition B via prior applies
      // at other times then a final apply at tstar (absolute overwrites).
      ApplyGradualOverstressAbsolute(a, p, T_nuc, tstar);
      ApplyGradualOverstressAbsolute(b, p, T_nuc, 0.3);
      ApplyGradualOverstressAbsolute(b, p, T_nuc, 1.2);
      ApplyGradualOverstressAbsolute(b, p, T_nuc, tstar);
      for (int i = 0; i < n; ++i)
      {
         CHECK(a[i].tau1_nuc == b[i].tau1_nuc && a[i].tau2_nuc == b[i].tau2_nuc,
               "gradual: partition-dependent tau_nuc at common stage time");
      }
   }

   // ---------------------------------------------------------------------
   // T4 — range restriction (gradual): disjoint ranges == whole-vector, and
   // out-of-range indices are byte-untouched.  Split at k.
   // ---------------------------------------------------------------------
   {
      const GradualOverstressPerDOFParams p = make_gradual(n);
      const real_t tstar = 0.77;
      const std::size_t k = 3;

      std::vector<DOFData> whole = zeros(n);
      ApplyGradualOverstressAbsolute(whole, p, T_nuc, tstar);

      std::vector<DOFData> ranged = zeros(n);
      ApplyGradualOverstressAbsolute(ranged, p, T_nuc, tstar, 0, k);
      // Sentinel: after the first range, [k, n) must still be exactly zero.
      for (std::size_t i = k; i < static_cast<std::size_t>(n); ++i)
      {
         CHECK(ranged[i].tau1_nuc == 0.0 && ranged[i].tau2_nuc == 0.0,
               "gradual range: out-of-range index was written");
      }
      ApplyGradualOverstressAbsolute(ranged, p, T_nuc, tstar, k, n);

      for (int i = 0; i < n; ++i)
      {
         CHECK(ranged[i].tau1_nuc == whole[i].tau1_nuc
               && ranged[i].tau2_nuc == whole[i].tau2_nuc,
               "gradual range: disjoint-range apply != whole-vector apply");
      }
   }

   // ---------------------------------------------------------------------
   // T5 — range restriction (compact-circular): same bit-for-bit property.
   // ---------------------------------------------------------------------
   {
      const CompactCircularPerDOFParams p = make_compact(n);
      const real_t tstar = 0.61;
      const std::size_t k = 5;

      std::vector<DOFData> whole = zeros(n);
      ApplyGradualOverstressCompactCircularAbsolute(whole, p, T_nuc, tstar);

      std::vector<DOFData> ranged = zeros(n);
      ApplyGradualOverstressCompactCircularAbsolute(ranged, p, T_nuc, tstar, 0, k);
      for (std::size_t i = k; i < static_cast<std::size_t>(n); ++i)
      {
         CHECK(ranged[i].tau2_nuc == 0.0,
               "compact range: out-of-range index was written");
      }
      ApplyGradualOverstressCompactCircularAbsolute(ranged, p, T_nuc, tstar, k, n);
      for (int i = 0; i < n; ++i)
      {
         CHECK(ranged[i].tau2_nuc == whole[i].tau2_nuc,
               "compact range: disjoint-range apply != whole-vector apply");
      }
   }

   // ---------------------------------------------------------------------
   // T6 — empty range is a no-op; a single mid-vector range touches only it.
   // ---------------------------------------------------------------------
   {
      const GradualOverstressPerDOFParams p = make_gradual(n);
      std::vector<DOFData> d = zeros(n);
      ApplyGradualOverstressAbsolute(d, p, T_nuc, 0.5, 2, 2);  // empty
      for (int i = 0; i < n; ++i)
      {
         CHECK(d[i].tau1_nuc == 0.0 && d[i].tau2_nuc == 0.0,
               "gradual range: empty range wrote something");
      }
      ApplyGradualOverstressAbsolute(d, p, T_nuc, 0.5, 2, 4);  // [2,4)
      const real_t S = SmoothStep(0.5, T_nuc);
      for (int i = 0; i < n; ++i)
      {
         const bool in = (i >= 2 && i < 4);
         const real_t exp1 = in ? S * p.amplitude_dip(i)    : 0.0;
         const real_t exp2 = in ? S * p.amplitude_strike(i) : 0.0;
         CHECK(d[i].tau1_nuc == exp1 && d[i].tau2_nuc == exp2,
               "gradual range: mid-vector range wrote the wrong indices");
      }
   }

   // ---------------------------------------------------------------------
   // T7 (R-002) — D-3 contract: the range ApplyAbsolute routes through
   // INucleationMethod.  The gradual strategy objects' range ApplyAbsolute must
   // equal the free-function range result bit-for-bit, and the base-class
   // full-range default must delegate to the whole-vector form.
   // ---------------------------------------------------------------------
   {
      const real_t tstar = 0.83;
      const std::size_t b = 2, e = 5;

      // Gaussian gradual.
      {
         const GradualOverstressPerDOFParams p = make_gradual(n);
         GaussianGradualOverstress method(p, T_nuc);   // copies params
         std::vector<DOFData> via_iface = zeros(n), via_free = zeros(n);
         method.ApplyAbsolute(via_iface, tstar, b, e);
         ApplyGradualOverstressAbsolute(via_free, p, T_nuc, tstar, b, e);
         real_t md = 0.0;
         for (int i = 0; i < n; ++i)
         {
            md = std::max(md, std::abs(via_iface[i].tau1_nuc - via_free[i].tau1_nuc));
            md = std::max(md, std::abs(via_iface[i].tau2_nuc - via_free[i].tau2_nuc));
         }
         CHECK(md == 0.0, "gaussian: INucleationMethod range == free-function range");

         // Base full-range default delegates to the whole-vector form.
         std::vector<DOFData> full_range = zeros(n), whole = zeros(n);
         method.ApplyAbsolute(full_range, tstar, 0, n);
         method.ApplyAbsolute(whole, tstar);
         real_t md2 = 0.0;
         for (int i = 0; i < n; ++i)
         {
            md2 = std::max(md2, std::abs(full_range[i].tau1_nuc - whole[i].tau1_nuc));
            md2 = std::max(md2, std::abs(full_range[i].tau2_nuc - whole[i].tau2_nuc));
         }
         CHECK(md2 == 0.0, "gaussian: range (0,N) == whole-vector ApplyAbsolute");
      }

      // Compact-circular gradual.
      {
         const CompactCircularPerDOFParams p = make_compact(n);
         CompactCircularGradualOverstress method(p, T_nuc);
         std::vector<DOFData> via_iface = zeros(n), via_free = zeros(n);
         method.ApplyAbsolute(via_iface, tstar, b, e);
         ApplyGradualOverstressCompactCircularAbsolute(via_free, p, T_nuc, tstar, b, e);
         real_t md = 0.0;
         for (int i = 0; i < n; ++i)
         { md = std::max(md, std::abs(via_iface[i].tau2_nuc - via_free[i].tau2_nuc)); }
         CHECK(md == 0.0, "compact: INucleationMethod range == free-function range");
      }
   }

   std::printf("test_lts_nucleation_absolute: %d checks, %d failures\n",
               g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}
