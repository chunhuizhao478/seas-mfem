// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_thermal_pressurization.cpp — 2D in-plane design rev
// 3, phase D (part 1): the thermal pressurization module
// fault/thermal_pressurization.hpp alone, against the solutions of the 1-D
// problem for a shear heating tau V that is constant in time (section 5.9,
// derivation 10.5).
//
// Analytic solutions (whole line, Gaussian source of standard deviation h
// and total power q = tau V per unit fault area, from t = 0):
//   theta(z, t; alpha) = q / (rho_c alpha sqrt(2 pi)) [F(s_t, z) - F(h, z)],
//      s_t = sqrt(h^2 + 2 alpha t),
//      F(s, z) = s exp(-z^2 / (2 s^2)) - |z| sqrt(pi / 2) erfc(|z| / (sqrt(2) s))
//   (the time integral of the spreading Gaussian; at z = 0 and h -> 0 it is
//   Rice's (2006) q / rho_c sqrt(t / (pi alpha)));
//   pi(z, t) = Lambda / (alpha_hy - alpha_th)
//              [alpha_hy theta(z, t; alpha_hy) - alpha_th theta(z, t; alpha_th)]
//   (Noda and Lapusta 2010 for a constant source; the undrained limit
//   alpha_hy -> 0 gives pi = Lambda theta).  A source switched on at t_1 and
//   off at t_2 is the difference of two of these.
//
//   P1  the grid: z_0 = 0, spacings dz0 r^k, z_{M-1} < z_max <= z_M, the
//       control volumes sum to z_{M-1/2}, the source weights to
//       erf(z_{M-1/2} / (sqrt(2) h)) / 2
//   P2  zero heating leaves T and p exactly at T_ini and p_ini; Lambda = 0
//       leaves p exactly at p_ini whatever the heating (the Lambda = 0
//       consistency run of gate G2 rests on it)
//   P3  heat balance: rho_c times the integral of T - T_ini equals the heat
//       put into the grid, q t times the source fraction, to 1e-12
//   P4  the analytic formulas: h -> 0 gives Rice's sqrt law; alpha_hy -> 0
//       gives pi = Lambda theta
//   P5  TPV105 parameters, constant tau V = 30 MPa x 1 m/s, on a fine grid and
//       step (Crank-Nicolson): T and p on the fault and at the nodes nearest
//       z = h, 2h, 5h, against the formulas to 1e-4 of the value on the fault,
//       at t = 0.01 .. 15 s; the same with the source switched off at 2 s
//   P6  convergence: halving the spacings (dz0 / 2, growth sqrt(r)) divides
//       the error by about 4; in time, the differences of runs with dt, dt/2
//       and dt/4 shrink by about 2 with backward Euler and 4 with
//       Crank-Nicolson
//   P7  production settings (dz0 = h / 10 and h / 20, growth 1.05, z_max
//       0.4 m; the 2D 100 m steps 7.2e-4 s at p = 1 and 3.1e-4 s at p = 2):
//       the largest error of T and p on the fault over 0 < t <= 15 s, both
//       schemes, a constant source and a 1 s pulse; printed, and held below
//       1e-3 at dz0 = h / 10
//   P8  points are independent (three points with different parameters and
//       heating equal three one-point runs bit for bit); GetState / SetState
//       continue a run bit for bit
//
// Modes: (default) P1-P8;  --refuse N (N = 1..7) must abort with:
//    1 dt = 0                                   "ThermalPressurization::Advance: dt = 0"
//    2 tau V < 0                                "the shear heating must be finite and >= 0"
//    3 V = NaN                                  "the shear heating must be finite and >= 0"
//    4 growth 0.9                               "need dz0 > 0, growth >= 1 and z_max > dz0"
//    5 alpha_th = 0                             "need alpha_th, alpha_hy, rho_c, h > 0"
//    6 Advance before SetParameters             "has no parameters (SetParameters)"
//    7 SetState of the wrong size               "ThermalPressurization::SetState: size"
//
//   make -f tests/2d_inplane/2d_inplane.mk test-thermal-pressurization
//
// Scope (Chunhui Zhao, 2026-10-26).

#include "mfem.hpp"
#include "../../fault/thermal_pressurization.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

// TPV105 (description, thermal pressurization table and initial conditions)
static TPParameters TPV105()
{
   TPParameters P;
   P.alpha_th = 1.0e-6;
   P.alpha_hy = 1.0e-4;
   P.rho_c = 2.7e6;
   P.Lambda = 0.1e6;
   P.h = 0.01;
   P.T_ini = 483.15;
   P.p_ini = 80.0e6;
   return P;
}

// ---------------------------------------------------------------- analytic
static real_t F(real_t s, real_t z)
{
   const real_t a = std::fabs(z);
   return s * std::exp(-z * z / (2.0 * s * s))
          - a * std::sqrt(M_PI / 2.0) * std::erfc(a / (std::sqrt(2.0) * s));
}

/// T - T_ini at (z, t) for a constant source q from t = 0 (0 for t <= 0).
static real_t Theta(real_t z, real_t t, real_t alpha, const TPParameters &P, real_t q)
{
   if (t <= 0.0) { return 0.0; }
   return q / (P.rho_c * alpha * std::sqrt(2.0 * M_PI))
          * (F(std::sqrt(P.h * P.h + 2.0 * alpha * t), z) - F(P.h, z));
}

/// p - p_ini at (z, t) for a constant source q from t = 0.
static real_t Pi(real_t z, real_t t, const TPParameters &P, real_t q)
{
   return P.Lambda / (P.alpha_hy - P.alpha_th)
          * (P.alpha_hy * Theta(z, t, P.alpha_hy, P, q) - P.alpha_th * Theta(z, t, P.alpha_th, P, q));
}

/// q on over [t1, t2): the difference of two constant-source solutions.
struct Pulse
{
   real_t q, t1, t2;
   real_t Q(real_t t) const { return (t >= t1 && t < t2) ? q : 0.0; }
   real_t T(real_t z, real_t t, const TPParameters &P) const
   { return Theta(z, t - t1, P.alpha_th, P, q) - Theta(z, t - t2, P.alpha_th, P, q); }
   real_t p(real_t z, real_t t, const TPParameters &P) const
   { return Pi(z, t - t1, P, q) - Pi(z, t - t2, P, q); }
};

/// The node nearest to z (the comparison is made at the node's own z, so no
/// interpolation error enters it).
static int NearestNode(const ThermalPressurization &tp, real_t z)
{
   const std::vector<real_t> &Z = tp.Z();
   int k = 0;
   for (int i = 1; i < tp.NumZ(); i++) { if (std::fabs(Z[i] - z) < std::fabs(Z[k] - z)) { k = i; } }
   return k;
}

/// Runs one point with source pulse s and step dt to the times in `ts` (multiples
/// of dt); calls `at(tp, t)` at each.
static void Run(const TPGridSpec &g, TPScheme sch, const TPParameters &P, const Pulse &s,
                real_t dt, const std::vector<real_t> &ts,
                const std::function<void(const ThermalPressurization &, real_t)> &at)
{
   ThermalPressurization tp(1, g, sch);
   tp.SetParameters(P);
   long n = 0;
   for (real_t te : ts)
   {
      const long ne = std::lround(te / dt);
      MFEM_VERIFY(std::fabs(ne * dt - te) < 1e-9 * te, "Run: t = " << te << " is not a multiple of dt");
      for (; n < ne; n++)
      {
         const real_t tau = s.Q((n + 0.5) * dt), V = 1.0;
         tp.Advance(dt, &tau, &V);
      }
      at(tp, te);
   }
}

/// Largest error of T and p on the fault, relative to the largest |analytic|,
/// over every step to t_end.
static void FaultError(const TPGridSpec &g, TPScheme sch, const TPParameters &P, const Pulse &s,
                       real_t dt, real_t t_end, real_t &eT, real_t &ep)
{
   ThermalPressurization tp(1, g, sch);
   tp.SetParameters(P);
   const long ns = std::lround(t_end / dt);
   real_t dTm = 0.0, dpm = 0.0, Tm = 0.0, pm = 0.0;
   for (long n = 0; n < ns; n++)
   {
      const real_t tau = s.Q((n + 0.5) * dt), V = 1.0;
      tp.Advance(dt, &tau, &V);
      const real_t t = (n + 1) * dt;
      const real_t Ta = s.T(0.0, t, P), pa = s.p(0.0, t, P);
      dTm = std::max(dTm, std::fabs(tp.TemperatureChange(0) - Ta));
      dpm = std::max(dpm, std::fabs(tp.PorePressureChange(0) - pa));
      Tm = std::max(Tm, std::fabs(Ta));
      pm = std::max(pm, std::fabs(pa));
   }
   eT = dTm / Tm;
   ep = dpm / pm;
}

// ---------------------------------------------------------------- P1
static void TestGrid()
{
   std::cout << "\nP1: the grid\n";
   const TPGridSpec g{1.0e-3, 1.05, 0.4};
   ThermalPressurization tp(2, g);
   tp.SetParameters(TPV105());
   const std::vector<real_t> &z = tp.Z();
   const int M = tp.NumZ();
   real_t e = 0.0;
   for (int k = 0; k < M; k++) { e = std::max(e, std::fabs((z[k + 1] - z[k]) / (g.dz0 * std::pow(g.growth, k)) - 1.0)); }
   std::cout << "  dz0 1 mm, growth 1.05, z_max 0.4 m: " << M << " unknowns, z_M = " << z[M]
             << " m, z_{M-1} = " << z[M - 1] << " m\n";
   TEST_ASSERT(z[0] == 0.0 && e < 1e-12, "P1 z_0 = 0 and z_{k+1} - z_k = dz0 r^k");
   TEST_ASSERT(z[M - 1] < g.z_max && z[M] >= g.z_max, "P1 z_{M-1} < z_max <= z_M");
   real_t sv = 0.0;
   for (real_t v : tp.Volumes()) { sv += v; }
   const real_t zh = 0.5 * (z[M - 1] + z[M]);
   TEST_ASSERT(std::fabs(sv - zh) < 1e-14 * zh, "P1 the control volumes sum to z_{M-1/2}");
   const real_t sf = 0.5 * std::erf(zh / (std::sqrt(2.0) * 0.01));
   TEST_ASSERT(std::fabs(tp.SourceFraction(1) - sf) < 1e-14,
               "P1 the source weights sum to erf(z_{M-1/2} / (sqrt(2) h)) / 2 = " << sf);
   TEST_ASSERT(tp.Temperature(0) == 483.15 && tp.PorePressure(1) == 80.0e6,
               "P1 T = T_ini and p = p_ini after SetParameters");
}

// ---------------------------------------------------------------- P2
static void TestZeros()
{
   std::cout << "\nP2: no heating, and Lambda = 0\n";
   const TPGridSpec g{1.0e-3, 1.05, 0.4};
   ThermalPressurization tp(2, g);
   TPParameters P = TPV105();
   tp.SetParameters(0, P);
   P.Lambda = 0.0;
   tp.SetParameters(1, P);
   const real_t tau0[2] = {0.0, 0.0}, V0[2] = {0.0, 0.0};
   for (int n = 0; n < 1000; n++) { tp.Advance(7.2e-4, tau0, V0); }
   bool zero = true;
   for (int j = 0; j < 2; j++)
      for (int k = 0; k < tp.NumZ(); k++)
      { zero = zero && tp.TemperatureChange(j, k) == 0.0 && tp.PorePressureChange(j, k) == 0.0; }
   TEST_ASSERT(zero, "P2 tau V = 0 for 1000 steps: T and p exactly unchanged at every node");
   const real_t tau[2] = {30.0e6, 30.0e6}, V[2] = {2.0, 2.0};
   for (int n = 0; n < 1000; n++) { tp.Advance(7.2e-4, tau, V); }
   bool p0 = true;
   for (int k = 0; k < tp.NumZ(); k++) { p0 = p0 && tp.PorePressureChange(1, k) == 0.0; }
   std::cout << "  after 1000 heated steps: T - T_ini on the fault " << tp.TemperatureChange(0)
             << " K; p - p_ini " << tp.PorePressureChange(0) << " Pa (Lambda 0.1 MPa/K), "
             << tp.PorePressureChange(1) << " Pa (Lambda = 0)\n";
   TEST_ASSERT(p0 && tp.TemperatureChange(1) == tp.TemperatureChange(0),
               "P2 Lambda = 0: p exactly p_ini at every node; T the same as with Lambda > 0");
}

// ---------------------------------------------------------------- P3
static void TestHeatBalance()
{
   std::cout << "\nP3: heat balance\n";
   const TPParameters P = TPV105();
   for (TPScheme sch : {TPScheme::BackwardEuler, TPScheme::CrankNicolson})
   {
      const TPGridSpec g{1.0e-3, 1.05, 2.0};   // far boundary well outside the heated zone
      ThermalPressurization tp(1, g, sch);
      tp.SetParameters(P);
      const real_t tau = 30.0e6, V = 1.0, dt = 7.2e-4;
      const int ns = 20000;   // 14.4 s
      for (int n = 0; n < ns; n++) { tp.Advance(dt, &tau, &V); }
      const real_t in = tau * V * ns * dt * tp.SourceFraction(0);
      const real_t e = std::fabs(tp.HeatContent(0) - in) / in;
      std::cout << "  " << (sch == TPScheme::BackwardEuler ? "backward Euler" : "Crank-Nicolson")
                << ": heat " << tp.HeatContent(0) << " J/m^2, put in " << in
                << ", relative difference " << e << "\n";
      TEST_ASSERT(e < 1e-12, "P3 " << (sch == TPScheme::BackwardEuler ? "backward Euler" : "Crank-Nicolson")
                  << ": the heat content equals the heat put in to 1e-12");
   }
}

// ---------------------------------------------------------------- P4
static void TestFormulas()
{
   std::cout << "\nP4: the analytic formulas\n";
   TPParameters P = TPV105();
   const real_t q = 3.0e7, t = 2.0;
   P.h = 1.0e-9;
   const real_t rice = q / P.rho_c * std::sqrt(t / (M_PI * P.alpha_th));
   const real_t th = Theta(0.0, t, P.alpha_th, P, q);
   std::cout << "  h = 1 nm: theta(0, 2 s) = " << th << " K, Rice " << rice << " K\n";
   TEST_ASSERT(std::fabs(th / rice - 1.0) < 1e-6, "P4 h -> 0: the formula is Rice's sqrt(t) law to 1e-6");
   P = TPV105();
   P.alpha_hy = 1.0e-12;   // nearly undrained
   real_t e = 0.0;
   for (real_t z : {0.0, 0.005, 0.02})
   { e = std::max(e, std::fabs(Pi(z, t, P, q) - P.Lambda * Theta(z, t, P.alpha_th, P, q)) / (P.Lambda * Theta(0.0, t, P.alpha_th, P, q))); }
   TEST_ASSERT(e < 1e-5, "P4 alpha_hy -> 0: pi = Lambda theta (" << e << ")");
}

// ---------------------------------------------------------------- P5
static void TestAnalytic()
{
   std::cout << "\nP5: TPV105 parameters against the analytic solution (dz0 0.1 mm, growth 1.01, dt 0.1 ms, Crank-Nicolson)\n";
   const TPParameters P = TPV105();
   const TPGridSpec g{1.0e-4, 1.01, 1.0};
   const real_t dt = 1.0e-4;
   const std::vector<real_t> ts = {0.01, 0.1, 1.0, 5.0, 15.0};
   for (int off = 0; off < 2; off++)
   {
      const Pulse s{3.0e7, 0.0, off ? 2.0 : 1.0e30};
      std::vector<real_t> tts = ts;
      if (off) { tts = {1.0, 2.5, 5.0, 15.0}; }
      real_t worst = 0.0;
      Run(g, TPScheme::CrankNicolson, P, s, dt, tts, [&](const ThermalPressurization &tp, real_t t) {
         const real_t T0 = s.T(0.0, t, P), p0 = s.p(0.0, t, P);
         real_t eT = 0.0, ep = 0.0;
         for (real_t z : {0.0, P.h, 2.0 * P.h, 5.0 * P.h})
         {
            const int k = NearestNode(tp, z);
            const real_t zk = tp.Z()[k];
            eT = std::max(eT, std::fabs(tp.TemperatureChange(0, k) - s.T(zk, t, P)) / std::fabs(T0));
            ep = std::max(ep, std::fabs(tp.PorePressureChange(0, k) - s.p(zk, t, P)) / std::fabs(p0));
         }
         std::cout << "  t = " << std::setw(5) << t << " s: T - T_ini " << std::setw(12) << tp.TemperatureChange(0)
                   << " K (analytic " << std::setw(12) << T0 << "), p - p_ini " << std::setw(12)
                   << tp.PorePressureChange(0) << " Pa (" << std::setw(12) << p0 << "); worst over z = 0, h, 2h, 5h: T "
                   << eT << ", p " << ep << "\n";
         worst = std::max(worst, std::max(eT, ep));
      });
      TEST_ASSERT(worst < 1e-4, "P5 " << (off ? "source off at 2 s" : "constant source")
                  << ": T and p at z = 0, h, 2h, 5h within 1e-4 of the fault value (" << worst << ")");
   }
}

// ---------------------------------------------------------------- P6
static void TestConvergence()
{
   std::cout << "\nP6: convergence (error of T and p on the fault at t = 1 s)\n";
   const TPParameters P = TPV105();
   const Pulse s{3.0e7, 0.0, 1.0e30};
   auto err = [&](const TPGridSpec &g, TPScheme sch, real_t dt, real_t &eT, real_t &ep) {
      Run(g, sch, P, s, dt, {1.0}, [&](const ThermalPressurization &tp, real_t t) {
         eT = std::fabs(tp.TemperatureChange(0) - s.T(0.0, t, P)) / s.T(0.0, t, P);
         ep = std::fabs(tp.PorePressureChange(0) - s.p(0.0, t, P)) / s.p(0.0, t, P);
      });
   };
   // space: a coarse grid, step small enough that the space error dominates
   real_t e1T, e1p, e2T, e2p;
   err({2.0e-3, 1.10, 1.0}, TPScheme::CrankNicolson, 1.0e-4, e1T, e1p);
   err({1.0e-3, std::sqrt(1.10), 1.0}, TPScheme::CrankNicolson, 1.0e-4, e2T, e2p);
   std::cout << "  space: dz0 2 mm, growth 1.10: T " << e1T << ", p " << e1p << "; dz0 1 mm, growth "
             << std::sqrt(1.10) << ": T " << e2T << ", p " << e2p << "; ratios " << e1T / e2T << ", "
             << e1p / e2p << "\n";
   TEST_ASSERT(e1T / e2T > 3.0 && e1p / e2p > 3.0, "P6 space: halving the spacings divides the error by > 3");
   // time: the differences of three runs with dt, dt/2, dt/4 on one grid
   // (the space error cancels): their ratio is 2 for a first-order scheme
   // and 4 for a second-order one
   for (TPScheme sch : {TPScheme::BackwardEuler, TPScheme::CrankNicolson})
   {
      real_t u[3][2];
      for (int i = 0; i < 3; i++)
      {
         Run({1.0e-3, 1.05, 1.0}, sch, P, s, 1.0e-2 / (1 << i), {1.0},
             [&](const ThermalPressurization &tp, real_t) {
                u[i][0] = tp.TemperatureChange(0);
                u[i][1] = tp.PorePressureChange(0);
             });
      }
      const real_t rT = (u[0][0] - u[1][0]) / (u[1][0] - u[2][0]);
      const real_t rp = (u[0][1] - u[1][1]) / (u[1][1] - u[2][1]);
      const bool be = sch == TPScheme::BackwardEuler;
      std::cout << "  time, " << (be ? "backward Euler" : "Crank-Nicolson") << ", dt 10, 5, 2.5 ms: "
                << "differences of p " << u[0][1] - u[1][1] << ", " << u[1][1] - u[2][1]
                << " Pa, ratio " << rp << "; of T ratio " << rT << "\n";
      if (be)
      {
         TEST_ASSERT(rp > 1.8 && rp < 2.2 && rT > 1.8 && rT < 2.2,
                     "P6 time, backward Euler: first order (ratios of T and p near 2)");
      }
      else
      {
         TEST_ASSERT(rp > 3.5 && rp < 4.5 && rT > 3.5 && rT < 4.5,
                     "P6 time, Crank-Nicolson: second order (ratios of T and p near 4)");
      }
   }
}

// ---------------------------------------------------------------- P7
static void TestProduction()
{
   std::cout << "\nP7: production settings (growth 1.05, z_max 0.4 m), 0 < t <= 15 s\n";
   const TPParameters P = TPV105();
   real_t worst = 0.0;
   for (const TPGridSpec &g : {TPGridSpec{1.0e-3, 1.05, 0.4}, TPGridSpec{0.5e-3, 1.05, 0.4}})
   {
      const ThermalPressurization probe(1, g);
      std::cout << "  dz0 " << g.dz0 * 1e3 << " mm (h / " << P.h / g.dz0 << "): " << probe.NumZ()
                << " unknowns per fault point\n";
      for (const Pulse &s : {Pulse{3.0e7, 0.0, 1.0e30}, Pulse{6.0e7, 0.5, 1.5}})
      {
         for (real_t dt : {7.2e-4, 3.1e-4})
         {
            for (TPScheme sch : {TPScheme::BackwardEuler, TPScheme::CrankNicolson})
            {
               // the pulse edges on step boundaries
               const real_t dtr = (s.t2 < 1e29) ? 0.5 / std::ceil(0.5 / dt) : dt;
               real_t eT, ep;
               FaultError(g, sch, P, s, dtr, 15.0, eT, ep);
               std::cout << "    " << (s.t2 < 1e29 ? "pulse 60 MW/m^2, 0.5-1.5 s" : "constant 30 MW/m^2  ")
                         << ", dt " << std::setw(9) << dtr << ", "
                         << (sch == TPScheme::BackwardEuler ? "backward Euler" : "Crank-Nicolson")
                         << ": T " << eT << ", p " << ep << "\n";
               if (g.dz0 == 1.0e-3) { worst = std::max(worst, std::max(eT, ep)); }
            }
         }
      }
   }
   TEST_ASSERT(worst < 1e-3, "P7 dz0 = h/10: T and p on the fault within 1e-3 (" << worst << ")");
}

// ---------------------------------------------------------------- P8
static void TestIndependence()
{
   std::cout << "\nP8: independent points, state round trip\n";
   const TPGridSpec g{1.0e-3, 1.05, 0.4};
   TPParameters P[3] = {TPV105(), TPV105(), TPV105()};
   P[1].h = 0.02; P[1].Lambda = 0.05e6; P[1].alpha_hy = 3.0e-5;
   P[2].alpha_th = 2.0e-6; P[2].rho_c = 3.0e6; P[2].T_ini = 400.0; P[2].p_ini = 50.0e6;
   const real_t tau[3] = {20.0e6, 35.0e6, 10.0e6};
   ThermalPressurization all(3, g, TPScheme::CrankNicolson);
   for (int j = 0; j < 3; j++) { all.SetParameters(j, P[j]); }
   std::vector<ThermalPressurization> one;
   for (int j = 0; j < 3; j++) { one.emplace_back(1, g, TPScheme::CrankNicolson); one[j].SetParameters(P[j]); }
   auto V = [](int n, int j) { return std::max(0.0, std::sin(0.01 * n + j)) * 3.0; };
   for (int n = 0; n < 500; n++)
   {
      const real_t v[3] = {V(n, 0), V(n, 1), V(n, 2)};
      all.Advance(7.2e-4, tau, v);
      for (int j = 0; j < 3; j++) { one[j].Advance(7.2e-4, &tau[j], &v[j]); }
   }
   bool same = true;
   for (int j = 0; j < 3; j++)
      for (int k = 0; k < all.NumZ(); k++)
      {
         same = same && all.TemperatureChange(j, k) == one[j].TemperatureChange(0, k)
                && all.PorePressureChange(j, k) == one[j].PorePressureChange(0, k);
      }
   TEST_ASSERT(same, "P8 three points equal three one-point runs bit for bit");
   Vector s;
   all.GetState(s);
   ThermalPressurization re(3, g, TPScheme::CrankNicolson);
   for (int j = 0; j < 3; j++) { re.SetParameters(j, P[j]); }
   re.SetState(s);
   for (int n = 500; n < 800; n++)
   {
      const real_t v[3] = {V(n, 0), V(n, 1), V(n, 2)};
      all.Advance(7.2e-4, tau, v);
      re.Advance(7.2e-4, tau, v);
   }
   Vector a, b;
   all.GetState(a);
   re.GetState(b);
   bool eq = a.Size() == b.Size();
   for (int i = 0; eq && i < a.Size(); i++) { eq = a(i) == b(i); }
   TEST_ASSERT(eq && re.PorePressure(2) == all.PorePressure(2),
               "P8 GetState / SetState: the restarted run equals the continuous one bit for bit");
}

// ---------------------------------------------------------------- refusals
static int Refuse(int n)
{
   std::cout << "test_thermal_pressurization --refuse " << n << " (must abort)\n";
   const TPGridSpec g{1.0e-3, 1.05, 0.4};
   const real_t tau = 1.0e6, V = 1.0;
   if (n == 4) { ThermalPressurization tp(1, {1.0e-3, 0.9, 0.4}); return 0; }
   ThermalPressurization tp(1, g);
   if (n == 6) { tp.Advance(1.0e-3, &tau, &V); return 0; }
   TPParameters P = TPV105();
   if (n == 5) { P.alpha_th = 0.0; }
   tp.SetParameters(P);
   if (n == 1) { tp.Advance(0.0, &tau, &V); }
   if (n == 2) { const real_t neg = -1.0e6; tp.Advance(1.0e-3, &neg, &V); }
   if (n == 3) { const real_t nan = std::numeric_limits<real_t>::quiet_NaN(); tp.Advance(1.0e-3, &tau, &nan); }
   if (n == 7) { Vector s(5); tp.SetState(s); }
   return 0;   // not refused: the make target reports it
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   int refuse = 0;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--refuse") == 0 && i + 1 < argc) { refuse = std::atoi(argv[++i]); }
   }
   if (refuse) { return Refuse(refuse); }
   std::cout << std::setprecision(6)
             << "test_thermal_pressurization (2D in-plane design rev 3, phase D)\n";
   TestGrid();
   TestZeros();
   TestHeatBalance();
   TestFormulas();
   TestAnalytic();
   TestConvergence();
   TestProduction();
   TestIndependence();
   std::cout << "\ntest_thermal_pressurization: " << (num_tests - num_failed) << " / " << num_tests
             << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
   return num_failed ? 1 : 0;
}
